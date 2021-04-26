#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "emmintrin.h"
#include "immintrin.h"
#include "smmintrin.h"
#include <sys/time.h> /* gettimeofday */

#include "config.h"            /* autoconf header */
#include "configs/base_configs.h"
#include "configs/eth_configs.h"

#include "utils/eth_data_structures.h"
#include "utils/data_generation.h"
#include "utils/io.h"

#include "utils/base_utils.h"
#include "utils/math.h"
#include "utils/barrier.h"
#include "utils/memory.h"
#include "utils/lock.h" 
#include "utils/learned_sort_for_sort_merge.h"


#ifndef KeyType
#define KeyType RELATION_KEY_TYPE
#define PayloadType RELATION_PAYLOAD_TYPE
#define TaskType Task<RELATION_KEY_TYPE, RELATION_PAYLOAD_TYPE>
#define NUM_THREADS NUM_THREADS_FOR_EVALUATION
#endif

#define RUN_NUMS 10 
#define INLJ_MORSE_SIZE 0 //100000

#define PREFETCH_SLOPES_AND_INTERCEPTS_MAJOR_BCKTS_UNIQUE_KEYS
#define SINGLE_TUPLE_PER_BUCKET

#ifdef INLJ_WITH_LEARNED_INDEX
#include "rmi/all_rmis.h"
using namespace INLJ_RMI_NAMESPACE;
#endif

#ifdef INLJ_WITH_CSS_TREE_INDEX
#include "utils/CC_CSSTree.h"
#endif

using namespace std;
using namespace learned_sort_for_sort_merge;

/*typedef struct StateSIMDForHashINLJ StateSIMDForHashINLJ;
struct StateSIMDForHashINLJ {
  __m512i key;
  __m512i ht_off;
  __mmask8 m_have_tuple;
  char stage;
};
*/

volatile static char inlj_g_lock = 0, inlj_g_lock_morse = 0;
volatile static uint64_t inlj_total_num = 0, inlj_global_curse = 0, inlj_global_morse_size;

typedef void (*INLJBuildFun)(IndexedNestedLoopJoinBuild<KeyType, PayloadType> *build_input, Relation<KeyType, PayloadType> * rel_r_partition, Relation<KeyType, PayloadType> * tmp_r);
volatile static struct Fun {
  INLJBuildFun fun_ptr;
  char fun_name[16];
} inlj_pfun[4];
volatile static int inlj_pf_num = 0;

typedef uint64_t (*INLJProbeFun)(Relation<KeyType, PayloadType> * rel_r_partition, Relation<KeyType, PayloadType> * rel_s_partition, IndexedNestedLoopJoinBuild<KeyType, PayloadType> *build_output);
volatile static struct Fun1 {
  INLJProbeFun fun_ptr;
  char fun_name[16];
} inlj_pfun1[4];

#ifdef INLJ_WITH_HASH_INDEX
Hashtable<KeyType, PayloadType> * ht;
#endif


#ifdef INLJ_WITH_HASH_INDEX
void inlj_with_hash_build_rel_r_partition(IndexedNestedLoopJoinBuild<KeyType, PayloadType> *build_input, Relation<KeyType, PayloadType> * rel_r_partition, Relation<KeyType, PayloadType> * tmp_r)
{
    Hashtable<KeyType, PayloadType>* ht = ((IndexedNestedLoopJoinBuild<KeyType, PayloadType> *)build_input)->ht;  
    BucketBuffer<KeyType, PayloadType>** overflowbuf = ((IndexedNestedLoopJoinBuild<KeyType, PayloadType> *)build_input)->overflowbuf;

    uint64_t i;
#ifndef USE_MURMUR3_HASH
    const uint64_t hashmask = ht->hash_mask;
    const uint64_t skipbits = ht->skip_bits;
#endif
#ifdef PREFETCH_INLJ
    size_t prefetch_index = PREFETCH_DISTANCE;
#endif
    
    for(i=0; i < rel_r_partition->num_tuples; i++){
        Tuple<KeyType, PayloadType> * dest;
        Bucket<KeyType, PayloadType> * curr, * nxt;

#ifdef PREFETCH_INLJ
        if (prefetch_index < rel_r_partition->num_tuples) {
#ifndef USE_MURMUR3_HASH
            uint64_t idx_prefetch = HASH(rel_r_partition->tuples[prefetch_index++].key,
                                         hashmask, skipbits);
#else
            uint64_t idx_prefetch_hash = murmur_hash_32(rel_r_partition->tuples[prefetch_index++].key);
            uint64_t idx_prefetch = alt_mod(idx_prefetch_hash, ht->num_buckets);
#endif
			__builtin_prefetch(ht->buckets + idx_prefetch, 1, 1);
        }
#endif

#ifndef USE_MURMUR3_HASH
        uint64_t idx = HASH(rel_r_partition->tuples[i].key, hashmask, skipbits);
#else
        uint64_t idx_hash = murmur_hash_32(rel_r_partition->tuples[i].key);
        uint64_t idx = alt_mod(idx_hash, ht->num_buckets);
#endif        
        /* copy the tuple to appropriate hash bucket */
        /* if full, follow nxt pointer to find correct place */
        curr = ht->buckets+idx;
        lock(&curr->latch);

        nxt = curr->next;

        if(curr->count == BUCKET_SIZE) {
            if(!nxt || nxt->count == BUCKET_SIZE) {
                Bucket<KeyType, PayloadType> * b;
                /* b = (bucket_t*) calloc(1, sizeof(bucket_t)); */
                /* instead of calloc() everytime, we pre-allocate */
                get_new_bucket(&b, overflowbuf);
                curr->next = b;
                b->next    = nxt;
                b->count   = 1;
                dest       = b->tuples;
            }
            else {
                dest = nxt->tuples + nxt->count;
                nxt->count ++;
            }
        }
        else 
        {
            dest = curr->tuples + curr->count;
            curr->count ++;
        }

        *dest = rel_r_partition->tuples[i];
        unlock(&curr->latch);
    }
}
#endif

#ifdef INLJ_WITH_HASH_INDEX
uint64_t inlj_with_hash_probe_rel_s_partition(Relation<KeyType, PayloadType> * rel_r_partition, Relation<KeyType, PayloadType> * rel_s_partition, IndexedNestedLoopJoinBuild<KeyType, PayloadType> *build_output)
{
    Hashtable<KeyType, PayloadType>* ht = ((IndexedNestedLoopJoinBuild<KeyType, PayloadType> *)build_output)->ht;  
    
    uint64_t i, j;
    uint64_t matches;

#ifndef USE_MURMUR3_HASH
    const uint64_t hashmask = ht->hash_mask;
    const uint64_t skipbits = ht->skip_bits;
#endif
#ifdef PREFETCH_INLJ    
    size_t prefetch_index = PREFETCH_DISTANCE;
#endif
    
    matches = 0; /*int curr_buckts_num;
    for(i=0; i < ht->num_buckets; i++)
    {
        Bucket<KeyType, PayloadType> * b = ht->buckets+i;
        //if((i < 5) && b && (b->count > 0))
        //    printf("naive i %ld key %ld \n", i, b->tuples[0].key);
        curr_buckts_num = 0;
        do {
            b = b->next;
            //if((i < 5) && b && (b->count > 0))
            //    printf("naive i %ld key %ld \n", i, b->tuples[0].key);
            curr_buckts_num++;
        } while(b);
        if((curr_buckts_num > 2) && (i < 100))
            printf("naive i %ld curr_buckets_num %d nbuckets %ld \n", i, curr_buckts_num, ht->num_buckets);
    }*/

    for (i = 0; i < rel_s_partition->num_tuples; i++)
    {
#ifdef PREFETCH_INLJ        
        if (prefetch_index < rel_s_partition->num_tuples) {
#ifndef USE_MURMUR3_HASH
			uint64_t idx_prefetch = HASH(rel_s_partition->tuples[prefetch_index++].key,
                                        hashmask, skipbits);
#else
            uint64_t idx_prefetch_hash = murmur_hash_32(rel_s_partition->tuples[prefetch_index++].key);
            uint64_t idx_prefetch = alt_mod(idx_prefetch_hash, ht->num_buckets);
#endif
			__builtin_prefetch(ht->buckets + idx_prefetch, 0, 1);
        }
#endif

#ifndef USE_MURMUR3_HASH        
        uint64_t idx = HASH(rel_s_partition->tuples[i].key, hashmask, skipbits);        
#else
        uint64_t idx_hash = murmur_hash_32(rel_s_partition->tuples[i].key);
        uint64_t idx = alt_mod(idx_hash, ht->num_buckets);
#endif
        Bucket<KeyType, PayloadType> * b = ht->buckets+idx;

        do {
        #ifdef SINGLE_TUPLE_PER_BUCKET    
            if(rel_s_partition->tuples[i].key == b->tuples[0].key){
                    matches ++;
            }
        #else
            for(j = 0; j < b->count; j++) {
                if(rel_s_partition->tuples[i].key == b->tuples[j].key){
                    matches ++;
                }
            }
        #endif

            b = b->next;/* follow overflow pointer */
        } while(b);

    }

    return matches;
}
#endif

#ifdef INLJ_WITH_LEARNED_INDEX
uint64_t inlj_with_rmi_probe_rel_s_partition(Relation<KeyType, PayloadType> * rel_r_partition, Relation<KeyType, PayloadType> * rel_s_partition, IndexedNestedLoopJoinBuild<KeyType, PayloadType> *build_output)
{
    uint64_t i;
    uint64_t matches = 0; 
    size_t err;
    uint64_t rmi_guess;
    uint64_t bound_start, bound_end;
    int n; uint64_t lower;
    KeyType * sorted_relation_r_keys_only = build_output->sorted_relation_r_keys_only;
    uint64_t original_relR_num_tuples = build_output->original_relR->num_tuples;

    for (i = 0; i < rel_s_partition->num_tuples; i++)
    {        
        rmi_guess = INLJ_RMI_NAMESPACE::lookup(rel_s_partition->tuples[i].key, &err);
        bound_start = rmi_guess - err; 
        bound_start = (bound_start < 0)? 0 : bound_start;
        bound_end = rmi_guess + err;
        bound_end = (bound_end > original_relR_num_tuples - 1)? original_relR_num_tuples - 1 : bound_end;

        n = bound_end - bound_start + 1; // `end` is inclusive.
        lower = bound_start;

        // Function adapted from https://github.com/gvinciguerra/rmi_pgm/blob/357acf668c22f927660d6ed11a15408f722ea348/main.cpp#L29.
        // Authored by Giorgio Vinciguerra.
        while (const int half = n / 2) {
            const uint64_t middle = lower + half;
            // Prefetch next two middles.
            __builtin_prefetch(&(sorted_relation_r_keys_only[lower + half / 2]), 0, 0);
            __builtin_prefetch(&(sorted_relation_r_keys_only[middle + half / 2]), 0, 0);
            lower = (sorted_relation_r_keys_only[middle] <= rel_s_partition->tuples[i].key) ? middle : lower;
            n -= half;
        }

        // Scroll back to the first occurrence.
        while (lower > 0 && sorted_relation_r_keys_only[lower - 1] == rel_s_partition->tuples[i].key) --lower;

        if (sorted_relation_r_keys_only[lower] == rel_s_partition->tuples[i].key) 
        {
            // Sum over all values with that key.
            for (unsigned int k = lower; sorted_relation_r_keys_only[k] == rel_s_partition->tuples[i].key && k < original_relR_num_tuples; ++k) {
                matches ++;
            }

        }
    }

    return matches;
}
#endif


#ifdef INLJ_WITH_CSS_TREE_INDEX
uint64_t inlj_with_csstree_probe_rel_s_partition(Relation<KeyType, PayloadType> * rel_r_partition, Relation<KeyType, PayloadType> * rel_s_partition, IndexedNestedLoopJoinBuild<KeyType, PayloadType> *build_output)
{
    uint64_t i, j;
    uint64_t matches = 0; 
    uint64_t curIndex=0;
	KeyType keyForSearch;

    CC_CSSTree<KeyType, PayloadType> *tree = build_output->tree;
    KeyType * sorted_relation_r_keys_only = build_output->sorted_relation_r_keys_only;
    uint64_t original_relR_num_tuples = build_output->original_relR->num_tuples;

    for (i = 0; i < rel_s_partition->num_tuples; i++)
    {
		keyForSearch = rel_s_partition->tuples[i].key;
		curIndex = tree->search(keyForSearch);
        
        for(j = curIndex-1; j > 0; j--)
		{
            if (sorted_relation_r_keys_only[j] == keyForSearch)
                matches++;
            else 
                if(sorted_relation_r_keys_only[j] < keyForSearch)
			    	break;
        }
        
        for(j = curIndex; j < original_relR_num_tuples; j++)
        {
            if (sorted_relation_r_keys_only[j] == keyForSearch)
                matches++;
            else 
                if(sorted_relation_r_keys_only[j] > keyForSearch)
			    	break;
        }
    }

    return matches;
}
#endif

void * inlj_join_thread(void * param)
{
    IndexedNestedLoopJoinThread<KeyType, PayloadType, TaskType> * args   = (IndexedNestedLoopJoinThread<KeyType, PayloadType, TaskType> *) param;
    int rv;   int deltaT = 0; struct timeval t1, t2;
#ifdef INLJ_WITH_LEARNED_INDEX        
    /*for (int rp = 0; rp < RUN_NUMS; ++rp) 
    {
        if(args->tid == 0)
            init_models_training_data_and_sample_counts<KeyType, PayloadType>(args->training_data, args->p.arch, 
                    args->sample_count, args->sample_count_R, args->sample_count_S, NUM_THREADS);

        BARRIER_ARRIVE(args->barrier, rv);
        if(args->tid == 0){
            gettimeofday(&t1, NULL);
        }

        sample_and_train_models_threaded(args);

        BARRIER_ARRIVE(args->barrier, rv);

        if(args->tid == 0){
            gettimeofday(&t2, NULL);

            deltaT = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
            printf("---- Sampling and training models time (ms) = %10.4lf\n",  deltaT * 1.0 / 1000);

#ifndef RUN_LEARNED_TECHNIQUES_WITH_FIRST_LEVEL_ONLY
            if(rp == RUN_NUMS - 1)
            {   
                for (unsigned int j = 0; j < args->rmi->hp.arch[1]; ++j) 
                {
                    args->slopes->push_back(args->rmi->models[1][j].slope);
                    args->intercepts->push_back(args->rmi->models[1][j].intercept);
                }
            } 
#endif            
        }        
    }*/
#endif
    
#ifdef INLJ_WITH_HASH_INDEX
    BucketBuffer<KeyType, PayloadType> * overflowbuf; // allocate overflow buffer for each thread

#if INPUT_HASH_TABLE_SIZE       
    uint32_t nbuckets = INPUT_HASH_TABLE_SIZE;
#else
    uint32_t nbuckets = (args->original_relR->num_tuples / BUCKET_SIZE / NUM_THREADS);
#endif

#endif    
    if (args->tid == 0) {

#ifdef INLJ_WITH_HASH_INDEX        
        strcpy(inlj_pfun[0].fun_name, "Hashing");

        inlj_pfun[0].fun_ptr = inlj_with_hash_build_rel_r_partition;

        strcpy(inlj_pfun1[0].fun_name, "Hashing");

        inlj_pfun1[0].fun_ptr = inlj_with_hash_probe_rel_s_partition;

        inlj_pf_num = 1;
#endif

#ifdef INLJ_WITH_LEARNED_INDEX        
        //strcpy(inlj_pfun[0].fun_name, "Learned");
        //inlj_pfun[0].fun_ptr = inlj_with_rmi_build_rel_r_partition;

        strcpy(inlj_pfun1[0].fun_name, "Learned");
        
        inlj_pfun1[0].fun_ptr = inlj_with_rmi_probe_rel_s_partition;

        inlj_pf_num = 1;
#endif

#ifdef INLJ_WITH_CSS_TREE_INDEX        
        //strcpy(inlj_pfun[0].fun_name, "CSSTree");
        //inlj_pfun[0].fun_ptr = inlj_with_csstree_build_rel_r_partition;

        strcpy(inlj_pfun1[0].fun_name, "CSSTree");
        
        inlj_pfun1[0].fun_ptr = inlj_with_csstree_probe_rel_s_partition;

        inlj_pf_num = 1;
#endif
    }
    BARRIER_ARRIVE(args->barrier, rv);
    
    
    IndexedNestedLoopJoinBuild<KeyType, PayloadType> build_data; 
    build_data.original_relR = args->original_relR;
    build_data.original_relS = args->original_relS;
    build_data.sorted_relation_r_keys_only = args->sorted_relation_r_keys_only;
#ifdef INLJ_WITH_CSS_TREE_INDEX        
    build_data.tree = args->tree;
#endif

    for (int fid = 0; fid < inlj_pf_num; ++fid) 
    {
        for (int rp = 0; rp < RUN_NUMS; ++rp) 
        {

        #ifdef  INLJ_WITH_HASH_INDEX          
            init_bucket_buffer(&overflowbuf);
            if(args->tid == 0)
                allocate_hashtable(&ht, nbuckets);
                //allocate_hashtable(&args->ht, nbuckets);
        #endif
            BARRIER_ARRIVE(args->barrier, rv);

        #ifdef INLJ_WITH_HASH_INDEX    
            args->ht = ht;
            build_data.ht = ht;
            build_data.overflowbuf = &overflowbuf;
        #endif

        #ifdef INLJ_WITH_LEARNED_INDEX   
            /*build_data.rmi = args->rmi;
            build_data.slopes = args->slopes;
            build_data.intercepts = args->intercepts;*/
        #endif
        #ifdef PERF_COUNTERS
            if(args->tid == 0){
                //TODO: performance counters to be implemented
            }
        #endif

            BARRIER_ARRIVE(args->barrier, rv);

#ifdef INLJ_WITH_HASH_INDEX 
            if(args->tid == 0){
                gettimeofday(&args->start_time, NULL);
            #ifndef DEVELOPMENT_MODE
                //args->e_start_to_partition.startCounters();
            #endif
            }  

        #if INLJ_MORSE_SIZE
            morse_driven(param, inlj_pfun[fid].fun_ptr, &overflowbuf);
        #else
            inlj_pfun[fid].fun_ptr(&build_data, &args->relR, NULL);
        #endif

            BARRIER_ARRIVE(args->barrier, rv);

        #ifdef PERF_COUNTERS
            if(args->tid == 0)
            {
                //TODO: performance counters to be implemented
            }
            BARRIER_ARRIVE(args->barrier, rv);
        #endif

            if(args->tid == 0){
                gettimeofday(&args->partition_end_time, NULL);

            #ifndef DEVELOPMENT_MODE
                //args->e_start_to_partition.stopCounters();
                //args->e_partition_to_end.startCounters();
            #endif
                deltaT = (args->partition_end_time.tv_sec - args->start_time.tv_sec) * 1000000 + args->partition_end_time.tv_usec - args->start_time.tv_usec;
                printf("---- %5s Build costs time (ms) = %10.4lf\n", inlj_pfun[fid].fun_name, deltaT * 1.0 / 1000);
                inlj_total_num = 0;
                inlj_global_curse = 0;
            }
#endif        
            if(!((fid == (inlj_pf_num - 1)) && (rp == (RUN_NUMS - 1)))){
            
            #ifdef INLJ_WITH_HASH_INDEX  
                if(args->tid == 0)
                    destroy_hashtable(args->ht);

                free_bucket_buffer(overflowbuf);
            #endif

                BARRIER_ARRIVE(args->barrier, rv);
            } 
        }
    }

    BARRIER_ARRIVE(args->barrier, rv);
    
    //Probe phase
    for (int fid = 0; fid < inlj_pf_num; ++fid) 
    {
        for (int rp = 0; rp < RUN_NUMS; ++rp) 
        {
            BARRIER_ARRIVE(args->barrier, rv);

            if(args->tid == 0){
                gettimeofday(&args->partition_end_time, NULL);
            }

            #if INLJ_MORSE_SIZE
                //TODO: to be done
            #else
                args->num_results = inlj_pfun1[fid].fun_ptr(NULL, &args->relS, &build_data);
            #endif
            
            BARRIER_ARRIVE(args->barrier, rv);
            // probe phase finished, thread-0 checkpoints the time
            if(args->tid == 0){
                gettimeofday(&args->end_time, NULL);

                deltaT = (args->end_time.tv_sec - args->partition_end_time.tv_sec) * 1000000 + args->end_time.tv_usec - args->partition_end_time.tv_usec;
                printf("---- %5s Probe costs time (ms) = %10.4lf\n", inlj_pfun1[fid].fun_name, deltaT * 1.0 / 1000);
            }
        }
    }

    return 0;
}

void initialize_inlj_join_thread_args(Relation<KeyType, PayloadType> * rel_r, 
                                 Relation<KeyType, PayloadType> * rel_s,
                        #ifdef INLJ_WITH_HASH_INDEX
                                 Hashtable<KeyType, PayloadType> * ht, 
                        #endif
                                 KeyType * sorted_relation_r_keys_only,
                        #ifdef INLJ_WITH_LEARNED_INDEX
                                 /*learned_sort_for_sort_merge::RMI<KeyType, PayloadType> * rmi,
                                 learned_sort_for_sort_merge::RMI<KeyType, PayloadType>::Params p,
                                 unsigned int SAMPLE_SZ_R, unsigned int SAMPLE_SZ_S,
                                 Tuple<KeyType, PayloadType> * tmp_training_sample_in,
                                 Tuple<KeyType, PayloadType> * sorted_training_sample_in,
                                 Tuple<KeyType, PayloadType> * r_tmp_training_sample_in,
                                 Tuple<KeyType, PayloadType> * r_sorted_training_sample_in,
                                 Tuple<KeyType, PayloadType> * s_tmp_training_sample_in,
                                 Tuple<KeyType, PayloadType> * s_sorted_training_sample_in,
                                 vector<vector<vector<training_point<KeyType, PayloadType>>>> * training_data,
                                 uint32_t * sample_count, uint32_t * sample_count_R, uint32_t * sample_count_S,
                                 vector<double>* slopes, vector<double>* intercepts,*/
                        #endif
                        #ifdef INLJ_WITH_CSS_TREE_INDEX
                                CC_CSSTree<KeyType, PayloadType> *tree,
                        #endif                        
                                 pthread_barrier_t* barrier_ptr,
                                 Result * joinresult,
                                 IndexedNestedLoopJoinThread<KeyType, PayloadType, TaskType> * args){
    int i;
    uint64_t numR, numS, numRthr, numSthr; /* total and per thread num */

#ifdef INLJ_WITH_LEARNED_INDEX
    //unsigned int SAMPLE_SZ_Rthr, SAMPLE_SZ_Sthr;
#endif
    numR = rel_r->num_tuples;
    numS = rel_s->num_tuples;
    numRthr = numR / NUM_THREADS;
    numSthr = numS / NUM_THREADS;
#ifdef INLJ_WITH_LEARNED_INDEX
    /*SAMPLE_SZ_Rthr = SAMPLE_SZ_R / NUM_THREADS;
    SAMPLE_SZ_Sthr = SAMPLE_SZ_S / NUM_THREADS;*/
#endif

    for(i = 0; i < NUM_THREADS; i++)
    {
        (*(args + i)).tid = i;
    #ifdef INLJ_WITH_HASH_INDEX
        (*(args + i)).ht = ht; 
    #endif
        /* assing part of the relR for next thread */
        (*(args + i)).relR.num_tuples = (i == (NUM_THREADS-1)) ? numR : numRthr;
        (*(args + i)).relR.tuples = rel_r->tuples + numRthr * i;
        numR -= numRthr;

        /* assing part of the relS for next thread */
        (*(args + i)).relS.num_tuples = (i == (NUM_THREADS-1)) ? numS : numSthr;
        (*(args + i)).relS.tuples = rel_s->tuples + numSthr * i;
        numS -= numSthr;

        (*(args + i)).original_relR = rel_r;
        (*(args + i)).original_relS = rel_s;

        (*(args + i)).sorted_relation_r_keys_only = sorted_relation_r_keys_only;

    #ifdef INLJ_WITH_LEARNED_INDEX
        /**** start stuff for learning RMI models ****/
        /*(*(args + i)).rmi = rmi;
        (*(args + i)).p = p;
        (*(args + i)).tmp_training_sample_in = tmp_training_sample_in;
        (*(args + i)).sorted_training_sample_in = sorted_training_sample_in;
        (*(args + i)).r_tmp_training_sample_in = r_tmp_training_sample_in;
        (*(args + i)).r_sorted_training_sample_in = r_sorted_training_sample_in;
        (*(args + i)).s_tmp_training_sample_in = s_tmp_training_sample_in;
        (*(args + i)).s_sorted_training_sample_in = s_sorted_training_sample_in;
        (*(args + i)).training_data = training_data;
        (*(args + i)).tmp_training_sample_R_offset = SAMPLE_SZ_Rthr * i;
        (*(args + i)).tmp_training_sample_S_offset = SAMPLE_SZ_Sthr * i;
        (*(args + i)).tmp_training_sample_offset = (SAMPLE_SZ_Rthr + SAMPLE_SZ_Sthr) * i;
        (*(args + i)).sample_count = sample_count;
        (*(args + i)).sample_count_R = sample_count_R;
        (*(args + i)).sample_count_S = sample_count_S;
        (*(args + i)).slopes = slopes;
        (*(args + i)).intercepts = intercepts;*/
        /**** end stuff for learning RMI models ****/
    #endif
    #ifdef INLJ_WITH_CSS_TREE_INDEX
        (*(args + i)).tree = tree;
    #endif
        (*(args + i)).barrier = barrier_ptr;
        (*(args + i)).threadresult  = &(joinresult->resultlist[i]);
    }
}


int main(int argc, char **argv) 
{
    Relation<KeyType, PayloadType> rel_r;
    Relation<KeyType, PayloadType> rel_s;
    
    int64_t result = 0;
    uint64_t curr_num_tuples_r = RELATION_R_NUM_TUPLES;
    uint64_t curr_num_tuples_s = RELATION_S_NUM_TUPLES; 

#ifdef LOAD_RELATIONS_FOR_EVALUATION
    // loading pre-built datasets
    string curr_rel_r_path = RELATION_R_PATH;
    string curr_rel_s_path = RELATION_S_PATH;

    load_relation<KeyType, PayloadType>(&rel_r, curr_rel_r_path.c_str(), curr_num_tuples_r);
    load_relation<KeyType, PayloadType>(&rel_s, curr_rel_s_path.c_str(), curr_num_tuples_s);    
#else
    // creating new datasets on-the-flay 
    result = create_eth_workload_relation_pk<KeyType, PayloadType>(&rel_r, curr_num_tuples_r, 0);
    //ASSERT_EQ(result, 0);
    #ifdef PERSIST_RELATIONS_FOR_EVALUATION
    write_relation<KeyType, PayloadType>(&rel_r, rel_r_path.c_str());
    #endif
    
    result = create_eth_workload_relation_pk<KeyType, PayloadType>(&rel_s, curr_num_tuples_s, 0);
    //ASSERT_EQ(result, 0);
    #ifdef PERSIST_RELATIONS_FOR_EVALUATION
    write_relation<KeyType, PayloadType>(&rel_s, rel_s_path.c_str());
    #endif
#endif

    int i, rv;
    pthread_barrier_t barrier;
    Result * joinresult;
    pthread_t tid[NUM_THREADS];
    pthread_attr_t attr;
    cpu_set_t set;
    
    joinresult = (Result *) malloc(sizeof(Result));

    IndexedNestedLoopJoinThread<KeyType, PayloadType, TaskType> args[NUM_THREADS];
    IndexedNestedLoopJoinThread<KeyType, PayloadType, TaskType> * args_ptr = args;

    rv = pthread_barrier_init(&barrier, NULL, NUM_THREADS);
    if(rv != 0){
        printf("[ERROR] Couldn't create the barrier\n");
        exit(EXIT_FAILURE);
    }

    pthread_attr_init(&attr);

#ifdef INLJ_WITH_HASH_INDEX
#if INPUT_HASH_TABLE_SIZE       
    uint32_t nbuckets = INPUT_HASH_TABLE_SIZE;
#else
    uint32_t nbuckets = (rel_r.num_tuples / BUCKET_SIZE / NUM_THREADS);

#endif        
    allocate_hashtable(&ht, nbuckets);
#endif

    KeyType * sorted_relation_r_keys_only = (KeyType *) alloc_aligned(rel_r.num_tuples  * sizeof(KeyType));

    for(int j = 0; j < rel_r.num_tuples; j++)
        sorted_relation_r_keys_only[j] = rel_r.tuples[j].key;
    
    std::sort((KeyType *)(sorted_relation_r_keys_only), (KeyType *)(sorted_relation_r_keys_only) + rel_r.num_tuples);

#ifdef INLJ_WITH_LEARNED_INDEX
  
    std::cout << "RMI status: " << INLJ_RMI_NAMESPACE::load(INLJ_RMI_DATA_PATH) << std::endl;    
/*
    //////////////////////////////////////////////////////////////////////////////
    // start stuff for sampling and building RMI models for both relations R and S
    //////////////////////////////////////////////////////////////////////////////
    Tuple<KeyType, PayloadType>* r_tmp_training_sample_in;
    Tuple<KeyType, PayloadType>* r_sorted_training_sample_in;
    Tuple<KeyType, PayloadType>* s_tmp_training_sample_in;
    Tuple<KeyType, PayloadType>* s_sorted_training_sample_in;
    unsigned int SAMPLE_SZ_R, SAMPLE_SZ_S;
    
    // Sampling and building RMI models for relations R and S together
    typename learned_sort_for_sort_merge::RMI<KeyType, PayloadType>::Params rmi_params;
    learned_sort_for_sort_merge::validate_params<KeyType, PayloadType>(rmi_params, rel_r.num_tuples);
    learned_sort_for_sort_merge::validate_params<KeyType, PayloadType>(rmi_params, rel_s.num_tuples);
    SAMPLE_SZ_R = std::min<unsigned int>(
        rel_r.num_tuples, std::max<unsigned int>(rmi_params.sampling_rate * rel_r.num_tuples,
                                        learned_sort_for_sort_merge::RMI<KeyType, PayloadType>::Params::MIN_SORTING_SIZE)) + 1;
    r_tmp_training_sample_in = (Tuple<KeyType, PayloadType>*) alloc_aligned(SAMPLE_SZ_R * sizeof(Tuple<KeyType, PayloadType>));
    #ifdef USE_AVXSORT_AS_STD_SORT
    r_sorted_training_sample_in = (Tuple<KeyType, PayloadType>*) alloc_aligned(SAMPLE_SZ_R * sizeof(Tuple<KeyType, PayloadType>));
    #endif
    SAMPLE_SZ_S = std::min<unsigned int>(
        rel_s.num_tuples, std::max<unsigned int>(rmi_params.sampling_rate * rel_s.num_tuples,
                                        learned_sort_for_sort_merge::RMI<KeyType, PayloadType>::Params::MIN_SORTING_SIZE)) + 1;
    s_tmp_training_sample_in = (Tuple<KeyType, PayloadType>*) alloc_aligned(SAMPLE_SZ_S * sizeof(Tuple<KeyType, PayloadType>));
    #ifdef USE_AVXSORT_AS_STD_SORT
    s_sorted_training_sample_in = (Tuple<KeyType, PayloadType>*) alloc_aligned(SAMPLE_SZ_S * sizeof(Tuple<KeyType, PayloadType>));
    #endif
    
    Tuple<KeyType, PayloadType>* tmp_training_sample_in;
    Tuple<KeyType, PayloadType>* sorted_training_sample_in;
    tmp_training_sample_in = (Tuple<KeyType, PayloadType>*) alloc_aligned((SAMPLE_SZ_R + SAMPLE_SZ_S) * sizeof(Tuple<KeyType, PayloadType>));
    #ifdef USE_AVXSORT_AS_STD_SORT
    sorted_training_sample_in = (Tuple<KeyType, PayloadType>*) alloc_aligned((SAMPLE_SZ_R + SAMPLE_SZ_S) * sizeof(Tuple<KeyType, PayloadType>));
    #endif

    RMI<KeyType, PayloadType> rmi(rmi_params, tmp_training_sample_in, sorted_training_sample_in,
                                              r_tmp_training_sample_in, r_sorted_training_sample_in,
                                              s_tmp_training_sample_in, s_sorted_training_sample_in);
    vector<vector<vector<training_point<KeyType, PayloadType>>>> training_data(rmi_params.arch.size());
    for (unsigned int layer_idx = 0; layer_idx < rmi_params.arch.size(); ++layer_idx) {
        training_data[layer_idx].resize(rmi_params.arch[layer_idx]);
    }

    uint32_t * sample_count = (uint32_t *) calloc(NUM_THREADS, sizeof(uint32_t)); 
    uint32_t * sample_count_R = (uint32_t *) calloc(NUM_THREADS, sizeof(uint32_t)); 
    uint32_t * sample_count_S = (uint32_t *) calloc(NUM_THREADS, sizeof(uint32_t));

    vector<double>* slopes = new vector<double>;                 
    vector<double>* intercepts = new vector<double>;
    //////////////////////////////////////////////////////////////////////////////
    // End stuff for sampling and building RMI models for both relations R and S
    //////////////////////////////////////////////////////////////////////////////
*/
#endif

#ifdef INLJ_WITH_CSS_TREE_INDEX
    for(int j = 0; j < rel_r.num_tuples; j++)
        rel_r.tuples[j].key = sorted_relation_r_keys_only[j];

	CC_CSSTree<KeyType, PayloadType> *tree=new CC_CSSTree<KeyType, PayloadType>(&rel_r, rel_r.num_tuples, INLJ_CSS_TREE_FANOUT);

#endif
    initialize_inlj_join_thread_args(&rel_r, &rel_s, 
                                #ifdef INLJ_WITH_HASH_INDEX
                                    ht, 
                                #endif
                                    sorted_relation_r_keys_only,
                                #ifdef INLJ_WITH_LEARNED_INDEX                                    
                                    /*&rmi, rmi_params,
                                    SAMPLE_SZ_R, SAMPLE_SZ_S,
                                    tmp_training_sample_in, sorted_training_sample_in, r_tmp_training_sample_in,
                                    r_sorted_training_sample_in, s_tmp_training_sample_in, s_sorted_training_sample_in,
                                    &training_data, sample_count, sample_count_R, sample_count_S,
                                    slopes, intercepts,*/
                                #endif
                                #ifdef INLJ_WITH_CSS_TREE_INDEX
                                    tree,
                                #endif
                                    &barrier, joinresult, args_ptr);

    inlj_global_curse = 0;
    if(NUM_THREADS==1){
        inlj_global_morse_size= rel_r.num_tuples;
    }else{
        inlj_global_morse_size = INLJ_MORSE_SIZE;
    }

    for(i = 0; i < NUM_THREADS; i++)
    {
        #ifdef DEVELOPMENT_MODE
        int cpu_idx = get_cpu_id_develop(i);
        #else
        int cpu_idx = get_cpu_id(i);
        #endif

        CPU_ZERO(&set);
        CPU_SET(cpu_idx, &set);
        pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);

        rv = pthread_create(&tid[i], &attr, inlj_join_thread, (void*)&args[i]);
        if (rv){
            printf("[ERROR] return code from pthread_create() is %d\n", rv);
            exit(-1);
        }
    }

    // wait for threads to finish
    for(i = 0; i < NUM_THREADS; i++){
        pthread_join(tid[i], NULL);
        result += args[i].num_results;
    }        

    printf("join results: %ld \n", result);

  return 0;
}






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
#define JoinThreadType ETHRadixJoinThread<KeyType, PayloadType, TaskType>
#define PartitionType ETHPartition<KeyType, PayloadType, TaskType, ETHRadixJoinThread<KeyType, PayloadType, TaskType>>
#define NUM_THREADS NUM_THREADS_FOR_EVALUATION
#endif

#define RUN_NUMS 1 // Fix it to be 1 for now 

#define PREFETCH_SLOPES_AND_INTERCEPTS_MAJOR_BCKTS_UNIQUE_KEYS

using namespace std;
using namespace learned_sort_for_sort_merge;

    int32_t ** histR;
    int32_t ** histS;
    Tuple<KeyType, PayloadType> * tmpRelR;
    Tuple<KeyType, PayloadType> * tmpRelS;
    TaskQueue<KeyType, PayloadType, TaskType> ** part_queue_ptr;
    TaskQueue<KeyType, PayloadType, TaskType> ** join_queue_ptr;
#if SKEW_HANDLING
    TaskQueue<KeyType, PayloadType, TaskType> * skew_queue;
#endif
    int small_padding_tupples = 0; 
    int padding_tuples =  small_padding_tupples * (FANOUT_PASS1 + 1);
    int relation_padding = padding_tuples * FANOUT_PASS1 * sizeof(Tuple<KeyType, PayloadType>); 
    int l1_cache_tuples = L1_CACHE_SIZE/sizeof(Tuple<KeyType, PayloadType>);

typedef void (*PJPartitionFun)(PartitionType * part);
volatile static struct PartitionFun {
  PJPartitionFun fun_ptr;
  char fun_name[16];
} pj_partition_pfun[4];
volatile static int pj_partition_pf_num = 0;

typedef uint32_t (*PJBuildFun)(ETHBucketChainingBuild *build_output, Relation<KeyType, PayloadType> * rel_r_partition, Relation<KeyType, PayloadType> * tmp_r);
volatile static struct BuildFun {
  PJBuildFun fun_ptr;
  char fun_name[16];
} pj_build_pfun[4];
volatile static int pj_build_pf_num = 0;

typedef uint64_t (*PJProbeFun)(Relation<KeyType, PayloadType> * rel_r_partition, uint32_t numR_from_build, Relation<KeyType, PayloadType> * rel_s_partition, ETHBucketChainingBuild *build_output);
volatile static struct ProbeFun {
  PJProbeFun fun_ptr;
  char fun_name[16];
} pj_probe_pfun[4];
volatile static int pj_probe_pf_num = 0;

void pj_partition_rel_segment_pass1(PartitionType * part) {
    const Tuple<KeyType, PayloadType> * restrict rel    = part->rel;
    int32_t **               hist   = part->hist;
    int64_t *       restrict output = part->output;

    const uint32_t my_tid     = part->thrargs->my_tid;
    const uint32_t nthreads   = part->thrargs->nthreads;
    const uint32_t num_tuples = part->num_tuples;

    const int32_t  R       = part->R;
    const int32_t  D       = part->D;
    const uint32_t fanOut  = 1 << D;
    const uint32_t MASK    = (fanOut - 1) << R;
    const uint32_t padding = part->padding;

    int64_t sum = 0;
    int64_t i, j;
    int rv;

    // compute local histogram for the assigned region of rel
    // compute histogram 
    int32_t * my_hist = hist[my_tid];
#ifndef USE_VECTORIZED_MURMUR3_HASH_FOR_RADIX_JOIN
    for(i = 0; i < num_tuples; i++) 
    {
        //printf("inside partition segment pass 1 %d key %d\n", my_tid, rel[i].key);
    #ifndef USE_MURMUR3_HASH_FOR_RADIX_JOIN
        uint32_t idx = HASH_BIT_MODULO(rel[i].key, MASK, R);
        //printf("inside partition segment pass 1 %d key %d idx %ld \n", my_tid, rel[i].key, idx);
    #else
        uint32_t idx_hash = murmur_hash_32(rel[i].key);
        uint32_t idx = HASH_BIT_MODULO(idx_hash, MASK, R);
        //printf("inside partition segment pass 1 %d key %d idx_hash %ld idx %ld \n", my_tid, rel[i].key, idx_hash, idx);
    #endif

        my_hist[idx] ++;
    }

#else
    printf("USE_VECTORIZED_MURMUR3_HASH_FOR_RADIX_JOIN needs to be double checked for correctness \n");
    /*
    int64_t num_tuples_batches = num_tuples / 8;
    uint32_t num_tuples_reminders = num_tuples % 8;
    i = 0;

    for(j = 0; j < num_tuples_batches; j++)
    {
        uint64_t idx_hash[8];
        avx_murmur_hash_32((uint64_t*)(&(rel[j * 8])), idx_hash);        
        for(uint32_t k = 0; k < 8; k++)
        {
            uint32_t idx = HASH_BIT_MODULO(idx_hash[k], MASK, R);

            my_hist[idx] ++;
        }   
    }
    
    for(uint32_t l = num_tuples_batches * 8; l < num_tuples_reminders; l++)
    {
        uint32_t idx_hash = murmur_hash_32(rel[l].key);
        uint32_t idx = HASH_BIT_MODULO(idx_hash, MASK, R);

        my_hist[idx] ++;
    }*/
#endif
    // compute local prefix sum on hist 
    for(i = 0; i < fanOut; i++){
        sum += my_hist[i];
        my_hist[i] = sum;
    }

    BARRIER_ARRIVE(part->thrargs->barrier, rv);

    // determine the start and end of each cluster 
    for(i = 0; i < my_tid; i++) {
        for(j = 0; j < fanOut; j++)
            output[j] += hist[i][j];
    }
    for(i = my_tid; i < nthreads; i++) {
        for(j = 1; j < fanOut; j++)
            output[j] += hist[i][j-1];
    }

    Tuple<KeyType, PayloadType> * restrict tmp = part->tmp;
    CacheLine<KeyType, PayloadType> buffer[fanOut] __attribute__((aligned(CACHE_LINE_SIZE)));

    for(i = 0; i < fanOut; i++ ) {
        uint64_t off = output[i] + i * padding;
        // pre        = (off + TUPLESPERCACHELINE) & ~(TUPLESPERCACHELINE-1); 
        // pre       -= off;
        output[i]  = off;
        buffer[i].data.slot = off;
    }
    output[fanOut] = part->total_tuples + fanOut * padding;

#ifndef USE_VECTORIZED_MURMUR3_HASH_FOR_RADIX_JOIN
    for(i = 0; i < num_tuples; i++ )
    {
    #ifndef USE_MURMUR3_HASH_FOR_RADIX_JOIN
        uint32_t  idx     = HASH_BIT_MODULO(rel[i].key, MASK, R);
    #else
        uint32_t idx_hash = murmur_hash_32(rel[i].key);
        uint32_t idx = HASH_BIT_MODULO(idx_hash, MASK, R);
    #endif
        uint64_t  slot    = buffer[idx].data.slot;
        Tuple<KeyType, PayloadType> * tup     = (Tuple<KeyType, PayloadType> *)(buffer + idx);
        uint32_t  slotMod = (slot) & ((CACHE_LINE_SIZE/sizeof(Tuple<KeyType, PayloadType>)) - 1);
        tup[slotMod]      = rel[i];

        if(slotMod == ((CACHE_LINE_SIZE/sizeof(Tuple<KeyType, PayloadType>))-1)){
            // write out 64-Bytes with non-temporal store
            store_nontemp_64B<KeyType, PayloadType>((tmp+slot-((CACHE_LINE_SIZE/sizeof(Tuple<KeyType, PayloadType>))-1)), (buffer+idx));
        }
        
        buffer[idx].data.slot = slot+1;
    }
#else
    printf("USE_VECTORIZED_MURMUR3_HASH_FOR_RADIX_JOIN needs to be double checked for correctness \n");
/*    num_tuples_batches = num_tuples / 8;
    num_tuples_reminders = num_tuples % 8;
    i = 0;

    for(j = 0; j < num_tuples_batches; j++)
    {
        uint64_t idx_hash[8];
        avx_murmur_hash_32((uint64_t*)(&(rel[j * 8])), idx_hash);
        for(uint32_t k = 0; k < 8; k++)
        {
            uint32_t idx = HASH_BIT_MODULO(idx_hash[k], MASK, R);

            uint64_t  slot    = buffer[idx].data.slot;
            Tuple<KeyType, PayloadType> * tup     = (Tuple<KeyType, PayloadType> *)(buffer + idx);
            uint32_t  slotMod = (slot) & ((CACHE_LINE_SIZE/sizeof(Tuple<KeyType, PayloadType>)) - 1);
            tup[slotMod]      = rel[j * 8 + k];

            if(slotMod == ((CACHE_LINE_SIZE/sizeof(Tuple<KeyType, PayloadType>))-1)){
                // write out 64-Bytes with non-temporal store 
                store_nontemp_64B<KeyType, PayloadType>((tmp+slot-((CACHE_LINE_SIZE/sizeof(Tuple<KeyType, PayloadType>))-1)), (buffer+idx));
            }

            buffer[idx].data.slot = slot+1;
        }   
    }
    
    for(uint32_t l = num_tuples_batches * 8; l < num_tuples_reminders; l++)
    {
        uint32_t idx_hash = murmur_hash_32(rel[l].key);
        uint32_t idx = HASH_BIT_MODULO(idx_hash, MASK, R);

        uint64_t  slot    = buffer[idx].data.slot;
        Tuple<KeyType, PayloadType> * tup     = (Tuple<KeyType, PayloadType> *)(buffer + idx);
        uint32_t  slotMod = (slot) & ((CACHE_LINE_SIZE/sizeof(Tuple<KeyType, PayloadType>)) - 1);
        tup[slotMod]      = rel[l];

        if(slotMod == ((CACHE_LINE_SIZE/sizeof(Tuple<KeyType, PayloadType>))-1)){
            // write out 64-Bytes with non-temporal store
            store_nontemp_64B<KeyType, PayloadType>((tmp+slot-((CACHE_LINE_SIZE/sizeof(Tuple<KeyType, PayloadType>))-1)), (buffer+idx));
        }
        
        buffer[idx].data.slot = slot+1;
    }
*/  
#endif
    // _mm_sfence (); 
    // write out the remainders in the buffer
    for(i = 0; i < fanOut; i++ ) {
        uint64_t slot  = buffer[i].data.slot;
        uint32_t sz    = (slot) & ((CACHE_LINE_SIZE/sizeof(Tuple<KeyType, PayloadType>)) - 1);
        slot          -= sz;
        for(uint32_t j = 0; j < sz; j++) {
            tmp[slot]  = buffer[i].data.tuples[j];
            slot ++;
        }
    }
}

uint32_t pj_build_rel_r_partition(ETHBucketChainingBuild *build_output, Relation<KeyType, PayloadType> * rel_r_partition, Relation<KeyType, PayloadType> * tmp_r)
{
    uint32_t * next, * bucket;
    const uint32_t numR = rel_r_partition->num_tuples;
    uint32_t N = numR;

    NEXT_POW_2(N);

    uint32_t numR_from_build = N;

    const uint32_t MASK = (N-1) << (NUM_RADIX_BITS);
    next   = (uint32_t*) malloc(sizeof(uint32_t) * numR);
    bucket = (uint32_t*) calloc(N, sizeof(uint32_t));

//#ifndef USE_VECTORIZED_MURMUR3_HASH_FOR_RADIX_JOIN
    for(uint32_t i=0; i < numR; )
    {
    //#ifndef USE_MURMUR3_HASH_FOR_RADIX_JOIN
        uint32_t idx = HASH_BIT_MODULO(rel_r_partition->tuples[i].key, MASK, NUM_RADIX_BITS);
    /*#else
        uint32_t idx_hash = murmur_hash_32(rel_r_partition->tuples[i].key);
        uint32_t idx = HASH_BIT_MODULO(idx_hash, MASK, NUM_RADIX_BITS);
    #endif*/    
        next[i]      = bucket[idx];
        bucket[idx]  = ++i;     /* we start pos's from 1 instead of 0 */

        /* Enable the following tO avoid the code elimination
           when running probe only for the time break-down experiment */
        /* matches += idx; */
    }
/*#else
    int64_t num_R_batches = numR / 8;
    uint32_t num_R_reminders = numR % 8;
    int64_t i = 0;

    for(int64_t j = 0; j < num_R_batches; j++)
    {
        uint64_t idx_hash[8];
        avx_murmur_hash_32((uint64_t*)(&(rel_r_partition->tuples[j * 8])), idx_hash);
        for(uint32_t k = 0; k < 8; k++)
        {
            uint32_t idx = HASH_BIT_MODULO(idx_hash[k], MASK, NUM_RADIX_BITS);
        
            next[i]      = bucket[idx];
            bucket[idx]  = ++i;     //we start pos's from 1 instead of 0
        }   
    }
    
    for(uint32_t l = num_R_batches * 8; l < num_R_reminders; l++)
    {
        uint32_t idx_hash = murmur_hash_32(rel_r_partition->tuples[l].key);
        uint32_t idx = HASH_BIT_MODULO(idx_hash, MASK, NUM_RADIX_BITS);
        
        next[i]      = bucket[idx];
        bucket[idx]  = ++i;     // we start pos's from 1 instead of 0 
    }
#endif*/

    ((ETHBucketChainingBuild *) build_output)->hist = bucket;
    ((ETHBucketChainingBuild *) build_output)->next = next;
    
    return numR_from_build;
}

uint64_t pj_probe_rel_s_partition(Relation<KeyType, PayloadType> * rel_r_partition, uint32_t numR_from_build, Relation<KeyType, PayloadType> * rel_s_partition, ETHBucketChainingBuild *build_output)
{
    uint32_t N = numR_from_build;
    const uint32_t MASK = (N-1) << (NUM_RADIX_BITS);                               
    uint64_t matches = 0;
    uint32_t * next, * bucket;
    const Tuple<KeyType, PayloadType> * const Rtuples = rel_r_partition->tuples;
    const Tuple<KeyType, PayloadType> * const Stuples = rel_s_partition->tuples;
    const uint32_t        numS    = rel_s_partition->num_tuples;   
    bucket = build_output->hist;    
    next = build_output->next;

    /* Disable the following loop for no-probe for the break-down experiments */
    /* PROBE- LOOP */
 //#ifndef USE_VECTORIZED_MURMUR3_HASH_FOR_RADIX_JOIN   
    for(uint32_t i=0; i < numS; i++ )
    {
    //#ifndef USE_MURMUR3_HASH_FOR_RADIX_JOIN
        uint32_t idx = HASH_BIT_MODULO(Stuples[i].key, MASK, NUM_RADIX_BITS);
    /*#else
        uint32_t idx_hash = murmur_hash_32(Stuples[i].key);
        uint32_t idx = HASH_BIT_MODULO(idx_hash, MASK, NUM_RADIX_BITS);
    #endif*/
        for(int hit = bucket[idx]; hit > 0; hit = next[hit-1]){
            if(Stuples[i].key == Rtuples[hit-1].key){
                matches ++;
            }
        }
    }
/*#else
    int64_t num_S_batches = numS / 8;
    uint32_t num_S_reminders = numS % 8;
    int64_t i = 0;

    for(int64_t j = 0; j < num_S_batches; j++)
    {
        uint64_t idx_hash[8];
        avx_murmur_hash_32((uint64_t*)(&(Stuples[j * 8])), idx_hash);
        for(uint32_t k = 0; k < 8; k++)
        {
            uint32_t idx = HASH_BIT_MODULO(idx_hash[k], MASK, NUM_RADIX_BITS);
            i = j * 8 + k;

            for(int hit = bucket[idx]; hit > 0; hit = next[hit-1]){
                if(Stuples[i].key == Rtuples[hit-1].key){
                    matches ++;
                }
            }
        }   
    }
    
    for(uint32_t l = num_S_batches * 8; l < num_S_reminders; l++)
    {
        uint32_t idx_hash = murmur_hash_32(Stuples[l].key);
        uint32_t idx = HASH_BIT_MODULO(idx_hash, MASK, NUM_RADIX_BITS);
        
        for(int hit = bucket[idx]; hit > 0; hit = next[hit-1]){
            if(Stuples[l].key == Rtuples[hit-1].key){
                matches ++;
            }
        }
    }
#endif*/

    /* PROBE-LOOP END  */
    
    /* clean up temp */
    free(bucket);
    free(next);

    return matches;
}

// Based on the implementation of radix hash join from ETH
void radix_cluster(Relation<KeyType, PayloadType> * restrict outRel, 
                   Relation<KeyType, PayloadType> * restrict inRel,
                   int32_t * restrict hist, int R, int D)
{
    int64_t i;
    uint32_t M = ((1 << D) - 1) << R;
    uint32_t offset;
    uint32_t fanOut = 1 << D;

    /* the following are fixed size when D is same for all the passes,
        and can be re-used from call to call. Allocating in this function 
        just in case D differs from call to call. */
    uint32_t dst[fanOut];

    /* count tuples per cluster */
    #ifndef USE_VECTORIZED_MURMUR3_HASH_FOR_RADIX_JOIN
    for( i=0; i < inRel->num_tuples; i++ )
    {
    #ifndef USE_MURMUR3_HASH_FOR_RADIX_JOIN
        uint32_t idx = HASH_BIT_MODULO(inRel->tuples[i].key, M, R);
    #else
        uint32_t idx_hash = murmur_hash_32(inRel->tuples[i].key);
        uint32_t idx = HASH_BIT_MODULO(idx_hash, M, R);
    #endif
        hist[idx]++;
    }
    #else

    int64_t num_tuples_batches = inRel->num_tuples / 8;
    uint32_t num_tuples_reminders = inRel->num_tuples % 8;
    i = 0;

    for(int64_t j = 0; j < num_tuples_batches; j++)
    {
        uint64_t idx_hash[8];
        avx_murmur_hash_32((uint64_t*)(&(inRel->tuples[j * 8])), idx_hash);
        for(uint32_t k = 0; k < 8; k++)
        {
            uint32_t idx = HASH_BIT_MODULO(idx_hash[k], M, R);

            hist[idx] ++;
        }   
    }

    for(uint32_t l = num_tuples_batches * 8; l < num_tuples_reminders; l++)
    {
        uint32_t idx_hash = murmur_hash_32(inRel->tuples[l].key);
        uint32_t idx = HASH_BIT_MODULO(idx_hash, M, R);

        hist[idx] ++;
    }

    #endif
    offset = 0;
    int small_padding_tupples = 0;
    //int small_padding_tupples = SMALL_PADDING_TUPLES_MULTIPLIER * CACHE_LINE_SIZE/sizeof(Tuple<KeyType, PayloadType>);
    /* determine the start and end of each cluster depending on the counts. */
    for ( i=0; i < fanOut; i++ ) {
        /* dst[i]      = outRel->tuples + offset; */
        /* determine the beginning of each partitioning by adding some
            padding to avoid L1 conflict misses during scatter. */
        dst[i] = offset + i * small_padding_tupples;
        offset += hist[i];
    }

    /* copy tuples to their corresponding clusters at appropriate offsets */
    #ifndef USE_VECTORIZED_MURMUR3_HASH_FOR_RADIX_JOIN   
    for( i=0; i < inRel->num_tuples; i++ )
    {
    #ifndef USE_MURMUR3_HASH_FOR_RADIX_JOIN
        uint32_t idx   = HASH_BIT_MODULO(inRel->tuples[i].key, M, R);
    #else
        uint32_t idx_hash = murmur_hash_32(inRel->tuples[i].key);
        uint32_t idx = HASH_BIT_MODULO(idx_hash, M, R);
    #endif
        outRel->tuples[ dst[idx] ] = inRel->tuples[i];
        ++dst[idx];
    }
    #else
    num_tuples_batches = inRel->num_tuples / 8;
    num_tuples_reminders = inRel->num_tuples % 8;
    i = 0;

    for(int64_t j = 0; j < num_tuples_batches; j++)
    {
        uint64_t idx_hash[8];
        avx_murmur_hash_32((uint64_t*)(&(inRel->tuples[j * 8])), idx_hash);
        for(uint32_t k = 0; k < 8; k++)
        {
            uint32_t idx = HASH_BIT_MODULO(idx_hash[k], M, R);

            outRel->tuples[ dst[idx] ] = inRel->tuples[j * 8 + k];
            ++dst[idx];
        }   
    }

    for(uint32_t l = num_tuples_batches * 8; l < num_tuples_reminders; l++)
    {
        uint32_t idx_hash = murmur_hash_32(inRel->tuples[l].key);
        uint32_t idx = HASH_BIT_MODULO(idx_hash, M, R);

        outRel->tuples[ dst[idx] ] = inRel->tuples[l];
        ++dst[idx];
    }
    #endif
}

void serial_radix_partition(TaskType * const task, 
                            TaskQueue<KeyType, PayloadType, TaskType> * join_queue, 
                            const int R, const int D) 
{
    int i;
    uint32_t offsetR = 0, offsetS = 0;
    const int fanOut = 1 << D;  /*(NUM_RADIX_BITS / NUM_PASSES);*/
    int32_t * outputR, * outputS;

    int small_padding_tupples = 0;
    //int small_padding_tupples = SMALL_PADDING_TUPLES_MULTIPLIER * CACHE_LINE_SIZE/sizeof(Tuple<KeyType, PayloadType>);

    outputR = (int32_t*)calloc(fanOut+1, sizeof(int32_t));
    outputS = (int32_t*)calloc(fanOut+1, sizeof(int32_t));
    /* TODO: measure the effect of memset() */
    /* memset(outputR, 0, fanOut * sizeof(int32_t)); */
    radix_cluster(&task->tmpR, &task->relR, outputR, R, D);

    /* memset(outputS, 0, fanOut * sizeof(int32_t)); */
    radix_cluster(&task->tmpS, &task->relS, outputS, R, D);

    /* TaskType t; */
    for(i = 0; i < fanOut; i++) {
        if(outputR[i] > 0 && outputS[i] > 0) {
            TaskType * t = task_queue_get_slot_atomic(join_queue);
            t->relR.num_tuples = outputR[i];
            t->relR.tuples = task->tmpR.tuples + offsetR 
                             + i * small_padding_tupples;
            t->tmpR.tuples = task->relR.tuples + offsetR 
                             + i * small_padding_tupples;
            offsetR += outputR[i];

            t->relS.num_tuples = outputS[i];
            t->relS.tuples = task->tmpS.tuples + offsetS 
                             + i * small_padding_tupples;
            t->tmpS.tuples = task->relS.tuples + offsetS 
                             + i * small_padding_tupples;
            offsetS += outputS[i];

            /* task_queue_copy_atomic(join_queue, &t); */
            task_queue_add_atomic(join_queue, t);
        } 
        else {
            offsetR += outputR[i];
            offsetS += outputS[i];
        }
    }
    free(outputR);
    free(outputS);
}

#ifdef BUILD_RMI_FROM_TWO_DATASETS
void sample_and_train_models_threaded(ETHRadixJoinThread<KeyType, PayloadType, TaskType> * args)
{
    int rv;
    int tid = args->my_tid;

    //----------------------------------------------------------//
    //                           SAMPLE                         //
    //----------------------------------------------------------//

    // Determine sample size
    unsigned int INPUT_SZ_R = args->original_relR->num_tuples;
    unsigned int INPUT_SZ_S = args->original_relS->num_tuples;
    const unsigned int SAMPLE_SZ_R = std::min<unsigned int>(
        INPUT_SZ_R, std::max<unsigned int>(args->p.sampling_rate * INPUT_SZ_R,
                                        RMI<KeyType, PayloadType>::Params::MIN_SORTING_SIZE));
    const unsigned int SAMPLE_SZ_S = std::min<unsigned int>(
        INPUT_SZ_S, std::max<unsigned int>(args->p.sampling_rate * INPUT_SZ_S,
                                        RMI<KeyType, PayloadType>::Params::MIN_SORTING_SIZE));
    //const unsigned int SAMPLE_SZ = SAMPLE_SZ_R + SAMPLE_SZ_S;        
    
    // Start sampling
    Tuple<KeyType, PayloadType> *   relR_start_sampling_ptr = args->relR;
    Tuple<KeyType, PayloadType> *   relR_end_sampling_ptr = args->relR + (args->numR - 1);
    Tuple<KeyType, PayloadType> *   relS_start_sampling_ptr = args->relS;
    Tuple<KeyType, PayloadType> *   relS_end_sampling_ptr = args->relS + (args->numS - 1);

    uint32_t * sample_count = args->sample_count; 
    Tuple<KeyType, PayloadType> * tmp_training_sample = args->rmi->tmp_training_sample + args->tmp_training_sample_offset;
    
    uint32_t * sample_count_R = args->sample_count_R; 
    Tuple<KeyType, PayloadType> * tmp_training_sample_R = args->rmi->tmp_training_sample_R + args->tmp_training_sample_R_offset;
    unsigned int offset_R = static_cast<unsigned int>(1. * INPUT_SZ_R / SAMPLE_SZ_R);
    for (auto i = relR_start_sampling_ptr; i <= relR_end_sampling_ptr; i += offset_R) {
      // NOTE:  We don't directly assign SAMPLE_SZ to rmi.training_sample_sz
      //        to avoid issues with divisibility
      tmp_training_sample[sample_count[tid]] = *i;
      ++sample_count[tid];

      //tmp_training_sample_R[sample_count_R[tid]] = *i;
      //++sample_count_R[tid];
    }
    
    uint32_t * sample_count_S = args->sample_count_S;
    Tuple<KeyType, PayloadType> * tmp_training_sample_S = args->rmi->tmp_training_sample_S + args->tmp_training_sample_S_offset;
    unsigned int offset_S = static_cast<unsigned int>(1. * INPUT_SZ_S / SAMPLE_SZ_S);
    for (auto i = relS_start_sampling_ptr; i <= relS_end_sampling_ptr; i += offset_S) {
      // NOTE:  We don't directly assign SAMPLE_SZ to rmi.training_sample_sz
      //        to avoid issues with divisibility
      tmp_training_sample[sample_count[tid]] = *i;
      ++sample_count[tid];

      //tmp_training_sample_S[sample_count_S[tid]] = *i;
      //++sample_count_S[tid];
    }

    BARRIER_ARRIVE(args->barrier, rv);

    #ifdef USE_AVXSORT_AS_STD_SORT          

    uint32_t total_sample_count = 0; 

    if(tid == 0)
    {
      for(int i = 0; i < NUM_THREADS; i++)
        total_sample_count += sample_count[i];

      Tuple<KeyType, PayloadType> * sorted_training_sample = args->rmi->sorted_training_sample;      
      int64_t * inputptr =  (int64_t *)(args->rmi->tmp_training_sample);
      int64_t * outputptr = (int64_t *)(sorted_training_sample);
      avxsort_int64(&inputptr, &outputptr, total_sample_count);
      Tuple<KeyType, PayloadType>* tmp_outputptr = (Tuple<KeyType, PayloadType>*) outputptr;
      for(unsigned int k = 0; k < total_sample_count; k++){
        sorted_training_sample[k] = tmp_outputptr[k];
      } 
      args->rmi->training_sample = &(sorted_training_sample);
      args->rmi->training_sample_size = total_sample_count;
    }

    /*if(tid == 1)
    {
      uint32_t total_sample_count_R = 0; 
      for(int i = 0; i < NUM_THREADS; i++)
        total_sample_count_R += sample_count_R[i];

      Tuple<KeyType, PayloadType> * sorted_training_sample_R = args->rmi->sorted_training_sample_R;
      int64_t * inputptr_R =  (int64_t *)(args->rmi->tmp_training_sample_R);
      int64_t * outputptr_R = (int64_t *)(sorted_training_sample_R);
      avxsort_int64(&inputptr_R, &outputptr_R, total_sample_count_R);
      Tuple<KeyType, PayloadType>* tmp_outputptr_R = (Tuple<KeyType, PayloadType>*) outputptr_R;
      for(unsigned int k = 0; k < total_sample_count_R; k++){
        sorted_training_sample_R[k] = tmp_outputptr_R[k];
      } 
      args->rmi->training_sample_R = &(sorted_training_sample_R);
      args->rmi->training_sample_size_R = total_sample_count_R;
    }

    if(tid == 2)
    {
      uint32_t total_sample_count_S = 0; 
      for(int i = 0; i < NUM_THREADS; i++)
        total_sample_count_S += sample_count_S[i];

      Tuple<KeyType, PayloadType> * sorted_training_sample_S = args->rmi->sorted_training_sample_S;
      int64_t * inputptr_S =  (int64_t *)(args->rmi->tmp_training_sample_S);
      int64_t * outputptr_S = (int64_t *)(args->rmi->sorted_training_sample_S);
      avxsort_int64(&inputptr_S, &outputptr_S, total_sample_count_S);
      Tuple<KeyType, PayloadType>* tmp_outputptr_S = (Tuple<KeyType, PayloadType>*) outputptr_S;
      for(unsigned int k = 0; k < total_sample_count_S; k++){
        sorted_training_sample_S[k] = tmp_outputptr_S[k];
      } 
      args->rmi->training_sample_S = &(sorted_training_sample_S);
      args->rmi->training_sample_size_S = total_sample_count_S;
    }*/
    #else
    uint32_t total_sample_count = 0;
    if(tid == 0)
    {
      for(int i = 0; i < NUM_THREADS; i++)
        total_sample_count += sample_count[i];

      std::sort((int64_t *)(args->rmi->tmp_training_sample), (int64_t *)(args->rmi->tmp_training_sample) + total_sample_count);
      args->rmi->training_sample = &(args->rmi->tmp_training_sample);
      args->rmi->training_sample_size = total_sample_count;
    }

    /*if(tid == 1)
    {
      uint32_t total_sample_count_R = 0; 
      for(int i = 0; i < NUM_THREADS; i++)
        total_sample_count_R += sample_count_R[i];

      std::sort((int64_t *)(args->rmi->tmp_training_sample_R), (int64_t *)(args->rmi->tmp_training_sample_R) + total_sample_count_R);
      args->rmi->training_sample_R = &(args->rmi->tmp_training_sample_R);
      args->rmi->training_sample_size_R = total_sample_count_R;
    }

    if(tid == 2)
    {
      uint32_t total_sample_count_S = 0; 
      for(int i = 0; i < NUM_THREADS; i++)
        total_sample_count_S += sample_count_S[i];

      std::sort((int64_t *)(args->rmi->tmp_training_sample_S), (int64_t *)(args->rmi->tmp_training_sample_S) + total_sample_count_S );
      args->rmi->training_sample_S = &(args->rmi->tmp_training_sample_S);
      args->rmi->training_sample_size_S = total_sample_count_S;
    }*/
    #endif

    //----------------------------------------------------------//
    //                     TRAIN THE MODELS                     //
    //----------------------------------------------------------//

    if(tid == 0)
    {
       // Stop early if the array is identical
      if (((*(args->rmi->training_sample))[0]).key == ((*(args->rmi->training_sample))[total_sample_count - 1]).key) 
      {
        return;
      }
         
      // Populate the training data for the root model
      vector<vector<vector<training_point<KeyType, PayloadType>>>> * training_data = args->training_data;
      for (unsigned int i = 0; i < total_sample_count; ++i) {
        (*training_data)[0][0].push_back({((*(args->rmi->training_sample))[i]), 1. * i / total_sample_count});
      }

      // Train the root model using linear interpolation
      auto *current_training_data = &(*training_data)[0][0];
      typename RMI<KeyType, PayloadType>::linear_model *current_model = &args->rmi->models[0][0];

      // Find the min and max values in the training set
      training_point<KeyType, PayloadType> min = current_training_data->front();
      training_point<KeyType, PayloadType> max = current_training_data->back();

      // Calculate the slope and intercept terms
      current_model->slope =
          1. / (max.x.key - min.x.key);  // Assuming min.y = 0 and max.y = 1
      current_model->intercept = -current_model->slope * min.x.key;

      // Extrapolate for the number of models in the next layer
      current_model->slope *= args->p.arch[1] - 1;
      current_model->intercept *= args->p.arch[1] - 1;
#ifndef RUN_LEARNED_TECHNIQUES_WITH_FIRST_LEVEL_ONLY
      // Populate the training data for the next layer
      for (const auto &d : *current_training_data) {
        // Predict the model index in next layer
        //unsigned int rank = current_model->slope * d.x.key + current_model->intercept;
        unsigned int rank = round(current_model->slope * d.x.key*1.00 + current_model->intercept);

        // Normalize the rank between 0 and the number of models in the next layer
        rank =
            std::max(static_cast<unsigned int>(0), std::min(args->p.arch[1] - 1, rank));
        
        //if(d.x.key >299 && d.x.key < 400)
            //printf("training: key %ld rank %ld \n", d.x.key, rank);
        // Place the data in the predicted training bucket
        (*training_data)[1][rank].push_back(d);
      }

      // Train the leaf models
      for (unsigned int model_idx = 0; model_idx < args->p.arch[1]; ++model_idx) {
        // Update iterator variables
        current_training_data = &(*training_data)[1][model_idx];
        current_model = &args->rmi->models[1][model_idx];

        // Interpolate the min points in the training buckets
        if (model_idx ==
            0) {  // The current model is the first model in the current layer

          if (current_training_data->size() <
              2) {  // Case 1: The first model in this layer is empty
            current_model->slope = 0;
            current_model->intercept = 0;

            // Insert a fictive training point to avoid propagating more than one
            // empty initial models.
            training_point<KeyType, PayloadType> tp;
            tp.x.key = 0;
            tp.x.payload = 0;
            tp.y = 0;
            current_training_data->push_back(tp);
          } else {  // Case 2: The first model in this layer is not empty

            min = current_training_data->front();
            max = current_training_data->back();

            current_model->slope =
                (max.y)*1.00 / (max.x.key - min.x.key);  // Hallucinating as if min.y = 0
            current_model->intercept = min.y - current_model->slope * min.x.key;
          }
        } else if (model_idx == args->p.arch[1] - 1) {
          if (current_training_data
                  ->empty()) {  // Case 3: The final model in this layer is empty

            current_model->slope = 0;
            current_model->intercept = 1;
          } else {  // Case 4: The last model in this layer is not empty

            min = (*training_data)[1][model_idx - 1].back();
            max = current_training_data->back();

            current_model->slope =
                (1 - min.y) * 1.00 / (max.x.key - min.x.key);  // Hallucinating as if max.y = 1
            current_model->intercept = min.y - current_model->slope * min.x.key;
          }
        } else {  // The current model is not the first model in the current layer

          if (current_training_data
                  ->empty()) {  // Case 5: The intermediate model in
            // this layer is empty
            current_model->slope = 0;
            current_model->intercept =
                (*training_data)[1][model_idx - 1].back().y;  // If the previous model
                                                          // was empty too, it will
                                                          // use the fictive
                                                          // training points

            // Insert a fictive training point to avoid propagating more than one
            // empty initial models.
            // NOTE: This will _NOT_ throw to DIV/0 due to identical x's and y's
            // because it is working backwards.
            training_point<KeyType, PayloadType> tp;
            tp.x = (*training_data)[1][model_idx - 1].back().x;
            tp.y = (*training_data)[1][model_idx - 1].back().y;
            current_training_data->push_back(tp);
          } else {  // Case 6: The intermediate leaf model is not empty

            min = (*training_data)[1][model_idx - 1].back();
            max = current_training_data->back();

            current_model->slope = (max.y - min.y) * 1.00 / (max.x.key - min.x.key);
            current_model->intercept = min.y - current_model->slope * min.x.key;
            //if(model_idx == 5)
            //printf("min.y %lf max.y %lf min.x.key %ld max.x.key %ld current_model->slope %lf current_model->intercept %lf\n", min.y, max.y, min.x.key, max.x.key, current_model->slope, current_model->intercept);
          }
        }
      }
#endif
      // NOTE:
      // The last stage (layer) of this model contains weights that predict the CDF
      // of the keys (i.e. Range is [0-1])
      // When using this model to predict the position of the keys in the sorted
      // order, you MUST scale the
      // weights of the last layer to whatever range you are predicting for. The
      // inner layers of the model have
      // already been extrapolated to the length of the stage.git
      //
      // This is a design choice to help with the portability of the model.
      //
      args->rmi->trained = true;         
    }

}
#endif

void * pj_join_thread(void * param)
{
    ETHRadixJoinThread<KeyType, PayloadType, TaskType> * args   = (ETHRadixJoinThread<KeyType, PayloadType, TaskType> *) param;
    int rv;   uint64_t results = 0; int i; 
    int deltaT = 0; struct timeval t1, t2;

    const int fanOut = 1 << (NUM_RADIX_BITS / NUM_PASSES);//PASS1RADIXBITS;
    const int R = (NUM_RADIX_BITS / NUM_PASSES);//PASS1RADIXBITS;
    const int D = (NUM_RADIX_BITS - (NUM_RADIX_BITS / NUM_PASSES));//PASS2RADIXBITS;

    int32_t my_tid = args->my_tid;
    
    TaskType * task;
    TaskQueue<KeyType, PayloadType, TaskType> * part_queue;
    TaskQueue<KeyType, PayloadType, TaskType> * join_queue;
#if SKEW_HANDLING
    TaskQueue<KeyType, PayloadType, TaskType> * skew_queue;
#endif

    int64_t * outputR = (int64_t *) calloc((fanOut+1), sizeof(int64_t));
    int64_t * outputS = (int64_t *) calloc((fanOut+1), sizeof(int64_t));
    MALLOC_CHECK((outputR && outputS));

    #ifdef DEVELOPMENT_MODE
    int numaid = get_numa_id_develop(my_tid);
    #else
    int numaid = get_numa_id(my_tid);
    #endif

    part_queue = args->part_queue[numaid];
    join_queue = args->join_queue[numaid];

#if SKEW_HANDLING
    skew_queue = args->skew_queue;
#endif

    args->histR[my_tid] = (int32_t *) calloc(fanOut, sizeof(int32_t));
    args->histS[my_tid] = (int32_t *) calloc(fanOut, sizeof(int32_t));

    /* in the first pass, partitioning is done together by all threads */

    args->parts_processed = 0;


#ifdef RUN_LEARNED_TECHNIQUES        
    for (int rp = 0; rp < RUN_NUMS; ++rp) 
    {
        if(args->my_tid == 0)
            init_models_training_data_and_sample_counts<KeyType, PayloadType>(args->training_data, args->p.arch, 
                    args->sample_count, args->sample_count_R, args->sample_count_S, NUM_THREADS);

        BARRIER_ARRIVE(args->barrier, rv);
        if(args->my_tid == 0){
            gettimeofday(&t1, NULL);
        }

        sample_and_train_models_threaded(args);

        BARRIER_ARRIVE(args->barrier, rv);

        if(args->my_tid == 0){
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
    }
#endif



    if (my_tid == 0) {

#ifndef RUN_LEARNED_TECHNIQUES        
        strcpy(pj_partition_pfun[0].fun_name, "Naive");

        pj_partition_pfun[0].fun_ptr = pj_partition_rel_segment_pass1;

        pj_partition_pf_num = 1;

        strcpy(pj_build_pfun[0].fun_name, "Naive");

        pj_build_pfun[0].fun_ptr = pj_build_rel_r_partition;

        pj_build_pf_num = 1;

        strcpy(pj_probe_pfun[0].fun_name, "Naive");

        pj_probe_pfun[0].fun_ptr = pj_probe_rel_s_partition;        

        pj_probe_pf_num = 1;
#else
        printf("Learned versions are not supported yet! \n");
        //TODO: to be done
        /*strcpy(npj_pfun[0].fun_name, "Learned");
        strcpy(npj_pfun[1].fun_name, "Learned IMV");

        npj_pfun[0].fun_ptr = npj_build_rel_r_partition_learned;
        npj_pfun[1].fun_ptr = npj_build_rel_r_partition_learned_imv;

        strcpy(npj_pfun1[0].fun_name, "Learned");
        strcpy(npj_pfun1[1].fun_name, "Learned IMV");
        
        npj_pfun1[0].fun_ptr = npj_probe_rel_s_partition_learned;
        npj_pfun1[1].fun_ptr = npj_probe_rel_s_partition_learned_imv;

        npj_pf_num = 2;*/
#endif        
    }
    
    #ifdef PERF_COUNTERS
        if(my_tid == 0){
            // TODO: performance counters to be implemented
        }
    #endif

    /* wait at a barrier until each thread starts and then start the timer */
    BARRIER_ARRIVE(args->barrier, rv);


    ETHPartition<KeyType, PayloadType, TaskType, ETHRadixJoinThread<KeyType, PayloadType, TaskType>> part;

    for (int fid = 0; fid < pj_partition_pf_num; ++fid) 
    {
        for (int rp = 0; rp < RUN_NUMS; ++rp) 
        {
            if(args->my_tid == 0){
                gettimeofday(&args->start_time, NULL);
            #ifndef DEVELOPMENT_MODE
                //args->e_start_to_partition.startCounters();
            #endif
            }

            /********** 1st pass of multi-pass partitioning ************/
            part.R       = 0;
            part.D       = NUM_RADIX_BITS / NUM_PASSES; //PASS1RADIXBITS
            part.thrargs = args;
            part.padding = padding_tuples;

            /* 1. partitioning for relation R */
            part.rel          = args->relR;
            part.tmp          = args->tmpR;
            part.hist         = args->histR;
            part.output       = outputR;
            part.num_tuples   = args->numR;
            part.total_tuples = args->totalR;
            part.relidx       = 0;
        
            pj_partition_pfun[fid].fun_ptr(&part);

            /* 2. partitioning for relation S */
            part.rel          = args->relS;
            part.tmp          = args->tmpS;
            part.hist         = args->histS;
            part.output       = outputS;
            part.num_tuples   = args->numS;
            part.total_tuples = args->totalS;
            part.relidx       = 1;
            
            pj_partition_pfun[fid].fun_ptr(&part);

            BARRIER_ARRIVE(args->barrier, rv);
            // end of 1st partitioning phase 

            #if SKEW_HANDLING
                // experimental skew threshold
                // const int thresh1 = MAX((1<<D), (1<<R)) * THRESHOLD1(args->nthreads);
                // const int thresh1 = MAX(args->totalR, arg->totalS)/MAX((1<<D),(1<<R));
                const int thresh1 = SKEWNESS_THRESHOLD_MULTIPLIER * l1_cache_tuples * args->nthreads;
            #endif

            // 3. first thread creates partitioning tasks for 2nd pass
            if(my_tid == 0) {
                for(i = 0; i < fanOut; i++) {
                    int32_t ntupR = outputR[i+1] - outputR[i] - padding_tuples;
                    int32_t ntupS = outputS[i+1] - outputS[i] - padding_tuples;

        #if SKEW_HANDLING

                    if(ntupR > thresh1 || ntupS > thresh1){

                        TaskType * t = task_queue_get_slot<KeyType, PayloadType, TaskType>(skew_queue);

                        t->relR.num_tuples = t->tmpR.num_tuples = ntupR;
                        t->relR.tuples = args->tmpR + outputR[i];
                        t->tmpR.tuples = args->relR + outputR[i];

                        t->relS.num_tuples = t->tmpS.num_tuples = ntupS;
                        t->relS.tuples = args->tmpS + outputS[i];
                        t->tmpS.tuples = args->relS + outputS[i];

                        task_queue_add<KeyType, PayloadType, TaskType>(skew_queue, t);
                    } 
                    else
        #endif
                    if(ntupR > 0 && ntupS > 0) {
                        // Determine the NUMA node of each partition:

                        #ifdef DEVELOPMENT_MODE
                        int pq_idx = get_numa_node_of_address_develop(NULL);
                        #else
                        void * ptr = (void*)&((args->tmpR + outputR[i])[0]);
                        int pq_idx = get_numa_node_of_address(ptr);
                        #endif

                        TaskQueue<KeyType, PayloadType, TaskType> * numalocal_part_queue = args->part_queue[pq_idx];

                        TaskType * t = task_queue_get_slot<KeyType, PayloadType, TaskType>(numalocal_part_queue);

                        t->relR.num_tuples = t->tmpR.num_tuples = ntupR;
                        t->relR.tuples = args->tmpR + outputR[i];
                        t->tmpR.tuples = args->relR + outputR[i];

                        t->relS.num_tuples = t->tmpS.num_tuples = ntupS;
                        t->relS.tuples = args->tmpS + outputS[i];
                        t->tmpS.tuples = args->relS + outputS[i];

                        task_queue_add<KeyType, PayloadType, TaskType>(numalocal_part_queue, t);
                    }
                }

                // debug partitioning task queue
                //DEBUGMSG(1, "Pass-2: # partitioning tasks = %d\n", part_queue->count);

                // DEBUG NUMA MAPPINGS
                // printf("Correct NUMA-mappings = %d, Wrong = %d\n",
                //        correct_numa_mapping, wrong_numa_mapping); 
                // printf("Counts -- 0=%d, 1=%d, 2=%d, 3=%d\n",  
                //        counts[0], counts[1], counts[2], counts[3]);
            }

            BARRIER_ARRIVE(args->barrier, rv);
            
            // 2nd pass of multi-pass partitioning
            // 4. now each thread further partitions and add to join task queue

        #if NUM_PASSES==1
            TaskQueue<KeyType, PayloadType, TaskType> * swap = join_queue;
            join_queue = part_queue;
            part_queue = swap;
            
        #elif NUM_PASSES==2

            while((task = task_queue_get_atomic<KeyType, PayloadType, TaskType>(part_queue))){
                serial_radix_partition(task, join_queue, R, D);
            }

        #else
        #warning Only 2-pass partitioning is implemented, set NUM_PASSES to 2!
        #endif

        #if SKEW_HANDLING
            // Partitioning pass-2 for skewed relations
            part.R         = R;
            part.D         = D;
            part.thrargs   = args;
            part.padding   = small_padding_tupples;

            while(1) {
                if(my_tid == 0) {
                    *args->skewtask = task_queue_get_atomic<KeyType, PayloadType, TaskType>(skew_queue);
                }
                BARRIER_ARRIVE(args->barrier, rv);
                if( *args->skewtask == NULL)
                    break;

                //DEBUGMSG((my_tid==0), "Got skew task = R: %d, S: %d\n", 
                //         (*args->skewtask)->relR.num_tuples,
                //         (*args->skewtask)->relS.num_tuples);

                int32_t numperthr = (*args->skewtask)->relR.num_tuples / args->nthreads;
                const int fanOut2 = (1 << D);

                free(outputR);
                free(outputS);

                outputR = (int64_t*) calloc(fanOut2 + 1, sizeof(int64_t));
                outputS = (int64_t*) calloc(fanOut2 + 1, sizeof(int64_t));

                free(args->histR[my_tid]);
                free(args->histS[my_tid]);

                args->histR[my_tid] = (int32_t*) calloc(fanOut2, sizeof(int32_t));
                args->histS[my_tid] = (int32_t*) calloc(fanOut2, sizeof(int32_t));

                BARRIER_ARRIVE(args->barrier, rv);

                // 1. partitioning for relation R
                part.rel          = (*args->skewtask)->relR.tuples + my_tid * numperthr;
                part.tmp          = (*args->skewtask)->tmpR.tuples;
                part.hist         = args->histR;
                part.output       = outputR;
                part.num_tuples   = (my_tid == (args->nthreads-1)) ? 
                                    ((*args->skewtask)->relR.num_tuples - my_tid * numperthr) 
                                    : numperthr;
                part.total_tuples = (*args->skewtask)->relR.num_tuples;
                part.relidx       = 2; // meaning this is pass-2, no syncstats 
                pj_partition_pfun[fid].fun_ptr(&part);

                numperthr = (*args->skewtask)->relS.num_tuples / args->nthreads;
                // 2. partitioning for relation S
                part.rel          = (*args->skewtask)->relS.tuples + my_tid * numperthr;
                part.tmp          = (*args->skewtask)->tmpS.tuples;
                part.hist         = args->histS;
                part.output       = outputS;
                part.num_tuples   = (my_tid == (args->nthreads-1)) ? 
                                    ((*args->skewtask)->relS.num_tuples - my_tid * numperthr)
                                    : numperthr;
                part.total_tuples = (*args->skewtask)->relS.num_tuples;
                part.relidx       = 2; // meaning this is pass-2, no syncstats
                pj_partition_pfun[fid].fun_ptr(&part);

                BARRIER_ARRIVE(args->barrier, rv);

                // first thread adds join tasks
                if(my_tid == 0) {
                    const int THR1 = l1_cache_tuples * args->nthreads;

                    for(i = 0; i < fanOut2; i++) {
                        int32_t ntupR = outputR[i+1] - outputR[i] - small_padding_tupples;
                        int32_t ntupS = outputS[i+1] - outputS[i] - small_padding_tupples;
                        if(ntupR > THR1 || ntupS > THR1){

                            //DEBUGMSG(1, "Large join task = R: %d, S: %d\n", ntupR, ntupS);

                            // use part_queue temporarily
                            for(int k=0; k < args->nthreads; k++) {
                                int ns = (k == args->nthreads-1)
                                        ? (ntupS - k*(ntupS/args->nthreads))
                                        : (ntupS/args->nthreads);
                                TaskType * t = task_queue_get_slot<KeyType, PayloadType, TaskType>(part_queue);

                                t->relR.num_tuples = t->tmpR.num_tuples = ntupR;
                                t->relR.tuples = (*args->skewtask)->tmpR.tuples + outputR[i];
                                t->tmpR.tuples = (*args->skewtask)->relR.tuples + outputR[i];

                                t->relS.num_tuples = t->tmpS.num_tuples = ns; //ntupS;
                                t->relS.tuples = (*args->skewtask)->tmpS.tuples + outputS[i] //;
                                                + k*(ntupS/args->nthreads);
                                t->tmpS.tuples = (*args->skewtask)->relS.tuples + outputS[i] //;
                                                + k*(ntupS/args->nthreads);

                                task_queue_add<KeyType, PayloadType, TaskType>(part_queue, t);
                            }
                        } 
                        else
                        if(ntupR > 0 && ntupS > 0) {
                            TaskType * t = task_queue_get_slot<KeyType, PayloadType, TaskType>(join_queue);

                            t->relR.num_tuples = t->tmpR.num_tuples = ntupR;
                            t->relR.tuples = (*args->skewtask)->tmpR.tuples + outputR[i];
                            t->tmpR.tuples = (*args->skewtask)->relR.tuples + outputR[i];

                            t->relS.num_tuples = t->tmpS.num_tuples = ntupS;
                            t->relS.tuples = (*args->skewtask)->tmpS.tuples + outputS[i];
                            t->tmpS.tuples = (*args->skewtask)->relS.tuples + outputS[i];

                            task_queue_add<KeyType, PayloadType, TaskType>(join_queue, t);

                            //DEBUGMSG(1, "Join added = R: %d, S: %d\n", 
                            //       t->relR.num_tuples, t->relS.num_tuples);
                        }
                    }

                }
            }

            // add large join tasks in part_queue to the front of the join queue 
            if(my_tid == 0) {
                while((task = task_queue_get_atomic<KeyType, PayloadType, TaskType>(part_queue)))
                    task_queue_add<KeyType, PayloadType, TaskType>(join_queue, task);
            }

        #endif
        
            free(outputR);
            free(outputS);

            BARRIER_ARRIVE(args->barrier, rv);

            if(args->my_tid == 0){
                gettimeofday(&args->partition_end_time, NULL);

            #ifndef DEVELOPMENT_MODE
                //args->e_start_to_partition.stopCounters();
                //args->e_partition_to_end.startCounters();
            #endif
                deltaT = (args->partition_end_time.tv_sec - args->start_time.tv_sec) * 1000000 + args->partition_end_time.tv_usec - args->start_time.tv_usec;
                printf("---- %5s Partition costs time (ms) = %10.4lf\n", pj_partition_pfun[fid].fun_name, deltaT * 1.0 / 1000);
            }

        
            if(!((fid == (pj_partition_pf_num - 1)) && (rp == (RUN_NUMS - 1)))){

                //free(outputR);
                //free(outputS);

                outputR = (int64_t *) calloc((fanOut+1), sizeof(int64_t));
                outputS = (int64_t *) calloc((fanOut+1), sizeof(int64_t));

                //free(args->histR[my_tid]);
                //free(args->histS[my_tid]);

                args->histR[my_tid] = (int32_t *) calloc(fanOut, sizeof(int32_t));
                args->histS[my_tid] = (int32_t *) calloc(fanOut, sizeof(int32_t));

                BARRIER_ARRIVE(args->barrier, rv);
            } 

        }
    }

    BARRIER_ARRIVE(args->barrier, rv);
    
    //DEBUGMSG((my_tid == 0), "Number of join tasks = %d\n", join_queue->count);

    //Join phase
    uint32_t numR_from_build = 0;

    for (int fid = 0; fid < pj_build_pf_num; ++fid) 
    {
        for (int rp = 0; rp < RUN_NUMS; ++rp) 
        {
            BARRIER_ARRIVE(args->barrier, rv);

            if(args->my_tid == 0){
                gettimeofday(&args->partition_end_time, NULL);
            }

            while((task = task_queue_get_atomic<KeyType, PayloadType, TaskType>(join_queue))){
                // do the actual join. join method differs for different algorithms,
                // i.e. bucket chaining, histogram-based, histogram-based with simd &
                // prefetching
                ETHBucketChainingBuild build_output; 
                //printf("thread %d start build \n", my_tid);       
                numR_from_build = pj_build_pfun[fid].fun_ptr(&build_output, &task->relR, NULL);    

                //printf("thread %d start probe \n", my_tid); 
                results += pj_probe_pfun[fid].fun_ptr(&task->relR, numR_from_build, &task->relS, &build_output);     
                //printf("thread %d end probe curr_results %ld parts_processed %d\n", my_tid, results, args->parts_processed); 
                args->parts_processed ++;
            }

            args->result = results;


            BARRIER_ARRIVE(args->barrier, rv);
            // probe phase finished, thread-0 checkpoints the time
            if(args->my_tid == 0){
                gettimeofday(&args->end_time, NULL);

                deltaT = (args->end_time.tv_sec - args->partition_end_time.tv_sec) * 1000000 + args->end_time.tv_usec - args->partition_end_time.tv_usec;
                printf("---- %5s Join costs time (ms) = %10.4lf\n", pj_probe_pfun[fid].fun_name, deltaT * 1.0 / 1000);
            }
        }
    }

    return 0;
}


void initialize_pj_join_thread_args(Relation<KeyType, PayloadType> * rel_r, 
                                    Relation<KeyType, PayloadType> * rel_s,
                                    learned_sort_for_sort_merge::RMI<KeyType, PayloadType> * rmi,
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
                                    vector<double>* slopes, vector<double>* intercepts,
                                    int32_t ** histR, int32_t ** histS,
                                    Tuple<KeyType, PayloadType> * tmpRelR, Tuple<KeyType, PayloadType> * tmpRelS, TaskQueue<KeyType, PayloadType, TaskType> ** part_queue_ptr, TaskQueue<KeyType, PayloadType, TaskType> ** join_queue_ptr,
                                #if SKEW_HANDLING
                                    TaskQueue<KeyType, PayloadType, TaskType> * skew_queue,
                                #endif
                                    pthread_barrier_t* barrier_ptr,
                                    Result * joinresult,
                                    ETHRadixJoinThread<KeyType, PayloadType, TaskType> * args){
    int i;
    int32_t numperthr[2];
    unsigned int SAMPLE_SZ_Rthr, SAMPLE_SZ_Sthr;
    SAMPLE_SZ_Rthr = SAMPLE_SZ_R / NUM_THREADS;
    SAMPLE_SZ_Sthr = SAMPLE_SZ_S / NUM_THREADS;

    #ifdef DEVELOPMENT_MODE
    int numnuma = get_num_numa_regions_develop();
    #else
    int numnuma = get_num_numa_regions();
    #endif

    #if SKEW_HANDLING
        TaskType * skewtask = NULL;
        skew_queue = task_queue_init<KeyType, PayloadType, TaskType>(FANOUT_PASS1);
    #endif

    for(i = 0; i < numnuma; i++){
        part_queue_ptr[i] = task_queue_init<KeyType, PayloadType, TaskType>(FANOUT_PASS1);
        join_queue_ptr[i] = task_queue_init<KeyType, PayloadType, TaskType>((1<<NUM_RADIX_BITS));
    }

    /* allocate temporary space for partitioning */
    tmpRelR = (Tuple<KeyType, PayloadType>*) alloc_aligned(rel_r->num_tuples * sizeof(Tuple<KeyType, PayloadType>) +
                                    relation_padding);
    tmpRelS = (Tuple<KeyType, PayloadType>*) alloc_aligned(rel_s->num_tuples * sizeof(Tuple<KeyType, PayloadType>) +
                                    relation_padding);
    MALLOC_CHECK((tmpRelR && tmpRelS));

    /* allocate histograms arrays, actual allocation is local to threads */
    histR = (int32_t**) alloc_aligned(NUM_THREADS * sizeof(int32_t*));
    histS = (int32_t**) alloc_aligned(NUM_THREADS * sizeof(int32_t*));
    MALLOC_CHECK((histR && histS));

    /* first assign chunks of relR & relS for each thread */
    numperthr[0] = rel_r->num_tuples / NUM_THREADS;
    numperthr[1] = rel_s->num_tuples / NUM_THREADS;

    for(i = 0; i < NUM_THREADS; i++){
        (*(args + i)).relR = rel_r->tuples + i * numperthr[0];
        (*(args + i)).tmpR = tmpRelR;
        (*(args + i)).histR = histR;

        (*(args + i)).relS = rel_s->tuples + i * numperthr[1];
        (*(args + i)).tmpS = tmpRelS;
        (*(args + i)).histS = histS;

        (*(args + i)).numR = (i == (NUM_THREADS-1)) ? 
            (rel_r->num_tuples - i * numperthr[0]) : numperthr[0];
        (*(args + i)).numS = (i == (NUM_THREADS-1)) ? 
            (rel_s->num_tuples - i * numperthr[1]) : numperthr[1];
        (*(args + i)).totalR = rel_r->num_tuples;
        (*(args + i)).totalS = rel_s->num_tuples;

        (*(args + i)).my_tid = i;
        (*(args + i)).part_queue = part_queue_ptr;
        (*(args + i)).join_queue = join_queue_ptr;

        /**** start stuff for learning RMI models ****/
        (*(args + i)).rmi = rmi;
        (*(args + i)).p = p;
        (*(args + i)).original_relR = rel_r;
        (*(args + i)).original_relS = rel_s;
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
        (*(args + i)).intercepts = intercepts;
        /**** end stuff for learning RMI models ****/

#if SKEW_HANDLING
        (*(args + i)).skew_queue = skew_queue;
        (*(args + i)).skewtask   = &skewtask;
#endif
        (*(args + i)).barrier       = barrier_ptr;
        (*(args + i)).nthreads      = NUM_THREADS;
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

    ETHRadixJoinThread<KeyType, PayloadType, TaskType> args[NUM_THREADS];
    ETHRadixJoinThread<KeyType, PayloadType, TaskType> * args_ptr = args;
    

    #ifdef DEVELOPMENT_MODE
    int numnuma = get_num_numa_regions_develop();
    #else
    int numnuma = get_num_numa_regions();
    #endif

    TaskQueue<KeyType, PayloadType, TaskType> * part_queue[numnuma];
    TaskQueue<KeyType, PayloadType, TaskType> * join_queue[numnuma];
    part_queue_ptr = part_queue;
    join_queue_ptr = join_queue;

    rv = pthread_barrier_init(&barrier, NULL, NUM_THREADS);
    if(rv != 0){
        printf("[ERROR] Couldn't create the barrier\n");
        exit(EXIT_FAILURE);
    }

    pthread_attr_init(&attr);

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

    initialize_pj_join_thread_args(&rel_r, &rel_s, &rmi, rmi_params,
                                 SAMPLE_SZ_R, SAMPLE_SZ_S,
                                 tmp_training_sample_in, sorted_training_sample_in, r_tmp_training_sample_in,
                                 r_sorted_training_sample_in, s_tmp_training_sample_in, s_sorted_training_sample_in,
                                 &training_data, sample_count, sample_count_R, sample_count_S,
                                 slopes, intercepts,
                                 histR, histS, tmpRelR, tmpRelS,
                                 part_queue_ptr, join_queue_ptr, 
                        #if SKEW_HANDLING         
                                 skew_queue,
                        #endif    
                                 &barrier, joinresult, args_ptr);

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

        rv = pthread_create(&tid[i], &attr, pj_join_thread, (void*)&args[i]);
        if (rv){
            printf("[ERROR] return code from pthread_create() is %d\n", rv);
            exit(-1);
        }
    }

    // wait for threads to finish
    for(i = 0; i < NUM_THREADS; i++){
        pthread_join(tid[i], NULL);
        result += args[i].result;
    }        

    printf("join results: %ld \n", result);

    return 0;
}



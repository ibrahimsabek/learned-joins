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

#define KeyType RELATION_KEY_TYPE
#define PayloadType RELATION_PAYLOAD_TYPE
#define TaskType Task<RELATION_KEY_TYPE, RELATION_PAYLOAD_TYPE>
#define NUM_THREADS NUM_THREADS_FOR_EVALUATION

#define RUN_NUMS 10 
#define NPJ_MORSE_SIZE 0 //100000

#define SINGLE_TUPLE_PER_BUCKET

using namespace std;

typedef struct StateSIMDForETHNPJ StateSIMDForETHNPJ;
struct StateSIMDForETHNPJ {
  __m512i key;
  __m512i tb_off;
  __m512i ht_off;
  __mmask8 m_have_tuple;
  char stage;
};

volatile static char npj_g_lock = 0, npj_g_lock_morse = 0;
volatile static uint64_t npj_total_num = 0, npj_global_curse = 0, npj_global_upper, npj_global_morse_size;
    
typedef void (*NPJBuildFun)(ETHNonPartitionJoinBuild<KeyType, PayloadType> *build_input, Relation<KeyType, PayloadType> * rel_r_partition, Relation<KeyType, PayloadType> * tmp_r);
volatile static struct Fun {
  NPJBuildFun fun_ptr;
  char fun_name[8];
} npj_pfun[2];
volatile static int npj_pf_num = 0;

//TOOD: to put the morse-driven implementation here !!!!!
void npj_build_rel_r_partition(ETHNonPartitionJoinBuild<KeyType, PayloadType> *build_input, Relation<KeyType, PayloadType> * rel_r_partition, Relation<KeyType, PayloadType> * tmp_r)
{
    Hashtable<KeyType, PayloadType>* ht = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_input)->ht;  
    BucketBuffer<KeyType, PayloadType>** overflowbuf = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_input)->overflowbuf;
#ifdef DEVELOPMENT_MODE
     /*unordered_map<uint64_t, uint64_t>* build_visits_map = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_input)->build_hash_bucket_visits;        
     vector<KeyType>* build_keys_list = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_input)->build_keys_list;  
     vector<uint64_t>* build_keys_hash_list = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_input)->build_keys_hash_list;        
     volatile char * keys_hash_latch = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_input)->keys_hash_latch;*/        
#endif

    uint64_t i;
#ifndef USE_MURMUR3_HASH
    const uint64_t hashmask = ht->hash_mask;
    const uint64_t skipbits = ht->skip_bits;
#endif
#ifdef PREFETCH_NPJ
    size_t prefetch_index = PREFETCH_DISTANCE;
#endif
    
    for(i=0; i < rel_r_partition->num_tuples; i++){
        Tuple<KeyType, PayloadType> * dest;
        Bucket<KeyType, PayloadType> * curr, * nxt;

#ifdef PREFETCH_NPJ
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

#ifdef DEVELOPMENT_MODE
        //(*build_visits_map)[idx]++;

        //lock(keys_hash_latch);
        //(*build_keys_list).push_back(rel_r_partition->tuples[i].key);
        //(*build_keys_hash_list).push_back(idx);
        //unlock(keys_hash_latch);
#endif

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

void npj_build_rel_r_partition_imv(ETHNonPartitionJoinBuild<KeyType, PayloadType> *build_input, Relation<KeyType, PayloadType> * rel_r_partition, Relation<KeyType, PayloadType> * tmp_r)
{
    Hashtable<KeyType, PayloadType>* ht = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_input)->ht;  
    BucketBuffer<KeyType, PayloadType>** overflowbuf = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_input)->overflowbuf;

    int32_t k = 0, done = 0, num, num_temp;
    __attribute__((aligned(CACHE_LINE_SIZE)))  __mmask8 mask[NPJ_VECTOR_SCALE + 1], m_to_insert = 0, m_no_conflict;
    __m512i v_offset = _mm512_set1_epi64(0), v_base_offset_upper = _mm512_set1_epi64(rel_r_partition->num_tuples * sizeof(Tuple<KeyType, PayloadType>)), v_base_offset, v_ht_cell, 
      v_cell_hash, v_neg_one512 = _mm512_set1_epi64(-1), v_zero512 = _mm512_set1_epi64(0), v_ht_addr = _mm512_set1_epi64((uint64_t) ht->buckets), 
      v_bucket_size = _mm512_set1_epi64(sizeof(Bucket<KeyType, PayloadType>)), v_next_off = _mm512_set1_epi64(8), v_addr, v_all_ones = _mm512_set1_epi64(-1),
      v_conflict, v_new_bucket, v_next, v_key_off = _mm512_set1_epi64(16), v_count_off = _mm512_set1_epi64(4);

    __m256i v256_one = _mm256_set1_epi32(1);
    uint64_t *new_bucket = (uint64_t*) &v_new_bucket;
    Bucket<KeyType, PayloadType> * bucket;
    __attribute__((aligned(CACHE_LINE_SIZE))) uint64_t cur_offset = 0, base_off[NPJ_MAX_VECTOR_SCALE], *ht_pos;

    #ifndef USE_MURMUR3_HASH
    __m512i v_factor = _mm512_set1_epi64(ht->hash_mask), v_shift = _mm512_set1_epi64(ht->skip_bits); 
    #endif

    for (int i = 0; i <= NPJ_VECTOR_SCALE; ++i) 
    {
        base_off[i] = i * sizeof(Tuple<KeyType, PayloadType>);
        mask[i] = (1 << i) - 1;
    }
    v_base_offset = _mm512_load_epi64(base_off);
    __attribute__((aligned(CACHE_LINE_SIZE)))   StateSIMDForETHNPJ state[NPJ_SIMDStateSize + 1];

    // init # of the state
    for (int i = 0; i <= NPJ_SIMDStateSize; ++i) {
        state[i].stage = 1;
        state[i].m_have_tuple = 0;
        state[i].ht_off = _mm512_set1_epi64(0);
        state[i].key = _mm512_set1_epi64(0);
    }

    for (uint64_t cur = 0; 1;) 
    {
        k = (k >= NPJ_SIMDStateSize) ? 0 : k;
        if ((cur >= rel_r_partition->num_tuples)) {
            if (state[k].m_have_tuple == 0 && state[k].stage != 3) 
            {
                ++done;
                state[k].stage = 3;
                ++k;
                continue;
            }
            if ((done >= NPJ_SIMDStateSize)) 
            {
                if (state[NPJ_SIMDStateSize].m_have_tuple > 0) 
                {
                    k = NPJ_SIMDStateSize;
                    state[NPJ_SIMDStateSize].stage = 0;
                } else {
                    break;
                }
            }
        }

        switch (state[k].stage) 
        {
            case 1: {
                _mm_prefetch((char *)(((void *)rel_r_partition->tuples) + cur_offset + NPJ_PDIS), _MM_HINT_T0);
                _mm_prefetch((char *)(((void *)rel_r_partition->tuples) + cur_offset + NPJ_PDIS + CACHE_LINE_SIZE), _MM_HINT_T0);
                _mm_prefetch((char *)(((void *)rel_r_partition->tuples) + cur_offset + NPJ_PDIS + 2 * CACHE_LINE_SIZE), _MM_HINT_T0);

                v_offset = _mm512_add_epi64(_mm512_set1_epi64(cur_offset), v_base_offset);
                // count the number of empty tuples
                cur_offset = cur_offset + base_off[NPJ_VECTOR_SCALE];
                state[k].m_have_tuple = _mm512_cmpgt_epi64_mask(v_base_offset_upper, v_offset);
                cur = cur + NPJ_VECTOR_SCALE;
                
                state[k].key = _mm512_mask_i64gather_epi64(state[k].key, state[k].m_have_tuple, v_offset, ((void * )rel_r_partition->tuples), 1);
                

                // Perform hashing for the keys
            #ifndef USE_MURMUR3_HASH
                v_cell_hash = _mm512_and_epi64(state[k].key, v_factor);
                v_cell_hash = _mm512_srlv_epi64(v_cell_hash, v_shift);
                v_cell_hash = _mm512_mullo_epi64(v_cell_hash, v_bucket_size);
            #else
                //TODO: murmur3 hashing to be implemented
                printf("Murmur3 hashing is not implemented in the probe phase yet!! \n");
            #endif

                state[k].ht_off = _mm512_mask_add_epi64(state[k].ht_off, state[k].m_have_tuple, v_cell_hash, v_ht_addr);
                state[k].stage = 2;

                ht_pos = (uint64_t *) &state[k].ht_off;
                for (int i = 0; i < NPJ_VECTOR_SCALE; ++i) 
                {
                    _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
                }

            }
            break;

            case 2: 
            {
                v_ht_cell = _mm512_mask_i64gather_epi64(v_neg_one512, state[k].m_have_tuple, state[k].ht_off, 0, 1);
                // inset new nodes
                m_to_insert = _mm512_cmpeq_epi64_mask(v_ht_cell, v_zero512);
                m_to_insert = _mm512_kand(m_to_insert, state[k].m_have_tuple);
                if (m_to_insert == 0) 
                {
                    state[k].stage = 0;
                    --k;
                    break;
                }
                v_addr = _mm512_mask_blend_epi64(m_to_insert, v_all_ones, state[k].ht_off);
                v_conflict = _mm512_conflict_epi64(v_addr);
                m_no_conflict = _mm512_testn_epi64_mask(v_conflict, v_all_ones);
                m_no_conflict = _mm512_kand(m_no_conflict, m_to_insert);

                if (m_no_conflict) 
                {
                    // write the key, payload, count, next to the nodes
                    _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_addr, v_key_off), state[k].key, 1);
                    _mm512_mask_i64scatter_epi32(0, m_no_conflict, _mm512_add_epi64(v_addr, v_count_off), v256_one, 1);
                    state[k].m_have_tuple = _mm512_kandn(m_no_conflict, state[k].m_have_tuple);
                    //found += _mm_popcnt_u32(m_no_conflict);
                }

                num = _mm_popcnt_u32(state[k].m_have_tuple);
                if (num == NPJ_VECTOR_SCALE || done >= NPJ_SIMDStateSize) 
                {
                    state[k].stage = 0;
                    --k;
                } 
                else if (num == 0) 
                {
                    state[k].stage = 1;
                    --k;
                } 
                else 
                {
                    if ((done < NPJ_SIMDStateSize)) 
                    {
                        num_temp = _mm_popcnt_u32(state[NPJ_SIMDStateSize].m_have_tuple);
                        if (num + num_temp < NPJ_VECTOR_SCALE) 
                        {
                            // compress v
                            state[k].ht_off = _mm512_maskz_compress_epi64(state[k].m_have_tuple, state[k].ht_off);
                            state[k].key = _mm512_maskz_compress_epi64(state[k].m_have_tuple, state[k].key);
                            // expand v -> temp
                            state[NPJ_SIMDStateSize].ht_off = _mm512_mask_expand_epi64(state[NPJ_SIMDStateSize].ht_off, _mm512_knot(state[NPJ_SIMDStateSize].m_have_tuple), state[k].ht_off);
                            state[NPJ_SIMDStateSize].key = _mm512_mask_expand_epi64(state[NPJ_SIMDStateSize].key, _mm512_knot(state[NPJ_SIMDStateSize].m_have_tuple), state[k].key);
                            state[NPJ_SIMDStateSize].m_have_tuple = mask[num + num_temp];
                            state[k].m_have_tuple = 0;
                            state[k].stage = 1;
                            --k;
                        } 
                        else 
                        {
                            // expand temp -> v
                            state[k].ht_off = _mm512_mask_expand_epi64(state[k].ht_off, _mm512_knot(state[k].m_have_tuple), state[NPJ_SIMDStateSize].ht_off);
                            state[k].key = _mm512_mask_expand_epi64(state[k].key, _mm512_knot(state[k].m_have_tuple), state[NPJ_SIMDStateSize].key);
                            // compress temp
                            state[NPJ_SIMDStateSize].m_have_tuple = _mm512_kand(state[NPJ_SIMDStateSize].m_have_tuple, _mm512_knot(mask[NPJ_VECTOR_SCALE - num]));
                            state[NPJ_SIMDStateSize].ht_off = _mm512_maskz_compress_epi64(state[NPJ_SIMDStateSize].m_have_tuple, state[NPJ_SIMDStateSize].ht_off);
                            state[NPJ_SIMDStateSize].key = _mm512_maskz_compress_epi64(state[NPJ_SIMDStateSize].m_have_tuple, state[NPJ_SIMDStateSize].key);
                            state[k].m_have_tuple = mask[NPJ_VECTOR_SCALE];
                            state[NPJ_SIMDStateSize].m_have_tuple = (state[NPJ_SIMDStateSize].m_have_tuple >> (NPJ_VECTOR_SCALE - num));
                            state[k].stage = 0;
                            --k;
                        }
                    }
                }
            }
            break;

            case 0: 
            {
                // insert new buckets after valid buckets
                m_to_insert = state[k].m_have_tuple;
                v_addr = _mm512_mask_blend_epi64(m_to_insert, v_all_ones, state[k].ht_off);
                v_conflict = _mm512_conflict_epi64(v_addr);
                m_no_conflict = _mm512_testn_epi64_mask(v_conflict, v_all_ones);
                m_no_conflict = _mm512_kand(m_no_conflict, m_to_insert);
                
                for (int i = 0; i < NPJ_VECTOR_SCALE; ++i) 
                {
                    new_bucket[i] = 0;
                    if (m_no_conflict & (1 << i)) 
                    {
                        get_new_bucket(&bucket, overflowbuf);
                        new_bucket[i] = (uint64_t)bucket;
                    }
                }
                v_next = _mm512_mask_i64gather_epi64(v_all_ones, state[k].m_have_tuple, _mm512_add_epi64(state[k].ht_off, v_next_off), 0, 1);
                _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_new_bucket, v_key_off), state[k].key, 1);
                //_mm512_mask_i64scatter_epi32(0, m_no_conflict, (v_new_bucket), v256_one, 1);
                _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_new_bucket, v_next_off), v_next, 1);
                _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(state[k].ht_off, v_next_off), v_new_bucket, 1);

                //found += _mm_popcnt_u32(m_no_conflict);
                state[k].m_have_tuple = _mm512_kandn(m_no_conflict, state[k].m_have_tuple);

                num = _mm_popcnt_u32(state[k].m_have_tuple);

                if (num == NPJ_VECTOR_SCALE || done >= NPJ_SIMDStateSize) 
                {
                    ht_pos = (uint64_t *) &state[k].ht_off;
                    for (int i = 0; i < NPJ_VECTOR_SCALE; ++i) {
                        _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
                        //_mm_prefetch((char * )(ht_pos[i] + CACHE_LINE_SIZE), _MM_HINT_T0);
                    }
                } 
                else if (num == 0) 
                {
                    state[k].stage = 1;
                    --k;
                    break;
                } 
                else
                {
                    if ((done < NPJ_SIMDStateSize)) 
                    {
                        num_temp = _mm_popcnt_u32(state[NPJ_SIMDStateSize].m_have_tuple);
                        if (num + num_temp < NPJ_VECTOR_SCALE) 
                        {
                            // compress v
                            state[k].ht_off = _mm512_maskz_compress_epi64(state[k].m_have_tuple, state[k].ht_off);
                            state[k].key = _mm512_maskz_compress_epi64(state[k].m_have_tuple, state[k].key);
                            // expand v -> temp
                            state[NPJ_SIMDStateSize].ht_off = _mm512_mask_expand_epi64(state[NPJ_SIMDStateSize].ht_off, _mm512_knot(state[NPJ_SIMDStateSize].m_have_tuple), state[k].ht_off);
                            state[NPJ_SIMDStateSize].key = _mm512_mask_expand_epi64(state[NPJ_SIMDStateSize].key, _mm512_knot(state[NPJ_SIMDStateSize].m_have_tuple), state[k].key);
                            state[NPJ_SIMDStateSize].m_have_tuple = mask[num + num_temp];
                            state[k].m_have_tuple = 0;
                            state[k].stage = 1;
                            --k;
                            break;
                        } 
                        else 
                        {
                            // expand temp -> v
                            state[k].ht_off = _mm512_mask_expand_epi64(state[k].ht_off, _mm512_knot(state[k].m_have_tuple), state[NPJ_SIMDStateSize].ht_off);
                            state[k].key = _mm512_mask_expand_epi64(state[k].key, _mm512_knot(state[k].m_have_tuple), state[NPJ_SIMDStateSize].key);

                            // compress temp
                            state[NPJ_SIMDStateSize].m_have_tuple = _mm512_kand(state[NPJ_SIMDStateSize].m_have_tuple, _mm512_knot(mask[NPJ_VECTOR_SCALE - num]));
                            state[NPJ_SIMDStateSize].ht_off = _mm512_maskz_compress_epi64(state[NPJ_SIMDStateSize].m_have_tuple, state[NPJ_SIMDStateSize].ht_off);
                            state[NPJ_SIMDStateSize].key = _mm512_maskz_compress_epi64(state[NPJ_SIMDStateSize].m_have_tuple, state[NPJ_SIMDStateSize].key);
                            state[k].m_have_tuple = mask[NPJ_VECTOR_SCALE];
                            state[NPJ_SIMDStateSize].m_have_tuple = (state[NPJ_SIMDStateSize].m_have_tuple >> (NPJ_VECTOR_SCALE - num));
                            state[k].stage = 0;

                            ht_pos = (uint64_t *) &state[k].ht_off;
                            for (int i = 0; i < NPJ_VECTOR_SCALE; ++i) {
                                _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
                                //_mm_prefetch((char * )(ht_pos[i] + CACHE_LINE_SIZE), _MM_HINT_T0);
                            }
                        }
                    }
                }
            }
            break;
        }
        ++k;   
    }
}

void * npj_join_thread(void * param)
{
    ETHNonPartitionJoinThread<KeyType, PayloadType, TaskType> * args   = (ETHNonPartitionJoinThread<KeyType, PayloadType, TaskType> *) param;
    int rv;   int deltaT = 0;
    BucketBuffer<KeyType, PayloadType> * overflowbuf; // allocate overflow buffer for each thread

    uint32_t nbuckets = (args->relR.num_tuples / BUCKET_SIZE / NUM_THREADS);

    if (args->tid == 0) {
        strcpy(npj_pfun[1].fun_name, "IMV");
        strcpy(npj_pfun[0].fun_name, "Naive");

        npj_pfun[1].fun_ptr = npj_build_rel_r_partition_imv;
        npj_pfun[0].fun_ptr = npj_build_rel_r_partition;

        npj_pf_num = 2;
    }
    BARRIER_ARRIVE(args->barrier, rv);
    
    ETHNonPartitionJoinBuild<KeyType, PayloadType> build_data; 
    for (int fid = 0; fid < npj_pf_num; ++fid) 
    {
        for (int rp = 0; rp < RUN_NUMS; ++rp) 
        {
            init_bucket_buffer(&overflowbuf);
            allocate_hashtable(&args->ht, nbuckets);
            build_data.ht = args->ht;
            build_data.overflowbuf = &overflowbuf;

        #ifdef PERF_COUNTERS
            if(args->tid == 0){
                //TODO: performance counters to be implemented
            }
        #endif

            /* wait at a barrier until each thread starts and start timer */
            BARRIER_ARRIVE(args->barrier, rv);

            /* the first thread checkpoints the start time */
            if(args->tid == 0){
                /* no partitionig phase, but we use partition variables to store building stats */
                gettimeofday(&args->start_time, NULL);
            #ifndef DEVELOPMENT_MODE
                //args->e_start_to_partition.startCounters();
            #endif
            }

        #ifdef DEVELOPMENT_MODE
            /*build_data.build_hash_bucket_visits = args->build_hash_bucket_visits;
            build_data.probe_hash_bucket_visits = args->probe_hash_bucket_visits; 
            build_data.keys_hash_latch = args->keys_hash_latch; 
            build_data.build_keys_list = args->build_keys_list;
            build_data.build_keys_hash_list = args->build_keys_hash_list;        
            build_data.probe_keys_list = args->probe_keys_list;
            build_data.probe_keys_hash_list = args->probe_keys_hash_list;*/        
        #endif        

        #if NPJ_MORSE_SIZE
            morse_driven(param, npj_pfun[fid].fun_ptr, &overflowbuf);
        #else
            npj_pfun[fid].fun_ptr(&build_data, &args->relR, NULL);
        #endif

            /* wait at a barrier until each thread completes build phase */
            BARRIER_ARRIVE(args->barrier, rv);

        #ifdef PERF_COUNTERS
            if(args->tid == 0)
            {
                //TODO: performance counters to be implemented
            }
            /* Just to make sure we get consistent performance numbers */
            BARRIER_ARRIVE(args->barrier, rv);
        #endif

            /* build phase finished, thread-0 checkpoints the time */
            if(args->tid == 0){
                gettimeofday(&args->partition_end_time, NULL);

            #ifndef DEVELOPMENT_MODE
                //args->e_start_to_partition.stopCounters();
                //args->e_partition_to_end.startCounters();
            #endif
                deltaT = (args->partition_end_time.tv_sec - args->start_time.tv_sec) * 1000000 + args->partition_end_time.tv_usec - args->start_time.tv_usec;
                printf("---- %5s BUILD costs time (ms) = %10.4lf\n", npj_pfun[fid].fun_name, deltaT * 1.0 / 1000);
                npj_total_num = 0;
                npj_global_curse = 0;
            }
            
            destroy_hashtable(args->ht);
            free_bucket_buffer(overflowbuf);
            
        }
    }

    //TODO: loop here for the prob function as well
    // probe for matching tuples from the assigned part of relS  
    /*#ifdef NPJ_ETH_AVX_IMV      
        args->num_results = join_steps.template probe_rel_s_partition_imv<ETHNonPartitionJoinBuild<KeyType, PayloadType>>(NULL, &args->relS, &build_data, args->tid); 
    #else
        args->num_results = join_steps.template probe_rel_s_partition<ETHNonPartitionJoinBuild<KeyType, PayloadType>>(NULL, &args->relS, &build_data);
    #endif

        // for a reliable timing we have to wait until all finishes
        BARRIER_ARRIVE(args->barrier, rv);

        // probe phase finished, thread-0 checkpoints the time
        if(args->tid == 0){
            gettimeofday(&args->end_time, NULL);
        #ifndef DEVELOPMENT_MODE
            //args->e_partition_to_end.stopCounters();
        #endif
        }

    #ifdef PERF_COUNTERS
        if(args->tid == 0) 
        {
            //TODO: performance counters to be implemented
        }
        // Just to make sure we get consistent performance numbers
        BARRIER_ARRIVE(args->barrier, rv);
    #endif
    */


    return 0;
}

void initialize_npj_join_thread_args(Relation<KeyType, PayloadType> * rel_r, 
                                 Relation<KeyType, PayloadType> * rel_s,
                                 Hashtable<KeyType, PayloadType> * ht, 
                                 pthread_barrier_t* barrier_ptr,
                                 Result * joinresult,
                                 ETHNonPartitionJoinThread<KeyType, PayloadType, TaskType> * args){
    int i;
    uint64_t numR, numS, numRthr, numSthr; /* total and per thread num */

    numR = rel_r->num_tuples;
    numS = rel_s->num_tuples;
    numRthr = numR / NUM_THREADS;
    numSthr = numS / NUM_THREADS;

#ifdef DEVELOPMENT_MODE
    /*for(uint32_t j = 0; j < nbuckets; j++)
    {
        build_visits_map[j] = 0;
        probe_visits_map[j] = 0;            
    }*/        
#endif

    for(i = 0; i < NUM_THREADS; i++)
    {
        (*(args + i)).tid = i;
        (*(args + i)).ht = ht; 
#ifdef DEVELOPMENT_MODE
    /*    (*(args + i)).build_hash_bucket_visits = &build_visits_map;
        (*(args + i)).probe_hash_bucket_visits = &probe_visits_map;
        (*(args + i)).keys_hash_latch = &keys_hash_latch;            
        (*(args + i)).build_keys_list = &build_keys_list;
        (*(args + i)).build_keys_hash_list = &build_keys_hash_list;
        (*(args + i)).probe_keys_list = &probe_keys_list;
        (*(args + i)).probe_keys_hash_list = &probe_keys_hash_list;*/
#endif
        /* assing part of the relR for next thread */
        (*(args + i)).relR.num_tuples = (i == (NUM_THREADS-1)) ? numR : numRthr;
        (*(args + i)).relR.tuples = rel_r->tuples + numRthr * i;
        numR -= numRthr;

        /* assing part of the relS for next thread */
        (*(args + i)).relS.num_tuples = (i == (NUM_THREADS-1)) ? numS : numSthr;
        (*(args + i)).relS.tuples = rel_s->tuples + numSthr * i;
        numS -= numSthr;

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

    ETHNonPartitionJoinThread<KeyType, PayloadType, TaskType> args[NUM_THREADS];
    ETHNonPartitionJoinThread<KeyType, PayloadType, TaskType> * args_ptr = args;

    rv = pthread_barrier_init(&barrier, NULL, NUM_THREADS);
    if(rv != 0){
        printf("[ERROR] Couldn't create the barrier\n");
        exit(EXIT_FAILURE);
    }

    pthread_attr_init(&attr);

    Hashtable<KeyType, PayloadType> * ht;
#ifdef INPUT_HASH_TABLE_SIZE       
    uint32_t nbuckets = hash_table_size;
#else
    uint32_t nbuckets = (rel_r.num_tuples / BUCKET_SIZE / NUM_THREADS);
#endif        
    allocate_hashtable(&ht, nbuckets);

    initialize_npj_join_thread_args(&rel_r, &rel_s, ht, &barrier, joinresult, args_ptr);

    npj_global_curse = 0;
    npj_global_upper = rel_r.num_tuples;
    if(NUM_THREADS==1){
        npj_global_morse_size= rel_r.num_tuples;
    }else{
        npj_global_morse_size = NPJ_MORSE_SIZE;
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

        rv = pthread_create(&tid[i], &attr, npj_join_thread, (void*)&args[i]);
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


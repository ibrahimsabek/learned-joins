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
#define NPJ_MORSE_SIZE 0 //100000

#define SINGLE_TUPLE_PER_BUCKET

using namespace std;
using namespace learned_sort_for_sort_merge;

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

typedef uint64_t (*NPJProbeFun)(Relation<KeyType, PayloadType> * rel_r_partition, Relation<KeyType, PayloadType> * rel_s_partition, ETHNonPartitionJoinBuild<KeyType, PayloadType> *build_output);
volatile static struct Fun1 {
  NPJProbeFun fun_ptr;
  char fun_name[8];
} npj_pfun1[2];

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
      v_conflict, v_new_bucket, v_next, v_key_off = _mm512_set1_epi64(16), v_count_off = _mm512_set1_epi64(4), 
      c1 = _mm512_set1_epi64(0x85ebca6b), c2 = _mm512_set1_epi64(0xc2b2ae35), s1, x1, s2, x2, s3, 
      v_factor = _mm512_set1_epi64(ht->hash_mask), v_shift = _mm512_set1_epi64(ht->skip_bits);

    __m256i v256_one = _mm256_set1_epi32(1);
    uint64_t *new_bucket = (uint64_t*) &v_new_bucket;
    Bucket<KeyType, PayloadType> * bucket;
    __attribute__((aligned(CACHE_LINE_SIZE))) uint64_t cur_offset = 0, base_off[NPJ_MAX_VECTOR_SCALE], *ht_pos;


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
                s1 = _mm512_srli_epi64(state[k].key, 16);
                x1 = _mm512_xor_epi64(state[k].key, s1);
                s1 = _mm512_mullo_epi64(x1, c1); 

                s2 = _mm512_srli_epi64(s1, 13);
                x2 = _mm512_xor_epi64(s1, s2);
                s2 = _mm512_mullo_epi64(x2, c2); 

                s3 = _mm512_srli_epi64(s2, 16);
                v_cell_hash = _mm512_xor_epi64(s2, s3);

                v_cell_hash = _mm512_and_epi64(v_cell_hash, v_factor);
                v_cell_hash = _mm512_srlv_epi64(v_cell_hash, v_shift);
                v_cell_hash = _mm512_mullo_epi64(v_cell_hash, v_bucket_size);           
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

uint64_t npj_probe_rel_s_partition(Relation<KeyType, PayloadType> * rel_r_partition, Relation<KeyType, PayloadType> * rel_s_partition, ETHNonPartitionJoinBuild<KeyType, PayloadType> *build_output)
{
    Hashtable<KeyType, PayloadType>* ht = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_output)->ht;  
#ifdef DEVELOPMENT_MODE
    /*unordered_map<uint64_t, uint64_t>* probe_visits_map = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_output)->probe_hash_bucket_visits;        
    vector<KeyType>* probe_keys_list = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_output)->probe_keys_list;  
    vector<uint64_t>* probe_keys_hash_list = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_output)->probe_keys_hash_list;        
    volatile char* keys_hash_latch = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_output)->keys_hash_latch;*/       
#endif
    
    uint64_t i, j;
    uint64_t matches;

#ifndef USE_MURMUR3_HASH
    const uint64_t hashmask = ht->hash_mask;
    const uint64_t skipbits = ht->skip_bits;
#endif
#ifdef PREFETCH_NPJ    
    size_t prefetch_index = PREFETCH_DISTANCE;
#endif
    
    matches = 0;

    for (i = 0; i < rel_s_partition->num_tuples; i++)
    {
#ifdef PREFETCH_NPJ        
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

#ifdef DEVELOPMENT_MODE
        /*lock(&b->latch);
        (*probe_visits_map)[idx]++;
        unlock(&b->latch);*/

        //lock(keys_hash_latch);
        //(*probe_keys_list).push_back(rel_s_partition->tuples[i].key);
        //(*probe_keys_hash_list).push_back(idx);
        //unlock(keys_hash_latch);
#endif

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

uint64_t npj_probe_rel_s_partition_imv(Relation<KeyType, PayloadType> * rel_r_partition, Relation<KeyType, PayloadType> * rel_s_partition, ETHNonPartitionJoinBuild<KeyType, PayloadType> *build_output)
{
    Hashtable<KeyType, PayloadType>* ht = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_output)->ht;  
    uint64_t matches = 0;

    int32_t k = 0, num, num_temp, done = 0, new_add = 0;
    __m512i v_ht_cell, v_base_offset_upper = _mm512_set1_epi64(rel_s_partition->num_tuples * sizeof(Tuple<KeyType, PayloadType>)), v_offset, v_base_offset, v_cell_hash,
            v_bucket_size = _mm512_set1_epi64(sizeof(Bucket<KeyType, PayloadType>)), v_ht_addr = _mm512_set1_epi64((uint64_t)ht->buckets), v_neg_one512 = _mm512_set1_epi64(-1),
            v_zero512 = _mm512_set1_epi64(0), v_tuple_size = _mm512_set1_epi64(16), v_next_off = _mm512_set1_epi64(8),
            c1 = _mm512_set1_epi64(0x85ebca6b), c2 = _mm512_set1_epi64(0xc2b2ae35), s1, x1, s2, x2, s3, 
            v_factor = _mm512_set1_epi64(ht->hash_mask), v_shift = _mm512_set1_epi64(ht->skip_bits);
    __attribute__((aligned(CACHE_LINE_SIZE))) __mmask8 mask[NPJ_VECTOR_SCALE + 1], m_valid_bucket = 0, m_match = 0;
    __attribute__((aligned(CACHE_LINE_SIZE))) uint64_t cur_offset = 0, base_off[NPJ_MAX_VECTOR_SCALE], *ht_pos;

    for (int i = 0; i <= NPJ_VECTOR_SCALE; ++i) 
    {
        base_off[i] = i * sizeof(Tuple<KeyType, PayloadType>);
        mask[i] = (1 << i) - 1;
    }
    v_base_offset = _mm512_load_epi64(base_off);

    __attribute__((aligned(CACHE_LINE_SIZE)))  StateSIMDForETHNPJ state[NPJ_SIMDStateSize + 1];
    // init # of the state
    for (int i = 0; i <= NPJ_SIMDStateSize; ++i) {
        state[i].stage = 1;
        state[i].m_have_tuple = 0;
    } 
   

    for (uint64_t cur = 0; 1;) 
    {
        k = (k >= NPJ_SIMDStateSize) ? 0 : k;

        if (cur >= rel_s_partition->num_tuples) {
            if (state[k].m_have_tuple == 0 && state[k].stage != 3) {
                ++done;
                state[k].stage = 3;
                ++k;
                continue;
            }
            if ((done >= NPJ_SIMDStateSize)) {
                if (state[NPJ_SIMDStateSize].m_have_tuple > 0) {
                    k = NPJ_SIMDStateSize;
                    state[NPJ_SIMDStateSize].stage = 0;
                } else {
                    break;
                }
            }
        }
        
        switch (state[k].stage) 
        {
            case 1: 
            {
                _mm_prefetch((char *)(((void *)rel_s_partition->tuples) + cur_offset + NPJ_PDIS), _MM_HINT_T0);
                _mm_prefetch((char *)(((void *)rel_s_partition->tuples) + cur_offset + NPJ_PDIS + CACHE_LINE_SIZE), _MM_HINT_T0);
                //_mm_prefetch((char *)(((void *)rel_s_partition->tuples) + cur_offset + NPJ_PDIS + 2 * CACHE_LINE_SIZE), _MM_HINT_T0);

                v_offset = _mm512_add_epi64(_mm512_set1_epi64(cur_offset), v_base_offset);
                // count the number of empty tuples
                cur_offset = cur_offset + base_off[NPJ_VECTOR_SCALE];
                cur = cur + NPJ_VECTOR_SCALE;

                state[k].m_have_tuple = _mm512_cmpgt_epi64_mask(v_base_offset_upper, v_offset);
                state[k].key = _mm512_mask_i64gather_epi64(state[k].key, state[k].m_have_tuple, v_offset, ((void *)rel_s_partition->tuples), 1);
                
                // Perform hashing for the keys
            #ifndef USE_MURMUR3_HASH
                v_cell_hash = _mm512_and_epi64(state[k].key, v_factor);
                v_cell_hash = _mm512_srlv_epi64(v_cell_hash, v_shift);
                v_cell_hash = _mm512_mullo_epi64(v_cell_hash, v_bucket_size);
            #else
                s1 = _mm512_srli_epi64(state[k].key, 16);
                x1 = _mm512_xor_epi64(state[k].key, s1);
                s1 = _mm512_mullo_epi64(x1, c1); 

                s2 = _mm512_srli_epi64(s1, 13);
                x2 = _mm512_xor_epi64(s1, s2);
                s2 = _mm512_mullo_epi64(x2, c2); 

                s3 = _mm512_srli_epi64(s2, 16);
                v_cell_hash = _mm512_xor_epi64(s2, s3);

                v_cell_hash = _mm512_and_epi64(v_cell_hash, v_factor);
                v_cell_hash = _mm512_srlv_epi64(v_cell_hash, v_shift);
                v_cell_hash = _mm512_mullo_epi64(v_cell_hash, v_bucket_size);  
            #endif            
                state[k].ht_off = _mm512_mask_add_epi64(state[k].ht_off, state[k].m_have_tuple, v_cell_hash, v_ht_addr);
                state[k].stage = 2;

                ht_pos = (uint64_t *) &state[k].ht_off;
                for (int i = 0; i < NPJ_VECTOR_SCALE; ++i) {
                    _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
                }
            }
            break;

            case 2: 
            {
                v_ht_cell = _mm512_mask_i64gather_epi64(v_neg_one512, state[k].m_have_tuple, state[k].ht_off, 0, 1);
                m_valid_bucket = _mm512_cmpneq_epi64_mask(v_ht_cell, v_zero512);
                state[k].m_have_tuple = _mm512_kand(m_valid_bucket, state[k].m_have_tuple);
                state[k].stage = 0;
                num = _mm_popcnt_u32(state[k].m_have_tuple);

                if (num == NPJ_VECTOR_SCALE) 
                {
                    ht_pos = (uint64_t *) &state[k].ht_off;
                    for (int i = 0; i < NPJ_VECTOR_SCALE; ++i) {
                        _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
                    }

                } else if(num ==0) 
                {
                    state[k].stage = 1;
                    --k;
                    break;
                } 
                else
                {
                     if (done < NPJ_SIMDStateSize) {
                        num_temp = _mm_popcnt_u32(state[NPJ_SIMDStateSize].m_have_tuple);
                        if (num + num_temp < NPJ_VECTOR_SCALE) {
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
                        } else {
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
                            }
                        }
                    }
                }
            }
            break;

            case 0: 
            {
                v_ht_cell = _mm512_mask_i64gather_epi64(v_neg_one512, state[k].m_have_tuple, _mm512_add_epi64(state[k].ht_off, v_tuple_size), 0, 1); 

                //TODO: handle the case of multiple tuples in the same bucket
                m_match = _mm512_cmpeq_epi64_mask(state[k].key, v_ht_cell);
                m_match = _mm512_kand(m_match, state[k].m_have_tuple);
                new_add = _mm_popcnt_u32(m_match);
                matches += new_add;

                // update next
                state[k].ht_off = _mm512_mask_i64gather_epi64(v_zero512, state[k].m_have_tuple, _mm512_add_epi64(state[k].ht_off, v_next_off), 0, 1);
                state[k].m_have_tuple = _mm512_kand(_mm512_cmpneq_epi64_mask(state[k].ht_off, v_zero512), state[k].m_have_tuple);

                num = _mm_popcnt_u32(state[k].m_have_tuple);

                if (num == NPJ_VECTOR_SCALE) {

                    ht_pos = (uint64_t *) &state[k].ht_off;
                    for (int i = 0; i < NPJ_VECTOR_SCALE; ++i) {
                        _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
                    }

                } else if(num ==0) {
                    state[k].stage = 1;
                    --k;
                    break;
                } 
                else
                {
                    if (done < NPJ_SIMDStateSize) 
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
                            }
                
                        }
                    }
                }
            }
            break;

        }
        ++k;
    }

    return matches;
}

#ifdef BUILD_RMI_FROM_TWO_DATASETS
void * sample_and_train_models_threaded(ETHNonPartitionJoinThread<KeyType, PayloadType, TaskType> * args)
{
    int rv;
    int tid = args->tid;

    vector<vector<vector<training_point<KeyType, PayloadType>>>> * training_data = args->training_data;
    static const unsigned int NUM_LAYERS = args->p.arch.size();
    if(tid == 0)
    {
      for (unsigned int layer_idx = 0; layer_idx < NUM_LAYERS; ++layer_idx) {
        (*training_data)[layer_idx].resize(args->p.arch[layer_idx]);
      }
    }

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
    Tuple<KeyType, PayloadType> *   relR_start_sampling_ptr = args->relR.tuples;
    Tuple<KeyType, PayloadType> *   relR_end_sampling_ptr = args->relR.tuples + (args->relR.num_tuples - 1);
    Tuple<KeyType, PayloadType> *   relS_start_sampling_ptr = args->relS.tuples;
    Tuple<KeyType, PayloadType> *   relS_end_sampling_ptr = args->relS.tuples + (args->relS.num_tuples - 1);

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

      tmp_training_sample_R[sample_count_R[tid]] = *i;
      ++sample_count_R[tid];
    }
    
    uint32_t * sample_count_S = args->sample_count_S;
    Tuple<KeyType, PayloadType> * tmp_training_sample_S = args->rmi->tmp_training_sample_S + args->tmp_training_sample_S_offset;
    unsigned int offset_S = static_cast<unsigned int>(1. * INPUT_SZ_S / SAMPLE_SZ_S);
    for (auto i = relS_start_sampling_ptr; i <= relS_end_sampling_ptr; i += offset_S) {
      // NOTE:  We don't directly assign SAMPLE_SZ to rmi.training_sample_sz
      //        to avoid issues with divisibility
      tmp_training_sample[sample_count[tid]] = *i;
      ++sample_count[tid];

      tmp_training_sample_S[sample_count_S[tid]] = *i;
      ++sample_count_S[tid];
    }
    BARRIER_ARRIVE(args->barrier, rv);
/*
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
        sorted_training_sample[k].key = tmp_outputptr[k].key;
        sorted_training_sample[k].payload = tmp_outputptr[k].payload;
      } 
      args->rmi->training_sample = &(sorted_training_sample);
      args->rmi->training_sample_size = total_sample_count;
    }

    if(tid == 1)
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
        sorted_training_sample_R[k].key = tmp_outputptr_R[k].key;
        sorted_training_sample_R[k].payload = tmp_outputptr_R[k].payload;
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
        sorted_training_sample_S[k].key = tmp_outputptr_S[k].key;
        sorted_training_sample_S[k].payload = tmp_outputptr_S[k].payload;
      } 
      args->rmi->training_sample_S = &(sorted_training_sample_S);
      args->rmi->training_sample_size_S = total_sample_count_S;
    }
    #else
    uint32_t total_sample_count = 0;
    if(tid == 0)
    {
      for(int i = 0; i < NUM_THREADS; i++)
        total_sample_count += sample_count[i];

      std::sort((int64_t *)(args->rmi->tmp_training_sample), (int64_t *)(args->rmi->tmp_training_sample) + total_sample_count - 1);
      args->rmi->training_sample = &(args->rmi->tmp_training_sample);
      args->rmi->training_sample_size = total_sample_count;
    }

    if(tid == 1)
    {
      uint32_t total_sample_count_R = 0; 
      for(int i = 0; i < NUM_THREADS; i++)
        total_sample_count_R += sample_count_R[i];

      std::sort((int64_t *)(args->rmi->tmp_training_sample_R), (int64_t *)(args->rmi->tmp_training_sample_R) + total_sample_count_R - 1);
      args->rmi->training_sample_R = &(args->rmi->tmp_training_sample_R);
      args->rmi->training_sample_size_R = total_sample_count_R;
    }

    if(tid == 2)
    {
      uint32_t total_sample_count_S = 0; 
      for(int i = 0; i < NUM_THREADS; i++)
        total_sample_count_S += sample_count_S[i];

      std::sort((int64_t *)(args->rmi->tmp_training_sample_S), (int64_t *)(args->rmi->tmp_training_sample_S) + total_sample_count_S - 1);
      args->rmi->training_sample_S = &(args->rmi->tmp_training_sample_S);
      args->rmi->training_sample_size_S = total_sample_count_S;
    }
    #endif

    //BARRIER_ARRIVE(args->barrier, rv);

    // Stop early if the array is identical
    if (((*(args->rmi->training_sample))[0]).key == ((*(args->rmi->training_sample))[total_sample_count - 1]).key) {
      return;
    }

    //----------------------------------------------------------//
    //                     TRAIN THE MODELS                     //
    //----------------------------------------------------------//

    if(tid==0)
    {
    printf("key1 %ld, key2 %ld \n", ((*(args->rmi->training_sample))[0]).key, ((*(args->rmi->training_sample))[total_sample_count - 1]).key);
    printf("key1 %ld, key2 %ld \n", ((*(args->rmi->training_sample))[0]).key, ((*(args->rmi->training_sample))[total_sample_count - 1]).key);        
      // Populate the training data for the root model
      //for (unsigned int i = 0; i < total_sample_count; ++i) {
        //printf("i %f \n", 1. * i / total_sample_count);
        //printf("key %ld, payload %ld, value %lf \n", ((*(args->rmi->training_sample))[i]).key, ((*(args->rmi->training_sample))[i]).payload, 1. * i / total_sample_count);
        //(*training_data)[0][0].push_back({((*(args->rmi->training_sample))[i]), 1. * i / total_sample_count});
      //}

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

      // Populate the training data for the next layer
      for (const auto &d : *current_training_data) {
        // Predict the model index in next layer
        unsigned int rank = current_model->slope * d.x.key + current_model->intercept;

        // Normalize the rank between 0 and the number of models in the next layer
        rank =
            std::max(static_cast<unsigned int>(0), std::min(args->p.arch[1] - 1, rank));

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
                (max.y) / (max.x.key - min.x.key);  // Hallucinating as if min.y = 0
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
                (min.y - 1) / (min.x.key - max.x.key);  // Hallucinating as if max.y = 1
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

            current_model->slope = (min.y - max.y) / (min.x.key - max.x.key);
            current_model->intercept = min.y - current_model->slope * min.x.key;
          }
        }
      }

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
*/
    //printf("here before last barrier tid %d \n", tid);
    //BARRIER_ARRIVE(args->barrier, rv);
    //printf("here after last barrier tid %d \n", tid);
}
#endif

void * npj_join_thread(void * param)
{
    ETHNonPartitionJoinThread<KeyType, PayloadType, TaskType> * args   = (ETHNonPartitionJoinThread<KeyType, PayloadType, TaskType> *) param;
    int rv;   int deltaT = 0;

    sample_and_train_models_threaded(args);

    BARRIER_ARRIVE(args->barrier, rv);

/*    BucketBuffer<KeyType, PayloadType> * overflowbuf; // allocate overflow buffer for each thread
    uint32_t nbuckets = (args->relR.num_tuples / BUCKET_SIZE / NUM_THREADS);

    if (args->tid == 0) {
        strcpy(npj_pfun[1].fun_name, "IMV");
        strcpy(npj_pfun[0].fun_name, "Naive");

        npj_pfun[1].fun_ptr = npj_build_rel_r_partition_imv;
        npj_pfun[0].fun_ptr = npj_build_rel_r_partition;

        strcpy(npj_pfun1[1].fun_name, "IMV");
        strcpy(npj_pfun1[0].fun_name, "Naive");

        npj_pfun1[1].fun_ptr = npj_probe_rel_s_partition_imv;
        npj_pfun1[0].fun_ptr = npj_probe_rel_s_partition;

        npj_pf_num = 2;
    }
    BARRIER_ARRIVE(args->barrier, rv);
    
    ETHNonPartitionJoinBuild<KeyType, PayloadType> build_data; 
    for (int fid = 0; fid < npj_pf_num; ++fid) 
    {
        for (int rp = 0; rp < RUN_NUMS; ++rp) 
        {
            init_bucket_buffer(&overflowbuf);
            if(args->tid == 0)
                allocate_hashtable(&args->ht, nbuckets);
            //BARRIER_ARRIVE(args->barrier, rv);

            build_data.ht = args->ht;
            build_data.overflowbuf = &overflowbuf;

        #ifdef PERF_COUNTERS
            if(args->tid == 0){
                //TODO: performance counters to be implemented
            }
        #endif

            BARRIER_ARRIVE(args->barrier, rv);

            if(args->tid == 0){
                gettimeofday(&args->start_time, NULL);
            #ifndef DEVELOPMENT_MODE
                //args->e_start_to_partition.startCounters();
            #endif
            }

        #ifdef DEVELOPMENT_MODE
            //build_data.build_hash_bucket_visits = args->build_hash_bucket_visits;
            //build_data.probe_hash_bucket_visits = args->probe_hash_bucket_visits; 
            //build_data.keys_hash_latch = args->keys_hash_latch; 
            //build_data.build_keys_list = args->build_keys_list;
            //build_data.build_keys_hash_list = args->build_keys_hash_list;        
            //build_data.probe_keys_list = args->probe_keys_list;
            //build_data.probe_keys_hash_list = args->probe_keys_hash_list;       
        #endif        

        #if NPJ_MORSE_SIZE
            morse_driven(param, npj_pfun[fid].fun_ptr, &overflowbuf);
        #else
            npj_pfun[fid].fun_ptr(&build_data, &args->relR, NULL);
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
                printf("---- %5s Build costs time (ms) = %10.4lf\n", npj_pfun[fid].fun_name, deltaT * 1.0 / 1000);
                npj_total_num = 0;
                npj_global_curse = 0;
            }
        
            if(!((fid == (npj_pf_num - 1)) && (rp == (RUN_NUMS - 1)))){
                if(args->tid == 0)
                    destroy_hashtable(args->ht);

                free_bucket_buffer(overflowbuf);
                
                BARRIER_ARRIVE(args->barrier, rv);
            } 
        }
    }

    BARRIER_ARRIVE(args->barrier, rv);

    //Probe phase
    for (int fid = 0; fid < npj_pf_num; ++fid) 
    {
        for (int rp = 0; rp < RUN_NUMS; ++rp) 
        {
            BARRIER_ARRIVE(args->barrier, rv);
            if(args->tid == 0){
                gettimeofday(&args->partition_end_time, NULL);

            #if NPJ_MORSE_SIZE
                //TODO: to be done
            #else
                args->num_results = npj_pfun1[fid].fun_ptr(NULL, &args->relS, &build_data);
            #endif
            }

            BARRIER_ARRIVE(args->barrier, rv);

            // probe phase finished, thread-0 checkpoints the time
            if(args->tid == 0){
                gettimeofday(&args->end_time, NULL);

                deltaT = (args->end_time.tv_sec - args->partition_end_time.tv_sec) * 1000000 + args->end_time.tv_usec - args->partition_end_time.tv_usec;
                printf("---- %5s Probe costs time (ms) = %10.4lf\n", npj_pfun1[fid].fun_name, deltaT * 1.0 / 1000);
            }
        }
    }
*/
    return 0;
}

void initialize_npj_join_thread_args(Relation<KeyType, PayloadType> * rel_r, 
                                 Relation<KeyType, PayloadType> * rel_s,
                                 Hashtable<KeyType, PayloadType> * ht, 
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
                                 pthread_barrier_t* barrier_ptr,
                                 Result * joinresult,
                                 ETHNonPartitionJoinThread<KeyType, PayloadType, TaskType> * args){
    int i;
    uint64_t numR, numS, numRthr, numSthr; /* total and per thread num */
    unsigned int SAMPLE_SZ_Rthr, SAMPLE_SZ_Sthr;

    numR = rel_r->num_tuples;
    numS = rel_s->num_tuples;
    numRthr = numR / NUM_THREADS;
    numSthr = numS / NUM_THREADS;
    SAMPLE_SZ_Rthr = SAMPLE_SZ_R / NUM_THREADS;
    SAMPLE_SZ_Sthr = SAMPLE_SZ_S / NUM_THREADS;

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
        /**** end stuff for learning RMI models ****/

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

    uint32_t * sample_count = (uint32_t *) calloc(NUM_THREADS, sizeof(uint32_t)); 
    uint32_t * sample_count_R = (uint32_t *) calloc(NUM_THREADS, sizeof(uint32_t)); 
    uint32_t * sample_count_S = (uint32_t *) calloc(NUM_THREADS, sizeof(uint32_t));

    //////////////////////////////////////////////////////////////////////////////
    // End stuff for sampling and building RMI models for both relations R and S
    //////////////////////////////////////////////////////////////////////////////

    //initialize_npj_join_thread_args(&rel_r, &rel_s, ht, &barrier, joinresult, args_ptr);
    initialize_npj_join_thread_args(&rel_r, &rel_s, ht, &rmi, rmi_params,
                                 SAMPLE_SZ_R, SAMPLE_SZ_S,
                                 tmp_training_sample_in, sorted_training_sample_in, r_tmp_training_sample_in,
                                 r_sorted_training_sample_in, s_tmp_training_sample_in, s_sorted_training_sample_in,
                                 &training_data, sample_count, sample_count_R, sample_count_S,
                                 &barrier, joinresult, args_ptr);

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


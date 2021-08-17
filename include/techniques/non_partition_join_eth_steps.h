#pragma once

/* An implementation for non-partitioned hash join (vectorized and non-vectorized) (NOT USED NOW) */

#include "emmintrin.h"
#include "immintrin.h"
#include "smmintrin.h"

#include "techniques/hash_join_steps_base.h"
#include "utils/base_utils.h"
#include "utils/math.h"
#include "utils/barrier.h"
#include "utils/eth_data_structures.h"
#include "utils/memory.h"
#include "utils/lock.h" 

#include "configs/base_configs.h"
#include "configs/eth_configs.h"

#define SINGLE_TUPLE_PER_BUCKET

typedef struct StateSIMDForETHNPJ StateSIMDForETHNPJ;
struct StateSIMDForETHNPJ {
  __m512i key;
  __m512i tb_off;
  __m512i ht_off;
  __mmask8 m_have_tuple;
  char stage;
};

// Based on the implementation of non partition join from ETH
template<typename KeyType, typename PayloadType, typename TaskType, typename JoinThreadType, typename PartitionType>
class ETHNonPartitionJoinSteps : public HashJoinSteps<KeyType, PayloadType, TaskType, JoinThreadType, PartitionType>{
 public:
  
  template<typename BuildType>
  void build_rel_r_partition(BuildType * build_input, const Relation<KeyType, PayloadType> * const rel_r_partition, Relation<KeyType, PayloadType> * const tmp_r) 
  {
    Hashtable<KeyType, PayloadType>* ht = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_input)->ht;  
    BucketBuffer<KeyType, PayloadType>** overflowbuf = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_input)->overflowbuf;
#ifdef DEVELOPMENT_MODE
     unordered_map<uint64_t, uint64_t>* build_visits_map = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_input)->build_hash_bucket_visits;        
     vector<KeyType>* build_keys_list = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_input)->build_keys_list;  
     vector<uint64_t>* build_keys_hash_list = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_input)->build_keys_hash_list;        
     volatile char * keys_hash_latch = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_input)->keys_hash_latch;        
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
        (*build_visits_map)[idx]++;

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
        else {
            dest = curr->tuples + curr->count;
            curr->count ++;
        }

        *dest = rel_r_partition->tuples[i];
        unlock(&curr->latch);
    }
  }

  template<typename BuildType>
  void build_rel_r_partition_imv(BuildType * build_input, const Relation<KeyType, PayloadType> * const rel_r_partition, Relation<KeyType, PayloadType> * const tmp_r) 
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

  template<typename BuildType>
  uint64_t probe_rel_s_partition(const Relation<KeyType, PayloadType> * const rel_r_partition, const Relation<KeyType, PayloadType> * const rel_s_partition,
                                 const BuildType * build_output) {
    Hashtable<KeyType, PayloadType>* ht = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_output)->ht;  
#ifdef DEVELOPMENT_MODE
    unordered_map<uint64_t, uint64_t>* probe_visits_map = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_output)->probe_hash_bucket_visits;        
    vector<KeyType>* probe_keys_list = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_output)->probe_keys_list;  
    vector<uint64_t>* probe_keys_hash_list = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_output)->probe_keys_hash_list;        
    volatile char* keys_hash_latch = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_output)->keys_hash_latch;        
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
        lock(&b->latch);
        (*probe_visits_map)[idx]++;
        unlock(&b->latch);

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

  // NOTE: v_tuple_size and v_next_off are calculated based on the implementation of Bucket<Keytype, PayloadType> 
  template<typename BuildType>
  uint64_t probe_rel_s_partition_imv(const Relation<KeyType, PayloadType> * const rel_r_partition, const Relation<KeyType, PayloadType> * const rel_s_partition,
                                     const BuildType * build_output, int thread_id) {

    Hashtable<KeyType, PayloadType>* ht = ((ETHNonPartitionJoinBuild<KeyType, PayloadType> *)build_output)->ht;  
    uint64_t matches = 0;

    int32_t k = 0, num, num_temp, done = 0, new_add = 0;
    __m512i v_ht_cell, v_base_offset_upper = _mm512_set1_epi64(rel_s_partition->num_tuples * sizeof(Tuple<KeyType, PayloadType>)), v_offset, v_base_offset, v_cell_hash,
            v_bucket_size = _mm512_set1_epi64(sizeof(Bucket<KeyType, PayloadType>)), v_ht_addr = _mm512_set1_epi64((uint64_t)ht->buckets), v_neg_one512 = _mm512_set1_epi64(-1),
            v_zero512 = _mm512_set1_epi64(0), v_tuple_size = _mm512_set1_epi64(16), v_next_off = _mm512_set1_epi64(8);
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

    #ifndef USE_MURMUR3_HASH
    __m512i v_factor = _mm512_set1_epi64(ht->hash_mask), v_shift = _mm512_set1_epi64(ht->skip_bits); 
    #endif   

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
                //TODO: murmur3 hashing to be implemented
                printf("Murmur3 hashing is not implemented in the probe phase yet!! \n");
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

  ~ETHNonPartitionJoinSteps() {
  }
 
};
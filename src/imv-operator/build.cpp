#include "imv-operator/build.hpp"

#include <map>
using namespace std;

#define SERIAL_BUILD 1
int64_t build_raw(hashtable_t *ht, relation_t *rel, bucket_buffer_t **overflowbuf) {
  uint64_t found = 0;
  uint32_t i;
  const uint32_t hashmask = ht->hash_mask;
  const uint32_t skipbits = ht->skip_bits;
  bucket_t *b;
  for (i = 0; i < rel->num_tuples; i++) {
    tuple_t *dest;
    bucket_t *curr, *nxt;

    int32_t idx = HASH(rel->tuples[i].key, hashmask, skipbits);
    /* copy the tuple to appropriate hash bucket */
    /* if full, follow nxt pointer to find correct place */
    curr = ht->buckets + idx;
    nxt = curr->next;
    if (curr->count == 0) {
      dest = curr->tuples;
      curr->count = 1;
      curr->lenth = 1;
    } else {
      get_new_bucket(&b, overflowbuf);
      b->next = curr->next;
      curr->next = b;
      b->count = 1;
      b->lenth++;
      dest = b->tuples;
    }
    ++found;
    *dest = rel->tuples[i];
  }
  return found;
}
int64_t build_AMAC(hashtable_t *ht, relation_t *rel, bucket_buffer_t **overflowbuf) {
  uint64_t found = 0;
  uint64_t key, hash_value, new_add = 0, cur = 0, end = rel->num_tuples;
  tuple_t *dest;
  bucket_t *cur_node, *next_node, *bucket;
  scalar_state_t state[ScalarStateSize];
  for (int i = 0; i < ScalarStateSize; ++i) {
    state[i].stage = 1;
  }
  int done = 0, k = 0;

  while (done < ScalarStateSize) {
    k = (k >= ScalarStateSize) ? 0 : k;
    switch (state[k].stage) {
      case 1: {
        if (cur >= end) {
          ++done;
          state[k].stage = 3;
          break;
        }
#if SEQPREFETCH
        _mm_prefetch(((char *)(rel->tuples + cur) + PDIS), _MM_HINT_T0);
#endif
        // step 1: get key
        key = rel->tuples[cur].key;
        // step 2: hash
        hash_value = HASH(key, ht->hash_mask, ht->skip_bits);
        state[k].b = ht->buckets + hash_value;
        state[k].stage = 0;
        state[k].key = key;
        state[k].payload = rel->tuples[cur].payload;
        ++cur;
        _mm_prefetch((char * )(state[k].b), _MM_HINT_T0);

      }
        break;
      case 0: {
        cur_node = state[k].b;
        if (cur_node->count == 0) {
          cur_node->count = 1;
          cur_node->tuples->key = state[k].key;
          cur_node->tuples->payload = state[k].payload;
          cur_node->next = NULL;
        } else {
          // step 4.3: insert new bucket at the tail
          get_new_bucket(&bucket, overflowbuf);
          bucket->next = cur_node->next;
          cur_node->next = bucket;
          bucket->tuples->key = state[k].key;
          bucket->tuples->payload = state[k].payload;
          bucket->count = 1;
        }
        state[k].stage = 1;
        --k;
        ++found;
      }
    }
    ++k;
  }

  return found;
}
#define WORDSIZE 8
size_t build_imv(hashtable_t *ht, relation_t *rel, bucket_buffer_t **overflowbuf) {
  int32_t found = 0, k = 0, done = 0, num, num_temp;
  __attribute__((aligned(64)))        __mmask8 to_scatt = 0, m_match = 0, m_new_cells = -1, m_valid_bucket = 0, mask[VECTOR_SCALE + 1], m_to_insert = 0, m_no_conflict;
  __m512i v_offset = _mm512_set1_epi64(0), v_base_offset_upper = _mm512_set1_epi64(rel->num_tuples * sizeof(tuple_t)), v_base_offset, v_ht_cell, v_factor = _mm512_set1_epi64(
      ht->hash_mask), v_shift = _mm512_set1_epi64(ht->skip_bits), v_cell_hash, v_neg_one512 = _mm512_set1_epi64(-1), v_zero512 = _mm512_set1_epi64(0), v_write_index =
      _mm512_set1_epi64(0), v_ht_addr = _mm512_set1_epi64((uint64_t) ht->buckets), v_word_size = _mm512_set1_epi64(WORDSIZE), v_tuple_size = _mm512_set1_epi64(sizeof(tuple_t)),
      v_bucket_size = _mm512_set1_epi64(sizeof(bucket_t)), v_next_off = _mm512_set1_epi64(offsetof(bucket_t, next)), v_right_payload, v_addr, v_all_ones = _mm512_set1_epi64(-1),
      v_conflict, v_one = _mm512_set1_epi64(1), v_new_bucket, v_next, v_63 = _mm512_set1_epi64(63), v_key_off = _mm512_set1_epi64(offsetof(bucket_t, tuples[0].key)), v_lzeros,
      v_previous, v_payload_off = _mm512_set1_epi64(offsetof(bucket_t, tuples[0].payload)), v_count_off = _mm512_set1_epi64(offsetof(bucket_t, count));
  __m256i v256_one = _mm256_set1_epi32(1);
  tuple_t *join_res = NULL;
  uint64_t *pos = NULL, *new_bucket = (uint64_t*) &v_new_bucket;
  bucket_t * bucket;
  int tail_add = 0;
  __attribute__((aligned(64)))     uint64_t cur_offset = 0, base_off[16], *ht_pos;
  for (int i = 0; i <= VECTOR_SCALE; ++i) {
    base_off[i] = i * sizeof(tuple_t);
    mask[i] = (1 << i) - 1;
  }
  v_base_offset = _mm512_load_epi64(base_off);
  __attribute__((aligned(64)))      StateSIMD state[SIMDStateSize + 1];
  // init # of the state
  for (int i = 0; i <= SIMDStateSize; ++i) {
    state[i].stage = 1;
    state[i].m_have_tuple = 0;
    state[i].ht_off = _mm512_set1_epi64(0);
    state[i].payload = _mm512_set1_epi64(0);
    state[i].key = _mm512_set1_epi64(0);
  }

  for (uint64_t cur = 0; 1;) {
    k = (k >= SIMDStateSize) ? 0 : k;
    if (UNLIKELY(cur >= rel->num_tuples)) {
      if (state[k].m_have_tuple == 0 && state[k].stage != 3) {
        ++done;
        state[k].stage = 3;
        ++k;
        continue;
      }
      if ((done >= SIMDStateSize)) {
        if (state[SIMDStateSize].m_have_tuple > 0) {
          k = SIMDStateSize;
          state[SIMDStateSize].stage = 0;
        } else {
          break;
        }
      }
    }
    switch (state[k].stage) {
      case 1: {
///////// step 1: load new tuples' address offsets
// the offset should be within MAX_32INT_
// the tail depends on the number of joins and tuples in each bucket
#if SEQPREFETCH
        _mm_prefetch((char *)(((void *)rel->tuples) + cur_offset + PDIS), _MM_HINT_T0);
        _mm_prefetch((char *)(((void *)rel->tuples) + cur_offset + PDIS + 64), _MM_HINT_T0);
        _mm_prefetch((char *)(((void *)rel->tuples) + cur_offset + PDIS + 128), _MM_HINT_T0);
#endif
        // directly use cur, instead of cur_offset to control the offset to rel.
        // In this case, using step = 16 to gather data, but step is larger
        // than the scale 1,2,4 or 8
        v_offset = _mm512_add_epi64(_mm512_set1_epi64(cur_offset), v_base_offset);
        // count the number of empty tuples
        cur_offset = cur_offset + base_off[VECTOR_SCALE];
        state[k].m_have_tuple = _mm512_cmpgt_epi64_mask(v_base_offset_upper, v_offset);
        cur = cur + VECTOR_SCALE;
        ///// step 2: load new cells from right tuples;
        // maybe need offset within a tuple
        state[k].key = _mm512_mask_i64gather_epi64(state[k].key, state[k].m_have_tuple, v_offset, ((void * )rel->tuples), 1);
        state[k].payload = _mm512_mask_i64gather_epi64(state[k].payload, state[k].m_have_tuple, _mm512_add_epi64(v_offset, v_word_size), ((void * )rel->tuples), 1);
        ///// step 3: load new values from hash tables;
        // hash the cell values
        v_cell_hash = _mm512_and_epi64(state[k].key, v_factor);
        v_cell_hash = _mm512_srlv_epi64(v_cell_hash, v_shift);
        v_cell_hash = _mm512_mullo_epi64(v_cell_hash, v_bucket_size);
        state[k].ht_off = _mm512_mask_add_epi64(state[k].ht_off, state[k].m_have_tuple, v_cell_hash, v_ht_addr);
        state[k].stage = 2;
#if KNL
        _mm512_mask_prefetch_i64gather_pd(
            state[k].ht_off, state[k].m_have_tuple, 0, 1, _MM_HINT_T0);
#else
        ht_pos = (uint64_t *) &state[k].ht_off;
        for (int i = 0; i < VECTOR_SCALE; ++i) {
          _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
        }
#endif
      }
        break;
      case 2: {
        /////////////////// random access
        // check valid bucket, insert directly for invalid buckets
        v_ht_cell = _mm512_mask_i64gather_epi64(v_neg_one512, state[k].m_have_tuple, state[k].ht_off, 0, 1);
        // inset new nodes
        m_to_insert = _mm512_cmpeq_epi64_mask(v_ht_cell, v_zero512);
        m_to_insert = _mm512_kand(m_to_insert, state[k].m_have_tuple);
        if (m_to_insert == 0) {
          state[k].stage = 0;
          --k;
          break;
        }
        v_addr = _mm512_mask_blend_epi64(m_to_insert, v_all_ones, state[k].ht_off);
        v_conflict = _mm512_conflict_epi64(v_addr);
        m_no_conflict = _mm512_testn_epi64_mask(v_conflict, v_all_ones);
        m_no_conflict = _mm512_kand(m_no_conflict, m_to_insert);
        if (m_no_conflict) {
          // write the key , payload, count, next to the nodes
          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_addr, v_key_off), state[k].key, 1);
          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_addr, v_payload_off), state[k].payload, 1);
          _mm512_mask_i64scatter_epi32(0, m_no_conflict, _mm512_add_epi64(v_addr, v_count_off), v256_one, 1);
//          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_addr,_mm512_set1_epi64(offsetof(bucket_t,next))), v_zero512, 1);
          state[k].m_have_tuple = _mm512_kandn(m_no_conflict, state[k].m_have_tuple);
          found += _mm_popcnt_u32(m_no_conflict);
        }
        num = _mm_popcnt_u32(state[k].m_have_tuple);
        if (num == VECTOR_SCALE || done >= SIMDStateSize) {
          state[k].stage = 0;
          --k;
        } else if (num == 0) {
          state[k].stage = 1;
          --k;
        } else {
          if (LIKELY(done < SIMDStateSize)) {
            num_temp = _mm_popcnt_u32(state[SIMDStateSize].m_have_tuple);
            if (num + num_temp < VECTOR_SCALE) {
              // compress v
              state[k].ht_off = _mm512_maskz_compress_epi64(state[k].m_have_tuple, state[k].ht_off);
              state[k].key = _mm512_maskz_compress_epi64(state[k].m_have_tuple, state[k].key);
              state[k].payload = _mm512_maskz_compress_epi64(state[k].m_have_tuple, state[k].payload);
              // expand v -> temp
              state[SIMDStateSize].ht_off = _mm512_mask_expand_epi64(state[SIMDStateSize].ht_off, _mm512_knot(state[SIMDStateSize].m_have_tuple), state[k].ht_off);
              state[SIMDStateSize].key = _mm512_mask_expand_epi64(state[SIMDStateSize].key, _mm512_knot(state[SIMDStateSize].m_have_tuple), state[k].key);
              state[SIMDStateSize].payload = _mm512_mask_expand_epi64(state[SIMDStateSize].payload, _mm512_knot(state[SIMDStateSize].m_have_tuple), state[k].payload);
              state[SIMDStateSize].m_have_tuple = mask[num + num_temp];
              state[k].m_have_tuple = 0;
              state[k].stage = 1;
              --k;
            } else {
              // expand temp -> v
              state[k].ht_off = _mm512_mask_expand_epi64(state[k].ht_off, _mm512_knot(state[k].m_have_tuple), state[SIMDStateSize].ht_off);
              state[k].key = _mm512_mask_expand_epi64(state[k].key, _mm512_knot(state[k].m_have_tuple), state[SIMDStateSize].key);
              state[k].payload = _mm512_mask_expand_epi64(state[k].payload, _mm512_knot(state[k].m_have_tuple), state[SIMDStateSize].payload);
              // compress temp
              state[SIMDStateSize].m_have_tuple = _mm512_kand(state[SIMDStateSize].m_have_tuple, _mm512_knot(mask[VECTOR_SCALE - num]));
              state[SIMDStateSize].ht_off = _mm512_maskz_compress_epi64(state[SIMDStateSize].m_have_tuple, state[SIMDStateSize].ht_off);
              state[SIMDStateSize].key = _mm512_maskz_compress_epi64(state[SIMDStateSize].m_have_tuple, state[SIMDStateSize].key);
              state[SIMDStateSize].payload = _mm512_maskz_compress_epi64(state[SIMDStateSize].m_have_tuple, state[SIMDStateSize].payload);
              state[k].m_have_tuple = mask[VECTOR_SCALE];
              state[SIMDStateSize].m_have_tuple = (state[SIMDStateSize].m_have_tuple >> (VECTOR_SCALE - num));
              state[k].stage = 0;
              --k;
            }
          }
        }

      }
        break;
      case 0: {
        // insert new buckets after valid buckets
        m_to_insert = state[k].m_have_tuple;
        v_addr = _mm512_mask_blend_epi64(m_to_insert, v_all_ones, state[k].ht_off);
        v_conflict = _mm512_conflict_epi64(v_addr);
        m_no_conflict = _mm512_testn_epi64_mask(v_conflict, v_all_ones);
        m_no_conflict = _mm512_kand(m_no_conflict, m_to_insert);
        if (true) {
#if SERIAL_BUILD
          for (int i = 0; i < VECTOR_SCALE; ++i) {
            new_bucket[i] = 0;
            if (m_no_conflict & (1 << i)) {
              get_new_bucket(&bucket, overflowbuf);
              new_bucket[i] = bucket;
            }
          }
          v_next = _mm512_mask_i64gather_epi64(v_all_ones, state[k].m_have_tuple, _mm512_add_epi64(state[k].ht_off, v_next_off), 0, 1);
          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_new_bucket, v_key_off), state[k].key, 1);
          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_new_bucket, v_payload_off), state[k].payload, 1);
//          _mm512_mask_i64scatter_epi32(0, m_no_conflict, (v_new_bucket), v256_one, 1);
          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_new_bucket, v_next_off), v_next, 1);

          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(state[k].ht_off, v_next_off), v_new_bucket, 1);

          found += _mm_popcnt_u32(m_no_conflict);
          state[k].m_have_tuple = _mm512_kandn(m_no_conflict, state[k].m_have_tuple);
#else
          for (int i = 0; i < VECTOR_SCALE; ++i) {
            new_bucket[i] = 0;
            if (m_to_insert & (1 << i)) {
              get_new_bucket(&bucket, overflowbuf);
              new_bucket[i] = bucket;
            }
          }
          v_next = _mm512_mask_i64gather_epi64(v_all_ones, state[k].m_have_tuple, _mm512_add_epi64(state[k].ht_off, v_next_off), 0, 1);
          _mm512_mask_i64scatter_epi64(0, m_to_insert, _mm512_add_epi64(v_new_bucket, v_key_off), state[k].key, 1);
          _mm512_mask_i64scatter_epi64(0, m_to_insert, _mm512_add_epi64(v_new_bucket, v_payload_off), state[k].payload, 1);
//          _mm512_mask_i64scatter_epi32(0, m_to_insert, (v_new_bucket), v256_one, 1);
          _mm512_mask_i64scatter_epi64(0, m_to_insert, _mm512_add_epi64(v_new_bucket, v_next_off), v_next, 1);

          // conflict-solved insert
          v_lzeros = _mm512_lzcnt_epi64(v_conflict);
          v_lzeros = _mm512_sub_epi64(v_63, v_lzeros);
          to_scatt = _mm512_kandn(m_no_conflict, m_to_insert);
          v_previous = _mm512_maskz_permutexvar_epi64(to_scatt, v_lzeros, v_new_bucket);
          _mm512_mask_i64scatter_epi64(0, to_scatt, _mm512_add_epi64(v_previous, v_next_off), v_new_bucket, 1);

          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(state[k].ht_off, v_next_off), v_new_bucket, 1);
          found += _mm_popcnt_u32(m_to_insert);
          state[k].m_have_tuple = _mm512_kandn(m_to_insert, state[k].m_have_tuple);
#endif
        }

        num = _mm_popcnt_u32(state[k].m_have_tuple);
#if 1

        if (num == VECTOR_SCALE || done >= SIMDStateSize) {
#if KNL
          _mm512_mask_prefetch_i64gather_pd(
              state[k].ht_off, state[k].m_have_tuple, 0, 1, _MM_HINT_T0);
#else
          ht_pos = (uint64_t *) &state[k].ht_off;
          for (int i = 0; i < VECTOR_SCALE; ++i) {
            _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
//            _mm_prefetch((char * )(ht_pos[i] + 64), _MM_HINT_T0);
          }
#endif
        } else if (num == 0) {
          state[k].stage = 1;
          --k;
          break;
        } else
#endif
        {
          if (LIKELY(done < SIMDStateSize)) {
            num_temp = _mm_popcnt_u32(state[SIMDStateSize].m_have_tuple);
            if (num + num_temp < VECTOR_SCALE) {
              // compress v
              state[k].ht_off = _mm512_maskz_compress_epi64(state[k].m_have_tuple, state[k].ht_off);
              state[k].key = _mm512_maskz_compress_epi64(state[k].m_have_tuple, state[k].key);
              state[k].payload = _mm512_maskz_compress_epi64(state[k].m_have_tuple, state[k].payload);
              // expand v -> temp
              state[SIMDStateSize].ht_off = _mm512_mask_expand_epi64(state[SIMDStateSize].ht_off, _mm512_knot(state[SIMDStateSize].m_have_tuple), state[k].ht_off);
              state[SIMDStateSize].key = _mm512_mask_expand_epi64(state[SIMDStateSize].key, _mm512_knot(state[SIMDStateSize].m_have_tuple), state[k].key);
              state[SIMDStateSize].payload = _mm512_mask_expand_epi64(state[SIMDStateSize].payload, _mm512_knot(state[SIMDStateSize].m_have_tuple), state[k].payload);
              state[SIMDStateSize].m_have_tuple = mask[num + num_temp];
              state[k].m_have_tuple = 0;
              state[k].stage = 1;
              --k;
              break;
            } else {
              // expand temp -> v
              state[k].ht_off = _mm512_mask_expand_epi64(state[k].ht_off, _mm512_knot(state[k].m_have_tuple), state[SIMDStateSize].ht_off);
              state[k].key = _mm512_mask_expand_epi64(state[k].key, _mm512_knot(state[k].m_have_tuple), state[SIMDStateSize].key);

              state[k].payload = _mm512_mask_expand_epi64(state[k].payload, _mm512_knot(state[k].m_have_tuple), state[SIMDStateSize].payload);
              // compress temp
              state[SIMDStateSize].m_have_tuple = _mm512_kand(state[SIMDStateSize].m_have_tuple, _mm512_knot(mask[VECTOR_SCALE - num]));
              state[SIMDStateSize].ht_off = _mm512_maskz_compress_epi64(state[SIMDStateSize].m_have_tuple, state[SIMDStateSize].ht_off);
              state[SIMDStateSize].key = _mm512_maskz_compress_epi64(state[SIMDStateSize].m_have_tuple, state[SIMDStateSize].key);
              state[SIMDStateSize].payload = _mm512_maskz_compress_epi64(state[SIMDStateSize].m_have_tuple, state[SIMDStateSize].payload);
              state[k].m_have_tuple = mask[VECTOR_SCALE];
              state[SIMDStateSize].m_have_tuple = (state[SIMDStateSize].m_have_tuple >> (VECTOR_SCALE - num));
              state[k].stage = 0;
#if KNL
              _mm512_mask_prefetch_i64gather_pd(
                  state[k].ht_off, state[k].m_have_tuple, 0, 1, _MM_HINT_T0);
#else
              ht_pos = (uint64_t *) &state[k].ht_off;
              for (int i = 0; i < VECTOR_SCALE; ++i) {
                _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
//                _mm_prefetch((char * )(ht_pos[i] + 64), _MM_HINT_T0);
              }
#endif
            }
          }
        }
      }
        break;
    }
    ++k;
  }

  return found;
}

size_t build_DVA(hashtable_t *ht, relation_t *rel, bucket_buffer_t **overflowbuf) {
  int32_t found = 0, k = 0, done = 0, num, num_temp;
  __attribute__((aligned(64)))     __mmask8 to_scatt = 0, m_match = 0, m_new_cells = -1, m_valid_bucket = 0, mask[VECTOR_SCALE + 1], m_to_insert = 0, m_no_conflict;
  __m512i v_offset = _mm512_set1_epi64(0), v_base_offset_upper = _mm512_set1_epi64(rel->num_tuples * sizeof(tuple_t)), v_base_offset, v_ht_cell, v_factor = _mm512_set1_epi64(
      ht->hash_mask), v_shift = _mm512_set1_epi64(ht->skip_bits), v_cell_hash, v_neg_one512 = _mm512_set1_epi64(-1), v_zero512 = _mm512_set1_epi64(0), v_write_index =
      _mm512_set1_epi64(0), v_ht_addr = _mm512_set1_epi64((uint64_t) ht->buckets), v_word_size = _mm512_set1_epi64(WORDSIZE), v_tuple_size = _mm512_set1_epi64(sizeof(tuple_t)),
      v_bucket_size = _mm512_set1_epi64(sizeof(bucket_t)), v_next_off = _mm512_set1_epi64(offsetof(bucket_t, next)), v_right_payload, v_addr, v_all_ones = _mm512_set1_epi64(-1),
      v_conflict, v_one = _mm512_set1_epi64(1), v_new_bucket, v_next, v_63 = _mm512_set1_epi64(63), v_key_off = _mm512_set1_epi64(offsetof(bucket_t, tuples[0].key)), v_lzeros,
      v_previous, v_payload_off = _mm512_set1_epi64(offsetof(bucket_t, tuples[0].payload)), v_count_off = _mm512_set1_epi64(offsetof(bucket_t, count));
  __m256i v256_one = _mm256_set1_epi32(1);
  tuple_t *join_res = NULL;
  uint64_t *pos = NULL, *new_bucket = (uint64_t*) &v_new_bucket;
  bucket_t * bucket;
  int tail_add = 0;
  __attribute__((aligned(64)))                       uint64_t cur_offset = 0, base_off[16], *ht_pos;
  for (int i = 0; i <= VECTOR_SCALE; ++i) {
    base_off[i] = i * sizeof(tuple_t);
    mask[i] = (1 << i) - 1;
  }
  v_base_offset = _mm512_load_epi64(base_off);
  __attribute__((aligned(64)))                       StateSIMD state[SIMDStateSize + 1];
  // init # of the state
  for (int i = 0; i <= SIMDStateSize; ++i) {
    state[i].stage = 1;
    state[i].m_have_tuple = 0;
    state[i].ht_off = _mm512_set1_epi64(0);
    state[i].payload = _mm512_set1_epi64(0);
    state[i].key = _mm512_set1_epi64(0);
  }
  for (uint64_t cur = 0; (cur < rel->num_tuples) || (done < SIMDStateSize);) {
    k = (k >= SIMDStateSize) ? 0 : k;
    if (UNLIKELY(cur >= rel->num_tuples)) {
      if (state[k].m_have_tuple == 0 && state[k].stage != 3) {
        ++done;
        state[k].stage = 3;
        ++k;
        continue;
      }
    }
    switch (state[k].stage) {
      case 1: {
///////// step 1: load new tuples' address offsets
// the offset should be within MAX_32INT_
// the tail depends on the number of joins and tuples in each bucket
#if SEQPREFETCH
        _mm_prefetch((char *)(((void *)rel->tuples) + cur_offset + PDIS), _MM_HINT_T0);
        _mm_prefetch((char *)(((void *)rel->tuples) + cur_offset + PDIS + 64), _MM_HINT_T0);
        _mm_prefetch((char *)(((void *)rel->tuples) + cur_offset + PDIS + 128), _MM_HINT_T0);
#endif
        // directly use cur, instead of cur_offset to control the offset to rel.
        // In this case, using step = 16 to gather data, but step is larger
        // than the scale 1,2,4 or 8
        v_offset = _mm512_add_epi64(_mm512_set1_epi64(cur_offset), v_base_offset);
        // count the number of empty tuples
        cur_offset = cur_offset + base_off[VECTOR_SCALE];
        cur = cur + VECTOR_SCALE;
        state[k].m_have_tuple = _mm512_cmpgt_epi64_mask(v_base_offset_upper, v_offset);
        ///// step 2: load new cells from right tuples;
        // maybe need offset within a tuple
        state[k].key = _mm512_mask_i64gather_epi64(state[k].key, state[k].m_have_tuple, v_offset, ((void * )rel->tuples), 1);
        state[k].payload = _mm512_mask_i64gather_epi64(state[k].payload, state[k].m_have_tuple, _mm512_add_epi64(v_offset, v_word_size), ((void * )rel->tuples), 1);
        ///// step 3: load new values from hash tables;
        // hash the cell values
        v_cell_hash = _mm512_and_epi64(state[k].key, v_factor);
        v_cell_hash = _mm512_srlv_epi64(v_cell_hash, v_shift);
        v_cell_hash = _mm512_mullo_epi64(v_cell_hash, v_bucket_size);
        state[k].ht_off = _mm512_mask_add_epi64(state[k].ht_off, state[k].m_have_tuple, v_cell_hash, v_ht_addr);
        state[k].stage = 2;
#if KNL
        _mm512_mask_prefetch_i64gather_pd(
            state[k].ht_off, state[k].m_have_tuple, 0, 1, _MM_HINT_T0);
#else
        ht_pos = (uint64_t *) &state[k].ht_off;
        for (int i = 0; i < VECTOR_SCALE; ++i) {
          _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
        }
#endif
        _mm_prefetch((char * )((*overflowbuf)->buf+(*overflowbuf)->count)+PDIS, _MM_HINT_T0);
        _mm_prefetch((char * )((*overflowbuf)->buf+(*overflowbuf)->count)+PDIS+64, _MM_HINT_T0);

      }
        break;
      case 2: {
        /////////////////// random access
        // check valid bucket, insert directly for invalid buckets
        v_ht_cell = _mm512_mask_i64gather_epi64(v_neg_one512, state[k].m_have_tuple, state[k].ht_off, 0, 1);
        // inset new nodes
        m_to_insert = _mm512_cmpeq_epi64_mask(v_ht_cell, v_zero512);
        m_to_insert = _mm512_kand(m_to_insert, state[k].m_have_tuple);
        v_addr = _mm512_mask_blend_epi64(m_to_insert, v_all_ones, state[k].ht_off);
        v_conflict = _mm512_conflict_epi64(v_addr);
        m_no_conflict = _mm512_testn_epi64_mask(v_conflict, v_all_ones);
        m_no_conflict = _mm512_kand(m_no_conflict, m_to_insert);
        if (m_no_conflict) {
          // write the key , payload, count, next to the nodes
          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_addr, v_key_off), state[k].key, 1);
          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_addr, v_payload_off), state[k].payload, 1);
          _mm512_mask_i64scatter_epi32(0, m_no_conflict, _mm512_add_epi64(v_addr, v_count_off), v256_one, 1);
//          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_addr,_mm512_set1_epi64(offsetof(bucket_t,next))), v_zero512, 1);

          state[k].m_have_tuple = _mm512_kandn(m_no_conflict, state[k].m_have_tuple);
          found += _mm_popcnt_u32(m_no_conflict);
        }
        state[k].stage = 0;
        --k;
      }
        break;
      case 0: {
        // insert new buckets after valid buckets
        m_to_insert = state[k].m_have_tuple;
        v_addr = _mm512_mask_blend_epi64(m_to_insert, v_all_ones, state[k].ht_off);
        v_conflict = _mm512_conflict_epi64(v_addr);
        m_no_conflict = _mm512_testn_epi64_mask(v_conflict, v_all_ones);
        m_no_conflict = _mm512_kand(m_no_conflict, m_to_insert);
        if (true) {
#if SERIAL_BUILD
          for (int i = 0; i < VECTOR_SCALE; ++i) {
            new_bucket[i] = 0;
            if (m_no_conflict & (1 << i)) {
              get_new_bucket(&bucket, overflowbuf);
              new_bucket[i] = bucket;
            }
          }
          v_next = _mm512_mask_i64gather_epi64(v_all_ones, state[k].m_have_tuple, _mm512_add_epi64(state[k].ht_off, v_next_off), 0, 1);
          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_new_bucket, v_key_off), state[k].key, 1);
          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_new_bucket, v_payload_off), state[k].payload, 1);
//          _mm512_mask_i64scatter_epi32(0, m_no_conflict, (v_new_bucket), v256_one, 1);
          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_new_bucket, v_next_off), v_next, 1);

          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(state[k].ht_off, v_next_off), v_new_bucket, 1);

          found += _mm_popcnt_u32(m_no_conflict);
          state[k].m_have_tuple = _mm512_kandn(m_no_conflict, state[k].m_have_tuple);
#else
          for (int i = 0; i < VECTOR_SCALE; ++i) {
            new_bucket[i] = 0;
            if (m_to_insert & (1 << i)) {
              get_new_bucket(&bucket, overflowbuf);
              new_bucket[i] = bucket;
            }
          }
          v_next = _mm512_mask_i64gather_epi64(v_all_ones, state[k].m_have_tuple, _mm512_add_epi64(state[k].ht_off, v_next_off), 0, 1);
          _mm512_mask_i64scatter_epi64(0, m_to_insert, _mm512_add_epi64(v_new_bucket, v_key_off), state[k].key, 1);
          _mm512_mask_i64scatter_epi64(0, m_to_insert, _mm512_add_epi64(v_new_bucket, v_payload_off), state[k].payload, 1);
//          _mm512_mask_i64scatter_epi32(0, m_to_insert, (v_new_bucket), v256_one, 1);
          _mm512_mask_i64scatter_epi64(0, m_to_insert, _mm512_add_epi64(v_new_bucket, v_next_off), v_next, 1);

          // conflict-solved insert
          v_lzeros = _mm512_lzcnt_epi64(v_conflict);
          v_lzeros = _mm512_sub_epi64(v_63, v_lzeros);
          to_scatt = _mm512_kandn(m_no_conflict, m_to_insert);
          v_previous = _mm512_maskz_permutexvar_epi64(to_scatt, v_lzeros, v_new_bucket);
          _mm512_mask_i64scatter_epi64(0, to_scatt, _mm512_add_epi64(v_previous, v_next_off), v_new_bucket, 1);

          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(state[k].ht_off, v_next_off), v_new_bucket, 1);
          found += _mm_popcnt_u32(m_to_insert);
          state[k].m_have_tuple = _mm512_kandn(m_to_insert, state[k].m_have_tuple);
#endif	
        }

        if (state[k].m_have_tuple) {
#if KNL
          _mm512_mask_prefetch_i64gather_pd(
              state[k].ht_off, state[k].m_have_tuple, 0, 1, _MM_HINT_T0);
#else
          ht_pos = (uint64_t *) &state[k].ht_off;
          for (int i = 0; i < VECTOR_SCALE; ++i) {
            _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
          }
#endif
        } else {
          state[k].stage = 1;
          --k;
        }
      }
        break;
    }
    ++k;
  }

  return found;
}
size_t build_FVA(hashtable_t *ht, relation_t *rel, bucket_buffer_t **overflowbuf) {
  int32_t found = 0, k = 0, done = 0, num, num_temp, new_add;
  __attribute__((aligned(64)))             __mmask8 to_scatt = 0, m_match = 0, m_new_cells = -1, m_valid_bucket = 0, mask[VECTOR_SCALE + 1], m_to_insert = 0, m_no_conflict, m_full = -1;
  __m512i v_offset = _mm512_set1_epi64(0), v_base_offset_upper = _mm512_set1_epi64(rel->num_tuples * sizeof(tuple_t)), v_base_offset, v_ht_cell, v_factor = _mm512_set1_epi64(
      ht->hash_mask), v_shift = _mm512_set1_epi64(ht->skip_bits), v_cell_hash, v_neg_one512 = _mm512_set1_epi64(-1), v_zero512 = _mm512_set1_epi64(0), v_write_index =
      _mm512_set1_epi64(0), v_ht_addr = _mm512_set1_epi64((uint64_t) ht->buckets), v_word_size = _mm512_set1_epi64(WORDSIZE), v_tuple_size = _mm512_set1_epi64(sizeof(tuple_t)),
      v_bucket_size = _mm512_set1_epi64(sizeof(bucket_t)), v_next_off = _mm512_set1_epi64(offsetof(bucket_t, next)), v_right_payload, v_addr, v_all_ones = _mm512_set1_epi64(-1),
      v_conflict, v_one = _mm512_set1_epi64(1), v_new_bucket, v_next, v_63 = _mm512_set1_epi64(63), v_key_off = _mm512_set1_epi64(offsetof(bucket_t, tuples[0].key)), v_lzeros,
      v_previous, v_payload_off = _mm512_set1_epi64(offsetof(bucket_t, tuples[0].payload)), v_count_off = _mm512_set1_epi64(offsetof(bucket_t, count));
  __m256i v256_one = _mm256_set1_epi32(1);
  tuple_t *join_res = NULL;
  uint64_t *pos = NULL, *new_bucket = (uint64_t*) &v_new_bucket;
  bucket_t * bucket;
  int tail_add = 0;
  __attribute__((aligned(64)))                       uint64_t cur_offset = 0, base_off[16], *ht_pos;
  for (int i = 0; i <= VECTOR_SCALE; ++i) {
    base_off[i] = i * sizeof(tuple_t);
    mask[i] = (1 << i) - 1;
  }
  v_base_offset = _mm512_load_epi64(base_off);
  __attribute__((aligned(64)))                       StateSIMD state[SIMDStateSize + 1];
  // init # of the state
  for (int i = 0; i <= SIMDStateSize; ++i) {
    state[i].stage = 1;
    state[i].m_have_tuple = 0;
    state[i].ht_off = _mm512_set1_epi64(0);
    state[i].payload = _mm512_set1_epi64(0);
    state[i].key = _mm512_set1_epi64(0);
  }
  for (uint64_t cur = 0; (cur < rel->num_tuples) || (done < SIMDStateSize);) {
    k = (k >= SIMDStateSize) ? 0 : k;
    if (UNLIKELY(cur >= rel->num_tuples)) {
      if (state[k].m_have_tuple == 0 && state[k].stage != 3) {
        ++done;
        state[k].stage = 3;
        ++k;
        continue;
      }
    }
    switch (state[k].stage) {
      case 1: {
        ///////// step 1: load new tuples' address offsets
        // the offset should be within MAX_32INT_
        // the tail depends on the number of joins and tuples in each bucket
#if SEQPREFETCH
        _mm_prefetch((char *)(((void *)rel->tuples) + cur_offset + PDIS), _MM_HINT_T0);
        _mm_prefetch((char *)(((void *)rel->tuples) + cur_offset + PDIS + 64), _MM_HINT_T0);
        _mm_prefetch((char *)(((void *)rel->tuples) + cur_offset + PDIS + 128), _MM_HINT_T0);
#endif
        // directly use cur, instead of cur_offset to control the offset to rel.
        // In this case, using step = 16 to gather data, but step is larger
        // than the scale 1,2,4 or 8
        v_offset = _mm512_add_epi64(_mm512_set1_epi64(cur_offset), v_base_offset);
        state[k].tb_off = _mm512_mask_expand_epi64(state[k].tb_off, _mm512_knot(state[k].m_have_tuple), v_offset);
        // count the number of empty tuples
        m_new_cells = _mm512_knot(state[k].m_have_tuple);
        new_add = _mm_popcnt_u32(m_new_cells);
        cur_offset = cur_offset + base_off[new_add];
        cur = cur + new_add;
        state[k].m_have_tuple = _mm512_cmpgt_epi64_mask(v_base_offset_upper, state[k].tb_off);
        ///// step 2: load new cells from right tuples;
        m_new_cells = _mm512_kand(m_new_cells, state[k].m_have_tuple);
        // maybe need offset within a tuple
        state[k].key = _mm512_mask_i64gather_epi64(state[k].key, m_new_cells, state[k].tb_off, ((void * )rel->tuples), 1);
        state[k].payload = _mm512_mask_i64gather_epi64(state[k].payload, state[k].m_have_tuple, _mm512_add_epi64(state[k].tb_off, v_word_size), ((void * )rel->tuples), 1);
        ///// step 3: load new values from hash tables;
        // hash the cell values
        v_cell_hash = _mm512_and_epi64(state[k].key, v_factor);
        v_cell_hash = _mm512_srlv_epi64(v_cell_hash, v_shift);
        v_cell_hash = _mm512_mullo_epi64(v_cell_hash, v_bucket_size);
        state[k].ht_off = _mm512_mask_add_epi64(state[k].ht_off, m_new_cells, v_cell_hash, v_ht_addr);
        state[k].stage = 2;
#if KNL
        _mm512_mask_prefetch_i64gather_pd(
            state[k].ht_off, state[k].m_have_tuple, 0, 1, _MM_HINT_T0);
#else
        ht_pos = (uint64_t *) &state[k].ht_off;
        for (int i = 0; i < VECTOR_SCALE; ++i) {
          _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
        }
#endif
        _mm_prefetch((char * )((*overflowbuf)->buf+(*overflowbuf)->count)+PDIS, _MM_HINT_T0);
        _mm_prefetch((char * )((*overflowbuf)->buf+(*overflowbuf)->count)+PDIS+64, _MM_HINT_T0);

      }
        break;
      case 2: {
        /////////////////// random access
        // check valid bucket, insert directly for invalid buckets
        v_ht_cell = _mm512_mask_i64gather_epi64(v_neg_one512, state[k].m_have_tuple, state[k].ht_off, 0, 1);
        // inset new nodes
        m_to_insert = _mm512_cmpeq_epi64_mask(v_ht_cell, v_zero512);
        m_to_insert = _mm512_kand(m_to_insert, state[k].m_have_tuple);
        v_addr = _mm512_mask_blend_epi64(m_to_insert, v_all_ones, state[k].ht_off);
        v_conflict = _mm512_conflict_epi64(v_addr);
        m_no_conflict = _mm512_testn_epi64_mask(v_conflict, v_all_ones);
        m_no_conflict = _mm512_kand(m_no_conflict, m_to_insert);
        if (m_no_conflict) {
          // write the key , payload, count, next to the nodes
          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_addr, v_key_off), state[k].key, 1);
          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_addr, v_payload_off), state[k].payload, 1);
          _mm512_mask_i64scatter_epi32(0, m_no_conflict, _mm512_add_epi64(v_addr, v_count_off), v256_one, 1);
//          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_addr,_mm512_set1_epi64(offsetof(bucket_t,next))), v_zero512, 1);
          state[k].m_have_tuple = _mm512_kandn(m_no_conflict, state[k].m_have_tuple);
          found += _mm_popcnt_u32(m_no_conflict);
        }
        if (m_full == state[k].m_have_tuple || cur >= rel->num_tuples) {
          state[k].stage = 0;
        } else {
          state[k].stage = 1;
        }
        --k;
      }
        break;
      case 0: {
        // insert new buckets after valid buckets
        m_to_insert = state[k].m_have_tuple;
        v_addr = _mm512_mask_blend_epi64(m_to_insert, v_all_ones, state[k].ht_off);
        v_conflict = _mm512_conflict_epi64(v_addr);
        m_no_conflict = _mm512_testn_epi64_mask(v_conflict, v_all_ones);
        m_no_conflict = _mm512_kand(m_no_conflict, m_to_insert);
        if (true) {
#if SERIAL_BUILD
          for (int i = 0; i < VECTOR_SCALE; ++i) {
            new_bucket[i] = 0;
            if (m_no_conflict & (1 << i)) {
              get_new_bucket(&bucket, overflowbuf);
              new_bucket[i] = bucket;
            }
          }
          v_next = _mm512_mask_i64gather_epi64(v_all_ones, state[k].m_have_tuple, _mm512_add_epi64(state[k].ht_off, v_next_off), 0, 1);
          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_new_bucket, v_key_off), state[k].key, 1);
          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_new_bucket, v_payload_off), state[k].payload, 1);
//          _mm512_mask_i64scatter_epi32(0, m_no_conflict, (v_new_bucket), v256_one, 1);
          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_new_bucket, v_next_off), v_next, 1);

          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(state[k].ht_off, v_next_off), v_new_bucket, 1);

          found += _mm_popcnt_u32(m_no_conflict);
          state[k].m_have_tuple = _mm512_kandn(m_no_conflict, state[k].m_have_tuple);
#else
          for (int i = 0; i < VECTOR_SCALE; ++i) {
            new_bucket[i] = 0;
            if (m_to_insert & (1 << i)) {
              get_new_bucket(&bucket, overflowbuf);
              new_bucket[i] = bucket;
            }
          }
          v_next = _mm512_mask_i64gather_epi64(v_all_ones, state[k].m_have_tuple, _mm512_add_epi64(state[k].ht_off, v_next_off), 0, 1);
          _mm512_mask_i64scatter_epi64(0, m_to_insert, _mm512_add_epi64(v_new_bucket, v_key_off), state[k].key, 1);
          _mm512_mask_i64scatter_epi64(0, m_to_insert, _mm512_add_epi64(v_new_bucket, v_payload_off), state[k].payload, 1);
//          _mm512_mask_i64scatter_epi32(0, m_to_insert, (v_new_bucket), v256_one, 1);
          _mm512_mask_i64scatter_epi64(0, m_to_insert, _mm512_add_epi64(v_new_bucket, v_next_off), v_next, 1);

          // conflict-solved insert
          v_lzeros = _mm512_lzcnt_epi64(v_conflict);
          v_lzeros = _mm512_sub_epi64(v_63, v_lzeros);
          to_scatt = _mm512_kandn(m_no_conflict, m_to_insert);
          v_previous = _mm512_maskz_permutexvar_epi64(to_scatt, v_lzeros, v_new_bucket);
          _mm512_mask_i64scatter_epi64(0, to_scatt, _mm512_add_epi64(v_previous, v_next_off), v_new_bucket, 1);

          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(state[k].ht_off, v_next_off), v_new_bucket, 1);
          found += _mm_popcnt_u32(m_to_insert);
          state[k].m_have_tuple = _mm512_kandn(m_to_insert, state[k].m_have_tuple);
#endif
        }
        if (m_full == state[k].m_have_tuple) {
          ht_pos = (uint64_t *) &state[k].ht_off;
          for (int i = 0; i < VECTOR_SCALE; ++i) {
            _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
          }
        } else {
          state[k].stage = 1;
          --k;
        }

      }
        break;
    }
    ++k;
  }

  return found;
}
size_t build_SIMD(hashtable_t *ht, relation_t *rel, bucket_buffer_t **overflowbuf) {
  int32_t found = 0, k = 0, done = 0, num, num_temp, new_add;
  __attribute__((aligned(64)))     __mmask8 to_scatt = 0, m_match = 0, m_new_cells = -1, m_valid_bucket = 0, mask[VECTOR_SCALE + 1], m_to_insert = 0, m_no_conflict, m_full = -1;
  __m512i v_offset = _mm512_set1_epi64(0), v_base_offset_upper = _mm512_set1_epi64(rel->num_tuples * sizeof(tuple_t)), v_base_offset, v_ht_cell, v_factor = _mm512_set1_epi64(
      ht->hash_mask), v_shift = _mm512_set1_epi64(ht->skip_bits), v_cell_hash, v_neg_one512 = _mm512_set1_epi64(-1), v_zero512 = _mm512_set1_epi64(0), v_write_index =
      _mm512_set1_epi64(0), v_ht_addr = _mm512_set1_epi64((uint64_t) ht->buckets), v_word_size = _mm512_set1_epi64(WORDSIZE), v_tuple_size = _mm512_set1_epi64(sizeof(tuple_t)),
      v_bucket_size = _mm512_set1_epi64(sizeof(bucket_t)), v_next_off = _mm512_set1_epi64(offsetof(bucket_t, next)), v_right_payload, v_addr, v_all_ones = _mm512_set1_epi64(-1),
      v_conflict, v_one = _mm512_set1_epi64(1), v_new_bucket, v_next, v_63 = _mm512_set1_epi64(63), v_key_off = _mm512_set1_epi64(offsetof(bucket_t, tuples[0].key)), v_lzeros,
      v_previous, v_payload_off = _mm512_set1_epi64(offsetof(bucket_t, tuples[0].payload)), v_count_off = _mm512_set1_epi64(offsetof(bucket_t, count));
  ;
  __m256i v256_one = _mm256_set1_epi32(1);
  tuple_t *join_res = NULL;
  uint64_t *pos = NULL, *new_bucket = (uint64_t*) &v_new_bucket;
  bucket_t * bucket;
  int tail_add = 0;
  __attribute__((aligned(64)))           uint64_t cur_offset = 0, base_off[16], *ht_pos;
  for (int i = 0; i <= VECTOR_SCALE; ++i) {
    base_off[i] = i * sizeof(tuple_t);
    mask[i] = (1 << i) - 1;
  }
  v_base_offset = _mm512_load_epi64(base_off);
  __attribute__((aligned(64)))          StateSIMD state[1];
  // init # of the state
  for (int i = 0; i < 1; ++i) {
    state[i].stage = 1;
    state[i].m_have_tuple = 0;
    state[i].ht_off = _mm512_set1_epi64(0);
    state[i].payload = _mm512_set1_epi64(0);
    state[i].key = _mm512_set1_epi64(0);
  }
  k = 0;
  for (uint64_t cur = 0; cur < rel->num_tuples || state[k].m_have_tuple;) {

    ///////// step 1: load new tuples' address offsets
    // the offset should be within MAX_32INT_
    // the tail depends on the number of joins and tuples in each bucket
#if SEQPREFETCH
    _mm_prefetch((char *)(((void *)rel->tuples) + cur_offset + PDIS), _MM_HINT_T0);
    _mm_prefetch((char *)(((void *)rel->tuples) + cur_offset + PDIS + 64), _MM_HINT_T0);
    _mm_prefetch((char *)(((void *)rel->tuples) + cur_offset + PDIS + 128), _MM_HINT_T0);
#endif
    // directly use cur, instead of cur_offset to control the offset to rel.
    // In this case, using step = 16 to gather data, but step is larger
    // than the scale 1,2,4 or 8
    v_offset = _mm512_add_epi64(_mm512_set1_epi64(cur_offset), v_base_offset);
    state[k].tb_off = _mm512_mask_expand_epi64(state[k].tb_off, _mm512_knot(state[k].m_have_tuple), v_offset);
    // count the number of empty tuples
    m_new_cells = _mm512_knot(state[k].m_have_tuple);
    new_add = _mm_popcnt_u32(m_new_cells);
    cur_offset = cur_offset + base_off[new_add];
    cur = cur + new_add;
    state[k].m_have_tuple = _mm512_cmpgt_epi64_mask(v_base_offset_upper, state[k].tb_off);
    ///// step 2: load new cells from right tuples;
    m_new_cells = _mm512_kand(m_new_cells, state[k].m_have_tuple);
    // maybe need offset within a tuple
    state[k].key = _mm512_mask_i64gather_epi64(state[k].key, m_new_cells, state[k].tb_off, ((void * )rel->tuples), 1);
    state[k].payload = _mm512_mask_i64gather_epi64(state[k].payload, state[k].m_have_tuple, _mm512_add_epi64(state[k].tb_off, v_word_size), ((void * )rel->tuples), 1);
    ///// step 3: load new values from hash tables;
    // hash the cell values
    v_cell_hash = _mm512_and_epi64(state[k].key, v_factor);
    v_cell_hash = _mm512_srlv_epi64(v_cell_hash, v_shift);
    v_cell_hash = _mm512_mullo_epi64(v_cell_hash, v_bucket_size);
    state[k].ht_off = _mm512_mask_add_epi64(state[k].ht_off, m_new_cells, v_cell_hash, v_ht_addr);

    /////////////////// random access
    // check valid bucket, insert directly for invalid buckets
    v_ht_cell = _mm512_mask_i64gather_epi64(v_neg_one512, state[k].m_have_tuple, state[k].ht_off, 0, 1);
    // inset new nodes
    m_to_insert = _mm512_cmpeq_epi64_mask(v_ht_cell, v_zero512);
    m_to_insert = _mm512_kand(m_to_insert, state[k].m_have_tuple);
    v_addr = _mm512_mask_blend_epi64(m_to_insert, v_all_ones, state[k].ht_off);
    v_conflict = _mm512_conflict_epi64(v_addr);
    m_no_conflict = _mm512_testn_epi64_mask(v_conflict, v_all_ones);
    m_no_conflict = _mm512_kand(m_no_conflict, m_to_insert);
    if (m_no_conflict) {
      // write the key , payload, count, next to the nodes
      _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_addr, v_key_off), state[k].key, 1);
      _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_addr, v_payload_off), state[k].payload, 1);
      _mm512_mask_i64scatter_epi32(0, m_no_conflict, _mm512_add_epi64(v_addr, v_count_off), v256_one, 1);
//          _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_addr,_mm512_set1_epi64(offsetof(bucket_t,next))), v_zero512, 1);
      state[k].m_have_tuple = _mm512_kandn(m_no_conflict, state[k].m_have_tuple);
      found += _mm_popcnt_u32(m_no_conflict);
    }

    // insert new buckets after valid buckets
    m_to_insert = state[k].m_have_tuple;
    v_addr = _mm512_mask_blend_epi64(m_to_insert, v_all_ones, state[k].ht_off);
    v_conflict = _mm512_conflict_epi64(v_addr);
    m_no_conflict = _mm512_testn_epi64_mask(v_conflict, v_all_ones);
    m_no_conflict = _mm512_kand(m_no_conflict, m_to_insert);
    if (true) {
#if SERIAL_BUILD
      for (int i = 0; i < VECTOR_SCALE; ++i) {
        new_bucket[i] = 0;
        if (m_no_conflict & (1 << i)) {
          get_new_bucket(&bucket, overflowbuf);
          new_bucket[i] = bucket;
        }
      }
      v_next = _mm512_mask_i64gather_epi64(v_all_ones, state[k].m_have_tuple, _mm512_add_epi64(state[k].ht_off, v_next_off), 0, 1);
      _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_new_bucket, v_key_off), state[k].key, 1);
      _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_new_bucket, v_payload_off), state[k].payload, 1);
//          _mm512_mask_i64scatter_epi32(0, m_no_conflict, (v_new_bucket), v256_one, 1);
      _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(v_new_bucket, v_next_off), v_next, 1);

      _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(state[k].ht_off, v_next_off), v_new_bucket, 1);

      found += _mm_popcnt_u32(m_no_conflict);
      state[k].m_have_tuple = _mm512_kandn(m_no_conflict, state[k].m_have_tuple);
#else
      for (int i = 0; i < VECTOR_SCALE; ++i) {
        new_bucket[i] = 0;
        if (m_to_insert & (1 << i)) {
          get_new_bucket(&bucket, overflowbuf);
          new_bucket[i] = bucket;
        }
      }
      v_next = _mm512_mask_i64gather_epi64(v_all_ones, state[k].m_have_tuple, _mm512_add_epi64(state[k].ht_off, v_next_off), 0, 1);
      _mm512_mask_i64scatter_epi64(0, m_to_insert, _mm512_add_epi64(v_new_bucket, v_key_off), state[k].key, 1);
      _mm512_mask_i64scatter_epi64(0, m_to_insert, _mm512_add_epi64(v_new_bucket, v_payload_off), state[k].payload, 1);
//          _mm512_mask_i64scatter_epi32(0, m_to_insert, (v_new_bucket), v256_one, 1);
      _mm512_mask_i64scatter_epi64(0, m_to_insert, _mm512_add_epi64(v_new_bucket, v_next_off), v_next, 1);

      // conflict-solved insert
      v_lzeros = _mm512_lzcnt_epi64(v_conflict);
      v_lzeros = _mm512_sub_epi64(v_63, v_lzeros);
      to_scatt = _mm512_kandn(m_no_conflict, m_to_insert);
      v_previous = _mm512_maskz_permutexvar_epi64(to_scatt, v_lzeros, v_new_bucket);
      _mm512_mask_i64scatter_epi64(0, to_scatt, _mm512_add_epi64(v_previous, v_next_off), v_new_bucket, 1);

      _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(state[k].ht_off, v_next_off), v_new_bucket, 1);
      found += _mm_popcnt_u32(m_to_insert);
      state[k].m_have_tuple = _mm512_kandn(m_to_insert, state[k].m_have_tuple);
#endif
    }
  }

  return found;
}

volatile static char g_lock = 0, g_lock_morse = 0;
volatile static uint64_t total_num = 0, global_curse = 0, global_upper, thread_num = 1,global_morse_size;
typedef int64_t (*BuildFun)(hashtable_t *ht, relation_t *rel, bucket_buffer_t **overflowbuf);
volatile static struct Fun {
  BuildFun fun_ptr;
  char fun_name[8];
} pfun[10];
volatile static int pf_num = 0;
static map<int64_t, int64_t> len2num;
static void search_ht(hashtable_t *ht) {
  int64_t len = 0;
  for (int64_t i = 0; i < ht->num_buckets; ++i) {
    len = 0;
    for (auto it = ht->buckets + i; it; it = it->next) {
      len++;
    }
    len2num[len]++;
  }
}
static void morse_driven(void*param, BuildFun fun, bucket_buffer_t **overflowbuf) {
  arg_t *args = (arg_t *) param;
  uint64_t base = 0, num = 0;
  args->num_results = 0;
  relation_t relS;
  relS.tuples = args->relS.tuples;
  relS.num_tuples = 0;
  while (1) {
    lock(&g_lock_morse);
    base = global_curse;
    global_curse += global_morse_size;
    unlock(&g_lock_morse);
    if (base >= global_upper) {
      break;
    }
    num = (global_upper - base) < global_morse_size ? (global_upper - base) : global_morse_size;
    relS.tuples = args->relS.tuples + base;
    relS.num_tuples = num;
    args->num_results += fun(args->ht, &relS, overflowbuf);
  }
}

void *build_thread(void *param) {
  int rv;
  arg_t *args = (arg_t *) param;
  struct timeval t1, t2;
  int deltaT = 0;
  bucket_buffer_t *overflowbuf;
  hashtable_t *ht;
  uint32_t nbuckets = (args->relR.num_tuples / BUCKET_SIZE / thread_num);
  if (args->tid == 0) {
    strcpy(pfun[5].fun_name, "IMV");
    strcpy(pfun[4].fun_name, "AMAC");
    strcpy(pfun[3].fun_name, "FVA");
    strcpy(pfun[2].fun_name, "DVA");
    strcpy(pfun[1].fun_name, "SIMD");
    strcpy(pfun[0].fun_name, "Naive");

    pfun[5].fun_ptr = build_imv;
    pfun[4].fun_ptr = build_AMAC;
    pfun[3].fun_ptr = build_FVA;
    pfun[2].fun_ptr = build_DVA;
    pfun[1].fun_ptr = build_SIMD;
    pfun[0].fun_ptr = build_raw;

    pf_num = 6;
  }
  BARRIER_ARRIVE(args->barrier, rv);

  for (int fid = 0; fid < pf_num; ++fid) {
    for (int rp = 0; rp < REPEAT_PROBE; ++rp) {
      init_bucket_buffer(&overflowbuf);
      allocate_hashtable(&ht, nbuckets);
      BARRIER_ARRIVE(args->barrier, rv);
      gettimeofday(&t1, NULL);

#if MORSE_SIZE
      args->ht = ht;
      morse_driven(param, pfun[fid].fun_ptr, &overflowbuf);
#else
      args->num_results = pfun[fid].fun_ptr(ht, &args->relS, &overflowbuf);
#endif
      lock(&g_lock);
#if DIVIDE
      total_num += args->num_results;
#elif MORSE_SIZE
      total_num += args->num_results;
#else
      total_num = args->num_results;
#endif
#if PRINT_HT
      search_ht(ht);
#endif
      unlock(&g_lock);
      BARRIER_ARRIVE(args->barrier, rv);
      if (args->tid == 0) {
        gettimeofday(&t2, NULL);
        printf("total result num = %lld\t", total_num);
        deltaT = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
        printf("---- %5s BUILD costs time (ms) = %10.4lf\n", pfun[fid].fun_name, deltaT * 1.0 / 1000);
        total_num = 0;
        global_curse = 0;
        for (auto iter : len2num) {
          cout << "ht  cnum = " << iter.first << " , times = " << iter.second << endl;
          if (total_num++ > 20)
            break;
        }
        total_num = 0;
        len2num.clear();
      }
      destroy_hashtable(ht);
      free_bucket_buffer(overflowbuf);
    }
  }
  return nullptr;
}

result_t *BUILD(relation_t *relR, relation_t *relS, int nthreads) {
  hashtable_t *ht;
  int64_t result = 0;
  int32_t numR, numS, numRthr, numSthr; /* total and per thread num */
  int i, rv;
  cpu_set_t set;
  arg_t args[nthreads];
  pthread_t tid[nthreads];
  pthread_attr_t attr;
  pthread_barrier_t barrier;

  result_t *joinresult = 0;
  joinresult = (result_t *) malloc(sizeof(result_t));

#ifdef JOIN_RESULT_MATERIALIZE
  joinresult->resultlist = (threadresult_t *) alloc_aligned(sizeof(threadresult_t) * nthreads);
#endif
//  thread_num = nthreads;

#if USE_TBB && 0
  pthread_attr_init(&attr);
  for (i = 0; i < nthreads; i++) {
    int cpu_idx = get_cpu_id(i);

    DEBUGMSG(1, "Assigning thread-%d to CPU-%d\n", i, cpu_idx);
#if AFFINITY
    CPU_ZERO(&set);
    CPU_SET(cpu_idx, &set);
    pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);
#endif
  }

  tbb_run(relR, relS, nthreads);
  joinresult->totalresults = result;
  joinresult->nthreads = nthreads;
  return joinresult;
#endif
  uint32_t nbuckets = (relS->num_tuples / BUCKET_SIZE / thread_num);
  allocate_hashtable(&ht, nbuckets);

  numR = relR->num_tuples;
  numS = relS->num_tuples;
  numRthr = numR / nthreads;
  numSthr = numS / nthreads;

  rv = pthread_barrier_init(&barrier, NULL, nthreads);
  if (rv != 0) {
    printf("Couldn't create the barrier\n");
    exit(EXIT_FAILURE);
  }
  global_curse = 0;
  global_upper = relS->num_tuples;
  if(nthreads==1){
    global_morse_size= relS->num_tuples;
  }else{
    global_morse_size = MORSE_SIZE;
  }
  pthread_attr_init(&attr);
  for (i = 0; i < nthreads; i++) {
    int cpu_idx = get_cpu_id(i);

    DEBUGMSG(1, "Assigning thread-%d to CPU-%d\n", i, cpu_idx);
#if AFFINITY
    CPU_ZERO(&set);
    CPU_SET(cpu_idx, &set);
    pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);
#endif
    args[i].tid = i;
    args[i].ht = ht;
    args[i].barrier = &barrier;
#if DIVIDE
    /* assing part of the relR for next thread */
    args[i].relR.num_tuples = (i == (nthreads - 1)) ? numR : numRthr;
    args[i].relR.tuples = relR->tuples + numRthr * i;
    numR -= numRthr;
#else
    args[i].relR.num_tuples = relR->num_tuples;
    args[i].relR.tuples = relR->tuples;
#endif
#if DIVIDE
    /* assing part of the relS for next thread */
    args[i].relS.num_tuples = (i == (nthreads - 1)) ? numS : numSthr;
    args[i].relS.tuples = relS->tuples + numSthr * i;
    numS -= numSthr;
#else
    args[i].relS.num_tuples = relS->num_tuples;
    args[i].relS.tuples = relS->tuples;
#endif
    args[i].threadresult = &(joinresult->resultlist[i]);

    rv = pthread_create(&tid[i], &attr, build_thread, (void *) &args[i]);
    if (rv) {
      printf("ERROR; return code from pthread_create() is %d\n", rv);
      exit(-1);
    }
  }

  for (i = 0; i < nthreads; i++) {
    pthread_join(tid[i], NULL);
    /* sum up results */
#if DIVIDE
    result += args[i].num_results;
#else
    result = args[i].num_results;
#endif
  }
  joinresult->totalresults = result;
  joinresult->nthreads = nthreads;

  destroy_hashtable(ht);

  return joinresult;
}

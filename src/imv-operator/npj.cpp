#ifndef _SIMD_PREFETCHING
#define _SIMD_PREFETCHING

#include "imv-operator/npj.hpp"

#include "imv-operator/tuple_buffer.hpp"
#define WORDSIZE 8
// target for 8B keys and 8B payload
int64_t probe_simd(hashtable_t *ht, relation_t *rel, void *output) {
  int64_t matches = 0;
  int32_t new_add = 0;
  __mmask8 m_match = 0, m_have_tuple = 0, m_new_cells = -1, m_valid_bucket = 0;
  __m512i v_offset = _mm512_set1_epi64(0), v_addr_offset = _mm512_set1_epi64(0), v_base_offset_upper = _mm512_set1_epi64(rel->num_tuples * sizeof(tuple_t)), v_tuple_cell =
      _mm512_set1_epi64(0), v_base_offset, v_ht_cell, v_factor = _mm512_set1_epi64(ht->hash_mask), v_shift = _mm512_set1_epi64(ht->skip_bits), v_cell_hash, v_ht_pos =
      _mm512_set1_epi64(0), v_neg_one512 = _mm512_set1_epi64(-1), v_zero512 = _mm512_set1_epi64(0), v_next_addr = _mm512_set1_epi64(0), v_write_index = _mm512_set1_epi64(0),
      v_ht_addr = _mm512_set1_epi64((uint64_t)(uint64_t)ht->buckets), v_word_size = _mm512_set1_epi64(WORDSIZE), v_tuple_size = _mm512_set1_epi64(sizeof(tuple_t)), v_bucket_size = _mm512_set1_epi64(
          sizeof(bucket_t)), v_next_off = _mm512_set1_epi64(8), v_left_payload, v_right_payload, v_payload_off = _mm512_set1_epi64(24);
  chainedtuplebuffer_t *chainedbuf = (chainedtuplebuffer_t *) output;
  tuple_t *join_res = NULL;
  uint64_t cur_offset = 0, base_off[16];

  for (int i = 0; i <= VECTOR_SCALE; ++i) {
    base_off[i] = i * sizeof(tuple_t);
  }
  v_base_offset = _mm512_load_epi64(base_off);
  for (uint64_t cur = 0; cur < rel->num_tuples || m_have_tuple;) {
    ///////// step 1: load new tuples' address offsets
    // the offset should be within MAX_32INT_
    // the tail depends on the number of joins and tuples in each bucket
    /*#if SEQPREFETCH
     _mm_prefetch((char *)(((void *)rel->tuples) + cur_offset + PDIS),
     _MM_HINT_T0);
     _mm_prefetch((char *)(((void *)rel->tuples) + cur_offset + PDIS + 64),
     _MM_HINT_T0);
     _mm_prefetch((char *)(((void *)rel->tuples) + cur_offset + PDIS + 128),
     _MM_HINT_T0);
     _mm_prefetch((char *)(((void *)rel->tuples) + cur_offset + PDIS + 192),
     _MM_HINT_T0);
     _mm_prefetch((char *)(((void *)rel->tuples) + cur_offset + PDIS + 256),
     _MM_HINT_T0);
     #endif
     */
    // directly use cur, instead of cur_offset to control the offset to rel.
    // In this case, using step = 16 to gather data, but step is larger
    // than the scale 1,2,4 or 8
    v_offset = _mm512_add_epi64(_mm512_set1_epi64(cur_offset), v_base_offset);
    v_addr_offset = _mm512_mask_expand_epi64(v_addr_offset, _mm512_knot(m_have_tuple), v_offset);
    // count the number of empty tuples
    m_new_cells = _mm512_knot(m_have_tuple);
    new_add = _mm_popcnt_u32(m_new_cells);
    cur_offset = cur_offset + base_off[new_add];
    cur = cur + new_add;
    m_have_tuple = _mm512_cmpgt_epi64_mask(v_base_offset_upper, v_addr_offset);
    ///// step 2: load new cells from right tuples;
    m_new_cells = _mm512_kand(m_new_cells, m_have_tuple);
    // maybe need offset within a tuple
    v_tuple_cell = _mm512_mask_i64gather_epi64(v_tuple_cell, m_new_cells, v_addr_offset, ((void * )rel->tuples), 1);
    ///// step 3: load new values from hash tables;
    // hash the cell values
    v_cell_hash = _mm512_and_epi64(v_tuple_cell, v_factor);
    v_cell_hash = _mm512_srlv_epi64(v_cell_hash, v_shift);
    v_cell_hash = _mm512_mullo_epi64(v_cell_hash, v_bucket_size);
    v_ht_pos = _mm512_mask_add_epi64(v_next_addr, m_new_cells, v_cell_hash, v_ht_addr);

    /////////////////// random access
    // check valid bucket
    v_ht_cell = _mm512_mask_i64gather_epi64(v_neg_one512, m_have_tuple, v_ht_pos, 0, 1);
    m_valid_bucket = _mm512_cmpneq_epi64_mask(v_ht_cell, v_zero512);
    m_have_tuple = _mm512_kand(m_valid_bucket, m_have_tuple);

    v_ht_cell = _mm512_mask_i64gather_epi64(v_neg_one512, m_have_tuple, _mm512_add_epi64(v_ht_pos, v_tuple_size), 0, 1);  // note the offset of the tuple in %bucket_t%

    ///// step 4: compare;
    m_match = _mm512_cmpeq_epi64_mask(v_tuple_cell, v_ht_cell);
    m_match = _mm512_kand(m_match, m_have_tuple);
    new_add = _mm_popcnt_u32(m_match);
    matches += new_add;
    v_next_addr = _mm512_mask_i64gather_epi64(v_zero512, m_have_tuple, _mm512_add_epi64(v_ht_pos, v_next_off), 0, 1);
    m_have_tuple = _mm512_kand(_mm512_cmpneq_epi64_mask(v_next_addr, v_zero512), m_have_tuple);
    // gather payloads
    v_left_payload = _mm512_mask_i64gather_epi64(v_neg_one512, m_match, _mm512_add_epi64(v_addr_offset, v_word_size), ((void * )rel->tuples), 1);
    v_right_payload = _mm512_mask_i64gather_epi64(v_neg_one512, m_match, _mm512_add_epi64(v_ht_pos, v_payload_off), 0, 1);

    // to scatter join results
    join_res = cb_next_n_writepos(chainedbuf, new_add);
    v_write_index = _mm512_mask_expand_epi64(v_zero512, m_match, v_base_offset);
    _mm512_mask_i64scatter_epi64((void * )join_res, m_match, v_write_index, v_left_payload, 1);
    v_write_index = _mm512_add_epi64(v_write_index, v_word_size);
    _mm512_mask_i64scatter_epi64((void * )join_res, m_match, v_write_index, v_right_payload, 1);
  }
  return matches;
}
int64_t probe_simd_amac(hashtable_t *ht, relation_t *rel, void *output) {
  int64_t matches = 0;
  int32_t new_add = 0, k = 0, done = 0;
  __mmask8 m_match = 0, m_new_cells = -1, m_valid_bucket = 0,m_full=-1;
  __m512i v_offset = _mm512_set1_epi64(0), v_base_offset_upper = _mm512_set1_epi64(rel->num_tuples * sizeof(tuple_t)), v_base_offset, v_ht_cell, v_factor = _mm512_set1_epi64(
      ht->hash_mask), v_shift = _mm512_set1_epi64(ht->skip_bits), v_cell_hash, v_neg_one512 = _mm512_set1_epi64(-1), v_zero512 = _mm512_set1_epi64(0), v_write_index =
      _mm512_set1_epi64(0), v_ht_addr = _mm512_set1_epi64((uint64_t)ht->buckets), v_word_size = _mm512_set1_epi64(WORDSIZE), v_tuple_size = _mm512_set1_epi64(sizeof(tuple_t)),
      v_bucket_size = _mm512_set1_epi64(sizeof(bucket_t)), v_next_off = _mm512_set1_epi64(8), v_left_payload, v_right_payload, v_payload_off = _mm512_set1_epi64(24);
  chainedtuplebuffer_t *chainedbuf = (chainedtuplebuffer_t *) output;
  tuple_t *join_res = NULL;
  __attribute__((aligned(64)))  uint64_t cur_offset = 0, base_off[16], *ht_pos;

  for (int i = 0; i <= VECTOR_SCALE; ++i) {
    base_off[i] = i * sizeof(tuple_t);
  }
  v_base_offset = _mm512_load_epi64(base_off);
  StateSIMD state[SIMDStateSize];
  // init # of the state
  for (int i = 0; i < SIMDStateSize; ++i) {
    state[i].stage = 1;
    state[i].m_have_tuple = 0;
    state[i].ht_off = _mm512_set1_epi64(0);
    state[i].tb_off = _mm512_set1_epi64(0);
    state[i].payload = _mm512_set1_epi64(0);
    state[i].key = _mm512_set1_epi64(0);
  }
  for (uint64_t cur = 0; (cur < rel->num_tuples) || (done < SIMDStateSize);) {
    k = (k >= SIMDStateSize) ? 0 : k;
    if (cur >= rel->num_tuples) {
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
        state[k].stage = 0;
#if KNL
        _mm512_mask_prefetch_i64gather_pd(
            state[k].ht_off, state[k].m_have_tuple, 0, 1, _MM_HINT_T0);
#elif DIR_PREFETCH
        ht_pos = (uint64_t *) &state[k].ht_off;
        for (int i = 0; i < VECTOR_SCALE; ++i) {
          _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
        }
#else
        m_have_tuple = state[k].m_have_tuple;
        ht_pos = (uint64_t *)&state[k].ht_off;
        for (int i = 0; (i < VECTOR_SCALE) & (m_have_tuple);
            ++i, (m_have_tuple >> 1)) {
          if (m_have_tuple & 1) {
            _mm_prefetch((char *)(ht_pos[i]), _MM_HINT_T0);
          }
        }
#endif
      }
        break;
      case 0: {
        /////////////////// random access
        // check valid bucket
        v_ht_cell = _mm512_mask_i64gather_epi64(v_neg_one512, state[k].m_have_tuple, state[k].ht_off, 0, 1);
        m_valid_bucket = _mm512_cmpneq_epi64_mask(v_ht_cell, v_zero512);
        state[k].m_have_tuple = _mm512_kand(m_valid_bucket, state[k].m_have_tuple);

        v_ht_cell = _mm512_mask_i64gather_epi64(v_neg_one512, state[k].m_have_tuple, _mm512_add_epi64(state[k].ht_off, v_tuple_size), 0, 1);  // note the offset of the tuple in %bucket_t%

        ///// step 4: compare;
        m_match = _mm512_cmpeq_epi64_mask(state[k].key, v_ht_cell);
        m_match = _mm512_kand(m_match, state[k].m_have_tuple);
        new_add = _mm_popcnt_u32(m_match);
        matches += new_add;

        // gather payloads
        v_right_payload = _mm512_mask_i64gather_epi64(v_neg_one512, m_match, _mm512_add_epi64(state[k].ht_off, v_payload_off), 0, 1);

        // update next
        state[k].ht_off = _mm512_mask_i64gather_epi64(v_zero512, state[k].m_have_tuple, _mm512_add_epi64(state[k].ht_off, v_next_off), 0, 1);
        state[k].m_have_tuple = _mm512_kand(_mm512_cmpneq_epi64_mask(state[k].ht_off, v_zero512), state[k].m_have_tuple);

        // to scatter join results
        join_res = cb_next_n_writepos(chainedbuf, new_add);
#if SEQPREFETCH
        _mm_prefetch((char *)(((void *)join_res) + PDIS), _MM_HINT_T0);
        _mm_prefetch((char *)(((void *)join_res) + PDIS + 64), _MM_HINT_T0);
        _mm_prefetch((char *)(((void *)join_res) + PDIS + 128), _MM_HINT_T0);
#endif
        v_write_index = _mm512_mask_expand_epi64(v_zero512, m_match, v_base_offset);
        _mm512_mask_i64scatter_epi64((void * )join_res, m_match, v_write_index, state[k].payload, 1);
        v_write_index = _mm512_add_epi64(v_write_index, v_word_size);
        _mm512_mask_i64scatter_epi64((void * )join_res, m_match, v_write_index, v_right_payload, 1);
        // return back every time
#if 0
        if (m_full == state[k].m_have_tuple) {
          ht_pos = (uint64_t *) &state[k].ht_off;
          for (int i = 0; i < VECTOR_SCALE; ++i) {
            _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
          }
        } else {
          state[k].stage = 1;
          --k;
        }
#else
        state[k].stage = 1;
        --k;
#endif
      }
        break;
    }
    ++k;
  }
  return matches;
}
int64_t probe_simd_amac_raw(hashtable_t *ht, relation_t *rel, void *output) {
  int64_t matches = 0;
  int32_t new_add = 0, k = 0, done = 0;
  __mmask8 m_match = 0, m_new_cells = -1, m_valid_bucket = 0;
  __m512i v_offset = _mm512_set1_epi64(0), v_base_offset_upper = _mm512_set1_epi64(rel->num_tuples * sizeof(tuple_t)), v_base_offset, v_left_size = _mm512_set1_epi64(8),
      v_bucket_offset = _mm512_set1_epi64(0), v_ht_cell, v_factor = _mm512_set1_epi64(ht->hash_mask), v_shift = _mm512_set1_epi64(ht->skip_bits), v_cell_hash, v_neg_one512 =
          _mm512_set1_epi64(-1), v_ht_upper, v_zero512 = _mm512_set1_epi64(0), v_write_index = _mm512_set1_epi64(0), v_ht_addr = _mm512_set1_epi64((uint64_t)ht->buckets), v_word_size =
          _mm512_set1_epi64(WORDSIZE), v_tuple_size = _mm512_set1_epi64(sizeof(tuple_t)), v_bucket_size = _mm512_set1_epi64(sizeof(bucket_t)), v_next_off = _mm512_set1_epi64(8),
      v_left_payload, v_right_payload, v_payload_off = _mm512_set1_epi64(24);
  chainedtuplebuffer_t *chainedbuf = (chainedtuplebuffer_t *) output;
  tuple_t *join_res = NULL;
  __attribute__((aligned(64)))  uint64_t cur_offset = 0, base_off[16], *left_payload, *right_payload, *ht_pos;
  left_payload = (uint64_t *) &v_left_payload;
  right_payload = (uint64_t *) &v_right_payload;

  for (int i = 0; i <= VECTOR_SCALE; ++i) {
    base_off[i] = i * sizeof(tuple_t);
  }
  v_base_offset = _mm512_load_epi64(base_off);
  StateSIMD state[SIMDStateSize];
  // init # of the state
  for (int i = 0; i < SIMDStateSize; ++i) {
    state[i].stage = 1;
    state[i].m_have_tuple = 0;
    state[i].ht_off = _mm512_set1_epi64(0);
    state[i].tb_off = _mm512_set1_epi64(0);
    state[i].key = _mm512_set1_epi64(0);
  }
  for (uint64_t cur = 0; (cur < rel->num_tuples) || (done < SIMDStateSize);) {
    k = (k >= SIMDStateSize) ? 0 : k;
    if (cur >= rel->num_tuples) {
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
        state[k].stage = 0;
#if KNL
        _mm512_mask_prefetch_i64gather_pd(
            state[k].ht_off, state[k].m_have_tuple, 0, 1, _MM_HINT_T0);
#elif DIR_PREFETCH
        ht_pos = (uint64_t *) &state[k].ht_off;
        for (int i = 0; i < VECTOR_SCALE; ++i) {
          _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
        }
#else
        m_have_tuple = state[k].m_have_tuple;
        ht_pos = (uint64_t *)&state[k].ht_off;
        for (int i = 0; (i < VECTOR_SCALE) & (m_have_tuple);
            ++i, (m_have_tuple >> 1)) {
          if (m_have_tuple & 1) {
            _mm_prefetch((char *)(ht_pos[i]), _MM_HINT_T0);
          }
        }
#endif
      }
        break;
      case 0: {
        /////////////////// random access
        // check valid bucket
        v_ht_cell = _mm512_mask_i64gather_epi64(v_neg_one512, state[k].m_have_tuple, state[k].ht_off, 0, 1);
        m_valid_bucket = _mm512_cmpneq_epi64_mask(v_ht_cell, v_zero512);
        state[k].m_have_tuple = _mm512_kand(m_valid_bucket, state[k].m_have_tuple);

        v_ht_cell = _mm512_mask_i64gather_epi64(v_neg_one512, state[k].m_have_tuple, _mm512_add_epi64(state[k].ht_off, v_tuple_size), 0, 1);  // note the offset of the tuple in %bucket_t%

        ///// step 4: compare;
        m_match = _mm512_cmpeq_epi64_mask(state[k].key, v_ht_cell);
        m_match = _mm512_kand(m_match, state[k].m_have_tuple);
        new_add = _mm_popcnt_u32(m_match);
        matches += new_add;

        // gather payloads
        v_right_payload = _mm512_mask_i64gather_epi64(v_neg_one512, m_match, _mm512_add_epi64(state[k].ht_off, v_payload_off), 0, 1);

        // update next
        state[k].ht_off = _mm512_mask_i64gather_epi64(v_zero512, state[k].m_have_tuple, _mm512_add_epi64(state[k].ht_off, v_next_off), 0, 1);
        state[k].m_have_tuple = _mm512_kand(_mm512_cmpneq_epi64_mask(state[k].ht_off, v_zero512), state[k].m_have_tuple);

        // to scatter join results
        join_res = cb_next_n_writepos(chainedbuf, new_add);
#if SEQPREFETCH
        _mm_prefetch((char *)(((void *)join_res) + PDIS), _MM_HINT_T0);
        _mm_prefetch((char *)(((void *)join_res) + PDIS + 64), _MM_HINT_T0);
        _mm_prefetch((char *)(((void *)join_res) + PDIS + 128), _MM_HINT_T0);
#endif
        v_write_index = _mm512_mask_expand_epi64(v_zero512, m_match, v_base_offset);
        _mm512_mask_i64scatter_epi64((void * )join_res, m_match, v_write_index, state[k].payload, 1);
        v_write_index = _mm512_add_epi64(v_write_index, v_word_size);
        _mm512_mask_i64scatter_epi64((void * )join_res, m_match, v_write_index, v_right_payload, 1);
        new_add = _mm_popcnt_u32(state[k].m_have_tuple);
        if (new_add == 0) {
          state[k].stage = 1;
        } else {
#if KNL
          _mm512_mask_prefetch_i64gather_pd(
              state[k].ht_off, state[k].m_have_tuple, 0, 1, _MM_HINT_T0);
#elif DIR_PREFETCH
          ht_pos = (uint64_t *) &state[k].ht_off;
          for (int i = 0; i < VECTOR_SCALE; ++i) {
            _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
          }
#else
          m_have_tuple = state[k].m_have_tuple;
          ht_pos = (uint64_t *)&state[k].ht_off;
          for (int i = 0; (i < VECTOR_SCALE) & (m_have_tuple);
              ++i, (m_have_tuple >> 1)) {
            if (m_have_tuple & 1) {
              _mm_prefetch((char *)(ht_pos[i]), _MM_HINT_T0);
            }
          }
#endif
        }
      }
        break;
    }
    ++k;
  }
  return matches;
}

int64_t probe_simd_amac_compact(hashtable_t *ht, relation_t *rel, void *output) {
  int64_t matches = 0;
  int32_t new_add = 0, k = 0, done = 0, num, num_temp;
  __mmask8 m_match = 0, m_new_cells = -1, m_valid_bucket = 0, mask[VECTOR_SCALE + 1];
  __m512i v_offset = _mm512_set1_epi64(0), v_base_offset_upper = _mm512_set1_epi64(rel->num_tuples * sizeof(tuple_t)), v_base_offset, v_left_size = _mm512_set1_epi64(8),
      v_bucket_offset = _mm512_set1_epi64(0), v_ht_cell, v_factor = _mm512_set1_epi64(ht->hash_mask), v_shift = _mm512_set1_epi64(ht->skip_bits), v_cell_hash, v_neg_one512 =
          _mm512_set1_epi64(-1), v_ht_upper, v_zero512 = _mm512_set1_epi64(0), v_write_index = _mm512_set1_epi64(0), v_ht_addr = _mm512_set1_epi64((uint64_t)ht->buckets), v_word_size =
          _mm512_set1_epi64(WORDSIZE), v_tuple_size = _mm512_set1_epi64(sizeof(tuple_t)), v_bucket_size = _mm512_set1_epi64(sizeof(bucket_t)), v_next_off = _mm512_set1_epi64(8),
      v_left_payload, v_right_payload, v_payload_off = _mm512_set1_epi64(24);
  chainedtuplebuffer_t *chainedbuf = (chainedtuplebuffer_t *) output;
  tuple_t *join_res = NULL;
  __attribute__((aligned(64)))  uint64_t cur_offset = 0, base_off[16], *left_payload, *right_payload, *ht_pos;
  left_payload = (uint64_t *) &v_left_payload;
  right_payload = (uint64_t *) &v_right_payload;

  for (int i = 0; i <= VECTOR_SCALE; ++i) {
    base_off[i] = i * sizeof(tuple_t);
    mask[i] = (1 << i) - 1;
  }
  v_base_offset = _mm512_load_epi64(base_off);
  StateSIMD state[SIMDStateSize + 1];
  // init # of the state
  for (int i = 0; i <= SIMDStateSize; ++i) {
    state[i].stage = 1;
    state[i].m_have_tuple = 0;
    state[i].ht_off = _mm512_set1_epi64(0);
    state[i].tb_off = _mm512_set1_epi64(0);
    state[i].key = _mm512_set1_epi64(0);
  }
  for (uint64_t cur = 0; 1;) {
    k = (k >= SIMDStateSize) ? 0 : k;
    if (cur >= rel->num_tuples) {
      if (state[k].m_have_tuple == 0 && state[k].stage != 3) {
        ++done;
        state[k].stage = 3;
        ++k;
        continue;
      }
      if (done >= SIMDStateSize) {
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
        state[k].stage = 0;
#if KNL
        _mm512_mask_prefetch_i64gather_pd(
            state[k].ht_off, state[k].m_have_tuple, 0, 1, _MM_HINT_T0);
#elif DIR_PREFETCH
        ht_pos = (uint64_t *) &state[k].ht_off;
        for (int i = 0; i < VECTOR_SCALE; ++i) {
          _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
        }
#else
        m_have_tuple = state[k].m_have_tuple;
        ht_pos = (uint64_t *)&state[k].ht_off;
        for (int i = 0; (i < VECTOR_SCALE) & (m_have_tuple);
            ++i, (m_have_tuple >> 1)) {
          if (m_have_tuple & 1) {
            _mm_prefetch((char *)(ht_pos[i]), _MM_HINT_T0);
          }
        }
#endif
      }
        break;
      case 0: {
        /////////////////// random access
        // check valid bucket
        v_ht_cell = _mm512_mask_i64gather_epi64(v_neg_one512, state[k].m_have_tuple, state[k].ht_off, 0, 1);
        m_valid_bucket = _mm512_cmpneq_epi64_mask(v_ht_cell, v_zero512);
        state[k].m_have_tuple = _mm512_kand(m_valid_bucket, state[k].m_have_tuple);

        v_ht_cell = _mm512_mask_i64gather_epi64(v_neg_one512, state[k].m_have_tuple, _mm512_add_epi64(state[k].ht_off, v_tuple_size), 0, 1);  // note the offset of the tuple in %bucket_t%

        ///// step 4: compare;
        m_match = _mm512_cmpeq_epi64_mask(state[k].key, v_ht_cell);
        m_match = _mm512_kand(m_match, state[k].m_have_tuple);
        new_add = _mm_popcnt_u32(m_match);
        matches += new_add;

        // gather payloads
        v_right_payload = _mm512_mask_i64gather_epi64(v_neg_one512, m_match, _mm512_add_epi64(state[k].ht_off, v_payload_off), 0, 1);

        // update next
        state[k].ht_off = _mm512_mask_i64gather_epi64(v_zero512, state[k].m_have_tuple, _mm512_add_epi64(state[k].ht_off, v_next_off), 0, 1);
        state[k].m_have_tuple = _mm512_kand(_mm512_cmpneq_epi64_mask(state[k].ht_off, v_zero512), state[k].m_have_tuple);

        // to scatter join results
        join_res = cb_next_n_writepos(chainedbuf, new_add);
#if SEQPREFETCH
        _mm_prefetch((char *)(((void *)join_res) + PDIS), _MM_HINT_T0);
        _mm_prefetch((char *)(((void *)join_res) + PDIS + 64), _MM_HINT_T0);
        _mm_prefetch((char *)(((void *)join_res) + PDIS + 128), _MM_HINT_T0);
#endif
        v_write_index = _mm512_mask_expand_epi64(v_zero512, m_match, v_base_offset);
        _mm512_mask_i64scatter_epi64((void * )join_res, m_match, v_write_index, state[k].payload, 1);
        v_write_index = _mm512_add_epi64(v_write_index, v_word_size);
        _mm512_mask_i64scatter_epi64((void * )join_res, m_match, v_write_index, v_right_payload, 1);
        num = _mm_popcnt_u32(state[k].m_have_tuple);
        if ((num == VECTOR_SCALE) || (done >= SIMDStateSize)) {
          state[k].stage = 0;
#if KNL
          _mm512_mask_prefetch_i64gather_pd(
              state[k].ht_off, state[k].m_have_tuple, 0, 1, _MM_HINT_T0);
#elif DIR_PREFETCH
          ht_pos = (uint64_t *) &state[k].ht_off;
          for (int i = 0; i < VECTOR_SCALE; ++i) {
            _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
          }
#else
          m_have_tuple = state[k].m_have_tuple;
          ht_pos = (uint64_t *)&state[k].ht_off;
          for (int i = 0; (i < VECTOR_SCALE) & (m_have_tuple);
              ++i, (m_have_tuple >> 1)) {
            if (m_have_tuple & 1) {
              _mm_prefetch((char *)(ht_pos[i]), _MM_HINT_T0);
            }
          }
#endif
        } else {
          if ((done < SIMDStateSize)) {
            state[k].stage = 2;
          } else {
            assert(0);
          }
          --k;
        }
      }
        break;
      case 2: {
        num = _mm_popcnt_u32(state[k].m_have_tuple);
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
#elif DIR_PREFETCH
          ht_pos = (uint64_t *) &state[k].ht_off;
          for (int i = 0; i < VECTOR_SCALE; ++i) {
            _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
          }
#else
          m_have_tuple = state[k].m_have_tuple;
          ht_pos = (uint64_t *)&state[k].ht_off;
          for (int i = 0; (i < VECTOR_SCALE) & (m_have_tuple);
              ++i, (m_have_tuple >> 1)) {
            if (m_have_tuple & 1) {
              _mm_prefetch((char *)(ht_pos[i]), _MM_HINT_T0);
            }
          }
#endif
        }
      }
        break;
    }
    ++k;
  }
  return matches;
}
int64_t smv_probe(hashtable_t *ht, relation_t *rel, void *output) {
  int64_t matches = 0;
  int32_t new_add = 0, k = 0, done = 0, num, num_temp;
  __attribute__((aligned(64)))  __mmask8 m_match = 0, m_new_cells = -1, m_valid_bucket = 0, mask[VECTOR_SCALE + 1],m_have_tuple;
  __m512i v_offset = _mm512_set1_epi64(0), v_base_offset_upper = _mm512_set1_epi64(rel->num_tuples * sizeof(tuple_t)), v_base_offset, v_ht_cell, v_factor = _mm512_set1_epi64(
      ht->hash_mask), v_shift = _mm512_set1_epi64(ht->skip_bits), v_cell_hash, v_neg_one512 = _mm512_set1_epi64(-1), v_zero512 = _mm512_set1_epi64(0), v_write_index =
      _mm512_set1_epi64(0), v_ht_addr = _mm512_set1_epi64((uint64_t)ht->buckets), v_word_size = _mm512_set1_epi64(WORDSIZE), v_tuple_size = _mm512_set1_epi64(sizeof(tuple_t)),
      v_bucket_size = _mm512_set1_epi64(sizeof(bucket_t)), v_next_off = _mm512_set1_epi64(8), v_right_payload, v_payload_off = _mm512_set1_epi64(24);
  chainedtuplebuffer_t *chainedbuf = (chainedtuplebuffer_t *) output;
  tuple_t *join_res = NULL;
  __attribute__((aligned(64)))  uint64_t cur_offset = 0, base_off[16], *ht_pos;
  for (int i = 0; i <= VECTOR_SCALE; ++i) {
    base_off[i] = i * sizeof(tuple_t);
    mask[i] = (1 << i) - 1;
  }
  v_base_offset = _mm512_load_epi64(base_off);
  __attribute__((aligned(64)))  StateSIMD state[SIMDStateSize + 1];
  // init # of the state
  for (int i = 0; i <= SIMDStateSize; ++i) {
    state[i].stage = 1;
    state[i].m_have_tuple = 0;
//    state[i].ht_off = _mm512_set1_epi64(0);
//    state[i].payload = _mm512_set1_epi64(0);
//    state[i].key = _mm512_set1_epi64(0);
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
//        _mm_prefetch((char *)(((void *)rel->tuples) + cur_offset + PDIS + 128), _MM_HINT_T0);
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
#elif DIR_PREFETCH
        ht_pos = (uint64_t *) &state[k].ht_off;
        for (int i = 0; i < VECTOR_SCALE; ++i) {
          _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
        }
#else
        m_have_tuple = state[k].m_have_tuple;
        ht_pos = (uint64_t *)&state[k].ht_off;
        for (int i = 0; (i < VECTOR_SCALE) & (m_have_tuple);
            ++i, (m_have_tuple >> 1)) {
          if (m_have_tuple & 1) {
            _mm_prefetch((char *)(ht_pos[i]), _MM_HINT_T0);
          }
        }
#endif
      }
        break;
      case 2: {
        /////////////////// random access
        // check valid bucket
        v_ht_cell = _mm512_mask_i64gather_epi64(v_neg_one512, state[k].m_have_tuple, state[k].ht_off, 0, 1);
        m_valid_bucket = _mm512_cmpneq_epi64_mask(v_ht_cell, v_zero512);
        state[k].m_have_tuple = _mm512_kand(m_valid_bucket, state[k].m_have_tuple);
        state[k].stage = 0;
        num = _mm_popcnt_u32(state[k].m_have_tuple);
#if 1
        if (num == VECTOR_SCALE) {
#if KNL
          _mm512_mask_prefetch_i64gather_pd(
              state[k].ht_off, state[k].m_have_tuple, 0, 1, _MM_HINT_T0);
#else
          ht_pos = (uint64_t *) &state[k].ht_off;
          for (int i = 0; i < VECTOR_SCALE; ++i) {
            _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
          }
#endif
        } else if(num ==0) {
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
#elif DIR_PREFETCH
              ht_pos = (uint64_t *) &state[k].ht_off;
              for (int i = 0; i < VECTOR_SCALE; ++i) {
                _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
              }
#else
              m_have_tuple = state[k].m_have_tuple;
              ht_pos = (uint64_t *)&state[k].ht_off;
              for (int i = 0; (i < VECTOR_SCALE) & (m_have_tuple);
                  ++i, (m_have_tuple >> 1)) {
                if (m_have_tuple & 1) {
                  _mm_prefetch((char *)(ht_pos[i]), _MM_HINT_T0);
                }
              }
#endif
            }
          }
        }

      }
        break;
      case 0: {
        v_ht_cell = _mm512_mask_i64gather_epi64(v_neg_one512, state[k].m_have_tuple, _mm512_add_epi64(state[k].ht_off, v_tuple_size), 0, 1);  // note the offset of the tuple in %bucket_t%

        ///// step 4: compare;
        m_match = _mm512_cmpeq_epi64_mask(state[k].key, v_ht_cell);
        m_match = _mm512_kand(m_match, state[k].m_have_tuple);

        new_add = _mm_popcnt_u32(m_match);
        matches += new_add;

        // gather payloads
        v_right_payload = _mm512_mask_i64gather_epi64(v_neg_one512, m_match, _mm512_add_epi64(state[k].ht_off, v_payload_off), 0, 1);

        // update next
        state[k].ht_off = _mm512_mask_i64gather_epi64(v_zero512, state[k].m_have_tuple, _mm512_add_epi64(state[k].ht_off, v_next_off), 0, 1);
        state[k].m_have_tuple = _mm512_kand(_mm512_cmpneq_epi64_mask(state[k].ht_off, v_zero512), state[k].m_have_tuple);
#if WRITE_RESULTS
        // to scatter join results
        join_res = cb_next_n_writepos(chainedbuf, new_add);
#if SEQPREFETCH
        _mm_prefetch((char *)(((void *)join_res) + PDIS), _MM_HINT_T0);
//        _mm_prefetch((char *)(((void *)join_res) + PDIS + 64), _MM_HINT_T0);
//        _mm_prefetch((char *)(((void *)join_res) + PDIS + 128), _MM_HINT_T0);
#endif
        v_write_index = _mm512_mask_expand_epi64(v_zero512, m_match, v_base_offset);
        _mm512_mask_i64scatter_epi64((void * )join_res, m_match, v_write_index, state[k].payload, 1);
        v_write_index = _mm512_add_epi64(v_write_index, v_word_size);
        _mm512_mask_i64scatter_epi64(((void * )join_res), m_match, v_write_index, v_right_payload, 1);
#endif
        num = _mm_popcnt_u32(state[k].m_have_tuple);
#if 1

        if (num == VECTOR_SCALE) {
#if KNL
          _mm512_mask_prefetch_i64gather_pd(
              state[k].ht_off, state[k].m_have_tuple, 0, 1, _MM_HINT_T0);
#else
          ht_pos = (uint64_t *) &state[k].ht_off;
          for (int i = 0; i < VECTOR_SCALE; ++i) {
            _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
          }
#endif
        } else if(num ==0) {
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
#elif DIR_PREFETCH
              ht_pos = (uint64_t *) &state[k].ht_off;
              for (int i = 0; i < VECTOR_SCALE; ++i) {
                _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
              }
#else
              m_have_tuple = state[k].m_have_tuple;
              ht_pos = (uint64_t *)&state[k].ht_off;
              for (int i = 0; (i < VECTOR_SCALE) & (m_have_tuple);
                  ++i, (m_have_tuple >> 1)) {
                if (m_have_tuple & 1) {
                  _mm_prefetch((char *)(ht_pos[i]), _MM_HINT_T0);
                }
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
  return matches;
}

int64_t probe_simd_amac_compact2(hashtable_t *ht, relation_t *rel, void *output) {
  int64_t matches = 0;
  int32_t new_add = 0, k = 0, done = 0, num, num_temp;
  __mmask8 m_match = 0, m_new_cells = -1, m_valid_bucket = 0, mask[VECTOR_SCALE + 1], flag, m_have_tuple;
  __m512i v_offset = _mm512_set1_epi64(0), v_base_offset_upper = _mm512_set1_epi64(rel->num_tuples * sizeof(tuple_t)), v_base_offset, v_left_size = _mm512_set1_epi64(8),
      v_bucket_offset = _mm512_set1_epi64(0), v_ht_cell, v_factor = _mm512_set1_epi64(ht->hash_mask), v_shift = _mm512_set1_epi64(ht->skip_bits), v_cell_hash, v_neg_one512 =
          _mm512_set1_epi64(-1), v_ht_upper, v_zero512 = _mm512_set1_epi64(0), v_write_index = _mm512_set1_epi64(0), v_ht_addr = _mm512_set1_epi64((uint64_t)ht->buckets), v_word_size =
          _mm512_set1_epi64(WORDSIZE), v_tuple_size = _mm512_set1_epi64(sizeof(tuple_t)), v_bucket_size = _mm512_set1_epi64(sizeof(bucket_t)), v_next_off = _mm512_set1_epi64(8),
      v_left_payload, v_right_payload, v_payload_off = _mm512_set1_epi64(24);
  chainedtuplebuffer_t *chainedbuf = (chainedtuplebuffer_t *) output;
  tuple_t *join_res = NULL;
  __attribute__((aligned(64)))  uint64_t cur_offset = 0, base_off[16], *left_payload, *right_payload, *ht_pos;
  left_payload = (uint64_t *) &v_left_payload;
  right_payload = (uint64_t *) &v_right_payload;

  for (int i = 0; i <= VECTOR_SCALE; ++i) {
    base_off[i] = i * sizeof(tuple_t);
    mask[i] = (1 << i) - 1;
  }
  v_base_offset = _mm512_load_epi64(base_off);
  StateSIMD state[SIMDStateSize + 1], state_temp[2];
  // init # of the state
  for (int i = 0; i <= SIMDStateSize; ++i) {
    state[i].stage = 1;
    state[i].m_have_tuple = 0;
    state[i].ht_off = _mm512_set1_epi64(0);
    state[i].tb_off = _mm512_set1_epi64(0);
    state[i].key = _mm512_set1_epi64(0);
  }
  for (uint64_t cur = 0; 1;) {
    k = (k >= SIMDStateSize) ? 0 : k;
    if (cur >= rel->num_tuples) {
      if (state[k].m_have_tuple == 0 && state[k].stage != 3) {
        ++done;
        state[k].stage = 3;
        ++k;
        continue;
      }
      if (done >= SIMDStateSize) {
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
        // directly use cur, instead of cur_offset to control the offset to
        // rel.
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
        state[k].stage = 0;
#if KNL
        _mm512_mask_prefetch_i64gather_pd(
            state[k].ht_off, state[k].m_have_tuple, 0, 1, _MM_HINT_T0);
#elif DIR_PREFETCH
        ht_pos = (uint64_t *) &state[k].ht_off;
        for (int i = 0; i < VECTOR_SCALE; ++i) {
          _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
        }
#else
        m_have_tuple = state[k].m_have_tuple;
        ht_pos = (uint64_t *)&state[k].ht_off;
        for (int i = 0; (i < VECTOR_SCALE) & (m_have_tuple);
            ++i, (m_have_tuple >> 1)) {
          if (m_have_tuple & 1) {
            _mm_prefetch((char *)(ht_pos[i]), _MM_HINT_T0);
          }
        }
#endif
      }
        break;
      case 0: {
        /////////////////// random access
        // check valid bucket
        v_ht_cell = _mm512_mask_i64gather_epi64(v_neg_one512, state[k].m_have_tuple, state[k].ht_off, 0, 1);
        m_valid_bucket = _mm512_cmpneq_epi64_mask(v_ht_cell, v_zero512);
        state[k].m_have_tuple = _mm512_kand(m_valid_bucket, state[k].m_have_tuple);

        v_ht_cell = _mm512_mask_i64gather_epi64(v_neg_one512, state[k].m_have_tuple, _mm512_add_epi64(state[k].ht_off, v_tuple_size), 0, 1);  // note the offset of the tuple in %bucket_t%

        ///// step 4: compare;
        m_match = _mm512_cmpeq_epi64_mask(state[k].key, v_ht_cell);
        m_match = _mm512_kand(m_match, state[k].m_have_tuple);
        new_add = _mm_popcnt_u32(m_match);
        matches += new_add;

        // gather payloads
        v_right_payload = _mm512_mask_i64gather_epi64(v_neg_one512, m_match, _mm512_add_epi64(state[k].ht_off, v_payload_off), 0, 1);

        state[k].ht_off = _mm512_mask_i64gather_epi64(v_zero512, state[k].m_have_tuple, _mm512_add_epi64(state[k].ht_off, v_next_off), 0, 1);
        state[k].m_have_tuple = _mm512_kand(_mm512_cmpneq_epi64_mask(state[k].ht_off, v_zero512), state[k].m_have_tuple);

        // to scatter join results
        join_res = cb_next_n_writepos(chainedbuf, new_add);
#if SEQPREFETCH
        _mm_prefetch((char *)(((void *)join_res) + PDIS), _MM_HINT_T0);
        _mm_prefetch((char *)(((void *)join_res) + PDIS + 64), _MM_HINT_T0);
        _mm_prefetch((char *)(((void *)join_res) + PDIS + 128), _MM_HINT_T0);
#endif
        v_write_index = _mm512_mask_expand_epi64(v_zero512, m_match, v_base_offset);
        _mm512_mask_i64scatter_epi64((void * )join_res, m_match, v_write_index, state[k].payload, 1);
        v_write_index = _mm512_add_epi64(v_write_index, v_word_size);
        _mm512_mask_i64scatter_epi64((void * )join_res, m_match, v_write_index, v_right_payload, 1);
#if 0
        if (num == VECTOR_SCALE) {
          state[k].stage = 0;
#if KNL
          _mm512_mask_prefetch_i64gather_pd(
              state[k].ht_off, state[k].m_have_tuple, 0, 1, _MM_HINT_T0);
#else
          ht_pos = (uint64_t *)&state[k].ht_off;
          for (int i = 0; i < VECTOR_SCALE; ++i) {
            _mm_prefetch((char *)(ht_pos[i]), _MM_HINT_T0);
          }
#endif
        } else
#else
        if (done < SIMDStateSize)
#endif
        {
          num = _mm_popcnt_u32(state[k].m_have_tuple);
          num_temp = _mm_popcnt_u32(state[SIMDStateSize].m_have_tuple);
          flag = (num + num_temp >= VECTOR_SCALE);
          state[k].stage = !flag;
          flag = flag - 1;  // denote <
          // case: num + num_temp < VECTOR_SCALE
          // compress v
          state_temp[0].ht_off = _mm512_maskz_compress_epi64(state[k].m_have_tuple, state[k].ht_off);
          state_temp[0].key = _mm512_maskz_compress_epi64(state[k].m_have_tuple, state[k].key);
          state_temp[0].payload = _mm512_maskz_compress_epi64(state[k].m_have_tuple, state[k].payload);
          // expand v -> temp
          state_temp[1].ht_off = _mm512_mask_expand_epi64(state[SIMDStateSize].ht_off, _mm512_knot(state[SIMDStateSize].m_have_tuple), state_temp[0].ht_off);
          state_temp[1].key = _mm512_mask_expand_epi64(state[SIMDStateSize].key, _mm512_knot(state[SIMDStateSize].m_have_tuple), state_temp[0].key);
          state_temp[1].payload = _mm512_mask_expand_epi64(state[SIMDStateSize].payload, _mm512_knot(state[SIMDStateSize].m_have_tuple), state_temp[0].payload);

          {
            // case: num + num_temp >= VECTOR_SCALE
            // expand temp -> v
            state[k].ht_off = _mm512_mask_expand_epi64(state[k].ht_off, _mm512_knot(state[k].m_have_tuple), state[SIMDStateSize].ht_off);
            state[k].key = _mm512_mask_expand_epi64(state[k].key, _mm512_knot(state[k].m_have_tuple), state[SIMDStateSize].key);
            state[k].payload = _mm512_mask_expand_epi64(state[k].payload, _mm512_knot(state[k].m_have_tuple), state[SIMDStateSize].payload);
            // compress temp
            state[SIMDStateSize].m_have_tuple = _mm512_kand(state[SIMDStateSize].m_have_tuple, _mm512_knot(mask[VECTOR_SCALE - num]));
            state[SIMDStateSize].ht_off = _mm512_maskz_compress_epi64(state[SIMDStateSize].m_have_tuple, state[SIMDStateSize].ht_off);
            state[SIMDStateSize].key = _mm512_maskz_compress_epi64(state[SIMDStateSize].m_have_tuple, state[SIMDStateSize].key);
            state[SIMDStateSize].payload = _mm512_maskz_compress_epi64(state[SIMDStateSize].m_have_tuple, state[SIMDStateSize].payload);
          }

          {  // merge the results of the two cases
            state[k].ht_off = _mm512_mask_mov_epi64(state[k].ht_off, flag, state_temp[0].ht_off);
            state[k].key = _mm512_mask_mov_epi64(state[k].key, flag, state_temp[0].key);
            state[k].payload = _mm512_mask_mov_epi64(state[k].payload, flag, state_temp[0].payload);

            state[SIMDStateSize].ht_off = _mm512_mask_mov_epi64(state[SIMDStateSize].ht_off, flag, state_temp[1].ht_off);
            state[SIMDStateSize].key = _mm512_mask_mov_epi64(state[SIMDStateSize].key, flag, state_temp[1].key);
            state[SIMDStateSize].payload = _mm512_mask_mov_epi64(state[SIMDStateSize].payload, flag, state_temp[1].payload);
            state[k].m_have_tuple = _mm512_kor(0, _mm512_kand(mask[VECTOR_SCALE], _mm512_knot(flag)));

            state[SIMDStateSize].m_have_tuple = _mm512_kor(_mm512_kand(mask[num + num_temp], flag),
                                                           _mm512_kand((state[SIMDStateSize].m_have_tuple >> (VECTOR_SCALE - num)), _mm512_knot(flag)));

#if KNL
            _mm512_mask_prefetch_i64gather_pd(
                state[k].ht_off, state[k].m_have_tuple, 0, 1, _MM_HINT_T0);
#elif DIR_PREFETCH
            ht_pos = (uint64_t *) &state[k].ht_off;
            for (int i = 0; i < VECTOR_SCALE; ++i) {
              _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
            }
#else
            m_have_tuple = state[k].m_have_tuple;
            ht_pos = (uint64_t *)&state[k].ht_off;
            for (int i = 0; (i < VECTOR_SCALE) & (m_have_tuple);
                ++i, (m_have_tuple >> 1)) {
              if (m_have_tuple & 1) {
                _mm_prefetch((char *)(ht_pos[i]), _MM_HINT_T0);
              }
            }
#endif
          }
        }
      }
        break;
    }
    ++k;
  }
  return matches;
}

int64_t probe_simd_gp(hashtable_t *ht, relation_t *rel, void *output) {
  int64_t matches = 0;
  int32_t new_add = 0, k = 0, done = 0, stage2_size = 0;
  __mmask8 m_match = 0, m_new_cells = -1, m_valid_bucket = 0;
  __m512i v_offset = _mm512_set1_epi64(0), v_base_offset_upper = _mm512_set1_epi64(rel->num_tuples * sizeof(tuple_t)), v_base_offset, v_left_size = _mm512_set1_epi64(8),
      v_bucket_offset = _mm512_set1_epi64(0), v_ht_cell, v_factor = _mm512_set1_epi64(ht->hash_mask), v_shift = _mm512_set1_epi64(ht->skip_bits), v_cell_hash, v_neg_one512 =
          _mm512_set1_epi64(-1), v_ht_upper, v_zero512 = _mm512_set1_epi64(0), v_write_index = _mm512_set1_epi64(0), v_ht_addr = _mm512_set1_epi64((uint64_t)ht->buckets), v_word_size =
          _mm512_set1_epi64(WORDSIZE), v_tuple_size = _mm512_set1_epi64(sizeof(tuple_t)), v_bucket_size = _mm512_set1_epi64(sizeof(bucket_t)), v_next_off = _mm512_set1_epi64(8),
      v_left_payload, v_right_payload, v_payload_off = _mm512_set1_epi64(24);
  chainedtuplebuffer_t *chainedbuf = (chainedtuplebuffer_t *) output;
  tuple_t *join_res = NULL;
  __attribute__((aligned(64)))  uint64_t cur_offset = 0, base_off[16], *left_payload, *right_payload, *ht_pos;

  left_payload = (uint64_t *) &v_left_payload;
  right_payload = (uint64_t *) &v_right_payload;

  for (int i = 0; i <= VECTOR_SCALE; ++i) {
    base_off[i] = i * sizeof(tuple_t);
  }
  v_base_offset = _mm512_load_epi64(base_off);
  StateSIMD state[SIMDStateSize];
  // init # of the state
  for (int i = 0; i < SIMDStateSize; ++i) {
    state[i].stage = 1;
    state[i].m_have_tuple = 0;
    state[i].ht_off = _mm512_set1_epi64(0);
    state[i].tb_off = _mm512_set1_epi64(0);
    state[i].key = _mm512_set1_epi64(0);
  }
  for (uint64_t cur = 0; (cur < rel->num_tuples);) {
    for (k = 0; (k < SIMDStateSize) && (cur < rel->num_tuples); ++k) {
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

      state[k].tb_off = _mm512_add_epi64(_mm512_set1_epi64(cur_offset), v_base_offset);
      new_add = VECTOR_SCALE;
      cur_offset = cur_offset + base_off[new_add];
      cur = cur + new_add;
      state[k].m_have_tuple = _mm512_cmpgt_epi64_mask(v_base_offset_upper, state[k].tb_off);
      ///// step 2: load new cells from right tuples;
      m_new_cells = state[k].m_have_tuple;
      // maybe need offset within a tuple
      state[k].key = _mm512_mask_i64gather_epi64(state[k].key, m_new_cells, state[k].tb_off, ((void * )rel->tuples), 1);
      state[k].payload = _mm512_mask_i64gather_epi64(state[k].payload, state[k].m_have_tuple, _mm512_add_epi64(v_offset, v_word_size), ((void * )rel->tuples), 1);
      ///// step 3: load new values from hash tables;
      // hash the cell values
      v_cell_hash = _mm512_and_epi64(state[k].key, v_factor);
      v_cell_hash = _mm512_srlv_epi64(v_cell_hash, v_shift);
      v_cell_hash = _mm512_mullo_epi64(v_cell_hash, v_bucket_size);
      state[k].ht_off = _mm512_mask_add_epi64(v_ht_addr, m_new_cells, v_cell_hash, v_ht_addr);
      state[k].stage = 0;
#if KNL
      _mm512_mask_prefetch_i64gather_pd(state[k].ht_off, state[k].m_have_tuple,
          0, 1, _MM_HINT_T0);
#elif DIR_PREFETCH
      ht_pos = (uint64_t *) &state[k].ht_off;
      for (int i = 0; i < VECTOR_SCALE; ++i) {
        _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
      }
#else
      m_have_tuple = state[k].m_have_tuple;
      ht_pos = (uint64_t *)&state[k].ht_off;
      for (int i = 0; (i < VECTOR_SCALE) & (m_have_tuple);
          ++i, (m_have_tuple >> 1)) {
        if (m_have_tuple & 1) {
          _mm_prefetch((char *)(ht_pos[i]), _MM_HINT_T0);
        }
      }
#endif
    }
    stage2_size = k;
    done = 0;
    while (done < stage2_size) {
      for (k = 0; k < stage2_size; ++k) {
        if (state[k].m_have_tuple == 0) {
          if (state[k].stage == 0) {
            ++done;
          }
          state[k].stage = 1;
          continue;
        }
        /////////////////// random access
        // check valid bucket
        v_ht_cell = _mm512_mask_i64gather_epi64(v_neg_one512, state[k].m_have_tuple, state[k].ht_off, 0, 1);
        m_valid_bucket = _mm512_cmpneq_epi64_mask(v_ht_cell, v_zero512);
        state[k].m_have_tuple = _mm512_kand(m_valid_bucket, state[k].m_have_tuple);

        v_ht_cell = _mm512_mask_i64gather_epi64(v_neg_one512, state[k].m_have_tuple, _mm512_add_epi64(state[k].ht_off, v_tuple_size), 0, 1);  // note the offset of the tuple in %bucket_t%

        ///// step 4: compare;
        m_match = _mm512_cmpeq_epi64_mask(state[k].key, v_ht_cell);
        m_match = _mm512_kand(m_match, state[k].m_have_tuple);
        new_add = _mm_popcnt_u32(m_match);
        matches += new_add;

        // gather payloads
        v_right_payload = _mm512_mask_i64gather_epi64(v_neg_one512, m_match, _mm512_add_epi64(state[k].ht_off, v_payload_off), 0, 1);

        state[k].ht_off = _mm512_mask_i64gather_epi64(v_zero512, state[k].m_have_tuple, _mm512_add_epi64(state[k].ht_off, v_next_off), 0, 1);
        state[k].m_have_tuple = _mm512_kand(_mm512_cmpneq_epi64_mask(state[k].ht_off, v_zero512), state[k].m_have_tuple);

        // to scatter join results
        join_res = cb_next_n_writepos(chainedbuf, new_add);
#if SEQPREFETCH
        _mm_prefetch((char *)(((void *)join_res) + PDIS), _MM_HINT_T0);
        _mm_prefetch((char *)(((void *)join_res) + PDIS + 64), _MM_HINT_T0);
        _mm_prefetch((char *)(((void *)join_res) + PDIS + 128), _MM_HINT_T0);
#endif
        v_write_index = _mm512_mask_expand_epi64(v_zero512, m_match, v_base_offset);
        _mm512_mask_i64scatter_epi64((void * )join_res, m_match, v_write_index, state[k].payload, 1);
        v_write_index = _mm512_add_epi64(v_write_index, v_word_size);
        _mm512_mask_i64scatter_epi64((void * )join_res, m_match, v_write_index, v_right_payload, 1);
#if KNL
        _mm512_mask_prefetch_i64gather_pd(
            state[k].ht_off, state[k].m_have_tuple, 0, 1, _MM_HINT_T0);
#elif DIR_PREFETCH
        ht_pos = (uint64_t *) &state[k].ht_off;
        for (int i = 0; i < VECTOR_SCALE; ++i) {
          _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
        }
#else
        m_have_tuple = state[k].m_have_tuple;
        ht_pos = (uint64_t *)&state[k].ht_off;
        for (int i = 0; (i < VECTOR_SCALE) & (m_have_tuple);
            ++i, (m_have_tuple >> 1)) {
          if (m_have_tuple & 1) {
            _mm_prefetch((char *)(ht_pos[i]), _MM_HINT_T0);
          }
        }
#endif
      }
    }
  }
  return matches;
}
#endif

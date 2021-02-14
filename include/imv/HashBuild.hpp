#pragma once

#include "head.hpp"
struct __attribute__((aligned(64))) BuildSIMDState {
  __m512i v_entry_addr;
  __m512i v_hash_value;
  __m512i v_build_key;
  __mmask8 m_valid;
  uint8_t valid_size, stage,mask[VECTORSIZE+1];
  BuildSIMDState()
      : v_entry_addr(_mm512_set1_epi64(0)),
        v_hash_value(_mm512_set1_epi64(0)),
        valid_size(VECTORSIZE),
        stage(1),
        m_valid(0){
    for (int i = 0; i <= VECTORSIZE; ++i) {
      mask[i] = (1 << i) - 1;
    }
  }
  inline void reset() {
    v_entry_addr = _mm512_set1_epi64(0);
    v_hash_value = _mm512_set1_epi64(0);
    valid_size = VECTORSIZE;
    stage = 1;
    m_valid=0;
  }
  void* operator new(size_t size) {
    return memalign(64, size);
  }
  void operator delete(void* mem) {
    return free(mem);
  }
  void* operator new[](size_t size) {
    return memalign(64, size);
  }
  void operator delete[](void* mem) {
    return free(mem);
  }
  inline void compress() {
    v_entry_addr = _mm512_maskz_compress_epi64(m_valid, v_entry_addr);
    v_hash_value = _mm512_maskz_compress_epi64(m_valid, v_hash_value);
    v_build_key = _mm512_maskz_compress_epi64(m_valid, v_build_key);
  }
  inline void expand(BuildSIMDState& src_state) {
    v_entry_addr = _mm512_mask_expand_epi64(v_entry_addr, _mm512_knot(m_valid), src_state.v_entry_addr);
    v_hash_value = _mm512_mask_expand_epi64(v_hash_value, _mm512_knot(m_valid), src_state.v_hash_value);
    v_build_key = _mm512_mask_expand_epi64(v_build_key, _mm512_knot(m_valid), src_state.v_build_key);


  }
  inline void compact(BuildSIMDState& RVS, uint8_t done, uint8_t imvNum,uint8_t& k,uint8_t next_stage,uint8_t src_stage){
   auto num = _mm_popcnt_u32(m_valid);
    if (num == VECTORSIZE || done >= imvNum) {
      stage = next_stage;
    } else {
       auto num_temp = _mm_popcnt_u32(RVS.m_valid);
        if (num + num_temp < VECTORSIZE) {
          // compress imv_state[k]
          compress();
          // expand imv_state[k] -> imv_state[imvNum1]
          RVS.expand(*this);
          RVS.m_valid = mask[num + num_temp];
          m_valid = 0;
          stage = src_stage;
          RVS.stage = next_stage;
          --k;
        } else {
          // expand imv_state[imvNum1] -> expand imv_state[k]
          expand(RVS);
          RVS.m_valid = _mm512_kand(RVS.m_valid, _mm512_knot(mask[VECTORSIZE - num]));
          // compress imv_state[imvNum]
          RVS.compress();
          RVS.m_valid = RVS.m_valid >> (VECTORSIZE - num);
          m_valid = mask[VECTORSIZE];
          RVS.stage = next_stage;
          stage = next_stage;
        }
    }
  }
};
struct BuildState{
  Hashmap::EntryHeader* ptr;
  hash_t hash_value;
  uint8_t stage;
};
size_t build_raw(size_t begin, size_t end, Database& db, runtime::Hashmap* hash_table, Allocator*allo, int entry_size);
size_t build_simd(size_t begin, size_t end, Database& db, runtime::Hashmap* hash_table, Allocator*allo, int entry_size);
size_t build_imv(size_t begin, size_t end, Database& db, runtime::Hashmap* hash_table, Allocator*allo, int entry_size);
size_t build_gp(size_t begin, size_t end, Database& db, runtime::Hashmap* hash_table, Allocator*allo, int entry_size);
size_t build_amac(size_t begin, size_t end, Database& db, runtime::Hashmap* hash_table, Allocator*allo, int entry_size);

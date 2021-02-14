#pragma once
#include "head.hpp"
#define ORDERKEY 1
#if ORDERKEY
#define HT_SIZE 150000
#else
#define HT_SIZE 100000
#endif
#define PARTITION_SIZE 4096
struct AMACState {
  uint8_t stage;
  int probeKey = 0;
  types::Numeric<12, 2> probeValue;
  pos_t tuple_id = 0;
  runtime::Hashmap::hash_t probeHash;
  runtime::Hashmap::EntryHeader* buildMatch;
};
struct __attribute__((aligned(64))) AggState {
  __m512i v_probe_keys;
  __m512i v_bucket_addrs;
  __m512i v_probe_offset, v_probe_hash, v_probe_value;
  __mmask8 m_valid_probe;
  uint8_t stage;
  AggState()
      : v_probe_keys(_mm512_set1_epi64(0)),
        v_bucket_addrs(_mm512_set1_epi64(0)),
        v_probe_offset(_mm512_set1_epi64(0)),
        m_valid_probe(0),
        stage(1) {
  }
  void reset() {
    v_probe_keys = _mm512_set1_epi64(0);
    v_bucket_addrs = _mm512_set1_epi64(0);
    v_probe_offset = _mm512_set1_epi64(0);
    m_valid_probe = 0;
    stage = 1;
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
    v_bucket_addrs = _mm512_maskz_compress_epi64(m_valid_probe, v_bucket_addrs);
    v_probe_keys = _mm512_maskz_compress_epi64(m_valid_probe, v_probe_keys);
    v_probe_value = _mm512_maskz_compress_epi64(m_valid_probe, v_probe_value);
    v_probe_offset = _mm512_maskz_compress_epi64(m_valid_probe, v_probe_offset);
    v_probe_hash = _mm512_maskz_compress_epi64(m_valid_probe, v_probe_hash);

  }
  inline void expand(AggState& src_state) {
    v_bucket_addrs = _mm512_mask_expand_epi64(v_bucket_addrs, _mm512_knot(m_valid_probe), src_state.v_bucket_addrs);
    v_probe_keys = _mm512_mask_expand_epi64(v_probe_keys, _mm512_knot(m_valid_probe), src_state.v_probe_keys);
    v_probe_value = _mm512_mask_expand_epi64(v_probe_value, _mm512_knot(m_valid_probe), src_state.v_probe_value);
    v_probe_offset = _mm512_mask_expand_epi64(v_probe_offset, _mm512_knot(m_valid_probe), src_state.v_probe_offset);
    v_probe_hash = _mm512_mask_expand_epi64(v_probe_hash, _mm512_knot(m_valid_probe), src_state.v_probe_hash);

  }
};
inline void write_results(int* o_key, uint64_t* o_value, void** results, size_t found) {
  __m512i v_values, address;
  __m256i v_keys;
  __mmask8 m_valid = -1;
  uint64_t key_off = 16, value_off = 24;
  for (int i = 0; i < found; i += VECTORSIZE) {
    if (i + VECTORSIZE >= found) {
      m_valid = (m_valid >> (i + VECTORSIZE - found));
    }
    address = _mm512_maskz_loadu_epi64(m_valid, (results + i));
    v_keys = _mm512_mask_i64gather_epi32(v_keys, m_valid, address+Vec8u(key_off), nullptr, 1);
    v_values = _mm512_mask_i64gather_epi64(v_values, m_valid, address+Vec8u(value_off), nullptr, 1);
    _mm256_mask_storeu_epi32((o_key + i), m_valid, v_keys);
    _mm512_mask_storeu_epi64((o_value + i), m_valid, v_values);
  }
}

size_t agg_raw(size_t begin, size_t end, Database& db, Hashmapx<types::Integer, types::Numeric<12, 2>, hashFun, false>* hash_table,
               PartitionedDeque<PARTITION_SIZE>* partition, void** entry_addrs = nullptr, void** results_entry = nullptr);
size_t agg_amac(size_t begin, size_t end, Database& db, Hashmapx<types::Integer, types::Numeric<12, 2>, hashFun, false>* hash_table,
                PartitionedDeque<PARTITION_SIZE>* partition, void** entry_addrs = nullptr, void** results_entry = nullptr);
size_t agg_gp(size_t begin, size_t end, Database& db, Hashmapx<types::Integer, types::Numeric<12, 2>, hashFun, false>* hash_table,
              PartitionedDeque<PARTITION_SIZE>* partition, void** entry_addrs = nullptr, void** results_entry = nullptr);
size_t agg_simd(size_t begin, size_t end, Database& db, Hashmapx<types::Integer, types::Numeric<12, 2>, hashFun, false>* hash_table,
                PartitionedDeque<PARTITION_SIZE>* partition, void** entry_addrs = nullptr, void** results_entry = nullptr);
size_t agg_imv_hybrid(size_t begin, size_t end, Database& db, Hashmapx<types::Integer, types::Numeric<12, 2>, hashFun, false>* hash_table,
               PartitionedDeque<PARTITION_SIZE>* partition, void** entry_addrs = nullptr, void** results_entry = nullptr);
size_t agg_imv1(size_t begin, size_t end, Database& db, Hashmapx<types::Integer, types::Numeric<12, 2>, hashFun, false>* hash_table,
               PartitionedDeque<PARTITION_SIZE>* partition, void** entry_addrs = nullptr, void** results_entry = nullptr);
size_t agg_imv_serial(size_t begin, size_t end, Database& db, Hashmapx<types::Integer, types::Numeric<12, 2>, hashFun, false>* hash_table,
               PartitionedDeque<PARTITION_SIZE>* partition, void** entry_addrs = nullptr, void** results_entry = nullptr);
size_t agg_imv_merged(size_t begin, size_t end, Database& db, Hashmapx<types::Integer, types::Numeric<12, 2>, hashFun, false>* hash_table,
               PartitionedDeque<PARTITION_SIZE>* partition, void** entry_addrs = nullptr, void** results_entry = nullptr);

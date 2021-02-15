#ifndef SRC_PREFETCH_H_
#define SRC_PREFETCH_H_
#include <assert.h>
#include <immintrin.h>
#include<algorithm>
#include "npj_types.hpp"
#include "types.hpp"

using std::sort;
typedef struct amac_state_t scalar_state_t;
typedef struct StateSIMD StateSIMD;
#define UNLIKELY(expr) expr //__builtin_expect(!!(expr), 0)
#define LIKELY(expr) expr // __builtin_expect(!!(expr), 1)

#define ScalarStateSize 20
#define PDIS 320
#define SIMDStateSize 5
#define WRITE_RESULTS 1
#define LOAD_FACTOR 1
#define MULTI_TUPLE (BUCKET_SIZE - 1)
#define REPEAT_PROBE 10
#define SLEEP_TIME 0
#define VECTOR_SCALE 8
#define DIR_PREFETCH 1
#define SEQPREFETCH PDIS
#define DIVIDE 0
#define USE_TBB 0
#define AFFINITY 0
//#define MORSE_SIZE 100000
#define SORTED 0
#define TEST_NUMA 0
#if KNL
#define _mm512_mullo_epi64(a, b) _mm512_mullo_epi32(a, b)
#endif

//#define _mm512_mask_i64scatter_epi64(addr, mask, idx, v, scale) \
  _mm512_mask_compressstoreu_epi64(addr, mask, v);
struct amac_state_t {
  int64_t tuple_id,key,payload;
  bucket_t *b;
  int16_t stage;
};
struct StateSIMD {
  __m512i key;
  __m512i payload;
  __m512i tb_off;
  __m512i ht_off;
  __mmask8 m_have_tuple;
  char stage;
};
//struct Vec8u{
//  union {
//     __m512i reg;
//     uint64_t entry[8];
//  };
//  // constructor
//  Vec8u(uint64_t x) { reg = _mm512_set1_epi64(x); };
//  Vec8u(void* p) { reg = _mm512_loadu_si512(p); };
//  Vec8u(__m512i x) { reg = x; };
//  Vec8u(uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7) { reg = _mm512_set_epi64(x0, x1, x2, x3, x4, x5, x6, x7); };
//
//  // implicit conversion to register
//  operator __m512i() { return reg; }
//
//  // print vector (for debugging)
//  friend std::ostream& operator<< (std::ostream& stream, const Vec8u& v) {
//     for (auto& e : v.entry)
//        stream << e << " ";
//     return stream;
//  }
//};
static bool cmp_tuples(tuple_t& a,tuple_t& b){
  return a.key<b.key;
}
inline void sort_rel(relation_t * rel){
  sort(rel->tuples,rel->tuples+rel->num_tuples,cmp_tuples);
}
inline void v_prefetch(__m512i& vec){
  uint64_t * ht_pos = (uint64_t*)&vec;
  for (int i = 0; i < VECTOR_SCALE; ++i) {
    _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
//    _mm_prefetch(((char * )(ht_pos[i]) + 64), _MM_HINT_T0);
  }
}
#endif /* SRC_PREFETCH_H_ */

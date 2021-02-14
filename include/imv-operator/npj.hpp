#ifndef _SIMD_PREFETCHING_HEAD
#define _SIMD_PREFETCHING_HEAD
#include "prefetch.hpp"

int64_t probe_simd(hashtable_t *ht, relation_t *rel, void *output);
int64_t probe_simd_amac(hashtable_t *ht, relation_t *rel, void *output);
int64_t smv_probe(hashtable_t *ht, relation_t *rel, void *output);

int64_t probe_simd(hashtable_t *ht, relation_t *rel, void *output);

int64_t probe_simd_amac_raw(hashtable_t *ht, relation_t *rel, void *output);


#endif

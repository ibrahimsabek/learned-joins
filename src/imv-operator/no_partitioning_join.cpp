/**
 * @file    no_partitioning_join.c
 * @author  Cagri Balkesen <cagri.balkesen@inf.ethz.ch>
 * @date    Sun Feb  5 20:16:58 2012
 * @version $Id: no_partitioning_join.c 4419 2013-10-21 16:24:35Z bcagri $
 *
 * @brief  The implementation of NPO, No Partitioning Optimized join algortihm.
 *
 * (c) 2012, ETH Zurich, Systems Group
 *
 */

#include "imv-operator/no_partitioning_join.h"

#include <string.h>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "imv-operator/npj.hpp"
//#include "tbb/tbb.h"
//#include "profile.h"

//using std::__cxx11::string;
using std::string;
using std::make_pair;
using std::pair;
using std::vector;

/**
 * @defgroup OverflowBuckets Buffer management for overflowing buckets.
 * Simple buffer management for overflow-buckets organized as a
 * linked-list of bucket_buffer_t.
 * @{
 */

/**
 * Initializes a new bucket_buffer_t for later use in allocating
 * buckets when overflow occurs.
 *
 * @param ppbuf [in,out] bucket buffer to be initialized
 */
void init_bucket_buffer(bucket_buffer_t **ppbuf) {
  bucket_buffer_t *overflowbuf;
  overflowbuf = (bucket_buffer_t *) malloc(sizeof(bucket_buffer_t));
  if (posix_memalign((void **) &(overflowbuf->buf), PAGE_SIZE, sizeof(bucket_t) * OVERFLOW_BUF_SIZE)) {
    perror("overflow buffer : Aligned allocation failed!\n");
    exit(EXIT_FAILURE);
  }
  overflowbuf->count = 0;
  overflowbuf->next = NULL;

  *ppbuf = overflowbuf;
}

/**
 * Returns a new bucket_t from the given bucket_buffer_t.
 * If the bucket_buffer_t does not have enough space, then allocates
 * a new bucket_buffer_t and adds to the list.
 *
 * @param result [out] the new bucket
 * @param buf [in,out] the pointer to the bucket_buffer_t pointer
 */
void get_new_bucket(bucket_t **result, bucket_buffer_t **buf) {
  if ((*buf)->count < OVERFLOW_BUF_SIZE) {
    *result = (*buf)->buf + (*buf)->count;
    (*buf)->count++;
  } else {
    /* need to allocate new buffer */
    bucket_buffer_t *new_buf = (bucket_buffer_t *) malloc(sizeof(bucket_buffer_t));
    if (posix_memalign((void **) &(new_buf->buf), PAGE_SIZE, sizeof(bucket_t) * OVERFLOW_BUF_SIZE)) {
      perror("overflow buffer : Aligned allocation failed!\n");
      exit(EXIT_FAILURE);
    }

    new_buf->count = 1;
    new_buf->next = *buf;
    *buf = new_buf;
    *result = new_buf->buf;
  }
}

/** De-allocates all the bucket_buffer_t */
void free_bucket_buffer(bucket_buffer_t *buf) {
  do {
    bucket_buffer_t *tmp = buf->next;
    free(buf->buf);
    free(buf);
    buf = tmp;
  } while (buf);
}

/** @} */

/**
 * @defgroup NPO The No Partitioning Optimized Join Implementation
 * @{
 */

/**
 * Allocates a hashtable of NUM_BUCKETS and inits everything to 0.
 *
 * @param ht pointer to a hashtable_t pointer
 */
void allocate_hashtable(hashtable_t **ppht, uint32_t nbuckets) {
  hashtable_t *ht;

  ht = (hashtable_t *) malloc(sizeof(hashtable_t));
  ht->num_buckets = nbuckets;
  NEXT_POW_2((ht->num_buckets));
  ht->num_buckets = ht->num_buckets / LOAD_FACTOR;

  /* allocate hashtable buckets cache line aligned */
  if (posix_memalign((void **) &ht->buckets, CACHE_LINE_SIZE, ht->num_buckets * sizeof(bucket_t))) {
    perror("Aligned allocation failed!\n");
    exit(EXIT_FAILURE);
  }

  /** Not an elegant way of passing whether we will numa-localize, but this
   feature is experimental anyway. */
  if (numalocalize) {
    tuple_t *mem = (tuple_t *) ht->buckets;
    uint32_t ntuples = (ht->num_buckets * sizeof(bucket_t)) / sizeof(tuple_t);
    numa_localize(mem, ntuples, nthreads);
  }

  memset(ht->buckets, 0, ht->num_buckets * sizeof(bucket_t));
  ht->skip_bits = 0; /* the default for modulo hash */
  ht->hash_mask = (ht->num_buckets - 1) << ht->skip_bits;
  *ppht = ht;
}

/**
 * Releases memory allocated for the hashtable.
 *
 * @param ht pointer to hashtable
 */
void destroy_hashtable(hashtable_t *ht) {
  free(ht->buckets);
  free(ht);
}

/**
 * print the distribution of the hash table
 *
 * @param ht pointer to hashtable
 */
#define LENGTH 100
void print_hashtable(hashtable_t * const ht) {
  uint32_t num[LENGTH], max = 0, len;
  uint64_t total = 0;
  bucket_t *curr = NULL;
  memset(num, 0, sizeof(num));
  for (uint32_t i = 0; i < ht->num_buckets; ++i) {
    curr = ht->buckets + i;
    len = (curr->lenth + 9) / 10;
    ++num[(len > LENGTH ? LENGTH - 1 : len)];
    if (curr->lenth > max) {
      max = curr->lenth;
    }
    total += curr->lenth;
  }
  // assert(max < LENGTH);
  printf("max = %d\t total = %lld\n", max, total);
  puts("======the statics of a hash table ====");
  for (uint32_t i = 0; i < LENGTH; ++i) {
    if (num[i] > 0) {
      printf("len= %5d, num= %10d\n", i * 10, num[i]);
    }
  }
  puts("==END=the statics of a hash table ====");
}
/**
 * Single-thread hashtable build method, ht is pre-allocated.
 *
 * @param ht hastable to be built
 * @param rel the build relation
 */
void build_hashtable_st(hashtable_t *ht, relation_t *rel) {
  uint32_t i;
  const uint32_t hashmask = ht->hash_mask;
  const uint32_t skipbits = ht->skip_bits;

  for (i = 0; i < rel->num_tuples; i++) {
    tuple_t *dest;
    bucket_t *curr, *nxt;
    int32_t idx = HASH(rel->tuples[i].key, hashmask, skipbits);

    /* copy the tuple to appropriate hash bucket */
    /* if full, follow nxt pointer to find correct place */
    curr = ht->buckets + idx;
    nxt = curr->next;

    if (curr->count == BUCKET_SIZE) {
      if (!nxt || nxt->count == BUCKET_SIZE) {
        bucket_t *b;
        b = (bucket_t *) calloc(1, sizeof(bucket_t));
        curr->next = b;
        b->next = nxt;
        b->count = 1;
        dest = b->tuples;
      } else {
        dest = nxt->tuples + nxt->count;
        nxt->count++;
      }
    } else {
      dest = curr->tuples + curr->count;
      curr->count++;
    }
    *dest = rel->tuples[i];
  }
}

/**
 * Probes the hashtable for the given outer relation, returns num results.
 * This probing method is used for both single and multi-threaded version.
 *
 * @param ht hashtable to be probed
 * @param rel the probing outer relation
 * @param output chained tuple buffer to write join results, i.e. rid pairs.
 *
 * @return number of matching tuples
 */
int64_t probe_hashtable(hashtable_t *ht, relation_t *rel, void *output) {
  uint32_t i, j;
  int64_t matches;

  const uint64_t hashmask = ht->hash_mask;
  const uint64_t skipbits = ht->skip_bits;
  matches = 0;
  j = 0;
#if WRITE_RESULTS
  chainedtuplebuffer_t *chainedbuf = (chainedtuplebuffer_t *) output;
#endif

  for (i = 0; i < rel->num_tuples; i++) {
    intkey_t idx = HASH(rel->tuples[i].key, hashmask, skipbits);
    bucket_t *b = ht->buckets + idx;
    do {
      if (b->count == 0) {
        break;
      }
      if (rel->tuples[i].key == b->tuples[0].key) {
        matches++;

#if WRITE_RESULTS
        /* copy to the result buffer */
        tuple_t *joinres = cb_next_writepos(chainedbuf);
        joinres->key = b->tuples[0].payload; /* R-rid */
        joinres->payload = rel->tuples[i].payload; /* S-rid */
#endif
      }
      b = b->next; /* follow overflow pointer */
    } while (b);
  }

  return matches;
}

/**
 * the raw prefetching in original ETH codes
 * Probes the hashtable for the given outer relation, returns num results.
 * This probing method is used for both single and multi-threaded version.
 *
 * @param ht hashtable to be probed
 * @param rel the probing outer relation
 * @param output chained tuple buffer to write join results, i.e. rid pairs.
 *
 * @return number of matching tuples
 */
#define PREFETCH_NPJ
int64_t probe_hashtable_raw_prefetch(hashtable_t *ht, relation_t *rel, void *output) {
  uint32_t i, j;
  int64_t matches;

  const uint32_t hashmask = ht->hash_mask;
  const uint32_t skipbits = ht->skip_bits;
#ifdef PREFETCH_NPJ
  size_t prefetch_index = PREFETCH_DISTANCE;
#endif

  matches = 0;

#ifdef JOIN_RESULT_MATERIALIZE
  chainedtuplebuffer_t *chainedbuf = (chainedtuplebuffer_t *) output;
#endif

  for (i = 0; i < rel->num_tuples; i++) {
#ifdef PREFETCH_NPJ
    _mm_prefetch((char *)(rel->tuples + prefetch_index) + PDIS, _MM_HINT_T0);
    if (prefetch_index < rel->num_tuples) {
      intkey_t idx_prefetch = HASH(rel->tuples[prefetch_index++].key, hashmask, skipbits);
      //__builtin_prefetch(ht->buckets + idx_prefetch, 0, 1);
      _mm_prefetch((char * )(ht->buckets + idx_prefetch), _MM_HINT_T0);
    }
#endif

    intkey_t idx = HASH(rel->tuples[i].key, hashmask, skipbits);
    bucket_t *b = ht->buckets + idx;
    do {
#if MULTI_TUPLE
      for (j = 0; j < b->count; j++) {
#else
      if (b->count == 0) {
        break;
      }
      j = 0;
#endif
      if (rel->tuples[i].key == b->tuples[j].key) {
        matches++;

#ifdef JOIN_RESULT_MATERIALIZE
        /* copy to the result buffer */
        tuple_t *joinres = cb_next_writepos(chainedbuf);
        joinres->key = b->tuples[j].payload; /* R-rid */
        joinres->payload = rel->tuples[i].payload; /* S-rid */
#endif
      }
#if MULTI_TUPLE
    }
#endif

      b = b->next; /* follow overflow pointer */
    } while (b);
  }

  return matches;
}

int64_t probe_gp(hashtable_t *ht, relation_t *rel, void *output) {
  int64_t matches = 0;
  int16_t k = 0, done = 0, j = 0, stage2_size = 0;
  scalar_state_t state[ScalarStateSize];
  const uint32_t hashmask = ht->hash_mask;
  const uint32_t skipbits = ht->skip_bits;
  chainedtuplebuffer_t *chainedbuf = (chainedtuplebuffer_t *) output;

  // init # of the state
  for (int i = 0; i < ScalarStateSize; ++i) {
    state[i].stage = 1;
  }
  for (uint64_t cur = 0; cur < rel->num_tuples;) {
    // step 1: load tuples from probing table
    for (k = 0; (k < ScalarStateSize) && (cur < rel->num_tuples); ++k) {
      _mm_prefetch((char *)(rel->tuples + cur) + PDIS, _MM_HINT_T0);

      intkey_t idx = HASH(rel->tuples[cur].key, hashmask, skipbits);
      state[k].b = ht->buckets + idx;
      //__builtin_prefetch(state[k].b, 0, 1);
      _mm_prefetch((char * )(state[k].b), _MM_HINT_T0);

      state[k].tuple_id = cur;
      state[k].stage = 0;
      ++cur;
    }
    // step 2: repeating step 2
    stage2_size = k;
    done = 0;
    while (done < stage2_size) {
      for (k = 0; k < stage2_size; ++k) {
        bucket_t *b = state[k].b;
        if (state[k].stage == 1) {
          continue;
        }
        if (b->count == 0) {
          if (state[k].stage == 0) {
            ++done;
          }
          state[k].stage = 1;
          continue;
        }
        if (rel->tuples[state[k].tuple_id].key == b->tuples[0].key) {
          ++matches;
#if WRITE_RESULTS
          /* copy to the result buffer */
          tuple_t *joinres = cb_next_writepos(chainedbuf);
          joinres->key = b->tuples[0].payload; /* R-rid */
          joinres->payload = rel->tuples[state[k].tuple_id].payload; /* S-rid */
#endif
        }
        b = b->next; /* follow overflow pointer */
        if (b) {
          state[k].b = b;
          // __builtin_prefetch(state[k].b, 0, 1);
          _mm_prefetch((char * )(state[k].b), _MM_HINT_T0);

        } else {
          ++done;
          state[k].stage = 1;
        }
      }
    }
  }
  return matches;
}
int64_t probe_AMAC(hashtable_t *ht, relation_t *rel, void *output) {
  int64_t matches = 0;
  int16_t k = 0, done = 0, j = 0;
  scalar_state_t state[ScalarStateSize];
  const uint32_t hashmask = ht->hash_mask;
  const uint32_t skipbits = ht->skip_bits;
  chainedtuplebuffer_t *chainedbuf = (chainedtuplebuffer_t *) output;

  // init # of the state
  for (int i = 0; i < ScalarStateSize; ++i) {
    state[i].stage = 1;
  }
  for (uint64_t cur = 0; (done < ScalarStateSize);) {
    k = (k >= ScalarStateSize) ? 0 : k;

    switch (state[k].stage) {
      case 1: {
        if (cur >= rel->num_tuples) {
          ++done;
          state[k].stage = 3;
          break;
        }
#if SEQPREFETCH
        _mm_prefetch(((char *)(rel->tuples + cur) + PDIS), _MM_HINT_T0);
#endif
        intkey_t idx = HASH(rel->tuples[cur].key, hashmask, skipbits);
        state[k].b = ht->buckets + idx;
        //__builtin_prefetch(state[k].b, 0, 1);
        _mm_prefetch((char * )(state[k].b), _MM_HINT_T0);

        state[k].tuple_id = cur;
        state[k].stage = 0;
        ++cur;
      }
        break;
      case 0: {
        bucket_t *b = state[k].b;
        //  _mm_lfence();
        //#pragma unroll(2)
        if (b->count == 0) {
          state[k].stage = 1;
          --k;
          break;
        }
        if (rel->tuples[state[k].tuple_id].key == b->tuples[0].key) {
          ++matches;
#if WRITE_RESULTS

          /* copy to the result buffer */
          tuple_t *joinres = cb_next_writepos(chainedbuf);
#if SEQPREFETCH
          _mm_prefetch((char *)(((void *)joinres) + PDIS), _MM_HINT_T0);
#endif
          joinres->key = b->tuples[0].payload; /* R-rid */
          joinres->payload = rel->tuples[state[k].tuple_id].payload; /* S-rid */
#endif
        }
        b = b->next; /* follow overflow pointer */
        if (b) {
          state[k].b = b;
          // __builtin_prefetch(state[k].b, 0, 1);
          _mm_prefetch((char * )(state[k].b), _MM_HINT_T0);
        } else {
          state[k].stage = 1;
          --k;
        }
      }
        break;
    }
    ++k;
  }

  return matches;
}

/** print out the execution time statistics of the join */
static void print_timing(uint64_t total, uint64_t build, uint64_t part, uint64_t numtuples, int64_t result, struct timeval *start, struct timeval *end) {
  double diff_usec = (((*end).tv_sec * 1000000L + (*end).tv_usec) - ((*start).tv_sec * 1000000L + (*start).tv_usec));
  double cyclestuple = total;
  cyclestuple /= numtuples;
  fprintf(stdout, "RUNTIME TOTAL, BUILD, PART (cycles): \n");
  fprintf(stderr, "%llu \t %llu \t %llu ", total, build, part);
  fprintf(stdout, "\n");
  fprintf(stdout, "TOTAL-TIME-USECS, TOTAL-TUPLES, CYCLES-PER-TUPLE: \n");
  fprintf(stdout, "%.4lf \t %llu \t ", diff_usec, result);
  fflush(stdout);
  fprintf(stderr, "%.4lf ", cyclestuple);
  fflush(stderr);
  fprintf(stdout, "\n");
}

/** \copydoc NPO_st */
result_t *NPO_st(relation_t *relR, relation_t *relS, int nthreads) {
  hashtable_t *ht;
  int64_t result = 0;
  result_t *joinresult;

#ifndef NO_TIMING
  struct timeval start, end;
  uint64_t timer1, timer2, timer3;
#endif
  uint32_t nbuckets = (relR->num_tuples / BUCKET_SIZE);
  allocate_hashtable(&ht, nbuckets);

  joinresult = (result_t *) malloc(sizeof(result_t));
#ifdef JOIN_RESULT_MATERIALIZE
  joinresult->resultlist = (threadresult_t *) malloc(sizeof(threadresult_t));
#endif

#ifndef NO_TIMING
  gettimeofday(&start, NULL);
  startTimer(&timer1);
  startTimer(&timer2);
  timer3 = 0; /* no partitioning */
#endif

  build_hashtable_st(ht, relR);

#ifndef NO_TIMING
  stopTimer(&timer2); /* for build */
#endif

#ifdef JOIN_RESULT_MATERIALIZE
  chainedtuplebuffer_t *chainedbuf = chainedtuplebuffer_init();
#else
  void *chainedbuf = NULL;
#endif

  result = probe_hashtable(ht, relS, chainedbuf);

#ifdef JOIN_RESULT_MATERIALIZE
  threadresult_t *thrres = &(joinresult->resultlist[0]); /* single-thread */
  thrres->nresults = result;
  thrres->threadid = 0;
  thrres->results = (void *) chainedbuf;
#endif

#ifndef NO_TIMING
  stopTimer(&timer1); /* over all */
  gettimeofday(&end, NULL);
  /* now print the timing results: */
  print_timing(timer1, timer2, timer3, relS->num_tuples, result, &start, &end);
#endif

  destroy_hashtable(ht);

  joinresult->totalresults = result;
  joinresult->nthreads = 1;

  return joinresult;
}

/**
 * Multi-thread hashtable build method, ht is pre-allocated.
 * Writes to buckets are synchronized via latches.
 *
 * @param ht hastable to be built
 * @param rel the build relation
 * @param overflowbuf pre-allocated chunk of buckets for overflow use.
 */
void build_hashtable_mt(hashtable_t *ht, relation_t *rel, bucket_buffer_t **overflowbuf) {
  uint32_t i;
  const uint32_t hashmask = ht->hash_mask;
  const uint32_t skipbits = ht->skip_bits;

#ifdef PREFETCH_NPJ
  size_t prefetch_index = PREFETCH_DISTANCE;
#endif

  for (i = 0; i < rel->num_tuples; i++) {
    tuple_t *dest;
    bucket_t *curr, *nxt;

#ifdef PREFETCH_NPJ
    if (prefetch_index < rel->num_tuples) {
      intkey_t idx_prefetch = HASH(rel->tuples[prefetch_index++].key, hashmask, skipbits);
      __builtin_prefetch(ht->buckets + idx_prefetch, 1, 1);
    }
#endif

    int32_t idx = HASH(rel->tuples[i].key, hashmask, skipbits);
    /* copy the tuple to appropriate hash bucket */
    /* if full, follow nxt pointer to find correct place */
    curr = ht->buckets + idx;
    lock(&curr->latch);
    nxt = curr->next;

    if (curr->count == BUCKET_SIZE) {
      // add new bucket after the first
      if (!nxt || nxt->count == BUCKET_SIZE) {
        bucket_t *b;
        /* b = (bucket_t*) calloc(1, sizeof(bucket_t)); */
        /* instead of calloc() everytime, we pre-allocate */
        get_new_bucket(&b, overflowbuf);
        curr->next = b;
        b->next = nxt;
        b->count = 1;
        dest = b->tuples;
        curr->lenth++;
      } else {
        dest = nxt->tuples + nxt->count;
        nxt->count++;
      }
    } else {
      dest = curr->tuples + curr->count;
      curr->count++;
      curr->lenth = 1;
    }

    *dest = rel->tuples[i];
    unlock(&curr->latch);
  }
}

/**
 * Just a wrapper to call the build and probe for each thread.
 *
 * @param param the parameters of the thread, i.e. tid, ht, reln, ...
 *
 * @return
 */
volatile static char g_lock = 0, g_lock_morse = 0;
volatile static uint64_t total_num = 0, global_curse = 0, global_upper, global_morse_size;
typedef int64_t (*ProbeFun)(hashtable_t *ht, relation_t *rel, void *output);
volatile struct Fun {
  ProbeFun fun_ptr;
  char fun_name[8];
} pfun[10];
volatile static int pf_num = 0;
void morse_driven(void*param, ProbeFun fun, void*output) {
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
    args->num_results += fun(args->ht, &relS, output);
  }
}
bucket_buffer_t *overflowbuf[2];
volatile static hashtable_t* ht[2];
volatile static relation_t rel[2];
void *npo_thread(void *param) {
  int rv;
  arg_t *args = (arg_t *) param;
  struct timeval t1, t2;
  int deltaT = 0, thread_num = 0;
  /* allocate overflow buffer for each thread */

#if TEST_NUMA
  if (args->tid == 0) {
    uint32_t nbuckets = (args->relR.num_tuples / BUCKET_SIZE);
    allocate_hashtable(&ht[args->tid], nbuckets);
    init_bucket_buffer(&overflowbuf[args->tid]);
    build_hashtable_mt(ht[args->tid], &args->relR, &overflowbuf[args->tid]);
    rel[args->tid].num_tuples= args->relS.num_tuples;
    rel[args->tid].tuples = (tuple_t*)alloc_aligned(rel[args->tid].num_tuples*sizeof(tuple_t));
    memcpy(rel[args->tid].tuples,args->relS.tuples,rel[args->tid].num_tuples*sizeof(tuple_t));
  }
  BARRIER_ARRIVE(args->barrier, rv);
  if (args->tid == 1) {
    uint32_t nbuckets = (args->relR.num_tuples / BUCKET_SIZE);
    allocate_hashtable(&ht[args->tid], nbuckets);
    init_bucket_buffer(&overflowbuf[args->tid]);
    build_hashtable_mt(ht[args->tid], &args->relR, &overflowbuf[args->tid]);
    rel[args->tid].num_tuples= args->relS.num_tuples;
    rel[args->tid].tuples = (tuple_t*)alloc_aligned(rel[args->tid].num_tuples*sizeof(tuple_t));
    memcpy(rel[args->tid].tuples,args->relS.tuples,rel[args->tid].num_tuples*sizeof(tuple_t));
  }
  BARRIER_ARRIVE(args->barrier, rv);
  args->ht = ht[(args->tid) % 2];
  args->relS.tuples=rel[(args->tid) % 2].tuples;
  BARRIER_ARRIVE(args->barrier, rv);

#else
  bucket_buffer_t *overflowbuf;
  init_bucket_buffer(&overflowbuf);
  gettimeofday(&t1, NULL);
  /* insert tuples from the assigned part of relR to the ht */
  build_hashtable_mt(args->ht, &args->relR, &overflowbuf);

  /* wait at a barrier until each thread completes build phase */
  BARRIER_ARRIVE(args->barrier, rv);
  if (args->tid == 0) {
    gettimeofday(&t2, NULL);
    deltaT = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
    printf("--------build costs time (ms) = %lf\n", deltaT * 1.0 / 1000);
    print_hashtable(args->ht);
    printf("size of bucket_t = %d\n", sizeof(bucket_t));
  }
#endif

  if (args->tid == 0) {
    puts("+++++sleep begin+++++");
  }
  sleep(SLEEP_TIME);
  if (args->tid == 0) {
    puts("+++++sleep end  +++++");
  }
  if (args->tid == 0) {
    strcpy(pfun[5].fun_name, "IMV");
    strcpy(pfun[4].fun_name, "AMAC");
    strcpy(pfun[3].fun_name, "FVA");
    strcpy(pfun[2].fun_name, "DVA");
    strcpy(pfun[1].fun_name, "SIMD");
    strcpy(pfun[0].fun_name, "Naive");

    pfun[5].fun_ptr = smv_probe;
    pfun[4].fun_ptr = probe_AMAC;
    pfun[3].fun_ptr = probe_simd_amac;
    pfun[2].fun_ptr = probe_simd_amac_raw;
    pfun[1].fun_ptr = probe_simd;
    pfun[0].fun_ptr = probe_hashtable;

    pf_num = 6;
  }
  BARRIER_ARRIVE(args->barrier, rv);

  int warm_up_times = 2;
  ////////// compact, do two branches in the integration

  chainedtuplebuffer_t *chainedbuf_compact = chainedtuplebuffer_init();
  for (int rp = 0; rp < warm_up_times; ++rp) {
    BARRIER_ARRIVE(args->barrier, rv);
    gettimeofday(&t1, NULL);
    args->num_results = probe_hashtable(args->ht, &args->relS, chainedbuf_compact);
    lock(&g_lock);
#if DIVIDE
    total_num += args->num_results;
#else
    total_num = args->num_results;
#endif
    unlock(&g_lock);
    BARRIER_ARRIVE(args->barrier, rv);
    if (args->tid == 0) {
      printf("WARNNP result num = %lld\t", total_num);
      gettimeofday(&t2, NULL);
      deltaT = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
      printf("---- Naive probe costs time (ms) = %lf\n", deltaT * 1.0 / 1000);
      total_num = 0;
    }
  }
  chainedtuplebuffer_free(chainedbuf_compact);
  if (args->tid == 0) {
    puts("+++++sleep begin+++++");
  }
  sleep(SLEEP_TIME);
  if (args->tid == 0) {
    puts("+++++sleep end  +++++");
  }

  BARRIER_ARRIVE(args->barrier, rv);

  for (int fid = 0; fid < pf_num; ++fid) {
    chainedtuplebuffer_t *chainedbuf_compact = chainedtuplebuffer_init();
    for (int rp = 0; rp < REPEAT_PROBE; ++rp) {
      BARRIER_ARRIVE(args->barrier, rv);
      gettimeofday(&t1, NULL);
      deltaT = 0;
      thread_num = 0;
#if MORSE_SIZE
      morse_driven(param, pfun[fid].fun_ptr, chainedbuf_compact);
#else
      args->num_results = pfun[fid].fun_ptr(args->ht, &args->relS, chainedbuf_compact);
#endif
      lock(&g_lock);
#if DIVIDE
      total_num += args->num_results;
#elif MORSE_SIZE
      total_num += args->num_results;
#else
      total_num = args->num_results;
#endif
      unlock(&g_lock);
      BARRIER_ARRIVE(args->barrier, rv);
      if (args->tid == 0) {
        gettimeofday(&t2, NULL);
        printf("total result num = %lld\t", total_num);
        deltaT = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
        printf("---- %5s probe costs time (ms) = %10.4lf\n", pfun[fid].fun_name, deltaT * 1.0 / 1000);
        total_num = 0;
        global_curse = 0;

      }
    }
    chainedtuplebuffer_free(chainedbuf_compact);
    if (args->tid == 0) {
      puts("+++++sleep begin+++++");
    }
    sleep(SLEEP_TIME);
    if (args->tid == 0) {
      puts("+++++sleep end  +++++");
    }
  }

//------------------------------------
#ifdef JOIN_RESULT_MATERIALIZE
  args->threadresult->nresults = args->num_results;
  args->threadresult->threadid = args->tid;
// args->threadresult->results = (void *)chainedbuf;
#endif
  BARRIER_ARRIVE(args->barrier, rv);

#if TEST_NUMA
  if (args->tid == 0) {
    destroy_hashtable(ht[0]);
    destroy_hashtable(ht[1]);
    free_bucket_buffer(overflowbuf[0]);
    free_bucket_buffer(overflowbuf[1]);
    delete_relation(&rel[0]);
    delete_relation(&rel[1]);
  }
#else
  /* clean-up the overflow buffers */
  free_bucket_buffer(overflowbuf);
#endif
  return 0;
}

/*void tbb_run(relation_t *relR, relation_t *relS, int nthreads) {
  int nrThreads = nthreads;

  tbb::task_scheduler_init scheduler(nrThreads);
//  auto resources = initQuery(nrThreads);
  using range = tbb::blocked_range<size_t>;
  const auto add = [](const size_t& a, const size_t& b) {return a + b;};
  bucket_buffer_t *overflowbuf;
  hashtable_t *ht;
  double start, end;
  init_bucket_buffer(&overflowbuf);
  uint32_t nbuckets = (relR->num_tuples / BUCKET_SIZE);
  allocate_hashtable(&ht, nbuckets);
  build_hashtable_mt(ht, relR, &overflowbuf);
  uint64_t found2 = 0;
  int morsesize = MORSE_SIZE;
  if (nrThreads == 1) {
    morsesize = relS->num_tuples;
  }
  auto probe = [&](ProbeFun probe_fun) {
    found2 = tbb::parallel_reduce(range(0, relS->num_tuples, morsesize), 0, [&](const tbb::blocked_range<size_t>& r, const size_t& f) {
          auto found = f;
          relation_t rel;
          chainedtuplebuffer_t *chainedbuf = chainedtuplebuffer_init();

          rel.num_tuples=r.end()-r.begin();
          rel.tuples= relS->tuples+r.begin();
          found+=probe_fun(ht,&rel,chainedbuf);
          chainedtuplebuffer_free(chainedbuf);

          return found;
        },
        add);
  };
  PerfEvents event;
  vector<pair<string, ProbeFun> > agg_name2fun;
  int repetitions = 5;
  agg_name2fun.push_back(make_pair("Naive", probe_hashtable));
  agg_name2fun.push_back(make_pair("SIMD", probe_simd));
  agg_name2fun.push_back(make_pair("DVA", probe_simd_amac_raw));
  agg_name2fun.push_back(make_pair("FVA", probe_simd_amac));
  agg_name2fun.push_back(make_pair("AMAC", probe_AMAC));
  agg_name2fun.push_back(make_pair("IMV", smv_probe));

  for (auto name2fun : agg_name2fun) {
#if 0
    event.timeAndProfile(name2fun.first, 10000,[&]() {probe(name2fun.second);}, repetitions);
#else
    probe(name2fun.second);
    start = gettime();
    for (int j = 0; j < REPEAT_PROBE; j++) {
      probe(name2fun.second);
    }
    end = gettime();
    printf("total result num = %lld\t", found2);
    printf("---- %5s probe costs time (ms) = %10.4lf\n", name2fun.first.c_str(), (end - start) * 1000 / REPEAT_PROBE);
    sleep(5);
#endif
  }
//  event.timeAndProfile("imv_probe",10000,probe,2,0);

  std::cout << "num of results = " << found2 << " threads = " << nrThreads << std::endl;
//  leaveQuery(nrThreads);
  destroy_hashtable(ht);

  scheduler.terminate();
}*/

/** \copydoc NPO */
result_t *NPO(relation_t *relR, relation_t *relS, int nthreads) {
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

#if USE_TBB
//  pthread_attr_init(&attr);
//  for (i = 0; i < nthreads; i++) {
//    int cpu_idx = get_cpu_id(i);
//
//    DEBUGMSG(1, "Assigning thread-%d to CPU-%d\n", i, cpu_idx);
//#if AFFINITY
//    CPU_ZERO(&set);
//    CPU_SET(cpu_idx, &set);
//    pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);
//#endif
//  }

  tbb_run(relR, relS, nthreads);
  joinresult->totalresults = result;
  joinresult->nthreads = nthreads;
  return joinresult;
#endif
  uint32_t nbuckets = (relR->num_tuples / BUCKET_SIZE);
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
  if (nthreads == 1) {
    global_morse_size = relS->num_tuples;
  } else {
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
#if TEST_NUMA
    args[i].relR.num_tuples = relR->num_tuples;
    args[i].relR.tuples = relR->tuples;
#else
    /* assing part of the relR for next thread */
    args[i].relR.num_tuples = (i == (nthreads - 1)) ? numR : numRthr;
    args[i].relR.tuples = relR->tuples + numRthr * i;
    numR -= numRthr;
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

    rv = pthread_create(&tid[i], &attr, npo_thread, (void *) &args[i]);
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

#ifndef NO_TIMING
  /* now print the timing results: */
  print_timing(args[0].timer1, args[0].timer2, args[0].timer3, relS->num_tuples, result, &args[0].start, &args[0].end);
#endif

  destroy_hashtable(ht);

  return joinresult;
}

/** @}*/

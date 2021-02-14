/**
 * @file    no_partitioning_join.h
 * @author  Cagri Balkesen <cagri.balkesen@inf.ethz.ch>
 * @date    Sun Feb  5 20:12:56 2012
 * @version $Id: no_partitioning_join.h 4419 2013-10-21 16:24:35Z bcagri $
 *
 * @brief  The interface of No partitioning optimized (NPO) join algorithm.
 *
 * (c) 2012, ETH Zurich, Systems Group
 *
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#ifndef NO_PARTITIONING_JOIN_H
#define NO_PARTITIONING_JOIN_H
#include "types.hpp" /* relation_t */
#include "prefetch.hpp"
#include <sys/types.h>
#include <unistd.h>
#include <sched.h>    /* CPU_ZERO, CPU_SET */
#include <pthread.h>  /* pthread_* */
#include <string.h>   /* memset */
#include <stdio.h>    /* printf */
#include <stdlib.h>   /* memalign */
#include <sys/time.h> /* gettimeofday */

#include <smmintrin.h>
#include "npj_params.hpp"  /* constant parameters */
#include "npj_types.hpp"   /* bucket_t, hashtable_t, bucket_buffer_t */
#include "rdtsc.h"       /* startTimer, stopTimer */
#include "lock.h"        /* lock, unlock */
#include "cpu_mapping.h" /* get_cpu_id */
#ifdef PERF_COUNTERS
#include "perf_counters.h" /* PCM_x */
#endif

#include "barrier.h"   /* pthread_barrier_* */
#include "affinity.h"  /* pthread_attr_setaffinity_np */
#include "generator.h" /* numa_localize() */

#ifdef JOIN_RESULT_MATERIALIZE
#include "tuple_buffer.hpp" /* for materialization */
#endif
#define NO_TIMING 1
#ifndef BARRIER_ARRIVE
/** barrier wait macro */
#define BARRIER_ARRIVE(B, RV)                           \
  RV = pthread_barrier_wait(B);                         \
  if (RV != 0 && RV != PTHREAD_BARRIER_SERIAL_THREAD) { \
    printf("Couldn't wait on barrier\n");               \
    exit(EXIT_FAILURE);                                 \
  }
#endif

#ifndef NEXT_POW_2
/**
 *  compute the next number, greater than or equal to 32-bit unsigned v.
 *  taken from "bit twiddling hacks":
 *  http://graphics.stanford.edu/~seander/bithacks.html
 */
#define NEXT_POW_2(V) \
  do {                \
    V--;              \
    V |= V >> 1;      \
    V |= V >> 2;      \
    V |= V >> 4;      \
    V |= V >> 8;      \
    V |= V >> 16;     \
    V++;              \
  } while (0)
#endif

#ifndef HASH
#define HASH(X, MASK, SKIP) (((X)&MASK) >> SKIP)
#endif

/** Debug msg logging method */
#ifndef DEBUG
#define DEBUGMSG(COND, MSG, ...)                    \
  if (COND) {                                       \
    fprintf(stdout, "[DEBUG] " MSG, ##__VA_ARGS__); \
  }
#else
#define DEBUGMSG(COND, MSG, ...)
#endif

/** An experimental feature to allocate input relations numa-local */
extern int numalocalize; /* defined in generator.c */
extern int nthreads;     /* defined in generator.c */

/**
 * \ingroup NPO arguments to the threads
 */
typedef struct npj_arg_t arg_t;

struct npj_arg_t {
  int32_t tid;
  hashtable_t *ht;
  relation_t relR;
  relation_t relS;
  pthread_barrier_t *barrier;
  int64_t num_results;

  /* results of the thread */
  threadresult_t *threadresult;

#ifndef NO_TIMING
  /* stats about the thread */
  uint64_t timer1, timer2, timer3;
  struct timeval start, end;
#endif
};

/**
 * NPO: No Partitioning Join Optimized.
 *
 * The "No Partitioning Join Optimized" implementation denoted as NPO
 * which was originally proposed by Blanas et al. in SIGMOD 2011.
 *
 * The following is a multi-threaded implementation. Just returns the
 * number of result tuples.
 *
 * @param relR input relation R - inner relation
 * @param relS input relation S - outer relation
 *
 * @return number of result tuples
 */
result_t *NPO(relation_t *relR, relation_t *relS, int nthreads);

/**
 * The No Partitioning Join Optimized (NPO) as a single-threaded
 * implementation. Just returns the number of result tuples.
 *
 * @param relR input relation R - inner relation
 * @param relS input relation S - outer relation
 *
 * @return number of result tuples
 */
result_t *NPO_st(relation_t *relR, relation_t *relS, int nthreads);
result_t *PIPELINE(relation_t *relR, relation_t *relS, int nthreads);
result_t *BTS(relation_t *relR, relation_t *relS, int nthreads);
void allocate_hashtable(hashtable_t **ppht, uint32_t nbuckets);
void destroy_hashtable(hashtable_t *ht);
void get_new_bucket(bucket_t **result, bucket_buffer_t **buf);
void init_bucket_buffer(bucket_buffer_t **ppbuf);
void free_bucket_buffer(bucket_buffer_t *buf);
void print_hashtable(hashtable_t * const ht);
#endif /* NO_PARTITIONING_JOIN_H */

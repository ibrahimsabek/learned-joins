#ifndef TREE_H
#define TREE_H
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
#include "barrier.h"     /* pthread_barrier_* */
#include "affinity.h"    /* pthread_attr_setaffinity_np */
#include "generator.h"   /* numa_localize() */
#include "../profile.hpp"

typedef struct tnode_t tnode_t;
typedef struct tnodebuffer_t tnodebuffer_t;
typedef struct chainedtnodebuffer_t chainedtnodebuffer_t;
typedef struct tree_t tree_t;
typedef struct tree_arg_t tree_arg_t;
typedef struct tree_state_t tree_state_t;
#define TNODEBUFF_NUMTUPLESPERBUF (1024 * 1024)
struct tnode_t {
  intkey_t key;
  value_t payload;
  tnode_t *lnext, *rnext;
};
struct tnodebuffer_t {
  tnode_t *tnode;
  tnodebuffer_t *next;
};
struct tree_t {
  tnode_t *first_node;
  uint64_t num;
  chainedtnodebuffer_t *buffer;
};
struct chainedtnodebuffer_t {
  tnodebuffer_t *buf;
  tnodebuffer_t *readcursor;
  tnodebuffer_t *writecursor;
  uint32_t writepos;
  uint32_t readpos;
  uint32_t readlen;
  uint32_t numbufs;
};

struct tree_arg_t {
  int32_t tid;
  tree_t *tree;
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
struct tree_state_t {
  int64_t tuple_id;
  tnode_t *b;
  int16_t stage;
};
static inline tnode_t *nb_next_writepos(chainedtnodebuffer_t *cb) {
  if (cb->writepos == TNODEBUFF_NUMTUPLESPERBUF) {
    tnodebuffer_t *newbuf = (tnodebuffer_t *)malloc(sizeof(tnodebuffer_t));
    posix_memalign((void **)&newbuf->tnode, CACHE_LINE_SIZE,
                   sizeof(tnode_t) * TNODEBUFF_NUMTUPLESPERBUF);
    memset(newbuf->tnode, 0, sizeof(tnode_t) * TNODEBUFF_NUMTUPLESPERBUF);
    newbuf->next = cb->buf;
    cb->buf = newbuf;
    cb->numbufs++;
    cb->writepos = 0;
  }
  return (cb->buf->tnode + cb->writepos++);
}

static chainedtnodebuffer_t *chainedtnodebuffer_init(void) {
  chainedtnodebuffer_t *newcb =
      (chainedtnodebuffer_t *)malloc(sizeof(chainedtnodebuffer_t));
  tnodebuffer_t *newbuf = (tnodebuffer_t *)malloc(sizeof(tnodebuffer_t));

  if (posix_memalign((void **)&(newbuf->tnode), CACHE_LINE_SIZE,
                     sizeof(tnode_t) * TNODEBUFF_NUMTUPLESPERBUF)) {
    perror("Aligned allocation failed!\n");
    exit(EXIT_FAILURE);
  }
  memset(newbuf->tnode, 0, sizeof(tnode_t) * TNODEBUFF_NUMTUPLESPERBUF);
  newbuf->next = NULL;
  newcb->buf = newcb->readcursor = newcb->writecursor = newbuf;
  newcb->writepos = newcb->readpos = 0;
  newcb->numbufs = 1;

  return newcb;
}

static void chainedtnodebuffer_free(chainedtnodebuffer_t *cb) {
  tnodebuffer_t *tmp = cb->buf;
  while (tmp) {
    tnodebuffer_t *tmp2 = tmp->next;
    free(tmp->tnode);
    free(tmp);
    tmp = tmp2;
  }

  free(cb);
  cb = NULL;
}
int64_t bts_smv(tree_t *tree, relation_t *rel, void *output);
int64_t bts_simd_amac(tree_t *tree, relation_t *rel, void *output);
int64_t bts_simd_amac_raw(tree_t *tree, relation_t *rel, void *output);
int64_t bts_simd(tree_t *tree, relation_t *rel, void *output);
int64_t bts_smv(tree_t *tree, relation_t *rel, void *output);


#endif

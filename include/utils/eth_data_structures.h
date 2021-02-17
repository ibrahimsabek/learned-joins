#pragma once

#include <unordered_map> 

#include "configs/eth_configs.h"
#include "utils/math.h"
#include "utils/memory.h"
#include "utils/data_structures.h"
#include "utils/eth_generic_task_queue.h"
#include "utils/learned_sort_for_sort_merge.h"

using namespace learned_sort_for_sort_merge;

/*********** Common data structures for ETH radix join ***********/
template<typename KeyType, typename PayloadType, typename TaskType>
struct ETHRadixJoinThread : JoinThreadBase<KeyType, PayloadType, TaskType> {
    int32_t ** histR;
    Tuple<KeyType, PayloadType> *  tmpR;
    int32_t ** histS;
    Tuple<KeyType, PayloadType> *  tmpS;
#ifdef SKEW_HANDLING
    TaskQueue<KeyType, PayloadType, TaskType> *      skew_queue;
    Task<KeyType, PayloadType> **        skewtask;
#endif
    /* stats about the thread */
    int32_t        parts_processed;
} __attribute__((aligned(CACHE_LINE_SIZE)));


template<typename KeyType, typename PayloadType, typename TaskType, typename JoinThreadType>
struct ETHPartition : PartitionBase<KeyType, PayloadType, TaskType, JoinThreadType> {
    Tuple<KeyType, PayloadType> *  tmp;
    int32_t ** hist;
    int64_t *  output;
    int32_t    R;
    uint32_t   D;
    uint32_t   padding;
} __attribute__((aligned(CACHE_LINE_SIZE)));

struct ETHBucketChainingBuild : BuildBase {
    uint32_t * next;
};

/**************** Common data structures for ETH non partition hash join ******************/

/*********** Hashtable and its buckets ******************************/
/********************************************************************/
#if PADDED_BUCKET==0
/** 
 * Normal hashtable buckets.
 *
 * if KEY_8B then key is 8B and sizeof(bucket_t) = 48B
 * else key is 16B and sizeof(bucket_t) = 32B
 */
template<typename KeyType, typename PayloadType>
struct Bucket {
    volatile char     latch;
    /* 3B hole */
    uint32_t          count;
    Bucket<KeyType, PayloadType> * next;
    Tuple<KeyType, PayloadType>  tuples[BUCKET_SIZE];
};
#else /* PADDED_BUCKET: bucket is padded to cache line size */
/** 
 * Cache-sized bucket where size of the bucket is padded
 * to cache line size (64B). 
 */
template<typename KeyType, typename PayloadType>
struct Bucket {
    volatile char     latch;
    /* 3B hole */
    uint32_t          count;
    Tuple<KeyType, PayloadType>  tuples[BUCKET_SIZE];
    Bucket<KeyType, PayloadType> * next;
} __attribute__ ((aligned(CACHE_LINE_SIZE)));
#endif /* PADDED_BUCKET */

/** Hashtable structure for NPO. */
template<typename KeyType, typename PayloadType>
struct Hashtable {
    Bucket<KeyType, PayloadType> * buckets;
    int32_t    num_buckets;
    uint32_t   hash_mask;
    uint32_t   skip_bits;
};

/** Pre-allocated bucket buffers are used for overflow-buckets. */
template<typename KeyType, typename PayloadType>
struct BucketBuffer {
    BucketBuffer<KeyType, PayloadType> * next;
    uint32_t count;
    Bucket<KeyType, PayloadType> buf[OVERFLOW_BUF_SIZE];
};

template<typename KeyType, typename PayloadType>
void 
init_bucket_buffer(BucketBuffer<KeyType, PayloadType> ** ppbuf)
{
    BucketBuffer<KeyType, PayloadType> * overflowbuf;
    overflowbuf = (BucketBuffer<KeyType, PayloadType>*) malloc(sizeof(BucketBuffer<KeyType, PayloadType>));
    overflowbuf->count = 0;
    overflowbuf->next  = NULL;

    *ppbuf = overflowbuf;
}

template<typename KeyType, typename PayloadType>
static inline void 
get_new_bucket(Bucket<KeyType, PayloadType> ** result, BucketBuffer<KeyType, PayloadType> ** buf)
{
    if((*buf)->count < OVERFLOW_BUF_SIZE) {
        *result = (*buf)->buf + (*buf)->count;
        (*buf)->count ++;
    }
    else {
        /* need to allocate new buffer */
        BucketBuffer<KeyType, PayloadType> * new_buf = (BucketBuffer<KeyType, PayloadType>*) 
                                                        malloc(sizeof(BucketBuffer<KeyType, PayloadType>));
        new_buf->count = 1;
        new_buf->next  = *buf;
        *buf    = new_buf;
        *result = new_buf->buf;
    }
}

template<typename KeyType, typename PayloadType>
void
free_bucket_buffer(BucketBuffer<KeyType, PayloadType> * buf)
{
    do {
        BucketBuffer<KeyType, PayloadType> * tmp = buf->next;
        free(buf);
        buf = tmp;
    } while(buf);
}

template<typename KeyType, typename PayloadType>
void 
allocate_hashtable(Hashtable<KeyType, PayloadType> ** ppht, uint32_t nbuckets)
{
    Hashtable<KeyType, PayloadType> * ht;

    ht              = (Hashtable<KeyType, PayloadType>*)malloc(sizeof(Hashtable<KeyType, PayloadType>));
    ht->num_buckets = nbuckets;
    NEXT_POW_2((ht->num_buckets));

    /* allocate hashtable buckets cache line aligned */
    if (posix_memalign((void**)&ht->buckets, CACHE_LINE_SIZE,
                       ht->num_buckets * sizeof(Bucket<KeyType, PayloadType>))){
        perror("Aligned allocation failed!\n");
        exit(EXIT_FAILURE);
    }

    memset(ht->buckets, 0, ht->num_buckets * sizeof(Bucket<KeyType, PayloadType>));
    ht->skip_bits = 0; /* the default for modulo hash */
    ht->hash_mask = (ht->num_buckets - 1) << ht->skip_bits;
    *ppht = ht;
}

template<typename KeyType, typename PayloadType>
void 
destroy_hashtable(Hashtable<KeyType, PayloadType> * ht)
{
    free(ht->buckets);
    free(ht);
}

/***************************************************************/

template<typename KeyType, typename PayloadType, typename TaskType>
struct ETHNonPartitionJoinThread {
    int32_t             tid;
    Hashtable<KeyType, PayloadType> *  ht;
    Relation<KeyType, PayloadType>     relR;
    Relation<KeyType, PayloadType>     relS;
    pthread_barrier_t * barrier;
    int64_t             num_results = 0;

    /* results of the thread */
    ThreadResult * threadresult;
    
    /**** start stuff for learning RMI models ****/
    learned_sort_for_sort_merge::RMI<KeyType, PayloadType> * rmi;
    typename learned_sort_for_sort_merge::RMI<KeyType, PayloadType>::Params p;
    Relation<KeyType, PayloadType> *     original_relR;
    Relation<KeyType, PayloadType> *     original_relS;
    Tuple<KeyType, PayloadType> *   relR_start_sampling_ptr;
    Tuple<KeyType, PayloadType> *   relR_end_sampling_ptr;
    Tuple<KeyType, PayloadType> *   relS_start_sampling_ptr;
    Tuple<KeyType, PayloadType> *   relS_end_sampling_ptr;
    Tuple<KeyType, PayloadType> * tmp_training_sample_in;
    Tuple<KeyType, PayloadType> * sorted_training_sample_in;
    Tuple<KeyType, PayloadType> * r_tmp_training_sample_in;
    Tuple<KeyType, PayloadType> * r_sorted_training_sample_in;
    Tuple<KeyType, PayloadType> * s_tmp_training_sample_in;
    Tuple<KeyType, PayloadType> * s_sorted_training_sample_in;
    vector<vector<vector<training_point<KeyType, PayloadType>>>> * training_data;
    uint32_t tmp_training_sample_R_offset, tmp_training_sample_S_offset, tmp_training_sample_offset;
    uint32_t * sample_count, * sample_count_R, * sample_count_S;
    /**** end stuff for learning RMI models ****/

    /* stats about the thread */
    struct timeval start_time, partition_end_time, end_time;
#ifndef DEVELOPMENT_MODE
    PerfEvent e_start_to_partition, e_partition_to_end;
#endif
#ifdef DEVELOPMENT_MODE
    unordered_map<uint64_t, uint64_t> * build_hash_bucket_visits;
    unordered_map<uint64_t, uint64_t> * probe_hash_bucket_visits;
    volatile char *    keys_hash_latch;
    vector<KeyType> * build_keys_list;
    vector<uint64_t> * build_keys_hash_list;
    vector<KeyType> * probe_keys_list;
    vector<uint64_t> * probe_keys_hash_list;    
#endif     
};



template<typename KeyType, typename PayloadType>
struct ETHNonPartitionJoinBuild {
    Hashtable<KeyType, PayloadType> * ht;
    BucketBuffer<KeyType, PayloadType> ** overflowbuf;
#ifdef DEVELOPMENT_MODE
    unordered_map<uint64_t, uint64_t> * build_hash_bucket_visits;
    unordered_map<uint64_t, uint64_t> * probe_hash_bucket_visits;
    volatile char *    keys_hash_latch;    
    vector<KeyType> * build_keys_list;
    vector<uint64_t> * build_keys_hash_list;
    vector<KeyType> * probe_keys_list;
    vector<uint64_t> * probe_keys_hash_list; 
#endif    
};

/**************** Common data structures for ETH non partition hash join ******************/

/**
 * Various NUMA shuffling strategies for data shuffling phase of join
 * algorithms as also described by NUMA-aware data shuffling paper [CIDR'13].
 *
 * NUMA_SHUFFLE_RANDOM, NUMA_SHUFFLE_RING, NUMA_SHUFFLE_NEXT
 */
enum numa_strategy_t {RANDOM, RING, NEXT};

// Join configuration parameters.
struct joinconfig_t {
    int NTHREADS;
    int PARTFANOUT;
    int LEARNEDSORT;
    int SCALARMERGE;
    int MWAYMERGEBUFFERSIZE;
    enum numa_strategy_t NUMASTRATEGY;
};

template<typename KeyType, typename PayloadType>
struct ETHSortMergeMultiwayJoinThread {
    Tuple<KeyType, PayloadType> *  relR;
    Tuple<KeyType, PayloadType> *  relS;

    // temporary relations for partitioning output 
    Tuple<KeyType, PayloadType> *  tmp_partR;
    Tuple<KeyType, PayloadType> *  tmp_partS;

    // temporary relations for sorting output
    Tuple<KeyType, PayloadType> *  tmp_sortR;
    Tuple<KeyType, PayloadType> *  tmp_sortS;

    int32_t numR;
    int32_t numS;

    int32_t my_tid;
    int     nthreads;

     // join configuration parameters:
    joinconfig_t * joincfg;

    pthread_barrier_t * barrier;
    int64_t result;

    RelationPair<KeyType, PayloadType> ** threadrelchunks;

    // used for multi-way merging, shared by active threads in each NUMA.
    Tuple<KeyType, PayloadType> ** sharedmergebuffer;

    // arguments specific to mpsm-join:
    uint32_t ** histR;
    //Tuple<KeyType, PayloadType> * tmpRglobal;
    //uint64_t totalR;

    ThreadResult * threadresult;

#ifdef SKEW_HANDLING
    // skew handling task queues (1 per NUMA region).
    taskqueue_t ** numa_taskqueues;
    pthread_mutex_t* numa_taskqueues_locks;
    int* is_numa_taskqueues_created;
#endif

    struct timeval start_time, partition_end_time, sort_end_time, tmp_sort_end_time, multiwaymerge_end_time, mergejoin_end_time, tmp_mergejoin_end_time;
#ifndef DEVELOPMENT_MODE
    PerfEvent e_start_to_partition, e_partition_to_sort, e_sort_to_multiwaymerge, e_multiwaymerge_to_mergejoin;
#endif     

}__attribute__((aligned(CACHE_LINE_SIZE)));

template<typename KeyType, typename PayloadType>
struct MergeNode {
    Tuple<KeyType, PayloadType> * buffer;
    volatile uint32_t count;
    volatile uint32_t head;
    volatile uint32_t tail;
} __attribute__((packed));

// This is a struct used for representing merge tasks when skew handling 
//    mechanism is enabled. Essentially, large merge tasks are decomposed into
//    smaller merge tasks and placed into a task queue.
template<typename KeyType, typename PayloadType>
struct MergeTask {
    Tuple<KeyType, PayloadType> * output;
    Relation<KeyType, PayloadType> ** runstomerge;
    int numruns;
    unsigned int totaltuples;
    // if heavy-hitter then not merged, directly copied
    int isheavyhitter; 
};

template<typename KeyType, typename PayloadType>
struct LearnedSortMergeMultiwayJoinThread : ETHSortMergeMultiwayJoinThread<KeyType, PayloadType> 
{
    Tuple<KeyType, PayloadType> * tmp_minor_bckts_r;
    Tuple<KeyType, PayloadType> * tmp_minor_bckts_s;
    int64_t * tmp_minor_bckt_sizes_r;
    int64_t * tmp_minor_bckt_sizes_s;
    Tuple<KeyType, PayloadType> * tmp_spill_bucket_r;
    Tuple<KeyType, PayloadType> * sorted_spill_bucket_r;
    Tuple<KeyType, PayloadType> * tmp_spill_bucket_s;
    Tuple<KeyType, PayloadType> * sorted_spill_bucket_s;
    Tuple<KeyType, PayloadType> * tmp_repeatedKeysPredictedRanksR;
    Tuple<KeyType, PayloadType> * tmp_repeatedKeysPredictedRanksS;
    int64_t * tmp_repeatedKeysPredictedRanksCountsR;
    int64_t * tmp_repeatedKeysPredictedRanksCountsS;
    int64_t tmp_repeatedKeysCountsR;
    int64_t tmp_repeatedKeysCountsS;
    int64_t tmp_total_repeatedKeysCountsR;
    int64_t tmp_total_repeatedKeysCountsS;
    unsigned int NUM_MINOR_BCKT_PER_MAJOR_BCKT_r;
    unsigned int NUM_MINOR_BCKT_PER_MAJOR_BCKT_s;
    unsigned int MINOR_BCKTS_OFFSET_r;
    unsigned int MINOR_BCKTS_OFFSET_s;
    unsigned int TOT_NUM_MINOR_BCKTS_r;
    unsigned int TOT_NUM_MINOR_BCKTS_s;
    unsigned int INPUT_SZ_r;
    unsigned int INPUT_SZ_s;

    learned_sort_for_sort_merge::RMI<KeyType, PayloadType> * rmi_r;
    learned_sort_for_sort_merge::RMI<KeyType, PayloadType> * rmi_s;

    struct timeval sample_end_time;
#ifndef DEVELOPMENT_MODE
    PerfEvent e_start_to_sample, e_sample_to_partition, e_sort_to_mergejoin;
#endif
} __attribute__((aligned(CACHE_LINE_SIZE)));

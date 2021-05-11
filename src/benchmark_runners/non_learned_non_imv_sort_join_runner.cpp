#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "emmintrin.h"
#include "immintrin.h"
#include "smmintrin.h"
#include <sys/time.h> /* gettimeofday */

#include "config.h"            /* autoconf header */
#include "configs/base_configs.h"
#include "configs/eth_configs.h"

#include "utils/eth_data_structures.h"
#include "utils/data_generation.h"
#include "utils/io.h"
#include "utils/eth_generic_task_queue.h"
#include "utils/cpu_mapping.h"
#include "utils/base_utils.h"
#include "utils/math.h"
#include "utils/barrier.h"
#include "utils/memory.h"
#include "utils/lock.h" 
#include "utils/learned_sort_for_sort_merge.h"

#include "techniques/sortmerge_multiway_join_eth_steps.h"

#include <chrono>
using namespace chrono;

#ifndef KeyType
#define KeyType RELATION_KEY_TYPE
#define PayloadType RELATION_PAYLOAD_TYPE
#define NUM_THREADS NUM_THREADS_FOR_EVALUATION
#endif


//#define RUN_NUMS 1 //10 

using namespace std;
using namespace learned_sort_for_sort_merge;

ETHSortMergeMultiwayJoinSteps<KeyType, PayloadType, 
                                ETHSortMergeMultiwayJoinThread<KeyType, PayloadType>> join_steps;
joinconfig_t joincfg;
int CACHELINEPADDING;
int RELATION_PADDING;
Tuple<KeyType, PayloadType> * tmpRelpartR;
Tuple<KeyType, PayloadType> * tmpRelpartS;
Tuple<KeyType, PayloadType> * tmpRelsortR;
Tuple<KeyType, PayloadType> * tmpRelsortS;
RelationPair<KeyType, PayloadType> ** threadrelchunks;
uint32_t ** histR;
Tuple<KeyType, PayloadType> ** ptrs_to_sharedmergebufs;
#ifdef SKEW_HANDLING
taskqueue_t ** ptrs_to_taskqueues;
pthread_mutex_t* ptrs_to_taskqueues_locks;
int* ptrs_to_is_numa_taskqueues_created;
#endif


void initialize_non_learned_non_imv_sort_join_thread_args(Relation<KeyType, PayloadType> * rel_r, 
                        Relation<KeyType, PayloadType> * rel_s, 
                        pthread_barrier_t* barrier_ptr,
                        Result * joinresult,
                        ETHSortMergeMultiwayJoinThread<KeyType, PayloadType> * args)
{
    int32_t numperthr[2];
    numperthr[0] = rel_r->num_tuples / NUM_THREADS;
    numperthr[1] = rel_s->num_tuples / NUM_THREADS;

    int i;
    for(i = 0; i < NUM_THREADS; i++)
    {
        (*(args + i)).relR = rel_r->tuples + i * (numperthr[0]);
        (*(args + i)).relS = rel_s->tuples + i * (numperthr[1]);

        /* temporary relations */
        (*(args + i)).tmp_partR = tmpRelpartR + i * (numperthr[0] + CACHELINEPADDING);
        (*(args + i)).tmp_partS = tmpRelpartS + i * (numperthr[1] + CACHELINEPADDING);
        (*(args + i)).tmp_sortR = tmpRelsortR + i * (numperthr[0]);
        (*(args + i)).tmp_sortS = tmpRelsortS + i * (numperthr[1]);

        (*(args + i)).numR = (i == (NUM_THREADS-1)) ?
            (rel_r->num_tuples - i * numperthr[0]) : numperthr[0];
        (*(args + i)).numS = (i == (NUM_THREADS-1)) ?
            (rel_s->num_tuples - i * numperthr[1]) : numperthr[1];

        (*(args + i)).my_tid        = i;/* this is the logical CPU-ID */
        (*(args + i)).nthreads      = NUM_THREADS;
        (*(args + i)).joincfg       = &joincfg;
        (*(args + i)).barrier       = barrier_ptr;
        (*(args + i)).threadrelchunks = threadrelchunks;
        (*(args + i)).sharedmergebuffer = ptrs_to_sharedmergebufs;

        /** information specific to mpsm-join */
        (*(args + i)).histR         = histR;
        //(*(args + i)).tmpRglobal    = tmpRelpartR;
        //(*(args + i)).totalR        = relR->num_tuples;
        (*(args + i)).threadresult  = &(joinresult->resultlist[i]);

#ifdef SKEW_HANDLING
        /** skew handling task queue ptrs. */
        (*(args + i)).numa_taskqueues     = ptrs_to_taskqueues;
        (*(args + i)).numa_taskqueues_locks = ptrs_to_taskqueues_locks;
        (*(args + i)).is_numa_taskqueues_created = ptrs_to_is_numa_taskqueues_created;
#endif
    }
}


void * non_learned_non_imv_sort_join_thread(void * param)
{
    ETHSortMergeMultiwayJoinThread<KeyType, PayloadType> * args   = (ETHSortMergeMultiwayJoinThread<KeyType, PayloadType> *) param;
    int32_t my_tid = args->my_tid;
    int rv;
    int deltaT = 0; struct timeval t1, t2;
    
    /*************************************************************************
    *
    *   Phase.1) NUMA-local partitioning.
    *
    *************************************************************************/
    Relation<KeyType, PayloadType> ** partsR = NULL;
    Relation<KeyType, PayloadType> ** partsS = NULL;
    
    auto partition_t1 = high_resolution_clock::now();
    auto partition_t2 = high_resolution_clock::now();
    vector<uint32_t> curr_partition_timings_in_ms;
    vector<uint32_t> final_partition_timings_in_ms;
    for (int rp = 0; rp < RUN_NUMS; ++rp) 
    {
        //DEBUGMSG(1, "Thread-%d started running ... \n", my_tid);

        // wait at a barrier until each thread starts and start timer
        BARRIER_ARRIVE(args->barrier, rv);

        // the first thread checkpoints the start time
        if(my_tid == 0){ 
            //gettimeofday(&args->start_time, NULL);
            partition_t1 = high_resolution_clock::now();
        }
        
        join_steps.partition_phase(&partsR, &partsS, args);


        // wait at a barrier until each thread completes the partition phase
        BARRIER_ARRIVE(args->barrier, rv);

        // partition phase finished, thread-0 checkpoints the time
        if(my_tid == 0){
            //gettimeofday(&args->partition_end_time, NULL);
            partition_t2 = high_resolution_clock::now();
            deltaT = std::chrono::duration_cast<std::chrono::microseconds>(partition_t2 - partition_t1).count();
            //deltaT = (args->partition_end_time.tv_sec - args->start_time.tv_sec) * 1000000 + args->partition_end_time.tv_usec - args->start_time.tv_usec;
            printf("---- %5s partitioning costs time (ms) = %10.4lf\n", "Non-learned sort join", deltaT * 1.0 / 1000);
            curr_partition_timings_in_ms.push_back((uint32_t)(deltaT * 1.0 / 1000)); //ms
        }
    
        if(!(rp == (RUN_NUMS - 1))){
            //TODO: make sure you can rerun here
            //if(args->tid == 0)
            //    destroy_hashtable(args->ht);

            //free_bucket_buffer(overflowbuf);

            BARRIER_ARRIVE(args->barrier, rv);
        } 
    }

    if(args->tid == 0){
        std::sort(curr_partition_timings_in_ms.begin(), curr_partition_timings_in_ms.end());
        final_partition_timings_in_ms.push_back(curr_partition_timings_in_ms[(int)(curr_partition_timings_in_ms.size()/2)]);
    }

    BARRIER_ARRIVE(args->barrier, rv);

    /*************************************************************************
    *
    *   Phase.2) NUMA-local sorting of cache-sized chunks
    *
    *************************************************************************/
    auto sorting_t1 = high_resolution_clock::now();
    auto sorting_t2 = high_resolution_clock::now();
    vector<uint32_t> curr_sorting_timings_in_ms;
    vector<uint32_t> final_sorting_timings_in_ms;
    for (int rp = 0; rp < RUN_NUMS; ++rp) 
    {
        BARRIER_ARRIVE(args->barrier, rv);

        // the first thread checkpoints the start time
        if(my_tid == 0){ 
            //gettimeofday(&args->partition_end_time, NULL);
            sorting_t1 = high_resolution_clock::now();
        }

        join_steps.sorting_phase(partsR, partsS, args);

        BARRIER_ARRIVE(args->barrier, rv);

        // partition phase finished, thread-0 checkpoints the time
        if(my_tid == 0){
            //gettimeofday(&args->sort_end_time, NULL);
            sorting_t2 = high_resolution_clock::now();
            deltaT = std::chrono::duration_cast<std::chrono::microseconds>(sorting_t2 - sorting_t1).count();
            //deltaT = (args->sort_end_time.tv_sec - args->partition_end_time.tv_sec) * 1000000 + args->sort_end_time.tv_usec - args->partition_end_time.tv_usec;
            printf("---- %5s sorting costs time (ms) = %10.4lf\n", "Non-learned sort join", deltaT * 1.0 / 1000);
            curr_sorting_timings_in_ms.push_back((uint32_t)(deltaT * 1.0 / 1000)); //ms
        }
    
        if(!(rp == (RUN_NUMS - 1))){
            //TODO: make sure you can rerun here
            //if(args->tid == 0)
            //    destroy_hashtable(args->ht);

            //free_bucket_buffer(overflowbuf);

            BARRIER_ARRIVE(args->barrier, rv);
        }
    }

    if(args->tid == 0){
        std::sort(curr_sorting_timings_in_ms.begin(), curr_sorting_timings_in_ms.end());
        final_sorting_timings_in_ms.push_back(curr_sorting_timings_in_ms[(int)(curr_sorting_timings_in_ms.size()/2)]);
    }

    /**
     * Allocate shared merge buffer for multi-way merge tree.
     * This buffer is further divided into given number of threads
     * active in the same NUMA-region.
     *
     * @note the first thread in each NUMA region allocates the shared L3 buffer.
     */
    int numaregionid = get_numa_region_id(my_tid);
#ifdef SKEW_HANDLING    
    pthread_mutex_lock(&(args->numa_taskqueues_locks[numaregionid]));
    if(!args->is_numa_taskqueues_created[numaregionid]) {
#endif
    //if(is_first_thread_in_numa_region(my_tid)) {
        /* TODO: make buffer size runtime parameter */
        Tuple<KeyType, PayloadType> * sharedmergebuffer = (Tuple<KeyType, PayloadType> *)
                alloc_aligned(args->joincfg->MWAYMERGEBUFFERSIZE);
        args->sharedmergebuffer[numaregionid] = sharedmergebuffer;
        /*    
        DEBUGMSG(1, "Thread-%d allocated %.3lf KiB merge buffer in "\
                "NUMA-region-%d to be used by %d active threads.\n",
                my_tid, (double)(args->joincfg->MWAYMERGEBUFFERSIZE/1024.0), 
                numaregionid, get_num_active_threads_in_numa(numaregionid));
        */
#ifdef SKEW_HANDLING
        /* initialize skew handling ; mwaytask_t taskqueues */
        args->numa_taskqueues[numaregionid] = taskqueue_init(32);
#endif
    //}
#ifdef SKEW_HANDLING    
        args->is_numa_taskqueues_created[numaregionid] = 1;
    }
    pthread_mutex_unlock(&(args->numa_taskqueues_locks[numaregionid]));
#endif

    BARRIER_ARRIVE(args->barrier, rv);

    /*************************************************************************
     *
     *   Phase.3) Apply multi-way merging with in-cache resident buffers +
     *   Detect skew & handle skew in a fine-grained manner for a robust perf.
     *
     *************************************************************************/
    Relation<KeyType, PayloadType> mergedRelR;
    Relation<KeyType, PayloadType> mergedRelS;
    auto merging_t1 = high_resolution_clock::now();
    auto merging_t2 = high_resolution_clock::now();
    vector<uint32_t> curr_merging_timings_in_ms;
    vector<uint32_t> final_merging_timings_in_ms;
    for (int rp = 0; rp < RUN_NUMS; ++rp) 
    {
        BARRIER_ARRIVE(args->barrier, rv);

        // the first thread checkpoints the start time
        if(my_tid == 0){ 
            //gettimeofday(&args->sort_end_time, NULL);
            merging_t1 = high_resolution_clock::now();
        }

        join_steps.multiwaymerge_phase(numaregionid, partsR, partsS, args,
            &mergedRelR, &mergedRelS);

        BARRIER_ARRIVE(args->barrier, rv);

        // partition phase finished, thread-0 checkpoints the time
        if(my_tid == 0){
            //gettimeofday(&args->multiwaymerge_end_time, NULL);
            merging_t2 = high_resolution_clock::now();
            deltaT = std::chrono::duration_cast<std::chrono::microseconds>(merging_t2 - merging_t1).count();
            //deltaT = (args->multiwaymerge_end_time.tv_sec - args->sort_end_time.tv_sec) * 1000000 + args->multiwaymerge_end_time.tv_usec - args->sort_end_time.tv_usec;
            printf("---- %5s merging costs time (ms) = %10.4lf\n", "Non-learned sort join", deltaT * 1.0 / 1000);
            curr_merging_timings_in_ms.push_back((uint32_t)(deltaT * 1.0 / 1000)); //ms
        }
    
        if(!(rp == (RUN_NUMS - 1))){
            //TODO: make sure you can rerun here
            //if(args->tid == 0)
            //    destroy_hashtable(args->ht);

            //free_bucket_buffer(overflowbuf);

            BARRIER_ARRIVE(args->barrier, rv);
        }
    }

    if(args->tid == 0){
        std::sort(curr_merging_timings_in_ms.begin(), curr_merging_timings_in_ms.end());
        final_merging_timings_in_ms.push_back(curr_merging_timings_in_ms[(int)(curr_merging_timings_in_ms.size()/2)]);
    }

    /* the thread that allocated the merge buffer releases it. */
    if(is_first_thread_in_numa_region(my_tid)) {
        free(args->sharedmergebuffer[numaregionid]);
        //free_threadlocal(args->sharedmergebuffer[numaregionid],
        //MWAY_MERGE_BUFFER_SIZE);
    #ifdef SKEW_HANDLING
        /* initialize skew handling ; mwaytask_t taskqueues */
        taskqueue_free(args->numa_taskqueues[numaregionid]);
    #endif
    }

    BARRIER_ARRIVE(args->barrier, rv);

    /*************************************************************************
     *
     *   Phase.4) NUMA-local merge-join on local sorted runs.
     *
     *************************************************************************/
    auto join_t1 = high_resolution_clock::now();
    auto join_t2 = high_resolution_clock::now();
    vector<uint32_t> curr_join_timings_in_ms;
    vector<uint32_t> final_join_timings_in_ms;
    for (int rp = 0; rp < RUN_NUMS; ++rp) 
    {
        BARRIER_ARRIVE(args->barrier, rv);

        // the first thread checkpoints the start time
        if(my_tid == 0){ 
            //gettimeofday(&args->multiwaymerge_end_time, NULL);
            join_t1 = high_resolution_clock::now();
        }

        join_steps.mergejoin_phase(partsR, partsS, &mergedRelR, &mergedRelS, args);

        BARRIER_ARRIVE(args->barrier, rv);

        // partition phase finished, thread-0 checkpoints the time
        if(my_tid == 0){
            //gettimeofday(&args->mergejoin_end_time, NULL);
            join_t2 = high_resolution_clock::now();

            deltaT = std::chrono::duration_cast<std::chrono::microseconds>(join_t2 - join_t1).count();
            //deltaT = (args->mergejoin_end_time.tv_sec - args->multiwaymerge_end_time.tv_sec) * 1000000 + args->mergejoin_end_time.tv_usec - args->multiwaymerge_end_time.tv_usec;
            printf("---- %5s joining costs time (ms) = %10.4lf\n", "Non-learned sort join", deltaT * 1.0 / 1000);
            curr_join_timings_in_ms.push_back((uint32_t)(deltaT * 1.0 / 1000)); //ms
        }
    
        if(!(rp == (RUN_NUMS - 1))){
            //TODO: make sure you can rerun here
            //if(args->tid == 0)
            //    destroy_hashtable(args->ht);

            //free_bucket_buffer(overflowbuf);

            BARRIER_ARRIVE(args->barrier, rv);
        }

    }

    if(args->tid == 0){
        std::sort(curr_join_timings_in_ms.begin(), curr_join_timings_in_ms.end());
        final_join_timings_in_ms.push_back(curr_join_timings_in_ms[(int)(curr_join_timings_in_ms.size()/2)]);
    }

    if(args->tid == 0){
        std::vector<std::pair<std::string, std::vector<uint32_t>>> final_results =
         {{"Partition_in_ms", final_partition_timings_in_ms},
          {"Sort_in_ms", final_sorting_timings_in_ms},
          {"Merge_in_ms", final_merging_timings_in_ms},
          {"Join_in_ms", final_join_timings_in_ms}};
        write_csv(BENCHMARK_RESULTS_PATH, final_results);
    }

    /* clean-up */
    join_steps.partitioning_cleanup(partsR, partsS);

    free(args->threadrelchunks[my_tid]);
    /* clean-up temporary relations */
    if(args->nthreads > 1){
        free_threadlocal(mergedRelR.tuples, mergedRelR.num_tuples * sizeof(Tuple<KeyType, PayloadType>));
        free_threadlocal(mergedRelS.tuples, mergedRelS.num_tuples * sizeof(Tuple<KeyType, PayloadType>));
    }

    return 0;
}

int main(int argc, char **argv) 
{
    Relation<KeyType, PayloadType> rel_r;
    Relation<KeyType, PayloadType> rel_s;
    
    int64_t result = 0;
    uint64_t curr_num_tuples_r = RELATION_R_NUM_TUPLES;
    uint64_t curr_num_tuples_s = RELATION_S_NUM_TUPLES; 

#ifdef LOAD_RELATIONS_FOR_EVALUATION
    // loading pre-built datasets
    string curr_rel_r_folder_path = RELATION_R_FOLDER_PATH;
    string curr_rel_s_folder_path = RELATION_S_FOLDER_PATH;

    string curr_rel_r_file_name = RELATION_R_FILE_NAME;
    string curr_rel_s_file_name = RELATION_S_FILE_NAME;

    string curr_rel_r_file_extension = RELATION_R_FILE_EXTENSION;
    string curr_rel_s_file_extension = RELATION_S_FILE_EXTENSION;

    load_relation_threaded<KeyType, PayloadType>(&rel_r, RELATION_R_FILE_NUM_PARTITIONS, curr_rel_r_folder_path.c_str(), curr_rel_r_file_name.c_str(), curr_rel_r_file_extension.c_str(), curr_num_tuples_r);
    load_relation_threaded<KeyType, PayloadType>(&rel_s, RELATION_S_FILE_NUM_PARTITIONS, curr_rel_s_folder_path.c_str(), curr_rel_s_file_name.c_str(), curr_rel_s_file_extension.c_str(), curr_num_tuples_s);
#else

    string curr_rel_r_path = RELATION_R_PATH;
    string curr_rel_s_path = RELATION_S_PATH;

    string curr_rel_r_folder_path = RELATION_R_FOLDER_PATH;
    string curr_rel_s_folder_path = RELATION_S_FOLDER_PATH;

    string curr_rel_r_file_name = RELATION_R_FILE_NAME;
    string curr_rel_s_file_name = RELATION_S_FILE_NAME;

    string curr_rel_r_file_extension = RELATION_R_FILE_EXTENSION;
    string curr_rel_s_file_extension = RELATION_S_FILE_EXTENSION;

    // creating new datasets on-the-flay 
    result = create_eth_workload_relation_pk<KeyType, PayloadType>(&rel_r, curr_num_tuples_r, 0);
    //ASSERT_EQ(result, 0);
    #ifdef PERSIST_RELATIONS_FOR_EVALUATION
    write_relation_threaded<KeyType, PayloadType>(&rel_r, RELATION_R_FILE_NUM_PARTITIONS, curr_rel_r_folder_path.c_str(), curr_rel_r_file_name.c_str(), curr_rel_r_file_extension.c_str());
    write_relation<KeyType, PayloadType>(&rel_r, curr_rel_r_path.c_str());    
    #endif
    
    result = create_eth_workload_relation_pk<KeyType, PayloadType>(&rel_s, curr_num_tuples_s, 0);
    //ASSERT_EQ(result, 0);
    #ifdef PERSIST_RELATIONS_FOR_EVALUATION
    write_relation_threaded<KeyType, PayloadType>(&rel_s, RELATION_S_FILE_NUM_PARTITIONS, curr_rel_s_folder_path.c_str(), curr_rel_s_file_name.c_str(), curr_rel_s_file_extension.c_str());
    write_relation<KeyType, PayloadType>(&rel_s, curr_rel_s_path.c_str());    
    #endif
#endif

/*#ifdef LOAD_RELATIONS_FOR_EVALUATION
    // loading pre-built datasets
    string curr_rel_r_path = RELATION_R_PATH;
    string curr_rel_s_path = RELATION_S_PATH;

    load_relation<KeyType, PayloadType>(&rel_r, curr_rel_r_path.c_str(), curr_num_tuples_r);
    load_relation<KeyType, PayloadType>(&rel_s, curr_rel_s_path.c_str(), curr_num_tuples_s);    
#else
    // creating new datasets on-the-flay 
    result = create_eth_workload_relation_pk<KeyType, PayloadType>(&rel_r, curr_num_tuples_r, 0);
    //ASSERT_EQ(result, 0);
    #ifdef PERSIST_RELATIONS_FOR_EVALUATION
    write_relation<KeyType, PayloadType>(&rel_r, rel_r_path.c_str());
    #endif
    
    result = create_eth_workload_relation_pk<KeyType, PayloadType>(&rel_s, curr_num_tuples_s, 0);
    //ASSERT_EQ(result, 0);
    #ifdef PERSIST_RELATIONS_FOR_EVALUATION
    write_relation<KeyType, PayloadType>(&rel_s, rel_s_path.c_str());
    #endif
#endif*/


    int i, rv;
    pthread_barrier_t barrier;
    Result * joinresult;
    pthread_t tid[NUM_THREADS];
    pthread_attr_t attr;
    cpu_set_t set;
    
    joinresult = (Result *) malloc(sizeof(Result));
    joinresult->resultlist = (ThreadResult *) malloc(sizeof(ThreadResult) * NUM_THREADS);


    joincfg.NTHREADS = NUM_THREADS;
    joincfg.PARTFANOUT = ETH_SORT_MERGE_PARTITION_FANOUT_PER_THREAD;
    joincfg.SCALARMERGE = ETH_SORT_MERGE_IS_SCALAR_MERGE;
    joincfg.LEARNEDSORT = USE_LEARNED_SORT;
    joincfg.MWAYMERGEBUFFERSIZE = MWAY_MERGE_BUFFER_SIZE_DEFAULT;
    joincfg.NUMASTRATEGY = ETH_SORT_MERGE_NUMA_STRATEGY;
    numa_shuffle_init(joincfg.NUMASTRATEGY, joincfg.NTHREADS);

    /* check whether nr. of threads is a power of 2 */
    if((joincfg.NTHREADS & (joincfg.NTHREADS-1)) != 0)
    {
        perror("[ERROR] m-way sort-merge join runs with a power of 2 #threads.");
        exit(EXIT_FAILURE);
    }

    CACHELINEPADDING = ((ETH_SORT_MERGE_PARTITION_FANOUT_PER_THREAD) * CACHE_LINE_SIZE/sizeof(Tuple<KeyType, PayloadType>));
    RELATION_PADDING = ((NUM_THREADS) * CACHELINEPADDING * sizeof(Tuple<KeyType, PayloadType>));
    

    /**** allocate temporary space for partitioning ****/
    tmpRelpartR = NULL; tmpRelpartS = NULL;
    tmpRelpartR = (Tuple<KeyType, PayloadType>*) alloc_aligned(rel_r.num_tuples * sizeof(Tuple<KeyType, PayloadType>)
                                    + RELATION_PADDING);
    tmpRelpartS = (Tuple<KeyType, PayloadType>*) alloc_aligned(rel_s.num_tuples * sizeof(Tuple<KeyType, PayloadType>)
                                    + RELATION_PADDING);

    /* this is just to make sure that chunks of the temporary memory
    will be numa local to threads. */
    numa_localize<KeyType, PayloadType>(tmpRelpartR, rel_r.num_tuples, NUM_THREADS);
    numa_localize<KeyType, PayloadType>(tmpRelpartS, rel_s.num_tuples, NUM_THREADS);

    /**** allocate temporary space for sorting ****/
    tmpRelsortR = NULL; tmpRelsortS = NULL;
    tmpRelsortR = (Tuple<KeyType, PayloadType>*) alloc_aligned(rel_r.num_tuples * sizeof(Tuple<KeyType, PayloadType>)
                                    + RELATION_PADDING);
    tmpRelsortS = (Tuple<KeyType, PayloadType>*) alloc_aligned(rel_s.num_tuples * sizeof(Tuple<KeyType, PayloadType>)
                                    + RELATION_PADDING);

    /* this is just to make sure that chunks of the temporary memory
    will be numa local to threads. */
    numa_localize<KeyType, PayloadType>(tmpRelsortR, rel_r.num_tuples, NUM_THREADS);
    numa_localize<KeyType, PayloadType>(tmpRelsortS, rel_s.num_tuples, NUM_THREADS);

    threadrelchunks = (RelationPair<KeyType, PayloadType> **) malloc(NUM_THREADS * sizeof(RelationPair<KeyType, PayloadType>*));

    /* allocate histograms arrays, actual allocation is local to threads */
    histR = (uint32_t**) malloc(NUM_THREADS * sizeof(uint32_t*));


    ETHSortMergeMultiwayJoinThread<KeyType, PayloadType> args[NUM_THREADS];
    ETHSortMergeMultiwayJoinThread<KeyType, PayloadType> * args_ptr = args;

    rv = pthread_barrier_init(&barrier, NULL, NUM_THREADS);
    if(rv != 0){
        printf("[ERROR] Couldn't create the barrier\n");
        exit(EXIT_FAILURE);
    }

    pthread_attr_init(&attr);

    int  err = 0;
    size_t  stackSize = 0;

    /* Get the default value */
    err = pthread_attr_getstacksize(&attr, &stackSize);
    if (err) {
        perror("[ERROR] pthread stack size could not be get!");
        exit(0);
    }

    /* set the attribute with our required value */
    if (stackSize < REQUIRED_STACK_SIZE) {
        err = pthread_attr_setstacksize (&attr, REQUIRED_STACK_SIZE);
        if (err) {
            perror("[ERROR] pthread stack size could not be set!");
            exit(0);
        }
    }


    #ifdef DEVELOPMENT_MODE
    int num_numa = get_num_numa_regions_develop();
    #else
    int num_numa = get_num_numa_regions_v2();
    #endif  
    ptrs_to_sharedmergebufs = (Tuple<KeyType, PayloadType> **)
        alloc_aligned(num_numa*sizeof(Tuple<KeyType, PayloadType>*));

    #ifdef SKEW_HANDLING
    ptrs_to_taskqueues = (taskqueue_t **) alloc_aligned(num_numa*sizeof(taskqueue_t*));
    ptrs_to_taskqueues_locks = (pthread_mutex_t*) malloc(num_numa*sizeof(pthread_mutex_t));
    ptrs_to_is_numa_taskqueues_created = (int*) malloc(num_numa*sizeof(int));
    for(i = 0; i < num_numa; i++)
    {
        ptrs_to_is_numa_taskqueues_created[i] = 0;
        ptrs_to_taskqueues_locks[i] = PTHREAD_MUTEX_INITIALIZER; 
    }
    #endif


    initialize_non_learned_non_imv_sort_join_thread_args(&rel_r, &rel_s, &barrier, joinresult, args_ptr);

    for(i = 0; i < NUM_THREADS; i++){
        #ifdef DEVELOPMENT_MODE
        int cpu_idx = get_cpu_id_develop(i);
        #else
        int cpu_idx = get_cpu_id_v2(i);
        numa_thread_mark_active(cpu_idx);
        #endif

        CPU_ZERO(&set);
        CPU_SET(cpu_idx, &set);
        pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);

        rv = pthread_create(&tid[i], &attr, non_learned_non_imv_sort_join_thread, (void*)&args[i]);
        if (rv){
            printf("[ERROR] return code from pthread_create() is %d\n", rv);
            exit(-1);
        }

    }

    // wait for threads to finish
    for(i = 0; i < NUM_THREADS; i++){
        pthread_join(tid[i], NULL);
        result += args[i].result;
    }        
    joinresult->totalresults = result;
    joinresult->nthreads     = NUM_THREADS;

    printf("join results: %ld \n", result);

    return 0;
}

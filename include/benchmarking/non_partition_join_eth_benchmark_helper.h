#pragma once

#include "benchmarking/benchmark_helper_base.h"
#include "techniques/non_partition_join_eth_steps.h"

#include "configs/eth_configs.h"

template<typename KeyType, typename PayloadType, class TaskType, int NUM_THREADS = 2>
class ETHNonPartitionJoinBenchmarkHelper : public BenchmarkHelper<KeyType, PayloadType, TaskType, 
                                                      ETHNonPartitionJoinThread<KeyType, PayloadType, TaskType>,
                                                      EmptyParition<KeyType, PayloadType, TaskType, ETHNonPartitionJoinThread<KeyType, PayloadType, TaskType>>, 
                                                      NUM_THREADS> {
    
 public:

    static void * non_partition_join_thread_helper(void * args) {
        ETHNonPartitionJoinThreadWrapper * wrapped_args = (ETHNonPartitionJoinThreadWrapper *)args;
        return wrapped_args->instance->join_thread(wrapped_args->thread_args);
    }
    
    ETHNonPartitionJoinBenchmarkHelper(uint64_t hash_table_size_in = 1000) : BenchmarkHelper<KeyType, PayloadType, TaskType, 
                                                      ETHNonPartitionJoinThread<KeyType, PayloadType, TaskType>,
                                                      EmptyParition<KeyType, PayloadType, TaskType, ETHNonPartitionJoinThread<KeyType, PayloadType, TaskType>>, 
                                                      NUM_THREADS>()
    {
        hash_table_size = hash_table_size_in;
    }

    void initialize_join_thread_args(Relation<KeyType, PayloadType> * rel_r, 
                                          Relation<KeyType, PayloadType> * rel_s, 
                                          ETHNonPartitionJoinThread<KeyType, PayloadType, TaskType> * args){
        int i;
        int32_t numR, numS, numRthr, numSthr; /* total and per thread num */

        numR = rel_r->num_tuples;
        numS = rel_s->num_tuples;
        numRthr = numR / NUM_THREADS;
        numSthr = numS / NUM_THREADS;
#if INPUT_HASH_TABLE_SIZE       
        uint32_t nbuckets = hash_table_size;
#else
        uint32_t nbuckets = (rel_r->num_tuples / BUCKET_SIZE);
#endif        
        allocate_hashtable(&ht, nbuckets);

#ifdef DEVELOPMENT_MODE
        /*for(uint32_t j = 0; j < nbuckets; j++)
        {
            build_visits_map[j] = 0;
            probe_visits_map[j] = 0;            
        }*/        
#endif

        for(i = 0; i < NUM_THREADS; i++)
        {
            (*(args + i)).tid = i;
            (*(args + i)).ht = ht; 
#ifdef DEVELOPMENT_MODE
        /*    (*(args + i)).build_hash_bucket_visits = &build_visits_map;
            (*(args + i)).probe_hash_bucket_visits = &probe_visits_map;
            (*(args + i)).keys_hash_latch = &keys_hash_latch;            
            (*(args + i)).build_keys_list = &build_keys_list;
            (*(args + i)).build_keys_hash_list = &build_keys_hash_list;
            (*(args + i)).probe_keys_list = &probe_keys_list;
            (*(args + i)).probe_keys_hash_list = &probe_keys_hash_list;*/
#endif
            /* assing part of the relR for next thread */
            (*(args + i)).relR.num_tuples = (i == (NUM_THREADS-1)) ? numR : numRthr;
            (*(args + i)).relR.tuples = rel_r->tuples + numRthr * i;
            numR -= numRthr;

            /* assing part of the relS for next thread */
            (*(args + i)).relS.num_tuples = (i == (NUM_THREADS-1)) ? numS : numSthr;
            (*(args + i)).relS.tuples = rel_s->tuples + numSthr * i;
            numS -= numSthr;

            (*(args + i)).barrier = &barrier;
            (*(args + i)).threadresult  = &(joinresult->resultlist[i]);
        }
    }

    void * join_thread(void * param){
        ETHNonPartitionJoinThread<KeyType, PayloadType, TaskType> * args   = (ETHNonPartitionJoinThread<KeyType, PayloadType, TaskType> *) param;

        int rv;

        /* allocate overflow buffer for each thread */
        BucketBuffer<KeyType, PayloadType> * overflowbuf;
        init_bucket_buffer(&overflowbuf);

    #ifdef PERF_COUNTERS
        if(args->tid == 0){
            //TODO: performance counters to be implemented
        }
    #endif

        /* wait at a barrier until each thread starts and start timer */
        BARRIER_ARRIVE(args->barrier, rv);

        /* the first thread checkpoints the start time */
        if(args->tid == 0){
            /* no partitionig phase, but we use partition variables to store building stats */
            gettimeofday(&args->start_time, NULL);
        #ifndef DEVELOPMENT_MODE
            //args->e_start_to_partition.startCounters();
        #endif
        }


        /* insert tuples from the assigned part of relR to the ht */
        ETHNonPartitionJoinBuild<KeyType, PayloadType> build_data; 
        build_data.ht = args->ht;
        build_data.overflowbuf = &overflowbuf;
#ifdef DEVELOPMENT_MODE
        /*build_data.build_hash_bucket_visits = args->build_hash_bucket_visits;
        build_data.probe_hash_bucket_visits = args->probe_hash_bucket_visits; 
        build_data.keys_hash_latch = args->keys_hash_latch; 
        build_data.build_keys_list = args->build_keys_list;
        build_data.build_keys_hash_list = args->build_keys_hash_list;        
        build_data.probe_keys_list = args->probe_keys_list;
        build_data.probe_keys_hash_list = args->probe_keys_hash_list;*/        
#endif        

    #ifdef NPJ_ETH_AVX_IMV 
        join_steps.template build_rel_r_partition_imv<ETHNonPartitionJoinBuild<KeyType, PayloadType>>(&build_data, &args->relR, NULL);
    #else
        join_steps.template build_rel_r_partition<ETHNonPartitionJoinBuild<KeyType, PayloadType>>(&build_data, &args->relR, NULL);
    #endif            


        /* wait at a barrier until each thread completes build phase */
        BARRIER_ARRIVE(args->barrier, rv);

    #ifdef PERF_COUNTERS
        if(args->tid == 0)
        {
           //TODO: performance counters to be implemented
        }
        /* Just to make sure we get consistent performance numbers */
        BARRIER_ARRIVE(args->barrier, rv);
    #endif

        /* build phase finished, thread-0 checkpoints the time */
        if(args->tid == 0){
           gettimeofday(&args->partition_end_time, NULL);
        #ifndef DEVELOPMENT_MODE
           //args->e_start_to_partition.stopCounters();
           //args->e_partition_to_end.startCounters();
        #endif
        }

        /* probe for matching tuples from the assigned part of relS */  
    #ifdef NPJ_ETH_AVX_IMV      
        args->num_results = join_steps.template probe_rel_s_partition_imv<ETHNonPartitionJoinBuild<KeyType, PayloadType>>(NULL, &args->relS, &build_data, args->tid); 
    #else
        args->num_results = join_steps.template probe_rel_s_partition<ETHNonPartitionJoinBuild<KeyType, PayloadType>>(NULL, &args->relS, &build_data);
    #endif

        /* for a reliable timing we have to wait until all finishes */
        BARRIER_ARRIVE(args->barrier, rv);

        /* probe phase finished, thread-0 checkpoints the time */
        if(args->tid == 0){
            gettimeofday(&args->end_time, NULL);
        #ifndef DEVELOPMENT_MODE
            //args->e_partition_to_end.stopCounters();
        #endif
        }

    #ifdef PERF_COUNTERS
        if(args->tid == 0) 
        {
            //TODO: performance counters to be implemented
        }
        /* Just to make sure we get consistent performance numbers */
        BARRIER_ARRIVE(args->barrier, rv);
    #endif

        /* clean-up the overflow buffers */
        free_bucket_buffer(overflowbuf);

        return 0;
    }

    Result * run_join(Relation<KeyType, PayloadType> * rel_r, Relation<KeyType, PayloadType> * rel_s) 
    {
        int i, rv;
        pthread_t tid[NUM_THREADS];
        pthread_attr_t attr;
        cpu_set_t set;

        int64_t result = 0;
        joinresult = (Result *) malloc(sizeof(Result));

        ETHNonPartitionJoinThread<KeyType, PayloadType, TaskType> args[NUM_THREADS];
        ETHNonPartitionJoinThread<KeyType, PayloadType, TaskType> * args_ptr = args;

        rv = pthread_barrier_init(&barrier, NULL, NUM_THREADS);
        if(rv != 0){
            printf("[ERROR] Couldn't create the barrier\n");
            exit(EXIT_FAILURE);
        }

        pthread_attr_init(&attr);

        initialize_join_thread_args(rel_r, rel_s, args_ptr);

        for(i = 0; i < NUM_THREADS; i++){
            #ifdef DEVELOPMENT_MODE
            int cpu_idx = get_cpu_id_develop(i);
            #else
            int cpu_idx = get_cpu_id(i);
            #endif

            CPU_ZERO(&set);
            CPU_SET(cpu_idx, &set);
            pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);

            wrappers[i].instance = this;
            wrappers[i].thread_args = (void*)&args[i];

            rv = pthread_create(&tid[i], &attr, (THREADFUNCPTR) &non_partition_join_thread_helper, (void*)&wrappers[i]);
            if (rv){
                printf("[ERROR] return code from pthread_create() is %d\n", rv);
                exit(-1);
            }
        }

        // wait for threads to finish
        for(i = 0; i < NUM_THREADS; i++){
            pthread_join(tid[i], NULL);
            result += args[i].num_results;
        }        
        joinresult->totalresults = result;
        joinresult->nthreads     = NUM_THREADS;

        this->fill_timing_stats(&args[0].start_time, &args[0].partition_end_time, &args[0].end_time);
    #ifndef DEVELOPMENT_MODE
        //this->fill_perf_events_stats(&args[0].e_start_to_partition, &args[0].e_partition_to_end);            
    #endif
    
        return joinresult;
    }

    void clean_up(){
        BenchmarkHelper<KeyType, PayloadType, TaskType, 
                    ETHNonPartitionJoinThread<KeyType, PayloadType, TaskType>,
                    EmptyParition<KeyType, PayloadType, TaskType, ETHNonPartitionJoinThread<KeyType, PayloadType, TaskType>>, 
                    NUM_THREADS>::clean_up();

        destroy_hashtable(ht);
        free(joinresult);
    }

    void set_hash_table_size (uint64_t hash_table_size_input)
    {
        hash_table_size = hash_table_size_input;
    }
    
    ~ETHNonPartitionJoinBenchmarkHelper() {
    }

    struct ETHNonPartitionJoinThreadWrapper {
        ETHNonPartitionJoinBenchmarkHelper<KeyType, PayloadType, TaskType, NUM_THREADS> * instance;
        void * thread_args;
    };

#ifdef DEVELOPMENT_MODE
    /*unordered_map<uint64_t, uint64_t> build_visits_map;
    unordered_map<uint64_t, uint64_t> probe_visits_map;
    volatile char keys_hash_latch;
    vector<KeyType> build_keys_list;
    vector<uint64_t> build_keys_hash_list;
    vector<KeyType> probe_keys_list;
    vector<uint64_t> probe_keys_hash_list;*/ 
#endif

 protected:
    pthread_barrier_t barrier;
    Hashtable<KeyType, PayloadType> * ht;    
    Result * joinresult;

 private:
    ETHNonPartitionJoinSteps<KeyType, PayloadType, TaskType, 
                            ETHNonPartitionJoinThread<KeyType, PayloadType, TaskType>,
                            EmptyParition<KeyType, PayloadType, TaskType, ETHNonPartitionJoinThread<KeyType, PayloadType, TaskType>>> join_steps;
    ETHNonPartitionJoinThreadWrapper wrappers[NUM_THREADS];
    uint64_t hash_table_size;

};



#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "config.h"            /* autoconf header */
#include "configs/eth_configs.h"

#include "utils/eth_data_structures.h"
#include "utils/data_generation.h"
#include "utils/io.h"

#ifdef USE_ETH_NON_PARTITION_JOIN_GOOGLE_BENCHMARK
#include "benchmarking/non_partition_join_eth_benchmark_helper.h"
#endif

#define KeyType RELATION_KEY_TYPE
#define PayloadType RELATION_PAYLOAD_TYPE
#define TaskType Task<RELATION_KEY_TYPE, RELATION_PAYLOAD_TYPE>
#define NUM_THREADS NUM_THREADS_FOR_EVALUATION

int main(int argc, char **argv) 
{
    Relation<KeyType, PayloadType> rel_r;
    Relation<KeyType, PayloadType> rel_s;
    
    int result = 0;
    uint64_t curr_num_tuples_r = RELATION_R_NUM_TUPLES;
    uint64_t curr_num_tuples_s = RELATION_S_NUM_TUPLES; 

#ifdef LOAD_RELATIONS_FOR_EVALUATION
    // loading pre-built datasets
    string curr_rel_r_path = RELATION_R_PATH;
    string curr_rel_s_path = RELATION_S_PATH;

    result = load_relation<KeyType, PayloadType>(&rel_r, curr_rel_r_path.c_str(), curr_num_tuples_r);
    ASSERT_EQ(result, 0);
    result = load_relation<KeyType, PayloadType>(&rel_s, curr_rel_s_path.c_str(), curr_num_tuples_s);
    ASSERT_EQ(result, 0);
#else
    // creating new datasets on-the-flay 
    result = create_eth_workload_relation_pk<KeyType, PayloadType>(&rel_r, curr_num_tuples_r, 0);
    ASSERT_EQ(result, 0);
    #ifdef PERSIST_RELATIONS_FOR_EVALUATION
    write_relation<KeyType, PayloadType>(&rel_r, rel_r_path.c_str());
    #endif
    
    result = create_eth_workload_relation_pk<KeyType, PayloadType>(&rel_s, curr_num_tuples_s, 0);
    ASSERT_EQ(result, 0);
    #ifdef PERSIST_RELATIONS_FOR_EVALUATION
    write_relation<KeyType, PayloadType>(&rel_s, rel_s_path.c_str());
    #endif
#endif

  Result * results;
  TimingStats timingStats;
#ifndef DEVELOPMENT_MODE
  PerfEventStats perfEventStats;
#endif

#ifdef USE_ETH_NON_PARTITION_JOIN_GOOGLE_BENCHMARK

  ETHNonPartitionJoinBenchmarkHelper<KeyType, PayloadType, TaskType, NUM_THREADS> benchmark_helper;
  #ifdef INPUT_HASH_TABLE_SIZE 
  benchmark_helper.set_hash_table_size(1000001);
  #endif

  results = benchmark_helper.run_join(&rel_r, &rel_s);

  timingStats.total_partitioning_time_usec += benchmark_helper.timingStats->total_partitioning_time_usec;
  timingStats.total_joining_time_usec += benchmark_helper.timingStats->total_joining_time_usec;
  timingStats.total_algorithm_time_usec += benchmark_helper.timingStats->total_algorithm_time_usec;

  #ifndef DEVELOPMENT_MODE
  perfEventStats.total_partitioning_cycles += benchmark_helper.perfEventStats->total_partitioning_cycles;
  perfEventStats.total_partitioning_instructions += benchmark_helper.perfEventStats->total_partitioning_instructions;
  perfEventStats.total_partitioning_l1_misses += benchmark_helper.perfEventStats->total_partitioning_l1_misses;
  perfEventStats.total_partitioning_llc_misses += benchmark_helper.perfEventStats->total_partitioning_llc_misses;
  perfEventStats.total_partitioning_branch_misses += benchmark_helper.perfEventStats->total_partitioning_branch_misses;
  perfEventStats.total_partitioning_task_clock += benchmark_helper.perfEventStats->total_partitioning_task_clock;
  perfEventStats.total_partitioning_instructions_per_cycle += benchmark_helper.perfEventStats->total_partitioning_instructions_per_cycle;
  perfEventStats.total_partitioning_cpus += benchmark_helper.perfEventStats->total_partitioning_cpus;
  perfEventStats.total_partitioning_ghz += benchmark_helper.perfEventStats->total_partitioning_ghz;

  perfEventStats.total_joining_cycles += benchmark_helper.perfEventStats->total_joining_cycles;
  perfEventStats.total_joining_instructions += benchmark_helper.perfEventStats->total_joining_instructions;
  perfEventStats.total_joining_l1_misses += benchmark_helper.perfEventStats->total_joining_l1_misses;
  perfEventStats.total_joining_llc_misses += benchmark_helper.perfEventStats->total_joining_llc_misses;
  perfEventStats.total_joining_branch_misses += benchmark_helper.perfEventStats->total_joining_branch_misses;
  perfEventStats.total_joining_task_clock += benchmark_helper.perfEventStats->total_joining_task_clock;
  perfEventStats.total_joining_instructions_per_cycle += benchmark_helper.perfEventStats->total_joining_instructions_per_cycle;
  perfEventStats.total_joining_cpus += benchmark_helper.perfEventStats->total_joining_cpus;
  perfEventStats.total_joining_ghz += benchmark_helper.perfEventStats->total_joining_ghz;
  #endif

  printf("partitioning_time in usec: %ld, joining_time in usec: %ld, total_time in usec: %ld \n", timingStats.total_partitioning_time_usec, timingStats.total_joining_time_usec, timingStats.total_algorithm_time_usec);
#endif

  return 0;
}
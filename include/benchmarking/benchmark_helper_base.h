#pragma once


/* A tool to collect benchmarking statistics (NOT USED NOW) */

#include <sched.h>              /* CPU_ZERO, CPU_SET */
#include <sys/time.h>           /* gettimeofday */

#include "utils/data_structures.h"

#include "utils/base_utils.h"
#include "utils/barrier.h"
#include "utils/rdtsc.h"
#include "utils/affinity.h"
#include "utils/memory.h"
#include "utils/io.h"
#include "utils/perf_event.h"
#include "utils/cpu_mapping.h"
#include "utils/cpu_mapping_one_numa.h"

#include "configs/base_configs.h"


template<typename KeyType, typename PayloadType, typename TaskType, typename JoinThreadType, typename PartitionType, int NUM_THREADS = 2>
class BenchmarkHelper {
 public:
 
  BenchmarkHelper()
  {
    timingStats = (TimingStats *) malloc(sizeof(TimingStats));
  #ifndef DEVELOPMENT_MODE
    //perfEventStats = (PerfEventStats *) malloc(sizeof(PerfEventStats));
  #endif
  }
  
  void initialize_join_thread_args(Relation<KeyType, PayloadType> * rel_r, 
                                        Relation<KeyType, PayloadType> * rel_s, 
                                        JoinThreadType * args){}

  void * join_thread(void * param)
  {
    return 0;
  }

  Result * run_join(Relation<KeyType, PayloadType> * rel_r, Relation<KeyType, PayloadType> * rel_s) 
  {
    return 0;
  }

  void fill_timing_stats(struct timeval * start_time, struct timeval * partition_end_time, struct timeval * end_time)
  {
    timingStats->total_partitioning_time_usec = (((*partition_end_time).tv_sec*1000000L + (*partition_end_time).tv_usec)
                        - ((*start_time).tv_sec*1000000L+(*start_time).tv_usec));
    timingStats->total_joining_time_usec = (((*end_time).tv_sec*1000000L + (*end_time).tv_usec)
                        - ((*partition_end_time).tv_sec*1000000L+(*partition_end_time).tv_usec));
    timingStats->total_algorithm_time_usec = (((*end_time).tv_sec*1000000L + (*end_time).tv_usec)
                        - ((*start_time).tv_sec*1000000L+(*start_time).tv_usec));
  }

#ifndef DEVELOPMENT_MODE
  /*void fill_perf_events_stats(PerfEvent * start_to_partition_event, PerfEvent * partition_to_end_event)
  {
    unsigned i;

    // Partitioning PerfEvent stats
    for (i = 0; i < start_to_partition_event->events.size(); i++)
    {
      if(start_to_partition_event->names[i] == "cycles")
        perfEventStats->total_partitioning_cycles = start_to_partition_event->events[i].readCounter();
      else if (start_to_partition_event->names[i] == "instructions")
        perfEventStats->total_partitioning_instructions = start_to_partition_event->events[i].readCounter();
      else if (start_to_partition_event->names[i] == "L1-misses")
        perfEventStats->total_partitioning_l1_misses = start_to_partition_event->events[i].readCounter();
      else if (start_to_partition_event->names[i] == "LLC-misses")
        perfEventStats->total_partitioning_llc_misses = start_to_partition_event->events[i].readCounter();
      else if (start_to_partition_event->names[i] == "branch-misses")
        perfEventStats->total_partitioning_branch_misses = start_to_partition_event->events[i].readCounter();
      else if (start_to_partition_event->names[i] == "task-clock")
        perfEventStats->total_partitioning_task_clock = start_to_partition_event->events[i].readCounter();
    }

    perfEventStats->total_partitioning_instructions_per_cycle = start_to_partition_event->getIPC();
    perfEventStats->total_partitioning_cpus = start_to_partition_event->getCPUs();
    perfEventStats->total_partitioning_ghz = start_to_partition_event->getGHz();

    // Joining PerfEvent stats
    for (i = 0; i < partition_to_end_event->events.size(); i++)
    {
      if(partition_to_end_event->names[i] == "cycles")
        perfEventStats->total_joining_cycles = partition_to_end_event->events[i].readCounter();
      else if (partition_to_end_event->names[i] == "instructions")
        perfEventStats->total_joining_instructions = partition_to_end_event->events[i].readCounter();
      else if (partition_to_end_event->names[i] == "L1-misses")
        perfEventStats->total_joining_l1_misses = partition_to_end_event->events[i].readCounter();
      else if (partition_to_end_event->names[i] == "LLC-misses")
        perfEventStats->total_joining_llc_misses = partition_to_end_event->events[i].readCounter();
      else if (partition_to_end_event->names[i] == "branch-misses")
        perfEventStats->total_joining_branch_misses = partition_to_end_event->events[i].readCounter();
      else if (partition_to_end_event->names[i] == "task-clock")
        perfEventStats->total_joining_task_clock = partition_to_end_event->events[i].readCounter();
    }

    perfEventStats->total_joining_instructions_per_cycle = partition_to_end_event->getIPC();
    perfEventStats->total_joining_cpus = partition_to_end_event->getCPUs();
    perfEventStats->total_joining_ghz = partition_to_end_event->getGHz();
  }*/
#endif

  void clean_up()
  {
    free(timingStats);
  #ifndef DEVELOPMENT_MODE
    //free(perfEventStats);
  #endif
  }

  TimingStats * timingStats;
#ifndef DEVELOPMENT_MODE
  //PerfEventStats * perfEventStats;
#endif
};

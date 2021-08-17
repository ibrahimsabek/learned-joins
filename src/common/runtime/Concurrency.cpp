/* An adapted implementation of concurrency for benchmarking with TPCH and SSB join queries based on https://github.com/fzhedu/db-imv */

#include "common/runtime/Concurrency.hpp"

#include <sys/time.h>
namespace runtime {

thread_local Worker* this_worker;
thread_local bool currentBarrier = false;

WorkerGroup mainGroup(1);
HierarchicBarrier mainBarrier(1, nullptr);
GlobalPool defaultPool;
Worker mainWorker(&mainGroup, &mainBarrier, defaultPool);

double gettime() {
   struct timeval now_tv;
   gettimeofday(&now_tv, NULL);
   return ((double)now_tv.tv_sec) + ((double)now_tv.tv_usec) / 1000000.0;
}
void Worker::log_time(std::string pipeline_name){
  if(group->main_worker==this_worker){
    group->pipeline_cost_time.push_back(make_pair(pipeline_name,gettime()));
  }
}

} // namespace runtime

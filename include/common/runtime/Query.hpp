#pragma once

/* An adapted implementation of Query class for benchmarking with TPCH and SSB join queries based on https://github.com/fzhedu/db-imv */

#include "common/runtime/MemoryPool.hpp"
#include <memory>

namespace runtime {

class Query {
 public:
   GlobalPool pool;
   std::unique_ptr<BlockRelation> result;
   Query() { result = std::make_unique<BlockRelation>(); }
   GlobalPool* participate() { return this_worker->allocator.setSource(&pool); }
   void leave(GlobalPool* prev) { this_worker->allocator.setSource(prev); }
};

} // namespace runtime

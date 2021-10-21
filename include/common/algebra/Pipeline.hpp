#pragma once

/* An adapted implementation of query pipelineing used for benchmarking with TPCH and SSB join queries based on https://github.com/fzhedu/db-imv */

#include "common/algebra/Operators.hpp"
#include <ostream>


namespace algebra{
   class Pipeline {
   public:
      std::unique_ptr<Operator> root;
   public:
      Pipeline& map();
      Pipeline& scan(runtime::Relation& r);
      Pipeline& print(std::vector<std::string> attrs);
      Pipeline& hashjoin(Pipeline& materialized);
      void print(std::ostream& s) const;
   };
}

#pragma once
#include "vectorwise/Operators.hpp"

/* An interface to implement queries for benchmarking with TPCH and SSB join queries based on https://github.com/fzhedu/db-imv */

class Query {
 protected:
   vectorwise::SharedStateManager shared;

 public:
   /// To be called by every worker thread to set up thread local resources
   void setup();
   /// To be called by every worker thread in order to execute query
   void run();
};

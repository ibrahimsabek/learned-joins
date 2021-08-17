#pragma once

/* An interface to import data for benchmarking with TPCH and SSB join queries based on https://github.com/fzhedu/db-imv */

#include "Database.hpp"
#include <string>

namespace runtime {
   /// imports tpch relations from CSVs in dir into db
   void importTPCH(std::string dir, Database& db);

   /// imports star schema benchmark from CSVs in dir into db
   void importSSB(std::string dir, Database& db);

   void importSSB_modified(std::string dir, Database& db);
}

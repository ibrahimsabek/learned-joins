#pragma once

#include "head.hpp"
#include "imv/HashBuild.hpp"
size_t amac_build_q11_date(size_t begin, size_t end, Database& db, runtime::Hashmap* hash_table, Allocator*allo, int entry_size, uint64_t* pos_buf);
size_t imv_build_q11_date(size_t begin, size_t end, Database& db, runtime::Hashmap* hash_table, Allocator*allo, int entry_size, uint64_t* pos_buf);

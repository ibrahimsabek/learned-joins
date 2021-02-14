#pragma once
#include "head.hpp"
#include "imv/HashProbe.hpp"
#include "common/runtime/Database.hpp"

size_t simd_filter_qa_build(size_t& begin, size_t end, Database& db, uint64_t* pos_buff, types::Date constrant);
size_t build_gp_qa(size_t begin, size_t end, Database& db, runtime::Hashmap* hash_table, Allocator*allo, int entry_size, uint64_t* pos_buff);
size_t build_imv_qa(size_t begin, size_t end, Database& db, runtime::Hashmap* hash_table, Allocator*allo, int entry_size, uint64_t* pos_buff);
size_t build_pipeline_imv_qa(size_t begin, size_t end, Database& db, runtime::Hashmap* hash_table, Allocator*allo, int entry_size, types::Date constrant);

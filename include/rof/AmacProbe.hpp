#pragma once

#include "head.hpp"
size_t amac_probe_q11(size_t begin, size_t end, Database& db, runtime::Hashmap* hash_table,uint64_t & results, uint64_t* pos_buff);

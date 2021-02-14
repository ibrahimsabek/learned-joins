#pragma once

#include "imv/HashBuild.hpp"
size_t build_imv_q1x(size_t begin, size_t end, Database& db, runtime::Hashmap* hash_table, Allocator*allo, int entry_size);

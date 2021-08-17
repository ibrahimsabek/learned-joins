#pragma once

/* An interface of implementing Q1 in TPCH benchmark based on https://github.com/fzhedu/db-imv */


#include "head.hpp"
#define PARTITION_SIZE 4096
using namespace types;

size_t agg_raw_q1(size_t begin, size_t end, Database& db,
                  Hashmapx<types::Integer, tuple<Numeric<12, 2>, Numeric<12, 2>, Numeric<12, 4>, Numeric<12, 6>, int64_t>, hashFun, false>* hash_table,
                  PartitionedDeque<PARTITION_SIZE>* partition, void** entry_addrs, void** results_entry);
size_t agg_gp_q1(size_t begin, size_t end, Database& db,
                  Hashmapx<types::Integer, tuple<Numeric<12, 2>, Numeric<12, 2>, Numeric<12, 4>, Numeric<12, 6>, int64_t>, hashFun, false>* hash_table,
                  PartitionedDeque<PARTITION_SIZE>* partition, void** entry_addrs, void** results_entry);
size_t agg_amac_q1(size_t begin, size_t end, Database& db,
                  Hashmapx<types::Integer, tuple<Numeric<12, 2>, Numeric<12, 2>, Numeric<12, 4>, Numeric<12, 6>, int64_t>, hashFun, false>* hash_table,
                  PartitionedDeque<PARTITION_SIZE>* partition, void** entry_addrs, void** results_entry);
size_t agg_imv_q1(size_t begin, size_t end, Database& db,
                  Hashmapx<types::Integer, tuple<Numeric<12, 2>, Numeric<12, 2>, Numeric<12, 4>, Numeric<12, 6>, int64_t>, hashFun, false>* hash_table,
                  PartitionedDeque<PARTITION_SIZE>* partition, void** entry_addrs, void** results_entry);
size_t agg_simd_q1(size_t begin, size_t end, Database& db,
                  Hashmapx<types::Integer, tuple<Numeric<12, 2>, Numeric<12, 2>, Numeric<12, 4>, Numeric<12, 6>, int64_t>, hashFun, false>* hash_table,
                  PartitionedDeque<PARTITION_SIZE>* partition, void** entry_addrs, void** results_entry);

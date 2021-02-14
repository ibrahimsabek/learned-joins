#pragma once
#include <iostream>
#include "head.hpp"
size_t simd_filter_q11_build_date(size_t& begin, size_t end, Database& db, uint64_t* pos_buff);
size_t simd_filter_q11_probe(size_t& begin, size_t end, Database& db, uint64_t* pos_buff);

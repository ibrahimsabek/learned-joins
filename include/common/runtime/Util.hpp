#pragma once

/* An adapted implementation of utilities for benchmarking with TPCH and SSB join queries based on https://github.com/fzhedu/db-imv */

#include "stddef.h"
#include "stdint.h"

#define CACHELINE_SIZE 64

template <typename T>
T* addBytes(T* t, size_t bytes)
/// Increments t by `bytes `bytes
{
   return reinterpret_cast<T*>(reinterpret_cast<uint8_t*>(t) + bytes);
}


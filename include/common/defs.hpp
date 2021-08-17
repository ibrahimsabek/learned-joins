/* An adapted implementation of some definitions for benchmarking with TPCH and SSB join queries based on https://github.com/fzhedu/db-imv */

#include <cstdint>

namespace defs{

#if HASH_SIZE == 32
  using hash_t = uint32_t;
#else
  using hash_t = uint64_t;
#endif

};

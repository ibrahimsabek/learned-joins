#pragma once

#include <cassert>
#include <vector>

#include "../convenience/builtins.hpp"
#include "bitvector.hpp"

namespace exotic_hashing::support {
   template<class T, class BitStream = FixedBitvector<sizeof(T) * 8, T>>
   struct FixedBitConverter {
      forceinline BitStream operator()(const T& data) const {
         BitStream bs(data);
         assert(bs.size() == sizeof(T) * 8);
#ifndef NDEBUG
         T reconstructed = 0;
         for (size_t i = 0; i < bs.size(); i++) {
            const uint64_t bit = bs[i];
            reconstructed |= bit << (sizeof(T) * 8 - i - 1);
         }
         assert(reconstructed == data);
#endif
         return bs;
      }
   };
} // namespace exotic_hashing::support

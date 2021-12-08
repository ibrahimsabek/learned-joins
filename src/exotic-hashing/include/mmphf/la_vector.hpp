#pragma once

#include <la_vector.hpp>

#include "include/convenience/builtins.hpp"

namespace exotic_hashing {
   template<class Data>
   class LAVector {
      la_vector<Data> vec{};

     public:
      LAVector() noexcept = default;

      template<class RandomIt>
      LAVector(const RandomIt& begin, const RandomIt& end) {
         vec = decltype(vec)(begin, end);
      }

      explicit LAVector(std::vector<Data> dataset) {
         // TODO(dominik): read paper to find out if data really has to be presorted!
         std::sort(dataset.begin(), dataset.end());
         vec = decltype(vec)(dataset);
      }

      forceinline size_t operator()(const Data& key) const {
         return vec.rank(key);
      }

      static std::string name() {
         return "LAVector";
      }

      size_t byte_size() const {
         return vec.size_in_bytes();
      }
   };
} // namespace exotic_hashing

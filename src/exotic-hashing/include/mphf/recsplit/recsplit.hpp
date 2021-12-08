/*
 * This is a lightweight wrapper for recsplit, reexporting symbols
 * from suxcpp (https://github.com/vigna/sux.git) to make them:
 * 1) usable within this codebase
 * 2) lightweight in the sense that most unused functionality was plainly removed
 *
 * Please see LICENSE.md for more information
 */

#include "../../convenience/builtins.hpp"
#include "src/RecSplit.hpp"

#include <string>
#include <vector>

namespace exotic_hashing {
   template<class Data, size_t BucketSize = 9, size_t LeafSize = 12,
            sux::util::AllocType AllocType = sux::util::AllocType::MALLOC>
   struct RecSplit : protected sux::function::RecSplit<LeafSize, AllocType> {
      explicit RecSplit(const std::vector<Data>& d)
         : sux::function::RecSplit<LeafSize, AllocType>(to_string(d), BucketSize) {}

      static std::string name() {
         return "RecSplit_leaf" + std::to_string(LeafSize) + "_bucket" + std::to_string(BucketSize);
      }

      forceinline size_t operator()(const Data& key) const {
         return sux::function::RecSplit<LeafSize, AllocType>::operator()(reinterpret_to_string(key));
      }

      size_t byte_size() const {
         return this->_totalBitSize / 8;
      };

     private:
      template<class T>
      static std::vector<std::string> to_string(const std::vector<T>& v) {
         std::vector<std::string> res;
         for (const auto k : v)
            res.emplace_back(reinterpret_to_string(k));

         return res;
      }

      template<class T>
      static std::string reinterpret_to_string(const T& v) {
         return std::string(reinterpret_cast<const char*>(&v), sizeof(T));
      }
   };
}; // namespace exotic_hashing

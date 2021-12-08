#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <fst.hpp>

#include "../convenience/builtins.hpp"

namespace exotic_hashing {
   /**
    * Lightweight wrapper around FST:
    * https://github.com/christophanneser/FST
    */
   template<class Key>
   class FastSuccinctTrie {
      static forceinline std::string convert(const Key& key) {
         Key endian_swapped_word = 0;
         switch (sizeof(Key)) {
            case 8:
               endian_swapped_word = __builtin_bswap64(key);
               break;
            case 4:
               endian_swapped_word = __builtin_bswap32(key);
               break;
            default:
               break;
         }
         return std::string(reinterpret_cast<const char*>(&endian_swapped_word), sizeof(Key));
      }

      std::unique_ptr<fst::FST> _fst{};
      Key _min_key, _max_key;

      // FST internally stores a pointer to the string keys
      // which is required during lookup. Therefore, we must
      // keep the string keys in memory!
      std::vector<std::string> converted_keys;

     public:
      FastSuccinctTrie() noexcept = default;

      /**
       * Builds a new fast succinct trie from a dataset
       */
      explicit FastSuccinctTrie(std::vector<Key> dataset) {
         std::sort(dataset.begin(), dataset.end());
         construct(dataset.begin(), dataset.end());
      }

      /**
       * Constructs on already sorted range of keys
       */
      template<class ForwardIt>
      FastSuccinctTrie(const ForwardIt& begin, const ForwardIt& end) {
         construct(begin, end);
      }

      /**
       * Lazily constructs on *already sorted* range of keys
       */
      template<class RandomIt>
      void construct(const RandomIt& begin, const RandomIt& end) {
         // find min and max elements. Assume dataset is not sorted!
         const auto min_it = std::min_element(begin, end);
         if (min_it < end)
            _min_key = *min_it;
         const auto max_it = std::max_element(begin, end);
         if (max_it < end)
            _max_key = *max_it;

         // convert keys to string. IMPORTANT: they must be sorted for fst to work
         converted_keys.reserve(std::distance(begin, end));
         for (auto it = begin; it < end; it++)
            converted_keys.emplace_back(convert(*it));

         // construct fst. Note that values is never used during construction, hence it is simply left empty
         _fst = std::make_unique<fst::FST>(converted_keys);
      }

      forceinline size_t operator()(const Key& key) const {
         const auto str_key = convert(key);

         std::uint64_t rank = 0;
         _fst->lookupKey(str_key, rank);
         return rank;
      }

      static std::string name() {
         return "FastSuccinctTrie";
      }

      size_t byte_size() const {
         size_t base_size = sizeof(FastSuccinctTrie<Key>) + _fst->getMemoryUsage();

         // ignore converted keys byte usage, since this additional overhead does not exist
         // in real world implementations.
         return base_size - sizeof(std::vector<std::string>);
      }
   };
}; // namespace exotic_hashing

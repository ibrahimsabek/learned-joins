#pragma once

#include <algorithm>
#include <string>
#include <vector>

#include "include/convenience/builtins.hpp"
#include "include/support/support.hpp"

namespace exotic_hashing {
   /**
    * Most basic minimal perfect hash function, i.e.,
    * mapping keys to their rank within the keyset.
    *
    * Meant as a space & speed baseline, i.e., slower than
    * this is only justified if space is also smaller.
    *
    * Using more space than this function is not desirable.
    */
   template<class Data>
   class RankHash {
      std::vector<Data> dataset;

      /// assumes sorted *full* data is in dataset
      void construct() {
         // omit every second element, deleting junk and ensuring the final dataset
         // vector is minimal, i.e., does not waste any additional space
         for (size_t i = 1, j = 2; j < dataset.size(); i++, j += 2)
            dataset[i] = dataset[j];
         const size_t middle = (dataset.size() / 2) + (dataset.size() & 0x1);
         dataset.erase(dataset.begin() + middle, dataset.end());
         dataset.resize(dataset.size());
      }

     public:
      RankHash() noexcept = default;

      /**
       * Constructs on already sorted range of keys
       */
      template<class ForwardIt>
      RankHash(const ForwardIt& begin, const ForwardIt& end) {
         construct(begin, end);
      }

      /**
       * Constructs on arbitrarily ordered keyset
       */
      explicit RankHash(const std::vector<Data>& d) : dataset(d) {
         // sort the dataset
         std::sort(dataset.begin(), dataset.end());

         // construct on internal dataset
         construct();
      }

      /**
       * Constructs on *already sorted* range of keys
       */
      template<class RandomIt>
      void construct(const RandomIt& begin, const RandomIt& end) {
         dataset = decltype(dataset)(begin, end);
         construct();
      }

      static std::string name() {
         return "RankHash";
      }

      forceinline size_t operator()(const Data& key) const {
         // primitively compute rank of key by:
         // 1. binary searching it in the sorted, compressed dataset
         const auto iter = std::lower_bound(dataset.begin(), dataset.end(), key);

         // 2. computing its rank based on the compressed keyset. If key does
         //    not match iter assume it was removed during compression and its
         //    index therefore is 2*iter_pos-1
         const size_t iter_pos = std::distance(dataset.begin(), iter);

         if (unlikely(iter == dataset.end()))
            return 2 * iter_pos - 1;
         return 2 * iter_pos - (*iter == key ? 0 : 1);
      }

      size_t byte_size() const {
         return dataset.size() * sizeof(Data) + sizeof(std::vector<Data>);
      };
   };

   /**
    * Most basic minimal perfect hash function, i.e.,
    * mapping keys to their rank within the keyset.
    *
    * Storage optimization compared to regular RankHash
    *
    * Meant as a space & speed baseline, i.e., slower than
    * this is only justified if space is also smaller.
    *
    * Using more space than this function is not desirable.
    */
   template<class Data>
   class CompressedRankHash {
      support::EliasFanoList<Data> efl{};

      /// assumes sorted *full* data is in dataset & requires capability to make modifications
      void construct(std::vector<Data>& dataset) {
         // omit every second element, deleting junk and ensuring the final dataset
         // vector is minimal, i.e., does not waste any additional space
         for (size_t i = 1, j = 2; j < dataset.size(); i++, j += 2)
            dataset[i] = dataset[j];
         const size_t middle = (dataset.size() / 2) + (dataset.size() & 0x1);
         dataset.erase(dataset.begin() + middle, dataset.end());
         dataset.resize(dataset.size());

         // store in compressed form
         efl = decltype(efl)(dataset.begin(), dataset.end());
      }

     public:
      CompressedRankHash() noexcept = default;

      /**
       * Constructs on *already sorted* range of keys
       */
      template<class ForwardIt>
      CompressedRankHash(const ForwardIt& begin, const ForwardIt& end) {
         construct(begin, end);
      }

      /**
       * Constructs on arbitrarily ordered keyset
       */
      explicit CompressedRankHash(std::vector<Data> dataset) {
         // sort the dataset
         std::sort(dataset.begin(), dataset.end());

         // construct on internal dataset
         construct(dataset);
      }

      /**
       * Constructs on *already sorted* range of keys
       */
      template<class RandomIt>
      void construct(const RandomIt& begin, const RandomIt& end) {
         std::vector<Data> dataset(begin, end);
         construct(dataset);
      }

      static std::string name() {
         return "CompressedRankHash";
      }

      constexpr forceinline size_t operator()(const Data& key) const {
         // compute rank of key by binary searching in sorted dataset
         const auto index = support::lower_bound(0, efl.size(), key, efl);

         if (index == efl.size())
            return 2 * index - 1;

         // account for sorted dataset omitting every second key
         return 2 * index - (efl[index] == key ? 0 : 1);
      }

      size_t byte_size() const {
         return efl.byte_size();
      };
   };
} // namespace exotic_hashing

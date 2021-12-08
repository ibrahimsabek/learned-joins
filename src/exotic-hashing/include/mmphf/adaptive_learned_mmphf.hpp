#pragma once

#include <algorithm>
#include <iostream>
#include <memory>
#include <sdsl/vectors.hpp>
#include <string>
#include <vector>

#include "../omphf/mwhc.hpp"
#include "../support/bitconverter.hpp"
#include "../support/clustering.hpp"
#include "../support/elias_fano_list.hpp"
#include "../support/support.hpp"
#include "hollow_trie.hpp"
#include "learned_linear.hpp"

// Ordering important
#include "../convenience/builtins.hpp"

namespace exotic_hashing {
   // turns out for compressed mwhc 832, 1665 are optimal sizes to obtain min bits/key.
   template<class Data, size_t MinRegionSize = 1665, size_t MaxRegionCount = 10000>
   class AdaptiveLearnedMMPHF {
      /**
       * Certain model types only work on specific dataset regions. To
       * compensate for this, dynamically choose the correct building block for
       * the job
       */
      struct BuildingBlock {
         enum Type
         {
            MWHC = 0x0,
            LearnedLinear = 0x1
         };

         template<class... ConstructorArgs>
         explicit BuildingBlock(Type type, ConstructorArgs... args) : type(type) {
            switch (type) {
               case MWHC: {
                  const auto* ptr = new exotic_hashing::CompressedMWHC<Data>(args...);
                  m = (decltype(m)) (ptr);
                  break;
               }
               case LearnedLinear: {
                  const auto* ptr = new exotic_hashing::LearnedLinear<Data>(args...);
                  m = (decltype(m)) (ptr);
                  break;
               }
            }
         }

         forceinline size_t operator()(const Data& key) const {
            switch (type) {
               case MWHC: {
                  const auto* ptr = reinterpret_cast<exotic_hashing::CompressedMWHC<Data>*>(m);
                  return ptr->operator()(key);
               }
               case LearnedLinear: {
                  const auto ptr = reinterpret_cast<exotic_hashing::LearnedLinear<Data>*>(m);
                  return ptr->operator()(key);
               }
            }
         }

         size_t byte_size() const {
            switch (type) {
               case MWHC: {
                  const auto* ptr = reinterpret_cast<exotic_hashing::CompressedMWHC<Data>*>(m);
                  return ptr->byte_size();
               }
               case LearnedLinear: {
                  const auto ptr = reinterpret_cast<exotic_hashing::LearnedLinear<Data>*>(m);
                  return ptr->byte_size();
               }
            }
         }

         // Can't copy BuildingBlock (for now)
         BuildingBlock(const BuildingBlock& other) = delete;
         BuildingBlock& operator=(const BuildingBlock& other) = delete;

         /// Custom copy constructor is necessary since sdsl's rank support contains a pointer to bitvec
         BuildingBlock(BuildingBlock&& other) noexcept {
            type = other.type;
            m = other.m;
            other.m = 0x0;
         }

         /// Custom copy constructor is necessary since sdsl's select support contains a pointer to upper
         BuildingBlock& operator=(BuildingBlock&& other) noexcept {
            type = other.type;
            m = other.m;
            other.m = 0x0;

            return *this;
         }

         ~BuildingBlock() noexcept {
            if (m != 0x0) {
               switch (type) {
                  case MWHC:
                     delete reinterpret_cast<exotic_hashing::CompressedMWHC<Data>*>(m);
                     break;
                  case LearnedLinear:
                     delete reinterpret_cast<exotic_hashing::LearnedLinear<Data>*>(m);
                     break;
               }
            }

            m = 0x0;
         }

        private:
         Type type : 1;
         std::uint64_t m : sizeof(void*) * 8 - 1;
      };

      /// root sorts incoming keys into one of the ranges that is then served by a leaf model
      support::EliasFanoList<Data> delimiters{};

      /// each leaf is a mmphf over a region (e.g., a dense cluster) of the dataset
      std::vector<BuildingBlock> leafs{};

      /// to obtain a global mmphf over the entire data, it is necessary to
      /// store the min rank of each region.
      support::EliasFanoList<size_t> region_offsets{};

      /**
       * Computes the required density threshold to match a certain bits per
       * key for LearnedLinear models
       */
      constexpr forceinline double learned_linear_density_threshold(double bits_per_key) {
         // bits_per_key >= (8 * sizeof(LearnedLinear) / |X|) + 1/density(X)
         // -> neglecting the first term (assumption: large |X|), we obtain
         // density(X) >= 1/bits_per_key
         return 1.0 / bits_per_key;
      }

     public:
      explicit AdaptiveLearnedMMPHF(std::vector<Data> data) {
         // ensure data is sorted before we start
         std::sort(data.begin(), data.end());

         // thresholds for clustering
         const auto density_threshold = learned_linear_density_threshold(5);
         const auto size_threshold = std::max(MinRegionSize, data.size() / MaxRegionCount);

         // 1. cluster based on density metric
         const auto clusters = support::cluster(data.begin(), data.end(), density_threshold);
         size_t clusters_count = clusters.size() - 1; // due to representation as iterators

         // 2. build one leaf model for each region
         std::vector<size_t> rs{0};
         std::vector<Data> d;
         for (size_t i = 1; i < clusters.size(); i++) {
            const auto a = clusters[i - 1];
            auto b = clusters[i];

            // merge small regions until region threshold is
            // reached or we're at the end of the dataset
            size_t region_size = std::distance(a, b);
            bool size_violation = false;
            while (region_size < size_threshold && i + 1 < clusters.size()) {
               size_violation = true;
               clusters_count--;
               b = clusters[++i];
               region_size = std::distance(a, b);
            }
            // ensure regions are now properly sized
            assert(region_size == std::distance(a, b));
            assert(i + 1 >= clusters.size() || region_size >= size_threshold);

            if (i + 1 < clusters.size())
               rs.emplace_back(region_size + rs.back());

            if (rs.size() > 2)
               d.emplace_back(*a);

            // choose model type based on size violation for now
            leafs.emplace_back(size_violation ? BuildingBlock::Type::MWHC : BuildingBlock::Type::LearnedLinear, a, b);
         }
         region_offsets = decltype(region_offsets)(rs.begin(), rs.end());
         delimiters = decltype(delimiters)(d.begin(), d.end());
      }

      forceinline size_t operator()(const Data& key) const {
         const auto lb = support::lower_bound(0, delimiters.size(), key, delimiters);
         const size_t region_ind = lb + ((lb < delimiters.size() && delimiters[lb] == key) & 0x1);

         const size_t offset = region_offsets[region_ind];
         const size_t local_rank = leafs[region_ind](key);

         return local_rank + offset;
      }

      forceinline size_t byte_size() const {
         size_t res = region_offsets.byte_size() + delimiters.byte_size();
         for (const auto& leaf : leafs)
            res += leaf.byte_size();

         return res;
      }

      static std::string name() {
         return "AdaptiveLearnedMMPHF";
      }
   };
} // namespace exotic_hashing

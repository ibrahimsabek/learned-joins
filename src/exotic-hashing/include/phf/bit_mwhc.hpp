#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <queue>
#include <random>
#include <sdsl/bit_vector_il.hpp>
#include <sdsl/io.hpp>
#include <sdsl/rrr_vector.hpp>
#include <sdsl/vectors.hpp>
#include <stack>
#include <string>
#include <tuple>
#include <vector>

#include <emmintrin.h>
#include <mmintrin.h>

#include <hashing.hpp>

#include "include/omphf/mwhc.hpp"

// order is important
#include "../convenience/builtins.hpp"

namespace exotic_hashing {
   template<class Data, class Hasher = support::Hasher<Data>, class HyperGraph = support::HyperGraph<Data, Hasher>>
   class BitMWHC {
      hashing::reduction::FastModulo<std::uint64_t> mod_N;
      Hasher hasher;

      sdsl::int_vector<2> vertex_values;

      static forceinline size_t vertices_count(const size_t& dataset_size, const long double& overalloc = 1.23) {
         return std::ceil(overalloc * dataset_size);
      }

     public:
      template<class RandomIt>
      BitMWHC(RandomIt begin, RandomIt end)
         : mod_N(vertices_count(std::distance(begin, end))), hasher(mod_N.N), vertex_values(mod_N.N, 3) {
         // Find suitable peel order
         std::vector<size_t> peel_order;
         while (peel_order.empty()) {
            // 1. Generate random Hypergraph
            hasher = Hasher(mod_N.N);
            HyperGraph g(begin, end, hasher, mod_N.N);

            // 2. Peel (i.e., check for acyclicity)
            peel_order = g.peel();
         }

         // 3. Assign values to vertices depending on reverse peel order
         for (auto it = peel_order.rbegin(); it != peel_order.rend(); it++) {
            // get next edge
            const auto edge_ind = *it;

            // lookup vertices of this edge
            const auto [h0, h1, h2] = hasher(*(begin + edge_ind));

            // compute vertex assign value. Handle the case that there are collisions, i.e. same h value twice
            const size_t curr_value =
               (vertex_values[h0] + (h1 != h0) * vertex_values[h1] + (h2 != h1 && h2 != h0) * vertex_values[h2]) % 3;

            // assign vertex values
            bool assigned = false;
            if (vertex_values[h0] == 3) {
               vertex_values[h0] = (3 + 0 - curr_value) % 3;
               assigned = true;
            }

            if (vertex_values[h1] == 3) {
               vertex_values[h1] = assigned ? 0 : (3 + 1 - curr_value) % 3;
               assigned = true;
            }

            if (vertex_values[h2] == 3)
               vertex_values[h2] = assigned ? 0 : (3 + 2 - curr_value) % 3;
         }
      }

      explicit BitMWHC(const std::vector<Data>& dataset) : BitMWHC(dataset.begin(), dataset.end()) {}

      static std::string name() {
         return "BitMWHC";
      }

      forceinline size_t operator()(const Data& key) const {
         const auto [h0, h1, h2] = hasher(key);
         size_t hash = vertex_values[h0];
         if (h1 != h0)
            hash += vertex_values[h1];
         if (h2 != h1 && h2 != h0)
            hash += vertex_values[h2];

         const std::array<size_t, 3> hashs{h0, h1, h2};
         return hashs[hash % 3];
      }

      size_t byte_size() const {
         return sizeof(hasher) + sizeof(mod_N) + sdsl::size_in_bytes(vertex_values);
      }
   };
} // namespace exotic_hashing

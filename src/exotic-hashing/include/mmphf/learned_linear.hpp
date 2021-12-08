#pragma once

#include <cstdint>
#include <sdsl/bit_vector_il.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>
#include <sdsl/util.hpp>
#include <string>

#include "../convenience/builtins.hpp"

namespace exotic_hashing {
   template<class Data>
   class LearnedLinear {
      Data neg_intercept = 0;
      sdsl::bit_vector_il<> bitvec{};
      decltype(bitvec)::rank_1_type rank{};

      forceinline size_t rank_index(const Data& key) const {
         // slope is always 1 since we don't want to produce collisions
         return key - neg_intercept;
      }

     public:
      LearnedLinear() noexcept = default;

      /**
       * Constructs from *already sorted* range [begin, end)
       */
      template<class It>
      LearnedLinear(const It& begin, const It& end) {
         construct(begin, end);
      }

      /**
       * Constructs given any dataset. Note: does not have
       * to be sorted
       */
      explicit LearnedLinear(std::vector<Data> dataset) {
         std::sort(dataset.begin(), dataset.end());
         construct(dataset.begin(), dataset.end());
      }

      /**
       * Constructs from *already sorted* range [begin, end)
       */
      template<class It>
      void construct(const It& begin, const It& end) {
         const size_t size = std::distance(begin, end);

         // Compute dataset properties
         const auto min = *begin;
         const auto max = *(end - 1);
         const size_t scale = (max - min + 1);

         // Set intercept
         neg_intercept = min;
         assert(neg_intercept == min);

         // build rank bitvector
         sdsl::bit_vector bv(scale, 0);
         for (size_t i = 0; i < size; i++) {
            const auto key = *(begin + i);
            const auto ind = rank_index(key);

            assert(ind < scale);
            bv[ind] = true;
         }
         bitvec = bv;
         sdsl::util::init_support(rank, &bitvec);
      }

      forceinline size_t operator()(const Data& key) const {
         const size_t ind = rank_index(key);
         assert(ind < bitvec.size());
         const size_t res = rank(ind);

         return res;
      }

      forceinline size_t byte_size() const {
         return sdsl::size_in_bytes(bitvec) + sdsl::size_in_bytes(rank) + sizeof(decltype(neg_intercept));
      }

      static std::string name() {
         return "LearnedLinear";
      }

      /// Custom copy constructor is necessary since sdsl's rank support contains a pointer to bitvec
      LearnedLinear(const LearnedLinear& other) noexcept {
         neg_intercept = other.neg_intercept;
         bitvec = other.bitvec;
         rank = other.rank;

         // reset vector as otherwise rank contains broken pointer
         rank.set_vector(&bitvec);
      }

      /// Custom copy constructor is necessary since sdsl's rank support contains a pointer to bitvec
      LearnedLinear(LearnedLinear&& other) noexcept {
         neg_intercept = other.neg_intercept;
         bitvec = other.bitvec;
         rank = other.rank;

         // reset vector as otherwise rank contains broken pointer
         rank.set_vector(&bitvec);
      }

      /// Custom copy constructor is necessary since sdsl's select support contains a pointer to upper
      LearnedLinear& operator=(const LearnedLinear& other) noexcept {
         if (this != &other) {
            neg_intercept = other.neg_intercept;
            bitvec = other.bitvec;
            rank = other.rank;

            // reset vector as otherwise rank contains broken pointer
            rank.set_vector(&bitvec);
         }

         return *this;
      }

      /// Custom copy constructor is necessary since sdsl's select support contains a pointer to upper
      LearnedLinear& operator=(LearnedLinear&& other) noexcept {
         neg_intercept = other.neg_intercept;
         bitvec = other.bitvec;
         rank = other.rank;

         // reset vector as otherwise rank contains broken pointer
         rank.set_vector(&bitvec);

         return *this;
      }

      // Destructor does not have to do any special work
      ~LearnedLinear() noexcept = default;
   };
} // namespace exotic_hashing

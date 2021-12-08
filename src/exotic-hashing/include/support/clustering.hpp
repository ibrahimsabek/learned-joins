#pragma once

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

#include "../convenience/builtins.hpp"

namespace exotic_hashing::support {
   struct DensityGauge {
      /**
       * Measures density of a given region [begin, end). Note that the region
       * *must be sorted*
       */
      template<class ConstIt, class Precision = double>
      forceinline Precision operator()(const ConstIt& begin, const ConstIt& end) const {
         assert(begin < end);
         assert(std::is_sorted(begin, end));

         // Since range is sorted we can simple fetch first & last element
         const auto min_elem = *begin;
         const auto max_elem = *(end - 1);

         const auto range_size = std::distance(begin, end);

         // Density defined as |X|/(max(X)-min(X))
         return static_cast<Precision>(range_size) / static_cast<Precision>(1 + max_elem - min_elem);
      }
   };

   /**
   * Clusters a given *pre sorted* range based on a certain Gauge
   *
   * @param begin iterator pointing to first element of the sorted range
   * @param end past the end iterator of the sorted range
   * @param score_threshold threshold until which merging is considered benificial
   * @param gauge Gauges are functors can assign a score \in [0, 1] to arbitraty,
   *        consecutive, sorted regions [a, b). The higher the score, the better
   *
   * @returns a vector of iterators pointing to the clustered ranges' boundaries, i.e.,
   *  first range is [vec[0], vec[1]), second range is [vec[1], vec[2]), ... until [vec[vec.size()-2], vec[vec.size()-1])
   */
   template<class RandomIt, class ClusterGauge = DensityGauge>
   std::vector<RandomIt> cluster(const RandomIt& begin, const RandomIt& end, const double& score_threshold,
                                 const ClusterGauge gauge = ClusterGauge()) {
      // Start out with one element per region.
      const size_t regions_count = std::distance(begin, end);
      size_t curr_size = regions_count + 1, next_i = 0;

      // Empty range -> zero clusters
      if (regions_count == 0)
         return {};
      if (regions_count == 1)
         return {begin, end};

      // Verify assumptions before proceeding
      assert(begin < end);
      assert(std::is_sorted(begin, end));
      assert(score_threshold >= 0);
      assert(score_threshold <= 1);

      // To eliminate additional alloc/deallocs use a double buffer technique
      const auto swap = [](auto& a, auto& b) {
         const auto tmp = a;
         a = b;
         b = tmp;
      };
      auto *curr_regions = new RandomIt[curr_size], *next_regions = new RandomIt[curr_size];

      // Initialize curr_regions
      for (size_t i = 0; i <= regions_count; i++) {
         assert(begin + i <= end);
         curr_regions[i] = begin + i;
      }

      // Iteratively merge the most promising neighboring regions as measured
      // by gauge, until no further progress is made
      for (size_t round = 0;; round++) {
         // Verify assumptions
         assert(curr_size > 0);

         // Merging pass where we only ever merge left should left exceed threshold
         // and should the right merge's score be lower
         for (size_t i = 1; i + 1 < curr_size; i++) {
            // Measure what would happen if we were to merge left
            const auto left_measurement = gauge(curr_regions[i - 1], curr_regions[i + 1]);

            // Commit left begin (won't change this round)
            if (next_i == 0 || next_regions[next_i - 1] != curr_regions[i - 1])
               next_regions[next_i++] = curr_regions[i - 1];

            // Don't merge at all if the left's score is lower than threshold or if the
            // right's score is higher than left's (only check if we have a right neighbor region at all)
            if (left_measurement >= score_threshold &&
                (i + 2 >= curr_size || gauge(curr_regions[i], curr_regions[i + 2]) <= left_measurement)) {
               next_regions[next_i++] = curr_regions[i + 1];
               i++;
            } else {
               next_regions[next_i++] = curr_regions[i];
            }
         }

         // Next must contain at least 2 entries, otherwise we're left with a
         // single (broken) range
         assert(next_i >= 2);

         // Don't forget to write out end iterator
         assert(curr_regions[curr_size - 1] == end);
         if (next_regions[next_i - 1] < end)
            next_regions[next_i++] = end;

         // Advance state (swap buffers etc)
         swap(curr_regions, next_regions);
         if (curr_size == next_i || next_i <= 2)
            break;
         curr_size = next_i;
         next_i = 0;
      }

      // copy data into vector for further use (TODO(dominik): eliminate this copy)
      std::vector<RandomIt> regions(curr_regions, curr_regions + next_i);

      // ensure we leak no memory
      delete[] curr_regions;
      delete[] next_regions;

      return regions;
   }

} // namespace exotic_hashing::support

#pragma once

#include <string>
#include <unordered_map>

#include "../convenience/builtins.hpp"

namespace exotic_hashing {
   /**
    * Naive order preserving mphf baseline using
    * as standard map container
    */
   template<class Data>
   struct MapOMPHF {
      explicit MapOMPHF(const std::vector<Data>& dataset) {
         for (size_t i = 0; i < dataset.size(); i++) {
            map[dataset[i]] = i;
         }
      }

      static std::string name() {
         return "MapOMPHF";
      }

      constexpr forceinline size_t operator()(const Data& key) const {
         return map.at(key);
      }

      size_t byte_size() const {
         size_t internal_map_size = 0;
         for (size_t i = 0; i < map.bucket_count(); i++)
            // This is a rough estimate and could vary from the actual size!
            internal_map_size += map.bucket_size(i) * (sizeof(Data) + sizeof(size_t)) + sizeof(void*);

         return sizeof(MapOMPHF<Data>) + internal_map_size;
      };

     private:
      std::unordered_map<Data, size_t> map;
   };
} // namespace exotic_hashing

#pragma once

#include "utils/data_structures.h"

template<typename KeyType, typename PayloadType, typename TaskType, typename JoinThreadType, typename PartitionType>
class HashJoinSteps {
 public:  
  void partition_rel_r_segment_pass1(PartitionType * const part_r) {}
  
  void partition_rel_s_segment_pass1(PartitionType * const part_s) {}

  template<typename BuildType>
  void build_rel_r_partition(BuildType * build_output, const Relation<KeyType, PayloadType> * const rel_r_partition, Relation<KeyType, PayloadType> * const tmp_r) {}
  
  template<typename BuildType>
  void build_rel_r_partition_imv(BuildType * build_output, const Relation<KeyType, PayloadType> * const rel_r_partition, Relation<KeyType, PayloadType> * const tmp_r) {}
  
  template<typename BuildType>
  uint64_t probe_rel_s_partition(const Relation<KeyType, PayloadType> * const rel_r_partition, const Relation<KeyType, PayloadType> * const rel_s_partition,
                                 const BuildType * build_output) {
    return 0;
  }

  template<typename BuildType>
  uint64_t probe_rel_s_partition_imv(const Relation<KeyType, PayloadType> * const rel_r_partition, const Relation<KeyType, PayloadType> * const rel_s_partition,
                                 const BuildType * build_output, int thread_id) {
    return 0;
  }

};

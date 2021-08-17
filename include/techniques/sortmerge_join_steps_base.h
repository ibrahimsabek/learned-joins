#pragma once

/* A generic base class for implementing sort-merge join steps  along with some utilities implementation based on https://systems.ethz.ch/research/data-processing-on-modern-hardware/projects/parallel-and-distributed-joins.html */

#include "utils/data_structures.h"
#include "utils/math.h"
#include "utils/memory.h"

template<typename KeyType, typename PayloadType, typename JoinThreadType>
class SortMergeJoinSteps {
 public:  
  
  void partition_phase(Relation<KeyType, PayloadType> *** relRparts, Relation<KeyType, PayloadType> *** relSparts, JoinThreadType * in_args) {}
  
  void sorting_phase(Relation<KeyType, PayloadType> ** relRparts, Relation<KeyType, PayloadType> ** relSparts, JoinThreadType * in_args) {}

  void multiwaymerge_phase(int numaregionid, Relation<KeyType, PayloadType> ** relRparts, Relation<KeyType, PayloadType> ** relSparts, JoinThreadType * in_args,
                                             Relation<KeyType, PayloadType> * mergedRelR, Relation<KeyType, PayloadType> * mergedRelS) {}

  void mergejoin_phase(Relation<KeyType, PayloadType> ** relRparts, Relation<KeyType, PayloadType> ** relSparts,
                       Relation<KeyType, PayloadType> * mergedRelR, Relation<KeyType, PayloadType> * mergedRelS, JoinThreadType * in_args) {}

  void partitioning_cleanup(Relation<KeyType, PayloadType> ** relRparts, Relation<KeyType, PayloadType> ** relSparts) {}

  // Based on the ETH implementation of sortmerge join
  void radix_cluster_optimized_V2(Relation<KeyType, PayloadType> * restrict outRel, 
                                  Relation<KeyType, PayloadType> * restrict inRel,
                                  int32_t * restrict hist, 
                                  int R, 
                                  int D)
  {
      uint32_t i;
      uint32_t offset = 0;
      const uint32_t M       = ((1 << D) - 1) << R;
      const uint32_t fanOut  = 1 << D;
      const uint32_t ntuples = inRel->num_tuples;

      Tuple<KeyType, PayloadType> * input  = inRel->tuples;
      Tuple<KeyType, PayloadType> * output = outRel->tuples;

      CacheLine<KeyType, PayloadType> buffer[fanOut] __attribute__((aligned(CACHE_LINE_SIZE)));

      for( i = 0; i < ntuples; i++ ){
          uint32_t idx = HASH_BIT_MODULO(input->key, M, R);
          hist[idx]++;
          input++;
      }

      /* determine the start and end of each cluster depending on the counts. */
      uint32_t TUPLESPERCACHELINE = CACHE_LINE_SIZE/sizeof(Tuple<KeyType, PayloadType>);
      for ( i = 0; i < fanOut; i++ ) {
          buffer[i].data.slot = offset;
          /* for aligning partition-outputs to cacheline: */
          /* hist[i] = (hist[i] + (64/sizeof(tuple_t))) & ~((64/sizeof(tuple_t))-1); */
          //offset += hist[i];
          offset += ALIGN_NUMTUPLES(TUPLESPERCACHELINE, hist[i]);
      }

      input = inRel->tuples;
      /* copy tuples to their corresponding clusters at appropriate offsets */
      for( i = 0; i < ntuples; i++ ){
          uint32_t  idx     = HASH_BIT_MODULO(input->key, M, R);
          /* store in the cache-resident buffer first */
          uint32_t  slot    = buffer[idx].data.slot;
          Tuple<KeyType, PayloadType> * tup     = (Tuple<KeyType, PayloadType> *)(buffer + idx);
          uint32_t  slotMod = (slot) & (TUPLESPERCACHELINE - 1); /* % operator */
          tup[slotMod] = *input;
          input ++;

          if(slotMod == (TUPLESPERCACHELINE-1)){
              /* uintptr_t p = (uintptr_t)(output+slot-TUPLESPERCACHELINE); */
              /* if(p % 64 != 0){ */
              /*     printf("there is a problem\n"); */
              /* } */
              store_nontemp_64B<KeyType, PayloadType>((output+slot-(TUPLESPERCACHELINE-1)), (buffer+idx));
          }
          buffer[idx].data.slot = slot+1;
      }

      /* flush the remainder tuples in the buffer */
      for ( i = 0; i < fanOut; i++ ) {
          uint32_t  slot = buffer[i].data.slot;
          uint32_t  num  = (slot) & (TUPLESPERCACHELINE - 1);
          if(num > 0){
              Tuple<KeyType, PayloadType> * dest = output + slot - num;

              for(uint32_t j = 0; j < num; j++) {
                  dest[j] = buffer[i].data.tuples[j];
              }
          }
      }
  }


  void partition_relation_optimized_V2(Relation<KeyType, PayloadType> ** partitions, 
                                       Relation<KeyType, PayloadType> * input, 
                                       Relation<KeyType, PayloadType> * output,
                                       uint32_t nbits,
                                       uint32_t shiftbits)
  {
      int i;
      uint32_t offset = 0;
      const int fanOut = 1 << nbits;

      int32_t * hist, * histAligned;
      hist = (int32_t*) calloc(fanOut + 16, sizeof(int32_t));
      histAligned = (int32_t *) ALIGNPTR(hist, CACHE_LINE_SIZE);

      /* int32_t histAligned[fanOut+1] __attribute__((aligned(CACHE_LINE_SIZE))); */
      /* memset(histAligned, 0, (fanOut+1)*sizeof(int32_t)); */

      radix_cluster_optimized_V2(output, input, histAligned, shiftbits, nbits);

      int TUPLESPERCACHELINE = CACHE_LINE_SIZE/sizeof(Tuple<KeyType, PayloadType>);
      for(i = 0; i < fanOut; i++) 
      {
        Relation<KeyType, PayloadType> * part = partitions[i];
        part->num_tuples = histAligned[i];
        part->tuples = output->tuples + offset;

        offset += ALIGN_NUMTUPLES(TUPLESPERCACHELINE, histAligned[i]);
      }

      free(hist);
  }

};

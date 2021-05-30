#include "vectorwise/Operators.hpp"

#include <assert.h>

#include "common/Compat.hpp"
#include "common/runtime/Concurrency.hpp"
#include "common/runtime/SIMD.hpp"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <tuple>
#include <x86intrin.h>
#include <string>

#include "common/runtime/Hash.hpp"

namespace vectorwise {
#define __AVX512F__ 1
using runtime::barrier;

size_t Select::next() {
   while (true) {
      auto n = child->next();
      if (n == EndOfStream) return EndOfStream;
      n = condition->evaluate(n);
      if (n > 0) return n;
   }
}

size_t Project::next() {
   auto n = child->next();
   if (n == EndOfStream) return EndOfStream;
   for (auto& expression : expressions) expression->evaluate(n);
   return n;
}

size_t FixedAggr::next() {
   if (!consumed) {
      size_t found = 0;
      for (auto n = child->next(); n != EndOfStream; n = child->next()) {
         found = aggregates.evaluate(n);
      }
      consumed = true;
      return found;
   } else {
      return EndOfStream;
   }
}

Scan::Scan(Shared& s, size_t n, size_t v)
    : shared(s), needsInit(true), currentChunk(0), lastOffset(0), nrTuples(n),
      vecSize(v) {
   scanChunkSize = 1;
   size_t scanMorselSize = 1024 * 10;
   if (vecSize < scanMorselSize) scanChunkSize = scanMorselSize / vecSize + 1;
   vecInChunk = scanChunkSize;
}

void Scan::addConsumer(void** colPtr, size_t typeSize) {
   consumers.emplace_back(colPtr, vecSize * typeSize);
}

size_t Scan::next() {
   auto step = 1;

   if (vecInChunk == scanChunkSize) {
      auto prevChunk = currentChunk;
      currentChunk = shared.pos.fetch_add(1);
      auto chunkSkip = currentChunk - prevChunk;
      if (needsInit) {
         step = chunkSkip * scanChunkSize;
         needsInit = false;
      } else {
         chunkSkip -= 1;
         step = chunkSkip * scanChunkSize + 1;
      }
      vecInChunk = 0;
   }

   auto nextBegin = lastOffset + step * vecSize;
   if (nextBegin >= nrTuples) return EndOfStream;
   auto nextBatchSize = std::min(nrTuples - nextBegin, vecSize);
   for (auto& cons : consumers)
      *cons.first = (void*)(*(uint8_t**)cons.first + step * cons.second);
   lastOffset = nextBegin;
   vecInChunk++;
   return nextBatchSize;
}

ResultWriter::Input::Input(void* d, size_t size,
                           runtime::BlockRelation::Attribute attr)
    : data(d), elementSize(size), attribute(attr) {}

ResultWriter::ResultWriter(Shared& s)
    : shared(s), currentBlock(nullptr, nullptr) {}

size_t ResultWriter::next() {
   size_t found = 0;
   for (pos_t n = child->next(); n != EndOfStream; n = child->next()) {
      found += n;
      // assure that enough space is available in current block to fit result of
      // all buffers
      if (currentBlock.spaceRemaining() < n)
         currentBlock = shared.result->result->createBlock(n);
      auto blockSize = currentBlock.size();
      for (const auto& input : inputs)
         // copy data from intermediate buffers into result relation
         std::memcpy(addBytes(currentBlock.data(input.attribute),
                              input.elementSize * blockSize),
                     input.data, n * input.elementSize);
      // update result relation size
      currentBlock.addedElements(n);
   }
   return found;
}
#define VECTORSIZE 8
int times=0;
void printVector(__m512i& vec, std::string str) {
  if(times>20) return;
  times++;
   uint64_t * ptr = (uint64_t*)&vec;
   std::cout<<str<<"  ";
   for(int i=0;i<VECTORSIZE;++i) {
     std::cout<<"vec["<< i <<"]=  "<<ptr[i]<<"    ";
   }
   std::cout<<std::endl;
}
void gather(__m512i& index, int* probe_keys) {
  if(times>20) return;
  uint64_t * ptr = (uint64_t*)&index;
  std::cout<<"gather keys ";
  for(int i=0;i<VECTORSIZE;++i) {
    std::cout<<"vec["<< ptr[i] <<"]=  "<<probe_keys[ptr[i]]<<"    ";
  }
  std::cout<<std::endl;
}
pos_t Hashjoin::joinSIMDAMAC() {
  //std::cout <<"JoinSIMDAMAC "<<"\n";
  size_t found=0;
  if(imv_cont.k==100) {
    for(int i=0;i<imvNum;++i) {
      imv_state[i]->reset();
    }
    imv_cont.k=0;
    imv_cont.done=0;
  }
  int* probeKeys = (reinterpret_cast<int*>(((F2_Op *)(probeHash.ops[0].get()))->param1));
  uint64_t keyOff = ((EqualityCheck*)(keyEquality.ops[0].get()))->offset;
  __mmask8 m_match = 0,  m_new_probes = -1;
  __m512i v_base_offset = _mm512_set_epi64(7,6,5,4,3,2,1,0);
  __m512i v_offset = _mm512_set1_epi64(0),v_new_build_key,v_build_keys;
  __m512i v_base_offset_upper = _mm512_set1_epi64(cont.numProbes);
  __m512i v_seed= _mm512_set1_epi64(primitives::seed),v_build_key_off = _mm512_set1_epi64(keyOff);
  __m512i v_probe_hash= _mm512_set1_epi64(0),v_zero=_mm512_set1_epi64(0);
  __m256i v256_zero= _mm256_set1_epi32(0),v256_probe_keys,v256_build_keys;

  while(imv_cont.done<imvNum) {
    imv_cont.k = (imv_cont.k>=imvNum)? 0 :imv_cont.k;
    if(cont.nextProbe >= cont.numProbes) {
      if(imv_state[imv_cont.k]->m_valid_probe==0 && imv_state[imv_cont.k]->stage!=3) {
        ++imv_cont.done;
        imv_state[imv_cont.k]->stage = 3;
        ++imv_cont.k;
        continue;
      }
    }
    switch(imv_state[imv_cont.k]->stage) {
      case 1: {
        /// step 1: load the offsets of probing tuples
          v_offset = _mm512_add_epi64(_mm512_set1_epi64(cont.nextProbe), v_base_offset);
          imv_state[imv_cont.k]->v_probe_offset = _mm512_mask_expand_epi64(imv_state[imv_cont.k]->v_probe_offset, _mm512_knot(imv_state[imv_cont.k]->m_valid_probe), v_offset);
          // count the number of empty tuples
          m_new_probes = _mm512_knot(imv_state[imv_cont.k]->m_valid_probe);
          cont.nextProbe = cont.nextProbe + _mm_popcnt_u32(m_new_probes);
          imv_state[imv_cont.k]->m_valid_probe = _mm512_cmpgt_epu64_mask(v_base_offset_upper, imv_state[imv_cont.k]->v_probe_offset);
          m_new_probes = _mm512_kand(m_new_probes, imv_state[imv_cont.k]->m_valid_probe);
        /// step 2: gather the probe keys
          v256_probe_keys= _mm512_mask_i64gather_epi32(v256_zero,m_new_probes, imv_state[imv_cont.k]->v_probe_offset, (void*)probeKeys, 4);
          imv_state[imv_cont.k]->v_probe_keys=_mm512_mask_blend_epi64(m_new_probes,imv_state[imv_cont.k]->v_probe_keys,_mm512_cvtepi32_epi64(v256_probe_keys));
        /// step 3: compute the hash values of probe keys
          v_probe_hash = runtime::MurMurHash()((imv_state[imv_cont.k]->v_probe_keys),(v_seed));
        /// step 4: find the addresses of corresponding buckets for new probes
          Vec8uM v_new_bucket_addrs = shared.ht.find_chain_tagged_sel((v_probe_hash),m_new_probes);
          // the addresses are null, then the corresponding probes are invalid
          imv_state[imv_cont.k]->m_valid_probe = _mm512_kand(_mm512_kor(_mm512_knot(m_new_probes),v_new_bucket_addrs.mask), imv_state[imv_cont.k]->m_valid_probe);
          imv_state[imv_cont.k]->v_bucket_addrs =  _mm512_mask_blend_epi64(v_new_bucket_addrs.mask,imv_state[imv_cont.k]->v_bucket_addrs,v_new_bucket_addrs.vec);
          imv_state[imv_cont.k]->stage = 0;
          uint64_t * ht_pos = (uint64_t *)&imv_state[imv_cont.k]->v_bucket_addrs;
          for (int i = 0; i < VECTORSIZE; ++i) {
            _mm_prefetch((char *)(ht_pos[i]), _MM_HINT_T0);
            _mm_prefetch(((char *)(ht_pos[i])+64), _MM_HINT_T0);
          }
      }break;
      case 0: {
        /// step 5: gather the all new build keys
          v256_build_keys= _mm512_mask_i64gather_epi32(v256_zero, imv_state[imv_cont.k]->m_valid_probe,
                                  _mm512_add_epi64(imv_state[imv_cont.k]->v_bucket_addrs,v_build_key_off),nullptr, 1);
          v_build_keys =_mm512_cvtepi32_epi64(v256_build_keys);
        /// step 6: compare the probe keys and build keys and write points
          m_match = _mm512_cmpeq_epi64_mask(imv_state[imv_cont.k]->v_probe_keys, v_build_keys);
          m_match = _mm512_kand(m_match, imv_state[imv_cont.k]->m_valid_probe);
          _mm512_mask_compressstoreu_epi64((buildMatches + found), m_match,imv_state[imv_cont.k]->v_bucket_addrs);
          _mm256_mask_compressstoreu_epi32((probeMatches + found), m_match, _mm512_cvtepi64_epi32(imv_state[imv_cont.k]->v_probe_offset));

          found+=_mm_popcnt_u32(m_match);
        /// step 7: move to the next bucket nodes
          imv_state[imv_cont.k]->v_bucket_addrs = _mm512_mask_i64gather_epi64(v_zero,imv_state[imv_cont.k]->m_valid_probe, imv_state[imv_cont.k]->v_bucket_addrs, nullptr, 1);
          imv_state[imv_cont.k]->m_valid_probe =_mm512_kand(imv_state[imv_cont.k]->m_valid_probe,_mm512_cmpneq_epi64_mask(imv_state[imv_cont.k]->v_bucket_addrs, v_zero));
          imv_state[imv_cont.k]->stage = 1;
          if(found+VECTORSIZE >= batchSize) {
            return found;
          }
      }break;
    }
    ++imv_cont.k;
  }
  imv_cont.k=100;
  return found;
}

pos_t Hashjoin::joinIMV() {
  //std::cout <<"joinIMV "<<"\n";
  size_t found = 0;
  if (imv_cont.k == 100) {
    for (int i = 0; i <= imvNum; ++i) {
      imv_state[i]->reset();
    }
    imv_cont.k = 0;
    imv_cont.done = 0;
  }
  int* probeKeys = (reinterpret_cast<int*>(((F2_Op *) (probeHash.ops[0].get()))->param1));
  uint64_t keyOff = ((EqualityCheck*) (keyEquality.ops[0].get()))->offset;
  __attribute__((aligned(64)))  __mmask8 m_match = 0, m_new_probes = -1, mask[VECTORSIZE + 1];

  __m512i v_base_offset = _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0);
  __m512i v_offset = _mm512_set1_epi64(0), v_new_build_key, v_build_keys;
  __m512i v_base_offset_upper = _mm512_set1_epi64(cont.numProbes);
  __m512i v_seed = _mm512_set1_epi64(primitives::seed), v_build_key_off = _mm512_set1_epi64(keyOff);
  __m512i  v_zero = _mm512_set1_epi64(0);
  __m256i v256_zero = _mm256_set1_epi32(0), v256_probe_keys, v256_build_keys;
  uint64_t * ht_pos = nullptr;
  uint8_t num, num_temp;
  for (int i = 0; i <= VECTORSIZE; ++i) {
    mask[i] = (1 << i) - 1;
  }
  while (true) {
    imv_cont.k = (imv_cont.k >= imvNum) ? 0 : imv_cont.k;
    if (cont.nextProbe >= cont.numProbes) {
      if (imv_state[imv_cont.k]->m_valid_probe == 0 && imv_state[imv_cont.k]->stage != 3) {
        ++imv_cont.done;
        imv_state[imv_cont.k]->stage = 3;
        ++imv_cont.k;
        continue;
      }
    }
    if (imv_cont.done >= imvNum) {
      if (imv_state[imvNum]->m_valid_probe > 0) {
        imv_cont.k = imvNum;
        imv_state[imvNum]->stage = 0;
      } else {
        break;
      }
    }
    switch (imv_state[imv_cont.k]->stage) {
      case 1: {
        /// step 1: load the offsets of probing tuples
        imv_state[imv_cont.k]->v_probe_offset = _mm512_add_epi64(_mm512_set1_epi64(cont.nextProbe), v_base_offset);
        imv_state[imv_cont.k]->m_valid_probe = _mm512_cmpgt_epu64_mask(v_base_offset_upper, imv_state[imv_cont.k]->v_probe_offset);
        /// step 2: gather the probe keys
        //v256_probe_keys = _mm512_mask_i64gather_epi32(v256_zero, imv_state[imv_cont.k]->m_valid_probe, imv_state[imv_cont.k]->v_probe_offset, (void* )probeKeys, 4);
        v256_probe_keys = _mm256_maskz_loadu_epi32(imv_state[imv_cont.k]->m_valid_probe, (char*)(probeKeys+cont.nextProbe));

        imv_state[imv_cont.k]->v_probe_keys = _mm512_cvtepi32_epi64(v256_probe_keys);
        /// step 3: compute the hash values of probe keys
        imv_state[imv_cont.k]->v_probe_hash = runtime::MurMurHash()((imv_state[imv_cont.k]->v_probe_keys), (v_seed));
        cont.nextProbe += VECTORSIZE;
        imv_state[imv_cont.k]->stage = 2;
        shared.ht.prefetchEntry((imv_state[imv_cont.k]->v_probe_hash));
      }break;
      case 2:{
        /// step 4: find the addresses of corresponding buckets for new probes
        Vec8uM v_new_bucket_addrs = shared.ht.find_chain_tagged((imv_state[imv_cont.k]->v_probe_hash));
        imv_state[imv_cont.k]->m_valid_probe = _mm512_kand(imv_state[imv_cont.k]->m_valid_probe, v_new_bucket_addrs.mask);
        imv_state[imv_cont.k]->v_bucket_addrs = v_new_bucket_addrs.vec;
        imv_state[imv_cont.k]->stage = 0;
        ht_pos = (uint64_t *) &imv_state[imv_cont.k]->v_bucket_addrs;
        for (int i = 0; i < VECTORSIZE; ++i) {
          _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
          _mm_prefetch(((char * )(ht_pos[i]) + 64), _MM_HINT_T0);

        }
      }
        break;
      case 0: {
        /// step 5: gather the all new build keys
        v256_build_keys = _mm512_mask_i64gather_epi32(v256_zero, imv_state[imv_cont.k]->m_valid_probe, _mm512_add_epi64(imv_state[imv_cont.k]->v_bucket_addrs, v_build_key_off), nullptr, 1);
        v_build_keys = _mm512_cvtepi32_epi64(v256_build_keys);
        /// step 6: compare the probe keys and build keys and write points
        m_match = _mm512_cmpeq_epi64_mask(imv_state[imv_cont.k]->v_probe_keys, v_build_keys);
        m_match = _mm512_kand(m_match, imv_state[imv_cont.k]->m_valid_probe);
        _mm512_mask_compressstoreu_epi64((buildMatches + found), m_match, imv_state[imv_cont.k]->v_bucket_addrs);
        _mm256_mask_compressstoreu_epi32((probeMatches + found), m_match, _mm512_cvtepi64_epi32(imv_state[imv_cont.k]->v_probe_offset));

        found += _mm_popcnt_u32(m_match);
        /// step 7: move to the next bucket nodes
        imv_state[imv_cont.k]->v_bucket_addrs = _mm512_mask_i64gather_epi64(v_zero, imv_state[imv_cont.k]->m_valid_probe, imv_state[imv_cont.k]->v_bucket_addrs, nullptr, 1);
        imv_state[imv_cont.k]->m_valid_probe = _mm512_kand(imv_state[imv_cont.k]->m_valid_probe, _mm512_cmpneq_epi64_mask(imv_state[imv_cont.k]->v_bucket_addrs, v_zero));

        num = _mm_popcnt_u32(imv_state[imv_cont.k]->m_valid_probe);
        if (num == VECTORSIZE) {
          ht_pos = (uint64_t *) &imv_state[imv_cont.k]->v_bucket_addrs;
          for (int i = 0; i < VECTORSIZE; ++i) {
            _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
            _mm_prefetch(((char * )(ht_pos[i]) + 64), _MM_HINT_T0);
          }
        } else {
          if (imv_cont.done < imvNum) {
            num_temp = _mm_popcnt_u32(imv_state[imvNum]->m_valid_probe);
            if (num + num_temp < VECTORSIZE) {
              // compress imv_state[imv_cont.k]
              compress(imv_state[imv_cont.k]);
              // expand imv_state[imv_cont.k] -> imv_state[imvNum]
              expand(imv_state[imv_cont.k],imv_state[imvNum]);
              imv_state[imvNum]->m_valid_probe = mask[num + num_temp];
              imv_state[imv_cont.k]->m_valid_probe = 0;
              imv_state[imv_cont.k]->stage = 1;
            } else {
              // expand imv_state[imvNum] -> expand imv_state[imv_cont.k]
              expand(imv_state[imvNum],imv_state[imv_cont.k]);
              imv_state[imvNum]->m_valid_probe = _mm512_kand(imv_state[imvNum]->m_valid_probe, _mm512_knot(mask[VECTORSIZE - num]));
              // compress imv_state[imvNum]
              compress(imv_state[imvNum]);
              imv_state[imvNum]->m_valid_probe = imv_state[imvNum]->m_valid_probe >> (VECTORSIZE - num);
              imv_state[imv_cont.k]->m_valid_probe = mask[VECTORSIZE];
              imv_state[imv_cont.k]->stage = 0;
              ht_pos = (uint64_t *) &imv_state[imv_cont.k]->v_bucket_addrs;
              for (int i = 0; i < VECTORSIZE; ++i) {
                _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
                _mm_prefetch(((char * )(ht_pos[i]) + 64), _MM_HINT_T0);
              }
            }
          }
        }
        if (found + VECTORSIZE >= batchSize) {
          return found;
        }
      }
        break;
    }
    ++imv_cont.k;
  }
  imv_cont.k = 100;
  return found;
}
#define PDISD 128
pos_t Hashjoin::joinSelIMV() {
    //std::cout <<"joinSelIMV "<<"\n";

  size_t found = 0;
  if (imv_cont.k == 100) {
    for (int i = 0; i <= imvNum; ++i) {
      imv_state[i]->reset();
    }
    imv_cont.k = 0;
    imv_cont.done = 0;
  }
  int* probeKeys = (reinterpret_cast<int*>(((F3_Op *) (probeHash.ops[0].get()))->param2));
  uint64_t keyOff = ((EqualityCheck*) (keyEquality.ops[0].get()))->offset;
  __attribute__((aligned(64)))  __mmask8 m_match = 0, m_new_probes = -1, mask[VECTORSIZE + 1];

  __m512i v_base_offset = _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0);
  __m512i v_offset = _mm512_set1_epi64(0), v_new_build_key, v_build_keys;
  __m512i v_base_offset_upper = _mm512_set1_epi64(cont.numProbes);
  __m512i v_seed = _mm512_set1_epi64(primitives::seed), v_build_key_off = _mm512_set1_epi64(keyOff);
  __m512i  v_zero = _mm512_set1_epi64(0);
  __m256i v256_zero = _mm256_set1_epi32(0), v256_probe_keys, v256_build_keys;
  uint64_t * ht_pos = nullptr;
  uint8_t num, num_temp;
  for (int i = 0; i <= VECTORSIZE; ++i) {
    mask[i] = (1 << i) - 1;
  }
  while (true) {
    imv_cont.k = (imv_cont.k >= imvNum) ? 0 : imv_cont.k;
    if (cont.nextProbe >= cont.numProbes) {
      if (imv_state[imv_cont.k]->m_valid_probe == 0 && imv_state[imv_cont.k]->stage != 3) {
        ++imv_cont.done;
        imv_state[imv_cont.k]->stage = 3;
        ++imv_cont.k;
        continue;
      }
    }
    if (imv_cont.done >= imvNum) {
      if (imv_state[imvNum]->m_valid_probe > 0) {
        imv_cont.k = imvNum;
        imv_state[imvNum]->stage = 0;
      } else {
        break;
      }
    }
    switch (imv_state[imv_cont.k]->stage) {
      case 1: {
        /// step 1: load the offsets of probing tuples
        imv_state[imv_cont.k]->v_probe_offset = _mm512_add_epi64(_mm512_set1_epi64(cont.nextProbe), v_base_offset);
        imv_state[imv_cont.k]->m_valid_probe = _mm512_cmpgt_epu64_mask(v_base_offset_upper, imv_state[imv_cont.k]->v_probe_offset);
        imv_state[imv_cont.k]->v_probe_offset = _mm512_cvtepi32_epi64(_mm256_maskz_loadu_epi32(imv_state[imv_cont.k]->m_valid_probe, (char*)(probeSel+cont.nextProbe)));
#if SEQ_PREFETCH ||1
          _mm_prefetch((char*)(probeKeys+ probeSel[cont.nextProbe+PDISD >= cont.numProbes? cont.nextProbe:cont.nextProbe+PDISD]), _MM_HINT_T0);
          _mm_prefetch(((char*)(probeKeys+ probeSel[cont.nextProbe+PDISD >= cont.numProbes? cont.nextProbe:cont.nextProbe+PDISD])+64), _MM_HINT_T0);
#endif
        /// step 2: gather the probe keys
        v256_probe_keys = _mm512_mask_i64gather_epi32(v256_zero, imv_state[imv_cont.k]->m_valid_probe, imv_state[imv_cont.k]->v_probe_offset, (void* )probeKeys, 4);
        //v256_probe_keys = _mm256_maskz_loadu_epi32(imv_state[imv_cont.k]->m_valid_probe, (char*)(probeKeys+cont.nextProbe));

        imv_state[imv_cont.k]->v_probe_keys = _mm512_cvtepi32_epi64(v256_probe_keys);
        /// step 3: compute the hash values of probe keys
        imv_state[imv_cont.k]->v_probe_hash = runtime::MurMurHash()((imv_state[imv_cont.k]->v_probe_keys), (v_seed));
        cont.nextProbe += VECTORSIZE;
        imv_state[imv_cont.k]->stage = 2;
        shared.ht.prefetchEntry((imv_state[imv_cont.k]->v_probe_hash));
      }break;
      case 2:{
        /// step 4: find the addresses of corresponding buckets for new probes
        Vec8uM v_new_bucket_addrs = shared.ht.find_chain_tagged((imv_state[imv_cont.k]->v_probe_hash));
        imv_state[imv_cont.k]->m_valid_probe = _mm512_kand(imv_state[imv_cont.k]->m_valid_probe, v_new_bucket_addrs.mask);
        imv_state[imv_cont.k]->v_bucket_addrs = v_new_bucket_addrs.vec;
        imv_state[imv_cont.k]->stage = 0;
        ht_pos = (uint64_t *) &imv_state[imv_cont.k]->v_bucket_addrs;
        for (int i = 0; i < VECTORSIZE; ++i) {
          _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
          _mm_prefetch(((char * )(ht_pos[i]) + 64), _MM_HINT_T0);

        }
      }
        break;
      case 0: {
        /// step 5: gather the all new build keys
        v256_build_keys = _mm512_mask_i64gather_epi32(v256_zero, imv_state[imv_cont.k]->m_valid_probe, _mm512_add_epi64(imv_state[imv_cont.k]->v_bucket_addrs, v_build_key_off), nullptr, 1);
        v_build_keys = _mm512_cvtepi32_epi64(v256_build_keys);
        /// step 6: compare the probe keys and build keys and write points
        m_match = _mm512_cmpeq_epi64_mask(imv_state[imv_cont.k]->v_probe_keys, v_build_keys);
        m_match = _mm512_kand(m_match, imv_state[imv_cont.k]->m_valid_probe);
        _mm512_mask_compressstoreu_epi64((buildMatches + found), m_match, imv_state[imv_cont.k]->v_bucket_addrs);
        _mm256_mask_compressstoreu_epi32((probeMatches + found), m_match, _mm512_cvtepi64_epi32(imv_state[imv_cont.k]->v_probe_offset));

        found += _mm_popcnt_u32(m_match);
        /// step 7: move to the next bucket nodes
        imv_state[imv_cont.k]->v_bucket_addrs = _mm512_mask_i64gather_epi64(v_zero, imv_state[imv_cont.k]->m_valid_probe, imv_state[imv_cont.k]->v_bucket_addrs, nullptr, 1);
        imv_state[imv_cont.k]->m_valid_probe = _mm512_kand(imv_state[imv_cont.k]->m_valid_probe, _mm512_cmpneq_epi64_mask(imv_state[imv_cont.k]->v_bucket_addrs, v_zero));

        num = _mm_popcnt_u32(imv_state[imv_cont.k]->m_valid_probe);
        if (num == VECTORSIZE) {
          ht_pos = (uint64_t *) &imv_state[imv_cont.k]->v_bucket_addrs;
          for (int i = 0; i < VECTORSIZE; ++i) {
            _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
            _mm_prefetch(((char * )(ht_pos[i]) + 64), _MM_HINT_T0);
          }
        } else {
          if (imv_cont.done < imvNum) {
            num_temp = _mm_popcnt_u32(imv_state[imvNum]->m_valid_probe);
            if (num + num_temp < VECTORSIZE) {
              // compress imv_state[imv_cont.k]
              compress(imv_state[imv_cont.k]);
              // expand imv_state[imv_cont.k] -> imv_state[imvNum]
              expand(imv_state[imv_cont.k],imv_state[imvNum]);
              imv_state[imvNum]->m_valid_probe = mask[num + num_temp];
              imv_state[imv_cont.k]->m_valid_probe = 0;
              imv_state[imv_cont.k]->stage = 1;
            } else {
              // expand imv_state[imvNum] -> expand imv_state[imv_cont.k]
              expand(imv_state[imvNum],imv_state[imv_cont.k]);
              imv_state[imvNum]->m_valid_probe = _mm512_kand(imv_state[imvNum]->m_valid_probe, _mm512_knot(mask[VECTORSIZE - num]));
              // compress imv_state[imvNum]
              compress(imv_state[imvNum]);
              imv_state[imvNum]->m_valid_probe = imv_state[imvNum]->m_valid_probe >> (VECTORSIZE - num);
              imv_state[imv_cont.k]->m_valid_probe = mask[VECTORSIZE];
              imv_state[imv_cont.k]->stage = 0;
              ht_pos = (uint64_t *) &imv_state[imv_cont.k]->v_bucket_addrs;
              for (int i = 0; i < VECTORSIZE; ++i) {
                _mm_prefetch((char * )(ht_pos[i]), _MM_HINT_T0);
                _mm_prefetch(((char * )(ht_pos[i]) + 64), _MM_HINT_T0);
              }
            }
          }
        }
        if (found + VECTORSIZE >= batchSize) {
          return found;
        }
      }
        break;
    }
    ++imv_cont.k;
  }
  imv_cont.k = 100;
  return found;
}

pos_t Hashjoin::joinFullSIMD() {
  //std::cout <<"joinFullSIMD "<<"\n";

  size_t found = 0;
  int* probeKeys = (reinterpret_cast<int*>(((F2_Op *) (probeHash.ops[0].get()))->param1));
  uint64_t keyOff = ((EqualityCheck*) (keyEquality.ops[0].get()))->offset;
  __mmask8 m_match = 0, m_new_probes = -1;
  __m512i v_base_offset = _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0);
  __m512i v_offset = _mm512_set1_epi64(0), v_new_build_key, v_build_keys;
  __m512i v_base_offset_upper = _mm512_set1_epi64(cont.numProbes);
  __m512i v_seed = _mm512_set1_epi64(primitives::seed), v_build_key_off = _mm512_set1_epi64(keyOff);
  __m512i v_probe_hash = _mm512_set1_epi64(0), v_zero = _mm512_set1_epi64(0);
  __m256i v256_zero = _mm256_set1_epi32(0), v256_probe_keys, v256_build_keys;

  for (; cont.nextProbe < cont.numProbes || SIMDcon->m_valid_probe;) {
    /// step 1: load the offsets of probing tuples
    v_offset = _mm512_add_epi64(_mm512_set1_epi64(cont.nextProbe), v_base_offset);
    SIMDcon->v_probe_offset = _mm512_mask_expand_epi64(SIMDcon->v_probe_offset, _mm512_knot(SIMDcon->m_valid_probe), v_offset);
    // count the number of empty tuples
    m_new_probes = _mm512_knot(SIMDcon->m_valid_probe);
    cont.nextProbe = cont.nextProbe + _mm_popcnt_u32(m_new_probes);
    SIMDcon->m_valid_probe = _mm512_cmpgt_epu64_mask(v_base_offset_upper, SIMDcon->v_probe_offset);
    m_new_probes = _mm512_kand(m_new_probes, SIMDcon->m_valid_probe);
#if 1
    /// step 2: gather the probe keys
    v256_probe_keys = _mm512_mask_i64gather_epi32(v256_zero, m_new_probes, SIMDcon->v_probe_offset, (void* )probeKeys, 4);
    SIMDcon->v_probe_keys = _mm512_mask_blend_epi64(m_new_probes, SIMDcon->v_probe_keys, _mm512_cvtepi32_epi64(v256_probe_keys));
    /// step 3: compute the hash values of probe keys
    v_probe_hash = runtime::MurMurHash()((SIMDcon->v_probe_keys), (v_seed));
    /// step 4: find the addresses of corresponding buckets for new probes
    Vec8uM v_new_bucket_addrs = shared.ht.find_chain_tagged_sel((v_probe_hash), m_new_probes);
    // the addresses are null, then the corresponding probes are invalid
    SIMDcon->m_valid_probe = _mm512_kand(_mm512_kor(_mm512_knot(m_new_probes), v_new_bucket_addrs.mask), SIMDcon->m_valid_probe);
    SIMDcon->v_bucket_addrs = _mm512_mask_blend_epi64(v_new_bucket_addrs.mask, SIMDcon->v_bucket_addrs, v_new_bucket_addrs.vec);
    /// step 5: gather the all new build keys
    v256_build_keys = _mm512_mask_i64gather_epi32(v256_zero, SIMDcon->m_valid_probe, _mm512_add_epi64(SIMDcon->v_bucket_addrs, v_build_key_off), nullptr, 1);
    v_build_keys = _mm512_cvtepi32_epi64(v256_build_keys);
    /// step 6: compare the probe keys and build keys and write points
    m_match = _mm512_cmpeq_epi64_mask(SIMDcon->v_probe_keys, v_build_keys);
    m_match = _mm512_kand(m_match, SIMDcon->m_valid_probe);
    _mm512_mask_compressstoreu_epi64((buildMatches + found), m_match, SIMDcon->v_bucket_addrs);
    _mm256_mask_compressstoreu_epi32((probeMatches + found), m_match, _mm512_cvtepi64_epi32(SIMDcon->v_probe_offset));
#else
    // just write the pair(buildEntries, probeOffset), note to open the probeHashes and Equality
    // fetch the hashes of probing tuples
    Vec8u v_probe_hash = _mm512_mask_i64gather_epi64(v_zero,SIMDcon->m_valid_probe,SIMDcon->v_probe_offset,(void*)probeHashes,sizeof(runtime::Hashmap::hash_t));
    Vec8uM entries = shared.ht.find_chain_tagged(v_probe_hash);
    SIMDcon->v_bucket_addrs = _mm512_mask_blend_epi64(m_new_probes,SIMDcon->v_bucket_addrs,entries.vec);
    SIMDcon->m_valid_probe= SIMDcon->m_valid_probe&entries.mask;
    m_match=SIMDcon->m_valid_probe;
    _mm512_mask_compressstoreu_epi64((buildMatches + found), m_match,SIMDcon->v_bucket_addrs);
    _mm256_mask_compressstoreu_epi32((probeMatches + found), m_match, _mm512_cvtepi64_epi32(SIMDcon->v_probe_offset));
    // gather the addresses of building entries
#endif
    found += _mm_popcnt_u32(m_match);
    /// step 7: move to the next bucket nodes
    SIMDcon->v_bucket_addrs = _mm512_mask_i64gather_epi64(v_zero, SIMDcon->m_valid_probe, SIMDcon->v_bucket_addrs, nullptr, 1);
    SIMDcon->m_valid_probe = _mm512_kand(SIMDcon->m_valid_probe, _mm512_cmpneq_epi64_mask(SIMDcon->v_bucket_addrs, v_zero));
    if (found + VECTORSIZE >= batchSize) {
      return found;
    }
  }
  SIMDcon->m_valid_probe = 0;
  cont.nextProbe = cont.numProbes;
  return found;
}
pos_t Hashjoin::joinAMAC() {
  //std::cout <<"joinAMAC "<<"\n";

  size_t found = 0;
  // initialization
  if (amac_cont.k == 100) {
    for (int i = 0; i < stateNum; ++i) {
      amac_state[i].stage = 1;
    }
    amac_cont.done = 0;
    amac_cont.k = 0;
  }
  auto probeKeys = (reinterpret_cast<int*>(((F2_Op *) (probeHash.ops[0].get()))->param1));
  auto keyOff = ((EqualityCheck*) (keyEquality.ops[0].get()))->offset;
  while (amac_cont.done < stateNum) {
    amac_cont.k = (amac_cont.k >= stateNum) ? 0 : amac_cont.k;
    switch (amac_state[amac_cont.k].stage) {

      case 1: {
        if (cont.nextProbe >= cont.numProbes) {
          ++amac_cont.done;
          amac_state[amac_cont.k].stage = 3;
          //   std::cout<<"amac done one "<<cont.numProbes<<" , "<<cont.nextProbe<<std::endl;
          break;
        }
        cont.probeKey = probeKeys[cont.nextProbe];
        cont.probeHash = (runtime::MurMurHash()(cont.probeKey, primitives::seed));
        amac_state[amac_cont.k].tuple_id = cont.nextProbe;
        ++cont.nextProbe;
        amac_state[amac_cont.k].probeKey = cont.probeKey;
        amac_state[amac_cont.k].probeHash = cont.probeHash;
        shared.ht.PrefetchEntry(cont.probeHash);
        amac_state[amac_cont.k].stage = 2;
      }
        break;
      case 2: {
        amac_state[amac_cont.k].buildMatch = shared.ht.find_chain_tagged(amac_state[amac_cont.k].probeHash);
        if (nullptr == amac_state[amac_cont.k].buildMatch) {
          amac_state[amac_cont.k].stage = 1;
          --amac_cont.k;
        } else {
          _mm_prefetch((char * )(amac_state[amac_cont.k].buildMatch), _MM_HINT_T0);
          _mm_prefetch((char * )(amac_state[amac_cont.k].buildMatch) + 64, _MM_HINT_T0);
          amac_state[amac_cont.k].stage = 0;
        }
      }
        break;
      case 0: {
        auto entry = amac_state[amac_cont.k].buildMatch;
        if (nullptr == entry) {
          amac_state[amac_cont.k].stage = 1;
          --amac_cont.k;
          break;
        } else {
          auto buildkey = *((addBytes((reinterpret_cast<int*>(entry)), keyOff)));
          if ((buildkey == amac_state[amac_cont.k].probeKey)) {
            buildMatches[found] = entry;
            probeMatches[found++] = amac_state[amac_cont.k].tuple_id;
          }
          entry = entry->next;
          if (nullptr == entry) {
            amac_state[amac_cont.k].stage = 1;
            --amac_cont.k;
          } else {
            amac_state[amac_cont.k].buildMatch = entry;
            _mm_prefetch((char * )(entry), _MM_HINT_T0);
            _mm_prefetch((char * )(entry) + 64, _MM_HINT_T0);
          }
          if (found == batchSize) {
            return batchSize;
          }
        }
      }
        break;
    }
    ++amac_cont.k;
  }
  amac_cont.k = 100;
  return found;
}
// Note the way to get buildKey and probeKey: multi-joinkey or selective vectors
pos_t Hashjoin::joinRow() {
    //std::cout <<"joinRow "<<"\n";
  // std::cout<<"use join row "<<std::endl;
  size_t found = 0;
  do {
    for (auto entry = cont.buildMatch; entry != shared.ht.end(); entry = entry->next) {
      auto buildkey = *((addBytes((reinterpret_cast<int*>(entry)), ((EqualityCheck*) (keyEquality.ops[0].get()))->offset)));
      if ((buildkey == cont.probeKey)) {
//      if (entry->hash == cont.probeHash) {
        buildMatches[found] = entry;
        probeMatches[found++] = cont.nextProbe;
        if (found == batchSize) {
          // output buffers are full, save state for continuation
          cont.buildMatch = entry->next;
          if (cont.buildMatch == shared.ht.end()) {
            ++cont.nextProbe;
          }
          return batchSize;
        }
      }
    }
    if (cont.buildMatch != shared.ht.end())
      ++cont.nextProbe;

    if (cont.nextProbe < cont.numProbes) {
//      cont.probeHash = probeHashes[cont.nextProbe];
      cont.probeKey = (reinterpret_cast<int*>(((F2_Op *) (probeHash.ops[0].get()))->param1))[cont.nextProbe];
      cont.probeHash = runtime::MurMurHash()(cont.probeKey, primitives::seed);
      cont.buildMatch = shared.ht.find_chain_tagged(cont.probeHash);
      if (cont.buildMatch == shared.ht.end()) {
        ++cont.nextProbe;
      }
    } else {
      break;
    }
  } while (true);
  cont.buildMatch = shared.ht.end();
  cont.nextProbe = cont.numProbes;
  return found;
}
#define OLD 1
pos_t Hashjoin::joinAll() {
      //std::cout <<"joinAll "<<"\n";

   size_t found = 0;
   // perform continuation
   for (auto entry = cont.buildMatch; entry != shared.ht.end();
        entry = entry->next) {
#if OLD
      if (entry->hash == cont.probeHash) {
#else
        auto buildkey = *((addBytes((reinterpret_cast<int*>(entry)),((EqualityCheck*)(keyEquality.ops[0].get()))->offset)));
        if ((buildkey==cont.probeKey)) {
#endif
        buildMatches[found] = entry;
         probeMatches[found++] = cont.nextProbe;
         if (found == batchSize) {
            // output buffers are full, save state for continuation
            cont.buildMatch = entry->next;
            return batchSize;
         }
      }
   }
   if (cont.buildMatch != shared.ht.end()) cont.nextProbe++;
   for (size_t i = cont.nextProbe, end = cont.numProbes; i < end; ++i) {
#if OLD
     auto hash = probeHashes[i];
#else
     auto probeKey =(reinterpret_cast<int*>(((F2_Op *)(probeHash.ops[0].get()))->param1))[i];
     auto hash = runtime::MurMurHash()(probeKey,primitives::seed);
#endif
      for (auto entry = shared.ht.find_chain_tagged(hash);
           entry != shared.ht.end(); entry = entry->next) {
#if OLD
        if (entry->hash == hash) {
#else
          auto buildkey = *((addBytes((reinterpret_cast<int*>(entry)),((EqualityCheck*)(keyEquality.ops[0].get()))->offset)));
          if ((buildkey==probeKey)) {
               cont.probeKey = probeKey;
#endif
            buildMatches[found] = entry;
            probeMatches[found++] = i;
            if (found == batchSize && (entry->next || i + 1 < end)) {
               // output buffers are full, save state for continuation
               cont.buildMatch = entry->next;
               cont.probeHash = hash;
               cont.nextProbe = i;
               return batchSize;
            }
         }
      }
   }
   cont.buildMatch = shared.ht.end();
   cont.nextProbe = cont.numProbes;
   return found;
}

pos_t Hashjoin::joinAllParallel() {
   //std::cout <<"joinAllParallel "<<"\n";

   size_t found = 0;
   auto followup = contCon.followup;
   auto followupWrite = contCon.followupWrite;

   if (followup == followupWrite) {
      for (size_t i = 0, end = cont.numProbes; i < end; ++i) {
         auto hash = probeHashes[i];
         auto entry = shared.ht.find_chain_tagged(hash);
         if (entry != shared.ht.end()) {
            if (entry->hash == hash) {
               buildMatches[found] = entry;
               probeMatches[found] = i;
               found += 1;
            }
            if (entry->next != shared.ht.end()) {
               followupIds[followupWrite] = i;
               followupEntries[followupWrite] = entry->next;
               followupWrite += 1;
            }
         }
      }
   }

   followupWrite %= followupBufferSize;

   while (followup != followupWrite) {
      auto remainingSpace = batchSize - found;
      auto nrFollowups = followup <= followupWrite
                             ? followupWrite - followup
                             : followupBufferSize - (followup - followupWrite);
      // std::cout << "nrFollowups: " << nrFollowups << "\n";
      auto fittingElements = std::min((size_t)nrFollowups, remainingSpace);
      for (size_t j = 0; j < fittingElements; ++j) {
         size_t i = followupIds[followup];
         auto entry = followupEntries[followup];
         // followup = (followup + 1) % followupBufferSize;
         followup = (followup + 1);
         if (followup == followupBufferSize) followup = 0;
         auto hash = probeHashes[i];
         if (entry->hash == hash) {
            buildMatches[found] = entry;
            probeMatches[found++] = i;
         }
         if (entry->next != shared.ht.end()) {
            followupIds[followupWrite] = i;
            followupEntries[followupWrite] = entry->next;
            followupWrite = (followupWrite + 1) % followupBufferSize;
         }
      }
      if (fittingElements < nrFollowups) {
         // continuation
         contCon.followupWrite = followupWrite;
         contCon.followup = followup;
         return found;
      }
   }
   cont.nextProbe = cont.numProbes;
   contCon.followup = 0;
   contCon.followupWrite = 0;
   return found;
}

pos_t Hashjoin::joinAllSIMD() {
     //std::cout <<"joinAllSIMD "<<"\n";

   size_t found = 0;
   auto followup = contCon.followup;
   auto followupWrite = contCon.followupWrite;

   if (followup == followupWrite) {

#ifdef __AVX512F__ // if AVX 512 available, use it!
#if HASH_SIZE == 32
      size_t rest = cont.numProbes % 8;
      auto ids =
          _mm512_set_epi32(0, 0, 0, 0, 0, 0, 0, 0, 7, 6, 5, 4, 3, 2, 1, 0);
      for (size_t i = 0, end = cont.numProbes - rest; i < end; i += 8) {

         // load hashes
         // auto hash = probeHashes[i];
         // Vec8u hashes(probeHashes + i);
         auto hashDense = _mm256_loadu_si256((const __m256i*)(probeHashes + i));
         Vec8u hashes = _mm512_cvtepu32_epi64(hashDense);
         // find entry pointers in ht
         Vec8uM entries = shared.ht.find_chain_tagged(hashes);
         // load entry hashes
         auto entryHashes = _mm512_mask_i64gather_epi32(
             hashDense, entries.mask,
             entries.vec +
                 Vec8u(offsetof(decltype(shared.ht)::EntryHeader, hash)),
             nullptr, 1);
         {
            // Check if hashes match
            __mmask8 hashesEq = _mm512_mask_cmpeq_epi32_mask(
                entries.mask, _mm512_castsi256_si512(entryHashes),
                _mm512_castsi256_si512(hashDense));
            // write pointers
            _mm512_mask_compressstoreu_epi64(buildMatches + found, hashesEq,
                                             entries.vec);
            static_assert(sizeof(pos_t) == 4,
                          "SIMD join assumes sizeof(pos_t) is 4"); // change the
                                                                   // types for
                                                                   // probeSels
                                                                   // if this
                                                                   // fails
            // write selection
            _mm512_mask_compressstoreu_epi32(probeMatches + found, hashesEq,
                                             ids);
            found += __builtin_popcount(hashesEq);
         }

         {
            // write continuations
            static_assert(offsetof(decltype(shared.ht)::EntryHeader, next) == 0,
                          "Hash is expected to be in first position");
            Vec8u nextPtrs = _mm512_mask_i64gather_epi64(
                entries.vec, entries.mask, entries.vec, nullptr, 1);
            __mmask8 hasNext = _mm512_mask_cmpneq_epi64_mask(
                entries.mask, nextPtrs, Vec8u(uint64_t(shared.ht.end())));
            if (hasNext) {
               // write pointers
               _mm512_mask_compressstoreu_epi64(followupEntries + followupWrite,
                                                hasNext, nextPtrs);
               static_assert(
                   sizeof(pos_t) == 4,
                   "SIMD join assumes sizeof(pos_t) is 4"); // change the types
                                                            // for probeSels if
                                                            // this fails
               // write selection
               _mm512_mask_compressstoreu_epi32(followupIds + followupWrite,
                                                hasNext, ids);
               followupWrite += __builtin_popcount(hasNext);
            }
            ids = _mm512_add_epi32(ids, _mm512_set1_epi32(8));
         }
      }
#else
      size_t rest = cont.numProbes % 8;
      // auto ids = _mm256_set_epi32(7,6,5,4,3,2,1,0);
      auto ids =
          _mm512_set_epi32(0, 0, 0, 0, 0, 0, 0, 0, 7, 6, 5, 4, 3, 2, 1, 0);
      for (size_t i = 0, end = cont.numProbes - rest; i < end; i += 8) {

         // load hashes
         // auto hash = probeHashes[i];
         Vec8u hashes(probeHashes + i);
         // find entry pointers in ht
         Vec8uM entries = shared.ht.find_chain_tagged(hashes);
         // load entry hashes
         Vec8u entryHashes = _mm512_mask_i64gather_epi64(
             entries.vec, entries.mask,
             entries.vec +
                 Vec8u(offsetof(decltype(shared.ht)::EntryHeader, hash)),
             nullptr, 1);
         {
            // Check if hashes match
            __mmask8 hashesEq =
                _mm512_mask_cmpeq_epi64_mask(entries.mask, entryHashes, hashes);
            // write pointers
            _mm512_mask_compressstoreu_epi64(buildMatches + found, hashesEq,
                                             entries.vec);
            static_assert(sizeof(pos_t) == 4,
                          "SIMD join assumes sizeof(pos_t) is 4"); // change the
                                                                   // types for
                                                                   // probeSels
                                                                   // if this
                                                                   // fails
            // write selection
            _mm512_mask_compressstoreu_epi32(probeMatches + found, hashesEq,
                                             ids);
            found += __builtin_popcount(hashesEq);
         }

         {
            // write continuations
            static_assert(offsetof(decltype(shared.ht)::EntryHeader, next) == 0,
                          "Hash is expected to be in first position");
            Vec8u nextPtrs = _mm512_mask_i64gather_epi64(
                entries.vec, entries.mask, entries.vec, nullptr, 1);
            __mmask8 hasNext = _mm512_mask_cmpneq_epi64_mask(
                entries.mask, nextPtrs, Vec8u(uint64_t(shared.ht.end())));
            if (hasNext) {
               // write pointers
               _mm512_mask_compressstoreu_epi64(followupEntries + followupWrite,
                                                hasNext, nextPtrs);
               static_assert(
                   sizeof(pos_t) == 4,
                   "SIMD join assumes sizeof(pos_t) is 4"); // change the types
                                                            // for probeSels if
                                                            // this fails
               // write selection
               _mm512_mask_compressstoreu_epi32(followupIds + followupWrite,
                                                hasNext, ids);
               followupWrite += __builtin_popcount(hasNext);
            }
            ids = _mm512_add_epi32(ids, _mm512_set1_epi32(8));
         }
      }
#endif // hash size
#else
      const size_t rest = cont.numProbes;
#endif
      for (size_t i = cont.numProbes - rest, end = cont.numProbes; i < end;
           ++i) {
         auto hash = probeHashes[i];
         auto entry = shared.ht.find_chain_tagged(hash);
         if (entry != shared.ht.end()) {
            if (entry->hash == hash) {
               buildMatches[found] = entry;
               probeMatches[found] = i;
               found += 1;
            }
            if (entry->next != shared.ht.end()) {
               followupIds[followupWrite] = i;
               followupEntries[followupWrite] = entry->next;
               followupWrite += 1;
            }
         }
      }
   }

   followupWrite %= followupBufferSize;

   while (followup != followupWrite) {
      auto remainingSpace = batchSize - found;
      auto nrFollowups = followup <= followupWrite
                             ? followupWrite - followup
                             : followupBufferSize - (followup - followupWrite);
      auto fittingElements = std::min((size_t)nrFollowups, remainingSpace);
      for (size_t j = 0; j < fittingElements; ++j) {
         size_t i = followupIds[followup];
         auto entry = followupEntries[followup];
         followup = (followup + 1);
         if (followup == followupBufferSize) followup = 0;
         auto hash = probeHashes[i];
         if (entry->hash == hash) {
            buildMatches[found] = entry;
            probeMatches[found++] = i;
         }
         if (entry->next != shared.ht.end()) {
            followupIds[followupWrite] = i;
            followupEntries[followupWrite] = entry->next;
            followupWrite = (followupWrite + 1) % followupBufferSize;
         }
      }
      if (fittingElements < nrFollowups) {
         // continuation
         contCon.followupWrite = followupWrite;
         contCon.followup = followup;
         return found;
      }
   }
   cont.nextProbe = cont.numProbes;
   contCon.followup = 0;
   contCon.followupWrite = 0;
   return found;
}

pos_t Hashjoin::joinSel() {
       //std::cout <<"joinSel "<<"\n";

   size_t found = 0;
   // perform continuation
   for (auto entry = cont.buildMatch; entry != shared.ht.end();
        entry = entry->next) {
      if (entry->hash == cont.probeHash) {
         buildMatches[found] = entry;
         probeMatches[found++] = probeSel[cont.nextProbe];
         if (found == batchSize) {
            // output buffers are full, save state for continuation
            cont.buildMatch = entry->next;
            return batchSize;
         }
      }
   }
   if (cont.buildMatch != shared.ht.end()) cont.nextProbe++;
   for (size_t i = cont.nextProbe, end = cont.numProbes; i < end; ++i) {
      auto hash = probeHashes[i];
      for (auto entry = shared.ht.find_chain_tagged(hash);
           entry != shared.ht.end(); entry = entry->next) {
         if (entry->hash == hash) {
            buildMatches[found] = entry;
            probeMatches[found++] = probeSel[i];
            if (found == batchSize && (entry->next || i + 1 < end)) {
               // output buffers are full, save state for continuation
               cont.buildMatch = entry->next;
               cont.probeHash = hash;
               cont.nextProbe = i;
               return batchSize;
            }
         }
      }
   }
   cont.buildMatch = shared.ht.end();
   cont.nextProbe = cont.numProbes;
   return found;
}

pos_t Hashjoin::joinSelParallel() {
    //std::cout <<"joinSelParallel "<<"\n";

   size_t found = 0;
   auto followup = contCon.followup;
   auto followupWrite = contCon.followupWrite;

   if (followup == followupWrite) {
      for (size_t i = 0, end = cont.numProbes; i < end; ++i) {
         auto hash = probeHashes[i];
         auto entry = shared.ht.find_chain_tagged(hash);
         if (entry != shared.ht.end()) {
            if (entry->hash == hash) {
               buildMatches[found] = entry;
               probeMatches[found] = probeSel[i];
               found += 1;
            }
            if (entry->next != shared.ht.end()) {
               followupIds[followupWrite] = i;
               followupEntries[followupWrite] = entry->next;
               followupWrite += 1;
            }
         }
      }
   }

   followupWrite %= followupBufferSize;

   while (followup != followupWrite) {
      auto remainingSpace = batchSize - found;
      auto nrFollowups = followup <= followupWrite
                             ? followupWrite - followup
                             : followupBufferSize - (followup - followupWrite);
      auto fittingElements = std::min((size_t)nrFollowups, remainingSpace);
      for (size_t j = 0; j < fittingElements; ++j) {
         size_t i = followupIds[followup];
         auto entry = followupEntries[followup];
         followup = (followup + 1) % followupBufferSize;
         auto hash = probeHashes[i];
         if (entry->hash == hash) {
            buildMatches[found] = entry;
            probeMatches[found++] = probeSel[i];
         }
         if (entry->next != shared.ht.end()) {
            followupIds[followupWrite] = i;
            followupEntries[followupWrite] = entry->next;
            followupWrite = (followupWrite + 1) % followupBufferSize;
         }
      }
      if (fittingElements < nrFollowups) {
         // continuation
         contCon.followupWrite = followupWrite;
         contCon.followup = followup;
         return found;
      }
   }
   cont.nextProbe = cont.numProbes;
   contCon.followup = 0;
   contCon.followupWrite = 0;
   return found;
}

pos_t Hashjoin::joinSelSIMD() {
      //std::cout <<"joinSelSIMD "<<"\n";

   size_t found = 0;
   auto followup = contCon.followup;
   auto followupWrite = contCon.followupWrite;

   if (followup == followupWrite) {

#ifdef __AVX512F__ // if AVX 512 available, use it!
#if HASH_SIZE == 32
      size_t rest = cont.numProbes % 8;
      auto ids =
          _mm512_set_epi32(0, 0, 0, 0, 0, 0, 0, 0, 7, 6, 5, 4, 3, 2, 1, 0);
      for (size_t i = 0, end = cont.numProbes - rest; i < end; i += 8) {

         // load hashes
         // auto hash = probeHashes[i];
         auto hashDense = _mm256_loadu_si256((const __m256i*)(probeHashes + i));
         Vec8u hashes = _mm512_cvtepu32_epi64(hashDense);
         // Vec8u hashes(probeHashes + i);
         // find entry pointers in ht
         Vec8uM entries = shared.ht.find_chain_tagged(hashes);
         // load entry hashes
         Vec8u hashPtrs =
             entries.vec +
             Vec8u(offsetof(decltype(shared.ht)::EntryHeader, hash));
         auto entryHashes = _mm512_mask_i64gather_epi32(hashDense, entries.mask,
                                                        hashPtrs, nullptr, 1);
         {
            // Check if hashes match
            __mmask8 hashesEq = _mm512_mask_cmpeq_epi32_mask(
                entries.mask, _mm512_castsi256_si512(entryHashes),
                _mm512_castsi256_si512(hashDense));
            // write pointers
            _mm512_mask_compressstoreu_epi64(buildMatches + found, hashesEq,
                                             entries.vec);
            static_assert(sizeof(pos_t) == 4,
                          "SIMD join assumes sizeof(pos_t) is 4"); // change the
                                                                   // types for
                                                                   // probeSels
                                                                   // if this
                                                                   // fails
            // write selection
            __m512i probeSels = _mm512_loadu_si512(probeSel + i);
            _mm512_mask_compressstoreu_epi32(probeMatches + found, hashesEq,
                                             probeSels);
            found += __builtin_popcount(hashesEq);
         }

         {
            // write continuations
            static_assert(offsetof(decltype(shared.ht)::EntryHeader, next) == 0,
                          "Hash is expected to be in first position");
            Vec8u nextPtrs = _mm512_mask_i64gather_epi64(
                entries.vec, entries.mask, entries.vec, nullptr, 1);
            __mmask8 hasNext = _mm512_mask_cmpneq_epi64_mask(
                entries.mask, nextPtrs, Vec8u(uint64_t(shared.ht.end())));
            if (hasNext) {
               // write pointers
               _mm512_mask_compressstoreu_epi64(followupEntries + followupWrite,
                                                hasNext, nextPtrs);
               static_assert(
                   sizeof(pos_t) == 4,
                   "SIMD join assumes sizeof(pos_t) is 4"); // change the types
                                                            // for probeSels if
                                                            // this fails
               // write selection
               _mm512_mask_compressstoreu_epi32(followupIds + followupWrite,
                                                hasNext, ids);
               followupWrite += __builtin_popcount(hasNext);
            }
            ids = _mm512_add_epi32(ids, _mm512_set1_epi32(8));
         }
      }
#else
      size_t rest = cont.numProbes % 8;
      auto ids =
          _mm512_set_epi32(0, 0, 0, 0, 0, 0, 0, 0, 7, 6, 5, 4, 3, 2, 1, 0);
      for (size_t i = 0, end = cont.numProbes - rest; i < end; i += 8) {

         // load hashes
         // auto hash = probeHashes[i];
         Vec8u hashes(probeHashes + i);
         // find entry pointers in ht
         Vec8uM entries = shared.ht.find_chain_tagged(hashes);
         // load entry hashes
         Vec8u hashPtrs =
             entries.vec +
             Vec8u(offsetof(decltype(shared.ht)::EntryHeader, hash));
         Vec8u entryHashes = _mm512_mask_i64gather_epi64(hashPtrs, entries.mask,
                                                         hashPtrs, nullptr, 1);
         {
            // Check if hashes match
            __mmask8 hashesEq =
                _mm512_mask_cmpeq_epi64_mask(entries.mask, entryHashes, hashes);
            // write pointers
            _mm512_mask_compressstoreu_epi64(buildMatches + found, hashesEq,
                                             entries.vec);
            static_assert(sizeof(pos_t) == 4,
                          "SIMD join assumes sizeof(pos_t) is 4"); // change the
                                                                   // types for
                                                                   // probeSels
                                                                   // if this
                                                                   // fails
            // write selection
            __m512i probeSels = _mm512_loadu_si512(probeSel + i);
            _mm512_mask_compressstoreu_epi32(probeMatches + found, hashesEq,
                                             probeSels);
            found += __builtin_popcount(hashesEq);
         }

         {
            // write continuations
            static_assert(offsetof(decltype(shared.ht)::EntryHeader, next) == 0,
                          "Hash is expected to be in first position");
            Vec8u nextPtrs = _mm512_mask_i64gather_epi64(
                entries.vec, entries.mask, entries.vec, nullptr, 1);
            __mmask8 hasNext = _mm512_mask_cmpneq_epi64_mask(
                entries.mask, nextPtrs, Vec8u(uint64_t(shared.ht.end())));
            if (hasNext) {
               // write pointers
               _mm512_mask_compressstoreu_epi64(followupEntries + followupWrite,
                                                hasNext, nextPtrs);
               static_assert(
                   sizeof(pos_t) == 4,
                   "SIMD join assumes sizeof(pos_t) is 4"); // change the types
                                                            // for probeSels if
                                                            // this fails
               // write selection
               _mm512_mask_compressstoreu_epi32(followupIds + followupWrite,
                                                hasNext, ids);
               followupWrite += __builtin_popcount(hasNext);
            }
            ids = _mm512_add_epi32(ids, _mm512_set1_epi32(8));
         }
      }
#endif // hash size
#else
      const size_t rest = cont.numProbes;
#endif
      for (size_t i = cont.numProbes - rest, end = cont.numProbes; i < end;
           ++i) {
         auto hash = probeHashes[i];
         auto entry = shared.ht.find_chain_tagged(hash);
         if (entry != shared.ht.end()) {
            if (entry->hash == hash) {
               buildMatches[found] = entry;
               probeMatches[found] = probeSel[i];
               found += 1;
            }
            if (entry->next != shared.ht.end()) {
               followupIds[followupWrite] = i;
               followupEntries[followupWrite] = entry->next;
               followupWrite += 1;
            }
         }
      }
   }

   followupWrite %= followupBufferSize;

   while (followup != followupWrite) {
      auto remainingSpace = batchSize - found;
      auto nrFollowups = followup <= followupWrite
                             ? followupWrite - followup
                             : followupBufferSize - (followup - followupWrite);
      auto fittingElements = std::min((size_t)nrFollowups, remainingSpace);
      for (size_t j = 0; j < fittingElements; ++j) {
         size_t i = followupIds[followup];
         auto entry = followupEntries[followup];
         followup = (followup + 1) % followupBufferSize;
         auto hash = probeHashes[i];
         if (entry->hash == hash) {
            buildMatches[found] = entry;
            probeMatches[found++] = probeSel[i];
         }
         if (entry->next != shared.ht.end()) {
            followupIds[followupWrite] = i;
            followupEntries[followupWrite] = entry->next;
            followupWrite = (followupWrite + 1) % followupBufferSize;
         }
      }
      if (fittingElements < nrFollowups) {
         // continuation
         contCon.followupWrite = followupWrite;
         contCon.followup = followup;
         return found;
      }
   }
   cont.nextProbe = cont.numProbes;
   contCon.followup = 0;
   contCon.followupWrite = 0;
   return found;
}

template <typename T, typename HT>
void INTERPRET_SEPARATE insertAllEntries(T& allocations, HT& ht,
                                         size_t ht_entry_size) {
   for (auto& block : allocations) {
      auto start =
          reinterpret_cast<runtime::Hashmap::EntryHeader*>(block.first);
      ht.insertAll_tagged(start, block.second, ht_entry_size);
   }
}

pos_t Hashjoin::joinBoncz() {
   //std::cout <<"joinBoncz "<<"\n";

   size_t followupWrite = contCon.followupWrite;
   size_t found = 0;
   if (followupWrite == 0)
      for (size_t i = 0, end = cont.numProbes; i < end; ++i) {
         auto hash = probeHashes[i];
         auto entry = shared.ht.find_chain_tagged(hash);
         if (entry != shared.ht.end()) {
            followupIds[followupWrite] = i;
            followupEntries[followupWrite] = entry;
            followupWrite += 1;
         }
      }

   while (followupWrite > 0) {
      size_t e = followupWrite;
      followupWrite = 0;
      for (size_t j = 0; j < e; j++) {
         auto i = followupIds[j];
         auto entry = followupEntries[j];
         // auto hash = probeHashes[i];
         // if (entry->hash == hash) {
         buildMatches[found] = entry;
         probeMatches[found++] = i;
         // }
         if (entry->next) {
            followupIds[followupWrite] = i;
            followupEntries[followupWrite++] = entry->next;
         }
      }
      if (followupWrite == 0) {
         cont.nextProbe = cont.numProbes;
         contCon.followupWrite = followupWrite;
         return found;
      } else if(found + followupWrite >= batchSize){
         contCon.followupWrite = followupWrite;
         assert(found);
         return found;
     }
   }
   cont.nextProbe = cont.numProbes;
   contCon.followupWrite = followupWrite;
   return 0;
}
struct __attribute__((aligned(64))) BuildSIMDState {
  __m512i v_entry_addr;
  __m512i v_hash_value;
  __m512i v_build_key;
  __mmask8 m_valid;
  uint8_t valid_size, stage,mask[VECTORSIZE+1];
  BuildSIMDState()
      : v_entry_addr(_mm512_set1_epi64(0)),
        v_hash_value(_mm512_set1_epi64(0)),
        valid_size(VECTORSIZE),
        stage(1),
        m_valid(0){
    for (int i = 0; i <= VECTORSIZE; ++i) {
      mask[i] = (1 << i) - 1;
    }
  }
  inline void reset() {
    v_entry_addr = _mm512_set1_epi64(0);
    v_hash_value = _mm512_set1_epi64(0);
    valid_size = VECTORSIZE;
    stage = 1;
    m_valid=0;
  }
  void* operator new(size_t size) {
    return memalign(64, size);
  }
  void operator delete(void* mem) {
    return free(mem);
  }
  void* operator new[](size_t size) {
    return memalign(64, size);
  }
  void operator delete[](void* mem) {
    return free(mem);
  }
  inline void compress() {
    v_entry_addr = _mm512_maskz_compress_epi64(m_valid, v_entry_addr);
    v_hash_value = _mm512_maskz_compress_epi64(m_valid, v_hash_value);
    v_build_key = _mm512_maskz_compress_epi64(m_valid, v_build_key);
  }
  inline void expand(BuildSIMDState& src_state) {
    v_entry_addr = _mm512_mask_expand_epi64(v_entry_addr, _mm512_knot(m_valid), src_state.v_entry_addr);
    v_hash_value = _mm512_mask_expand_epi64(v_hash_value, _mm512_knot(m_valid), src_state.v_hash_value);
    v_build_key = _mm512_mask_expand_epi64(v_build_key, _mm512_knot(m_valid), src_state.v_build_key);


  }
  inline void compact(BuildSIMDState& RVS, uint8_t done, uint8_t imvNum,uint8_t& k,uint8_t next_stage,uint8_t src_stage){
   auto num = _mm_popcnt_u32(m_valid);
    if (num == VECTORSIZE || done >= imvNum) {
      stage = next_stage;
    } else {
       auto num_temp = _mm_popcnt_u32(RVS.m_valid);
        if (num + num_temp < VECTORSIZE) {
          // compress imv_state[k]
          compress();
          // expand imv_state[k] -> imv_state[imvNum1]
          RVS.expand(*this);
          RVS.m_valid = mask[num + num_temp];
          m_valid = 0;
          stage = src_stage;
          RVS.stage = next_stage;
          --k;
        } else {
          // expand imv_state[imvNum1] -> expand imv_state[k]
          expand(RVS);
          RVS.m_valid = _mm512_kand(RVS.m_valid, _mm512_knot(mask[VECTORSIZE - num]));
          // compress imv_state[imvNum]
          RVS.compress();
          RVS.m_valid = RVS.m_valid >> (VECTORSIZE - num);
          m_valid = mask[VECTORSIZE];
          RVS.stage = next_stage;
          stage = next_stage;
        }
    }
  }
};
static int stateNumSIMD=vectorwise::Hashjoin::imvNum;
size_t HashBuild(size_t begin, size_t end, uint32_t* buildKey, runtime::Hashmap* hash_table, runtime::Allocator*allo, int entry_size, uint32_t* pos_buff) {
  size_t found = 0, cur = begin;
  uint8_t valid_size = VECTORSIZE, done = 0, k = 0;

  int build_key_off = sizeof(runtime::Hashmap::EntryHeader);
  __m512i v_build_key, v_offset, v_base_entry_off, v_key_off = _mm512_set1_epi64(build_key_off), v_build_hash_mask, v_zero = _mm512_set1_epi64(0), v_all_ones = _mm512_set1_epi64(
      -1), v_conflict, v_base_offset = _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0), v_seed = _mm512_set1_epi64(primitives::seed), v_offset_upper = _mm512_set1_epi64(end);

  __mmask8 m_no_conflict, m_rest;
  __m256i v256_zero = _mm256_set1_epi32(0), v256_build_key;
  v_base_entry_off = _mm512_mullo_epi64(v_base_offset, _mm512_set1_epi64(entry_size));
  uint64_t* hash_value = nullptr;
  BuildSIMDState state[stateNumSIMD];
  while (done < stateNumSIMD) {
    k = (k >= stateNumSIMD) ? 0 : k;
    if (cur >= end) {
      if (state[k].m_valid == 0 && state[k].stage != 3) {
        ++done;
        state[k].stage = 3;
        ++k;
        continue;
      }
    }
    switch (state[k].stage) {
      case 1: {
        /// step 1: gather build keys (using gather to compilate with loading discontinuous values)
        v_offset = _mm512_add_epi64(_mm512_set1_epi64(cur), v_base_offset);
        state[k].m_valid = _mm512_cmpgt_epu64_mask(v_offset_upper, v_offset);
        if (pos_buff) {
          v_offset = _mm512_cvtepi32_epi64(_mm512_mask_i64gather_epi32(v256_zero, state[k].m_valid, v_offset, (const int* )pos_buff, 4));
        }
        v256_build_key = _mm512_mask_i64gather_epi32(v256_zero, state[k].m_valid, v_offset, (const int* )buildKey, 4);
        v_build_key = _mm512_cvtepi32_epi64(v256_build_key);
        cur += VECTORSIZE;

        /// step 2: allocate new entries
        runtime::Hashmap::EntryHeader* ptr = (runtime::Hashmap::EntryHeader*) allo->allocate(VECTORSIZE * entry_size);

        /// step 3: write build keys to new entries
        state[k].v_entry_addr = _mm512_add_epi64(_mm512_set1_epi64((uint64_t) ptr), v_base_entry_off);
        _mm512_i64scatter_epi64(0, _mm512_add_epi64(state[k].v_entry_addr, v_key_off), v_build_key, 1);

        /// step 4: hashing the build keys (note the hash value cannot be used to directly fetch buckets)
        state[k].v_hash_value = runtime::MurMurHash()(v_build_key, v_seed);
        state[k].stage = 0;

        hash_table->prefetchEntry(state[k].v_hash_value);
      }
        break;
      case 0: {
        /// scalar codes due to writhe conflicts among multi-threads
        hash_table->insert_tagged_sel((Vec8u*) (&state[k].v_entry_addr), (Vec8u*) (&state[k].v_hash_value), state[k].m_valid);
        found += _mm_popcnt_u32(state[k].m_valid);
        state[k].stage = 1;
        state[k].m_valid = 0;
      }
        break;
    }
    ++k;
  }

  return found;
}

size_t Hashjoin::next() {
   using runtime::Hashmap;
   // --- build
   if (!consumed) {
      size_t found = 0;

    if (join == &vectorwise::Hashjoin::joinSelIMV) {
      barrier([&]() {
        if (!shared.sizeIsSet.load()) {
          shared.sizeIsSet.store(false);
          if (shared.ht.capacity < ht_size) {
            shared.ht.setSize(ht_size);
//            std::cout<<"set size "<<ht_size<<std::endl;
          }
        }

      });

      for (auto n = left->next(); n != EndOfStream; n = left->next()) {
        found += n;
        uint32_t* buildKeys = (reinterpret_cast<int*>(((F3_Op *) (buildHash.ops[0].get()))->param2));
        pos_t * sel = (reinterpret_cast<int*>(((F3_Op *) (buildHash.ops[0].get()))->outputSelectionV));
        HashBuild(0, n, buildKeys, &shared.ht, &runtime::this_worker->allocator, ht_entry_size, sel);
      }
      shared.found.fetch_add(found);
      auto globalFound = shared.found.load();
      if (globalFound == 0) {
        consumed = true;
        return EndOfStream;
      }
    } else {
      // --- build phase 1: materialize ht entries
      for (auto n = left->next(); n != EndOfStream; n = left->next()) {
        found += n;
        // build hashes
        buildHash.evaluate(n);
        // scatter hash, keys and values into ht entries
        auto alloc = runtime::this_worker->allocator.allocate(n * ht_entry_size);
        if (!alloc)
          throw std::runtime_error("malloc failed");
        allocations.push_back(std::make_pair(alloc, n));
        scatterStart = reinterpret_cast<decltype(scatterStart)>(alloc);
        buildScatter.evaluate(n);
      }

      // --- build phase 2: insert ht entries
      shared.found.fetch_add(found);
      barrier([&]() {
        auto globalFound = shared.found.load();
        if (globalFound) shared.ht.setSize(globalFound);
      });
      auto globalFound = shared.found.load();
      if (globalFound == 0) {
        consumed = true;
        return EndOfStream;
      }
      insertAllEntries(allocations, shared.ht, ht_entry_size);
    }

    consumed = true;
    barrier();  // wait for all threads to finish build phase
    runtime::this_worker->log_time("build");
//   if(!shared.printed.load()) {
//     shared.printed.store(true);
//   shared.ht.printStaTag();
//   }
  }
   // --- lookup
   if(join == &vectorwise::Hashjoin::joinRow) {
     while (true) {
         if (cont.nextProbe >= cont.numProbes) {
           cont.numProbes = right->next();
            cont.nextProbe = 0;
            if (cont.numProbes == EndOfStream) return EndOfStream;
         }
         // create join pair vectors with matching hashes (Entry*, pos), where
         // Entry* is for the build side, pos a selection index to the right side
         auto n = (this->*join)();
         if (n == 0) continue;
         // materialize build side
         buildGather.evaluate(n);
         return n;
      }
   }else if(join == &vectorwise::Hashjoin::joinAMAC) {
     while (true) {
       if(amac_cont.k==100) {
          cont.numProbes = right->next();
           cont.nextProbe = 0;
           if (cont.numProbes == EndOfStream) return EndOfStream;
        }
        // create join pair vectors with matching hashes (Entry*, pos), where
        // Entry* is for the build side, pos a selection index to the right side
        auto n = (this->*join)();
        if (n == 0) continue;
        // materialize build side
        buildGather.evaluate(n);
        return n;
     }
   }else if(join == &vectorwise::Hashjoin::joinFullSIMD){
     while (true) {
        if (cont.nextProbe >= cont.numProbes && SIMDcon->m_valid_probe==0) {
           cont.numProbes = right->next();
           cont.nextProbe = 0;
           SIMDcon->reset();
           if (cont.numProbes == EndOfStream) return EndOfStream;
     //      probeHash.evaluate(cont.numProbes);
        }
        // create join pair vectors with matching hashes (Entry*, pos), where
        // Entry* is for the build side, pos a selection index to the right side
        auto n = (this->*join)();
        // check key equality and remove non equal keys from join result
     //   n = keyEquality.evaluate(n);
        if (n == 0) continue;
        // materialize build side
        buildGather.evaluate(n);
        return n;
     }

  }else if(join == &vectorwise::Hashjoin::joinSIMDAMAC || join ==&vectorwise::Hashjoin::joinIMV || join ==&vectorwise::Hashjoin::joinSelIMV){
    while (true) {
       if (imv_cont.k==100) {
          cont.numProbes = right->next();
          cont.nextProbe = 0;
          if (cont.numProbes == EndOfStream) return EndOfStream;
    //      probeHash.evaluate(cont.numProbes);
       }
       // create join pair vectors with matching hashes (Entry*, pos), where
       // Entry* is for the build side, pos a selection index to the right side
       auto n = (this->*join)();
       // check key equality and remove non equal keys from join result
    //   n = keyEquality.evaluate(n);
       if (n == 0) continue;
       // materialize build side
       buildGather.evaluate(n);
       return n;
    }

 }else {
    while (true) {
       if (cont.nextProbe >= cont.numProbes) {
         cont.numProbes = right->next();
          cont.nextProbe = 0;
          if (cont.numProbes == EndOfStream) return EndOfStream;
          probeHash.evaluate(cont.numProbes);
       }
       // create join pair vectors with matching hashes (Entry*, pos), where
       // Entry* is for the build side, pos a selection index to the right side
       auto n = (this->*join)();
       // check key equality and remove non equal keys from join result
      n = keyEquality.evaluate(n);
       if (n == 0) continue;
       // materialize build side
       buildGather.evaluate(n);
       return n;
    }
 }
}

Hashjoin::Hashjoin(Shared& sm) : shared(sm) { SIMDcon = new SIMDContinuation();

  for(int i=0;i<=imvNum;++i)
  imv_state[i] = new IMVState();
}

Hashjoin::~Hashjoin() {
   // for (auto& block : allocations) free(block.first);
  if(SIMDcon) {
    delete SIMDcon;
    SIMDcon = nullptr;
  }
  for(int i=0;i<=imvNum;++i){
    if(imv_state[i]) {
      delete imv_state[i];
      imv_state[i] = nullptr;
      }
   }
}

HashGroup::HashGroup(Shared& s)
    : shared(s), preAggregation(*this), globalAggregation(*this) {
   maxFill = ht.setSize(initialMapSize);
}
HashGroup::~HashGroup() {
   // for (auto& alloc : preAggregation.allocations) free(alloc.first);
   // for (auto& alloc : globalAggregation.allocations) free(alloc.first);
}

pos_t HashGroup::findGroupsFromPartition(void* data, size_t n) {
   globalAggregation.groupHashes = reinterpret_cast<hash_t*>(data);
   return globalAggregation.findGroups(n, ht);
}

size_t HashGroup::next() {
   using header_t = decltype(ht)::EntryHeader;
   if (!cont.consumed) {
      /// ------ phase 1: local preaggregation
      /// aggregate all incoming tuples into local hashtable
      size_t groups = 0;
      auto& spill = shared.spillStorage.local();
      auto entry_size = preAggregation.ht_entry_size;

      auto flushAndClear = [&]() INTERPRET_SEPARATE {
         assert(offsetof(header_t, next) + sizeof(header_t::next) ==
                offsetof(header_t, hash));
         // flush ht entries into spillStorage
         for (auto& alloc : preAggregation.allocations) {
            for (auto entry = reinterpret_cast<header_t*>(alloc.first),
                      end = addBytes(entry, alloc.second * entry_size);
                 entry < end; entry = addBytes(entry, entry_size))
               spill.push_back(&entry->hash, entry->hash);
         }
         preAggregation.allocations.clear();
         preAggregation.clearHashtable(ht);
      };

      for (pos_t n = child->next(); n != EndOfStream; n = child->next()) {
         groupHash.evaluate(n);
         preAggregation.findGroups(n, ht);
         auto groupsCreated = preAggregation.createMissingGroups(ht, false);
         updateGroups.evaluate(n);
         groups += groupsCreated;
         if (groups >= maxFill) flushAndClear();
      }
      flushAndClear(); // flush remaining entries into spillStorage
      barrier();       // Wait until all workers have finished phase 1

      cont.consumed = true;
      cont.partition = shared.partition.fetch_add(1);
      cont.partitionNeedsAggregation = true;
   }

   /// ------ phase 2: global aggregation
   // get a partition
   for (; cont.partition < nrPartitions;) {
      if (cont.partitionNeedsAggregation) {
         auto partNr = cont.partition;
         // for all thread local partitions
         for (auto& threadPartitions : shared.spillStorage.threadData) {
            // aggregate data from thread local partition
            auto& partition = threadPartitions.second.getPartitions()[partNr];
            for (auto chunk = partition.first; chunk; chunk = chunk->next) {
               auto elementSize = threadPartitions.second.entrySize;
               auto nPart = partition.size(chunk, elementSize);
               for (size_t n = std::min(nPart, vecSize), pos = 0; n;
                    nPart -= n, pos += n, n = std::min(nPart, vecSize)) {

                  // communicate data position of current chunk to primitives
                  // for group lookup and creation
                  auto data = addBytes(chunk->data<void>(), pos * elementSize);
                  globalAggregation.rowData = data;
                  findGroupsFromPartition(data, n);
                  auto cGroups = [&]() INTERPRET_SEPARATE {
                     globalAggregation.createMissingGroups(ht, true);
                  };
                  cGroups();
                  updateGroupsFromPartition.evaluate(n);
               }
            }
         }
         cont.partitionNeedsAggregation = false;
         cont.iter = globalAggregation.allocations.begin();
      }
      // send aggregation result to parent operator
      // TODO: refill result instead of sending all allocations separately
      if (cont.iter != globalAggregation.allocations.end()) {
         auto& block = *cont.iter;
         // write current block start to htMatches so the the gather primitives
         // can read the offset from there
         *globalAggregation.htMatches =
             reinterpret_cast<header_t*>(block.first);
         auto n = block.second;
         gatherGroups.evaluate(n);
         cont.iter++;
         return n;
      } else {
         auto htClear = [&]() INTERPRET_SEPARATE {
            globalAggregation.clearHashtable(ht);
         };
         htClear();
         cont.partitionNeedsAggregation = true;
         cont.partition = shared.partition.fetch_add(1);
      }
   }
   return EndOfStream;
}
} // namespace vectorwise

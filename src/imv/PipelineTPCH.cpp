#include "imv/PipelineTPCH.hpp"
using namespace std;
__m512i v_all_ones = _mm512_set1_epi64(-1), v_zero = _mm512_set1_epi64(0), v_one = _mm512_set1_epi64(1), v_63 = _mm512_set1_epi64(63);

struct AMACStateQ1 {
  uint64_t probeValue[5];
  uint8_t stage;
  int probeKey = 0;
  pos_t tuple_id = 0;
  runtime::Hashmap::hash_t probeHash;
  runtime::Hashmap::EntryHeader* buildMatch;
};
size_t agg_raw_q1(size_t begin, size_t end, Database& db,
                  Hashmapx<types::Integer, tuple<Numeric<12, 2>, Numeric<12, 2>, Numeric<12, 4>, Numeric<12, 6>, int64_t>, hashFun, false>* hash_table,
                  PartitionedDeque<PARTITION_SIZE>* partition, void** entry_addrs, void** results_entry) {
  size_t found = 0, pos = 0, cur = begin;
  types::Date c1 = types::Date::castString("1998-09-02");
  types::Numeric<12, 2> one = types::Numeric<12, 2>::castString("1.00");
  auto& li = db["lineitem"];
  auto l_returnflag = li["l_returnflag"].data<types::Char<1>>();
  auto l_linestatus = li["l_linestatus"].data<types::Char<1>>();
  auto l_extendedprice = li["l_extendedprice"].data<types::Numeric<12, 2>>();
  auto l_discount = li["l_discount"].data<types::Numeric<12, 2>>();
  auto l_tax = li["l_tax"].data<types::Numeric<12, 2>>();
  auto l_quantity = li["l_quantity"].data<types::Numeric<12, 2>>();
  auto l_shipdate = li["l_shipdate"].data<types::Date>();
  using group_t = Hashmapx<types::Integer, tuple<Numeric<12, 2>, Numeric<12, 2>,
  Numeric<12, 4>, Numeric<12, 6>, int64_t>, hashFun, false>::Entry;
  hash_t hash_value;
  uint32_t probeKey = 0;
  group_t* entry = nullptr, *old_entry = nullptr;
  uint64_t* values, *values_old;
  for (size_t cur = begin; cur < end; ++cur) {
    if (nullptr == entry_addrs) {
      *(char*) (((char*) (&probeKey)) + 0) = l_returnflag[cur].value;
      *(char*) (((char*) (&probeKey)) + 1) = l_linestatus[cur].value;
      hash_value = hashFun()(probeKey, primitives::seed);
      entry = hash_table->findOneEntry(probeKey, hash_value);
      if (!entry) {
        entry = (group_t*) partition->partition_allocate(hash_value);
        entry->h.hash = hash_value;
        entry->h.next = nullptr;
        entry->k = probeKey;
        values = (uint64_t*) (((char*) entry) + offsetof(group_t, v));
        // initialize values
        *values = types::Numeric<12, 2>().value;
        *(values + 1) = types::Numeric<12, 2>().value;
        *(values + 2) = types::Numeric<12, 4>().value;
        *(values + 3) = types::Numeric<12, 6>().value;
        *(values + 4) = 0;
        hash_table->insert<false>(*entry);
        ++found;
      }
      values = (uint64_t*) (((char*) entry) + offsetof(group_t, v));
      // update aggregators
      *values += l_quantity[cur].value;
      *(values + 1) += l_extendedprice[cur].value;
      auto disc_price = l_extendedprice[cur] * (one - l_discount[cur]);
      *(values + 2) += disc_price.value;
      auto charge = disc_price * (one + l_tax[cur]);
      *(values + 3) += charge.value;
      *(values + 4) += 1;
    } else {
      old_entry = (group_t*) entry_addrs[cur];
      entry = hash_table->findOneEntry(old_entry->k, old_entry->h.hash);
      if (!entry) {
        old_entry->h.next = nullptr;
        hash_table->insert<false>(*old_entry);
        results_entry[found++] = entry_addrs[cur];
      } else {
        // update aggregators with old's
        values = (uint64_t*) (((char*) entry) + offsetof(group_t, v));
        values_old = (uint64_t*) (((char*) old_entry) + offsetof(group_t, v));
        *values += *values_old;
        *(values + 1) += *(values_old + 1);
        *(values + 2) += *(values_old + 2);
        *(values + 3) += *(values_old + 3);
        *(values + 4) += *(values_old + 4);
      }
    }
  }

  return found;
}

size_t agg_gp_q1(size_t begin, size_t end, Database& db,
                 Hashmapx<types::Integer, tuple<Numeric<12, 2>, Numeric<12, 2>, Numeric<12, 4>, Numeric<12, 6>, int64_t>, hashFun, false>* hash_table,
                 PartitionedDeque<PARTITION_SIZE>* partition, void** entry_addrs, void** results_entry) {
  size_t found = 0, pos = 0, cur = begin;
  int k = 0, done = 0, buildkey, probeKey = 0, valid_size;
  AMACStateQ1 state[stateNum];
  hash_t probeHash;
  types::Date c1 = types::Date::castString("1998-09-02");
  types::Numeric<12, 2> one = types::Numeric<12, 2>::castString("1.00");
  auto& li = db["lineitem"];
  auto l_returnflag = li["l_returnflag"].data<types::Char<1>>();
  auto l_linestatus = li["l_linestatus"].data<types::Char<1>>();
  auto l_extendedprice = li["l_extendedprice"].data<types::Numeric<12, 2>>();
  auto l_discount = li["l_discount"].data<types::Numeric<12, 2>>();
  auto l_tax = li["l_tax"].data<types::Numeric<12, 2>>();
  auto l_quantity = li["l_quantity"].data<types::Numeric<12, 2>>();
  auto l_shipdate = li["l_shipdate"].data<types::Date>();
  using group_t = Hashmapx<types::Integer, tuple<Numeric<12, 2>, Numeric<12, 2>,
  Numeric<12, 4>, Numeric<12, 6>, int64_t>, hashFun, false>::Entry;
  hash_t hash_value;
  group_t* entry = nullptr, *old_entry = nullptr;
  uint64_t* values, *values_old;

  auto insetNewEntry = [&](AMACStateQ1& state) {
    if(nullptr == entry_addrs) {
      entry = (group_t*) partition->partition_allocate(state.probeHash);
      entry->h.hash = state.probeHash;
      entry->h.next = nullptr;
      entry->k = types::Integer(state.probeKey);
      values = (uint64_t*) (((char*) entry) + offsetof(group_t, v));
      values[0]=state.probeValue[0];
      values[1]=state.probeValue[1];
      values[2]=state.probeValue[2];
      values[3]=state.probeValue[3];
      values[4]=state.probeValue[4];
    } else {
      entry = (group_t*) entry_addrs[state.tuple_id];
      entry->h.next = nullptr;
      results_entry[found] = entry_addrs[state.tuple_id];
    }
    auto lastEntry = (group_t*) state.buildMatch;
    if(lastEntry == nullptr) { /* the bucket is empty*/
      hash_table->insert<false>(*entry);
    } else {
      lastEntry->h.next = (decltype(lastEntry->h.next))entry;
    }

    ++found;
  };
  while (cur < end) {
    /// step 1: get the hash key and compute hash value
    for (k = 0; (k < stateNum) && (cur < end); ++k, ++cur) {
      if (nullptr == entry_addrs) {
#if SEQ_PREFETCH
        _mm_prefetch((((char* )(l_returnflag+cur))+PDIS), _MM_HINT_T0);
        _mm_prefetch((((char* )(l_returnflag+cur))+PDIS+64), _MM_HINT_T0);
        _mm_prefetch((((char* )(l_linestatus+cur))+PDIS), _MM_HINT_T0);
        _mm_prefetch((((char* )(l_linestatus+cur))+PDIS+64), _MM_HINT_T0);
        _mm_prefetch((((char* )(l_quantity+cur))+PDIS), _MM_HINT_T0);
        _mm_prefetch((((char* )(l_quantity+cur))+PDIS+64), _MM_HINT_T0);
        _mm_prefetch((((char* )(l_extendedprice+cur))+PDIS), _MM_HINT_T0);
        _mm_prefetch((((char* )(l_extendedprice+cur))+PDIS+64), _MM_HINT_T0);
        _mm_prefetch((((char* )(l_discount+cur))+PDIS), _MM_HINT_T0);
        _mm_prefetch((((char* )(l_discount+cur))+PDIS+64), _MM_HINT_T0);
        _mm_prefetch((((char* )(l_tax+cur))+PDIS), _MM_HINT_T0);
        _mm_prefetch((((char* )(l_tax+cur))+PDIS+64), _MM_HINT_T0);
#endif
        *(char*) (((char*) (&probeKey)) + 0) = l_returnflag[cur].value;
        *(char*) (((char*) (&probeKey)) + 1) = l_linestatus[cur].value;
        probeHash = (hashFun()(probeKey, primitives::seed));
        state[k].probeValue[0] = l_quantity[cur].value;
        state[k].probeValue[1] = l_extendedprice[cur].value;
        auto disc_price = l_extendedprice[cur] * (one - l_discount[cur]);
        state[k].probeValue[2] = disc_price.value;
        auto charge = disc_price * (one + l_tax[cur]);
        state[k].probeValue[3] = charge.value;
        state[k].probeValue[4] = 1;
        state[k].tuple_id = cur;
        state[k].probeKey = probeKey;
        state[k].probeHash = probeHash;
      } else {
#if SEQ_PREFETCH
        _mm_prefetch((((char* )(entry_addrs[cur+PDISD]))), _MM_HINT_T0);
        _mm_prefetch((((char* )(entry_addrs[cur+PDISD]))+64), _MM_HINT_T0);
#endif
        entry = (group_t*) entry_addrs[cur];
        values = (uint64_t*) (((char*) entry) + offsetof(group_t, v));
        state[k].probeValue[0] = values[0];
        state[k].probeValue[1] = values[1];
        state[k].probeValue[2] = values[2];
        state[k].probeValue[3] = values[3];
        state[k].probeValue[4] = values[4];
        state[k].probeHash = entry->h.hash;
        state[k].tuple_id = cur;
        state[k].probeKey = entry->k.value;
      }
      state[k].stage = 0;
      if (stateNum > 1)
        hash_table->PrefetchEntry(state[k].probeHash);
    }
    valid_size = k;
    done = 0;
    /// step 2: fetch the first node in the hash table bucket
    for (k = 0; k < valid_size; ++k) {
      state[k].buildMatch = hash_table->find_chain(state[k].probeHash);
      if (nullptr == state[k].buildMatch) {
        //// must immediately write a new entry
        insetNewEntry(state[k]);
        state[k].stage = 4;
        ++done;
      } else {
        if (stateNum > 1) {
          _mm_prefetch((char * )(state[k].buildMatch), _MM_HINT_T0);
          _mm_prefetch((char * )(state[k].buildMatch) + 64, _MM_HINT_T0);
        }
      }
    }
    /// step 3: repeating probing the hash buckets
    while (done < valid_size) {
      for (k = 0; k < valid_size; ++k) {
        // done or need to insert
        if (state[k].stage >= 3) {
          continue;
        }
        entry = (group_t*) state[k].buildMatch;
        buildkey = entry->k.value;
        // found, then update the aggregators
        if ((buildkey == state[k].probeKey)) {
          values = (uint64_t*) (((char*) entry) + offsetof(group_t, v));
          values[0] += state[k].probeValue[0];
          values[1] += state[k].probeValue[1];
          values[2] += state[k].probeValue[2];
          values[3] += state[k].probeValue[3];
          values[4] += state[k].probeValue[4];
          state[k].stage = 3;
          ++done;
          continue;
        }
        auto entryHeader = entry->h.next;
        // not found, to insert
        if (nullptr == entryHeader) {
          //// must immediately write a new entry
          insetNewEntry(state[k]);
          state[k].stage = 4;
          ++done;
          continue;
        } else {
          // not found, then continue
          state[k].buildMatch = entryHeader;
          if (stateNum > 1) {
            _mm_prefetch((char * )(entryHeader), _MM_HINT_T0);
            _mm_prefetch((char * )(entryHeader) + 64, _MM_HINT_T0);
          }
        }
      }
    }
  }

  return found;
}
size_t agg_amac_q1(size_t begin, size_t end, Database& db,
                   Hashmapx<types::Integer, tuple<Numeric<12, 2>, Numeric<12, 2>, Numeric<12, 4>, Numeric<12, 6>, int64_t>, hashFun, false>* hash_table,
                   PartitionedDeque<PARTITION_SIZE>* partition, void** entry_addrs, void** results_entry) {
  size_t found = 0, pos = 0, cur = begin;
  int k = 0, done = 0, buildkey, probeKey = 0, valid_size;
  AMACStateQ1 state[stateNum];
  hash_t probeHash;
  types::Date c1 = types::Date::castString("1998-09-02");
  types::Numeric<12, 2> one = types::Numeric<12, 2>::castString("1.00");
  auto& li = db["lineitem"];
  auto l_returnflag = li["l_returnflag"].data<types::Char<1>>();
  auto l_linestatus = li["l_linestatus"].data<types::Char<1>>();
  auto l_extendedprice = li["l_extendedprice"].data<types::Numeric<12, 2>>();
  auto l_discount = li["l_discount"].data<types::Numeric<12, 2>>();
  auto l_tax = li["l_tax"].data<types::Numeric<12, 2>>();
  auto l_quantity = li["l_quantity"].data<types::Numeric<12, 2>>();
  auto l_shipdate = li["l_shipdate"].data<types::Date>();
  using group_t = Hashmapx<types::Integer, tuple<Numeric<12, 2>, Numeric<12, 2>,
  Numeric<12, 4>, Numeric<12, 6>, int64_t>, hashFun, false>::Entry;
  hash_t hash_value;
  group_t* entry = nullptr, *old_entry = nullptr;
  uint64_t* values, *values_old;

  // initialization
  for (int i = 0; i < stateNum; ++i) {
    state[i].stage = 1;
  }

  while (done < stateNum) {
    k = (k >= stateNum) ? 0 : k;
    switch (state[k].stage) {
      case 1: {
        if (cur >= end) {
          ++done;
          state[k].stage = 3;
          break;
        }
        if (nullptr == entry_addrs) {
#if SEQ_PREFETCH
          _mm_prefetch((((char* )(l_returnflag+cur))+PDIS), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_returnflag+cur))+PDIS+64), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_linestatus+cur))+PDIS), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_linestatus+cur))+PDIS+64), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_quantity+cur))+PDIS), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_quantity+cur))+PDIS+64), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_extendedprice+cur))+PDIS), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_extendedprice+cur))+PDIS+64), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_discount+cur))+PDIS), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_discount+cur))+PDIS+64), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_tax+cur))+PDIS), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_tax+cur))+PDIS+64), _MM_HINT_T0);
#endif
          *(char*) (((char*) (&probeKey)) + 0) = l_returnflag[cur].value;
          *(char*) (((char*) (&probeKey)) + 1) = l_linestatus[cur].value;
          probeHash = (hashFun()(probeKey, primitives::seed));
          state[k].probeValue[0] = l_quantity[cur].value;
          state[k].probeValue[1] = l_extendedprice[cur].value;
          auto disc_price = l_extendedprice[cur] * (one - l_discount[cur]);
          state[k].probeValue[2] = disc_price.value;
          auto charge = disc_price * (one + l_tax[cur]);
          state[k].probeValue[3] = charge.value;
          state[k].probeValue[4] = 1;
          state[k].tuple_id = cur;
          state[k].probeKey = probeKey;
          state[k].probeHash = probeHash;
        } else {
#if SEQ_PREFETCH
          _mm_prefetch((((char* )(entry_addrs[cur+PDISD]))), _MM_HINT_T0);
          _mm_prefetch((((char* )(entry_addrs[cur+PDISD]))+64), _MM_HINT_T0);
#endif
          entry = (group_t*) entry_addrs[cur];
          values = (uint64_t*) (((char*) entry) + offsetof(group_t, v));
          state[k].probeValue[0] = values[0];
          state[k].probeValue[1] = values[1];
          state[k].probeValue[2] = values[2];
          state[k].probeValue[3] = values[3];
          state[k].probeValue[4] = values[4];
          state[k].probeHash = entry->h.hash;
          state[k].tuple_id = cur;
          state[k].probeKey = entry->k.value;
        }
        ++cur;
        state[k].stage = 2;
        if (stateNum > 1)
          hash_table->PrefetchEntry(state[k].probeHash);
      }
        break;
      case 2: {
        state[k].buildMatch = hash_table->find_chain(state[k].probeHash);
        if (nullptr == state[k].buildMatch) {
          state[k].stage = 4;
          --k;  // must immediately shift to case 4
        } else {
          if (stateNum > 1) {
            _mm_prefetch((char * )(state[k].buildMatch), _MM_HINT_T0);
            _mm_prefetch((char * )(state[k].buildMatch) + 64, _MM_HINT_T0);
          }
          state[k].stage = 0;
        }
      }
        break;
      case 0: {
        entry = (group_t*) state[k].buildMatch;
        buildkey = entry->k.value;
        if ((buildkey == state[k].probeKey)) {
          values = (uint64_t*) (((char*) entry) + offsetof(group_t, v));
          values[0] += state[k].probeValue[0];
          values[1] += state[k].probeValue[1];
          values[2] += state[k].probeValue[2];
          values[3] += state[k].probeValue[3];
          values[4] += state[k].probeValue[4];
          state[k].stage = 1;
          --k;
          break;
        }
        auto entryHeader = entry->h.next;
        if (nullptr == entryHeader) {
          state[k].stage = 4;
          --k;  // must immediately shift to case 4
        } else {
          state[k].buildMatch = entryHeader;
          if (stateNum > 1) {
            _mm_prefetch((char * )(entryHeader), _MM_HINT_T0);
            _mm_prefetch((char * )(entryHeader) + 64, _MM_HINT_T0);
          }
        }
      }
        break;
      case 4: {
        if (nullptr == entry_addrs) {
          entry = (group_t*) partition->partition_allocate(state[k].probeHash);
          entry->h.hash = state[k].probeHash;
          entry->h.next = nullptr;
          entry->k = types::Integer(state[k].probeKey);
          values = (uint64_t*) (((char*) entry) + offsetof(group_t, v));
          values[0] = state[k].probeValue[0];
          values[1] = state[k].probeValue[1];
          values[2] = state[k].probeValue[2];
          values[3] = state[k].probeValue[3];
          values[4] = state[k].probeValue[4];
        } else {
          entry = (group_t*) entry_addrs[state[k].tuple_id];
          entry->h.next = nullptr;
          results_entry[found] = entry_addrs[state[k].tuple_id];
        }
        auto lastEntry = (group_t*) state[k].buildMatch;
        if (lastEntry == nullptr) { /* the bucket is empty*/
          hash_table->insert<false>(*entry);
        } else {
          lastEntry->h.next = (decltype(lastEntry->h.next)) entry;
        }
        state[k].stage = 1;
        ++found;
        --k;
      }
        break;
    }
    ++k;
  }
  return found;
}
struct __attribute__((aligned(64))) AggStateQ1 {
  __m512i v_probe_keys;
  __m512i v_bucket_addrs;
  __m512i v_probe_offset, v_probe_hash, v_probe_value[5];
  __mmask8 m_valid_probe;
  uint8_t stage;
  AggStateQ1()
      : v_probe_keys(_mm512_set1_epi64(0)),
        v_bucket_addrs(_mm512_set1_epi64(0)),
        v_probe_offset(_mm512_set1_epi64(0)),
        m_valid_probe(0),
        stage(1) {
  }
  void reset() {
    v_probe_keys = _mm512_set1_epi64(0);
    v_bucket_addrs = _mm512_set1_epi64(0);
    v_probe_offset = _mm512_set1_epi64(0);
    m_valid_probe = 0;
    stage = 1;
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
    v_bucket_addrs = _mm512_maskz_compress_epi64(m_valid_probe, v_bucket_addrs);
    v_probe_keys = _mm512_maskz_compress_epi64(m_valid_probe, v_probe_keys);
    v_probe_value[0] = _mm512_maskz_compress_epi64(m_valid_probe, v_probe_value[0]);
    v_probe_value[1] = _mm512_maskz_compress_epi64(m_valid_probe, v_probe_value[1]);
    v_probe_value[2] = _mm512_maskz_compress_epi64(m_valid_probe, v_probe_value[2]);
    v_probe_value[3] = _mm512_maskz_compress_epi64(m_valid_probe, v_probe_value[3]);
    v_probe_value[4] = _mm512_maskz_compress_epi64(m_valid_probe, v_probe_value[4]);
    v_probe_offset = _mm512_maskz_compress_epi64(m_valid_probe, v_probe_offset);
    v_probe_hash = _mm512_maskz_compress_epi64(m_valid_probe, v_probe_hash);

  }
  inline void expand(AggStateQ1& src_state) {
    v_bucket_addrs = _mm512_mask_expand_epi64(v_bucket_addrs, _mm512_knot(m_valid_probe), src_state.v_bucket_addrs);
    v_probe_keys = _mm512_mask_expand_epi64(v_probe_keys, _mm512_knot(m_valid_probe), src_state.v_probe_keys);
    v_probe_value[0] = _mm512_mask_expand_epi64(v_probe_value[0], _mm512_knot(m_valid_probe), src_state.v_probe_value[0]);
    v_probe_value[1] = _mm512_mask_expand_epi64(v_probe_value[1], _mm512_knot(m_valid_probe), src_state.v_probe_value[1]);
    v_probe_value[2] = _mm512_mask_expand_epi64(v_probe_value[2], _mm512_knot(m_valid_probe), src_state.v_probe_value[2]);
    v_probe_value[3] = _mm512_mask_expand_epi64(v_probe_value[3], _mm512_knot(m_valid_probe), src_state.v_probe_value[3]);
    v_probe_value[4] = _mm512_mask_expand_epi64(v_probe_value[4], _mm512_knot(m_valid_probe), src_state.v_probe_value[4]);
    v_probe_offset = _mm512_mask_expand_epi64(v_probe_offset, _mm512_knot(m_valid_probe), src_state.v_probe_offset);
    v_probe_hash = _mm512_mask_expand_epi64(v_probe_hash, _mm512_knot(m_valid_probe), src_state.v_probe_hash);

  }
};
inline void insertNewEntry(AggStateQ1& state, Vec8u& u_new_addrs, __mmask8 m_no_conflict, PartitionedDeque<PARTITION_SIZE>* partition, __m512i& u_offset_hash, __m512i& u_offset_k, __m512i& u_offset_v) {
  Vec8u u_probe_hash(state.v_probe_hash);
  for(int i=0;i<VECTORSIZE;++i) {
    u_new_addrs.entry[i] =0;
    if(m_no_conflict & (1<<i)) {
      u_new_addrs.entry[i] = (uint64_t)partition->partition_allocate(u_probe_hash.entry[i]);
    }
  }
  // write entry->next
  _mm512_mask_i64scatter_epi64(0,m_no_conflict,u_new_addrs.reg,v_zero,1);
  // write entry->hash
  _mm512_mask_i64scatter_epi64(0,m_no_conflict,u_new_addrs + u_offset_hash,state.v_probe_hash,1);
  // write entry->k , NOTE it is 32 bits
  _mm512_mask_i64scatter_epi32(0,m_no_conflict,u_new_addrs + u_offset_k,_mm512_cvtepi64_epi32(state.v_probe_keys),1);
  // write entry->v
  _mm512_mask_i64scatter_epi64(0,m_no_conflict,u_new_addrs + u_offset_v,state.v_probe_value[0],1);
  _mm512_mask_i64scatter_epi64(0,m_no_conflict,u_new_addrs + u_offset_v+Vec8u(uint64_t(8)),state.v_probe_value[1],1);
  _mm512_mask_i64scatter_epi64(0,m_no_conflict,u_new_addrs + u_offset_v+Vec8u(uint64_t(16)),state.v_probe_value[2],1);
  _mm512_mask_i64scatter_epi64(0,m_no_conflict,u_new_addrs + u_offset_v+Vec8u(uint64_t(24)),state.v_probe_value[3],1);
  _mm512_mask_i64scatter_epi64(0,m_no_conflict,u_new_addrs + u_offset_v+Vec8u(uint64_t(32)),state.v_probe_value[4],1);

}
inline void mergeKeys(AggStateQ1& state) {
  auto v_keys = _mm512_mask_blend_epi64(state.m_valid_probe, v_all_ones, state.v_probe_keys);
  auto v_conflict = _mm512_conflict_epi64(v_keys);
  auto m_no_conflict = _mm512_testn_epi64_mask(v_conflict, v_all_ones);
  m_no_conflict = _mm512_kand(m_no_conflict, state.m_valid_probe);
  if (m_no_conflict == 0)
    return;
  uint64_t* pos_v[5];
  pos_v[0] = (uint64_t*) &state.v_probe_value[0];
  pos_v[1] = (uint64_t*) &state.v_probe_value[1];
  pos_v[2] = (uint64_t*) &state.v_probe_value[2];
  pos_v[3] = (uint64_t*) &state.v_probe_value[3];
  pos_v[4] = (uint64_t*) &state.v_probe_value[4];
  auto m_conflict = _mm512_kandn(m_no_conflict, state.m_valid_probe);

  auto v_lzeros = _mm512_lzcnt_epi64(v_conflict);
  v_lzeros = _mm512_sub_epi64(v_63, v_lzeros);
  uint64_t* pos_lz = (uint64_t*) &v_lzeros;
  for (int i = VECTORSIZE - 1; i >= 0; --i) {
    if ((m_conflict & (1 << i))) {
      pos_v[0][pos_lz[i]] += pos_v[0][i];
      pos_v[1][pos_lz[i]] += pos_v[1][i];
      pos_v[2][pos_lz[i]] += pos_v[2][i];
      pos_v[3][pos_lz[i]] += pos_v[3][i];
      pos_v[4][pos_lz[i]] += pos_v[4][i];
    }
  }
  state.m_valid_probe = m_no_conflict;
  state.v_probe_keys = _mm512_mask_blend_epi64(m_no_conflict, v_all_ones, state.v_probe_keys);
  return;
}
inline void mergeKeys_(AggStateQ1& state) {
  auto v_keys = _mm512_mask_blend_epi64(state.m_valid_probe, v_all_ones, state.v_probe_keys);
  auto v_conflict = _mm512_conflict_epi64(v_keys);
  auto m_no_conflict = _mm512_testn_epi64_mask(v_conflict, v_all_ones);
  m_no_conflict = _mm512_kand(m_no_conflict, state.m_valid_probe);
  if (m_no_conflict == 0)
    return;
  uint64_t* pos_v[5];
  pos_v[0] = (uint64_t*) &state.v_probe_value[0];
  pos_v[1] = (uint64_t*) &state.v_probe_value[1];
  pos_v[2] = (uint64_t*) &state.v_probe_value[2];
  pos_v[3] = (uint64_t*) &state.v_probe_value[3];
  pos_v[4] = (uint64_t*) &state.v_probe_value[4];
  auto m_conflict = _mm512_kandn(m_no_conflict, state.m_valid_probe);

  auto v_lzeros = _mm512_lzcnt_epi64(v_conflict);
  v_lzeros = _mm512_sub_epi64(v_63, v_lzeros);
  uint64_t* pos_lz = (uint64_t*) &v_lzeros;
  for (int i = VECTORSIZE - 1; i >= 0; --i) {
    if ((m_conflict & (1 << i))) {
      pos_v[0][pos_lz[i]] += pos_v[0][i];
      pos_v[1][pos_lz[i]] += pos_v[1][i];
      pos_v[2][pos_lz[i]] += pos_v[2][i];
      pos_v[3][pos_lz[i]] += pos_v[3][i];
      pos_v[4][pos_lz[i]] += pos_v[4][i];
    }
  }
  state.m_valid_probe = m_no_conflict;
  state.v_probe_keys = _mm512_mask_blend_epi64(m_no_conflict, v_all_ones, state.v_probe_keys);
  return;
}
size_t agg_imv_q1(size_t begin, size_t end, Database& db,
                  Hashmapx<types::Integer, tuple<Numeric<12, 2>, Numeric<12, 2>, Numeric<12, 4>, Numeric<12, 6>, int64_t>, hashFun, false>* hash_table,
                  PartitionedDeque<PARTITION_SIZE>* partition, void** entry_addrs, void** results_entry) {
  size_t found = 0, pos = 0, cur = begin;
  int k = 0, done = 0, buildkey, probeKey, valid_size, imvNum = stateNumSIMD;
  AggStateQ1 state[stateNumSIMD + 1];
  hash_t probeHash;
  types::Date c1 = types::Date::castString("1998-09-02");
  types::Numeric<12, 2> one = types::Numeric<12, 2>::castString("1.00");
  auto& li = db["lineitem"];
  auto l_returnflag = li["l_returnflag"].data<types::Char<1>>();
  auto l_linestatus = li["l_linestatus"].data<types::Char<1>>();
  auto l_extendedprice = li["l_extendedprice"].data<types::Numeric<12, 2>>();
  auto l_discount = li["l_discount"].data<types::Numeric<12, 2>>();
  auto l_tax = li["l_tax"].data<types::Numeric<12, 2>>();
  auto l_quantity = li["l_quantity"].data<types::Numeric<12, 2>>();
  auto l_shipdate = li["l_shipdate"].data<types::Date>();
  using group_t = Hashmapx<types::Integer, tuple<Numeric<12, 2>, Numeric<12, 2>,
  Numeric<12, 4>, Numeric<12, 6>, int64_t>, hashFun, false>::Entry;
  hash_t hash_value;
  group_t* entry = nullptr, *old_entry = nullptr;
  uint64_t* values, *values_old;

  __m512i v_base_offset = _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0), v_zero = _mm512_set1_epi64(0), v_lzeros, v_63 = _mm512_set1_epi64(63);
  __m512i v_offset = _mm512_set1_epi64(0), v_base_offset_upper = _mm512_set1_epi64(end - begin), v_seed = _mm512_set1_epi64(vectorwise::primitives::seed), v_all_ones =
      _mm512_set1_epi64(-1), v_conflict, v_ht_keys, v_hash_mask, v_ht_value, v_next, v_one_num = _mm512_set1_epi64(one.value);
  Vec8u u_new_addrs(uint64_t(0)), u_offset_hash(offsetof(group_t, h.hash)), u_offset_k(offsetof(group_t, k)), u_offset_v(offsetof(group_t, v)), u_8(uint64_t(8)), u_16(
      uint64_t(16)), u_24(uint64_t(24)), u_32(uint64_t(32));
  __mmask8 m_no_conflict, m_rest, m_match, m_to_insert, mask[VECTORSIZE + 1];
  __m256i v256_zero = _mm256_set1_epi32(0), v256_probe_keys, v256_probe_value, v256_ht_keys;
  uint8_t num, num_temp;
  Vec8u u_keys(uint64_t(0));
  for (int i = 0; i <= VECTORSIZE; ++i) {
    mask[i] = (1 << i) - 1;
  }

  while (true) {
    k = (k >= imvNum) ? 0 : k;
    if ((cur >= end)) {
      if (state[k].m_valid_probe == 0 && state[k].stage != 3) {
        ++done;
        state[k].stage = 3;
        ++k;
        continue;
      }
    }
    if (done >= imvNum) {

      if (state[imvNum].m_valid_probe > 0) {
        k = imvNum;
      } else {
        break;
      }

    }
    switch (state[k].stage) {
      case 1: {
        /// step 1: get offsets
        state[k].v_probe_offset = _mm512_add_epi64(_mm512_set1_epi64(cur), v_base_offset);
        state[k].m_valid_probe = -1;
        if (cur + VECTORSIZE >= end) {
          state[k].m_valid_probe = (state[k].m_valid_probe >> (cur + VECTORSIZE - end));
        }
        if (nullptr == entry_addrs) {
          /// step 2: gather probe keys and values
#if SEQ_PREFETCH
          _mm_prefetch((((char* )(l_returnflag+cur))+PDIS), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_returnflag+cur))+PDIS+64), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_linestatus+cur))+PDIS), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_linestatus+cur))+PDIS+64), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_quantity+cur))+PDIS), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_quantity+cur))+PDIS+64), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_extendedprice+cur))+PDIS), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_extendedprice+cur))+PDIS+64), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_discount+cur))+PDIS), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_discount+cur))+PDIS+64), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_tax+cur))+PDIS), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_tax+cur))+PDIS+64), _MM_HINT_T0);
#endif
          for (int i = 0; (i < VECTORSIZE) && (cur + i < end); ++i) {
            *(char*) (((char*) (&u_keys.entry[i])) + 0) = l_returnflag[cur + i].value;
            *(char*) (((char*) (&u_keys.entry[i])) + 1) = l_linestatus[cur + i].value;
          }
          state[k].v_probe_keys = u_keys.reg;
          state[k].v_probe_value[0] = _mm512_mask_i64gather_epi64(v_zero, state[k].m_valid_probe, state[k].v_probe_offset, (const long long int* )l_quantity, 8);
          state[k].v_probe_value[1] = _mm512_mask_i64gather_epi64(v_zero, state[k].m_valid_probe, state[k].v_probe_offset, (const long long int* )l_extendedprice, 8);
          state[k].v_probe_value[2] = _mm512_mask_i64gather_epi64(v_zero, state[k].m_valid_probe, state[k].v_probe_offset, (const long long int* )l_discount, 8);
          state[k].v_probe_value[2] = _mm512_mullo_epi64(state[k].v_probe_value[1], _mm512_sub_epi64(v_one_num, state[k].v_probe_value[2]));
          state[k].v_probe_value[3] = _mm512_mask_i64gather_epi64(v_zero, state[k].m_valid_probe, state[k].v_probe_offset, (const long long int* )l_tax, 8);
          state[k].v_probe_value[3] = _mm512_mullo_epi64(state[k].v_probe_value[2], _mm512_add_epi64(v_one_num, state[k].v_probe_value[3]));
          state[k].v_probe_value[4] = v_one;
          /// step 3: compute hash values
          state[k].v_probe_hash = runtime::MurMurHash()((state[k].v_probe_keys), (v_seed));
        } else {
#if SEQ_PREFETCH
          _mm_prefetch((((char* )(entry_addrs[cur+PDISD]))), _MM_HINT_T0);
          _mm_prefetch((((char* )(entry_addrs[cur+PDISD]))+64), _MM_HINT_T0);
#endif
          // gather the addresses of entries
          state[k].v_probe_offset = _mm512_mask_i64gather_epi64(v_zero, state[k].m_valid_probe, state[k].v_probe_offset, (const long long int * )entry_addrs, 8);
          v256_probe_keys = _mm512_mask_i64gather_epi32(v256_zero, state[k].m_valid_probe, state[k].v_probe_offset + u_offset_k, nullptr, 1);
          state[k].v_probe_keys = _mm512_cvtepi32_epi64(v256_probe_keys);
          state[k].v_probe_hash = _mm512_mask_i64gather_epi64(v_zero, state[k].m_valid_probe, state[k].v_probe_offset+u_offset_hash, nullptr, 1);
          state[k].v_probe_value[0] = _mm512_mask_i64gather_epi64(v_zero, state[k].m_valid_probe, state[k].v_probe_offset+u_offset_v, nullptr, 1);
          state[k].v_probe_value[1] = _mm512_mask_i64gather_epi64(v_zero, state[k].m_valid_probe, state[k].v_probe_offset+u_offset_v+u_8, nullptr, 1);
          state[k].v_probe_value[2] = _mm512_mask_i64gather_epi64(v_zero, state[k].m_valid_probe, state[k].v_probe_offset+u_offset_v+u_16, nullptr, 1);
          state[k].v_probe_value[3] = _mm512_mask_i64gather_epi64(v_zero, state[k].m_valid_probe, state[k].v_probe_offset+u_offset_v+u_24, nullptr, 1);
          state[k].v_probe_value[4] = _mm512_mask_i64gather_epi64(v_zero, state[k].m_valid_probe, state[k].v_probe_offset+u_offset_v+u_32, nullptr, 1);
        }

        cur += VECTORSIZE;
        state[k].stage = 2;
        if (stateNumSIMD > 1)
          hash_table->prefetchEntry(state[k].v_probe_hash);
      }
        break;
      case 2: {

        mergeKeys(state[k]);
        /// step 4: find the addresses of corresponding buckets for new probes
        Vec8uM v_new_bucket_addrs = hash_table->find_chain(state[k].v_probe_hash);

        /// insert new nodes in the corresponding hash buckets
        m_to_insert = _mm512_kandn(v_new_bucket_addrs.mask, state[k].m_valid_probe);
        v_hash_mask = ((Vec8u(state[k].v_probe_hash) & Vec8u(hash_table->mask)));
        v_hash_mask = _mm512_mask_blend_epi64(state[k].m_valid_probe, v_all_ones, v_hash_mask);
        v_conflict = _mm512_conflict_epi64(v_hash_mask);
        m_no_conflict = _mm512_testn_epi64_mask(v_conflict, v_all_ones);
        m_no_conflict = _mm512_kand(m_no_conflict, m_to_insert);
        if (m_no_conflict) {
          if (nullptr == entry_addrs) {
            insertNewEntry(state[k], u_new_addrs, m_no_conflict, partition, u_offset_hash.reg, u_offset_k.reg, u_offset_v.reg);
            // insert the new addresses to the hash table
            _mm512_mask_i64scatter_epi64((long long int* )hash_table->entries, m_no_conflict, v_hash_mask, u_new_addrs.reg, 8);
          } else {
            // set the next of entries = 0
            _mm512_mask_i64scatter_epi64(nullptr, m_no_conflict, state[k].v_probe_offset, v_zero, 1);
            // must write back the merged values!!!!!!!!!
            _mm512_mask_i64scatter_epi64(nullptr, m_no_conflict, state[k].v_probe_offset+u_offset_v, state[k].v_probe_value[0], 1);
            _mm512_mask_i64scatter_epi64(nullptr, m_no_conflict, state[k].v_probe_offset+u_offset_v + u_8, state[k].v_probe_value[1], 1);
            _mm512_mask_i64scatter_epi64(nullptr, m_no_conflict, state[k].v_probe_offset+u_offset_v + u_16, state[k].v_probe_value[2], 1);
            _mm512_mask_i64scatter_epi64(nullptr, m_no_conflict, state[k].v_probe_offset+u_offset_v + u_24, state[k].v_probe_value[3], 1);
            _mm512_mask_i64scatter_epi64(nullptr, m_no_conflict, state[k].v_probe_offset+u_offset_v + u_32, state[k].v_probe_value[4], 1);

            _mm512_mask_i64scatter_epi64((long long int* )hash_table->entries, m_no_conflict, v_hash_mask, state[k].v_probe_offset, 8);
            _mm512_mask_compressstoreu_epi64((results_entry + found), m_no_conflict, state[k].v_probe_offset);
          }
          found += _mm_popcnt_u32(m_no_conflict);
          // get rid of no-conflict elements
          state[k].m_valid_probe = _mm512_kandn(m_no_conflict, state[k].m_valid_probe);
        }
        state[k].v_bucket_addrs = _mm512_mask_i64gather_epi64(v_all_ones, state[k].m_valid_probe, v_hash_mask, (const long long int* )hash_table->entries, 8);
#if 0
        v_prefetch(state[k].v_bucket_addrs);
        state[k].stage = 0;
#else
        num = _mm_popcnt_u32(state[k].m_valid_probe);
        if (num == VECTORSIZE) {
          if (stateNumSIMD > 1)
            v_prefetch(state[k].v_bucket_addrs);
          state[k].stage = 0;
        } else {
          if ((done < imvNum)) {
            num_temp = _mm_popcnt_u32(state[imvNum].m_valid_probe);
            if (num + num_temp < VECTORSIZE) {
              // compress imv_state[k]
              state[k].compress();
              // expand imv_state[k] -> imv_state[imvNum]
              state[imvNum].expand(state[k]);
              state[imvNum].m_valid_probe = mask[num + num_temp];
              state[k].m_valid_probe = 0;
              state[k].stage = 1;
              state[imvNum].stage = 0;
              --k;
              break;
            } else {
              // expand imv_state[imvNum] -> expand imv_state[k]
              state[k].expand(state[imvNum]);
              state[imvNum].m_valid_probe = _mm512_kand(state[imvNum].m_valid_probe, _mm512_knot(mask[VECTORSIZE - num]));
              // compress imv_state[imvNum]
              state[imvNum].compress();
              state[imvNum].m_valid_probe = state[imvNum].m_valid_probe >> (VECTORSIZE - num);
              state[k].m_valid_probe = mask[VECTORSIZE];
              state[k].stage = 0;
              state[imvNum].stage = 0;
              if (stateNumSIMD > 1)
                v_prefetch(state[k].v_bucket_addrs);
            }
          }
        }
#endif
      }
        break;
      case 0: {
        /// step 5: gather the all new build keys
        v256_ht_keys = _mm512_mask_i64gather_epi32(v256_zero, state[k].m_valid_probe, _mm512_add_epi64(state[k].v_bucket_addrs, u_offset_k.reg), nullptr, 1);
        v_ht_keys = _mm512_cvtepi32_epi64(v256_ht_keys);
        /// step 6: compare the probe keys and build keys and write points
        m_match = _mm512_cmpeq_epi64_mask(state[k].v_probe_keys, v_ht_keys);
        m_match = _mm512_kand(m_match, state[k].m_valid_probe);

        /// update the aggregators
        v_conflict = _mm512_conflict_epi64(state[k].v_bucket_addrs);
        m_no_conflict = _mm512_testn_epi64_mask(v_conflict, v_all_ones);
        m_no_conflict = _mm512_kand(m_no_conflict, m_match);

        v_ht_value = _mm512_mask_i64gather_epi64(v_zero, m_no_conflict, _mm512_add_epi64(state[k].v_bucket_addrs, u_offset_v.reg), nullptr, 1);
        _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(state[k].v_bucket_addrs, u_offset_v.reg), _mm512_add_epi64(state[k].v_probe_value[0], v_ht_value), 1);
        v_ht_value = _mm512_mask_i64gather_epi64(v_zero, m_no_conflict, _mm512_add_epi64(state[k].v_bucket_addrs, u_offset_v.reg) + u_8, nullptr, 1);
        _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(state[k].v_bucket_addrs, u_offset_v.reg) + u_8, _mm512_add_epi64(state[k].v_probe_value[1], v_ht_value), 1);
        v_ht_value = _mm512_mask_i64gather_epi64(v_zero, m_no_conflict, _mm512_add_epi64(state[k].v_bucket_addrs, u_offset_v.reg) + u_16, nullptr, 1);
        _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(state[k].v_bucket_addrs, u_offset_v.reg)+ u_16, _mm512_add_epi64(state[k].v_probe_value[2], v_ht_value), 1);
        v_ht_value = _mm512_mask_i64gather_epi64(v_zero, m_no_conflict, _mm512_add_epi64(state[k].v_bucket_addrs, u_offset_v.reg) + u_24, nullptr, 1);
        _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(state[k].v_bucket_addrs, u_offset_v.reg)+ u_24, _mm512_add_epi64(state[k].v_probe_value[3], v_ht_value), 1);
        v_ht_value = _mm512_mask_i64gather_epi64(v_zero, m_no_conflict, _mm512_add_epi64(state[k].v_bucket_addrs, u_offset_v.reg) + u_32, nullptr, 1);
        _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(state[k].v_bucket_addrs, u_offset_v.reg)+ u_32, _mm512_add_epi64(state[k].v_probe_value[4], v_ht_value), 1);

        state[k].m_valid_probe = _mm512_kandn(m_no_conflict, state[k].m_valid_probe);
        // the remaining matches, DO NOT get next
        m_match = _mm512_kandn(m_no_conflict, m_match);

        /// step 7: NOT found, then insert
        v_next = _mm512_mask_i64gather_epi64(v_all_ones, _mm512_kandn(m_match, state[k].m_valid_probe), state[k].v_bucket_addrs, nullptr, 1);
        m_to_insert = _mm512_kand(_mm512_kandn(m_match, state[k].m_valid_probe), _mm512_cmpeq_epi64_mask(v_next, v_zero));
        // get rid of bucket address of matched probes
        v_next = _mm512_mask_blend_epi64(_mm512_kandn(m_match, state[k].m_valid_probe), v_all_ones, state[k].v_bucket_addrs);
        v_conflict = _mm512_conflict_epi64(v_next);
        m_no_conflict = _mm512_testn_epi64_mask(v_conflict, v_all_ones);
        m_no_conflict = _mm512_kand(m_no_conflict, m_to_insert);
        if (m_no_conflict) {
          if (nullptr == entry_addrs) {
            insertNewEntry(state[k], u_new_addrs, m_no_conflict, partition, u_offset_hash.reg, u_offset_k.reg, u_offset_v.reg);
            // insert the new addresses to the hash table
            _mm512_mask_i64scatter_epi64(0, m_no_conflict, state[k].v_bucket_addrs, u_new_addrs.reg, 1);
          } else {
            // set the next of entries = 0
            _mm512_mask_i64scatter_epi64(nullptr, m_no_conflict, state[k].v_probe_offset, v_zero, 1);
            _mm512_mask_i64scatter_epi64(nullptr, m_no_conflict, state[k].v_bucket_addrs, state[k].v_probe_offset, 1);
            _mm512_mask_compressstoreu_epi64((results_entry + found), m_no_conflict, state[k].v_probe_offset);
          }
          found += _mm_popcnt_u32(m_no_conflict);

          state[k].m_valid_probe = _mm512_kandn(m_no_conflict, state[k].m_valid_probe);
        }
        v_next = _mm512_mask_i64gather_epi64(v_all_ones, state[k].m_valid_probe, state[k].v_bucket_addrs, nullptr, 1);
        // the remaining matches, DO NOT get next
        state[k].v_bucket_addrs = _mm512_mask_blend_epi64(m_match, v_next, state[k].v_bucket_addrs);

        num = _mm_popcnt_u32(state[k].m_valid_probe);
        if (num == VECTORSIZE || done >= imvNum) {
          if (stateNumSIMD > 1)
            v_prefetch(state[k].v_bucket_addrs);
        } else {
          if ((done < imvNum)) {
            num_temp = _mm_popcnt_u32(state[imvNum].m_valid_probe);
            if (num + num_temp < VECTORSIZE) {
              // compress imv_state[k]
              state[k].compress();
              // expand imv_state[k] -> imv_state[imvNum]
              state[imvNum].expand(state[k]);
              state[imvNum].m_valid_probe = mask[num + num_temp];
              state[k].m_valid_probe = 0;
              state[k].stage = 1;
              state[imvNum].stage = 0;
              --k;
              break;
            } else {
              // expand imv_state[imvNum] -> expand imv_state[k]
              state[k].expand(state[imvNum]);
              state[imvNum].m_valid_probe = _mm512_kand(state[imvNum].m_valid_probe, _mm512_knot(mask[VECTORSIZE - num]));
              // compress imv_state[imvNum]
              state[imvNum].compress();
              state[imvNum].m_valid_probe = state[imvNum].m_valid_probe >> (VECTORSIZE - num);
              state[k].m_valid_probe = mask[VECTORSIZE];
              state[k].stage = 0;
              state[imvNum].stage = 0;
              if (stateNumSIMD > 1)
                v_prefetch(state[k].v_bucket_addrs);
            }
          }
        }

      }
        break;
    }
    ++k;
  }

  return found;
}
size_t agg_simd_q1(size_t begin, size_t end, Database& db,
                   Hashmapx<types::Integer, tuple<Numeric<12, 2>, Numeric<12, 2>, Numeric<12, 4>, Numeric<12, 6>, int64_t>, hashFun, false>* hash_table,
                   PartitionedDeque<PARTITION_SIZE>* partition, void** entry_addrs, void** results_entry) {
  size_t found = 0, pos = 0, cur = begin;
  int stateNumSIMD = 1;
  int k = 0, done = 0, buildkey, probeKey, valid_size, imvNum = stateNumSIMD;
  AggStateQ1 state[stateNumSIMD + 1];
  hash_t probeHash;
  types::Date c1 = types::Date::castString("1998-09-02");
  types::Numeric<12, 2> one = types::Numeric<12, 2>::castString("1.00");
  auto& li = db["lineitem"];
  auto l_returnflag = li["l_returnflag"].data<types::Char<1>>();
  auto l_linestatus = li["l_linestatus"].data<types::Char<1>>();
  auto l_extendedprice = li["l_extendedprice"].data<types::Numeric<12, 2>>();
  auto l_discount = li["l_discount"].data<types::Numeric<12, 2>>();
  auto l_tax = li["l_tax"].data<types::Numeric<12, 2>>();
  auto l_quantity = li["l_quantity"].data<types::Numeric<12, 2>>();
  auto l_shipdate = li["l_shipdate"].data<types::Date>();
  using group_t = Hashmapx<types::Integer, tuple<Numeric<12, 2>, Numeric<12, 2>,
  Numeric<12, 4>, Numeric<12, 6>, int64_t>, hashFun, false>::Entry;
  hash_t hash_value;
  group_t* entry = nullptr, *old_entry = nullptr;
  uint64_t* values, *values_old;

  __m512i v_base_offset = _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0), v_zero = _mm512_set1_epi64(0), v_lzeros, v_63 = _mm512_set1_epi64(63);
  __m512i v_offset = _mm512_set1_epi64(0), v_base_offset_upper = _mm512_set1_epi64(end - begin), v_seed = _mm512_set1_epi64(vectorwise::primitives::seed), v_all_ones =
      _mm512_set1_epi64(-1), v_conflict, v_ht_keys, v_hash_mask, v_ht_value, v_next, v_one_num = _mm512_set1_epi64(one.value);
  Vec8u u_new_addrs(uint64_t(0)), u_offset_hash(offsetof(group_t, h.hash)), u_offset_k(offsetof(group_t, k)), u_offset_v(offsetof(group_t, v)), u_8(uint64_t(8)), u_16(
      uint64_t(16)), u_24(uint64_t(24)), u_32(uint64_t(32));
  __mmask8 m_no_conflict, m_rest, m_match, m_to_insert, mask[VECTORSIZE + 1];
  __m256i v256_zero = _mm256_set1_epi32(0), v256_probe_keys, v256_probe_value, v256_ht_keys;
  uint8_t num, num_temp;
  Vec8u u_keys(uint64_t(0));
  for (int i = 0; i <= VECTORSIZE; ++i) {
    mask[i] = (1 << i) - 1;
  }

  while (true) {
    k = (k >= stateNumSIMD) ? 0 : k;
    if ((cur >= end)) {
      if (state[k].m_valid_probe == 0 && state[k].stage != 3) {
        ++done;
        state[k].stage = 3;
        ++k;
        continue;
      }
    }
    if (done >= imvNum) {
      if (state[imvNum].m_valid_probe > 0) {
        k = imvNum;
      } else {
        break;
      }

    }
    switch (state[k].stage) {
      case 1: {
        /// step 1: get offsets
        state[k].v_probe_offset = _mm512_add_epi64(_mm512_set1_epi64(cur), v_base_offset);
        state[k].m_valid_probe = -1;
        if (cur + VECTORSIZE >= end) {
          state[k].m_valid_probe = (state[k].m_valid_probe >> (cur + VECTORSIZE - end));
        }
        if (nullptr == entry_addrs) {
          /// step 2: gather probe keys and values
#if SEQ_PREFETCH
          _mm_prefetch((((char* )(l_returnflag+cur))+PDIS), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_returnflag+cur))+PDIS+64), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_linestatus+cur))+PDIS), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_linestatus+cur))+PDIS+64), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_quantity+cur))+PDIS), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_quantity+cur))+PDIS+64), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_extendedprice+cur))+PDIS), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_extendedprice+cur))+PDIS+64), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_discount+cur))+PDIS), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_discount+cur))+PDIS+64), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_tax+cur))+PDIS), _MM_HINT_T0);
          _mm_prefetch((((char* )(l_tax+cur))+PDIS+64), _MM_HINT_T0);
#endif
          for (int i = 0; (i < VECTORSIZE) && (cur + i < end); ++i) {
            *(char*) (((char*) (&u_keys.entry[i])) + 0) = l_returnflag[cur + i].value;
            *(char*) (((char*) (&u_keys.entry[i])) + 1) = l_linestatus[cur + i].value;
          }
          state[k].v_probe_keys = u_keys.reg;
          state[k].v_probe_value[0] = _mm512_mask_i64gather_epi64(v_zero, state[k].m_valid_probe, state[k].v_probe_offset, (const long long int* )l_quantity, 8);
          state[k].v_probe_value[1] = _mm512_mask_i64gather_epi64(v_zero, state[k].m_valid_probe, state[k].v_probe_offset, (const long long int* )l_extendedprice, 8);
          state[k].v_probe_value[2] = _mm512_mask_i64gather_epi64(v_zero, state[k].m_valid_probe, state[k].v_probe_offset, (const long long int* )l_discount, 8);
          state[k].v_probe_value[2] = _mm512_mullo_epi64(state[k].v_probe_value[1], _mm512_sub_epi64(v_one_num, state[k].v_probe_value[2]));
          state[k].v_probe_value[3] = _mm512_mask_i64gather_epi64(v_zero, state[k].m_valid_probe, state[k].v_probe_offset, (const long long int* )l_tax, 8);
          state[k].v_probe_value[3] = _mm512_mullo_epi64(state[k].v_probe_value[2], _mm512_add_epi64(v_one_num, state[k].v_probe_value[3]));
          state[k].v_probe_value[4] = v_one;
          /// step 3: compute hash values
          state[k].v_probe_hash = runtime::MurMurHash()((state[k].v_probe_keys), (v_seed));
        } else {
#if SEQ_PREFETCH
          _mm_prefetch((((char* )(entry_addrs[cur+PDISD]))), _MM_HINT_T0);
          _mm_prefetch((((char* )(entry_addrs[cur+PDISD]))+64), _MM_HINT_T0);
#endif
          // gather the addresses of entries
          state[k].v_probe_offset = _mm512_mask_i64gather_epi64(v_zero, state[k].m_valid_probe, state[k].v_probe_offset, (const long long int * )entry_addrs, 8);
          v256_probe_keys = _mm512_mask_i64gather_epi32(v256_zero, state[k].m_valid_probe, state[k].v_probe_offset + u_offset_k, nullptr, 1);
          state[k].v_probe_keys = _mm512_cvtepi32_epi64(v256_probe_keys);
          state[k].v_probe_hash = _mm512_mask_i64gather_epi64(v_zero, state[k].m_valid_probe, state[k].v_probe_offset+u_offset_hash, nullptr, 1);
          state[k].v_probe_value[0] = _mm512_mask_i64gather_epi64(v_zero, state[k].m_valid_probe, state[k].v_probe_offset+u_offset_v, nullptr, 1);
          state[k].v_probe_value[1] = _mm512_mask_i64gather_epi64(v_zero, state[k].m_valid_probe, state[k].v_probe_offset+u_offset_v+u_8, nullptr, 1);
          state[k].v_probe_value[2] = _mm512_mask_i64gather_epi64(v_zero, state[k].m_valid_probe, state[k].v_probe_offset+u_offset_v+u_16, nullptr, 1);
          state[k].v_probe_value[3] = _mm512_mask_i64gather_epi64(v_zero, state[k].m_valid_probe, state[k].v_probe_offset+u_offset_v+u_24, nullptr, 1);
          state[k].v_probe_value[4] = _mm512_mask_i64gather_epi64(v_zero, state[k].m_valid_probe, state[k].v_probe_offset+u_offset_v+u_32, nullptr, 1);
        }

        cur += VECTORSIZE;
//        state[k].stage = 2;
//        if (stateNumSIMD > 1)
//          hash_table->prefetchEntry(state[k].v_probe_hash);
//      }
//        break;
//      case 2: {

        mergeKeys(state[k]);
        /// step 4: find the addresses of corresponding buckets for new probes
        Vec8uM v_new_bucket_addrs = hash_table->find_chain(state[k].v_probe_hash);

        /// insert new nodes in the corresponding hash buckets
        m_to_insert = _mm512_kandn(v_new_bucket_addrs.mask, state[k].m_valid_probe);
        v_hash_mask = ((Vec8u(state[k].v_probe_hash) & Vec8u(hash_table->mask)));
        v_hash_mask = _mm512_mask_blend_epi64(state[k].m_valid_probe, v_all_ones, v_hash_mask);
        v_conflict = _mm512_conflict_epi64(v_hash_mask);
        m_no_conflict = _mm512_testn_epi64_mask(v_conflict, v_all_ones);
        m_no_conflict = _mm512_kand(m_no_conflict, m_to_insert);
        if (m_no_conflict) {
          if (nullptr == entry_addrs) {
            insertNewEntry(state[k], u_new_addrs, m_no_conflict, partition, u_offset_hash.reg, u_offset_k.reg, u_offset_v.reg);
            // insert the new addresses to the hash table
            _mm512_mask_i64scatter_epi64((long long int* )hash_table->entries, m_no_conflict, v_hash_mask, u_new_addrs.reg, 8);
          } else {
            // set the next of entries = 0
            _mm512_mask_i64scatter_epi64(nullptr, m_no_conflict, state[k].v_probe_offset, v_zero, 1);
            // must write back the merged values!!!!!!!!!
            _mm512_mask_i64scatter_epi64(nullptr, m_no_conflict, state[k].v_probe_offset+u_offset_v, state[k].v_probe_value[0], 1);
            _mm512_mask_i64scatter_epi64(nullptr, m_no_conflict, state[k].v_probe_offset+u_offset_v + u_8, state[k].v_probe_value[1], 1);
            _mm512_mask_i64scatter_epi64(nullptr, m_no_conflict, state[k].v_probe_offset+u_offset_v + u_16, state[k].v_probe_value[2], 1);
            _mm512_mask_i64scatter_epi64(nullptr, m_no_conflict, state[k].v_probe_offset+u_offset_v + u_24, state[k].v_probe_value[3], 1);
            _mm512_mask_i64scatter_epi64(nullptr, m_no_conflict, state[k].v_probe_offset+u_offset_v + u_32, state[k].v_probe_value[4], 1);

            _mm512_mask_i64scatter_epi64((long long int* )hash_table->entries, m_no_conflict, v_hash_mask, state[k].v_probe_offset, 8);
            _mm512_mask_compressstoreu_epi64((results_entry + found), m_no_conflict, state[k].v_probe_offset);
          }
          found += _mm_popcnt_u32(m_no_conflict);
          // get rid of no-conflict elements
          state[k].m_valid_probe = _mm512_kandn(m_no_conflict, state[k].m_valid_probe);
        }
        state[k].v_bucket_addrs = _mm512_mask_i64gather_epi64(v_all_ones, state[k].m_valid_probe, v_hash_mask, (const long long int* )hash_table->entries, 8);
#if 0
        v_prefetch(state[k].v_bucket_addrs);
        state[k].stage = 0;
#else
        num = _mm_popcnt_u32(state[k].m_valid_probe);
        if (num == VECTORSIZE) {
//          if (stateNumSIMD > 1)
//            v_prefetch(state[k].v_bucket_addrs);
          state[k].stage = 0;
        } else {
          if ((done < imvNum)) {
            num_temp = _mm_popcnt_u32(state[imvNum].m_valid_probe);
            if (num + num_temp < VECTORSIZE) {
              // compress imv_state[k]
              state[k].compress();
              // expand imv_state[k] -> imv_state[imvNum]
              state[imvNum].expand(state[k]);
              state[imvNum].m_valid_probe = mask[num + num_temp];
              state[k].m_valid_probe = 0;
              state[k].stage = 1;
              state[imvNum].stage = 0;
              --k;
              break;
            } else {
              // expand imv_state[imvNum] -> expand imv_state[k]
              state[k].expand(state[imvNum]);
              state[imvNum].m_valid_probe = _mm512_kand(state[imvNum].m_valid_probe, _mm512_knot(mask[VECTORSIZE - num]));
              // compress imv_state[imvNum]
              state[imvNum].compress();
              state[imvNum].m_valid_probe = state[imvNum].m_valid_probe >> (VECTORSIZE - num);
              state[k].m_valid_probe = mask[VECTORSIZE];
              state[k].stage = 0;
              state[imvNum].stage = 0;
//              if (stateNumSIMD > 1)
//                v_prefetch(state[k].v_bucket_addrs);
            }
          }
        }
#endif
      }
        break;
      case 0: {
        /// step 5: gather the all new build keys
        v256_ht_keys = _mm512_mask_i64gather_epi32(v256_zero, state[k].m_valid_probe, _mm512_add_epi64(state[k].v_bucket_addrs, u_offset_k.reg), nullptr, 1);
        v_ht_keys = _mm512_cvtepi32_epi64(v256_ht_keys);
        /// step 6: compare the probe keys and build keys and write points
        m_match = _mm512_cmpeq_epi64_mask(state[k].v_probe_keys, v_ht_keys);
        m_match = _mm512_kand(m_match, state[k].m_valid_probe);

        /// update the aggregators
        v_conflict = _mm512_conflict_epi64(state[k].v_bucket_addrs);
        m_no_conflict = _mm512_testn_epi64_mask(v_conflict, v_all_ones);
        m_no_conflict = _mm512_kand(m_no_conflict, m_match);

        v_ht_value = _mm512_mask_i64gather_epi64(v_zero, m_no_conflict, _mm512_add_epi64(state[k].v_bucket_addrs, u_offset_v.reg), nullptr, 1);
        _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(state[k].v_bucket_addrs, u_offset_v.reg), _mm512_add_epi64(state[k].v_probe_value[0], v_ht_value), 1);
        v_ht_value = _mm512_mask_i64gather_epi64(v_zero, m_no_conflict, _mm512_add_epi64(state[k].v_bucket_addrs, u_offset_v.reg) + u_8, nullptr, 1);
        _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(state[k].v_bucket_addrs, u_offset_v.reg) + u_8, _mm512_add_epi64(state[k].v_probe_value[1], v_ht_value), 1);
        v_ht_value = _mm512_mask_i64gather_epi64(v_zero, m_no_conflict, _mm512_add_epi64(state[k].v_bucket_addrs, u_offset_v.reg) + u_16, nullptr, 1);
        _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(state[k].v_bucket_addrs, u_offset_v.reg)+ u_16, _mm512_add_epi64(state[k].v_probe_value[2], v_ht_value), 1);
        v_ht_value = _mm512_mask_i64gather_epi64(v_zero, m_no_conflict, _mm512_add_epi64(state[k].v_bucket_addrs, u_offset_v.reg) + u_24, nullptr, 1);
        _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(state[k].v_bucket_addrs, u_offset_v.reg)+ u_24, _mm512_add_epi64(state[k].v_probe_value[3], v_ht_value), 1);
        v_ht_value = _mm512_mask_i64gather_epi64(v_zero, m_no_conflict, _mm512_add_epi64(state[k].v_bucket_addrs, u_offset_v.reg) + u_32, nullptr, 1);
        _mm512_mask_i64scatter_epi64(0, m_no_conflict, _mm512_add_epi64(state[k].v_bucket_addrs, u_offset_v.reg)+ u_32, _mm512_add_epi64(state[k].v_probe_value[4], v_ht_value), 1);

        state[k].m_valid_probe = _mm512_kandn(m_no_conflict, state[k].m_valid_probe);
        // the remaining matches, DO NOT get next
        m_match = _mm512_kandn(m_no_conflict, m_match);

        /// step 7: NOT found, then insert
        v_next = _mm512_mask_i64gather_epi64(v_all_ones, _mm512_kandn(m_match, state[k].m_valid_probe), state[k].v_bucket_addrs, nullptr, 1);
        m_to_insert = _mm512_kand(_mm512_kandn(m_match, state[k].m_valid_probe), _mm512_cmpeq_epi64_mask(v_next, v_zero));
        // get rid of bucket address of matched probes
        v_next = _mm512_mask_blend_epi64(_mm512_kandn(m_match, state[k].m_valid_probe), v_all_ones, state[k].v_bucket_addrs);
        v_conflict = _mm512_conflict_epi64(v_next);
        m_no_conflict = _mm512_testn_epi64_mask(v_conflict, v_all_ones);
        m_no_conflict = _mm512_kand(m_no_conflict, m_to_insert);
        if (m_no_conflict) {
          if (nullptr == entry_addrs) {
            insertNewEntry(state[k], u_new_addrs, m_no_conflict, partition, u_offset_hash.reg, u_offset_k.reg, u_offset_v.reg);
            // insert the new addresses to the hash table
            _mm512_mask_i64scatter_epi64(0, m_no_conflict, state[k].v_bucket_addrs, u_new_addrs.reg, 1);
          } else {
            // set the next of entries = 0
            _mm512_mask_i64scatter_epi64(nullptr, m_no_conflict, state[k].v_probe_offset, v_zero, 1);
            _mm512_mask_i64scatter_epi64(nullptr, m_no_conflict, state[k].v_bucket_addrs, state[k].v_probe_offset, 1);
            _mm512_mask_compressstoreu_epi64((results_entry + found), m_no_conflict, state[k].v_probe_offset);
          }
          found += _mm_popcnt_u32(m_no_conflict);

          state[k].m_valid_probe = _mm512_kandn(m_no_conflict, state[k].m_valid_probe);
        }
        v_next = _mm512_mask_i64gather_epi64(v_all_ones, state[k].m_valid_probe, state[k].v_bucket_addrs, nullptr, 1);
        // the remaining matches, DO NOT get next
        state[k].v_bucket_addrs = _mm512_mask_blend_epi64(m_match, v_next, state[k].v_bucket_addrs);

        num = _mm_popcnt_u32(state[k].m_valid_probe);
        if (num == VECTORSIZE || done >= imvNum) {
//          if (stateNumSIMD > 1)
//            v_prefetch(state[k].v_bucket_addrs);
        } else {
          if ((done < imvNum)) {
            num_temp = _mm_popcnt_u32(state[imvNum].m_valid_probe);
            if (num + num_temp < VECTORSIZE) {
              // compress imv_state[k]
              state[k].compress();
              // expand imv_state[k] -> imv_state[imvNum]
              state[imvNum].expand(state[k]);
              state[imvNum].m_valid_probe = mask[num + num_temp];
              state[k].m_valid_probe = 0;
              state[k].stage = 1;
              state[imvNum].stage = 0;
              --k;
              break;
            } else {
              // expand imv_state[imvNum] -> expand imv_state[k]
              state[k].expand(state[imvNum]);
              state[imvNum].m_valid_probe = _mm512_kand(state[imvNum].m_valid_probe, _mm512_knot(mask[VECTORSIZE - num]));
              // compress imv_state[imvNum]
              state[imvNum].compress();
              state[imvNum].m_valid_probe = state[imvNum].m_valid_probe >> (VECTORSIZE - num);
              state[k].m_valid_probe = mask[VECTORSIZE];
              state[k].stage = 0;
              state[imvNum].stage = 0;
//              if (stateNumSIMD > 1)
//                v_prefetch(state[k].v_bucket_addrs);
            }
          }
        }

      }
        break;
    }
    ++k;
  }

  return found;
}

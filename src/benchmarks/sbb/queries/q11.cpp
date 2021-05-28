#include <deque>
#include <iostream>

#include "benchmarks/ssb/Queries.hpp"
#include "common/runtime/Hash.hpp"
#include "common/runtime/Types.hpp"
#include "hyper/GroupBy.hpp"
#include "hyper/ParallelHelper.hpp"
#include "tbb/tbb.h"
#include "vectorwise/Operations.hpp"
#include "vectorwise/Operators.hpp"
#include "vectorwise/Primitives.hpp"
#include "vectorwise/QueryBuilder.hpp"
#include "vectorwise/VectorAllocator.hpp"
#include "imv/HashBuildSSB.hpp"
#include "imv/PipelineSSB.hpp"
#include "rof/SimdFilter.hpp"
#include "rof/AmacBuild.hpp"
#include "rof/AmacProbe.hpp"
#include "rof/AmacAgg.hpp"

using namespace runtime;
using namespace std;

namespace ssb {
uint64_t ht_date_size=1024;
NOVECTORIZE std::unique_ptr<runtime::Query> q11_hyper(Database& db, size_t nrThreads) {
  // --- aggregates

  auto resources = initQuery(nrThreads);

  // --- constants
  const auto relevant_year = types::Integer(1993);
  const auto discount_min = types::Numeric<18, 2>::castString("1.00");
  const auto discount_max = types::Numeric<18, 2>::castString("3.00");
  const auto quantity_max = types::Integer(25);

  using hash = runtime::MurMurHash;
 // const size_t morselSize = 100000;

  // --- ht for join date-lineorder
  Hashset<types::Integer, hash> ht;
  tbb::enumerable_thread_specific<runtime::Stack<decltype(ht)::Entry>> entries1;
  auto& d = db["date"];
  auto d_year = d["d_year"].data<types::Integer>();
  auto d_datekey = d["d_datekey"].data<types::Integer>();
  // do selection on part and put selected elements into ht
  auto found = PARALLEL_SELECT(d.nrTuples, entries1, {
    auto& year = d_year[i];
    auto& datekey = d_datekey[i];
    if (year == relevant_year) {
      entries.emplace_back(ht.hash(datekey), datekey)
      ;
      found++
      ;
    }
  });
  ht.setSize(found);
  ht_date_size = found;
  parallel_insert(entries1, ht);
#if DEBUG
  cout << "q11_hyper build num = " << found << endl;
#endif
  // --- scan lineorder
  auto& lo = db["lineorder"];
  auto lo_orderdate = lo["lo_orderdate"].data<types::Integer>();
  auto lo_quantity = lo["lo_quantity"].data<types::Integer>();
  auto lo_discount = lo["lo_discount"].data<types::Numeric<18, 2>>();
  auto lo_extendedprice = lo["lo_extendedprice"].data<types::Numeric<18, 2>>();

  auto result_revenue = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, lo.nrTuples), types::Numeric<18, 4>(0), [&](const tbb::blocked_range<size_t>& r,
      const types::Numeric<18, 4>& s) {
    auto revenue = s;
    for (size_t i = r.begin(), end = r.end(); i != end; ++i) {
      auto& quantity = lo_quantity[i];
      auto& discount = lo_discount[i];
      auto& extendedprice = lo_extendedprice[i];
      if ((quantity < quantity_max) & (discount >= discount_min) &
          (discount <= discount_max)) {
        if (ht.contains(lo_orderdate[i])) {
          /*  --- aggregation*/
          revenue += extendedprice * discount;
        }

      }
    }
    return revenue;
  },
                                             [](const types::Numeric<18, 4>& x, const types::Numeric<18, 4>& y) {
                                               return x + y;
                                             });

  // --- output
  auto& result = resources.query->result;
  auto revAttr = result->addAttribute("revenue", sizeof(types::Numeric<18, 4>));
  auto block = result->createBlock(1);
  auto revenue = static_cast<types::Numeric<18, 4>*>(block.data(revAttr));
  *revenue = result_revenue;
  block.addedElements(1);
#if DEBUG
  cout << "q11_hyper results = " << result_revenue << endl;
#endif
  leaveQuery(nrThreads);
  return std::move(resources.query);
}
NOVECTORIZE std::unique_ptr<runtime::Query> q11_imv(Database& db, size_t nrThreads) {
  // --- aggregates

  auto resources = initQuery(nrThreads);
  const auto relevant_year = types::Integer(1993);

  using hash = runtime::MurMurHash;
  using range = tbb::blocked_range<size_t>;
  const auto add = [](const size_t& a, const size_t& b) {return a + b;};
  const auto numeric_add = [](const types::Numeric<18, 4>& x, const types::Numeric<18, 4>& y) {
    return x + y;
  };
  // --- ht for join date-lineorder
  Hashset<types::Integer, hash> ht;

  auto& d = db["date"];
  ht.setSize(ht_date_size);
  auto entry_size = sizeof(decltype(ht)::Entry);
  auto found = tbb::parallel_reduce(range(0, d.nrTuples, morselSize), 0, [&](const tbb::blocked_range<size_t>& r, const size_t& f) {
    auto found1 = f;
    found1 +=build_imv_q1x(r.begin(),r.end(),db,&ht,&this_worker->allocator,entry_size);
    return found1;
  },
                                    add);
#if DEBUG
  cout << "q11_imv build num = " << found << endl;
#endif
  // --- lineorder scan -> filter -> probe ->fixAgg
  auto& lo = db["lineorder"];
  auto result_revenue = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, lo.nrTuples), types::Numeric<18, 4>(0), [&](const tbb::blocked_range<size_t>& r,
      const types::Numeric<18, 4>& s) {
    auto revenue = s;
    uint64_t res=0;
    pipeline_imv_q1x(r.begin(),r.end(),db,&ht,res);
    revenue.value += res;
    return revenue;
  },
                                             numeric_add);

  // --- output
  auto& result = resources.query->result;
  auto revAttr = result->addAttribute("revenue", sizeof(types::Numeric<18, 4>));
  auto block = result->createBlock(1);
  auto revenue = static_cast<types::Numeric<18, 4>*>(block.data(revAttr));
  *revenue = result_revenue;
  block.addedElements(1);
#if DEBUG
  cout << "q11_imv results = " << result_revenue << endl;
#endif
  leaveQuery(nrThreads);
  return std::move(resources.query);
}
NOVECTORIZE std::unique_ptr<runtime::Query> q11_rof(Database& db, size_t nrThreads) {
  // --- aggregates

  auto resources = initQuery(nrThreads);
  const auto relevant_year = types::Integer(1993);

  using hash = runtime::MurMurHash;
  using range = tbb::blocked_range<size_t>;
  const auto add = [](const size_t& a, const size_t& b) {return a + b;};
  const auto numeric_add = [](const types::Numeric<18, 4>& x, const types::Numeric<18, 4>& y) {
    return x + y;
  };
  // --- ht for join date-lineorder
  Hashset<types::Integer, hash> ht;
  tbb::enumerable_thread_specific<vector<uint64_t>>position;
  auto& d = db["date"];
  ht.setSize(ht_date_size);
  auto entry_size = sizeof(decltype(ht)::Entry);
  auto found = tbb::parallel_reduce(range(0, d.nrTuples, morselSize), 0, [&](const tbb::blocked_range<size_t>& r, const size_t& f) {
    auto found1 = f;
    auto pos_buff = position.local();
    pos_buff.resize(ROF_VECTOR_SIZE);
    pos_buff.clear();
    for (size_t i = r.begin(), size = 0; i < r.end(); ) {
       size = simd_filter_q11_build_date(i,r.end(),db,&pos_buff[0]);
       found1 += amac_build_q11_date(0,size,db,&ht,&this_worker->allocator,entry_size, &pos_buff[0]);
    }
    return found1;
  },
                                    add);
#if DEBUG
  cout << "q11_imv build num = " << found << endl;
#endif
  // --- lineorder scan -> filter -> probe ->fixAgg
  auto& lo = db["lineorder"];
  auto result_revenue = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, lo.nrTuples), types::Numeric<18, 4>(0), [&](const tbb::blocked_range<size_t>& r,
      const types::Numeric<18, 4>& s) {
    auto revenue = s;
    uint64_t res=0;
    auto pos_buff = position.local();
    pos_buff.resize(ROF_VECTOR_SIZE);
    pos_buff.clear();
    for (size_t i = r.begin(), size = 0; i < r.end(); ) {
      size = simd_filter_q11_probe(i,r.end(),db,&pos_buff[0]);
      amac_probe_q11(0,size,db,&ht,res, &pos_buff[0]);
    }
    revenue.value += res;
    return revenue;
  },
                                             numeric_add);

  // --- output
  auto& result = resources.query->result;
  auto revAttr = result->addAttribute("revenue", sizeof(types::Numeric<18, 4>));
  auto block = result->createBlock(1);
  auto revenue = static_cast<types::Numeric<18, 4>*>(block.data(revAttr));
  *revenue = result_revenue;
  block.addedElements(1);
#if DEBUG
  cout << "q11_imv results = " << result_revenue << endl;
#endif
  leaveQuery(nrThreads);
  return std::move(resources.query);
}
std::unique_ptr<Q11Builder::Q11> Q11Builder::getQuery() {
  using namespace vectorwise;

  auto r = make_unique<Q11>();
  auto date = Scan("date");
  // select d_year = 1993
  Select(Expression().addOp(BF(primitives::sel_equal_to_int32_t_col_int32_t_val), Buffer(sel_year, sizeof(pos_t)), Column(date, "d_year"), Value(&r->year)));

  auto lineorder = Scan("lineorder");
  // select lo_discount between 1 and 3, lo_quantity < 25
  Select(
      Expression().addOp(BF(primitives::sel_less_int32_t_col_int32_t_val), Buffer(sel_qty, sizeof(pos_t)), Column(lineorder, "lo_quantity"), Value(&r->quantity_max)).addOp(
          BF(primitives::selsel_greater_equal_int64_t_col_int64_t_val), Buffer(sel_qty, sizeof(pos_t)), Buffer(sel_discount_low, sizeof(pos_t)), Column(lineorder, "lo_discount"),
          Value(&r->discount_min)).addOp(BF(primitives::selsel_less_equal_int64_t_col_int64_t_val), Buffer(sel_discount_low, sizeof(pos_t)),
                                         Buffer(sel_discount_high, sizeof(pos_t)), Column(lineorder, "lo_discount"), Value(&r->discount_max)));

  HashJoin(Buffer(join_result, sizeof(pos_t)), conf.joinAll()).setProbeSelVector(Buffer(sel_discount_high), conf.joinSel()).addBuildKey(Column(date, "d_datekey"), Buffer(sel_year),
                                                                                                                                        conf.hash_sel_int32_t_col(),
                                                                                                                                        primitives::scatter_sel_int32_t_col)
      .addProbeKey(Column(lineorder, "lo_orderdate"), Buffer(sel_discount_high), conf.hash_sel_int32_t_col(), primitives::keys_equal_int32_t_col);

  Project().addExpression(
      Expression().addOp(primitives::proj_sel_both_multiplies_int64_t_col_int64_t_col, Buffer(join_result), Buffer(result_project, sizeof(int64_t)),
                         Column(lineorder, "lo_extendedprice"), Column(lineorder, "lo_discount")));

  FixedAggregation(Expression().addOp(primitives::aggr_static_plus_int64_t_col, Value(&r->aggregator), Buffer(result_project)));
  r->rootOp = popOperator();
  assert(operatorStack.size() == 0);
  return r;
}

std::unique_ptr<runtime::Query> q11_vectorwise(Database& db, size_t nrThreads, size_t vectorSize) {
  using namespace vectorwise;

  // runtime::Relation result;
  auto result = std::make_unique<runtime::Query>();
  vectorwise::SharedStateManager shared;
  WorkerGroup workers(nrThreads);
  std::atomic<size_t> n;
  std::atomic<int64_t> aggr;
  aggr = 0;
  n = 0;
  workers.run([&]() {
    Q11Builder b(db, shared, vectorSize);
    b.previous = result->participate();
    auto query = b.getQuery();
    auto n_ = query->rootOp->next();
    if (n_) {
      aggr.fetch_add(query->aggregator);
      n.fetch_add(n_);
    }

    auto leader = barrier();
    if (leader) {
      auto attribute = result->result->addAttribute(
          "revenue", sizeof(types::Numeric<18, 4>));
      if (n.load()) {
        auto block = result->result->createBlock(1);
        int64_t* revenue = static_cast<int64_t*>(block.data(attribute));
        *revenue = aggr.load();
        block.addedElements(1);
#if DEBUG
        cout<<"q11_vectorwise results = "<<types::Numeric<18, 4>(*revenue)<<endl;
#endif
      }
    }
  });

  return result;
}

}  // namespace ssb

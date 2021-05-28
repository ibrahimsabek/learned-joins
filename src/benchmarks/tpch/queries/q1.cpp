#include "benchmarks/tpch/Queries.hpp"
#include "common/runtime/Hash.hpp"
#include "common/runtime/Types.hpp"
//#include "hyper/GroupBy.hpp"
//#include "hyper/ParallelHelper.hpp"
#include "tbb/tbb.h"
#include "vectorwise/Operations.hpp"
#include "vectorwise/Operators.hpp"
#include "vectorwise/Primitives.hpp"
#include "vectorwise/QueryBuilder.hpp"
#include "vectorwise/VectorAllocator.hpp"
#include <deque>
#include <iostream>
#include "imv/PipelineTPCH.hpp"

using namespace runtime;
using namespace std;
using vectorwise::primitives::Char_1;
using vectorwise::primitives::hash_t;

//  select
//    l_returnflag,
//    l_linestatus,
//    sum(l_quantity) as sum_qty,
//    sum(l_extendedprice) as sum_base_price,
//    sum(l_extendedprice * (1 - l_discount)) as sum_disc_price,
//    sum(l_extendedprice * (1 - l_discount) * (1 + l_tax)) as sum_charge,
//    avg(l_quantity) as avg_qty,
//    avg(l_extendedprice) as avg_price,
//    avg(l_discount) as avg_disc,
//    count(*) as count_order
//  from
//    lineitem
//  where
//    l_shipdate <= date '1998-12-01' - interval '90' day
//  group by
//    l_returnflag,
//    l_linestatus
void printResultQ1(BlockRelation* result) {
  size_t found = 0;
  auto retAttr = result->getAttribute("l_returnflag");
  auto statusAttr = result->getAttribute("l_linestatus");
  auto qtyAttr = result->getAttribute("sum_qty");
  auto base_priceAttr = result->getAttribute("sum_base_price");
  auto disc_priceAttr = result->getAttribute("sum_disc_price");
  auto chargeAttr = result->getAttribute("sum_charge");
  auto count_orderAttr = result->getAttribute("count_order");
  for (auto& block : *result) {
    auto elementsInBlock = block.size();
    found += elementsInBlock;
    auto ret = reinterpret_cast<Char<1>*>(block.data(retAttr));
    auto status = reinterpret_cast<Char<1>*>(block.data(statusAttr));
    auto qty = reinterpret_cast<types::Numeric<12, 2>*>(block.data(qtyAttr));
    auto base = reinterpret_cast<types::Numeric<12, 2>*>(block.data(base_priceAttr));
    auto disc = reinterpret_cast<types::Numeric<12, 2>*>(block.data(disc_priceAttr));
    auto charge = reinterpret_cast<types::Numeric<12, 2>*>(block.data(chargeAttr));
    auto order = reinterpret_cast<int64_t*>(block.data(count_orderAttr));

    for (size_t i = 0; i < elementsInBlock; ++i) {
      cout << ret[i] << "\t" << status[i] << "\t" << qty[i] << "\t" << base[i] << "\t" << disc[i] << "\t" << charge[i] << "\t" << order[i] << endl;
    }
  }
  cout << "total results number = " << found << endl;
}

/*NOVECTORIZE std::unique_ptr<runtime::Query> q1_hyper(Database& db, size_t nrThreads) {
  using namespace types;
  using namespace std;
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

  auto resources = initQuery(nrThreads);

//  using hash = runtime::CRC32Hash;
  using hash = runtime::MurMurHash;
  auto groupOp = make_GroupBy<tuple<Char<1>, Char<1>>, tuple<Numeric<12, 2>, Numeric<12, 2>, Numeric<12, 4>, Numeric<12, 6>, int64_t>, hash>(
      [](auto& acc, auto&& value) {
        get<0>(acc) += get<0>(value);
        get<1>(acc) += get<1>(value);
        get<2>(acc) += get<2>(value);
        get<3>(acc) += get<3>(value);
        get<4>(acc) += get<4>(value);
      },
      make_tuple(Numeric<12, 2>(), Numeric<12, 2>(), Numeric<12, 4>(), Numeric<12, 6>(), int64_t(0)), nrThreads);

  tbb::parallel_for(tbb::blocked_range<size_t>(0, li.nrTuples, morselSize), [&](const tbb::blocked_range<size_t>& r) {
    auto locals = groupOp.preAggLocals();
    for (size_t i = r.begin(), end = r.end(); i != end; ++i) {
//          if (l_shipdate[i] <= c1) {
                    auto& group = locals.getGroup(make_tuple(l_returnflag[i], l_linestatus[i]));

                    get<0>(group) += l_quantity[i];
                    get<1>(group) += l_extendedprice[i];
                    auto disc_price = l_extendedprice[i] * (one - l_discount[i]);
                    get<2>(group) += disc_price;
                    auto charge = disc_price * (one + l_tax[i]);
                    get<3>(group) += charge;
                    get<4>(group) += 1;
//                       }
    }
  });

  auto& result = resources.query->result;
  auto retAttr = result->addAttribute("l_returnflag", sizeof(Char<1> ));
  auto statusAttr = result->addAttribute("l_linestatus", sizeof(Char<1> ));
  auto qtyAttr = result->addAttribute("sum_qty", sizeof(Numeric<12, 2> ));
  auto base_priceAttr = result->addAttribute("sum_base_price", sizeof(Numeric<12, 2> ));
  auto disc_priceAttr = result->addAttribute("sum_disc_price", sizeof(Numeric<12, 2> ));
  auto chargeAttr = result->addAttribute("sum_charge", sizeof(Numeric<12, 2> ));
  auto count_orderAttr = result->addAttribute("count_order", sizeof(int64_t));

  //groupOp.forallGroups([&](runtime::Stack<decltype(groupOp)::group_t>& auto&entries) {
    auto n = entries.size();
    auto block = result->createBlock(n);
    auto ret = reinterpret_cast<Char<1>*>(block.data(retAttr));
    auto status = reinterpret_cast<Char<1>*>(block.data(statusAttr));
    auto qty = reinterpret_cast<Numeric<12, 2>*>(block.data(qtyAttr));
    auto base_price =
    reinterpret_cast<Numeric<12, 2>*>(block.data(base_priceAttr));
    auto disc_price =
    reinterpret_cast<Numeric<12, 4>*>(block.data(disc_priceAttr));
    auto charge = reinterpret_cast<Numeric<12, 6>*>(block.data(chargeAttr));
    auto count_order =
    reinterpret_cast<int64_t*>(block.data(count_orderAttr));
    for (auto block : entries)
    for (auto& entry : block) {
      *ret++ = get<0>(entry.k);
      *status++ = get<1>(entry.k);
      *qty++ = get<0>(entry.v);
      *base_price++ = get<1>(entry.v);
      *disc_price++ = get<2>(entry.v);
      *charge++ = get<3>(entry.v);
      *count_order++ = get<4>(entry.v);
    }
    block.addedElements(n);
  });
  printResultQ1(result.get());
  leaveQuery(nrThreads);
  return move(resources.query);
}
*/
std::unique_ptr<runtime::Query> q1_imv(Database& db, size_t nrThreads) {
  using namespace types;
  using namespace std;
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

  auto resources = initQuery(nrThreads);
  using hash = runtime::MurMurHash;
  using range = tbb::blocked_range<size_t>;
  const auto add = [](const size_t& a, const size_t& b) {return a + b;};

  tbb::enumerable_thread_specific<Hashmapx<types::Integer, tuple<Numeric<12, 2>, Numeric<12, 2>, Numeric<12, 4>, Numeric<12, 6>, int64_t>, hash, false>> hash_table;

  using group_t = typename decltype(hash_table)::value_type::Entry;

  /// Memory for materialized entries in hashmap
  tbb::enumerable_thread_specific<runtime::Stack<group_t>> entries;
  /// Memory for spilling hastable entries
  tbb::enumerable_thread_specific<runtime::PartitionedDeque<PARTITION_SIZE>> partitionedDeques;
  /// globally collect entry addresses
  tbb::enumerable_thread_specific<vector<group_t*>> entry_addrs;
  tbb::enumerable_thread_specific<vector<group_t*>> results_addrs;

  // local aggregation
  auto found2 = tbb::parallel_reduce(range(0, li.nrTuples, morselSize), 0, [&](const tbb::blocked_range<size_t>& r, const size_t& f) {
    auto found = f;
    bool exist=false;
    auto& ht = hash_table.local(exist);

    if(!exist) {
      ht.setSize(1024);

    }
    auto& partition = partitionedDeques.local(exist);
    if(!exist) {
      partition.postConstruct(nrThreads * 4, sizeof(group_t));
    }
    found += agg_imv_q1(r.begin(),r.end(),db,&ht,&partition,nullptr,nullptr);
    return found;
  },
                                     add);
//  cout << "local agg num = " << found2 << endl;
  auto& result = resources.query->result;
  auto retAttr = result->addAttribute("l_returnflag", sizeof(Char<1> ));
  auto statusAttr = result->addAttribute("l_linestatus", sizeof(Char<1> ));
  auto qtyAttr = result->addAttribute("sum_qty", sizeof(Numeric<12, 2> ));
  auto base_priceAttr = result->addAttribute("sum_base_price", sizeof(Numeric<12, 2> ));
  auto disc_priceAttr = result->addAttribute("sum_disc_price", sizeof(Numeric<12, 2> ));
  auto chargeAttr = result->addAttribute("sum_charge", sizeof(Numeric<12, 2> ));
  auto count_orderAttr = result->addAttribute("count_order", sizeof(int64_t));

  // global aggregation
  auto nrPartitions = partitionedDeques.begin()->getPartitions().size();

  tbb::parallel_for(0ul, nrPartitions, [&](auto partitionNr) {
    bool exist=false;
    auto& ht = hash_table.local(exist);
    if(!exist) {
      ht.setSize(1024);
    }
    ht.clear();
    auto& partition = partitionedDeques.local(exist);
    if(!exist) {
      partition.postConstruct(nrThreads * 4, sizeof(group_t));
    }

    /*  collect entry addresses from a partition */

    auto& entry_addrs_ = entry_addrs.local();
    entry_addrs_.clear();
    auto& results_addrs_ = results_addrs.local();
    results_addrs_.clear();
    for (auto& deque : partitionedDeques) {
      auto& partition = deque.getPartitions()[partitionNr];
      for (auto chunk = partition.first; chunk; chunk = chunk->next) {
        for (auto value = chunk->template data<group_t>(), end = value + partition.size(chunk, sizeof(group_t)); value < end; value++) {
          entry_addrs_.push_back(value);
        }
      }
    }
    results_addrs_.resize(entry_addrs_.size());
    auto found = agg_imv_q1(0,entry_addrs_.size(),db,&ht,nullptr,(void**)&entry_addrs_[0],(void**)&results_addrs_[0]);
    if(found>0) {
      auto block = result->createBlock(found);
      auto ret = reinterpret_cast<Char<1>*>(block.data(retAttr));
      auto status = reinterpret_cast<Char<1>*>(block.data(statusAttr));
      auto qty = reinterpret_cast<Numeric<12, 2>*>(block.data(qtyAttr));
      auto base_price =
      reinterpret_cast<Numeric<12, 2>*>(block.data(base_priceAttr));
      auto disc_price =
      reinterpret_cast<Numeric<12, 4>*>(block.data(disc_priceAttr));
      auto charge = reinterpret_cast<Numeric<12, 6>*>(block.data(chargeAttr));
      auto count_order =
      reinterpret_cast<int64_t*>(block.data(count_orderAttr));

      for(int i=0;i<found;++i) {
       auto values = (uint64_t*) (((char*) results_addrs_[i]) + offsetof(group_t, v));
        *ret++ = Char<1>::build((((char*)&(results_addrs_[i]->k))+0));
        *status++ = Char<1>::build(((char*)&(results_addrs_[i]->k))+1);
        qty->value = *values;
        base_price->value = *(values+1);
        disc_price->value = *(values+2);
        charge->value = *(values+3);
        *count_order = *(values+4);
        ++qty;
        ++base_price;
        ++disc_price;
        ++charge;
        ++count_order;
      }
      block.addedElements(found);
    }

  });
  printResultQ1(result.get());

  leaveQuery(nrThreads);
  return move(resources.query);
}
/*
std::unique_ptr<runtime::Query> q1_rof(Database& db, size_t nrThreads) {
  using namespace types;
  using namespace std;
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

  auto resources = initQuery(nrThreads);
  using hash = runtime::MurMurHash;
  using range = tbb::blocked_range<size_t>;
  const auto add = [](const size_t& a, const size_t& b) {return a + b;};

  tbb::enumerable_thread_specific<Hashmapx<types::Integer, tuple<Numeric<12, 2>, Numeric<12, 2>, Numeric<12, 4>, Numeric<12, 6>, int64_t>, hash, false>> hash_table;

  using group_t = typename decltype(hash_table)::value_type::Entry;

  /// Memory for materialized entries in hashmap
  tbb::enumerable_thread_specific<runtime::Stack<group_t>> entries;
  /// Memory for spilling hastable entries
  tbb::enumerable_thread_specific<runtime::PartitionedDeque<PARTITION_SIZE>> partitionedDeques;
  /// globally collect entry addresses
  tbb::enumerable_thread_specific<vector<group_t*>> entry_addrs;
  tbb::enumerable_thread_specific<vector<group_t*>> results_addrs;

  // local aggregation
  auto found2 = tbb::parallel_reduce(range(0, li.nrTuples, morselSize), 0, [&](const tbb::blocked_range<size_t>& r, const size_t& f) {
    auto found = f;
    bool exist=false;
    auto& ht = hash_table.local(exist);

    if(!exist) {
      ht.setSize(1024);

    }
    auto& partition = partitionedDeques.local(exist);
    if(!exist) {
      partition.postConstruct(nrThreads * 4, sizeof(group_t));
    }
    found += agg_gp_q1(r.begin(),r.end(),db,&ht,&partition,nullptr,nullptr);
    return found;
  },
                                     add);
//  cout << "local agg num = " << found2 << endl;
  auto& result = resources.query->result;
  auto retAttr = result->addAttribute("l_returnflag", sizeof(Char<1> ));
  auto statusAttr = result->addAttribute("l_linestatus", sizeof(Char<1> ));
  auto qtyAttr = result->addAttribute("sum_qty", sizeof(Numeric<12, 2> ));
  auto base_priceAttr = result->addAttribute("sum_base_price", sizeof(Numeric<12, 2> ));
  auto disc_priceAttr = result->addAttribute("sum_disc_price", sizeof(Numeric<12, 2> ));
  auto chargeAttr = result->addAttribute("sum_charge", sizeof(Numeric<12, 2> ));
  auto count_orderAttr = result->addAttribute("count_order", sizeof(int64_t));

  // global aggregation
  auto nrPartitions = partitionedDeques.begin()->getPartitions().size();

  tbb::parallel_for(0ul, nrPartitions, [&](auto partitionNr) {
    bool exist=false;
    auto& ht = hash_table.local(exist);
    if(!exist) {
      ht.setSize(1024);
    }
    ht.clear();
    auto& partition = partitionedDeques.local(exist);
    if(!exist) {
      partition.postConstruct(nrThreads * 4, sizeof(group_t));
    }

    //  collect entry addresses from a partition 

    auto& entry_addrs_ = entry_addrs.local();
    entry_addrs_.clear();
    auto& results_addrs_ = results_addrs.local();
    results_addrs_.clear();
    for (auto& deque : partitionedDeques) {
      auto& partition = deque.getPartitions()[partitionNr];
      for (auto chunk = partition.first; chunk; chunk = chunk->next) {
        for (auto value = chunk->template data<group_t>(), end = value + partition.size(chunk, sizeof(group_t)); value < end; value++) {
          entry_addrs_.push_back(value);
        }
      }
    }
    results_addrs_.resize(entry_addrs_.size());
    auto found = agg_gp_q1(0,entry_addrs_.size(),db,&ht,nullptr,(void**)&entry_addrs_[0],(void**)&results_addrs_[0]);
    if(found>0) {
      auto block = result->createBlock(found);
      auto ret = reinterpret_cast<Char<1>*>(block.data(retAttr));
      auto status = reinterpret_cast<Char<1>*>(block.data(statusAttr));
      auto qty = reinterpret_cast<Numeric<12, 2>*>(block.data(qtyAttr));
      auto base_price =
      reinterpret_cast<Numeric<12, 2>*>(block.data(base_priceAttr));
      auto disc_price =
      reinterpret_cast<Numeric<12, 4>*>(block.data(disc_priceAttr));
      auto charge = reinterpret_cast<Numeric<12, 6>*>(block.data(chargeAttr));
      auto count_order =
      reinterpret_cast<int64_t*>(block.data(count_orderAttr));

      for(int i=0;i<found;++i) {
       auto values = (uint64_t*) (((char*) results_addrs_[i]) + offsetof(group_t, v));
        *ret++ = Char<1>::build((((char*)&(results_addrs_[i]->k))+0));
        *status++ = Char<1>::build(((char*)&(results_addrs_[i]->k))+1);
        qty->value = *values;
        base_price->value = *(values+1);
        disc_price->value = *(values+2);
        charge->value = *(values+3);
        *count_order = *(values+4);
        ++qty;
        ++base_price;
        ++disc_price;
        ++charge;
        ++count_order;
      }
      block.addedElements(found);
    }

  });
  printResultQ1(result.get());

  leaveQuery(nrThreads);
  return move(resources.query);
}
*/
std::unique_ptr<Q1Builder::Q1> Q1Builder::getQuery() {
  using namespace vectorwise;
  auto result = Result();
  previous = result.resultWriter.shared.result->participate();

  auto r = make_unique<Q1>();
  auto lineitem = Scan("lineitem");
  Select(Expression().addOp(BF(primitives::sel_less_equal_Date_col_Date_val), Buffer(sel_date, sizeof(pos_t)), Column(lineitem, "l_shipdate"), Value(&r->c1)));
  Project().addExpression(
      Expression().addOp(conf.proj_sel_minus_int64_t_val_int64_t_col(), Buffer(sel_date), Buffer(result_proj_minus, sizeof(int64_t)), Value(&r->one),
                         Column(lineitem, "l_discount")).addOp(conf.proj_multiplies_sel_int64_t_col_int64_t_col(), Buffer(sel_date), Buffer(disc_price, sizeof(int64_t)),
                                                               Column(lineitem, "l_extendedprice"), Buffer(result_proj_minus, sizeof(int64_t)))).addExpression(
      Expression().addOp(conf.proj_sel_plus_int64_t_col_int64_t_val(), Buffer(sel_date), Buffer(result_proj_plus, sizeof(int64_t)), Column(lineitem, "l_tax"), Value(&r->one)).addOp(
          conf.proj_multiplies_int64_t_col_int64_t_col(), Buffer(charge, sizeof(int64_t)), Buffer(disc_price, sizeof(int64_t)), Buffer(result_proj_plus, sizeof(int64_t))));
  HashGroup().pushKeySelVec(Buffer(sel_date), Buffer(sel_date_grouped, sizeof(pos_t))).addKey(Column(lineitem, "l_returnflag"), Buffer(sel_date), primitives::hash_sel_Char_1_col,
                                                                                              primitives::keys_not_equal_sel_Char_1_col,
                                                                                              primitives::partition_by_key_sel_Char_1_col, Buffer(sel_date_grouped, sizeof(pos_t)),
                                                                                              primitives::scatter_sel_Char_1_col, primitives::keys_not_equal_row_Char_1_col,
                                                                                              primitives::partition_by_key_row_Char_1_col, primitives::scatter_sel_row_Char_1_col,
                                                                                              primitives::gather_val_Char_1_col, Buffer(returnflag, sizeof(Char_1))).addKey(
      Column(lineitem, "l_linestatus"), Buffer(sel_date), primitives::rehash_sel_Char_1_col, primitives::keys_not_equal_sel_Char_1_col, primitives::partition_by_key_sel_Char_1_col,
      Buffer(sel_date_grouped, sizeof(pos_t)), primitives::scatter_sel_Char_1_col, primitives::keys_not_equal_row_Char_1_col, primitives::partition_by_key_row_Char_1_col,
      primitives::scatter_sel_row_Char_1_col, primitives::gather_val_Char_1_col, Buffer(linestatus, sizeof(Char_1))).padToAlign(sizeof(types::Numeric<12, 4>)).addValue(
      Buffer(disc_price), primitives::aggr_init_plus_int64_t_col, primitives::aggr_plus_int64_t_col, primitives::aggr_row_plus_int64_t_col, primitives::gather_val_int64_t_col,
      Buffer(sum_disc_price, sizeof(types::Numeric<12, 4>))).addValue(Buffer(charge), primitives::aggr_init_plus_int64_t_col, primitives::aggr_plus_int64_t_col,
                                                                      primitives::aggr_row_plus_int64_t_col, primitives::gather_val_int64_t_col,
                                                                      Buffer(sum_charge, sizeof(types::Numeric<12, 4>))).addValue(Column(lineitem, "l_quantity"), Buffer(sel_date),
                                                                                                                                  primitives::aggr_init_plus_int64_t_col,
                                                                                                                                  primitives::aggr_sel_plus_int64_t_col,
                                                                                                                                  primitives::aggr_row_plus_int64_t_col,
                                                                                                                                  primitives::gather_val_int64_t_col,
                                                                                                                                  Buffer(sum_qty, sizeof(types::Numeric<12, 2>)))
      .addValue(Column(lineitem, "l_extendedprice"), Buffer(sel_date), primitives::aggr_init_plus_int64_t_col, primitives::aggr_sel_plus_int64_t_col,
                primitives::aggr_row_plus_int64_t_col, primitives::gather_val_int64_t_col, Buffer(sum_base_price, sizeof(types::Numeric<12, 2>))).addValue(
      Buffer(charge, sizeof(uint64_t)), primitives::aggr_init_plus_int64_t_col, primitives::aggr_count_star, primitives::aggr_row_plus_int64_t_col,
      primitives::gather_val_int64_t_col, Buffer(count_order, sizeof(uint64_t)));

  result.addValue("l_returnflag", Buffer(returnflag)).addValue("l_linestatus", Buffer(linestatus)).addValue("sum_qty", Buffer(sum_qty)).addValue("sum_base_price",
                                                                                                                                                 Buffer(sum_base_price)).addValue(
      "sum_disc_price", Buffer(sum_disc_price)).addValue("sum_charge", Buffer(sum_charge)).addValue("count_order", Buffer(count_order)).finalize();

  // TODO: add averages
  r->rootOp = popOperator();
  return r;
}

std::unique_ptr<runtime::Query> q1_vectorwise(Database& db, size_t nrThreads, size_t vectorSize) {
  using namespace vectorwise;
  WorkerGroup workers(nrThreads);
  vectorwise::SharedStateManager shared;

  std::unique_ptr<runtime::Query> result;
  workers.run([&]() {
    Q1Builder builder(db, shared, vectorSize);
    auto query = builder.getQuery();
    /* auto found = */query->rootOp->next();
    auto leader = barrier();
    if (leader)
    result = move(
        dynamic_cast<ResultWriter*>(query->rootOp.get())->shared.result);
  });
  printResultQ1(result.get()->result.get());

  return result;
}

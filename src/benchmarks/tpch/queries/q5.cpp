/* An adapted implementation of TPCH query q5 based on https://github.com/fzhedu/db-imv */

#include "benchmarks/tpch/Queries.hpp"
#include "common/runtime/Concurrency.hpp"
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
#include <iostream>

#include "config.h"            /* autoconf header */
#include "configs/base_configs.h"
#include "configs/eth_configs.h"
//#include "utils/io.h"


using namespace runtime;
using namespace std;
using vectorwise::primitives::Char_10;
using vectorwise::primitives::Char_25;
using vectorwise::primitives::hash_t;

// select
//  n_name,
//  sum(l_extendedprice * (1 - l_discount)) as revenue
// from
//   customer,
//   orders,
//   lineitem,
//   supplier,
//   nation,
//   region
// where
//   c_custkey = o_custkey
//   and l_orderkey = o_orderkey
//   and l_suppkey = s_suppkey
//   and c_nationkey = s_nationkey
//   and s_nationkey = n_nationkey
//   and n_regionkey = r_regionkey
//   and r_name = 'ASIA'
//   and o_orderdate >= date '1994-01-01'
//   and o_orderdate < date '1995-01-01'
// group by
//  n_name

using namespace runtime;
using namespace std;
/*NOVECTORIZE std::unique_ptr<runtime::Query> q5_hyper(Database& db,
                                                     size_t nrThreads) {

   const size_t morselSize = 10000;

   auto resources = initQuery(nrThreads);

   // --- constants
   auto c1 = types::Date::castString("1994-01-01");
   auto c2 = types::Date::castString("1995-01-01");
   string b = "ASIA";
   auto c3 = types::Char<25>::castString(b.data(), b.size());

   auto& su = db["supplier"];
   auto& re = db["region"];
   auto& na = db["nation"];
   auto& cu = db["customer"];
   auto& ord = db["orders"];
   auto& li = db["lineitem"];

   using hash = runtime::CRC32Hash;

   auto r_name = re["r_name"].data<types::Char<25>>();
   auto r_regionkey = re["r_regionkey"].data<types::Integer>();
   // --- select region and build ht
   Hashset<types::Integer, hash> ht1;
   tbb::enumerable_thread_specific<runtime::Stack<decltype(ht1)::Entry>>
       entries1;
   auto found1 = PARALLEL_SELECT(re.nrTuples, entries1, {
      if (r_name[i] == c3) {
         entries.emplace_back(ht1.hash(r_regionkey[i]), r_regionkey[i]);
         found++;
      }
   });
   ht1.setSize(found1);
   parallel_insert(entries1, ht1);

   // --- join on region and build ht
   auto n_regionkey = na["n_regionkey"].data<types::Integer>();
   auto n_nationkey = na["n_nationkey"].data<types::Integer>();
   auto n_name = na["n_name"].data<types::Char<25>>();
   Hashmapx<types::Integer, types::Char<25>, hash> ht2;
   tbb::enumerable_thread_specific<runtime::Stack<decltype(ht2)::Entry>>
       entries2;
   auto found2 = PARALLEL_SELECT(na.nrTuples, entries2, {
      if (ht1.contains(n_regionkey[i])) {
         entries.emplace_back(ht2.hash(n_nationkey[i]), n_nationkey[i],
                              n_name[i]);
         found++;
      }
   });
   ht2.setSize(found2);
   parallel_insert(entries2, ht2);

   // --- join on nation and build ht
   auto c_nationkey = cu["c_nationkey"].data<types::Integer>();
   auto c_custkey = cu["c_custkey"].data<types::Integer>();
   Hashmapx<types::Integer, std::tuple<types::Integer, types::Char<25>>, hash>
       ht3;
   tbb::enumerable_thread_specific<runtime::Stack<decltype(ht3)::Entry>>
       entries3;

   auto found3 = PARALLEL_SELECT(cu.nrTuples, entries3, {
      decltype(ht2)::value_type* v;
      if ((v = ht2.findOne(c_nationkey[i]))) {
         entries.emplace_back(ht3.hash(c_custkey[i]), c_custkey[i],
                              make_tuple(c_nationkey[i], *v));
         found++;
      }
   });
   ht3.setSize(found3);
   parallel_insert(entries3, ht3);

   // --- join on customer and build ht
   auto o_orderkey = ord["o_orderkey"].data<types::Integer>();
   auto o_orderdate = ord["o_orderdate"].data<types::Date>();
   auto o_custkey = ord["o_custkey"].data<types::Integer>();
   Hashmapx<types::Integer, std::tuple<types::Integer, types::Char<25>>, hash>
       ht4;
   tbb::enumerable_thread_specific<runtime::Stack<decltype(ht4)::Entry>>
       entries4;

   auto found4 = PARALLEL_SELECT(ord.nrTuples, entries4, {
      decltype(ht3)::value_type* v;
      if ((o_orderdate[i] < c2) & (o_orderdate[i] >= c1) &&
          (v = ht3.findOne(o_custkey[i]))) {
         entries.emplace_back(ht4.hash(o_orderkey[i]), o_orderkey[i], *v);
         found++;
      }
   });
   ht4.setSize(found4);
   parallel_insert(entries4, ht4);

   // --- build ht for supplier
   auto s_suppkey = su["s_suppkey"].data<types::Integer>();
   auto s_nationkey = su["s_nationkey"].data<types::Integer>();
   Hashset<std::tuple<types::Integer, types::Integer>, hash> ht5;
   tbb::enumerable_thread_specific<runtime::Stack<decltype(ht5)::Entry>>
       entries5;

   PARALLEL_SCAN(su.nrTuples, entries5, {
      auto key = make_tuple(s_suppkey[i], s_nationkey[i]);
      entries.emplace_back(ht5.hash(key), key);
   });

   ht5.setSize(su.nrTuples);
   parallel_insert(entries5, ht5);

   // --- join on customer and build ht
   auto l_orderkey = li["l_orderkey"].data<types::Integer>();
   auto l_suppkey = li["l_suppkey"].data<types::Integer>();
   auto l_extendedprice = li["l_extendedprice"].data<types::Numeric<12, 2>>();
   auto l_discount = li["l_discount"].data<types::Numeric<12, 2>>();

   const auto one = types::Numeric<12, 2>::castString("1.00");
   const auto zero = types::Numeric<12, 4>::castString("0.00");

   auto groupOp = make_GroupBy<types::Char<25>, types::Numeric<12, 4>, hash>(
       [](auto& acc, auto&& value) { acc += value; }, zero, nrThreads);

   // preaggregation
   tbb::parallel_for(
       tbb::blocked_range<size_t>(0, li.nrTuples, morselSize),
       [&](const tbb::blocked_range<size_t>& r) {
          auto groupLocals = groupOp.preAggLocals();

          for (size_t i = r.begin(), end = r.end(); i != end; ++i) {
             auto v = ht4.findOne(l_orderkey[i]);
             if (v) {
                auto suppkey = make_tuple(l_suppkey[i], get<0>(*v));
                if (ht5.contains(suppkey)) {
                   groupLocals.consume(get<1>(*v), l_extendedprice[i] *
                                                       (one - l_discount[i]));
                }
             }
          }
       });

   // --- output
   auto& result = resources.query->result;
   auto revAttr =
       result->addAttribute("revenue", sizeof(types::Numeric<12, 4>));
   auto nameAttr = result->addAttribute("n_name", sizeof(types::Char<25>));

   groupOp.forallGroups([&](auto& groups) {
      // write aggregates to result
      auto block = result->createBlock(groups.size());
      auto name = reinterpret_cast<types::Char<25>*>(block.data(nameAttr));
      auto rev = reinterpret_cast<types::Numeric<12, 4>*>(block.data(revAttr));
      for (auto block : groups)
         for (auto& group : block) {
            *name++ = group.k;
            *rev++ = group.v;
         }
      block.addedElements(groups.size());
   });

   leaveQuery(nrThreads);
   return move(resources.query);
}*/
//the full query
/*unique_ptr<Q5Builder::Q5> Q5Builder::getQuery() {
   using namespace vectorwise;
   auto result = Result();
   previous = result.resultWriter.shared.result->participate();
   auto r = make_unique<Q5>();
   auto supplier = Scan("supplier");
   auto region = Scan("region");
   Select(
       Expression().addOp(BF(primitives::sel_equal_to_Char_25_col_Char_25_val),
                          Buffer(sel_region, sizeof(pos_t)), //
                          Column(region, "r_name"),          //
                          Value(&r->c3)));
   auto nation = Scan("nation");
   HashJoin(Buffer(join_reg_nat, sizeof(pos_t)), conf.joinAll())
       .addBuildKey(Column(region, "r_regionkey"), Buffer(sel_region),
                    conf.hash_sel_int32_t_col(),
                    primitives::scatter_sel_int32_t_col)
       .addProbeKey(Column(nation, "n_regionkey"), conf.hash_int32_t_col(),
                    primitives::keys_equal_int32_t_col);
   auto customer = Scan("customer");
   HashJoin(Buffer(join_cust, sizeof(pos_t)), conf.joinAll())
       .addBuildKey(Column(nation, "n_nationkey"), Buffer(join_reg_nat),
                    conf.hash_sel_int32_t_col(),
                    primitives::scatter_sel_int32_t_col)
       .addProbeKey(Column(customer, "c_nationkey"),
                    conf.hash_int32_t_col(),
                    primitives::keys_equal_int32_t_col)
       .addBuildValue(Column(nation, "n_name"), Buffer(join_reg_nat),
                      primitives::scatter_sel_Char_25_col,
                      Buffer(n_name, sizeof(Char_25)),
                      primitives::gather_col_Char_25_col);
   auto orders = Scan("orders");
   Select(Expression()
              .addOp(BF(primitives::sel_less_Date_col_Date_val),
                     Buffer(sel_ord, sizeof(pos_t)),
                     Column(orders, "o_orderdate"), Value(&r->c2))
              .addOp(BF(primitives::selsel_greater_equal_Date_col_Date_val),
                     Buffer(sel_ord, sizeof(pos_t)),
                     Buffer(sel_ord2, sizeof(pos_t)),
                     Column(orders, "o_orderdate"), Value(&r->c1)));
   HashJoin(Buffer(join_ord, sizeof(pos_t)), conf.joinAll())
       .setProbeSelVector(Buffer(sel_ord2), conf.joinSel())
       .addBuildKey(Column(customer, "c_custkey"), Buffer(join_cust),
                    conf.hash_sel_int32_t_col(),
                    primitives::scatter_sel_int32_t_col)
       .addProbeKey(Column(orders, "o_custkey"), Buffer(sel_ord2),
                    conf.hash_sel_int32_t_col(),
                    primitives::keys_equal_int32_t_col)
       .addBuildValue(Column(customer, "c_nationkey"), Buffer(join_cust),
                      primitives::scatter_sel_int32_t_col,
                      Buffer(join_ord_nationkey, sizeof(int32_t)),
                      primitives::gather_col_int32_t_col)
       .addBuildValue(Buffer(n_name), primitives::scatter_Char_25_col,
                      Buffer(n_name2, sizeof(Char_25)),
                      primitives::gather_col_Char_25_col);
   auto lineitem = Scan("lineitem");
   HashJoin(Buffer(join_line, sizeof(pos_t)), conf.joinAll())
       .addBuildKey(Column(orders, "o_orderkey"), Buffer(join_ord),
                    conf.hash_sel_int32_t_col(),
                    primitives::scatter_sel_int32_t_col)
       .addProbeKey(Column(lineitem, "l_orderkey"),
                    conf.hash_int32_t_col(),
                    primitives::keys_equal_int32_t_col)
       .addBuildValue(Buffer(join_ord_nationkey),
                      primitives::scatter_int32_t_col,
                      Buffer(join_line_nationkey, sizeof(int32_t)),
                      primitives::gather_col_int32_t_col)
       .addBuildValue(Buffer(n_name2), primitives::scatter_Char_25_col,
                      Buffer(n_name, sizeof(Char_25)),
                      primitives::gather_col_Char_25_col);
   HashJoin(Buffer(join_supp, sizeof(pos_t)), conf.joinAll())
       .addBuildKey(Column(supplier, "s_nationkey"),
                     conf.hash_int32_t_col(),
                    primitives::scatter_int32_t_col)
       .addBuildKey(Column(supplier, "s_suppkey"),
                    conf.rehash_int32_t_col(),
                    primitives::scatter_int32_t_col)
       .addProbeKey(Buffer(join_line_nationkey),  conf.hash_int32_t_col(),
                    primitives::keys_equal_int32_t_col)
       .pushProbeSelVector(Buffer(join_line),
                           Buffer(join_supp_line, sizeof(pos_t)))
       .addProbeKey(Column(lineitem, "l_suppkey"), Buffer(join_line),
                    conf.rehash_sel_int32_t_col(),
                    Buffer(join_supp_line, sizeof(pos_t)),
                    primitives::keys_equal_int32_t_col);
   Project().addExpression(
       Expression()
           .addOp(primitives::proj_sel_minus_int64_t_val_int64_t_col,
                  Buffer(join_supp_line),
                  Buffer(result_proj_minus, sizeof(int64_t)), Value(&r->one),
                  Column(lineitem, "l_discount"))
           .addOp(primitives::proj_multiplies_sel_int64_t_col_int64_t_col,
                  Buffer(join_supp_line),
                  Buffer(result_project, sizeof(int64_t)),
                  Column(lineitem, "l_extendedprice"),
                  Buffer(result_proj_minus, sizeof(int64_t))));
   HashGroup()
       .pushKeySelVec(Buffer(join_supp), Buffer(selScat, sizeof(pos_t)))
       .addKey(
           Buffer(n_name), Buffer(join_supp), primitives::hash_sel_Char_25_col,
           primitives::keys_not_equal_sel_Char_25_col,
           primitives::partition_by_key_sel_Char_25_col,
           Buffer(selScat, sizeof(pos_t)), primitives::scatter_sel_Char_25_col,
           primitives::keys_not_equal_row_Char_25_col,
           primitives::partition_by_key_row_Char_25_col,
           primitives::scatter_sel_row_Char_25_col,
           primitives::gather_val_Char_25_col, Buffer(n_name))
       .addValue(Buffer(result_project), primitives::aggr_init_plus_int64_t_col,
                 primitives::aggr_plus_int64_t_col,
                 primitives::aggr_row_plus_int64_t_col,
                 primitives::gather_val_int64_t_col,
                 Buffer(sum, sizeof(int64_t)));

   result //
       .addValue("n_name", Buffer(n_name))
       .addValue("revenue", Buffer(sum))
       .finalize();

   r->rootOp = popOperator();
   assert(operatorStack.size() == 0);
   return r;
}*/
/*
unique_ptr<Q5Builder::Q5> Q5Builder::getQuery() {
   using namespace vectorwise;
   auto result = Result();
   previous = result.resultWriter.shared.result->participate();
   auto r = make_unique<Q5>();
   auto region = Scan("region");
   Select(
       Expression().addOp(BF(primitives::sel_equal_to_Char_25_col_Char_25_val),
                          Buffer(sel_region, sizeof(pos_t)), //
                          Column(region, "r_name"),          //
                          Value(&r->c3)));
   result //
       .addValue("sel_region", Buffer(sel_region))
       .finalize();

   r->rootOp = popOperator();
   assert(operatorStack.size() == 0);
   return r;
}

void printResultQ5(BlockRelation* result) {
  using namespace types;  
  size_t found = 0;
  auto selRegionAttr = result->getAttribute("sel_region");
  uint32_t* selRegion;
  for (auto& block : *result) {
    auto elementsInBlock = block.size();
    found += elementsInBlock;
    selRegion = reinterpret_cast<uint32_t*>(block.data(selRegionAttr));
    for (size_t i = 0; i < elementsInBlock; ++i) {
      cout << selRegion[i] << endl;
    }
  }
  cout << "sel_region total results number = " << found << endl;

  materialize_one_relation<RELATION_KEY_TYPE, RELATION_PAYLOAD_TYPE>(selRegion, found);   
}


unique_ptr<Q5Builder::Q5> Q5Builder::getQuery() {
   using namespace vectorwise;
   auto result = Result();
   previous = result.resultWriter.shared.result->participate();
   auto r = make_unique<Q5>();

   auto nation = Scan("nation");

   result //
       .addValue("n_regionkey", Column(nation, "n_regionkey"))
       .finalize();

   r->rootOp = popOperator();
   assert(operatorStack.size() == 0);
   return r;
}

void printResultQ5(BlockRelation* result) {
  using namespace types;  
  size_t found = 0;
  auto nRegionKeyAttr = result->getAttribute("n_regionkey");
  uint32_t* nRegionKey;
  for (auto& block : *result) {
    auto elementsInBlock = block.size();
    found += elementsInBlock;
    nRegionKey = reinterpret_cast<uint32_t*>(block.data(nRegionKeyAttr));
    //for (size_t i = 0; i < elementsInBlock; ++i) {
    //  cout << nRegionKey[i] << endl;
    //}
  }
  cout << "n_regionkey total results number = " << found << endl;

  materialize_one_relation<RELATION_KEY_TYPE, RELATION_PAYLOAD_TYPE>(nRegionKey, found);   
}
*/

/*
unique_ptr<Q5Builder::Q5> Q5Builder::getQuery() {
   using namespace vectorwise;
   auto result = Result();
   previous = result.resultWriter.shared.result->participate();
   auto r = make_unique<Q5>();
   auto region = Scan("region");
   Select(
       Expression().addOp(BF(primitives::sel_equal_to_Char_25_col_Char_25_val),
                          Buffer(sel_region, sizeof(pos_t)), //
                          Column(region, "r_name"),          //
                          Value(&r->c3)));
   auto nation = Scan("nation");
   HashJoin(Buffer(join_reg_nat, sizeof(pos_t)), conf.joinAll())
       .addBuildKey(Column(region, "r_regionkey"), Buffer(sel_region),
                    conf.hash_sel_int32_t_col(),
                    primitives::scatter_sel_int32_t_col)
       .addProbeKey(Column(nation, "n_regionkey"), conf.hash_int32_t_col(),
                    primitives::keys_equal_int32_t_col);

   result //
       .addValue("join_reg_nat", Buffer(join_reg_nat))
       .finalize();

   r->rootOp = popOperator();
   assert(operatorStack.size() == 0);
   return r;
}

void printResultQ5(BlockRelation* result) {
  using namespace types;  
  size_t found = 0;
  auto regNatAttr = result->getAttribute("join_reg_nat");
  uint32_t* regNat;
  for (auto& block : *result) {
    auto elementsInBlock = block.size();
    found += elementsInBlock;
    regNat = reinterpret_cast<uint32_t*>(block.data(regNatAttr));
    //for (size_t i = 0; i < elementsInBlock; ++i) {
    //  cout << regNat[i] << endl;
    //}
  }
  cout << "join_reg_nat total results number = " << found << endl;

  //materialize_one_relation<RELATION_KEY_TYPE, RELATION_PAYLOAD_TYPE>(regNat, found);   
}
*/


unique_ptr<Q5Builder::Q5> Q5Builder::getQuery() {
   using namespace vectorwise;
   auto result = Result();
   previous = result.resultWriter.shared.result->participate();
   auto r = make_unique<Q5>();
   auto customer = Scan("customer");
  
   result //
       .addValue("c_nationkey", Column(customer, "c_nationkey"))
       .finalize();

   r->rootOp = popOperator();
   assert(operatorStack.size() == 0);
   return r;
}

void printResultQ5(BlockRelation* result) {
  using namespace types;  
  size_t found = 0;
  auto nationKeyAttr = result->getAttribute("c_nationkey");
  uint32_t* nationKey;
  for (auto& block : *result) {
    auto elementsInBlock = block.size();
    found += elementsInBlock;
    nationKey = reinterpret_cast<uint32_t*>(block.data(nationKeyAttr));
    //for (size_t i = 0; i < elementsInBlock; ++i) {
    //  cout << nationKey[i] << endl;
    //}
  }
  cout << "c_nationkey total results number = " << found << endl;

  //materialize_one_relation<RELATION_KEY_TYPE, RELATION_PAYLOAD_TYPE>(nationKey, found);   
}

/*
unique_ptr<Q5Builder::Q5> Q5Builder::getQuery() {
   using namespace vectorwise;
   auto result = Result();
   previous = result.resultWriter.shared.result->participate();
   auto r = make_unique<Q5>();
   auto region = Scan("region");
   Select(
       Expression().addOp(BF(primitives::sel_equal_to_Char_25_col_Char_25_val),
                          Buffer(sel_region, sizeof(pos_t)), //
                          Column(region, "r_name"),          //
                          Value(&r->c3)));
   auto nation = Scan("nation");
   HashJoin(Buffer(join_reg_nat, sizeof(pos_t)), conf.joinAll())
       .addBuildKey(Column(region, "r_regionkey"), Buffer(sel_region),
                    conf.hash_sel_int32_t_col(),
                    primitives::scatter_sel_int32_t_col)
       .addProbeKey(Column(nation, "n_regionkey"), conf.hash_int32_t_col(),
                    primitives::keys_equal_int32_t_col);
   auto customer = Scan("customer");
   HashJoin(Buffer(join_cust, sizeof(pos_t)), conf.joinAll())
       .addBuildKey(Column(nation, "n_nationkey"), Buffer(join_reg_nat),
                    conf.hash_sel_int32_t_col(),
                    primitives::scatter_sel_int32_t_col)
       .addProbeKey(Column(customer, "c_nationkey"),
                    conf.hash_int32_t_col(),
                    primitives::keys_equal_int32_t_col)
       .addBuildValue(Column(nation, "n_name"), Buffer(join_reg_nat),
                      primitives::scatter_sel_Char_25_col,
                      Buffer(n_name, sizeof(Char_25)),
                      primitives::gather_col_Char_25_col);

   result //
       .addValue("join_cust", Buffer(join_cust))
       .finalize();

   r->rootOp = popOperator();
   assert(operatorStack.size() == 0);
   return r;
}

void printResultQ5(BlockRelation* result) {
  using namespace types;  
  size_t found = 0;
  auto joinCustAttr = result->getAttribute("join_cust");
  uint32_t* joinCust;
  for (auto& block : *result) {
    auto elementsInBlock = block.size();
    found += elementsInBlock;
    joinCust = reinterpret_cast<uint32_t*>(block.data(joinCustAttr));
    //for (size_t i = 0; i < elementsInBlock; ++i) {
    //  cout << joinCust[i] << endl;
    //}
  }
  cout << "join_cust total results number = " << found << endl;

  materialize_one_relation<RELATION_KEY_TYPE, RELATION_PAYLOAD_TYPE>(joinCust, found);   
}



unique_ptr<Q5Builder::Q5> Q5Builder::getQuery() {
   using namespace vectorwise;
   auto result = Result();
   previous = result.resultWriter.shared.result->participate();
   auto r = make_unique<Q5>();
   auto orders = Scan("orders");
   Select(Expression()
              .addOp(BF(primitives::sel_less_Date_col_Date_val),
                     Buffer(sel_ord, sizeof(pos_t)),
                     Column(orders, "o_orderdate"), Value(&r->c2))
              .addOp(BF(primitives::selsel_greater_equal_Date_col_Date_val),
                     Buffer(sel_ord, sizeof(pos_t)),
                     Buffer(sel_ord2, sizeof(pos_t)),
                     Column(orders, "o_orderdate"), Value(&r->c1)));
   result //
       .addValue("sel_ord2", Buffer(sel_ord2, sizeof(pos_t)))
       .finalize();

   r->rootOp = popOperator();
   assert(operatorStack.size() == 0);
   return r;
}

void printResultQ5(BlockRelation* result) {
  using namespace types;  
  size_t found = 0;
  auto selOrd2Attr = result->getAttribute("sel_ord2");
  uint32_t* selOrd2;
  for (auto& block : *result) {
    auto elementsInBlock = block.size();
    found += elementsInBlock;
    selOrd2 = reinterpret_cast<uint32_t*>(block.data(selOrd2Attr));
    //for (size_t i = 0; i < elementsInBlock; ++i) {
    //  cout << selOrd2[i] << endl;
    //}
  }
  cout << "sel_ord2 total results number = " << found << endl;

  materialize_one_relation<RELATION_KEY_TYPE, RELATION_PAYLOAD_TYPE>(selOrd2, found);   
}

unique_ptr<Q5Builder::Q5> Q5Builder::getQuery() {
   using namespace vectorwise;
   auto result = Result();
   previous = result.resultWriter.shared.result->participate();
   auto r = make_unique<Q5>();
   auto region = Scan("region");
   Select(
       Expression().addOp(BF(primitives::sel_equal_to_Char_25_col_Char_25_val),
                          Buffer(sel_region, sizeof(pos_t)), //
                          Column(region, "r_name"),          //
                          Value(&r->c3)));
   auto nation = Scan("nation");
   HashJoin(Buffer(join_reg_nat, sizeof(pos_t)), conf.joinAll())
       .addBuildKey(Column(region, "r_regionkey"), Buffer(sel_region),
                    conf.hash_sel_int32_t_col(),
                    primitives::scatter_sel_int32_t_col)
       .addProbeKey(Column(nation, "n_regionkey"), conf.hash_int32_t_col(),
                    primitives::keys_equal_int32_t_col);
   auto customer = Scan("customer");
   HashJoin(Buffer(join_cust, sizeof(pos_t)), conf.joinAll())
       .addBuildKey(Column(nation, "n_nationkey"), Buffer(join_reg_nat),
                    conf.hash_sel_int32_t_col(),
                    primitives::scatter_sel_int32_t_col)
       .addProbeKey(Column(customer, "c_nationkey"),
                    conf.hash_int32_t_col(),
                    primitives::keys_equal_int32_t_col)
       .addBuildValue(Column(nation, "n_name"), Buffer(join_reg_nat),
                      primitives::scatter_sel_Char_25_col,
                      Buffer(n_name, sizeof(Char_25)),
                      primitives::gather_col_Char_25_col);
   auto orders = Scan("orders");
   Select(Expression()
              .addOp(BF(primitives::sel_less_Date_col_Date_val),
                     Buffer(sel_ord, sizeof(pos_t)),
                     Column(orders, "o_orderdate"), Value(&r->c2))
              .addOp(BF(primitives::selsel_greater_equal_Date_col_Date_val),
                     Buffer(sel_ord, sizeof(pos_t)),
                     Buffer(sel_ord2, sizeof(pos_t)),
                     Column(orders, "o_orderdate"), Value(&r->c1)));
   HashJoin(Buffer(join_ord, sizeof(pos_t)), conf.joinAll())
       .setProbeSelVector(Buffer(sel_ord2), conf.joinSel())
       .addBuildKey(Column(customer, "c_custkey"), Buffer(join_cust),
                    conf.hash_sel_int32_t_col(),
                    primitives::scatter_sel_int32_t_col)
       .addProbeKey(Column(orders, "o_custkey"), Buffer(sel_ord2),
                    conf.hash_sel_int32_t_col(),
                    primitives::keys_equal_int32_t_col)
       .addBuildValue(Column(customer, "c_nationkey"), Buffer(join_cust),
                      primitives::scatter_sel_int32_t_col,
                      Buffer(join_ord_nationkey, sizeof(int32_t)),
                      primitives::gather_col_int32_t_col)
       .addBuildValue(Buffer(n_name), primitives::scatter_Char_25_col,
                      Buffer(n_name2, sizeof(Char_25)),
                      primitives::gather_col_Char_25_col);
   

   result //
       .addValue("join_ord", Buffer(join_ord))
       .finalize();

   r->rootOp = popOperator();
   assert(operatorStack.size() == 0);
   return r;
}

void printResultQ5(BlockRelation* result) {
  using namespace types;  
  size_t found = 0;
  auto joinOrdAttr = result->getAttribute("join_ord");
  uint32_t* joinOrd;
  for (auto& block : *result) {
    auto elementsInBlock = block.size();
    found += elementsInBlock;
    joinOrd = reinterpret_cast<uint32_t*>(block.data(joinOrdAttr));
    //for (size_t i = 0; i < elementsInBlock; ++i) {
    //  cout << joinOrd[i] << endl;
    //}
  }
  cout << "join_ord total results number = " << found << endl;

  materialize_one_relation<RELATION_KEY_TYPE, RELATION_PAYLOAD_TYPE>(joinOrd, found);   
}


unique_ptr<Q5Builder::Q5> Q5Builder::getQuery() {
   using namespace vectorwise;
   auto result = Result();
   previous = result.resultWriter.shared.result->participate();
   auto r = make_unique<Q5>();
   
   auto lineitem = Scan("lineitem");

   result //
       .addValue("l_orderkey", Column(lineitem, "l_orderkey"))
       .finalize();

   r->rootOp = popOperator();
   assert(operatorStack.size() == 0);
   return r;
}

void printResultQ5(BlockRelation* result) {
  using namespace types;  
  size_t found = 0;
  auto orderKeyAttr = result->getAttribute("l_orderkey");
  uint32_t* orderKey;
  for (auto& block : *result) {
    auto elementsInBlock = block.size();
    found += elementsInBlock;
    orderKey = reinterpret_cast<uint32_t*>(block.data(orderKeyAttr));
    //for (size_t i = 0; i < elementsInBlock; ++i) {
    //  cout << orderKey[i] << endl;
    //}
  }
  cout << "l_orderkey total results number = " << found << endl;

  materialize_one_relation<RELATION_KEY_TYPE, RELATION_PAYLOAD_TYPE>(orderKey, found);   
}

unique_ptr<Q5Builder::Q5> Q5Builder::getQuery() {
   using namespace vectorwise;
   auto result = Result();
   previous = result.resultWriter.shared.result->participate();
   auto r = make_unique<Q5>();
   auto region = Scan("region");
   Select(
       Expression().addOp(BF(primitives::sel_equal_to_Char_25_col_Char_25_val),
                          Buffer(sel_region, sizeof(pos_t)), //
                          Column(region, "r_name"),          //
                          Value(&r->c3)));
   auto nation = Scan("nation");
   HashJoin(Buffer(join_reg_nat, sizeof(pos_t)), conf.joinAll())
       .addBuildKey(Column(region, "r_regionkey"), Buffer(sel_region),
                    conf.hash_sel_int32_t_col(),
                    primitives::scatter_sel_int32_t_col)
       .addProbeKey(Column(nation, "n_regionkey"), conf.hash_int32_t_col(),
                    primitives::keys_equal_int32_t_col);
   auto customer = Scan("customer");
   HashJoin(Buffer(join_cust, sizeof(pos_t)), conf.joinAll())
       .addBuildKey(Column(nation, "n_nationkey"), Buffer(join_reg_nat),
                    conf.hash_sel_int32_t_col(),
                    primitives::scatter_sel_int32_t_col)
       .addProbeKey(Column(customer, "c_nationkey"),
                    conf.hash_int32_t_col(),
                    primitives::keys_equal_int32_t_col)
       .addBuildValue(Column(nation, "n_name"), Buffer(join_reg_nat),
                      primitives::scatter_sel_Char_25_col,
                      Buffer(n_name, sizeof(Char_25)),
                      primitives::gather_col_Char_25_col);
   auto orders = Scan("orders");
   Select(Expression()
              .addOp(BF(primitives::sel_less_Date_col_Date_val),
                     Buffer(sel_ord, sizeof(pos_t)),
                     Column(orders, "o_orderdate"), Value(&r->c2))
              .addOp(BF(primitives::selsel_greater_equal_Date_col_Date_val),
                     Buffer(sel_ord, sizeof(pos_t)),
                     Buffer(sel_ord2, sizeof(pos_t)),
                     Column(orders, "o_orderdate"), Value(&r->c1)));
   HashJoin(Buffer(join_ord, sizeof(pos_t)), conf.joinAll())
       .setProbeSelVector(Buffer(sel_ord2), conf.joinSel())
       .addBuildKey(Column(customer, "c_custkey"), Buffer(join_cust),
                    conf.hash_sel_int32_t_col(),
                    primitives::scatter_sel_int32_t_col)
       .addProbeKey(Column(orders, "o_custkey"), Buffer(sel_ord2),
                    conf.hash_sel_int32_t_col(),
                    primitives::keys_equal_int32_t_col)
       .addBuildValue(Column(customer, "c_nationkey"), Buffer(join_cust),
                      primitives::scatter_sel_int32_t_col,
                      Buffer(join_ord_nationkey, sizeof(int32_t)),
                      primitives::gather_col_int32_t_col)
       .addBuildValue(Buffer(n_name), primitives::scatter_Char_25_col,
                      Buffer(n_name2, sizeof(Char_25)),
                      primitives::gather_col_Char_25_col);
   auto lineitem = Scan("lineitem");
   HashJoin(Buffer(join_line, sizeof(pos_t)), conf.joinAll())
       .addBuildKey(Column(orders, "o_orderkey"), Buffer(join_ord),
                    conf.hash_sel_int32_t_col(),
                    primitives::scatter_sel_int32_t_col)
       .addProbeKey(Column(lineitem, "l_orderkey"),
                    conf.hash_int32_t_col(),
                    primitives::keys_equal_int32_t_col)
       .addBuildValue(Buffer(join_ord_nationkey),
                      primitives::scatter_int32_t_col,
                      Buffer(join_line_nationkey, sizeof(int32_t)),
                      primitives::gather_col_int32_t_col)
       .addBuildValue(Buffer(n_name2), primitives::scatter_Char_25_col,
                      Buffer(n_name, sizeof(Char_25)),
                      primitives::gather_col_Char_25_col);
   result //
       .addValue("join_line", Buffer(join_line))
       .finalize();

   r->rootOp = popOperator();
   assert(operatorStack.size() == 0);
   return r;
}

void printResultQ5(BlockRelation* result) {
  using namespace types;  
  size_t found = 0;
  auto joinLineAttr = result->getAttribute("join_line");
  uint32_t* joinLine;
  for (auto& block : *result) {
    auto elementsInBlock = block.size();
    found += elementsInBlock;
    joinLine = reinterpret_cast<uint32_t*>(block.data(joinLineAttr));
    //for (size_t i = 0; i < elementsInBlock; ++i) {
    //  cout << joinLine[i] << endl;
    //}
  }
  cout << "join_line total results number = " << found << endl;

  materialize_one_relation<RELATION_KEY_TYPE, RELATION_PAYLOAD_TYPE>(joinLine, found);   
}


unique_ptr<Q5Builder::Q5> Q5Builder::getQuery() {
   using namespace vectorwise;
   auto result = Result();
   previous = result.resultWriter.shared.result->participate();
   auto r = make_unique<Q5>();
   auto supplier = Scan("supplier");
   
   result //
       .addValue("s_nationkey", Column(supplier, "s_nationkey"))
       .finalize();

   r->rootOp = popOperator();
   assert(operatorStack.size() == 0);
   return r;
}

void printResultQ5(BlockRelation* result) {
  using namespace types;  
  size_t found = 0;
  auto nationKeyAttr = result->getAttribute("s_nationkey");
  uint32_t* nationKey;
  for (auto& block : *result) {
    auto elementsInBlock = block.size();
    found += elementsInBlock;
    nationKey = reinterpret_cast<uint32_t*>(block.data(nationKeyAttr));
    //for (size_t i = 0; i < elementsInBlock; ++i) {
    //  cout << nationKey[i] << endl;
    //}
  }
  cout << "s_nationkey total results number = " << found << endl;

  materialize_one_relation<RELATION_KEY_TYPE, RELATION_PAYLOAD_TYPE>(nationKey, found);   
}
*/


std::unique_ptr<runtime::Query> q5_vectorwise(Database& db, size_t nrThreads,
                                              size_t vectorSize) {
   using namespace vectorwise;

   WorkerGroup workers(nrThreads);
   vectorwise::SharedStateManager shared;
   std::unique_ptr<runtime::Query> result;
   workers.run([&]() {
      Q5Builder builder(db, shared, vectorSize);
      auto resources = builder.getQuery();
      /* auto found = */ resources->rootOp->next();
      auto leader = barrier();
      if (leader)
         result = move(dynamic_cast<ResultWriter*>(resources->rootOp.get())
                           ->shared.result);
   });

   printResultQ5(result.get()->result.get());
   return result;
}

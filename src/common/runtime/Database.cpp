#include "common/runtime/Database.hpp"

#include <malloc.h>
#include <string.h>

#include "common/runtime/Concurrency.hpp"
#include "common/runtime/Types.hpp"
#include <cstdlib>

namespace runtime {

Attribute::Attribute(std::string n, std::unique_ptr<Type> t)
    : name(n),
      type(move(t)) {
}

Attribute& Relation::operator[](std::string key) {
  auto att = attributes.find(key);
  if (att != attributes.end())
    return att->second;
  else
    throw std::range_error("Unknown attribute " + key + " in relation " + name);
}

Attribute& Relation::insert(std::string name, std::unique_ptr<Type> t) {
  return attributes.emplace(std::piecewise_construct, std::forward_as_tuple(name), std::forward_as_tuple(name, std::move(t))).first->second;
}

bool Database::hasRelation(std::string name) {
  return relations.find(name) != relations.end();
}

Relation& Database::operator[](std::string key) {
  return relations[key];
}
;

BlockRelation::Block BlockRelation::createBlock(size_t minNrElements) {
  auto elements = std::max(minBlockSize, minNrElements);
  auto a = this_worker->allocator.allocate(sizeof(BlockHeader) + elements * currentAttributeSize);
  auto header = new (a) BlockHeader(elements);
  {
    std::lock_guard<std::mutex> lock(insertMutex);
    blocks.push_back(header);
  }
  return Block(this, header);
}

BlockRelation::Attribute BlockRelation::addAttribute(std::string name, size_t elementSize) {
  Attribute attr = attributes.size();
  attributes.emplace_back(elementSize, currentAttributeSize);
  attributeNames.emplace(name, attr);
  currentAttributeSize += elementSize;
  return attr;
}

BlockRelation::~BlockRelation() {
  // Block storage is freed by queries own memory pool
  // Thus, no freeing needed here.
}
uint64_t runtime::ModifyCol::get_id(std::string str) {
  auto it = str2id.find(str);
  if (it == str2id.end()) {
    throw std::range_error("Unknown str " + str + " in string 2 id ");
  }
  return it->second;
}

std::string& runtime::ModifyCol::get_str(uint64_t id) {
  auto it = id2str.find(id);
  if (it == id2str.end()) {
    throw std::range_error("Unknown id  in id 2 string ");
  }
  return it->second;
}

uint64_t* runtime::ModifyCol::get_addr(std::string col) {
  auto it = col_addresses.find(col);
  if (it == col_addresses.end()) {
    throw std::range_error("Unknown column " + std::string(col) + " in modified");
  }
  return it->second;
}

void runtime::Database::modifyDB() {
  uint64_t id = 0;
  std::string col_name;
  std::string cell;
  char temp[32];

  auto& tb_cust = relations["customer"];
  col_name = "c_region";
  auto ptr = tb_cust[col_name].data<types::Char<12>>();
  uint64_t* col = (uint64_t*) memalign(64, tb_cust.nrTuples * sizeof(uint64_t));
  modify.col_addresses[col_name] = col;
  for (uint64_t i = 0; i < tb_cust.nrTuples; ++i) {
    memcpy(temp,ptr[i].value,ptr[i].length());
    temp[ptr[i].length()]='\0';
    cell = std::string(temp);
    auto iter = modify.str2id.find(cell);
    if (iter == modify.str2id.end()) {
      modify.str2id[cell] = id;
      modify.id2str[id] = cell;
      col[i] = id;
      ++id;
    } else {
      col[i] = iter->second;
    }
  }

  col_name = "c_city";
  auto cust_city = tb_cust[col_name].data<types::Char<10>>();
  col = (uint64_t*) memalign(64, tb_cust.nrTuples * sizeof(uint64_t));
  modify.col_addresses[col_name] = col;
  for (uint64_t i = 0; i < tb_cust.nrTuples; ++i) {
    memcpy(temp,cust_city[i].value,cust_city[i].length());
    temp[cust_city[i].length()]='\0';
    cell = std::string(temp);
    auto iter = modify.str2id.find(cell);
    if (iter == modify.str2id.end()) {
      modify.str2id[cell] = id;
      modify.id2str[id] = cell;
      col[i] = id;
      ++id;
    } else {
      col[i] = iter->second;
    }
  }
  col_name = "c_nation";
  auto cust_nation = tb_cust[col_name].data<types::Char<15>>();
  col = (uint64_t*) memalign(64, tb_cust.nrTuples * sizeof(uint64_t));
  modify.col_addresses[col_name] = col;
  for (uint64_t i = 0; i < tb_cust.nrTuples; ++i) {
    memcpy(temp,cust_nation[i].value,cust_nation[i].length());
    temp[cust_nation[i].length()]='\0';
    cell = std::string(temp);
    auto iter = modify.str2id.find(cell);
    if (iter == modify.str2id.end()) {
      modify.str2id[cell] = id;
      modify.id2str[id] = cell;
      col[i] = id;
      ++id;
    } else {
      col[i] = iter->second;
    }
  }
  std::cout<<"customer is ok-----------size = "<<modify.str2id.size()<<std::endl;

  auto& tb_supp = relations["supplier"];
  col_name = "s_region";
  auto supp_region = tb_supp[col_name].data<types::Char<12>>();
  col = (uint64_t*) memalign(64, tb_supp.nrTuples * sizeof(uint64_t));
  modify.col_addresses[col_name] = col;
  for (uint64_t i = 0; i < tb_supp.nrTuples; ++i) {
    memcpy(temp,supp_region[i].value,supp_region[i].length());
    temp[supp_region[i].length()]='\0';
    cell = std::string(temp);
    auto iter = modify.str2id.find(cell);
    if (iter == modify.str2id.end()) {
      modify.str2id[cell] = id;
      modify.id2str[id] = cell;
      col[i] = id;
      ++id;
    } else {
      col[i] = iter->second;
    }
  }
  col_name = "s_city";
  auto supp_city = tb_supp[col_name].data<types::Char<10>>();
  col = (uint64_t*) memalign(64, tb_supp.nrTuples * sizeof(uint64_t));
  modify.col_addresses[col_name] = col;
  for (uint64_t i = 0; i < tb_supp.nrTuples; ++i) {
    memcpy(temp,supp_city[i].value,supp_city[i].length());
    temp[supp_city[i].length()]='\0';
    cell = std::string(temp);
    auto iter = modify.str2id.find(cell);
    if (iter == modify.str2id.end()) {
      modify.str2id[cell] = id;
      modify.id2str[id] = cell;
      col[i] = id;
      ++id;
    } else {
      col[i] = iter->second;
    }
  }
  col_name = "s_nation";
  auto supp_nation = tb_supp[col_name].data<types::Char<15>>();
  col = (uint64_t*) memalign(64, tb_supp.nrTuples * sizeof(uint64_t));
  modify.col_addresses[col_name] = col;
  for (uint64_t i = 0; i < tb_supp.nrTuples; ++i) {
    memcpy(temp,supp_nation[i].value,supp_nation[i].length());
    temp[supp_nation[i].length()]='\0';
    cell = std::string(temp);
    auto iter = modify.str2id.find(cell);
    if (iter == modify.str2id.end()) {
      modify.str2id[cell] = id;
      modify.id2str[id] = cell;
      col[i] = id;
      ++id;
    } else {
      col[i] = iter->second;
    }
  }
  std::cout<<"supplier is ok-----------"<<"size = "<<modify.str2id.size()<<std::endl;

  ///////////
  return;
//  auto& tb_part = relations["part"];
//  col_name = "p_mfgr";
//  auto p_mfgr = tb_part[col_name].data<types::Char<6>>();
//  col = (uint64_t*) memalign(64, tb_part.nrTuples * sizeof(uint64_t));
//  modify.col_addresses[col_name] = col;
//  for(uint64_t i=0;i<tb_part.nrTuples;++i) {
//    cell = p_mfgr[i].value;
//    cell = cell.substr(5,10);
//    col[i] = stoi(cell);
//  }
//  col_name = "p_category";
//  auto p_category = tb_part[col_name].data<types::Char<7>>();
//  col = (uint64_t*) memalign(64, tb_part.nrTuples * sizeof(uint64_t));
//  modify.col_addresses[col_name] = col;
//  for(uint64_t i=0;i<tb_part.nrTuples;++i) {
//    cell = p_category[i].value;
//    cell = cell.substr(5,10);
//    col[i] = stoi(cell);
//  }
//
//  col_name = "p_brand1";
//  auto p_brand1 = tb_part[col_name].data<types::Char<9>>();
//  col = (uint64_t*) memalign(64, tb_part.nrTuples * sizeof(uint64_t));
//  modify.col_addresses[col_name] = col;
//  for(uint64_t i=0;i<tb_part.nrTuples;++i) {
//    cell = p_brand1[i].value;
//    cell = cell.substr(5,10);
//    col[i] = stoi(cell);
//  }


}

}  // namespace runtime

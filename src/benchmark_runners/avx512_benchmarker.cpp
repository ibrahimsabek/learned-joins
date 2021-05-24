#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "emmintrin.h"
#include "immintrin.h"
#include "smmintrin.h"
#include <sys/time.h> /* gettimeofday */

#include "config.h"            /* autoconf header */
#include "configs/base_configs.h"
#include "configs/eth_configs.h"

#include "utils/eth_data_structures.h"
#include "utils/data_generation.h"
#include "utils/io.h"
#include "utils/csv.h"

#include "utils/base_utils.h"
#include "utils/math.h"
#include "utils/barrier.h"
#include "utils/memory.h"
#include "utils/lock.h" 
#include "utils/learned_sort_for_sort_merge.h"

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

#include <chrono>
using namespace chrono;

#ifndef KeyType
#define KeyType RELATION_KEY_TYPE
#define PayloadType RELATION_PAYLOAD_TYPE
#endif

using namespace std;


void printDuration(std::chrono::steady_clock::time_point start, std::chrono::steady_clock::time_point end, const char *message) {
         auto diff = end - start;
         std::cout << message << ' ' << std::chrono::duration <double, std::milli> (diff).count() << " ms\n";
 }

 void test_vectorized_hashing(Relation<RELATION_KEY_TYPE, RELATION_PAYLOAD_TYPE>* rel_r, uint32_t MASK)
 {

    int64_t num_tuples_batches = rel_r->num_tuples / 8;
    uint32_t num_tuples_reminders = rel_r->num_tuples % 8;
    int64_t i = 0;    

    PerfEvent profiler;
 
    profiler.startCounters();
    
    const auto start = std::chrono::steady_clock::now();     
    for(int64_t j = 0; j < num_tuples_batches; j++)
    {
        uint64_t idx_hash[8];
        avx_murmur_hash_32((uint64_t*)(&(rel_r->tuples[j * 8])), idx_hash);
        
        for(uint32_t k = 0; k < 8; k++){
            uint32_t idx = HASH_BIT_MODULO(idx_hash[k], MASK, NUM_RADIX_BITS);
        }

    }

    for(uint32_t l = num_tuples_batches * 8; l < num_tuples_reminders; l++)
    {
        uint32_t idx_hash = murmur_hash_32(rel_r->tuples[l].key);
        uint32_t idx = HASH_BIT_MODULO(idx_hash, MASK, NUM_RADIX_BITS);
    }
    const auto end = std::chrono::steady_clock::now();

    profiler.stopCounters();


    for (i = 0; i < profiler.events.size(); i++)
    {
    if(profiler.names[i] == "cycles")
        std::cout << "cycles: " << profiler.events[i].readCounter() << "\n";
    else if (profiler.names[i] == "instructions")
        std::cout << "instructions: " << profiler.events[i].readCounter() << "\n";
    else if (profiler.names[i] == "L1-misses")
        std::cout << "L1-misses: " << profiler.events[i].readCounter() << "\n";
    else if (profiler.names[i] == "LLC-misses")
        std::cout << "LLC-misses: " << profiler.events[i].readCounter() << "\n";
    else if (profiler.names[i] == "branch-misses")
        std::cout << "branch-misses: " << profiler.events[i].readCounter() << "\n";
    else if (profiler.names[i] == "task-clock")
        std::cout << "task-clock: " << profiler.events[i].readCounter() << "\n";
    }

    std::cout << "IPC: " << profiler.getIPC() << "\n";
    std::cout << "CPUs: " << profiler.getCPUs() << "\n";
    std::cout << "GHz: " << profiler.getGHz() << "\n";

    printDuration(start, end, "vectorized_hashing: ");
 } 

void test_vectorized_hashing_v1(Relation<RELATION_KEY_TYPE, RELATION_PAYLOAD_TYPE>* rel_r, uint32_t MASK)
 {

    int64_t num_tuples_batches = rel_r->num_tuples / 8;
    int64_t i = 0;    
    uint64_t idx_hash[8];
    double total_time = 0;
    uint64_t total_sum = 0;

    __m512i k, s1, x1, m1, s2, x2, m2, s3, x3;
    const __m512i c1 = _mm512_set1_epi64(0x85ebca6b);
    const __m512i c2 = _mm512_set1_epi64(0xc2b2ae35);
    const __m512i mask_avx = _mm512_set1_epi64((uint64_t)MASK);

    //PerfEvent profiler;
 
    //profiler.startCounters();
    
    for(int64_t j = 0; j < num_tuples_batches; j++)
    {
        const auto start = std::chrono::steady_clock::now();     

        k = _mm512_load_epi64((uint64_t *)(&(rel_r->tuples[j * 8])));
        
        s1 = _mm512_srli_epi64(k, 16);
        x1 = _mm512_xor_epi64(k, s1);
        m1 = _mm512_mullo_epi64(x1, c1);

        s2 = _mm512_srli_epi64(m1, 13);
        x2 = _mm512_xor_epi64(m1, s2);
        m2 = _mm512_mullo_epi64(x2, c2);

        s3 = _mm512_srli_epi64(m2, 16);
        x3 = _mm512_xor_epi64(m2, s3);
        
        x3 = _mm512_and_epi64(x3, mask_avx);
        x3 = _mm512_srli_epi64(x3, NUM_RADIX_BITS);

        _mm512_store_epi64((uint64_t *) idx_hash, x3);

        const auto end = std::chrono::steady_clock::now(); 
        auto diff = end - start;
        total_time = total_time + std::chrono::duration <double, std::milli> (diff).count();

        for(int l = 0; l < 8; l++)
            total_sum = total_sum + idx_hash[l]; 
    }

    //profiler.stopCounters();


    /*for (i = 0; i < profiler.events.size(); i++)
    {
    if(profiler.names[i] == "cycles")
        std::cout << "cycles: " << profiler.events[i].readCounter() << "\n";
    else if (profiler.names[i] == "instructions")
        std::cout << "instructions: " << profiler.events[i].readCounter() << "\n";
    else if (profiler.names[i] == "L1-misses")
        std::cout << "L1-misses: " << profiler.events[i].readCounter() << "\n";
    else if (profiler.names[i] == "LLC-misses")
        std::cout << "LLC-misses: " << profiler.events[i].readCounter() << "\n";
    else if (profiler.names[i] == "branch-misses")
        std::cout << "branch-misses: " << profiler.events[i].readCounter() << "\n";
    else if (profiler.names[i] == "task-clock")
        std::cout << "task-clock: " << profiler.events[i].readCounter() << "\n";
    }

    std::cout << "IPC: " << profiler.getIPC() << "\n";
    std::cout << "CPUs: " << profiler.getCPUs() << "\n";
    std::cout << "GHz: " << profiler.getGHz() << "\n";
*/
    std::cout << "vectorized_hashing: " <<  total_time << " ms\n";
    std::cout << "total_sum: " << total_sum << endl;
 } 

 void test_vectorized_learned_models(Relation<RELATION_KEY_TYPE, RELATION_PAYLOAD_TYPE>* rel_r)
 {

  unsigned int TOT_NUM_MINOR_BCKTS = 10000;
  unsigned int MINOR_BCKTS_OFFSET = 100000;
  unsigned int NUM_MINOR_BCKT_PER_MAJOR_BCKT = 1000;

  auto root_slope = 0.5;
  auto root_intrcpt = 1.5;
  unsigned int num_models = 1000;
  vector<double> slopes, intercepts;
  for (unsigned int i = 0; i < num_models; ++i) {
    slopes.push_back(0.6);
    intercepts.push_back(2.2);
  }

    int64_t num_tuples_batches = rel_r->num_tuples / 8;
    uint32_t num_tuples_reminders = rel_r->num_tuples % 8;
    int64_t i = 0;

    __m512i  v_slopes_addr = _mm512_set1_epi64((uint64_t) (&slopes[0])), v_intercepts_addr = _mm512_set1_epi64((uint64_t) (&intercepts[0])),
             general_reg_1, general_reg_2;

    __m512d general_reg_1_double, general_reg_2_double, root_slope_avx = _mm512_set1_pd(root_slope), root_intrcpt_avx = _mm512_set1_pd(root_intrcpt),
          num_models_minus_one_avx = _mm512_set1_pd((double)num_models - 1.), v_zero512_double = _mm512_set1_pd(0.), 
          v_64bit_elem_size_double = _mm512_set1_pd(8.), minor_bckts_offset_avx = _mm512_set1_pd((double)MINOR_BCKTS_OFFSET),
          num_minor_bckt_per_major_bckt_minus_one_avx = _mm512_set1_pd((double)NUM_MINOR_BCKT_PER_MAJOR_BCKT - 1.),
          tot_num_minor_bckts_avx = _mm512_set1_pd((double)TOT_NUM_MINOR_BCKTS), 
          intercepts_avx, slopes_avx;
    __m512i k, pred_model_idx;

    //slopes_avx = _mm512_i64gather_pd(_mm512_add_epi64(pred_model_idx, v_slopes_addr), 0, 1);
    //intercepts_avx = _mm512_i64gather_pd(_mm512_add_epi64(pred_model_idx, v_intercepts_addr), 0, 1);
    slopes_avx = _mm512_set1_pd(1.);
    intercepts_avx = _mm512_set1_pd(8.);
    
    k = _mm512_set1_epi64(1);

    PerfEvent profiler;
 
    profiler.startCounters();

    const auto start = std::chrono::steady_clock::now();     
    for(int64_t j = 0; j < num_tuples_batches; j++)
    {

        k = _mm512_load_epi64((uint64_t*)(&(rel_r->tuples[j * 8])));

        //general_reg_2_double = _mm512_cvtepi64_pd(k);

        general_reg_1_double = _mm512_mul_pd(
                                    _mm512_floor_pd(
                                        _mm512_max_pd(
                                            _mm512_min_pd(
                                                _mm512_fmadd_pd(general_reg_2_double, root_slope_avx, root_intrcpt_avx), num_models_minus_one_avx), v_zero512_double)), v_64bit_elem_size_double);            

        pred_model_idx = _mm512_cvtpd_epi64(general_reg_1_double);

        general_reg_1_double = _mm512_fmadd_pd(general_reg_2_double, slopes_avx, intercepts_avx);

        general_reg_1_double = /*_mm512_floor_pd(*/
                                    _mm512_max_pd(
                                        _mm512_min_pd(
                                            _mm512_fmsub_pd(general_reg_1_double, tot_num_minor_bckts_avx, minor_bckts_offset_avx), num_minor_bckt_per_major_bckt_minus_one_avx), v_zero512_double)/*)*/;

        //pred_model_idx = _mm512_cvtpd_epi64(general_reg_1_double);
    }
    const auto end = std::chrono::steady_clock::now();

    profiler.stopCounters();


    for (i = 0; i < profiler.events.size(); i++)
    {
    if(profiler.names[i] == "cycles")
        std::cout << "cycles: " << profiler.events[i].readCounter() << "\n";
    else if (profiler.names[i] == "instructions")
        std::cout << "instructions: " << profiler.events[i].readCounter() << "\n";
    else if (profiler.names[i] == "L1-misses")
        std::cout << "L1-misses: " << profiler.events[i].readCounter() << "\n";
    else if (profiler.names[i] == "LLC-misses")
        std::cout << "LLC-misses: " << profiler.events[i].readCounter() << "\n";
    else if (profiler.names[i] == "branch-misses")
        std::cout << "branch-misses: " << profiler.events[i].readCounter() << "\n";
    else if (profiler.names[i] == "task-clock")
        std::cout << "task-clock: " << profiler.events[i].readCounter() << "\n";
    }

    std::cout << "IPC: " << profiler.getIPC() << "\n";
    std::cout << "CPUs: " << profiler.getCPUs() << "\n";
    std::cout << "GHz: " << profiler.getGHz() << "\n";

    printDuration(start, end, "vectorized_learned_models: ");
 } 

 void test_vectorized_learned_models_v1(Relation<RELATION_KEY_TYPE, RELATION_PAYLOAD_TYPE>* rel_r)
 {

  unsigned int TOT_NUM_MINOR_BCKTS = 10000;
  unsigned int MINOR_BCKTS_OFFSET = 100000;
  unsigned int NUM_MINOR_BCKT_PER_MAJOR_BCKT = 1000;

  auto root_slope = 0.5;
  auto root_intrcpt = 1.5;
  unsigned int num_models = 1000;
  vector<double> slopes, intercepts;
  for (unsigned int i = 0; i < num_models; ++i) {
    slopes.push_back(0.6);
    intercepts.push_back(2.2);
  }

    int64_t num_tuples_batches = rel_r->num_tuples / 8;
    int64_t i = 0;

    uint64_t idx_hash[8];
    double total_time = 0;
    uint64_t total_sum = 0;

    __m512i  v_slopes_addr = _mm512_set1_epi64((uint64_t) (&slopes[0])), v_intercepts_addr = _mm512_set1_epi64((uint64_t) (&intercepts[0])),
             general_reg_1, general_reg_2;

    __m512d general_reg_1_double, general_reg_2_double, root_slope_avx = _mm512_set1_pd(root_slope), root_intrcpt_avx = _mm512_set1_pd(root_intrcpt),
          num_models_minus_one_avx = _mm512_set1_pd((double)num_models - 1.), v_zero512_double = _mm512_set1_pd(0.), 
          v_64bit_elem_size_double = _mm512_set1_pd(8.), minor_bckts_offset_avx = _mm512_set1_pd((double)MINOR_BCKTS_OFFSET),
          num_minor_bckt_per_major_bckt_minus_one_avx = _mm512_set1_pd((double)NUM_MINOR_BCKT_PER_MAJOR_BCKT - 1.),
          tot_num_minor_bckts_avx = _mm512_set1_pd((double)TOT_NUM_MINOR_BCKTS), 
          intercepts_avx, slopes_avx;
    __m512i k, pred_model_idx;

    //slopes_avx = _mm512_i64gather_pd(_mm512_add_epi64(pred_model_idx, v_slopes_addr), 0, 1);
    //intercepts_avx = _mm512_i64gather_pd(_mm512_add_epi64(pred_model_idx, v_intercepts_addr), 0, 1);
    slopes_avx = _mm512_set1_pd(1.);
    intercepts_avx = _mm512_set1_pd(8.);
    
    //PerfEvent profiler;
 
    //profiler.startCounters();

    for(int64_t j = 0; j < num_tuples_batches; j++)
    {
        const auto start = std::chrono::steady_clock::now();     

        k = _mm512_load_epi64((uint64_t*)(&(rel_r->tuples[j * 8])));

        general_reg_2_double = _mm512_cvtepi64_pd(k);

        general_reg_1_double = _mm512_mul_pd(
                                    /*_mm512_floor_pd(*/
                                        _mm512_max_pd(
                                            _mm512_min_pd(
                                                _mm512_fmadd_pd(general_reg_2_double, root_slope_avx, root_intrcpt_avx), num_models_minus_one_avx), v_zero512_double)/*)*/, v_64bit_elem_size_double);            

        pred_model_idx = _mm512_cvtpd_epi64(general_reg_1_double);

        general_reg_1_double = _mm512_fmadd_pd(general_reg_2_double, slopes_avx, intercepts_avx);

        general_reg_1_double = /*_mm512_floor_pd(*/
                                    _mm512_max_pd(
                                        _mm512_min_pd(
                                            _mm512_fmsub_pd(general_reg_1_double, tot_num_minor_bckts_avx, minor_bckts_offset_avx), num_minor_bckt_per_major_bckt_minus_one_avx), v_zero512_double)/*)*/;

        pred_model_idx = _mm512_cvtpd_epi64(general_reg_1_double);

        _mm512_store_epi64((uint64_t *) idx_hash, pred_model_idx);

        const auto end = std::chrono::steady_clock::now();
        auto diff = end - start;
        total_time = total_time + std::chrono::duration <double, std::milli> (diff).count();

        for(int l = 0; l < 8; l++)
            total_sum = total_sum + idx_hash[l]; 
    }

    //profiler.stopCounters();


    /*for (i = 0; i < profiler.events.size(); i++)
    {
    if(profiler.names[i] == "cycles")
        std::cout << "cycles: " << profiler.events[i].readCounter() << "\n";
    else if (profiler.names[i] == "instructions")
        std::cout << "instructions: " << profiler.events[i].readCounter() << "\n";
    else if (profiler.names[i] == "L1-misses")
        std::cout << "L1-misses: " << profiler.events[i].readCounter() << "\n";
    else if (profiler.names[i] == "LLC-misses")
        std::cout << "LLC-misses: " << profiler.events[i].readCounter() << "\n";
    else if (profiler.names[i] == "branch-misses")
        std::cout << "branch-misses: " << profiler.events[i].readCounter() << "\n";
    else if (profiler.names[i] == "task-clock")
        std::cout << "task-clock: " << profiler.events[i].readCounter() << "\n";
    }

    std::cout << "IPC: " << profiler.getIPC() << "\n";
    std::cout << "CPUs: " << profiler.getCPUs() << "\n";
    std::cout << "GHz: " << profiler.getGHz() << "\n";
    */

    std::cout << "vectorized_learned_models: " <<  total_time << " ms\n";
    std::cout << "total_sum: " << total_sum << endl;
 } 



int main(int argc, char **argv) 
{
    Relation<KeyType, PayloadType> rel_r;
    Relation<KeyType, PayloadType> rel_s;
    
    int64_t result = 0;
    uint64_t curr_num_tuples_r = RELATION_R_NUM_TUPLES;
    uint64_t curr_num_tuples_s = RELATION_S_NUM_TUPLES; 

#ifdef LOAD_RELATIONS_FOR_EVALUATION
    // loading pre-built datasets
    string curr_rel_r_folder_path = RELATION_R_FOLDER_PATH;
    string curr_rel_s_folder_path = RELATION_S_FOLDER_PATH;

    string curr_rel_r_file_name = RELATION_R_FILE_NAME;
    string curr_rel_s_file_name = RELATION_S_FILE_NAME;

    string curr_rel_r_file_extension = RELATION_R_FILE_EXTENSION;
    string curr_rel_s_file_extension = RELATION_S_FILE_EXTENSION;

    load_relation_threaded<KeyType, PayloadType>(&rel_r, RELATION_R_FILE_NUM_PARTITIONS, curr_rel_r_folder_path.c_str(), curr_rel_r_file_name.c_str(), curr_rel_r_file_extension.c_str(), curr_num_tuples_r);
    load_relation_threaded<KeyType, PayloadType>(&rel_s, RELATION_S_FILE_NUM_PARTITIONS, curr_rel_s_folder_path.c_str(), curr_rel_s_file_name.c_str(), curr_rel_s_file_extension.c_str(), curr_num_tuples_s);
#else

    string curr_rel_r_path = RELATION_R_PATH;
    string curr_rel_s_path = RELATION_S_PATH;

    string curr_rel_r_folder_path = RELATION_R_FOLDER_PATH;
    string curr_rel_s_folder_path = RELATION_S_FOLDER_PATH;

    string curr_rel_r_file_name = RELATION_R_FILE_NAME;
    string curr_rel_s_file_name = RELATION_S_FILE_NAME;

    string curr_rel_r_file_extension = RELATION_R_FILE_EXTENSION;
    string curr_rel_s_file_extension = RELATION_S_FILE_EXTENSION;

    // creating new datasets on-the-flay 
    result = create_eth_workload_relation_pk<KeyType, PayloadType>(&rel_r, curr_num_tuples_r, 0);
    //ASSERT_EQ(result, 0);
    #ifdef PERSIST_RELATIONS_FOR_EVALUATION
    write_relation_threaded<KeyType, PayloadType>(&rel_r, RELATION_R_FILE_NUM_PARTITIONS, curr_rel_r_folder_path.c_str(), curr_rel_r_file_name.c_str(), curr_rel_r_file_extension.c_str());
    write_relation<KeyType, PayloadType>(&rel_r, curr_rel_r_path.c_str());    
    #endif
    
    result = create_eth_workload_relation_pk<KeyType, PayloadType>(&rel_s, curr_num_tuples_s, 0);
    //ASSERT_EQ(result, 0);
    #ifdef PERSIST_RELATIONS_FOR_EVALUATION
    write_relation_threaded<KeyType, PayloadType>(&rel_s, RELATION_S_FILE_NUM_PARTITIONS, curr_rel_s_folder_path.c_str(), curr_rel_s_file_name.c_str(), curr_rel_s_file_extension.c_str());
    write_relation<KeyType, PayloadType>(&rel_s, curr_rel_s_path.c_str());    
    #endif
#endif

    uint32_t N = rel_r.num_tuples;
    NEXT_POW_2(N);
    const uint32_t MASK = (N-1) << (NUM_RADIX_BITS);

    //test_vectorized_hashing(&rel_r, MASK);
    test_vectorized_hashing_v1(&rel_r, MASK);

    //test_vectorized_learned_models(&rel_r);
    test_vectorized_learned_models_v1(&rel_r);


    return 0;
}

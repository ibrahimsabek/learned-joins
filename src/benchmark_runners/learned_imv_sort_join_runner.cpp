#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "emmintrin.h"
#include "immintrin.h"
#include "smmintrin.h"
#include <sys/time.h> /* gettimeofday */
#include <map> 
#include <iterator> 

#include "config.h"            /* autoconf header */
#include "configs/base_configs.h"
#include "configs/eth_configs.h"

#include "utils/eth_data_structures.h"
#include "utils/data_generation.h"
#include "utils/io.h"
#include "utils/eth_generic_task_queue.h"
#include "utils/cpu_mapping.h"
#include "utils/base_utils.h"
#include "utils/math.h"
#include "utils/barrier.h"
#include "utils/memory.h"
#include "utils/lock.h" 
#include "utils/numa_shuffle.h"
#include "utils/learned_sort_for_sort_merge.h"


#include "techniques/sortmerge_multiway_join_learned_steps.h"

#ifndef KeyType
#define KeyType RELATION_KEY_TYPE
#define PayloadType RELATION_PAYLOAD_TYPE
#define NUM_THREADS NUM_THREADS_FOR_EVALUATION
#endif


#define RUN_NUMS 1 //10 

#define PREFETCH_SLOPES_AND_INTERCEPTS_MAJOR_BCKTS_UNIQUE_KEYS

using namespace std;
using namespace learned_sort_for_sort_merge;


LearnedSortMergeMultiwayJoinSteps<KeyType, PayloadType, 
                                LearnedSortMergeMultiwayJoinThread<KeyType, PayloadType>> join_steps;
joinconfig_t joincfg;
int CACHELINEPADDING;
int RELATION_PADDING;

learned_sort_for_sort_merge::RMI<KeyType, PayloadType>* rmi_ptr;
learned_sort_for_sort_merge::RMI<KeyType, PayloadType>* r_rmi_ptr;
learned_sort_for_sort_merge::RMI<KeyType, PayloadType>* s_rmi_ptr;

int64_t* r_initial_partition_sizes_for_threads;
int64_t* s_initial_partition_sizes_for_threads;
int64_t* r_partition_sizes_for_threads;
int64_t* s_partition_sizes_for_threads;
int64_t* r_partition_offsets;
int64_t* s_partition_offsets;
Tuple<KeyType, PayloadType> * tmpRelpartR;
Tuple<KeyType, PayloadType> * tmpRelpartS;


int64_t* r_initial_repeated_keys_sizes_for_threads;
int64_t* s_initial_repeated_keys_sizes_for_threads;
int64_t* r_repeated_keys_sizes_for_threads;
int64_t* s_repeated_keys_sizes_for_threads;
int64_t* r_total_repeated_keys_sizes_for_threads;
int64_t* s_total_repeated_keys_sizes_for_threads;    
int64_t* r_repeated_keys_offsets;
int64_t* s_repeated_keys_offsets;
Tuple<KeyType, PayloadType> * tmpRepeatedKeysPredictedRanksR;
Tuple<KeyType, PayloadType> * tmpRepeatedKeysPredictedRanksS;
int64_t* tmpRepeatedKeysCountsR;
int64_t* tmpRepeatedKeysCountsS;

unsigned int * r_NUM_MINOR_BCKT_PER_MAJOR_BCKT;
unsigned int * s_NUM_MINOR_BCKT_PER_MAJOR_BCKT;
unsigned int * r_MINOR_BCKTS_OFFSET;
unsigned int * s_MINOR_BCKTS_OFFSET;
unsigned int r_TOT_NUM_MINOR_BCKTS;
unsigned int s_TOT_NUM_MINOR_BCKTS;

Tuple<KeyType, PayloadType> * tmpRelsortR;
Tuple<KeyType, PayloadType> * tmpRelsortS;
Tuple<KeyType, PayloadType> * r_tmp_spill_bucket;
Tuple<KeyType, PayloadType> * s_tmp_spill_bucket;
Tuple<KeyType, PayloadType> * r_sorted_spill_bucket;
Tuple<KeyType, PayloadType> * s_sorted_spill_bucket;
Tuple<KeyType, PayloadType> * r_tmp_minor_bckts;
Tuple<KeyType, PayloadType> * s_tmp_minor_bckts;
int64_t* r_minor_bckts_offsets;
int64_t* s_minor_bckts_offsets;
int64_t* r_tmp_minor_bckt_sizes;
int64_t* s_tmp_minor_bckt_sizes;
int64_t* r_minor_bckt_sizes_offsets;
int64_t* s_minor_bckt_sizes_offsets;

RelationPair<KeyType, PayloadType> ** threadrelchunks;
//uint32_t ** histR;
Tuple<KeyType, PayloadType> ** ptrs_to_sharedmergebufs;
#ifdef SKEW_HANDLING
taskqueue_t ** ptrs_to_taskqueues;
pthread_mutex_t* ptrs_to_taskqueues_locks;
int* ptrs_to_is_numa_taskqueues_created;
#endif



void initialize_learned_imv_sort_join_thread_args(Relation<KeyType, PayloadType> * rel_r, 
                        Relation<KeyType, PayloadType> * rel_s, 
                        learned_sort_for_sort_merge::RMI<KeyType, PayloadType> * rmi,
                        learned_sort_for_sort_merge::RMI<KeyType, PayloadType> * r_rmi,
                        learned_sort_for_sort_merge::RMI<KeyType, PayloadType> * s_rmi,
                        learned_sort_for_sort_merge::RMI<KeyType, PayloadType>::Params p,
                        learned_sort_for_sort_merge::RMI<KeyType, PayloadType>::Params r_p,
                        learned_sort_for_sort_merge::RMI<KeyType, PayloadType>::Params s_p,                                 
                        unsigned int SAMPLE_SZ_R, unsigned int SAMPLE_SZ_S,
                        Tuple<KeyType, PayloadType> * tmp_training_sample_in,
                        Tuple<KeyType, PayloadType> * sorted_training_sample_in,
                        Tuple<KeyType, PayloadType> * r_tmp_training_sample_in,
                        Tuple<KeyType, PayloadType> * r_sorted_training_sample_in,
                        Tuple<KeyType, PayloadType> * s_tmp_training_sample_in,
                        Tuple<KeyType, PayloadType> * s_sorted_training_sample_in,
                        vector<vector<vector<learned_sort_for_sort_merge::training_point<KeyType, PayloadType>>>> * training_data,
                        vector<vector<vector<learned_sort_for_sort_merge::training_point<KeyType, PayloadType>>>> * r_training_data,
                        vector<vector<vector<learned_sort_for_sort_merge::training_point<KeyType, PayloadType>>>> * s_training_data,
                        uint32_t * sample_count, uint32_t * sample_count_R, uint32_t * sample_count_S,
                        vector<double>* slopes, vector<double>* intercepts,
                        vector<double>* r_slopes, vector<double>* r_intercepts,
                        vector<double>* s_slopes, vector<double>* s_intercepts,
                        pthread_barrier_t* barrier_ptr,
                        Result * joinresult,
                        LearnedSortMergeMultiwayJoinThread<KeyType, PayloadType> * args)
{
    int32_t numperthr[2];
    unsigned int SAMPLE_SZ_Rthr, SAMPLE_SZ_Sthr;

    numperthr[0] = rel_r->num_tuples / NUM_THREADS;
    numperthr[1] = rel_s->num_tuples / NUM_THREADS;
    SAMPLE_SZ_Rthr = SAMPLE_SZ_R / NUM_THREADS;
    SAMPLE_SZ_Sthr = SAMPLE_SZ_S / NUM_THREADS;

    int i;
    for(i = 0; i < NUM_THREADS; i++)
    {
        (*(args + i)).relR = rel_r->tuples + i * (numperthr[0]);
        (*(args + i)).relS = rel_s->tuples + i * (numperthr[1]);

        (*(args + i)).numR_to_be_partitioned = (i == (NUM_THREADS-1)) ?
            (rel_r->num_tuples - i * numperthr[0]) : numperthr[0];
        (*(args + i)).numS_to_be_partitioned = (i == (NUM_THREADS-1)) ?
            (rel_s->num_tuples - i * numperthr[1]) : numperthr[1];

        // temporary relations
        (*(args + i)).tmp_partR = tmpRelpartR + r_partition_offsets[i];
        (*(args + i)).tmp_partS = tmpRelpartS + s_partition_offsets[i];
        (*(args + i)).tmp_repeatedKeysPredictedRanksR = tmpRepeatedKeysPredictedRanksR + r_repeated_keys_offsets[i];
        (*(args + i)).tmp_repeatedKeysPredictedRanksS = tmpRepeatedKeysPredictedRanksS + s_repeated_keys_offsets[i];
        (*(args + i)).tmp_repeatedKeysPredictedRanksCountsR = tmpRepeatedKeysCountsR + r_repeated_keys_offsets[i];
        (*(args + i)).tmp_repeatedKeysPredictedRanksCountsS = tmpRepeatedKeysCountsS + s_repeated_keys_offsets[i];
        (*(args + i)).tmp_sortR = tmpRelsortR + r_partition_offsets[i];
        (*(args + i)).tmp_sortS = tmpRelsortS + s_partition_offsets[i];
        (*(args + i)).tmp_spill_bucket_r = r_tmp_spill_bucket + r_partition_offsets[i];
        (*(args + i)).sorted_spill_bucket_r = r_sorted_spill_bucket + r_partition_offsets[i];
        (*(args + i)).tmp_spill_bucket_s = s_tmp_spill_bucket + s_partition_offsets[i];
        (*(args + i)).sorted_spill_bucket_s = s_sorted_spill_bucket + s_partition_offsets[i];
        (*(args + i)).tmp_minor_bckts_r = r_tmp_minor_bckts + r_minor_bckts_offsets[i];
        (*(args + i)).tmp_minor_bckts_s = s_tmp_minor_bckts + s_minor_bckts_offsets[i];
        (*(args + i)).tmp_minor_bckt_sizes_r = r_tmp_minor_bckt_sizes + r_minor_bckt_sizes_offsets[i];            
        (*(args + i)).tmp_minor_bckt_sizes_s = s_tmp_minor_bckt_sizes + s_minor_bckt_sizes_offsets[i];

        (*(args + i)).numR = r_partition_sizes_for_threads[i] + r_total_repeated_keys_sizes_for_threads[i];
        (*(args + i)).numS = s_partition_sizes_for_threads[i] + s_total_repeated_keys_sizes_for_threads[i];
        (*(args + i)).tmp_repeatedKeysCountsR = r_repeated_keys_sizes_for_threads[i];
        (*(args + i)).tmp_repeatedKeysCountsS = s_repeated_keys_sizes_for_threads[i];
        (*(args + i)).tmp_total_repeatedKeysCountsR = r_total_repeated_keys_sizes_for_threads[i];
        (*(args + i)).tmp_total_repeatedKeysCountsS = s_total_repeated_keys_sizes_for_threads[i];

        (*(args + i)).NUM_MINOR_BCKT_PER_MAJOR_BCKT_r = r_NUM_MINOR_BCKT_PER_MAJOR_BCKT[i]; 
        (*(args + i)).NUM_MINOR_BCKT_PER_MAJOR_BCKT_s = s_NUM_MINOR_BCKT_PER_MAJOR_BCKT[i];
        (*(args + i)).MINOR_BCKTS_OFFSET_r = r_MINOR_BCKTS_OFFSET[i];
        (*(args + i)).MINOR_BCKTS_OFFSET_s = s_MINOR_BCKTS_OFFSET[i];
        (*(args + i)).TOT_NUM_MINOR_BCKTS_r = r_TOT_NUM_MINOR_BCKTS;
        (*(args + i)).TOT_NUM_MINOR_BCKTS_s = s_TOT_NUM_MINOR_BCKTS;
        (*(args + i)).INPUT_SZ_r = rel_r->num_tuples;
        (*(args + i)).INPUT_SZ_s = rel_s->num_tuples; 

        /**** start stuff for learning RMI models ****/
        (*(args + i)).rmi = rmi;
        (*(args + i)).rmi_r = r_rmi_ptr;
        (*(args + i)).rmi_s = s_rmi_ptr;
        (*(args + i)).p = p;
        (*(args + i)).r_p = r_p;
        (*(args + i)).s_p = s_p;
        (*(args + i)).original_relR = rel_r;
        (*(args + i)).original_relS = rel_s;
        (*(args + i)).tmp_training_sample_in = tmp_training_sample_in;
        (*(args + i)).sorted_training_sample_in = sorted_training_sample_in;
        (*(args + i)).r_tmp_training_sample_in = r_tmp_training_sample_in;
        (*(args + i)).r_sorted_training_sample_in = r_sorted_training_sample_in;
        (*(args + i)).s_tmp_training_sample_in = s_tmp_training_sample_in;
        (*(args + i)).s_sorted_training_sample_in = s_sorted_training_sample_in;
        (*(args + i)).training_data = training_data;
        (*(args + i)).r_training_data = r_training_data;
        (*(args + i)).s_training_data = s_training_data;
        (*(args + i)).tmp_training_sample_R_offset = SAMPLE_SZ_Rthr * i;
        (*(args + i)).tmp_training_sample_S_offset = SAMPLE_SZ_Sthr * i;
        (*(args + i)).tmp_training_sample_offset = (SAMPLE_SZ_Rthr + SAMPLE_SZ_Sthr) * i;
        (*(args + i)).sample_count = sample_count;
        (*(args + i)).sample_count_R = sample_count_R;
        (*(args + i)).sample_count_S = sample_count_S;
        (*(args + i)).slopes = slopes;
        (*(args + i)).intercepts = intercepts;
        (*(args + i)).r_slopes = r_slopes;
        (*(args + i)).r_intercepts = r_intercepts;
        (*(args + i)).s_slopes = s_slopes;
        (*(args + i)).s_intercepts = s_intercepts;
        /**** end stuff for learning RMI models ****/

        (*(args + i)).my_tid        = i;/* this is the logical CPU-ID */
        (*(args + i)).nthreads      = NUM_THREADS;
        (*(args + i)).joincfg       = &joincfg;
        (*(args + i)).barrier       = barrier_ptr;
        (*(args + i)).threadrelchunks = threadrelchunks;
        //(*(args + i)).sharedmergebuffer = ptrs_to_sharedmergebufs;

        // information specific to mpsm-join
        //(*(args + i)).histR         = histR;
        //(*(args + i)).tmpRglobal    = tmpRelpartR;
        //(*(args + i)).totalR        = relR->num_tuples;
        (*(args + i)).threadresult  = &(joinresult->resultlist[i]);
        
    #ifdef SKEW_HANDLING
        // skew handling task queue ptrs.
        (*(args + i)).numa_taskqueues     = ptrs_to_taskqueues;
        (*(args + i)).numa_taskqueues_locks = ptrs_to_taskqueues_locks;
        (*(args + i)).is_numa_taskqueues_created = ptrs_to_is_numa_taskqueues_created;
    #endif

    }
}

#ifdef BUILD_RMI_FROM_TWO_DATASETS
void sample_and_train_models_threaded(LearnedSortMergeMultiwayJoinThread<KeyType, PayloadType> * args)
{
    int rv;
    int tid = args->my_tid;

    //----------------------------------------------------------//
    //                           SAMPLE                         //
    //----------------------------------------------------------//

    // Determine sample size
    unsigned int INPUT_SZ_R = args->original_relR->num_tuples;
    unsigned int INPUT_SZ_S = args->original_relS->num_tuples;
    const unsigned int SAMPLE_SZ_R = std::min<unsigned int>(
        INPUT_SZ_R, std::max<unsigned int>(args->p.sampling_rate * INPUT_SZ_R,
                                        learned_sort_for_sort_merge::RMI<KeyType, PayloadType>::Params::MIN_SORTING_SIZE));
    const unsigned int SAMPLE_SZ_S = std::min<unsigned int>(
        INPUT_SZ_S, std::max<unsigned int>(args->p.sampling_rate * INPUT_SZ_S,
                                        learned_sort_for_sort_merge::RMI<KeyType, PayloadType>::Params::MIN_SORTING_SIZE));
    //const unsigned int SAMPLE_SZ = SAMPLE_SZ_R + SAMPLE_SZ_S;        
    
    // Start sampling
    Tuple<KeyType, PayloadType> *   relR_start_sampling_ptr = args->relR;
    Tuple<KeyType, PayloadType> *   relR_end_sampling_ptr = args->relR + (args->numR_to_be_partitioned - 1);
    Tuple<KeyType, PayloadType> *   relS_start_sampling_ptr = args->relS;
    Tuple<KeyType, PayloadType> *   relS_end_sampling_ptr = args->relS + (args->numS_to_be_partitioned - 1);

    uint32_t * sample_count = args->sample_count; 
    Tuple<KeyType, PayloadType> * tmp_training_sample = args->rmi->tmp_training_sample + args->tmp_training_sample_offset;
    
    uint32_t * sample_count_R = args->sample_count_R; 
    Tuple<KeyType, PayloadType> * tmp_training_sample_R = args->rmi->tmp_training_sample_R + args->tmp_training_sample_R_offset;
    unsigned int offset_R = static_cast<unsigned int>(1. * INPUT_SZ_R / SAMPLE_SZ_R);
    for (auto i = relR_start_sampling_ptr; i <= relR_end_sampling_ptr; i += offset_R) {
      // NOTE:  We don't directly assign SAMPLE_SZ to rmi.training_sample_sz
      //        to avoid issues with divisibility
      tmp_training_sample[sample_count[tid]] = *i;
      ++sample_count[tid];

      //tmp_training_sample_R[sample_count_R[tid]] = *i;
      //++sample_count_R[tid];
    }
    
    uint32_t * sample_count_S = args->sample_count_S;
    Tuple<KeyType, PayloadType> * tmp_training_sample_S = args->rmi->tmp_training_sample_S + args->tmp_training_sample_S_offset;
    unsigned int offset_S = static_cast<unsigned int>(1. * INPUT_SZ_S / SAMPLE_SZ_S);
    for (auto i = relS_start_sampling_ptr; i <= relS_end_sampling_ptr; i += offset_S) {
      // NOTE:  We don't directly assign SAMPLE_SZ to rmi.training_sample_sz
      //        to avoid issues with divisibility
      tmp_training_sample[sample_count[tid]] = *i;
      ++sample_count[tid];

      //tmp_training_sample_S[sample_count_S[tid]] = *i;
      //++sample_count_S[tid];
    }

    BARRIER_ARRIVE(args->barrier, rv);

    #ifdef USE_AVXSORT_AS_STD_SORT          
    uint32_t total_sample_count = 0; 
    if(tid == 0)
    {
      for(int i = 0; i < NUM_THREADS; i++)
        total_sample_count += sample_count[i];

      Tuple<KeyType, PayloadType> * sorted_training_sample = args->rmi->sorted_training_sample;      
      int64_t * inputptr =  (int64_t *)(args->rmi->tmp_training_sample);
      int64_t * outputptr = (int64_t *)(sorted_training_sample);
      avxsort_int64(&inputptr, &outputptr, total_sample_count);
      Tuple<KeyType, PayloadType>* tmp_outputptr = (Tuple<KeyType, PayloadType>*) outputptr;
      for(unsigned int k = 0; k < total_sample_count; k++){
        sorted_training_sample[k] = tmp_outputptr[k];
      } 
      args->rmi->training_sample = &(sorted_training_sample);
      args->rmi->training_sample_size = total_sample_count;
    }
    #else
    uint32_t total_sample_count = 0;
    if(tid == 0)
    {
      for(int i = 0; i < NUM_THREADS; i++)
        total_sample_count += sample_count[i];

      std::sort((int64_t *)(args->rmi->tmp_training_sample), (int64_t *)(args->rmi->tmp_training_sample) + total_sample_count - 1);
      args->rmi->training_sample = &(args->rmi->tmp_training_sample);
      args->rmi->training_sample_size = total_sample_count;
    }
    #endif

    //----------------------------------------------------------//
    //                     TRAIN THE MODELS                     //
    //----------------------------------------------------------//

    if(tid == 0)
    {
       // Stop early if the array is identical
      if (((*(args->rmi->training_sample))[0]).key == ((*(args->rmi->training_sample))[total_sample_count - 1]).key) 
      {
        return;
      }
         
      // Populate the training data for the root model
      vector<vector<vector<learned_sort_for_sort_merge::training_point<KeyType, PayloadType>>>> * training_data = args->training_data;
      for (unsigned int i = 0; i < total_sample_count; ++i) {
        (*training_data)[0][0].push_back({((*(args->rmi->training_sample))[i]), 1. * i / total_sample_count});
      }

      // Train the root model using linear interpolation
      auto *current_training_data = &(*training_data)[0][0];
      typename learned_sort_for_sort_merge::RMI<KeyType, PayloadType>::linear_model *current_model = &args->rmi->models[0][0];

      // Find the min and max values in the training set
      learned_sort_for_sort_merge::training_point<KeyType, PayloadType> min = current_training_data->front();
      learned_sort_for_sort_merge::training_point<KeyType, PayloadType> max = current_training_data->back();

      // Calculate the slope and intercept terms
      current_model->slope =
          1. / (max.x.key - min.x.key);  // Assuming min.y = 0 and max.y = 1
      current_model->intercept = -current_model->slope * min.x.key;

      // Extrapolate for the number of models in the next layer
      current_model->slope *= args->p.arch[1] - 1;
      current_model->intercept *= args->p.arch[1] - 1;
#ifndef RUN_LEARNED_TECHNIQUES_WITH_FIRST_LEVEL_ONLY
      // Populate the training data for the next layer
      for (const auto &d : *current_training_data) {
        // Predict the model index in next layer
        unsigned int rank = current_model->slope * d.x.key + current_model->intercept;

        // Normalize the rank between 0 and the number of models in the next layer
        rank =
            std::max(static_cast<unsigned int>(0), std::min(args->p.arch[1] - 1, rank));
        
        //if(d.x.key >299 && d.x.key < 400)
            //printf("training: key %ld rank %ld \n", d.x.key, rank);
        // Place the data in the predicted training bucket
        (*training_data)[1][rank].push_back(d);
      }

      // Train the leaf models
      for (unsigned int model_idx = 0; model_idx < args->p.arch[1]; ++model_idx) {
        // Update iterator variables
        current_training_data = &(*training_data)[1][model_idx];
        current_model = &args->rmi->models[1][model_idx];

        // Interpolate the min points in the training buckets
        if (model_idx ==
            0) {  // The current model is the first model in the current layer

          if (current_training_data->size() <
              2) {  // Case 1: The first model in this layer is empty
            current_model->slope = 0;
            current_model->intercept = 0;

            // Insert a fictive training point to avoid propagating more than one
            // empty initial models.
            learned_sort_for_sort_merge::training_point<KeyType, PayloadType> tp;
            tp.x.key = 0;
            tp.x.payload = 0;
            tp.y = 0;
            current_training_data->push_back(tp);
          } else {  // Case 2: The first model in this layer is not empty

            min = current_training_data->front();
            max = current_training_data->back();

            current_model->slope =
                (max.y) / (max.x.key - min.x.key);  // Hallucinating as if min.y = 0
            current_model->intercept = min.y - current_model->slope * min.x.key;
          }
        } else if (model_idx == args->p.arch[1] - 1) {
          if (current_training_data
                  ->empty()) {  // Case 3: The final model in this layer is empty

            current_model->slope = 0;
            current_model->intercept = 1;
          } else {  // Case 4: The last model in this layer is not empty

            min = (*training_data)[1][model_idx - 1].back();
            max = current_training_data->back();

            current_model->slope =
                (min.y - 1) / (min.x.key - max.x.key);  // Hallucinating as if max.y = 1
            current_model->intercept = min.y - current_model->slope * min.x.key;
          }
        } else {  // The current model is not the first model in the current layer

          if (current_training_data->empty()) {  // Case 5: The intermediate model in
            // this layer is empty
            current_model->slope = 0;
            current_model->intercept =
                (*training_data)[1][model_idx - 1].back().y;  // If the previous model
                                                          // was empty too, it will
                                                          // use the fictive
                                                          // training points

            // Insert a fictive training point to avoid propagating more than one
            // empty initial models.
            // NOTE: This will _NOT_ throw to DIV/0 due to identical x's and y's
            // because it is working backwards.
            learned_sort_for_sort_merge::training_point<KeyType, PayloadType> tp;
            tp.x = (*training_data)[1][model_idx - 1].back().x;
            tp.y = (*training_data)[1][model_idx - 1].back().y;
            current_training_data->push_back(tp);
          } else {  // Case 6: The intermediate leaf model is not empty

            min = (*training_data)[1][model_idx - 1].back();
            max = current_training_data->back();

            current_model->slope = (min.y - max.y) / (min.x.key - max.x.key);
            current_model->intercept = min.y - current_model->slope * min.x.key;
            //if(model_idx == 5)
            //printf("min.y %lf max.y %lf min.x.key %ld max.x.key %ld current_model->slope %lf current_model->intercept %lf\n", min.y, max.y, min.x.key, max.x.key, current_model->slope, current_model->intercept);
          }
        }
      }
#endif
      
      args->rmi->trained = true;         
    }

}
#else
void sample_and_train_models_threaded(LearnedSortMergeMultiwayJoinThread<KeyType, PayloadType> * args)
{
    int rv;
    int tid = args->my_tid;

    //----------------------------------------------------------//
    //                           SAMPLE                         //
    //----------------------------------------------------------//

    // Determine sample size
    unsigned int INPUT_SZ_R = args->original_relR->num_tuples;
    unsigned int INPUT_SZ_S = args->original_relS->num_tuples;
    const unsigned int SAMPLE_SZ_R = std::min<unsigned int>(
        INPUT_SZ_R, std::max<unsigned int>(args->r_p.sampling_rate * INPUT_SZ_R,
                                        learned_sort_for_sort_merge::RMI<KeyType, PayloadType>::Params::MIN_SORTING_SIZE));
    const unsigned int SAMPLE_SZ_S = std::min<unsigned int>(
        INPUT_SZ_S, std::max<unsigned int>(args->s_p.sampling_rate * INPUT_SZ_S,
                                        learned_sort_for_sort_merge::RMI<KeyType, PayloadType>::Params::MIN_SORTING_SIZE));
    //const unsigned int SAMPLE_SZ = SAMPLE_SZ_R + SAMPLE_SZ_S;        
    
    // Start sampling
    Tuple<KeyType, PayloadType> *   relR_start_sampling_ptr = args->relR;
    Tuple<KeyType, PayloadType> *   relR_end_sampling_ptr = args->relR + (args->numR_to_be_partitioned - 1);
    Tuple<KeyType, PayloadType> *   relS_start_sampling_ptr = args->relS;
    Tuple<KeyType, PayloadType> *   relS_end_sampling_ptr = args->relS + (args->numS_to_be_partitioned - 1);

    //uint32_t * sample_count = args->sample_count; 
    //Tuple<KeyType, PayloadType> * tmp_training_sample = args->rmi->tmp_training_sample + args->tmp_training_sample_offset;
    
    uint32_t * sample_count_R = args->sample_count_R; 
    Tuple<KeyType, PayloadType> * tmp_training_sample_R = args->r_rmi->tmp_training_sample + args->tmp_training_sample_R_offset;
    unsigned int offset_R = static_cast<unsigned int>(1. * INPUT_SZ_R / SAMPLE_SZ_R);
    for (auto i = relR_start_sampling_ptr; i <= relR_end_sampling_ptr; i += offset_R) {
      // NOTE:  We don't directly assign SAMPLE_SZ to rmi.training_sample_sz
      //        to avoid issues with divisibility
      //tmp_training_sample[sample_count[tid]] = *i;
      //++sample_count[tid];

      tmp_training_sample_R[sample_count_R[tid]] = *i;
      ++sample_count_R[tid];
    }
    
    uint32_t * sample_count_S = args->sample_count_S;
    Tuple<KeyType, PayloadType> * tmp_training_sample_S = args->s_rmi->tmp_training_sample + args->tmp_training_sample_S_offset;
    unsigned int offset_S = static_cast<unsigned int>(1. * INPUT_SZ_S / SAMPLE_SZ_S);
    for (auto i = relS_start_sampling_ptr; i <= relS_end_sampling_ptr; i += offset_S) {
      // NOTE:  We don't directly assign SAMPLE_SZ to rmi.training_sample_sz
      //        to avoid issues with divisibility
      //tmp_training_sample[sample_count[tid]] = *i;
      //++sample_count[tid];

      tmp_training_sample_S[sample_count_S[tid]] = *i;
      ++sample_count_S[tid];
    }

    BARRIER_ARRIVE(args->barrier, rv);

    #ifdef USE_AVXSORT_AS_STD_SORT     
    uint32_t total_sample_count_R = 0;  
    if(tid == 0)
    {
      for(int i = 0; i < NUM_THREADS; i++)
        total_sample_count_R += sample_count_R[i];

      Tuple<KeyType, PayloadType> * sorted_training_sample_R = args->r_rmi->sorted_training_sample;
      int64_t * inputptr_R =  (int64_t *)(args->r_rmi->tmp_training_sample);
      int64_t * outputptr_R = (int64_t *)(sorted_training_sample_R);
      avxsort_int64(&inputptr_R, &outputptr_R, total_sample_count_R);
      Tuple<KeyType, PayloadType>* tmp_outputptr_R = (Tuple<KeyType, PayloadType>*) outputptr_R;
      for(unsigned int k = 0; k < total_sample_count_R; k++){
        sorted_training_sample_R[k] = tmp_outputptr_R[k];
      } 
      args->r_rmi->training_sample = &(sorted_training_sample_R);
      args->r_rmi->training_sample_size = total_sample_count_R;
    }

    uint32_t total_sample_count_S = 0; 
    if(tid == 1)
    {
      for(int i = 0; i < NUM_THREADS; i++)
        total_sample_count_S += sample_count_S[i];

      Tuple<KeyType, PayloadType> * sorted_training_sample_S = args->s_rmi->sorted_training_sample;
      int64_t * inputptr_S =  (int64_t *)(args->s_rmi->tmp_training_sample);
      int64_t * outputptr_S = (int64_t *)(sorted_training_sample_S);
      avxsort_int64(&inputptr_S, &outputptr_S, total_sample_count_S);
      Tuple<KeyType, PayloadType>* tmp_outputptr_S = (Tuple<KeyType, PayloadType>*) outputptr_S;
      for(unsigned int k = 0; k < total_sample_count_S; k++){
        sorted_training_sample_S[k] = tmp_outputptr_S[k];
      } 
      args->s_rmi->training_sample = &(sorted_training_sample_S);
      args->s_rmi->training_sample_size = total_sample_count_S;
    }
    #else

    uint32_t total_sample_count_R = 0; 
    if(tid == 0)
    {
      for(int i = 0; i < NUM_THREADS; i++)
        total_sample_count_R += sample_count_R[i];

      std::sort((int64_t *)(args->r_rmi->tmp_training_sample), (int64_t *)(args->r_rmi->tmp_training_sample) + total_sample_count_R - 1);
      args->r_rmi->training_sample = &(args->r_rmi->tmp_training_sample);
      args->r_rmi->training_sample_size = total_sample_count_R;
    }

    uint32_t total_sample_count_S = 0; 
    if(tid == 1)
    {
      for(int i = 0; i < NUM_THREADS; i++)
        total_sample_count_S += sample_count_S[i];

      std::sort((int64_t *)(args->s_rmi->tmp_training_sample), (int64_t *)(args->s_rmi->tmp_training_sample) + total_sample_count_S - 1);
      args->s_rmi->training_sample = &(args->s_rmi->tmp_training_sample);
      args->s_rmi->training_sample_size = total_sample_count_S;
    }
    #endif

    //----------------------------------------------------------//
    //                     TRAIN THE MODELS                     //
    //----------------------------------------------------------//

    if(tid == 0)
    {
       // Stop early if the array is identical
      if (((*(args->r_rmi->training_sample))[0]).key == ((*(args->r_rmi->training_sample))[total_sample_count_R - 1]).key) 
      {
        return;
      }
         
      // Populate the training data for the root model
      vector<vector<vector<learned_sort_for_sort_merge::training_point<KeyType, PayloadType>>>> * r_training_data = args->r_training_data;
      for (unsigned int i = 0; i < total_sample_count_R; ++i) {
        (*r_training_data)[0][0].push_back({((*(args->r_rmi->training_sample))[i]), 1. * i / total_sample_count_R});
      }

      // Train the root model using linear interpolation
      auto *current_training_data = &(*r_training_data)[0][0];
      typename learned_sort_for_sort_merge::RMI<KeyType, PayloadType>::linear_model *current_model = &args->r_rmi->models[0][0];

      // Find the min and max values in the training set
      learned_sort_for_sort_merge::training_point<KeyType, PayloadType> min = current_training_data->front();
      learned_sort_for_sort_merge::training_point<KeyType, PayloadType> max = current_training_data->back();

      // Calculate the slope and intercept terms
      current_model->slope =
          1. / (max.x.key - min.x.key);  // Assuming min.y = 0 and max.y = 1
      current_model->intercept = -current_model->slope * min.x.key;

      // Extrapolate for the number of models in the next layer
      current_model->slope *= args->r_p.arch[1] - 1;
      current_model->intercept *= args->r_p.arch[1] - 1;
#ifndef RUN_LEARNED_TECHNIQUES_WITH_FIRST_LEVEL_ONLY
      // Populate the training data for the next layer
      for (const auto &d : *current_training_data) {
        // Predict the model index in next layer
        unsigned int rank = current_model->slope * d.x.key + current_model->intercept;

        // Normalize the rank between 0 and the number of models in the next layer
        rank =
            std::max(static_cast<unsigned int>(0), std::min(args->r_p.arch[1] - 1, rank));
        
        //if(d.x.key >299 && d.x.key < 400)
            //printf("training: key %ld rank %ld \n", d.x.key, rank);
        // Place the data in the predicted training bucket
        (*r_training_data)[1][rank].push_back(d);
      }

      // Train the leaf models
      for (unsigned int model_idx = 0; model_idx < args->r_p.arch[1]; ++model_idx) {
        // Update iterator variables
        current_training_data = &(*r_training_data)[1][model_idx];
        current_model = &args->r_rmi->models[1][model_idx];

        // Interpolate the min points in the training buckets
        if (model_idx ==
            0) {  // The current model is the first model in the current layer

          if (current_training_data->size() <
              2) {  // Case 1: The first model in this layer is empty
            current_model->slope = 0;
            current_model->intercept = 0;

            // Insert a fictive training point to avoid propagating more than one
            // empty initial models.
            learned_sort_for_sort_merge::training_point<KeyType, PayloadType> tp;
            tp.x.key = 0;
            tp.x.payload = 0;
            tp.y = 0;
            current_training_data->push_back(tp);
          } else {  // Case 2: The first model in this layer is not empty

            min = current_training_data->front();
            max = current_training_data->back();

            current_model->slope =
                (max.y) / (max.x.key - min.x.key);  // Hallucinating as if min.y = 0
            current_model->intercept = min.y - current_model->slope * min.x.key;
          }
        } else if (model_idx == args->r_p.arch[1] - 1) {
          if (current_training_data
                  ->empty()) {  // Case 3: The final model in this layer is empty

            current_model->slope = 0;
            current_model->intercept = 1;
          } else {  // Case 4: The last model in this layer is not empty

            min = (*r_training_data)[1][model_idx - 1].back();
            max = current_training_data->back();

            current_model->slope =
                (min.y - 1) / (min.x.key - max.x.key);  // Hallucinating as if max.y = 1
            current_model->intercept = min.y - current_model->slope * min.x.key;
          }
        } else {  // The current model is not the first model in the current layer

          if (current_training_data
                  ->empty()) {  // Case 5: The intermediate model in
            // this layer is empty
            current_model->slope = 0;
            current_model->intercept =
                (*r_training_data)[1][model_idx - 1].back().y;  // If the previous model
                                                          // was empty too, it will
                                                          // use the fictive
                                                          // training points

            // Insert a fictive training point to avoid propagating more than one
            // empty initial models.
            // NOTE: This will _NOT_ throw to DIV/0 due to identical x's and y's
            // because it is working backwards.
            learned_sort_for_sort_merge::training_point<KeyType, PayloadType> tp;
            tp.x = (*r_training_data)[1][model_idx - 1].back().x;
            tp.y = (*r_training_data)[1][model_idx - 1].back().y;
            current_training_data->push_back(tp);
          } else {  // Case 6: The intermediate leaf model is not empty

            min = (*r_training_data)[1][model_idx - 1].back();
            max = current_training_data->back();

            current_model->slope = (min.y - max.y) / (min.x.key - max.x.key);
            current_model->intercept = min.y - current_model->slope * min.x.key;
            //if(model_idx == 5)
            //printf("min.y %lf max.y %lf min.x.key %ld max.x.key %ld current_model->slope %lf current_model->intercept %lf\n", min.y, max.y, min.x.key, max.x.key, current_model->slope, current_model->intercept);
          }
        }
      }
#endif
      args->r_rmi->trained = true;         
    }


    if(tid == 1)
    {
       // Stop early if the array is identical
      if (((*(args->s_rmi->training_sample))[0]).key == ((*(args->s_rmi->training_sample))[total_sample_count_S - 1]).key) 
      {
        return;
      }
         
      // Populate the training data for the root model
      vector<vector<vector<learned_sort_for_sort_merge::training_point<KeyType, PayloadType>>>> * s_training_data = args->s_training_data;
      for (unsigned int i = 0; i < total_sample_count; ++i) {
        (*s_training_data)[0][0].push_back({((*(args->s_rmi->training_sample))[i]), 1. * i / total_sample_count_S});
      }

      // Train the root model using linear interpolation
      auto *current_training_data = &(*s_training_data)[0][0];
      typename learned_sort_for_sort_merge::RMI<KeyType, PayloadType>::linear_model *current_model = &args->s_rmi->models[0][0];

      // Find the min and max values in the training set
      learned_sort_for_sort_merge::training_point<KeyType, PayloadType> min = current_training_data->front();
      learned_sort_for_sort_merge::training_point<KeyType, PayloadType> max = current_training_data->back();

      // Calculate the slope and intercept terms
      current_model->slope =
          1. / (max.x.key - min.x.key);  // Assuming min.y = 0 and max.y = 1
      current_model->intercept = -current_model->slope * min.x.key;

      // Extrapolate for the number of models in the next layer
      current_model->slope *= args->s_p.arch[1] - 1;
      current_model->intercept *= args->s_p.arch[1] - 1;
#ifndef RUN_LEARNED_TECHNIQUES_WITH_FIRST_LEVEL_ONLY
      // Populate the training data for the next layer
      for (const auto &d : *current_training_data) {
        // Predict the model index in next layer
        unsigned int rank = current_model->slope * d.x.key + current_model->intercept;

        // Normalize the rank between 0 and the number of models in the next layer
        rank =
            std::max(static_cast<unsigned int>(0), std::min(args->s_p.arch[1] - 1, rank));
        
        //if(d.x.key >299 && d.x.key < 400)
            //printf("training: key %ld rank %ld \n", d.x.key, rank);
        // Place the data in the predicted training bucket
        (*s_training_data)[1][rank].push_back(d);
      }

      // Train the leaf models
      for (unsigned int model_idx = 0; model_idx < args->s_p.arch[1]; ++model_idx) {
        // Update iterator variables
        current_training_data = &(*s_training_data)[1][model_idx];
        current_model = &args->s_rmi->models[1][model_idx];

        // Interpolate the min points in the training buckets
        if (model_idx ==
            0) {  // The current model is the first model in the current layer

          if (current_training_data->size() <
              2) {  // Case 1: The first model in this layer is empty
            current_model->slope = 0;
            current_model->intercept = 0;

            // Insert a fictive training point to avoid propagating more than one
            // empty initial models.
            learned_sort_for_sort_merge::training_point<KeyType, PayloadType> tp;
            tp.x.key = 0;
            tp.x.payload = 0;
            tp.y = 0;
            current_training_data->push_back(tp);
          } else {  // Case 2: The first model in this layer is not empty

            min = current_training_data->front();
            max = current_training_data->back();

            current_model->slope =
                (max.y) / (max.x.key - min.x.key);  // Hallucinating as if min.y = 0
            current_model->intercept = min.y - current_model->slope * min.x.key;
          }
        } else if (model_idx == args->s_p.arch[1] - 1) {
          if (current_training_data
                  ->empty()) {  // Case 3: The final model in this layer is empty

            current_model->slope = 0;
            current_model->intercept = 1;
          } else {  // Case 4: The last model in this layer is not empty

            min = (*s_training_data)[1][model_idx - 1].back();
            max = current_training_data->back();

            current_model->slope =
                (min.y - 1) / (min.x.key - max.x.key);  // Hallucinating as if max.y = 1
            current_model->intercept = min.y - current_model->slope * min.x.key;
          }
        } else {  // The current model is not the first model in the current layer

          if (current_training_data
                  ->empty()) {  // Case 5: The intermediate model in
            // this layer is empty
            current_model->slope = 0;
            current_model->intercept =
                (*s_training_data)[1][model_idx - 1].back().y;  // If the previous model
                                                          // was empty too, it will
                                                          // use the fictive
                                                          // training points

            // Insert a fictive training point to avoid propagating more than one
            // empty initial models.
            // NOTE: This will _NOT_ throw to DIV/0 due to identical x's and y's
            // because it is working backwards.
            learned_sort_for_sort_merge::training_point<KeyType, PayloadType> tp;
            tp.x = (*training_data)[1][model_idx - 1].back().x;
            tp.y = (*training_data)[1][model_idx - 1].back().y;
            current_training_data->push_back(tp);
          } else {  // Case 6: The intermediate leaf model is not empty

            min = (*s_training_data)[1][model_idx - 1].back();
            max = current_training_data->back();

            current_model->slope = (min.y - max.y) / (min.x.key - max.x.key);
            current_model->intercept = min.y - current_model->slope * min.x.key;
            //if(model_idx == 5)
            //printf("min.y %lf max.y %lf min.x.key %ld max.x.key %ld current_model->slope %lf current_model->intercept %lf\n", min.y, max.y, min.x.key, max.x.key, current_model->slope, current_model->intercept);
          }
        }
      }
#endif
      args->s_rmi->trained = true;         
    }

}
#endif


void * learned_imv_sort_join_thread(void * param)
{
    LearnedSortMergeMultiwayJoinThread<KeyType, PayloadType> * args   = (LearnedSortMergeMultiwayJoinThread<KeyType, PayloadType> *) param;
    int32_t my_tid = args->my_tid;
    int rv;
    int deltaT = 0; struct timeval t1, t2;

    /*************************************************************************
    *
    *   Phase.0) Sampling and training RMI models.
    *
    *************************************************************************/
    for (int rp = 0; rp < RUN_NUMS; ++rp) 
    {
        if(my_tid == 0){
            #ifdef BUILD_RMI_FROM_TWO_DATASETS
                init_models_training_data_and_sample_counts<KeyType, PayloadType>(args->training_data, args->p.arch, 
                        args->sample_count, args->sample_count_R, args->sample_count_S, NUM_THREADS);
            #else
                init_models_training_data_and_sample_counts<KeyType, PayloadType>(args->r_training_data, args->r_p.arch, 
                        args->sample_count_R, args->sample_count_R, args->sample_count_R, NUM_THREADS);
                init_models_training_data_and_sample_counts<KeyType, PayloadType>(args->s_training_data, args->s_p.arch, 
                        args->sample_count_S, args->sample_count_S, args->sample_count_S, NUM_THREADS);
            #endif
        }

        BARRIER_ARRIVE(args->barrier, rv);

        if(my_tid == 0){
            gettimeofday(&t1, NULL);
        }

        sample_and_train_models_threaded(args);

        BARRIER_ARRIVE(args->barrier, rv);

        if(my_tid == 0){
            gettimeofday(&t2, NULL);

            deltaT = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
            printf("---- Learned sort join sampling and training models time (ms) = %10.4lf\n",  deltaT * 1.0 / 1000);

    #ifndef RUN_LEARNED_TECHNIQUES_WITH_FIRST_LEVEL_ONLY
            if(rp == RUN_NUMS - 1)
            {   
            #ifdef BUILD_RMI_FROM_TWO_DATASETS
                for (unsigned int j = 0; j < args->rmi->hp.arch[1]; ++j) 
                {
                    args->slopes->push_back(args->rmi->models[1][j].slope);
                    args->intercepts->push_back(args->rmi->models[1][j].intercept);
                }
            #else
                for (unsigned int j = 0; j < args->r_rmi->hp.arch[1]; ++j) 
                {
                    args->r_slopes->push_back(args->r_rmi->models[1][j].slope);
                    args->r_intercepts->push_back(args->r_rmi->models[1][j].intercept);
                }
                for (unsigned int j = 0; j < args->s_rmi->hp.arch[1]; ++j) 
                {
                    args->s_slopes->push_back(args->s_rmi->models[1][j].slope);
                    args->s_intercepts->push_back(args->s_rmi->models[1][j].intercept);
                }
            #endif
            } 
    #endif            
        }        
    }

    /*************************************************************************
    *
    *   Phase.1) global partitioning across different threads (i.e., major buckets)
    *
    *************************************************************************/
    for (int rp = 0; rp < RUN_NUMS; ++rp) 
    {
        //DEBUGMSG(1, "Thread-%d started running ... \n", my_tid);

        // wait at a barrier until each thread starts and start timer
        BARRIER_ARRIVE(args->barrier, rv);

        // the first thread checkpoints the start time
        if(my_tid == 0){ 
            gettimeofday(&args->start_time, NULL);
        }
         
        if(my_tid == 0)
        {
            learned_sort_for_sort_merge::partition_major_buckets<KeyType, PayloadType>(1, r_rmi_ptr, NUM_THREADS, args->original_relR->tuples, args->original_relR->num_tuples, 
                                                tmpRelpartR, r_partition_offsets, r_partition_sizes_for_threads, 
                                                tmpRepeatedKeysPredictedRanksR, tmpRepeatedKeysCountsR, r_repeated_keys_offsets, r_repeated_keys_sizes_for_threads, r_total_repeated_keys_sizes_for_threads, -1, -1);

            learned_sort_for_sort_merge::partition_major_buckets<KeyType, PayloadType>(0, s_rmi_ptr, NUM_THREADS, args->original_relS->tuples, args->original_relS->num_tuples, 
                                                tmpRelpartS, s_partition_offsets, s_partition_sizes_for_threads, 
                                                tmpRepeatedKeysPredictedRanksS, tmpRepeatedKeysCountsS, s_repeated_keys_offsets, s_repeated_keys_sizes_for_threads, s_total_repeated_keys_sizes_for_threads, -1, -1);


            for(int i = 0; i < NUM_THREADS; i++)
            {
               printf("r_partition_offsets for thread %d is %d\n", i, r_partition_offsets[i]);
               printf("r_partition_sizes_for_threads for thread %d is %d\n", i, r_partition_sizes_for_threads[i]);
               uint32_t min_key = args->original_relR->num_tuples; uint32_t max_key = 0;
               for(int j = 0; j < r_partition_sizes_for_threads[i]; j++)
               {
                  if(tmpRelpartR[r_partition_offsets[i] + j].key < min_key)
                     min_key = tmpRelpartR[r_partition_offsets[i] + j].key;
                  if(tmpRelpartR[r_partition_offsets[i] + j].key > max_key)
                     max_key = tmpRelpartR[r_partition_offsets[i] + j].key;
                  //if(i == 3)
                   //printf("tmpRelpartR[r_partition_offsets[i] + j].key %ld\n", tmpRelpartR[r_partition_offsets[i] + j].key);
               }
               printf("for thread %d: min_key %ld, max_key %ld \n", i, min_key, max_key);
            }
            
            for(int i = 0; i < NUM_THREADS; i++)
            {
               printf("s_partition_offsets for thread %d is %d\n", i, s_partition_offsets[i]);
               printf("s_partition_sizes_for_threads for thread %d is %d\n", i, s_partition_sizes_for_threads[i]);
               uint32_t min_key = args->original_relS->num_tuples; uint32_t max_key = 0;
               for(int j = 0; j < s_partition_sizes_for_threads[i]; j++)
               {
                  if(tmpRelpartS[s_partition_offsets[i] + j].key < min_key)
                     min_key = tmpRelpartS[s_partition_offsets[i] + j].key;
                  if(tmpRelpartS[s_partition_offsets[i] + j].key > max_key)
                     max_key = tmpRelpartS[s_partition_offsets[i] + j].key;
                  //if(i == 3)
                   //printf("tmpRelpartS[s_partition_offsets[i] + j].key %ld\n", tmpRelpartS[s_partition_offsets[i] + j].key);
               }
               printf("for thread %d: min_key %ld, max_key %ld \n", i, min_key, max_key);
            }
        }

        // wait at a barrier until each thread completes the partition phase
        BARRIER_ARRIVE(args->barrier, rv);

        // partition phase finished, thread-0 checkpoints the time
        if(my_tid == 0){
            gettimeofday(&args->partition_end_time, NULL);

            deltaT = (args->partition_end_time.tv_sec - args->start_time.tv_sec) * 1000000 + args->partition_end_time.tv_usec - args->start_time.tv_usec;
            printf("---- %5s partitioning costs time (ms) = %10.4lf\n", "Learned sort join", deltaT * 1.0 / 1000);
        }
    
        if(!(rp == (RUN_NUMS - 1))){
            //TODO: make sure you can rerun here
            //if(args->tid == 0)
            //    destroy_hashtable(args->ht);

            //free_bucket_buffer(overflowbuf);

            BARRIER_ARRIVE(args->barrier, rv);
        } 
    }

    BARRIER_ARRIVE(args->barrier, rv);

    /*************************************************************************
    *
    *   Phase.2) NUMA-local sorting of cache-sized chunks
    *
    *************************************************************************/
    for (int rp = 0; rp < RUN_NUMS; ++rp) 
    {
        BARRIER_ARRIVE(args->barrier, rv);

        // the first thread checkpoints the start time
        if(my_tid == 0){ 
            gettimeofday(&args->partition_end_time, NULL);
        }

        join_steps.sorting_phase(0, 0, args);

        BARRIER_ARRIVE(args->barrier, rv);

        // partition phase finished, thread-0 checkpoints the time
        if(my_tid == 0){
            gettimeofday(&args->sort_end_time, NULL);

            deltaT = (args->sort_end_time.tv_sec - args->partition_end_time.tv_sec) * 1000000 + args->sort_end_time.tv_usec - args->partition_end_time.tv_usec;
            printf("---- %5s sorting costs time (ms) = %10.4lf\n", "Learned sort join", deltaT * 1.0 / 1000);
        }
    
        if(!(rp == (RUN_NUMS - 1))){
            for(int j = 0; j < args->NUM_MINOR_BCKT_PER_MAJOR_BCKT_r; j++)
                args->tmp_minor_bckt_sizes_r[j] = 0;
            for(int j = 0; j < args->NUM_MINOR_BCKT_PER_MAJOR_BCKT_s; j++)
                args->tmp_minor_bckt_sizes_s[j] = 0;

            BARRIER_ARRIVE(args->barrier, rv);
        }
    }

    BARRIER_ARRIVE(args->barrier, rv);

    /*************************************************************************
     *
     *   Phase.4) NUMA-local merge-join on local sorted runs.
     *
     *************************************************************************/
    for (int rp = 0; rp < RUN_NUMS; ++rp) 
    {
        BARRIER_ARRIVE(args->barrier, rv);

        // the first thread checkpoints the start time
        if(my_tid == 0){ 
            gettimeofday(&args->tmp_mergejoin_end_time, NULL);
        }

        join_steps.mergejoin_phase(0, 0, 0, 0, args);

        BARRIER_ARRIVE(args->barrier, rv);

        // partition phase finished, thread-0 checkpoints the time
        if(my_tid == 0){
            gettimeofday(&args->mergejoin_end_time, NULL);

            deltaT = (args->mergejoin_end_time.tv_sec - args->tmp_mergejoin_end_time.tv_sec) * 1000000 + args->mergejoin_end_time.tv_usec - args->tmp_mergejoin_end_time.tv_usec;
            printf("---- %5s joining costs time (ms) = %10.4lf\n", "Learned sort join", deltaT * 1.0 / 1000);
        }
    }

    return 0;
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
    string curr_rel_r_path = RELATION_R_PATH;
    string curr_rel_s_path = RELATION_S_PATH;

    load_relation<KeyType, PayloadType>(&rel_r, curr_rel_r_path.c_str(), curr_num_tuples_r);
    load_relation<KeyType, PayloadType>(&rel_s, curr_rel_s_path.c_str(), curr_num_tuples_s);    
#else
    // creating new datasets on-the-flay 
    result = create_eth_workload_relation_pk<KeyType, PayloadType>(&rel_r, curr_num_tuples_r, 0);
    //ASSERT_EQ(result, 0);
    #ifdef PERSIST_RELATIONS_FOR_EVALUATION
    write_relation<KeyType, PayloadType>(&rel_r, rel_r_path.c_str());
    #endif
    
    result = create_eth_workload_relation_pk<KeyType, PayloadType>(&rel_s, curr_num_tuples_s, 0);
    //ASSERT_EQ(result, 0);
    #ifdef PERSIST_RELATIONS_FOR_EVALUATION
    write_relation<KeyType, PayloadType>(&rel_s, rel_s_path.c_str());
    #endif
#endif

    int i, rv;
    pthread_barrier_t barrier;
    Result * joinresult;
    pthread_t tid[NUM_THREADS];
    pthread_attr_t attr;
    cpu_set_t set;
    
    joinresult = (Result *) malloc(sizeof(Result));
    joinresult->resultlist = (ThreadResult *) malloc(sizeof(ThreadResult) * NUM_THREADS);

    joincfg.NTHREADS = NUM_THREADS;
    joincfg.LEARNEDSORT = USE_LEARNED_SORT;
    joincfg.PARTFANOUT = ETH_SORT_MERGE_PARTITION_FANOUT_PER_THREAD;
    joincfg.SCALARMERGE = ETH_SORT_MERGE_IS_SCALAR_MERGE;
    joincfg.MWAYMERGEBUFFERSIZE = MWAY_MERGE_BUFFER_SIZE_DEFAULT;
    joincfg.NUMASTRATEGY = ETH_SORT_MERGE_NUMA_STRATEGY;
    numa_shuffle_init(joincfg.NUMASTRATEGY, joincfg.NTHREADS);

    /* check whether nr. of threads is a power of 2 */
    if((joincfg.NTHREADS & (joincfg.NTHREADS-1)) != 0)
    {
        perror("[ERROR] Learned sort-merge join runs with a power of 2 #threads.");
        exit(EXIT_FAILURE);
    }

    //CACHELINEPADDING = ((ETH_SORT_MERGE_PARTITION_FANOUT_PER_THREAD) * CACHE_LINE_SIZE/sizeof(Tuple<KeyType, PayloadType>));
    CACHELINEPADDING = 0; //NOTE: It should be always zero in the case of sortmerge_multiway_join_learned
    RELATION_PADDING = ((NUM_THREADS) * CACHELINEPADDING * sizeof(Tuple<KeyType, PayloadType>));


    //////////////////////////////////////////////////////////////////////////////
    // start stuff for sampling and building RMI models for both relations R and S
    //////////////////////////////////////////////////////////////////////////////

    Tuple<KeyType, PayloadType>* r_tmp_training_sample_in;
    Tuple<KeyType, PayloadType>* r_sorted_training_sample_in;
    Tuple<KeyType, PayloadType>* s_tmp_training_sample_in;
    Tuple<KeyType, PayloadType>* s_sorted_training_sample_in;
    unsigned int SAMPLE_SZ_R, SAMPLE_SZ_S;

#ifdef BUILD_RMI_FROM_TWO_DATASETS

    // Sampling and building RMI models for relations R and S together
    typename learned_sort_for_sort_merge::RMI<KeyType, PayloadType>::Params rmi_params;
    learned_sort_for_sort_merge::validate_params<KeyType, PayloadType>(rmi_params, rel_r.num_tuples);
    learned_sort_for_sort_merge::validate_params<KeyType, PayloadType>(rmi_params, rel_s.num_tuples);
    SAMPLE_SZ_R = std::min<unsigned int>(
        rel_r.num_tuples, std::max<unsigned int>(rmi_params.sampling_rate * rel_r.num_tuples,
                                        learned_sort_for_sort_merge::RMI<KeyType, PayloadType>::Params::MIN_SORTING_SIZE)) + 1;
    r_tmp_training_sample_in = (Tuple<KeyType, PayloadType>*) alloc_aligned(SAMPLE_SZ_R * sizeof(Tuple<KeyType, PayloadType>));
    #ifdef USE_AVXSORT_AS_STD_SORT
    r_sorted_training_sample_in = (Tuple<KeyType, PayloadType>*) alloc_aligned(SAMPLE_SZ_R * sizeof(Tuple<KeyType, PayloadType>));
    #endif
    SAMPLE_SZ_S = std::min<unsigned int>(
        rel_s.num_tuples, std::max<unsigned int>(rmi_params.sampling_rate * rel_s.num_tuples,
                                        learned_sort_for_sort_merge::RMI<KeyType, PayloadType>::Params::MIN_SORTING_SIZE)) + 1;
    s_tmp_training_sample_in = (Tuple<KeyType, PayloadType>*) alloc_aligned(SAMPLE_SZ_S * sizeof(Tuple<KeyType, PayloadType>));
    #ifdef USE_AVXSORT_AS_STD_SORT
    s_sorted_training_sample_in = (Tuple<KeyType, PayloadType>*) alloc_aligned(SAMPLE_SZ_S * sizeof(Tuple<KeyType, PayloadType>));
    #endif
    
    Tuple<KeyType, PayloadType>* tmp_training_sample_in;
    Tuple<KeyType, PayloadType>* sorted_training_sample_in;
    tmp_training_sample_in = (Tuple<KeyType, PayloadType>*) alloc_aligned((SAMPLE_SZ_R + SAMPLE_SZ_S) * sizeof(Tuple<KeyType, PayloadType>));
    #ifdef USE_AVXSORT_AS_STD_SORT
    sorted_training_sample_in = (Tuple<KeyType, PayloadType>*) alloc_aligned((SAMPLE_SZ_R + SAMPLE_SZ_S) * sizeof(Tuple<KeyType, PayloadType>));
    #endif

    learned_sort_for_sort_merge::RMI<KeyType, PayloadType> rmi(rmi_params, tmp_training_sample_in, sorted_training_sample_in,
                                              r_tmp_training_sample_in, r_sorted_training_sample_in,
                                              s_tmp_training_sample_in, s_sorted_training_sample_in);

    vector<vector<vector<learned_sort_for_sort_merge::training_point<KeyType, PayloadType>>>> training_data(rmi_params.arch.size());
    for (unsigned int layer_idx = 0; layer_idx < rmi_params.arch.size(); ++layer_idx) {
        training_data[layer_idx].resize(rmi_params.arch[layer_idx]);
    }

    uint32_t * sample_count = (uint32_t *) calloc(NUM_THREADS, sizeof(uint32_t)); 
    uint32_t * sample_count_R = (uint32_t *) calloc(NUM_THREADS, sizeof(uint32_t)); 
    uint32_t * sample_count_S = (uint32_t *) calloc(NUM_THREADS, sizeof(uint32_t));

    vector<double>* slopes = new vector<double>;                 
    vector<double>* intercepts = new vector<double>;

    rmi_ptr = &rmi;
    r_rmi_ptr = &rmi;
    s_rmi_ptr = &rmi;

#else
    // Sampling and building RMI models for relation R and S seperately
    typename learned_sort_for_sort_merge::RMI<KeyType, PayloadType>::Params r_rmi_params;
    learned_sort_for_sort_merge::validate_params<KeyType, PayloadType>(r_rmi_params, rel_r.num_tuples);
    SAMPLE_SZ_R = std::min<unsigned int>(
        rel_r.num_tuples, std::max<unsigned int>(r_rmi_params.sampling_rate * rel_r.num_tuples,
                                        learned_sort_for_sort_merge::RMI<KeyType, PayloadType>::Params::MIN_SORTING_SIZE)) + 1;

    r_tmp_training_sample_in = (Tuple<KeyType, PayloadType>*) alloc_aligned(SAMPLE_SZ_R * sizeof(Tuple<KeyType, PayloadType>));
    #ifdef USE_AVXSORT_AS_STD_SORT
    r_sorted_training_sample_in = (Tuple<KeyType, PayloadType>*) alloc_aligned(SAMPLE_SZ_R * sizeof(Tuple<KeyType, PayloadType>));
    #endif

    learned_sort_for_sort_merge::RMI<KeyType, PayloadType> r_rmi(r_rmi_params, r_tmp_training_sample_in, r_sorted_training_sample_in);

    typename learned_sort_for_sort_merge::RMI<KeyType, PayloadType>::Params s_rmi_params;
    learned_sort_for_sort_merge::validate_params<KeyType, PayloadType>(s_rmi_params, rel_s.num_tuples);
    SAMPLE_SZ_S = std::min<unsigned int>(
        rel_s.num_tuples, std::max<unsigned int>(s_rmi_params.sampling_rate * rel_s.num_tuples,
                                        learned_sort_for_sort_merge::RMI<KeyType, PayloadType>::Params::MIN_SORTING_SIZE)) + 1;

    s_tmp_training_sample_in = (Tuple<KeyType, PayloadType>*) alloc_aligned(SAMPLE_SZ_S * sizeof(Tuple<KeyType, PayloadType>));
    #ifdef USE_AVXSORT_AS_STD_SORT
    s_sorted_training_sample_in = (Tuple<KeyType, PayloadType>*) alloc_aligned(SAMPLE_SZ_S * sizeof(Tuple<KeyType, PayloadType>));
    #endif

    learned_sort_for_sort_merge::RMI<KeyType, PayloadType> s_rmi(s_rmi_params, s_tmp_training_sample_in, s_sorted_training_sample_in);
    
    
    vector<vector<vector<learned_sort_for_sort_merge::training_point<KeyType, PayloadType>>>> r_training_data(r_rmi_params.arch.size());
    for (unsigned int layer_idx = 0; layer_idx < r_rmi_params.arch.size(); ++layer_idx) {
        r_training_data[layer_idx].resize(r_rmi_params.arch[layer_idx]);
    }
    vector<vector<vector<learned_sort_for_sort_merge::training_point<KeyType, PayloadType>>>> s_training_data(s_rmi_params.arch.size());
    for (unsigned int layer_idx = 0; layer_idx < s_rmi_params.arch.size(); ++layer_idx) {
        s_training_data[layer_idx].resize(s_rmi_params.arch[layer_idx]);
    }

    uint32_t * sample_count_R = (uint32_t *) calloc(NUM_THREADS, sizeof(uint32_t)); 
    uint32_t * sample_count_S = (uint32_t *) calloc(NUM_THREADS, sizeof(uint32_t));

    vector<double>* r_slopes = new vector<double>;                 
    vector<double>* r_intercepts = new vector<double>;
    vector<double>* s_slopes = new vector<double>;                 
    vector<double>* s_intercepts = new vector<double>;

    r_rmi_ptr = &r_rmi;
    s_rmi_ptr = &s_rmi;    

#endif        
    //////////////////////////////////////////////////////////////////////////////
    // End stuff for sampling and building RMI models for both relations R and S
    //////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Determining the size of partitions for the different threads annd allocating required data structures
    //////////////////////////////////////////////////////////////////////////////////////////////////////// 
    r_initial_partition_sizes_for_threads = (int64_t*) calloc(NUM_THREADS, sizeof(int64_t));
    s_initial_partition_sizes_for_threads = (int64_t*) calloc(NUM_THREADS, sizeof(int64_t));
    int64_t total_r_initial_partition_sizes_for_threads = 0; int64_t total_s_initial_partition_sizes_for_threads = 0;
    //learned_sort_for_sort_merge::histogram_and_get_max_capacity_for_major_buckets(r_initial_partition_sizes_for_threads, NUM_THREADS, r_rmi_ptr, rel_r.tuples, rel_r.num_tuples, -1, -1);
    //learned_sort_for_sort_merge::histogram_and_get_max_capacity_for_major_buckets(s_initial_partition_sizes_for_threads, NUM_THREADS, s_rmi_ptr, rel_s.tuples, rel_s.num_tuples, -1, -1);
    for(i=0; i<NUM_THREADS; i++)
    {
        r_initial_partition_sizes_for_threads[i] = (rel_r.num_tuples/NUM_THREADS) * OVERALLOCATION_SIZE_RATIO + CACHELINEPADDING;
        s_initial_partition_sizes_for_threads[i] = (rel_s.num_tuples/NUM_THREADS) * OVERALLOCATION_SIZE_RATIO + CACHELINEPADDING;
        total_r_initial_partition_sizes_for_threads += r_initial_partition_sizes_for_threads[i];
        total_s_initial_partition_sizes_for_threads += s_initial_partition_sizes_for_threads[i];
    }
    for(i=0; i<NUM_THREADS; i++)
    {
        printf("R partition for thread %d has %d elements\n", i, r_initial_partition_sizes_for_threads[i]);
        printf("S partition for thread %d has %d elements\n", i, s_initial_partition_sizes_for_threads[i]);
    } 

    r_NUM_MINOR_BCKT_PER_MAJOR_BCKT = (unsigned int *) calloc(NUM_THREADS, sizeof(unsigned int));
    s_NUM_MINOR_BCKT_PER_MAJOR_BCKT = (unsigned int *) calloc(NUM_THREADS, sizeof(unsigned int));
    r_MINOR_BCKTS_OFFSET = (unsigned int *) calloc(NUM_THREADS, sizeof(unsigned int)); 
    s_MINOR_BCKTS_OFFSET = (unsigned int *) calloc(NUM_THREADS, sizeof(unsigned int));
    r_TOT_NUM_MINOR_BCKTS = 0;
    s_TOT_NUM_MINOR_BCKTS = 0;

    static const double r_OA_RATIO = r_rmi_ptr->hp.overallocation_ratio;
    static const double s_OA_RATIO = s_rmi_ptr->hp.overallocation_ratio;
    static const unsigned int r_THRESHOLD = r_rmi_ptr->hp.threshold;
    static const unsigned int s_THRESHOLD = s_rmi_ptr->hp.threshold;
    unsigned int r_curr_minor_bckts_offset = 0; unsigned int s_curr_minor_bckts_offset = 0; 
    for(i = 0; i < NUM_THREADS; i++)
    {
        /*r_NUM_MINOR_BCKT_PER_MAJOR_BCKT[i] = std::max(1u, 
                            static_cast<unsigned>(r_initial_partition_sizes_for_threads[i] * r_OA_RATIO / r_THRESHOLD));
        s_NUM_MINOR_BCKT_PER_MAJOR_BCKT[i] = std::max(1u, 
                            static_cast<unsigned>(s_initial_partition_sizes_for_threads[i] * s_OA_RATIO / s_THRESHOLD));*/
        r_NUM_MINOR_BCKT_PER_MAJOR_BCKT[i] = std::max(1u, 
                            static_cast<unsigned>((rel_r.num_tuples/NUM_THREADS) * r_OA_RATIO / r_THRESHOLD));
        s_NUM_MINOR_BCKT_PER_MAJOR_BCKT[i] = std::max(1u, 
                            static_cast<unsigned>((rel_s.num_tuples/NUM_THREADS) * s_OA_RATIO / s_THRESHOLD));
        r_TOT_NUM_MINOR_BCKTS += r_NUM_MINOR_BCKT_PER_MAJOR_BCKT[i];
        s_TOT_NUM_MINOR_BCKTS += s_NUM_MINOR_BCKT_PER_MAJOR_BCKT[i];
        r_MINOR_BCKTS_OFFSET[i] = r_curr_minor_bckts_offset;
        s_MINOR_BCKTS_OFFSET[i] = s_curr_minor_bckts_offset;
        r_curr_minor_bckts_offset += r_NUM_MINOR_BCKT_PER_MAJOR_BCKT[i];
        s_curr_minor_bckts_offset += s_NUM_MINOR_BCKT_PER_MAJOR_BCKT[i];
    }

    //for(i = 0; i < NUM_THREADS; i++)
    //{
    //    printf("for %d, r_NUM_MINOR_BCKT_PER_MAJOR_BCKT[i] %ld s_NUM_MINOR_BCKT_PER_MAJOR_BCKT[i] %ld \n", i, r_NUM_MINOR_BCKT_PER_MAJOR_BCKT[i], s_NUM_MINOR_BCKT_PER_MAJOR_BCKT[i]);
    //    printf("for %d, r_MINOR_BCKTS_OFFSET[i] %ld s_MINOR_BCKTS_OFFSET[i] %ld \n", i, r_MINOR_BCKTS_OFFSET[i], s_MINOR_BCKTS_OFFSET[i]);
    //    printf("for %d, r_TOT_NUM_MINOR_BCKTS %ld s_TOT_NUM_MINOR_BCKTS %ld \n", i, r_TOT_NUM_MINOR_BCKTS, s_TOT_NUM_MINOR_BCKTS);
    //}

    // allocate temporary space for partitioning and just to make sure that chunks of the temporary memory will be numa local to threads.
    tmpRelpartR = NULL; tmpRelpartS = NULL;
    tmpRelpartR = (Tuple<KeyType, PayloadType>*) alloc_aligned(OVERALLOCATION_SIZE_RATIO * total_r_initial_partition_sizes_for_threads * sizeof(Tuple<KeyType, PayloadType>)
                                    + RELATION_PADDING);
    tmpRelpartS = (Tuple<KeyType, PayloadType>*) alloc_aligned(OVERALLOCATION_SIZE_RATIO * total_s_initial_partition_sizes_for_threads * sizeof(Tuple<KeyType, PayloadType>)
                                    + RELATION_PADDING);

    numa_localize_varlen<KeyType, PayloadType>(tmpRelpartR, r_initial_partition_sizes_for_threads, NUM_THREADS);
    numa_localize_varlen<KeyType, PayloadType>(tmpRelpartS, s_initial_partition_sizes_for_threads, NUM_THREADS);

    r_partition_offsets = (int64_t*) calloc(NUM_THREADS, sizeof(int64_t));
    s_partition_offsets = (int64_t*) calloc(NUM_THREADS, sizeof(int64_t));
    int64_t r_rel_offset = 0; int64_t s_rel_offset = 0; 
    for(i = 0; i < NUM_THREADS; i++)
    {
        r_partition_offsets[i] = r_rel_offset;
        s_partition_offsets[i] = s_rel_offset;
        r_rel_offset += r_initial_partition_sizes_for_threads[i];
        s_rel_offset += s_initial_partition_sizes_for_threads[i];            
    }
    r_partition_sizes_for_threads = (int64_t*) calloc(NUM_THREADS, sizeof(int64_t));
    s_partition_sizes_for_threads = (int64_t*) calloc(NUM_THREADS, sizeof(int64_t));

    // allocate temporary space for repeated keys themselves and just to make sure that chunks of the temporary memory will be numa local to threads.
    r_initial_repeated_keys_sizes_for_threads = (int64_t*) calloc(NUM_THREADS, sizeof(int64_t));
    s_initial_repeated_keys_sizes_for_threads = (int64_t*) calloc(NUM_THREADS, sizeof(int64_t));
    int64_t total_r_initial_repeated_keys_size = 0; int64_t total_s_initial_repeated_keys_size = 0;
    for(i=0; i<NUM_THREADS; i++)
    {
        r_initial_repeated_keys_sizes_for_threads[i] = (REPEATED_KEYS_SIZE_RATIO * r_initial_partition_sizes_for_threads[i]) + 1 ;    
        s_initial_repeated_keys_sizes_for_threads[i] = (REPEATED_KEYS_SIZE_RATIO * s_initial_partition_sizes_for_threads[i]) + 1 ;      
        total_r_initial_repeated_keys_size += r_initial_repeated_keys_sizes_for_threads[i];
        total_s_initial_repeated_keys_size += s_initial_repeated_keys_sizes_for_threads[i];                
    }
    for(i=0; i<NUM_THREADS; i++)
    {
        r_initial_repeated_keys_sizes_for_threads[i] = r_initial_repeated_keys_sizes_for_threads[i] + CACHELINEPADDING;    
        s_initial_repeated_keys_sizes_for_threads[i] = s_initial_repeated_keys_sizes_for_threads[i] + CACHELINEPADDING;      
    }
    tmpRepeatedKeysPredictedRanksR = NULL; tmpRepeatedKeysPredictedRanksS = NULL;
    tmpRepeatedKeysPredictedRanksR = (Tuple<KeyType, PayloadType>*) alloc_aligned(total_r_initial_repeated_keys_size * sizeof(Tuple<KeyType, PayloadType>)
                                    + RELATION_PADDING);
    tmpRepeatedKeysPredictedRanksS = (Tuple<KeyType, PayloadType>*) alloc_aligned(total_s_initial_repeated_keys_size * sizeof(Tuple<KeyType, PayloadType>)
                                    + RELATION_PADDING);

    numa_localize_varlen<KeyType, PayloadType>(tmpRepeatedKeysPredictedRanksR, r_initial_repeated_keys_sizes_for_threads, NUM_THREADS);
    numa_localize_varlen<KeyType, PayloadType>(tmpRepeatedKeysPredictedRanksS, s_initial_repeated_keys_sizes_for_threads, NUM_THREADS);

    tmpRepeatedKeysCountsR = NULL; tmpRepeatedKeysCountsS = NULL;
    tmpRepeatedKeysCountsR = (int64_t*) alloc_aligned(total_r_initial_repeated_keys_size * sizeof(int64_t)
                                    + RELATION_PADDING);
    tmpRepeatedKeysCountsS = (int64_t*) alloc_aligned(total_s_initial_repeated_keys_size * sizeof(int64_t)
                                    + RELATION_PADDING);

    numa_localize_generic_varlen<int64_t>(tmpRepeatedKeysCountsR, r_initial_repeated_keys_sizes_for_threads, NUM_THREADS);
    numa_localize_generic_varlen<int64_t>(tmpRepeatedKeysCountsS, s_initial_repeated_keys_sizes_for_threads, NUM_THREADS);
    
    r_repeated_keys_offsets = (int64_t *) calloc(NUM_THREADS, sizeof(int64_t));
    s_repeated_keys_offsets = (int64_t *) calloc(NUM_THREADS, sizeof(int64_t));
    int64_t r_repeated_keys_offset = 0; int64_t s_repeated_keys_offset = 0; 
    for(i = 0; i < NUM_THREADS; i++)
    {
        r_repeated_keys_offsets[i] = r_repeated_keys_offset;
        s_repeated_keys_offsets[i] = s_repeated_keys_offset;
        r_repeated_keys_offset += r_initial_repeated_keys_sizes_for_threads[i];
        s_repeated_keys_offset += s_initial_repeated_keys_sizes_for_threads[i];            
    }
    r_repeated_keys_sizes_for_threads = (int64_t*) calloc(NUM_THREADS, sizeof(int64_t));
    s_repeated_keys_sizes_for_threads = (int64_t*) calloc(NUM_THREADS, sizeof(int64_t));
    r_total_repeated_keys_sizes_for_threads = (int64_t*) calloc(NUM_THREADS, sizeof(int64_t));
    s_total_repeated_keys_sizes_for_threads = (int64_t*) calloc(NUM_THREADS, sizeof(int64_t));

    // allocate temporary space for sorting and just to make sure that chunks of the temporary memory will be numa local to threads.
    tmpRelsortR = NULL; tmpRelsortS = NULL;
    tmpRelsortR = (Tuple<KeyType, PayloadType>*) alloc_aligned(OVERALLOCATION_SIZE_RATIO * total_r_initial_partition_sizes_for_threads * sizeof(Tuple<KeyType, PayloadType>)
                                    + RELATION_PADDING);
    tmpRelsortS = (Tuple<KeyType, PayloadType>*) alloc_aligned(OVERALLOCATION_SIZE_RATIO * total_s_initial_partition_sizes_for_threads * sizeof(Tuple<KeyType, PayloadType>)
                                    + RELATION_PADDING);

    numa_localize_varlen<KeyType, PayloadType>(tmpRelsortR, r_initial_partition_sizes_for_threads, NUM_THREADS);
    numa_localize_varlen<KeyType, PayloadType>(tmpRelsortS, s_initial_partition_sizes_for_threads, NUM_THREADS);

    r_tmp_spill_bucket = NULL; r_sorted_spill_bucket = NULL;
    r_tmp_spill_bucket = (Tuple<KeyType, PayloadType>*) alloc_aligned(SPILL_BUCKET_SIZE_RATIO * total_r_initial_partition_sizes_for_threads * sizeof(Tuple<KeyType, PayloadType>)
                                    + RELATION_PADDING);
    r_sorted_spill_bucket = (Tuple<KeyType, PayloadType>*) alloc_aligned(SPILL_BUCKET_SIZE_RATIO * total_r_initial_partition_sizes_for_threads * sizeof(Tuple<KeyType, PayloadType>)
                                    + RELATION_PADDING);

    numa_localize_varlen<KeyType, PayloadType>(r_tmp_spill_bucket, r_initial_partition_sizes_for_threads, NUM_THREADS);
    numa_localize_varlen<KeyType, PayloadType>(r_sorted_spill_bucket, r_initial_partition_sizes_for_threads, NUM_THREADS);

    s_tmp_spill_bucket = NULL; s_sorted_spill_bucket = NULL;
    s_tmp_spill_bucket = (Tuple<KeyType, PayloadType>*) alloc_aligned(SPILL_BUCKET_SIZE_RATIO * total_s_initial_partition_sizes_for_threads * sizeof(Tuple<KeyType, PayloadType>)
                                    + RELATION_PADDING);
    s_sorted_spill_bucket = (Tuple<KeyType, PayloadType>*) alloc_aligned(SPILL_BUCKET_SIZE_RATIO * total_s_initial_partition_sizes_for_threads * sizeof(Tuple<KeyType, PayloadType>)
                                    + RELATION_PADDING);

    numa_localize_varlen<KeyType, PayloadType>(s_tmp_spill_bucket, s_initial_partition_sizes_for_threads, NUM_THREADS);
    numa_localize_varlen<KeyType, PayloadType>(s_sorted_spill_bucket, s_initial_partition_sizes_for_threads, NUM_THREADS);

    r_tmp_minor_bckts = NULL; s_tmp_minor_bckts = NULL; 
    int64_t r_minor_bckts_sizes_for_threads[NUM_THREADS]; int64_t s_minor_bckts_sizes_for_threads[NUM_THREADS];
    r_minor_bckts_offsets = (int64_t *) calloc(NUM_THREADS, sizeof(int64_t));
    s_minor_bckts_offsets = (int64_t *) calloc(NUM_THREADS, sizeof(int64_t));
    int64_t curr_r_minor_bckts_offset = 0; int64_t curr_s_minor_bckts_offset = 0;
    for(i = 0; i < NUM_THREADS; i++)
    {
        r_minor_bckts_sizes_for_threads[i] = r_NUM_MINOR_BCKT_PER_MAJOR_BCKT[i] * r_THRESHOLD;
        s_minor_bckts_sizes_for_threads[i] = s_NUM_MINOR_BCKT_PER_MAJOR_BCKT[i] * s_THRESHOLD;
        r_minor_bckts_offsets[i] = curr_r_minor_bckts_offset;
        s_minor_bckts_offsets[i] = curr_s_minor_bckts_offset;
        curr_r_minor_bckts_offset += r_minor_bckts_sizes_for_threads[i];
        curr_s_minor_bckts_offset += s_minor_bckts_sizes_for_threads[i];
    }
    r_tmp_minor_bckts = (Tuple<KeyType, PayloadType>*) alloc_aligned(r_TOT_NUM_MINOR_BCKTS * r_THRESHOLD * sizeof(Tuple<KeyType, PayloadType>)
                                    + RELATION_PADDING);
    s_tmp_minor_bckts = (Tuple<KeyType, PayloadType>*) alloc_aligned(s_TOT_NUM_MINOR_BCKTS * s_THRESHOLD * sizeof(Tuple<KeyType, PayloadType>)
                                    + RELATION_PADDING);
    
    numa_localize_varlen<KeyType, PayloadType>(r_tmp_minor_bckts, r_minor_bckts_sizes_for_threads, NUM_THREADS);
    numa_localize_varlen<KeyType, PayloadType>(s_tmp_minor_bckts, s_minor_bckts_sizes_for_threads, NUM_THREADS);

    r_tmp_minor_bckt_sizes = NULL; s_tmp_minor_bckt_sizes = NULL; 
    int64_t r_minor_bckt_sizes_for_numa_initialization[NUM_THREADS]; int64_t s_minor_bckt_sizes_for_numa_initialization[NUM_THREADS];
    r_minor_bckt_sizes_offsets = (int64_t *) calloc(NUM_THREADS, sizeof(int64_t));
    s_minor_bckt_sizes_offsets = (int64_t *) calloc(NUM_THREADS, sizeof(int64_t));
    int64_t curr_r_minor_bckt_sizes_offset = 0; int64_t curr_s_minor_bckt_sizes_offset = 0;
    for(i = 0; i < NUM_THREADS; i++)
    {
        r_minor_bckt_sizes_for_numa_initialization[i] = r_NUM_MINOR_BCKT_PER_MAJOR_BCKT[i];
        s_minor_bckt_sizes_for_numa_initialization[i] = s_NUM_MINOR_BCKT_PER_MAJOR_BCKT[i];
        r_minor_bckt_sizes_offsets[i] = curr_r_minor_bckt_sizes_offset;
        s_minor_bckt_sizes_offsets[i] = curr_s_minor_bckt_sizes_offset;
        curr_r_minor_bckt_sizes_offset += r_NUM_MINOR_BCKT_PER_MAJOR_BCKT[i];
        curr_s_minor_bckt_sizes_offset += s_NUM_MINOR_BCKT_PER_MAJOR_BCKT[i];
    }
    r_tmp_minor_bckt_sizes = (int64_t *) alloc_aligned(r_TOT_NUM_MINOR_BCKTS * sizeof(int64_t)
                                    + RELATION_PADDING);
    s_tmp_minor_bckt_sizes = (int64_t *) alloc_aligned(s_TOT_NUM_MINOR_BCKTS * sizeof(int64_t)
                                    + RELATION_PADDING);
    /*for(i = 0; i < NUM_THREADS; i++)
    {
        for(int j = 0; j < r_NUM_MINOR_BCKT_PER_MAJOR_BCKT[i]; j++)
            r_tmp_minor_bckt_sizes[r_minor_bckt_sizes_offsets[i] + j] = 0;
        for(int j = 0; j < s_NUM_MINOR_BCKT_PER_MAJOR_BCKT[i]; j++)
            s_tmp_minor_bckt_sizes[s_minor_bckt_sizes_offsets[i] + j] = 0;
    }*/

    numa_localize_generic_varlen<int64_t>(r_tmp_minor_bckt_sizes, r_minor_bckt_sizes_for_numa_initialization, NUM_THREADS);
    numa_localize_generic_varlen<int64_t>(s_tmp_minor_bckt_sizes, s_minor_bckt_sizes_for_numa_initialization, NUM_THREADS);

    //////////////////////////////////////////////////
    // Initialize threads and required data structures
    //////////////////////////////////////////////////
    threadrelchunks = (RelationPair<KeyType, PayloadType> **) malloc(NUM_THREADS *
                                                sizeof(RelationPair<KeyType, PayloadType>*));


    LearnedSortMergeMultiwayJoinThread<KeyType, PayloadType> args[NUM_THREADS];
    LearnedSortMergeMultiwayJoinThread<KeyType, PayloadType> * args_ptr = args;

    rv = pthread_barrier_init(&barrier, NULL, NUM_THREADS);
    if(rv != 0){
        printf("[ERROR] Couldn't create the barrier\n");
        exit(EXIT_FAILURE);
    }

    pthread_attr_init(&attr);

    int             err = 0;
    size_t          stackSize = 0;

    // Get the default value
    err = pthread_attr_getstacksize(&attr, &stackSize);
    if (err) {
        perror("[ERROR] pthread stack size could not be get!");
        exit(0);
    }

    // set the attribute with our required value
    if (stackSize < REQUIRED_STACK_SIZE) {
        err = pthread_attr_setstacksize (&attr, REQUIRED_STACK_SIZE);
        if (err) {
            perror("[ERROR] pthread stack size could not be set!");
            exit(0);
        }
    }

    #ifdef DEVELOPMENT_MODE
    int num_numa = get_num_numa_regions_develop();
    #else
    int num_numa = get_num_numa_regions_v2();
    #endif  

    ptrs_to_sharedmergebufs = (Tuple<KeyType, PayloadType> **)
        alloc_aligned(num_numa*sizeof(Tuple<KeyType, PayloadType>*));
    
    #ifdef SKEW_HANDLING
    ptrs_to_taskqueues = (taskqueue_t **) alloc_aligned(num_numa*sizeof(taskqueue_t*));
    ptrs_to_taskqueues_locks = (pthread_mutex_t*) malloc(num_numa*sizeof(pthread_mutex_t));
    ptrs_to_is_numa_taskqueues_created = (int*) malloc(num_numa*sizeof(int));
    for(i = 0; i < num_numa; i++)
        ptrs_to_is_numa_taskqueues_created[i] = 0;
    #endif


#ifdef BUILD_RMI_FROM_TWO_DATASETS
    initialize_learned_imv_sort_join_thread_args(&rel_r, &rel_s, 
                        rmi_ptr, r_rmi_ptr, s_rmi_ptr,
                        rmi_params, rmi_params, rmi_params,                              
                        SAMPLE_SZ_R, SAMPLE_SZ_S,
                        tmp_training_sample_in, sorted_training_sample_in,
                        r_tmp_training_sample_in, r_sorted_training_sample_in,
                        s_tmp_training_sample_in, s_sorted_training_sample_in,
                        &training_data, &training_data, &training_data,
                        sample_count, sample_count_R, sample_count_S,
                        slopes, intercepts,
                        slopes, intercepts,
                        slopes, intercepts,
                        &barrier, joinresult, args_ptr);
#else
    initialize_learned_imv_sort_join_thread_args(&rel_r, &rel_s, 
                        rmi_ptr, r_rmi_ptr, s_rmi_ptr,
                        r_rmi_params, r_rmi_params, s_rmi_params,                              
                        SAMPLE_SZ_R, SAMPLE_SZ_S,
                        r_tmp_training_sample_in, r_sorted_training_sample_in,
                        r_tmp_training_sample_in, r_sorted_training_sample_in,
                        s_tmp_training_sample_in, s_sorted_training_sample_in,
                        &r_training_data, &r_training_data, &s_training_data,
                        sample_count_R, sample_count_R, sample_count_S,
                        r_slopes, r_intercepts,
                        r_slopes, r_intercepts,
                        s_slopes, s_intercepts,
                        &barrier, joinresult, args_ptr);
#endif


    for(i = 0; i < NUM_THREADS; i++){
        #ifdef DEVELOPMENT_MODE
        int cpu_idx = get_cpu_id_develop(i);
        #else
        int cpu_idx = get_cpu_id_v2(i);
        numa_thread_mark_active(cpu_idx);
        #endif

        CPU_ZERO(&set);
        CPU_SET(cpu_idx, &set);
        pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);

        rv = pthread_create(&tid[i], &attr, learned_imv_sort_join_thread, (void*)&args[i]);
        if (rv){
            printf("[ERROR] return code from pthread_create() is %d\n", rv);
            exit(-1);
        }

    }


    // wait for threads to finish
    for(i = 0; i < NUM_THREADS; i++){
        pthread_join(tid[i], NULL);
        result += args[i].result;
    }        
    joinresult->totalresults = result;
    joinresult->nthreads     = NUM_THREADS;

    printf("join results: %ld \n", result);

    return 0;
}

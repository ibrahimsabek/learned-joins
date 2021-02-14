
/* data generation functionalities for input relations */

#pragma once

#include <stdio.h>              /* perror */
#include <time.h>               /* time() */
#include <random>
#include <type_traits>

#include "configs/base_configs.h"
#include "utils/data_structures.h"
#include "utils/memory.h"
#include "utils/math.h"
#include "utils/barrier.h"
#include "utils/cpu_mapping.h"

#define ZERO_PAYLOAD

using namespace std;

enum synthetic_workload_distr_t
{
    NORMAL,
    UNIFORM,
    EXPONENTIAL,
    LOGNORMAL
};

struct DataDistnParams
{
    double normal_mean = 0;
    double normal_stddev = 1;
    double uniform_a = 0;
    double uniform_b = 1;
    double exponential_lambda = 2;
    double exponential_scale = 1e6;
    double lognormal_mean = 0;
    double lognormal_stddev = 0.25;
    double lognormal_scale = 1e6;
};

template<class KeyType, class PayloadType>
struct create_arg_t {
    Relation<KeyType, PayloadType>          rel;
    int64_t             firstkey;
    int64_t             maxid;
    Relation<KeyType, PayloadType> *        fullrel;
    volatile void *     locks;
    pthread_barrier_t * barrier;
};

template<class Type>
struct create_arg_generic_t {
    Type*          rel;
    int64_t     num_elems;
    int64_t             firstkey;
    int64_t             maxid;
    Type**        fullrel;
    volatile void *     locks;
    pthread_barrier_t * barrier;
};

//based on the ETH implementation 
//create an ETH workload relation of primary keys, where the keys are unique, have no-gaps, and randomly shuffled
template<class KeyType, class PayloadType>
int create_eth_workload_relation_pk(Relation<KeyType, PayloadType> * relation, int64_t num_tuples, int relation_padding);

//based on the ETH implementation 
//create an ETH workload relation of foreign keys, where the keys are non-unique, randomly shuffled, and each key between 0 and maxid exists at least once
template<class KeyType, class PayloadType>
int create_eth_workload_relation_fk(Relation<KeyType, PayloadType> *relation, int64_t num_tuples, const int64_t maxid, int relation_padding);

//create a synthetic workload relation of foreign keys, where the keys are following a certain data distribution
template<class KeyType, class PayloadType>
int create_synthetic_workload_relation_fk(Relation<KeyType, PayloadType> *relation, int64_t num_tuples, synthetic_workload_distr_t data_distn_type, DataDistnParams* data_distn_params, int relation_padding);

template<class KeyType, class PayloadType>
int numa_localize(Tuple<KeyType, PayloadType> * relation, int64_t num_tuples, uint32_t nthreads);

template<class KeyType, class PayloadType>
int numa_localize_varlen(Tuple<KeyType, PayloadType> * relation, int64_t* num_tuples_for_threads, uint32_t nthreads);

template<class Type>
int numa_localize_generic_varlen(Type * relation, int64_t* num_tuples_for_threads, uint32_t nthreads);

//based on the ETH implementation 
static int seeded = 0;
static unsigned int seedValue;

//based on the ETH implementation 
static void check_seed()
{
    if(!seeded) {
        seedValue = time(NULL);
        srand(seedValue);
        seeded = 1;
    }
}

//based on the ETH implementation 
template<class KeyType, class PayloadType>
void knuth_shuffle(Relation<KeyType, PayloadType> * relation)
{
    uint64_t i;
    for (i = relation->num_tuples - 1; i > 0; i--) {
        int64_t  j             = RAND_RANGE(i);
        KeyType tmp            = relation->tuples[i].key;
        relation->tuples[i].key = relation->tuples[j].key;
        relation->tuples[j].key = tmp;

        PayloadType tmp1            = relation->tuples[i].payload;
        relation->tuples[i].payload = relation->tuples[j].payload;
        relation->tuples[j].payload = tmp1;
    }
}

//based on the ETH implementation 
template<class KeyType, class PayloadType>
void random_unique_gen(Relation<KeyType, PayloadType> * rel) 
{
    uint64_t i;

    for (i = 0; i < rel->num_tuples; i++) {
        rel->tuples[i].key = (KeyType)(i+1);
#ifdef ZERO_PAYLOAD
        rel->tuples[i].payload = (PayloadType) 0; 
#else
        rel->tuples[i].payload = (PayloadType)(i+1);
#endif
    }

    /* randomly shuffle elements */
    knuth_shuffle<KeyType, PayloadType>(rel);
}


template<class KeyType, class PayloadType>
int create_eth_workload_relation_pk(Relation<KeyType, PayloadType> *relation, int64_t num_tuples, int relation_padding) 
{
    check_seed();

    relation->num_tuples = num_tuples;
    relation->tuples = (Tuple<KeyType, PayloadType> *) alloc_aligned(num_tuples * sizeof(Tuple<KeyType, PayloadType>) + relation_padding);

    if (!relation->tuples) { 
        perror("out of memory");
        return -1; 
    }
  
    random_unique_gen<KeyType, PayloadType>(relation);

    return 0;
}

template<class KeyType, class PayloadType>
int create_eth_workload_relation_fk(Relation<KeyType, PayloadType> *relation, int64_t num_tuples, const int64_t maxid, int relation_padding)
{
    int32_t i, iters;
    int64_t remainder;
    Relation<KeyType, PayloadType> tmp;

    check_seed();

    relation->num_tuples = num_tuples;
    relation->tuples = (Tuple<KeyType, PayloadType>*) alloc_aligned(relation->num_tuples * sizeof(Tuple<KeyType, PayloadType>) + relation_padding);
      
    if (!relation->tuples) { 
        perror("out of memory");
        return -1; 
    }
  
    /* alternative generation method */
    iters = num_tuples / maxid;
    for(i = 0; i < iters; i++){
        tmp.num_tuples = maxid;
        tmp.tuples = relation->tuples + maxid * i;
        random_unique_gen<KeyType, PayloadType>(&tmp);
    }

    /* if num_tuples is not an exact multiple of maxid */
    remainder = num_tuples % maxid;
    if(remainder > 0) {
        tmp.num_tuples = remainder;
        tmp.tuples = relation->tuples + maxid * iters;
        random_unique_gen<KeyType, PayloadType>(&tmp);
    }

    return 0;
}

template<class KeyType, class PayloadType>
void normal_distr_gen(Relation<KeyType, PayloadType> * relation, double mean = 0, double stddev = 1)
{
    random_device rd;
    mt19937 generator(rd());
    normal_distribution<> distribution(mean, stddev);

    uint64_t i; 

    for (i = 0; i < relation->num_tuples; i++)
    {
        relation->tuples[i].key = (KeyType)(distribution(generator) * relation->num_tuples);
#ifdef ZERO_PAYLOAD
        relation->tuples[i].payload = (PayloadType)0;
#else
        relation->tuples[i].payload = (PayloadType)(relation->tuples[i].key);
#endif
    }
}

template<class KeyType, class PayloadType>
void uniform_distr_gen(Relation<KeyType, PayloadType> * relation, double a = 0, double b = 1)
{
    random_device rd;
    mt19937 generator(rd());
    uniform_real_distribution<> distribution(a, b);

    uint64_t i; 

    for (i = 0; i < relation->num_tuples; i++)
    {
        relation->tuples[i].key = (KeyType)(distribution(generator) * relation->num_tuples);
#ifdef ZERO_PAYLOAD
        relation->tuples[i].payload = (PayloadType)0;
#else
        relation->tuples[i].payload = (PayloadType)(relation->tuples[i].key);
#endif
    }
}

template<class KeyType, class PayloadType>
void exponential_distr_gen(Relation<KeyType, PayloadType> * relation, double lambda = 2, double scale = 1e6)
{
    random_device rd;
    mt19937 generator(rd());
    exponential_distribution<> distribution(lambda);

    uint64_t i; 

    for (i = 0; i < relation->num_tuples; i++)
    {
        relation->tuples[i].key = (KeyType)(distribution(generator) * scale);
#ifdef ZERO_PAYLOAD
        relation->tuples[i].payload = (PayloadType)0;
#else
        relation->tuples[i].payload = (PayloadType)(relation->tuples[i].key);
#endif
    }
}

template<class KeyType, class PayloadType>
void lognormal_distr_gen(Relation<KeyType, PayloadType> * relation, double mean = 0, double stddev = 1, double scale = 1e6)
{
    random_device rd;
    mt19937 generator(rd());
    lognormal_distribution<> distribution(mean, stddev);

    uint64_t i; 

    for (i = 0; i < relation->num_tuples; i++)
    {
        relation->tuples[i].key = (KeyType)(distribution(generator) * scale);
#ifdef ZERO_PAYLOAD
        relation->tuples[i].payload = (PayloadType)0;
#else
        relation->tuples[i].payload = (PayloadType)(relation->tuples[i].key);
#endif
    }
}

template<class KeyType, class PayloadType>
int create_synthetic_workload_relation_fk(Relation<KeyType, PayloadType> *relation, int64_t num_tuples, synthetic_workload_distr_t data_distn_type, DataDistnParams* data_distn_params, int relation_padding) 
{
    relation->num_tuples = num_tuples;
    relation->tuples = (Tuple<KeyType, PayloadType> *) alloc_aligned(num_tuples * sizeof(Tuple<KeyType, PayloadType>) + relation_padding);

    if (!relation->tuples) { 
        perror("out of memory");
        return -1; 
    }

    switch (data_distn_type) {
      case NORMAL:
        normal_distr_gen<KeyType, PayloadType>(relation, data_distn_params->normal_mean, data_distn_params->normal_stddev);
        break;
      case UNIFORM:
        uniform_distr_gen<KeyType, PayloadType>(relation, data_distn_params->uniform_a, data_distn_params->uniform_b);
        break;
      case EXPONENTIAL:
        exponential_distr_gen<KeyType, PayloadType>(relation, data_distn_params->exponential_lambda, data_distn_params->exponential_scale);
        break;
      case LOGNORMAL:
        lognormal_distr_gen<KeyType, PayloadType>(relation, data_distn_params->lognormal_mean, data_distn_params->lognormal_stddev, data_distn_params->lognormal_scale);
        break;
      default:
        perror("undefined data distribution type");
        return -1;
    }

    return 0;
}

template<class KeyType, class PayloadType>
void * numa_localize_thread(void * args) 
{
    create_arg_t<KeyType, PayloadType> * arg = (create_arg_t<KeyType, PayloadType> *) args;
    Relation<KeyType, PayloadType> *   rel = & arg->rel;
    uint64_t i;
    
    for (i = 0; i < rel->num_tuples; i++) {
        rel->tuples[i].key = 0;
    }

    return 0;
}

template<class Type>
void * numa_localize_generic_thread(void * args) 
{
    create_arg_generic_t<Type> * arg = (create_arg_generic_t<Type> *) args;
    Type *   rel = arg->rel;
    uint64_t num_elems = arg->num_elems;
    uint64_t i;
    
    for (i = 0; i < num_elems; i++) {
        rel[i] = 0;
    }

    return 0;
}

template<class KeyType, class PayloadType>
int numa_localize(Tuple<KeyType, PayloadType> * relation, int64_t num_tuples, uint32_t nthreads)
{
    uint32_t i, rv;
    uint64_t offset = 0;

    /* we need aligned allocation of items */
    create_arg_t<KeyType, PayloadType> args[nthreads];
    pthread_t tid[nthreads];
    cpu_set_t set;
    pthread_attr_t attr;

    unsigned int pagesize;
    unsigned int npages;
    unsigned int npages_perthr;
    uint64_t ntuples_perthr;
    uint64_t ntuples_lastthr;

    pagesize        = getpagesize();
    npages          = (num_tuples * sizeof(Tuple<KeyType, PayloadType>)) / pagesize + 1;
    npages_perthr   = npages / nthreads;
    ntuples_perthr  = npages_perthr * (pagesize/sizeof(Tuple<KeyType, PayloadType>));
    ntuples_lastthr = num_tuples - ntuples_perthr * (nthreads-1);

    pthread_attr_init(&attr);

    for( i = 0; i < nthreads; i++ ) {
    #ifdef DEVELOPMENT_MODE
    int cpu_idx = get_cpu_id_develop(i);
    #else
    int cpu_idx = get_cpu_id_v2(i);
    #endif
       
        CPU_ZERO(&set);
        CPU_SET(cpu_idx, &set);
        pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);

        args[i].firstkey       = offset + 1;
        args[i].rel.tuples     = relation + offset;
        args[i].rel.num_tuples = (i == nthreads-1) ? ntuples_lastthr 
                                 : ntuples_perthr;
        offset += ntuples_perthr;

        rv = pthread_create(&tid[i], &attr, numa_localize_thread<KeyType, PayloadType>, (void*)&args[i]);
        if (rv){
            fprintf(stderr, "[ERROR] pthread_create() return code is %d\n", rv);
            exit(-1);
        }
    }

    for(i = 0; i < nthreads; i++){
        pthread_join(tid[i], NULL);
    }

    return 0;

}

template<class KeyType, class PayloadType>
int numa_localize_varlen(Tuple<KeyType, PayloadType> * relation, int64_t* num_tuples_for_threads, uint32_t nthreads)
{
    uint32_t i, rv;
    uint64_t offset = 0;

    /* we need aligned allocation of items */
    create_arg_t<KeyType, PayloadType> args[nthreads];
    pthread_t tid[nthreads];
    cpu_set_t set;
    pthread_attr_t attr;

/*
    unsigned int pagesize;
    unsigned int npages_perthr;
    uint64_t ntuples_perthr [nthreads];

    pagesize        = getpagesize();
    for(i = 0; i < nthreads; i++)
    {
        npages_perthr = (num_tuples_for_threads[i] * sizeof(Tuple<KeyType, PayloadType>)) / pagesize + 1;
        ntuples_perthr[i]  = npages_perthr * (pagesize/sizeof(Tuple<KeyType, PayloadType>));
    }
*/
    pthread_attr_init(&attr);

    for( i = 0; i < nthreads; i++ ) {
        #ifdef DEVELOPMENT_MODE
        int cpu_idx = get_cpu_id_develop(i);
        #else
        int cpu_idx = get_cpu_id_v2(i);
        #endif

        CPU_ZERO(&set);
        CPU_SET(cpu_idx, &set);
        pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);

        args[i].firstkey       = offset + 1;
        args[i].rel.tuples     = relation + offset;
        args[i].rel.num_tuples = num_tuples_for_threads[i];

        offset += num_tuples_for_threads[i];

        rv = pthread_create(&tid[i], &attr, numa_localize_thread<KeyType, PayloadType>, (void*)&args[i]);
        if (rv){
            fprintf(stderr, "[ERROR] pthread_create() return code is %d\n", rv);
            exit(-1);
        }
    }

    for(i = 0; i < nthreads; i++){
        pthread_join(tid[i], NULL);
    }

    return 0;

}


template<class Type>
int numa_localize_generic_varlen(Type * relation, int64_t* num_tuples_for_threads, uint32_t nthreads)
{
    uint32_t i, rv;
    uint64_t offset = 0;

    /* we need aligned allocation of items */
    create_arg_generic_t<Type> args[nthreads];
    pthread_t tid[nthreads];
    cpu_set_t set;
    pthread_attr_t attr;

    pthread_attr_init(&attr);

    for( i = 0; i < nthreads; i++ ) {
        #ifdef DEVELOPMENT_MODE
        int cpu_idx = get_cpu_id_develop(i);
        #else
        int cpu_idx = get_cpu_id_v2(i);
        #endif

        CPU_ZERO(&set);
        CPU_SET(cpu_idx, &set);
        pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);

        args[i].firstkey       = offset + 1;
        args[i].rel     = relation + offset;
        args[i].num_elems = num_tuples_for_threads[i];

        offset += num_tuples_for_threads[i];

        rv = pthread_create(&tid[i], &attr, numa_localize_generic_thread<Type>, (void*)&args[i]);
        if (rv){
            fprintf(stderr, "[ERROR] pthread_create() return code is %d\n", rv);
            exit(-1);
        }
    }

    for(i = 0; i < nthreads; i++){
        pthread_join(tid[i], NULL);
    }

    return 0;

}

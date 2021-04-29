/* IO read/write functionalities for input relations based on the ETH implementation */

#pragma once

#include <stdio.h>              /* perror */
#include <iostream>
#include <fstream>
#include <vector>

#include "utils/data_structures.h"
#include "utils/memory.h"


template<class KeyType, class PayloadType>
struct write_read_arg_t {
    Relation<KeyType, PayloadType>          rel;
    int                 thread_id;
    const char * folder_path;    
    const char * filename;
    const char * file_extension;
    Relation<KeyType, PayloadType> *        fullrel;
};


//based on the ETH implementation 
// persisting a relation to the disk
template<class KeyType, class PayloadType>
void write_relation(Relation<KeyType, PayloadType>* rel, const char * filename);

template<class KeyType, class PayloadType>
void write_relation_threaded(Relation<KeyType, PayloadType>* rel, int nthreads, const char * folder_path, const char * filename, const char * file_extension);

//based on the ETH implementation 
// reading a persisted relation from the disk
template<class KeyType, class PayloadType>
int load_relation(Relation<KeyType, PayloadType>* relation, const char * filename, uint64_t num_tuples);

template<class KeyType, class PayloadType>
int load_relation_threaded(Relation<KeyType, PayloadType>* relation, int nthreads, const char * folder_path, const char * filename, const char * file_extension, uint64_t num_tuples);

/**
 * Free memory allocated for only tuples.
 */
template<class KeyType, class PayloadType>
void delete_relation(Relation<KeyType, PayloadType> * rel);

template<typename T>
std::vector<T> load_binary_vector_data(std::string& filename);

template<class KeyType, class PayloadType>
void write_relation(Relation<KeyType, PayloadType>* rel, const char * filename)
{
    FILE * fp = fopen(filename, "w");
    uint64_t i;

    fprintf(fp, "#KEY, VAL\n");

    if (std::is_same<KeyType, int>::value)
    {
        if (std::is_same<PayloadType, int>::value)
            for (i = 0; i < rel->num_tuples; i++)
                fprintf(fp, "%d %d\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, unsigned int>::value)
            for (i = 0; i < rel->num_tuples; i++)                
                fprintf(fp, "%d %u\n", rel->tuples[i].key, rel->tuples[i].payload);            
        else if(std::is_same<PayloadType, long long int>::value)
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%d %ld\n", rel->tuples[i].key, rel->tuples[i].payload);            
        else if(std::is_same<PayloadType, unsigned long long int>::value)
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%d %lu\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, double>::value)
            for (i = 0; i < rel->num_tuples; i++)
                fprintf(fp, "%d %lf\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, long double>::value)
            for (i = 0; i < rel->num_tuples; i++)              
                fprintf(fp, "%d %Lf\n", rel->tuples[i].key, rel->tuples[i].payload);            
    }
    else if(std::is_same<KeyType, unsigned int>::value)
    {
        if (std::is_same<PayloadType, int>::value)
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%u %d\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, unsigned int>::value)
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%u %u\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, long long int>::value)
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%u %ld\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, unsigned long long int>::value)
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%u %lu\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, double>::value)
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%u %lf\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, long double>::value)
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%u %Lf\n", rel->tuples[i].key, rel->tuples[i].payload);            
    }
    else if(std::is_same<KeyType, long long int>::value)
    {
        if (std::is_same<PayloadType, int>::value)
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%ld %d\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, unsigned int>::value)
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%ld %u\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, long long int>::value)
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%ld %ld\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, unsigned long long int>::value)
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%ld %lu\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, double>::value)
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%ld %lf\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, long double>::value)
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%ld %Lf\n", rel->tuples[i].key, rel->tuples[i].payload);
    }
    else if(std::is_same<KeyType, unsigned long long int>::value)
    {
        if (std::is_same<PayloadType, int>::value)
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%lu %d\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, unsigned int>::value)
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%lu %u\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, long long int>::value)
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%lu %ld\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, unsigned long long int>::value)
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%lu %lu\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, double>::value)
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%lu %lf\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, long double>::value)  
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%lu %Lf\n", rel->tuples[i].key, rel->tuples[i].payload);                      
    }
    else if(std::is_same<KeyType, double>::value)
    {
        if (std::is_same<PayloadType, int>::value)
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%lf %d\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, unsigned int>::value)
            for (i = 0; i < rel->num_tuples; i++)
                fprintf(fp, "%lf %u\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, long long int>::value)
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%lf %ld\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, unsigned long long int>::value)
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%lf %lu\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, double>::value)
            for (i = 0; i < rel->num_tuples; i++)
                fprintf(fp, "%lf %lf\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, long double>::value)            
            for (i = 0; i < rel->num_tuples; i++)
                fprintf(fp, "%lf %Lf\n", rel->tuples[i].key, rel->tuples[i].payload);
    }
    else if(std::is_same<KeyType, long double>::value)
    {
        if (std::is_same<PayloadType, int>::value)
            for (i = 0; i < rel->num_tuples; i++)
                fprintf(fp, "%Lf %d\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, unsigned int>::value)
            for (i = 0; i < rel->num_tuples; i++)
                fprintf(fp, "%Lf %u\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, long long int>::value)
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%Lf %ld\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, unsigned long long int>::value)
            for (i = 0; i < rel->num_tuples; i++)            
                fprintf(fp, "%Lf %lu\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, double>::value)
            for (i = 0; i < rel->num_tuples; i++)
                fprintf(fp, "%Lf %lf\n", rel->tuples[i].key, rel->tuples[i].payload);
        else if(std::is_same<PayloadType, long double>::value)           
            for (i = 0; i < rel->num_tuples; i++)
                fprintf(fp, "%Lf %Lf\n", rel->tuples[i].key, rel->tuples[i].payload);        
    }

    fclose(fp);
}

template<class KeyType, class PayloadType>
void * write_relation_thread(void * args) 
{
    write_read_arg_t<KeyType, PayloadType> * arg = (write_read_arg_t<KeyType, PayloadType> *) args;

    stringstream full_filename;
    
    full_filename << arg->folder_path;
    full_filename << arg->filename;
    full_filename << "_";
    full_filename << arg->thread_id;
    full_filename << arg->file_extension;

    write_relation(&(arg->rel), full_filename.str().c_str());

    return 0;
}

template<class KeyType, class PayloadType>
void write_relation_threaded(Relation<KeyType, PayloadType>* rel, int nthreads, const char * folder_path, const char * filename, const char * file_extension)
{    
    uint32_t i, rv;
    uint64_t offset = 0;

    write_read_arg_t<KeyType, PayloadType> args[nthreads];
    pthread_t tid[nthreads];
    cpu_set_t set;
    pthread_attr_t attr;

    uint64_t ntuples_perthr;
    uint64_t ntuples_lastthr;

    ntuples_perthr  = rel->num_tuples / nthreads;
    ntuples_lastthr = rel->num_tuples - ntuples_perthr * (nthreads-1);

    pthread_attr_init(&attr);

    for(i = 0; i < nthreads; i++ ) 
    {
        #ifdef DEVELOPMENT_MODE
        int cpu_idx = get_cpu_id_develop(i);
        #else
        int cpu_idx = get_cpu_id_v2(i);
        #endif
    
        CPU_ZERO(&set);
        CPU_SET(cpu_idx, &set);
        pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);

        args[i].folder_path = folder_path;
        args[i].filename = filename;
        args[i].file_extension = file_extension;
        args[i].thread_id = i;
        args[i].rel.tuples     = rel->tuples + offset;
        args[i].rel.num_tuples = (i == nthreads-1) ? ntuples_lastthr 
                                 : ntuples_perthr;
        offset += ntuples_perthr;

        rv = pthread_create(&tid[i], &attr, write_relation_thread<KeyType, PayloadType>, (void*)&args[i]);
        if (rv){
            fprintf(stderr, "[ERROR] pthread_create() return code is %d\n", rv);
            exit(-1);
        }
    }

    for(i = 0; i < nthreads; i++){
        pthread_join(tid[i], NULL);
    }
}


// based on the ETH implementation
template<class KeyType, class PayloadType>
void read_relation(Relation<KeyType, PayloadType> * rel, const char * filename){

    FILE * fp = fopen(filename, "r");

    /* skip the header line */
    char c;
    do{
        c = fgetc(fp);
    } while (c != '\n');

    /* search for a whitespace for "key payload" format */
    int fmtspace = 0;
    //int fmtcomma = 0;
    do{
        c = fgetc(fp);
        if(c == ' '){
            fmtspace = 1;
            break;
        }
    //    if(c == ','){
    //        fmtcomma = 1;
    //        break;
    //    }
    } while (c != '\n');

    /* rewind back to the beginning and start parsing again */
    rewind(fp);
    /* skip the header line */
    do{
        c = fgetc(fp);
    } while (c != '\n');

    uint64_t ntuples = rel->num_tuples;
    KeyType key;
    PayloadType payload = 0;
    int warn = 1;
    for(uint64_t i = 0; i < ntuples; i++){
        if(fmtspace){
            if (std::is_same<KeyType, int>::value)
            {
                if (std::is_same<PayloadType, int>::value)
                        fscanf(fp, "%d %d\n", &key, &payload);
                else if(std::is_same<PayloadType, unsigned int>::value)
                        fscanf(fp, "%d %u\n", &key, &payload);            
                else if(std::is_same<PayloadType, long long int>::value)
                        fscanf(fp, "%d %ld\n", &key, &payload);            
                else if(std::is_same<PayloadType, unsigned long long int>::value)
                        fscanf(fp, "%d %lu\n", &key, &payload);
                else if(std::is_same<PayloadType, double>::value)
                        fscanf(fp, "%d %lf\n", &key, &payload);
                else if(std::is_same<PayloadType, long double>::value)
                        fscanf(fp, "%d %Lf\n", &key, &payload);            
            }
            else if(std::is_same<KeyType, unsigned int>::value)
            {
                if (std::is_same<PayloadType, int>::value)
                        fscanf(fp, "%u %d\n", &key, &payload);
                else if(std::is_same<PayloadType, unsigned int>::value)
                        fscanf(fp, "%u %u\n", &key, &payload);
                else if(std::is_same<PayloadType, long long int>::value)
                        fscanf(fp, "%u %ld\n", &key, &payload);
                else if(std::is_same<PayloadType, unsigned long long int>::value)
                        fscanf(fp, "%u %lu\n", &key, &payload);
                else if(std::is_same<PayloadType, double>::value)
                        fscanf(fp, "%u %lf\n", &key, &payload);
                else if(std::is_same<PayloadType, long double>::value)
                        fscanf(fp, "%u %Lf\n", &key, &payload);            
            }
            else if(std::is_same<KeyType, long long int>::value)
            {
                if (std::is_same<PayloadType, int>::value)
                        fscanf(fp, "%ld %d\n", &key, &payload);
                else if(std::is_same<PayloadType, unsigned int>::value)
                        fscanf(fp, "%ld %u\n", &key, &payload);
                else if(std::is_same<PayloadType, long long int>::value)
                        fscanf(fp, "%ld %ld\n", &key, &payload);
                else if(std::is_same<PayloadType, unsigned long long int>::value)
                        fscanf(fp, "%ld %lu\n", &key, &payload);
                else if(std::is_same<PayloadType, double>::value)
                        fscanf(fp, "%ld %lf\n", &key, &payload);
                else if(std::is_same<PayloadType, long double>::value)
                        fscanf(fp, "%ld %Lf\n", &key, &payload);
            }
            else if(std::is_same<KeyType, unsigned long long int>::value)
            {
                if (std::is_same<PayloadType, int>::value)
                        fscanf(fp, "%lu %d\n", &key, &payload);
                else if(std::is_same<PayloadType, unsigned int>::value)
                        fscanf(fp, "%lu %u\n", &key, &payload);
                else if(std::is_same<PayloadType, long long int>::value)
                        fscanf(fp, "%lu %ld\n", &key, &payload);
                else if(std::is_same<PayloadType, unsigned long long int>::value)
                        fscanf(fp, "%lu %lu\n", &key, &payload);
                else if(std::is_same<PayloadType, double>::value)
                        fscanf(fp, "%lu %lf\n", &key, &payload);
                else if(std::is_same<PayloadType, long double>::value)  
                        fscanf(fp, "%lu %Lf\n", &key, &payload);                      
            }
            else if(std::is_same<KeyType, double>::value)
            {
                if (std::is_same<PayloadType, int>::value)
                        fscanf(fp, "%lf %d\n", &key, &payload);
                else if(std::is_same<PayloadType, unsigned int>::value)
                        fscanf(fp, "%lf %u\n", &key, &payload);
                else if(std::is_same<PayloadType, long long int>::value)
                        fscanf(fp, "%lf %ld\n", &key, &payload);
                else if(std::is_same<PayloadType, unsigned long long int>::value)
                        fscanf(fp, "%lf %lu\n", &key, &payload);
                else if(std::is_same<PayloadType, double>::value)
                        fscanf(fp, "%lf %lf\n", &key, &payload);
                else if(std::is_same<PayloadType, long double>::value)            
                        fscanf(fp, "%lf %Lf\n", &key, &payload);
            }
            else if(std::is_same<KeyType, long double>::value)
            {
                if (std::is_same<PayloadType, int>::value)
                        fscanf(fp, "%Lf %d\n", &key, &payload);
                else if(std::is_same<PayloadType, unsigned int>::value)
                        fscanf(fp, "%Lf %u\n", &key, &payload);
                else if(std::is_same<PayloadType, long long int>::value)
                        fscanf(fp, "%Lf %ld\n", &key, &payload);
                else if(std::is_same<PayloadType, unsigned long long int>::value)
                        fscanf(fp, "%Lf %lu\n", &key, &payload);
                else if(std::is_same<PayloadType, double>::value)
                        fscanf(fp, "%Lf %lf\n", &key, &payload);
                else if(std::is_same<PayloadType, long double>::value)           
                        fscanf(fp, "%Lf %Lf\n", &key, &payload);        
            }
        }
        else
        {
            perror("Incorrect format");
        }

        if(warn && key < 0){
            warn = 0;
            printf("[WARN ] key=%d, payload=%d\n", key, payload);
        }
        rel->tuples[i].key = key;
        rel->tuples[i].payload = payload;
    }

    fclose(fp);
}

template<class KeyType, class PayloadType>
int load_relation(Relation<KeyType, PayloadType>* relation, const char * filename, uint64_t num_tuples)
{
    relation->num_tuples = num_tuples;

    /* we need aligned allocation of items */
    relation->tuples = (Tuple<KeyType, PayloadType> *) alloc_aligned(num_tuples * sizeof(Tuple<KeyType, PayloadType>));

    if (!relation->tuples) { 
        perror("out of memory");
        return -1; 
    }

    /* load from the given input file */
    read_relation(relation, filename);

    return 0;    
}

template<class KeyType, class PayloadType>
void * read_relation_thread(void * args) 
{
    write_read_arg_t<KeyType, PayloadType> * arg = (write_read_arg_t<KeyType, PayloadType> *) args;

    stringstream full_filename;
    
    full_filename << arg->folder_path;
    full_filename << arg->filename;
    full_filename << "_";
    full_filename << arg->thread_id;
    full_filename << arg->file_extension;

    read_relation(&(arg->rel), full_filename.str().c_str());

    return 0;
}

template<class KeyType, class PayloadType>
void read_relation_threaded(Relation<KeyType, PayloadType>* rel, int nthreads, const char * folder_path, const char * filename, const char * file_extension)
{    
    uint32_t i, rv;
    uint64_t offset = 0;

    write_read_arg_t<KeyType, PayloadType> args[nthreads];
    pthread_t tid[nthreads];
    cpu_set_t set;
    pthread_attr_t attr;

    uint64_t ntuples_perthr;
    uint64_t ntuples_lastthr;

    ntuples_perthr  = rel->num_tuples / nthreads;
    ntuples_lastthr = rel->num_tuples - ntuples_perthr * (nthreads-1);

    pthread_attr_init(&attr);

    for(i = 0; i < nthreads; i++ ) 
    {
        #ifdef DEVELOPMENT_MODE
        int cpu_idx = get_cpu_id_develop(i);
        #else
        int cpu_idx = get_cpu_id_v2(i);
        #endif
    
        CPU_ZERO(&set);
        CPU_SET(cpu_idx, &set);
        pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &set);

        args[i].folder_path = folder_path;
        args[i].filename = filename;
        args[i].file_extension = file_extension;
        args[i].thread_id = i;
        args[i].rel.tuples     = rel->tuples + offset;
        args[i].rel.num_tuples = (i == nthreads-1) ? ntuples_lastthr 
                                 : ntuples_perthr;
        offset += ntuples_perthr;

        rv = pthread_create(&tid[i], &attr, read_relation_thread<KeyType, PayloadType>, (void*)&args[i]);
        if (rv){
            fprintf(stderr, "[ERROR] pthread_create() return code is %d\n", rv);
            exit(-1);
        }
    }

    for(i = 0; i < nthreads; i++){
        pthread_join(tid[i], NULL);
    }
}

template<class KeyType, class PayloadType>
int load_relation_threaded(Relation<KeyType, PayloadType>* relation, int nthreads, const char * folder_path, const char * filename, const char * file_extension, uint64_t num_tuples)
{
    relation->num_tuples = num_tuples;

    /* we need aligned allocation of items */
    relation->tuples = (Tuple<KeyType, PayloadType> *) alloc_aligned(num_tuples * sizeof(Tuple<KeyType, PayloadType>));

    if (!relation->tuples) { 
        perror("out of memory");
        return -1; 
    }

    read_relation_threaded(relation, nthreads, folder_path, filename, file_extension);

    return 0;
}


template<class KeyType, class PayloadType>
void delete_relation(Relation<KeyType, PayloadType> * rel)
{
    /* clean up */
    free(rel->tuples);
}

// Adapted implementation from the SOSD benchmark
// Loads values from binary file into vector.
template<typename T>
std::vector<T> load_binary_vector_data(std::string& filename) {
    std::vector<T> data;

    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open()) {
        std::cerr << "unable to open " << filename << std::endl;
        exit(EXIT_FAILURE);
    }
    // Read size.
    uint64_t size;
    in.read(reinterpret_cast<char*>(&size), sizeof(uint64_t));
    data.resize(size);
    // Read values.
    in.read(reinterpret_cast<char*>(data.data()), size*sizeof(T));
    in.close();

    return data;
}
#!/bin/sh

process_radix_join()
{
    threads=(16) #(2 4 8 16 32 64)
    num_radix_bits=(12 14) #(9 10 11 12 14 16)
    num_passes=(1) #(1 2)
    use_murmur3_hash_for_radix_join=(1) #(0 1)
    use_vectorized_murmur3_hash_for_radix_join=(1) #(0 1)

    dataset_folder_path=/spinning/sabek/learned_join_datasets/

    r_datasets=$1
    r_datasets_sizes=$2
    s_datasets=$3
    s_datasets_sizes=$4

    output_folder_path=$5

    mkdir $output_folder_path

    for ds in ${!r_datasets[@]}
    do
        curr_r_dataset='"'$dataset_folder_path${r_datasets[$ds]}.txt'"'
        curr_r_dataset_size=${r_datasets_sizes[$ds]}
        curr_s_dataset='"'$dataset_folder_path${s_datasets[$ds]}.txt'"'
        curr_s_dataset_size=${s_datasets_sizes[$ds]}

        echo 'Joining '$curr_r_dataset' '$curr_r_dataset_size' '$curr_s_dataset' '$curr_s_dataset_size'...'
        
        for th in ${!threads[@]}
        do
            curr_threads=${threads[$th]}
            for rb in ${!num_radix_bits[@]}
            do
                curr_num_radix_bits=${num_radix_bits[$rb]}
                for np in ${!num_passes[@]}
                do
                    curr_num_passes=${num_passes[$np]}

                    curr_output_file=$output_folder_path'radix_join_tuning_v'$((ds))'_th_'$curr_threads'_rb_'$curr_num_radix_bits'_np_'$curr_num_passes'.csv'

                    sh $(dirname "$0")/base_configs_maker.sh -USE_LEARNED_SORT_MERGE_JOIN_GOOGLE_BENCHMARK 0 \
                                                            -USE_ETH_SORT_MERGE_JOIN_GOOGLE_BENCHMARK 0 \
                                                            -USE_ETH_NON_PARTITION_JOIN_GOOGLE_BENCHMARK 0 \
                                                            -USE_ETH_RADIX_JOIN_GOOGLE_BENCHMARK 1 \
                                                            -NUM_THREADS_FOR_EVALUATION $curr_threads \
                                                            -RELATION_R_PATH $curr_r_dataset \
                                                            -RELATION_R_NUM_TUPLES $curr_r_dataset_size \
                                                            -RELATION_S_PATH $curr_s_dataset \
                                                            -RELATION_S_NUM_TUPLES $curr_s_dataset_size \
                                                            -SKEW_HANDLING 1
                                                      
                    sh $(dirname "$0")/eth_configs_maker.sh -NUM_RADIX_BITS $curr_num_radix_bits \
                                                            -NUM_PASSES $curr_num_passes \
                                                            -SKEWNESS_THRESHOLD_MULTIPLIER 64 \
                                                            -USE_MURMUR3_HASH_FOR_RADIX_JOIN $use_murmur3_hash_for_radix_join \
                                                            -USE_VECTORIZED_MURMUR3_HASH_FOR_RADIX_JOIN $USE_VECTORIZED_MURMUR3_HASH_FOR_RADIX_JOIN

                    cd $(dirname "$0")/../../
                    sh makefile_am_maker.sh && aclocal && autoconf && automake --add-missing && ./configure --prefix=/spinning/sabek/learned_join_binaries --enable-syncstats --enable-swwc-part > /dev/null && make clean -s && make -s && make install -s > /dev/null

                    ./benchmarkhashjoins --mode=1 --benchmark_out=$curr_output_file --benchmark_out_format=csv
                    
                    echo "$(sed 1,10d $curr_output_file)" > $curr_output_file #remove the header of google benchmark output
                done
            done   
        done
    done
}

#unique datasets
#r_unique_datasets=(r_UNIQUE_v1_uint32_uint32_16000000) #(r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v4_uint32_uint32_384000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v6_uint32_uint32_896000000 r_UNIQUE_v7_uint32_uint32_1152000000 r_UNIQUE_v8_uint32_uint32_1664000000 r_UNIQUE_v9_uint32_uint32_1920000000)
#s_unique_datasets=(s_UNIQUE_v1_uint32_uint32_16000000) #(s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v4_uint32_uint32_384000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v6_uint32_uint32_896000000 s_UNIQUE_v7_uint32_uint32_1152000000 s_UNIQUE_v8_uint32_uint32_1664000000 s_UNIQUE_v9_uint32_uint32_1920000000)
#r_unique_datasets_sizes=(16E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
#s_unique_datasets_sizes=(16E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
#process_radix_join $r_unique_datasets $r_unique_datasets_sizes $s_unique_datasets $s_unique_datasets_sizes /spinning/sabek/learned_join_results/radix_join_tuning_unique/


#lognormal datasets
r_lognormal_datasets=(r_LOGNORMAL_v1_uint32_uint32_segma_1_16000000) #(r_LOGNORMAL_v1_uint32_uint32_segma_1_16000000 r_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 r_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 r_LOGNORMAL_v4_uint32_uint32_segma_1_384000000 r_LOGNORMAL_v5_uint32_uint32_segma_1_640000000 r_LOGNORMAL_v6_uint32_uint32_segma_1_896000000 r_LOGNORMAL_v7_uint32_uint32_segma_1_1152000000 r_LOGNORMAL_v8_uint32_uint32_segma_1_1664000000 r_LOGNORMAL_v9_uint32_uint32_segma_1_1920000000)
s_lognormal_datasets=(s_LOGNORMAL_v1_uint32_uint32_segma_1_16000000) #(s_LOGNORMAL_v1_uint32_uint32_segma_1_16000000 s_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 s_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 s_LOGNORMAL_v4_uint32_uint32_segma_1_384000000 s_LOGNORMAL_v5_uint32_uint32_segma_1_640000000 s_LOGNORMAL_v6_uint32_uint32_segma_1_896000000 s_LOGNORMAL_v7_uint32_uint32_segma_1_1152000000 s_LOGNORMAL_v8_uint32_uint32_segma_1_1664000000 s_LOGNORMAL_v9_uint32_uint32_segma_1_1920000000)
r_lognormal_datasets_sizes=(16E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_lognormal_datasets_sizes=(16E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
process_radix_join $r_lognormal_datasets $r_lognormal_datasets_sizes $s_lognormal_datasets $s_lognormal_datasets_sizes /spinning/sabek/learned_join_results/radix_join_tuning_lognormal/

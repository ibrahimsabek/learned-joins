# A script to run the experiments for MWAY sort merge join algorithm

#!/bin/sh

process_sort_merge_join()
{
    threads=(32) #(2 4 8 16 32 64)
    fanout_per_thread=(128) #(4 16 128 1024 2048 4096)
    use_learned_sort=(0) #(0 1)
    use_avxsort_as_std_sort=(0) #(0 1)

    #dataset_folder_path=/spinning/sabek/learned_join_datasets/
    #dataset_folder_path=/spinning/sabek/learned_join_datasets_sosd/
    dataset_folder_path=/spinning/sabek/learned_join_datasets_tpch/

    #r_datasets=$1
    #r_datasets_sizes=$2
    #r_datasets_file_num_partitions=$3
    r_datasets_file_extension='"'.txt'"'
    #s_datasets=$4
    #s_datasets_sizes=$5
    #s_datasets_file_num_partitions=$6
    s_datasets_file_extension='"'.txt'"'

    output_folder_path=$7
    mkdir $output_folder_path

    run_nums=$8

    load_relations_for_evaluation=$9
    persist_relations_for_evaluation=${10}


    for ds in ${!r_datasets[@]}
    do
        curr_r_dataset='"'$dataset_folder_path${r_datasets[$ds]}.txt'"'
        curr_r_dataset_size=${r_datasets_sizes[$ds]}
        curr_r_dataset_file_num_partitions=${r_datasets_file_num_partitions[$ds]}
        curr_s_dataset='"'$dataset_folder_path${s_datasets[$ds]}.txt'"'
        curr_s_dataset_size=${s_datasets_sizes[$ds]}
        curr_s_dataset_file_num_partitions=${s_datasets_file_num_partitions[$ds]}

        echo 'Joining '$curr_r_dataset' '$curr_r_dataset_size' '$curr_s_dataset' '$curr_s_dataset_size'... \n'
        
        for th in ${!threads[@]}
        do
            curr_threads=${threads[$th]}
            for ft in ${!fanout_per_thread[@]}
            do
                curr_fanout_per_thread=${fanout_per_thread[$ft]}
                for ls in ${!use_learned_sort[@]}
                do
                    curr_use_learned_sort=${use_learned_sort[$ls]}

                    for lsavx in ${!use_avxsort_as_std_sort[@]}
                    do
                        curr_use_avxsort_as_std_sort=${use_avxsort_as_std_sort[$lsavx]}

                        curr_output_file=$output_folder_path'sj_with_eth_tuning_'$curr_r_dataset_size'_'$curr_s_dataset_size'_th_'$curr_threads'_ls_'$curr_use_learned_sort'_ft_'$curr_fanout_per_thread'_lsavx_'$curr_use_avxsort_as_std_sort'.csv'

                        sh $(dirname "$0")/base_configs_maker.sh -USE_LEARNED_SORT_FOR_SORT_MERGE 0 \
                                                                -NUM_THREADS_FOR_EVALUATION $curr_threads \
                                                                -RELATION_R_PATH $curr_r_dataset \
                                                                -RELATION_R_FOLDER_PATH '"'$dataset_folder_path'"' \
                                                                -RELATION_R_FILE_NAME '"'${r_datasets[$ds]}'"' \
                                                                -RELATION_R_FILE_EXTENSION ${r_datasets_file_extension} \
                                                                -RELATION_R_NUM_TUPLES $curr_r_dataset_size \
                                                                -RELATION_R_FILE_NUM_PARTITIONS $curr_r_dataset_file_num_partitions \
                                                                -RELATION_S_PATH $curr_s_dataset \
                                                                -RELATION_S_FOLDER_PATH '"'$dataset_folder_path'"' \
                                                                -RELATION_S_FILE_NAME '"'${s_datasets[$ds]}'"' \
                                                                -RELATION_S_FILE_EXTENSION ${s_datasets_file_extension} \
                                                                -RELATION_S_NUM_TUPLES ${curr_s_dataset_size} \
                                                                -RELATION_S_FILE_NUM_PARTITIONS ${curr_s_dataset_file_num_partitions} \
                                                                -BENCHMARK_RESULTS_PATH '"'${curr_output_file}'"' \
                                                                -RUN_NUMS ${run_nums} -LOAD_RELATIONS_FOR_EVALUATION ${load_relations_for_evaluation} \
                                                                -PERSIST_RELATIONS_FOR_EVALUATION ${persist_relations_for_evaluation} \
                                                                -ETH_SORT_MERGE_PARTITION_FANOUT_PER_THREAD $curr_fanout_per_thread \
                                                                -USE_LEARNED_SORT $curr_use_learned_sort \
                                                                -USE_AVXSORT_AS_STD_SORT $curr_use_avxsort_as_std_sort \
                                                                -USE_AVXSORT_FOR_SORTING_MINOR_BCKTS 0 \
                                                                -ETH_SORT_MERGE_IS_SCALAR_MERGE 0 #\
                                                                #-CUSTOM_CPU_MAPPING '"'../../include/configs/cpu-mapping_berners_lee.txt'"' \
                                                                #-CUSTOM_CPU_MAPPING_V2 '"'../../include/configs/cpu-mapping-v2_berners_lee.txt'"'

                        cmake -DCMAKE_BUILD_TYPE=Release -DVECTORWISE_BRANCHING=on $(dirname "$0")/../.. > /dev/null

                        cd $(dirname "$0")/../../build/release

                        make > /dev/null

                        ./non_learned_non_imv_sort_join_runner

                        cd ../../scripts/evaluation/
                    done    
                done
            done   
        done
    done
}



run_nums=1 #10
load_relations_for_evaluation=1
persist_relations_for_evaluation=0

#VIP NOTE: number of partitions should be multiples of the number of threads
#unique datasets
################

r_datasets=(r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v5_uint32_uint32_640000000) #(r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v4_uint32_uint32_384000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v6_uint32_uint32_896000000 r_UNIQUE_v7_uint32_uint32_1152000000 r_UNIQUE_v8_uint32_uint32_1664000000 r_UNIQUE_v9_uint32_uint32_1920000000)
s_datasets=(s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000) #(s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v4_uint32_uint32_384000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v6_uint32_uint32_896000000 s_UNIQUE_v7_uint32_uint32_1152000000 s_UNIQUE_v8_uint32_uint32_1664000000 s_UNIQUE_v9_uint32_uint32_1920000000)
r_datasets_sizes=(16E6 16E6 16E6 32E6 32E6 32E6 128E6 128E6 128E6 640E6 640E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(32E6 128E6 640E6 16E6 128E6 640E6 16E6 32E6 640E6 16E6 32E6 128E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32 32 32 32 32 32 32 32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32 32 32 32 32 32 32 32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)

#output_folder_path=/spinning/sabek/learned_join_results/sj_with_eth_unique/
#process_sort_merge_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation

r_datasets=(r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v5_uint32_uint32_640000000) #(r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v4_uint32_uint32_384000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v6_uint32_uint32_896000000 r_UNIQUE_v7_uint32_uint32_1152000000 r_UNIQUE_v8_uint32_uint32_1664000000 r_UNIQUE_v9_uint32_uint32_1920000000)
s_datasets=(s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v5_uint32_uint32_640000000) #(s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v4_uint32_uint32_384000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v6_uint32_uint32_896000000 s_UNIQUE_v7_uint32_uint32_1152000000 s_UNIQUE_v8_uint32_uint32_1664000000 s_UNIQUE_v9_uint32_uint32_1920000000)
r_datasets_sizes=(16E6 32E6 128E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(16E6 32E6 128E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32 32 32 32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32 32 32 32) #(64 64 64 64 64 64 64 64 64)

r_datasets=(r_UNIQUE_v5_uint32_uint32_640000000) #(r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v4_uint32_uint32_384000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v6_uint32_uint32_896000000 r_UNIQUE_v7_uint32_uint32_1152000000 r_UNIQUE_v8_uint32_uint32_1664000000 r_UNIQUE_v9_uint32_uint32_1920000000)
s_datasets=(s_UNIQUE_v5_uint32_uint32_640000000) #(s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v4_uint32_uint32_384000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v6_uint32_uint32_896000000 s_UNIQUE_v7_uint32_uint32_1152000000 s_UNIQUE_v8_uint32_uint32_1664000000 s_UNIQUE_v9_uint32_uint32_1920000000)
r_datasets_sizes=(640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32 32 32 32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32 32 32 32) #(64 64 64 64 64 64 64 64 64)

output_folder_path=/spinning/sabek/learned_join_results/sj_with_eth_unique/
process_sort_merge_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation


#lognormal datasets
###################

r_datasets=(r_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 r_LOGNORMAL_v2_uint32_uint32_segma_1_32000000  r_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 r_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 r_LOGNORMAL_v5_uint32_uint32_segma_1_640000000 r_LOGNORMAL_v5_uint32_uint32_segma_1_640000000) #(r_LOGNORMAL_v1_uint32_uint32_segma_1_16000000 r_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 r_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 r_LOGNORMAL_v4_uint32_uint32_segma_1_384000000 r_LOGNORMAL_v5_uint32_uint32_segma_1_640000000 r_LOGNORMAL_v6_uint32_uint32_segma_1_896000000 r_LOGNORMAL_v7_uint32_uint32_segma_1_1152000000 r_LOGNORMAL_v8_uint32_uint32_segma_1_1664000000 r_LOGNORMAL_v9_uint32_uint32_segma_1_1920000000)
s_datasets=(s_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 s_LOGNORMAL_v5_uint32_uint32_segma_1_640000000 s_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 s_LOGNORMAL_v5_uint32_uint32_segma_1_640000000 s_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 s_LOGNORMAL_v3_uint32_uint32_segma_1_128000000) #(s_LOGNORMAL_v1_uint32_uint32_segma_1_16000000 s_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 s_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 s_LOGNORMAL_v4_uint32_uint32_segma_1_384000000 s_LOGNORMAL_v5_uint32_uint32_segma_1_640000000 s_LOGNORMAL_v6_uint32_uint32_segma_1_896000000 s_LOGNORMAL_v7_uint32_uint32_segma_1_1152000000 s_LOGNORMAL_v8_uint32_uint32_segma_1_1664000000 s_LOGNORMAL_v9_uint32_uint32_segma_1_1920000000)
r_datasets_sizes=(32E6 32E6 128E6 128E6 640E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(128E6 640E6 32E6 640E6 32E6 128E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32 32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32 32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)

#output_folder_path=/spinning/sabek/learned_join_results/sj_with_eth_lognormal/
#process_sort_merge_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation


r_datasets=(r_LOGNORMAL_v1_uint32_uint32_segma_1_16000000 r_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 r_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 r_LOGNORMAL_v5_uint32_uint32_segma_1_640000000) #(r_LOGNORMAL_v1_uint32_uint32_segma_1_16000000 r_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 r_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 r_LOGNORMAL_v4_uint32_uint32_segma_1_384000000 r_LOGNORMAL_v5_uint32_uint32_segma_1_640000000 r_LOGNORMAL_v6_uint32_uint32_segma_1_896000000 r_LOGNORMAL_v7_uint32_uint32_segma_1_1152000000 r_LOGNORMAL_v8_uint32_uint32_segma_1_1664000000 r_LOGNORMAL_v9_uint32_uint32_segma_1_1920000000)
s_datasets=(s_LOGNORMAL_v1_uint32_uint32_segma_1_16000000 s_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 s_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 s_LOGNORMAL_v5_uint32_uint32_segma_1_640000000) #(s_LOGNORMAL_v1_uint32_uint32_segma_1_16000000 s_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 s_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 s_LOGNORMAL_v4_uint32_uint32_segma_1_384000000 s_LOGNORMAL_v5_uint32_uint32_segma_1_640000000 s_LOGNORMAL_v6_uint32_uint32_segma_1_896000000 s_LOGNORMAL_v7_uint32_uint32_segma_1_1152000000 s_LOGNORMAL_v8_uint32_uint32_segma_1_1664000000 s_LOGNORMAL_v9_uint32_uint32_segma_1_1920000000)
r_datasets_sizes=(16E6 32E6 128E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(16E6 32E6 128E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32 32 32 32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32 32 32 32) #(64 64 64 64 64 64 64 64 64)

#output_folder_path=/spinning/sabek/learned_join_results/sj_with_eth_lognormal/
#process_sort_merge_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation

#seq_hole datasets
################

r_datasets=(r_SEQ_HOLE_v1_uint32_uint32_16000000 r_SEQ_HOLE_v1_uint32_uint32_16000000 r_SEQ_HOLE_v1_uint32_uint32_16000000 r_SEQ_HOLE_v2_uint32_uint32_32000000 r_SEQ_HOLE_v2_uint32_uint32_32000000 r_SEQ_HOLE_v2_uint32_uint32_32000000 r_SEQ_HOLE_v3_uint32_uint32_128000000 r_SEQ_HOLE_v3_uint32_uint32_128000000 r_SEQ_HOLE_v3_uint32_uint32_128000000 r_SEQ_HOLE_v5_uint32_uint32_640000000 r_SEQ_HOLE_v5_uint32_uint32_640000000 r_SEQ_HOLE_v5_uint32_uint32_640000000) #(r_SEQ_HOLE_v1_uint32_uint32_16000000 r_SEQ_HOLE_v2_uint32_uint32_32000000 r_SEQ_HOLE_v3_uint32_uint32_128000000 r_SEQ_HOLE_v4_uint32_uint32_384000000 r_SEQ_HOLE_v5_uint32_uint32_640000000 r_SEQ_HOLE_v6_uint32_uint32_896000000 r_SEQ_HOLE_v7_uint32_uint32_1152000000 r_SEQ_HOLE_v8_uint32_uint32_1664000000 r_SEQ_HOLE_v9_uint32_uint32_1920000000)
s_datasets=(s_SEQ_HOLE_v2_uint32_uint32_32000000 s_SEQ_HOLE_v3_uint32_uint32_128000000 s_SEQ_HOLE_v5_uint32_uint32_640000000 s_SEQ_HOLE_v1_uint32_uint32_16000000 s_SEQ_HOLE_v3_uint32_uint32_128000000 s_SEQ_HOLE_v5_uint32_uint32_640000000 s_SEQ_HOLE_v1_uint32_uint32_16000000 s_SEQ_HOLE_v2_uint32_uint32_32000000 s_SEQ_HOLE_v5_uint32_uint32_640000000 s_SEQ_HOLE_v1_uint32_uint32_16000000 s_SEQ_HOLE_v2_uint32_uint32_32000000 s_SEQ_HOLE_v3_uint32_uint32_128000000) #(s_SEQ_HOLE_v1_uint32_uint32_16000000 s_SEQ_HOLE_v2_uint32_uint32_32000000 s_SEQ_HOLE_v3_uint32_uint32_128000000 s_SEQ_HOLE_v4_uint32_uint32_384000000 s_SEQ_HOLE_v5_uint32_uint32_640000000 s_SEQ_HOLE_v6_uint32_uint32_896000000 s_SEQ_HOLE_v7_uint32_uint32_1152000000 s_SEQ_HOLE_v8_uint32_uint32_1664000000 s_SEQ_HOLE_v9_uint32_uint32_1920000000)
r_datasets_sizes=(16E6 16E6 16E6 32E6 32E6 32E6 128E6 128E6 128E6 640E6 640E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(32E6 128E6 640E6 16E6 128E6 640E6 16E6 32E6 640E6 16E6 32E6 128E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32 32 32 32 32 32 32 32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32 32 32 32 32 32 32 32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)

#output_folder_path=/spinning/sabek/learned_join_results/sj_with_eth_seq_hole/
#process_sort_merge_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation


r_datasets=(r_SEQ_HOLE_v1_uint32_uint32_16000000 r_SEQ_HOLE_v2_uint32_uint32_32000000 r_SEQ_HOLE_v3_uint32_uint32_128000000 r_SEQ_HOLE_v5_uint32_uint32_640000000) #(r_SEQ_HOLE_v1_uint32_uint32_16000000 r_SEQ_HOLE_v2_uint32_uint32_32000000 r_SEQ_HOLE_v3_uint32_uint32_128000000 r_SEQ_HOLE_v4_uint32_uint32_384000000 r_SEQ_HOLE_v5_uint32_uint32_640000000 r_SEQ_HOLE_v6_uint32_uint32_896000000 r_SEQ_HOLE_v7_uint32_uint32_1152000000 r_SEQ_HOLE_v8_uint32_uint32_1664000000 r_SEQ_HOLE_v9_uint32_uint32_1920000000)
s_datasets=(s_SEQ_HOLE_v1_uint32_uint32_16000000 s_SEQ_HOLE_v2_uint32_uint32_32000000 s_SEQ_HOLE_v3_uint32_uint32_128000000 s_SEQ_HOLE_v5_uint32_uint32_640000000) #(s_SEQ_HOLE_v1_uint32_uint32_16000000 s_SEQ_HOLE_v2_uint32_uint32_32000000 s_SEQ_HOLE_v3_uint32_uint32_128000000 s_SEQ_HOLE_v4_uint32_uint32_384000000 s_SEQ_HOLE_v5_uint32_uint32_640000000 s_SEQ_HOLE_v6_uint32_uint32_896000000 s_SEQ_HOLE_v7_uint32_uint32_1152000000 s_SEQ_HOLE_v8_uint32_uint32_1664000000 s_SEQ_HOLE_v9_uint32_uint32_1920000000)
r_datasets_sizes=(16E6 32E6 128E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(16E6 32E6 128E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32 32 32 32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32 32 32 32) #(64 64 64 64 64 64 64 64 64)

#output_folder_path=/spinning/sabek/learned_join_results/sj_with_eth_seq_hole/
#process_sort_merge_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation


#uniform datasets
################

r_datasets=(r_UNIFORM_v2_uint32_uint32_32000000 r_UNIFORM_v2_uint32_uint32_32000000  r_UNIFORM_v3_uint32_uint32_128000000 r_UNIFORM_v3_uint32_uint32_128000000 r_UNIFORM_v5_uint32_uint32_640000000 r_UNIFORM_v5_uint32_uint32_640000000) #(r_UNIFORM_v1_uint32_uint32_16000000 r_UNIFORM_v2_uint32_uint32_32000000 r_UNIFORM_v3_uint32_uint32_128000000 r_UNIFORM_v4_uint32_uint32_384000000 r_UNIFORM_v5_uint32_uint32_640000000 r_UNIFORM_v6_uint32_uint32_896000000 r_UNIFORM_v7_uint32_uint32_1152000000 r_UNIFORM_v8_uint32_uint32_1664000000 r_UNIFORM_v9_uint32_uint32_1920000000)
s_datasets=(s_UNIFORM_v3_uint32_uint32_128000000 s_UNIFORM_v5_uint32_uint32_640000000 s_UNIFORM_v2_uint32_uint32_32000000 s_UNIFORM_v5_uint32_uint32_640000000 s_UNIFORM_v2_uint32_uint32_32000000 s_UNIFORM_v3_uint32_uint32_128000000) #(s_UNIFORM_v1_uint32_uint32_16000000 s_UNIFORM_v2_uint32_uint32_32000000 s_UNIFORM_v3_uint32_uint32_128000000 s_UNIFORM_v4_uint32_uint32_384000000 s_UNIFORM_v5_uint32_uint32_640000000 s_UNIFORM_v6_uint32_uint32_896000000 s_UNIFORM_v7_uint32_uint32_1152000000 s_UNIFORM_v8_uint32_uint32_1664000000 s_UNIFORM_v9_uint32_uint32_1920000000)
r_datasets_sizes=(32E6 32E6 128E6 128E6 640E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(128E6 640E6 32E6 640E6 32E6 128E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32 32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32 32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)

#output_folder_path=/spinning/sabek/learned_join_results/sj_with_eth_uniform/
#process_sort_merge_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation


r_datasets=(r_UNIFORM_v1_uint32_uint32_16000000 r_UNIFORM_v2_uint32_uint32_32000000 r_UNIFORM_v3_uint32_uint32_128000000 r_UNIFORM_v5_uint32_uint32_640000000) #(r_UNIFORM_v1_uint32_uint32_16000000 r_UNIFORM_v2_uint32_uint32_32000000 r_UNIFORM_v3_uint32_uint32_128000000 r_UNIFORM_v4_uint32_uint32_384000000 r_UNIFORM_v5_uint32_uint32_640000000 r_UNIFORM_v6_uint32_uint32_896000000 r_UNIFORM_v7_uint32_uint32_1152000000 r_UNIFORM_v8_uint32_uint32_1664000000 r_UNIFORM_v9_uint32_uint32_1920000000)
s_datasets=(s_UNIFORM_v1_uint32_uint32_16000000 s_UNIFORM_v2_uint32_uint32_32000000 s_UNIFORM_v3_uint32_uint32_128000000 s_UNIFORM_v5_uint32_uint32_640000000) #(s_UNIFORM_v1_uint32_uint32_16000000 s_UNIFORM_v2_uint32_uint32_32000000 s_UNIFORM_v3_uint32_uint32_128000000 s_UNIFORM_v4_uint32_uint32_384000000 s_UNIFORM_v5_uint32_uint32_640000000 s_UNIFORM_v6_uint32_uint32_896000000 s_UNIFORM_v7_uint32_uint32_1152000000 s_UNIFORM_v8_uint32_uint32_1664000000 s_UNIFORM_v9_uint32_uint32_1920000000)
r_datasets_sizes=(16E6 32E6 128E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(16E6 32E6 128E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32 32 32 32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32 32 32 32) #(64 64 64 64 64 64 64 64 64)

#output_folder_path=/spinning/sabek/learned_join_results/sj_with_eth_uniform/
#process_sort_merge_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation


#sosd datasets
################

r_datasets=(books_200M_uint32) #(books_200M_uint32 books_800M_uint64 fb_200M_uint64 osm_cellids_800M_uint64 wiki_ts_200M_uint64) 
s_datasets=(books_200M_uint32) #(books_200M_uint32 books_800M_uint64 fb_200M_uint64 osm_cellids_800M_uint64 wiki_ts_200M_uint64)
r_datasets_sizes=(200E6) #(200E6 800E6 200E6 800E6 200E6)
s_datasets_sizes=(200E6) #(200E6 800E6 200E6 800E6 200E6)
r_datasets_file_num_partitions=(32) #(32 32 32 32 32)
s_datasets_file_num_partitions=(32) #(32 32 32 32 32)

#output_folder_path=/spinning/sabek/learned_join_results/sj_with_eth_sosd_books_200M_uint32/
#process_sort_merge_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation

r_datasets=(books_800M_uint64) #(books_200M_uint32 books_800M_uint64 fb_200M_uint64 osm_cellids_800M_uint64 wiki_ts_200M_uint64) 
s_datasets=(books_800M_uint64) #(books_200M_uint32 books_800M_uint64 fb_200M_uint64 osm_cellids_800M_uint64 wiki_ts_200M_uint64)
r_datasets_sizes=(800E6) #(200E6 800E6 200E6 800E6 200E6)
s_datasets_sizes=(800E6) #(200E6 800E6 200E6 800E6 200E6)
r_datasets_file_num_partitions=(32) #(32 32 32 32 32)
s_datasets_file_num_partitions=(32) #(32 32 32 32 32)

#output_folder_path=/spinning/sabek/learned_join_results/sj_with_eth_sosd_books_800M_uint64/
#process_sort_merge_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation


r_datasets=(fb_200M_uint64) #(books_200M_uint32 books_800M_uint64 fb_200M_uint64 osm_cellids_800M_uint64 wiki_ts_200M_uint64) 
s_datasets=(fb_200M_uint64) #(books_200M_uint32 books_800M_uint64 fb_200M_uint64 osm_cellids_800M_uint64 wiki_ts_200M_uint64)
r_datasets_sizes=(200E6) #(200E6 800E6 200E6 800E6 200E6)
s_datasets_sizes=(200E6) #(200E6 800E6 200E6 800E6 200E6)
r_datasets_file_num_partitions=(32) #(32 32 32 32 32)
s_datasets_file_num_partitions=(32) #(32 32 32 32 32)

#output_folder_path=/spinning/sabek/learned_join_results/sj_with_eth_sosd_fb_200M_uint64/
#process_sort_merge_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation


r_datasets=(osm_cellids_800M_uint64) #(books_200M_uint32 books_800M_uint64 fb_200M_uint64 osm_cellids_800M_uint64 wiki_ts_200M_uint64) 
s_datasets=(osm_cellids_800M_uint64) #(books_200M_uint32 books_800M_uint64 fb_200M_uint64 osm_cellids_800M_uint64 wiki_ts_200M_uint64)
r_datasets_sizes=(800E6) #(200E6 800E6 200E6 800E6 200E6)
s_datasets_sizes=(800E6) #(200E6 800E6 200E6 800E6 200E6)
r_datasets_file_num_partitions=(32) #(32 32 32 32 32)
s_datasets_file_num_partitions=(32) #(32 32 32 32 32)

#output_folder_path=/spinning/sabek/learned_join_results/sj_with_eth_sosd_osm_cellids_800M_uint64/
#process_sort_merge_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation


r_datasets=(wiki_ts_200M_uint64) #(books_200M_uint32 books_800M_uint64 fb_200M_uint64 osm_cellids_800M_uint64 wiki_ts_200M_uint64) 
s_datasets=(wiki_ts_200M_uint64) #(books_200M_uint32 books_800M_uint64 fb_200M_uint64 osm_cellids_800M_uint64 wiki_ts_200M_uint64)
r_datasets_sizes=(200E6) #(200E6 800E6 200E6 800E6 200E6)
s_datasets_sizes=(200E6) #(200E6 800E6 200E6 800E6 200E6)
r_datasets_file_num_partitions=(32) #(32 32 32 32 32)
s_datasets_file_num_partitions=(32) #(32 32 32 32 32)

#output_folder_path=/spinning/sabek/learned_join_results/sj_with_eth_sosd_wiki_ts_200M_uint64/
#process_sort_merge_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation

#tpch workloads
################

r_datasets=(r_q3_v1_uint32_uint32_178 r_q3_v2_uint32_uint32_87)
s_datasets=(s_q3_v1_uint32_uint32_5787 s_q3_v2_uint32_uint32_233105) 
r_datasets_sizes=(178 87) 
s_datasets_sizes=(5787 233105)
r_datasets_file_num_partitions=(32 32) 
s_datasets_file_num_partitions=(32 32)
input_hash_table_size=(16777216 16777216) #(16777216(for_16E6) 33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

#output_folder_path=/spinning/sabek/learned_join_results/sj_with_eth_tpch_q3/
#process_sort_merge_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation


r_datasets=(r_q5_v1_uint32_uint32_1 r_q5_v2_uint32_uint32_5 r_q5_v3_uint32_uint32_191 r_q5_v4_uint32_uint32_88 r_q5_v5_uint32_uint32_38)
s_datasets=(s_q5_v1_uint32_uint32_20 s_q5_v2_uint32_uint32_885 s_q5_v3_uint32_uint32_131 s_q5_v4_uint32_uint32_865032 s_q5_v5_uint32_uint32_690) 
r_datasets_sizes=(1 5 191 88 38) 
s_datasets_sizes=(20 885 131 865032 690)
r_datasets_file_num_partitions=(32 32 32 32 32) 
s_datasets_file_num_partitions=(32 32 32 32 32)
input_hash_table_size=(16777216 16777216 16777216 16777216 16777216) #(16777216(for_16E6) 33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

#output_folder_path=/spinning/sabek/learned_join_results/sj_with_eth_tpch_q5/
#process_sort_merge_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation


r_datasets=(r_q9_v1_uint32_uint32_22 r_q9_v2_uint32_uint32_47 r_q9_v3_uint32_uint32_705 r_q9_v4_uint32_uint32_73)
s_datasets=(s_q9_v1_uint32_uint32_690 s_q9_v2_uint32_uint32_11078 s_q9_v3_uint32_uint32_38 s_q9_v4_uint32_uint32_742875) 
r_datasets_sizes=(22 47 705 73) 
s_datasets_sizes=(690 11078 38 742875)
r_datasets_file_num_partitions=(32 32 32 32) 
s_datasets_file_num_partitions=(32 32 32 32)
input_hash_table_size=(16777216 16777216 16777216 16777216) #(16777216(for_16E6) 33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

#output_folder_path=/spinning/sabek/learned_join_results/sj_with_eth_tpch_q9/
#process_sort_merge_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation


r_datasets=(r_q18_v1_uint32_uint32_1 r_q18_v2_uint32_uint32_922 r_q18_v3_uint32_uint32_1)
s_datasets=(s_q18_v1_uint32_uint32_41501 s_q18_v2_uint32_uint32_1 s_q18_v3_uint32_uint32_752845) 
r_datasets_sizes=(1 922 1) 
s_datasets_sizes=(41501 1 752845)
r_datasets_file_num_partitions=(32 32 32) 
s_datasets_file_num_partitions=(32 32 32)
input_hash_table_size=(16777216 16777216 16777216) #(16777216(for_16E6) 33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

#output_folder_path=/spinning/sabek/learned_join_results/sj_with_eth_tpch_q18/
#process_sort_merge_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation

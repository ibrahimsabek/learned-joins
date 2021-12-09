# A script to run the experiments for the non-partitioned hash join (including learned) algorithms

#!/bin/sh

process_non_partition_join()
{
    threads=(32) #(2 4 8 16 32 64)
    bucket_size=(2) #2 4 8 10
    prefetch_distance=(32) #(4 16 32)
    #npj_eth_avx_imv=1 #(0 1)
    use_murmur3_hash=1 #(0 1)
    curr_ls_default_threshold=2097152 #100 1000 2000 10000 20000 30000 40000 60000 80000 
    curr_ls_default_arch=2097152 #1000 2000 5000 10000 20000 50000 100000 200000 500000 800000 1000000

    dataset_folder_path=/spinning/sabek/learned_join_datasets/

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

    run_non_learned=${11}
    run_learned=${12}

    if [ ${run_non_learned} == 1 ]
    then
        echo "Running NPJ with non-learned ..."

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
                for bs in ${!bucket_size[@]}
                do
                    curr_bucket_size=${bucket_size[$bs]}
                    for pd in ${!prefetch_distance[@]}
                    do
                        curr_prefetch_distance=${prefetch_distance[$pd]}

                        curr_output_file=$output_folder_path'npj_tuning_'$curr_r_dataset_size'_'$curr_s_dataset_size'_th_'$curr_threads'_bs_'$curr_bucket_size'_pd_'$curr_prefetch_distance'.csv'
        
                        sh $(dirname "$0")/base_configs_maker.sh -NUM_THREADS_FOR_EVALUATION $curr_threads \
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
                                                                -RUN_LEARNED_TECHNIQUES 0 \
                                                                -RUN_LEARNED_TECHNIQUES_WITH_FIRST_LEVEL_ONLY 0 \
                                                                -LS_FOR_SORT_MERGE_DEFAULT_FANOUT $curr_ls_default_threshold \
                                                                -LS_FOR_SORT_MERGE_DEFAULT_ARCH_SECOND_LEVEL $curr_ls_default_arch \
                                                                -LS_FOR_SORT_MERGE_DEFAULT_THRESHOLD $curr_ls_default_threshold \
                                                                -LS_FOR_SORT_MERGE_DEFAULT_SAMPLING_RATE 0.01 \
                                                                -CUSTOM_CPU_MAPPING '"'../../include/configs/cpu-mapping_berners_lee.txt'"' \
                                                                -CUSTOM_CPU_MAPPING_V2 '"'../../include/configs/cpu-mapping-v2_berners_lee.txt'"'

                        sh $(dirname "$0")/eth_configs_maker.sh -BUCKET_SIZE $curr_bucket_size \
                                                                -PREFETCH_DISTANCE $curr_prefetch_distance \
                                                                -USE_MURMUR3_HASH $use_murmur3_hash \
                                                                -NPJ_ETH_AVX_IMV 1 \
                                                                -NPJ_SIMDStateSize 5


                        cmake -DCMAKE_BUILD_TYPE=Release -DVECTORWISE_BRANCHING=on $(dirname "$0")/../.. > /dev/null

                        cd $(dirname "$0")/../../build/release

                        make > /dev/null

                        ./non_partition_join_runner

                        cd ../../scripts/evaluation/

                    done
                done   
            done
        done    
    fi

    if [ ${run_learned} == 1 ]
    then
        echo "Running NPJ with learned ..."

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
                for bs in ${!bucket_size[@]}
                do
                    curr_bucket_size=${bucket_size[$bs]}
                    for pd in ${!prefetch_distance[@]}
                    do
                        curr_prefetch_distance=${prefetch_distance[$pd]}

                        curr_output_file=$output_folder_path'npj_tuning_'$curr_r_dataset_size'_'$curr_s_dataset_size'_th_'$curr_threads'_bs_'$curr_bucket_size'_pd_'$curr_prefetch_distance'.csv'
        
                        sh $(dirname "$0")/base_configs_maker.sh -NUM_THREADS_FOR_EVALUATION $curr_threads \
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
                                                                -RUN_LEARNED_TECHNIQUES 1 \
                                                                -RUN_LEARNED_TECHNIQUES_WITH_FIRST_LEVEL_ONLY 0 \
                                                                -LS_FOR_SORT_MERGE_DEFAULT_FANOUT $curr_ls_default_threshold \
                                                                -LS_FOR_SORT_MERGE_DEFAULT_ARCH_SECOND_LEVEL $curr_ls_default_arch \
                                                                -LS_FOR_SORT_MERGE_DEFAULT_THRESHOLD $curr_ls_default_threshold \
                                                                -LS_FOR_SORT_MERGE_DEFAULT_SAMPLING_RATE 0.01 \
                                                                -CUSTOM_CPU_MAPPING '"'../../include/configs/cpu-mapping_berners_lee.txt'"' \
                                                                -CUSTOM_CPU_MAPPING_V2 '"'../../include/configs/cpu-mapping-v2_berners_lee.txt'"'

                        sh $(dirname "$0")/eth_configs_maker.sh -BUCKET_SIZE $curr_bucket_size \
                                                                -PREFETCH_DISTANCE $curr_prefetch_distance \
                                                                -USE_MURMUR3_HASH $use_murmur3_hash \
                                                                -NPJ_ETH_AVX_IMV 1 \
                                                                -NPJ_SIMDStateSize 5


                        cmake -DCMAKE_BUILD_TYPE=Release -DVECTORWISE_BRANCHING=on $(dirname "$0")/../.. > /dev/null

                        cd $(dirname "$0")/../../build/release

                        make > /dev/null

                        ./non_partition_join_runner

                        cd ../../scripts/evaluation/

                    done
                done   
            done
        done        
    fi
}


run_nums=1 #5 10
load_relations_for_evaluation=0 #1
persist_relations_for_evaluation=0

#unique datasets
################

r_datasets=(r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v5_uint32_uint32_640000000) #(r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v4_uint32_uint32_384000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v6_uint32_uint32_896000000 r_UNIQUE_v7_uint32_uint32_1152000000 r_UNIQUE_v8_uint32_uint32_1664000000 r_UNIQUE_v9_uint32_uint32_1920000000)
s_datasets=(s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000) #(s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v4_uint32_uint32_384000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v6_uint32_uint32_896000000 s_UNIQUE_v7_uint32_uint32_1152000000 s_UNIQUE_v8_uint32_uint32_1664000000 s_UNIQUE_v9_uint32_uint32_1920000000)
r_datasets_sizes=(16E6 16E6 16E6 32E6 32E6 32E6 128E6 128E6 128E6 640E6 640E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(32E6 128E6 640E6 16E6 128E6 640E6 16E6 32E6 640E6 16E6 32E6 128E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32 32 32 32 32 32 32 32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32 32 32 32 32 32 32 32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)

#output_folder_path=/spinning/sabek/learned_join_results/npj_with_learned_unique/
#process_non_partition_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 1
#output_folder_path=/spinning/sabek/learned_join_results/npj_with_non_learned_unique/
#process_non_partition_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0


#r_datasets=(r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v8_uint32_uint32_1664000000) #(r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v4_uint32_uint32_384000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v6_uint32_uint32_896000000 r_UNIQUE_v7_uint32_uint32_1152000000 r_UNIQUE_v8_uint32_uint32_1664000000 r_UNIQUE_v9_uint32_uint32_1920000000)
#s_datasets=(s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v8_uint32_uint32_1664000000) #(s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v4_uint32_uint32_384000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v6_uint32_uint32_896000000 s_UNIQUE_v7_uint32_uint32_1152000000 s_UNIQUE_v8_uint32_uint32_1664000000 s_UNIQUE_v9_uint32_uint32_1920000000)
#r_datasets_sizes=(32E6 128E6 640E6 1664E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
#s_datasets_sizes=(32E6 128E6 640E6 1664E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
#r_datasets_file_num_partitions=(32 32 32 32) #(64 64 64 64 64 64 64 64 64)
#s_datasets_file_num_partitions=(32 32 32 32) #(64 64 64 64 64 64 64 64 64)

r_datasets=(r_UNIQUE_v3_uint32_uint32_128000000) #(r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v4_uint32_uint32_384000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v6_uint32_uint32_896000000 r_UNIQUE_v7_uint32_uint32_1152000000 r_UNIQUE_v8_uint32_uint32_1664000000 r_UNIQUE_v9_uint32_uint32_1920000000)
s_datasets=(s_UNIQUE_v3_uint32_uint32_128000000) #(s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v4_uint32_uint32_384000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v6_uint32_uint32_896000000 s_UNIQUE_v7_uint32_uint32_1152000000 s_UNIQUE_v8_uint32_uint32_1664000000 s_UNIQUE_v9_uint32_uint32_1920000000)
r_datasets_sizes=(128E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(128E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32) #(64 64 64 64 64 64 64 64 64)

#output_folder_path=/spinning/sabek/learned_join_results/npj_with_learned_unique/
#process_non_partition_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 1
output_folder_path=/spinning/sabek/learned_join_results/npj_with_non_learned_unique/
process_non_partition_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0


#seq_hole datasets
################

r_datasets=(r_SEQ_HOLE_v1_uint32_uint32_16000000 r_SEQ_HOLE_v1_uint32_uint32_16000000 r_SEQ_HOLE_v1_uint32_uint32_16000000 r_SEQ_HOLE_v2_uint32_uint32_32000000 r_SEQ_HOLE_v2_uint32_uint32_32000000 r_SEQ_HOLE_v2_uint32_uint32_32000000 r_SEQ_HOLE_v3_uint32_uint32_128000000 r_SEQ_HOLE_v3_uint32_uint32_128000000 r_SEQ_HOLE_v3_uint32_uint32_128000000 r_SEQ_HOLE_v5_uint32_uint32_640000000 r_SEQ_HOLE_v5_uint32_uint32_640000000 r_SEQ_HOLE_v5_uint32_uint32_640000000) #(r_SEQ_HOLE_v1_uint32_uint32_16000000 r_SEQ_HOLE_v2_uint32_uint32_32000000 r_SEQ_HOLE_v3_uint32_uint32_128000000 r_SEQ_HOLE_v4_uint32_uint32_384000000 r_SEQ_HOLE_v5_uint32_uint32_640000000 r_SEQ_HOLE_v6_uint32_uint32_896000000 r_SEQ_HOLE_v7_uint32_uint32_1152000000 r_SEQ_HOLE_v8_uint32_uint32_1664000000 r_SEQ_HOLE_v9_uint32_uint32_1920000000)
s_datasets=(s_SEQ_HOLE_v2_uint32_uint32_32000000 s_SEQ_HOLE_v3_uint32_uint32_128000000 s_SEQ_HOLE_v5_uint32_uint32_640000000 s_SEQ_HOLE_v1_uint32_uint32_16000000 s_SEQ_HOLE_v3_uint32_uint32_128000000 s_SEQ_HOLE_v5_uint32_uint32_640000000 s_SEQ_HOLE_v1_uint32_uint32_16000000 s_SEQ_HOLE_v2_uint32_uint32_32000000 s_SEQ_HOLE_v5_uint32_uint32_640000000 s_SEQ_HOLE_v1_uint32_uint32_16000000 s_SEQ_HOLE_v2_uint32_uint32_32000000 s_SEQ_HOLE_v3_uint32_uint32_128000000) #(s_SEQ_HOLE_v1_uint32_uint32_16000000 s_SEQ_HOLE_v2_uint32_uint32_32000000 s_SEQ_HOLE_v3_uint32_uint32_128000000 s_SEQ_HOLE_v4_uint32_uint32_384000000 s_SEQ_HOLE_v5_uint32_uint32_640000000 s_SEQ_HOLE_v6_uint32_uint32_896000000 s_SEQ_HOLE_v7_uint32_uint32_1152000000 s_SEQ_HOLE_v8_uint32_uint32_1664000000 s_SEQ_HOLE_v9_uint32_uint32_1920000000)
r_datasets_sizes=(16E6 16E6 16E6 32E6 32E6 32E6 128E6 128E6 128E6 640E6 640E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(32E6 128E6 640E6 16E6 128E6 640E6 16E6 32E6 640E6 16E6 32E6 128E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32 32 32 32 32 32 32 32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32 32 32 32 32 32 32 32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)

#output_folder_path=/spinning/sabek/learned_join_results/npj_with_learned_seq_hole/
#process_non_partition_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 1
#output_folder_path=/spinning/sabek/learned_join_results/npj_with_non_learned_seq_hole/
#process_non_partition_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0


r_datasets=(r_SEQ_HOLE_v2_uint32_uint32_32000000 r_SEQ_HOLE_v3_uint32_uint32_128000000 r_SEQ_HOLE_v5_uint32_uint32_640000000) #(r_SEQ_HOLE_v1_uint32_uint32_16000000 r_SEQ_HOLE_v2_uint32_uint32_32000000 r_SEQ_HOLE_v3_uint32_uint32_128000000 r_SEQ_HOLE_v4_uint32_uint32_384000000 r_SEQ_HOLE_v5_uint32_uint32_640000000 r_SEQ_HOLE_v6_uint32_uint32_896000000 r_SEQ_HOLE_v7_uint32_uint32_1152000000 r_SEQ_HOLE_v8_uint32_uint32_1664000000 r_SEQ_HOLE_v9_uint32_uint32_1920000000)
s_datasets=(s_SEQ_HOLE_v2_uint32_uint32_32000000 s_SEQ_HOLE_v3_uint32_uint32_128000000 s_SEQ_HOLE_v5_uint32_uint32_640000000) #(s_SEQ_HOLE_v1_uint32_uint32_16000000 s_SEQ_HOLE_v2_uint32_uint32_32000000 s_SEQ_HOLE_v3_uint32_uint32_128000000 s_SEQ_HOLE_v4_uint32_uint32_384000000 s_SEQ_HOLE_v5_uint32_uint32_640000000 s_SEQ_HOLE_v6_uint32_uint32_896000000 s_SEQ_HOLE_v7_uint32_uint32_1152000000 s_SEQ_HOLE_v8_uint32_uint32_1664000000 s_SEQ_HOLE_v9_uint32_uint32_1920000000)
r_datasets_sizes=(32E6 128E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(32E6 128E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32 32 32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32 32 32) #(64 64 64 64 64 64 64 64 64)

#output_folder_path=/spinning/sabek/learned_join_results/npj_with_learned_seq_hole/
#process_non_partition_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 1
#output_folder_path=/spinning/sabek/learned_join_results/npj_with_non_learned_seq_hole/
#process_non_partition_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0

#uniform datasets
################

r_datasets=(r_UNIFORM_v2_uint32_uint32_32000000 r_UNIFORM_v2_uint32_uint32_32000000  r_UNIFORM_v3_uint32_uint32_128000000 r_UNIFORM_v3_uint32_uint32_128000000 r_UNIFORM_v5_uint32_uint32_640000000 r_UNIFORM_v5_uint32_uint32_640000000) #(r_UNIFORM_v1_uint32_uint32_16000000 r_UNIFORM_v2_uint32_uint32_32000000 r_UNIFORM_v3_uint32_uint32_128000000 r_UNIFORM_v4_uint32_uint32_384000000 r_UNIFORM_v5_uint32_uint32_640000000 r_UNIFORM_v6_uint32_uint32_896000000 r_UNIFORM_v7_uint32_uint32_1152000000 r_UNIFORM_v8_uint32_uint32_1664000000 r_UNIFORM_v9_uint32_uint32_1920000000)
s_datasets=(s_UNIFORM_v3_uint32_uint32_128000000 s_UNIFORM_v5_uint32_uint32_640000000 s_UNIFORM_v2_uint32_uint32_32000000 s_UNIFORM_v5_uint32_uint32_640000000 s_UNIFORM_v2_uint32_uint32_32000000 s_UNIFORM_v3_uint32_uint32_128000000) #(s_UNIFORM_v1_uint32_uint32_16000000 s_UNIFORM_v2_uint32_uint32_32000000 s_UNIFORM_v3_uint32_uint32_128000000 s_UNIFORM_v4_uint32_uint32_384000000 s_UNIFORM_v5_uint32_uint32_640000000 s_UNIFORM_v6_uint32_uint32_896000000 s_UNIFORM_v7_uint32_uint32_1152000000 s_UNIFORM_v8_uint32_uint32_1664000000 s_UNIFORM_v9_uint32_uint32_1920000000)
r_datasets_sizes=(32E6 32E6 128E6 128E6 640E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(128E6 640E6 32E6 640E6 32E6 128E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32 32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32 32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)

#output_folder_path=/spinning/sabek/learned_join_results/npj_with_learned_uniform/
#process_non_partition_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 1
#output_folder_path=/spinning/sabek/learned_join_results/npj_with_non_learned_uniform/
#process_non_partition_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0


r_datasets=(r_UNIFORM_v2_uint32_uint32_32000000 r_UNIFORM_v3_uint32_uint32_128000000 r_UNIFORM_v5_uint32_uint32_640000000) #(r_UNIFORM_v1_uint32_uint32_16000000 r_UNIFORM_v2_uint32_uint32_32000000 r_UNIFORM_v3_uint32_uint32_128000000 r_UNIFORM_v4_uint32_uint32_384000000 r_UNIFORM_v5_uint32_uint32_640000000 r_UNIFORM_v6_uint32_uint32_896000000 r_UNIFORM_v7_uint32_uint32_1152000000 r_UNIFORM_v8_uint32_uint32_1664000000 r_UNIFORM_v9_uint32_uint32_1920000000)
s_datasets=(s_UNIFORM_v2_uint32_uint32_32000000 s_UNIFORM_v3_uint32_uint32_128000000 s_UNIFORM_v5_uint32_uint32_640000000) #(s_UNIFORM_v1_uint32_uint32_16000000 s_UNIFORM_v2_uint32_uint32_32000000 s_UNIFORM_v3_uint32_uint32_128000000 s_UNIFORM_v4_uint32_uint32_384000000 s_UNIFORM_v5_uint32_uint32_640000000 s_UNIFORM_v6_uint32_uint32_896000000 s_UNIFORM_v7_uint32_uint32_1152000000 s_UNIFORM_v8_uint32_uint32_1664000000 s_UNIFORM_v9_uint32_uint32_1920000000)
r_datasets_sizes=(32E6 128E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(32E6 128E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32 32 32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32 32 32) #(64 64 64 64 64 64 64 64 64)

#output_folder_path=/spinning/sabek/learned_join_results/npj_with_learned_uniform/
#process_non_partition_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 1
#output_folder_path=/spinning/sabek/learned_join_results/npj_with_non_learned_uniform/
#process_non_partition_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0

#lognormal datasets
###################

r_datasets=(r_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 r_LOGNORMAL_v2_uint32_uint32_segma_1_32000000  r_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 r_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 r_LOGNORMAL_v5_uint32_uint32_segma_1_640000000 r_LOGNORMAL_v5_uint32_uint32_segma_1_640000000) #(r_LOGNORMAL_v1_uint32_uint32_segma_1_16000000 r_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 r_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 r_LOGNORMAL_v4_uint32_uint32_segma_1_384000000 r_LOGNORMAL_v5_uint32_uint32_segma_1_640000000 r_LOGNORMAL_v6_uint32_uint32_segma_1_896000000 r_LOGNORMAL_v7_uint32_uint32_segma_1_1152000000 r_LOGNORMAL_v8_uint32_uint32_segma_1_1664000000 r_LOGNORMAL_v9_uint32_uint32_segma_1_1920000000)
s_datasets=(s_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 s_LOGNORMAL_v5_uint32_uint32_segma_1_640000000 s_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 s_LOGNORMAL_v5_uint32_uint32_segma_1_640000000 s_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 s_LOGNORMAL_v3_uint32_uint32_segma_1_128000000) #(s_LOGNORMAL_v1_uint32_uint32_segma_1_16000000 s_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 s_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 s_LOGNORMAL_v4_uint32_uint32_segma_1_384000000 s_LOGNORMAL_v5_uint32_uint32_segma_1_640000000 s_LOGNORMAL_v6_uint32_uint32_segma_1_896000000 s_LOGNORMAL_v7_uint32_uint32_segma_1_1152000000 s_LOGNORMAL_v8_uint32_uint32_segma_1_1664000000 s_LOGNORMAL_v9_uint32_uint32_segma_1_1920000000)
r_datasets_sizes=(32E6 32E6 128E6 128E6 640E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(128E6 640E6 32E6 640E6 32E6 128E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32 32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32 32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)

#output_folder_path=/spinning/sabek/learned_join_results/npj_with_learned_lognormal/
#process_non_partition_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 1
#output_folder_path=/spinning/sabek/learned_join_results/npj_with_non_learned_lognormal/
#process_non_partition_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0


r_datasets=(r_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 r_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 r_LOGNORMAL_v5_uint32_uint32_segma_1_640000000) #(r_LOGNORMAL_v1_uint32_uint32_segma_1_16000000 r_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 r_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 r_LOGNORMAL_v4_uint32_uint32_segma_1_384000000 r_LOGNORMAL_v5_uint32_uint32_segma_1_640000000 r_LOGNORMAL_v6_uint32_uint32_segma_1_896000000 r_LOGNORMAL_v7_uint32_uint32_segma_1_1152000000 r_LOGNORMAL_v8_uint32_uint32_segma_1_1664000000 r_LOGNORMAL_v9_uint32_uint32_segma_1_1920000000)
s_datasets=(s_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 s_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 s_LOGNORMAL_v5_uint32_uint32_segma_1_640000000) #(s_LOGNORMAL_v1_uint32_uint32_segma_1_16000000 s_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 s_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 s_LOGNORMAL_v4_uint32_uint32_segma_1_384000000 s_LOGNORMAL_v5_uint32_uint32_segma_1_640000000 s_LOGNORMAL_v6_uint32_uint32_segma_1_896000000 s_LOGNORMAL_v7_uint32_uint32_segma_1_1152000000 s_LOGNORMAL_v8_uint32_uint32_segma_1_1664000000 s_LOGNORMAL_v9_uint32_uint32_segma_1_1920000000)
r_datasets_sizes=(32E6 128E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(32E6 128E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32 32 32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32 32 32) #(64 64 64 64 64 64 64 64 64)

#output_folder_path=/spinning/sabek/learned_join_results/npj_with_learned_lognormal/
#process_non_partition_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 1
#output_folder_path=/spinning/sabek/learned_join_results/npj_with_non_learned_lognormal/
#process_non_partition_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0

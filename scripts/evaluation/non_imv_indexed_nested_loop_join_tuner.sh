#!/bin/sh

process_non_imv_indexed_nested_loop_join()
{
    threads=(32) #(2 4 8 16 32 64)
    rmi_models=(0 1 2 3 4 5 6 7 8 9)
    css_fanouts=(33) #(10 33 40)

    dataset_folder_path=/spinning/sabek/learned_join_datasets_tpch/
    #dataset_folder_path=/spinning/sabek/learned_join_datasets_sosd/
    #dataset_folder_path=/spinning/sabek/learned_join_datasets/

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

    run_inlj_with_hash_index=${11}
    run_inlj_with_learned_index=${12}
    run_inlj_with_csstree_index=${13}
    run_inlj_with_art32tree_index=${14}

    #input_hash_table_size=${15}

    if [ ${run_inlj_with_hash_index} == 1 ]
    then
        echo "Running INLJ with hash index ..."
        
        for ds in ${!r_datasets[@]}
        do
            curr_r_dataset='"'$dataset_folder_path${r_datasets[$ds]}.txt'"'
            curr_r_dataset_size=${r_datasets_sizes[$ds]}
            curr_r_dataset_file_num_partitions=${r_datasets_file_num_partitions[$ds]}
            curr_s_dataset='"'$dataset_folder_path${s_datasets[$ds]}.txt'"'
            curr_s_dataset_size=${s_datasets_sizes[$ds]}
            curr_s_dataset_file_num_partitions=${s_datasets_file_num_partitions[$ds]}

            curr_input_hash_table_size=${input_hash_table_size[$ds]}

            echo 'Joining '$curr_r_dataset' '$curr_r_dataset_size' '$curr_s_dataset' '$curr_s_dataset_size'... \n'
            
            for th in ${!threads[@]}
            do
                curr_threads=${threads[$th]}

                curr_output_file=$output_folder_path'non_imv_inlj_with_hash_index_tuning_'$curr_r_dataset_size'_'$curr_s_dataset_size'_th_'$curr_threads'_hts_'$curr_input_hash_table_size'.csv'
                
                sh $(dirname "$0")/base_configs_maker.sh -INLJ_WITH_HASH_INDEX 1 \
                                                -INLJ_WITH_LEARNED_INDEX 0 \
                                                -INLJ_WITH_CSS_TREE_INDEX 0 \
                                                -INLJ_WITH_ART32_TREE_INDEX 0 \
                                                -PREFETCH_INLJ 0 \
                                                -RUN_LEARNED_TECHNIQUES_WITH_FIRST_LEVEL_ONLY 1 \
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
                                                -PERSIST_RELATIONS_FOR_EVALUATION ${persist_relations_for_evaluation} #\
                                                #-CUSTOM_CPU_MAPPING '"'../../include/configs/cpu-mapping_berners_lee.txt'"' \
                                                #-CUSTOM_CPU_MAPPING_V2 '"'../../include/configs/cpu-mapping-v2_berners_lee.txt'"'

                sh $(dirname "$0")/eth_configs_maker.sh   -BUCKET_SIZE 1 \
                                                -PREFETCH_DISTANCE 128 \
                                                -USE_MURMUR3_HASH 1 \
                                                -INPUT_HASH_TABLE_SIZE $curr_input_hash_table_size

                cmake -DCMAKE_BUILD_TYPE=Release -DVECTORWISE_BRANCHING=on $(dirname "$0")/../.. > /dev/null

                cd $(dirname "$0")/../../build/release

                make > /dev/null

                ./non_imv_indexed_nested_loop_join_runner

                cd ../../scripts/evaluation/
            done
        done

    fi

    if [ ${run_inlj_with_learned_index} == 1 ]
    then
        echo "Running INLJ with learned index ..."

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

                for model in  ${!rmi_models[@]}
                do
                    #curr_rmi_model=${r_datasets[$ds]}'_key_uint32_'${rmi_models[$model]}
                    curr_rmi_model=${r_datasets[$ds]}'_'${rmi_models[$model]}
                    
                    curr_output_file=$output_folder_path'non_imv_inlj_with_learned_index_tuning_'$curr_r_dataset_size'_'$curr_s_dataset_size'_th_'$curr_threads'_rmi_'${rmi_models[$model]}'.csv'
                
                    sh $(dirname "$0")/base_configs_maker.sh -INLJ_WITH_HASH_INDEX 0 \
                                                    -INLJ_WITH_LEARNED_INDEX 1 \
                                                    -INLJ_WITH_CSS_TREE_INDEX 0 \
                                                    -INLJ_WITH_ART32_TREE_INDEX 0 \
                                                    -RUN_LEARNED_TECHNIQUES_WITH_FIRST_LEVEL_ONLY 1 \
                                                    -INLJ_RMI_DATA_PATH '"'/spinning/sabek/rmi_data'"' \
                                                    -INLJ_RMI_NAMESPACE ${curr_rmi_model} \
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
                                                    -PERSIST_RELATIONS_FOR_EVALUATION ${persist_relations_for_evaluation} #\
                                                    #-CUSTOM_CPU_MAPPING '"'../../include/configs/cpu-mapping_berners_lee.txt'"' \
                                                    #-CUSTOM_CPU_MAPPING_V2 '"'../../include/configs/cpu-mapping-v2_berners_lee.txt'"'

                    cmake -DCMAKE_BUILD_TYPE=Release -DVECTORWISE_BRANCHING=on $(dirname "$0")/../.. > /dev/null

                    cd $(dirname "$0")/../../build/release

                    make > /dev/null

                    ./non_imv_indexed_nested_loop_join_runner

                    cd ../../scripts/evaluation/
                done
            done
        done

    fi

    if [ ${run_inlj_with_csstree_index} == 1 ]
    then
        echo "Running INLJ with css tree index ..."

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

                for f in ${!css_fanouts[@]}
                do
                    curr_fanout=${css_fanouts[$f]}

                    curr_output_file=$output_folder_path'non_imv_inlj_with_csstree_index_tuning_'$curr_r_dataset_size'_'$curr_s_dataset_size'_th_'$curr_threads'_f_'$curr_fanout'.csv'

                    sh $(dirname "$0")/base_configs_maker.sh -INLJ_WITH_HASH_INDEX 0 \
                                                -INLJ_WITH_LEARNED_INDEX 0 \
                                                -INLJ_WITH_CSS_TREE_INDEX 1 \
                                                -INLJ_WITH_ART32_TREE_INDEX 0 \
                                                -INLJ_CSS_TREE_FANOUT $curr_fanout \
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
                                                -PERSIST_RELATIONS_FOR_EVALUATION ${persist_relations_for_evaluation} #\
                                                #-CUSTOM_CPU_MAPPING '"'../../include/configs/cpu-mapping_berners_lee.txt'"' \
                                                #-CUSTOM_CPU_MAPPING_V2 '"'../../include/configs/cpu-mapping-v2_berners_lee.txt'"'

                    cmake -DCMAKE_BUILD_TYPE=Release -DVECTORWISE_BRANCHING=on $(dirname "$0")/../.. > /dev/null

                    cd $(dirname "$0")/../../build/release

                    make > /dev/null

                    ./non_imv_indexed_nested_loop_join_runner

                    cd ../../scripts/evaluation/

                done
                
            done
        done

    fi

    if [ ${run_inlj_with_art32tree_index} == 1 ]
    then
        echo "Running INLJ with ART32 tree index ..."

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

                curr_output_file=$output_folder_path'non_imv_inlj_with_art32tree_index_tuning_'$curr_r_dataset_size'_'$curr_s_dataset_size'_th_'$curr_threads'.csv'

                sh $(dirname "$0")/base_configs_maker.sh -INLJ_WITH_HASH_INDEX 0 \
                                            -INLJ_WITH_LEARNED_INDEX 0 \
                                            -INLJ_WITH_CSS_TREE_INDEX 0 \
                                            -INLJ_WITH_ART32_TREE_INDEX 1 \
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
                                            -PERSIST_RELATIONS_FOR_EVALUATION ${persist_relations_for_evaluation} #\
                                            #-CUSTOM_CPU_MAPPING '"'../../include/configs/cpu-mapping_berners_lee.txt'"' \
                                            #-CUSTOM_CPU_MAPPING_V2 '"'../../include/configs/cpu-mapping-v2_berners_lee.txt'"'

                cmake -DCMAKE_BUILD_TYPE=Release -DVECTORWISE_BRANCHING=on $(dirname "$0")/../.. > /dev/null

                cd $(dirname "$0")/../../build/release

                make > /dev/null

                ./non_imv_indexed_nested_loop_join_runner

                cd ../../scripts/evaluation/
            done
        done

    fi
}

run_nums=10 #1 10
load_relations_for_evaluation=1 #0
persist_relations_for_evaluation=0 #1

#unique datasets
################

r_datasets=(r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v5_uint32_uint32_640000000) #(r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v4_uint32_uint32_384000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v6_uint32_uint32_896000000 r_UNIQUE_v7_uint32_uint32_1152000000 r_UNIQUE_v8_uint32_uint32_1664000000 r_UNIQUE_v9_uint32_uint32_1920000000)
s_datasets=(s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000) #(s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v4_uint32_uint32_384000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v6_uint32_uint32_896000000 s_UNIQUE_v7_uint32_uint32_1152000000 s_UNIQUE_v8_uint32_uint32_1664000000 s_UNIQUE_v9_uint32_uint32_1920000000)
r_datasets_sizes=(16E6 16E6 16E6 32E6 32E6 32E6 128E6 128E6 128E6 640E6 640E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(32E6 128E6 640E6 16E6 128E6 640E6 16E6 32E6 640E6 16E6 32E6 128E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32 32 32 32 32 32 32 32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32 32 32 32 32 32 32 32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)
input_hash_table_size=(16777216 16777216 16777216 33554432 33554432 33554432 134217728 134217728 134217728 536870912 536870912 536870912) #(16777216(for_16E6) 33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_unique/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_learned_index_unique/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 1 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_csstree_index_unique/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 1 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_art32tree_index_unique/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 0 1 $input_hash_table_size

r_datasets=(r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v8_uint32_uint32_1664000000) #(r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v4_uint32_uint32_384000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v6_uint32_uint32_896000000 r_UNIQUE_v7_uint32_uint32_1152000000 r_UNIQUE_v8_uint32_uint32_1664000000 r_UNIQUE_v9_uint32_uint32_1920000000)
s_datasets=(s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v8_uint32_uint32_1664000000) #(s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v4_uint32_uint32_384000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v6_uint32_uint32_896000000 s_UNIQUE_v7_uint32_uint32_1152000000 s_UNIQUE_v8_uint32_uint32_1664000000 s_UNIQUE_v9_uint32_uint32_1920000000)
r_datasets_sizes=(16E6 32E6 128E6 640E6 1664E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(16E6 32E6 128E6 640E6 1664E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)
input_hash_table_size=(16777216 33554432 134217728 536870912 1073741824) #(33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_unique/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_learned_index_unique/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 1 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_csstree_index_unique/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 1 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_art32tree_index_unique/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 0 1 $input_hash_table_size


#lognormal datasets
###################

r_datasets=(r_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 r_LOGNORMAL_v2_uint32_uint32_segma_1_32000000  r_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 r_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 r_LOGNORMAL_v5_uint32_uint32_segma_1_640000000 r_LOGNORMAL_v5_uint32_uint32_segma_1_640000000) #(r_LOGNORMAL_v1_uint32_uint32_segma_1_16000000 r_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 r_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 r_LOGNORMAL_v4_uint32_uint32_segma_1_384000000 r_LOGNORMAL_v5_uint32_uint32_segma_1_640000000 r_LOGNORMAL_v6_uint32_uint32_segma_1_896000000 r_LOGNORMAL_v7_uint32_uint32_segma_1_1152000000 r_LOGNORMAL_v8_uint32_uint32_segma_1_1664000000 r_LOGNORMAL_v9_uint32_uint32_segma_1_1920000000)
s_datasets=(s_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 s_LOGNORMAL_v5_uint32_uint32_segma_1_640000000 s_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 s_LOGNORMAL_v5_uint32_uint32_segma_1_640000000 s_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 s_LOGNORMAL_v3_uint32_uint32_segma_1_128000000) #(s_LOGNORMAL_v1_uint32_uint32_segma_1_16000000 s_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 s_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 s_LOGNORMAL_v4_uint32_uint32_segma_1_384000000 s_LOGNORMAL_v5_uint32_uint32_segma_1_640000000 s_LOGNORMAL_v6_uint32_uint32_segma_1_896000000 s_LOGNORMAL_v7_uint32_uint32_segma_1_1152000000 s_LOGNORMAL_v8_uint32_uint32_segma_1_1664000000 s_LOGNORMAL_v9_uint32_uint32_segma_1_1920000000)
r_datasets_sizes=(32E6 32E6 128E6 128E6 640E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(128E6 640E6 32E6 640E6 32E6 128E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32 32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32 32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)
input_hash_table_size=(33554432 33554432 134217728 134217728 536870912 536870912) #(16777216(for_16E6) 33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_lognormal/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_learned_index_lognormal/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 1 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_csstree_index_lognormal/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 1 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_art32tree_index_lognormal/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 0 1 $input_hash_table_size


r_datasets=(r_LOGNORMAL_v1_uint32_uint32_segma_1_16000000 r_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 r_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 r_LOGNORMAL_v5_uint32_uint32_segma_1_640000000) #(r_LOGNORMAL_v1_uint32_uint32_segma_1_16000000 r_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 r_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 r_LOGNORMAL_v4_uint32_uint32_segma_1_384000000 r_LOGNORMAL_v5_uint32_uint32_segma_1_640000000 r_LOGNORMAL_v6_uint32_uint32_segma_1_896000000 r_LOGNORMAL_v7_uint32_uint32_segma_1_1152000000 r_LOGNORMAL_v8_uint32_uint32_segma_1_1664000000 r_LOGNORMAL_v9_uint32_uint32_segma_1_1920000000)
s_datasets=(s_LOGNORMAL_v1_uint32_uint32_segma_1_16000000 s_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 s_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 s_LOGNORMAL_v5_uint32_uint32_segma_1_640000000) #(s_LOGNORMAL_v1_uint32_uint32_segma_1_16000000 s_LOGNORMAL_v2_uint32_uint32_segma_1_32000000 s_LOGNORMAL_v3_uint32_uint32_segma_1_128000000 s_LOGNORMAL_v4_uint32_uint32_segma_1_384000000 s_LOGNORMAL_v5_uint32_uint32_segma_1_640000000 s_LOGNORMAL_v6_uint32_uint32_segma_1_896000000 s_LOGNORMAL_v7_uint32_uint32_segma_1_1152000000 s_LOGNORMAL_v8_uint32_uint32_segma_1_1664000000 s_LOGNORMAL_v9_uint32_uint32_segma_1_1920000000)
r_datasets_sizes=(16E6 32E6 128E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(16E6 32E6 128E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32 32 32 32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32 32 32 32) #(64 64 64 64 64 64 64 64 64)
input_hash_table_size=(16777216 33554432 134217728 536870912) #(33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_lognormal/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_learned_index_lognormal/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 1 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_csstree_index_lognormal/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 1 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_art32tree_index_lognormal/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 0 1 $input_hash_table_size

#seq_hole datasets
################

r_datasets=(r_SEQ_HOLE_v1_uint32_uint32_16000000 r_SEQ_HOLE_v1_uint32_uint32_16000000 r_SEQ_HOLE_v1_uint32_uint32_16000000 r_SEQ_HOLE_v2_uint32_uint32_32000000 r_SEQ_HOLE_v2_uint32_uint32_32000000 r_SEQ_HOLE_v2_uint32_uint32_32000000 r_SEQ_HOLE_v3_uint32_uint32_128000000 r_SEQ_HOLE_v3_uint32_uint32_128000000 r_SEQ_HOLE_v3_uint32_uint32_128000000 r_SEQ_HOLE_v5_uint32_uint32_640000000 r_SEQ_HOLE_v5_uint32_uint32_640000000 r_SEQ_HOLE_v5_uint32_uint32_640000000) #(r_SEQ_HOLE_v1_uint32_uint32_16000000 r_SEQ_HOLE_v2_uint32_uint32_32000000 r_SEQ_HOLE_v3_uint32_uint32_128000000 r_SEQ_HOLE_v4_uint32_uint32_384000000 r_SEQ_HOLE_v5_uint32_uint32_640000000 r_SEQ_HOLE_v6_uint32_uint32_896000000 r_SEQ_HOLE_v7_uint32_uint32_1152000000 r_SEQ_HOLE_v8_uint32_uint32_1664000000 r_SEQ_HOLE_v9_uint32_uint32_1920000000)
s_datasets=(s_SEQ_HOLE_v2_uint32_uint32_32000000 s_SEQ_HOLE_v3_uint32_uint32_128000000 s_SEQ_HOLE_v5_uint32_uint32_640000000 s_SEQ_HOLE_v1_uint32_uint32_16000000 s_SEQ_HOLE_v3_uint32_uint32_128000000 s_SEQ_HOLE_v5_uint32_uint32_640000000 s_SEQ_HOLE_v1_uint32_uint32_16000000 s_SEQ_HOLE_v2_uint32_uint32_32000000 s_SEQ_HOLE_v5_uint32_uint32_640000000 s_SEQ_HOLE_v1_uint32_uint32_16000000 s_SEQ_HOLE_v2_uint32_uint32_32000000 s_SEQ_HOLE_v3_uint32_uint32_128000000) #(s_SEQ_HOLE_v1_uint32_uint32_16000000 s_SEQ_HOLE_v2_uint32_uint32_32000000 s_SEQ_HOLE_v3_uint32_uint32_128000000 s_SEQ_HOLE_v4_uint32_uint32_384000000 s_SEQ_HOLE_v5_uint32_uint32_640000000 s_SEQ_HOLE_v6_uint32_uint32_896000000 s_SEQ_HOLE_v7_uint32_uint32_1152000000 s_SEQ_HOLE_v8_uint32_uint32_1664000000 s_SEQ_HOLE_v9_uint32_uint32_1920000000)
r_datasets_sizes=(16E6 16E6 16E6 32E6 32E6 32E6 128E6 128E6 128E6 640E6 640E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(32E6 128E6 640E6 16E6 128E6 640E6 16E6 32E6 640E6 16E6 32E6 128E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32 32 32 32 32 32 32 32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32 32 32 32 32 32 32 32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)
input_hash_table_size=(16777216 16777216 16777216 33554432 33554432 33554432 134217728 134217728 134217728 536870912 536870912 536870912) #(16777216(for_16E6) 33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_seq_hole/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_learned_index_seq_hole/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 1 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_csstree_index_seq_hole/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 1 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_art32tree_index_seq_hole/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 0 1 $input_hash_table_size


r_datasets=(r_SEQ_HOLE_v1_uint32_uint32_16000000 r_SEQ_HOLE_v2_uint32_uint32_32000000 r_SEQ_HOLE_v3_uint32_uint32_128000000 r_SEQ_HOLE_v5_uint32_uint32_640000000) #(r_SEQ_HOLE_v1_uint32_uint32_16000000 r_SEQ_HOLE_v2_uint32_uint32_32000000 r_SEQ_HOLE_v3_uint32_uint32_128000000 r_SEQ_HOLE_v4_uint32_uint32_384000000 r_SEQ_HOLE_v5_uint32_uint32_640000000 r_SEQ_HOLE_v6_uint32_uint32_896000000 r_SEQ_HOLE_v7_uint32_uint32_1152000000 r_SEQ_HOLE_v8_uint32_uint32_1664000000 r_SEQ_HOLE_v9_uint32_uint32_1920000000)
s_datasets=(s_SEQ_HOLE_v1_uint32_uint32_16000000 s_SEQ_HOLE_v2_uint32_uint32_32000000 s_SEQ_HOLE_v3_uint32_uint32_128000000 s_SEQ_HOLE_v5_uint32_uint32_640000000) #(s_SEQ_HOLE_v1_uint32_uint32_16000000 s_SEQ_HOLE_v2_uint32_uint32_32000000 s_SEQ_HOLE_v3_uint32_uint32_128000000 s_SEQ_HOLE_v4_uint32_uint32_384000000 s_SEQ_HOLE_v5_uint32_uint32_640000000 s_SEQ_HOLE_v6_uint32_uint32_896000000 s_SEQ_HOLE_v7_uint32_uint32_1152000000 s_SEQ_HOLE_v8_uint32_uint32_1664000000 s_SEQ_HOLE_v9_uint32_uint32_1920000000)
r_datasets_sizes=(16E6 32E6 128E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(16E6 32E6 128E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32 32 32 32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32 32 32 32) #(64 64 64 64 64 64 64 64 64)
input_hash_table_size=(16777216 33554432 134217728 536870912) #(33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_seq_hole/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_learned_index_seq_hole/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 1 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_csstree_index_seq_hole/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 1 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_art32tree_index_seq_hole/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 0 1 $input_hash_table_size


#uniform datasets
################

r_datasets=(r_UNIFORM_v2_uint32_uint32_32000000 r_UNIFORM_v2_uint32_uint32_32000000  r_UNIFORM_v3_uint32_uint32_128000000 r_UNIFORM_v3_uint32_uint32_128000000 r_UNIFORM_v5_uint32_uint32_640000000 r_UNIFORM_v5_uint32_uint32_640000000) #(r_UNIFORM_v1_uint32_uint32_16000000 r_UNIFORM_v2_uint32_uint32_32000000 r_UNIFORM_v3_uint32_uint32_128000000 r_UNIFORM_v4_uint32_uint32_384000000 r_UNIFORM_v5_uint32_uint32_640000000 r_UNIFORM_v6_uint32_uint32_896000000 r_UNIFORM_v7_uint32_uint32_1152000000 r_UNIFORM_v8_uint32_uint32_1664000000 r_UNIFORM_v9_uint32_uint32_1920000000)
s_datasets=(s_UNIFORM_v3_uint32_uint32_128000000 s_UNIFORM_v5_uint32_uint32_640000000 s_UNIFORM_v2_uint32_uint32_32000000 s_UNIFORM_v5_uint32_uint32_640000000 s_UNIFORM_v2_uint32_uint32_32000000 s_UNIFORM_v3_uint32_uint32_128000000) #(s_UNIFORM_v1_uint32_uint32_16000000 s_UNIFORM_v2_uint32_uint32_32000000 s_UNIFORM_v3_uint32_uint32_128000000 s_UNIFORM_v4_uint32_uint32_384000000 s_UNIFORM_v5_uint32_uint32_640000000 s_UNIFORM_v6_uint32_uint32_896000000 s_UNIFORM_v7_uint32_uint32_1152000000 s_UNIFORM_v8_uint32_uint32_1664000000 s_UNIFORM_v9_uint32_uint32_1920000000)
r_datasets_sizes=(32E6 32E6 128E6 128E6 640E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(128E6 640E6 32E6 640E6 32E6 128E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32 32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32 32 32 32 32 32) #(64 64 64 64 64 64 64 64 64)
input_hash_table_size=(33554432 33554432 134217728 134217728 536870912 536870912) #(16777216(for_16E6) 33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_uniform/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_learned_index_uniform/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 1 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_csstree_index_uniform/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 1 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_art32tree_index_uniform/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 0 1 $input_hash_table_size


r_datasets=(r_UNIFORM_v1_uint32_uint32_16000000 r_UNIFORM_v2_uint32_uint32_32000000 r_UNIFORM_v3_uint32_uint32_128000000 r_UNIFORM_v5_uint32_uint32_640000000) #(r_UNIFORM_v1_uint32_uint32_16000000 r_UNIFORM_v2_uint32_uint32_32000000 r_UNIFORM_v3_uint32_uint32_128000000 r_UNIFORM_v4_uint32_uint32_384000000 r_UNIFORM_v5_uint32_uint32_640000000 r_UNIFORM_v6_uint32_uint32_896000000 r_UNIFORM_v7_uint32_uint32_1152000000 r_UNIFORM_v8_uint32_uint32_1664000000 r_UNIFORM_v9_uint32_uint32_1920000000)
s_datasets=(s_UNIFORM_v1_uint32_uint32_16000000 s_UNIFORM_v2_uint32_uint32_32000000 s_UNIFORM_v3_uint32_uint32_128000000 s_UNIFORM_v5_uint32_uint32_640000000) #(s_UNIFORM_v1_uint32_uint32_16000000 s_UNIFORM_v2_uint32_uint32_32000000 s_UNIFORM_v3_uint32_uint32_128000000 s_UNIFORM_v4_uint32_uint32_384000000 s_UNIFORM_v5_uint32_uint32_640000000 s_UNIFORM_v6_uint32_uint32_896000000 s_UNIFORM_v7_uint32_uint32_1152000000 s_UNIFORM_v8_uint32_uint32_1664000000 s_UNIFORM_v9_uint32_uint32_1920000000)
r_datasets_sizes=(16E6 32E6 128E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_datasets_sizes=(16E6 32E6 128E6 640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_datasets_file_num_partitions=(32 32 32 32) #(64 64 64 64 64 64 64 64 64)
s_datasets_file_num_partitions=(32 32 32 32) #(64 64 64 64 64 64 64 64 64)
input_hash_table_size=(16777216 33554432 134217728 536870912) #(33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_uniform/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_learned_index_uniform/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 1 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_csstree_index_uniform/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 1 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_art32tree_index_uniform/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 0 1 $input_hash_table_size

#sosd datasets
################

r_datasets=(books_200M_uint32) #(books_200M_uint32 books_800M_uint64 fb_200M_uint64 osm_cellids_800M_uint64 wiki_ts_200M_uint64) 
s_datasets=(books_200M_uint32) #(books_200M_uint32 books_800M_uint64 fb_200M_uint64 osm_cellids_800M_uint64 wiki_ts_200M_uint64)
r_datasets_sizes=(200E6) #(200E6 800E6 200E6 800E6 200E6)
s_datasets_sizes=(200E6) #(200E6 800E6 200E6 800E6 200E6)
r_datasets_file_num_partitions=(32) #(32 32 32 32 32)
s_datasets_file_num_partitions=(32) #(32 32 32 32 32)
input_hash_table_size=(536870912) #(16777216(for_16E6) 33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_sosd_books_200M_uint32/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_learned_index_sosd_books_200M_uint32/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 1 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_csstree_index_sosd_books_200M_uint32/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 1 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_art32tree_index_sosd_books_200M_uint32/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 0 1 $input_hash_table_size

r_datasets=(books_800M_uint64) #(books_200M_uint32 books_800M_uint64 fb_200M_uint64 osm_cellids_800M_uint64 wiki_ts_200M_uint64) 
s_datasets=(books_800M_uint64) #(books_200M_uint32 books_800M_uint64 fb_200M_uint64 osm_cellids_800M_uint64 wiki_ts_200M_uint64)
r_datasets_sizes=(800E6) #(200E6 800E6 200E6 800E6 200E6)
s_datasets_sizes=(800E6) #(200E6 800E6 200E6 800E6 200E6)
r_datasets_file_num_partitions=(32) #(32 32 32 32 32)
s_datasets_file_num_partitions=(32) #(32 32 32 32 32)
input_hash_table_size=(1073741824) #(16777216(for_16E6) 33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_sosd_books_800M_uint64/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_learned_index_sosd_books_800M_uint64/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 1 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_csstree_index_sosd_books_800M_uint64/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 1 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_art32tree_index_sosd_books_800M_uint64/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 0 1 $input_hash_table_size


r_datasets=(fb_200M_uint64) #(books_200M_uint32 books_800M_uint64 fb_200M_uint64 osm_cellids_800M_uint64 wiki_ts_200M_uint64) 
s_datasets=(fb_200M_uint64) #(books_200M_uint32 books_800M_uint64 fb_200M_uint64 osm_cellids_800M_uint64 wiki_ts_200M_uint64)
r_datasets_sizes=(200E6) #(200E6 800E6 200E6 800E6 200E6)
s_datasets_sizes=(200E6) #(200E6 800E6 200E6 800E6 200E6)
r_datasets_file_num_partitions=(32) #(32 32 32 32 32)
s_datasets_file_num_partitions=(32) #(32 32 32 32 32)
input_hash_table_size=(536870912) #(16777216(for_16E6) 33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_sosd_fb_200M_uint64/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_learned_index_sosd_fb_200M_uint64/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 1 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_csstree_index_sosd_fb_200M_uint64/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 1 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_art32tree_index_sosd_fb_200M_uint64/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 0 1 $input_hash_table_size


r_datasets=(osm_cellids_800M_uint64) #(books_200M_uint32 books_800M_uint64 fb_200M_uint64 osm_cellids_800M_uint64 wiki_ts_200M_uint64) 
s_datasets=(osm_cellids_800M_uint64) #(books_200M_uint32 books_800M_uint64 fb_200M_uint64 osm_cellids_800M_uint64 wiki_ts_200M_uint64)
r_datasets_sizes=(800E6) #(200E6 800E6 200E6 800E6 200E6)
s_datasets_sizes=(800E6) #(200E6 800E6 200E6 800E6 200E6)
r_datasets_file_num_partitions=(32) #(32 32 32 32 32)
s_datasets_file_num_partitions=(32) #(32 32 32 32 32)
input_hash_table_size=(1073741824) #(16777216(for_16E6) 33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_sosd_osm_cellids_800M_uint64/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_learned_index_sosd_osm_cellids_800M_uint64/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 1 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_csstree_index_sosd_osm_cellids_800M_uint64/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 1 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_art32tree_index_sosd_osm_cellids_800M_uint64/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 0 1 $input_hash_table_size

r_datasets=(wiki_ts_200M_uint64) #(books_200M_uint32 books_800M_uint64 fb_200M_uint64 osm_cellids_800M_uint64 wiki_ts_200M_uint64) 
s_datasets=(wiki_ts_200M_uint64) #(books_200M_uint32 books_800M_uint64 fb_200M_uint64 osm_cellids_800M_uint64 wiki_ts_200M_uint64)
r_datasets_sizes=(200E6) #(200E6 800E6 200E6 800E6 200E6)
s_datasets_sizes=(200E6) #(200E6 800E6 200E6 800E6 200E6)
r_datasets_file_num_partitions=(32) #(32 32 32 32 32)
s_datasets_file_num_partitions=(32) #(32 32 32 32 32)
input_hash_table_size=(536870912) #(16777216(for_16E6) 33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_sosd_wiki_ts_200M_uint64/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_learned_index_sosd_wiki_ts_200M_uint64/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 1 0 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_csstree_index_sosd_wiki_ts_200M_uint64/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 1 0 $input_hash_table_size
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_art32tree_index_sosd_wiki_ts_200M_uint64/
#process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 0 1 $input_hash_table_size


#tpch workloads
################

r_datasets=(r_q18_v1_uint32_uint32_57 r_q18_v2_uint32_uint32_150000 r_q18_v3_uint32_uint32_57)
s_datasets=(s_q18_v1_uint32_uint32_1500000 s_q18_v2_uint32_uint32_57 s_q18_v3_uint32_uint32_6001215) 
r_datasets_sizes=(57 150000 57) 
s_datasets_sizes=(1500000 57 6001215)
r_datasets_file_num_partitions=(32 32 32) 
s_datasets_file_num_partitions=(32 32 32)
input_hash_table_size=(16777216 16777216 16777216) #(16777216(for_16E6) 33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_tpch_q18/
process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0 0 0 $input_hash_table_size
output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_learned_index_tpch_q18/
process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 1 0 0 $input_hash_table_size
output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_csstree_index_tpch_q18/
process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 1 0 $input_hash_table_size
output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_art32tree_index_tpch_q18/
process_non_imv_indexed_nested_loop_join $r_datasets $r_datasets_sizes $r_datasets_file_num_partitions $s_datasets $s_datasets_sizes $s_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 0 1 $input_hash_table_size



#-rw-r--r-- 1 sabek sabek 3.7K May 30 18:22 r_q3_v1_uint32_uint32_30142_sorted_0.txt
#-rw-r--r-- 1 sabek sabek  18K May 30 18:47 r_q3_v2_uint32_uint32_147126_sorted_0.txt
#-rw-r--r-- 1 sabek sabek   10 May 30 23:51 r_q5_v1_uint32_uint32_1_sorted_0.txt
#-rw-r--r-- 1 sabek sabek   10 May 31 01:38 r_q5_v2_uint32_uint32_5_sorted_0.txt
#-rw-r--r-- 1 sabek sabek 3.7K May 31 00:22 r_q5_v3_uint32_uint32_30183_sorted_0.txt
#-rw-r--r-- 1 sabek sabek 5.7K May 31 00:38 r_q5_v4_uint32_uint32_46008_sorted_0.txt
#-rw-r--r-- 1 sabek sabek  23K May 31 01:00 r_q5_v5_uint32_uint32_184082_sorted_0.txt
#-rw-r--r-- 1 sabek sabek   10 May 31 12:29 r_q9_v1_uint32_uint32_25_sorted_0.txt
#-rw-r--r-- 1 sabek sabek 1.4K May 31 12:41 r_q9_v2_uint32_uint32_10664_sorted_0.txt
#-rw-r--r-- 1 sabek sabek 1.3K May 31 12:47 r_q9_v3_uint32_uint32_10000_sorted_0.txt
#-rw-r--r-- 1 sabek sabek 5.3K May 31 12:53 r_q9_v4_uint32_uint32_42656_sorted_0.txt

#-rw-r--r-- 1 sabek sabek  89K May 31 01:55 s_q3_v1_uint32_uint32_727305_sorted_0.txt
#-rw-r--r-- 1 sabek sabek 396K May 30 19:04 s_q3_v2_uint32_uint32_3241776_sorted_0.txt
#-rw-r--r-- 1 sabek sabek   10 May 31 01:33 s_q5_v1_uint32_uint32_25_sorted_0.txt
#-rw-r--r-- 1 sabek sabek  19K May 31 00:09 s_q5_v2_uint32_uint32_150000_sorted_0.txt
#-rw-r--r-- 1 sabek sabek  28K May 31 00:30 s_q5_v3_uint32_uint32_227597_sorted_0.txt
#-rw-r--r-- 1 sabek sabek 733K May 31 00:51 s_q5_v4_uint32_uint32_6001215_sorted_0.txt
#-rw-r--r-- 1 sabek sabek 1.3K May 31 01:11 s_q5_v5_uint32_uint32_10000_sorted_0.txt
#-rw-r--r-- 1 sabek sabek 1.3K May 31 12:37 s_q9_v1_uint32_uint32_10000_sorted_0.txt
#-rw-r--r-- 1 sabek sabek  98K May 31 12:44 s_q9_v2_uint32_uint32_800000_sorted_0.txt
#-rw-r--r-- 1 sabek sabek 5.3K May 31 12:49 s_q9_v3_uint32_uint32_42656_sorted_0.txt
#-rw-r--r-- 1 sabek sabek 733K May 31 12:57 s_q9_v4_uint32_uint32_6001215_sorted_0.txt
#!/bin/sh

process_non_imv_indexed_nested_loop_join()
{
    threads=(4 16 32 64) #(2 4 8 16 32 64)
    rmi_models=(0 1 2 3 4 5 6 7 8 9)
    css_fanouts=(10 33 50)

    dataset_folder_path=/spinning/sabek/learned_join_datasets/

    r_datasets=$1
    r_datasets_sizes=$2
    r_datasets_file_num_partitions=$3
    r_datasets_file_extension='"'.txt'"'
    s_datasets=$4
    s_datasets_sizes=$5
    s_datasets_file_num_partitions=$6
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

    input_hash_table_size=${15}

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
												-RUN_NUMS ${run_nums} \
												-LOAD_RELATIONS_FOR_EVALUATION ${load_relations_for_evaluation} \
												-PERSIST_RELATIONS_FOR_EVALUATION ${persist_relations_for_evaluation}

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
                    curr_rmi_model=${r_datasets[$ds]}'_key_uint32_'${rmi_models[$model]}

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
                                                    -RUN_NUMS ${run_nums} \
                                                    -LOAD_RELATIONS_FOR_EVALUATION ${load_relations_for_evaluation} \
                                                    -PERSIST_RELATIONS_FOR_EVALUATION ${persist_relations_for_evaluation}

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
												-RUN_NUMS ${run_nums} \
												-LOAD_RELATIONS_FOR_EVALUATION ${load_relations_for_evaluation} \
												-PERSIST_RELATIONS_FOR_EVALUATION ${persist_relations_for_evaluation}

                done
                
                cmake -DCMAKE_BUILD_TYPE=Release -DVECTORWISE_BRANCHING=on $(dirname "$0")/../.. > /dev/null

                cd $(dirname "$0")/../../build/release

                make > /dev/null

                ./non_imv_indexed_nested_loop_join_runner

                cd ../../scripts/evaluation/
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
                                            -RUN_NUMS ${run_nums} \
                                            -LOAD_RELATIONS_FOR_EVALUATION ${load_relations_for_evaluation} \
                                            -PERSIST_RELATIONS_FOR_EVALUATION ${persist_relations_for_evaluation}

                
                cmake -DCMAKE_BUILD_TYPE=Release -DVECTORWISE_BRANCHING=on $(dirname "$0")/../.. > /dev/null

                cd $(dirname "$0")/../../build/release

                make > /dev/null

                ./non_imv_indexed_nested_loop_join_runner

                cd ../../scripts/evaluation/
            done
        done

    fi
}

run_nums=10
load_relations_for_evaluation=1
persist_relations_for_evaluation=0

#unique datasets
################
r_unique_datasets=(r_UNIQUE_v2_uint32_uint32_32000000) #(r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v4_uint32_uint32_384000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v6_uint32_uint32_896000000 r_UNIQUE_v7_uint32_uint32_1152000000 r_UNIQUE_v8_uint32_uint32_1664000000 r_UNIQUE_v9_uint32_uint32_1920000000)
s_unique_datasets=(s_UNIQUE_v2_uint32_uint32_32000000) #(s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v4_uint32_uint32_384000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v6_uint32_uint32_896000000 s_UNIQUE_v7_uint32_uint32_1152000000 s_UNIQUE_v8_uint32_uint32_1664000000 s_UNIQUE_v9_uint32_uint32_1920000000)
r_unique_datasets_sizes=(32E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_unique_datasets_sizes=(32E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_unique_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
s_unique_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_unique/
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_learned_index_unique/
output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_csstree_index_unique/
input_hash_table_size=(33554432) #(33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

process_non_imv_indexed_nested_loop_join $r_unique_datasets $r_unique_datasets_sizes $r_unique_datasets_file_num_partitions $s_unique_datasets $s_unique_datasets_sizes $s_unique_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 1 0 $input_hash_table_size

r_unique_datasets=(r_UNIQUE_v3_uint32_uint32_128000000) #(r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v4_uint32_uint32_384000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v6_uint32_uint32_896000000 r_UNIQUE_v7_uint32_uint32_1152000000 r_UNIQUE_v8_uint32_uint32_1664000000 r_UNIQUE_v9_uint32_uint32_1920000000)
s_unique_datasets=(s_UNIQUE_v3_uint32_uint32_128000000) #(s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v4_uint32_uint32_384000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v6_uint32_uint32_896000000 s_UNIQUE_v7_uint32_uint32_1152000000 s_UNIQUE_v8_uint32_uint32_1664000000 s_UNIQUE_v9_uint32_uint32_1920000000)
r_unique_datasets_sizes=(128E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_unique_datasets_sizes=(128E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_unique_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
s_unique_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_unique/
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_learned_index_unique/
output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_csstree_index_unique/
input_hash_table_size=(134217728) #(33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

process_non_imv_indexed_nested_loop_join $r_unique_datasets $r_unique_datasets_sizes $r_unique_datasets_file_num_partitions $s_unique_datasets $s_unique_datasets_sizes $s_unique_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 1 0 $input_hash_table_size

r_unique_datasets=(r_UNIQUE_v5_uint32_uint32_640000000) #(r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v4_uint32_uint32_384000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v6_uint32_uint32_896000000 r_UNIQUE_v7_uint32_uint32_1152000000 r_UNIQUE_v8_uint32_uint32_1664000000 r_UNIQUE_v9_uint32_uint32_1920000000)
s_unique_datasets=(s_UNIQUE_v5_uint32_uint32_640000000) #(s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v4_uint32_uint32_384000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v6_uint32_uint32_896000000 s_UNIQUE_v7_uint32_uint32_1152000000 s_UNIQUE_v8_uint32_uint32_1664000000 s_UNIQUE_v9_uint32_uint32_1920000000)
r_unique_datasets_sizes=(640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_unique_datasets_sizes=(640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_unique_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
s_unique_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_unique/
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_learned_index_unique/
output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_csstree_index_unique/
input_hash_table_size=(536870912) #(33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

process_non_imv_indexed_nested_loop_join $r_unique_datasets $r_unique_datasets_sizes $r_unique_datasets_file_num_partitions $s_unique_datasets $s_unique_datasets_sizes $s_unique_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 1 0 $input_hash_table_size


r_unique_datasets=(r_UNIQUE_v8_uint32_uint32_1664000000) #(r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v4_uint32_uint32_384000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v6_uint32_uint32_896000000 r_UNIQUE_v7_uint32_uint32_1152000000 r_UNIQUE_v8_uint32_uint32_1664000000 r_UNIQUE_v9_uint32_uint32_1920000000)
s_unique_datasets=(s_UNIQUE_v8_uint32_uint32_1664000000) #(s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v4_uint32_uint32_384000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v6_uint32_uint32_896000000 s_UNIQUE_v7_uint32_uint32_1152000000 s_UNIQUE_v8_uint32_uint32_1664000000 s_UNIQUE_v9_uint32_uint32_1920000000)
r_unique_datasets_sizes=(1664E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_unique_datasets_sizes=(1664E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_unique_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
s_unique_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_unique/
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_learned_index_unique/
output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_csstree_index_unique/
input_hash_table_size=(1073741824) #(134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

process_non_imv_indexed_nested_loop_join $r_unique_datasets $r_unique_datasets_sizes $r_unique_datasets_file_num_partitions $s_unique_datasets $s_unique_datasets_sizes $s_unique_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 1 0 $input_hash_table_size




r_unique_datasets=(r_UNIQUE_v2_uint32_uint32_32000000) #(r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v4_uint32_uint32_384000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v6_uint32_uint32_896000000 r_UNIQUE_v7_uint32_uint32_1152000000 r_UNIQUE_v8_uint32_uint32_1664000000 r_UNIQUE_v9_uint32_uint32_1920000000)
s_unique_datasets=(s_UNIQUE_v2_uint32_uint32_32000000) #(s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v4_uint32_uint32_384000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v6_uint32_uint32_896000000 s_UNIQUE_v7_uint32_uint32_1152000000 s_UNIQUE_v8_uint32_uint32_1664000000 s_UNIQUE_v9_uint32_uint32_1920000000)
r_unique_datasets_sizes=(32E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_unique_datasets_sizes=(32E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_unique_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
s_unique_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_unique/
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_learned_index_unique/
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_csstree_index_unique/
output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_art32tree_index_unique/
input_hash_table_size=(33554432) #(33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

process_non_imv_indexed_nested_loop_join $r_unique_datasets $r_unique_datasets_sizes $r_unique_datasets_file_num_partitions $s_unique_datasets $s_unique_datasets_sizes $s_unique_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 0 1 $input_hash_table_size

r_unique_datasets=(r_UNIQUE_v3_uint32_uint32_128000000) #(r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v4_uint32_uint32_384000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v6_uint32_uint32_896000000 r_UNIQUE_v7_uint32_uint32_1152000000 r_UNIQUE_v8_uint32_uint32_1664000000 r_UNIQUE_v9_uint32_uint32_1920000000)
s_unique_datasets=(s_UNIQUE_v3_uint32_uint32_128000000) #(s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v4_uint32_uint32_384000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v6_uint32_uint32_896000000 s_UNIQUE_v7_uint32_uint32_1152000000 s_UNIQUE_v8_uint32_uint32_1664000000 s_UNIQUE_v9_uint32_uint32_1920000000)
r_unique_datasets_sizes=(128E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_unique_datasets_sizes=(128E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_unique_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
s_unique_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_unique/
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_learned_index_unique/
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_csstree_index_unique/
output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_art32tree_index_unique/
input_hash_table_size=(134217728) #(33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

process_non_imv_indexed_nested_loop_join $r_unique_datasets $r_unique_datasets_sizes $r_unique_datasets_file_num_partitions $s_unique_datasets $s_unique_datasets_sizes $s_unique_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 0 1 $input_hash_table_size

r_unique_datasets=(r_UNIQUE_v5_uint32_uint32_640000000) #(r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v4_uint32_uint32_384000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v6_uint32_uint32_896000000 r_UNIQUE_v7_uint32_uint32_1152000000 r_UNIQUE_v8_uint32_uint32_1664000000 r_UNIQUE_v9_uint32_uint32_1920000000)
s_unique_datasets=(s_UNIQUE_v5_uint32_uint32_640000000) #(s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v4_uint32_uint32_384000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v6_uint32_uint32_896000000 s_UNIQUE_v7_uint32_uint32_1152000000 s_UNIQUE_v8_uint32_uint32_1664000000 s_UNIQUE_v9_uint32_uint32_1920000000)
r_unique_datasets_sizes=(640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_unique_datasets_sizes=(640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_unique_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
s_unique_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_unique/
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_learned_index_unique/
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_csstree_index_unique/
output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_art32tree_index_unique/
input_hash_table_size=(536870912) #(33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

process_non_imv_indexed_nested_loop_join $r_unique_datasets $r_unique_datasets_sizes $r_unique_datasets_file_num_partitions $s_unique_datasets $s_unique_datasets_sizes $s_unique_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 0 1 $input_hash_table_size


r_unique_datasets=(r_UNIQUE_v8_uint32_uint32_1664000000) #(r_UNIQUE_v1_uint32_uint32_16000000 r_UNIQUE_v2_uint32_uint32_32000000 r_UNIQUE_v3_uint32_uint32_128000000 r_UNIQUE_v4_uint32_uint32_384000000 r_UNIQUE_v5_uint32_uint32_640000000 r_UNIQUE_v6_uint32_uint32_896000000 r_UNIQUE_v7_uint32_uint32_1152000000 r_UNIQUE_v8_uint32_uint32_1664000000 r_UNIQUE_v9_uint32_uint32_1920000000)
s_unique_datasets=(s_UNIQUE_v8_uint32_uint32_1664000000) #(s_UNIQUE_v1_uint32_uint32_16000000 s_UNIQUE_v2_uint32_uint32_32000000 s_UNIQUE_v3_uint32_uint32_128000000 s_UNIQUE_v4_uint32_uint32_384000000 s_UNIQUE_v5_uint32_uint32_640000000 s_UNIQUE_v6_uint32_uint32_896000000 s_UNIQUE_v7_uint32_uint32_1152000000 s_UNIQUE_v8_uint32_uint32_1664000000 s_UNIQUE_v9_uint32_uint32_1920000000)
r_unique_datasets_sizes=(1664E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
s_unique_datasets_sizes=(1664E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
r_unique_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
s_unique_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_unique/
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_learned_index_unique/
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_csstree_index_unique/
output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_art32tree_index_unique/
input_hash_table_size=(1073741824) #(134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

process_non_imv_indexed_nested_loop_join $r_unique_datasets $r_unique_datasets_sizes $r_unique_datasets_file_num_partitions $s_unique_datasets $s_unique_datasets_sizes $s_unique_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 0 0 0 1 $input_hash_table_size


#lognormal datasets
###################
#r_lognormal_datasets=(r_LOGNORMAL_v2_uint32_uint32_segma_1.000000_32000000) #(r_LOGNORMAL_v1_uint32_uint32_segma_1.000000_16000000 r_LOGNORMAL_v2_uint32_uint32_segma_1.000000_32000000 r_LOGNORMAL_v3_uint32_uint32_segma_1.000000_128000000 r_LOGNORMAL_v4_uint32_uint32_segma_1.000000_384000000 r_LOGNORMAL_v5_uint32_uint32_segma_1.000000_640000000 r_LOGNORMAL_v6_uint32_uint32_segma_1.000000_896000000 r_LOGNORMAL_v7_uint32_uint32_segma_1.000000_1152000000 r_LOGNORMAL_v8_uint32_uint32_segma_1.000000_1664000000 r_LOGNORMAL_v9_uint32_uint32_segma_1.000000_1920000000)
#s_lognormal_datasets=(s_LOGNORMAL_v2_uint32_uint32_segma_1.000000_32000000) #(s_LOGNORMAL_v1_uint32_uint32_segma_1.000000_16000000 s_LOGNORMAL_v2_uint32_uint32_segma_1.000000_32000000 s_LOGNORMAL_v3_uint32_uint32_segma_1.000000_128000000 s_LOGNORMAL_v4_uint32_uint32_segma_1.000000_384000000 s_LOGNORMAL_v5_uint32_uint32_segma_1.000000_640000000 s_LOGNORMAL_v6_uint32_uint32_segma_1.000000_896000000 s_LOGNORMAL_v7_uint32_uint32_segma_1.000000_1152000000 s_LOGNORMAL_v8_uint32_uint32_segma_1.000000_1664000000 s_LOGNORMAL_v9_uint32_uint32_segma_1.000000_1920000000)
#r_lognormal_datasets_sizes=(32E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
#s_lognormal_datasets_sizes=(32E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
#r_lognormal_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
#s_lognormal_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_lognormal/
#input_hash_table_size=(33554432) #(33554432(for_32E6) 134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

#process_non_imv_indexed_nested_loop_join $r_lognormal_datasets $r_lognormal_datasets_sizes $r_lognormal_datasets_file_num_partitions $s_lognormal_datasets $s_lognormal_datasets_sizes $s_lognormal_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0 0 0 $input_hash_table_size


#r_lognormal_datasets=(r_LOGNORMAL_v3_uint32_uint32_segma_1.000000_128000000) #(r_LOGNORMAL_v1_uint32_uint32_segma_1.000000_16000000 r_LOGNORMAL_v2_uint32_uint32_segma_1.000000_32000000 r_LOGNORMAL_v3_uint32_uint32_segma_1.000000_128000000 r_LOGNORMAL_v4_uint32_uint32_segma_1.000000_384000000 r_LOGNORMAL_v5_uint32_uint32_segma_1.000000_640000000 r_LOGNORMAL_v6_uint32_uint32_segma_1.000000_896000000 r_LOGNORMAL_v7_uint32_uint32_segma_1.000000_1152000000 r_LOGNORMAL_v8_uint32_uint32_segma_1.000000_1664000000 r_LOGNORMAL_v9_uint32_uint32_segma_1.000000_1920000000)
#s_lognormal_datasets=(s_LOGNORMAL_v3_uint32_uint32_segma_1.000000_128000000) #(s_LOGNORMAL_v1_uint32_uint32_segma_1.000000_16000000 s_LOGNORMAL_v2_uint32_uint32_segma_1.000000_32000000 s_LOGNORMAL_v3_uint32_uint32_segma_1.000000_128000000 s_LOGNORMAL_v4_uint32_uint32_segma_1.000000_384000000 s_LOGNORMAL_v5_uint32_uint32_segma_1.000000_640000000 s_LOGNORMAL_v6_uint32_uint32_segma_1.000000_896000000 s_LOGNORMAL_v7_uint32_uint32_segma_1.000000_1152000000 s_LOGNORMAL_v8_uint32_uint32_segma_1.000000_1664000000 s_LOGNORMAL_v9_uint32_uint32_segma_1.000000_1920000000)
#r_lognormal_datasets_sizes=(128E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
#s_lognormal_datasets_sizes=(128E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)#
#r_lognormal_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
#s_lognormal_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_lognormal/
#input_hash_table_size=(134217728) #(134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

#process_non_imv_indexed_nested_loop_join $r_lognormal_datasets $r_lognormal_datasets_sizes $r_lognormal_datasets_file_num_partitions $s_lognormal_datasets $s_lognormal_datasets_sizes $s_lognormal_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0 0 0 $input_hash_table_size

#r_lognormal_datasets=(r_LOGNORMAL_v5_uint32_uint32_segma_1.000000_640000000) #(r_LOGNORMAL_v1_uint32_uint32_segma_1.000000_16000000 r_LOGNORMAL_v2_uint32_uint32_segma_1.000000_32000000 r_LOGNORMAL_v3_uint32_uint32_segma_1.000000_128000000 r_LOGNORMAL_v4_uint32_uint32_segma_1.000000_384000000 r_LOGNORMAL_v5_uint32_uint32_segma_1.000000_640000000 r_LOGNORMAL_v6_uint32_uint32_segma_1.000000_896000000 r_LOGNORMAL_v7_uint32_uint32_segma_1.000000_1152000000 r_LOGNORMAL_v8_uint32_uint32_segma_1.000000_1664000000 r_LOGNORMAL_v9_uint32_uint32_segma_1.000000_1920000000)
#s_lognormal_datasets=(s_LOGNORMAL_v5_uint32_uint32_segma_1.000000_640000000) #(s_LOGNORMAL_v1_uint32_uint32_segma_1.000000_16000000 s_LOGNORMAL_v2_uint32_uint32_segma_1.000000_32000000 s_LOGNORMAL_v3_uint32_uint32_segma_1.000000_128000000 s_LOGNORMAL_v4_uint32_uint32_segma_1.000000_384000000 s_LOGNORMAL_v5_uint32_uint32_segma_1.000000_640000000 s_LOGNORMAL_v6_uint32_uint32_segma_1.000000_896000000 s_LOGNORMAL_v7_uint32_uint32_segma_1.000000_1152000000 s_LOGNORMAL_v8_uint32_uint32_segma_1.000000_1664000000 s_LOGNORMAL_v9_uint32_uint32_segma_1.000000_1920000000)
#r_lognormal_datasets_sizes=(640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
#s_lognormal_datasets_sizes=(640E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
#r_lognormal_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
#s_lognormal_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_lognormal/
#input_hash_table_size=(536870912) #(134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

#process_non_imv_indexed_nested_loop_join $r_lognormal_datasets $r_lognormal_datasets_sizes $r_lognormal_datasets_file_num_partitions $s_lognormal_datasets $s_lognormal_datasets_sizes $s_lognormal_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0 0 0 $input_hash_table_size


#r_lognormal_datasets=(r_LOGNORMAL_v8_uint32_uint32_segma_1.000000_1664000000) #(r_LOGNORMAL_v1_uint32_uint32_segma_1.000000_16000000 r_LOGNORMAL_v2_uint32_uint32_segma_1.000000_32000000 r_LOGNORMAL_v3_uint32_uint32_segma_1.000000_128000000 r_LOGNORMAL_v4_uint32_uint32_segma_1.000000_384000000 r_LOGNORMAL_v5_uint32_uint32_segma_1.000000_640000000 r_LOGNORMAL_v6_uint32_uint32_segma_1.000000_896000000 r_LOGNORMAL_v7_uint32_uint32_segma_1.000000_1152000000 r_LOGNORMAL_v8_uint32_uint32_segma_1.000000_1664000000 r_LOGNORMAL_v9_uint32_uint32_segma_1.000000_1920000000)
#s_lognormal_datasets=(s_LOGNORMAL_v8_uint32_uint32_segma_1.000000_1664000000) #(s_LOGNORMAL_v1_uint32_uint32_segma_1.000000_16000000 s_LOGNORMAL_v2_uint32_uint32_segma_1.000000_32000000 s_LOGNORMAL_v3_uint32_uint32_segma_1.000000_128000000 s_LOGNORMAL_v4_uint32_uint32_segma_1.000000_384000000 s_LOGNORMAL_v5_uint32_uint32_segma_1.000000_640000000 s_LOGNORMAL_v6_uint32_uint32_segma_1.000000_896000000 s_LOGNORMAL_v7_uint32_uint32_segma_1.000000_1152000000 s_LOGNORMAL_v8_uint32_uint32_segma_1.000000_1664000000 s_LOGNORMAL_v9_uint32_uint32_segma_1.000000_1920000000)
#r_lognormal_datasets_sizes=(1664E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
#s_lognormal_datasets_sizes=(1664E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
#r_lognormal_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
#s_lognormal_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_lognormal/
#input_hash_table_size=(1073741824) #(134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

#process_non_imv_indexed_nested_loop_join $r_lognormal_datasets $r_lognormal_datasets_sizes $r_lognormal_datasets_file_num_partitions $s_lognormal_datasets $s_lognormal_datasets_sizes $s_lognormal_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0 0 0 $input_hash_table_size


#r_lognormal_datasets=(r_LOGNORMAL_v9_uint32_uint32_segma_1.000000_1920000000) #(r_LOGNORMAL_v1_uint32_uint32_segma_1.000000_16000000 r_LOGNORMAL_v2_uint32_uint32_segma_1.000000_32000000 r_LOGNORMAL_v3_uint32_uint32_segma_1.000000_128000000 r_LOGNORMAL_v4_uint32_uint32_segma_1.000000_384000000 r_LOGNORMAL_v5_uint32_uint32_segma_1.000000_640000000 r_LOGNORMAL_v6_uint32_uint32_segma_1.000000_896000000 r_LOGNORMAL_v7_uint32_uint32_segma_1.000000_1152000000 r_LOGNORMAL_v8_uint32_uint32_segma_1.000000_1664000000 r_LOGNORMAL_v9_uint32_uint32_segma_1.000000_1920000000)
#s_lognormal_datasets=(s_LOGNORMAL_v9_uint32_uint32_segma_1.000000_1920000000) #(s_LOGNORMAL_v1_uint32_uint32_segma_1.000000_16000000 s_LOGNORMAL_v2_uint32_uint32_segma_1.000000_32000000 s_LOGNORMAL_v3_uint32_uint32_segma_1.000000_128000000 s_LOGNORMAL_v4_uint32_uint32_segma_1.000000_384000000 s_LOGNORMAL_v5_uint32_uint32_segma_1.000000_640000000 s_LOGNORMAL_v6_uint32_uint32_segma_1.000000_896000000 s_LOGNORMAL_v7_uint32_uint32_segma_1.000000_1152000000 s_LOGNORMAL_v8_uint32_uint32_segma_1.000000_1664000000 s_LOGNORMAL_v9_uint32_uint32_segma_1.000000_1920000000)
#r_lognormal_datasets_sizes=(1920E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
#s_lognormal_datasets_sizes=(1920E6) #(16E6 32E6 128E6 384E6 640E6 896E6 1152E6 1664E6 1920E6)
#r_lognormal_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
#s_lognormal_datasets_file_num_partitions=(64) #(64 64 64 64 64 64 64 64 64)
#output_folder_path=/spinning/sabek/learned_join_results/non_imv_inlj_with_hash_index_lognormal/
#input_hash_table_size=(2147483648) #(134217728(for_128E6) 536870912(for_640E6) 1073741824(for_1664E6) 2147483648(for_1920E6))

#process_non_imv_indexed_nested_loop_join $r_lognormal_datasets $r_lognormal_datasets_sizes $r_lognormal_datasets_file_num_partitions $s_lognormal_datasets $s_lognormal_datasets_sizes $s_lognormal_datasets_file_num_partitions $output_folder_path $run_nums $load_relations_for_evaluation $persist_relations_for_evaluation 1 0 0 0 $input_hash_table_size

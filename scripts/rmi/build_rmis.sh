#! /usr/bin/env bash
# An adapted version of build_rmis.sh in the SOSD benchmark 
#git submodule update --init --recursive

function optimize_rmis() {
    DATA_NAME=$1
    DATA_PATH=$2
    RMI_SPECS_FOLDER=$3
    THREADS_NUM=$4

    if [ ! -f ${RMI_SPECS_FOLDER}/${DATA_NAME}.json ]; then
        echo "Optimizing the RMI configs for $DATA_NAME"
        $(dirname "$0")/../../RMI/target/release/rmi --optimize ${RMI_SPECS_FOLDER}/${DATA_NAME}.json --threads ${THREADS_NUM} ${DATA_PATH}
    fi
}


function build_rmi_set() {
    DATA_NAME=$1
    DATA_PATH=$2
    HEADER_FOLDER=$3
    RMI_SPECS_FOLDER=$4
    RMI_DATA_FOLDER=$5
    THREADS_NUM=$6

    JSON_PATH=$(dirname "$0")/../../include/rmi_specs/${DATA_NAME}.json

    shift 1
    if [ ! -f ${HEADER_FOLDER}/${DATA_NAME}_0.h ]; then
        echo "Building RMI set for $DATA_NAME"
        $(dirname "$0")/../../RMI/target/release/rmi ${DATA_PATH} --param-grid ${RMI_SPECS_FOLDER}/${DATA_NAME}.json -d ${RMI_DATA_FOLDER} --threads ${THREADS_NUM} --zero-build-time
        mv ${DATA_NAME}_* ${HEADER_FOLDER}/
    fi
}

mkdir -p /spinning/sabek/rmi_data
mkdir -p $(dirname "$0")/../../include/rmi_specs
mkdir -p $(dirname "$0")/../../include/rmi

rm -r $(dirname "$0")/../../include/rmi/*
rm -r /spinning/sabek/include/rmi_data/*
rm -r $(dirname "$0")/../../include/rmi_specs/*

cd $(dirname "$0")/../../RMI && cargo build --release && cd ../build/release

#optimize_rmis r_LOGNORMAL_v2_int_int_100000000 /spinning/sabek/learned_join_datasets/r_LOGNORMAL_v2_int_int_100000000_key_uint32 $(dirname "$0")/../../include/rmi_specs 8
#build_rmi_set r_LOGNORMAL_v2_int_int_100000000 /spinning/sabek/learned_join_datasets/r_LOGNORMAL_v2_int_int_100000000_key_uint32 $(dirname "$0")/../../include/rmi $(dirname "$0")/../../include/rmi_specs /spinning/sabek/rmi_data 8

optimize_rmis r_UNIQUE_v3_uint32_uint32_128000000 /spinning/sabek/learned_join_datasets/r_UNIQUE_v3_uint32_uint32_128000000_key_uint32 $(dirname "$0")/../../include/rmi_specs 8
build_rmi_set r_UNIQUE_v3_uint32_uint32_128000000 /spinning/sabek/learned_join_datasets/r_UNIQUE_v3_uint32_uint32_128000000_key_uint32 $(dirname "$0")/../../include/rmi $(dirname "$0")/../../include/rmi_specs /spinning/sabek/rmi_data 8


/bin/bash $(dirname "$0")/generate_all_rmis.sh
#!/bin/bash

## Script taken from the blcok partitioning project

SCALE_FACTOR=${1:-1.0}  # Optional: first parameter of script is SCALE_FACTOR, defaults to 1.0
TPCH_TOOL_URL="https://github.com/electrum/tpch-dbgen.git"
TPCH_TOOL_DIR="./cache/tpch-tool"
TPCH_DATA_DIR="./cache/tpch-data-sf=${SCALE_FACTOR}"
CSV_FILE="${TPCH_DATA_DIR}/denormalized-tpch.csv"
DENORMALIZE_SCRIPT_URL="https://raw.githubusercontent.com/Jibbow/denormalized-tpch/master/denormalize.py"
DENORMALIZE_SCRIPT_PATH="./cache/denormalize-tpch.py"

mkdir -p "${TPCH_DATA_DIR}"

if [ ! -d "${TPCH_TOOL_DIR}" ]; then
    echo "Downloading tpch-tool..."
    git clone "${TPCH_TOOL_URL}" "${TPCH_TOOL_DIR}"
fi

if [ ! -f "${TPCH_TOOL_DIR}/dbgen" ]; then
    echo "Compiling tpch-tool..."
    cd "${TPCH_TOOL_DIR}"
    make
    cd -
fi

if [ ! "$(ls -A ${TPCH_DATA_DIR})" ]; then
    echo "Generating TPC-H data..."
    mkdir -p "${TPCH_DATA_DIR}"
    DSS_PATH="${TPCH_DATA_DIR}" \
    DSS_CONFIG="${TPCH_TOOL_DIR}" \
    "${TPCH_TOOL_DIR}/dbgen" -v -f -s ${SCALE_FACTOR}
fi

if [ ! -f ${CSV_FILE} ]; then
    echo "Denormalizing TPC-H data..."
    wget --output-document "${DENORMALIZE_SCRIPT_PATH}" "${DENORMALIZE_SCRIPT_URL}"
    python3 "${DENORMALIZE_SCRIPT_PATH}" --tpch-dir "${TPCH_DATA_DIR}" --out "${CSV_FILE}"
fi

echo "blockpartitioning:data-file=${CSV_FILE}"

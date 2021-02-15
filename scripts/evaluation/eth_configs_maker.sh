#!/bin/sh

# Set some default values:
PROBE_BUFFER_SIZE=4
PADDED_BUCKET=0 # default case: not padded
OVERFLOW_BUF_SIZE=1024 
NRADIXBITS_DEFAULT=7
MERGEBITONICWIDTH=16 #4,8,16
HAVE_AVX=1

# Specific Radix join parameters
NUM_RADIX_BITS=10 #9 10 11 12 14 16
NUM_PASSES=1  # 1 2
SKEWNESS_THRESHOLD_MULTIPLIER=64 #fixed
USE_MURMUR3_HASH_FOR_RADIX_JOIN=1 #fixed
USE_VECTORIZED_MURMUR3_HASH_FOR_RADIX_JOIN=1 #fixed

# Specific non-partition join parameters
BUCKET_SIZE=1
PREFETCH_DISTANCE=16  #10
PREFETCH_NPJ=1 #fixed
USE_MURMUR3_HASH=1 #fixed
INPUT_HASH_TABLE_SIZE=0 #fixed
NPJ_VECTOR_SCALE=8
NPJ_MAX_VECTOR_SCALE=16
NPJ_PDIS=320
NPJ_SIMDStateSize=5
NPJ_ETH_AVX_IMV=1


usage()
{
  echo "Usage: sh base_configs_maker.sh [ --NUM_RADIX_BITS 10] [ --NUM_PASSES 1] [ --USE_MURMUR3_HASH_FOR_RADIX_JOIN 1] [ --USE_VECTORIZED_MURMUR3_HASH_FOR_RADIX_JOIN 1] [ --PROBE_BUFFER_SIZE 4] [ --BUCKET_SIZE 2]
                        [ --PADDED_BUCKET 0] [ --OVERFLOW_BUF_SIZE  1024] [ --SKEWNESS_THRESHOLD_MULTIPLIER 64] 
                        [ --PREFETCH_DISTANCE 16] [ --NRADIXBITS_DEFAULT 7] [ --MERGEBITONICWIDTH 16]
                        [ --PREFETCH_NPJ 1] [ --NPJ_VECTOR_SCALE 8] [ --NPJ_MAX_VECTOR_SCALE 16] [ --USE_MURMUR3_HASH 1] [ --HAVE_AVX 1] [ --INPUT_HASH_TABLE_SIZE 0] [ --NPJ_PDIS 320] [ --NPJ_SIMDStateSize 5] [ --NPJ_ETH_AVX_IMV 1]"
  exit 2
}

PARSED_ARGUMENTS=$(getopt -a -n base_configs -o t --long NUM_RADIX_BITS:,NUM_PASSES:,USE_MURMUR3_HASH_FOR_RADIX_JOIN:,USE_VECTORIZED_MURMUR3_HASH_FOR_RADIX_JOIN:,PROBE_BUFFER_SIZE:,BUCKET_SIZE:,PADDED_BUCKET:,OVERFLOW_BUF_SIZE:,PREFETCH_DISTANCE:,NRADIXBITS_DEFAULT:,MERGEBITONICWIDTH:,PREFETCH_NPJ:,NPJ_VECTOR_SCALE:,NPJ_MAX_VECTOR_SCALE:,USE_MURMUR3_HASH:,HAVE_AVX:,INPUT_HASH_TABLE_SIZE:,NPJ_PDIS:,NPJ_SIMDStateSize:,NPJ_ETH_AVX_IMV:,SKEWNESS_THRESHOLD_MULTIPLIER: -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
  usage
fi

echo "PARSED_ARGUMENTS is $PARSED_ARGUMENTS"
eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    --NUM_RADIX_BITS)   NUM_RADIX_BITS="$2"      ; shift 2  ;;
    --NUM_PASSES)   NUM_PASSES="$2"      ; shift 2  ;;
    --USE_MURMUR3_HASH_FOR_RADIX_JOIN)   USE_MURMUR3_HASH_FOR_RADIX_JOIN="$2"   ; shift 2 ;;
    --USE_VECTORIZED_MURMUR3_HASH_FOR_RADIX_JOIN)   USE_VECTORIZED_MURMUR3_HASH_FOR_RADIX_JOIN="$2"   ; shift 2 ;;
    --PROBE_BUFFER_SIZE) PROBE_BUFFER_SIZE="$2" ;  shift 2 ;;
    --BUCKET_SIZE)   BUCKET_SIZE="$2"   ; shift 2 ;;
    --PADDED_BUCKET)   PADDED_BUCKET="$2"   ; shift 2 ;;
    --OVERFLOW_BUF_SIZE)   OVERFLOW_BUF_SIZE="$2"   ; shift 2 ;;
    --PREFETCH_DISTANCE)   PREFETCH_DISTANCE="$2"   ; shift 2 ;;
    --NRADIXBITS_DEFAULT)   NRADIXBITS_DEFAULT="$2"   ; shift 2 ;;
    --MERGEBITONICWIDTH)   MERGEBITONICWIDTH="$2"   ; shift 2 ;;
    --PREFETCH_NPJ)   PREFETCH_NPJ="$2"   ; shift 2 ;;
    --NPJ_VECTOR_SCALE)   NPJ_VECTOR_SCALE="$2"   ; shift 2 ;;    
    --NPJ_MAX_VECTOR_SCALE)   NPJ_MAX_VECTOR_SCALE="$2"   ; shift 2 ;;    
    --USE_MURMUR3_HASH)   USE_MURMUR3_HASH="$2"   ; shift 2 ;;
    --HAVE_AVX)   HAVE_AVX="$2"  ; shift 2 ;;
    --INPUT_HASH_TABLE_SIZE)   INPUT_HASH_TABLE_SIZE="$2"   ; shift 2 ;;
    --SKEWNESS_THRESHOLD_MULTIPLIER)   SKEWNESS_THRESHOLD_MULTIPLIER="$2"   ; shift 2 ;;
    --NPJ_PDIS)   NPJ_PDIS="$2"   ; shift 2 ;;    
    --NPJ_SIMDStateSize)   NPJ_SIMDStateSize="$2"   ; shift 2 ;;    
    --NPJ_ETH_AVX_IMV)   NPJ_ETH_AVX_IMV="$2"   ; shift 2 ;;    
    # -- means the end of the arguments; drop this, and break out of the while loop
    --) shift; break ;;
    # If invalid options were passed, then getopt should have reported an error,
    # which we checked as VALID_ARGUMENTS when getopt was called...
    *) echo "Unexpected option: $1 - this should not happen."
       usage ;;
  esac
done

PASS1RADIXBITS='('$NUM_RADIX_BITS/$NUM_PASSES')' 
PASS2RADIXBITS='('$NUM_RADIX_BITS-'('$NUM_RADIX_BITS/$NUM_PASSES')'')'

FANOUT_PASS1='('1'<<''('$NUM_RADIX_BITS/$NUM_PASSES'))' #(1<<PASS1RADIXBITS)
FANOUT_PASS2='('1'<<''('$NUM_RADIX_BITS-'('$NUM_RADIX_BITS/$NUM_PASSES')))' #(1<<PASS2RADIXBITS)

PARTFANOUT_DEFAULT='('1'<<'$NRADIXBITS_DEFAULT')'

MWAY_MERGE_BUFFER_SIZE_DEFAULT=L3_CACHE_SIZE #20MB L3 cache as default value

: '
echo "NUM_RADIX_BITS   : $NUM_RADIX_BITS"
echo "NUM_PASSES   : $NUM_PASSES "
echo "USE_MURMUR3_HASH_FOR_RADIX_JOIN   : $USE_MURMUR3_HASH_FOR_RADIX_JOIN"
echo "USE_VECTORIZED_MURMUR3_HASH_FOR_RADIX_JOIN   : $USE_VECTORIZED_MURMUR3_HASH_FOR_RADIX_JOIN"
echo "PROBE_BUFFER_SIZE : $PROBE_BUFFER_SIZE"
echo "BUCKET_SIZE   : $BUCKET_SIZE"
echo "PADDED_BUCKET : $PADDED_BUCKET"
echo "OVERFLOW_BUF_SIZE   : $OVERFLOW_BUF_SIZE"
echo "PREFETCH_DISTANCE   : $PREFETCH_DISTANCE"
echo "NRADIXBITS_DEFAULT   : $NRADIXBITS_DEFAULT"
echo "MERGEBITONICWIDTH   : $MERGEBITONICWIDTH"
echo "PREFETCH_NPJ   : $PREFETCH_NPJ"
echo "NPJ_VECTOR_SCALE   : $NPJ_VECTOR_SCALE"
echo "NPJ_MAX_VECTOR_SCALE   : $NPJ_MAX_VECTOR_SCALE"
echo "USE_MURMUR3_HASH   : $USE_MURMUR3_HASH"
echo "HAVE_AVX   : $HAVE_AVX"
echo "INPUT_HASH_TABLE_SIZE   : $INPUT_HASH_TABLE_SIZE"
echo "PASS1RADIXBITS   : $PASS1RADIXBITS"
echo "PASS2RADIXBITS   : $PASS2RADIXBITS"
echo "FANOUT_PASS1   : $FANOUT_PASS1"
echo "FANOUT_PASS2   : $FANOUT_PASS2"
echo "PARTFANOUT_DEFAULT   : $PARTFANOUT_DEFAULT"
echo "MWAY_MERGE_BUFFER_SIZE_DEFAULT   : $MWAY_MERGE_BUFFER_SIZE_DEFAULT"
echo "SKEWNESS_THRESHOLD_MULTIPLIER   : $SKEWNESS_THRESHOLD_MULTIPLIER"
echo "NPJ_PDIS   : $NPJ_PDIS"
echo "NPJ_SIMDStateSize   : $NPJ_SIMDStateSize"
echo "NPJ_ETH_AVX_IMV   : $NPJ_ETH_AVX_IMV"
echo "Parameters remaining are: $@"
'


echo "$(echo -n '#pragma once'; echo -n $'\n\n';
        
        echo -n '#include "base_configs.h"'; echo -n $'\n\n';

        echo -n $'#define NUM_RADIX_BITS '$NUM_RADIX_BITS; echo -n $'\n\n';

        echo -n $'#define NUM_PASSES '$NUM_PASSES; echo -n $'\n\n';

        if [ "$USE_MURMUR3_HASH_FOR_RADIX_JOIN" = 1 ]; then
           echo -n $'#define USE_MURMUR3_HASH_FOR_RADIX_JOIN '$USE_MURMUR3_HASH_FOR_RADIX_JOIN; echo -n $'\n\n';
        fi

        if [ "$USE_VECTORIZED_MURMUR3_HASH_FOR_RADIX_JOIN" = 1 ]; then
           echo -n $'#define USE_VECTORIZED_MURMUR3_HASH_FOR_RADIX_JOIN '$USE_VECTORIZED_MURMUR3_HASH_FOR_RADIX_JOIN; echo -n $'\n\n';
        fi

        echo -n $'#define PROBE_BUFFER_SIZE '$PROBE_BUFFER_SIZE; echo -n $'\n\n';
        
        echo -n $'#define BUCKET_SIZE '$BUCKET_SIZE; echo -n $'\n\n';

        echo -n $'#define PADDED_BUCKET '$PADDED_BUCKET; echo -n $'\n\n';

        echo -n $'#define OVERFLOW_BUF_SIZE '$OVERFLOW_BUF_SIZE; echo -n $'\n\n';

        echo -n $'#define PREFETCH_DISTANCE '$PREFETCH_DISTANCE; echo -n $'\n\n';

        echo -n $'#define NRADIXBITS_DEFAULT '$NRADIXBITS_DEFAULT; echo -n $'\n\n';

        echo -n $'#define MERGEBITONICWIDTH '$MERGEBITONICWIDTH; echo -n $'\n\n';

        if [ "$PREFETCH_NPJ" = 1 ]; then
             echo -n $'#define PREFETCH_NPJ '$PREFETCH_NPJ; echo -n $'\n\n';
        fi

        if [ "$USE_MURMUR3_HASH" = 1 ]; then
           echo -n $'#define USE_MURMUR3_HASH '$USE_MURMUR3_HASH; echo -n $'\n\n';
        fi

        if [ "$HAVE_AVX" = 1 ]; then
           echo -n $'#define HAVE_AVX '$HAVE_AVX; echo -n $'\n\n';
        fi

        if [ "$INPUT_HASH_TABLE_SIZE" = 1 ]; then
           echo -n $'#define INPUT_HASH_TABLE_SIZE '$INPUT_HASH_TABLE_SIZE; echo -n $'\n\n';
        fi

        echo -n $'#define PASS1RADIXBITS '$PASS1RADIXBITS; echo -n $'\n\n';

        echo -n $'#define PASS2RADIXBITS '$PASS2RADIXBITS; echo -n $'\n\n';

        echo -n $'#define FANOUT_PASS1 '$FANOUT_PASS1; echo -n $'\n\n';

        echo -n $'#define FANOUT_PASS2 '$FANOUT_PASS2; echo -n $'\n\n';

        echo -n $'#define PARTFANOUT_DEFAULT '$PARTFANOUT_DEFAULT; echo -n $'\n\n';

        echo -n $'#define MWAY_MERGE_BUFFER_SIZE_DEFAULT '$MWAY_MERGE_BUFFER_SIZE_DEFAULT; echo -n $'\n\n';

        echo -n $'#define SKEWNESS_THRESHOLD_MULTIPLIER '$SKEWNESS_THRESHOLD_MULTIPLIER; echo -n $'\n\n';

        if [ "$NPJ_ETH_AVX_IMV" = 1 ]; then
           echo -n $'#define NPJ_ETH_AVX_IMV '$NPJ_ETH_AVX_IMV; echo -n $'\n\n';
        fi

        echo -n $'#define NPJ_VECTOR_SCALE '$NPJ_VECTOR_SCALE; echo -n $'\n\n';

        echo -n $'#define NPJ_MAX_VECTOR_SCALE '$NPJ_MAX_VECTOR_SCALE; echo -n $'\n\n';

        echo -n $'#define NPJ_PDIS '$NPJ_PDIS; echo -n $'\n\n';

        echo -n $'#define NPJ_SIMDStateSize '$NPJ_SIMDStateSize; echo -n $'\n\n';

        )" > $(dirname "$0")/../../include/configs/eth_configs.h
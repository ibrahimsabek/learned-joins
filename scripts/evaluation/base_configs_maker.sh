#!/bin/sh

# Some default values:
INTEL_E5=0
INTEL_XEON=1
MAX_CPU_NODES=512
HAVE_LIBNUMA=1
CUSTOM_CPU_MAPPING='"'../../include/configs/cpu-mapping.txt'"'
CUSTOM_CPU_MAPPING_V2='"'../../include/configs/cpu-mapping-v2.txt'"'
MAX_THREADS=$MAX_CPU_NODES 
CACHE_LINE_SIZE=64
SMALL_PADDING_TUPLES_MULTIPLIER=3 
L1_CACHE_SIZE='('1024*1024')' #32768
L1_ASSOCIATIVITY=8
L2_CACHE_SIZE='('2*16*1024*1024')' #(256*1024)
BLOCKSIZE='('$L2_CACHE_SIZE'/(2*sizeof(int64_t)))'
L3_CACHE_SIZE='('2*16*1024*1024')' #(20*1024*1024)
SKEW_HANDLING=1
DEVELOPMENT_MODE=0
USE_SWWC_OPTIMIZED_PART=1
TYPICAL_COPY_MODE=0
USE_AVX_512=1
PERF_COUNTERS=0
ETH_SORT_MERGE_NUMA_STRATEGY=NEXT

# For experimental evaluation
USE_LEARNED_SORT_MERGE_JOIN_GOOGLE_BENCHMARK=1
USE_ETH_SORT_MERGE_JOIN_GOOGLE_BENCHMARK=0
USE_ETH_NON_PARTITION_JOIN_GOOGLE_BENCHMARK=0
USE_ETH_RADIX_JOIN_GOOGLE_BENCHMARK=0
RUN_LEARNED_TECHNIQUES=0
LOAD_RELATIONS_FOR_EVALUATION=1
PERSIST_RELATIONS_FOR_EVALUATION=0
RELATION_KEY_TYPE=uint32_t #fixed
RELATION_PAYLOAD_TYPE=uint32_t #fixed
RELATION_R_PATH='"'/spinning/sabek/learned_join_datasets/r_UNIQUE_v1_uint32_uint32_16000000.txt'"'
RELATION_R_NUM_TUPLES=16E6
RELATION_S_PATH='"'/spinning/sabek/learned_join_datasets/s_UNIQUE_v1_uint32_uint32_16000000.txt'"'
RELATION_S_NUM_TUPLES=16E6
NUM_THREADS_FOR_EVALUATION=4 #2 4 8 16 32 

# General Sort merge join parameters (whether ETH Multi-way or learned)
USE_LEARNED_SORT=1            # enabling the "learned sort" mode, otherwise will use the avx_sort as a black box for the whole sorting phase
USE_LEARNED_SORT_AVX=0  #fixed  # uses the the "original learned sort with AVX extensions" as a black box for the whole sorting phase, if our proposed algorithm not enabled
USE_AVXSORT_AS_STD_SORT=1     # uses avx_sort instead of any "std::sort" call
LS_DEFAULT_BATCH_SZ=10
LS_DEFAULT_FANOUT=1e3
LS_DEFAULT_OVERALLOCATION_RATIO=1.1 #fixed
LS_DEFAULT_SAMPLING_RATE=.01 #fixed
LS_DEFAULT_THRESHOLD=100
LS_DEFAULT_ARCH_SECOND_LEVEL=1000
LS_MIN_SORTING_SIZE=1e4 #fixed

# Specific ETH Multi-way sort merge join parameters
REQUIRED_STACK_SIZE='('2*32*1024*1024')' #fixed #(32*1024*1024)
SKEW_DECOMPOSE_MARGIN=(1.10)                    #10% margin
SKEW_DECOMPOSE_SAMPLES=64                       #nr. of samples for range partitioning
SKEW_MAX_HEAVY_HITTERS=16                       #max nr. of heavy hitters to detect
SKEW_HEAVY_HITTER_THR=0.5                       #heavy hitter threshold freq
ETH_SORT_MERGE_IS_SCALAR_MERGE=0   #fixed
ETH_SORT_MERGE_PARTITION_FANOUT_PER_THREAD=4

# Specific Learned sort merge join parameters
USE_LEARNED_SORT_FOR_SORT_MERGE=1
USE_AVXSORT_FOR_SORTING_MINOR_BCKTS=1
BUILD_RMI_FROM_TWO_DATASETS=1
OVERALLOCATION_SIZE_RATIO=1
REPEATED_KEYS_SIZE_RATIO=0.5
SPILL_BUCKET_SIZE_RATIO=1
USE_FIXED_PARTITION_SIZES=0
LEARNED_SORT_MERGE_PARTITION_FANOUT_PER_THREAD=1 #fixed
LS_FOR_SORT_MERGE_DEFAULT_BATCH_SZ=10
LS_FOR_SORT_MERGE_DEFAULT_FANOUT=1e3
LS_FOR_SORT_MERGE_DEFAULT_OVERALLOCATION_RATIO=1.1 #fixed
LS_FOR_SORT_MERGE_DEFAULT_SAMPLING_RATE=.01 #fixed
LS_FOR_SORT_MERGE_DEFAULT_THRESHOLD=1000
LS_FOR_SORT_MERGE_DEFAULT_ARCH_SECOND_LEVEL=20000
LS_FOR_SORT_MERGE_MIN_SORTING_SIZE=1e4
LS_FOR_SORT_MERGE_IMV_AVX=0
LS_FOR_SORT_MERGE_IMV_AVX_MINOR_BCKTS=0
LS_FOR_SORT_MERGE_PREFETCH_INPUT_FOR_MINOR_BCKTS=0
LS_FOR_SORT_MERGE_PREFETCH_MINOR_BCKT_SIZES_OFF=0
LS_FOR_SORT_MERGE_PREFETCH_SLOPES_AND_INTERCEPTS_MINOR_BCKTS=0
LS_FOR_SORT_MERGE_WORDSIZE=8
LS_FOR_SORT_MERGE_VECTOR_SCALE=8
LS_FOR_SORT_MERGE_MAX_VECTOR_SCALE=16
LS_FOR_SORT_MERGE_SIMDStateSize=5
LS_FOR_SORT_MERGE_PDIS=320


usage()
{
  echo "Usage: sh base_configs_maker.sh [ --INTEL_E5 0] [ --INTEL_XEON 1] [ --MAX_CPU_NODES 512] [ --HAVE_LIBNUMA  1]
                        [ --CUSTOM_CPU_MAPPING src/configs/cpu-mapping.txt] [ --CUSTOM_CPU_MAPPING_V2  src/configs/cpu-mapping-v2.txt] 
                        [ --MAX_THREADS 512] [ --REQUIRED_STACK_SIZE (2*32*1024*1024)]
                        [ --CACHE_LINE_SIZE 64] [ --L1_CACHE_SIZE (1024*1024)] [ --L1_ASSOCIATIVITY 8] [ --L2_CACHE_SIZE (2*16*1024*1024)]
                        [ --BLOCKSIZE (L2_CACHE_SIZE/(2*sizeof(int64_t)))]
                        [ --L3_CACHE_SIZE (2*16*1024*1024)] [ --SKEW_HANDLING 1] [ --SKEW_DECOMPOSE_MARGIN 1.10]
                        [ --SKEW_DECOMPOSE_SAMPLES 64] [ --SKEW_MAX_HEAVY_HITTERS 16] [ --SKEW_HEAVY_HITTER_THR 0.5]
                        [ --PERSIST_RELATIONS 0] [ --DEVELOPMENT_MODE 0] [ --USE_SWWC_OPTIMIZED_PART 1] [--TYPICAL_COPY_MODE 0]                        
                        [ --SMALL_PADDING_TUPLES_MULTIPLIER 3] [ --USE_AVX_512 1] [--USE_LEARNED_SORT_AVX 0]
                        [ --USE_AVXSORT_AS_STD_SORT 1] [ --USE_LEARNED_SORT_FOR_SORT_MERGE 1] [ --USE_AVXSORT_FOR_SORTING_MINOR_BCKTS 1] [ --BUILD_RMI_FROM_TWO_DATASETS 1]
                        [ --OVERALLOCATION_SIZE_RATIO 1] [ --REPEATED_KEYS_SIZE_RATIO 0.5] [ --SPILL_BUCKET_SIZE_RATIO 1]
                        [ --USE_FIXED_PARTITION_SIZES 0] [ --PERF_COUNTERS 0]
                        [ --USE_LEARNED_SORT_MERGE_JOIN_GOOGLE_BENCHMARK 1] [ --USE_ETH_SORT_MERGE_JOIN_GOOGLE_BENCHMARK 0]
                        [ --USE_ETH_NON_PARTITION_JOIN_GOOGLE_BENCHMARK 0] [ --USE_ETH_RADIX_JOIN_GOOGLE_BENCHMARK 0]
                        [ --RUN_LEARNED_TECHNIQUES 0] [ --RUN_LEARNED_TECHNIQUES 0]
                        [ --LOAD_RELATIONS_FOR_EVALUATION 1] [ --PERSIST_RELATIONS_FOR_EVALUATION 0] 
                        [ --RELATION_KEY_TYPE uint32_t] [ --RELATION_PAYLOAD_TYPE uint32_t]
                        [ --RELATION_R_PATH /spinning/sabek/learned_join_datasets/r_UNIQUE_v1_uint32_uint32_16000000.txt] [ --RELATION_R_NUM_TUPLES 16E6]
                        [ --RELATION_S_PATH /spinning/sabek/learned_join_datasets/s_UNIQUE_v1_uint32_uint32_16000000.txt] [ --RELATION_S_NUM_TUPLES 16E6]
                        [ --NUM_THREADS_FOR_EVALUATION 4] [ --ETH_SORT_MERGE_IS_SCALAR_MERGE 0] [ --ETH_SORT_MERGE_NUMA_STRATEGY NEXT] 
                        [ --ETH_SORT_MERGE_PARTITION_FANOUT_PER_THREAD 4] [ --LEARNED_SORT_MERGE_PARTITION_FANOUT_PER_THREAD 1] [-- USE_LEARNED_SORT 1]
                        [ --LS_DEFAULT_BATCH_SZ 10] [ --LS_DEFAULT_FANOUT 1e3] [ --LS_DEFAULT_OVERALLOCATION_RATIO 1.1] [ --LS_DEFAULT_SAMPLING_RATE .01]
                        [ --LS_DEFAULT_THRESHOLD 100] [ --LS_DEFAULT_ARCH_SECOND_LEVEL 1000] [ --LS_MIN_SORTING_SIZE 1e4]
                        [ --LS_FOR_SORT_MERGE_DEFAULT_BATCH_SZ 10] [ --LS_FOR_SORT_MERGE_DEFAULT_FANOUT 1e3] [ --LS_FOR_SORT_MERGE_DEFAULT_OVERALLOCATION_RATIO 1.1] [ --LS_FOR_SORT_MERGE_DEFAULT_SAMPLING_RATE .01]
                        [ --LS_FOR_SORT_MERGE_DEFAULT_THRESHOLD 1000] [ --LS_FOR_SORT_MERGE_DEFAULT_ARCH_SECOND_LEVEL 20000] [ --LS_FOR_SORT_MERGE_MIN_SORTING_SIZE 1e4]
                        [ --LS_FOR_SORT_MERGE_IMV_AVX 0] [ --LS_FOR_SORT_MERGE_IMV_AVX_MINOR_BCKTS 0] [ --LS_FOR_SORT_MERGE_PREFETCH_INPUT_FOR_MINOR_BCKTS 0]
                        [ --LS_FOR_SORT_MERGE_PREFETCH_MINOR_BCKT_SIZES_OFF 0] [ --LS_FOR_SORT_MERGE_PREFETCH_SLOPES_AND_INTERCEPTS_MINOR_BCKTS 0] [ --LS_FOR_SORT_MERGE_WORDSIZE 8]
                        [ --LS_FOR_SORT_MERGE_VECTOR_SCALE 8] [ --LS_FOR_SORT_MERGE_MAX_VECTOR_SCALE 16] [ --LS_FOR_SORT_MERGE_SIMDStateSize 5] [ --LS_FOR_SORT_MERGE_PDIS 320]"
  exit 2
}

PARSED_ARGUMENTS=$(getopt -a -n base_configs -o t --long INTEL_E5:,INTEL_XEON:,MAX_CPU_NODES:,HAVE_LIBNUMA:,CUSTOM_CPU_MAPPING:,CUSTOM_CPU_MAPPING_V2:,MAX_THREADS:,REQUIRED_STACK_SIZE:,CACHE_LINE_SIZE:,L1_CACHE_SIZE:,L1_ASSOCIATIVITY:,L2_CACHE_SIZE:,BLOCKSIZE:,L3_CACHE_SIZE:,SKEW_HANDLING:,SKEW_DECOMPOSE_MARGIN:,SKEW_DECOMPOSE_SAMPLES:,SKEW_MAX_HEAVY_HITTERS:,SKEW_HEAVY_HITTER_THR:,PERSIST_RELATIONS:,DEVELOPMENT_MODE:,USE_SWWC_OPTIMIZED_PART:,TYPICAL_COPY_MODE:,SMALL_PADDING_TUPLES_MULTIPLIER:,USE_AVX_512:,USE_LEARNED_SORT_AVX:,USE_AVXSORT_AS_STD_SORT:,USE_LEARNED_SORT_FOR_SORT_MERGE:,USE_AVXSORT_FOR_SORTING_MINOR_BCKTS:,BUILD_RMI_FROM_TWO_DATASETS:,OVERALLOCATION_SIZE_RATIO:,REPEATED_KEYS_SIZE_RATIO:,SPILL_BUCKET_SIZE_RATIO:,USE_FIXED_PARTITION_SIZES:,PERF_COUNTERS:,USE_ETH_SORT_MERGE_JOIN_GOOGLE_BENCHMARK:,USE_LEARNED_SORT_MERGE_JOIN_GOOGLE_BENCHMARK:,USE_ETH_NON_PARTITION_JOIN_GOOGLE_BENCHMARK:,USE_ETH_RADIX_JOIN_GOOGLE_BENCHMARK:,RUN_LEARNED_TECHNIQUES:,LOAD_RELATIONS_FOR_EVALUATION:,PERSIST_RELATIONS_FOR_EVALUATION:,RELATION_R_PATH:,RELATION_R_NUM_TUPLES:,RELATION_S_PATH:,RELATION_S_NUM_TUPLES:,NUM_THREADS_FOR_EVALUATION:,ETH_SORT_MERGE_IS_SCALAR_MERGE:,ETH_SORT_MERGE_NUMA_STRATEGY:,ETH_SORT_MERGE_PARTITION_FANOUT_PER_THREAD:,LEARNED_SORT_MERGE_PARTITION_FANOUT_PER_THREAD:,USE_LEARNED_SORT:,LS_DEFAULT_BATCH_SZ:,LS_DEFAULT_FANOUT:,LS_DEFAULT_OVERALLOCATION_RATIO:,LS_DEFAULT_SAMPLING_RATE:,LS_DEFAULT_THRESHOLD:,LS_DEFAULT_ARCH_SECOND_LEVEL:,LS_MIN_SORTING_SIZE:,LS_FOR_SORT_MERGE_DEFAULT_BATCH_SZ:,LS_FOR_SORT_MERGE_DEFAULT_FANOUT:,LS_FOR_SORT_MERGE_DEFAULT_OVERALLOCATION_RATIO:,LS_FOR_SORT_MERGE_DEFAULT_SAMPLING_RATE:,LS_FOR_SORT_MERGE_DEFAULT_THRESHOLD:,LS_FOR_SORT_MERGE_DEFAULT_ARCH_SECOND_LEVEL:,LS_FOR_SORT_MERGE_MIN_SORTING_SIZE:,LS_FOR_SORT_MERGE_IMV_AVX:,LS_FOR_SORT_MERGE_IMV_AVX_MINOR_BCKTS:,LS_FOR_SORT_MERGE_PREFETCH_INPUT_FOR_MINOR_BCKTS:,LS_FOR_SORT_MERGE_PREFETCH_MINOR_BCKT_SIZES_OFF:,LS_FOR_SORT_MERGE_PREFETCH_SLOPES_AND_INTERCEPTS_MINOR_BCKTS:,LS_FOR_SORT_MERGE_WORDSIZE:,LS_FOR_SORT_MERGE_VECTOR_SCALE:,LS_FOR_SORT_MERGE_MAX_VECTOR_SCALE:,LS_FOR_SORT_MERGE_SIMDStateSize:,LS_FOR_SORT_MERGE_PDIS: -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
  usage
fi

echo "PARSED_ARGUMENTS is $PARSED_ARGUMENTS"
eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    --INTEL_E5)   INTEL_E5="$2"      ; shift 2  ;;
    --INTEL_XEON)   INTEL_XEON="$2"      ; shift 2  ;;
    --MAX_CPU_NODES) MAX_CPU_NODES="$2" ; MAX_THREADS=$MAX_CPU_NODES ; shift 2 ;;
    --HAVE_LIBNUMA)   HAVE_LIBNUMA="$2"   ; shift 2 ;;
    --CUSTOM_CPU_MAPPING)   CUSTOM_CPU_MAPPING="$2"   ; shift 2 ;;
    --CUSTOM_CPU_MAPPING_V2)   CUSTOM_CPU_MAPPING_V2="$2"   ; shift 2 ;;
    --MAX_THREADS)   MAX_THREADS="$2"   ; shift 2 ;;
    --REQUIRED_STACK_SIZE)   REQUIRED_STACK_SIZE="$2"   ; shift 2 ;;
    --CACHE_LINE_SIZE)   CACHE_LINE_SIZE="$2"   ; shift 2 ;;
    --L1_CACHE_SIZE)   L1_CACHE_SIZE="$2"   ; shift 2 ;;
    --L1_ASSOCIATIVITY)   L1_ASSOCIATIVITY="$2"   ; shift 2 ;;
    --L2_CACHE_SIZE)   L2_CACHE_SIZE="$2"   ; BLOCKSIZE='('$L2_CACHE_SIZE'/(2*sizeof(int64_t)))' ; shift 2 ;;
    --BLOCKSIZE)   BLOCKSIZE="$2"   ; shift 2 ;;
    --L3_CACHE_SIZE)   L3_CACHE_SIZE="$2"   ; shift 2 ;;
    --SKEW_HANDLING)   SKEW_HANDLING="$2"   ; shift 2 ;;
    --SKEW_DECOMPOSE_MARGIN)   SKEW_DECOMPOSE_MARGIN="$2"   ; shift 2 ;;
    --SKEW_DECOMPOSE_SAMPLES)   SKEW_DECOMPOSE_SAMPLES="$2"   ; shift 2 ;;
    --SKEW_MAX_HEAVY_HITTERS)   SKEW_MAX_HEAVY_HITTERS="$2"   ; shift 2 ;;
    --SKEW_HEAVY_HITTER_THR)   SKEW_HEAVY_HITTER_THR="$2"   ; shift 2 ;;
    --PERSIST_RELATIONS)   PERSIST_RELATIONS="$2"   ; shift 2 ;;
    --DEVELOPMENT_MODE)   DEVELOPMENT_MODE="$2"   ; shift 2 ;;
    --USE_SWWC_OPTIMIZED_PART)   USE_SWWC_OPTIMIZED_PART="$2"   ; shift 2 ;;
    --TYPICAL_COPY_MODE)   TYPICAL_COPY_MODE="$2"   ; shift 2 ;;
    --SMALL_PADDING_TUPLES_MULTIPLIER)   SMALL_PADDING_TUPLES_MULTIPLIER="$2"   ; shift 2 ;;
    --USE_AVX_512)   USE_AVX_512="$2"   ; shift 2 ;;
    --USE_LEARNED_SORT_AVX)   USE_LEARNED_SORT_AVX="$2"   ; shift 2 ;;
    --USE_AVXSORT_AS_STD_SORT)   USE_AVXSORT_AS_STD_SORT="$2"   ; shift 2 ;;
    --USE_LEARNED_SORT_FOR_SORT_MERGE)   USE_LEARNED_SORT_FOR_SORT_MERGE="$2"   ; shift 2 ;;
    --USE_AVXSORT_FOR_SORTING_MINOR_BCKTS)   USE_AVXSORT_FOR_SORTING_MINOR_BCKTS="$2"   ; shift 2 ;;
    --BUILD_RMI_FROM_TWO_DATASETS)   BUILD_RMI_FROM_TWO_DATASETS="$2"   ; shift 2 ;;
    --OVERALLOCATION_SIZE_RATIO)   OVERALLOCATION_SIZE_RATIO="$2"   ; shift 2 ;;
    --REPEATED_KEYS_SIZE_RATIO)   REPEATED_KEYS_SIZE_RATIO="$2"   ; shift 2 ;;
    --SPILL_BUCKET_SIZE_RATIO)   SPILL_BUCKET_SIZE_RATIO="$2"   ; shift 2 ;;    
    --USE_FIXED_PARTITION_SIZES)   USE_FIXED_PARTITION_SIZES="$2"   ; shift 2 ;;    
    --PERF_COUNTERS)   PERF_COUNTERS="$2"   ; shift 2 ;;    
    --USE_LEARNED_SORT_MERGE_JOIN_GOOGLE_BENCHMARK)   USE_LEARNED_SORT_MERGE_JOIN_GOOGLE_BENCHMARK="$2"   ; shift 2 ;;    
    --USE_ETH_SORT_MERGE_JOIN_GOOGLE_BENCHMARK)   USE_ETH_SORT_MERGE_JOIN_GOOGLE_BENCHMARK="$2"   ; shift 2 ;;
    --USE_ETH_NON_PARTITION_JOIN_GOOGLE_BENCHMARK)   USE_ETH_NON_PARTITION_JOIN_GOOGLE_BENCHMARK="$2"   ; shift 2 ;;
    --USE_ETH_RADIX_JOIN_GOOGLE_BENCHMARK)   USE_ETH_RADIX_JOIN_GOOGLE_BENCHMARK="$2"   ; shift 2 ;;    
    --RUN_LEARNED_TECHNIQUES)   RUN_LEARNED_TECHNIQUES="$2"   ; shift 2 ;;    
    --LOAD_RELATIONS_FOR_EVALUATION)   LOAD_RELATIONS_FOR_EVALUATION="$2"   ; shift 2 ;;    
    --PERSIST_RELATIONS_FOR_EVALUATION)   PERSIST_RELATIONS_FOR_EVALUATION="$2"   ; shift 2 ;;    
    --RELATION_KEY_TYPE)   RELATION_KEY_TYPE="$2"   ; shift 2 ;;    
    --RELATION_PAYLOAD_TYPE)   RELATION_PAYLOAD_TYPE="$2"   ; shift 2 ;;    
    --RELATION_R_PATH)   RELATION_R_PATH="$2"   ; shift 2 ;;    
    --RELATION_R_NUM_TUPLES)   RELATION_R_NUM_TUPLES="$2"   ; shift 2 ;;    
    --RELATION_S_PATH)   RELATION_S_PATH="$2"   ; shift 2 ;;    
    --RELATION_S_NUM_TUPLES)   RELATION_S_NUM_TUPLES="$2"   ; shift 2 ;;    
    --NUM_THREADS_FOR_EVALUATION)   NUM_THREADS_FOR_EVALUATION="$2"   ; shift 2 ;;    
    --ETH_SORT_MERGE_IS_SCALAR_MERGE)   ETH_SORT_MERGE_IS_SCALAR_MERGE="$2"   ; shift 2 ;;    
    --ETH_SORT_MERGE_NUMA_STRATEGY)   ETH_SORT_MERGE_NUMA_STRATEGY="$2"   ; shift 2 ;;    
    --ETH_SORT_MERGE_PARTITION_FANOUT_PER_THREAD)   ETH_SORT_MERGE_PARTITION_FANOUT_PER_THREAD="$2"   ; shift 2 ;;    
    --LEARNED_SORT_MERGE_PARTITION_FANOUT_PER_THREAD)   LEARNED_SORT_MERGE_PARTITION_FANOUT_PER_THREAD="$2"   ; shift 2 ;;    
    --USE_LEARNED_SORT)   USE_LEARNED_SORT="$2"   ; shift 2 ;;    
    --LS_DEFAULT_BATCH_SZ)   LS_DEFAULT_BATCH_SZ="$2"   ; shift 2 ;;    
    --LS_DEFAULT_FANOUT)   LS_DEFAULT_FANOUT="$2"   ; shift 2 ;;    
    --LS_DEFAULT_OVERALLOCATION_RATIO)   LS_DEFAULT_OVERALLOCATION_RATIO="$2"   ; shift 2 ;;    
    --LS_DEFAULT_SAMPLING_RATE)   LS_DEFAULT_SAMPLING_RATE="$2"   ; shift 2 ;;    
    --LS_DEFAULT_THRESHOLD)   LS_DEFAULT_THRESHOLD="$2"   ; shift 2 ;;    
    --LS_DEFAULT_ARCH_SECOND_LEVEL)   LS_DEFAULT_ARCH_SECOND_LEVEL="$2"   ; shift 2 ;;    
    --LS_MIN_SORTING_SIZE)   LS_MIN_SORTING_SIZE="$2"   ; shift 2 ;;    
    --LS_FOR_SORT_MERGE_DEFAULT_BATCH_SZ)   LS_FOR_SORT_MERGE_DEFAULT_BATCH_SZ="$2"   ; shift 2 ;;    
    --LS_FOR_SORT_MERGE_DEFAULT_FANOUT)   LS_FOR_SORT_MERGE_DEFAULT_FANOUT="$2"   ; shift 2 ;;    
    --LS_FOR_SORT_MERGE_DEFAULT_OVERALLOCATION_RATIO)   LS_FOR_SORT_MERGE_DEFAULT_OVERALLOCATION_RATIO="$2"   ; shift 2 ;;    
    --LS_FOR_SORT_MERGE_DEFAULT_SAMPLING_RATE)   LS_FOR_SORT_MERGE_DEFAULT_SAMPLING_RATE="$2"   ; shift 2 ;;    
    --LS_FOR_SORT_MERGE_DEFAULT_THRESHOLD)   LS_FOR_SORT_MERGE_DEFAULT_THRESHOLD="$2"   ; shift 2 ;;    
    --LS_FOR_SORT_MERGE_DEFAULT_ARCH_SECOND_LEVEL)   LS_FOR_SORT_MERGE_DEFAULT_ARCH_SECOND_LEVEL="$2"   ; shift 2 ;;    
    --LS_FOR_SORT_MERGE_MIN_SORTING_SIZE)   LS_FOR_SORT_MERGE_MIN_SORTING_SIZE="$2"   ; shift 2 ;;    
    --LS_FOR_SORT_MERGE_IMV_AVX)   LS_FOR_SORT_MERGE_IMV_AVX="$2"   ; shift 2 ;;    
    --LS_FOR_SORT_MERGE_IMV_AVX_MINOR_BCKTS)   LS_FOR_SORT_MERGE_IMV_AVX_MINOR_BCKTS="$2"   ; shift 2 ;;    
    --LS_FOR_SORT_MERGE_PREFETCH_INPUT_FOR_MINOR_BCKTS)   LS_FOR_SORT_MERGE_PREFETCH_INPUT_FOR_MINOR_BCKTS="$2"   ; shift 2 ;;    
    --LS_FOR_SORT_MERGE_PREFETCH_MINOR_BCKT_SIZES_OFF)   LS_FOR_SORT_MERGE_PREFETCH_MINOR_BCKT_SIZES_OFF="$2"   ; shift 2 ;;    
    --LS_FOR_SORT_MERGE_PREFETCH_SLOPES_AND_INTERCEPTS_MINOR_BCKTS)   LS_FOR_SORT_MERGE_PREFETCH_SLOPES_AND_INTERCEPTS_MINOR_BCKTS="$2"   ; shift 2 ;;    
    --LS_FOR_SORT_MERGE_WORDSIZE)   LS_FOR_SORT_MERGE_WORDSIZE="$2"   ; shift 2 ;;    
    --LS_FOR_SORT_MERGE_VECTOR_SCALE)   LS_FOR_SORT_MERGE_VECTOR_SCALE="$2"   ; shift 2 ;;    
    --LS_FOR_SORT_MERGE_MAX_VECTOR_SCALE)   LS_FOR_SORT_MERGE_MAX_VECTOR_SCALE="$2"   ; shift 2 ;;    
    --LS_FOR_SORT_MERGE_SIMDStateSize)   LS_FOR_SORT_MERGE_SIMDStateSize="$2"   ; shift 2 ;;    
    --LS_FOR_SORT_MERGE_PDIS)   LS_FOR_SORT_MERGE_PDIS="$2"   ; shift 2 ;;    
    # -- means the end of the arguments; drop this, and break out of the while loop
    --) shift; break ;;
    # If invalid options were passed, then getopt should have reported an error,
    # which we checked as VALID_ARGUMENTS when getopt was called...
    *) echo "Unexpected option: $1 - this should not happen."
       usage ;;
  esac
done

LOG2_BLOCKSIZE='(log2('$BLOCKSIZE'))'
USE_AVXMERGE_AS_STD_MERGE=$USE_AVXSORT_AS_STD_SORT

: '
echo "INTEL_E5   : $INTEL_E5"
echo "INTEL_XEON   : $INTEL_XEON "
echo "MAX_CPU_NODES : $MAX_CPU_NODES"
echo "HAVE_LIBNUMA   : $HAVE_LIBNUMA"
echo "CUSTOM_CPU_MAPPING : $CUSTOM_CPU_MAPPING"
echo "CUSTOM_CPU_MAPPING_V2   : $CUSTOM_CPU_MAPPING_V2"
echo "MAX_THREADS   : $MAX_THREADS"
echo "REQUIRED_STACK_SIZE   : $REQUIRED_STACK_SIZE"
echo "CACHE_LINE_SIZE   : $CACHE_LINE_SIZE"
echo "L1_CACHE_SIZE   : $L1_CACHE_SIZE"
echo "L1_ASSOCIATIVITY   : $L1_ASSOCIATIVITY"
echo "L2_CACHE_SIZE   : $L2_CACHE_SIZE"
echo "BLOCKSIZE   : $BLOCKSIZE"
echo "LOG2_BLOCKSIZE   : $LOG2_BLOCKSIZE"
echo "L3_CACHE_SIZE   : $L3_CACHE_SIZE"
echo "SKEW_HANDLING   : $SKEW_HANDLING"
echo "SKEW_DECOMPOSE_MARGIN   : $SKEW_DECOMPOSE_MARGIN"
echo "SKEW_DECOMPOSE_SAMPLES   : $SKEW_DECOMPOSE_SAMPLES"
echo "SKEW_MAX_HEAVY_HITTERS   : $SKEW_MAX_HEAVY_HITTERS"
echo "SKEW_HEAVY_HITTER_THR   : $SKEW_HEAVY_HITTER_THR"
echo "PERSIST_RELATIONS   : $PERSIST_RELATIONS"
echo "DEVELOPMENT_MODE   : $DEVELOPMENT_MODE"
echo "USE_SWWC_OPTIMIZED_PART   : $USE_SWWC_OPTIMIZED_PART"
echo "TYPICAL_COPY_MODE   : $TYPICAL_COPY_MODE"
echo "SMALL_PADDING_TUPLES_MULTIPLIER   : $SMALL_PADDING_TUPLES_MULTIPLIER"
echo "USE_AVX_512   : $USE_AVX_512"
echo "USE_LEARNED_SORT   : $USE_LEARNED_SORT"
echo "USE_LEARNED_SORT_AVX   : $USE_LEARNED_SORT_AVX"
echo "USE_AVXSORT_AS_STD_SORT   : $USE_AVXSORT_AS_STD_SORT"
echo "USE_AVXMERGE_AS_STD_MERGE   : $USE_AVXMERGE_AS_STD_MERGE"
echo "USE_LEARNED_SORT_FOR_SORT_MERGE   : $USE_LEARNED_SORT_FOR_SORT_MERGE"
echo "USE_AVXSORT_FOR_SORTING_MINOR_BCKTS   : $USE_AVXSORT_FOR_SORTING_MINOR_BCKTS"
echo "BUILD_RMI_FROM_TWO_DATASETS   : $BUILD_RMI_FROM_TWO_DATASETS"
echo "OVERALLOCATION_SIZE_RATIO   : $OVERALLOCATION_SIZE_RATIO"
echo "REPEATED_KEYS_SIZE_RATIO   : $REPEATED_KEYS_SIZE_RATIO"
echo "SPILL_BUCKET_SIZE_RATIO   : $SPILL_BUCKET_SIZE_RATIO"
echo "USE_FIXED_PARTITION_SIZES   : $USE_FIXED_PARTITION_SIZES"
echo "PERF_COUNTERS   : $PERF_COUNTERS"
echo "USE_LEARNED_SORT_MERGE_JOIN_GOOGLE_BENCHMARK   : $USE_LEARNED_SORT_MERGE_JOIN_GOOGLE_BENCHMARK"
echo "USE_ETH_SORT_MERGE_JOIN_GOOGLE_BENCHMARK   : $USE_ETH_SORT_MERGE_JOIN_GOOGLE_BENCHMARK"
echo "USE_ETH_NON_PARTITION_JOIN_GOOGLE_BENCHMARK   : $USE_ETH_NON_PARTITION_JOIN_GOOGLE_BENCHMARK"
echo "USE_ETH_RADIX_JOIN_GOOGLE_BENCHMARK   : $USE_ETH_RADIX_JOIN_GOOGLE_BENCHMARK"
echo "RUN_LEARNED_TECHNIQUES   : $RUN_LEARNED_TECHNIQUES"
echo "LOAD_RELATIONS_FOR_EVALUATION   : $LOAD_RELATIONS_FOR_EVALUATION"
echo "PERSIST_RELATIONS_FOR_EVALUATION   : $PERSIST_RELATIONS_FOR_EVALUATION"
echo "RELATION_KEY_TYPE   : $RELATION_KEY_TYPE"
echo "RELATION_PAYLOAD_TYPE   : $RELATION_PAYLOAD_TYPE"
echo "RELATION_R_PATH   : $RELATION_R_PATH"
echo "RELATION_R_NUM_TUPLES   : $RELATION_R_NUM_TUPLES"
echo "RELATION_S_PATH   : $RELATION_S_PATH"
echo "RELATION_S_NUM_TUPLES   : $RELATION_S_NUM_TUPLES"
echo "NUM_THREADS_FOR_EVALUATION   : $NUM_THREADS_FOR_EVALUATION"
echo "ETH_SORT_MERGE_IS_SCALAR_MERGE   : $ETH_SORT_MERGE_IS_SCALAR_MERGE"
echo "ETH_SORT_MERGE_NUMA_STRATEGY   : $ETH_SORT_MERGE_NUMA_STRATEGY"
echo "ETH_SORT_MERGE_PARTITION_FANOUT_PER_THREAD   : $ETH_SORT_MERGE_PARTITION_FANOUT_PER_THREAD"
echo "LEARNED_SORT_MERGE_PARTITION_FANOUT_PER_THREAD   : $LEARNED_SORT_MERGE_PARTITION_FANOUT_PER_THREAD"
echo "LS_DEFAULT_BATCH_SZ   : $LS_DEFAULT_BATCH_SZ"
echo "LS_DEFAULT_FANOUT   : $LS_DEFAULT_FANOUT"
echo "LS_DEFAULT_OVERALLOCATION_RATIO   : $LS_DEFAULT_OVERALLOCATION_RATIO"
echo "LS_DEFAULT_SAMPLING_RATE   : $LS_DEFAULT_SAMPLING_RATE"
echo "LS_DEFAULT_THRESHOLD   : $LS_DEFAULT_THRESHOLD"
echo "LS_DEFAULT_ARCH_SECOND_LEVEL   : $LS_DEFAULT_ARCH_SECOND_LEVEL"
echo "LS_MIN_SORTING_SIZE   : $LS_MIN_SORTING_SIZE"
echo "LS_FOR_SORT_MERGE_DEFAULT_BATCH_SZ   : $LS_FOR_SORT_MERGE_DEFAULT_BATCH_SZ"
echo "LS_FOR_SORT_MERGE_DEFAULT_FANOUT   : $LS_FOR_SORT_MERGE_DEFAULT_FANOUT"
echo "LS_FOR_SORT_MERGE_DEFAULT_OVERALLOCATION_RATIO   : $LS_FOR_SORT_MERGE_DEFAULT_OVERALLOCATION_RATIO"
echo "LS_FOR_SORT_MERGE_DEFAULT_SAMPLING_RATE   : $LS_FOR_SORT_MERGE_DEFAULT_SAMPLING_RATE"
echo "LS_FOR_SORT_MERGE_DEFAULT_THRESHOLD   : $LS_FOR_SORT_MERGE_DEFAULT_THRESHOLD"
echo "LS_FOR_SORT_MERGE_DEFAULT_ARCH_SECOND_LEVEL   : $LS_FOR_SORT_MERGE_DEFAULT_ARCH_SECOND_LEVEL"
echo "LS_FOR_SORT_MERGE_MIN_SORTING_SIZE   : $LS_FOR_SORT_MERGE_MIN_SORTING_SIZE"
echo "LS_FOR_SORT_MERGE_IMV_AVX   : $LS_FOR_SORT_MERGE_IMV_AVX"
echo "LS_FOR_SORT_MERGE_IMV_AVX_MINOR_BCKTS   : $LS_FOR_SORT_MERGE_IMV_AVX_MINOR_BCKTS"
echo "LS_FOR_SORT_MERGE_PREFETCH_INPUT_FOR_MINOR_BCKTS   : $LS_FOR_SORT_MERGE_PREFETCH_INPUT_FOR_MINOR_BCKTS"
echo "LS_FOR_SORT_MERGE_PREFETCH_MINOR_BCKT_SIZES_OFF   : $LS_FOR_SORT_MERGE_PREFETCH_MINOR_BCKT_SIZES_OFF"
echo "LS_FOR_SORT_MERGE_PREFETCH_SLOPES_AND_INTERCEPTS_MINOR_BCKTS   : $LS_FOR_SORT_MERGE_PREFETCH_SLOPES_AND_INTERCEPTS_MINOR_BCKTS"
echo "LS_FOR_SORT_MERGE_WORDSIZE   : $LS_FOR_SORT_MERGE_WORDSIZE"
echo "LS_FOR_SORT_MERGE_VECTOR_SCALE   : $LS_FOR_SORT_MERGE_VECTOR_SCALE"
echo "LS_FOR_SORT_MERGE_MAX_VECTOR_SCALE   : $LS_FOR_SORT_MERGE_MAX_VECTOR_SCALE"
echo "LS_FOR_SORT_MERGE_SIMDStateSize   : $LS_FOR_SORT_MERGE_SIMDStateSize"
echo "LS_FOR_SORT_MERGE_PDIS   : $LS_FOR_SORT_MERGE_PDIS"
echo "Parameters remaining are: $@" 
'

echo "$(echo -n '#pragma once'; echo -n $'\n\n';
        
        echo -n $'#define INTEL_E5 '$INTEL_E5; echo -n $'\n\n';

        echo -n $'#define INTEL_XEON '$INTEL_XEON; echo -n $'\n\n';
        
        echo -n $'#define MAX_CPU_NODES '$MAX_CPU_NODES; echo -n $'\n\n';

        echo -n $'#define CUSTOM_CPU_MAPPING '$CUSTOM_CPU_MAPPING; echo -n $'\n\n';

        echo -n $'#define CUSTOM_CPU_MAPPING_V2 '$CUSTOM_CPU_MAPPING_V2; echo -n $'\n\n';

        echo -n $'#define HAVE_LIBNUMA '$HAVE_LIBNUMA; echo -n $'\n\n';

        echo -n $'#define MAX_THREADS '$MAX_THREADS; echo -n $'\n\n';

        echo -n $'#define REQUIRED_STACK_SIZE '$REQUIRED_STACK_SIZE; echo -n $'\n\n';

        echo -n $'#define CACHE_LINE_SIZE '$CACHE_LINE_SIZE; echo -n $'\n\n';

        echo -n $'#define L1_CACHE_SIZE '$L1_CACHE_SIZE; echo -n $'\n\n';

        echo -n $'#define L1_ASSOCIATIVITY '$L1_ASSOCIATIVITY; echo -n $'\n\n';

        echo -n $'#define L2_CACHE_SIZE '$L2_CACHE_SIZE; echo -n $'\n\n';

        echo -n $'#define BLOCKSIZE '$BLOCKSIZE; echo -n $'\n\n';

        echo -n $'#define LOG2_BLOCKSIZE '$LOG2_BLOCKSIZE; echo -n $'\n\n';

        echo -n $'#define L3_CACHE_SIZE '$L3_CACHE_SIZE; echo -n $'\n\n';

        echo -n $'#define SKEW_HANDLING '$SKEW_HANDLING; echo -n $'\n\n';

        echo -n $'#define SKEW_DECOMPOSE_MARGIN '$SKEW_DECOMPOSE_MARGIN; echo -n $'\n\n';

        echo -n $'#define SKEW_DECOMPOSE_SAMPLES '$SKEW_DECOMPOSE_SAMPLES; echo -n $'\n\n';

        echo -n $'#define SKEW_MAX_HEAVY_HITTERS '$SKEW_MAX_HEAVY_HITTERS; echo -n $'\n\n';

        echo -n $'#define SKEW_HEAVY_HITTER_THR '$SKEW_HEAVY_HITTER_THR; echo -n $'\n\n';

	     if [ "$PERSIST_RELATIONS" = 1 ]; then
 	         echo -n $'#define PERSIST_RELATIONS '$PERSIST_RELATIONS; echo -n $'\n\n';
	     fi

	     if [ "$DEVELOPMENT_MODE" = 1 ]; then
 	         echo -n $'#define DEVELOPMENT_MODE '$DEVELOPMENT_MODE; echo -n $'\n\n';
	     fi

        echo -n $'#define USE_SWWC_OPTIMIZED_PART '$USE_SWWC_OPTIMIZED_PART; echo -n $'\n\n';

        if [ "$TYPICAL_COPY_MODE" = 1 ]; then
           echo -n $'#define TYPICAL_COPY_MODE '$TYPICAL_COPY_MODE; echo -n $'\n\n';
        fi

        echo -n $'#define SMALL_PADDING_TUPLES_MULTIPLIER '$SMALL_PADDING_TUPLES_MULTIPLIER; echo -n $'\n\n';

        if [ "$USE_AVX_512" = 1 ]; then
           echo -n $'#define USE_AVX_512 '$USE_AVX_512; echo -n $'\n\n';
        fi

        echo -n $'#define USE_LEARNED_SORT '$USE_LEARNED_SORT; echo -n $'\n\n';

        if [ "$USE_LEARNED_SORT_AVX" = 1 ]; then
           echo -n $'#define USE_LEARNED_SORT_AVX '$USE_LEARNED_SORT_AVX; echo -n $'\n\n';
        fi

        if [ "$USE_AVXSORT_AS_STD_SORT" = 1 ]; then
           echo -n $'#define USE_AVXSORT_AS_STD_SORT '$USE_AVXSORT_AS_STD_SORT; echo -n $'\n\n';
        fi

        if [ "$USE_AVXSORT_AS_STD_SORT" = 1 ]; then
           echo -n $'#define USE_AVXMERGE_AS_STD_MERGE '$USE_AVXMERGE_AS_STD_MERGE; echo -n $'\n\n';
        fi


        if [ "$USE_LEARNED_SORT_FOR_SORT_MERGE" = 1 ]; then
           echo -n $'#define USE_LEARNED_SORT_FOR_SORT_MERGE '$USE_LEARNED_SORT_FOR_SORT_MERGE; echo -n $'\n\n';
        fi

        if [ "$USE_AVXSORT_FOR_SORTING_MINOR_BCKTS" = 1 ]; then
           echo -n $'#define USE_AVXSORT_FOR_SORTING_MINOR_BCKTS '$USE_AVXSORT_FOR_SORTING_MINOR_BCKTS; echo -n $'\n\n';
        fi

        if [ "$BUILD_RMI_FROM_TWO_DATASETS" = 1 ]; then
           echo -n $'#define BUILD_RMI_FROM_TWO_DATASETS '$BUILD_RMI_FROM_TWO_DATASETS; echo -n $'\n\n';
        fi

        echo -n $'#define OVERALLOCATION_SIZE_RATIO '$OVERALLOCATION_SIZE_RATIO; echo -n $'\n\n';

        echo -n $'#define REPEATED_KEYS_SIZE_RATIO '$REPEATED_KEYS_SIZE_RATIO; echo -n $'\n\n';

        echo -n $'#define SPILL_BUCKET_SIZE_RATIO '$SPILL_BUCKET_SIZE_RATIO; echo -n $'\n\n';

        if [ "$USE_FIXED_PARTITION_SIZES" = 1 ]; then
           echo -n $'#define USE_FIXED_PARTITION_SIZES '$USE_FIXED_PARTITION_SIZES; echo -n $'\n\n';
        fi

        if [ "$PERF_COUNTERS" = 1 ]; then
           echo -n $'#define PERF_COUNTERS '$PERF_COUNTERS; echo -n $'\n\n';
        fi

        if [ "$USE_LEARNED_SORT_MERGE_JOIN_GOOGLE_BENCHMARK" = 1 ]; then
           echo -n $'#define USE_LEARNED_SORT_MERGE_JOIN_GOOGLE_BENCHMARK '$USE_LEARNED_SORT_MERGE_JOIN_GOOGLE_BENCHMARK; echo -n $'\n\n';
        fi

        if [ "$USE_ETH_SORT_MERGE_JOIN_GOOGLE_BENCHMARK" = 1 ]; then
           echo -n $'#define USE_ETH_SORT_MERGE_JOIN_GOOGLE_BENCHMARK '$USE_ETH_SORT_MERGE_JOIN_GOOGLE_BENCHMARK; echo -n $'\n\n';
        fi

        if [ "$USE_ETH_NON_PARTITION_JOIN_GOOGLE_BENCHMARK" = 1 ]; then
           echo -n $'#define USE_ETH_NON_PARTITION_JOIN_GOOGLE_BENCHMARK '$USE_ETH_NON_PARTITION_JOIN_GOOGLE_BENCHMARK; echo -n $'\n\n';
        fi

        if [ "$USE_ETH_RADIX_JOIN_GOOGLE_BENCHMARK" = 1 ]; then
           echo -n $'#define USE_ETH_RADIX_JOIN_GOOGLE_BENCHMARK '$USE_ETH_RADIX_JOIN_GOOGLE_BENCHMARK; echo -n $'\n\n';
        fi

        if [ "$RUN_LEARNED_TECHNIQUES" = 1 ]; then
           echo -n $'#define RUN_LEARNED_TECHNIQUES '$RUN_LEARNED_TECHNIQUES; echo -n $'\n\n';
        fi

	     if [ "$LOAD_RELATIONS_FOR_EVALUATION" = 1 ]; then
 	         echo -n $'#define LOAD_RELATIONS_FOR_EVALUATION '$LOAD_RELATIONS_FOR_EVALUATION; echo -n $'\n\n';
	     fi

	     if [ "$PERSIST_RELATIONS_FOR_EVALUATION" = 1 ]; then
 	         echo -n $'#define PERSIST_RELATIONS_FOR_EVALUATION '$PERSIST_RELATIONS_FOR_EVALUATION; echo -n $'\n\n';
	     fi

        echo -n $'#define RELATION_KEY_TYPE '$RELATION_KEY_TYPE; echo -n $'\n\n';

        echo -n $'#define RELATION_PAYLOAD_TYPE '$RELATION_PAYLOAD_TYPE; echo -n $'\n\n';

        echo -n $'#define RELATION_R_PATH '$RELATION_R_PATH; echo -n $'\n\n';

        echo -n $'#define RELATION_R_NUM_TUPLES '$RELATION_R_NUM_TUPLES; echo -n $'\n\n';

        echo -n $'#define RELATION_S_PATH '$RELATION_S_PATH; echo -n $'\n\n';

        echo -n $'#define RELATION_S_NUM_TUPLES '$RELATION_S_NUM_TUPLES; echo -n $'\n\n';

        echo -n $'#define NUM_THREADS_FOR_EVALUATION '$NUM_THREADS_FOR_EVALUATION; echo -n $'\n\n';

        echo -n $'#define ETH_SORT_MERGE_IS_SCALAR_MERGE '$ETH_SORT_MERGE_IS_SCALAR_MERGE; echo -n $'\n\n';

        echo -n $'#define ETH_SORT_MERGE_NUMA_STRATEGY '$ETH_SORT_MERGE_NUMA_STRATEGY; echo -n $'\n\n';

        echo -n $'#define ETH_SORT_MERGE_PARTITION_FANOUT_PER_THREAD '$ETH_SORT_MERGE_PARTITION_FANOUT_PER_THREAD; echo -n $'\n\n';

        echo -n $'#define LEARNED_SORT_MERGE_PARTITION_FANOUT_PER_THREAD '$LEARNED_SORT_MERGE_PARTITION_FANOUT_PER_THREAD; echo -n $'\n\n';

        echo -n $'#define LS_DEFAULT_BATCH_SZ '$LS_DEFAULT_BATCH_SZ; echo -n $'\n\n';

        echo -n $'#define LS_DEFAULT_FANOUT '$LS_DEFAULT_FANOUT; echo -n $'\n\n';

        echo -n $'#define LS_DEFAULT_OVERALLOCATION_RATIO '$LS_DEFAULT_OVERALLOCATION_RATIO; echo -n $'\n\n';

        echo -n $'#define LS_DEFAULT_SAMPLING_RATE '$LS_DEFAULT_SAMPLING_RATE; echo -n $'\n\n';

        echo -n $'#define LS_DEFAULT_THRESHOLD '$LS_DEFAULT_THRESHOLD; echo -n $'\n\n';

        echo -n $'#define LS_DEFAULT_ARCH_SECOND_LEVEL '$LS_DEFAULT_ARCH_SECOND_LEVEL; echo -n $'\n\n';

        echo -n $'#define LS_MIN_SORTING_SIZE '$LS_MIN_SORTING_SIZE; echo -n $'\n\n';

        echo -n $'#define LS_FOR_SORT_MERGE_DEFAULT_BATCH_SZ '$LS_FOR_SORT_MERGE_DEFAULT_BATCH_SZ; echo -n $'\n\n';

        echo -n $'#define LS_FOR_SORT_MERGE_DEFAULT_FANOUT '$LS_FOR_SORT_MERGE_DEFAULT_FANOUT; echo -n $'\n\n';

        echo -n $'#define LS_FOR_SORT_MERGE_DEFAULT_OVERALLOCATION_RATIO '$LS_FOR_SORT_MERGE_DEFAULT_OVERALLOCATION_RATIO; echo -n $'\n\n';

        echo -n $'#define LS_FOR_SORT_MERGE_DEFAULT_SAMPLING_RATE '$LS_FOR_SORT_MERGE_DEFAULT_SAMPLING_RATE; echo -n $'\n\n';

        echo -n $'#define LS_FOR_SORT_MERGE_DEFAULT_THRESHOLD '$LS_FOR_SORT_MERGE_DEFAULT_THRESHOLD; echo -n $'\n\n';

        echo -n $'#define LS_FOR_SORT_MERGE_DEFAULT_ARCH_SECOND_LEVEL '$LS_FOR_SORT_MERGE_DEFAULT_ARCH_SECOND_LEVEL; echo -n $'\n\n';

        echo -n $'#define LS_FOR_SORT_MERGE_MIN_SORTING_SIZE '$LS_FOR_SORT_MERGE_MIN_SORTING_SIZE; echo -n $'\n\n';

	     if [ "$LS_FOR_SORT_MERGE_IMV_AVX" = 1 ]; then
 	         echo -n $'#define LS_FOR_SORT_MERGE_IMV_AVX '$LS_FOR_SORT_MERGE_IMV_AVX; echo -n $'\n\n';
	     fi
	     if [ "$LS_FOR_SORT_MERGE_IMV_AVX_MINOR_BCKTS" = 1 ]; then
 	         echo -n $'#define LS_FOR_SORT_MERGE_IMV_AVX_MINOR_BCKTS '$LS_FOR_SORT_MERGE_IMV_AVX_MINOR_BCKTS; echo -n $'\n\n';
	     fi
	     if [ "$LS_FOR_SORT_MERGE_PREFETCH_INPUT_FOR_MINOR_BCKTS" = 1 ]; then
 	         echo -n $'#define LS_FOR_SORT_MERGE_PREFETCH_INPUT_FOR_MINOR_BCKTS '$LS_FOR_SORT_MERGE_PREFETCH_INPUT_FOR_MINOR_BCKTS; echo -n $'\n\n';
	     fi
        if [ "$LS_FOR_SORT_MERGE_PREFETCH_MINOR_BCKT_SIZES_OFF" = 1 ]; then
 	         echo -n $'#define LS_FOR_SORT_MERGE_PREFETCH_MINOR_BCKT_SIZES_OFF '$LS_FOR_SORT_MERGE_PREFETCH_MINOR_BCKT_SIZES_OFF; echo -n $'\n\n';
	     fi
	     if [ "$LS_FOR_SORT_MERGE_PREFETCH_SLOPES_AND_INTERCEPTS_MINOR_BCKTS" = 1 ]; then
 	         echo -n $'#define LS_FOR_SORT_MERGE_PREFETCH_SLOPES_AND_INTERCEPTS_MINOR_BCKTS '$LS_FOR_SORT_MERGE_PREFETCH_SLOPES_AND_INTERCEPTS_MINOR_BCKTS; echo -n $'\n\n';
	     fi

        echo -n $'#define LS_FOR_SORT_MERGE_WORDSIZE '$LS_FOR_SORT_MERGE_WORDSIZE; echo -n $'\n\n';

        echo -n $'#define LS_FOR_SORT_MERGE_VECTOR_SCALE '$LS_FOR_SORT_MERGE_VECTOR_SCALE; echo -n $'\n\n';

        echo -n $'#define LS_FOR_SORT_MERGE_MAX_VECTOR_SCALE '$LS_FOR_SORT_MERGE_MAX_VECTOR_SCALE; echo -n $'\n\n';

        echo -n $'#define LS_FOR_SORT_MERGE_SIMDStateSize '$LS_FOR_SORT_MERGE_SIMDStateSize; echo -n $'\n\n';

        echo -n $'#define LS_FOR_SORT_MERGE_PDIS '$LS_FOR_SORT_MERGE_PDIS; echo -n $'\n\n';

        )" > $(dirname "$0")/../../include/configs/base_configs.h

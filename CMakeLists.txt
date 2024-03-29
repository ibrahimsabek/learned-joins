#cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
#cmake_minimum_required(VERSION 2.8...3.20.2)
cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

if(UNIX AND NOT APPLE)
  set(LINUX true)
else()
  set(LINUX )
endif()

if(LINUX)
set(CMAKE_C_COMPILER gcc) #gcc 
set(CMAKE_CXX_COMPILER g++) #g++
ADD_DEFINITIONS(-D_LIBCPP_VERSION=1) # tbb clang compatibility
else()
set(CMAKE_C_COMPILER gcc) #gcc
set(CMAKE_CXX_COMPILER g++) #g++
endif()
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

OPTION(VECTORWISE_BRANCHING "Use branching vectorwise primitives" OFF)
OPTION(AUTOVECTORIZE "Allow the compiler to autovectorize" OFF)
OPTION(DATADIR "Directory containing testdata" "")
#OPTION(AVX512EXPERIMENTS ON)
OPTION(HARDWARE_BENCHMARKS ON)
OPTION(INTERPRET_SEPARATE OFF)
OPTION(HASH_SIZE_32 OFF)
#OPTION(KNL OFF)


#-mtune=native -g 
# Compiler flags for the different targets
#if(KNL)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -mtune=native -std=c++17 -g -fPIC -w -Wextra -Wno-psabi  -lnuma -fno-omit-frame-pointer -Wno-unknown-pragmas -mavx512f -mavx512pf -mavx512er -mavx512cd -mavx2 -mbmi2 -mavx -fpermissive")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -w -Wextra -lnuma -fno-omit-frame-pointer -march=native -fdiagnostics-color -mavx512f -mavx512pf -mavx512er -mavx512cd -mavx2 -mbmi2 -mavx -fpermissive")
#else()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=skylake-avx512 -g -mtune=native -std=c++17 -fPIC -w -Wextra -Wno-psabi -lnuma -fno-omit-frame-pointer -Wno-unknown-pragmas -mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx2 -mbmi2 -mavx -fpermissive -lpthread -ffast-math")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -w -Wextra -lnuma -fno-omit-frame-pointer -march=native -fdiagnostics-color -mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx2 -mbmi2 -mavx -fpermissive -lpthread")
#endif(KNL)

if(LINUX)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -static-libasan")
  # set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ")
else()
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -fsanitize=address")
  # set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ")
endif(LINUX)


string(REPLACE "-O2" "" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
string(REPLACE "-O2" "" CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}  -O3")

string(REPLACE "-O2" "" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3")
#set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -fopt-info -fopt-info-missed -O3")
string(REPLACE "-O2" "" CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
set(CMAKE_CFLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -O3")

set(CMAKE_MACOSX_RPATH 1)

IF(VECTORWISE_BRANCHING)
    ADD_DEFINITIONS(-DBRANCHING)
ENDIF(VECTORWISE_BRANCHING)
IF(INTERPRET_SEPARATE)
  ADD_DEFINITIONS(-DINTERPRET_SEPARATE)
ENDIF(INTERPRET_SEPARATE)
IF(HASH_SIZE_32)
  ADD_DEFINITIONS(-DHASH_SIZE=32)
ELSE()
  ADD_DEFINITIONS(-DHASH_SIZE=64)
ENDIF()

#IF(KNL)
#	ADD_DEFINITIONS(-DKNL=1)
#ELSE()
#	ADD_DEFINITIONS(-DKNL=0)
#ENDIF()



IF(AUTOVECTORIZE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize")
ELSE(AUTOVECTORIZE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-tree-vectorize")
ENDIF(AUTOVECTORIZE)


# Download and unpack googletest at configure time
#configure_file(3rdparty/GTest.CMakeLists.txt googletest-download/CMakeLists.txt)
#execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
#    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download" )
#execute_process(COMMAND "${CMAKE_COMMAND}" --build .
#    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download" )

# Prevent GoogleTest from overriding our compiler/linker options
# when building with Visual Studio
#set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This adds
# the following targets: gtest, gtest_main, gmock
# and gmock_main
#add_subdirectory("${CMAKE_BINARY_DIR}/googletest-src"
#                 "${CMAKE_BINARY_DIR}/googletest-build")

# The gtest/gmock targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
#if(CMAKE_VERSION VERSION_LESS 2.8.11)
#    include_directories("${gtest_SOURCE_DIR}/include"
#                        "${gmock_SOURCE_DIR}/include")
#endif()

# Now simply link your own targets against gtest, gmock,
# etc. as appropriate

project("Learned In-memory Joins")

IF(DATADIR)
    message(STATUS "setting datadir to user defined")
    message(STATUS "${DATADIR}")

    ADD_DEFINITIONS(-DDATADIR="${DATADIR}")
ELSE(DATADIR)
    message(STATUS "setting datadir to proj dir ")
    message(STATUS "${PROJECT_SOURCE_DIR}/data")
    ADD_DEFINITIONS(-DDATADIR="${PROJECT_SOURCE_DIR}/data")
ENDIF(DATADIR)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/3rdparty/")
find_package(TBB)
include_directories(${TBB_INCLUDE_DIRS})


IF(LINUX)
  add_subdirectory(3rdparty/jevents)
  set(JEVENTSLIB "jevents")
ELSE(LINUX)
  set(JEVENTSLIB "")
ENDIF(LINUX)

# ==== General Dependencies ====
include(${PROJECT_SOURCE_DIR}/thirdparty/hashing.cmake)
include(${PROJECT_SOURCE_DIR}/thirdparty/learned_hashing.cmake)
include(${PROJECT_SOURCE_DIR}/thirdparty/exotic_hashing.cmake)
include(${PROJECT_SOURCE_DIR}/thirdparty/hashtable.cmake)

add_library(common
#  src/common/algebra/Operators.cpp
#  src/common/algebra/Types.cpp
#  src/common/runtime/Database.cpp
#  src/common/runtime/MemoryPool.cpp
#  src/common/runtime/Types.cpp
#  src/common/runtime/String.cpp
#  src/common/runtime/Import.cpp
#  src/common/runtime/Hashmap.cpp
#  src/common/runtime/Concurrency.cpp
  src/common/runtime/profile.cpp
  )
target_include_directories(common PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src)

#FILE(GLOB PRIMITIVES ${CMAKE_CURRENT_SOURCE_DIR}/src/vectorwise/primitives/*.cpp)
#add_library(vectorwise
#  src/vectorwise/Operations.cpp
#  src/vectorwise/Operators.cpp
#  ${PRIMITIVES}
#  src/vectorwise/QueryBuilder.cpp
#  )
#target_include_directories(vectorwise PUBLIC
#    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#    $<INSTALL_INTERFACE:include>
#    PRIVATE src)

#SET_TARGET_PROPERTIES(vectorwise PROPERTIES
#  COMPILE_FLAGS "-march=skylake-avx512"
#)

#target_link_libraries(vectorwise common)

#add_library(benchmark_config
#  src/benchmarks/config.cpp)
#target_include_directories(benchmark_config PUBLIC
#    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#    $<INSTALL_INTERFACE:include>
#    PRIVATE src)
#target_link_libraries(benchmark_config vectorwise common)


# imv labrary
#add_library(imv
#    src/imv/PipelineTPCH.cpp
#  )
#target_include_directories(imv PUBLIC
#    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#    $<INSTALL_INTERFACE:include>
#    PRIVATE src)
#target_link_libraries(imv common)


#file(GLOB TPCHQUERIES ${CMAKE_CURRENT_SOURCE_DIR}/src/benchmarks/tpch/queries/*.cpp)
#add_library(tpch
#  ${TPCHQUERIES}
#  )
#target_include_directories(tpch PUBLIC
#    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#    $<INSTALL_INTERFACE:include>
#    PRIVATE src)
#target_link_libraries(tpch imv benchmark_config pthread ${TBB_LIBRARIES} vectorwise common)

#add_executable(run_tpch
#  ${HYPER_TRANSLATORS}
#  src/benchmarks/tpch/run.cpp
#  )
#target_include_directories(run_tpch PUBLIC
#    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#    $<INSTALL_INTERFACE:include>
#    PRIVATE src)
#target_link_libraries(run_tpch imv tpch vectorwise common ${JEVENTSLIB})


# Enable tests
enable_testing()
set(CTEST_OUTPUT_ON_FAILURE "1")

FILE(GLOB RMICPPSources ${CMAKE_CURRENT_SOURCE_DIR}/include/rmi/*.cpp)
add_library(rmi_lib ${RMICPPSources})

target_include_directories(rmi_lib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src)

add_executable(non_partition_join_runner src/benchmark_runners/non_partition_join_runner.cpp)
target_include_directories(non_partition_join_runner PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src)
  target_link_libraries(non_partition_join_runner common ${JEVENTSLIB} ${HASHING_LIBRARY} ${LEARNED_HASHING_LIBRARY} ${EXOTIC_HASHING_LIBRARY} ${HASHTABLE_LIBRARY})


add_executable(partition_join_runner src/benchmark_runners/partition_join_runner.cpp)
target_include_directories(partition_join_runner PUBLIC
$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
$<INSTALL_INTERFACE:include>
PRIVATE src)
target_link_libraries(partition_join_runner common ${JEVENTSLIB} ${HASHING_LIBRARY} ${LEARNED_HASHING_LIBRARY} ${EXOTIC_HASHING_LIBRARY} ${HASHTABLE_LIBRARY})

add_executable(non_learned_non_imv_sort_join_runner src/benchmark_runners/non_learned_non_imv_sort_join_runner.cpp)
target_include_directories(non_learned_non_imv_sort_join_runner PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src)
  target_link_libraries(non_learned_non_imv_sort_join_runner common ${JEVENTSLIB} ${HASHING_LIBRARY} ${LEARNED_HASHING_LIBRARY} ${EXOTIC_HASHING_LIBRARY} ${HASHTABLE_LIBRARY})

add_executable(learned_imv_sort_join_runner src/benchmark_runners/learned_imv_sort_join_runner.cpp)
target_include_directories(learned_imv_sort_join_runner PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src)
  target_link_libraries(learned_imv_sort_join_runner common ${JEVENTSLIB} ${HASHING_LIBRARY} ${LEARNED_HASHING_LIBRARY} ${EXOTIC_HASHING_LIBRARY} ${HASHTABLE_LIBRARY})

add_executable(non_imv_indexed_nested_loop_join_runner src/benchmark_runners/non_imv_indexed_nested_loop_join_runner.cpp)
target_include_directories(non_imv_indexed_nested_loop_join_runner PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src)  
  target_link_libraries(non_imv_indexed_nested_loop_join_runner common ${JEVENTSLIB} rmi_lib ${HASHING_LIBRARY} ${LEARNED_HASHING_LIBRARY} ${EXOTIC_HASHING_LIBRARY} ${HASHTABLE_LIBRARY})

#add_executable(avx512_benchmarker src/benchmark_runners/avx512_benchmarker.cpp)
#target_include_directories(avx512_benchmarker PUBLIC
#  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#  $<INSTALL_INTERFACE:include>
#  PRIVATE src)  
#  target_link_libraries(avx512_benchmarker common ${JEVENTSLIB})

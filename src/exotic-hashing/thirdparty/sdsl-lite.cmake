include(ExternalProject)
find_package(Git REQUIRED)

# library name
set(SDSL_LIBRARY sdsl-lite)

ExternalProject_Add(
        ${SDSL_LIBRARY}_src
        PREFIX external/${SDSL_LIBRARY}
        GIT_REPOSITORY "https://github.com/DominikHorn/sdsl-lite.git"
        GIT_TAG 3b18ab2
        TIMEOUT 10
        CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external/${SDSL_LIBRARY}
        -DCMAKE_INSTALL_LIBDIR=${PROJECT_BINARY_DIR}/external/${SDSL_LIBRARY}/lib
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
        -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
        -DBENCHMARK_ENABLE_GTEST_TESTS=0
        UPDATE_COMMAND ""
)

# path to installed artifacts
ExternalProject_Get_Property(${SDSL_LIBRARY}_src install_dir)
set(SDSL_INCLUDE_DIR ${install_dir}/include)
set(SDSL_LIBRARY_PATH ${install_dir}/lib)

# build library from external project
file(MAKE_DIRECTORY ${SDSL_INCLUDE_DIR})
add_library(${SDSL_LIBRARY}::main UNKNOWN IMPORTED)
set_target_properties(${SDSL_LIBRARY}::main PROPERTIES
  IMPORTED_LOCATION ${SDSL_LIBRARY_PATH}/libsdsl.a
  INTERFACE_INCLUDE_DIRECTORIES ${SDSL_INCLUDE_DIR}
  )
add_library(${SDSL_LIBRARY}::divsufsort UNKNOWN IMPORTED)
set_target_properties(${SDSL_LIBRARY}::divsufsort PROPERTIES
  IMPORTED_LOCATION ${SDSL_LIBRARY_PATH}/libdivsufsort.a
  INTERFACE_INCLUDE_DIRECTORIES ${SDSL_INCLUDE_DIR}
  )
add_library(${SDSL_LIBRARY}::divsufsort64 UNKNOWN IMPORTED)
set_target_properties(${SDSL_LIBRARY}::divsufsort64 PROPERTIES
  IMPORTED_LOCATION ${SDSL_LIBRARY_PATH}/libdivsufsort64.a
  INTERFACE_INCLUDE_DIRECTORIES ${SDSL_INCLUDE_DIR}
  )

add_library(${SDSL_LIBRARY} INTERFACE IMPORTED)
set_property(TARGET ${SDSL_LIBRARY} PROPERTY
  INTERFACE_LINK_LIBRARIES ${SDSL_LIBRARY}::main ${SDSL_LIBRARY}::divsufsort ${SDSL_LIBRARY}::divsufsort64
  )
add_dependencies(${SDSL_LIBRARY} ${SDSL_LIBRARY}_src)

message(STATUS "[SDSL] settings")
message(STATUS "    SDSL_LIBRARY = ${SDSL_LIBRARY}")
message(STATUS "    SDSL_INCLUDE_DIR = ${SDSL_INCLUDE_DIR}")
message(STATUS "    SDSL_LIBRARY_PATH = ${SDSL_LIBRARY_PATH}")

#==============================================================================
# Preamble
#==============================================================================
cmake_minimum_required(VERSION 3.24)
project(Numerixx)

# Add the cmake folder so the FindSphinx module is found
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

#==============================================================================
# Set compiler settings
#==============================================================================
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_definitions(NOMINMAX)
message(STATUS "Using toolchain file: ${CMAKE_TOOLCHAIN_FILE}")
#==============================================================================
# Add individual subdirectories
#==============================================================================
add_subdirectory(numerixx)

option(NUMERIXX_ENABLE_TESTS "Enable unit tests" ${PROJECT_IS_TOP_LEVEL})
if(NUMERIXX_ENABLE_TESTS)
    add_subdirectory(test)
endif()

option(NUMERIXX_ENABLE_DOCS "Enable documentation" ${PROJECT_IS_TOP_LEVEL})
if(NUMERIXX_ENABLE_DOCS)
    add_subdirectory(docs)
endif()

option(NUMERIXX_ENABLE_DEMOS "Enable demo programs" ${PROJECT_IS_TOP_LEVEL})
if(NUMERIXX_ENABLE_DEMOS)
    add_subdirectory(demo)
endif()

#option(NUMERIXX_ENABLE_BENCHMARKS "Enable benchmarks" ${PROJECT_IS_TOP_LEVEL})
#if(NUMERIXX_ENABLE_BENCHMARKS)
#    add_subdirectory(benchmark)
#endif()



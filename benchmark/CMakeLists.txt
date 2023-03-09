#==================================================================================================================
# Preamble
#==================================================================================================================
cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(numerixx.benchmark)

#==============================================================================
# Set compiler settings
#==============================================================================
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(BENCHMARK_ENABLE_GTEST_TESTS OFF)
set(BENCHMARK_ENABLE_TESTING OFF)
add_subdirectory(gbench)


add_executable(NumerixxBench)
target_sources(NumerixxBench
        PRIVATE
        benchmark.cpp
        )

target_link_libraries(NumerixxBench PRIVATE benchmark::benchmark numerixx)
cmake_minimum_required(VERSION 3.25)
project(Numerixx VERSION 2.0.0 LANGUAGES CXX)

add_library(numerixx INTERFACE)
add_library(numerixx::numerixx ALIAS numerixx)
target_include_directories(numerixx INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>)
target_compile_features(numerixx INTERFACE cxx_std_23)

option(NUMERIXX_BUILD_TESTS "Build the Numerixx tests" ${PROJECT_IS_TOP_LEVEL})
if (NUMERIXX_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()

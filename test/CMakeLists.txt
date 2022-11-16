#=======================================================================================================================
# Preamble
#=======================================================================================================================
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(Numerix.Tests)

#=======================================================================================================================
# Set up Catch2 library
#=======================================================================================================================
#add_library(Catch INTERFACE IMPORTED)
#target_include_directories(Catch SYSTEM INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/catch/>)
add_subdirectory(Catch2)

#=======================================================================================================================
# Define TEST target
#=======================================================================================================================
add_executable(NumerixTests EXCLUDE_FROM_ALL "")
target_sources(NumerixTests
        PRIVATE
        testMatrix.cpp
        )

target_link_libraries(NumerixTests
        PRIVATE
        numerix
        Catch2WithMain
        )
#==============================================================================
# Preamble
#==============================================================================
cmake_minimum_required(VERSION 3.20)
project(Expected)

#==============================================================================
# Set compiler settings
#==============================================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(expected INTERFACE)
target_include_directories(expected INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/expected>)
#==============================================================================
# Preamble
#==============================================================================
cmake_minimum_required(VERSION 3.20)
project(Numerix)

#==============================================================================
# Set compiler settings
#==============================================================================
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(numerix INTERFACE IMPORTED)
target_include_directories(numerix SYSTEM INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/numerix>)

#add_subdirectory(numerix)
add_subdirectory(demo)
add_subdirectory(test)


#==============================================================================
# Preamble
#==============================================================================
cmake_minimum_required(VERSION 3.24)
project(Numerixx.Lib)

#==============================================================================
# Set compiler settings
#==============================================================================
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library( blaze INTERFACE IMPORTED)
target_include_directories( blaze SYSTEM INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/.dependencies>)

add_library(numerixx INTERFACE)
target_include_directories(numerixx INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)
target_link_libraries(numerixx INTERFACE blaze)
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

find_path(GCEM_INCLUDE_DIR NAMES gcem.hpp)
find_path(EXPECTED_INCLUDE_DIR NAMES expected.hpp PATH_SUFFIXES tl)

#add_library( blaze INTERFACE IMPORTED)
#target_include_directories( blaze SYSTEM INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/.dependencies>)

find_package(blaze REQUIRED)

add_library(numerixx INTERFACE)
target_include_directories(numerixx INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)
target_include_directories(numerixx SYSTEM INTERFACE ${GCEM_INCLUDE_DIR} ${EXPECTED_INCLUDE_DIR})
target_link_libraries(numerixx INTERFACE blaze::blaze)

find_package(tl-expected CONFIG REQUIRED)
target_link_libraries(numerixx INTERFACE tl::expected)
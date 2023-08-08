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

#==============================================================================
# Find dependencies
#==============================================================================
find_path(GCEM_INCLUDE_DIR NAMES gcem.hpp)
find_path(EXPECTED_INCLUDE_DIR NAMES expected.hpp PATH_SUFFIXES tl)

find_package(blaze REQUIRED)
find_package(LAPACK)

#==============================================================================
# Add libraries
#==============================================================================
add_library(nxx_func INTERFACE)
target_include_directories(nxx_func INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/func)
target_include_directories(nxx_func SYSTEM INTERFACE ${GCEM_INCLUDE_DIR} ${EXPECTED_INCLUDE_DIR})
add_library(numerixx::func ALIAS nxx_func)

add_library(nxx_poly INTERFACE)
target_include_directories(nxx_poly INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/poly)
target_include_directories(nxx_poly SYSTEM INTERFACE ${GCEM_INCLUDE_DIR} ${EXPECTED_INCLUDE_DIR})
add_library(numerixx::poly ALIAS nxx_poly)

add_library(nxx_deriv INTERFACE)
target_include_directories(nxx_deriv INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/deriv)
target_include_directories(nxx_deriv SYSTEM INTERFACE ${GCEM_INCLUDE_DIR} ${EXPECTED_INCLUDE_DIR})
target_link_libraries(nxx_deriv INTERFACE nxx_poly)
add_library(numerixx::deriv ALIAS nxx_deriv)

add_library(nxx_roots INTERFACE)
target_include_directories(nxx_roots INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/roots)
target_include_directories(nxx_roots SYSTEM INTERFACE ${GCEM_INCLUDE_DIR} ${EXPECTED_INCLUDE_DIR})
target_link_libraries(nxx_roots INTERFACE nxx_deriv)
add_library(numerixx::roots ALIAS nxx_roots)

add_library(nxx_multiroots INTERFACE)
target_include_directories(nxx_multiroots INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/multiroots)
target_include_directories(nxx_multiroots SYSTEM INTERFACE ${GCEM_INCLUDE_DIR} ${EXPECTED_INCLUDE_DIR})
target_link_libraries(nxx_multiroots INTERFACE LAPACK::LAPACK nxx_deriv blaze::blaze)
add_library(numerixx::multiroots ALIAS nxx_multiroots)



#add_library(numerixx INTERFACE)
#target_include_directories(numerixx INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)
#target_include_directories(numerixx SYSTEM INTERFACE ${GCEM_INCLUDE_DIR} ${EXPECTED_INCLUDE_DIR})
#target_link_libraries(numerixx INTERFACE blaze::blaze)

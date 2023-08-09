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
find_package(gcem REQUIRED)
find_package(tl-expected REQUIRED)
find_package(blaze REQUIRED)
find_package(LAPACK)

#==============================================================================
# Add libraries
#==============================================================================
add_library(nxx_func INTERFACE)
target_include_directories(nxx_func INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/func)
target_link_libraries(nxx_func INTERFACE gcem tl::expected)
add_library(numerixx::func ALIAS nxx_func)

add_library(nxx_poly INTERFACE)
target_include_directories(nxx_poly INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/poly)
target_link_libraries(nxx_poly INTERFACE gcem tl::expected)
add_library(numerixx::poly ALIAS nxx_poly)

add_library(nxx_deriv INTERFACE)
target_include_directories(nxx_deriv INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/deriv)
target_link_libraries(nxx_deriv INTERFACE gcem tl::expected)
target_link_libraries(nxx_deriv INTERFACE nxx_poly)
add_library(numerixx::deriv ALIAS nxx_deriv)

add_library(nxx_roots INTERFACE)
target_include_directories(nxx_roots INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/roots)
target_link_libraries(nxx_roots INTERFACE gcem tl::expected)
target_link_libraries(nxx_roots INTERFACE nxx_deriv)
add_library(numerixx::roots ALIAS nxx_roots)

add_library(nxx_multiroots INTERFACE)
target_include_directories(nxx_multiroots INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/multiroots)
target_link_libraries(nxx_multiroots INTERFACE gcem tl::expected)
target_link_libraries(nxx_multiroots INTERFACE LAPACK::LAPACK nxx_deriv blaze::blaze)
add_library(numerixx::multiroots ALIAS nxx_multiroots)

add_library(nxx_all INTERFACE)
target_link_libraries(nxx_multiroots INTERFACE nxx_func nxx_poly nxx_deriv nxx_roots nxx_multiroots)
add_library(numerixx::all ALIAS nxx_all)

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
find_package(Boost REQUIRED)

set(FETCHCONTENT_SOURCE_DIR_HWINFO ${CMAKE_CURRENT_LIST_DIR}/../../hwinfo)
include(FetchContent)
FetchContent_Declare(hwinfo
        GIT_REPOSITORY https://github.com/troldal/hwinfo.git
        GIT_TAG        header_only_dev
)
FetchContent_MakeAvailable(hwinfo)

#==============================================================================
# Define utility interface library
#==============================================================================
add_library(nxx_utility INTERFACE)
target_include_directories(nxx_utility INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/.utils)

#==============================================================================
# nxx::func library
#==============================================================================
add_library(nxx_func INTERFACE)
target_include_directories(nxx_func INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/func)
add_library(numerixx::func ALIAS nxx_func)

#==============================================================================
# nxx::poly library
#==============================================================================
add_library(nxx_poly INTERFACE)
target_include_directories(nxx_poly INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/poly)
add_library(numerixx::poly ALIAS nxx_poly)

#==============================================================================
# nxx::deriv library
#==============================================================================
add_library(nxx_deriv INTERFACE)
target_include_directories(nxx_deriv INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/deriv)
add_library(numerixx::deriv ALIAS nxx_deriv)

#==============================================================================
# nxx::roots library
#==============================================================================
add_library(nxx_roots INTERFACE)
target_include_directories(nxx_roots INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/roots)
add_library(numerixx::roots ALIAS nxx_roots)

#==============================================================================
# nxx::multiroots library
#==============================================================================
add_library(nxx_multiroots INTERFACE)
target_include_directories(nxx_multiroots INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/multiroots)
add_library(numerixx::multiroots ALIAS nxx_multiroots)

#==============================================================================
# nxx::all library
#==============================================================================
add_library(nxx_all INTERFACE)
add_library(numerixx::all ALIAS nxx_all)

#==============================================================================
# Link libraries
#==============================================================================
target_link_libraries(nxx_func INTERFACE nxx_utility gcem tl::expected)
target_link_libraries(nxx_deriv INTERFACE nxx_utility gcem tl::expected)
target_link_libraries(nxx_poly INTERFACE nxx_utility nxx_deriv nxx_roots gcem tl::expected)
target_link_libraries(nxx_roots INTERFACE nxx_utility nxx_poly nxx_deriv gcem tl::expected)
target_link_libraries(nxx_multiroots INTERFACE nxx_utility nxx_deriv gcem tl::expected)
target_link_libraries(nxx_multiroots INTERFACE LAPACK::LAPACK blaze::blaze)
target_link_libraries(nxx_all INTERFACE nxx_func nxx_poly nxx_deriv nxx_roots nxx_multiroots)

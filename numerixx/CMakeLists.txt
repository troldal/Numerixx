#==============================================================================
# PREAMBLE
#==============================================================================
cmake_minimum_required(VERSION 3.26)
project(Numerixx.Lib CXX)

#==============================================================================
# SET COMPILER SETTINGS
#==============================================================================
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#==============================================================================
# FIND DEPENDENCIES
#==============================================================================
add_subdirectory(.external/blaze SYSTEM)
find_package(Boost REQUIRED)
find_package(tl-expected CONFIG REQUIRED)

#==============================================================================
# DEFINE NUMERIXX LIBRARIES
#==============================================================================

# Define utility interface library
add_library(nxx_utility INTERFACE)
target_include_directories(nxx_utility INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/.utils)

# nxx::poly library
add_library(nxx_poly INTERFACE)
target_include_directories(nxx_poly INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/poly)
add_library(numerixx::poly ALIAS nxx_poly)

# nxx::deriv library
add_library(nxx_deriv INTERFACE)
target_include_directories(nxx_deriv INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/deriv)
add_library(numerixx::deriv ALIAS nxx_deriv)

# nxx::integrate library
add_library(nxx_integrate INTERFACE)
target_include_directories(nxx_integrate INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/integrate)
add_library(numerixx::integrate ALIAS nxx_integrate)

# nxx::interpolate library
add_library(nxx_interpolate INTERFACE)
target_include_directories(nxx_interpolate INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/interpolate)
add_library(numerixx::interpolate ALIAS nxx_interpolate)

# nxx::optimize library
add_library(nxx_optimize INTERFACE)
target_include_directories(nxx_optimize INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/optimize)
add_library(numerixx::optimize ALIAS nxx_optimize)

# nxx::roots library
add_library(nxx_roots INTERFACE)
target_include_directories(nxx_roots INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/roots)
add_library(numerixx::roots ALIAS nxx_roots)

# nxx::multiroots library
add_library(nxx_multiroots INTERFACE)
target_include_directories(nxx_multiroots INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/multiroots)
add_library(numerixx::multiroots ALIAS nxx_multiroots)

# nxx::all library (all libraries)
add_library(nxx_all INTERFACE)
add_library(numerixx::all ALIAS nxx_all)

#==============================================================================
# LINK LIBRARIES
#==============================================================================
target_link_libraries(nxx_utility
        INTERFACE
        Boost::boost
        tl::expected)

target_link_libraries(nxx_deriv
        INTERFACE
        nxx_utility )

target_link_libraries(nxx_integrate
        INTERFACE
        nxx_utility )

target_link_libraries(nxx_interpolate
        INTERFACE
        nxx_utility
        nxx_poly
        nxx_deriv
        blaze)

target_link_libraries(nxx_optimize
        INTERFACE
        nxx_utility
        nxx_deriv
        nxx_poly
        nxx_interpolate
        nxx_roots)

target_link_libraries(nxx_poly
        INTERFACE
        nxx_utility
        nxx_deriv
        nxx_roots)

target_link_libraries(nxx_roots
        INTERFACE
        nxx_utility
        nxx_poly
        nxx_deriv)

target_link_libraries(nxx_multiroots
        INTERFACE
        nxx_utility
        nxx_deriv
        nxx_roots
        nxx_poly
        nxx_interpolate
        blaze)

#if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_SIMULATE_ID STREQUAL "MSVC")
#    find_package(OpenMP REQUIRED)
#    target_link_libraries(nxx_multiroots INTERFACE OpenMP::OpenMP_CXX)
#    target_link_libraries(nxx_interpolate INTERFACE OpenMP::OpenMP_CXX)
#endif()

target_link_libraries(nxx_all
        INTERFACE
        nxx_poly
        nxx_deriv
        nxx_integrate
        nxx_interpolate
        nxx_optimize
        nxx_roots
        nxx_multiroots)

target_compile_definitions(nxx_all
        INTERFACE
        NOMINMAX)

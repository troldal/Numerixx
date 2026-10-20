#==============================================================================
# PREAMBLE
#==============================================================================
cmake_minimum_required(VERSION 3.26)
project(Numerixx.Lib CXX)

#==============================================================================
# SET COMPILER SETTINGS
#==============================================================================
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#==============================================================================
# FIND DEPENDENCIES
#==============================================================================
find_package(Boost REQUIRED)

#==============================================================================
# DEFINE NUMERIXX LIBRARIES
#==============================================================================

if (NUMERIXX_MULTIROOTS OR NUMERIXX_OPTIMIZE)
    set(NUMERIXX_INTERPOLATE ON)
endif ()

if (NUMERIXX_INTERPOLATE OR NUMERIXX_MULTIROOTS OR NUMERIXX_ALL)
    add_subdirectory(.external/blaze SYSTEM)
endif ()


# Define utility interface library
add_library(nxx_utility INTERFACE)
target_include_directories(nxx_utility INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/.utils)

# nxx::poly library
if (NUMERIXX_POLY OR NUMERIXX_ALL)
    include(.cmake/AddPoly.cmake)
    include(.cmake/AddDeriv.cmake)
    include(.cmake/AddRoots.cmake)
    #    message(STATUS "[Numerixx]: Including numerixx::poly library")
    #    add_library(nxx_poly INTERFACE)
    #    target_include_directories(nxx_poly INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/poly)
    #    add_library(numerixx::poly ALIAS nxx_poly)
endif ()

# nxx::deriv library
if (NUMERIXX_DERIV OR NUMERIXX_ALL)
    include(.cmake/AddDeriv.cmake)
    #    message(STATUS "[Numerixx]: Including numerixx::deriv library")
    #    add_library(nxx_deriv INTERFACE)
    #    target_include_directories(nxx_deriv INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/deriv)
    #    add_library(numerixx::deriv ALIAS nxx_deriv)
endif ()

# nxx::integrate library
if (NUMERIXX_INTEGRATE OR NUMERIXX_ALL)
    include(.cmake/AddIntegrate.cmake)
    #    message(STATUS "[Numerixx]: Including numerixx::integrate library")
    #    add_library(nxx_integrate INTERFACE)
    #    target_include_directories(nxx_integrate INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/integrate)
    #    add_library(numerixx::integrate ALIAS nxx_integrate)
endif ()

# nxx::interpolate library
if (NUMERIXX_INTERPOLATE OR NUMERIXX_ALL)
    include(.cmake/AddInterpolate.cmake)
    include(.cmake/AddPoly.cmake)
    include(.cmake/AddDeriv.cmake)
    #    message(STATUS "[Numerixx]: Including numerixx::interpolate library")
    #    add_library(nxx_interpolate INTERFACE)
    #    target_include_directories(nxx_interpolate INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/interpolate)
    #    add_library(numerixx::interpolate ALIAS nxx_interpolate)
endif ()

# nxx::optimize library
if (NUMERIXX_OPTIMIZE OR NUMERIXX_ALL)
    include(.cmake/AddOptimize.cmake)
    include(.cmake/AddPoly.cmake)
    include(.cmake/AddDeriv.cmake)
    include(.cmake/AddInterpolate.cmake)
    include(.cmake/AddRoots.cmake)
    #    message(STATUS "[Numerixx]: Including numerixx::optimize library")
    #    add_library(nxx_optimize INTERFACE)
    #    target_include_directories(nxx_optimize INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/optimize)
    #    add_library(numerixx::optimize ALIAS nxx_optimize)
endif ()

# nxx::roots library
if (NUMERIXX_ROOTS OR NUMERIXX_ALL)
    include(.cmake/AddRoots.cmake)
    include(.cmake/AddPoly.cmake)
    include(.cmake/AddDeriv.cmake)
    #    message(STATUS "[Numerixx]: Including numerixx::roots library")
    #    add_library(nxx_roots INTERFACE)
    #    target_include_directories(nxx_roots INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/roots)
    #    add_library(numerixx::roots ALIAS nxx_roots)
endif ()

# nxx::multiroots library
if (NUMERIXX_MULTIROOTS OR NUMERIXX_ALL)
    include(.cmake/AddMultiroots.cmake)
    include(.cmake/AddRoots.cmake)
    include(.cmake/AddPoly.cmake)
    include(.cmake/AddDeriv.cmake)
    include(.cmake/AddInterpolate.cmake)
    #    message(STATUS "[Numerixx]: Including numerixx::multiroots library")
    #    add_library(nxx_multiroots INTERFACE)
    #    target_include_directories(nxx_multiroots INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/multiroots)
    #    add_library(numerixx::multiroots ALIAS nxx_multiroots)
endif ()

# nxx::all library (all libraries)
add_library(nxx_all INTERFACE)
add_library(numerixx::all ALIAS nxx_all)

#==============================================================================
# LINK LIBRARIES
#==============================================================================
target_link_libraries(nxx_utility
        INTERFACE
        Boost::boost)

if (NUMERIXX_DERIV OR NUMERIXX_ALL)
    target_link_libraries(nxx_deriv
            INTERFACE
            nxx_utility)
endif ()

if (NUMERIXX_INTEGRATE OR NUMERIXX_ALL)
    target_link_libraries(nxx_integrate
            INTERFACE
            nxx_utility)
endif ()

if (NUMERIXX_INTERPOLATE OR NUMERIXX_ALL)
    target_link_libraries(nxx_interpolate
            INTERFACE
            nxx_utility
            nxx_poly
            nxx_deriv
            blaze)
endif ()

if (NUMERIXX_OPTIMIZE OR NUMERIXX_ALL)
    target_link_libraries(nxx_optimize
            INTERFACE
            nxx_utility
            nxx_deriv
            nxx_poly
            nxx_interpolate
            nxx_roots)
endif ()

if (NUMERIXX_POLY OR NUMERIXX_ALL)
    target_link_libraries(nxx_poly
            INTERFACE
            nxx_utility
            nxx_deriv
            nxx_roots)
endif ()

if (NUMERIXX_ROOTS OR NUMERIXX_ALL)
    target_link_libraries(nxx_roots
            INTERFACE
            nxx_utility
            nxx_poly
            nxx_deriv)
endif ()


if (NUMERIXX_MULTIROOTS OR NUMERIXX_ALL)
    target_link_libraries(nxx_multiroots
            INTERFACE
            nxx_utility
            nxx_deriv
            nxx_roots
            nxx_poly
            nxx_interpolate
            blaze)
endif ()


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

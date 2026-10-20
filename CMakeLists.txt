#==============================================================================
# PREAMBLE
#==============================================================================
cmake_minimum_required(VERSION 3.26)
project(Numerixx)

# Add the cmake folder so the FindSphinx module is found
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

#==============================================================================
# SET COMPILER SETTINGS
#==============================================================================
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
add_compile_definitions(NOMINMAX)
message(STATUS "Using toolchain file: ${CMAKE_TOOLCHAIN_FILE}")

#==============================================================================
#
#==============================================================================
option(NUMERIXX_DERIV "Enable the numerixx::deriv library" OFF)
option(NUMERIXX_INTEGRATE "Enable the numerixx::integrate library" OFF)
option(NUMERIXX_OPTIMIZE "Enable the numerixx::optimize library" OFF)
option(NUMERIXX_POLY "Enable the numerixx::poly library" OFF)
option(NUMERIXX_MULTIROOTS "Enable the numerixx::multiroots library" OFF)
option(NUMERIXX_ROOTS "Enable the numerixx::roots library" OFF)
option(NUMERIXX_INTERPOLATE "Enable the numerixx::interpolate library" OFF)
option(NUMERIXX_ALL "Enable all numerixx libraries" OFF)

if (PROJECT_IS_TOP_LEVEL)
    set(NUMERIXX_ALL ON)
endif ()


#==============================================================================
# Add individual subdirectories
#==============================================================================
add_subdirectory(numerixx)

#option(NUMERIXX_ENABLE_TESTS "Enable unit tests" ${PROJECT_IS_TOP_LEVEL})
#if(NUMERIXX_ENABLE_TESTS)
#    add_subdirectory(test)
#endif()

option(NUMERIXX_ENABLE_DOCS "Enable documentation" ${PROJECT_IS_TOP_LEVEL})
if(NUMERIXX_ENABLE_DOCS)
    add_subdirectory(docs)
endif()

option(NUMERIXX_ENABLE_DEMOS "Enable demo programs" ${PROJECT_IS_TOP_LEVEL})
if(NUMERIXX_ENABLE_DEMOS)
    add_subdirectory(demo)
endif()

#option(NUMERIXX_ENABLE_BENCHMARKS "Enable benchmarks" ${PROJECT_IS_TOP_LEVEL})
#if(NUMERIXX_ENABLE_BENCHMARKS)
#    add_subdirectory(benchmark)
#endif()



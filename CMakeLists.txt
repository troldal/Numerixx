#==============================================================================
# Preamble
#==============================================================================
cmake_minimum_required(VERSION 3.20)
project(Numerixx)

# Add the cmake folder so the FindSphinx module is found
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

#==============================================================================
# Set compiler settings
#==============================================================================
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library( blaze INTERFACE IMPORTED)
target_include_directories( blaze SYSTEM INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/numerixx/.dependencies> )

add_library(numerixx INTERFACE)
target_include_directories(numerixx INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/numerixx>)
target_link_libraries(numerixx INTERFACE blaze)

add_subdirectory(docs)
add_subdirectory(demo)
add_subdirectory(test)
add_subdirectory(benchmark)


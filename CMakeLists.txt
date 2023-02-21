#==============================================================================
# Preamble
#==============================================================================
cmake_minimum_required(VERSION 3.20)
project(Numerix)

# Add the cmake folder so the FindSphinx module is found
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

#==============================================================================
# Set compiler settings
#==============================================================================
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(deps)

add_library(numerix INTERFACE)
target_include_directories(numerix INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/numerix>)
target_link_libraries(numerix INTERFACE expected)

#add_subdirectory(numerix)
add_subdirectory(docs)
add_subdirectory(demo)
add_subdirectory(test)


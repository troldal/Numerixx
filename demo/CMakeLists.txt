#==================================================================================================================
# Preamble
#==================================================================================================================
cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(Numerix.Samples)

#==============================================================================
# Set compiler settings
#==============================================================================
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#=================================================================================================================
# Define TEST target
#=================================================================================================================
add_executable(DemoRoot)
target_sources(DemoRoot
        PRIVATE
        DemoRoot.cpp
        )

target_link_libraries(DemoRoot PRIVATE fmt numerix)

add_executable(DemoPoly)
target_sources(DemoPoly
        PRIVATE
        DemoPoly.cpp
        )

target_link_libraries(DemoPoly PRIVATE fmt numerix)

add_executable(DemoMatrix)
target_sources(DemoMatrix
        PRIVATE
        DemoMatrix.cpp
        )

target_link_libraries(DemoMatrix PRIVATE numerix)

add_executable(DemoDeriv)
target_sources(DemoDeriv
        PRIVATE
        DemoDeriv.cpp
        )

target_link_libraries(DemoDeriv PRIVATE numerix)
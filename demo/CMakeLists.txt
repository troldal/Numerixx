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

add_subdirectory(fmt)
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

#=======================================================================================================================
# COMPILER FLAGS
#   Set compiler debug flags for GCC, Clang and MSVC.
#=======================================================================================================================
list(APPEND OPENXLSX_DEBUG_FLAGS_GNU
        "-Wmisleading-indentation"
        "-Wduplicated-cond"
        "-Wduplicated-branches"
        "-Wlogical-op"
        "-Wnull-dereference")
list(APPEND OPENXLSX_DEBUG_FLAGS_GNUCLANG
        "-Wall"
        "-Wextra"
        "-Wshadow"
        "-Wnon-virtual-dtor"
        "-Wold-style-cast"
        "-Wcast-align"
        "-Wunused"
        "-Woverloaded-virtual"
        "-Wpedantic"
        "-Wconversion"
        "-Wdouble-promotion"
        "-Wformat=2"
        "-Weffc++"
        "-Wno-unknown-pragmas")
list(APPEND OPENXLSX_DEBUG_FLAGS_MSVC
        "/permissive"
        "/W4"
        "/w14242"
        "/w14254"
        "/w14263"
        "/w14265"
        "/w14287"
        "/we4289"
        "/w14296"
        "/w14311"
        "/w14545"
        "/w14546"
        "/w14547"
        "/w14549"
        "/w14555"
        "/w14619"
        "/w14640"
        "/w14826"
        "/w14905"
        "/w14906"
        "/w14928"
        "/wd4251"
        "/wd4275")
list(APPEND OPENXLSX_RELEASE_FLAGS_MSVC
        "/wd4251"
        "/wd4275")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options(DemoMatrix PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_GNU}>)

elseif (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR
("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang") OR
("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"))
    target_compile_options(DemoMatrix PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_GNUCLANG}>)

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    string(REGEX REPLACE "/W[3|4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    target_compile_options(DemoMatrix PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_MSVC}>)
    target_compile_options(DemoMatrix PRIVATE $<$<CONFIG:Release>:${OPENXLSX_RELEASE_FLAGS_MSVC}>)
endif ()
#=======================================================================================================================
# Preamble
#=======================================================================================================================
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(Numerixx.Tests)

#=======================================================================================================================
# Set up Catch2 library
#=======================================================================================================================
#add_subdirectory(Catch2)
find_package(Catch2 REQUIRED)

#=======================================================================================================================
# Define TEST target
#=======================================================================================================================
add_executable(NumerixxTests EXCLUDE_FROM_ALL "")
target_sources(NumerixxTests
        PRIVATE
        testPolynomials.cpp
#        testDerivatives.cpp
#        testMatrix.cpp
#        testRootBracketing.cpp
#        testRootPolishing.cpp
        )

target_link_libraries(NumerixxTests
        PRIVATE
        numerixx::poly
        Catch2::Catch2WithMain
        )

list(APPEND NUMERIXX_TEST_DEBUG_FLAGS_GNU
        "-Wmisleading-indentation"
        "-Wduplicated-cond"
        "-Wduplicated-branches"
        "-Wlogical-op"
        "-Wnull-dereference")
list(APPEND NUMERIXX_TEST_DEBUG_FLAGS_GNUCLANG
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
list(APPEND NUMERIXX_TEST_DEBUG_FLAGS_MSVC
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
list(APPEND NUMERIXX_TEST_RELEASE_FLAGS_MSVC
        "/wd4251"
        "/wd4275")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options(NumerixxTests PRIVATE $<$<CONFIG:Debug>:${NUMERIXX_TEST_DEBUG_FLAGS_GNU}>)

elseif (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR
("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang") OR
("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"))
    target_compile_options(NumerixxTests PRIVATE $<$<CONFIG:Debug>:${NUMERIXX_TEST_DEBUG_FLAGS_GNUCLANG}>)

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    string(REGEX REPLACE "/W[3|4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    target_compile_options(NumerixxTests PRIVATE $<$<CONFIG:Debug>:${NUMERIXX_TEST_DEBUG_FLAGS_MSVC}>)
    target_compile_options(NumerixxTests PRIVATE $<$<CONFIG:Release>:${NUMERIXX_TEST_RELEASE_FLAGS_MSVC}>)
endif ()
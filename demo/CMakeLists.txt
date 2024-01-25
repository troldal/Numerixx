#==================================================================================================================
# Preamble
#==================================================================================================================
cmake_minimum_required(VERSION 3.26 FATAL_ERROR)
project(numerixx.samples)

#==============================================================================
# Set compiler settings
#==============================================================================
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#==================================================================================================================
# Dependencies
#==================================================================================================================
include(SetDefaultWarnings)
include(EnableLTO)
find_package(Boost REQUIRED)
find_package(sciplot CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)

#=================================================================================================================
# Define targets
#=================================================================================================================

# Target: Root finding
add_executable(DemoRootFinding DemoRootFinding.cpp)
target_link_libraries(DemoRootFinding
        PRIVATE
        numerixx::roots
        fmt::fmt-header-only)
target_set_warnings(DemoRootFinding)

# Target: Root searching
add_executable(DemoRootSearching DemoRootSearching.cpp)
target_link_libraries(DemoRootSearching
        PRIVATE
        numerixx::roots
        fmt::fmt-header-only)
target_set_warnings(DemoRootSearching)

# Target: Polynomials
add_executable(DemoPoly DemoPoly.cpp)
target_link_libraries(DemoPoly
        PRIVATE
        numerixx::poly
        fmt::fmt-header-only)
target_set_warnings(DemoPoly)

# Target: Numerical derivatives
add_executable(DemoDeriv DemoDeriv.cpp)
target_link_libraries(DemoDeriv
        PRIVATE
        numerixx::deriv
        fmt::fmt-header-only)
target_set_warnings(DemoDeriv)

# Target: Numerical integration
add_executable(DemoIntegrate DemoIntegrate.cpp)
target_link_libraries(DemoIntegrate
        PRIVATE
        numerixx::integrate
        fmt::fmt-header-only)
target_set_warnings(DemoIntegrate)

# Target: Numerical interpolation
add_executable(DemoInterpolate DemoInterpolate.cpp)
target_link_libraries(DemoInterpolate
        PRIVATE
        numerixx::interpolate
        fmt::fmt-header-only
        sciplot::sciplot)
target_set_warnings(DemoInterpolate)
target_enable_lto(DemoInterpolate)

# Target: Numerical optimization
add_executable(DemoOptimize DemoOptimize.cpp)
target_link_libraries(DemoOptimize
        PRIVATE
        numerixx::optimize
        fmt::fmt-header-only
        sciplot::sciplot)
target_set_warnings(DemoOptimize)

# Target: Multiple root-finding
add_executable(DemoMultiroot DemoMultiroot.cpp)
target_link_libraries(DemoMultiroot
        PRIVATE
        numerixx::multiroots
        fmt::fmt-header-only)
target_set_warnings(DemoMultiroot)

#=======================================================================================================================
# COMPILER FLAGS
#   Set compiler debug flags for GCC, Clang and MSVC.
#=======================================================================================================================
#list(APPEND OPENXLSX_DEBUG_FLAGS_GNU
#        "-Wmisleading-indentation"
#        "-Wduplicated-cond"
#        "-Wduplicated-branches"
#        "-Wlogical-op"
#        "-Wnull-dereference"
#        "-Wall"
#        "-Wextra"
#        "-Wshadow"
#        "-Wnon-virtual-dtor"
#        "-Wold-style-cast"
#        "-Wuseless-cast"
#        "-Wcast-align"
#        "-Wunused"
#        "-Woverloaded-virtual"
#        "-Wpedantic"
#        "-Wconversion"
#        "-Wsign-conversion"
#        "-Wdouble-promotion"
#        "-Wformat=2"
#        "-Wimplicit-fallthrough"
#        "-Weffc++"
#        "-Wno-unknown-pragmas"
#)
#list(APPEND OPENXLSX_DEBUG_FLAGS_GNUCLANG
#        "-Wall"
#        "-Wextra"
#        "-Wshadow"
#        "-Wnon-virtual-dtor"
#        "-Wold-style-cast"
#        "-Wuseless-cast"
#        "-Wcast-align"
#        "-Wunused"
#        "-Woverloaded-virtual"
#        "-Wpedantic"
#        "-Wconversion"
#        "-Wsign-conversion"
#        "-Wmisleading-indentation"
#        "-Wdouble-promotion"
#        "-Wformat=2"
#        "-Wlifetime"
#        "-Wimplicit-fallthrough"
#        "-Weverything"
#        "-Weffc++"
#        "-Wno-unknown-pragmas"
#        "-Wno-c++98-compat"
#        "-Wno-c++98-compat-pedantic"
#        "-Wno-documentation"
#        "-Wno-documentation-unknown-command"
#)
#list(APPEND OPENXLSX_DEBUG_FLAGS_MSVC
#        "/permissive-"
#        "/W4"
#        "/w14242"
#        "/w14254"
#        "/w14263"
#        "/w14265"
#        "/w14287"
#        "/we4289"
#        "/w14296"
#        "/w14311"
#        "/w14545"
#        "/w14546"
#        "/w14547"
#        "/w14549"
#        "/w14555"
#        "/w14619"
#        "/w14640"
#        "/w14826"
#        "/w14905"
#        "/w14906"
#        "/w14928"
#        "/wd4251"
#        "/wd4275"
#        "/guard:cf")

#if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#    target_compile_options(DemoRootFinding PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_GNU}>)
#    target_compile_options(DemoRootSearching PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_GNU}>)
#    target_compile_options(DemoPoly PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_GNU}>)
#    target_compile_options(DemoDeriv PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_GNU}>)
#    target_compile_options(DemoMultiroot PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_GNU}>)
#
#
#elseif (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR
#("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang"))
#    target_compile_options(DemoRootFinding PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_GNUCLANG}>)
#    target_compile_options(DemoRootSearching PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_GNUCLANG}>)
#    target_compile_options(DemoPoly PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_GNUCLANG}>)
#    target_compile_options(DemoDeriv PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_GNUCLANG}>)
#    target_compile_options(DemoMultiroot PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_GNUCLANG}>)
#
#elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
#    string(REGEX REPLACE "/W[3|4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
#    target_compile_options(DemoRootFinding PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_MSVC}>)
#    target_compile_options(DemoRootSearching PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_MSVC}>)
#    target_compile_options(DemoPoly PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_MSVC}>)
#    target_compile_options(DemoDeriv PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_MSVC}>)
#    target_compile_options(DemoMultiroot PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_MSVC}>)
#endif ()
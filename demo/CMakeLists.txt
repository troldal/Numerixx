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

find_package(Boost REQUIRED)
#find_package(fmt REQUIRED)

include(FetchContent)
FetchContent_Declare(fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG master
)
FetchContent_MakeAvailable(fmt)



#=================================================================================================================
# Define targets
#=================================================================================================================

# Target: Root finding
add_executable(DemoRootFinding)
target_sources(DemoRootFinding
        PRIVATE
        DemoRootFinding.cpp
        )

target_link_libraries(DemoRootFinding PRIVATE numerixx::roots fmt::fmt)

# Target: Root searching
add_executable(DemoRootSearching)
target_sources(DemoRootSearching
        PRIVATE
        DemoRootSearching.cpp
        )

target_link_libraries(DemoRootSearching PRIVATE numerixx::roots fmt::fmt)

# Target: Function
add_executable(DemoFunction)
target_sources(DemoFunction
        PRIVATE
        DemoFunction.cpp
        )

target_link_libraries(DemoFunction PRIVATE numerixx::func fmt::fmt)


# Target: Polynomials
add_executable(DemoPoly)
target_sources(DemoPoly
        PRIVATE
        DemoPoly.cpp
        )

target_link_libraries(DemoPoly PRIVATE numerixx::poly fmt::fmt)

# Target: Numerical derivatives
add_executable(DemoDeriv)
target_sources(DemoDeriv
        PRIVATE
        DemoDeriv.cpp
        )

target_link_libraries(DemoDeriv PRIVATE numerixx::deriv fmt::fmt)

# Target: Numerical integration
add_executable(DemoIntegrate)
target_sources(DemoIntegrate
        PRIVATE
        DemoIntegrate.cpp
)

target_link_libraries(DemoIntegrate PRIVATE numerixx::integrate fmt::fmt)

# Target: Numerical interpolation
add_executable(DemoInterpolate)
target_sources(DemoInterpolate
        PRIVATE
        DemoInterpolate.cpp
)

target_link_libraries(DemoInterpolate PRIVATE numerixx::interpolate fmt::fmt)

# Target: Numerical optimization
add_executable(DemoOptimize)
target_sources(DemoOptimize
        PRIVATE
        DemoOptimize.cpp
)

target_link_libraries(DemoOptimize PRIVATE numerixx::optimize fmt::fmt)


# Target: Multiple root-finding
add_executable(DemoMultiroot)
target_sources(DemoMultiroot
        PRIVATE
        DemoMultiroot.cpp
        )

target_link_libraries(DemoMultiroot PRIVATE fmt numerixx::multiroots fmt::fmt)

#=======================================================================================================================
# COMPILER FLAGS
#   Set compiler debug flags for GCC, Clang and MSVC.
#=======================================================================================================================
list(APPEND OPENXLSX_DEBUG_FLAGS_GNU
        "-Wmisleading-indentation"
        "-Wduplicated-cond"
        "-Wduplicated-branches"
        "-Wlogical-op"
        "-Wnull-dereference"
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
        "-Wno-unknown-pragmas"
)
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
        "-Wno-unknown-pragmas"
        "-Wno-c++98-compat"
        "-Wno-c++98-compat-pedantic"
        "-Wno-documentation"
        "-Wno-documentation-unknown-command"
)
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
        "/wd4275"
        "/guard:cf")
list(APPEND OPENXLSX_RELEASE_FLAGS_MSVC
        "/wd4251"
        "/wd4275")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options(DemoRootFinding PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_GNU}>)
    target_compile_options(DemoRootSearching PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_GNU}>)
    target_compile_options(DemoPoly PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_GNU}>)
    target_compile_options(DemoDeriv PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_GNU}>)
    target_compile_options(DemoMultiroot PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_GNU}>)


elseif (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR
("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang"))
    target_compile_options(DemoRootFinding PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_GNUCLANG}>)
    target_compile_options(DemoRootSearching PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_GNUCLANG}>)
    target_compile_options(DemoPoly PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_GNUCLANG}>)
    target_compile_options(DemoDeriv PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_GNUCLANG}>)
    target_compile_options(DemoMultiroot PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_GNUCLANG}>)

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    string(REGEX REPLACE "/W[3|4]" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    target_compile_options(DemoRootFinding PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_MSVC}>)
    target_compile_options(DemoRootFinding PRIVATE $<$<CONFIG:Release>:${OPENXLSX_RELEASE_FLAGS_MSVC}>)

    target_compile_options(DemoRootSearching PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_MSVC}>)
    target_compile_options(DemoRootSearching PRIVATE $<$<CONFIG:Release>:${OPENXLSX_RELEASE_FLAGS_MSVC}>)

    target_compile_options(DemoPoly PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_MSVC}>)
    target_compile_options(DemoPoly PRIVATE $<$<CONFIG:Release>:${OPENXLSX_RELEASE_FLAGS_MSVC}>)

    target_compile_options(DemoDeriv PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_MSVC}>)
    target_compile_options(DemoDeriv PRIVATE $<$<CONFIG:Release>:${OPENXLSX_RELEASE_FLAGS_MSVC}>)

    target_compile_options(DemoMultiroot PRIVATE $<$<CONFIG:Debug>:${OPENXLSX_DEBUG_FLAGS_MSVC}>)
    target_compile_options(DemoMultiroot PRIVATE $<$<CONFIG:Release>:${OPENXLSX_RELEASE_FLAGS_MSVC}>)
endif ()
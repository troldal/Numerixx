cmake_minimum_required(VERSION 3.25)
project(Numerixx VERSION 2.0.0 LANGUAGES CXX)

add_library(numerixx INTERFACE)
add_library(numerixx::numerixx ALIAS numerixx)
target_include_directories(numerixx INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>)
target_compile_features(numerixx INTERFACE cxx_std_23)

# The tests, demos and examples compile through ccache where it is found: a file whose source, headers and flags are
# unchanged is taken from the cache instead of compiled again, locally and in CI. Top level only; a project that adds
# Numerixx keeps its own choice. MSVC's objects are cacheable only with their debug information inside them (/Z7),
# not in a .pdb shared by the target (/Zi, CMake's default). Not clang-cl: CMake asks it for a dependency file with
# -clang:-MD -clang:-MF, which ccache 4.14 refuses and 4.9 does not restore on a hit, so that Ninja would lose the
# object's headers and not rebuild it when one changes.
option(NUMERIXX_CCACHE "Compile the tests, demos and examples through ccache where it is found" ${PROJECT_IS_TOP_LEVEL})
if (NUMERIXX_CCACHE AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
    message(STATUS "ccache: not with clang-cl, whose dependency files it does not keep")
elseif (NUMERIXX_CCACHE)
    find_program(NUMERIXX_CCACHE_PROGRAM ccache)
    if (NUMERIXX_CCACHE_PROGRAM)
        set(CMAKE_CXX_COMPILER_LAUNCHER ${NUMERIXX_CCACHE_PROGRAM})
        set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<$<CONFIG:Debug,RelWithDebInfo>:Embedded>")
    else ()
        message(STATUS "ccache not found: compiling without it")
    endif ()
endif ()

option(NUMERIXX_BUILD_TESTS "Build the Numerixx tests" ${PROJECT_IS_TOP_LEVEL})
if (NUMERIXX_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()

option(NUMERIXX_BUILD_DEMOS "Build the Numerixx demos" ${PROJECT_IS_TOP_LEVEL})
if (NUMERIXX_BUILD_DEMOS)
    add_subdirectory(demo)
endif ()

# The book on nxx::deriv (docs/deriv-book): its code examples build and run as tests; its data and PDF are
# on-demand targets (numerixx_book_data, numerixx_book).
option(NUMERIXX_BUILD_BOOK "Build the book's examples and add its targets" ${PROJECT_IS_TOP_LEVEL})
if (NUMERIXX_BUILD_BOOK)
    add_subdirectory(docs/deriv-book)
endif ()

option(NUMERIXX_BUILD_DOCS "Add the Numerixx documentation targets (needs Doxygen; the site also uv)"
        ${PROJECT_IS_TOP_LEVEL})
if (NUMERIXX_BUILD_DOCS)
    find_package(Doxygen 1.10)
    if (DOXYGEN_FOUND)
        add_subdirectory(docs)
    else ()
        message(STATUS "Doxygen not found: no documentation targets")
    endif ()
endif ()

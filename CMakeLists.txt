cmake_minimum_required(VERSION 3.25)
project(Numerixx VERSION 2.0.0 LANGUAGES CXX)

add_library(numerixx INTERFACE)
add_library(numerixx::numerixx ALIAS numerixx)
target_include_directories(numerixx INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>)
target_compile_features(numerixx INTERFACE cxx_std_23)

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

option(NUMERIXX_BUILD_DOCS "Add the Numerixx documentation targets (needs Doxygen)" ${PROJECT_IS_TOP_LEVEL})
if (NUMERIXX_BUILD_DOCS)
    find_package(Doxygen 1.10)
    if (DOXYGEN_FOUND)
        add_subdirectory(docs)
    else ()
        message(STATUS "Doxygen not found: no documentation targets")
    endif ()
endif ()

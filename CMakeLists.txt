cmake_minimum_required(VERSION 3.25...4.4)

# CPM runs cmake_minimum_required(3.14) internally, so a plain cmake_policy(SET CMP0168 NEW) would not reach its
# FetchContent calls; the policy *default* does. NEW populates without a sub-build (faster, shorter paths).
set(CMAKE_POLICY_DEFAULT_CMP0168 NEW)

project(Numerixx
        VERSION 2.0.0
        DESCRIPTION "A general-purpose, header-only numerical library for C++23"
        HOMEPAGE_URL "https://github.com/troldal/Numerixx"
        LANGUAGES CXX)

# ---- Options (DESIGN §4.3). Development options are OFF when Numerixx is not the top-level project. -------------
option(NUMERIXX_BUILD_TESTS             "Build the test suite"                                   ${PROJECT_IS_TOP_LEVEL})
option(NUMERIXX_BUILD_EXAMPLES          "Build the examples (registered as smoke tests)"          ON)
option(NUMERIXX_BUILD_BENCHMARKS        "Build the google/benchmark benchmarks"                   OFF)
option(NUMERIXX_BUILD_DOCS              "Build the API documentation (needs Doxygen)"             OFF)
option(NUMERIXX_BUILD_INTEGRATION_TESTS "Register the consumer-build scenarios as tests"          OFF)
option(NUMERIXX_INSTALL                 "Generate install and export rules"                       ${PROJECT_IS_TOP_LEVEL})
option(NUMERIXX_WITH_FXT                "Fetch FXT and provide numerixx::pipes"                   ON)
option(NUMERIXX_WITH_LINALG             "Fetch Eigen and provide numerixx::linalg and numerixx::multiroots" ON)
option(NUMERIXX_WITH_MULTIPRECISION     "Fetch standalone Boost.Multiprecision and provide numerixx::multiprecision" OFF)
option(NUMERIXX_TEST_ORACLES            "Fetch standalone Boost.Math for the oracle tests"        OFF)
option(NUMERIXX_NO_EXCEPTIONS           "Build tests, examples and benchmarks with exceptions disabled" OFF)
option(NUMERIXX_WARNINGS_AS_ERRORS      "Treat warnings as errors in Numerixx-owned targets"      OFF)
set(NUMERIXX_SANITIZE "" CACHE STRING "Sanitizers for the tests, e.g. address;undefined")

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(GNUInstallDirs)
include("${PROJECT_SOURCE_DIR}/cmake/CPM.cmake")   # by path: a parent's own CPM.cmake on CMAKE_MODULE_PATH must not win
include(NumerixxDependencies)
include(NumerixxTargets)

if(NUMERIXX_INSTALL)
  include(NumerixxInstall)
endif()

if(NUMERIXX_BUILD_TESTS OR NUMERIXX_BUILD_INTEGRATION_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
if(NUMERIXX_BUILD_EXAMPLES)
  enable_testing()
  add_subdirectory(examples)
endif()
if(NUMERIXX_BUILD_BENCHMARKS)
  add_subdirectory(benchmarks)
endif()
if(NUMERIXX_BUILD_DOCS)
  add_subdirectory(docs)
endif()

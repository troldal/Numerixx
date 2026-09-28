# Numerixx test suite (DESIGN §9).
include(NumerixxBuildMode)

if(NUMERIXX_BUILD_INTEGRATION_TESTS)
  add_subdirectory(integration)
endif()
if(NOT NUMERIXX_BUILD_TESTS)
  return()
endif()

include(NumerixxWarnings)
include(NumerixxCompileFail)

# ---- doctest (fetched here, after NumerixxBuildMode) -------------------------------------------------------------
# doctest's own static main library is compiled as C++11, so the test main is built below instead, as C++23 and with
# the tests' exception model. doctest detects -fno-exceptions (and /EHs-c-) by itself, but then drops REQUIRE.
CPMAddPackage(
  NAME doctest
  VERSION 2.5.3
  URL https://github.com/doctest/doctest/archive/refs/tags/v2.5.3.tar.gz
  URL_HASH SHA256=174ebc4e769928959614789c5b4e9c3d0a0f81a62bb608756b127bfebfb21331
  EXCLUDE_FROM_ALL YES
  SYSTEM YES
  OPTIONS "DOCTEST_WITH_TESTS OFF" "DOCTEST_WITH_MAIN_IN_STATIC_LIB OFF" "DOCTEST_NO_INSTALL ON"
)
if(doctest_SOURCE_DIR)   # empty when CPM used an installed doctest, whose package config provides doctest.cmake
  list(APPEND CMAKE_MODULE_PATH "${doctest_SOURCE_DIR}/scripts/cmake")
endif()
include(doctest)

add_library(numerixx_doctest_main STATIC doctest_main.cpp)
target_link_libraries(numerixx_doctest_main PUBLIC doctest::doctest)
if(NUMERIXX_NO_EXCEPTIONS)
  # Keeps REQUIRE compiling without exceptions, where a failed REQUIRE reports but does not end the test case. Only
  # here: the macro also switches doctest to its no-exceptions mode (DESIGN §9.1).
  target_compile_definitions(numerixx_doctest_main PUBLIC DOCTEST_CONFIG_NO_EXCEPTIONS_BUT_WITH_ALL_ASSERTS)
endif()

# ---- Boost.Math oracles (standalone; tests only) -----------------------------------------------------------------
if(NUMERIXX_TEST_ORACLES AND NOT TARGET Boost::math)
  CPMAddPackage(
    NAME boost_math
    VERSION 1.92.0
    URL https://github.com/boostorg/math/archive/refs/tags/${NUMERIXX_BOOST_TAG}.tar.gz
    URL_HASH SHA256=aa84eec67402a51cb98172ece53b82f0500e95a80656de109074a71920cd2c48
    EXCLUDE_FROM_ALL YES
    SYSTEM YES
    OPTIONS "BOOST_MATH_STANDALONE ON"
  )
endif()

# numerixx_add_test(<name> SOURCES <files...> LINK <targets...>)
# One doctest executable per module, linking only that module (DESIGN §9.1). Every test case becomes a CTest test
# named <name>.<test case> and labelled <name> (ctest -L roots). Discovery runs the executable once, after it is built;
# doctest's ADD_LABELS would run it once more per test case, which is slow under Emscripten's node.
function(numerixx_add_test name)
  cmake_parse_arguments(ARG "" "" "SOURCES;LINK" ${ARGN})
  set(tgt numerixx_test_${name})
  add_executable(${tgt} ${ARG_SOURCES})
  target_link_libraries(${tgt} PRIVATE ${ARG_LINK} numerixx_doctest_main)
  numerixx_target_warnings(${tgt})
  doctest_discover_tests(${tgt}
    TEST_PREFIX "${name}."
    PROPERTIES LABELS ${name})
endfunction()

# ---- Unit tests, one executable per module -----------------------------------------------------------------------
numerixx_add_test(core SOURCES core/test_version.cpp core/test_refined.cpp core/test_criteria.cpp LINK numerixx::core)
target_compile_definitions(numerixx_test_core PRIVATE NUMERIXX_EXPECTED_VERSION="${PROJECT_VERSION}")

numerixx_add_test(roots
  SOURCES roots/test_solvers.cpp roots/test_soundness.cpp roots/test_combinators.cpp roots/test_steps.cpp
          roots/test_any_solver.cpp roots/test_determinism.cpp
  LINK numerixx::roots)
numerixx_add_test(deriv SOURCES deriv/test_deriv.cpp LINK numerixx::deriv)

# ---- Usage: the canonical calls (DESIGN §6.14) and cross-module properties ---------------------------------------
numerixx_add_test(usage
  SOURCES usage/canonical_calls.cpp usage/test_composition.cpp usage/test_regularity.cpp
  LINK numerixx::roots numerixx::deriv)

if(TARGET numerixx::linalg)
  numerixx_add_test(linalg SOURCES linalg/test_linalg.cpp LINK numerixx::linalg)
endif()
if(TARGET numerixx::pipes)
  numerixx_add_test(pipes SOURCES pipes/test_pipes.cpp LINK numerixx::pipes numerixx::roots numerixx::deriv)
endif()
if(TARGET numerixx::multiprecision)
  # Scalar tests: Boost.Multiprecision without the adapter, because the scalar traits need none (DESIGN D16).
  # Phases 1-8 add their multiprecision instantiations here.
  numerixx_add_test(multiprecision SOURCES multiprecision/test_multiprecision.cpp
                    LINK numerixx::roots numerixx::deriv Boost::multiprecision Boost::config)
  numerixx_add_test(multiprecision_adapter SOURCES multiprecision/test_adapter.cpp LINK numerixx::multiprecision)
endif()
if(NUMERIXX_TEST_ORACLES)
  numerixx_add_test(oracles SOURCES oracles/test_oracles.cpp LINK numerixx::core Boost::math)
endif()

# ---- Structural: header self-containment (DESIGN §5.3) ----------------------------------------------------------
# One translation unit per public header, which includes it twice (#pragma once). Each is compiled against its own
# module's target only, so a header that needs something its module does not provide fails here. Headers whose
# module is not configured are skipped.
file(GLOB_RECURSE _nxx_headers CONFIGURE_DEPENDS RELATIVE "${PROJECT_SOURCE_DIR}/include"
     "${PROJECT_SOURCE_DIR}/include/numerixx/*.hpp")
set(_nxx_header_targets "")
foreach(header IN LISTS _nxx_headers)
  if(header STREQUAL "numerixx/numerixx.hpp")
    set(module numerixx)                                        # the umbrella header: numerixx::numerixx
  elseif(header MATCHES "^numerixx/adapters/")
    set(module multiprecision)
  elseif(header MATCHES "^numerixx/([^/.]+)")
    set(module "${CMAKE_MATCH_1}")
    if(module STREQUAL "version" OR module STREQUAL "config")
      set(module core)
    endif()
  endif()
  if(NOT TARGET numerixx::${module})
    continue()
  endif()
  if(header MATCHES "_linalg\\.hpp$" AND NOT TARGET numerixx::linalg)
    continue()
  endif()
  string(MAKE_C_IDENTIFIER "${header}" id)
  set(source "${CMAKE_CURRENT_BINARY_DIR}/header_check/${id}.cpp")
  file(CONFIGURE OUTPUT "${source}" CONTENT "#include <${header}>\n#include <${header}>\n")
  set(tgt numerixx_header_check_${module})
  if(NOT TARGET ${tgt})
    add_library(${tgt} OBJECT)
    target_link_libraries(${tgt} PRIVATE numerixx::${module})
    numerixx_target_warnings(${tgt})
    list(APPEND _nxx_header_targets ${tgt})
  endif()
  target_sources(${tgt} PRIVATE "${source}")
endforeach()

# ---- Structural: a consumer with strict warnings (DESIGN §9.1) ---------------------------------------------------
# Compiled with warnings as errors whatever NUMERIXX_WARNINGS_AS_ERRORS says, without /Zc:__cplusplus on MSVC, and
# with a global `f`, which is what makes MSVC's C4459 fire in headers that name a parameter `f`.
add_library(numerixx_consumer_warnings OBJECT structural/consumer_warnings.cpp)
target_link_libraries(numerixx_consumer_warnings PRIVATE numerixx::numerixx)
numerixx_target_warnings(numerixx_consumer_warnings CONSUMER)

# ---- Structural: layering against the module DAG (DESIGN §5.2) ---------------------------------------------------
add_test(NAME structural.layering
         COMMAND ${CMAKE_COMMAND} -DROOT=${PROJECT_SOURCE_DIR}/include/numerixx
                 -P ${PROJECT_SOURCE_DIR}/cmake/CheckLayering.cmake)
set_tests_properties(structural.layering PROPERTIES LABELS structural)

# ---- Compile-fail tests (DESIGN §9.1) ----------------------------------------------------------------------------
include(compile_fail/compile_fail.cmake)

# ---- Structural: compile time of the umbrella header and of a linalg TU (DESIGN §3.6) ----------------------------
include(structural/compile_time.cmake)

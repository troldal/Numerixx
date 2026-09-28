# Consumer-build scenarios (DESIGN §4.1 principle 2, §9.4 "consumers"): Numerixx consumed by a CPM parent and by a
# FetchContent parent, each in both declaration orders; next to a parent's own "Boost" package; as a scalar-only
# subproject that must not download FXT or Eigen; and as an installed package through find_package.
#
# Each scenario configures, builds and runs a small parent project in its own build tree; they run only for host
# builds. They need network access: most of them are served by a warm CPM_SOURCE_CACHE, but
# fetchcontent_parent_deps_first fetches FXT and Eigen through plain FetchContent every time, and
# cpm_parent_deps_first downloads CPM 0.42.1 on its first run.
if(CMAKE_CROSSCOMPILING)
  message(STATUS "Numerixx: integration scenarios are skipped when cross-compiling")
  return()
endif()

if(DEFINED ENV{CPM_SOURCE_CACHE})
  set(_nxx_cache "$ENV{CPM_SOURCE_CACHE}")
elseif(CPM_SOURCE_CACHE)
  set(_nxx_cache "${CPM_SOURCE_CACHE}")
else()
  set(_nxx_cache "${CMAKE_CURRENT_BINARY_DIR}/cpm-cache")      # shared by all scenarios
endif()

get_property(_nxx_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(_nxx_multi_config)
  set(_nxx_build_type Release)
else()
  set(_nxx_build_type ${CMAKE_BUILD_TYPE})
endif()

# numerixx_add_scenario(<name> [INSTALLED])
function(numerixx_add_scenario name)
  cmake_parse_arguments(ARG "INSTALLED" "" "" ${ARGN})
  set(mode build)
  if(ARG_INSTALLED)
    set(mode installed)
  endif()
  add_test(NAME integration.${name}
           COMMAND ${CMAKE_COMMAND}
                   -DSCENARIO=${name}
                   -DMODE=${mode}
                   -DSCENARIO_DIR=${CMAKE_CURRENT_SOURCE_DIR}/${name}
                   -DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}/${name}
                   -DNUMERIXX_SOURCE_DIR=${PROJECT_SOURCE_DIR}
                   -DGENERATOR=${CMAKE_GENERATOR}
                   -DMAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
                   -DCXX_COMPILER=${CMAKE_CXX_COMPILER}
                   -DBUILD_TYPE=${_nxx_build_type}
                   -DCPM_SOURCE_CACHE=${_nxx_cache}
                   -DFXT_REF=${NUMERIXX_FXT_REF}
                   -DFXT_SHA256=${NUMERIXX_FXT_SHA256}
                   -DEIGEN_VERSION=${NUMERIXX_EIGEN_VERSION}
                   -DEIGEN_SHA256=${NUMERIXX_EIGEN_SHA256}
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/run_scenario.cmake)
  set_tests_properties(integration.${name} PROPERTIES LABELS integration TIMEOUT 1800)
endfunction()

numerixx_add_scenario(cpm_parent_deps_first)
numerixx_add_scenario(cpm_parent_numerixx_first)
numerixx_add_scenario(fetchcontent_parent_deps_first)
numerixx_add_scenario(fetchcontent_parent_numerixx_first)
numerixx_add_scenario(boost_named_parent)
numerixx_add_scenario(boost_named_parent_numerixx_first)
numerixx_add_scenario(scalar_only)
numerixx_add_scenario(installed INSTALLED)

# Scenarios share one CPM cache; run them one at a time so that no two populate it at once.
set_tests_properties(
  integration.cpm_parent_deps_first integration.cpm_parent_numerixx_first
  integration.fetchcontent_parent_deps_first integration.fetchcontent_parent_numerixx_first
  integration.boost_named_parent integration.boost_named_parent_numerixx_first integration.scalar_only
  integration.installed
  PROPERTIES RESOURCE_LOCK numerixx_cpm_cache)

# A parent that uses only the scalar modules (NUMERIXX_WITH_FXT=OFF, NUMERIXX_WITH_LINALG=OFF). Numerixx must download
# nothing at all, and the FXT and Eigen targets must not exist.
cmake_minimum_required(VERSION 3.25)
project(scalar_only LANGUAGES CXX)

include("${NUMERIXX_SOURCE_DIR}/cmake/get_cpm.cmake")
CPMAddPackage(NAME Numerixx SOURCE_DIR "${NUMERIXX_SOURCE_DIR}"
              OPTIONS "NUMERIXX_WITH_FXT OFF" "NUMERIXX_WITH_LINALG OFF")

# CPM records every package it adds in the cache variable CPM_PACKAGES, which (unlike <name>_SOURCE_DIR, set in
# Numerixx's directory scope) is visible here. Numerixx itself must be the only entry.
set(extra_packages ${CPM_PACKAGES})
list(REMOVE_ITEM extra_packages Numerixx)
if(extra_packages)
  message(FATAL_ERROR "a scalar-only Numerixx must not download anything (CPM packages: ${CPM_PACKAGES})")
endif()
foreach(target IN ITEMS fxt::fxt Eigen3::Eigen numerixx::pipes numerixx::linalg numerixx::multiroots)
  if(TARGET ${target})
    message(FATAL_ERROR "a scalar-only Numerixx must not define ${target}")
  endif()
endforeach()

add_executable(consumer "${NUMERIXX_SOURCE_DIR}/tests/integration/common/scalar.cpp")
target_link_libraries(consumer PRIVATE numerixx::numerixx)

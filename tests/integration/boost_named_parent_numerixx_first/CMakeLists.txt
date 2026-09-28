# A parent that adds Numerixx (with the multiprecision adapter) first and declares its own CPM package named "Boost"
# afterwards. This is the order in which a Numerixx package called "Boost" would take the parent's place: CPM would
# deduplicate the parent's declaration away. Numerixx must not claim the name, so the parent still gets its own.
cmake_minimum_required(VERSION 3.25)
project(boost_named_parent_numerixx_first LANGUAGES CXX)

include("${NUMERIXX_SOURCE_DIR}/cmake/get_cpm.cmake")

CPMAddPackage(NAME Numerixx SOURCE_DIR "${NUMERIXX_SOURCE_DIR}" OPTIONS "NUMERIXX_WITH_MULTIPRECISION ON")

# The parent's stand-in for the full Boost distribution (shared with boost_named_parent).
set(parent_boost_dir "${NUMERIXX_SOURCE_DIR}/tests/integration/boost_named_parent/fake_boost")
CPMAddPackage(NAME Boost SOURCE_DIR "${parent_boost_dir}" DOWNLOAD_ONLY YES)
if(NOT Boost_ADDED)
  message(FATAL_ERROR "the parent's Boost package was deduplicated against a package Numerixx declared")
endif()

include("${NUMERIXX_SOURCE_DIR}/tests/integration/common/check_boost.cmake")
nxx_check_parent_boost("${parent_boost_dir}")

add_executable(consumer "${NUMERIXX_SOURCE_DIR}/tests/integration/common/multiprecision.cpp")
target_link_libraries(consumer PRIVATE numerixx::numerixx numerixx::multiprecision)

# A parent that declares its own CPM package named "Boost" (as projects that fetch the full Boost distribution do),
# before it adds Numerixx. Numerixx's standalone Boost packages use other names, so the parent's Boost stays the
# only package of that name, and it is untouched. boost_named_parent_numerixx_first covers the other order.
cmake_minimum_required(VERSION 3.25)
project(boost_named_parent LANGUAGES CXX)

include("${NUMERIXX_SOURCE_DIR}/cmake/get_cpm.cmake")

# A stand-in for the full Boost distribution: a local directory, so that the scenario downloads nothing extra.
set(parent_boost_dir "${CMAKE_CURRENT_SOURCE_DIR}/fake_boost")
CPMAddPackage(NAME Boost SOURCE_DIR "${parent_boost_dir}" DOWNLOAD_ONLY YES)

CPMAddPackage(NAME Numerixx SOURCE_DIR "${NUMERIXX_SOURCE_DIR}" OPTIONS "NUMERIXX_WITH_MULTIPRECISION ON")

include("${NUMERIXX_SOURCE_DIR}/tests/integration/common/check_boost.cmake")
nxx_check_parent_boost("${parent_boost_dir}")

add_executable(consumer "${NUMERIXX_SOURCE_DIR}/tests/integration/common/multiprecision.cpp")
target_link_libraries(consumer PRIVATE numerixx::numerixx numerixx::multiprecision)

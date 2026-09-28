# A CPM parent that adds Numerixx first and declares FXT and Eigen afterwards under the same package names. CPM must
# deduplicate both, so the parent uses the copies Numerixx fetched.
cmake_minimum_required(VERSION 3.25)
project(cpm_parent_numerixx_first LANGUAGES CXX)

include("${NUMERIXX_SOURCE_DIR}/cmake/get_cpm.cmake")

CPMAddPackage(NAME Numerixx SOURCE_DIR "${NUMERIXX_SOURCE_DIR}")

CPMAddPackage(NAME FXT URL https://github.com/troldal/FXT/archive/${FXT_REF}.tar.gz URL_HASH SHA256=${FXT_SHA256}
              DOWNLOAD_ONLY YES)
CPMAddPackage(NAME Eigen3 VERSION ${EIGEN_VERSION}
              URL https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_VERSION}/eigen-${EIGEN_VERSION}.tar.gz
              URL_HASH SHA256=${EIGEN_SHA256} DOWNLOAD_ONLY YES)

if(FXT_ADDED OR Eigen3_ADDED)
  message(FATAL_ERROR "CPM fetched FXT or Eigen a second time instead of reusing Numerixx's copies")
endif()
# The parent's own declarations were deduplicated, so it relies on the usual target names that Numerixx defines.
if(NOT TARGET fxt::fxt OR NOT TARGET Eigen3::Eigen OR NOT TARGET numerixx::linalg)
  message(FATAL_ERROR "Numerixx did not provide fxt::fxt, Eigen3::Eigen and numerixx::linalg")
endif()

add_executable(consumer "${NUMERIXX_SOURCE_DIR}/tests/integration/common/full.cpp")
target_link_libraries(consumer PRIVATE numerixx::numerixx numerixx::linalg fxt::fxt Eigen3::Eigen)

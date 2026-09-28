# A CPM parent on an older CPM (0.42.1) that declares FXT and Eigen itself, with its own targets, before it adds
# Numerixx. Numerixx must reuse the parent's CPM, fxt::fxt and Eigen3::Eigen instead of creating its own.
cmake_minimum_required(VERSION 3.25)
project(cpm_parent_deps_first LANGUAGES CXX)

set(CPM_DOWNLOAD_VERSION 0.42.1)
set(CPM_DOWNLOAD_LOCATION "${CPM_SOURCE_CACHE}/cpm/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
file(DOWNLOAD https://github.com/cpm-cmake/CPM.cmake/releases/download/v${CPM_DOWNLOAD_VERSION}/CPM.cmake
     "${CPM_DOWNLOAD_LOCATION}"
     EXPECTED_HASH SHA256=f3a6dcc6a04ce9e7f51a127307fa4f699fb2bade357a8eb4c5b45df76e1dc6a5)
include("${CPM_DOWNLOAD_LOCATION}")

CPMAddPackage(NAME FXT URL https://github.com/troldal/FXT/archive/${FXT_REF}.tar.gz URL_HASH SHA256=${FXT_SHA256}
              DOWNLOAD_ONLY YES)
add_library(parent_fxt INTERFACE)
target_include_directories(parent_fxt SYSTEM INTERFACE "${FXT_SOURCE_DIR}/include")
add_library(fxt::fxt ALIAS parent_fxt)

CPMAddPackage(NAME Eigen3 VERSION ${EIGEN_VERSION}
              URL https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_VERSION}/eigen-${EIGEN_VERSION}.tar.gz
              URL_HASH SHA256=${EIGEN_SHA256} DOWNLOAD_ONLY YES)
add_library(parent_eigen INTERFACE)
target_include_directories(parent_eigen SYSTEM INTERFACE "${Eigen3_SOURCE_DIR}")
add_library(Eigen3::Eigen ALIAS parent_eigen)

CPMAddPackage(NAME Numerixx SOURCE_DIR "${NUMERIXX_SOURCE_DIR}")

if(TARGET numerixx_fxt OR TARGET numerixx_eigen)
  message(FATAL_ERROR "Numerixx created its own FXT or Eigen target instead of reusing the parent's")
endif()
if(TARGET numerixx_test_core)
  message(FATAL_ERROR "Numerixx built its tests although it is a subproject")
endif()

add_executable(consumer "${NUMERIXX_SOURCE_DIR}/tests/integration/common/full.cpp")
target_link_libraries(consumer PRIVATE numerixx::numerixx numerixx::linalg)

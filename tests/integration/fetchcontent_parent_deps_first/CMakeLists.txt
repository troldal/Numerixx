# A FetchContent parent (no CPM of its own) that fetches FXT and Eigen and defines fxt::fxt and Eigen3::Eigen before
# it adds Numerixx. Numerixx must reuse both targets.
cmake_minimum_required(VERSION 3.25)
project(fetchcontent_parent_deps_first LANGUAGES CXX)

include(FetchContent)
# SOURCE_SUBDIR points at a directory that does not exist, so FetchContent populates without running the projects'
# own CMake; the parent defines the targets itself.
FetchContent_Declare(fxt URL https://github.com/troldal/FXT/archive/${FXT_REF}.tar.gz URL_HASH SHA256=${FXT_SHA256}
                     SOURCE_SUBDIR no-cmake)
FetchContent_Declare(eigen URL https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_VERSION}/eigen-${EIGEN_VERSION}.tar.gz
                     URL_HASH SHA256=${EIGEN_SHA256} SOURCE_SUBDIR no-cmake)
FetchContent_MakeAvailable(fxt eigen)

add_library(parent_fxt INTERFACE)
target_include_directories(parent_fxt SYSTEM INTERFACE "${fxt_SOURCE_DIR}/include")
add_library(fxt::fxt ALIAS parent_fxt)
add_library(parent_eigen INTERFACE)
target_include_directories(parent_eigen SYSTEM INTERFACE "${eigen_SOURCE_DIR}")
add_library(Eigen3::Eigen ALIAS parent_eigen)

FetchContent_Declare(numerixx SOURCE_DIR "${NUMERIXX_SOURCE_DIR}")
FetchContent_MakeAvailable(numerixx)

if(TARGET numerixx_fxt OR TARGET numerixx_eigen)
  message(FATAL_ERROR "Numerixx created its own FXT or Eigen target instead of reusing the parent's")
endif()

add_executable(consumer "${NUMERIXX_SOURCE_DIR}/tests/integration/common/full.cpp")
target_link_libraries(consumer PRIVATE numerixx::numerixx numerixx::linalg)

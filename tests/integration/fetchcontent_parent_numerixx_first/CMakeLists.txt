# A FetchContent parent (no CPM of its own) that adds Numerixx first and then uses the fxt::fxt and Eigen3::Eigen
# targets that Numerixx provides, instead of fetching FXT or Eigen a second time.
cmake_minimum_required(VERSION 3.25)
project(fetchcontent_parent_numerixx_first LANGUAGES CXX)

include(FetchContent)
FetchContent_Declare(numerixx SOURCE_DIR "${NUMERIXX_SOURCE_DIR}")
FetchContent_MakeAvailable(numerixx)

foreach(target IN ITEMS fxt::fxt Eigen3::Eigen numerixx::pipes numerixx::linalg)
  if(NOT TARGET ${target})
    message(FATAL_ERROR "Numerixx did not provide ${target}")
  endif()
endforeach()

add_executable(consumer "${NUMERIXX_SOURCE_DIR}/tests/integration/common/full.cpp")
target_link_libraries(consumer PRIVATE numerixx::numerixx numerixx::linalg fxt::fxt Eigen3::Eigen)

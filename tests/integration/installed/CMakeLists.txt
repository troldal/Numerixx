# A consumer of an installed Numerixx (find_package in config mode). run_scenario.cmake installs Numerixx into a
# staging prefix first and passes it as CMAKE_PREFIX_PATH.
cmake_minimum_required(VERSION 3.25)
project(installed LANGUAGES CXX)

find_package(numerixx 2.0 CONFIG REQUIRED)
foreach(target IN ITEMS numerixx::numerixx numerixx::pipes numerixx::linalg numerixx::multiroots fxt::fxt)
  if(NOT TARGET ${target})
    message(FATAL_ERROR "the installed package does not provide ${target}")
  endif()
endforeach()

# The FXT and Eigen headers Numerixx fetched are installed below a private directory, never into the shared
# include/eigen3 or include/fxt, where they could overwrite another installation; all their licences come along.
set(stage "${CMAKE_PREFIX_PATH}")
foreach(required IN ITEMS include/numerixx/numerixx.hpp include/numerixx-deps/fxt/fxt/monads/Expected.hpp
                          include/numerixx-deps/eigen3/Eigen/Core share/doc/Numerixx/LICENSE
                          share/doc/Numerixx/fxt/LICENSE share/doc/Numerixx/eigen/COPYING.MPL2
                          share/doc/Numerixx/eigen/COPYING.BSD)
  if(NOT EXISTS "${stage}/${required}")
    message(FATAL_ERROR "the installed package lacks ${required}")
  endif()
endforeach()
foreach(forbidden IN ITEMS include/eigen3 include/fxt include/fxt.hpp)
  if(EXISTS "${stage}/${forbidden}")
    message(FATAL_ERROR "the installed package writes into the shared location ${forbidden}")
  endif()
endforeach()

add_executable(consumer "${NUMERIXX_SOURCE_DIR}/tests/integration/common/full.cpp")
target_link_libraries(consumer PRIVATE numerixx::numerixx numerixx::linalg)

# The same consumer with the imported include directories as plain -I paths: they must be absolute and exist.
add_executable(consumer_plain_includes "${NUMERIXX_SOURCE_DIR}/tests/integration/common/full.cpp")
target_link_libraries(consumer_plain_includes PRIVATE numerixx::numerixx numerixx::linalg)
set_target_properties(consumer_plain_includes PROPERTIES NO_SYSTEM_FROM_IMPORTED ON)

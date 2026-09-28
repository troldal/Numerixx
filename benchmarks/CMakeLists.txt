# Benchmarks (google/benchmark). They report function evaluations next to time; CI builds them on the gcc leg but
# does not time them (DESIGN §9.1).
if(EMSCRIPTEN)
  message(STATUS "Numerixx: benchmarks are not built under Emscripten")
  return()
endif()

include(NumerixxBuildMode)
include(NumerixxWarnings)

# google/benchmark derives its version from `git describe`. Its archive has no .git, so git would walk up into an
# enclosing checkout (Numerixx's own, when the CPM cache lives inside it) and report that repository's tag. With
# GIT_DIR pointing at nothing, git fails and benchmark falls back to its project version, 1.9.5.
if(DEFINED ENV{GIT_DIR})
  set(_nxx_saved_git_dir "$ENV{GIT_DIR}")
endif()
set(ENV{GIT_DIR} "${CMAKE_CURRENT_BINARY_DIR}/no-git")
CPMAddPackage(
  NAME benchmark
  VERSION 1.9.5
  URL https://github.com/google/benchmark/archive/refs/tags/v1.9.5.tar.gz
  URL_HASH SHA256=9631341c82bac4a288bef951f8b26b41f69021794184ece969f8473977eaa340
  EXCLUDE_FROM_ALL YES
  SYSTEM YES
  OPTIONS "BENCHMARK_ENABLE_TESTING OFF" "BENCHMARK_ENABLE_GTEST_TESTS OFF" "BENCHMARK_ENABLE_INSTALL OFF"
          "BENCHMARK_ENABLE_WERROR OFF" "BENCHMARK_INSTALL_DOCS OFF"
)
if(DEFINED _nxx_saved_git_dir)
  set(ENV{GIT_DIR} "${_nxx_saved_git_dir}")
else()
  unset(ENV{GIT_DIR})
endif()

add_executable(numerixx_benchmarks bench_smoke.cpp)
target_link_libraries(numerixx_benchmarks PRIVATE numerixx::numerixx benchmark::benchmark_main)
numerixx_target_warnings(numerixx_benchmarks)

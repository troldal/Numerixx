# Warning flags for Numerixx-owned targets only (tests, examples, benchmarks). Consumers never inherit them
# (DESIGN §4.1, principle 4).
#
#   numerixx_target_warnings(<target> [CONSUMER])
#
# CONSUMER is for the structural consumer test: warnings are errors whatever NUMERIXX_WARNINGS_AS_ERRORS says, and
# MSVC's /Zc:__cplusplus is left out, so the headers are compiled the way a default MSVC project sees them
# (__cplusplus == 199711L).
function(numerixx_target_warnings tgt)
  cmake_parse_arguments(ARG "CONSUMER" "" "" ${ARGN})
  set(_werror OFF)
  if(NUMERIXX_WARNINGS_AS_ERRORS OR ARG_CONSUMER)
    set(_werror ON)
  endif()

  if(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")      # cl and clang-cl
    target_compile_options(${tgt} PRIVATE /W4 /permissive- /utf-8)
    if(NOT ARG_CONSUMER)
      target_compile_options(${tgt} PRIVATE /Zc:__cplusplus)     # cl reports 199711L otherwise
    endif()
    if(_werror)
      target_compile_options(${tgt} PRIVATE /WX)
    endif()
  else()
    target_compile_options(${tgt} PRIVATE
      -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -Wnon-virtual-dtor -Wold-style-cast
      -Wcast-align -Wunused -Woverloaded-virtual -Wnull-dereference -Wdouble-promotion -Wformat=2
      $<$<CXX_COMPILER_ID:GNU>:-Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wuseless-cast>)
    if(_werror)
      target_compile_options(${tgt} PRIVATE -Werror)
    endif()
  endif()
endfunction()

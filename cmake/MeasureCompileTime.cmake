# cmake -DBUILD_DIR=<dir> -DTARGET=<target> -DCONFIG=<config> -DOBJECTS=<objects> -DNAME=<name>
#       [-DCOMPILER_ID=<id> -DCOMPILER_VERSION=<v> -DCOMPILER_FRONTEND=<GNU|MSVC>]
#       [-DGUARD_SECONDS=<s> -DGUARD_COMPILER=<id>] [-DRUNS=<n>] [-DREPORT=<file>]
#       -P MeasureCompileTime.cmake
#
# Compile-time measurement (DESIGN §3.6): removes TARGET's object files (OBJECTS, from $<TARGET_OBJECTS:...>), rebuilds
# TARGET with `cmake --build`, RUNS times (default 3), and takes the best wall time. That time includes the build
# tool's own start-up, a few tens of milliseconds. It prints "compile time: <name> <seconds> s (best of <n>)" and
# appends a line to REPORT. When COMPILER_ID equals GUARD_COMPILER, a best time above GUARD_SECONDS fails the test
# (the guard); otherwise the time is only recorded.
cmake_minimum_required(VERSION 3.25)   # script mode starts without policies; %f in string(TIMESTAMP) needs 3.23

foreach(var IN ITEMS BUILD_DIR TARGET CONFIG OBJECTS NAME)
  if(NOT DEFINED ${var} OR "${${var}}" STREQUAL "")
    message(FATAL_ERROR "MeasureCompileTime: ${var} is not set")
  endif()
endforeach()
if(NOT DEFINED RUNS)
  set(RUNS 3)
endif()

# Microseconds since the epoch: %s (seconds) followed by %f (microseconds, six digits). 64-bit math(EXPR) holds it.
function(_nxx_now out)
  string(TIMESTAMP now "%s%f" UTC)
  set(${out} "${now}" PARENT_SCOPE)
endfunction()

# "1.25" -> 1250000 microseconds (at most six decimals are used).
function(_nxx_seconds_to_us seconds out)
  if(NOT seconds MATCHES "^([0-9]+)(\\.([0-9]*))?$")
    message(FATAL_ERROR "MeasureCompileTime: '${seconds}' is not a number of seconds")
  endif()
  set(whole "${CMAKE_MATCH_1}")
  string(SUBSTRING "${CMAKE_MATCH_3}000000" 0 6 frac)
  string(REGEX REPLACE "^0+([0-9])" "\\1" frac "${frac}")    # no leading zeros in math(EXPR)
  string(REGEX REPLACE "^0+([0-9])" "\\1" whole "${whole}")
  math(EXPR us "${whole} * 1000000 + ${frac}")
  set(${out} "${us}" PARENT_SCOPE)
endfunction()

# 1234567 microseconds -> "1.235" (seconds, rounded to milliseconds).
function(_nxx_us_to_seconds us out)
  math(EXPR ms "(${us} + 500) / 1000")
  math(EXPR whole "${ms} / 1000")
  math(EXPR frac "${ms} % 1000")
  if(frac LESS 10)
    set(frac "00${frac}")
  elseif(frac LESS 100)
    set(frac "0${frac}")
  endif()
  set(${out} "${whole}.${frac}" PARENT_SCOPE)
endfunction()

set(best "")
set(all "")
foreach(run RANGE 1 ${RUNS})
  file(REMOVE ${OBJECTS})
  foreach(object IN LISTS OBJECTS)
    if(EXISTS "${object}")
      message(FATAL_ERROR "MeasureCompileTime: could not remove '${object}'")
    endif()
  endforeach()

  _nxx_now(start)
  execute_process(
    COMMAND "${CMAKE_COMMAND}" --build "${BUILD_DIR}" --target "${TARGET}" --config "${CONFIG}"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE output
    ERROR_VARIABLE  output)
  _nxx_now(stop)

  if(NOT result EQUAL 0)
    message(FATAL_ERROR "MeasureCompileTime: building '${TARGET}' failed:\n${output}")
  endif()
  foreach(object IN LISTS OBJECTS)
    if(NOT EXISTS "${object}")
      message(FATAL_ERROR "MeasureCompileTime: building '${TARGET}' did not produce '${object}'")
    endif()
  endforeach()

  math(EXPR elapsed "${stop} - ${start}")
  _nxx_us_to_seconds(${elapsed} seconds)
  list(APPEND all "${seconds}")
  if(best STREQUAL "" OR elapsed LESS best)
    set(best ${elapsed})
  endif()
endforeach()

_nxx_us_to_seconds(${best} best_seconds)
list(JOIN all " " all_text)
set(compiler "${COMPILER_ID}")
if(COMPILER_VERSION)
  string(APPEND compiler " ${COMPILER_VERSION}")
endif()
if(COMPILER_ID STREQUAL "Clang" AND COMPILER_FRONTEND STREQUAL "MSVC")
  string(APPEND compiler " (clang-cl)")
endif()
message(STATUS "compile time: ${NAME} ${best_seconds} s (best of ${RUNS})")
message(STATUS "  runs: ${all_text} s; ${compiler}, config ${CONFIG}")

if(REPORT)
  string(TIMESTAMP when "%Y-%m-%d %H:%M:%S")
  file(APPEND "${REPORT}" "${when} ${NAME} ${compiler} ${CONFIG}: ${best_seconds} s (best of ${RUNS}: ${all_text})\n")
endif()

if(DEFINED GUARD_SECONDS AND NOT GUARD_SECONDS STREQUAL "" AND COMPILER_ID STREQUAL GUARD_COMPILER)
  _nxx_seconds_to_us(${GUARD_SECONDS} guard)
  if(best GREATER guard)
    message(FATAL_ERROR "compile time: ${NAME} took ${best_seconds} s on ${COMPILER_ID}, above the ${GUARD_SECONDS} s guard "
                        "(DESIGN §3.6)")
  endif()
  message(STATUS "  within the ${GUARD_SECONDS} s guard for ${GUARD_COMPILER}")
endif()

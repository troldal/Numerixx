# cmake -DBUILD_DIR=<dir> -DTARGET=<target> -DCONFIG=<config> [-DEXPECT=<regex> -DCHECK_REASON=ON]
#       -P RunCompileFail.cmake
#
# Builds TARGET, which must fail to compile: the build must fail with a compiler diagnostic ("error:"), not merely
# fail (a missing target or a build-tool error does not count). With CHECK_REASON, EXPECT must match the first
# error: the first diagnostic line and the lines after it, up to the next error or the end of the compiler output.
cmake_minimum_required(VERSION 3.25)   # script mode starts without policies

execute_process(
  COMMAND "${CMAKE_COMMAND}" --build "${BUILD_DIR}" --target "${TARGET}" --config "${CONFIG}"
  RESULT_VARIABLE result
  OUTPUT_VARIABLE output
  ERROR_VARIABLE  output)

if(result EQUAL 0)
  message(FATAL_ERROR "compile-fail: '${TARGET}' compiled, but it must not")
endif()

# Strip ANSI colour sequences (coloured diagnostics end the sequence with a letter right before "error:").
string(ASCII 27 esc)
string(REGEX REPLACE "${esc}\\[[0-9;]*[A-Za-z]" "" output "${output}")

# Split into lines. Brackets and semicolons are masked first: in a CMake list, a ';' after an unbalanced '[' is not
# a separator, which would merge the rest of the output into one element.
string(REPLACE "[" "<nxx-lb>" masked "${output}")
string(REPLACE "]" "<nxx-rb>" masked "${masked}")
string(REPLACE ";" "<nxx-sc>" masked "${masked}")
string(REGEX REPLACE "\r?\n" ";" lines "${masked}")

# A diagnostic looks like "file:line:col: error: ..." (GCC, Clang), "file(line,col): error: ..." (clang-cl) or
# "file(line): error C1234: ..." (cl). The compile command echoed by the build tool contains "-Werror" but never
# "error:", and lines of the build tool itself ("ninja: ...", "FAILED: ...") are not diagnostics.
set(first_error "")
set(in_first OFF)
foreach(line IN LISTS lines)
  if(line MATCHES "^(ninja|FAILED|make|gmake|MSBuild)")
    if(in_first)
      break()
    endif()
    continue()
  endif()
  if(line MATCHES "^error( [A-Z]+[0-9]+)?:" OR line MATCHES "[^A-Za-z]error( [A-Z]+[0-9]+)?:")
    if(in_first)
      break()
    endif()
    set(in_first ON)
  endif()
  if(in_first)
    string(APPEND first_error "${line}\n")
  endif()
endforeach()

string(REPLACE "<nxx-lb>" "[" first_error "${first_error}")
string(REPLACE "<nxx-rb>" "]" first_error "${first_error}")
string(REPLACE "<nxx-sc>" ";" first_error "${first_error}")

if(first_error STREQUAL "")
  message(FATAL_ERROR "compile-fail: building '${TARGET}' failed without a compiler error:\n${output}")
endif()
if(CHECK_REASON AND NOT first_error MATCHES "${EXPECT}")
  message(FATAL_ERROR "compile-fail: the first error of '${TARGET}' does not match '${EXPECT}':\n${first_error}")
endif()
message(STATUS "compile-fail: '${TARGET}' failed to compile, as required:\n${first_error}")

# cmake -DBUILD_DIR=<dir> -DTARGET=<target> -DCONFIG=<config> [-DEXPECT=<regex> -DCHECK_REASON=ON [-DREASON_IN_MESSAGE=ON]]
#       [-DREPORT=<file>] -P RunCompileFail.cmake
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
# GCC reports a failed consteval call as "call to consteval function ... is not a constant expression" and explains it
# in a nested error that follows its "in 'constexpr' expansion of ..." context lines; that nested error belongs to the
# first one.
set(first_error "")
set(in_first OFF)
set(nested OFF)
foreach(line IN LISTS lines)
  if(line MATCHES "^(ninja|FAILED|make|gmake|MSBuild)")
    if(in_first)
      break()
    endif()
    continue()
  endif()
  if(line MATCHES "^error( [A-Z]+[0-9]+)?:" OR line MATCHES "[^A-Za-z]error( [A-Z]+[0-9]+)?:")
    if(in_first AND NOT nested)
      break()
    endif()
    set(in_first ON)
    set(nested OFF)
  elseif(in_first AND line MATCHES "in 'constexpr' expansion of")
    set(nested ON)
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

# Line counts, recorded (not gated): the diagnostic lines of the whole build output and of the first error.
set(diagnostic_lines 0)
foreach(line IN LISTS lines)
  if(NOT line MATCHES "^(ninja|FAILED|make|gmake|MSBuild|<nxx-lb>[0-9]+/[0-9]+<nxx-rb>)" AND NOT line STREQUAL "")
    math(EXPR diagnostic_lines "${diagnostic_lines} + 1")
  endif()
endforeach()
string(REGEX MATCHALL "\n" first_error_newlines "${first_error}")
list(LENGTH first_error_newlines first_error_lines)
if(REPORT)
  file(APPEND "${REPORT}" "${TARGET}: ${diagnostic_lines} diagnostic lines, first error ${first_error_lines} lines\n")
endif()

# With REASON_IN_MESSAGE (deletion reasons, which GCC 15+ and Clang 19+ print in the message), the reason must be in
# the compiler's own text, not in a source line it quotes: a note such as "declared here" echoes the NXX_DELETE("...")
# line, which would match even if the compiler printed no reason. GCC and Clang quote source as "  123 | code",
# followed by a caret line "      |   ^~~"; both are left out of the match. A rejected consteval literal is different:
# GCC names the failing Tag::reject() and quotes its source line, which holds the reason, so there the quote counts.
string(REPLACE "[" "<nxx-lb>" first_masked "${first_error}")
string(REPLACE "]" "<nxx-rb>" first_masked "${first_masked}")
string(REPLACE ";" "<nxx-sc>" first_masked "${first_masked}")
string(REGEX REPLACE "\n" ";" first_lines "${first_masked}")
set(first_error_text "")
foreach(line IN LISTS first_lines)
  if(NOT line MATCHES "^ *[0-9]* \\|")
    string(APPEND first_error_text "${line}\n")
  endif()
endforeach()
string(REPLACE "<nxx-lb>" "[" first_error_text "${first_error_text}")
string(REPLACE "<nxx-rb>" "]" first_error_text "${first_error_text}")
string(REPLACE "<nxx-sc>" ";" first_error_text "${first_error_text}")

if(NOT REASON_IN_MESSAGE)
  set(first_error_text "${first_error}")
endif()
if(CHECK_REASON AND NOT first_error_text MATCHES "${EXPECT}")
  message(FATAL_ERROR "compile-fail: the first error of '${TARGET}' does not match '${EXPECT}':\n${first_error}")
endif()
message(STATUS "compile-fail: '${TARGET}' failed to compile, as required (${diagnostic_lines} diagnostic lines, first error "
               "${first_error_lines} lines):\n${first_error}")

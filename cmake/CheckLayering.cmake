# cmake -DROOT=<path to include/numerixx> -P CheckLayering.cmake
#
# Enforces the module DAG of DESIGN §5.2 on the public headers:
#   - a module may include itself, core (incl. config.hpp and version.hpp) and the modules listed in ALLOW_<module>;
#   - Eigen may be included only by linalg, multiroots and adapters/multiprecision_linalg.hpp, and Eigen's unsupported
#     modules by no public header (they are not installed);
#   - FXT only by pipes; Boost only by the adapters;
#   - the umbrella numerixx.hpp includes the scalar modules and pipes, never linalg, multiroots, the adapters or
#     core/any_solver.hpp;
#   - public headers include with angle brackets only (<numerixx/...>, <Eigen/...>), so that every include is
#     classified here; a quoted include ("linalg.hpp") is itself a violation.
# A header's module is its top-level directory (numerixx/roots/brent.hpp -> roots) or, for a top-level header, its
# name (numerixx/roots.hpp -> roots).
cmake_minimum_required(VERSION 3.25)   # script mode starts without policies (IN_LIST needs CMP0057)

set(MODULES core deriv roots optimize poly integrate interpolate linalg multiroots pipes adapters)
set(ALLOW_core        "")
set(ALLOW_deriv       "")
set(ALLOW_roots       "")
set(ALLOW_optimize    "")
set(ALLOW_poly        "")
set(ALLOW_integrate   "")
set(ALLOW_interpolate "")
set(ALLOW_linalg      "")
set(ALLOW_multiroots  "linalg;deriv")
set(ALLOW_pipes       "")
set(ALLOW_adapters    "${MODULES}")
set(UMBRELLA_FORBIDDEN linalg multiroots adapters)

if(NOT IS_DIRECTORY "${ROOT}")
  message(FATAL_ERROR "CheckLayering: ROOT='${ROOT}' is not a directory")
endif()

# Maps a path below include/numerixx ("roots/brent.hpp", "roots.hpp", "config.hpp") to its module.
function(_nxx_module_of path out)
  if(path MATCHES "^([^/]+)/")
    set(mod "${CMAKE_MATCH_1}")
  else()
    string(REGEX REPLACE "\\.hpp$" "" mod "${path}")
  endif()
  if(mod STREQUAL "config" OR mod STREQUAL "version")
    set(mod core)
  endif()
  set(${out} "${mod}" PARENT_SCOPE)
endfunction()

file(GLOB_RECURSE files RELATIVE "${ROOT}" "${ROOT}/*.hpp")
list(SORT files)
set(violations 0)
# A function, not a macro: a macro would substitute the message text and evaluate it again, so a backslash in an
# offending include line would be read as an escape sequence.
function(_nxx_violation msg)
  string(REPLACE "<nxx-lb>" "[" msg "${msg}")
  string(REPLACE "<nxx-rb>" "]" msg "${msg}")
  string(REPLACE "<nxx-sc>" ";" msg "${msg}")
  message(SEND_ERROR "layering: ${msg}")
  math(EXPR count "${violations} + 1")
  set(violations ${count} PARENT_SCOPE)
endfunction()

foreach(file IN LISTS files)
  # Read line by line. Brackets and semicolons are masked first: in a CMake list, a ';' after an unbalanced '['
  # (say, in a trailing comment) is not a separator, which would merge the following lines into one element.
  file(READ "${ROOT}/${file}" content)
  string(REPLACE "[" "<nxx-lb>" content "${content}")
  string(REPLACE "]" "<nxx-rb>" content "${content}")
  string(REPLACE ";" "<nxx-sc>" content "${content}")
  string(REGEX REPLACE "\r?\n" ";" includes "${content}")
  list(FILTER includes INCLUDE REGEX "^[ \t]*#[ \t]*include")

  # Only canonical spellings are classified below, so reject the ones a case-insensitive file system or a lenient
  # compiler would still accept (<numerixx\linalg.hpp>, <eigen/Core>) instead of letting them pass unclassified.
  foreach(line IN LISTS includes)
    string(TOLOWER "${line}" lower)
    if(NOT line MATCHES "^[ \t]*#[ \t]*include[ \t]*<[^>]+>")
      _nxx_violation("numerixx/${file}: use an angle-bracket include such as <numerixx/...>: ${line}")
    elseif(line MATCHES "<[^>]*\\\\")
      _nxx_violation("numerixx/${file}: use '/' in include paths: ${line}")
    elseif(lower MATCHES "<(numerixx|eigen|eigen3|unsupported|fxt|boost)[/.]"
           AND NOT line MATCHES "<(numerixx|Eigen|eigen3|unsupported|fxt|boost)[/.]")
      _nxx_violation("numerixx/${file}: include path has the wrong case: ${line}")
    endif()
  endforeach()

  if(file STREQUAL "numerixx.hpp")
    foreach(line IN LISTS includes)
      if(line MATCHES "<numerixx/([^>]+)>")
        set(included "${CMAKE_MATCH_1}")
        _nxx_module_of("${included}" dep)
        if(dep IN_LIST UMBRELLA_FORBIDDEN OR included STREQUAL "core/any_solver.hpp")
          _nxx_violation("the umbrella numerixx.hpp must not include numerixx/${included}")
        endif()
      elseif(line MATCHES "<(Eigen|unsupported|eigen3|fxt|boost)[/.]")
        _nxx_violation("the umbrella numerixx.hpp must include only Numerixx headers: ${line}")
      endif()
    endforeach()
    continue()
  endif()

  _nxx_module_of("${file}" mod)
  if(NOT mod IN_LIST MODULES)
    _nxx_violation("numerixx/${file} belongs to unknown module '${mod}'; add it to cmake/CheckLayering.cmake")
    continue()
  endif()

  foreach(line IN LISTS includes)
    if(line MATCHES "<numerixx/([^>]+)>")
      _nxx_module_of("${CMAKE_MATCH_1}" dep)
      if(NOT (dep STREQUAL mod OR dep STREQUAL "core" OR dep IN_LIST ALLOW_${mod}))
        _nxx_violation("numerixx/${file} (module ${mod}) includes numerixx/${CMAKE_MATCH_1} (module ${dep})")
      endif()
    elseif(line MATCHES "<unsupported/")
      # Eigen's unsupported modules are not installed with Numerixx; they may be used only by tests (as oracles).
      _nxx_violation("numerixx/${file} includes Eigen's unsupported modules: ${line}")
    elseif(line MATCHES "<(Eigen|eigen3)/")
      if(NOT (mod STREQUAL "linalg" OR mod STREQUAL "multiroots" OR file STREQUAL "adapters/multiprecision_linalg.hpp"))
        _nxx_violation("numerixx/${file} includes Eigen: ${line}")
      endif()
    elseif(line MATCHES "<fxt[/.]")
      if(NOT mod STREQUAL "pipes")
        _nxx_violation("numerixx/${file} includes FXT: ${line}")
      endif()
    elseif(line MATCHES "<boost/")
      if(NOT mod STREQUAL "adapters")
        _nxx_violation("numerixx/${file} includes Boost: ${line}")
      endif()
    endif()
  endforeach()
endforeach()

list(LENGTH files count)
if(violations GREATER 0)
  message(FATAL_ERROR "layering: ${violations} violation(s) in ${count} headers")
endif()
message(STATUS "layering: OK (${count} headers)")

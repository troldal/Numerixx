# Examples, each registered as a smoke test (DESIGN §9.1). Every module phase adds its examples here.
include(NumerixxBuildMode)
include(NumerixxWarnings)

# The examples print with std::println. With MinGW's libstdc++ (still in GCC 16.1), std::print needs libstdc++exp,
# which holds std::__open_terminal and std::__write_to_terminal; the other toolchains need no library. This is a probe
# rather than a test for MINGW, so that Clang on MinGW's libstdc++ is covered too. The library never uses <print>.
function(numerixx_print_libraries out)
  include(CheckCXXSourceCompiles)
  include(CMakePushCheckState)
  set(CMAKE_CXX_STANDARD 23)    # try_compile honours it (CMP0067); the library asks for C++23 only as a target feature
  set(src "#include <print>\nint main() { std::println(\"{}\", 1); }\n")
  cmake_push_check_state(RESET)
  check_cxx_source_compiles("${src}" NUMERIXX_PRINT_LINKS)
  if(NOT NUMERIXX_PRINT_LINKS)
    set(CMAKE_REQUIRED_LIBRARIES stdc++exp)
    check_cxx_source_compiles("${src}" NUMERIXX_PRINT_LINKS_WITH_STDCXXEXP)
  endif()
  cmake_pop_check_state()
  if(NUMERIXX_PRINT_LINKS)
    set(${out} "" PARENT_SCOPE)
  elseif(NUMERIXX_PRINT_LINKS_WITH_STDCXXEXP)
    set(${out} stdc++exp PARENT_SCOPE)
  else()
    message(FATAL_ERROR "Numerixx examples: std::println (C++23 <print>) does not link with this toolchain, with or "
                        "without stdc++exp. Numerixx needs GCC 14, or Clang 19 with libc++ 19 or libstdc++ 14.3 "
                        "(DESIGN D2); a Clang-family compiler on an older libstdc++ stops here too.")
  endif()
endfunction()
numerixx_print_libraries(_nxx_print_libs)

# numerixx_add_example(<name> <source> LINK <targets...>)
function(numerixx_add_example name source)
  cmake_parse_arguments(ARG "" "" "LINK" ${ARGN})
  add_executable(numerixx_example_${name} ${source})
  target_link_libraries(numerixx_example_${name} PRIVATE ${ARG_LINK} ${_nxx_print_libs})
  numerixx_target_warnings(numerixx_example_${name})
  add_test(NAME example.${name} COMMAND numerixx_example_${name})
  set_tests_properties(example.${name} PROPERTIES LABELS example)
endfunction()

numerixx_add_example(version version.cpp LINK numerixx::numerixx)
numerixx_add_example(quick_tour quick_tour.cpp LINK numerixx::roots numerixx::deriv)

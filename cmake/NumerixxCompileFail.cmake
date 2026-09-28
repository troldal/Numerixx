# Compile-fail tests (DESIGN §9.1): illegal states must not compile, and the diagnostic must say why.
#
#   numerixx_add_compile_fail_test(<file.cpp> LINK <targets...> [EXPECT <regex>] [DELETE_REASON] [MSVC_ONLY_FAILS])
#
# Each file must contain an `#ifdef NUMERIXX_CF_CONTROL` branch that is legal code. Two tests are added:
#   cf.<name>.control  builds with NUMERIXX_CF_CONTROL defined and must succeed (guards against typos);
#   cf.<name>          builds without it and must fail with a compiler error. On GCC and Clang (clang-cl too),
#                      EXPECT must also match the first error (the error line and its notes); cl does not print
#                      deletion reasons, so there only the compiler error itself is checked.
# DELETE_REASON: the reason comes from a deleted function (NXX_DELETE), which only GCC 15+ and Clang 19+ print, so
# older compilers (the nightly floor legs) check only that the build fails.
# The number of diagnostic lines of each case is recorded in <build>/compile_fail_report.txt (DESIGN §9.1: recorded,
# not gated).
set(_NXX_RUN_COMPILE_FAIL "${CMAKE_CURRENT_LIST_DIR}/RunCompileFail.cmake")

function(numerixx_add_compile_fail_test src)
  cmake_parse_arguments(ARG "DELETE_REASON" "EXPECT" "LINK" ${ARGN})
  get_filename_component(name "${src}" NAME_WE)

  set(check_reason OFF)
  if(ARG_EXPECT AND NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(check_reason ON)
    if(ARG_DELETE_REASON)
      if((CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 15)
         OR (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 19)
         OR NOT CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang)$")
        set(check_reason OFF)
      endif()
    endif()
  endif()

  foreach(kind IN ITEMS control negative)
    set(tgt "nxx_cf_${name}_${kind}")
    add_library(${tgt} OBJECT EXCLUDE_FROM_ALL "${src}")
    target_link_libraries(${tgt} PRIVATE ${ARG_LINK})
    if(kind STREQUAL "control")
      target_compile_definitions(${tgt} PRIVATE NUMERIXX_CF_CONTROL=1)
      add_test(NAME cf.${name}.control
               COMMAND ${CMAKE_COMMAND} --build "${CMAKE_BINARY_DIR}" --target ${tgt} --config $<CONFIG>)
      set(test_name cf.${name}.control)
    else()
      add_test(NAME cf.${name}
               COMMAND ${CMAKE_COMMAND}
                       -DBUILD_DIR=${CMAKE_BINARY_DIR} -DTARGET=${tgt} -DCONFIG=$<CONFIG>
                       -DEXPECT=${ARG_EXPECT} -DCHECK_REASON=${check_reason}
                       -DREPORT=${CMAKE_BINARY_DIR}/compile_fail_report.txt
                       -P "${_NXX_RUN_COMPILE_FAIL}")
      set(test_name cf.${name})
    endif()
    set_tests_properties(${test_name} PROPERTIES LABELS "compile-fail" RESOURCE_LOCK numerixx_build_tree)
  endforeach()
endfunction()

# numerixx_add_msvc_escalation_probe(<file.cpp> LINK <targets...>)
# The P2564 (consteval escalation) probe of DESIGN §6.2: legal C++23, which GCC, Clang and clang-cl compile and cl
# rejects (C7595). The test builds it and is WILL_FAIL on cl only, documenting the MSVC gap.
function(numerixx_add_msvc_escalation_probe src)
  cmake_parse_arguments(ARG "" "" "LINK" ${ARGN})
  get_filename_component(name "${src}" NAME_WE)
  set(tgt "nxx_probe_${name}")
  add_library(${tgt} OBJECT EXCLUDE_FROM_ALL "${src}")
  target_link_libraries(${tgt} PRIVATE ${ARG_LINK})
  add_test(NAME probe.${name} COMMAND ${CMAKE_COMMAND} --build "${CMAKE_BINARY_DIR}" --target ${tgt} --config $<CONFIG>)
  set_tests_properties(probe.${name} PROPERTIES LABELS "compile-fail" RESOURCE_LOCK numerixx_build_tree)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set_tests_properties(probe.${name} PROPERTIES WILL_FAIL TRUE)
  endif()
endfunction()

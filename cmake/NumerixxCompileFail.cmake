# Compile-fail tests (DESIGN §9.1): illegal states must not compile, and the diagnostic must say why.
#
#   numerixx_add_compile_fail_test(<file.cpp> LINK <targets...> [EXPECT <regex>])
#
# Each file must contain an `#ifdef NUMERIXX_CF_CONTROL` branch that is legal code. Two tests are added:
#   cf.<name>.control  builds with NUMERIXX_CF_CONTROL defined and must succeed (guards against typos);
#   cf.<name>          builds without it and must fail with a compiler error. On GCC and Clang (clang-cl too),
#                      EXPECT must also match the first error (the error line and its notes); cl does not print
#                      deletion reasons, so there only the compiler error itself is checked.
set(_NXX_RUN_COMPILE_FAIL "${CMAKE_CURRENT_LIST_DIR}/RunCompileFail.cmake")

function(numerixx_add_compile_fail_test src)
  cmake_parse_arguments(ARG "" "EXPECT" "LINK" ${ARGN})
  get_filename_component(name "${src}" NAME_WE)

  set(check_reason OFF)
  if(ARG_EXPECT AND NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(check_reason ON)
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
                       -P "${_NXX_RUN_COMPILE_FAIL}")
      set(test_name cf.${name})
    endif()
    set_tests_properties(${test_name} PROPERTIES LABELS "compile-fail" RESOURCE_LOCK numerixx_build_tree)
  endforeach()
endfunction()

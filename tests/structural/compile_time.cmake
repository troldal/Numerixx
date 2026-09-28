# Compile-time measurements (DESIGN §3.6, spike exit criterion 5). Included by tests/CMakeLists.txt.
#
#   structural.compile_time.umbrella  a TU that includes <numerixx/numerixx.hpp> and instantiates nothing. On GCC the
#                                     test fails when its best time exceeds 2 s (the guard); elsewhere it is recorded.
#   structural.compile_time.linalg    a TU that instantiates Eigen solves through <numerixx/linalg.hpp>: recorded, not
#                                     gated.
#
# Each target is EXCLUDE_FROM_ALL, so a normal build never compiles it; the test removes its object file and rebuilds
# it three times (cmake/MeasureCompileTime.cmake), and appends the best time to <build>/compile_time_report.txt. The
# tests rebuild targets of this build tree, so they share the build-tree lock with the compile-fail tests.
set(_nxx_measure_compile_time "${PROJECT_SOURCE_DIR}/cmake/MeasureCompileTime.cmake")
set(_nxx_compile_time_report "${CMAKE_BINARY_DIR}/compile_time_report.txt")

# numerixx_add_compile_time_test(<name> <source> LINK <targets...> [GUARD_SECONDS <s>])
function(numerixx_add_compile_time_test name src)
  cmake_parse_arguments(ARG "" "GUARD_SECONDS" "LINK" ${ARGN})
  set(tgt numerixx_compile_time_${name})
  add_library(${tgt} OBJECT EXCLUDE_FROM_ALL "${src}")
  target_link_libraries(${tgt} PRIVATE ${ARG_LINK})
  add_test(NAME structural.compile_time.${name}
           COMMAND ${CMAKE_COMMAND}
                   -DBUILD_DIR=${CMAKE_BINARY_DIR} -DTARGET=${tgt} -DCONFIG=$<CONFIG>
                   "-DOBJECTS=$<TARGET_OBJECTS:${tgt}>" -DNAME=${name}
                   -DCOMPILER_ID=${CMAKE_CXX_COMPILER_ID} -DCOMPILER_VERSION=${CMAKE_CXX_COMPILER_VERSION}
                   -DCOMPILER_FRONTEND=${CMAKE_CXX_COMPILER_FRONTEND_VARIANT}
                   -DGUARD_SECONDS=${ARG_GUARD_SECONDS} -DGUARD_COMPILER=GNU
                   -DREPORT=${_nxx_compile_time_report}
                   -P "${_nxx_measure_compile_time}")
  set_tests_properties(structural.compile_time.${name} PROPERTIES
                       LABELS "compile-time" RESOURCE_LOCK numerixx_build_tree)
endfunction()

numerixx_add_compile_time_test(umbrella "${CMAKE_CURRENT_LIST_DIR}/compile_time_umbrella.cpp"
                               LINK numerixx::numerixx GUARD_SECONDS 2.0)

if(TARGET numerixx::linalg)
  numerixx_add_compile_time_test(linalg "${CMAKE_CURRENT_LIST_DIR}/compile_time_linalg.cpp" LINK numerixx::linalg)
endif()

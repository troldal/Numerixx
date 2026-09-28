# Examples, each registered as a smoke test (DESIGN §9.1). Every module phase adds its examples here.
include(NumerixxBuildMode)
include(NumerixxWarnings)

# numerixx_add_example(<name> <source> LINK <targets...>)
function(numerixx_add_example name source)
  cmake_parse_arguments(ARG "" "" "LINK" ${ARGN})
  add_executable(numerixx_example_${name} ${source})
  target_link_libraries(numerixx_example_${name} PRIVATE ${ARG_LINK})
  numerixx_target_warnings(numerixx_example_${name})
  add_test(NAME example.${name} COMMAND numerixx_example_${name})
  set_tests_properties(example.${name} PROPERTIES LABELS example)
endfunction()

numerixx_add_example(version version.cpp LINK numerixx::numerixx)

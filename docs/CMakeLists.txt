# The documentation targets, run on demand and never part of the default build. Each runs docs/Doxyfile from the
# project root, through a small Doxyfile here that includes it and moves its output into this build directory.
function(numerixx_add_docs target output sections)
    set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/${target}.doxyfile)
    file(WRITE ${doxyfile}
            "@INCLUDE = \"${PROJECT_SOURCE_DIR}/docs/Doxyfile\"\n"
            "OUTPUT_DIRECTORY = \"${output}\"\n"
            "ENABLED_SECTIONS = ${sections}\n")
    add_custom_target(${target}
            COMMAND Doxygen::doxygen ${doxyfile}
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            COMMENT "Generating ${output}/html/index.html"
            VERBATIM)
endfunction()

numerixx_add_docs(numerixx_docs ${CMAKE_CURRENT_BINARY_DIR}/public "")
numerixx_add_docs(numerixx_docs_internal ${CMAKE_CURRENT_BINARY_DIR}/internal INTERNAL)

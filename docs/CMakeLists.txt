# The documentation targets, run on demand and never part of the default build (docs/sphinx/PHASE2_PLAN.md):
# - numerixx_docs: the public site. Doxygen writes XML without the detail namespaces, then Sphinx renders it with
#   Breathe and MyST (docs/sphinx), through uv with the pinned requirements, and fails on any warning or check;
# - numerixx_docs_internal: Doxygen's HTML with the detail namespaces, for developers.
# Each runs docs/Doxyfile from the project root, through a small Doxyfile here that includes it.
function(numerixx_write_doxyfile target)
    set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/${target}.doxyfile)
    list(JOIN ARGN "\n" settings)
    file(WRITE ${doxyfile} "@INCLUDE = \"${PROJECT_SOURCE_DIR}/docs/Doxyfile\"\n${settings}\n")
    set(${target}_doxyfile ${doxyfile} PARENT_SCOPE)
endfunction()

numerixx_write_doxyfile(numerixx_docs_internal
        "OUTPUT_DIRECTORY = \"${CMAKE_CURRENT_BINARY_DIR}/internal\""
        "ENABLED_SECTIONS = INTERNAL")
add_custom_target(numerixx_docs_internal
        COMMAND Doxygen::doxygen ${numerixx_docs_internal_doxyfile}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMENT "Generating ${CMAKE_CURRENT_BINARY_DIR}/internal/html/index.html"
        VERBATIM)

# The Breathe patch in docs/sphinx/conf.py was verified on Doxygen 1.18's XML.
find_program(NUMERIXX_UV uv)
if (DOXYGEN_VERSION VERSION_LESS 1.18 OR NOT NUMERIXX_UV)
    message(STATUS "numerixx_docs needs Doxygen 1.18 or later (found ${DOXYGEN_VERSION}) and uv: not added")
    return()
endif ()

set(xml ${CMAKE_CURRENT_BINARY_DIR}/public/xml)
set(html ${CMAKE_CURRENT_BINARY_DIR}/public/html)
# The only [[deprecated]] declarations are derivative's low-precision warnings, not deprecations: no list of them.
numerixx_write_doxyfile(numerixx_docs
        "OUTPUT_DIRECTORY = \"${CMAKE_CURRENT_BINARY_DIR}/public\""
        "GENERATE_HTML = NO"
        "GENERATE_XML = YES"
        "XML_PROGRAMLISTING = NO"
        "EXCLUDE_SYMBOLS = *::detail *::detail::*"
        "GENERATE_DEPRECATEDLIST = NO")
# -E re-reads every page, so that -W sees every warning on a rebuild too.
add_custom_target(numerixx_docs
        COMMAND Doxygen::doxygen ${numerixx_docs_doxyfile}
        COMMAND ${CMAKE_COMMAND} -E env NUMERIXX_DOXYGEN_XML=${xml}
                ${NUMERIXX_UV} run --no-project --python 3.14
                --with-requirements ${PROJECT_SOURCE_DIR}/docs/sphinx/requirements.txt
                python -m sphinx -E -W -b html ${PROJECT_SOURCE_DIR}/docs/sphinx ${html}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMENT "Generating ${html}/index.html"
        VERBATIM)

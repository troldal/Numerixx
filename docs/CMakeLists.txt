# API documentation (NUMERIXX_BUILD_DOCS, OFF by default). Skipped with a message when Doxygen is missing.
# The user documentation is rewritten in phase 9 (DESIGN §10.3); docs/redesign holds the design documents.
find_package(Doxygen)
if(NOT DOXYGEN_FOUND)
  message(STATUS "Numerixx: Doxygen not found; the numerixx_docs target is not available")
  return()
endif()

set(DOXYGEN_PROJECT_NAME "Numerixx")
set(DOXYGEN_PROJECT_NUMBER "${PROJECT_VERSION}")
set(DOXYGEN_EXTRACT_ALL YES)
set(DOXYGEN_RECURSIVE YES)
set(DOXYGEN_GENERATE_LATEX NO)
set(DOXYGEN_STRIP_FROM_INC_PATH "${PROJECT_SOURCE_DIR}/include")
doxygen_add_docs(numerixx_docs "${PROJECT_SOURCE_DIR}/include/numerixx"
                 COMMENT "Generating the Numerixx API documentation")

find_package(Doxygen REQUIRED)
find_package(Sphinx REQUIRED)

# Find all the public headers
get_target_property(NUMERIXX_PUBLIC_HEADER_DIR nxx_all INTERFACE_INCLUDE_DIRECTORIES)
file(GLOB_RECURSE NUMERIXX_PUBLIC_HEADERS ${NUMERIXX_PUBLIC_HEADER_DIR}/*.hpp)

set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/numerixx)
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/docs/doxygen)
set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/html/index.html)
set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

#Replace variables inside @@ with the current values
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

# Create the Doxygen output directory, and define a custom command that will run the Doxygen application
# every time the index.html file is accessed by subsequent target builds.
file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR}) #Doxygen won't create this for us
add_custom_command(
        OUTPUT ${DOXYGEN_INDEX_FILE}
        DEPENDS ${NUMERIXX_PUBLIC_HEADERS}
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
        MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
        COMMENT "Generating docs")

# Target for running the Doxygen application. When run, it will also run the previous custom command.
add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})

set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/docs/sphinx)

add_custom_target(Sphinx ALL DEPENDS ${DOXYGEN_INDEX_FILE}
        # Tell Breathe where to find the Doxygen output
        COMMAND ${SPHINX_EXECUTABLE} -b html -Dbreathe_projects.Numerixx=${DOXYGEN_OUTPUT_DIR} ${SPHINX_SOURCE} ${SPHINX_BUILD}
        COMMAND ${SPHINX_EXECUTABLE} -b latex -Dbreathe_projects.numerix=${DOXYGEN_OUTPUT_DIR} ${SPHINX_SOURCE} ${SPHINX_BUILD}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating documentation with Sphinx")

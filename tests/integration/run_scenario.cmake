# Runs one consumer-build scenario: configure, build and run the parent project in SCENARIO_DIR. With
# MODE=installed, Numerixx is first configured, built and installed into a staging prefix, and the parent finds it
# with find_package. Invoked by tests/integration/CMakeLists.txt; see there for the variables.
cmake_minimum_required(VERSION 3.25)

function(nxx_run)
  string(REPLACE ";" " " printable "${ARGN}")
  message(STATUS "[${SCENARIO}] ${printable}")
  execute_process(COMMAND ${ARGN} RESULT_VARIABLE result)
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "[${SCENARIO}] command failed (${result}): ${printable}")
  endif()
endfunction()

file(REMOVE_RECURSE "${BINARY_DIR}")

# --no-warn-unused-cli: every scenario receives the same variables, and not every scenario uses all of them.
set(common -G "${GENERATOR}" --no-warn-unused-cli "-DCMAKE_CXX_COMPILER=${CXX_COMPILER}" "-DCMAKE_BUILD_TYPE=${BUILD_TYPE}"
           "-DCPM_SOURCE_CACHE=${CPM_SOURCE_CACHE}")
if(MAKE_PROGRAM)
  list(APPEND common "-DCMAKE_MAKE_PROGRAM=${MAKE_PROGRAM}")
endif()

set(parent_args "-DNUMERIXX_SOURCE_DIR=${NUMERIXX_SOURCE_DIR}"
                "-DFXT_REF=${FXT_REF}" "-DFXT_SHA256=${FXT_SHA256}"
                "-DEIGEN_VERSION=${EIGEN_VERSION}" "-DEIGEN_SHA256=${EIGEN_SHA256}")
set(parent_binary "${BINARY_DIR}")

if(MODE STREQUAL "installed")
  set(stage "${BINARY_DIR}/stage")
  nxx_run(${CMAKE_COMMAND} -S "${NUMERIXX_SOURCE_DIR}" -B "${BINARY_DIR}/numerixx" ${common}
          -DNUMERIXX_BUILD_TESTS=OFF -DNUMERIXX_INSTALL=ON "-DCMAKE_INSTALL_PREFIX=${stage}")
  nxx_run(${CMAKE_COMMAND} --build "${BINARY_DIR}/numerixx" --config "${BUILD_TYPE}")
  nxx_run(${CMAKE_COMMAND} --install "${BINARY_DIR}/numerixx" --config "${BUILD_TYPE}")
  list(APPEND parent_args "-DCMAKE_PREFIX_PATH=${stage}")
  set(parent_binary "${BINARY_DIR}/consumer")
endif()

nxx_run(${CMAKE_COMMAND} -S "${SCENARIO_DIR}" -B "${parent_binary}" ${common} ${parent_args})
nxx_run(${CMAKE_COMMAND} --build "${parent_binary}" --config "${BUILD_TYPE}")

file(GLOB_RECURSE executables LIST_DIRECTORIES false "${parent_binary}/consumer" "${parent_binary}/consumer.exe")
list(FILTER executables EXCLUDE REGEX "/CMakeFiles/")
if(NOT executables)
  message(FATAL_ERROR "[${SCENARIO}] the parent built no 'consumer' executable")
endif()
list(GET executables 0 executable)
nxx_run("${executable}")
message(STATUS "[${SCENARIO}] OK")

# Numerixx CPM bootstrap (DESIGN §4.2, D25).
#
# Reuse a CPM that a parent project already loaded (the first loaded CPM wins; parents on CPM 0.42.x work),
# otherwise load the pinned v0.43.2 through the committed get_cpm.cmake, which checks the download's SHA256.
if(COMMAND CPMAddPackage)
  message(STATUS "Numerixx: reusing the parent's CPM ${CURRENT_CPM_VERSION}")
  return()
endif()
include("${CMAKE_CURRENT_LIST_DIR}/get_cpm.cmake")

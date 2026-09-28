# Install and export rules (DESIGN §4.3, NUMERIXX_INSTALL): find_package(numerixx 2.0 CONFIG).
#
# FXT and Eigen that Numerixx fetched itself are installed with it, below include/numerixx-deps (see
# NumerixxDependencies.cmake), together with their licences. A dependency that a parent provided must be findable
# through its own package (Eigen3). Every check runs before the first install() rule is registered: if the package
# cannot be exported completely, no rule is registered at all and a warning says why.
include(CMakePackageConfigHelpers)

set(_nxx_config_dir "${CMAKE_INSTALL_DATADIR}/numerixx")
get_property(_nxx_export GLOBAL PROPERTY NUMERIXX_MODULES)
list(PREPEND _nxx_export numerixx)
set(NUMERIXX_CONFIG_FIND_EIGEN OFF)
set(NUMERIXX_CONFIG_FIND_BOOST OFF)

# ---- Checks ------------------------------------------------------------------------------------------------------
if(TARGET numerixx::pipes AND NOT TARGET numerixx_fxt)
  message(WARNING "Numerixx: install rules disabled: numerixx::pipes uses a parent-provided fxt::fxt, which "
                  "cannot be exported. Set NUMERIXX_INSTALL=OFF, or let Numerixx fetch FXT.")
  return()
endif()

if(TARGET numerixx::linalg AND NOT TARGET numerixx_eigen)
  get_target_property(_nxx_eigen_imported Eigen3::Eigen IMPORTED)
  if(NOT _nxx_eigen_imported)
    message(WARNING "Numerixx: install rules disabled: numerixx::linalg uses a parent-provided, non-imported "
                    "Eigen3::Eigen, which cannot be exported. Set NUMERIXX_INSTALL=OFF.")
    return()
  endif()
  set(NUMERIXX_CONFIG_FIND_EIGEN ON)   # numerixx-config.cmake calls find_dependency(Eigen3)
endif()

# The adapter is exported only when Boost comes from an installed package; over the standalone repositories that
# Numerixx fetched itself it stays a build-tree target.
if(TARGET numerixx::multiprecision)
  get_target_property(_nxx_boost_imported Boost::multiprecision IMPORTED)
  if(NOT _nxx_boost_imported)
    list(REMOVE_ITEM _nxx_export numerixx_multiprecision)
    message(STATUS "Numerixx: numerixx::multiprecision is not installed (Boost was fetched, not found)")
  else()
    set(NUMERIXX_CONFIG_FIND_BOOST ON)  # numerixx-config.cmake finds Boost.Config and Boost.Multiprecision again
  endif()
endif()

# ---- Rules -------------------------------------------------------------------------------------------------------
if(TARGET numerixx::pipes)
  list(APPEND _nxx_export numerixx_fxt)            # exported as numerixx::fxt; numerixx-config.cmake adds fxt::fxt
  install(DIRECTORY "${FXT_SOURCE_DIR}/include/" DESTINATION ${NUMERIXX_INSTALL_DEPS_INCLUDEDIR}/fxt)
  install(FILES "${FXT_SOURCE_DIR}/LICENSE" DESTINATION ${CMAKE_INSTALL_DOCDIR}/fxt)
endif()

if(TARGET numerixx_eigen)
  list(APPEND _nxx_export numerixx_eigen)          # exported as numerixx::eigen
  install(DIRECTORY "${Eigen3_SOURCE_DIR}/Eigen" DESTINATION ${NUMERIXX_INSTALL_DEPS_INCLUDEDIR}/eigen3)
  # Eigen is MPL-2.0; some of its files carry BSD, Apache or MINPACK notices (see COPYING.README).
  file(GLOB _nxx_eigen_licences "${Eigen3_SOURCE_DIR}/COPYING.*")
  install(FILES ${_nxx_eigen_licences} DESTINATION ${CMAKE_INSTALL_DOCDIR}/eigen)
endif()

install(TARGETS ${_nxx_export} EXPORT numerixx-targets)
install(DIRECTORY "${PROJECT_SOURCE_DIR}/include/numerixx" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES "${PROJECT_SOURCE_DIR}/LICENSE" DESTINATION ${CMAKE_INSTALL_DOCDIR})
install(EXPORT numerixx-targets NAMESPACE numerixx:: DESTINATION ${_nxx_config_dir})

configure_package_config_file("${CMAKE_CURRENT_LIST_DIR}/numerixx-config.cmake.in"
  "${PROJECT_BINARY_DIR}/numerixx-config.cmake"
  INSTALL_DESTINATION ${_nxx_config_dir})
write_basic_package_version_file("${PROJECT_BINARY_DIR}/numerixx-config-version.cmake"
  COMPATIBILITY SameMinorVersion
  ARCH_INDEPENDENT)
install(FILES "${PROJECT_BINARY_DIR}/numerixx-config.cmake" "${PROJECT_BINARY_DIR}/numerixx-config-version.cmake"
        DESTINATION ${_nxx_config_dir})

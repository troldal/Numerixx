# Numerixx library dependencies (DESIGN §4.2). Development-only dependencies (doctest, google/benchmark,
# the Boost.Math oracles) are fetched where they are used, under tests/ and benchmarks/.
#
# Every package is pinned to an exact commit or version plus a SHA256, and declared DOWNLOAD_ONLY with an own
# INTERFACE target where the upstream CMake would do more than we need. A target that a parent project already
# provides is reused. A local checkout can replace any package through CPM_<Name>_SOURCE, e.g.
#   -DCPM_FXT_SOURCE=C:/Dev/XLThermo/FXT
#
# When Numerixx is installed, the FXT and Eigen headers it fetched are installed with it, below a private directory
# (not include/eigen3 or include/fxt), so they can never overwrite another installation of those libraries.
set(NUMERIXX_INSTALL_DEPS_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/numerixx-deps")

# ---------------------------------------------------------------------------------------------------------------
# FXT: only for numerixx::pipes (D18, D26).
# DOWNLOAD_ONLY, because FXT's own CMake downloads an unhashed CPM and the TartanLlama expected/optional
# repositories. Parents that use FXT themselves declare it as `NAME FXT` too, so CPM deduplicates it.
# ---------------------------------------------------------------------------------------------------------------
set(NUMERIXX_FXT_REF    "9208e597156b196dd5cbf55ab8460ffb2cead30b"
    CACHE STRING "FXT commit fetched by Numerixx (FXT has no release tags yet)")
set(NUMERIXX_FXT_SHA256 "dcd46f5cce7973ecdda69088cab4f2e364fbdbbc285e03bc7eee5799b689397d"
    CACHE STRING "SHA256 of the GitHub archive of NUMERIXX_FXT_REF")

if(NUMERIXX_WITH_FXT)
  if(TARGET fxt::fxt)
    message(STATUS "Numerixx: reusing the parent's fxt::fxt")
    get_target_property(_nxx_fxt_defs fxt::fxt INTERFACE_COMPILE_DEFINITIONS)
    if(_nxx_fxt_defs MATCHES "FXT_USE_TL_(EXPECTED|OPTIONAL)")
      message(WARNING "Numerixx: the parent's fxt::fxt uses tl::expected/tl::optional; Numerixx's public "
                      "signatures use std::expected, so the FXT pipes may not apply to Numerixx results.")
    endif()
  else()
    CPMAddPackage(
      NAME FXT
      URL https://github.com/troldal/FXT/archive/${NUMERIXX_FXT_REF}.tar.gz
      URL_HASH SHA256=${NUMERIXX_FXT_SHA256}
      DOWNLOAD_ONLY YES
    )
    if(NOT TARGET fxt::fxt)
      add_library(numerixx_fxt INTERFACE)
      set_target_properties(numerixx_fxt PROPERTIES EXPORT_NAME fxt)
      target_include_directories(numerixx_fxt SYSTEM INTERFACE $<BUILD_INTERFACE:${FXT_SOURCE_DIR}/include>)
      # Not SYSTEM: imported targets are system headers anyway, and CMake 3.x exports a relative (broken) path for
      # SYSTEM install interfaces.
      target_include_directories(numerixx_fxt INTERFACE $<INSTALL_INTERFACE:${NUMERIXX_INSTALL_DEPS_INCLUDEDIR}/fxt>)
      target_compile_features(numerixx_fxt INTERFACE cxx_std_23)
      add_library(fxt::fxt ALIAS numerixx_fxt)
    endif()
  endif()
endif()

# ---------------------------------------------------------------------------------------------------------------
# Eigen 5.0.1 (MPL-2.0, header-only, no BLAS/LAPACK): the backend of numerixx::linalg and numerixx::multiroots
# (D22). Named Eigen3, like Eigen's own CMake package, so that CPM deduplicates it with parents that fetch Eigen.
# If GitLab ever regenerates the archive and the hash stops matching, use GIT_TAG 5.0.1 instead of the URL.
# ---------------------------------------------------------------------------------------------------------------
set(NUMERIXX_EIGEN_VERSION "5.0.1")
set(NUMERIXX_EIGEN_SHA256  "e9c326dc8c05cd1e044c71f30f1b2e34a6161a3b6ecf445d56b53ff1669e3dec"
    CACHE STRING "SHA256 of the Eigen ${NUMERIXX_EIGEN_VERSION} archive")

set(NUMERIXX_CONFIG_EIGEN_VERSION "")   # the version numerixx-config.cmake asks find_dependency(Eigen3) for
if(NUMERIXX_WITH_LINALG)
  if(TARGET Eigen3::Eigen)
    message(STATUS "Numerixx: reusing the parent's Eigen3::Eigen")
    set(NUMERIXX_EIGEN_TARGET Eigen3::Eigen)
  else()
    CPMAddPackage(
      NAME Eigen3
      VERSION ${NUMERIXX_EIGEN_VERSION}
      URL https://gitlab.com/libeigen/eigen/-/archive/${NUMERIXX_EIGEN_VERSION}/eigen-${NUMERIXX_EIGEN_VERSION}.tar.gz
      URL_HASH SHA256=${NUMERIXX_EIGEN_SHA256}
      DOWNLOAD_ONLY YES
    )
    if(TARGET Eigen3::Eigen)
      # CPM used an installed Eigen (CPM_USE_LOCAL_PACKAGES or CPM_LOCAL_PACKAGES_ONLY): use its imported target.
      # CPM now considers Eigen3 added and deduplicates a parent's later declaration, so the parent must see this
      # target too; imported targets are directory-scoped unless made global.
      get_target_property(_nxx_eigen_aliased Eigen3::Eigen ALIASED_TARGET)
      get_target_property(_nxx_eigen_global Eigen3::Eigen IMPORTED_GLOBAL)
      if(NOT _nxx_eigen_aliased AND NOT _nxx_eigen_global)
        set_target_properties(Eigen3::Eigen PROPERTIES IMPORTED_GLOBAL TRUE)
      endif()
      set(NUMERIXX_EIGEN_TARGET Eigen3::Eigen)
      set(NUMERIXX_CONFIG_EIGEN_VERSION ${NUMERIXX_EIGEN_VERSION})
    else()
      add_library(numerixx_eigen INTERFACE)
      set_target_properties(numerixx_eigen PROPERTIES EXPORT_NAME eigen)
      target_include_directories(numerixx_eigen SYSTEM INTERFACE $<BUILD_INTERFACE:${Eigen3_SOURCE_DIR}>)
      target_include_directories(numerixx_eigen INTERFACE
        $<INSTALL_INTERFACE:${NUMERIXX_INSTALL_DEPS_INCLUDEDIR}/eigen3>)   # not SYSTEM: see numerixx_fxt above
      # A parent that declares Eigen3 after Numerixx is deduplicated by CPM; it finds Eigen's usual target here.
      add_library(Eigen3::Eigen ALIAS numerixx_eigen)
      set(NUMERIXX_EIGEN_TARGET numerixx_eigen)
    endif()
  endif()
endif()

# ---------------------------------------------------------------------------------------------------------------
# Boost.Config + Boost.Multiprecision 1.92.0, standalone boostorg repositories: only for the optional
# numerixx::multiprecision adapter (D16, D27). Never declared under the package name "Boost", which parent
# projects commonly use for the full Boost distribution.
# ---------------------------------------------------------------------------------------------------------------
set(NUMERIXX_BOOST_TAG "boost-1.92.0")

if(NUMERIXX_WITH_MULTIPRECISION AND NOT TARGET Boost::multiprecision)
  if(NOT TARGET Boost::config)
    CPMAddPackage(
      NAME boost_config
      VERSION 1.92.0
      URL https://github.com/boostorg/config/archive/refs/tags/${NUMERIXX_BOOST_TAG}.tar.gz
      URL_HASH SHA256=b4171037f13373203ba79cbc141d612982052283e696a315185ab5bea46102a0
      EXCLUDE_FROM_ALL YES
      SYSTEM YES
    )
  endif()
  CPMAddPackage(
    NAME boost_multiprecision
    VERSION 1.92.0
    URL https://github.com/boostorg/multiprecision/archive/refs/tags/${NUMERIXX_BOOST_TAG}.tar.gz
    URL_HASH SHA256=9da997843edd802f9d5e25d53872a7676330f31a52a5e47a418d1c94ba3b01d0
    EXCLUDE_FROM_ALL YES
    SYSTEM YES
    OPTIONS "BOOST_MP_STANDALONE ON"
  )
endif()

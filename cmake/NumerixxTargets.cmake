# Numerixx targets (DESIGN §5.2, D24): one include root, one thin INTERFACE target per module, named
# numerixx::<module> (namespace = header = target). The module DAG below is also enforced on the headers by
# cmake/CheckLayering.cmake.

# Consumers see Numerixx's headers as system headers, so warnings inside them never reach a consumer's build.
# In Numerixx's own top-level build they are ordinary include paths, so our tests see every warning. Installed
# (imported) targets are treated as system headers by CMake anyway.
if(PROJECT_IS_TOP_LEVEL)
  set(_nxx_include_kind "")
else()
  set(_nxx_include_kind SYSTEM)
endif()

set_property(GLOBAL PROPERTY NUMERIXX_MODULES "")

# numerixx_add_module(<name> [DEPS <targets...>])
function(numerixx_add_module name)
  cmake_parse_arguments(ARG "" "" "DEPS" ${ARGN})
  add_library(numerixx_${name} INTERFACE)
  add_library(numerixx::${name} ALIAS numerixx_${name})
  set_target_properties(numerixx_${name} PROPERTIES EXPORT_NAME ${name})
  target_link_libraries(numerixx_${name} INTERFACE ${ARG_DEPS})
  set_property(GLOBAL APPEND PROPERTY NUMERIXX_MODULES numerixx_${name})
endfunction()

# core: the standard library only.
numerixx_add_module(core)
target_include_directories(numerixx_core ${_nxx_include_kind} INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>)
# Not SYSTEM: imported targets are system headers anyway, and CMake 3.x exports a relative path for SYSTEM install
# interfaces.
target_include_directories(numerixx_core INTERFACE $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_compile_features(numerixx_core INTERFACE cxx_std_23)

# Scalar modules: core only, so a user of these never downloads Eigen or FXT.
numerixx_add_module(deriv       DEPS numerixx::core)
numerixx_add_module(roots       DEPS numerixx::core)
numerixx_add_module(optimize    DEPS numerixx::core)
numerixx_add_module(poly        DEPS numerixx::core)
numerixx_add_module(integrate   DEPS numerixx::core)
numerixx_add_module(interpolate DEPS numerixx::core)   # in-house O(n) tridiagonal solvers

# Eigen-backed modules.
if(NUMERIXX_WITH_LINALG)
  numerixx_add_module(linalg     DEPS numerixx::core ${NUMERIXX_EIGEN_TARGET})
  numerixx_add_module(multiroots DEPS numerixx::linalg numerixx::deriv)
endif()

# The only FXT consumer.
if(NUMERIXX_WITH_FXT)
  numerixx_add_module(pipes DEPS numerixx::core fxt::fxt)
endif()

# Optional adapter; never part of numerixx::numerixx.
if(NUMERIXX_WITH_MULTIPRECISION)
  set(_nxx_mp_deps numerixx::core Boost::multiprecision Boost::config)
  if(TARGET numerixx::linalg)
    list(APPEND _nxx_mp_deps numerixx::linalg)   # for adapters/multiprecision_linalg.hpp
  endif()
  numerixx_add_module(multiprecision DEPS ${_nxx_mp_deps})
endif()

# Everything but the adapters. Linking it adds no include: the umbrella header <numerixx/numerixx.hpp> never
# includes linalg or multiroots, so Eigen's compile cost is paid only where they are included explicitly.
set(_nxx_all numerixx::deriv numerixx::roots numerixx::optimize numerixx::poly numerixx::integrate
             numerixx::interpolate)
foreach(_nxx_optional IN ITEMS multiroots pipes)
  if(TARGET numerixx::${_nxx_optional})
    list(APPEND _nxx_all numerixx::${_nxx_optional})
  endif()
endforeach()
add_library(numerixx INTERFACE)
add_library(numerixx::numerixx ALIAS numerixx)
target_link_libraries(numerixx INTERFACE ${_nxx_all})

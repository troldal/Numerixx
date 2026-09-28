# Build mode of Numerixx-owned directories (tests/, examples/, benchmarks/). Included at the top of each of them,
# before any dependency is fetched there, so that fetched dependencies (doctest, google/benchmark) are compiled the
# same way: every translation unit of one program must use the same exception model (DESIGN §4.1, principle 5).
# The library targets themselves add no such flags; the exception model belongs to the consumer.
include_guard(DIRECTORY)

# Numerixx uses no C++20 modules. CMake >= 3.28 scans C++20-and-later sources for them by default, which costs time
# and fails with GCC on Windows (its module mapper cannot read CMake's backslash paths).
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

# Test the library in strict ISO C++23 (-std=c++23, not gnu++23): it must also build with MSVC and clang-cl.
set(CMAKE_CXX_EXTENSIONS OFF)

if(NUMERIXX_NO_EXCEPTIONS)
  if(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
    string(REGEX REPLACE "/EH[a-z]+" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    add_compile_options(/EHs-c-)
    add_compile_definitions(_HAS_EXCEPTIONS=0)
  else()
    add_compile_options(-fno-exceptions)
    if(EMSCRIPTEN)
      add_link_options(-fno-exceptions)
    endif()
  endif()
endif()

if(NUMERIXX_SANITIZE)
  if(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
    message(WARNING "Numerixx: NUMERIXX_SANITIZE is not supported with MSVC-style compilers; ignored")
  else()
    list(JOIN NUMERIXX_SANITIZE "," _nxx_sanitizers)
    add_compile_options(-fsanitize=${_nxx_sanitizers} -fno-omit-frame-pointer -fno-sanitize-recover=all)
    add_link_options(-fsanitize=${_nxx_sanitizers})
  endif()
endif()

if(EMSCRIPTEN)
  # Tests and examples run under node. NODERAWFS: tests see the host file system (test data, report files).
  # EXIT_RUNTIME: the exit code reaches CTest.
  # -Wno-pthreads-mem-growth: with -pthread, memory growth makes non-wasm access slower, which the tests accept.
  add_link_options(-sALLOW_MEMORY_GROWTH=1 -sSTACK_SIZE=1MB -sEXIT_RUNTIME=1 -sNODERAWFS=1 -Wno-pthreads-mem-growth)
endif()

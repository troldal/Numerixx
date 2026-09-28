# Numerixx

A general-purpose, header-only numerical library for C++23, in the spirit of GSL but with a smaller scope:
numerical differentiation, root finding, minimisation, polynomials, systems of nonlinear equations, quadrature and
interpolation, done carefully, generically and composably.

- **Solvers are values.** A configured solver is a small immutable object; solvers compose (`first_of` tries the
  next solver if one fails, `then` stages search, solve and polish).
- **No exceptions.** Every result is a `std::expected` that carries the estimate, or the best estimate so far with
  the reason for the failure and the iteration and evaluation counts.
- **Illegal states do not compile** where the type system can express it, and are validated once at the boundary
  where it cannot.
- **Portable.** GCC, Clang, MSVC, clang-cl and Emscripten (WebAssembly). MIT licence.

> **Status: Numerixx 2 is being rebuilt.** This branch contains the new build system and an empty library skeleton
> (roadmap phase 0). The modules are ported phase by phase; see [the plan](docs/redesign/PLAN.md) and
> [the design reference](docs/redesign/DESIGN.md). The previous API is preserved at the tags `v1.0.0` (master) and
> `v1.1.0-legacy` (the last development branch); [MIGRATION.md](MIGRATION.md) maps it to the new one.

## Requirements

- A C++23 compiler: GCC ≥ 14, Clang ≥ 18 with libc++, MSVC 19.51 (`/std:c++latest`), clang-cl, or Emscripten ≥ 6.0.8.
- CMake ≥ 3.25 (≥ 3.30 recommended on Windows) and Ninja.

Dependencies are fetched by [CPM.cmake](https://github.com/cpm-cmake/CPM.cmake), pinned by version and SHA256.
The scalar modules need only the standard library.

| Dependency | Needed by | Option (default) |
|---|---|---|
| [FXT](https://github.com/troldal/FXT) | `numerixx::pipes` (FXT pipe syntax for results) | `NUMERIXX_WITH_FXT` (ON) |
| [Eigen](https://eigen.tuxfamily.org) 5.0.1 | `numerixx::linalg`, `numerixx::multiroots` | `NUMERIXX_WITH_LINALG` (ON) |
| Boost.Multiprecision 1.92 (standalone) | `numerixx::multiprecision` adapter | `NUMERIXX_WITH_MULTIPRECISION` (OFF) |
| doctest 2.5.3, google/benchmark 1.9.5, Boost.Math 1.92 | tests, benchmarks, test oracles | development only |

## Using Numerixx

With CPM (declare it as `NAME Numerixx`, so that `-DCPM_Numerixx_SOURCE=<path>` can substitute a local checkout):

```cmake
CPMAddPackage(NAME Numerixx GITHUB_REPOSITORY troldal/Numerixx GIT_TAG <commit-or-tag>)
target_link_libraries(my_target PRIVATE numerixx::numerixx)   # every module, or link single ones: numerixx::roots
```

With FetchContent:

```cmake
include(FetchContent)
FetchContent_Declare(numerixx GIT_REPOSITORY https://github.com/troldal/Numerixx.git GIT_TAG <commit-or-tag>)
FetchContent_MakeAvailable(numerixx)
```

As an installed package:

```cmake
find_package(numerixx 2.0 CONFIG REQUIRED)
```

Each module has its own target and header: `numerixx::roots` and `<numerixx/roots.hpp>`, and so on.
`numerixx::numerixx` links every module, including the Eigen-backed `linalg` and `multiroots` when
`NUMERIXX_WITH_LINALG` is ON (the default), and `numerixx::pipes` when `NUMERIXX_WITH_FXT` is ON. The umbrella
header `<numerixx/numerixx.hpp>` includes the scalar modules (and the pipes) but never the Eigen-backed ones, so
Eigen's compile cost is paid only where `<numerixx/linalg.hpp>` or `<numerixx/multiroots.hpp>` is included. A user
of the scalar modules only can set `NUMERIXX_WITH_FXT=OFF` and `NUMERIXX_WITH_LINALG=OFF` and downloads neither
dependency.

If the parent project uses FXT or Eigen itself:

- **CPM parents:** declare FXT as `NAME FXT` (with `FXT_USE_TL_EXPECTED` and `FXT_USE_TL_OPTIONAL` OFF) and Eigen as
  `NAME Eigen3`. CPM deduplicates them in either order: declared first, Numerixx reuses the parent's `fxt::fxt` and
  `Eigen3::Eigen`; declared after Numerixx, the parent gets those targets from Numerixx.
- **FetchContent parents:** either fetch FXT and Eigen before Numerixx and define `fxt::fxt` and `Eigen3::Eigen`,
  which Numerixx then reuses, or add Numerixx first and use the `fxt::fxt` and `Eigen3::Eigen` targets it defines.
  After Numerixx, do not fetch FXT or Eigen again, define those targets again, or add FXT's or Eigen's own CMake
  project: that is either an error (the target already exists) or a second copy of the library.

## Building and testing

Every configuration is a CMake preset; `cmake --workflow --preset <name>` configures, builds and tests it.

| Preset | Toolchain |
|---|---|
| `gcc`, `gcc-noexcept`, `gcc-multiprecision` | GCC + libstdc++ with library assertions; without exceptions; with the multiprecision adapter and Boost.Math oracles |
| `clang`, `clang-asan` | Clang + libc++; with AddressSanitizer, UndefinedBehaviorSanitizer and libc++ debug hardening |
| `msvc`, `clang-cl` | MSVC and clang-cl (run from a Developer PowerShell) |
| `emscripten`, `emscripten-jsexcept`, `emscripten-noexcept`, `emscripten-pthread` | Emscripten with wasm, JavaScript or no exceptions, and with `-pthread`; tests run under node (activate emsdk first) |
| `integration` | the consumer-build scenarios: CPM and FetchContent parents, a scalar-only parent, an installed package |

On Windows, keep the CPM cache and build directories short, because long paths lose files silently. Workflow
presets accept no `-D` options, so set the cache location through the environment (or a `CMakeUserPresets.json`):

```bash
CPM_SOURCE_CACHE=C:/cpm cmake --workflow --preset gcc
```

## Licence

MIT. See [LICENSE](LICENSE).

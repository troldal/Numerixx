# Feasibility prototype (throwaway)

This directory is **evidence for the redesign plan, not library code**. It was written from an
earlier draft of the design to prove that the core abstractions compile and behave as intended on
every target toolchain. Some names in [`../DESIGN.md`](../DESIGN.md) have since changed (those
places are marked *[prototyped mechanism]* there); the mechanisms are the same.

The prototype also still contains an in-house constexpr LU (`nxx/linalg.hpp`) and zero-heap
choices that the design no longer requires: the library uses Eigen 5.0.1 behind `nxx::linalg`
(DESIGN D22), and heap allocation is allowed (DESIGN §3.6).

Everything here is header-only C++23 and uses only the standard library, except `nxx/pipes.hpp`
(FXT) and `nxx/linalg_eigen.hpp` (Eigen 5.0.1).

## Contents

| Path | What it shows |
|---|---|
| `nxx/core.hpp` | scalar traits, refined types with `consteval` literal checks and `make()`, `errc`/`fault`/`failure`/`solution`, stop criteria algebra, the single bounded-iteration driver, `first_of`/`then`/`warm_fallback`, `steps_view` |
| `nxx/roots.hpp` | bisection, Brent, secant, Newton (with derivative policy and projection), `expand` through the driver |
| `nxx/deriv.hpp` | integer stencils, relative steps, `derivative_of`, the `numeric` derivative policy |
| `nxx/linalg.hpp` | fixed-size `vec`/`mat`, constexpr pivoted LU, `vector_traits` |
| `nxx/linalg_eigen.hpp` | Eigen-backed `lu_solve` returning `std::expected` (checks dimensions, `rcond`, finiteness) |
| `nxx/multiroots.hpp` | damped Newton with Armijo backtracking, box projection, best-iterate tracking |
| `nxx/pipes.hpp` | the only FXT include: `using fxt::operator|` for consumers |
| `test_core.cpp` | 1-D core: refined types, headline solver chain (also inside `static_assert`), criteria, projection, warm start, fallible callbacks, FXT pipes, sizes |
| `test_nd.cpp` | 2×2 constexpr damped Newton on in-house LU; with `-DNXX_WITH_EIGEN`, the same solver on `Eigen::Vector2d`/`VectorXd` |
| `neg/` | 16 compile-fail tests (each must **fail** to compile, ideally with a readable reason) and 4 probes |
| `base.cpp`, `eigen_only.cpp` | compile-time baselines for `time_*.{sh,ps1}` |
| `probe_features.cpp` | prints the C++23 feature-test macros of the current toolchain |
| `fxt-1.patch` | the 2+2-line FXT fix (`throw 0;` → `std::unreachable();`) needed for `-fno-exceptions` builds |
| `nxx/runtime.hpp`, `test_runtime.cpp` | run-time chains: `any_solver<F, Est, UE>` (copyable, never-empty `std::function` wrapper for a curried solver) and `first_of` over a run-time range (same semantics as the static one); the test builds a chain from strings, checks it equals the static chain bit for bit, and counts heap allocations. No FXT (`./build_gnu.sh test_runtime.cpp rt`, `./build_msvc.ps1 test_runtime.cpp rt`) |

## Building

The scripts assume the toolchains of the development machine (`C:/Toolchains/GCC16`,
`C:/Toolchains/LLVM22`, `C:/Toolchains/emsdk`, Visual Studio 18 for MSVC/clang-cl) and FXT at
`C:/Dev/XLThermo/FXT`. Adjust the paths at the top of each script if yours differ.

Environment variables:

- `EIGEN_DIR`: an extracted Eigen 5.0.1 source tree (default `../eigen-5.0.1`). Only needed for
  `test_nd.cpp -DNXX_WITH_EIGEN`.
- `FXT_PATCHED_DIR`: the `include/` directory of an FXT checkout with `fxt-1.patch` applied
  (default `fxt_patched/include`). Only needed for the `-fno-exceptions` legs.

```bash
./build_gnu.sh test_core.cpp core -DNXX_FIX_UNWRAP_FAILURE     # GCC, Clang+libc++, em++ in 3 EH modes
./build_gnu.sh test_nd.cpp nd -DNXX_WITH_EIGEN
./neg_gnu.sh                                                   # compile-fail suite on GCC and Clang
```

```powershell
./build_msvc.ps1 test_core.cpp core NXX_FIX_UNWRAP_FAILURE      # MSVC 19.51 and clang-cl 22
./neg_msvc.ps1
```

Quick single-compiler check:

```bash
g++ -std=c++23 -O2 -DNXX_FIX_UNWRAP_FAILURE -I. -IC:/Dev/XLThermo/FXT/include test_core.cpp -o out/core && out/core
```

Both `test_core.cpp` and `test_nd.cpp` end with `PASS (0 check failures)`. The build output goes
to `out/`.

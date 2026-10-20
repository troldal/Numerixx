# Numerixx

Numerixx 2 is a header-only C++23 numerical library (MIT) for programs that need numerical algorithms on built-in or
user-defined real types. It is inspired by the composability of GSL: small algorithms and building blocks that combine,
as in the STL. Numerixx 1.x is at the tag `v1.1.0-legacy`.

## Goal and success criteria

Current milestone: numerical derivatives of a scalar function (`nxx::deriv`). Its phases and checkboxes are in
`docs/plan.md`; a phase is done when its `static_assert` and run-time tests pass on every toolchain below, locally and
in CI.

## Layout

- `include/numerixx/` public headers; the library is header-only, so there is no `src/`
- `tests/` doctest, with standalone Boost.Multiprecision for the multiprecision cases
- `demo/` small demo programs (`NUMERIXX_BUILD_DEMOS`, on when top level), e.g. `build/<preset>/demo/numerixx_deriv_demo`
- `cmake/CPM.cmake` fetches the test dependencies; the library itself has none
- `CMakePresets.json` one preset per toolchain, used locally and in CI (`.github/workflows/ci.yml`)
- `docs/plan.md` current plan; `docs/design.md` design decisions; also a technique survey, a tutorial and API
  proposals (`docs/PHASE2_API.md`)

## Build and test

Run from the project root, with `CPM_SOURCE_CACHE=C:/cpm` in the environment (a long path loses files to MAX_PATH).
Presets: `clang-cl`, `cl`, `gcc`, `clang`, `emscripten`; each builds Debug into `build/<preset>`.
- Configure: `cmake --preset <preset>`
- Build: `cmake --build --preset <preset>`
- Run affected tests: `build/<preset>/tests/numerixx_tests --test-case="deriv:*"` (one doctest executable)
- Full suite: `ctest --preset <preset>`
- `clang-cl` and `cl` need the MSVC environment, which the agent's shell lacks: run them from a batch file that first
  calls `"C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat"`.
- `g++` on PATH is Strawberry Perl's; configure `gcc` with `-DCMAKE_CXX_COMPILER=C:/Toolchains/GCC16/bin/g++.exe` and
  put `C:/Toolchains/GCC16/bin` first on PATH.
- `emscripten` needs `EMSDK` set: source `C:/Toolchains/emsdk/emsdk_env.sh` first. The tests run under node.
- WSL (openSUSE Tumbleweed) has GCC 16 and Clang 23. Build there from a copy of the tracked and new files (a tar of
  `git ls-files` and `git ls-files --others --exclude-standard`), not from `/mnt/c`: `build/<preset>` would clash
  with the Windows builds, and the untracked CLion build folders are gigabytes.
- Primary toolchain: clang-cl + MSVC STL, x64. Also must build with MSVC `cl`, GCC and Clang (MinGW and WSL), and
  Emscripten.

## Project-specific constraints

- **Workflow.** Work happens on `dev`. Claude commits and pushes to `dev` after each step the user approves; a PR to
  `master` is opened only when the user asks. One small step at a time, and the user decides the next one. Since the
  library is header-only, every step changes the public API: propose its API before writing it. Before pushing, build
  and test on all five toolchains; after pushing, check CI (`gh run list --workflow ci.yml --branch dev`).
- **Licensing.** Code may come only from Numerixx 1.x and Boost.Math (BSL-1.0; keep its notice). Never use GSL or
  other GPL, LGPL or AGPL code, Numerical Recipes listings, or code without a licence. Implement algorithms from the
  papers and textbooks, and do not read GSL source.
- **Errors.** Nothing throws, not even for exceptional conditions: run-time failures return `std::expected`.
- **User types.** Algorithms take any `nxx::real`, such as Boost.Multiprecision types. Algorithm code uses only
  `+ - * /`, comparisons and the `nxx::math` helpers; write `T(1)`, never a `double` literal; never `auto` on an
  arithmetic expression (Boost's expression templates).
- **Emscripten.** Everything must build there: no threads, stack traces or other platform facilities.
- **Callbacks** are named `fn`. Qualify internal calls (`nxx::deriv::...`) so user types cannot hijack them by ADL.
- **Tests.** Reference values come from closed forms or independent high-precision sources, never from Numerixx
  itself; each tolerance is justified by an error estimate in a comment.
- **Generality.** Justify features on general numerical grounds; do not name downstream projects in code, docs or
  commits.

## Lessons learned

- None yet.

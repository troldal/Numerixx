# Numerixx

Numerixx 2 is a header-only C++23 numerical library (MIT) for programs that need numerical algorithms on built-in or
user-defined real types. It is inspired by the composability of GSL: small algorithms and building blocks that combine,
as in the STL. Numerixx 1.x is at the tag `v1.1.0-legacy`.

## Goal and success criteria

Current milestone: numerical derivatives of a scalar function (`nxx::deriv`). An item is done when its `static_assert`
and run-time tests pass on every toolchain below, locally and in CI.
- [x] Central difference with a default step, on `double` and Boost.Multiprecision (29ba464, CI run 39437081443)
- [ ] Higher-order central stencils, and forward and backward stencils
- [ ] Ridders' extrapolation with an error estimate
- [ ] Second derivatives

## Layout

- `include/numerixx/` public headers; the library is header-only, so there is no `src/`
- `tests/` doctest, with standalone Boost.Multiprecision for the multiprecision cases
- `cmake/CPM.cmake` fetches the test dependencies; the library itself has none
- `.github/workflows/ci.yml` CI on GCC, Clang with libc++, MSVC, clang-cl and Emscripten

## Build and test

Run from the project root. There are no presets yet; `<tc>` names the build directory of one toolchain. The agent's
shell lacks the MSVC environment, so wrap the MSVC and clang-cl commands:
`cmd /c '"C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat" >nul && <command>'`.
- Configure: `cmake -S . -B build/<tc> -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER=<compiler> -DCPM_SOURCE_CACHE=C:/cpm`
- Build: `cmake --build build/<tc>`
- Run affected tests: `build/<tc>/tests/numerixx_tests --test-case="deriv:*"`
- Full suite: `ctest --test-dir build/<tc> --output-on-failure`
- Emscripten: run `emsdk_env` first, then `emcmake cmake ...`; the tests run under node.
- Primary toolchain: clang-cl + MSVC STL, x64. Also must build with MSVC `cl`, GCC 16 and Clang 22 (MinGW and WSL),
  and Emscripten.

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

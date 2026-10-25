# Numerixx

Numerixx 2 is a header-only C++23 numerical library (MIT) for programs that need numerical algorithms on built-in or
user-defined real types. It is inspired by the composability of GSL: small algorithms and building blocks that combine,
as in the STL. Numerixx 1.x is at the tag `v1.1.0-legacy`.

## Modules and milestones

One module per namespace and header. `docs/plan.md` holds a milestone per module with its phases and checkboxes, and
the deferred features; a phase is done when its `static_assert` and run-time tests pass on every toolchain below,
locally and in CI.

| Module | Header | Documents | Test cases | Status |
|---|---|---|---|---|
| `nxx::deriv`: derivatives of a scalar function | `deriv.hpp` | `docs/deriv/` | `deriv:*`, `orders:*`, `ridders:*`, `onesided:*` | Milestone 1, done |
| `nxx::roots`: roots of a real function of one variable | `roots.hpp` | `docs/roots/` | `roots:*` | Milestone 2, current |

A module's folder holds its design record (`design.md`: the decisions in force for the module, its primary sources
and its test rules), its phase specifications (`PHASE<n>_API.md`) and its briefs and surveys. Read the module's
`design.md` and the current phase document before working on it; `docs/design.md` holds what applies to every module.
A module may use another module's public API, never its `detail`.

## Layout

- `include/numerixx/` public headers: one per module, plus `math.hpp` for what the modules share. Header-only, no
  `src/`
- `tests/` one doctest executable (`numerixx_tests`), with standalone Boost.Multiprecision for the multiprecision
  cases; a new module's test files and test cases carry its prefix (`test_roots_*.cpp`, `roots:*`)
- `demo/` one demo program per module (`NUMERIXX_BUILD_DEMOS`, on when top level), e.g.
  `build/<preset>/demo/numerixx_deriv_demo`
- `cmake/CPM.cmake` fetches the test dependencies; the library itself has none
- `CMakePresets.json` one preset per toolchain, used locally and in CI (`.github/workflows/ci.yml`)
- `docs/plan.md` the plan; `docs/design.md` the decisions for every module; `docs/<module>/` per module (above)
- `docs/Doxyfile` the API documentation settings, used by the docs targets in `docs/CMakeLists.txt`
- `docs/deriv-book/` the book on `nxx::deriv` (LaTeX): its own `CLAUDE.md` holds its rules, roles and process; its
  examples build as tests (`NUMERIXX_BUILD_BOOK`), its PDF is the target `numerixx_book`

## Build and test

Run from the project root, with `CPM_SOURCE_CACHE=C:/cpm` in the environment (a long path loses files to MAX_PATH).
Presets: `clang-cl`, `cl`, `gcc`, `clang`, `emscripten`; each builds Debug into `build/<preset>`.
- Configure: `cmake --preset <preset>`
- Build: `cmake --build --preset <preset>`
- Run affected tests: `build/<preset>/tests/numerixx_tests --test-case="roots:*"` (the module's prefixes, above)
- Full suite: `ctest --preset <preset>`
- Docs: `cmake --build --preset <preset> --target numerixx_docs` (or `numerixx_docs_internal`, which adds the `detail`
  namespaces) writes `build/<preset>/docs/public/html`; a Doxygen warning fails it. The targets exist when Doxygen
  1.10 or later is found (Strawberry Perl's 1.10 is first on PATH; 1.18 is in `C:/Program Files/doxygen/bin`).
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
  library is header-only, every step changes the public API: propose its API before writing it, as
  `docs/<module>/PHASE<n>_API.md`, reviewed by `numerics-reviewer`, `simplicity-critic` and `usage-sketcher` before it
  goes to the user. Before pushing, build and test on all five toolchains; after pushing, check CI
  (`gh run list --workflow ci.yml --branch dev`).
- **Independent review at the end of each phase.** After the phase's implementation is committed and CI passes, run a
  second-opinion review by another model with the OpenAI Codex CLI, read-only, on the phase's commit:
  `codex exec -s read-only --ephemeral --color never -o <scratchpad>/codex_review.md "<prompt>"` from the project
  root, with a prompt that names the commit, the spec (`docs/<module>/PHASE<n>_API.md`), `docs/design.md`, the
  module's `design.md` and the rules here, and asks for ranked findings, each with a file, line, counterexample or
  derivation, plus what was found correct. Reproduce every finding before acting on it, propose the fixes to
  `include/` to the user before making them, and record the adopted ones in the phase's API document. (deriv Phase 4:
  three edge cases our own reviews had missed.)
- **Licensing.** Code may come only from Numerixx 1.x and Boost.Math (BSL-1.0; keep its notice). Never use GSL or
  other GPL, LGPL or AGPL code, Numerical Recipes listings, or code without a licence. Implement algorithms from the
  papers and textbooks, and do not read GSL source.
- **Errors.** Nothing throws, not even for exceptional conditions: run-time failures return `std::expected`.
- **User types.** Algorithms take any `nxx::real`, such as Boost.Multiprecision types. Algorithm code uses only
  `+ - * /`, comparisons and the `nxx::math` helpers; write `T(1)`, never a `double` literal; never `auto` on an
  arithmetic expression (Boost's expression templates).
- **Emscripten.** Everything must build there: no threads, stack traces or other platform facilities.
- **Naming.** Functions, concepts and enumerations follow the standard library (`derivative`, `isfinite`, `real`,
  `error::non_finite`). Classes and class templates, meaning types that can have member functions, are PascalCase
  (`Central<P>`, `Estimate<T>`). No digits in a class name where a template parameter belongs (`Central<4>`, not
  `Central4`).
- **Callbacks** are named `fn`. Qualify internal calls (`nxx::math::...`, `nxx::deriv::...`, `nxx::roots::...`) so
  user types cannot hijack them by ADL.
- **Tests.** Reference values come from closed forms or independent high-precision sources, never from Numerixx
  itself; each tolerance is justified by an error estimate in a comment. One exception: a bound on the quality of a
  result for which no independent value exists, such as Ridders' error estimate in a multiprecision type whose guard
  digits no prototype models, may be a regression baseline from a recorded run of Numerixx. It is labelled as such in
  the test, with the commit it came from, and never replaces the correctness check against a closed form.
- **Generality.** Justify features on general numerical grounds; do not name downstream projects in code, docs or
  commits.

## Lessons learned

- None yet.

# Numerixx: instructions for Claude Code

Numerixx 2 is a header-only C++23 numerical library, MIT-licensed, general-purpose in the spirit of GSL with a
smaller scope. Nothing throws: every result is a `std::expected`. The library is being rewritten (Numerixx 1.x is kept
at the tags `v1.0.0` and `v1.1.0-legacy`).

- `docs/redesign/PLAN.md`: the short plan. Its **Status** line says where the work is.
- `docs/redesign/DESIGN.md`: the reference. Its § numbers are stable; cite them in code comments, commits and docs.
  §10.3 is the roadmap: phases, scope, acceptance criteria and what is left of each. §2 and §12 list the decisions.
- Machine-specific setup (toolchain paths, environments) belongs in `CLAUDE.local.md`, which is not committed.

## Rules

1. **Stay inside the current roadmap phase** (DESIGN §10.3). Build only what the phase's scope and criteria name.
   When a task seems to need a later phase's code, or the plan is unclear on scope, ask before building.
2. **The DESIGN §2 and §12 decisions are settled.** Do not reopen them.
3. **No GPL code.** Never port, paraphrase or copy code under the GPL, LGPL or AGPL (GSL and MPSolve among them),
   including code found on the web. Port from Numerixx 1.x, or implement from the papers.
   - The default 1.x source is `v1.1.0-legacy` (dev-reorg). Use `v1.0.0` (master) where DESIGN §7 or §10.1 says
     dev-reorg regressed. Read files with `git show <tag>:<path>`.
   - Cite the source in the header and in the commit message, and keep the BSL-1.0 notice on Boost-derived code.
   - Each "must not port" bug in DESIGN §7 gets a regression test, in the PR that ports its module.
4. **Justify everything on general numerical grounds.** Do not name downstream or consumer projects in code, docs,
   tests or commit messages. Use general examples: box constraints, expensive functions, fallible callbacks, and
   published test suites (Alefeld–Potra–Shi, Moré–Garbow–Hillstrom).
5. **No silent failure.** A success with `stop_reason::criterion` must meet that criterion's guarantee (DESIGN §9.3).
   A failure is a value with an error code, the cost and the best estimate. Never report an endpoint, a pole or the
   last iterate as a success.
6. **Claims must be verified.** Numbers in docs, commits and PR text are measured, not estimated. Before writing
   "tested", "every" or "all", check it. When local and hosted results differ, the hosted CI result counts.
7. **Check hosted CI before reporting work as done.**
   - `ci.yml` runs on pull requests, on pushes to `master` and on manual dispatch. A plain branch push starts
     nothing, so the branch needs an open PR.
   - Use the run whose `headSha` is `git rev-parse HEAD`:
     `gh run list --workflow ci.yml --branch <branch> --json headSha,status,conclusion,databaseId`.
   - Hosted compilers differ from local ones: the floating `gcc:16` container (GCC 16.2 at the time of writing),
     Clang 22 with libc++ on Ubuntu, MSVC 19.51 with its bundled clang-cl 20 on `windows-2025-vs2026`, and the emsdk
     version pinned in `ci.yml`. A job log's "The CXX compiler identification is" line gives the exact version.
8. **Git.**
   - Work on a branch and open PRs against `master`. The user merges PRs and deletes branches. Do not merge,
     force-push, delete branches or rewrite published history.
   - Only the main session commits and pushes, and only after the preset runs below. Subagents leave their changes
     uncommitted.

## Build and test

Everything goes through the presets in `CMakePresets.json`; build trees are `build/<preset>`.

| Preset | What it checks |
|---|---|
| `gcc` | GCC + libstdc++ with `_GLIBCXX_ASSERTIONS`, benchmarks built |
| `gcc-noexcept` | GCC with `-fno-exceptions` |
| `gcc-multiprecision` | GCC with Boost.Multiprecision (`cpp_bin_float_50`) and the Boost.Math oracles |
| `clang` | Clang + libc++ |
| `clang-asan` | Clang + libc++, ASan + UBSan + libc++ debug hardening |
| `msvc`, `clang-cl` | Windows; configure from a Developer PowerShell (or after `vcvars64.bat`) |
| `emscripten`, `emscripten-jsexcept`, `emscripten-noexcept`, `emscripten-pthread` | em++ with wasm, JavaScript-based and no exceptions, and wasm exceptions with `-pthread`; tests run under node; need `EMSDK` |
| `integration` | the consumer-build scenarios only (CPM and FetchContent parents, scalar-only, installed package); builds no unit tests |

```bash
cmake --workflow --preset gcc --fresh      # configure, build, test from scratch
cmake --build --preset gcc                 # incremental build
ctest --preset gcc -R "^roots\."           # one module's cases
ctest --preset gcc -L compile-fail         # labels: compile-fail, compile-time, example, structural, <module>
ctest --preset gcc -I 57,57 --output-on-failure   # one test by its number (names contain regex characters)
build/gcc/tests/numerixx_test_roots -tc="solvers: brent*"   # a doctest binary directly
```

- **Before a PR, run all 12 presets from a fresh configure.** Clang-only or GCC-only greens are not enough: MSVC
  overload resolution, em++ and the no-exceptions modes each caught real bugs.
- **Read the build output before trusting ctest.** If the build fails, ctest runs the previous binaries and reports
  a stale pass.
- **"could not load cache"** from every compile-fail and compile-time test means `build/<preset>/CMakeCache.txt` is
  gone: reconfigure with `cmake --preset <preset>`.
- **Do not run em++ while an Emscripten preset is building.** A different emsdk configuration clears the shared cache
  and breaks the running build.
- **Compiler floor** (DESIGN D2): GCC 14, and Clang 18 + libc++ 18.
  - Neither the 12 presets nor branch CI check it; they use GCC 16 and Clang 22. Only `.github/workflows/nightly.yml`
    does, together with MinGW g++, Intel ICX and clang-cl without exceptions. It runs on a schedule, on `master`.
  - Check `gh run list --workflow nightly.yml --limit 3` when a phase starts and before opening a PR.
  - When a change adds a language or library feature, say so in the PR. Ask the user before running
    `gh workflow run nightly.yml --ref <branch>`.
- **Format** with clang-format 22 (`.clang-format`). CI runs `clang-format-22 --dry-run -Werror` over the `.hpp` and
  `.cpp` files in `include/`, `tests/`, `examples/` and `benchmarks/`.
- Tests and examples compile with `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion ... -Werror`
  (`/W4 /WX` on Windows). A library header may never add a warning to a consumer's build, including a compiler's
  false positive.

## Code conventions (DESIGN §3, §5, §6)

- **Layout:** `include/numerixx/<module>.hpp` includes `include/numerixx/<module>/*.hpp`. Each module is one
  INTERFACE target (`numerixx::<module>`), and the module graph in DESIGN §5.2 has no cycles.
- **Every header with arithmetic** wraps its body in `NXX_BEGIN_HEADER` … `NXX_END_HEADER` (`config.hpp`).
  - On Clang, clang-cl and em++ that turns floating-point contraction off, so that, given the same values of f, a
    solver takes the same path on every preset. The golden table in `tests/roots/test_determinism.cpp` checks this.
  - GCC still contracts on FMA targets. Whether to add `-ffp-contract=off` is an open phase-1 question (DESIGN
    §5.3); do not decide it without the user.
- **Maths:**
  - Call the `nxx::math` helpers (`abs`, `isfinite`, `isnan`, `sqrt`, `midpoint`, `pow2`, `root_eps`), not `<cmath>`,
    on constexpr and multiprecision paths.
  - Stop tests, step rules and midpoints use only `+ - * /` and the exact or correctly rounded helpers (DESIGN §6.1).
  - No `auto` on arithmetic expressions of Eigen or multiprecision types, and no Eigen expression returned through a
    deduced type.
- **Reasoned deletions:** `NXX_DELETE("reason")`, starting on the declarator's own line (cl reports only that line).
  Invalid calls must make `std::is_invocable_v` false, not a hard error: constrain with `bool` variable templates and
  add a deleted sibling that carries the reason.
- **Qualify internal calls** (`nxx::evaluate(...)`, `nxx::detail::...`) so user types cannot hijack them through ADL.
  - Write `(std::min)`, `(std::max)` and `(std::numeric_limits<T>::max)()`, so that `<windows.h>` macros cannot break
    the headers.
  - Name callback parameters `fn`, never `f`.
  - Mark parameters that only some `if constexpr` branches use `[[maybe_unused]]`.
- **Refined types** validate literals in consteval constructors; run-time values go through `make()` and return a
  `std::expected`. The library validates run-time inputs in-band and never throws. Keep solvers constexpr where
  portable.
- **Iterative methods** follow the protocol in DESIGN §6.6: `prepare`, `init`, `step`, `estimate`, `best`, `view`,
  `intrinsic`, and optionally `finish`. They derive from a family facade and keep their configuration in
  `options<Stop, Deriv, Proj, Obs>`. Evaluation counts must equal the calls of the user's function.
- **Determinism:** when a solver's path changes on purpose, regenerate the golden table in
  `tests/roots/test_determinism.cpp`. Compile it (with `tests/doctest_main.cpp`) with `-DNXX_PRINT_GOLDEN`, run the
  case `determinism: golden values`, and paste the printed rows. Then check that the rows pass on every preset that
  builds the tests. A row that differs on one toolchain is a determinism bug; do not regenerate it there.

## Tests (DESIGN §9)

- doctest, one executable per module. The CTest name of a case is `<name>.<case name>`, where `<name>` is the first
  argument of `numerixx_add_test` in `tests/CMakeLists.txt`.
- **Sources are listed, not globbed:** add a new test file to the `SOURCES` of its `numerixx_add_test(...)` call.
- **Never guard a dereference with `REQUIRE`:** under `-fno-exceptions` it reports but does not stop. Compare whole
  results (`CHECK(r == expected)`) or dereference inside `if (r)`.
- Property tests draw from a `std::mt19937` with a fixed seed.
- **Reference values** (DESIGN §9.2):
  - at least 21 significant digits (40 for multiprecision), generated once at high precision and committed with the
    command that generated them;
  - tolerances go through `tol<T>(k, ref_eps)`;
  - Numerixx output is never its own reference.
- **Compile-fail cases** live in `tests/compile_fail/<case>.cpp`: the `#ifdef NUMERIXX_CF_CONTROL` branch must compile
  and the other must not.
  - Register each case in `tests/compile_fail/compile_fail.cmake` with an `EXPECT` regex. On GCC and Clang
    (clang-cl too) it must match the first error; quoted source lines count, and the consteval literal cases rely on
    that.
  - Add `DELETE_REASON` for an `NXX_DELETE` reason. The reason must then appear in the compiler's own message, and
    only GCC 15+ and Clang 19+ check it. Older compilers and cl only check that the build fails.
  - Line counts go to `build/<preset>/compile_fail_report.txt`.
- Each module phase adds its canonical calls (DESIGN §6.14) to `tests/usage/canonical_calls.cpp`.
- `examples/` are built and run as smoke tests; `examples/quick_tour.cpp` shows the current API.

## Docs and changes

- Every user-visible change goes into `CHANGELOG.md`.
- `MIGRATION.md` maps 1.x to Numerixx 2 (DESIGN §10.4), and updating it is part of every phase's acceptance
  criteria. When a change replaces a 1.x API, or changes a result because 1.x was wrong, add or correct its row.
- Keep DESIGN in step with the code: the status marks (`[sketch]`, `[prototyped]`, `[spike]`), the §10.2 and §10.3
  status, and Appendix D, whose numbers come from the build reports and from measurements.
- Commit messages: an imperative summary line, then why the change was made.

## Subagents (`.claude/agents/`)

| Agent | Use it to |
|---|---|
| `architect` | design a new family, or answer a core question, as a note with 2-3 options for the user to approve |
| `api-ergonomics-reviewer` | review an architect note from the caller's side: user code, cross-family consistency, misuse |
| `phase-scope-checker` | check a plan, design note or diff against the current roadmap phase |
| `algorithm-implementer` | implement or port one algorithm, with tests and docs |
| `test-author` | write doctest cases, property tests, corpus cases, compile-fail cases and canonical calls |
| `numerics-reviewer` | review numerical correctness adversarially, with probes |
| `cpp-reviewer` | review C++ mechanics and portability across GCC, Clang, MSVC, clang-cl and em++ |
| `matrix-runner` | build and test presets from scratch and report the results |
| `ci-investigator` | find out why hosted CI failed, and reproduce it |
| `docs-auditor` | check that the docs' claims match the code and the measurements |

**A new family** (optimize, multiroots with the linalg facade, integrate, interpolate, poly) is designed before it is
built:
1. `phase-scope-checker` settles what the phase builds and what it only accommodates.
2. `architect` writes the note.
3. The note is reviewed in parallel: always by `api-ergonomics-reviewer` and by `phase-scope-checker` (its
   `build now` and `accommodate, do not build` marks), and by `cpp-reviewer` and `numerics-reviewer` in design-note
   mode when the options differ in their area.
4. The main session calls `architect` again with its note and every reviewer's findings, verbatim. It revises once,
   or records each disagreement.
5. The user approves an option and settles every `needs a decision` item. The main session writes the approved option
   into the family's DESIGN §7 section: its types, canonical calls, misuse catalogue, and the 1.x entry points that
   `MIGRATION.md` must map.
6. `algorithm-implementer` builds the first member from that section, including every rejection in the misuse
   catalogue. `test-author` then turns the catalogue into compile-fail and doctest cases and adds the canonical calls.
7. The usual reviews, `docs-auditor` (for the `MIGRATION.md` rows), `matrix-runner` and the commit follow.

A core note (a change to `core/` types, error codes or the criteria algebra) follows steps 1-5, and its result goes
into the DESIGN §6 section it changes. In step 3, `api-ergonomics-reviewer` reviews it only when it changes a
user-visible spelling, result field or error code.

`matrix-runner`, `algorithm-implementer`, `test-author` and `ci-investigator` build in `build/<preset>`. A `--fresh`
configure deletes the cache under any build running in the same tree, so run only one of them at a time. The other
agents never build in or write to `build/`: the architect compiles nothing, the reviewers and `docs-auditor` compile
probes only in scratch directories outside the repository, and `docs-auditor` may read the build reports and run
`ctest -N`.

# Milestone 3, Phase 2: the documentation skeleton

Status: approved by the user (2026-10-27) after the reviews (`simplicity-critic`, `numerics-reviewer`; the record is
at the end), with the comment fix to `deriv.hpp:14`, and with uv instead of a hand-made virtual environment.
Milestone 3 is in `docs/plan.md`; Phase 1 (the spike) chose Sphinx, Breathe and MyST, with a patch to Breathe for
`requires` clauses and concepts.

## Scope

The machinery and the pages that need no prose of their own: the build and its checks, CI, a getting-started page,
the reference page of each namespace, and tests for the claims of the headers' `@code` blocks. The user guides of
`nxx::roots` and `nxx::deriv` are Phases 3 and 4. `include/` changes only if the user approves the comment fix below.

## Layout

```
docs/sphinx/
  conf.py              Sphinx settings, the Breathe patch and the checks
  requirements.txt     pip freeze of the venv: sphinx 9.2.0, breathe 5.1.0, myst-parser 5.1.0, furo 2025.12.19 and
                       every transitive package (docutils 0.22.4, Pygments 2.21.0, ...)
  index.md             what the library is, the modules, links
  getting-started.md   the first program (no install section until the first release: there is no tag to pin)
  reference/roots.md   doxygennamespace nxx::roots, doxygenstruct std::formatter< nxx::roots::error >
  reference/deriv.md   doxygennamespace nxx::deriv, doxygenstruct std::formatter< nxx::deriv::error >
  reference/core.md    doxygenconcept nxx::real, doxygennamespace nxx::math
  examples/getting_started.cpp
  PHASE2_PLAN.md       this document, excluded from the site
```

## Build

`docs/CMakeLists.txt` keeps its two targets; `numerixx_add_docs` only writes the small Doxyfile both need.

- `numerixx_docs`: the public site. Doxygen writes XML only (`GENERATE_XML = YES`, `GENERATE_HTML = NO`,
  `EXCLUDE_SYMBOLS = *::detail *::detail::*`, `GENERATE_DEPRECATEDLIST = NO`), then `python -m sphinx -E -W -b html`,
  writing `build/<preset>/docs/public/html`. The settings go in the generated Doxyfile, so `docs/Doxyfile` is
  unchanged.
  - `EXCLUDE_SYMBOLS` because `deriv.hpp` has `detail` blocks outside `@cond INTERNAL`, which the current public
    Doxygen HTML shows (about 40 entries); this fixes it without touching the header.
  - `GENERATE_DEPRECATEDLIST = NO` because the only `[[deprecated]]` declarations are the low-precision warnings of
    `derivative` (`detail::noisy`), not deprecations; the list made empty "Deprecated" sections and broken links.
  - `-E` so that a local rebuild re-reads every page and `-W` sees every warning (the site is small); `-n` is
    `nitpicky = True` in `conf.py`; `--keep-going` is the default in Sphinx 9.2.
- `numerixx_docs_internal`: unchanged, Doxygen HTML with the `detail` namespaces, for developers.

`numerixx_docs` exists when uv is found and Doxygen is 1.18 or later, the version whose XML the patch was verified
on; otherwise a status message and only the internal target. The target runs Sphinx as `uv run --no-project
--python 3.14 --with-requirements docs/sphinx/requirements.txt python -m sphinx ...`: uv builds and caches the pinned
environment itself, so there is no virtual environment to make or select, locally or in CI. `conf.py` first checks
that Breathe is 5.1.0 and stops with a message naming `requirements.txt` otherwise, so a run outside that
environment fails rather than rendering without the patch.

## conf.py

- `breathe_projects` points at the XML that CMake passes in (`-D breathe_projects.Numerixx=<dir>`), so `conf.py`
  holds no build path.
- MyST, Furo, `primary_domain = "cpp"`, `highlight_language = "cpp"`, `breathe_show_include = False`,
  `nitpicky = True`, `exclude_patterns = ["PHASE*.md"]`, and no `suppress_warnings`. No math extension until a page
  has a formula.
- `nitpick_ignore_regex`: names in `detail` and `std`, and `formatter<string_view>` and `string_view` (the base of
  the formatters, written without `std::` by Doxygen), which have no target.
- The patch, about 30 lines, each part with a comment that names breathe#482 and the function it relies on:
  1. `SphinxRenderer.create_template_prefix` appends ` requires <clause>` from `<requiresclause>`; Sphinx accepts the
     clause after the template parameter list since 5.2.
  2. The `Node_compounddefType` entry of `SphinxRenderer.node_handlers` puts a concept's `<initializer>`, its
     definition, as a C++ code block before its description.
  Each part counts what it emits.
- The checks, in a `build-finished` handler that fails the build on a mismatch:
  1. Clauses and concepts: the clauses part 1 appended equal the `<requiresclause>` elements in the XML, and the
     definitions part 2 emitted equal the concepts in the XML, so that a Breathe change that bypasses either part
     fails the build rather than dropping them silently.
  2. Coverage: every symbol in the XML (`EXCLUDE_SYMBOLS` has removed `detail`; the `std::formatter`
     specialisations included) is reached by a directive on some page: a `doxygennamespace` covers its namespace's
     members, the other directives one entity each. From the spike's `coverage.py`, simplified to the directive forms
     the pages use.
- A Python warnings filter for Breathe's `unexpected element "conceptparts"` (a Doxygen 1.18 element Breathe 5.1 does
  not know; it reads the rest of the concept), for a clean log.

## Examples

A page shows code only by `literalinclude` from `docs/sphinx/examples/*.cpp`, between `// <name>` and `// </name>`
markers, as the book does. Phase 2 has one, `getting_started.cpp`: a complete program with Brent on x^2 - 2 over
[1, 2] and the derivative of exp at 1 by Ridders. It is one `add_executable` and `add_test`
(`numerixx_docs_example_getting_started`) in `tests/CMakeLists.txt`, with the tests' warning flags (and `stdc++exp`
on MinGW GCC for `std::print`), so it builds on all five toolchains and in every CI job, Doxygen or not. `main`
returns non-zero unless the page's claims hold, by the tests' rules: Brent's value against sqrt(2) within
`abs_error` + 8e-17 (fn's rounding over |f'|), Ridders' against e with the bounds of `deriv.hpp:12` below, each
derivation in a comment. A folder with its own CMakeLists.txt waits for the second example (Phase 3).

## The headers' examples

The 16 `@code` blocks in the headers are fragments (`use`, `fn`, `psat` are not defined), so they cannot be compiled
as they are. Each claim a block makes is asserted by a test in the module's own test files (`roots:*`, `deriv:*`),
with a comment naming the block, so that the affected-tests command of `CLAUDE.md` catches a stale comment. Claims
already asserted get that comment; the blocks with no checkable claim (`roots.hpp:142`, `1922`, `2033`) are compiled
as written, with the free names defined. Double only: every claim is about double. Each reference is a closed form or
an existing independent source, and each tolerance is the bound an existing test already derives for the same case.

| Block | Claim | Checked against | Now |
|---|---|---|---|
| `deriv.hpp:12` | e to about 1e-15, estimate about 6e-14 (wrong: below) | e; `test_ridders.cpp`'s rule on the prototype's 9.8e-16 and 6.0e-14 relative (`docs/deriv/PHASE3_API.md`) | new |
| `deriv.hpp:68` | `abs_error` present for Ridders | `CHECK(r->abs_error)` | new |
| `deriv.hpp:283` | 100 to about 1e-13 relative, after the restart | 1/0.01; `test_ridders.cpp`'s rule on the prototype's 5.0e-14 and 1.8e-13 | new |
| `deriv.hpp:310` | `Ridders` gives `non_finite`, `ForwardRidders` 2 to about 1e-13 | the error value; 2, with the bound from a run of the independent one-sided prototype at 0.5 (`judged()` of `test_onesided.cpp`) | new |
| `deriv.hpp:1654` | the default step is what `derivative(fn, x)` uses | equality with the call at `default_step` (an API identity) | new |
| `deriv.hpp:1722` | the nodes x +- h, x +- 2h, and x - 5h ... x | a recording fn, as in `test_onesided.cpp` | new |
| `deriv.hpp:1895` | `derivative_of(fn)` is f' by `Central<2>` | equality with `derivative(fn, 1.0)` | new |
| `deriv.hpp:1922` | `no estimate: non_finite` | `std::format` | `test_deriv.cpp:289` |
| `roots.hpp:22` | sqrt(2) to about 1e-16 after 9 evaluations; Newton (pair) from 1 after 7 | `brackets(r, sqrt2, 8e-17)`; the prototype counts | Brent: `test_roots_bracketing.cpp:134`; Newton in pair form new |
| `roots.hpp:86` | deduces `Bracket<double>` | a `static_assert` | new |
| `roots.hpp:1980` | Newton 7 evaluations, Secant 9 | the prototype counts of `test_roots_open.cpp:155` | Secant new |
| `roots.hpp:2270` | bracket about [231.7, 254.4] after 14 evaluations; root 250 | the exact ends 231.65625 and 254.4375 and the count 14, derived by hand from the search's steps (6 * 1.5^k from 300, all exact; the smallest \|f\| at a point is 0.42, so no rounding can change a sign); 5000/20 within c = 2.8e-14, fn's rounding over \|f'(250)\| = 0.08 | new |
| `roots.hpp:2299` | `no root: no_sign_change` | `std::format` | `test_roots_bracketing.cpp:393` |

The counts are independent values (the Python prototypes, and a hand derivation for `search`), under the existing
tests' assumption of IEEE double without FMA contraction; if a change alters one, the count is derived again or the
prototype rerun, never copied from Numerixx.

**A header comment to fix (needs approval, it is in `include/`).** `deriv.hpp:14` says Ridders gives "e to about
1e-15, with an error estimate of about 6e-14" for exp at 1, and shows them as `d->value` and `*d->abs_error`.
`abs_error` is absolute, and these are the relative figures. Measured (GCC 16, -O2): error 2.7e-15 absolute
(9.8e-16 relative), `*d->abs_error` 1.63e-13 (6.0e-14 relative). Proposed: "e to about 3e-15, with an error estimate
of about 2e-13". The test asserts the corrected claim.

## CI

A `docs` job on `ubuntu-26.04`: Doxygen 1.18.0 from its release tarball (`doxygen-1.18.0.linux.bin.tar.gz`, SHA-256
pinned), `ninja-build`, and uv from `astral-sh/setup-uv` (the action pinned by commit, as `msvc-dev-cmd` is), then
configure with tests, demos and the book off and build `numerixx_docs`, which installs Python 3.14 and the pinned
packages through uv. The example and the claims tests build in the existing five jobs.

## Changes to `docs/plan.md` and `CLAUDE.md`

- Plan, Phase 2: the target is `numerixx_docs`, not `numerixx_sphinx`; `-W` without `--keep-going`; the checks run in
  `conf.py`; standalone example programs as tests (the book's pattern: a page shows a whole program a reader can
  copy) instead of examples in `numerixx_tests`; the claims tests in the module test files; linkcheck deferred to
  hosting (below); the program that writes tables moves to Phase 3, where the first table appears.
- `CLAUDE.md`, Build and test: `numerixx_docs` builds the Sphinx site and needs uv and Doxygen 1.18. Layout:
  `docs/sphinx/`.

## Success criteria

- [ ] `numerixx_docs` builds with no warning from Doxygen or Sphinx, checks included, on the `clang-cl` preset
- [ ] Each check fails when it should: a directive removed from a page, either patch part disabled (shown once, by
  hand, and recorded here)
- [ ] `numerixx_docs_example_getting_started` and the new `roots:*` and `deriv:*` cases pass on all five toolchains,
  locally and in CI; the `docs` CI job passes
- [ ] The rendered pages read by eye: signatures with their clauses, `nxx::real` with its definition, the formatters,
  no `detail`
- [ ] The phase's reviews and the independent review (Codex), with the adopted findings recorded here

## Results

- The site builds with no warning from Doxygen 1.18.0 or Sphinx 9.2.0 on `clang-cl`; the patch rendered 13 requires
  clauses and 1 concept definition, as many as the XML holds.
- The checks fail when they should (2026-10-27, on copies of `docs/sphinx`):
  - the formatter directive removed from the roots page: `struct std::formatter< nxx::roots::error > is on no page`;
  - part 1 of the patch disabled: `0 requires clauses rendered, 13 in the XML`;
  - part 2 disabled: `0 concept definitions rendered, 1 in the XML`;
  - the math page naming a wrong namespace: each function of `nxx::math` is on no page.
- First draft of part 2 rendered the concept's text through Breathe, which dropped the spaces and line breaks between
  its parts (`conceptnxx::real=`); it now reads the initializer from the XML file.
- 31 of 31 tests pass on `clang-cl`, `cl`, GCC 16 (MinGW), Clang 22 (llvm-mingw) and Emscripten, with no warning:
  `numerixx_tests` with the new `roots: the examples of the header comments` (15 assertions) and the cases added to
  `onesided: domain edges` and `ridders: accuracy ...`, and `numerixx_docs_example_getting_started`.
- `ForwardRidders` at the edge 0.5 (`deriv.hpp:310`), from a run of the independent prototype
  `.scratch/onesided/ridders5.py`: relative error 9.4e-15, estimate 1.7e-12, 10 evaluations. The comment's "about
  1e-13" is ten times the error: conservative, not wrong.

## Deferred

- `linkcheck` (a target and a CI step): with hosting. Until then only internal references exist that matter, and
  `-W` with `nitpicky` checks them; a network step on every push would make CI flaky.
- Hosting and versioning (GitHub Pages or Read the Docs; per release or `dev` only), and the install section of
  getting-started, pinned to a release tag: before the first release.
- An examples folder with a glob, and the math extensions: Phase 3, at the second example and the first formula.
- A separate concept check: if a second public concept appears.
- Offering the clause patch to Breathe upstream: the user's decision.
- `@snippet` in the headers instead of `@code`: a change to `include/`, after Phase 4 shows whether the claims tests
  are enough.

## Reviews of this proposal

- `numerics-reviewer`: the `deriv.hpp:12` comment gives relative figures as absolute (reproduced, above); the counts
  are independent, not baselines (sources cited); the search bracket is exact; tolerances from the existing tests'
  derivations instead of "the claim times 10"; five blocks with claims the first draft missed. All adopted.
- `simplicity-critic`: no configure-time package probe (a Breathe version check in `conf.py`); the example in
  `tests/CMakeLists.txt`; the checks as counters in `conf.py`, not HTML scraping; claims in the module tests, not
  `docs:*`; no math extensions; linkcheck deferred; `-E`; a full `pip freeze`; Doxygen 1.18 required; the formatters
  covered (rendering verified); Ninja in CI; no install section; the reviews in the criteria. All adopted. Python
  comes from uv at the user's request, which also settles the reviewer's `actions/setup-python` point.

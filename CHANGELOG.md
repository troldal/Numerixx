# Changelog

## Unreleased: Numerixx 2.0.0

Numerixx 2 is a rewrite; see [docs/redesign/PLAN.md](docs/redesign/PLAN.md). The previous API is preserved at the
tags `v1.0.0` (master) and `v1.1.0-legacy` (the dev-reorg branch). [MIGRATION.md](MIGRATION.md) maps it to the new
one.

### Phase 0: build skeleton

- New CMake build: C++23, one include root (`<numerixx/...>`), one INTERFACE target per module
  (`numerixx::core`, `::deriv`, `::roots`, `::optimize`, `::poly`, `::integrate`, `::interpolate`, `::linalg`,
  `::multiroots`, `::pipes`, `::multiprecision`) and the umbrella `numerixx::numerixx`; install and export rules for
  `find_package(numerixx 2.0 CONFIG)`, which install the fetched FXT and Eigen headers (and their licences) below
  `include/numerixx-deps`, never over another installation.
- Dependencies through CPM 0.43.2, each pinned by version and SHA256: FXT (only for `numerixx::pipes`), Eigen 5.0.1
  (only for `numerixx::linalg` and `numerixx::multiroots`), standalone Boost.Multiprecision (only for the optional
  adapter). doctest, google/benchmark and Boost.Math are fetched for development only.
- Removed: vcpkg, Blaze, LAPACK, OpenBLAS, OpenMP, gcem, tl-expected, nlohmann-json, fmt, Boost.Stacktrace,
  Boost.MultiArray, the vendored Google Benchmark copy, and the old library, tests, demos and documentation.
- CMake presets for GCC (with library assertions), Clang + libc++ (with sanitizers and hardening), MSVC, clang-cl,
  and Emscripten with wasm, JavaScript or no exceptions and with `-pthread`; GitHub Actions CI with every preset, the
  consumer-build scenarios and a format check; nightly floor-compiler legs.
- Test infrastructure: doctest 2.5.3 with test discovery (one CTest test per test case); header self-containment against each module's own target; a
  strict-warnings consumer translation unit; a layering check against the module DAG; compile-fail tests that check
  the diagnostic; and the consumer-build scenarios (CPM and FetchContent parents in both declaration orders, a
  parent with its own `Boost` package in both orders, a scalar-only parent, an installed package).
- Deferred: the MSVC P2564 (consteval escalation) probe arrives with the refined types in phase 1.

### Spike (roadmap phase S; DESIGN §10.2)

- FXT is pinned to a commit with the FXT-1 probe fix (troldal/FXT#1), so the FXT pipes build without exceptions:
  `gcc-noexcept` now includes `numerixx::pipes`, and the `gcc-noexcept-pipes` canary is gone.
- First cut of the core vocabulary (`<numerixx/core.hpp>`): scalar traits and maths helpers, refined types
  (`tolerance`, `abs_tolerance`, `rel_tolerance`, `evaluation_budget`, `max_iterations`, `bracket`) with consteval
  literals and `make()`, errors and results (`errc`, `fault`, `failure`, `solution`, `result`, `best_x`, `is_fatal`),
  evaluation (`evaluate`, `evaluate_sample`, `cost_of`, the unwrapping and common-cause rules), stop criteria with view
  kinds (`x_tol`, `step_tol`, `width_tol`, `floored_width`, `f_tol`, `max_evaluations`, `min_iterations`, `never`,
  `||`, `&&`), the driver `nxx::iterate`, `steps_view`, the family facades with generic builders (`with_stop`,
  `with_budget`, `with_derivative`, `with_projection`, `with_observer`), `.on()`, the named combinators (`first_of`,
  `first_of_with`, `then`, `warm_fallback`), `fn::counted`, and the opt-in run-time chains
  (`<numerixx/core/any_solver.hpp>`).
- First cut of 1-D root finding (`<numerixx/roots.hpp>`): `bisection`, `brent`, `secant`, `newton` (with a derivative
  source: a callable, the `deriv::numeric` policy, or a structural `.derivative()`), `expand` and
  `solve(f, {lo, hi})`; run-time brackets (`{lo, hi}`, `std::pair`, `bracket<T>::make`) validated in-band; per-iterate
  projection (`clamp_to`); the pole check.
- First cut of numerical differentiation (`<numerixx/deriv.hpp>`): stencils as integer data, `optimal`/`relative`/
  `absolute` steps, `diff`, `central`, `derivative_of` and the `numeric` policy.
- Tests: 156 doctest cases (criterion soundness, determinism with a golden table of 22 solves that is bit-identical on
  GCC, Clang, MSVC, clang-cl and em++, run-time chains equal to static chains, evaluation counts equal to instrumented
  calls, poles, extreme brackets, the canonical calls with run-time inputs, regularity, composition), 29 compile-fail
  cases, each with a control that must compile, whose reason must appear in the first error on GCC and Clang
  (diagnostic line counts are written to `compile_fail_report.txt` and recorded in DESIGN Appendix D), two harness
  self-tests (one of which the harness must reject), the P2564 probe (it must fail on MSVC with C7595), a
  strict-warnings consumer TU that instantiates the library with a global `f`, and compile-time measurements (the
  umbrella header is guarded at 2 s on GCC; measured 0.58 s on GCC 16, 0.46 s on Clang 22, 0.55 s on MSVC, 0.53 s on
  clang-cl). All 12 presets pass: 225 CTest tests on most, 223 on MSVC, 234 with multiprecision, 8 on integration.
- Found and fixed by the spike (recorded in DESIGN): Clang's default floating-point contraction made solver paths
  platform-dependent, so the headers turn it off for library code; `better_than` was not transitive, so the best
  estimate of a chain depended on how it was grouped; the solver constant for the view kind is spelled `views`;
  bisection's default budget is 200 (enough for `cpp_bin_float_50`); the driver's `finish` hook returns an optional
  failure (the sketched shape triggered a false GCC warning); builders that a solver would ignore are deleted with a
  reason.
- Found by the spike's review and fixed: brent reported `stop_reason::criterion` for enclosures up to 4ε|b| wider than
  the tolerance, and measured the relative part at b; it now reports `criterion` only when the width criterion holds,
  and `resolution_limit` when the tolerance is below its floor. The pole check was switched off by an infinite end
  sample (1/x on [-1, 0] was reported as a root); it now uses the finite samples and rejects a non-finite fx.
  `expand` could grow to an infinite endpoint, after which bisection reported a false `resolution_limit` success; it
  now saturates at ±max and stalls there. secant's second point falls back to x0 − h when x0 + h overflows. The
  facades now check that the function can be called with the input's scalar type, so `std::is_invocable_v` is
  `false` instead of a hard error; so does `bound` (what `.on()` returns). The combinators keep a `static_assert`,
  which gives the one-line "did you forget .on(input)?" diagnostic. `any_solver` is no
  longer a conversion target for every type. `evaluation_budget` rejects a negative literal instead of wrapping it.
  MSVC: C-array brackets are no longer ambiguous, and plain functions no longer trigger warning C4180.
  `root_estimate<double>{1.0}` gave an `exact_zero` success without evaluating f; a `root_estimate` now needs x and
  f(x). `min_iterations` alone (or under `||`) reported success after n iterations; it is now a guard that a solver
  accepts only as `test && min_iterations{n}`, and anything else is deleted with a reason. `best_x` no longer
  hard-errors on a search result (it is constrained away). The best estimate of a failure with neither an enclosure
  nor a step reports an unknown (infinite) uncertainty instead of the window width. The compile-fail harness now
  requires a deletion reason in the compiler's own message, not in a source line it quotes, and the MSVC P2564 probe
  must fail with C7595. The minor findings left open are listed in DESIGN Appendix D.
- Found by an audit of the spike's documentation against the code, and fixed: `expand`'s public `rebuild(options)`
  accepted a stop criterion, which could forge a sign bracket without a sign change (it now takes only `never{}`,
  and every solver's `rebuild` rejects a bare `min_iterations`); secant started at an exact root reported an
  uncertainty of 0 (now inf, as Newton does); a braced `{lo, hi}` with a function of the wrong signature now gets the
  facade's reason; a one-element braced list `{x}` bound to the two-element array overload as `{x, 0}` and was solved
  on [0, x], so a list or array is now a bracket only with two ends of a real type (anything else is deleted with a
  reason); `.on(C array)` is no longer ambiguous on MSVC, and `.on(pointer)` keeps its reason; `NXX_DELETE` now starts on the declarator's own
  line, where cl's note points; the Clang header brackets now save and restore the floating-point setting on RISC-V,
  PowerPC and SystemZ too. The audit also found that GCC contracts `a * b + c` by default in C++, ISO mode included,
  and that the headers cannot stop it; DESIGN §5.3 documents this (build with `-ffp-contract=off` for bit-identical
  results on FMA targets), and phase 1 decides whether the GCC interface flags should add it.
- `examples/quick_tour.cpp`: a tour of the library as it is now, built with the strict warning flags and run as a
  smoke test.
- Known limits, documented in DESIGN §7.2: `expand` from a window on one side of 0 cannot cross 0 (phase 3), and a
  large but finite initial sample can hide a pole from the pole check (add `&& f_tol{…}` for a residual guarantee).
- Found by hosted CI and fixed: GCC 16.2 (the `gcc:16` container) reported a false `-Wmaybe-uninitialized` when a
  solver holding a lambda that captures a `std::vector` was copy-assigned, which failed every `-Werror` GCC build of
  the spike; GCC 16.1 locally did not warn. `copyable_box` now copies such a capture into a temporary and moves it in
  place, without `std::optional`; a throwing copy now leaves the box unchanged instead of empty.
- Found by the nightly legs, red on master since their first run on 2026-09-29. GCC 14 (the `gcc:14` container)
  reports 8 `-Wnull-dereference` false positives inside Eigen 5.0.1's LU kernel, which `-isystem` does not hide;
  `tests/linalg/test_linalg.cpp` now silences them with a pragma region around `<Eigen/LU>` (DESIGN §5.3; checked
  with GCC 16.1, which reports the same 8 without `NDEBUG`, and with GCC 14.4 in the nightly). The Intel
  ICX leg compiled against Ubuntu 24.04's default libstdc++ 13, which is below the floor and lacks
  `std::forward_like`; Ubuntu's libstdc++ 14.2 is not enough either (see the floor entry below). The leg now installs
  g++-14 14.3 from Ubuntu's toolchain PPA, selects it with `--gcc-install-dir`, stops early if the library or the
  floating-point model is not what it tests, and pins its image by digest. All five nightly legs passed on the spike
  branch on 2026-10-02.
- Compiler floor: Clang 19 + libc++ 19, raised from 18 on 2026-10-01. The first nightly run on the spike branch
  showed that Clang 18 rejects the refined literals: `nxx::tolerance t{1e-8}` needs class template argument deduction
  for alias templates (P1814), and the consteval literal checks need P2448. A Clang-family compiler on libstdc++
  needs libstdc++ 14.3 or newer (DESIGN D2).
- Documented: Numerixx does not support `-ffast-math`, `-ffinite-math-only` or icpx's default fast floating-point
  model, because its failure checks rely on infinities and NaN (DESIGN §5.3). The ICX leg builds with
  `-fp-model=precise`.
- The spike's code is kept: phases 1-3 continue from it, with unchanged scope and acceptance criteria (DESIGN
  §10.3 lists what the spike built and what is left).
- Deferred to phases 1-3: the progress window and step-length cap of the open methods, the representation-space
  bisection midpoint, `illinois`, `ridders`, `rtsafe`, `scan`, `subdivide`, `inverse_of`, `.from_enclosure()`,
  `with_evaluation_budget`, `solve(f, x0)` and `solve(f, df, x0)`, noise steps, `diff_with_error`, Ridders
  differentiation and the mixed partials.

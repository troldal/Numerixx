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
- Tests: 154 doctest cases (criterion soundness, determinism with a golden table of 22 solves that is bit-identical on
  GCC, Clang, MSVC, clang-cl and em++, run-time chains equal to static chains, evaluation counts equal to instrumented
  calls, poles, extreme brackets, the canonical calls with run-time inputs, regularity, composition), 24 compile-fail
  cases, each with a control that must compile, whose reason must appear in the first error on GCC and Clang
  (diagnostic line counts are written to `compile_fail_report.txt` and recorded in DESIGN Appendix D), the MSVC-only
  P2564 probe, a strict-warnings consumer TU that instantiates the library with a global `f`, and compile-time
  measurements (the umbrella header is guarded at 2 s on GCC; measured 0.58 s on GCC 16, 0.46 s on Clang 22, 0.55 s
  on MSVC, 0.53 s on clang-cl). All 12 presets pass.
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
- Known limits, documented in DESIGN §7.2: `expand` from a window on one side of 0 cannot cross 0 (phase 3), and a
  large but finite initial sample can hide a pole from the pole check (add `&& f_tol{…}` for a residual guarantee).
- Deferred to phases 1-3: the progress window and step-length cap of the open methods, the representation-space
  bisection midpoint, `illinois`, `ridders`, `rtsafe`, `scan`, `subdivide`, `inverse_of`, `.from_enclosure()`,
  `with_evaluation_budget`, `solve(f, x0)` and `solve(f, df, x0)`, noise steps, `diff_with_error`, Ridders
  differentiation and the mixed partials.

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
- Tests: 164 doctest cases (criterion soundness, determinism with a golden table of 22 solves that is bit-identical on
  GCC, Clang, MSVC, clang-cl and em++, run-time chains equal to static chains, evaluation counts equal to instrumented
  calls, poles, extreme brackets, the canonical calls with run-time inputs, regularity, composition), 31 compile-fail
  cases, each with a control that must compile, whose reason must appear in the first error on GCC and Clang
  (diagnostic line counts are written to `compile_fail_report.txt` and recorded in DESIGN Appendix D), two harness
  self-tests (one of which the harness must reject), the P2564 probe (it must fail on MSVC with C7595), a
  strict-warnings consumer TU that instantiates the library with a global `f`, and compile-time measurements (the
  umbrella header is guarded at 2 s on GCC; measured 0.58 s on GCC 16, 0.46 s on Clang 22, 0.55 s on MSVC, 0.53 s on
  clang-cl). All 12 presets pass: 237 CTest tests on most, 235 on MSVC, 246 with multiprecision, 8 on integration.
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
  results on FMA targets); the phase-1 design (DESIGN §12.20) leaves the flag to the consumer and does not add it to
  the GCC interface flags.
- `examples/quick_tour.cpp`: a tour of the library as it is now, built with the strict warning flags and run as a
  smoke test.
- Known limits, documented in DESIGN §7.2: `expand` from a window on one side of 0 cannot cross 0 (phase 3), and a
  large but finite initial sample can hide a pole from the pole check (for a residual guarantee, check |fx| of the
  result, or use bisection with `&& f_tol{…}` and accept only `stop_reason::criterion` or `exact_zero`; brent's
  `with_stop(f_tol{…})` is only an early exit).
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
- Found by the review of PR #3 and fixed (DESIGN Appendix D): secant and newton reject a projected iterate that is
  not finite, which gave false `exact_zero` and `criterion` successes at x = ±inf or NaN; `width_tol` and `x_tol`
  thresholds saturate at the largest finite value, so brent no longer reports `criterion` at 0 iterations when
  abs + rel·min(|lo|, |hi|) overflows; `real` requires a specialised `std::numeric_limits`, because a type with only a
  `scalar_traits` specialisation got collapsed defaults and false successes; `diff` returns `invalid_input` before any
  evaluation when a stencil point is not finite; `solve(f, {1.0})` no longer solves on [0, 1]; with
  `<numerixx/core/any_solver.hpp>` included, the one-solver `first_of` and `first_of_with` compile again;
  `NUMERIXX_BUILD_EXAMPLES` defaults to `PROJECT_IS_TOP_LEVEL`, so projects that pull Numerixx in through CPM or
  FetchContent no longer build its examples. A follow-up review of these fixes found more of the same: the criteria of
  secant and newton now see the larger of the proposed and the projected step, so a projection to a far finite value
  no longer gives a `criterion` success there; secant's second point falls back to x0 - h when x0 + h is projected off
  the reals; `solve` is constrained on the function, so `std::is_invocable_v` is false for one that cannot take the
  bracket's scalar type, and `solve(f, 1.0)` states that it needs a bracket; `is_real_v` of an array or function type is
  false instead of a hard error.
- Documented: Numerixx does not support `-ffast-math`, `-ffinite-math-only` or icpx's default fast floating-point
  model, because its failure checks rely on infinities and NaN (DESIGN §5.3). The ICX leg builds with
  `-fp-model=precise`.
- The spike's code is kept: phases 1-3 continue from it, with unchanged scope and acceptance criteria (DESIGN
  §10.3 lists what the spike built and what is left).
- Deferred to phases 1-3: the progress window and step-length cap of the open methods, the representation-space
  bisection midpoint, `illinois`, `ridders`, `rtsafe`, `scan`, `subdivide`, `inverse_of`, `.from_enclosure()`,
  `with_evaluation_budget`, `solve(f, x0)` and `solve(f, df, x0)`, noise steps, `diff_with_error`, Ridders
  differentiation and the mixed partials.

### Phase 1: core vocabulary (DESIGN §10.3)

- Fixed: `brent` took a width criterion through `with_stop` and through its public `rebuild(options)`, and reported
  `stop_reason::criterion` where its own tolerance held instead of the given one.
  `brent{}.with_stop(width_tol{1e-20, 0})`, alone or with `|| f_tol{…}` or `&& f_tol{…}`, stopped at width 6.66e-16
  on x² − 2 over [1, 2]. Its intrinsic test runs before the stop criterion, so `with_stop` can only add an early exit
  (`f_tol`) or a failure (`max_evaluations`). `with_stop` and `rebuild` now reject any stop criterion that contains a
  width criterion (`width_tol`, `floored_width`, at any depth of `||` and `&&`), with a reason that points to the
  constructor: `brent{nxx::width_tol{1e-20, 0}}` reports `resolution_limit` at that width. The rule is keyed on a
  facade trait, `internal_tolerance`, for the later solvers with their own tolerance, and `with_stop` and every
  solver's `rebuild` share one predicate, `detail::stop_allowed_v` (DESIGN §6.8). brent has no spelling that
  guarantees a residual, and never had one: `with_stop(floored_width{} && f_tol{1e-12})` reported `criterion` at
  |fx| = 1 on a jump from −1 to +1, and `with_stop(f_tol{…})` is OR-ed with its tolerance. The docs that advised
  `&& f_tol{…}` for a residual guarantee now restrict it to bisection (DESIGN §7.2, §9.3).
- "newton needs a derivative" is now deleted in `open_facade`, keyed on the solver's `ready_v`, so `newton` declares no
  call operator and no `using open_facade::operator();`. Calls, `std::is_invocable_v` and the reason text are
  unchanged. GCC and cl now name `nxx::open_facade::operator()` in the error, and Clang lists the candidate in the
  notes of secant's misuses. CLion's ReSharper C++ engine (2026.2) took Newton's deleted overload to hide the facade's,
  and marked every valid Newton call as an error (DESIGN §6.6).
- The examples print with `std::println`, not `std::printf`. Their format strings are checked at compile time, and
  doubles print in the shortest form that reads back to the same value, not with 17 significant digits. The
  quick tour's error and stop-reason names are `constexpr std::string_view` functions with `using enum`, and its
  stepping loop counts with a range-for initializer. With MinGW's libstdc++ (GCC 16.1 still), `std::print` needs
  libstdc++exp. `examples/CMakeLists.txt` probes for it and links it to the examples only where the toolchain needs it.
  The library headers do not use `<print>`, so consumers link nothing extra.
- A bare number as a solver's tolerance gets a reason. `brent{1e-10}` was a hard error inside `brent.hpp` ("'applies_to'
  is not a member of double"), and so were `std::is_constructible_v<brent<double>, double>` and CTAD probes: its
  `width_tolerance_v` used `&&` in a variable template's initializer, which formed `W::applies_to` for every `W`. The
  trait is now false for a non-criterion, and `brent`, `bisection`, `secant` and `newton` delete a number with "a
  tolerance is a criterion, not a number: write brent{nxx::width_tol{1e-10}}" (each solver names its own criterion).
  `bisection{1e-10}`, `secant{1e-10}` and `newton{1e-10}` failed class template argument deduction without a reason.
  A number that the solver's own criterion type converts from is still accepted, so a user criterion with a
  converting constructor keeps `brent<my_width>{tol}`. A deduction guide sends a bare number to `brent<>`: Clang 19.1
  deduced `brent<double>` from the implicit guide of `brent(Tol)` despite its constraint, and gave a bare "no matching
  constructor" (found by the nightly clang-floor leg). Each solver has a compile-fail case for the reason.
- An open method given a braced list or a C array (`secant{}(f, {1.0, 2.0})`, `newton{}.with_derivative(df).on({1.0,
  2.0})`) is deleted with "open methods take one guess of a real type or a root estimate, not a braced list: write 1.0,
  or pass {lo, hi} to a bracketing solver", not a bare "no matching function". `open_facade`'s `.on` catch-all takes a
  forwarding reference, as the other facades' do, so a named array reaches the array overload on cl.
- The facades check the solver protocol (DESIGN §6.6). `std::is_invocable_v` is now `false` for a solver without
  `accepts_v`, or without a member that `prepare(f, in)` or `nxx::iterate` needs (`id`, `options`, `init`, `step`,
  `view`, `estimate`, `best`, `intrinsic`), or whose `prepare` does not return a `std::expected`, where it was a hard
  error inside `detail::run`; a call names the protocol in a deleted overload. A member of the wrong type (a `prepare`
  error that `init`'s failure type cannot hold, an `options()` that is not an options aggregate, an `init` error that is
  not a `nxx::failure`) is still a hard error inside `detail::run` or `nxx::iterate`. `with_stop` on a facade-derived
  type without `views` is `false` too, where `detail::stop_allowed_v` made it a hard error, and so is one whose
  `views` is not a `view_kind` constant, on GCC, Clang and clang-cl (cl 19.51 still rejects that case with a hard
  error, as before). On Clang 22.1.8, misuses of
  the call operators list the new candidates in their notes (`bisection_given_guess` 35 / 27 to 44 / 36 lines,
  `open_int_guess` 38 / 30 to 47 / 39 with the braced-list deletion, `newton_mixed_errors` 62 / 24 to 84 / 34); the
  error line and its reason are unchanged. On GCC 16.1 only `newton_mixed_errors` changed (27 / 4 to 33 / 4).
- DESIGN §8.1 and `pipes.hpp` no longer say that no other `operator|` is declared in `nxx`: `operator|(view_kind,
  view_kind)` is, and it never competes with the pipe.
- Tests for these fixes: 6 compile-fail cases (`brent_number_tolerance`, `bisection_number_tolerance`,
  `secant_number_tolerance`, `open_braced_bracket`, `newton_on_braced_bracket`, `solver_incomplete`), static asserts in
  `tests/roots/test_solvers.cpp` and one doctest case, a client solver run under its own id. `gcc`, `gcc-noexcept`
  and `clang` run 255 CTest tests (242 before), `gcc-multiprecision` 264 (251 before).

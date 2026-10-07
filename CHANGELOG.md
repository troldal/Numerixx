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
- Tests for these fixes: 7 compile-fail cases (`brent_number_tolerance`, `bisection_number_tolerance`,
  `secant_number_tolerance`, `newton_number_tolerance`, `open_braced_bracket`, `newton_on_braced_bracket`,
  `solver_incomplete`), static asserts in `tests/roots/test_solvers.cpp` and one doctest case, a client solver run
  under its own id. `gcc`, `gcc-noexcept` and `clang` run 257 CTest tests (242 before), `gcc-multiprecision` 266 (251
  before), as hosted CI measured on `master` at e6b8e44 (ci.yml run 37463362893); this entry said 255 and 264, the
  counts before `newton_number_tolerance` was added.
- Renamed: `failure::where` is now `failure::by`, the name `solution` already used, so a result's algorithm is
  `res ? res->by : res.error().by`; and `fault::evals` is now `fault::evaluations`, the name `counters` uses. Both are
  Numerixx 2 alpha spellings: the layout and the positional construction (`failure{code, id, used, best, cause}`) are
  unchanged, so only code that reads `.where` or `.evals` changes (DESIGN §6.3, D7, §12.20).
- Added `nxx::best(r)`: the solution's estimate, or the failure's best estimate, as one `std::optional<Est>`
  (`std::nullopt` only when no evaluation succeeded), for every result whose success and failure carry the same estimate,
  through `first_of` and `any_solver` too. On a search result, which succeeds with a `sign_bracket` and fails with a
  `root_estimate`, it is deleted with "nxx::best: this result succeeds and fails with different estimates (a search
  result: a sign_bracket, then a root_estimate): read *r and r.error().best separately". A one-shot derivative's
  `std::expected<T, fault<UE>>` is not a result, and `best` is not invocable on it. `best_x` is unchanged; its comment
  now names the remedy for a search result: read `r->lo()` and `r->hi()` (DESIGN §6.3, §6.14 call 15, §12.21).
- Tests for these changes: `tests/roots/test_results.cpp` (4 doctest cases, and static asserts for `best`, `best_x`
  and `detail::is_result_v`, with `best` in a constant expression), the compile-fail case `best_search_result` with its
  reason, and canonical call 15 in `tests/usage/canonical_calls.cpp`. `gcc` runs 264 CTest tests (7 added),
  `gcc-multiprecision` 273 (7 added).
- Changed: a NaN or infinite input value fails with `errc::non_finite_input`, cost {0, 0} and no best estimate, in
  every form. A bracket end (`bracket<T>::make`, and so braced `{lo, hi}`, `std::pair`, `.on(...)`, `solve` and
  `expand`'s window) and the x of `deriv::diff`, `central` and `derivative_of(f)(x)` gave `invalid_input`, while guesses,
  root estimates and projected starts already gave `non_finite_input`. Equal ends, a derivative step that vanishes or
  overflows, and a stencil point that overflows at a finite x keep `invalid_input`, and so does every bad value given
  to a refined type's `make()` (DESIGN §6.3, §6.5, §7.1, §12.20 decision 11).
- Fixed: an input code from a nested callable no longer surfaces as the solver's own input error. When a callback
  returns a Numerixx fault (`derivative_of(g)` used as f, the `diff` inside Newton's numeric derivative, or a fallible
  callback that returns `fault<UE>`), `nxx::evaluate` and `evaluate_sample` now turn its `invalid_input` and
  `non_finite_input` into `non_finite_value`, keeping the fault's evaluations and cause, at every evaluation: in
  `prepare`, in `init` and in a step. Before, Newton with `deriv::numeric{}` on log(x) − log(1.79769e308) from 1e308
  failed with `invalid_input` after 5 iterations and 13 evaluations, when a stencil point overflowed, and brent on
  `derivative_of(g)` over [−max, 1] failed with `invalid_input` at cost {0, 0}, exactly like a rejected bracket; so
  `first_of_with` with a policy that stops on input errors stopped there, or not, depending on the order of its
  alternatives. Direct calls keep their codes: `diff(f, nan)` gives `non_finite_input`, and `derivative_of(f)(x)` with
  an overflowing stencil gives `invalid_input`. `nxx::iterate` and `steps_view` also step through
  `detail::checked_step`, which maps the same two codes from a user-written step that returns one directly. The
  guarantee covers these two codes, the ones Numerixx's own callables produce: a fallible callback or a nested solve
  can still pass another input code (`no_sign_change`, `out_of_domain`, ...) through, so `is_input_error(code)` alone
  does not mean "rejected before iterating" (DESIGN §6.3, §6.4, §6.7, §12.20 choice 1, §12 item 22).
- Tests for these changes: the `bracket<T>::make` rows in `tests/core/test_refined.cpp`, every bracket input form in
  `tests/roots/test_solvers.cpp`, a new case for a non-finite x in `tests/deriv/test_deriv.cpp`, `step_fault` and a
  user-written step that returns `invalid_input` directly in `tests/roots/test_steps.cpp`, four cases in
  `tests/usage/test_composition.cpp` (the Newton regression row and the secant with `derivative_of(g)` as f, the
  `first_of_with` fall-through, `steps_view` against the driver, and one row per evaluation site: brent's first and
  second sample, `solve`, `expand`'s two samples, the secant's x0, x1 and a seeded x1, Newton's x0, `steps_view`'s
  element 0, `nxx::evaluate` and `diff` over a nested derivative, the direct calls, a fallible callback's fault, and the
  stop-on-input chain in both orders), and a `cpp_bin_float_50` case. `gcc` runs 270 CTest tests (6 added),
  `gcc-multiprecision` 280 (7 added).
- Added `nxx::better_than(a, b)`, the customisation point for the order of failure estimates, which picks the best
  estimate of a solve and of a chain. It calls the `better_than(const Est&, const Est&)` that ADL finds next to the
  estimate type, and is not invocable when there is none; a member function or a function in another namespace does
  not count. The driver's fallback to `merit_of(e) < merit_of(best)` is removed, and `iterative_solver_for` requires
  an order for the failure estimate type, so a solver whose estimate has none makes `std::is_invocable_v` false, and
  the facades' protocol reason now ends "…, and better_than(const Est&, const Est&) for its estimate type". roots'
  order is a hidden friend of `root_estimate`: `nxx::roots::better_than`, a namespace-scope function, is removed (call
  `nxx::better_than`), so the unqualified name also works under `using namespace nxx;` and `using namespace
  nxx::roots;`. A side effect: a solver whose `init` error has no `estimate_type` is now `false` through the same
  reason, where it was a hard error inside `nxx::iterate` (DESIGN §6.6, §12.20 decision 9, §12.21 item 7).
- Fixed: roots' order of failure estimates was not a strict weak order once a NaN |f(x)| entered, because it compared
  |f(x)| with `<`, so the best estimate of a chain could depend on how its alternatives were grouped or folded. Two
  enclosures whose widths both overflow to inf also tied, whatever their size. The order is now: an estimate with an
  enclosure first; the smaller width, and the smaller hi/2 − lo/2 only when both widths overflow; then the smaller
  |f(x)|, with a NaN last. It still never prefers a strictly wider enclosure, which a half-width key would do near the
  subnormal range ([2d, 5d] before [d, 3d], d = `denorm_min()`). `sign_bracket`'s constructor gains the precondition
  that both ends are finite, checked in assert builds (DESIGN §6.7, §7.2, §12.21 item 8).
- Fixed: a pole failure (`sign_change_not_root`) carried the final enclosure as its best estimate. The enclosure holds
  the pole, not a root, and it outranked every estimate without an enclosure, so a chain reported the pole as its best
  estimate: `first_of(brent{}.on({1.0, 2.0}), secant{}.with_budget(2).on(3.0))` on tan failed with best x =
  1.5707963267948974 and |f| = 1.21e15. The failure now carries the estimate without its enclosure (x, f(x) and the
  uncertainty unchanged, so the pole is at `best->x`), and that chain's best is the secant's x = 3.1415807758403682
  with |f| = 1.19e-5. Known limits until phase 3: an enclosure around a pole that the solver did not detect (a
  bracketing failure that ends before its pole check: out of budget, or a step fault next to the pole) still ranks
  first, and after a `first_of`, whose failure keeps the last alternative's code
  with the best alternative's estimate, a `sign_change_not_root` failure can carry another alternative's enclosure
  (DESIGN §6.7, §6.10, §7.2, §12.20 choice 4).
- Tests for these changes: `tests/roots/test_order.cpp` (7 doctest cases, each over `float`, `double` and
  `long double`: NaN against finite in both orders, NaN against NaN, two overflowing widths, a finite width before an
  overflowing one, [d, 3d] before [2d, 5d], a fixed-seed property test of the four strict-weak-order axioms (the two
  transitivity axioms over every third of its 160 estimates, a pool drawn the same on every standard library and
  checked to hold an enclosure whose width overflows), of agreement with the key (e, o, k, n, a) and
  of nested enclosures, and `first_of` over three extreme estimates in every
  order and grouping, static and run-time), concept tests for estimates without an order in
  `tests/core/test_criteria.cpp`, the compile-fail case `solver_without_better_than` with its reason, the pole-payload
  row in `tests/roots/test_combinators.cpp`, the pole cases of `tests/roots/test_solvers.cpp` (now checking that the
  estimate has no enclosure), and a `cpp_bin_float_50` case. `gcc` runs 294 CTest tests (24 added),
  `gcc-multiprecision` 305 (25 added).
- Fixed: `nxx::better_than` was declared `noexcept` whenever the user's `better_than` was, although it also converts
  that result to `bool`, which may throw; it is now `noexcept` only when both are. A static assert in
  `tests/core/test_criteria.cpp` checks an order whose result converts through a conversion that is not `noexcept`
  (DESIGN §6.6).
- Fixed: on MSVC (cl), a `first_of` nested as the first alternative of another, `first_of(first_of(a, b), c)` or
  `first_of_with` with the same policy type at both levels, gave wrong results. `first_of_t` declared its leading empty
  policy `[[msvc::no_unique_address]]`, and cl 19.51 then overlapped the inner chain with the outer chain's next
  alternative: over newton, secant and bisection on x² + 1, the 40-byte chain held a 32-byte inner chain and a 24-byte
  bisection, and it failed after 5 iterations and 10 evaluations instead of 6 and 11. The policy no longer has the
  attribute, on any compiler, so that the type is correct on cl; a chain is a byte plus padding larger per link. cl
  and clang-cl still lay out `nxx::options`, and so every solver and chain, differently (a nested empty
  `[[msvc::no_unique_address]]` box: a curried secant is 16 bytes on cl and 24 on clang-cl), so mixing cl and clang-cl
  translation units that share Numerixx types is unsupported (DESIGN §5.3). Tests: a case in
  `tests/roots/test_combinators.cpp` with a static assert on the size of a left-nested `first_of` and a run-time
  comparison with the flat chain. `gcc` runs 295 CTest tests (1 added), `gcc-multiprecision` 306 (1 added).
- Fixed: stage 2 of `then` and `warm_fallback` starts from stage 1's output, not from the caller's input, but its
  `invalid_input` or `non_finite_input` reached the chain's caller as an input error: `then(secant{}.on(3.0),
  newton{}.with_derivative(df).with_projection(p))`, with a projection p that sends stage 1's root off the reals,
  failed with `non_finite_input` after 8 iterations and 10 evaluations, so a `first_of_with` chain that stops on input
  errors stopped or succeeded depending on the order of its alternatives. Such a failure is now `non_finite_value`,
  with stage 2's cost and cause. And `then` returned stage 2's failure as it was, so a stage 2 that failed at its own
  start left `best`, and `nxx::best(r)`, empty after stage 1's success and its evaluations; the failure now carries the
  better of stage 2's best and stage 1's estimate (a search's `sign_bracket` gives its `best()`), except when stage 2
  fails with `sign_change_not_root`: stage 2 has then shown that stage 1's bracket holds a pole, and the failure keeps
  stage 2's pole estimate without an enclosure, so a solver's pole payload holds through `then` (on tan,
  `then(expand{}.on({1.0, 2.0}), brent{})` keeps brent's x = 1.5707963267948974, not the search's window end x = 1
  with the enclosure [1, 2], and a `first_of` of it and a starved secant reports the secant's estimate). The rule
  reads stage 2's code, so a `first_of` stage 2 whose last alternative found a pole can still report another
  alternative's estimate, with its enclosure, under that code, a known limit until phase 3 (DESIGN §6.10). Like
  `first_of` and `warm_fallback`, `then` now needs an order (`better_than`) for its failure estimate type: a `then` over
  two hand-written stages whose estimate has none no longer compiles. The result types and the costs are unchanged
  (DESIGN §3.4, §6.7, §6.10, §12 item 23). Tests: two cases in `tests/roots/test_combinators.cpp` (the projected chains
  under `then` and `warm_fallback`, `first_of_with` with a stop-on-input policy in both orders, a clamp onto a point
  where f is NaN under a Newton and a secant stage 2, with every call of f counted, and a search stage 1), and a third
  for the pole (`then` over `expand` and brent or bisection on tan, with the solver's own cost and best estimate, and a
  `first_of` of the brent chain and the starved secant). `gcc`, `gcc-noexcept` and `clang` run 297 CTest tests (2
  added), `msvc` 295, `gcc-multiprecision` 308 (2 added), before the pole case and the next fix.
- Fixed: `warm_fallback` restarted its open method from a pole failure's `best->x` (`sign_change_not_root`), where
  Newton's step and the secant's are tiny, so the step criterion accepted the pole and the chain reported it as a
  success with `stop_reason::criterion`: on tan, `warm_fallback(brent{}.on({1.0, 2.0}), newton{}.with_derivative(df))`
  succeeded at |f| = 5.8e14, and `warm_fallback(bisection{}.on({1.0, 2.0}), secant{})` at |f| = 652. `warm_fallback`
  now returns stage 1's `sign_change_not_root` failure as it is, with its estimate and cost. Known limits until phase
  3, where the open-method safeguards and the pole check address them: the rule reads stage 1's code, so
  `warm_fallback` still restarts, and can succeed at the pole, from a bracketing failure that ended before its pole check next to a pole it had not detected (out of budget, or a
  step fault such as a callback that refuses x near its singularity; on tan, over `bisection{}.with_budget(40).on({1.0, 2.0})` and Newton, at |f| =
  6.7e11) and from a pole estimate that a `first_of` merge put under another code (over `first_of(brent{}.on({1.0,
  2.0}), bisection{}.on({1.0, 1.0}))` and Newton, at |f| = 5.8e14); it skips a valid restart when the last
  alternative of a `first_of` stage 1 found a pole; and an open method that the caller starts at a pole can still
  meet its step criterion there (measured on GCC 16.1 and Clang 22.1.8; DESIGN §3.4, §6.10, §7.2, §12 item 23).
  Tests: a case in `tests/roots/test_combinators.cpp` (brent and Newton, bisection and the secant, and
  `then(expand, brent)` and Newton, on tan, with stage 2 never run).
  With the `then` pole case above, `gcc`, `gcc-noexcept` and `clang` run 299 CTest tests (2 added), `msvc` 297,
  `gcc-multiprecision` 310 (2 added).
- Changed: the mixed criteria `x_tol` and `width_tol` name their relative part, so the two numbers can no longer be
  swapped. One number is absolute, `width_tol{1e-10}`; mixed is `width_tol{1e-10, nxx::rel_tolerance{1e-8}}` and
  purely relative `width_tol{0.0, nxx::rel_tolerance{1e-8}}` (the same for `x_tol`), each the same (abs, rel) pair
  as the old positional `width_tol{abs, rel}`, so no threshold changes. The two-number literal (`width_tol{1e-10, 1e-8}`,
  also with run-time numbers) is deleted with "say which number is relative: …", and a part alone
  (`width_tol{nxx::rel_tolerance{1e-8}}`, `width_tol{nxx::abs_tolerance{1e-10}}`) with "a part alone is not a
  criterion: …", which names both spellings and their run-time paths; the second was a deduction failure with no
  reason. At run time `make(T, T)` is deleted with a reason that names the direct remedy, "say which number is
  relative: make(a, *rel) with rel = rel_tolerance<T>::make(r); make(a) for an absolute tolerance"; `make(abs)` is
  added (finite and > 0, as `tolerance<T>::make`), `make(abs, rel_tolerance<T>)` mirrors the mixed literal and checks
  the absolute part in-band (finite and >= 0), and `make(abs_tolerance<T>, rel_tolerance<T>)` delegates to it; every
  failure is `invalid_input`. The mixed literal constructor is consteval, and cl 19.51 lacks P2564, so generic code
  that forwards the parts calls `make` (DESIGN §6.2 FLAG). bisection's reason for `x_tol` and brent's for a criterion
  that is not a width criterion now quote `width_tol{abs}, width_tol{abs, nxx::rel_tolerance{rel}} or floored_width{}`
  (DESIGN §3.3, §6.2, §6.8, §6.14 calls 13 and 14, §12.20 decisions 3 and 4, §12.21 items 4 and 5, §12.24).
- Added: a validated tolerance takes a relative part, `width_tol{*tol, *rel}` and `x_tol{*tol, *rel}` with `tol` a
  `tolerance<T>`, as a literal or with run-time values, and `make(tolerance<T>, rel_tolerance<T>)`. A `tolerance<T>`
  is finite and > 0, so the constructor checks nothing and is constexpr: it forwards through a constexpr template on
  cl too, and `make` cannot fail. Before, every way of adding a relative part to a validated tolerance failed with
  only the compiler's error (106 lines on GCC 16.1 for the literal). A literal with a `float` relative part,
  `width_tol{1e-10, nxx::rel_tolerance{1e-8f}}`, builds a `float` criterion and narrows the absolute part (to 0 below
  about 7e-46), with a warning only under `-Wconversion` on GCC and Clang: write both parts in one type (DESIGN §6.2,
  §12.24).
- Tests for these changes: the compile-fail cases `width_tol_two_numbers`, `x_tol_two_numbers`,
  `width_tol_make_two_numbers` (whose EXPECT names `make(a, *rel)`) and `width_tol_relative_alone`, each with its
  reason; `x_tol_zero_zero` rewritten to `x_tol{0.0, nxx::rel_tolerance{0.0}}`, and `rel_tolerance_as_tolerance`,
  which now reaches the part-alone reason (35 / 33 diagnostic lines to 15 / 13 on GCC 16.1, 20 / 19 to 8 / 7 on
  Clang 22.1.8); doctest cases in `tests/core/test_refined.cpp` (`make(a, *rel)` for a valid and a negative `a`,
  generic code that forwards the parts through `make`, a validated tolerance with a relative part as a literal, with
  explicit T, through `make`, forwarded and at run time, and concept and `is_constructible` tests for two numbers, a
  part alone, swapped roles, two tolerances and `make`) besides the rewritten literal and `make` case; the
  refined-forwarding case in the P2564 probe; canonical calls 13 and 14 in `tests/usage/canonical_calls.cpp`, call 13
  also with `width_tol{*abs, *rel}`; and two `cpp_bin_float_50` cases. The tests, the quick tour and the determinism
  table use the new spellings; the table's labels and values are unchanged.
- A validated tolerance, or one of its parts, given to a solver's constructor gets a reason. `brent{*tol}` and
  `bisection{*tol}`, with `tol` a `tolerance<double>` from `make()`, failed class template argument deduction with no
  reason (104 lines on GCC 16 and 59 on Clang 22 for bisection); they now report it in the first error (13 / 11 lines
  on GCC 16.1, 8 / 7 on Clang 22.1.8). No new deletion: the bare-number deletions of `brent`, `bisection`, `secant`
  and `newton`, and brent's deduction guide, now take `tolerance<T>`, `abs_tolerance<T>` and `rel_tolerance<T>` too,
  read from one trait, and the text gains the remedies: "…; a validated tolerance is not a criterion either; wrap it:
  brent{nxx::width_tol{*tol}}; a part (abs_tolerance, rel_tolerance) is not one either: build the criterion with
  nxx::width_tol<T>::make" (secant and newton name `x_tol`). A solver whose criterion type converts from a tolerance
  still takes one, `brent<width_tol<double>>{*tol}` (DESIGN §6.6, §7.2, §12.20 decision 12, §12.21 item 2).
- `with_stop` given something that is not a criterion, a number (`with_stop(1e-10)`), a validated tolerance or a part,
  got the false reason "this criterion does not apply to this solver …". A deleted sibling now says "with_stop takes a
  criterion, not a number or a validated tolerance: wrap it in the test you mean: width_tol{*tol} (bisection; brent
  takes its width in its constructor) bounds the error in x; x_tol{*tol} (open methods) bounds only the last step in
  x; f_tol{*tol} bounds only |f(x)|; …", and the criterion catch-all keeps its reason for criteria that do not apply.
  On Clang 22.1.8 the `with_stop` misuse cases list the new overload among their candidates (26 / 25 to 32 / 31
  lines); GCC's counts are unchanged (DESIGN §6.6, §12.21 item 3).
- A curried solver given a function it cannot take, `brent{}.on({1.0, 2.0})(g)`, gets a reason: "this solver cannot
  take this function with its bound input: call solver(f, input) for the reason". It was Clang's bare "no matching
  function for call to object of type 'bound<…>'" (GCC showed the facade's reason at line 33 of 41); `std::is_invocable_v`
  stays `false` (DESIGN §6.6).
- Tests for these changes: the compile-fail cases `brent_validated_tolerance` (through brent's guide),
  `bisection_validated_tolerance` (through the implicit guide), `bisection_with_stop_tolerance` and
  `bound_wrong_function`, each with its reason; static asserts in `tests/roots/test_solvers.cpp` that the four solvers
  reject a tolerance and both parts, that a solver over `width_tol` or `x_tol` still takes a tolerance but not a part,
  that `with_stop` on bisection, brent, secant, newton and expand rejects a number, a tolerance and a part while each
  remedy the text names compiles where it names it, and that a curried solver is not invocable with a wrong function;
  a doctest case that solves with each remedy from run-time tolerances; and `cpp_bin_float_50` static asserts (each of
  the four solvers rejects a tolerance and both parts, each over `width_tol` or `x_tol` takes a tolerance, bisection
  and secant over those reject a part, and `with_stop` on bisection, secant and brent) and a case. The existing
  `*_number_tolerance` and `with_stop` cases pass with unchanged regexes.
- Tests added by the test audit of this PR (DESIGN §6.2, §7.2): the compile-fail cases `x_tol_abs_part_alone` (the
  `x_tol` copy of the part-alone text, and the only case that reaches its `abs_tolerance` path, §12.21 item 5),
  `width_tol_zero_zero`, `width_tol_negative_abs_literal` and `width_tol_runtime_parts`, each with a control; the
  EXPECT of the two-number and part-alone cases now includes the criterion's name, and those four deletion cases spell
  the failing line as a plain variable, so that cl reports C2280 at the declaration rather than only C2131; and
  `is_constructible` checks that the legal mixed spellings stay constructible.
- The defaults are now checked achievable as the solvers use them, and in `cpp_bin_float_50` too (DESIGN §3.5, A1).
  A `TEST_CASE_TEMPLATE` in `tests/roots/test_solvers.cpp` reads each solver's own default criterion
  (`bisection<>{}.options().stop`, `brent<>{}.tolerance()`, `newton<>{}.options().stop`, `secant<>{}.options().stop`)
  and `static_assert`s that its thresholds are >= 4 eps for `float`, `double` and `long double`; the checks on
  `floored_width` and `step_tol` themselves passed with an unachievable solver default (with Newton's default changed to
  `step_tol<99, 100>`, the new asserts fail). `cpp_bin_float_50` is not a literal type, so a run-time case in the
  multiprecision test binary reads the same defaults at 0, 1 and 10^6 and compares them with 4 eps computed at run
  time (4 eps · 10^6 at 10^6). It also checks that the comparison rejects a threshold below 4 eps (`step_tol<99, 100>`,
  2^-167 at 168 digits); without the 4 eps floor in `floored_width::factor`, 6 of its 14 checks fail (DESIGN §12.21
  item 11).
- With the changes of this PR, all 12 presets pass from a fresh configure: `gcc`, `gcc-noexcept`, `clang`,
  `clang-asan`, `clang-cl` and the four Emscripten presets run 333 CTest tests (34 added), `msvc` 331,
  `gcc-multiprecision` 348 (38 added) and `integration` 8.

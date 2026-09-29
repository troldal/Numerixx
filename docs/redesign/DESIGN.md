# Numerixx 2 — Design Reference

- **Date:** 2026-09-27
- **Status:** Approved on 2026-09-28, with every §12 default accepted
- **Scope:** every module of Numerixx (core, deriv, poly, roots, optimize, linalg + multiroots, integrate, interpolate), the build system, tests and CI, upstream work in FXT, migration from Numerixx 1.x, and the families planned after v2.0 (§1.1, §10.5).

**Purpose.** This is the detailed companion to [`docs/redesign/PLAN.md`](PLAN.md), the short plan. It records the library's scope and positioning, every design decision with its rationale and the alternatives considered, the principles behind them, the build and dependency design, the core abstractions with code sketches, the per-module algorithm tables (including the bugs that must not be ported), the testing and CI design, the roadmap with phase sizes and acceptance criteria, the risks, and the decisions that were open until 2026-09-28, with their defaults (all accepted). Section numbers are stable, so other documents and code comments can cite them (for example "DESIGN §6.8"). Requirement IDs (R-B1, R-E3, R-A5, …) are listed in Appendix B.

**Evidence markers.**
- **[prototyped]**: compiled and run in the feasibility prototype preserved at [`docs/redesign/prototype/`](prototype/), on 9 configurations with bit-identical output: GCC 16.1 and Clang 22 + libc++, each with and without `-fno-exceptions`; em++ 6.0.8 with `-fexceptions`, `-fno-exceptions` and `-fwasm-exceptions`; MSVC 19.51; clang-cl 22. MSVC and clang-cl also passed with `/EHs-c-`.
- **[prototyped mechanism]**: the prototype verified the mechanism, but the spelling here differs (the prototype predates some renames).
- **[sketch]**: not yet compiled. Spike exit criteria 7–11 (§10.2) make the important sketches concrete first.

The prototype is throwaway evidence, not library code. It still contains an in-house constexpr LU (`nxx/linalg.hpp`) and zero-heap design choices that this design no longer requires. Claims that the prose calls *verified* or *measured* without a bracketed marker come from exploratory work that is not preserved: build and CMake experiments, the linear-algebra research, exploratory solver prototypes, and numerical and usability probes of a simpler predecessor API. Where the text says "a naive design would …", the figure quoted was measured on that predecessor.

---

## 1. Executive summary

**What changes.** Everything below the CMake target names. Numerixx 2 is a header-only C++23 library built from four things:
- **immutable solver values**, each a pure `init`/`step` pair;
- **one bounded iteration driver**;
- **stop criteria typed by what they can soundly judge**;
- **validated input types**.

There is one error model, built on `std::expected`. Any solver can be curried into a callable `f → result`, and named, assignable combinators compose solvers:
- `first_of`: if one fails, the next one tries;
- `then`: search, then solve, then polish;
- `warm_fallback`: restart the next solver from the last one's best estimate;
- `with_evaluation_budget`: one evaluation budget for a whole chain.

Chains are templates built at compile time by default. An opt-in `nxx::any_solver` type-erases a curried solver, so chains can also be assembled at run time (for example from configuration). Manual stepping and tracing go through one lazy range, `nxx::steps_view`. The library never throws; it builds and is CI-tested with `-fno-exceptions`.

Removed: vcpkg, Boost in the library, Blaze, LAPACK, OpenMP, gcem, tl-expected, the five duplicated driver loops, the CRTP state machines and the 21 `throw` sites. CPM 0.43.2 fetches FXT (for `<numerixx/pipes.hpp>`, default ON), Eigen 5.0.1 (the backend of `numerixx::linalg` and `numerixx::multiroots`, `NUMERIXX_WITH_LINALG`, default ON), and, optionally, standalone Boost for the multiprecision adapter and the test oracles.

**Safeguards a naive design lacks:**
- **Sound criteria.** A successive-iterate test reports false convergence on a bracketing solver: the prototype, which used one criterion type for all solvers, returned x = 1.5 for `bisection{x_tol{1e-9}}` on x²−2 **[prototyped]**. Here `x_tol` does not compile on a bracketing solver; bracketing solvers use `width_tol` or `floored_width`.
- **Assignable values.** Combinators and function-returning APIs are named class templates, so they can be copy-assigned (closures cannot).
- **No error nesting.** `derivative_of` returns a `fault`, so its errors do not nest inside Newton's.
- **Numeric-derivative Newton in curried chains.** A derivative *policy* makes this possible.
- **Run-time inputs are first-class.** Solvers accept braced `{lo, hi}`, `std::pair`, or the result of `make()`.
- **Safer numerical defaults:**
  - overflow-safe and representation-space bisection;
  - pole detection;
  - relative derivative steps (a floor-1 step would be 46 % of x at x = 1e-3);
  - rounding-aware open-method tests;
  - cycle and divergence detection;
  - honest evaluation counts;
  - N-D scaling;
  - QUADPACK error rules;
  - defaults that are functions of `T`.
- **Dependencies first.** deriv, a small module, lands before the driver and combinators: it needs only core, and the numeric-derivative Newton (phase 3) and the N-D Jacobians (phase 5) build on it.
- **Scope kept to what users need.** `first_of_all` is deferred; `retry` is dropped (static inputs: `first_of` over `.on(…)`; run-time inputs: `any_solver` ranges); TOMS748, ITP, Halley and Steffensen form an optional late phase; arclength continuation is dropped; algorithm ids are an open enum.

**Headline decisions**
1. **Backbone.** The public API is a value-and-composition model: of the candidate models, it is the only one that chains heterogeneous solvers without adapters. A pragmatic engineering discipline (plain loops inside steps, measured guardrails) governs kernels, scope and the one-call facade. A types-first vocabulary supplies the input types and the error design. The build, test and migration design was verified in separate build experiments. The prototype shows that the whole core is feasible on all 9 configurations **[prototyped]**.
2. **Solvers are values.** For example, `roots::brent{nxx::width_tol{1e-12}}.with_budget(60)`. They are called as `solver(f, input)`, curried as `solver.on(input)`, and driven by `nxx::iterate`, which returns the best iterate on every exit.
3. **Results.** `std::expected<solution<Est>, failure<Est, UE>>`. The failure is cheap to copy and holds no heap memory of its own: cause code, algorithm id, counters, best estimate, and the user's own callback error.
4. **Illegal states.** Invalid *literals* are compile errors (`bracket{2.0, 1.0}`, `max_iterations m = 0` or `= true`), and so are solver/input and criterion/solver mismatches. *Run-time* values are validated once, in-band (`make()` or the solver's `invalid_input`), and a validated value cannot become invalid.
5. **Scalars are open** (`nxx::scalar_traits<T>`, user-specialisable, defaulted from `std::numeric_limits`). Every default tolerance is an expression in `T`, so `float`, `long double` and multiprecision types work; multiprecision gets an optional adapter and CI leg.
6. **Linear algebra is Eigen 5.0.1**, as you asked: no BLAS/LAPACK, verified under Emscripten, fetched by CPM and wrapped in a thin `nxx::linalg` facade that returns `std::expected` and concrete types. Only the O(n) tridiagonal solvers for splines stay in-house.
7. **Order: dependencies first, then by breadth of use.** Core vocabulary → deriv → driver + 1-D roots (with `steps_view` and `any_solver`) → optimize → linalg + multiroots → poly → integrate + interpolate → optional roots → multiprecision, docs and release. Multidimensional minimisation and nonlinear least squares follow in v2.1, ODE initial-value solvers in v2.2 (§10.5).

### 1.1 Scope and positioning

**What Numerixx is.** A general-purpose, header-only C++23 numerical library in the spirit of GSL, with a smaller scope: the core numerical methods that most scientific and engineering code needs, done carefully, generically and composably. No application domain is a design driver. The author uses Numerixx for thermodynamics, among other things, but every feature in this document is justified on general grounds.

**Positioning.**

| Library | Licence and language | Model | Overlapping scope |
|---|---|---|---|
| GSL | GPL-3.0-or-later; C | double-only algorithms; mutable workspace objects (`alloc`/`set`/`iterate`/`free`); errors through return codes and a global error handler (by default, abort) | very broad: roots, minimisation, multiroots, quadrature, interpolation, ODEs, least squares, and much more (special functions, statistics, FFT, …) |
| Boost.Math tools | BSL-1.0; header-only C++ | generic scalar types; errors through policies (by default mostly exceptions) | roots (bisection, TOMS748, Newton, Halley, Schröder), minima (Brent), quadrature (Gauss, Gauss–Kronrod, tanh-sinh family, trapezoidal) |
| Eigen, with its unsupported `NonLinearOptimization` module | MPL-2.0; header-only C++. The unsupported nonlinear solvers are ports of MINPACK, and Eigen ships the Minpack licence (`COPYING.MINPACK`, BSD-style with an acknowledgement clause) for that code; its `LevenbergMarquardt` module states it explicitly. Numerixx uses this code only as a test oracle | expression templates; fixed and dynamic sizes | linear algebra; hybrid Powell (MINPACK `hybrd`/`hybrj`) and Levenberg–Marquardt for systems |

**Numerixx's niche:**
- value semantics and composable solvers: fallback chains, staging, warm restarts and shared evaluation budgets, static or assembled at run time (§6.10);
- `std::expected` results that carry the best estimate, the counters and the user's own callback error on every exit, with no exceptions and no global state (§3.4, §6.3);
- generic scalar types, from `float` to multiprecision, with defaults that are expressions in `T` (§3.5);
- functions as values: derivative, inverse, integral, antiderivative and interpolant as callables (§6.12);
- the MIT licence, which, unlike GSL's GPL, places no conditions on the licence of the code that uses it;
- Emscripten support in every exception mode (§4.5).

**Scope, mapped onto GSL's areas.**

| GSL area | Numerixx | Where, or what to use instead |
|---|---|---|
| Numerical differentiation | **v2.0** | `deriv` (§7.1); gradient, Jacobian and Hessian in `multiroots/derivatives.hpp` |
| One-dimensional root finding, bracket search | **v2.0** | `roots` (§7.2): bracketing and open methods; `expand`, `scan`, `subdivide` |
| One-dimensional minimisation | **v2.0** | `optimize` (§7.3): golden section, Brent, `bracket_minimum` |
| Polynomials (evaluation, algebra, roots) | **v2.0** | `poly` (§7.4): Horner, ring operations, closed forms, Aberth–Ehrlich |
| Multidimensional root finding | **v2.0** | `multiroots` (§7.5): damped Newton, Broyden, dogleg |
| Numerical integration (1-D): adaptive integration of general integrands over finite and infinite ranges, and fixed Gauss–Legendre | **v2.0** | `integrate` (§7.6): adaptive Gauss–Kronrod (QUADPACK QAG-style), tanh-sinh family (exp-sinh and sinh-sinh for infinite ranges), Romberg, Gauss–Legendre |
| Numerical integration (1-D): QUADPACK's weighted and oscillatory rules (QAWC, QAWS, QAWO, QAWF), QAGS's extrapolation, QAGP's user-supplied break points, CQUAD, and the fixed-point quadratures (fixed Gauss rules) other than Gauss–Legendre (Gauss–Laguerre, Gauss–Hermite, Gauss–Jacobi, …) | candidate | endpoint singularities: tanh-sinh; semi-infinite and infinite ranges: exp-sinh and sinh-sinh; known interior singularities: split the range there and sum; oscillatory integrands: split at the zeros and sum the pieces; QAGS-style extrapolation would reuse series acceleration (§10.5) |
| Interpolation (1-D) | **v2.0** | `interpolate` (§7.7): linear, cubic splines, monotone (PCHIP, Steffen), rational (Floater–Hormann) |
| 2-D interpolation (bilinear, bicubic on rectilinear grids) | candidate | a later extension of `interpolate`; it needs only core |
| Basis splines (B-splines) | candidate | revisit with linear least-squares fitting (§10.5): its basis is a user callable, so a B-spline basis plugs into it |
| Linear algebra | **v2.0, small** | `linalg` (§7.5): a dense LU/QR/Cholesky facade over Eigen |
| Multidimensional minimisation | planned, v2.1 | new module `multimin`: Nelder–Mead, BFGS/L-BFGS, nonlinear conjugate gradient (§10.5) |
| Nonlinear least-squares fitting | planned, v2.1 | new module `fit`: Levenberg–Marquardt (geodesic acceleration optional) on the Eigen QR facade (§10.5) |
| Ordinary differential equations (initial-value problems) | planned, v2.2 | new module `ode`: Dormand–Prince RK45 with dense output first; a stiff solver (Rosenbrock or BDF) later (§10.5) |
| Chebyshev approximation | planned, later | new module `chebyshev` (§10.5) |
| Series acceleration | planned, later | new module `series`: Richardson, Wynn epsilon, Levin u (§10.5) |
| Linear least-squares fitting | planned, later | in `fit`, on `qr_solve` (§10.5) |
| Special functions | out of scope | Boost.Math (portable); the C++17 special maths in `<cmath>` where the standard library provides it (libstdc++ and MSVC's library do; libc++, and so Emscripten, only partially) |
| Random numbers, quasi-random sequences, distributions | out of scope | `<random>`; specialised libraries |
| Statistics (including running and moving-window statistics), histograms, N-tuples | out of scope | specialised libraries |
| FFT, digital filtering | out of scope | specialised signal-processing libraries |
| Sparse matrices and sparse linear algebra | out of scope | Eigen's sparse modules |
| Eigensystems | out of scope beyond what Eigen provides | Eigen's eigensolvers |
| Monte Carlo integration, simulated annealing | out of scope | specialised libraries |
| Wavelet transforms, discrete Hankel transforms | out of scope | specialised libraries |
| Physical constants | out of scope | `std::numbers` for mathematical constants only |
| Vectors and matrices, BLAS, permutations, combinations and multisets, sorting, complex numbers, elementary maths functions, IEEE utilities | not needed | the standard library (`std::ranges`, `std::next_permutation`, `std::complex`, `<cmath>`, `<limits>`) and Eigen |

Status legend: **v2.0** = in the first release; **planned** = scheduled in §10.5; **candidate** = not scheduled, and can be added after v2.0 without changing the core; **out of scope** = use another library; **not needed** = the standard library or Eigen already covers it.

**How the core extends to the planned families [sketch].** The core abstractions (§6) were chosen so that the planned families reuse them instead of adding parallel machinery. A **minimiser** over ℝⁿ is a solver whose state holds an N-D point (plus, for BFGS, an inverse-Hessian approximation) and whose estimate is extremum-like: x, f(x) and a gradient norm. **Nonlinear least squares** reuses the systems machinery of `multiroots` with a residual vector r: ℝⁿ → ℝᵐ (m ≥ n): the FD Jacobian, the `typical` and `project` hooks, and `qr_solve`. An **ODE integrator** is a solver whose step advances t with an error-controlled step size; its stop criteria include reaching t_end; `steps_view` yields the trajectory; and dense output is a function-returning API that gives the solution as a callable t → `expected<y, fault>`. The driver, the criteria algebra, the `solution`/`failure`/`fault` types, cost accounting and the combinators apply unchanged, for example `first_of(nonstiff, stiff)` for an integration, or `with_evaluation_budget` over a minimisation. Each family is a new module placed downstream in the module DAG (§5.2), so no v2.0 module gains a dependency.

**FLAG, licensing.** GSL is GPL-3.0-or-later and Numerixx is MIT, so no GSL source may be ported or paraphrased. Algorithms are implemented from the literature (for example Brent 1973; Alefeld, Potra and Shi 1995; the MINPACK and QUADPACK reports; Dormand and Prince 1980; Nocedal and Wright), and each header cites its references. Values from GSL's documentation may be used as reference facts, but the test suite's oracles are Boost.Math, Eigen and high-precision reference tables (§9.2). Code derived from Boost (Brent, TOMS748) keeps its BSL-1.0 notice (§10.1).

### 1.2 Flags: where your wishes need nuance

- **FP.** Purity holds at function boundaries. Inside a step (and in Horner's rule or a tridiagonal solve) the code is an ordinary loop on locals. Chains are templates by default: zero overhead and usable in `constexpr`. Run-time chains are opt-in through `nxx::any_solver` **[prototyped]**, which is built on `std::function` because `std::move_only_function` is missing on libc++ 22 and hence on Emscripten. It may allocate when a solver is wrapped or copied (never on a call with an existing `F`), costs an indirect call per alternative and (with a `std::function` callable) per evaluation, is 24–64 bytes depending on the standard library, and is not `constexpr`. It wraps only copyable solvers, and a move is a copy, so that it can never be empty. `first_of` needs alternatives with the same result type; `.on(input)` is what unifies different inputs.
- **FXT.** Numerixx 2's internals use only `std::expected` members. FXT is the user-facing pipe vocabulary (`<numerixx/pipes.hpp>`, target `numerixx::pipes`) and, once FXT-4..7 land, the home of the generic parts of `first_of`/`then`. Today FXT lacks every numerics combinator (first-success, bounded iteration, Kleisli, refined types) and breaks `-fno-exceptions` (`throw 0;`; the fix is FXT-1).
- **Functions as return values.** They recompute on every call (no memoisation), can fail (they return `expected`), capture by value (`std::ref` for shared state), and are assignable even when they hold a capturing lambda.
- **Immutability.** It is enforced by the interface (private members, `with_*`, no mutators), **not** by `const` members, which break assignment, `std::expected` and the driver. Owning types deep-copy. There is **no exception**: N-D solvers carry value states over Eigen vectors, fixed-size when N is known at compile time and dynamic otherwise. A dynamic state allocates on every step; that is negligible next to evaluating any nontrivial system, and allocation is allowed.
- **Illegal states.** Compile-time rejection covers literals (tests, examples, constants). In real applications nearly every bracket, tolerance and budget is a run-time value (read from input, computed from data, or set by an outer algorithm), so in practice the guarantee is "validated once at the boundary, never invalid afterwards", plus in-band numerical errors. Numerical failures (NaN mid-iteration, a one-ulp bracket, a singular Jacobian, a stall, an exhausted budget) are not states a type can exclude.
- **Exceptions.** The library never throws, but user callbacks may (third-party code, or a model that reports failure by throwing). It is exception-neutral with conditional `noexcept`. `bad_alloc` from owning types (polynomials, interpolants, Eigen dynamic storage, `any_solver`) is not converted. Programmer errors are assertions, not `expected`. `-fno-exceptions` is a supported, CI-tested mode: it proves that the library never throws, serves codebases built without exceptions, and shrinks Emscripten builds.
- **Boost.** "Fetch Boost via CPM" becomes **no Boost in the library at all**; today Boost serves one trait, a stack trace and a 2-D table. Standalone boostorg `config` + `multiprecision` (+ `math`) come via CPM only for the optional multiprecision adapter and the test oracles, never under the package name `Boost`.
- **Linear algebra.** Eigen, as you asked, behind a facade. Eigen 5.0.1 has no BLAS/LAPACK dependency and was verified under em++ 6.0.8 (including `-fno-exceptions`) and with Boost `cpp_bin_float_50` scalars. Two costs are handled explicitly:
  - **compile time**, measured at +2.5–6.4 s per TU once LU and N-D Newton are instantiated **[prototyped]**. It stays inside TUs that include linalg or multiroots: scalar modules and the umbrella header never include Eigen, and a user who needs only the scalar modules sets `NUMERIXX_WITH_LINALG=OFF` and never downloads it;
  - **expression templates dangle under `auto`**, which clashes with an `auto`-heavy functional style. The facade returns concrete types only.

  Eigen is not `constexpr`, so N-D solves are not either.
- **Complex numbers** leave the root solvers and live only in `poly`, where complex roots arise routinely; complex open methods for general analytic functions are a specialised need that can be added later without changing the core (open decision 14).

---

## 2. Decisions

| # | Decision | Choice | Alternatives considered | Rationale |
|---|---|---|---|---|
| D1 | Starting point | Fresh tree on the redesign branch. Port algorithm bodies by hand, mostly from dev-reorg. Tag `v1.0.0` → 5de1e07 and `v1.1.0-legacy` → 8528e94. Never merge dev-reorg. | Merge dev-reorg; continue from master | Every file is rewritten anyway. A merge would put 38.8 MB of Blaze into history for good, and dev-reorg carries regressions (§7). |
| D2 | Language and toolchains | C++23. Floors: GCC 14, Clang 18 + libc++, MSVC with `/std:c++latest` (a toolset with deducing `this` and `std::expected`'s monadic operations), clang-cl, em++ ≥ 6.0.8. Tested: GCC 16, Clang 22, MSVC 19.51, clang-cl 22, em++ 6.0.8. A nightly floor-compiler job confirms the floors; CI runs the newest toolchains. **MSVC lacks P2564** (`__cpp_consteval` = 201811L), so user code must forward refined types, never raw scalars, into refined constructors (§6.2). | C++20 | FXT needs C++23. The design uses `std::expected` and deducing `this`, which set the GCC 14 / Clang 18 floor (FXT's README states GCC 13, Clang 17, MSVC 19.30). |
| D3 | Solver model | Immutable, configured solver **values** with `init`/`step`. `.on(input)` curries a solver into `f → result`. **[prototyped]** | Problem types that bind f; stateless function objects with per-call criteria | Only this model lets `first_of(newton…, secant…, brent…)` mix inputs without adapters. |
| D4 | Call order | `solver(f, input)`; curried form `solver.on(input)(f)`; `std::bind_front(solver, f)` fixes f **[sketch]** | `(input, f)` | Matches the facade, the old `fsolve(f, bounds)` and every candidate API but one. |
| D5 | Driver | One `nxx::iterate(alg, problem[, observer])`. Successes carry `estimate(s)`; failures carry the tracked `best(s)`, which may be a different type (searchers). A driver-internal `nxx::detail::advance` maps a step fault to a failure. There is an optional `finish()` post-condition. **[prototyped; `advance` and `finish` are sketch]** | Five loops (today); a `Next{state, done}` step result | Fixes the off-by-one and success-at-maxiter bugs once. The best iterate on every exit is what callers act on (R-A3). |
| D6 | Stop criteria | First-class values combined with `\|\|`/`&&`, returning a `verdict`, and **typed by view kind**: `x_tol` and `step_tol` judge successive iterates (open and N-D methods only); `width_tol` and `floored_width` judge enclosures (bracketing only). A mismatch is a compile error with a reason. Solvers with an internal convergence test (Brent, golden, Brent-min; TOMS748 and ITP when added) take a width criterion *as their tolerance*; their external stop defaults to `never{}`. The budget is a separate, mandatory field. | One criterion type for all views | With one criterion type, `x_tol` on a bracketing view compares repeated endpoints; the prototype reported success at x = 1.5 **[prototyped]**. |
| D7 | Result type | `std::expected<solution<Est>, failure<Est, UE>>`. `failure` holds `errc`, `algo`, `counters`, `optional<Est> best` and `cause<UE>`, with no strings, stack traces or `source_location`. **Guideline: errors are cheap to copy and hold no heap memory of their own** (the best estimate of a dynamically sized system owns its vectors). Sizes are recorded, not gated. `result<>` is **not** trivially copyable on libc++ **[prototyped]**. | Error codes of ≤ 16 bytes without an estimate; `source_location` in errors | One type per problem family makes chains type-check. The best estimate and counters are what callers act on (R-A3, `warm_fallback`). |
| D8 | Error codes | `errc` in ranges: input 1–31, numerical 32–63, callback 64+ **[prototyped]** | One flat enum | Cheap classification; stable numeric codes. |
| D9 | Error categories | Numerical failure → `expected`. Contract violation → `NXX_EXPECTS` (assert; C++26 `pre` later). Foreign exceptions propagate. | Everything in `expected` | Different owners, different fixes. |
| D10 | Exceptions | The library never throws. `noexcept` is conditional on the callbacks. Exception-neutral. Opt-in `nxx::fn::catching(f)` under `__cpp_exceptions`. `-fno-exceptions` is a supported, CI-tested mode. | Blanket `noexcept` | A blanket `noexcept` would `std::terminate` the program whenever a user callback throws, for example third-party code inside f (R-E2). |
| D11 | Illegal states | Three tiers: compile time, construction time, run time. `consteval` literal constructors plus `make() → expected` **[prototyped]**. `detail::refined<Tag, T>` for scalar refinements; hand-written classes for multi-field invariants. | A generic `fxt::refined` now | One mechanism with readable diagnostics. Can be upstreamed later. |
| D12 | Bracket types | `bracket<T>` (finite, lo < hi) and `sign_bracket<T>` (bracket + endpoint samples of opposite sign, or a zero; samples may be ±inf). Not bound to f's type. **[prototyped]** | `Bracketed<F,T>` bound to the closure type | Binding to the closure type breaks heterogeneous chains and wrapped functions. |
| D13 | Capability errors | A solver/input mismatch does not compile; family facades carry a reasoned deletion. **Newton needs an explicit derivative source**, which is one of: (a) a callable `df`; (b) a **derivative policy** with `bind(f)`, e.g. `deriv::numeric{}`; (c) a structural `.derivative()` on f. **[prototyped]** | A silent finite-difference fallback | An implicit FD silently changes cost and accuracy. The policy keeps it explicit and works in curried chains without a roots→deriv edge **[prototyped]**. |
| D14 | Scalars | Open `nxx::scalar_traits<T>` (user-specialisable; the default covers every type whose `numeric_limits` is specialised, non-integer and inexact); concept `real`; maths through the ADL idiom (`using std::abs; abs(x)`) or thin `nxx::math` helpers; no `double` round-trips; **every default tolerance is an expression in `T`** (`root_eps<T>`), with no run-time `pow` **[prototyped mechanism]** | `IsFloat` + Boost | `float`, `long double` and multiprecision work without an adapter; double-literal defaults are unattainable in `float`; power-of-two constants are exact on every platform. |
| D15 | Complex numbers | Only in `poly` | Complex open methods | Complex roots arise routinely for polynomials; complex open methods for general analytic functions are a specialised need and can be added later without changing the core. |
| D16 | Multiprecision | Works through the open trait for expression-template-off types. Optional `numerixx::multiprecision` adapter (its linalg header includes `<boost/multiprecision/eigen.hpp>`) and CI leg. | Drop the claim | Cheap once the core is Boost-free. |
| D17 | Composition | Named class templates `first_of_t`, `then_t`, `warm_fallback_t`, `bound`, `budgeted_t`; factories `first_of`, `first_of_with(policy, …)`, `then`, `warm_fallback`, `with_evaluation_budget`. User callables are stored in a copyable box. Written on `std::expected` members. **[prototyped mechanism; named types sketch]** | Closures (not assignable **[prototyped]**); `fxt::operator\|\|` (eager, same type); `std::function` everywhere | No allocation, constexpr, clang-cl-safe, assignable, readable symbols. |
| D18 | Where FXT is included | Kernels and combinators include only the standard library. `<numerixx/pipes.hpp>` is the only FXT include, in its own target `numerixx::pipes`; `numerixx::core` does not link FXT. **[prototyped]** | FXT inside the facade | The core and its `-fno-exceptions` legs do not wait for FXT-1. A user who does not want the pipes need not fetch FXT (`NUMERIXX_WITH_FXT=OFF`). One include to audit. |
| D19 | Manual stepping and tracing | **`steps_view`**, a lazy input range of `expected<state, fault>`, is *the* public manual-stepping and tracing API, delivered with the roots phase **[prototyped]**. The observer handles logging inside `iterate`. `init`/`step` remain the solver protocol; `advance` is driver-internal. | `std::generator` (missing on libc++ 22); a public `init`/`step` + `advance` loop | One public way to step; composes with `std::views`; a hand-written loop does not even type-check (`init` and `step` return different error types) **[prototyped]**. |
| D20 | Per-iterate projection | `.with_projection(p)` on open and N-D methods, applied to the proposed iterate *before* evaluation. An iterate pinned at the edge → `stalled`. A final projection is a `transform`. **[prototyped]** | Extending the objective beyond its domain (`fn::extend_linearly`) | Box and domain constraints (bounds from the application, or a mathematical domain such as x > 0) become a library feature instead of a change to the user's function, and f is never evaluated outside its domain (R-A8). |
| D21 | Functions as values | `derivative_of`, `inverse_of`, `integral_of`, `antiderivative`, interpolants, `minimizer_of`, `poly::derivative`. Named `_fn` types that capture by value and are **copy-assignable whenever the captured callables are copy-constructible** (semiregular when those are also default-constructible). Two return conventions (§6.4): **scalar-valued** callables (`derivative_of`, `second_derivative_of`, `reject_outside` interpolants) **return `expected<T, fault<UE>>`**, so they plug into solvers as callbacks without nesting; **estimate-valued** callables (`integral_of`, `antiderivative`, `inverse_of`, `minimizer_of`) return `result<Est, UE>` (value, error estimate, counters) and become callbacks through `fn::value_of` **[sketch]**. | `ResultProxy` (dev-reorg); a `failure`-returning `derivative_of` | Fixes today's broken `derivativeOf`/`integralOf`; a `failure` return would nest errors inside Newton's **[prototyped]**. |
| D22 | Linear algebra | **Eigen 5.0.1 via CPM is the backend** of `numerixx::linalg` and `numerixx::multiroots`, behind a thin `nxx::linalg` facade: `lu_solve`, `qr_solve` (`ColPivHouseholderQR`, rank-deficient and least squares) and `cholesky_solve` return `std::expected<…, errc>`; they check dimensions and `allFinite()`, `lu_solve` also checks `rcond` against n·ε (`PartialPivLU` never reports singularity itself), and all return concrete types only. Aliases `nxx::linalg::vector<T, N = Eigen::Dynamic>` and `matrix<T, R, C>`. Systems accept `std::array`, `std::vector` and Eigen vectors. In-house code is limited to the O(n) tridiagonal solvers inside `interpolate`. Option `NUMERIXX_WITH_LINALG` (default ON). | In-house LU/Cholesky/QR kernels over views (constexpr, allocation-free); keep Blaze (needs BLAS/LAPACK); Eigen only as an optional adapter | Your request: a BLAS/LAPACK-free library that compiles with Emscripten. Eigen is mature, handles fixed and dynamic sizes, and is scalar-generic (verified with `cpp_bin_float_50`). Measured cost +2.5–6.4 s per TU in use **[prototyped]**, confined to linalg/multiroots TUs. |
| D23 | Naming | snake_case for public names. Namespace = header = CMake target per module. The integrate facade is `integrate::quad` (not `integrate::integrate`). | PascalCase; `nxx::optim` | Uniform pipelines; one name per module; no namespace/function clash. |
| D24 | Headers and targets | A single include root `<numerixx/…>`. Thin per-module INTERFACE targets with **SYSTEM** include directories for consumers. The DAG is enforced by a layering test. | Per-module include roots | Fixes case and name collisions; consumers with `/W4 /WX` never see Numerixx warnings. |
| D25 | CMake and CPM | CMake ≥ 3.25 (3.30+ recommended on Windows). CPM 0.43.2 committed, with a wrapper that reuses a parent's CPM. `set(CMAKE_POLICY_DEFAULT_CMP0168 NEW)` before CPM. | vcpkg; plain FetchContent | Verified: plain `cmake_policy(SET …)` does not reach CPM. |
| D26 | FXT fetch | `NAME FXT`, pinned SHA + SHA256, `DOWNLOAD_ONLY`, own `fxt::fxt` shim only if no parent provides it; `-DCPM_FXT_SOURCE=` override. Fetched when `NUMERIXX_WITH_FXT` (default ON). | Running FXT's own CMake | FXT's CMake downloads unhashed CPM and TL repositories that fail under MAX_PATH. |
| D27 | Boost | None in the library. Standalone boostorg `config` + `multiprecision` (+ `math` for oracles) 1.92.0 under non-`Boost` package names, OFF by default. | Full Boost; vcpkg | The library needs none, so no user pays for Boost headers in every TU (measured: in Numerixx 1.x, the deriv header alone brought about 477 Boost headers into every TU that included it). Claiming the common package name `Boost` would clash with a parent project that declares its own (possibly partial) Boost, where the first declaration wins (R-B4). |
| D28 | Tests | doctest 2.5.3 via CPM (the author's choice). Corpus + property + compile-fail (control twins, reason-string checks) + `static_assert` suites; header self-containment; layering; EH-mode legs; scalar matrix (`float`, `double`, `long double`, `cpp_bin_float_50`); determinism; oracles (Boost.Math, Eigen `HybridNonLinearSolver`). | Catch2 3.16 via CPM; FetchContent Catch2 3.4 (today) | Only 2 of 4 test files compile today. doctest is a single header with a small compile-time footprint and detects `-fno-exceptions` by itself. |
| D29 | Defaults | `roots::solve(f, bracket)`: the corpus-chosen bracketing solver (Brent provisionally). `solve(f, x0)`: `then(expand, brent)`. `solve(f, df, x0)`: `then(expand, rtsafe)`. `optimize::minimize`: Brent-min. `integrate::quad`: adaptive G7K15. `deriv::central`: optimal relative step. `multiroots::solve`: dogleg with Broyden updates once it passes the phase-5 corpus; until then damped Newton with an FD Jacobian. | TOMS748 by fiat | Choose by mean and worst-case evaluation counts on the §9.2 corpus, recorded in phase 3 and re-run in phase 8 when TOMS748 and ITP arrive. |
| D30 | Polynomial roots | Deterministic Aberth–Ehrlich with polishing of every root; sum types for the closed forms; no dependency on `roots` | Laguerre with `random_device` | Reproducible results. Breaks the roots↔poly cycle. |
| D31 | Migration | Tags `v1.0.0` and `v1.1.0-legacy` let existing users pin the old API; `MIGRATION.md` maps old calls to new ones; no compatibility shim (§10.4). | A compatibility shim | A shim would reproduce removed bugs: unchecked results, success at maxiter, the O(1)-wrong default second and mixed derivatives. |
| D32 | Scale | One `typical` scale concept threads through defaults. Derivative steps are **relative** (typical defaults to \|x\|, and to 1 at x = 0). Stopping tolerances keep an absolute floor at `scale = 1`, visible and configurable, so roots at 0 terminate. | Hard-coded `max(1, \|x\|)` everywhere | A floor-1 step is 46 % of x at x = 1e-3 (a 27 % error in the derivative of 1/x, measured). A purely relative stopping test never terminates at a root of 0. |
| D33 | Evaluation accounting | `counters.evaluations` counts calls of the user's f. `cost_of(fn)` (e.g. a `derivative_fn` costs its stencil size); `fault` carries the failing step's evaluations; stages seeded with (x, fx) or a sampled enclosure do not re-evaluate. | 1 per callable call | Each evaluation can be expensive: a simulation, an inner iterative solve, a table look-up (R-A5). Counting 1 per callable call under-counts by 33 % with an FD derivative (measured). |
| D34 | Algorithm ids | `enum class algo : std::uint8_t` with no closed enumerator list. Each module defines constants in reserved ranges; user solvers use ≥ 200. | A closed core enum | A closed list is a layering inversion and contradicts the open solver set. |
| D35 | Run-time chains | Opt-in `nxx::any_solver<F, Est, UE = none>` (header `<numerixx/core/any_solver.hpp>`) type-erases a curried solver `const F& → result<Est, UE>` for a fixed callable type `F` (e.g. `std::function<double(double)>`). It is copyable and never empty (no default constructor, no move operations). `nxx::first_of` also accepts a run-time range of them and returns an `any_solver`. Static template chains stay the default. **[prototyped: all 9 configurations, bit-identical to the static chains; the policy overload is sketch]** | Static chains only; `std::move_only_function` | Allocation is allowed, and chains assembled from configuration are a real need. `std::move_only_function` is missing on libc++ 22/Emscripten, so it is `std::function`. |

---

## 3. Design principles

### 3.1 What "functional programming style" means here

1. **Values in, values out.** Problems, solvers, criteria, estimates, results and errors are copyable (and copy-assignable) values with no mutating API.
2. **Pure steps.** `init(problem) → expected<state, failure>` and `step(problem, state) → expected<state, fault>` are pure except for calling the user's function. The driver rebinds a *local* state variable.
3. **One loop.** Every iterative algorithm runs through `nxx::iterate`. Manual stepping goes through `steps_view`, which walks the same `init`/`step` protocol lazily.
4. **First-class solvers.** A configured solver is a value. `.on(input)` curries it to `f → result`. Combinators build new solver values from old ones; `any_solver` erases a curried solver's type when a chain must be assembled at run time.
5. **Functions as return values** wherever the mathematics is "an operator on functions": derivative, inverse, integral, antiderivative, interpolant, minimiser-of-a-family.
6. **Monads at the boundary only.** `std::expected` appears once per evaluation, once per step and once per result, never per arithmetic operation. An exploratory measurement put the whole abstraction at 1.18–1.22× a hand-written bisection loop on the cheapest possible f; with any nontrivial f (≥ 10 µs per evaluation, for example one that runs an inner iterative solve), that is noise.
7. **FXT at the edge.** Users compose results with FXT pipes (`transform`, `and_then`, `or_else`, `match`, `tap`, `value_or`) **[prototyped]**. The library's own code uses the members of `std::expected`.

### 3.2 Immutability discipline

- **No `const` data members, anywhere.** They delete copy and move assignment.
- Types with an invariant have **private members, `constexpr` getters and `with_*()` builders** that return a new value. If a change could break the invariant, `with_*` returns `expected`.
- **Records are public aggregates**: `counters`, `solution`, `failure`, `fault`, `extremum`, `integral`, `derivative_estimate`, `system_estimate`. `root_estimate` is a public struct too, but its constructor requires x and f(x) (§7.2), so it is not an aggregate **[spike]**.
- **Combinators and function-returning APIs are named class templates**, not closures, because closure copy assignment is deleted **[prototyped: the prototype's closure-based chain, `then` and capturing `derivative_of` were not copy-assignable]**. User callables are held in `detail::copyable_box<F>`, the `std::ranges` movable-box technique: assignment is implemented as destroy + construct, and `[[no_unique_address]]` applies for empty F. A `static_assert` suite, started in phase 1 and extended by every module phase, checks `std::copyable` for every solver, chain, `any_solver` and returned function, and `std::semiregular` when the parts are default-constructible (`any_solver` is deliberately not default-constructible, §6.10).
- **No `mutable` members, no lazy caches, no statics.** Derived data (spline coefficients, PCHIP slopes, Gauss nodes) is computed eagerly in the smart constructor. The reasons are thread safety (dev-reorg's `mutable std::optional` caches are data races) and reproducibility (bit-identical results on repeated calls).
- **Owning values deep-copy** (`polynomial`, interpolants, Eigen dynamic vectors, `any_solver`).
- **No exception to the rule.** N-D solvers use value states over Eigen vectors: `linalg::vector<T, N>` when N is known at compile time (`std::array` or fixed-size Eigen input), dynamic otherwise. Each step builds a new state, so a dynamic state allocates a few times per step. An exploratory measurement of value states over `std::vector` found 5.2 allocations per step and 1.7–1.8× the time of an in-place kernel on a cheap f; next to the evaluation of any nontrivial system, that is negligible.
- **User callables cannot be made pure.** They often carry state: a model object that caches intermediate results, an instrumentation counter, or a wrapper around a stateful third-party library. Numerixx guarantees sequential, deterministic, `const`-invoked calls and no gratuitous copies. One-shot calls hold f by `std::cref`. Returned callables capture by value; stateful objects go in through `std::ref`.

### 3.3 Illegal-state strategy: type level versus run time

| Tier | Mechanism | Examples | Run-time cost |
|---|---|---|---|
| **A. Compile time** | types, concepts, `consteval` literal constructors, deleted overloads with reasons | Invalid literals: `bracket{2.0, 1.0}`, `tolerance{-1e-8}`, `max_iterations m = 0` and `= true`, `x_tol{0.0, 0.0}`. Solver/input mismatches: bisection on a guess; Newton without a derivative source; a fixed-size guess of the wrong length; `quadratic{0.0, 1.0, 2.0}`. Criterion/solver mismatches: **`x_tol` or `step_tol` on a bracketing solver**; `width_tol` on an open method; **`f_tol` on a minimiser**. Composition errors: `first_of` over mismatched result types; `then(newton, brent)`, where a bracket solver cannot take a point estimate without `.from_enclosure()`; an `any_solver` built from a solver whose result type differs. Roles: `rel_tolerance` where `tolerance` is expected. **[prototyped: 16 compile-fail tests on 5 compilers]** | none |
| **B. Construction time** | private constructors + `make() → std::expected<T, errc>` | a tolerance or budget from a config file; a bracket from run-time values; `sign_bracket::make(f, b)`; strictly increasing knots; finite polynomial coefficients | one check, once, at the boundary |
| **C. Run time** | the error channel of every solver; unvalidated inputs accepted by solvers | braced `{lo, hi}`, `std::pair` and `expected<In, errc>` inputs validated in `prepare()`; NaN from the callback; the callback's own error; zero derivative; singular Jacobian; stall; cycling; divergence; budget exhausted; a sign change that is a pole; a run-time dimension mismatch between a system and its guess | once per evaluation or step |

Rules:
- **Normalise when every input has a canonical legal representative**: `polynomial{1.0, 2.0, 0.0}` becomes degree 1. **Validate when some inputs have none**: a reversed bracket from `make` is re-ordered, but equal endpoints are an error.
- **Iteration states change only through invariant-preserving transitions, where that is cheap.** `sign_bracket::narrowed(m, fm)` keeps the half that still changes sign. Heavier states are plain aggregates, reachable only through manual stepping.
- **Solver values never hold an invalid configuration.** Builders take only validated types. A run-time budget or tolerance goes through `make()` once, when configuration is parsed.
- **What stays at run time is documented:** non-finite values, callback errors and exceptions, singularity, stall, cycling, divergence, budget, poles and jump discontinuities (§7.2), a noisy f contradicting cached end signs, run-time dimension mismatch, wrong user derivatives (test helper `check_derivative`), allocation failure, and internal bugs (`NXX_ASSERT`).

### 3.4 Error model

- **`std::expected<V, E>`** for every operation that can fail for a reason: solvers, smart constructors, evaluations, linear solves, function-returning APIs, and interpolant evaluation outside the knots.
- **`std::optional<T>`** only where absence *is* the answer: `polynomial::degree()`, `failure::best` (nothing evaluated yet), `root_estimate::enclosure`, `system_estimate::rcond`, narrowing conversions, `best_x(r)`.
- **Errors are cheap to copy and hold no heap memory of their own**: no strings, `std::stacktrace`, `exception_ptr` or `source_location`. The best estimate is the family's estimate type, so a failure of a dynamically sized system owns its vectors.
- **Contracts.** Violated preconditions of the low-level protocol (calling `init`/`step` with a problem that `prepare()` did not produce, `narrowed()` with a point outside the bracket) are `NXX_EXPECTS(cond)`. Run-time dimension mismatches in linalg and multiroots return `errc::dimension_mismatch` instead.
- **Exceptions from user code** propagate untouched. `noexcept(std::is_nothrow_invocable_v<…>)` everywhere. Under `__cpp_exceptions`, `nxx::fn::catching(f)` turns a throwing callback into a fallible one.
- **Fallible callbacks** `x → std::expected<T, E>` are first-class. The failure's code is `callback_failed`, and `failure::cause` holds the user's `E` unchanged **[prototyped]**. Plain callbacks carry `none`, which takes no space. A user can mark an error fatal for fallback chains with the CPO `nxx::is_fatal(e)` (default `false`).
- **No silent failure.** Returning the last iterate, NaN, an endpoint or a pole as success is forbidden. **A success with `stop_reason::criterion` implies that the criterion's guarantee holds for the returned estimate** (§9.3). Every failure after the first evaluation carries the best estimate.

### 3.5 Genericity

- **Real scalars:** `float`, `double` and `long double` in every test; `cpp_bin_float_50` in the multiprecision leg.
- **Maths calls** use the ADL idiom (`using std::abs; abs(x)`) or the thin `nxx::math` helpers that wrap it, so the functions of a multiprecision type, which live in its own namespace, are found.
- **Defaults are functions of `T`:** `rel = math::root_eps<T>(1, 2)` and similar. A `static_assert` test instantiates every module default for `float`, `double`, `long double` and `cpp_bin_float_50` and checks that it is achievable (`default_rel >= 4·eps_T`).
- **Multiprecision:** `cpp_bin_float_50` satisfies `real` with no adapter. Only expression-template-off types are supported: algorithms write `T x = …`, never `auto x = a - b;`. The same rule protects against Eigen's expression templates (§5.3).
- **AD scalars are not a design goal** and are not tested. A user type that satisfies `real` (with a `scalar_traits` specialisation if it has no `numeric_limits`) is accepted, but control flow compares values of `T` directly, and nothing guarantees convergence of derivative parts.
- **Complex:** `poly` only.
- **Domain versus codomain** are separate type parameters.
- **No double literals inside algorithms**, except tabulated rule constants (Gauss–Kronrod nodes and weights, ≥ 36 significant digits; §7.6). Stencils are integer weights over a common denominator. Constants are `T(k)` or powers of two.

### 3.6 Performance guardrails

- **Allocation is allowed.** The 1-D solvers, criteria and static chains hold only scalars and small arrays, so they allocate nothing and run in constant expressions **[prototyped]**. Dynamically sized N-D solvers, polynomials, interpolants and `any_solver` may allocate (the prototype's allocation counts for `any_solver` are in §6.10). The library's tests do not count allocations.
- **Evaluation economy** (D33):
  - endpoint values cached in `sign_bracket`, so bisection costs 1 evaluation per step and Ridders 2;
  - seeded stages skip re-evaluation;
  - `inverse_of` caches its endpoint samples;
  - benchmarks report evaluations alongside time;
  - efficiency index per evaluation: secant 1.618; Newton with a same-cost analytic `df` 1.414; Newton with a central-difference `df` 1.26 (dominated by secant; documented at `deriv::numeric`).
- **Algorithm choice beats micro-optimisation.** Brent needs 7 iterations and 9 evaluations on x²−2 over [1, 2] **[prototyped]**, where bisection needs about 50.
- **State size.**
  - 1-D states are a few scalars (Brent's is 72 bytes).
  - Damped Newton's Jacobian and LU factors are step temporaries.
  - **Broyden and dogleg states carry O(n²) data** (a Jacobian or its QR factors, the trust radius), so each step builds a new O(n²) state. That is acceptable next to O(n³) factorisations and expensive evaluations of F, but documented.
  - Fixed-size Eigen storage lives on the stack and suits small systems; larger systems use dynamic size.
  - `std::vector` systems are converted to Eigen storage at the evaluation boundary, which costs an allocation per evaluation **[sketch]**.
- **Compile time (seconds per TU, best of 3, `-c`, -O2; shared machine, so noisy) [prototyped]:**

  | TU | GCC | Clang | em++ | MSVC | clang-cl |
  |---|---|---|---|---|---|
  | standard-library baseline | 0.61 | 0.51 | 0.75 | 0.39 | 0.38 |
  | core + 15 constexpr solves | 1.36 | 1.28 | 1.37 | 1.16 | 0.95 |
  | core + FXT pipes | 1.28 | 1.55 | 1.37 | 1.24 | 1.01 |
  | N-D Newton, in-house LU (prototype baseline) | 0.68 | 0.55 | 0.71 | 0.65 | 0.81 |
  | N-D Newton, Eigen (fixed + dynamic) | 5.46 | 3.45 | 3.54 | 7.09 | 3.31 |
  | `<Eigen/Core>` + `<Eigen/LU>`, include only | 0.96 | 0.74 | 0.95 | 1.40 | 0.82 |

  FXT pipes cost nothing measurable. The design uses Eigen, so a TU that instantiates linalg or multiroots pays the Eigen row; the in-house row shows what that choice costs. Scalar modules never include Eigen, and neither does the umbrella header (§5.2). CI warns when a TU that includes `<numerixx/numerixx.hpp>` and instantiates nothing exceeds 2 s on GCC; the instantiation cost of one linalg/multiroots TU is recorded per release, not gated.
- **Diagnostics** (claims limited to what tests enforce, §9.1):
  - Contract violations of combinators and solver inputs produce a `static_assert` or a reasoned deletion whose message is in the **first** error on GCC and Clang.
  - The `first_of` mismatch is a single error line on GCC and MSVC **[prototyped]**.
  - MSVC does not print deletion reasons (it rejects `= delete("…")`); it shows the deleted declaration, whose source line holds the reason.
  - Symbol length: the longest mangled name in the core TU is 381 characters, about 1.37k demangled **[prototyped]**.

---

## 4. Build system and dependencies

### 4.1 Principles (verified in build experiments)

1. The core needs only the standard library. FXT is needed only by `numerixx::pipes`, and Eigen only by `numerixx::linalg` and `numerixx::multiroots`. Everything else is optional (Boost) or development-only (doctest, benchmark, Doxygen).
2. **Be a good subproject.** Tests, examples, benchmarks and docs are OFF when not top-level. No global flags, no `export(PACKAGE)`, no REQUIRED `find_package`. Reuse any CPM, FXT, Eigen or Boost targets a parent provides. Never claim a common CPM package name such as `Boost` (D27). Tested with a CPM parent and a FetchContent parent, in both declaration orders (§9.4). Documented for parents:
   - declare Numerixx as `CPMAddPackage(NAME Numerixx …)`, so that `CPM_Numerixx_SOURCE` can replace it with a local checkout;
   - a parent that uses FXT itself declares it as `NAME FXT` with `FXT_USE_TL_EXPECTED` and `FXT_USE_TL_OPTIONAL` OFF (their defaults); CPM then deduplicates it with Numerixx's declaration, and Numerixx reuses the parent's `fxt::fxt` (D26).
3. **Pin everything** to an exact version or commit plus a SHA256. Local overrides use CPM's `CPM_<Name>_SOURCE`.
4. **Consumers get usage requirements only:** `cxx_std_23` and **SYSTEM** include directories (CMake ≥ 3.25 gives `/external:I` on MSVC). The library's own tests add the include root non-system, so its warnings stay visible to us. Warning, exception and sanitizer flags stay on Numerixx-owned targets. `NumerixxWarnings` adds `/Zc:__cplusplus` for MSVC-owned targets (cl reports 199711L otherwise **[prototyped]**).
5. **The exception model belongs to the consumer.** The library adds no EH flags. The test directory builds every EH mode.

### 4.2 Dependencies

| Dependency | Pin | How | Needed by | Default |
|---|---|---|---|---|
| CPM | 0.43.2, `get_cpm.cmake` SHA256 `49a3bef9…f232aa` | committed file + a 6-line wrapper that reuses a parent's CPM (tested with 0.42.1) | build | always |
| **FXT** | `9570f44d…`, which adds the FXT-1 probe fix (troldal/FXT#1; moves to the merge commit on FXT's main once merged); archive SHA256 | `CPMAddPackage(NAME FXT URL … URL_HASH … DOWNLOAD_ONLY YES)` + `fxt::fxt` shim if absent | `numerixx::pipes` only | `NUMERIXX_WITH_FXT=ON` |
| **Eigen** | 5.0.1, SHA256 `e9c326dc…3dec` (fallback `GIT_TAG 5.0.1`) | `DOWNLOAD_ONLY` + own INTERFACE target, unless a parent provides `Eigen3::Eigen` | `numerixx::linalg`, `numerixx::multiroots`; the `HybridNonLinearSolver` oracle tests | `NUMERIXX_WITH_LINALG=ON` |
| Boost.Config + Boost.Multiprecision | boost-1.92.0 (SHA256 `b4171037…a0`, `9da99784…01d0`) | standalone boostorg repos, `BOOST_MP_STANDALONE ON`; package names `boost_config`/`boost_multiprecision`; guarded by `if(NOT TARGET Boost::multiprecision)` | `numerixx::multiprecision` adapter | OFF |
| Boost.Math | boost-1.92.0 (SHA256 `aa84eec6…2c48`) | standalone, `BOOST_MATH_STANDALONE ON` | oracle tests only (roots, minima, quadrature; §9.1) | OFF (`NUMERIXX_TEST_ORACLES`) |
| doctest | v2.5.3, SHA256 `174ebc4e…1331` | URL + hash, fetched inside `tests/`; the test main is built by Numerixx (C++23, the tests' exception model) | tests | top-level only |
| google/benchmark | v1.9.5, SHA256 `9631341c…a340` | URL + hash | benchmarks | OFF |
| Doxygen / Sphinx / Breathe | system | `find_package` without REQUIRED | docs | OFF |

**Deleted:** `vcpkg.json` (and the `.idea` toolchain paths), gcem, tl-expected, Blaze, LAPACK, OpenBLAS, OpenMP, nlohmann-json, fmt, hwinfo, Boost.Stacktrace, Boost.MultiArray, and the 172-file vendored `benchmark/gbench`.

**FLAG, Boost.** Your request to "fetch Boost via CPM" is honoured in form but not in substance: the library needs no Boost. `IsFloat` becomes `scalar_traits`. Stacktrace is dropped: it has no Emscripten backend, is absent from libc++, and would be paid on every failed attempt in a fallback chain. MultiArray becomes two rows of an array. CPM fetches the standalone repositories (about 5 MB instead of 108 MB), and only when you opt in.

### 4.3 Options

| Option | Default | Meaning |
|---|---|---|
| `NUMERIXX_BUILD_TESTS` | `PROJECT_IS_TOP_LEVEL` | doctest suite, header check, compile-fail, layering |
| `NUMERIXX_BUILD_EXAMPLES` | OFF | examples, registered as smoke tests |
| `NUMERIXX_BUILD_BENCHMARKS` | OFF | ignored under Emscripten |
| `NUMERIXX_BUILD_DOCS` | OFF | skipped with a message if the tools are missing |
| `NUMERIXX_INSTALL` | `PROJECT_IS_TOP_LEVEL` | install/export; `find_package(numerixx 2.0 CONFIG)` |
| `NUMERIXX_WITH_FXT` | ON | fetch FXT; provide `numerixx::pipes` (set OFF when the pipes are not wanted; nothing else needs FXT) |
| `NUMERIXX_WITH_LINALG` | ON | fetch Eigen 5.0.1; provide `numerixx::linalg` and `numerixx::multiroots` and their oracle tests (a user of the scalar modules only can set it OFF and never download Eigen) |
| `NUMERIXX_WITH_MULTIPRECISION` | OFF | `numerixx::multiprecision` adapter and MP tests |
| `NUMERIXX_TEST_ORACLES` | OFF | Boost.Math oracles |
| `NUMERIXX_NO_EXCEPTIONS` | OFF | tests/examples with `-fno-exceptions` or `/EHs-c- /D_HAS_EXCEPTIONS=0` |
| `NUMERIXX_WARNINGS_AS_ERRORS` | OFF (ON in presets) | `-Werror` / `/WX` on owned targets |
| `NUMERIXX_SANITIZE` | "" | e.g. `address;undefined` (tests only) |
| `NUMERIXX_BUILD_INTEGRATION_TESTS` | OFF | register the consumer-build scenarios of §9.4 as CTest tests (`tests/integration`, preset `integration`); host builds only |

### 4.4 CMake sketch (condensed from verified build experiments; the linalg wiring is [sketch])

```cmake
# CMakeLists.txt
cmake_minimum_required(VERSION 3.25...4.4)
set(CMAKE_POLICY_DEFAULT_CMP0168 NEW)          # MUST be the *default*: CPM runs cmake_minimum_required(3.14) internally
project(Numerixx VERSION 2.0.0 LANGUAGES CXX)
option(NUMERIXX_BUILD_TESTS "Build tests" ${PROJECT_IS_TOP_LEVEL})
option(NUMERIXX_WITH_FXT "FXT pipes (numerixx::pipes)" ON)
option(NUMERIXX_WITH_LINALG "Eigen-backed numerixx::linalg + numerixx::multiroots" ON)
# ... remaining options from 4.3 ...
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(GNUInstallDirs)
include(CPM)                     # cmake/CPM.cmake: if(COMMAND CPMAddPackage) reuse parent; else include(get_cpm.cmake)
include(NumerixxDependencies)    # FXT (if NUMERIXX_WITH_FXT), Eigen (if NUMERIXX_WITH_LINALG), Boost (optional)
include(NumerixxTargets)
if(NUMERIXX_INSTALL)     include(NumerixxInstall) endif()
if(NUMERIXX_BUILD_TESTS) enable_testing() add_subdirectory(tests) endif()

# cmake/NumerixxDependencies.cmake (FXT part)
if(NUMERIXX_WITH_FXT AND NOT TARGET fxt::fxt)
  CPMAddPackage(NAME FXT URL https://github.com/troldal/FXT/archive/${NUMERIXX_FXT_REF}.tar.gz
                URL_HASH SHA256=${NUMERIXX_FXT_SHA256} DOWNLOAD_ONLY YES)
  if(NOT TARGET fxt::fxt)
    add_library(numerixx_fxt INTERFACE)
    target_include_directories(numerixx_fxt SYSTEM INTERFACE
      $<BUILD_INTERFACE:${FXT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
    target_compile_features(numerixx_fxt INTERFACE cxx_std_23)
    add_library(fxt::fxt ALIAS numerixx_fxt)
  endif()
endif()

# cmake/NumerixxDependencies.cmake (Eigen part: same DOWNLOAD_ONLY + own-target pattern)
if(NUMERIXX_WITH_LINALG)
  if(TARGET Eigen3::Eigen)                                               # a parent provides Eigen: reuse it
    set(NUMERIXX_EIGEN_TARGET Eigen3::Eigen)
  else()
    CPMAddPackage(NAME Eigen URL https://gitlab.com/libeigen/eigen/-/archive/5.0.1/eigen-5.0.1.tar.gz
                  URL_HASH SHA256=${NUMERIXX_EIGEN_SHA256} DOWNLOAD_ONLY YES)   # fallback: GIT_TAG 5.0.1
    add_library(numerixx_eigen INTERFACE)
    target_include_directories(numerixx_eigen SYSTEM INTERFACE
      $<BUILD_INTERFACE:${Eigen_SOURCE_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/numerixx-deps/eigen3>)
    # (the implementation installs fetched FXT and Eigen below include/numerixx-deps, never over another installation)
    set(NUMERIXX_EIGEN_TARGET numerixx_eigen)
  endif()
endif()

# cmake/NumerixxTargets.cmake
function(numerixx_add_module name)
  cmake_parse_arguments(ARG "" "" "DEPS" ${ARGN})
  add_library(numerixx_${name} INTERFACE)
  add_library(numerixx::${name} ALIAS numerixx_${name})
  set_target_properties(numerixx_${name} PROPERTIES EXPORT_NAME ${name})
  target_link_libraries(numerixx_${name} INTERFACE ${ARG_DEPS})
endfunction()
numerixx_add_module(core)                                              # standard library only
target_include_directories(numerixx_core SYSTEM INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_compile_features(numerixx_core INTERFACE cxx_std_23)
numerixx_add_module(deriv       DEPS numerixx::core)                   # scalar derivatives only: no Eigen
numerixx_add_module(roots       DEPS numerixx::core)
numerixx_add_module(optimize    DEPS numerixx::core)
numerixx_add_module(poly        DEPS numerixx::core)
numerixx_add_module(integrate   DEPS numerixx::core)
numerixx_add_module(interpolate DEPS numerixx::core)                   # in-house O(n) tridiagonal solvers
if(NUMERIXX_WITH_LINALG)
  numerixx_add_module(linalg     DEPS numerixx::core ${NUMERIXX_EIGEN_TARGET})
  numerixx_add_module(multiroots DEPS numerixx::linalg numerixx::deriv) # incl. gradient/Jacobian/Hessian
endif()
if(TARGET fxt::fxt)
  numerixx_add_module(pipes DEPS numerixx::core fxt::fxt)              # the only FXT consumer
endif()
add_library(numerixx INTERFACE)
add_library(numerixx::numerixx ALIAS numerixx)
target_link_libraries(numerixx INTERFACE numerixx::deriv numerixx::roots numerixx::optimize numerixx::poly
                      numerixx::integrate numerixx::interpolate
                      $<TARGET_NAME_IF_EXISTS:numerixx::multiroots> $<TARGET_NAME_IF_EXISTS:numerixx::pipes>)
if(NUMERIXX_WITH_MULTIPRECISION)
  numerixx_add_module(multiprecision DEPS numerixx::core Boost::multiprecision Boost::config
                      $<TARGET_NAME_IF_EXISTS:numerixx::linalg>)       # adapter, never in numerixx::numerixx
endif()
```

The install/export (`numerixx-config.cmake`, `SameMinorVersion`, a `find_dependency` or installed headers for the fetched Eigen), the warnings module, the compile-fail helper and the layering script follow the verified build experiments.

### 4.5 Presets, Emscripten and Windows

- **Presets** (`CMakePresets.json` v6, verified to parse on CMake 3.29 and 4.3): `msvc`, `clang-cl`, `gcc`, `clang` (libc++), `clang-asan`, `gcc-noexcept` (with the FXT pipes, now that the FXT pin includes the FXT-1 probe fix), `gcc-multiprecision`, `integration` (the consumer-build scenarios), `emscripten` (`-fwasm-exceptions`), `emscripten-jsexcept` (JavaScript-based `-fexceptions`), `emscripten-noexcept`, `emscripten-pthread` (`-fwasm-exceptions -pthread`), and workflow presets.
- **Emscripten.**
  - Test link flags: `-sALLOW_MEMORY_GROWTH=1 -sSTACK_SIZE=1MB -sEXIT_RUNTIME=1 -sNODERAWFS=1`.
  - `long double` on wasm32 is software quad.
  - Consumers may build with `-pthread`. The library creates no threads and has no statics or thread-local state (§3.2), so nothing changes; the `emscripten-pthread` leg confirms that it builds and runs under that flag.
  - The core never uses `std::stacktrace`, `std::generator`, `std::move_only_function`, `std::function_ref`, `std::copyable_function` or constexpr `<cmath>`. (`any_solver` uses `std::function`, which libc++ has.)
  - Windows: keep `EM_CACHE` short (the emsdk default location or `C:\emcache`); a long path breaks the sysroot install through MAX_PATH **[prototyped]**. Always activate through `emsdk_env`, because changing the `EM_CONFIG` path spelling clears the shared cache **[prototyped]**.
- **Windows MAX_PATH.** Short `CPM_SOURCE_CACHE` (`C:\cpm`); `CMAKE_POLICY_DEFAULT_CMP0168=NEW`; CMake ≥ 3.30; `CMAKE_INTERMEDIATE_DIR_STRATEGY=SHORT` on CMake ≥ 4.2; short binary directories.

---

## 5. Library layout

### 5.1 Directory tree

```
Numerixx/
├─ CMakeLists.txt  CMakePresets.json  README.md  CHANGELOG.md  MIGRATION.md  LICENSE
├─ cmake/  get_cpm.cmake CPM.cmake NumerixxDependencies.cmake NumerixxTargets.cmake NumerixxInstall.cmake
│          numerixx-config.cmake.in NumerixxWarnings.cmake NumerixxCompileFail.cmake CheckLayering.cmake
├─ include/numerixx/
│  ├─ numerixx.hpp                    umbrella: every scalar module; pipes.hpp if FXT is available; never
│  │                                  linalg.hpp/multiroots.hpp (Eigen), adapters or any_solver.hpp
│  ├─ config.hpp                      NXX_DELETE, NXX_NO_UNIQUE_ADDRESS, NXX_EXPECTS, NXX_ASSERT
│  ├─ core.hpp  core/{scalar,math,refined,interval,error,callable,criteria,iterate,steps,facade,compose,fn}.hpp
│  │            core/any_solver.hpp   opt-in runtime chains (std::function); not included by core.hpp
│  ├─ pipes.hpp                       the ONLY FXT include (target numerixx::pipes)
│  ├─ deriv.hpp       deriv/{stencil,step,diff,ridders,mixed,derivative_of}.hpp
│  ├─ roots.hpp       roots/{bracket,bisection,illinois,ridders,brent,rtsafe,secant,newton,search,inverse,solve}.hpp
│  │                  roots/{toms748,itp,halley,steffensen}.hpp   (optional phase 8)
│  ├─ optimize.hpp    optimize/{extremum,golden,brent_min,bracket_minimum,maximizing,newton_min,solve}.hpp
│  ├─ poly.hpp        poly/{polynomial,closed_form,aberth,format}.hpp
│  ├─ linalg.hpp      linalg/{types,solve,traits}.hpp          Eigen facade: aliases, lu/qr/cholesky_solve, vector_traits
│  ├─ multiroots.hpp  multiroots/{system,hooks,derivatives,newton,broyden,dogleg,solve}.hpp
│  ├─ integrate.hpp   integrate/{domain,romberg,gauss_legendre,gauss_kronrod,tanh_sinh,integral_of,quad}.hpp
│  ├─ interpolate.hpp interpolate/{knots,policy,tridiagonal,linear,cubic_spline,pchip,steffen,barycentric}.hpp
│  └─ adapters/{multiprecision,multiprecision_linalg}.hpp
├─ tests/  (§9; includes tests/usage/canonical_calls.cpp)   examples/   benchmarks/   docs/   tools/gen_reference.cpp
└─ .github/workflows/ci.yml   .clang-format  .clang-tidy (FXT's)
```

### 5.2 Targets and the module DAG (no cycles)

```
core ─┬─► deriv ─────────────────────┐
      ├─► roots                      │
      ├─► optimize                   │
      ├─► poly                       │
      ├─► integrate                  │
      ├─► interpolate                │
      ├─► linalg (+ Eigen 5.0.1) ────┴──► multiroots
      └─► pipes (+ fxt::fxt)
adapter (leaf): multiprecision ─► core (+ Boost::multiprecision, Boost::config; + linalg when it exists)
```

- **Scalar deriv never depends on linalg.** `deriv` holds every scalar derivative, including the two-variable mixed partial, and needs only core, so a user of the scalar modules never downloads Eigen. The vector-valued derivatives (gradient, Jacobian, Hessian) return Eigen types, so they live in `multiroots/derivatives.hpp` (namespace `nxx::multiroots`), which may depend on linalg.
- **`nxx::linalg::vector_traits<V>`** maps a user vector type to Eigen storage: `std::array<T, N>` → `vector<T, N>`, `std::vector<T>` → `vector<T>` (dynamic), Eigen column vectors → themselves **[prototyped mechanism: the same damped Newton ran unchanged on in-house and Eigen storage through the prototype's traits]**.
- **`interpolate` depends only on core**: its tridiagonal solvers are in-house O(n) code (§7.7).
- Today's back-edges disappear. `roots` does not include `deriv` or `poly`: derivative sources are recognised structurally (§6.6) **[prototyped]**. `poly` does its own Horner–Newton polish. `IsPolynomial` leaves the core. `multiroots` no longer calls `polysolve` for a parabola vertex.
- Layering allow-list: `roots:""`, `optimize:""`, `poly:""`, `deriv:""`, `integrate:""`, `interpolate:""`, `linalg:"Eigen"`, `multiroots:"linalg;deriv;Eigen"`, `pipes:"FXT"`, `adapters:*`. `core` is always allowed. multiroots includes Eigen directly because its states, Jacobians and `derivatives.hpp` results are Eigen types; its linear solves still go through the linalg facade.
- **Modules planned after v2.0 (§10.5) [sketch]** each get their own namespace, header and target (D23) and sit downstream of the v2.0 modules: `multimin` → optimize, multiroots; `fit` → multiroots (and poly once linear least-squares fitting lands); `ode` → multiroots; `chebyshev` → poly; `series` → core. No v2.0 allow-list entry changes, so `optimize` and the other scalar modules stay free of Eigen, and the Eigen-backed families (`multimin`, `fit`, `ode`) are built only when `NUMERIXX_WITH_LINALG` is ON.
- **The umbrella `<numerixx/numerixx.hpp>` never includes linalg or multiroots.** A TU that wants them includes `<numerixx/linalg.hpp>` or `<numerixx/multiroots.hpp>` explicitly, so Eigen's compile cost is paid only where it is used, and the 2 s umbrella guard (§3.6) measures the scalar modules only. The CMake target `numerixx::numerixx` still links `numerixx::multiroots` when it exists; linking adds no include.
- **Namespaces:** `nxx` (core vocabulary), `nxx::math`, `nxx::fn`, and one per module. Implementation details go in `<module>::detail`.

### 5.3 Header conventions

- `#pragma once`. Every header is self-contained; a CTest OBJECT library compiles each public header twice.
- No FXT include outside `pipes.hpp`, and no Eigen include outside `linalg/`, `multiroots/` and `adapters/multiprecision_linalg.hpp`. No `<cmath>` on constexpr paths: use `nxx::math` helpers.
- **No `auto` on arithmetic expressions** of Eigen or multiprecision types inside the library: write the concrete type (`vector<T, N> dx = …;`). Facade functions never return an Eigen expression through a deduced return type.
- **MSVC `/W4` hygiene [prototyped]:**
  - callback parameters are named `fn` or `func`, never `f`, because C4459 fires when a consumer has a global `f`;
  - parameters used only in some `if constexpr` branches are `[[maybe_unused]]` (C4100);
  - a consumer test TU compiled with `/W4 /WX` defines a global `f`.
- `NXX_NO_UNIQUE_ADDRESS` expands to `[[msvc::no_unique_address]]` under `_MSC_VER` **[prototyped]**.
- **`NXX_DELETE(reason)` [prototyped mechanism]:**
  - `#if defined(__clang__) && __clang_major__ >= 19`: `= delete(reason)`, wrapped in an in-macro `_Pragma` push/ignore(`-Wc++26-extensions`)/pop (verified: 0 warnings, reason printed). Clang 18, the floor, lacks the syntax and takes the `#else` branch.
  - `#elif defined(__GNUC__) && __GNUC__ >= 15`: `= delete(reason)`, with each public header's body wrapped in `#pragma GCC diagnostic push` / `ignored "-Wc++26-extensions"` / `pop`. GCC rejects `_Pragma` inside a declaration, and it does not define `__cpp_deleted_function` in C++23 mode even though it prints the reason.
  - `#else`: `= delete` (MSVC rejects the syntax; GCC 14 and Clang 18 lack it).
  - The macro starts on the declarator's own line **[spike]**: cl and the floor compilers print only the file and line of the deleted declaration, so that is where the reason must begin.
  - Trade-off, recorded: a constrained overload whose body is `static_assert(false, reason)` would show the reason on every compiler, but it makes `std::is_invocable_v` true, which generic code and `first_of` rely on being false.
- **`NXX_BEGIN_HEADER` / `NXX_END_HEADER` [verified in the spike]** bracket the body of every header with arithmetic.
  - **Clang** (also clang-cl and em++): no floating-point contraction in the library's code. Clang contracts `a * b + c` into a fused multiply-add by default; on FMA hardware (ARM64, x86 with `-march=haswell`) and in constant folding that rounds once instead of twice, so without the pragma a solver's path would depend on the platform (§6.1). Measured on Clang 22 with `-march=haswell`: 19 contracted FMAs in the library's code before, none after. Where Clang supports `#pragma float_control` (x86, x86-64, AArch64, RISC-V, PowerPC and SystemZ, checked with Clang 22), the setting is saved and restored (`float_control(push)`, `clang fp contract(off)`, `float_control(pop)`), so a user's `x * y + z` still fuses. Elsewhere (WebAssembly, 32-bit ARM, MIPS, LoongArch, …) the pragma is ignored with a warning, so only `clang fp contract(off)` is emitted, and contraction stays off for the rest of the translation unit; WebAssembly has no scalar FMA, so there it only makes constant folding agree with run time.
  - **GCC 15+**: they silence `-Wc++26-extensions` for `NXX_DELETE`, and do nothing about contraction. **GCC contracts by default in C++, ISO mode included** (`-ffp-contract=fast`; only ISO C defaults to off), and it has no pragma that turns contraction off for a region of code; `#pragma GCC optimize("fp-contract=off")` was not used, because GCC does not inline across differing optimisation attributes. So a GCC build for FMA hardware (AArch64; x86 with `-mfma` or `-march=haswell`) may fuse the library's arithmetic and differ from other platforms in the last ulp. Measured with GCC 16, `-std=c++23 -O2 -march=haswell`: `a * b + c` compiles to one `vfmadd`, and with `-ffp-contract=off` to a multiply and an add. The presets build baseline x86-64, which has no FMA, so the tested results are unaffected; a consumer who needs bit-identical results on an FMA target builds with `-ffp-contract=off`. **Open question for phase 1** (the `config.hpp` macros): leave this to the consumer, or add `-ffp-contract=off` to the GCC interface flags, which would change the consumer's own code as well.
  - **MSVC**: C4459 (a declaration hides a global) is silenced in the library's templates. MSVC does not contract by default (`/fp:contract` is off); a consumer that builds with `/fp:fast` or `/fp:contract` may see last-ulp differences.
- **clang-cl rules:** no folds over lambda-bearing concepts in variadic requires-clauses (use `bool` variable templates or unconstrained templates with `static_assert`); deduced return types on constrained factories. No clang-cl mangling problems were found with these shapes **[prototyped]**.

---

## 6. Core abstractions

### 6.1 Scalars and maths helpers [prototyped mechanism]

```cpp
namespace nxx {
template<class T> struct scalar_traits {};                               // primary: "not a scalar"; users specialise
template<class T>
  requires(std::numeric_limits<T>::is_specialized && !std::numeric_limits<T>::is_integer &&
           !std::numeric_limits<T>::is_exact)                           // rejects int, bool, cpp_int, cpp_rational
struct scalar_traits<T> {
    static constexpr bool is_real = true;
    static constexpr T epsilon() noexcept { return std::numeric_limits<T>::epsilon(); }
};
template<class T> inline constexpr bool is_real_v =                     // bool variable template: clang-cl safe
    requires { requires scalar_traits<std::remove_cvref_t<T>>::is_real; };
template<class T> concept real = is_real_v<T> && std::regular<T> && std::totally_ordered<T> &&
    requires(const T a, const T b) { {a+b}->std::convertible_to<T>; {a-b}->std::convertible_to<T>;
                                     {a*b}->std::convertible_to<T>; {a/b}->std::convertible_to<T>; {-a}->std::convertible_to<T>; };

namespace math {   // thin helpers over the ADL idiom: std:: serves built-in floats, ADL finds multiprecision's functions
template<real T> constexpr T abs(const T& x) noexcept {
    if consteval { return x < T(0) ? -x : x; }
    else { using std::abs; return abs(x); }                              // block-scope using: no recursion into nxx::math::abs
}
template<real T> constexpr bool isfinite(const T& x) noexcept;  // consteval branch: COMPARISONS ONLY (make(NaN), make(0, inf) rejected in static_assert everywhere)
// signbit, copysign, fma, frexp, ldexp, nextafter, floor, ceil, sqrt: the same pattern; exact or correctly rounded (IEEE 754)
// exp, log, pow: the same pattern; NOT correctly rounded, so results may differ in the last ulp between standard libraries
template<class P> constexpr P pow2(int e) noexcept;                // exact
template<class P> constexpr P root_eps(int num, int den) noexcept  // eps^(num/den) rounded to a power of two: 0x1p-26 = sqrt(eps)
{ return pow2<P>(-((std::numeric_limits<P>::digits - 1) * num) / den); }
template<class T> constexpr T midpoint(const T& a, const T& b) noexcept;   // std::midpoint for built-in floats (overflow-safe);
                                                                           // generic: a/2 + b/2 if signs differ, else a + (b-a)/2
}
}
```

- **Two groups, on purpose.** The first group is exact or correctly rounded, so its results are the same on every conforming platform. Stop tests, step-size rules, bisection midpoints and the power-of-two constants (`pow2`, `root_eps`, ITP's ⌈log₂⌉ via `frexp`) use only that group and the four basic operations, so, given the same values of f, a solver takes the same path on every platform. That extends the determinism rule (§3.2, §9.1) beyond repeated runs on one machine. It also needs the library's own arithmetic to be free of floating-point contraction. Clang and GCC both contract by default; the headers turn it off on Clang (§5.3, found and fixed in the spike), but GCC has no such pragma, so on GCC the guarantee needs a target without FMA or a build with `-ffp-contract=off` (§5.3). The second group is not correctly rounded: where the library needs it (tanh-sinh abscissae, for example) and where a user's f calls it, results may differ in the last ulp between standard libraries.
- `cpp_bin_float_50` and every other type with a suitable `numeric_limits` specialisation need nothing. A type without one specialises `scalar_traits` with `is_real = true` and `epsilon()`.
- The prototype verified the trait, the constexpr-safe `abs`/`isfinite` and the power-of-two constants; there they were spelled as CPO objects and the trait also carried an AD "primal" type, which this design drops.

### 6.2 Refined types and smart or consteval constructors [prototyped]

```cpp
namespace nxx {
namespace detail {
    inline void literal_violates_invariant(const char*) noexcept {}   // NOT constexpr: reaching it in a consteval ctor = compile error
    struct trust_me { explicit constexpr trust_me() = default; };      // library-internal construction key
    template<class Tag, class T> class refined {                       // [prototyped]
        T v_;
    public:
        constexpr refined(trust_me, T v) noexcept : v_(v) {}
        consteval refined(T v) : v_(v) { if (!Tag::check(v)) literal_violates_invariant(Tag::message); }
        template<class B> requires std::same_as<B, bool> refined(B) NXX_DELETE("a bool is not a numeric refinement");
        static constexpr auto make(T v) noexcept -> std::expected<refined, errc>;
        constexpr T value() const noexcept { return v_; }
        friend constexpr auto operator<=>(const refined&, const refined&) = default;
    };
}
template<real T> using tolerance     = detail::refined<tag::positive_tolerance, T>;  // finite, > 0: a criterion's only threshold
template<real T> using abs_tolerance = detail::refined<tag::abs_tolerance, T>;       // finite, >= 0: absolute PART of a mixed test
template<real T> using rel_tolerance = detail::refined<tag::rel_tolerance, T>;       // finite, 0 <= r < 1; roles not interchangeable
class evaluation_budget;   // 1 .. 2^32-1, the argument of max_evaluations; shaped like max_iterations below   [spike]
using progress_window   = detail::refined<tag::progress_window, std::uint8_t>;       // 2..16 (open methods, §7.2)
template<real T> using growth_factor = detail::refined<tag::growth, T>;              // > 1 (bracket expansion)
template<class T> inline constexpr bool is_refined_v = /* refined<>, max_iterations, bracket, criteria */;

class max_iterations {                                                  // 1 .. 2^32-1   [prototyped]
public:
    template<std::integral I> requires(!std::same_as<I, bool>) consteval max_iterations(I n);
    template<class B> requires std::same_as<B, bool> max_iterations(B) NXX_DELETE("a bool is not an iteration count");
    static constexpr auto make(long long n) noexcept -> std::expected<max_iterations, errc>;   // checked in the SOURCE type
    constexpr std::uint32_t value() const noexcept;
};
template<real T> class bracket {                                        // finite, lo < hi   [prototyped]
public:
    consteval bracket(T lo, T hi);                                      // bracket{2.0, 1.0}: compile error
    static constexpr auto make(T a, T b) noexcept -> std::expected<bracket, errc>;   // re-orders; a == b or non-finite: error
    constexpr T lo() const noexcept; constexpr T hi() const noexcept;
    constexpr T midpoint() const noexcept;                              // math::midpoint: no overflow on [-1.7e308, 1.7e308]
    constexpr T half_width() const noexcept;                            // hi/2 - lo/2: finite for every finite bracket
};
template<real T> class interval;      // integration domain: finite, ANY order; orientation() in {-1, 0, +1}
template<real T> class sign_bracket;  // §6.6
}
```

- **Literals** are checked when passed as arguments: `roots::brent{nxx::width_tol{1e-12}}` validates `1e-12` at compile time, and alias-template CTAD works (`constexpr nxx::tolerance t{1e-8};` **[prototyped]**).
- **Run-time values** cannot use the literal constructor. `bracket{lo, 2.0}` with a run-time `lo` fails with "not a constant expression", a diagnostic that does not mention `make()`, so the docs and the `.on()` overloads carry the run-time path:
  - `bracket<double>::make(lo, hi)`;
  - braced `{lo, hi}` or `std::pair` passed straight to a solver (§6.5);
  - configuration parsed into refined types once, with the config struct holding `tolerance<double>` and `max_iterations`.
- **Mixed criteria** take `abs_tolerance` (≥ 0) and `rel_tolerance`, with the invariant **`abs > 0 || rel > 0`** enforced in the consteval constructors and `make()`. `x_tol{0.0, 1e-8}` (purely relative, QUADPACK/GSL style) is legal; `x_tol{0.0, 0.0}` is not.
- A chain-wide budget is `with_evaluation_budget` (§6.10); there is no per-solver "remaining budget" arithmetic.
- **FLAG (portability).** MSVC 19.51 lacks P2564. A user template such as `template<class T> constexpr auto make_solver(T t) { return r::bisection{nxx::width_tol{t}}; }` escalates to an immediate function and compiles on GCC, Clang and clang-cl, but fails on MSVC with C7595 **[verified in the spike: tests/compile_fail/probe_p2564_escalation.cpp]**; only constexpr functions escalate (P2564), so without `constexpr` the call is an error on every compiler. Rule: *every function that forwards a tolerance or budget takes the refined type (`nxx::tolerance<T>`, `nxx::max_iterations`), never a raw scalar; wrap run-time values with `make()`*. On the MSVC leg the probe must fail with C7595 as its first error, which documents the gap. Also, `std::is_constructible_v<tolerance<double>, double>` is **true** on all five compilers, so generic code must detect validated inputs with `nxx::is_refined_v`, not `constructible_from` or `convertible_to` **[prototyped]**.

### 6.3 Error and result types [prototyped]

```cpp
namespace nxx {
enum class errc : std::uint8_t {
    invalid_input = 1, no_sign_change, non_finite_input, dimension_mismatch, not_increasing,
    leading_coefficient_zero, out_of_domain,                                          //  1..31 input (before iterating)
    budget_exhausted = 32, evaluations_exhausted, stalled, non_finite_value, zero_derivative,
    singular, line_search_failed, diverged, sign_change_not_root, local_minimum,      // 32..63 numerical (while iterating)
    callback_failed = 64 };                                                           // 64..   the user's callback said no
constexpr bool is_input_error(errc e) noexcept  { return std::to_underlying(e) < 32; }
constexpr bool is_budget_error(errc e) noexcept { return e == errc::budget_exhausted || e == errc::evaluations_exhausted; }

enum class algo : std::uint8_t { none = 0, user_first = 200 };   // open: roots 1-39, optimize 40-59, multiroots 60-79,
                                                                 // integrate 80-99, deriv 100-109, poly 110-119, user >= 200
namespace roots::algos { inline constexpr algo bisection{1}, brent{4}, secant{9}, newton{11}, expand{13} /* ... */; }
enum class stop_reason : std::uint8_t { exact_zero, criterion, resolution_limit, algorithm };
struct counters { std::uint32_t iterations = 0, evaluations = 0; /* operator+, == */ };
struct none { friend constexpr bool operator==(none, none) = default; };
template<class UE> using cause_slot = std::conditional_t<std::is_same_v<UE, none>, none, std::optional<UE>>;

template<class UE = none> struct fault {                // per-evaluation / per-step error
    errc code{}; std::uint32_t evals = 0;               // evaluations consumed by the failing step (D33)   [sketch: evals]
    NXX_NO_UNIQUE_ADDRESS cause_slot<UE> cause{};
};
template<class Est> struct solution : Est {             // r->x, r->used, r->by, r->how
    counters used{}; algo by = algo::none; stop_reason how{};
};
template<class Est, class UE = none> struct failure {   // trivially copyable when Est and UE are   [prototyped]
    using estimate_type = Est; using cause_type = UE;
    errc code{}; algo where = algo::none; counters used{};
    std::optional<Est> best{};                           // best-so-far; nullopt only if nothing was evaluated
    NXX_NO_UNIQUE_ADDRESS cause_slot<UE> cause{};        // the user's callback error, unchanged
};
template<class Est, class UE = none> using result = std::expected<solution<Est>, failure<Est, UE>>;
template<class R> requires requires(const R& r) { r->x; }   // not for a search result: a sign_bracket has two ends   [spike]
constexpr auto best_x(const R& r);    // optional: the solution's x, or the failure's best->x   [prototyped]
}
namespace nxx::roots {
template<real T> struct root_estimate {
    T x; T fx;
    T uncertainty;                                // enclosure width for bracketing methods (a bound); |x_k - x_{k-1}| for open methods (an indicator)
    std::optional<sign_bracket<T>> enclosure;     // WITH samples: a bracketing successor starts without re-evaluating
};
}
```

- **One result type per problem family:**

  | Family | Result type |
  |---|---|
  | roots | `result<root_estimate<T>, UE>` |
  | optimisation | `result<extremum<T>, UE>` |
  | quadrature | `result<integral<T>, UE>` |
  | systems | `result<system_estimate<V>, UE>` |
  | Ridders differentiation | `result<derivative_estimate<T>, UE>` |

  Every 1-D root solver applied to the same f returns exactly the same type, so chains type-check. One-shot differentiation (`diff`, `central`) returns `expected<T, fault<UE>>` (§6.12).
- **Search results** have a different success type (`solution<sign_bracket<T>>`) but the *same* failure type (`failure<root_estimate<T>, UE>`). The driver supports this through the `estimate`/`best` split (§6.7) **[prototyped]**. `then(search, bisection)` composes; `first_of(search, bisection)` is rejected at compile time **[prototyped]**.
- **Size.** The prototype's layout, without `uncertainty` and enclosure samples, measured `failure<root_estimate<double>>` = 64 bytes and `result` = 72 bytes on x64, MSVC and wasm32 **[prototyped]**; adding them projects to about 88 and 96 bytes. Sizes are recorded, not gated: the guideline is that errors are cheap to copy and hold no heap memory of their own (§3.4). `std::expected<solution, failure>` is not trivially copyable on libc++ (Clang, em++) although both members are **[prototyped]**, so no test asserts trivial copyability of `result<>`.

### 6.4 Callables and evaluation [prototyped]

```cpp
template<class F, class X> using callback_error_t =
    /* none for x -> T;  E for x -> expected<T, E>;  UE for x -> expected<T, fault<UE>> (unwrapped: no nesting) */;
template<class X, class F>
constexpr auto evaluate(const F& fn, const X& x) noexcept(std::is_nothrow_invocable_v<const F&, const X&>)
    -> std::expected<X, fault<callback_error_t<F, X>>>;   // NaN or ±inf -> non_finite_value   [prototyped]
template<class X, class F>
constexpr auto evaluate_sample(const F& fn, const X& x) noexcept(/*...*/)
    -> std::expected<X, fault<callback_error_t<F, X>>>;   // bracketing: only NaN is an error; ±inf is a signed sample   [sketch]
template<class F> constexpr std::uint32_t cost_of(const F&) noexcept;   // CPO: default 1; derivative_fn: its stencil's non-zero points
template<class A, class B> using common_cause_t = /* A if A == B; the non-none one; else static_assert naming transform_error */;  // [prototyped]
template<class E> constexpr bool is_fatal(const E&) noexcept;           // CPO, default false: stops first_of's fall-through
```

- **One unwrapping rule, two return conventions (D21).**
  - **Scalar-valued** function-returning APIs (`derivative_of`, `second_derivative_of`, interpolants with the `reject_outside` policy) return `x → std::expected<T, fault<UE>>`, which the rule already unwraps. Total ones (polynomials, clamping or extrapolating interpolants) return plain `T`.
  - So `r::newton{}.with_derivative(d::derivative_of(g))` has `callback_error_t == UE` and keeps the derivative's cause.
  - Declaring `failure<T, UE>` there instead would nest errors and break `first_of` **[prototyped: `neg/neg_deriv_unwrap.cpp`]**.
  - **Estimate-valued** APIs (`integral_of`, `antiderivative`, `inverse_of`, `minimizer_of`) return `result<Est, UE>` because their estimate carries more than a value (an error estimate, an enclosure, counters). To use one as a solver callback, adapt it with `nxx::fn::value_of(g[, projection])`. That maps `result<Est, UE>` to `expected<T, fault<UE>>`, keeping the cause and the evaluation count. The default projection is `&integral<T>::value` or `&root_estimate<T>::x`; for `extremum`, the projection must be named **[sketch]**.
- **Several callbacks** (f, f′, f″, J): the solver's `UE` is `common_cause_t<UE_f, UE_df, …>`. It is the same type, or the non-`none` one. Anything else is a compile error that says "map one with `.transform_error` so they agree" **[prototyped: GCC 20 lines, Clang 123, MSVC 12]**. A cause is never silently dropped.
- **Counting (D33).** States add `cost_of(fn)` per call, so `counters.evaluations` counts calls of the user's f. Counting one per callable call would count a central-difference `df` call as 1 while it makes 2 f-calls, a 33 % under-count (measured). `nxx::fn::counted(f, counter&)` remains available for instrumentation, and a property test compares the two for every solver (§9.3).

### 6.5 Problem types and accepted inputs

| Family | Accepted inputs (validated in `prepare()`, errors in-band) |
|---|---|
| bracketing roots | `bracket<T>`, `sign_bracket<T>`, `solution<sign_bracket<T>>` (search result), **braced `{lo, hi}`** (a `const T(&)[2]` overload), **`std::pair<T,T>`**, **`std::expected<bracket<T>, errc>`** (straight from `make()`); a `root_estimate` only through `.from_enclosure()` |
| open-method roots | a guess `T` (a real type: `1.0`, not `1`), or a `root_estimate<T>` whose `fx` seeds the first evaluation |
| searchers | `bracket<T>`/braced/pair start window, or a guess (+ optional limits) |
| optimisation | `bracket<T>`/braced/pair, or `min_bracket<T>` |
| quadrature | `interval<T>` (finite, any orientation), `semi_infinite<T>`, `whole_line<T>`, or `(a, b)` |
| systems | **`std::array<T,N>`** (fixed size: a wrong length is a compile error), `std::vector<T>`, Eigen column vectors (fixed or dynamic); internally Eigen storage via `linalg::vector_traits` (§5.2) |

- `.on(input)` accepts the same forms. A braced `.on({lo, hi})` needs its own `const T(&)[2]` overload; a braced list cannot pass through a generic template parameter.
- Unvalidated forms fail as `failure{errc::invalid_input or the forwarded errc, id, {0, 0}, nullopt}`, which flows through `first_of` and `then` **[sketch]**. A naive design that rejects these forms at compile time produced 18 to 161 lines of errors with no reason (measured).
- Internally, `solver.prepare(std::cref(f), input)` produces `problem<F, In>{f, in, nfev0}` or a failure. A domain's problem factories return that domain's failure type **[prototyped]**. Problems are exposed for manual stepping only.

### 6.6 The solver protocol [prototyped]

A solver type `A` over problem `P` provides, all `const`:

```cpp
static constexpr algo id;
template<class In> static constexpr bool accepts_v;       // which inputs prepare() takes (bool variable template: clang-cl safe)
static constexpr view_kind views;                         // point | enclosure | system: which criteria apply
                                                          //   (spelled views: view(s) below is the protocol function)
prepare(f, input) -> std::expected<P, failure<FEst, UE>>  // validate; evaluate what must be evaluated
init(p)           -> std::expected<S, failure<FEst, UE>>  // initial state
step(p, s)        -> std::expected<S, fault<UE>>          // one iteration: pure except for calling p.f
view(s)           -> V                                    // what stop criteria see
estimate(s)       -> SEst                                 // success payload
best(s)           -> FEst                                 // failure payload candidate (default: estimate(s))   [prototyped]
intrinsic(s)      -> std::optional<stop_reason>           // exact zero, unsplittable, Brent's tol1, ITP's bound
finish(p, sol)    -> std::optional<failure<FEst, UE>>     // optional post-condition (bracketing: pole check, §7.2): a failure if it does not hold
options()                                                 // stop, budget, derivative source, projection, observer
s.nfev                                                    // evaluations so far (in cost_of units)
```

- **Configuration is one aggregate, builders are generic.** Each solver holds `options<Stop, Deriv, Proj, Obs>` plus algorithm-specific parameters, and declares `template<class O2> using rebind = …`. The facade provides `with_stop(c)`, `with_budget(max_iterations)`, `with_derivative(d)`, `with_projection(p)` and `with_observer(o)` by rebinding that aggregate with deducing `this`. The alternative, hand-written builders per solver with a private all-members constructor and a friend declaration, was prototyped and is boilerplate **[prototyped]**.
  - Constructors take nothing, a criterion (`bisection{nxx::width_tol{1e-4}}`), or algorithm-specific parameters (`itp{itp_params{…}}`), with CTAD guides.
  - There are no positional budgets, which is what makes configuration order-independent.
  - `with_stop` is constrained on `criterion_for<C, Self::views>` and has a reasoned deleted sibling. A Tier-A test checks `S{}` and `S{criterion}` for every solver **[sketch]**.
  - `with_derivative` and `with_projection` exist only where they mean something (`uses_derivative`, `projects`), and are deleted with a reason elsewhere: "this solver does not use a derivative (newton does)"; "projection applies to open methods (secant, newton): a bracketing method keeps every iterate inside its bracket" **[spike]**.
- **Family facades carry the reasons.** Solvers derive from `bracketing_facade`, `open_facade`, `search_facade` or `system_facade`, all deducing-`this` bases with no CRTP **[prototyped]**. Each facade constrains `operator()` on `accepts_v<In>` and "F invocable on the scalar", and `.on()` on `accepts_v<In>` alone (F is not known yet), with a reasoned deletion for every rejected input and for a function that cannot take the input's scalar type. The curried solver (`bound`) is constrained on the solver being invocable with F and its input, without a reasoned deletion, so a wrong function there gets the compiler's generic error. `std::is_invocable_v<bisection<>, F, double>` is then `false`, not a hard error (as it is when the check sits inside the body). The reasons:
  - "bracketing solvers need a bracket: pass {lo, hi}, nxx::bracket<T>::make(a, b), or a search result";
  - "open methods take a guess of a real type: write 1.0, not 1";
  - "searchers take a start window or a guess";
  - "the function cannot be called with the scalar type of the bracket" (and of the guess, or of the window; for a braced `{lo, hi}` too) **[spike]**. The spike review found that the facades did not check "F invocable on the scalar", so `std::is_invocable_v` with a function of the wrong signature was a hard error inside the solver, not `false`. Each solver now declares `callable_v<F, In>`, and the facade asks it only once `accepts_v` holds (`detail::input_callable_v`, with F decayed so that a plain function does not form a `const` function type, which MSVC warns about).
  - **[spike]** A braced list or a C array binds to `const T (&)[N]` overloads with N deduced: N == 2 with a real T is a bracket (or a window), and any other length or element type is deleted with "a bracket has two ends of a real type: write {lo, hi}, for example {1.0, 2.0}". With a fixed `[2]`, a one-element `{x}` bound as `{x, 0}` and was solved on [0, x]. Arrays are left out of the call operators' catch-alls, and the `.on` catch-alls take a forwarding reference, so an array never decays to a pointer, which cl cannot order against the array overloads; a pointer keeps its reason.

  Solver-specific deletions (Newton's) stay in the solver, which must then repeat `using facade::operator();`. A test catches it if that line is forgotten.
- **Derivative sources (D13) [prototyped]:**
  - (a) A callable `df`, stored in a copyable box.
  - (b) **A derivative policy**: any `P` with `P::bind(f) -> df`, e.g. `with_derivative(deriv::numeric{deriv::central_1_2, deriv::noise{1e-10}})`. Recognised structurally, so roots does not include deriv. Each bound `df` call costs its stencil's points, which are counted.
  - (c) A structural `.derivative()` on f (polynomial, spline, `fdf(f, df)`), found by the CPO `nxx::derivative_source(f)`, named so to avoid clashing with `poly::derivative` and `deriv::derivative_estimate` **[prototyped mechanism]**.

  Without a source, the deleted overload fires with this message on Clang (GCC shows it too with the §5.3 macro):

  > `newton needs a derivative: .with_derivative(df), .with_derivative(deriv::numeric{}), a callable with .derivative(), or use secant`
- The only transition of `sign_bracket` is `narrowed(m, fm)`. Sign tests use comparisons, not products, and accept an exact zero at an endpoint **[prototyped]**.

### 6.7 The single bounded-iteration driver [prototyped; `detail::advance`, `finish` and the enclosure-aware `better` are implemented in the spike; the listing below is the pre-spike sketch]

```cpp
template<class A, class P, class Obs = no_observer>
    requires iterative_solver_for<A, P>
constexpr auto iterate(const A& alg, const P& p, const Obs& observe = {}) {
    using Fail = detail::init_failure_t<A, P>;                 // failure<FEst, UE>, from init()'s error type
    using FEst = typename Fail::estimate_type;                 // root_estimate, also for searchers
    using SEst = detail::estimate_t<A, P>;                     // sign_bracket for searchers, else FEst
    using R    = std::expected<solution<SEst>, Fail>;
    auto first = alg.init(p);
    if (!first) return R{std::unexpect, std::move(first).error()};
    auto s = *std::move(first);
    FEst best = alg.best(s);
    if (auto how = alg.intrinsic(s)) return detail::succeed<R>(alg, p, solution<SEst>{alg.estimate(s), {0, s.nfev}, A::id, *how});
    const auto& stop = alg.options().stop;
    const std::uint32_t n = alg.options().budget.value();
    for (std::uint32_t k = 1; k <= n; ++k) {
        auto next = detail::advance(alg, p, s, counters{k, s.nfev}, best);   // step(); a fault becomes a failure with best
        if (!next) return R{std::unexpect, std::move(next).error()};         //   and the failing step's evaluations
        if (const FEst e = alg.best(*next); detail::better(e, best)) best = e;   // best iterate on EVERY exit
        std::invoke(observe, alg.view(*next));                                     // logging lives here
        if (auto how = alg.intrinsic(*next))
            return detail::succeed<R>(alg, p, solution<SEst>{alg.estimate(*next), {k, next->nfev}, A::id, *how});
        switch (stop(alg.view(s), alg.view(*next), counters{k, next->nfev})) {
            case verdict::converged:
                return detail::succeed<R>(alg, p, solution<SEst>{alg.estimate(*next), {k, next->nfev}, A::id, stop_reason::criterion});
            case verdict::stalled:   return R{std::unexpect, Fail{errc::stalled, A::id, {k, next->nfev}, best, {}}};
            case verdict::exhausted: return R{std::unexpect, Fail{errc::evaluations_exhausted, A::id, {k, next->nfev}, best, {}}};
            case verdict::proceed:   break;
        }
        s = *std::move(next);                                              // the only mutation: a local
    }
    return R{std::unexpect, Fail{errc::budget_exhausted, A::id, {n, s.nfev}, best, {}}};
}
namespace detail {                                                        // internal: not public API
template<class A, class P, class S, class FEst>
constexpr auto advance(const A& alg, const P& p, const S& s, counters so_far, const FEst& best)
    -> std::expected<S, failure<FEst, /*UE*/>>;
}
```

**Implemented in the spike** (`core/iterate.hpp`), with one change of shape: a success goes through `detail::succeed<R>(alg, p, solution)`, which runs the solver's `finish(p, solution) -> std::optional<failure>` and builds the result once, in place. The sketched `finish(p, r) -> r`, which passes the `std::expected` through a by-value hook, made GCC 16 report a false `-Wmaybe-uninitialized` in every `-Werror` build of a chain.

This one loop replaces `fsolve_impl`, `fdfsolve_impl`, `search_impl`, `integrate` and `multisolve_impl`. The prototype runs `expand` through it **[prototyped]**. It removes, in one place:
- the off-by-one `maxiter`;
- max-iterations reported as success;
- the missing best estimate;
- printing to `std::cout`.

`detail::better` is per family:
- roots: **an estimate with a sign-changing enclosure beats one without; between two enclosures the narrower (nested, hence newer) wins; ties, and two estimates without enclosures, go by the smaller |fx|** — a strict weak order. Without the enclosure rule, a bisection failure could carry an older, wider enclosure. The earlier wording ("if both carry enclosures, the narrower wins; otherwise the smaller |fx|") is not transitive: the spike found three failures where each beat the next, so a static `first_of` (which merges right to left) and the run-time chain (left to right) returned different best estimates, and regrouping a static chain changed its answer. With a strict weak order, every fold order selects the same estimate;
- optimisation: fx under the optimisation sense;
- systems: the weighted merit.

Open methods can end on their worst iterate when they diverge, which is why the best iterate is tracked. The prototype's Newton on x²+1 returned `budget_exhausted` with best x = 0.0078, |f| = 1.00006 **[prototyped]**; with the progress window (§7.2) it fails earlier, as `stalled`.

### 6.8 Stop criteria [prototyped mechanism]

```cpp
enum class verdict : std::uint8_t { proceed, converged, stalled, exhausted };
enum class view_kind : std::uint8_t { point = 1, enclosure = 2, system = 4 };
struct criterion_base {   // hidden friends, bool-variable-template constraints (clang-cl safe)   [prototyped]
    template<class A, class B> requires(is_criterion_v<A> && is_criterion_v<B>)
    friend constexpr auto operator||(const A& a, const B& b) noexcept { return any_of_t<A, B>{{}, a, b}; }  // first non-proceed wins
    template<class A, class B> requires(is_criterion_v<A> && is_criterion_v<B>)
    friend constexpr auto operator&&(const A& a, const B& b) noexcept { return all_of_t<A, B>{{}, a, b}; }  // both converge; failure wins
};
// every criterion declares `static constexpr view_kind applies_to`; any_of_t / all_of_t take the intersection
```

| Criterion | Applies to | Test | Notes |
|---|---|---|---|
| `x_tol{abs[, rel]}` | point, system | \|x_k − x_{k−1}\| ≤ abs + rel·\|x_k\|. For systems it is componentwise: max_i \|dx_i\| / (abs_i + rel·\|x_i\|) ≤ 1, measured on the **full** Newton step | **Ill-formed on enclosure views**: "x_tol and step_tol compare successive iterates; bracketing methods converge on the enclosure: use width_tol{abs[, rel]} or floored_width{}" |
| `step_tol<Num, Den>{}` | point, system | \|dx\| ≤ 2^(−⌈p·Num/Den⌉)·max(\|x\|, typical), p = digits of `T`; returns the post-step iterate | Open-method default: Newton/Halley `step_tol<3,5>` (quadratic convergence puts the post-step error at O(ε)), secant `step_tol<7,10>`. The spike uses typical = 1; the `typical` hook arrives with the options in phase 3 |
| `width_tol{abs[, rel]}` | enclosure | hi − lo ≤ abs + rel·min(\|lo\|, \|hi\|): every point of the enclosure, including the returned x, is within tolerance | abs may be 0 when rel > 0 |
| `floored_width{bits = digits, scale = 1}` | enclosure | w ≤ max(2^(1−bits), 4ε)·max(scale, min(\|a\|, \|b\|)) | The default tolerance of bracketing methods. The absolute floor at `scale = 1` makes roots at 0 terminate; set `scale` for small-magnitude roots. The spike implements `bits` with scale = 1; `scale` is phase 3 |
| `f_tol{abs}` | all | \|f(best)\| ≤ abs; systems: the weighted norm with `hooks.weights` | opt-in; **ill-formed on minimisers** (an f test is meaningless for a minimum) |
| `max_evaluations{evaluation_budget}` | all | counter-aware → `exhausted` | budgets for expensive f (simulations, inner iterative solvers) |
| `min_iterations{n}` | all | counter-aware | ships with integrate (phase 7): false convergence on sin²(8πx). **A guard** **[spike]**: it delays convergence but never establishes it, so a solver accepts it only under `&&` with a convergence test (`stop_criterion_for_v`). Alone, or under `\|\|`, it would report success after n iterations without testing accuracy; the constructors and `with_stop` delete that with a reason |
| `never{}` | all | never stops | leaves the decision to `intrinsic` and the budget |
| `custom{λ}` | declared | `(prev, next[, counters]) → verdict` | user-defined tests, e.g. per-component tests on a system |

- **Implementation.** Views expose `x()`, `fx()`, `residual()` and `scale()`. Point and system views add `distance(prev)`; enclosure views add `enclosure()` but **not** `distance()`. The prototype put a `static_assert` inside the criterion (`width_tol` on an open method fails with its reason on all compilers **[prototyped]**). The design adds the solver-level check at construction (`criterion_for_v<C, S::views>`) so the error points at the user's line **[implemented in the spike: a constrained constructor and `with_stop`, each with a reasoned deleted sibling]**.
- **Windowed tests are not criteria.** A criterion is a pure function of (previous view, next view, counters), so it cannot see a window; a windowed `no_progress{factor}` criterion is unimplementable. Cycle and divergence detection live in open-method states (§7.2).
- **Internal-test solvers.** Brent, golden and Brent-min (and TOMS748 and ITP, phase 8) take a width criterion *as their tolerance*: `brent{nxx::width_tol{abs, rel}}`, default `floored_width{}`.
  - Brent's `tol1 = max(threshold/2, 2ε|b|)`, where `threshold` is the tolerance's enclosure form on the current enclosure [min(b, c), max(b, c)] (`width_tol`: abs + rel·min(|b|, |c|)). The intrinsic test stops when |c − b|/2 ≤ tol1 **[implemented in the spike]**.
  - It reports `stop_reason::criterion` only when |c − b| ≤ threshold. That is the width criterion's own guarantee (§9.3), so Brent is sound without a weaker Brent-specific bound. When the tolerance is below Brent's floor (threshold < 4ε|b|), Brent stops at the floor, with width ≤ 4ε|b|: it reports `stop_reason::criterion` if that width still meets the threshold, and `stop_reason::resolution_limit` otherwise.
  - The earlier form had two defects, and the spike's soundness test caught both. The form was `tol1 = 2ε|b| + threshold/2` with `threshold = abs + rel·|b|`. It reported `criterion` for widths up to threshold + 4ε|b|, and it measured the relative part at b instead of at the enclosure's smaller end.
  - Their external `Stop` defaults to `never{}`, so there are not two sources of truth. Extra external criteria are allowed only when combined explicitly.
- **Per-module defaults (all expressions in `T`):**

  | Module | Default criterion | Budget |
  |---|---|---|
  | roots, bracketing | `floored_width{}`, as internal tolerance or external stop | bisection 200 (enough for `cpp_bin_float_50`: 168 bits need about 165 halvings from a unit bracket), brent 100 |
  | roots, open | `step_tol` + progress window (§7.2) | Newton/Halley 30, secant 50 |
  | optimize | Brent's intrinsic test, `tol1 = rel·\|x\| + abs/3`, `rel = root_eps<T>(1, 2)`, `abs = rel·typical` | 200 |
  | integrate | QUADPACK acceptance with `rel = root_eps<T>(1, 2)` and a roundoff floor (§7.6) | per rule |
  | systems | `x_tol` on the full step \|\| weighted `f_tol`, plus the `local_minimum` test (§7.5) | per solver |

  A literal default such as `x_tol{1e-10, 1.5e-8}` can never be met in `float`.

### 6.9 Manual stepping, observation and per-iteration projection

```cpp
namespace nxx {
template<class A, class P>
class steps_view : public std::ranges::view_interface<steps_view<A, P>> {     // [prototyped]
    using UE = typename detail::init_failure_t<A, P>::cause_type;             // the user's callback error (none for plain f)
public:
    using element = std::expected<detail::state_t<A, P>, fault<UE>>;
    class iterator;                                       // input iterator; caches the current element
    constexpr steps_view(A alg, P problem);
    constexpr iterator begin() const;                     // element 0 is init(p), its failure mapped to a fault
    constexpr std::default_sentinel_t end() const noexcept;
};  // input range; ends after the first error or the first intrinsic stop, otherwise infinite
}

// Manual stepping and tracing.                                                                [prototyped mechanism]
namespace r = nxx::roots;
const auto solver = r::brent{};
auto p = solver.prepare(std::cref(fn), nxx::bracket{1.0, 2.0});          // expected<problem, failure>
if (p)
    for (const auto& st : nxx::steps_view{solver, *p} | std::views::take(8))
        if (st) trace(solver.estimate(*st).x);

// States to estimates, lazily.                                                                          [prototyped]
const auto nt = r::newton{}.with_derivative(df);
auto pn = nt.prepare(std::cref(fn), 3.0);
auto xs = nxx::steps_view{nt, *pn}
        | std::views::transform([&](const auto& st) { return st.transform([&](const auto& s) { return nt.estimate(s); }); })
        | std::views::take(4);

// Observer: logging without polluting stop criteria.                                                    [sketch]
auto logged = r::brent{}.with_observer([&](const auto& view) { log(view.x()); });

// Per-iterate projection: applied to the PROPOSED iterate before f is evaluated.                        [prototyped]
auto res = r::secant{}.with_projection(r::clamp_to{xmin, xmax})(objective, guess);   // box constraint on each iterate
// A projected iterate pinned at the edge -> errc::stalled with best at the edge (never "the root is the boundary").
auto final_x = res.transform([&](const auto& s) { return Length{std::clamp(s.x, xmin, xmax)}; });   // final projection into a user type
```

- **`steps_view` is the public manual-stepping and tracing API**, delivered in phase 3. It is an input range of `expected<state, fault>` **[prototyped: satisfies `input_range` and `view`, composes with `views::take` and `views::transform` on all 9 configurations]**.
  - The range ends after the first error or the first intrinsic stop, and is otherwise infinite: bound it with `views::take` or stop on your own test. It does not apply the solver's stop criterion or budget; those belong to `nxx::iterate`.
  - Like `std::generator`, its iterator caches the current element.
  - f must outlive the view, because the problem holds `std::cref(f)`.
  - A hand-written `s = step(p, *s)` loop does not type-check, because `init` and `step` return different error types **[prototyped: `neg/neg_manual_step.cpp`]**; the view hides that.
  - A hand-rolled iteration that inspects or adjusts each iterate (for example an outer loop around a nested solve) maps onto it, or onto `multiroots::newton` with a `project` hook.
- Projection cost: a clamped Newton from 10 on [1, 3] needs 8 iterations **[prototyped]**. Where a bracket is known, `rtsafe` is the better tool.

### 6.10 Combinators [prototyped mechanism; the named types and policies are sketch; `any_solver` is prototyped]

```cpp
namespace detail { template<class F> class copyable_box; }   // std::ranges movable-box technique: assignment = destroy + construct

struct continue_unless_fatal {                                // default first_of policy
    template<class Est, class UE> constexpr bool operator()(const failure<Est, UE>& e) const noexcept {
        if constexpr (std::is_same_v<UE, none>) return true;  // input and numerical errors: try the next alternative
        else return !(e.cause && nxx::is_fatal(*e.cause));    // the user's CPO decides which callback errors are final
    }
};

template<class Policy, class S1, class S2>
class first_of_t {                                            // a named value: copy-assignable whenever S1, S2 are copyable
    NXX_NO_UNIQUE_ADDRESS Policy policy_;
    detail::copyable_box<S1> s1_; detail::copyable_box<S2> s2_;
public:
    constexpr first_of_t(Policy p, S1 a, S2 b);
    template<class... A> constexpr auto operator()(const A&... a) const {
        if constexpr (!detail::callable_v<S1, A...> || !detail::callable_v<S2, A...>) {
            detail::report_not_callable<S1, S2, A...>();      // static_assert: "nxx::first_of: alternative #i is not callable
                                                              //  with these arguments; did you forget .on(input)?"
        } else if constexpr (!std::is_same_v<detail::result_t<S1, A...>, detail::result_t<S2, A...>>) {
            static_assert(detail::always_false<A...>, "nxx::first_of: every alternative must return the same "
                          "std::expected<solution<Est>, failure<Est, UE>>; adapt the odd one with .transform/.transform_error");
        } else {                                              // [prototyped: this body and the mismatch message]
            auto r1 = std::invoke(*s1_, a...);
            if (r1 || !policy_(r1.error())) return r1;        // lazy: later alternatives never run
            auto r2 = std::invoke(*s2_, a...);
            if (r2) { r2->used = r2->used + r1.error().used; return r2; }   // success pays for failed attempts
            return decltype(r1){std::unexpect, detail::merge(r1.error(), std::move(r2).error())};
        }                                                     // merge: last code/cause, BEST estimate (enclosure-aware), total cost
    }
};
template<class... S> constexpr auto first_of(S... s);                    // first_of_t<continue_unless_fatal, S1, first_of_t<...>>
template<class P, class... S> constexpr auto first_of_with(P policy, S... s);

template<class S1, class S2>
class then_t {                                                // Kleisli with environment: stage 2 gets f AND stage 1's value
    detail::copyable_box<S1> s1_; detail::copyable_box<S2> s2_;
public:
    template<class F> constexpr auto operator()(const F& fn) const {
        auto r1 = std::invoke(*s1_, fn);
        using V1 = typename decltype(r1)::value_type;
        if constexpr (!std::is_invocable_v<const S2&, const F&, const V1&>) {         // [prototyped: NXX_THEN_CONTRACT]
            static_assert(std::is_invocable_v<const S2&, const F&, const V1&>,
                "nxx::then: stage 2 cannot start from stage 1's result (a bracketing solver needs a bracket or a "
                "search result; use .from_enclosure() or put a searcher first)");
            return r1;
        } else {
            using R2 = decltype(std::invoke(*s2_, fn, *r1));
            if (!r1) return R2{std::unexpect, std::move(r1).error()};
            return detail::add_cost(std::invoke(*s2_, fn, *r1), r1->used);
        }
    }
};   // then(s1, s2, s3, ...) == then(then(s1, s2), s3, ...)
```

| Combinator | Semantics |
|---|---|
| `first_of(s…)` | First success. On total failure it reports the last code and cause, the best estimate over all attempts, and the total cost. Falls through unless the user's error `is_fatal` |
| `first_of_with(policy, s…)` | Same, with a custom `policy(const failure&) → bool` (continue?) |
| `then(s1, s2, …)` | Kleisli staging with static typing: search → bracket solver ✓, bracket solver → open method ✓ (seeded with x and fx, no re-evaluation), bracket solver → rtsafe or another bracket solver via `.from_enclosure()` ✓ (reuses the sampled enclosure; fails with `no_sign_change` if there is none), open method → bracket solver ✗ |
| `warm_fallback(s1, s2)` | Restart s2 from s1's `failure::best`, including its sampled enclosure |
| `with_evaluation_budget(evaluation_budget n, chain)` | One evaluation budget for a whole chain. Each call wraps f in a local counting guard; once n evaluations are spent, every later evaluation fails at zero cost with `evaluations_exhausted`, so the remaining alternatives fail instantly and the merged failure carries the best estimate **[sketch]** |
| `any_solver<F, Est, UE>` + `first_of(range)` | Opt-in run-time chains over a run-time list of type-erased curried solvers: the same laziness, merge and cost accounting as `first_of`, and the result is itself an `any_solver` **[prototyped]**; the `is_fatal` policy and `first_of_with(policy, range)` as for static chains **[sketch]** |

**Run-time chains (opt-in) [prototyped].** Static chains cannot be assembled from configuration. `any_solver<F, Est, UE>` erases the type of any curried solver (a `bound` solver, a `then_t`, a `first_of_t`, another run-time chain) whose call `const F& → result<Est, UE>` matches exactly, for one fixed callable type `F` such as `std::function<double(double)>`. In the library it is `<numerixx/core/any_solver.hpp>`, the only header that uses `std::function`, included by neither `core.hpp` nor `numerixx.hpp`. The code below is the prototype's `nxx/runtime.hpp`, verbatim as compiled; it passed on all 9 configurations, with results bit-identical to the equivalent static chains:

```cpp
// Prototype of run-time solver chains: any_solver (a copyable, type-erased curried solver) and first_of
// over a run-time range, with the same semantics as the static first_of in core.hpp.
// Standard library only (no FXT). Heap: std::function may allocate when a solver is wrapped or copied
// (target larger than its small buffer); a call never allocates by itself.
#pragma once
#include "core.hpp"
#include <functional>
#include <ranges>
#include <vector>

namespace nxx {

// A curried solver for callables of type F: a copyable value s with s(f) -> R for const F& f. The call is checked
// before copyability [spike; see Exact types below].
template<class S, class F, class R>
concept curried_solver_for = std::invocable<const S&, const F&> &&
                             std::same_as<std::invoke_result_t<const S&, const F&>, R> && std::copy_constructible<S>;

// any_solver<F, Est, UE>: "const F& -> result<Est, UE>" for ONE fixed callable type F
// (e.g. std::function<double(double)>). Built on std::function (std::move_only_function is missing on
// libc++ 22 / Emscripten, and copyability is the point here).
// Invariant: never empty. There is no default constructor and there are no move operations (a move is a
// copy), because a moved-from std::function is unspecified and may be empty.
template<class F, class Est, class UE = none>
class any_solver {
public:
    using function_type = F;
    using estimate_type = Est;
    using cause_type = UE;
    using result_type = result<Est, UE>;

    template<class S>
        requires(!std::same_as<std::remove_cvref_t<S>, any_solver> &&
                 curried_solver_for<std::remove_cvref_t<S>, F, result<Est, UE>>)
    any_solver(S&& s) : impl_(std::forward<S>(s)) {}   // implicit: every matching solver value IS an any_solver

    template<class S>
        requires(!std::same_as<std::remove_cvref_t<S>, any_solver> &&
                 !curried_solver_for<std::remove_cvref_t<S>, F, result<Est, UE>>)
    any_solver(S&&) NXX_DELETE("nxx::any_solver<F, Est, UE>: needs a copyable curried solver `const F& -> "
                               "result<Est, UE>` with exactly this F, Est and UE");

    any_solver(const any_solver&) = default;
    any_solver& operator=(const any_solver&) = default;   // replaces the whole value (strong guarantee)
    ~any_solver() = default;

    [[nodiscard]] result_type operator()(const F& f) const { return impl_(f); }

private:
    std::function<result_type(const F&)> impl_;
};

template<class T> inline constexpr bool is_any_solver_v = false;
template<class F, class Est, class UE> inline constexpr bool is_any_solver_v<any_solver<F, Est, UE>> = true;

namespace detail {
    template<class F, class Est, class UE>
    class runtime_first_of {   // owns its alternatives; never changed after construction
        std::vector<any_solver<F, Est, UE>> alts_;
    public:
        explicit runtime_first_of(std::vector<any_solver<F, Est, UE>> a) : alts_(std::move(a)) {}
        [[nodiscard]] result<Est, UE> operator()(const F& f) const {
            using R = result<Est, UE>;
            using Fail = failure<Est, UE>;
            if (alts_.empty()) return R{std::unexpect, Fail{errc::invalid_input, algo::none, {}, std::nullopt, {}}};
            std::optional<Fail> failed;   // merged so far: last code/cause, BEST estimate, total cost
            for (const auto& s : alts_) {
                R r = s(f);
                if (r) { if (failed) r->used = r->used + failed->used; return r; }   // lazy; success pays for failures
                failed = failed ? nxx::detail::merge(*failed, std::move(r).error()) : std::move(r).error();
            }
            return R{std::unexpect, *std::move(failed)};
        }
    };
}    // namespace detail

// Same shape as core.hpp's single-argument first_of(S) (one by-value parameter, `class` template head) so
// that this constrained overload is more specialised and wins for ranges of any_solver.
// Accepts std::vector, std::span, std::array, ... of any_solver; the chain owns a copy (it is curried,
// so borrowing the range would dangle). An empty range fails with errc::invalid_input when called.
template<class Rng>
    requires(std::ranges::input_range<Rng> &&
             is_any_solver_v<std::remove_cvref_t<std::ranges::range_reference_t<Rng>>>)
auto first_of(Rng alts) {
    using S = std::remove_cvref_t<std::ranges::range_reference_t<Rng>>;
    using Chain = detail::runtime_first_of<typename S::function_type, typename S::estimate_type, typename S::cause_type>;
    if constexpr (std::is_same_v<Rng, std::vector<S>>) return S{Chain{std::move(alts)}};
    else return S{Chain{std::ranges::to<std::vector<S>>(alts)}};
}

}    // namespace nxx
```

Usage, from the prototype's `test_runtime.cpp` (abridged; its solver spellings predate some renames):

```cpp
using fn_t     = std::function<double(double)>;
using solver_t = nxx::any_solver<fn_t, r::root_estimate<double>>;          // plain callbacks: UE = none
using result_t = solver_t::result_type;
std::optional<solver_t> make_solver(std::string_view name) {
    if (name == "newton")    return r::newton{}.with_derivative(df).on(0.0);
    if (name == "newton_fd") return r::newton{}.with_derivative(d::numeric{}).on(1.0);
    if (name == "secant")    return r::secant{nxx::default_step{}, 5}.on(0.0);
    if (name == "bisection") return r::bisection{}.on(nxx::bracket{0.0, 2.0});
    if (name == "brent")     return r::brent{}.on(nxx::bracket{0.0, 2.0});
    return std::nullopt;
}
std::expected<solver_t, std::string> make_chain(const std::vector<std::string>& names) {
    std::vector<solver_t> alts;
    for (const auto& n : names) {
        auto s = make_solver(n);
        if (!s) return std::unexpected("unknown solver '" + n + "'");
        alts.push_back(*s);
    }
    return nxx::first_of(std::move(alts));
}
auto chain = make_chain({"newton", "secant", "bisection"});
result_t rr = (*chain)(f);                                                  // f is a const fn_t
solver_t s  = static_chain;                                                 // a static chain converts too
auto mixed  = nxx::first_of(r::newton{}.with_derivative(df).on(0.0), *chain);   // static first_of over a run-time chain
```

- **Design deltas from the prototype [sketch].** The library adds the same `Policy` as static chains (default `continue_unless_fatal`) and `first_of_with(policy, range)`; the prototype's run-time chain, like its static `first_of`, falls through on every failure. The merge gets the enclosure-aware `better` of §6.7; the prototype keeps the lower |f| and lets the later attempt win ties, exactly as its static chain does.
- **Verified [prototyped], identical output on all 9 configurations:**
  - `{newton, secant, bisection}` on x²−2 equals the static chain bit for bit (x = 1.4142135623730949, 57 iterations, 61 evaluations, found by bisection; the cost includes every attempt), both when the static chain is given an `fn_t` and when it is given a plain lambda;
  - on x²+1 it fails like the static chain (`no_sign_change` from bisection, 6 iterations and 10 evaluations in total, best x = 0 with its bracket [0, 2]);
  - `{newton_fd, brent}` equals the static numeric-derivative chain (13 evaluations);
  - an unknown name is an error string; an empty range fails with `invalid_input`, zero cost and no best estimate; an alternative after a success is never called;
  - assignment, `std::move`, swaps inside a vector, a `span` overload, and nesting inside the static `first_of` all behave;
  - a fallible callback (`any_solver<gfn_t, est, E>`, with a user error enum `E`) keeps the user's error when everything fails;
  - `static_assert`s: copy-constructible and copy-assignable, not default-constructible, not nothrow-invocable; a searcher (`expand_out`, another estimate type) and an `int` are rejected; `first_of` over a `vector`, `span<const>` or `array` returns an `any_solver`.
- **Size and allocations [prototyped].** `sizeof(any_solver)` is that of the underlying `std::function`: 32 B (x64 libstdc++), 48 B (x64 libc++), 64 B (x64 MSVC STL), 24 B (wasm32 libc++). For comparison, `bisection.on(…)` is 24 B, `newton.on(…)` 16 B, the 3-solver static chain 56 B and `result<root_estimate<double>>` 72 B. Counted by replacing `operator new`: wrapping `bisection.on(…)` allocates once on libstdc++ and wasm32 and not at all on x64 libc++ and MSVC; wrapping the static chain allocates once everywhere; **a call with an existing `fn_t` never allocates**, but passing a lambda that is not already an `fn_t` converts it on every call (one allocation per call for a 64-byte capture). Timing was not measured.
- **Pitfalls, recorded:**
  - *Never empty.* A moved-from `std::function` may be empty, and calling it throws `bad_function_call` (or aborts under `-fno-exceptions`). So `any_solver` has no default constructor and no move operations: a move is a copy, which may allocate. Chains are built at configuration time, so that is acceptable.
  - *Copyable solvers only.* A solver with a move-only capture cannot be wrapped; `std::move_only_function` and `std::copyable_function` are missing on libc++ 22 / Emscripten.
  - *Small buffers differ*: 16 B and trivially copyable targets only (libstdc++), three pointers (libc++), 56 B (MSVC). Convert f to `F` once, outside the call.
  - *Indirection.* Each evaluation of f goes through the caller's `std::function`, and each alternative adds one indirect call; nothing inlines across the boundary. `operator()` cannot be `noexcept`, because f may throw. All three standard libraries compile `std::function` under `-fno-exceptions`.
  - *Exact types.* `F`, `Est` and `UE` must match exactly. A solver that is callable with F but returns another result type hits the reasoned deleted constructor (Clang prints the reason; GCC shows the deleted function and the source line holding the reason); anything else is simply not convertible, so `any_solver` is not a viable conversion target for every type, and overload sets that mention it stay unambiguous **[spike]**. The concept checks the call before copyability: copying a chain that holds an `any_solver` asks whether its copyable box converts to an `any_solver`, and checking copyability first makes that question depend on itself, which Clang rejects **[spike]**. A `static_assert` in the constructor body would instead make `std::is_constructible_v` true (§5.3).
  - *Overload shape.* To beat the one-argument `first_of(S)`, the range overload takes one by-value parameter with a plain `class` template head and a requires-clause; `template<std::ranges::input_range R>` or `R&&` made calls ambiguous. All 5 compilers resolve the compiled form.
  - *Ownership.* `first_of` copies the range into a vector (a curried chain must not borrow a `span`); an rvalue vector is moved in.
  - **FLAG (illegal states).** An empty run-time chain is representable and fails in-band with `invalid_input` when called, because run-time configuration can be empty. The alternatives, `first_of` returning `expected<any_solver, errc>` or taking a non-empty range type, were not chosen; `make_chain` above already rejects bad configuration before a chain exists.
- **Costs, in short:** a possible allocation when a solver is wrapped or copied, an indirect call per alternative, and an indirect call per evaluation through `F`. None of it matters next to an expensive f; static chains stay the default because they cost nothing and run in `constexpr`.
- **Dropped:** `retry(solver, inputs)`. For static inputs it is `first_of(s.on(a), s.on(b))`; over a run-time range it is a `first_of` over `any_solver`s (or FXT-7 `first_success` once it exists). **Deferred:** `first_of_all` (no demonstrated need; it can be added later without changing the core). **Moved:** `inverse_of` lives in roots (§6.12), because its default solver is Brent.
- `fxt::operator||` is deliberately not used: it is eager, requires both sides to have the same type, and silently becomes the built-in `bool` operator without its using-declaration.
- **Robust idioms (documented):**
  - `then(search, brent)` (or `then(search, toms748)` once phase 8 lands);
  - `then(bisection{width_tol}, rtsafe.from_enclosure())`;
  - `first_of` only over alternatives with different failure modes, cheapest and fastest-failing first.

  In a numerical probe, a cycling Newton spent 102 of a chain's 130 evaluations before the next alternative ran. The progress window (§7.2) now cuts that short.

### 6.11 The headline example, end to end

The spellings below are the design's. The prototype compiled and ran the same example with older spellings on all 9 configurations **[prototyped]**: `r::secant{nxx::default_step{}, 5}` and `r::bisection{nxx::x_tol{1e-4}}`. The spike compiles and runs this block as written (checked with GCC 16), and its tests cover each part on all 12 presets **[spike]**.

```cpp
#include <numerixx/roots.hpp>
#include <numerixx/pipes.hpp>                        // only for the FXT pipes at the end
namespace r = nxx::roots;

constexpr auto f  = [](double x) { return x * x - 2.0; };
constexpr auto df = [](double x) { return 2.0 * x; };

// "If one solver fails, the next one tries": three different solvers, three different inputs, one value.
constexpr auto chain = nxx::first_of(
    r::newton{}.with_derivative(df).on(0.0),         // fails: f'(0) = 0 -> errc::zero_derivative
    r::secant{}.with_budget(5).on(0.0),              // fails within its 5-iteration budget
    r::bisection{}.on(nxx::bracket{0.0, 2.0}));      // succeeds (literal bracket, validated at compile time)
static_assert(chain(f).has_value() && chain(f)->by == r::algos::bisection);   // the whole chain runs at compile time
auto res = chain(f);   // x = 1.4142135623730949; used = {57 iterations, 62 evaluations} incl. failed attempts [spike; the prototype counted 61: the derivative call in Newton's failing step is now counted]

// Numeric-derivative Newton inside a curried chain (f unknown when the chain is built).
constexpr auto chain_fd = nxx::first_of(r::newton{}.with_derivative(nxx::deriv::numeric{}).on(1.0),
                                        r::brent{}.on(nxx::bracket{0.0, 2.0}));   // [spike: newton wins, 5 it, 16 evals incl. the stencil's 2 per derivative; the prototype reported 6 it, 13 evals]

// Sequential refinement: grow a window until f changes sign, coarse enclosure, Newton polish.
constexpr auto pipeline = nxx::then(r::expand{}.on(nxx::bracket{2.0, 2.5}),
                                    r::bisection{nxx::width_tol{1e-4}},        // width, NOT x_tol (x_tol would falsely converge)
                                    r::newton{}.with_derivative(df));          // rtsafe{}.from_enclosure() stays inside it

// Run-time inputs are first-class.                                                                   [spike]
auto r1 = r::brent{}(f, {lo, hi});                                         // braced, validated in-band
auto r2 = r::brent{nxx::width_tol{1e-12}}.with_budget(60).on(nxx::bracket<double>::make(lo, hi))(f);

// A fallible callback keeps its own error through solver and combinator.                              [prototyped]
enum class eval_error { domain };
auto g = [](double x) -> std::expected<double, eval_error> {
    if (x < 0) return std::unexpected(eval_error::domain); return x * x - 2.0; };
auto bad = r::bisection{}(g, nxx::bracket{-1.0, 2.0});   // code == callback_failed, *cause == eval_error::domain
auto ok  = nxx::first_of(r::bisection{}.on(nxx::bracket{-1.0, 2.0}), r::brent{}.on(nxx::bracket{0.0, 2.0}))(g);

// Consumption with FXT pipes (or with std::expected members).                                         [prototyped]
using fxt::operator|;
double x  = chain(f) | fxt::transform([](const auto& s) { return s.x; })
                     | fxt::value_or(std::numeric_limits<double>::quiet_NaN());
auto len  = chain(f) | fxt::transform([](const auto& s) { return Length{s.x}; });        // into a user strong type
double x2 = nxx::best_x(chain(f)).value_or(fallback);                                   // accept the best estimate on failure
```

Also verified **[prototyped]**:
- no sign change is an error carrying the best endpoint;
- a warm start from a starved bisection into Newton (8 iterations, 16 evaluations in total);
- a clamped Newton, and a Newton pinned at the edge giving `stalled` with best x = 3;
- criteria algebra (`never || max_evaluations{10}` gives `evaluations_exhausted` after 10 evaluations);
- the `first_of` mismatch as a single diagnostic line on GCC and MSVC.

### 6.12 Function-returning APIs

| Factory | Returns | Notes |
|---|---|---|
| `deriv::derivative_of(fn, stencil = central_1_2, step = optimal{})` | `derivative_fn`: x → `expected<T, fault<UE>>`; `cost_of` = stencil points | holds f in a copyable box; `const` call; runs at compile time (d(x³)/dx at 2 within 1e-8) **[prototyped]** |
| `deriv::numeric{stencil, step}` | a derivative **policy**: `bind(f) → derivative_fn` | for Newton in curried chains **[prototyped]** |
| `deriv::second_derivative_of` | callable | per-order optimal step; core only |
| `multiroots::gradient_of`, `jacobian_of`, `hessian_of` | callables returning `linalg::vector`/`linalg::matrix` | in `multiroots/derivatives.hpp` (needs linalg); per-order optimal steps; a correct Hessian |
| `roots::inverse_of(fn, bracket, solver = brent{})` | **smart constructor** → `expected<inverse_fn, failure>`, evaluating the endpoints once; each call y → `result<root_estimate<T>>` reuses the cached samples shifted by y, and y outside [min(f_lo, f_hi), max(f_lo, f_hi)] → `out_of_domain` with zero evaluations | D33 |
| `optimize::minimizer_of(family, bracket, solver)` | θ → `result<extremum<T>>` | parametric optimisation: the minimiser (or, through `maximizing`, the maximiser) as a function of a parameter θ |
| `optimize::maximizing(solver)` | a solver value that minimises −f and restores fx on success and failure | composes with `.on`, `first_of`, `then` |
| `integrate::integral_of(fn, rule)` | (a, b) or domain → `result<integral<T>>` | fixes `integralOf` |
| `integrate::antiderivative(fn, a, rule)` | x → `result<integral<T>>` | |
| `interpolate::make_cubic_spline(xs, ys, boundary, outside)` etc. | `expected<interpolant, errc>`; tag arguments deduce the policies | validated knots; eager coefficients |
| `poly::derivative(p)`, `antiderivative(p)`, `compose(p, q)` | polynomials | `.derivative()` makes Newton exact on polynomials |
| `fn::negate`, `fn::shift`, `fn::extend_linearly(fn, bracket, slopes)`, `fn::counted`, `fn::catching`, `fn::value_of(g[, projection])` | callables | `value_of` turns an estimate-valued callable into a solver callback (§6.4) **[sketch]**; `extend_linearly` extends a function beyond its domain by its end tangents (C¹), for algorithms that must probe outside it; where the solver supports it, `.with_projection` (D20) is the better tool |

Every `_fn` type is a named class template, copy-assignable whenever its callables are copy-constructible. The prototype's closure-based `derivative_of` holding a capturing lambda was not assignable **[prototyped]**; the copyable box fixes that.

### 6.13 The one-call convenience facade

```cpp
nxx::roots::solve(f, {lo, hi});                        // corpus-chosen bracketing default (Brent provisionally)
nxx::roots::solve(f, x0);                              // then(expand{}.on(window(x0)), brent{}): window [x0(1-2^-5), x0(1+2^-5)], ±2^-5·typical at 0
nxx::roots::solve(f, df, x0);                          // then(expand{}.on(window(x0)), rtsafe{}.with_derivative(df)): safeguarded
nxx::optimize::minimize(f, {a, b});                    // brent_min -> extremum{x, fx}
nxx::optimize::maximize(f, {a, b}, solver = brent_min{});   // maximizing(solver); fx = f(x) on success AND failure
nxx::deriv::central(f, x);                             // central_1_2, optimal relative step
nxx::integrate::quad(f, a, b);                         // adaptive G7K15; a > b and a == b legal
nxx::multiroots::solve(F, x0);                         // dogleg + Broyden once it passes the phase-5 corpus; until then damped Newton, FD Jacobian
nxx::interpolate::make_cubic_spline(xs, ys);           // expected<cubic_spline<T>>
```

Each facade is a *composition of the core*, not a second implementation, and returns the core's result types.

### 6.14 Canonical calls (normative; compiled in `tests/usage/canonical_calls.cpp`)

Each module phase adds its calls to this file. Calls 1–3, 9 and 10 are spike exit criterion 8.

| # | Task | Spelling |
|---|---|---|
| 1 | Root in a run-time bracket | `nxx::roots::solve(f, {lo, hi})` or `r::brent{}(f, {lo, hi})` |
| 2 | Newton with an analytic derivative | `r::newton{}.with_derivative(df)(f, x0)`; safeguarded: `nxx::roots::solve(f, df, x0)` |
| 3 | Secant with a box clamp | `r::secant{}.with_projection(r::clamp_to{xmin, xmax})(objective, g).transform([](const auto& s) { return Length{s.x}; })`. A pinned iterate is `stalled`: use `nxx::best_x(r)` to accept the best estimate. |
| 4 | Golden-section maximum on an interval | `nxx::optimize::maximize(f, {a, b}, nxx::optimize::golden{})` → `r->x`, `r->fx` |
| 5 | Derivative at a point | `nxx::deriv::central(f, x).value_or(NaN)`; for an f with about 1e-10 relative noise: `nxx::deriv::diff(f, x, d::central_1_2, d::noise{1e-10})` |
| 6 | Integral | `nxx::integrate::quad(f, a, b)` → `->value`, `->error` |
| 7 | Spline | `auto s = nxx::interpolate::make_cubic_spline(xs, ys, natural{}, clamp_to_domain{}); double y = s ? (*s)(x) : NaN;` (a clamping policy returns `T`, not `expected`) |
| 8 | 2×2 system | `nxx::multiroots::solve([](const std::array<double, 2>& x) { return std::array{x[0]*x[0] + x[1]*x[1] - 4.0, x[0] - x[1]}; }, std::array{1.0, 0.5})` |
| 9 | Three-solver chain | `nxx::first_of(r::newton{}.with_derivative(df).on(x0), r::secant{}.on(x0), r::brent{}.on({lo, hi}))(f)` |
| 10 | Coarse bisection, then secant | `nxx::then(r::bisection{nxx::width_tol{1e-3}}.with_budget(100).on({lo, hi}), r::secant{})(f)` |
| 11 | Run-time chain | `using solver_t = nxx::any_solver<fn_t, r::root_estimate<double>>; nxx::first_of(std::vector<solver_t>{r::secant{}.on(x0), r::brent{}.on({lo, hi})})(fn_t{f})` **[prototyped mechanism]** (convert f to `fn_t` once) |
| 12 | Trace the iterates | `auto p = s.prepare(std::cref(f), nxx::bracket<double>::make(lo, hi)); for (const auto& st : nxx::steps_view{s, *p} \| std::views::take(20)) …` **[sketch]** |
| — | Fix f, vary the input | `auto solve_at = std::bind_front(r::brent{}, f); solve_at({lo, hi});` |

---

## 7. Module-by-module design

Legend: **Keep** = port the algorithm as a pure step. **Fix** = port and correct. **Add** = new. **Drop** = removed. "Must not port" lists bugs confirmed in master and dev-reorg; each needs a regression test in the PR that ports its module.

### 7.1 deriv

```cpp
namespace nxx::deriv {
template<int Order, int Accuracy, std::size_t Points> struct stencil {       // integer data: exact for any scalar   [prototyped]
    static constexpr int order = Order, accuracy = Accuracy;
    std::array<int, Points> offset, weight; int denominator; };
inline constexpr stencil<1, 2, 2> central_1_2{{-1, 1}, {-1, 1}, 2};
inline constexpr stencil<1, 4, 4> central_1_4{{-2, -1, 1, 2}, {1, -8, 8, -1}, 12};
inline constexpr stencil<2, 2, 3> central_2_2{{-1, 0, 1}, {1, -2, 1}, 1};
inline constexpr stencil<2, 4, 5> central_2_4{{-2, -1, 0, 1, 2}, {-1, 16, -30, 16, -1}, 12};
// forward_1_1/1_2/1_3, backward_1_1/1_2/1_3, forward_2_1/2_2, backward_2_1/2_2 (one-sided, for domain edges)

// Step specifications are scalar-agnostic values, resolved per stencil <O, A> and per x inside diff.  [prototyped mechanism]
// h = factor * max(|x|, typical); typical defaults to |x| (x != 0) and to 1 at x == 0 (D32); then h = (x + h) - x.
struct optimal {};                                                   // factor = eps^(1/(O+A)) as an exact power of two
template<class P> struct relative { tolerance<P> factor; std::optional<tolerance<P>> typical; };  // relative{1e-5, 1.0}: h = 1e-5·max(|x|, 1)
template<class P> struct noise    { rel_tolerance<P> eps_f; std::optional<tolerance<P>> typical; }; // factor = eps_f^(1/(O+A)) via frexp/ldexp
template<class P> struct absolute { tolerance<P> h; };
// (CTAD guides: relative{1e-5, 1.0} -> relative<double>, noise{1e-10} -> noise<double>)
// [spike] relative stores the typical value and a flag, and typical() returns the std::optional: MSVC 19.51 rejects a
// consteval constructor that initialises a std::optional member. noise follows the same shape when it lands.

template<class F, real T, int O = 1, int A = 2, std::size_t N = 2, class H = optimal>
constexpr auto diff(const F& fn, T x, const stencil<O, A, N>& s = central_1_2, H h = {})
    -> std::expected<T, fault<callback_error_t<F, T>>>;               // default template args: diff(f, x) deduces   [prototyped]
constexpr auto central(fn, x, h = optimal{}); forward(...); backward(...); second(...);   // named conveniences
template<real T> struct derivative_estimate { T value; T error; T step; };
template<class F, real T, class Lo = decltype(central_1_2), class Hi = decltype(central_1_4), class H = optimal>
constexpr auto diff_with_error(const F& fn, T x, Lo lo = central_1_2, Hi hi = central_1_4, H h = {})
    -> std::expected<derivative_estimate<T>, fault<callback_error_t<F, T>>>;   // error = |D_hi - D_lo|; 2 extra evaluations
template<class F, real T, class H = relative<T>>
constexpr auto ridders(const F& fn, T x, H h0 = relative{0x1p-3}, max_iterations levels = 10,
                       T con = T(7) / T(5), T safe = T(2))
    -> result<derivative_estimate<T>, callback_error_t<F, T>>;   // <= 2*levels evaluations; best tableau entry + its error;
                                                                 // stops when the error grows by `safe` (success, best-so-far)
constexpr auto mixed(const F2& fn, T x, T y, /* mixed_5 | mixed_9 */ H hx = {}, H hy = {})
    -> std::expected<T, fault<callback_error_t2<F2, T>>>;        // per-axis relative steps, eps^(1/4) / eps^(1/6)
template<class F, class S = decltype(central_1_2), class H = optimal>
constexpr auto derivative_of(F fn, S s = central_1_2, H h = {});  // -> derivative_fn   [prototyped]
template<class S = decltype(central_1_2), class H = optimal> struct numeric {   // policy: bind(f) -> derivative_fn   [prototyped]
    S s = central_1_2; H h = {};
    template<class F> constexpr auto bind(const F& fn) const { return derivative_of(fn, s, h); } };
// Vector-valued derivatives (gradient, Jacobian, Hessian) return Eigen types and live in multiroots/derivatives.hpp (§7.5).
}
```

| | Items |
|---|---|
| Keep | The stencil formulas, as data. The mixed 5- and 9-point stencils (dev-reorg). The fixed `DerivativeFunctor` idea, now `derivative_of`. |
| Fix | Duplicate stencils (`Order1CentralRichardson` ≡ `Order1Central5Point`; "ForwardRichardson" is an honestly named `forward_1_3`). The sign-blind step `max(h, h·x)` becomes `factor·max(\|x\|, typical)`. **Per-stencil optimal step** (v1's second-derivative and mixed overloads with default steps are wrong by O(1)). **Relative by default**: at x = 1e-3, `noise(1e-10)` with a floor of 1 gives h = 4.6e-4, a 27 % error in the derivative of f(x) = 1/x (measured); the relative default gives h = 4.6e-7. `validateStepSize` no longer throws. |
| Add | `diff_with_error` and `ridders` with error estimates (a low-confidence signal for the caller); `noise` steps for functions computed with limited precision or by inner iterative solvers; fallible callbacks; `derivative_of` + `numeric` policy. Complex-step is deferred. |
| Move | gradient, Jacobian and a correct Hessian go to `multiroots/derivatives.hpp`, because they return Eigen types; scalar deriv stays core-only. |
| Drop | Template-template `diff<ALGO>`; `IsDiffSolver`; the `requires(!poly::IsPolynomial)` back-edge; gcem. |
| Must not port | the discarded argument in `derivativeOf`; non-const `operator()`; the sign-blind step; dev-reorg's √ε step for every order (the sin′(1) error rises from 2e-12 to 1.2e-9); absolute √ε in `mdiff` (errors of 0.57); examples that dereference a result without checking it. |
| Genericity | real T (`float`, `double`, `long double`, multiprecision); steps in `T`; no `pow` at run time (factors are exact powers of two). Measured: sin′(1) error −1.86e-13 with `central_1_2`, −5.34e-14 with `central_1_4` **[prototyped]**. |

### 7.2 roots (1-D)

**Common rules for bracketing solvers:**
- The midpoint is overflow-safe (`math::midpoint`). The naive `lo + (hi−lo)/2` reported a false `resolution_limit` success at x = −1.7e308 on [−1.7e308, 1.7e308] (measured).
- Endpoint samples may be ±inf (`evaluate_sample`), and the solver bisects while an endpoint value is infinite, so log(x) on [0, 2] works; a design that rejects infinite samples fails there with `non_finite_value`.
- **Pole check in `finish()`.** After a criterion or `resolution_limit` stop, if min(|f(a)|, |f(b)|) > max(|f(lo₀)|, |f(hi₀)|) (the residual grew while the bracket shrank), the solve fails with `errc::sign_change_not_root`, carrying the final enclosure as best. Without it, tan on [1, 2] is reported as a root with |f| = 1.2e15 (measured on the spike's bisection and brent, which both stop at x = 1.5707963267948974; the prototype measured 6e15).
  - **Implemented in the spike** (`detail::pole_check`, shared by bisection and brent), with two refinements found by the review. First, a non-finite fx at the returned point is always a pole. Second, the reference max(|f(lo₀)|, |f(hi₀)|) is taken over the *finite* initial samples. An infinite end sample (1/x on [−1, 0]) would otherwise make the reference infinite and switch the check off, and 1/x on [−1, 0] was reported as a root with |f| = 1.1e15.
  - When both initial samples are infinite, only the finiteness test applies (1/x³ on [−1e-200, 1e-200] is caught by it; tested).
  - **Known limit:** the reference is the *larger* initial sample, so a large but finite one hides a pole near the other end. bisection on 1/x over [−1e-20, 1] (reference 1e20) returns a success at |f| = 1.1e15 (measured). Comparing with the smaller sample, or side by side, would reject legitimate roots when the bracket starts next to a near-root and the tolerance is loose. It is a heuristic by design (§11); callers who need a residual guarantee add `&& f_tol{…}`, as for jump discontinuities.
- Jump discontinuities (for example a step from −1 to +1) still return the sign-change location with its fx. Callers who need a residual guarantee add `&& f_tol{…}`.

| Solver | Input | Phase | Status | Notes |
|---|---|---|---|---|
| `bisection` | bracket-like | 3 | Keep/Fix | 1 evaluation per step. Exact zero → `exact_zero`; unsplittable → `resolution_limit` **[prototyped]**. **Representation-space midpoint** for `float`, `double`, and `long double` where it is binary64 or binary128 (x87's padding bits break constexpr `bit_cast`), when the bracket straddles 0 or spans more than 2 binades: bisection on the monotone integer image of the bit pattern gives at most 64 steps for double at full relative accuracy at any magnitude. Value-space bisection needs about 665 halvings for a root at 1e-200 and exhausts its budget. Other types use the arithmetic midpoint. |
| `illinois` | bracket-like | 3 | Fix | **Anderson–Björck is the default variant**; plain false position only as `variant::plain` |
| `ridders` | bracket-like | 3 | Keep/Fix | guard the `sqrt`, falling back to bisection; 2 evaluations per step |
| `brent` | bracket-like | 3 | Add, provisional **default** | zeroin; tolerance = a width criterion (default `floored_width{}`); 7 iterations and 9 evaluations on x²−2 over [1, 2]; constexpr **[prototyped]** |
| `rtsafe` | sign_bracket + derivative | 3 | Add | a Newton step is accepted only if it stays inside the enclosure and makes progress; otherwise bisect. Used by `solve(f, df, x0)` |
| `secant` | guess | 3 | **Fix: derivative-free** | second point x0 + 2⁻¹⁰·max(\|x0\|, typical), or x0 − h where x0 + h overflows or a projection pins it to x0 **[spike]**; a caller-supplied `x1` is phase 3; a flat secant → `stalled` **[prototyped]**; `step_tol<7,10>`; budget 50 |
| `newton` | guess + derivative source | 3 | Keep/Fix | zero or non-finite derivative → error; `step_tol<3,5>`; budget 30; optional projection **[prototyped]** |
| `expand` / `scan` / `subdivide` | window or guess (+ limits) | 3 | Merge the six searchers | success `solution<sign_bracket<T>>` **is** the input of every bracketing solver (R-A6) **[prototyped: `expand` through the driver]**. **No configurable stop criterion** (intrinsic stop + budget only), because `estimate()` on a non-bracketing state would build an invalid `sign_bracket`: `with_stop` is deleted, and the public `rebuild(options)` accepts only `never{}` **[spike]**. Rules for `expand`: **geometric growth** (lo/g, hi·g, g = 1.6) when lo > 0, mirrored when hi < 0, additive otherwise; **expand only the end with smaller \|f\|** (1 evaluation per step); a NaN at a trial point marks a domain limit on that side, so it backtracks halfway to the last finite point and stops expanding that side (a symmetric additive search walked sqrt(x)−10 into x < 0, measured); it returns the **tightest** sign-changing pair among cached samples. **Spike:** endpoints saturate at ±max (the width is computed as 2·(hi/2 − lo/2)), because an infinite endpoint gave a `sign_bracket` that bisection called unsplittable, a false `resolution_limit` success at \|f\| = 4.5e307. When both ends are at ±max, `expand` fails with `stalled`. **Known limit, for phase 3:** geometric growth moves the end nearer 0 toward 0 without ever crossing it, so x + 1 from [1, 2] exhausts its budget. Phase 3 switches that end to additive growth once it is within one window width of 0, together with the NaN backtrack, which is needed exactly where crossing 0 leaves a domain (log, sqrt). The NaN backtrack, a guess as start, `scan` and `subdivide` are phase 3. |
| `toms748` | bracket-like | 8 (optional) | Add | attributed port of Boost.Math `toms748_solve` (BSL-1.0 notice kept); tolerance `floored_width`. A strong general method, usually fewer evaluations than Brent on smooth f; checked against Boost.Math as an oracle |
| `itp` | bracket-like | 8 (optional) | Add | `itp_params{kappa1 = relative{0.2} /* kappa1 = 0.2/(b0-a0) at init */ or absolute{k}, kappa2 = 2 /* 1 <= kappa2 < 1+phi */, n0 = 1}`. ε_ITP = half the tolerance's threshold at (a₀, b₀), frozen at init; the ⌈log₂(w₀/2ε)⌉ + n₀ worst-case bound is relative to it; ⌈log₂⌉ via `frexp` (exact). Best worst-case bound among the bracketing methods |
| `halley` | guess + f′, f″ (+ optional sign_bracket) | 8 (optional) | Add | accepts one callable returning `{f, f′, f″}`; with a bracket it falls back to bisection when a step leaves it (as rtsafe does) |
| `steffensen` | guess | 8 (optional) | Fix (low priority) | no Newton first step, no throw; deleted if the corpus shows no benefit |

**Open-method safeguards** (newton, secant; halley and steffensen when added):
- **Progress window** in the state (default `progress_window{4}`, configurable, and can be disabled) **[sketch]**:
  - `stalled` when the best |f| has not improved for `window` steps **and** the step length has not shrunk. This catches cycling, such as Newton on x³−2x+2 from 0, which alternates 0 ↔ 1.
  - `diverged` when |x| and |dx| have both grown monotonically over the window.
  - The step-shrink guard prevents false stalls on Newton's long approach from a poor start: x²−2 from 1e-3 overshoots to about 1000, then needs about 10 steps of halving.
- **Before each evaluation**: a non-finite proposed x → `diverged`; |dx| is **capped** at `max_step` (default 2^10·max(|x|, typical)). A binding cap feeds the window's divergence rule; it is not an immediate failure.
- A `root_estimate` is built from x and f(x) at least **[spike]**. Its constructor has no default for fx, because an open method starts from an estimate's fx without evaluating f again. Before the spike's review, `root_estimate<double>{1.0}` gave an `exact_zero` success at x = 1 without a single evaluation.
- `root_estimate::uncertainty` = |x_k − x_{k−1}|, or inf for an estimate with neither an enclosure nor a step: the best endpoint of a no-sign-change or search failure, or an open method that stops at its start on an exact zero (secant and Newton agree) **[spike]**. `exact_zero` means "f evaluated to 0", not "full accuracy". Newton on the expanded (x−1)³ stops `exact_zero` at an error of 4.7e-6, the ε^(1/3) conditioning limit (measured).

Other roots items:
- **Default choice (D29):** the phase-3 corpus records mean and worst-case evaluation counts for bisection, illinois, ridders and brent, and fixes the bracketing default. Phase 8 re-runs it with toms748 and itp; the default changes only with a documented behaviour note. TOMS748 usually needs fewer evaluations on smooth f; ITP has the best worst-case bound.
- **Drop:**
  - `fsolve`/`fdfsolve`/`search` and their template-template drivers;
  - the CRTP bases, `*Traits`, `RootErrorImpl<T>`, `StopToken`/`ArgTypes`/`makeToken`, `ResultProxy`/`SearchResult`;
  - `PolishingIterData`'s history vector;
  - complex support.
- **Must not port:**
  - residual-only convergence; product sign tests;
  - the strict `< 0` test that rejects endpoint roots;
  - no sign-change check (x²−5 on [3, 4] "converged" to 3.99999997; dev-reorg's `fsolve<Bisection>` on a bracket [5, 6] without a sign change silently returned 5.998, measured);
  - non-convergence returned as a plain value (dev-reorg's `fdfsolve<Secant>` on x²+1, which has no real root, returned 3.89 from 1.5 through `.result()`, measured);
  - Secant and Steffensen needing df;
  - the `eps*x + eps/2` tolerance, negative for x < −0.5;
  - silently ignored single eps/maxiter arguments; swapped (maxiter, eps) returning a point outside the bracket;
  - `SearchStopToken` default-constructing f;
  - Newton returning `inf` or −0.87 on x²+1;
  - `validateBounds` throwing out of `fsolve`.
- **Genericity:** real T. The stale master tests' 8-function table becomes the corpus.

### 7.3 optimize (1-D; N-D minimisation is the v2.1 module `multimin`)

```cpp
namespace nxx::optimize {
template<real T> struct extremum { T x; T fx; std::optional<bracket<T>> enclosure; };   // named x AND f(x)
template<real T> class min_bracket;       // a < b < c, f(b) <= min(f(a), f(c)); from bracket_minimum (mnbrak)
class golden;                             // dev-reorg step, pure, 1 evaluation per iteration; Brent's intrinsic width test
class brent_min;                          // parabolic + golden; tol1 = rel*|x| + abs/3, test |x - m| <= 2*tol1 - (b - a)/2
class bracket_minimum;                    // dev-reorg AutoSearch, downhill-direction bug fixed -> solution<min_bracket<T>>
template<class S> constexpr auto maximizing(S s);   // solver adaptor: minimise -f; fx restored on success AND failure
template<class F, class In, class S = brent_min> constexpr auto maximize(const F& fn, const In& in, S s = {});   // == maximizing(s)(fn, in)
template<class F, class In, class S = brent_min> constexpr auto minimize(const F& fn, const In& in, S s = {});
}
```

- **Termination.** Golden and Brent-min end on Brent's intrinsic test, with `rel = root_eps<T>(1, 2)` (0x1p-26 for double) and `abs = rel·typical`. Never `x_tol` on the best point, which has the same repeated-endpoint defect as bracketing roots. `f_tol` is ill-formed on minimisers.
- **Keep:** golden section and the Brent minimiser (dev-reorg `OptimBracket.hpp` L136-334; rewrite Brent from Brent 1973 or keep the BSL-1.0 attribution); `AutoSearch` → `bracket_minimum`.
- **Add:** `maximizing` and `maximize`, which return `extremum{x, fx}` with fx = f(x) restored, so the argument of an extremum is never confused with its value, nor −f with f (R-A9); `minimizer_of`; `newton_min` with analytic f′ and f″ and a curvature check (later). Nelder–Mead, BFGS/L-BFGS and nonlinear CG are planned for v2.1 in a separate module, `multimin`, because they need Eigen and `optimize` stays core-only (§5.2, §10.5).
- **Drop:**
  - 1-D gradient descent;
  - Newton as a finite difference of a finite difference;
  - the Minimize/Maximize template-template mode;
  - a `maximize(solver, f, in)` signature, which breaks D4's call order and is not a solver value;
  - the Vandermonde "parabola vertex via polysolve".
- **Must not port:** the `eps*x` terminator (a minimum at −2 took 226 iterations against 40 at +2); `foptimize_impl`'s hard-coded `IterData<size_t, double>`; the non-self-contained `Optim.hpp`.

### 7.4 poly

```cpp
namespace nxx::poly {
template<class T> class polynomial {          // T real or std::complex<real>; c[i]*x^i; canonical: empty == zero, else back() != 0
public:
    polynomial() = default;                                             // the zero polynomial
    polynomial(std::initializer_list<T>);                               // NORMALISING: trims exact trailing zeros
    static auto make(std::span<const T>) -> std::expected<polynomial, errc>;   // also rejects non-finite coefficients
    std::optional<std::size_t> degree() const noexcept;                 // nullopt for the zero polynomial
    constexpr T operator()(const T& x) const noexcept;                  // Horner; total
    std::array<T, 3> eval_with_derivatives(const T& x) const noexcept;  // p, p', p''
    polynomial derivative() const;                                      // TOTAL: d/dx(c) is the zero polynomial
    polynomial antiderivative(const T& c0 = T(0)) const;
    friend polynomial operator+(const polynomial&, const polynomial&);
    friend polynomial operator-(const polynomial&, const polynomial&);  // correct when deg(lhs) < deg(rhs)
    friend polynomial operator*(const polynomial&, const polynomial&);
    friend polynomial operator*(const T&, const polynomial&);
};
template<class T> std::optional<nonzero<T>> nonzero_of(polynomial<T>);
template<class T> struct divmod_result { polynomial<T> quotient, remainder; };
template<class T> divmod_result<T> divmod(const polynomial<T>&, const nonzero<T>&);      // TOTAL: no throw, remainder reduced
template<real T> class linear; template<real T> class quadratic; template<real T> class cubic;   // leading coeff != 0
template<real T> using quadratic_roots = std::variant<real_pair<T>, complex_conjugates<T>>;      // never NaN-as-success
template<real T> using cubic_roots     = std::variant<three_real<T>, one_real_two_complex<T>>;
template<class T> auto roots(const nonzero<T>&, aberth<real_t<T>> = {})
    -> result<root_set<std::complex<real_t<T>>>>;                                          // (heap: std::vector)
}
```

- **Aberth–Ehrlich, specified:**
  - starting points on a circle (or on Newton-polygon radii, as in MPSolve) with an angular offset;
  - per-root stop when |p(z)| ≤ γ₂ₙ·Σ|aᵢ||z|ⁱ (Horner's running error bound), freezing converged roots;
  - for real coefficients, Im z := 0 when |Im z| is within that bound; real roots are polished in real arithmetic; conjugates are paired;
  - documented accuracy ε^(1/m) for m-fold roots.
- **Closed forms:** Kahan's fma-based discriminant and q = −(b + sign(b)·√disc)/2 for the quadratic, then one or two Newton steps on the original polynomial for every closed-form root.
- **Keep:** Horner, the O(nm) product, the stable quadratic, the cubic, `from_roots`, and formatting (as `std::formatter`).
- **Fix:**
  - `operator-` when deg(lhs) < deg(rhs), with the test that locks in `{-1,-1,-1,8}`;
  - `divide` (total, remainder reduced);
  - trimming by **exact zero** (not dev-reorg's hidden 1.2e-4 / 0.01 tolerances);
  - the derivative of a constant;
  - distinguishing zero from a constant;
  - `linear`/`quadratic`/`cubic` accepting higher degrees;
  - `sortRoots` (a strict weak ordering).
- **Add:** `antiderivative`, `compose`, scalar operations, `eval_with_derivatives`, Aberth–Ehrlich, and an optional closed-form quartic.
- **Not a library feature:** a companion-matrix eigenvalue solver. It would need Eigen, and poly stays core-only (`poly:""`). Instead, when `NUMERIXX_WITH_LINALG` is ON, the poly tests use Eigen's `EigenSolver` on the companion matrix as an oracle for Aberth (§9.1), the same role `HybridNonLinearSolver` has for multiroots.
- **Drop:** Laguerre with `random_device`; `polysolve`'s silent real/complex switch; the dependency on `roots::fdfsolve`; the compound mutators.
- **FLAG:** the invariant "leading coefficient ≠ 0" is exact. A leading coefficient of 1e-300 is legal; near-degeneracy is the algorithm's problem, not the type's.

### 7.5 linalg and multiroots

```cpp
namespace nxx::linalg {                     // Eigen 5.0.1 behind a thin facade
template<real T, int N = Eigen::Dynamic>              using vector = Eigen::Matrix<T, N, 1>;
template<real T, int R = Eigen::Dynamic, int C = R>   using matrix = Eigen::Matrix<T, R, C>;
// Every function returns a concrete Eigen type or std::expected of one: never an expression template through auto.
template<real T, int N>                                                   // [prototyped mechanism]
auto lu_solve(const matrix<T, N, N>& a, const vector<T, N>& b) -> std::expected<vector<T, N>, errc>;
    // PartialPivLU: dimension_mismatch; non_finite_input (allFinite); rcond() <= n*eps -> singular (PartialPivLU never
    // reports singularity itself); non-finite solution -> singular
template<real T> struct least_squares_solution { vector<T> x; Eigen::Index rank; };   // [sketch]
template<real T, int R, int C>
auto qr_solve(const matrix<T, R, C>& a, const vector<T, R>& b) -> std::expected<least_squares_solution<T>, errc>;
    // ColPivHouseholderQR: square, over-determined (least squares) and rank-deficient (basic solution + numerical rank)
template<real T, int N>                                                   // [sketch]
auto cholesky_solve(const matrix<T, N, N>& a, const vector<T, N>& b) -> std::expected<vector<T, N>, errc>;
    // LLT; info() != Success (not positive definite) -> singular
template<class V> struct vector_traits;   // std::array<T,N> -> vector<T,N>; std::vector<T> -> vector<T>; Eigen column vectors: identity
                                          // scalar, storage, to_storage(const V&), from_storage(const storage&)   [prototyped mechanism]
}
namespace nxx::multiroots {
template<class V> struct system_estimate { V x; V fx; scalar_of<V> merit; scalar_of<V> step_norm; std::optional<scalar_of<V>> rcond; };
struct default_hooks {                            // a template parameter: no type erasure
    template<class X> constexpr X project(X x) const noexcept { return x; }            // box or domain constraints
    template<class X> constexpr X typical(const X& x) const noexcept;                  // D_x scaling (MINPACK diag)
    template<class X> constexpr X weights(const X& x, const X& fx) const noexcept;     // optional (default 1): frozen per iterate;
};                                                                                     //   merit 0.5*sum (w_i f_i)^2 and f_tol's norm
template<class Opt = /* x_tol on the full step || weighted f_tol; g_tol */, class Jac = forward_difference, class Hooks = default_hooks>
class newton;          // damped Newton + Armijo backtracking; value state over linalg::vector<T, N>   [prototyped: in-house and Eigen storage]
class broyden; class dogleg;
// derivatives.hpp: gradient_of, jacobian_of (forward default), hessian_of (true symmetric mixed partials); Eigen results
}
```

- **The facade.**
  - Eigen 5.0.1 is fetched by CPM (`NUMERIXX_WITH_LINALG`, default ON). It needs no BLAS/LAPACK and was verified under em++ 6.0.8, including `-fno-exceptions`; the prototype's facade and N-D Newton on Eigen storage ran on all 9 configurations **[prototyped]**.
  - Measured on the prototype's `lu_solve` **[prototyped]**: [[1,2],[3,4]] → (−4, 4.5); the singular [[1,2],[2,4]] → `singular`; a 2×2 against a 3-vector → `dimension_mismatch`; the same damped Newton runs on `Vector2d` and `VectorXd`.
  - **Singularity.** A threshold of n·ε·max|a| falsely flags a well-conditioned Jacobian whose columns differ by 1e10 in scale (measured). The facade instead compares Eigen's `rcond()` estimate against n·ε and checks the solution's finiteness. `rcond` is not scale-invariant either; if the corpus (§9.2) shows badly scaled Jacobians flagged singular, the facade equilibrates (row/column scaling) before the rcond test. In Newton, the `typical` scaling (below) helps as well.
  - **Concrete types only.** Facade functions and solver code never bind `auto` to an Eigen expression (§5.3), because Eigen's expression templates dangle under `auto`.
  - **Multiprecision.** `cpp_bin_float_50` scalars were verified with Eigen. `adapters/multiprecision_linalg.hpp` includes `<boost/multiprecision/eigen.hpp>`. Eigen's `operator<<` with `cpp_bin_float` fails to compile on Boost 1.92, so tests never stream such matrices.
  - Eigen takes no part in constexpr tests, so N-D solves are not `constexpr`. Allocation failure in dynamic storage is not converted (as elsewhere).
  - `Eigen::HybridNonLinearSolver` is a test oracle, never a code path.
- **FLAG: linalg as requested.** Eigen, as you asked, behind a facade:
  - (i) **measured compile cost** **[prototyped]**: +0.2–1.0 s per TU to include `<Eigen/Core>` + `<Eigen/LU>`, and +2.5–6.4 s over the prototype's in-house TU once LU and N-D Newton are instantiated. It stays in TUs that include linalg or multiroots: scalar modules and the umbrella header never include Eigen (§5.2), and a user of the scalar modules only can skip the download (`NUMERIXX_WITH_LINALG=OFF`);
  - (ii) **the expression-template/`auto` hazard** is handled by concrete return types;
  - (iii) N-D solving loses `constexpr`, and dynamic sizes allocate, both acceptable;
  - (iv) the only in-house linear algebra is the O(n) tridiagonal solvers for splines (§7.7), which are algorithms rather than a library need.
- **multiroots, damped Newton (phase 5):**
  - best iterate on every exit;
  - `singular`, `line_search_failed`, `stalled` and `budget_exhausted` are distinct;
  - no success at maxiter; no `std::cout`;
  - rounding-level step acceptance (‖dx‖ ≤ 16ε·max(1, ‖x‖)) avoids a false stall at merit ≈ 1e-32 **[prototyped]**;
  - **D_x scaling** from `hooks.typical`;
  - componentwise `x_tol` measured on the **full** Newton step ‖J⁻¹F‖_D (a damped step λ ≪ 1 must not declare convergence);
  - λ < λ_min (1e-10) → `line_search_failed`;
  - **`errc::local_minimum`** when ‖JᵀF‖_D ≤ g_tol·max(merit, 1) while f_tol is unmet, so a chain can react (for example restart from another guess);
  - Armijo α = 1e-4 with quadratic/cubic backtracking, and `stpmax = 100·max(‖x‖_D, n)` (Dennis & Schnabel A6.3.1);
  - the FD Jacobian uses h_j = √ε_f·max(|x_j|, typical_j)·sign(x_j), rounded as (x_j + h) − x_j, and takes a **backward** difference when `hooks.project(x + h·e_j) ≠ x + h·e_j`. A floor-1 step perturbs x_j = 1e-10 by 1.5e-8 (a 99 % error for ln x) and can step out of the box;
  - state: a value over `linalg::vector<T, N>`, fixed-size when N is known at compile time, dynamic otherwise (§3.2).
- **multiroots, Broyden and dogleg (phase 5):**
  - Broyden keeps the QR factors in its state, with a rank-1 update (MINPACK-style Givens update written over Eigen matrices, since Eigen has no QR update), and refreshes the FD Jacobian after a failed line search or 2 steps with < 10 % merit decrease.
  - Powell dogleg follows the MINPACK `hybrd` trust-region rules, with the trust region in the scaled norm.
  - `system_of(f1, …, fN)` for fixed N.
  - **`multiroots::solve(F, x0)` defaults to dogleg (hybrid) with Broyden updates**, the robust default of MINPACK, GSL and Eigen, once it passes the corpus; until then it is damped Newton with an FD Jacobian. Damped Newton stays available as `multiroots::newton`.
  - Levenberg–Marquardt belongs to nonlinear least squares, planned for v2.1 in the module `fit` (§10.5).
- **Hooks.** `project`, `typical` and the optional `weights` are the whole hook set. A caller-defined acceptance test goes through a `custom{λ}` criterion (§6.8).
- **multiroots, Drop:**
  - `MultiFunction`/`MultiFunctionArray` (a system is now one callable `V → V` or `V → expected<V, E>`);
  - `SteepestDescent` as a solver (its −JᵀF survives inside dogleg);
  - `ContainerTraits`;
  - Blaze; LAPACK.
- **Must not port:** success at maxiter; `std::cout`; the float cast in `MultiFunction`; `hessian()` returning only the pure second partials; three residual evaluations per loop; the Vandermonde solve plus `polysolve` for a parabola vertex.

### 7.6 integrate

```cpp
namespace nxx::integrate {
template<real T> struct integral { T value; T error; };
using min_level = detail::refined<tag::quad_level, std::uint8_t>;    // 1..30: (1u << level) is defined BY TYPE
struct romberg { min_level min = 4; max_level max = 12; };            // 4097 evaluations max (a cap of 20 would allow ~1.05M)
template<std::size_t MaxSegments = 64> struct gauss_kronrod15 {};     // adaptive G7K15, fixed-capacity max-heap worklist; DEFAULT
struct tanh_sinh { max_level max = 8; };  struct exp_sinh { /*...*/ };  struct sinh_sinh { /*...*/ };
template<std::size_t N, class F, real T>
constexpr auto gauss_legendre(const F& fn, T a, T b) -> std::expected<T, fault<callback_error_t<F, T>>>;  // plain function: no error estimate
template<class F, real T, class Rule = gauss_kronrod15<>>
constexpr auto quad(const F& fn, T a, T b, Rule = {}, abs_tolerance<T> abs = /* 0 */, rel_tolerance<T> rel = /* root_eps<T>(1,2) */)
    -> result<integral<T>, callback_error_t<F, T>>;   // a == b -> {0, 0}; a > b -> negated; worklist full -> budget_exhausted + best
template<class F, class Rule = gauss_kronrod15<>> constexpr auto integral_of(F fn, Rule r = {});      // FIXED: keeps f
template<class F, real T, class Rule> constexpr auto antiderivative(F fn, T a, Rule r);
}
```

- **QUADPACK rules for G7K15:**
  - per segment `err = resasc·min(1, (200·|K−G|/resasc)^1.5)`, floored at `50·ε·resabs`;
  - accept when Σerr ≤ max(abs, rel·|ΣI|);
  - always bisect the segment with the largest error;
  - roundoff detection → `stalled` with best.
- **`abs` may be 0** (purely relative, the QUADPACK default). The roundoff floor lets zero-valued integrals (sin over a period) converge.
- **Tabulated Kronrod constants.** Nodes and weights (G7K15, optionally G10K21) are shipped with ≥ 36 significant digits as literal constants (the §3.5 exception), with `static_assert(digits10<T> <= 33)`. Higher precision routes to tanh-sinh or Romberg. Laurie's algorithm would need an eigensolver, and integrate stays core-only, so the nodes are not computed. Gauss–Legendre nodes *are* computed for T (Newton on Legendre polynomials).
- **tanh-sinh** computes abscissae as distances to the endpoints (the complement form; optionally `f(x, xc)`) and truncates the tail where the abscissa rounds to the endpoint or f is non-finite in the last level. Otherwise 1 − x == 1 samples the singular endpoint. exp-sinh and sinh-sinh handle `semi_infinite` and `whole_line`.
- **Chains:** `first_of(gauss_kronrod15<>{}, tanh_sinh{}, romberg{})` is meaningful because the three fail differently. Romberg comes last: placed first with a level cap of 20, it has a 1M-evaluation worst case.
- **Keep/Fix:** Trapezoid, Romberg and Simpson become one family over the shared trapezoid sequence. It is generic over the scalar, reuses samples, and uses a relative + absolute test with `min_iterations`. That fixes false convergence on sin²(8πx), removes the `1 << m_iter` UB by type, and makes a > b and a == b legal.
- **Drop:** `boost::multi_array`, `IntegrationBase`/`IntegrationTraits`, public mutable members, the throwing `IntegrationFunctor`.
- **Must not port:** `integralOf` discarding f; "lower ≥ upper throws"; re-running `init()` after iterating.

### 7.7 interpolate

```cpp
namespace nxx::interpolate {
template<real T> class knots;                       // size >= 2, finite, strictly increasing; segment(x) never reads out of bounds
struct reject_outside {}; struct clamp_to_domain {}; struct linear_extrapolation {};   // the out-of-range policy is a TYPE
template<real T, class Outside = reject_outside> class linear;
template<real T, class Boundary = natural, class Outside = reject_outside> class cubic_spline;   // natural | clamped{d0,dn} | not_a_knot | periodic
template<real T, class Outside = reject_outside> class pchip;          // monotone Hermite (dev-reorg's "Steffen", correctly named)
template<real T, class Outside = reject_outside> class steffen;        // true Steffen (1990)
template<real T, class Outside = reject_outside> class barycentric;    // Floater-Hormann rational, d = 3 default; pure polynomial
                                                                       // barycentric only for Chebyshev nodes
// make_linear / make_cubic_spline / make_pchip / make_steffen / make_barycentric(xs, ys [, boundary] [, outside])
//   -> std::expected<interp, errc>; policies deduced from tag arguments; coefficients computed EAGERLY; binary search
// operator()(x): reject_outside -> std::expected<T, errc{out_of_domain}>; clamp / extrapolate -> T
// linear_extrapolation = the C^1 continuation with the interpolant's end derivative
// cubic_spline, pchip, steffen expose .derivative() -> Newton on an interpolant uses exact slopes
namespace detail {   // tridiagonal.hpp: in-house O(n) solvers; interpolate depends only on core
    // thomas (diagonally dominant: natural, clamped); pivoted tridiagonal LU, dgtsv-style (not-a-knot);
    // Sherman-Morrison on top of thomas (periodic, cyclic tridiagonal)
}
}
```

- **Keep/Fix:**
  - Linear (dev-reorg `lower_bound`, fixing x == x0);
  - splines: natural and clamped via Thomas; **not-a-knot via a dgtsv-style pivoted tridiagonal LU**; **periodic via Sherman–Morrison**, since neither of the last two is diagonally dominant or tridiagonal;
  - the monotone Hermite, renamed and fixed at the last knot;
  - the `interpolationOf` idea that an interpolant is a callable.
- **Drop:** `makepoly` (Vandermonde + LAPACK), `InterpBase` CRTP, throwing constructors and `operator()`, `mutable` caches, Blaze storage. Plain barycentric Lagrange on equispaced knots is dropped too (Runge-unstable).
- **Must not port:** dev-reorg's out-of-bounds read of `Steffen(x_last)`; the uninitialised `d[n]`; master's ODR violation in `evaluateSpline`; `LinearInterp` rejecting x == x0; thread-unsafe lazy caches; no duplicate-x check.
- Storage is `std::vector` and copies deep.

---

## 8. FXT

### 8.1 What Numerixx uses

- **`<numerixx/pipes.hpp>` only**, in target `numerixx::pipes`:
  - `monads/Expected.hpp` (`operator|`) and the adaptors `Transform`, `AndThen`, `OrElse`, `TransformError`, `ValueOr`, `Tap`, `Match` **[prototyped]**;
  - `namespace nxx { using fxt::operator|; }`, with no other `operator|` declared in `nxx`.
- **Available to users:** `curry`, `compose`, `FXT_LIFT`, `zip`/`with` (validating several refined inputs at once), `sequence`/`traverse`.
- **Never in core paths:** `operator||`, `failure`/`result<T>`, `lazy`, `attempt` (only in `fn::catching` under `__cpp_exceptions`), `value()`, `immutable`, the enums, `Format.hpp`.
- **Public signatures use `std::expected`**, never the configuration-dependent `fxt::expected` alias. The build warns if a parent configured FXT for `tl::expected`.
- **ADL hazard:** once FXT-4 exists, both `fxt::first_of` and `nxx::first_of` exist. Numerixx code always qualifies its calls.

### 8.2 Upstream changes to FXT, prioritised

| Pri | ID | Change | Signature / location | When |
|---|---|---|---|---|
| **P0** | FXT-1 | `throw 0;` → `std::unreachable();` in the `expected_like`/`optional_like` probes; `#if __cpp_exceptions` guards for Attempt, Failure, Lazy, Format and the enums | `concepts/IsExpected.hpp:87-88`, `concepts/IsOptional.hpp:76-77` | before phase 1. P0 because it is a 2+2-line diff; it gates only the `-fno-exceptions` legs that use pipes. **Verified necessary and sufficient** for transform, and_then, value_or, match and tap under `-fno-exceptions` on GCC, Clang, em++ and clang-cl `/EHs-c-` **[prototyped]**; the patch is `docs/redesign/prototype/fxt-1.patch`. **Status:** the probe fix is troldal/FXT#1, pinned by Numerixx since the spike; the `__cpp_exceptions` guards are still to do (`numerixx::pipes` includes none of those headers) |
| **P1** | FXT-2 | CMake hygiene: `if(NOT COMMAND CPMAddPackage)` or a hash-pinned CPM; TL repos only when their options are ON; `target_compile_features(fxt INTERFACE cxx_std_23)`; `install(EXPORT)` + `fxtConfig.cmake`; semver tags | `CMakeLists.txt` | phases 0–1 (not blocking) |
| **P2** | FXT-3 | `std::expected<void, E>` in `expected_like` (or document `fxt::unit`) | `IsExpected.hpp` | phase 1 |
| **P3** | FXT-4 | Generic lazy alternatives as **named, assignable class templates** | `template<class Policy, class... Fs> constexpr auto first_of_with(Policy p, Fs... fs); template<class... Fs> constexpr auto first_of(Fs... fs);` | after phase 3 by default; during phase 3 if you choose (open decision 12) |
| P3 | FXT-5 | Error aggregation | `template<class... Fs> constexpr auto collect_errors(Fs... fs); // -> expected<V, std::array<E, N>>` | with FXT-4 |
| P3 | FXT-6 | Kleisli composition | `template<class F, class... Gs> constexpr auto kleisli(F f, Gs... gs);` | with FXT-4 |
| P3 | FXT-7 | First success over a range (statically typed elements; complements `nxx::first_of` over `any_solver` ranges) | `template<std::ranges::input_range R, class F, class Select = keep_last> constexpr auto first_success(R&& r, F f, Select s = {});` | with FXT-4 |
| P4 | FXT-8 | Lazy unfold view (`steps_view` could build on it); lazy fallback value | `template<class S, class F> constexpr auto iterate(S s0, F step); template<class F> constexpr auto value_or_else(F f);` | optional |
| P4 | FXT-9 | Self-contained headers (TupleAsArray, TuplePrepend, TupleReverse, TupleTransform, TypeValue, Lazy); fix `demos/LogicalOr.cpp`; document `operator\|\|`'s eager, same-type semantics; **macro include guards or angle-bracket includes**, so consumers can override single headers | various | any time |

Once FXT-4..7 exist, Numerixx's `first_of_t` can delegate to `fxt::first_of_with(nxx::detail::merge_and_charge, …)` with the same behaviour. `then` keeps its environment-threading wrapper.

### 8.3 What stays in Numerixx

- the driver (verdicts, intrinsic stop reasons, views versus estimates, counters, best-so-far), `detail::advance`, `finish`, `steps_view`;
- the stop criteria and their view-kind algebra;
- `solution`/`failure`/`fault`;
- cost accounting and `cost_of`;
- `.on()`;
- `then` (it threads f and accumulates cost, a Reader-like concern), `warm_fallback` (needs `failure::best`), `with_evaluation_budget` and `any_solver`;
- the refined *numeric* types;
- every algorithm.

---

## 9. Testing and CI

### 9.1 Layers (verified in build experiments on GCC 16, MSVC 19.51 and em++ in both EH modes; clang-cl by direct compile)

| Layer | Mechanism |
|---|---|
| Unit and property | doctest `TEST_CASE`/`TEST_CASE_TEMPLATE`/`SUBCASE` inside `TEST_SUITE`s; one executable per module, linking only that module; `doctest_discover_tests`, one CTest test per test case, labelled with its module. Without exceptions a failed `REQUIRE` reports but does not end the test case, so no test relies on `REQUIRE` to guard a dereference: results are compared whole (`CHECK(r == 3)` on an `expected`). doctest has no generators, so property tests loop over samples drawn from a `std::mt19937` with a fixed seed |
| Compile-time contracts | `static_assert` on concepts and traits; constexpr solves where portable (1-D); **regularity** (`std::copyable` for every solver, chain, `_fn` and `any_solver`; `semiregular` when parts are default-constructible); **defaults achievable for every T** (§3.5). Result sizes are printed, not asserted |
| Compile-fail | Each case is a CTest test that runs `cmake/RunCompileFail.cmake`, with a `NUMERIXX_CF_CONTROL` twin that must compile; `RESOURCE_LOCK`. The build must fail with a compiler error (a build that succeeds, or fails without one, fails the test). **On GCC and Clang (clang-cl too) the script extracts the first error (the error line and its notes, including GCC's nested error after "in 'constexpr' expansion of") and checks that the reason matches it**; for deletion reasons only on GCC 15+ and Clang 19+. On cl the case only has to fail with a compiler error. Line counts are recorded, not gated: GCC 16's nested explanations alone made 77 lines for the prototype's `then` contract **[prototyped]**. Cases: bisection given a guess; int guess; `first_of` without `.on`; open method → bracket solver in `then`; run-time `bracket{}` literal; wrong-length fixed-size guess; `x_tol` on a bracketing solver; `f_tol` on a minimiser; Newton without a derivative; mixed callback error types; an `any_solver` from a solver with a different result type. Plus the documented MSVC-only failure of the P2564 escalation probe. **[spike]** A deletion reason (`NXX_DELETE`) must appear in the compiler's own message: quoted source lines ("  123 \| code" and caret lines) are left out of the match, because a "declared here" note echoes the `NXX_DELETE("…")` line and would match even if the compiler printed no reason. The self-test `harness_quoted_reason`, which the harness must reject (`PASS_REGULAR_EXPRESSION` on its rejection), guards that rule. A rejected consteval literal is different: GCC reports the call to the non-constexpr `literal_violates_invariant("…")` and quotes the source line that makes it, which holds the reason, so there the quote counts. The P2564 probe must fail on cl with C7595 as its first error, so an unrelated build failure cannot pass for the documented gap. |
| Structural | header self-containment (every header compiled twice); layering against the DAG allow-list (Eigen only under linalg/multiroots); **consumer TU with `/W4 /WX` (MSVC) and `-Wall -Wextra -Werror` defining a global `f`** |
| Usage | `tests/usage/canonical_calls.cpp` (§6.14) on every toolchain |
| Determinism | bit-identical results on repeated runs and repeated calls; no statics or `mutable` members (review and clang-tidy) |
| EH modes | `gcc-noexcept`, `emscripten-noexcept`, `emscripten` (wasm EH), `emscripten-jsexcept` (JS EH); exception-neutrality tests under `#if __cpp_exceptions`. These legs prove that the library never throws and serve users who build without exceptions |
| Scalar matrix | `float`, `double`, `long double` everywhere; `cpp_bin_float_50` in the multiprecision leg, including linalg through `<boost/multiprecision/eigen.hpp>` |
| Oracles | Boost.Math behind `NUMERIXX_TEST_ORACLES`: its root finders (`bisect`, `toms748_solve`, `newton_raphson_iterate`), `brent_find_minima` and its quadrature (`gauss_kronrod`, `tanh_sinh`); Eigen: its own decompositions and residuals for the linalg facade, `HybridNonLinearSolver` for multiroots, `EigenSolver` on the companion matrix for poly roots (when `NUMERIXX_WITH_LINALG` is ON; test code only, so `poly` stays core-only) |
| Acceptance scenarios | generic usage scenarios (§9.3) |
| Hardening | ASan/UBSan + libc++ debug hardening; `_GLIBCXX_ASSERTIONS` |
| Examples, benchmarks | examples as smoke tests; benchmarks built on the gcc leg (not timed), reporting `fevals` |

### 9.2 Corpus and reference values

- **Published test suites**, run in addition to the targeted cases below:
  - **Alefeld–Potra–Shi** (ACM TOMS 21(3), 1995, the TOMS 748 paper): its 15 problem families of bracketing test problems, each over several parameter values (about 150 instances), for every bracketing root finder. Their mean and worst-case evaluation counts feed the default choice (D29);
  - **Moré–Garbow–Hillstrom** (ACM TOMS 7(1), 1981): the square (m = n) problems as systems of nonlinear equations (Rosenbrock, Freudenstein–Roth, Powell badly scaled, helical valley, Powell singular, extended Rosenbrock, extended Powell singular, trigonometric, Brown almost-linear, discrete boundary value, discrete integral equation, Broyden tridiagonal, Broyden banded, Chebyquad) from the standard starting point and from 10× and 100× that point, for multiroots. The least-squares and unconstrained-minimisation problems of the same collection serve v2.1 (§10.5);
  - **QUADPACK/Piessens** test integrals (the QUADPACK book's test families) and **Bailey–Borwein-style** high-precision test integrals (for example Bailey, Jeyabalan and Li, 2005), for integrate: algebraic and logarithmic endpoint singularities, interior peaks, oscillatory integrands, and semi-infinite and infinite ranges.
- **Roots:**
  - the 8 salvaged functions;
  - negative root; root at the first midpoint; root at an endpoint;
  - no sign change (x²−5 on [3, 4]); NaN inside the domain;
  - **poles** tan on [1, 2] and 1/(x−1/3) on [0, 1], which must fail `sign_change_not_root`; a **step function** (documented behaviour);
  - **extreme brackets** [−1.7e308, 1.7e308]; **roots at 1e-200, 1e-12, 1e-8 and 1e300**;
  - **a root exactly at 0**: f(x) = x + x³ on [−1, 2] (nonlinear, so an interpolating step does not land on 0 at once; the first midpoint is not 0) with the default `floored_width{}`, for every bracketing solver. It must succeed (a criterion or exact-zero stop) within the default budget. This tests the absolute floor at `scale = 1` (D32, §6.8): a purely relative width test cannot be met by an enclosure of 0;
  - **log(x) on [0, 2]** (an infinite endpoint);
  - a triple root; a flat function (x⁹);
  - **Newton cycle** x³−2x+2 from 0; **Newton from 1e-3** on x²−2 (must not false-stall);
  - **expand into a domain edge** (sqrt(x)−10 from [1, 2]);
  - a 1e-10-relative-noise f;
  - generic mixed scales: roots and brackets at magnitudes from 1e-8 to 1e8, with `scale` set and unset.
- **deriv:** x ∈ {−1e6, −2, 0, **1e-8, 1e-3**, 1, 1e3, 1e8}; second and mixed partials (x²y at (3, 2) = 6; xy at (1e6, 1e-3) = 1); an f with 1e-10 relative noise, as from an inner iterative solve.
- **optimize:** minima at ±2 and at 1000.5, a plateau, and `float` instantiations of every default.
- **integrate:** polynomials exact per rule; sin²(8πx); ∫ sin over a period (zero value); endpoint singularities; reversed limits.
- **multiroots:**
  - Rosenbrock, Powell singular, Freudenstein–Roth (local minimum), and the rest of the Moré–Garbow–Hillstrom systems above;
  - a small **badly scaled** system, with variables of magnitude 1e3 alongside variables of magnitude 1e-6 (plus Powell badly scaled);
  - a Jacobian with columns scaled 1e10 apart (must *not* be flagged singular);
  - a singular Jacobian;
  - fixed-size (`std::array`) and dynamic (`std::vector`, `VectorXd`) variants of each.
- **Reference rules** (a `long double` reference failure in the build experiments showed why):
  - references carry ≥ 21 significant digits (≥ 40 for multiprecision), generated once at high precision, by `tools/gen_reference.cpp` (with a multiprecision type) or offline with mpmath, and committed together with the command that generated them;
  - tolerances go through `tol<T>(k, ref_eps) = k·max(ε_T, ref_eps)·(1 + |x|)`;
  - Numerixx output is never its own reference.
- **Oracles** (§9.1): Boost.Math for roots, minima and quadrature; Eigen for the linalg facade and, through `HybridNonLinearSolver`, for multiroots. GSL is not an oracle and its code is not used (§1.1); values published in its documentation may serve as reference facts.

### 9.3 Properties and acceptance scenarios

- **Properties:**
  - **Criterion soundness**, run on every solver *standalone* (never behind a polishing stage): success with `stop_reason::criterion` implies the guarantee holds for the returned estimate. For width criteria, width ≤ abs + rel·min(|lo|, |hi|); for `f_tol`, |fx| ≤ tol; for `x_tol`, |x_k − x_{k−1}| ≤ tol. Regression: exp(x)−1.0001 on [0, 10] with a 1e-3 width tolerance.
  - A bracketing solver never reports convergence while its width exceeds the requested tolerance.
  - `resolution_limit` implies `hi == nextafter(lo, +inf)` (for built-in floats) for bisection. For Brent it means that the tolerance was below Brent's floor, and it implies width ≤ 4ε|x| (§6.8) **[both tested in the spike]**.
  - `first_of(s1, s2)(a) == s1(a)` when s1 succeeds, **and s2 is never invoked**; the same for `first_of` over an `any_solver` range.
  - A run-time `any_solver` chain returns exactly what the equivalent static chain returns, bit for bit, on success and on total failure **[prototyped]**.
  - `steps_view` yields the same iterates as the driver, up to the driver's stopping point.
  - Every intermediate bracket keeps a sign change and strictly shrinks.
  - **Evaluation counts equal instrumented f-call counts** (`fn::counted`) for every solver, combinator and `derivative_fn`.
  - The progress window never fires on a corpus run that converges within its budget without it.
  - Relative accuracy is scale-invariant for |x| ≥ typical (derivatives; roots with `scale` set).
  - **Error-estimate reliability:** `ridders`, `diff_with_error` and G7K15 estimates are ≥ the true error on at least 95 % of corpus cases.
  - `maximize(f) == minimize(−f)` with `fx == f(x)`.
  - Poly ring laws on integer coefficients.
  - ∫ₐᵇ = −∫ᵇₐ, and additivity.
  - Interpolants pass through their knots; monotone data stays monotone.
  - ‖Ax − b‖ ≤ c·ε·‖A‖‖x‖ for the facade's solves.
  - N-D merit decreases monotonically under line search.
  - stall ≠ budget.
- **Acceptance scenarios** (generic usage; each is a test under `tests/usage/`, except the last, which is a leg of the §9.4 consumers job):
  - coarse bisection (`width_tol{1e-3}`, budget 100), then derivative-free secant polish with a box clamp;
  - search, then solve (the search output type *is* the solver input type);
  - `first_of(expand-from-guess, subdivide)`;
  - the best estimate on budget exhaustion and on every other failure;
  - minimum and maximum of the same f, each returning `{x, fx}` with fx = f(x);
  - a final `transform` that clamps and converts into a user strong type;
  - damped N-D Newton with box constraints (`project` hook), returning the best iterate;
  - manual stepping through `steps_view` of a solve whose iterates the caller inspects;
  - first, second and mixed derivatives at mixed scales (1e-8 to 1e8) with relative, noise-aware steps;
  - a fallible callback differentiated with `relative{1e-5, 1.0}` steps and an error estimate, its error propagated unchanged;
  - a callable that mutates captured state is called sequentially and deterministically, and never copied gratuitously;
  - a scalar-only build with `NUMERIXX_WITH_FXT=OFF` and `NUMERIXX_WITH_LINALG=OFF` downloads neither FXT nor Eigen.

### 9.4 CI matrix (GitHub Actions)

| Job | Runner | Legs |
|---|---|---|
| windows | `windows-2025-vs2026` (VS 18.9, MSVC 14.51, clang-cl) | `msvc`, `clang-cl` (`CPM_SOURCE_CACHE=C:\cpm`); `/W4 /WX` consumer TU; P2564 probe must fail on `msvc` with C7595 as its first error |
| linux-clang | `ubuntu-26.04` | `clang` (libc++ 22), `clang-asan` |
| linux-gcc | container `gcc:16` | `gcc`, `gcc-noexcept` (including the FXT pipes), `gcc-multiprecision` |
| emscripten | `ubuntu-26.04` + emsdk 6.0.10 | `emscripten` (wasm EH), `emscripten-jsexcept` (JS EH), `emscripten-noexcept`, `emscripten-pthread` (wasm EH with `-pthread`); node runs the tests |
| consumers | `ubuntu-26.04` | a CPM parent with an older CPM and a parent `fxt::fxt`, and a FetchContent parent, each in both declaration orders (the parent declares FXT and Eigen first, or Numerixx first); a parent that declares its own `Boost` package; a scalar-only parent (`NUMERIXX_WITH_FXT=OFF`, `NUMERIXX_WITH_LINALG=OFF`, asserting that neither FXT nor Eigen is downloaded); install + `find_package` |
| format | `ubuntu-26.04` | clang-format 22 dry run |
| nightly | — | floor compilers (GCC 14, Clang 18 + libc++ 18), MinGW g++ (a common GCC distribution on Windows), Intel ICX (common in HPC), clang-cl `/EHs-c-` |

---

## 10. Migration roadmap

### 10.1 Starting point

Use a fresh tree on this redesign branch, which descends from master. Port algorithm *bodies* by hand, mostly from dev-reorg. It contains all of master's library code plus optimize, interpolate, `mdiff`, the fixed `DerivativeFunctor`, the `1LL` shift fix and the lowercase layout. Where dev-reorg regressed, take master's version (polynomial trimming).
- Cite provenance in each commit.
- Keep the BSL-1.0 notice for Boost-derived Brent or TOMS748 text.
- Never port or paraphrase GSL source (GPL-3.0-or-later, §1.1); implement new algorithms from the literature and cite it in each header.
- Add each module's regression tests in the same PR as its port.
- Start from the prototype (`docs/redesign/prototype/nxx/*.hpp`), which already implements the protocol, driver, criteria algebra, facade, bisection, Brent, secant, Newton, `expand`, `derivative_of`, `steps_view`, `any_solver` run-time chains, the Eigen `lu_solve` facade and damped N-D Newton. Its in-house LU is not carried over.

### 10.2 The first de-risking spike (3–4 days; right after phase 0)

Already proven:
- **Build experiments:** the CPM 0.43.2 bootstrap with parent reuse; FXT `DOWNLOAD_ONLY` + shim + override; Eigen and standalone Boost fetched `DOWNLOAD_ONLY` with own targets; test discovery (including under node; verified with Catch2 before the switch to doctest in phase 0); compile-fail with controls; header self-containment; layering; install + `find_package`; CPM and FetchContent parents; presets.
- **Linear-algebra research:** Eigen 5.0.1 under em++ 6.0.8, including `-fno-exceptions`, and with `cpp_bin_float_50` scalars.
- **Exploratory solver prototypes:** constexpr solves; ITP; N-D Newton.
- **Prototype [prototyped]:**
  - the full core on 9 configurations with bit-identical output;
  - five constexpr pipelines: the headline chain, expand → bisection → Newton, the numeric-derivative chain, Brent, `derivative_of`, plus a 2×2 damped Newton;
  - the `estimate`/`best` driver split; derivative policies; the `then` contract; the common-cause rule;
  - `steps_view` as an input range and view;
  - run-time chains (`any_solver` + `first_of` over a range), bit-identical to the static chains;
  - 16 compile-fail tests on 5 compilers;
  - FXT-1 sufficiency;
  - the Eigen facade and N-D Newton on Eigen storage;
  - compile times;
  - no clang-cl mangling problems.

Exit criteria:
1. Hosted CI green on all legs with the skeleton: MSVC and clang-cl configured **through CMake**, clang-22 + libc++, `gcc:16`, emsdk 6.0.10.
2. FXT-1 merged upstream and the pin bumped, or a patched fork commit pinned via `NUMERIXX_FXT_REF`/`CPM_FXT_SOURCE`, so `gcc-noexcept` turns green. (Met by pinning the troldal/FXT#1 commit; the pin moves to FXT's main once that is merged.)
3. `first_of`/`then` over solvers with different state types **and** fallible callbacks through `pipes.hpp` under clang-cl in CMake builds, including `.with_projection`/`.with_observer` values and the family facades.
4. CPM deduplication with a CPM parent and a FetchContent parent, in both declaration orders (the parent declares FXT and Eigen first; Numerixx first).
5. A compile-time measurement for the umbrella header (guard: 2 s on GCC), and one recorded for a linalg/multiroots TU.
6. Your written decision on §12 items 1–8 (done on 2026-09-28: every default accepted).
7. **Criterion soundness:** the §9.3 property holds on every spike solver run standalone; `x_tol` on a bracketing solver does not compile, with its reason.
8. **Canonical calls** 1, 2, 3, 9 and 10 (§6.14) compile and run on GCC, Clang + libc++, MSVC, clang-cl and em++ with **run-time** brackets, tolerances and budgets.
9. **Diagnostics:** the compile-fail cases in §9.1 put the reason string in the first error on GCC and Clang; line counts are recorded.
10. **Regularity:** combinator and function-returning results are copy-assignable (`static_assert`).
11. **Composition:** `newton.with_derivative(derivative_of(g))` with a fallible g preserves `UE` and the derivative's cause; the numeric policy works in a curried chain.

**Status (2026-09-28, branch `claude/numerixx-spike`).** Every criterion is met locally: the test-based criteria on the 11 presets that build the tests, and criterion 4 (the consumer-build scenarios) on `integration`, which builds only those. Hosted CI runs on the spike's pull request. The evidence is in Appendix D.

| # | Status | Where |
|---|---|---|
| 1 | met in phase 0 (hosted CI green on every leg) | `.github/workflows/ci.yml` |
| 2 | met: the troldal/FXT#1 commit is pinned; the pin moves to FXT's main once that PR is merged | `cmake/NumerixxDependencies.cmake` |
| 3 | met: combinators over bisection, brent, secant and newton with fallible callbacks, `.with_projection`, `.with_observer` and every facade, through `pipes.hpp`, on the `clang-cl` CMake preset (and on every other preset that builds the tests, including `-fno-exceptions`) | `tests/pipes/test_pipes.cpp` |
| 4 | met in phase 0 | `tests/integration/` |
| 5 | met: umbrella header 0.58 s on GCC (guard 2 s); a linalg TU 4.6 s on GCC (Eigen LU solves called directly: `<numerixx/linalg.hpp>` is still the phase-0 skeleton), both recorded on all four desktop compilers (Appendix D) | `tests/structural/compile_time.cmake` |
| 6 | met on 2026-09-28 | §12 |
| 7 | met: the property holds for bisection and brent (width criteria, `f_tol`, near the resolution limit, the exp(x)−1.0001 regression) and for secant and newton (`x_tol`, `step_tol`, `f_tol`), each run standalone; `x_tol` on bisection and brent (in the constructor and in `with_stop`) and `step_tol` on bisection do not compile, with the reason | `tests/roots/test_soundness.cpp`; `tests/compile_fail/{bisection_x_tol,brent_x_tol,bisection_with_stop_x_tol,bisection_step_tol}.cpp` |
| 8 | met: calls 1, 2, 3, 9 and 10 with run-time brackets, tolerances (`make()`) and budgets, on GCC, Clang + libc++, MSVC, clang-cl and em++ (4 EH/thread modes) | `tests/usage/canonical_calls.cpp` |
| 9 | met for the spike's modules: 29 cases, the reason in the first error on GCC and Clang (for deletion reasons, in the compiler's own message, not in a quoted source line), plus two harness self-tests; line counts recorded (Appendix D). §9.1's wrong-length fixed-size guess (multiroots, phase 5) and `f_tol` on a minimiser (optimize, phase 4) arrive with their modules | `tests/compile_fail/` |
| 10 | met: `static_assert` copy-assignability for solvers, criteria, refined types and results, and for solvers, curried solvers, `first_of`/`first_of_with`/`then`/`warm_fallback` chains, `derivative_of` and `fn::counted` holding capturing (hence non-assignable) lambdas, and `any_solver` | `tests/usage/test_regularity.cpp` |
| 11 | met: `UE` and the derivative's cause are preserved; the numeric policy works in a curried chain, including at compile time | `tests/usage/test_composition.cpp` |

### 10.3 Phases

Sizes are focused developer-days for one developer (rough). Every phase also has these acceptance criteria:
- all CI legs green;
- header, layering, regularity, determinism and `static_assert` suites pass;
- the phase's canonical calls compile;
- no new warnings;
- `CHANGELOG.md`/`MIGRATION.md` updated;
- tag `v2.0.0-alpha.<n>`.

| # | Phase | Scope | Deliverables | Acceptance criteria | Size |
|---|---|---|---|---|---|
| 0 | Skeleton | Tag `v1.0.0`/`v1.1.0-legacy` (§10.4); delete the old tree, vcpkg, gcem, Blaze, gbench, `.idea`; new CMake, presets, CI with every §9.4 leg, including `gcc-multiprecision` (standalone Boost.Multiprecision fetched for tests; `cpp_bin_float_50` satisfies `real` without the adapter, D16), so phases 1–8 add MP instantiations as they land | buildable empty library; CI | `cmake --workflow --preset gcc` from a clean clone; consumers job green, including the scalar-only leg with `NUMERIXX_WITH_FXT=OFF` and `NUMERIXX_WITH_LINALG=OFF`; `gcc-multiprecision` leg green | 2.5–3.5 |
| S | Spike | §10.2 | spike branch merged | exit criteria 1–11 | 3–4 |
| 1 | Core vocabulary | `config.hpp` macros; `scalar_traits`, `nxx::math` helpers (incl. `midpoint`); refined types (`abs_tolerance`, `evaluation_budget`, `is_refined_v`); `errc`/`fault`/`failure`/`solution`; `evaluate`/`evaluate_sample`; `cost_of`; unwrapping and common cause; `copyable_box`; `numerixx::pipes` | `numerixx::core`, `numerixx::pipes` | illegal-state and compile-fail suites; defaults-achievable `static_assert`s for `float`…`cpp_bin_float_50`; result sizes recorded | 2.5–3.5 |
| 2 | deriv | stencils, step specs (optimal, relative, noise, absolute), `diff` + conveniences, `diff_with_error`, `ridders`, `mixed`, `derivative_of`, `numeric` | `numerixx::deriv` (core only) | deriv corpus incl. x ∈ {1e-3, 1e-8} and mixed scales; stencil-order log-log property; error-estimate reliability; the §9.3 derivative scenarios | 3–5 |
| 3 | Driver + 1-D roots | criteria with view kinds; driver (`detail::advance`, `finish`); family facades, options, builders, input overloads; `first_of_t`/`then_t`/`warm_fallback_t`/`with_evaluation_budget`; **`any_solver` + `first_of` over a range**; **`steps_view`**; bisection, brent, illinois (Anderson–Björck), ridders, rtsafe, secant, newton; expand, scan, subdivide; `solve` facade; `inverse_of` | `numerixx::roots` | criterion soundness; roots corpus (`float`/`double`/`long double`, MP) incl. poles, extreme brackets, small roots, a root exactly at 0, cycles, and the Alefeld–Potra–Shi problems; counts match instrumented f; default bracketing solver chosen by corpus counts; the §9.3 root scenarios; `steps_view` iterates equal the driver's; run-time chain equals the static chain; canonical calls 1–3, 9–12 | 10–14 |
| 4 | optimize | golden, brent_min (intrinsic test), bracket_minimum, `maximizing`/`maximize`, `minimizer_of` | `numerixx::optimize` | optimize corpus incl. `float`; max = min(−f); the §9.3 optimisation scenarios; canonical call 4 | 2.5–4.5 |
| 5 | linalg + multiroots | Eigen 5.0.1 via CPM behind `NUMERIXX_WITH_LINALG`; the `nxx::linalg` facade (aliases, `lu_solve`, `qr_solve`, `cholesky_solve`, `vector_traits` for `std::array`/`std::vector`/Eigen); `multiroots/derivatives` (gradient, FD Jacobian with typical and projection-aware sides, true Hessian); damped Newton with `project`/`typical`/`weights` hooks, D_x scaling, full-step test, `local_minimum`, `stpmax`; Broyden (QR in state, refresh rule); dogleg (hybrd rules, scaled trust region); `system_of`; `multiroots::solve` default | `numerixx::linalg`, `numerixx::multiroots` | facade properties (‖Ax − b‖; singular → error; dimension mismatch); systems corpus incl. badly scaled systems and 1e10 column scaling, fixed and dynamic sizes; the Moré–Garbow–Hillstrom systems; Powell singular; Freudenstein–Roth → `local_minimum` handled; matches Eigen `HybridNonLinearSolver` on random systems; stall ≠ budget; the §9.3 N-D scenario (box constraints); canonical call 8; linalg TU compile time recorded | 8–15 |
| 6 | poly | polynomial, closed forms (Kahan discriminant + polish), Aberth (specified starts, running-error stop, conjugate pairing), formatting | `numerixx::poly` | ring laws; `roots ∘ from_roots`; agreement with the Eigen companion-matrix oracle when linalg is ON (test code only); regressions (operator−, trimming, divmod) | 4–6 |
| 7 | integrate + interpolate | G7K15 (QUADPACK rules, tabulated nodes), tanh-sinh (complement) + exp-sinh/sinh-sinh, Romberg (max_level 12, `min_iterations`), `gauss_legendre` (plain function), `quad`, `integral_of`, `antiderivative`; linear, splines with the in-house tridiagonal solvers (Thomas, pivoted tridiagonal, Sherman–Morrison), pchip, steffen, Floater–Hormann | `numerixx::integrate`, `numerixx::interpolate` (both core only) | corpora and properties; the QUADPACK/Piessens and Bailey–Borwein-style test integrals; sin²(8πx); zero-valued integral; last-knot regression; canonical calls 5–7 | 8–11 |
| 8 | Optional roots | toms748 (attributed Boost.Math port), itp (`itp_params`, frozen ε_ITP), halley (combined callable, bracket fallback), steffensen (deleted if the corpus shows no benefit); re-run the default-choice corpus | roots additions | agreement with the Boost.Math TOMS748 oracle; ITP within its worst-case bound on the corpus; criterion soundness for each; default-solver decision recorded | 4 |
| 9 | Multiprecision, docs, release | multiprecision adapter (`numerixx::multiprecision`: `adapters/multiprecision.hpp`, `adapters/multiprecision_linalg.hpp`) and the MP linalg tests on the existing `gcc-multiprecision` leg; docs rewritten (dev-reorg `docRoots.rst` structure, "no error handling" policy inverted); examples; benchmarks | `v2.0.0` | docs build; examples are smoke tests; MIGRATION.md complete; MP leg green including linalg | 5–7 |

**Phases 1–3 start from the spike's code** (decided on 2026-09-29). The spike built first cuts of phase 1–3 code in the library tree, more than its exit criteria needed (§10.2, Appendix D). That code is kept, and phases 1–3 continue from it. Their scope and acceptance criteria are unchanged, and a phase is done only when all of them are met. What the spike built, and what is left:

| Phase | In the spike (a first cut, tested on every preset) | Still to do | Size: planned → left |
|---|---|---|---|
| 1 Core vocabulary | the whole scope: `config.hpp` macros, scalar traits, `nxx::math`, refined types, errors and results, `evaluate`/`evaluate_sample`, `cost_of`, unwrapping and common cause, `copyable_box`, `numerixx::pipes` | the acceptance work: defaults-achievable `static_assert`s for `float`…`cpp_bin_float_50`, recorded result sizes; one error code for every non-finite input; the role types in the literal forms of `x_tol` and `width_tol`; the GCC `-ffp-contract` decision (§5.3) | 2.5–3.5 → 0.5–1 |
| 2 deriv | the stencils `central_1_2`, `central_1_4`, `central_2_2`, `central_2_4`, `forward_1_1` and `backward_1_1`; the steps `optimal`, `relative` and `absolute`; `diff`, `central`, `derivative_of`, `numeric` | `noise` steps, `diff_with_error`, `ridders`, mixed partials, `second_derivative_of`, the remaining stencils; the deriv corpus (x ∈ {1e-3, 1e-8}, mixed scales), the stencil-order property, error-estimate reliability and the §9.3 derivative scenarios. The default step needs that corpus: relative steps are scale-invariant for power laws (the derivative of 1/x has a relative error of about 6e-11 from x = 1e-8 to 1e8), but not for exp at large x (8.7e-7 at x = 300 with `central_1_2`, 2.5e-4 with `central_1_4`; measured) | 3–5 → 2–3.5 |
| 3 Driver + 1-D roots | criteria with view kinds; the driver; facades, options, builders and input overloads; `first_of`, `then`, `warm_fallback`; `any_solver` and `first_of` over a range; `steps_view`; bisection, brent, secant, newton, expand; `solve(f, bracket)`. Tested: criterion soundness, poles, extreme brackets, a root at 0, counts equal to instrumented calls, `steps_view` equal to the driver, run-time chains equal to static chains, canonical calls 1–3, 9 and 10 | `with_evaluation_budget`; illinois, ridders, rtsafe; scan, subdivide; `expand` from a guess, its NaN backtrack and crossing 0; `solve(f, x0)`, `solve(f, df, x0)`; `inverse_of`; the open-method safeguards (without the progress window, Newton on x³ − 2x + 2 from 0 runs its whole budget, 30 iterations and 61 evaluations, before failing); the representation-space bisection midpoint (without it, a root at 1e-200 with a relative tolerance exhausts bisection's 200-step budget); `floored_width`'s `scale`, `step_tol`'s `typical`, secant's `x1`; `custom` criteria and `fdf`; the roots corpus with the Alefeld–Potra–Shi problems, and the default-solver choice; canonical calls 11 and 12 | 10–14 → 6–9 |

This moves work into the spike rather than saving it: the spike did much more than its 3–4 days. Left for phases 1–9: 40–61 developer-days, or 36–57 without phase 8 (planned: 47–70, or 43–66).

**Totals, and how they follow from the baseline.** The baseline is 57.5–82.5 developer-days: the estimate for an earlier scope that also required zero heap allocation and AD scalars, both now out of scope (§3.5, §3.6), taken without its project-specific migration work. Each change below is a single number applied to both ends of the range.

| Change | Phase | Days |
|---|---|---|
| Spike: drop the AD-scalar exit criterion and the two-TU link test | S | −1 |
| Core: drop the AD primal type, the AD-aware routing of maths calls and the error-size gate (sizes are recorded, not gated, D7) | 1 | −0.5 |
| Roots: drop an externally driven single-step ITP stepper and its budget-in/budget-out tests (manual stepping is `steps_view`, D19) | 3 | −0.5 |
| Roots: drop the AD scalar tests and the extra Newton step for AD types | 3 | −0.5 |
| Roots: drop `roots::sensitivity` | 4 | −0.5 |
| Linalg: drop the in-house kernels and storage (vec/mat, views, equilibrated LU + rcond, Cholesky, QR) | 5 | −3.5 |
| Multiroots: drop `newton_inplace` and the span + workspace API | 5 | −1 |
| Drop the zero-heap tests (allocation counting, a run with a small default stack, and an `nm` check that no allocator symbol is linked) | 5 | −0.5 |
| Drop arclength continuation | 5 | −2 |
| CI: drop the zero-heap checks' small-stack preset and `nm` symbol-check job | 0 | −0.5 |
| CI: add the `NUMERIXX_WITH_LINALG=OFF` consumer leg | 0 | +0.5 |
| Add `any_solver` and `first_of` over a range (prototyped; productise, add the policy overload, test, document) | 3 | +1 |
| Add `steps_view` as a phase-3 deliverable (prototyped; productise, test, document) | 3 | +1 |
| Remove the optional `steps_view` slot from the release phase | 9 | −0.5 |
| Widen the Eigen facade to QR and Cholesky and three storage adapters | 5 | +1 |
| Re-run the default-choice corpus with TOMS748 and ITP | 8 | +0.5 |
| Multiprecision adapter: Eigen interop and MP linalg tests | 9 | +0.5 |
| Add the published test suites (§9.2): Alefeld–Potra–Shi (phase 3), Moré–Garbow–Hillstrom (phase 5), QUADPACK/Piessens and Bailey–Borwein-style integrals (phase 7), 0.5 each | 3, 5, 7 | +1.5 |
| **Net** | | **−5** |

Work also moves between phases (zero net): toms748 (1.5) and itp (1) to phase 8; illinois, ridders and scan (1 together) into phase 3; halley and steffensen (1 together) to phase 8; the tridiagonal solvers (0.5) to phase 7. The `gcc-multiprecision` CI leg also moves from phase 9 to phase 0, so that the multiprecision checks of phases 1 and 3 have a leg to run on; it is one preset and one CI leg next to the others phase 0 already sets up, well under the rounding of either range, so no number changes. The same holds for the `emscripten-pthread` and `emscripten-jsexcept` legs (§4.5). Per phase:

| New phase | Baseline phase (size) | Changes | New size |
|---|---|---|---|
| 0 Skeleton | skeleton (2.5–3.5) | −0.5 + 0.5 | 2.5–3.5 |
| S Spike | spike (4–5) | −1 | 3–4 |
| 1 Core vocabulary | core vocabulary (3–4) | −0.5 | 2.5–3.5 |
| 2 deriv | deriv (3–5) | 0 | 3–5 |
| 3 Driver + 1-D roots | driver + roots (10–14) | −1 removed; +2 `any_solver`, `steps_view`; −2.5 toms748, itp out; +1 illinois, ridders, scan in; +0.5 Alefeld–Potra–Shi suite | 10–14 |
| 4 optimize | optimize + remaining roots (5–7) | −0.5 `sensitivity`; −1 illinois, ridders, scan out; −1 halley, steffensen out | 2.5–4.5 |
| 5 linalg + multiroots | linalg + N-D Newton (8–12) plus N-D extensions (6–9) = 14–21 | −7 removed; +1 facade; −0.5 tridiagonal out; +0.5 Moré–Garbow–Hillstrom suite | 8–15 |
| 6 poly | poly (4–6) | 0 | 4–6 |
| 7 integrate + interpolate | integrate + interpolate (7–10) | +0.5 tridiagonal in; +0.5 QUADPACK/Bailey–Borwein integrals | 8–11 |
| 8 Optional roots | — | +2.5 toms748, itp in; +1 halley, steffensen in; +0.5 corpus re-run | 4 |
| 9 Multiprecision, docs, release | adapters, docs, release (5–7) | −0.5 + 0.5 | 5–7 |

- **Total:** 57.5–82.5 − 5 = **52.5–77.5 developer-days**. Check by phase: 2.5 + 3 + 2.5 + 3 + 10 + 2.5 + 8 + 4 + 8 + 4 + 5 = 52.5 and 3.5 + 4 + 3.5 + 5 + 14 + 4.5 + 15 + 6 + 11 + 4 + 7 = 77.5. Check of the baseline column: 2.5 + 4 + 3 + 3 + 10 + 5 + 14 + 4 + 7 + 5 = 57.5 and 3.5 + 5 + 4 + 5 + 14 + 7 + 21 + 6 + 10 + 7 = 82.5. Without the optional phase 8: 48.5–73.5.
- **The published test suites** add 1.5 days (0.5 each in phases 3, 5 and 7); without them the total would be 51–76.
- **Order:** phases 0 → S → 1 → 2 → 3 are sequential. deriv (phase 2) comes before the driver on general grounds: it is small, needs only core, exercises the evaluation vocabulary (faults, callback errors, `cost_of`) before the driver lands, and is used by phase 3 (the numeric-derivative Newton) and phase 5 (FD Jacobians). After phase 3, phases 4, 5, 6 and 7 depend only on core and the driver and follow breadth of use; with a second developer, 6 and 7 run in parallel with 4 and 5. Phase 8 is optional and may follow `v2.0.0` without changing the API of the other phases. The families planned after v2.0 are in §10.5.

### 10.4 Migration from Numerixx 1.x

1. **Tags.** Phase 0 tags `v1.0.0` (master, 5de1e07) and `v1.1.0-legacy` (the dev-reorg tip, 8528e94), so existing users can pin the old API (for example `GIT_TAG v1.0.0` in CPM or FetchContent) and migrate on their own schedule.
2. **`MIGRATION.md`** maps the old API to the new one and grows with each phase:

   | Numerixx 1.x | Numerixx 2 |
   |---|---|
   | `fsolve<Bisection>` | `roots::bisection{}(f, {lo, hi})`, or the facade `roots::solve(f, {lo, hi})` |
   | `fdfsolve<Newton>(f, df, x0)` | `roots::newton{}.with_derivative(df)(f, x0)`; safeguarded: `roots::solve(f, df, x0)` |
   | `fdfsolve<Secant>` | the derivative-free `roots::secant{}(f, x0)`: no derivative argument |
   | `.result()`, `.result<T>()` | the returned `std::expected` itself: check it (`if (r)`, then `r->x`) or use `value_or`; a failure never arrives as a plain value |
   | `.result(fn)` | `transform(fn)`, a member of `std::expected` or an FXT pipe |
   | `search<...>` | `roots::expand`, `roots::scan` or `roots::subdivide`; the result is the input of every bracketing solver |
   | `fminimize` / `fmaximize` | `optimize::minimize` / `optimize::maximize`, returning `extremum{x, fx}` |
   | `diff<ALGO>(f, x)` | `deriv::diff(f, x, stencil)`, for example `deriv::diff(f, x, deriv::central_1_4)` |
   | `derivativeOf(f)`, `integralOf(f)` | `deriv::derivative_of(f)`, `integrate::integral_of(f)`, both of which now keep f |
   | `multisolve<MultiNewton>` | `multiroots::newton{}(F, x0)`, or the facade `multiroots::solve(F, x0)` |
   | `polysolve` | `poly::roots`, with the closed forms `linear`, `quadratic` and `cubic` for degrees 1–3 |
   | per-module includes such as `<Deriv.hpp>` | one include root: `<numerixx/deriv.hpp>` and so on (§5.1) |

   Every module is `nxx::<module>` in `<numerixx/<module>.hpp>` with the CMake target `numerixx::<module>`, because namespace = header = target (D23).
3. **No compatibility shim.** A shim would have to reproduce semantics the redesign removes on purpose: unchecked results, a silent last iterate, `ResultProxy`, a bare-double argmin, success at maxiter.
4. **Results change where v1 was wrong**, and `MIGRATION.md` and the release notes say so:
   - the default second and mixed derivatives, which were wrong by O(1) because of their default steps (§7.1);
   - success reported at the iteration limit, now `budget_exhausted` carrying the best estimate (§6.7);
   - NaN, infinite or arbitrary values returned as roots (v1 Newton returned `inf` or −0.87 on x²+1), now in-band errors such as `zero_derivative`, `non_finite_value`, `stalled` or `diverged` (§7.2);
   - brackets without a sign change that "converged" (x²−5 on [3, 4] gave 3.99999997; dev-reorg's bisection on [5, 6] gave 5.998), now `no_sign_change`;
   - non-convergence returned as a plain value (dev-reorg's secant on x²+1 gave 3.89), now an in-band error carrying the best estimate (§7.2).

### 10.5 Beyond v2.0

**Rough estimates, not planned in detail.** Sizes use the §10.3 basis (focused developer-days for one developer), but they come from analogy with the v2.0 phases, not from a task breakdown: treat them as ±50 %. Each family adds a module with its own namespace, header and target (D23), placed downstream of the v2.0 modules in the DAG (§5.2), and reuses the core unchanged (§1.1). The Eigen-backed modules exist only when `NUMERIXX_WITH_LINALG` is ON; the umbrella header leaves them out, as it does linalg and multiroots.

| Release | Scope | Module → dependencies | Reused core pieces | New vocabulary | Rough size |
|---|---|---|---|---|---|
| **v2.1** | **Multidimensional minimisation**: Nelder–Mead (adaptive parameters); BFGS and L-BFGS with a line search that satisfies the strong Wolfe conditions (Moré–Thuente); nonlinear conjugate gradient (Polak–Ribière+ or Hager–Zhang) | `multimin` → optimize, multiroots (hence linalg, deriv, Eigen) | driver, criteria algebra, `solution`/`failure`, combinators; `maximizing` from `optimize`; `gradient_of` from `multiroots/derivatives.hpp`; linalg aliases and `vector_traits` (fixed and dynamic sizes); the `project` and `typical` hooks | `extremum_nd<V>{x, fx, gradient_norm}`; an N-D extremum view kind for the criteria; `g_tol`, today an option of multiroots' `local_minimum` test (§7.5), promoted to a stop criterion for minimisers, which adds a row to the §6.8 table; a line-search policy type; a validated initial `simplex` (n + 1 affinely independent points) | 8–12 |
| **v2.1** | **Nonlinear least squares**: Levenberg–Marquardt (Moré's trust-region form with scaling; geodesic acceleration optional) on the `qr_solve` facade | `fit` → multiroots (hence linalg, deriv, Eigen) | the systems machinery of `multiroots` (a callable V → W with m ≥ n; the FD Jacobian with typical scaling and projection-aware sides; the hooks; the weighted merit), driver, criteria, the `local_minimum` test, combinators | `least_squares_estimate<V>{x, residual, cost, rank, covariance}`, with the covariance σ²(JᵀJ)⁻¹ taken from the QR factors; relative-cost and gradient criteria | 6–9 |
| **v2.2** | **ODE initial-value problems**: Dormand–Prince RK45 (FSAL, embedded error estimate, PI step-size control, 4th-order dense output) first; then a stiff method (Rosenbrock or variable-order BDF) with stiffness detection | `ode` → multiroots (hence linalg, deriv, Eigen): vector states through `vector_traits`; `jacobian_of` and `lu_solve` for the stiff method | driver (a step advances t; the budget bounds the number of steps), criteria (reaching t_end is an intrinsic stop), `solution`/`failure` (the best estimate is the last accepted point), `steps_view` (the trajectory), function-returning APIs (dense output), `first_of` (non-stiff, then stiff); for the stiff method, the linalg facade and `jacobian_of` | `ode_problem` (f(t, y) and y₀); `time_interval` (finite, either orientation, like `interval<T>`); `step_tolerance{abs, rel}` (error per step, componentwise); `ode_point{t, y}`; `dense_solution` (a callable t → `expected<y, fault>`) | RK45 6–9; stiff 8–12 |
| later | **Chebyshev approximation** on an interval: fit at Chebyshev nodes, Clenshaw evaluation, derivative and integral as series | `chebyshev` → poly (no Eigen) | `interval<T>`; function-returning APIs (the series is a callable with `.derivative()`, so Newton is exact on it); conversion to `poly`; `roots` on the approximant (user-level composition, no module edge) | `chebyshev_series<T>`, with its degree and an error estimate from the tail coefficients | 3–5 |
| later | **Series acceleration**: Richardson, Wynn epsilon, Levin u | `series` → core (no Eigen) | driver and criteria over a sequence (a step consumes one term); `steps_view`; the Richardson extrapolation of `deriv::ridders`, moved into a `core/detail` header so that neither module depends on the other | a term source (a callable n → `expected<T, E>`); `accelerated_sum<T>{value, error, terms}` | 3–4 |
| later | **Linear least-squares fitting**: weighted fits over a general basis, and polynomial fits | in `fit`, which gains an edge to poly | `qr_solve`; `poly::polynomial`; `least_squares_estimate` from v2.1 | a basis: a callable x → one row of the design matrix | 2–3 |

- **Totals (rough):** v2.1 = (8–12) + (6–9) = 14–21 days; v2.2 = (6–9) + (8–12) = 14–21; later = (3–5) + (3–4) + (2–3) = 8–12; in all 36–54 developer-days beyond v2.0.
- **Tests.** v2.1 adds the Moré–Garbow–Hillstrom least-squares and minimisation problems (§9.2); v2.2 adds standard non-stiff and stiff test problems (for example the DETEST set, Robertson and Van der Pol) with high-precision references.

---

## 11. Risks and mitigations

| Risk | Likelihood / impact | Mitigation |
|---|---|---|
| Existing 1.x users cannot upgrade in place (the API breaks) | certain / medium | tags `v1.0.0` and `v1.1.0-legacy` to pin the old API; `MIGRATION.md` maps every old call (§10.4) |
| FXT-1 delayed | medium / low | only `numerixx::pipes` includes FXT; pin `NUMERIXX_FXT_REF` or `-DCPM_FXT_SOURCE` to a **patched fork commit**. Shadowing the two headers through the include path does not work, because FXT includes them by quoted relative path and guards only with `#pragma once` **[prototyped]**. FXT-9 adds guards. |
| clang-cl mangling on variadic constrained combinators | low / medium | none found with these shapes **[prototyped]**; unconstrained variadics with `static_assert`s; bool variable templates; clang-cl leg from day one |
| MSVC lacks P2564 (consteval escalation) | certain / low | rule: forward refined types, never raw scalars; P2564 probe on the MSVC leg; `is_refined_v` for generic detection |
| MSVC shows no deletion reasons | certain / low | deleted declarations sit on a line that holds the reason; `static_assert`-based contracts in combinators give messages everywhere |
| The floor compilers (GCC 14, Clang 18) miss something the design uses | medium / low | nightly floor job; `NXX_DELETE` falls back to plain `= delete` below GCC 15 and Clang 19; raise the floor if a needed feature is missing |
| Numerical defaults unsuitable at small scales | medium / medium | D32 `typical`; relative derivative steps; `floored_width{scale}`; small-scale corpus entries; scale-invariance property |
| False stall or divergence from the progress window | medium / medium | step-shrink guard; corpus property "never fires on a run that converges"; configurable window and an off switch |
| Pole check rejects a legitimate root | low / medium | applies only to criterion or `resolution_limit` stops and compares with the initial samples; corpus covers steep simple roots; documented |
| Learning curve (`.on`, `then`, views versus estimates) | medium / medium | facade (§6.13) and canonical calls (§6.14) are the front door; reasoned diagnostics |
| Compile time and symbol length grow | medium / low | per-module headers; no `fxt.hpp`; umbrella excludes linalg, multiroots, adapters and `any_solver.hpp`; CI compile-time guard; `extern template` for `double` in the facade if needed; named combinator types shorten symbols |
| Eigen compile cost spreads beyond linalg/multiroots TUs | medium / low | only linalg, multiroots and the MP-linalg adapter include Eigen (layering test), and the umbrella header includes none of them; `NUMERIXX_WITH_LINALG=OFF` for scalar-only consumers; per-TU time recorded per release |
| The facade misreports singularity (`PartialPivLU` never flags it; `rcond` is not scale-invariant) | medium / medium | rcond against n·ε plus a finiteness check of the solution; corpus with 1e10 column scaling; equilibrate before the rcond test if the corpus demands it; Eigen `HybridNonLinearSolver` oracle; MP leg |
| Eigen expression templates dangle under `auto` | medium / medium | concrete return types in the facade; no `auto` on Eigen expressions in library code (§5.3); documented for users |
| Eigen `operator<<` with `cpp_bin_float` fails on Boost 1.92 | certain / low | tests never stream MP matrices; revisit on the next Boost release |
| Dynamic N-D states allocate on every step | certain / low | negligible next to evaluations of a nontrivial F; fixed-size path when N is known at compile time |
| `any_solver` overhead or misuse (allocation on wrap or copy, per-call conversion of a non-`F` callable, indirect calls, not constexpr, copyable solvers only) | low / low | opt-in header outside the umbrella; static chains stay the default; never-empty invariant (no default constructor, no moves); costs and pitfalls documented (§6.10); equivalence property against static chains, bit-identical on 9 configurations **[prototyped]** |
| Result types grow (failure ≈ 88 B for roots) | certain / low | still cheap to copy and allocation-free for scalar families; sizes recorded |
| Behaviour change for 1.x users (default second and mixed derivatives, success at maxiter, NaN roots) | certain / low | intended, because v1 was wrong; listed in `MIGRATION.md` and the release notes (§10.4) |
| A noisy f contradicts cached endpoint signs | low / medium | documented; `noise` steps and relative tolerances for noisy f |
| Windows MAX_PATH (silent file loss; emsdk cache) | medium / high on developer machines | `CPM_SOURCE_CACHE=C:\cpm`; short `EM_CACHE`; CMP0168 default; CMake ≥ 3.30; SHORT intermediate directories |
| GitLab regenerates the Eigen archive and the hash breaks | low / low | `GIT_TAG 5.0.1` fallback; `CPM_Eigen_SOURCE` override |
| A parent's FXT is configured for `tl::expected` | low / medium | configure-time warning; public signatures use `std::expected` |
| Scope creep (8 modules; pressure toward GSL's breadth) | medium / medium | the §1.1 scope table, with an explicit out-of-scope list and unscheduled candidates; dependency-first order; poly, integrate and interpolate may slip past `v2.0.0-beta`; phase 8 may follow `v2.0.0`; the v2.1 and v2.2 families wait for `v2.0.0` (§10.5) |
| Provenance and licences (Brent, TOMS748 from Boost.Math; GSL is GPL-3.0-or-later) | low / high | BSL-1.0 attribution headers; no GSL source ported or paraphrased; algorithms implemented from the literature, with references in each header (§1.1, §10.1) |
| Reference-value precision makes wide-type tests flaky | medium / low | §9.2 reference rules; the generator tool |

---

## 12. Decisions

**Decided on 2026-09-28: the author accepted every default below.** Until then these were the open decisions; items 1–8 gated the spike (§10.2, criterion 6).

1. **Boost.** None in the library; standalone Boost.Config + Multiprecision (+ Math) via CPM only for the optional multiprecision adapter and test oracles. **Default: accept.**
2. **Linear algebra.** Eigen 5.0.1 as the backend behind the `nxx::linalg` facade (+2.5–6.4 s per linalg/multiroots TU in use; no `constexpr` N-D solves), or in-house kernels (constexpr and allocation-free, but more code to write and validate)? **Default: Eigen behind the facade.**
3. **Criterion mismatches.** `x_tol` on a bracketing solver (and `f_tol` on a minimiser) is a compile error with a reason, rather than being silently reinterpreted as a width test. **Default: compile error.**
4. **Reach of "illegal states unrepresentable".** Compile-time rejection for literals; run-time inputs (braced, pair, `make()` results) accepted by every solver and validated in-band; solver values never hold invalid configuration. **Default: accept.**
5. **API break and naming.** snake_case; namespace = header = CMake target for every module (`nxx::<module>`, `<numerixx/<module>.hpp>`, `numerixx::<module>`); builders instead of positional configuration (`secant{}.with_budget(5)`); `integrate::quad`. **Default: accept.**
6. **Scale defaults (D32).** Derivative steps are relative to |x| (1 at x = 0). Stopping tolerances keep an absolute floor at scale 1, so roots at 0 terminate. Both can be overridden through `typical`/`scale`. **Default: accept.**
7. **v2.0 scope.** The eight current modules (deriv, roots, optimize, poly, linalg + multiroots, integrate, interpolate, on top of core); multidimensional minimisation (module `multimin`) and nonlinear least squares (module `fit`) in v2.1, ODE initial-value solvers (module `ode`) in v2.2, each a new module downstream of the v2.0 ones (§5.2, §10.5); the "candidate" areas of §1.1 unscheduled; the out-of-scope areas of §1.1 stay out. **Default: accept.**
8. **Features within the v2.0 modules.** `steps_view` in the roots phase as the only public manual-stepping API (`advance` internal); defer `first_of_all`; drop `retry`, arclength continuation and `roots::sensitivity`; TOMS748, ITP, Halley and Steffensen in the optional phase 8; algorithm ids are an open enum. **Default: accept.**
9. **`first_of` fall-through.** Continue after input and numerical errors; stop only when the user marks a callback error fatal (`nxx::is_fatal`); `first_of_with(policy, …)` for anything else. **Default: accept.**
10. **Facade defaults.** `solve(f, x0)` = expand + Brent; `solve(f, df, x0)` = expand + rtsafe; the bracketing default is chosen by corpus evaluation counts in phase 3 (re-checked in phase 8); `multiroots::solve` = dogleg + Broyden once it passes the phase-5 corpus. **Default: accept.**
11. **FXT dependency shape.** `numerixx::core` has no FXT; `numerixx::pipes` links it; `NUMERIXX_WITH_FXT` is ON by default (users who do not want the pipes set it OFF). **Default: accept.**
12. **When combinators move to FXT.** Keep `first_of_t`/`then_t` in Numerixx until phase 3 has settled their semantics, then upstream the generic parts (FXT-4..7), or land FXT-4 during phase 3 so `nxx::first_of` is FXT code from day one? **Default: after phase 3.**
13. **Multiprecision.** Keep an optional adapter and CI leg, or drop the claim? **Default: keep.**
14. **Complex scope.** Poly only, or also complex open methods? **Default: poly only.**
15. **Compiler floor.** GCC 14, Clang 18 + libc++, MSVC with `/std:c++latest` (19.51 tested), clang-cl (22 tested), em++ ≥ 6.0.8; confirmed by the nightly floor job; MinGW g++ and Intel ICX nightly only. Note that this is one step above FXT's own stated floor (its README says GCC 13+, Clang 17+, MSVC 19.30+): Numerixx's use of deducing `this` sets GCC 14 and Clang 18 (D2). Lowering the floor to FXT's would mean replacing deducing `this` in the facades with CRTP. **Default: accept.**
16. **Low-value algorithms.** Keep Ridders and Illinois (Anderson–Björck by default); plain regula falsi only as an Illinois variant; Steffensen at low priority in phase 8, deleted if the corpus shows no benefit; Gauss–Legendre only as a plain function. **Default: accept.**
17. **`NUMERIXX_WITH_LINALG` default.** ON, so the target `numerixx::numerixx` is complete (the umbrella header still leaves linalg and multiroots out, §5.2); users who need only the scalar modules set it OFF and never download Eigen. **Default: accept.**
18. **Branch hygiene.** After tagging `v1.0.0` and `v1.1.0-legacy` (and any other branch tip worth keeping), delete `dev`, `dev-terminator` and `dev-reorg`. **Default: keep the tags, delete the branches.**
19. **Run-time chains.** Offer `nxx::any_solver` and `first_of` over a run-time range, in an opt-in header built on `std::function`, next to the default static chains (prototyped on all 9 configurations). An empty run-time chain fails in-band with `invalid_input` when called, rather than being made unrepresentable. **Default: yes, opt-in header.**

---

## Appendix A: Alternatives considered

| Question | Alternatives | Decision | Rationale |
|---|---|---|---|
| Call order | `(input, f)` vs `(f, input)` | `(f, input)` + `.on(input)` | matches the facade and the old APIs |
| Bracket type | `Bracketed<F,T>` vs `sign_bracket<T>` vs plain `Bracket<T>` | `sign_bracket<T>` (samples may be ±inf) | binding to the closure type breaks chains |
| Newton's derivative | implicit FD vs explicit | explicit: callable, policy or structural | hidden cost; the policy keeps curried chains possible |
| Where FXT lives | inside the facade vs at the edges | edges (`numerixx::pipes`) | the core and its `-fno-exceptions` legs do not wait for FXT-1; a user who does not want the pipes need not fetch FXT |
| Linear algebra | in-house kernels vs Eigen in the core vs Eigen as an optional adapter | Eigen 5.0.1 as the backend, behind a facade | your request; BLAS/LAPACK-free and verified under Emscripten; mature fixed and dynamic sizes; compile cost confined to linalg TUs. In-house kernels would give constexpr, allocation-free solves, which the design does not require (allocation is allowed, §3.6) |
| Run-time-sized N-D systems | bounded vectors vs span + workspace vs Eigen dynamic | Eigen dynamic vectors in value states | allocation is allowed; one storage model; no in-place exception to immutability |
| Chains | static templates only vs also type-erased | static by default; opt-in `any_solver` | configuration-driven chains are a real need; `std::move_only_function` is missing on libc++ 22, so `std::function`; prototyped bit-identical to the static chains |
| Manual stepping | public `init`/`step` + `advance` loop vs a lazy range | `steps_view` | one public API; composes with `std::views`; a hand loop does not type-check |
| Scalar genericity | AD primal-type machinery vs plain open trait | open `scalar_traits` + ADL maths | AD scalars are not required; `float`, `long double` and multiprecision still work |
| Criteria composition | operators vs named `any_of` | operators via hidden friends; `any_of`/`all_of` as spellings | verified on all compilers |
| `x_tol` on brackets | reinterpret as a width test vs reject | reject with a reason | a silent change of meaning is what caused the bug |
| Budget | inside criteria vs separate | separate and mandatory | no criterion can make the loop unbounded |
| Error payload | `source_location` vs algorithm id | open `algo` id, no `source_location` | a small, stable id names the failing stage of a chain; `source_location` points into the library |
| Error size | ≤ 16-byte code only vs best estimate | best estimate + enclosure in the failure; sizes recorded, not gated | callers act on the best estimate (R-A3, `warm_fallback`); libc++ `expected` is not trivially copyable anyway |
| Windowed `no_progress` | a criterion vs solver state | solver state | criteria are stateless |
| Divergence | fail on the first large step vs cap it | cap, and diverge via the window | avoids killing slow but convergent Newton runs |
| Run-time budgets | `with_budget(expected<…>)` vs validated only | validated only | a solver value must not hold an invalid configuration |
| Diagnostics | a ≤ 10-line budget vs reason in the first error | reason in the first error; lines recorded | GCC 16's nested explanations alone exceed 10 lines |
| N-D hooks | acceptance, merit, weights, projection, scaling vs a minimal set | `project`, `typical`, optional `weights` | the general needs (box and domain constraints, scaling, weighting); custom acceptance goes through `custom{λ}` |
| Quadrature rule | G7K15 vs G10K21 | G7K15 default, G10K21 optional | QUADPACK standard |
| Namespace | `nxx::optim` vs `nxx::optimize` | `nxx::optimize` | namespace = header = target |
| Phase order | roots first vs deriv first | dependencies first, then by breadth of use: vocabulary → deriv → driver + 1-D roots → optimize → linalg/multiroots → poly → integrate/interpolate → optional roots | deriv is small, needs only core, exercises the evaluation vocabulary before the driver lands, and is used by phase 3 (numeric-derivative Newton) and phase 5 (FD Jacobians); after that, the most widely used modules come first |

---

## Appendix B: Requirements

Derived from the author's brief (functional style, solver chaining, functions as first-class values, `std::expected` instead of exceptions, immutable objects, illegal states unrepresentable, Eigen for linear algebra, Boost via CPM only where something needs it, no vcpkg, Emscripten) and from the general-purpose goal (§1.1). Each line says where the design meets it.

**Not required**, by the author's decision; the design states the general reason where each topic appears: zero heap allocation (§3.6), AD scalars (§3.5) and `constexpr` N-D solves (§7.5).

**Build**

| ID | Requirement | Where met |
|---|---|---|
| R-B1 | Header-only C++23; one include root `<numerixx/…>`; one INTERFACE target per module, `numerixx::<module>` (namespace = header = target) | D2, D23, D24, §5 |
| R-B2 | Dependencies through CPM only, pinned to a version or commit plus a SHA256; no vcpkg | D25, D26, §4.2 |
| R-B3 | The library itself needs no Boost; standalone Boost only for the optional multiprecision adapter and the test oracles | D27, §4.2 |
| R-B4 | A good subproject: reuse a parent's CPM, FXT and Eigen targets; never claim a common package name (`Boost`); no global flags; overridable as CPM package `Numerixx` (`CPM_Numerixx_SOURCE`); tested with CPM and FetchContent parents in both declaration orders | D26, §4.1, §4.4, §9.4 |
| R-B5 | Optional parts stay optional: FXT only for `numerixx::pipes`, Eigen only for linalg and multiroots; a scalar-only build downloads neither | D18, D22, §4.3 |
| R-B6 | Consumers see no Numerixx warnings (SYSTEM includes) and stay clean under `/W4 /WX` and `-Wall -Wextra -Werror` | D24, §5.3, §9.1 |
| R-B7 | Linear algebra through Eigen behind the `nxx::linalg` facade; no BLAS or LAPACK | D22, §7.5 |
| R-B8 | MIT licence: no GPL-derived code (GSL); code derived from Boost keeps its BSL-1.0 notice | §1.1, §10.1 |

**Error model**

| ID | Requirement | Where met |
|---|---|---|
| R-E1 | The library never throws; `std::expected` for failures with a reason; `std::optional` only where absence is the answer | D9, D10, §3.4 |
| R-E2 | Exception neutrality: callbacks may throw (third-party code); conditional `noexcept`; `-fno-exceptions` supported and CI-tested | D10, §3.4, §9.1 |
| R-E3 | Fallible callbacks `x → expected<T, E>` are first-class, and the user's `E` reaches the caller unchanged | §3.4, §6.4 |
| R-E4 | No silent failure: no sign change, non-convergence, NaN, poles and the iteration limit are in-band errors, never successes | §3.4, §6.7, §7.2 |
| R-E5 | Errors are cheap to copy and hold no heap memory of their own | D7, §6.3 |

**API semantics**

| ID | Requirement | Where met |
|---|---|---|
| R-A1 | Functional style: immutable solver values, pure `init`/`step`, one driver | D3, D5, §3.1, §3.2 |
| R-A2 | Solver chaining (fallback, staging, warm restart, a shared evaluation budget), static or assembled at run time | D17, D35, §6.10 |
| R-A3 | The best estimate and the counters on every exit, success or failure | D5, D7, §6.7 |
| R-A4 | Functions as first-class values: derivative, inverse, integral, antiderivative, interpolant and minimiser as callables | D21, §6.12 |
| R-A5 | Honest evaluation accounting and evaluation budgets, because f may be expensive (a simulation, an inner iterative solve, a table look-up) | D33, §6.4, §6.10 |
| R-A6 | Run-time inputs are first-class (braced `{lo, hi}`, `std::pair`, `make()` results); the output of a bracket search is the input of a bracketing solver | §6.5, §7.2 |
| R-A7 | Illegal states unrepresentable where C++ allows: invalid literals and solver/input/criterion mismatches do not compile; run-time values are validated once | D6, D11, §3.3 |
| R-A8 | Box and domain constraints on the iterates (per-iterate projection), and a final clamp through `transform` | D20, §6.9, §7.5 |
| R-A9 | Minimum and maximum return a named `{x, fx}`, so the argument of an extremum is never confused with its value | §7.3 |
| R-A10 | Manual stepping and tracing through one lazy range | D19, §6.9 |
| R-A11 | Stateful user callables are invoked sequentially and deterministically, without gratuitous copies | §3.2, §9.3 |
| R-A12 | Derivative steps suited to noisy f (computed with limited precision or by inner iterative solvers), per-stencil optimal steps, one-sided stencils for the edges of a domain, and error estimates | D32, §7.1 |
| R-A13 | Deterministic results: bit-identical on repeated runs and calls; no statics or mutable caches | §3.2, §6.1, §9.1 |

**Genericity**

| ID | Requirement | Where met |
|---|---|---|
| R-S1 | Generic real scalars (`float`, `double`, `long double`, multiprecision) through an open `scalar_traits` | D14, §3.5, §6.1 |
| R-S2 | Every default tolerance is an expression in `T` and achievable for every supported `T` | D14, §3.5 |
| R-S3 | Scale awareness: relative steps and tolerances with an explicit absolute floor; mixed scales (1e-8 to 1e8) and badly scaled systems work | D32, §9.2 |
| R-S4 | N-D systems of fixed and dynamic size (`std::array`, `std::vector`, Eigen vectors) | §6.5, §7.5 |
| R-S5 | Complex arithmetic where it is routine: polynomial roots | D15, §7.4 |

**Portability**

| ID | Requirement | Where met |
|---|---|---|
| R-P1 | GCC ≥ 14, Clang ≥ 18 with libc++, MSVC (`/std:c++latest`), clang-cl, em++ ≥ 6.0.8; MinGW g++ and Intel ICX nightly | D2, §9.4 |
| R-P2 | Emscripten in every exception mode (`-fexceptions`, `-fwasm-exceptions`, `-fno-exceptions`) and with `-pthread`, also under a consumer's global flags | §4.5, §9.4 |
| R-P3 | The exception model belongs to the consumer: no EH flags in the usage requirements | §4.1 |
| R-P4 | Public signatures use `std::expected`, never a configuration-dependent alias | §8.1 |
| R-P5 | Builds survive Windows path-length limits (MAX_PATH) | §4.5 |

---

## Appendix C: Prototype evidence

**Prototype** (`docs/redesign/prototype/`, 2,027 lines):
- sources: `nxx/{core,roots,deriv,linalg,linalg_eigen,multiroots,pipes,runtime}.hpp`, `test_core.cpp`, `test_nd.cpp`, `test_runtime.cpp` (the run-time chain test uses no FXT);
- `neg/`: 16 compile-fail tests and 4 probes;
- build and timing scripts for GNU-style compilers and MSVC; `fxt-1.patch`.

An earlier synthesis spike (about 560 lines, not preserved) ran the same API on GCC 16.1, Clang 22 + libc++ (± `-fno-exceptions`), MSVC 19.51 `/W4`, clang-cl 22 `/W4` and em++ 6.0.8 `-fno-exceptions`; its headline chain took 59 iterations and 64 evaluations.

| Result | Value (identical on all 9 configurations unless stated) |
|---|---|
| Headline chain (Newton → secant(5) → bisection), constexpr | x = 1.4142135623730949, 57 iterations, 61 evaluations |
| expand → bisection{x_tol{1e-4}} → Newton, constexpr | 7 iterations, 15 evaluations (masks the x_tol defect) |
| Newton with `deriv::numeric{}` policy in a curried chain, constexpr | Newton wins: 6 iterations, 13 evaluations |
| Brent on x²−2 over [1, 2] | 7 iterations, 9 evaluations |
| `bisection{x_tol{…}}` on x²−2 over [0, 2] | **x = 1.5 after 2–3 iterations (false convergence)** |
| Newton on x²+1 from 0.5 | `budget_exhausted`, best x = 0.00784753, \|f\| = 1.00006 |
| Clamped Newton from 10 on [1, 3]; pinned at edge | 8 iterations; `stalled` with best x = 3 |
| warm_fallback (starved bisection → Newton) | 8 iterations, 16 evaluations in total |
| `derivative_of` sin′(1) error | −1.86e-13 (`central_1_2`), −5.34e-14 (`central_1_4`); constexpr d(x³)/dx(2) within 1e-8 |
| `steps_view` | satisfies `input_range` and `view`; `take(8)` over Brent ends at the intrinsic stop; composes with `views::transform` |
| Sizes | `failure<root_estimate<double>>` 64 B, `result` 72 B, chain 56 B; `result` trivially copyable only on GCC/MSVC |
| Regularity | copy-assignable: closure chain 0, closure `then` 0, solver 1, curried 1, closure `derivative_of` (capturing) 0 |
| N-D | constexpr 2×2 damped Newton (in-house LU) 4 iterations; Eigen `Vector2d` from (1, 1): x = (−1.81626406882515, 0.837367799891248), 8 iterations, 29 evaluations, merit 6.16e-33; from (10, 10): 14 iterations; `VectorXd` 3×3 → (1, 2, 3) in 5 iterations; singular → `errc::singular`; box-pinned → `stalled` |
| Run-time chains (`any_solver`) | `{newton, secant, bisection}` bit-identical to the static chain (57 iterations, 61 evaluations); all-fail and numeric-derivative chains identical too; `sizeof` 32 B (libstdc++), 48 B (libc++), 64 B (MSVC), 24 B (wasm32); no allocation on a call with an existing `fn_t`; test compile + link 1.4–2.4 s per configuration |
| Eigen TU compile + link (`test_nd -DNXX_WITH_EIGEN`) | 3.6–7.3 s per configuration (GCC 6.6, Clang + libc++ 3.6–4.4, em++ 4.2–4.3, MSVC 7.3, clang-cl 3.6); one MSVC `/W4` warning (C4100, an unused lambda parameter on the vector path of the prototype's `evaluate`, which the §5.3 `[[maybe_unused]]` rule covers) |
| FXT-1 | the 2+2-line patch is necessary and sufficient under `-fno-exceptions`; header shadowing does not work |
| Diagnostics | `first_of` mismatch: 1 line (GCC, MSVC); `then` contract: first error line is the message (GCC 77 lines total, Clang 10); `NXX_DELETE` reasons on Clang/clang-cl/em++ only as the prototype wrote it |
| Portability probes | MSVC lacks P2564 (C7595); `is_constructible_v<tolerance<double>, double>` true everywhere; MSVC `/W4` C4459/C4100 in consumer TUs |

**Numerical probes** (not preserved; GCC 16 -O2, on the synthesis spike):
- false x_tol convergence;
- bracket-overflow false success;
- budget exhaustion at 1e-200;
- tan pole reported as a root;
- Newton at noise level burning 102 evaluations;
- (x−1)³ linear rate;
- Newton 0 ↔ 1 cycle costing 102 of 130 chain evaluations;
- expand stepping into sqrt's invalid domain;
- log on [0, 2] rejected;
- 33 % evaluation under-count with an FD derivative.

**Usability probes** (not preserved; GCC 16 and Clang 22, on the synthesis spike):
- run-time braced brackets, budgets and tolerances rejected;
- `make()` inside `.on()`: 161 / 61 lines of errors;
- `is_invocable` a hard error;
- misuse diagnostics of 18–161 lines;
- closures not assignable;
- `failure`-returning callbacks nesting.

## Appendix D: Spike evidence

**Code** (`include/numerixx/`, branch `claude/numerixx-spike`): `core/` 2,167 lines, `roots/` 1,341, `deriv/` 293, and `config.hpp` plus the umbrella headers 121. The spike's tests (the `.cpp` files in `tests/{core,roots,deriv,usage,pipes,compile_fail}/`) are 6,920 lines. Every header with arithmetic in `core/`, `roots/` and `deriv/` is bracketed by `NXX_BEGIN_HEADER`/`NXX_END_HEADER` (§5.3).

**Test matrix** (2026-09-28, each preset configured `--fresh`): 224 CTest tests on `gcc`, `gcc-noexcept`, `clang` (libc++), `clang-asan`, `clang-cl`, `emscripten`, `emscripten-jsexcept`, `emscripten-noexcept` and `emscripten-pthread`; 222 on `msvc`, where the harness self-test `harness_quoted_reason` is not added because cl prints no deletion reasons; 233 on `gcc-multiprecision`; 8 on `integration`. All pass. The 224 are 155 doctest cases, 31 compile-fail cases with their controls (29 illegal states and 2 harness self-tests), the P2564 probe, 3 structural tests, 2 examples and 1 linalg smoke test from phase 0.

**Determinism.** `tests/roots/test_determinism.cpp` holds a golden table of 22 solves. Each row records success or error code, stop reason, iterations, evaluations, x, fx, uncertainty and enclosure, in hexadecimal floating point. The table covers bisection, brent, secant, newton, expand, chains and failures. It is bit-identical on GCC 16.1, Clang 22.1.8 + libc++, MSVC 19.51, clang-cl 22.1.3 and em++ 6.0.8 (CI: 6.0.10) in all four modes. It was bit-identical only after Clang's floating-point contraction was switched off inside the headers (§5.3).

| Result (measured on GCC 16.1; the headline chain and brent counts are asserted on every preset that builds the tests, all but `integration`) | Value |
|---|---|
| Headline chain (§6.11: Newton → secant(5) → bisection), constexpr | x = 1.4142135623730949, 57 iterations, 62 evaluations |
| Newton with the `deriv::numeric{}` policy, then brent, in a curried chain | Newton wins: 5 iterations, 16 evaluations |
| expand from [2, 2.5] → `bisection{width_tol{1e-4}}` → Newton | x = 1.4142135623730949, 17 iterations, 22 evaluations |
| brent on x² − 2 over [1, 2] | 7 iterations, 9 evaluations |

**Compile times** (best of 3, RelWithDebInfo, one TU, measured serially on an idle machine):

| Compiler | Umbrella header (`<numerixx/numerixx.hpp>`, nothing instantiated) | linalg TU (Eigen `partialPivLu` solves, fixed-size and dynamic, called directly; `<numerixx/linalg.hpp>` is still the phase-0 skeleton, so there is no facade cost yet) |
|---|---|---|
| GCC 16.1 | 0.58 s (guard 2 s) | 4.58 s |
| Clang 22.1.8 + libc++ | 0.46 s | 2.55 s |
| MSVC 19.51 | 0.55 s | 3.98 s |
| clang-cl 22.1.3 | 0.53 s | 2.55 s |

**Diagnostics.** On GCC and Clang, each compile-fail case checks that the reason is in the first error. Line counts, as total diagnostic lines / lines of the first error (the harness self-tests: `harness_selftest` 5 / 4 on GCC and 8 / 7 on Clang; `harness_quoted_reason`, which the harness must reject, 8 / 6 and 8 / 7):

| Case | GCC 16.1 | Clang 22.1.8 |
|---|---|---|
| `bisection_given_guess` | 14 / 12 | 35 / 27 |
| `open_int_guess` | 14 / 12 | 32 / 24 |
| `newton_no_derivative` | 13 / 11 | 38 / 30 |
| `newton_mixed_errors` | 27 / 4 | 62 / 24 |
| `bisection_x_tol` | 13 / 11 | 8 / 7 |
| `brent_x_tol` | 13 / 11 | 8 / 7 |
| `bisection_step_tol` | 13 / 11 | 8 / 7 |
| `bisection_with_stop_x_tol` (`with_stop`) | 14 / 12 | 20 / 19 |
| `newton_width_tol` | 13 / 11 | 8 / 7 |
| `secant_width_tol` | 13 / 11 | 8 / 7 |
| `bisection_with_projection` | 14 / 12 | 14 / 13 |
| `secant_with_derivative` | 14 / 12 | 14 / 13 |
| `secant_min_iterations_alone` (`with_stop`) | 14 / 12 | 20 / 19 |
| `bisection_min_iterations_or` (constructor) | 13 / 11 | 8 / 7 |
| `brent_braced_wrong_function` | 14 / 12 | 38 / 30 |
| `bracket_one_end` | 14 / 12 | 35 / 27 |
| `bisection_on_pointer` (`.on`) | 14 / 12 | 20 / 19 |
| `expand_with_stop` | 11 / 9 | 8 / 7 |
| `first_of_without_on` | 12 / 4 | 11 / 7 |
| `first_of_mismatch` | 12 / 4 | 11 / 7 |
| `then_open_to_bracket` | 12 / 4 | 11 / 7 |
| `any_solver_wrong_result` | 12 / 10 | 8 / 7 |
| `bracket_runtime_literal` | 11 / 3 | 11 / 10 |
| `bracket_reversed_literal` | 16 / 10 | 26 / 12 |
| `tolerance_negative_literal` | 21 / 6 | 32 / 15 |
| `max_iterations_zero` | 17 / 15 | 14 / 13 |
| `max_iterations_bool` | 16 / 14 | 8 / 7 |
| `x_tol_zero_zero` | 15 / 7 | 26 / 12 |
| `rel_tolerance_as_tolerance` | 35 / 33 | 20 / 19 |

The prototype's worst case was the `then` contract at 77 lines on GCC; the spike's worst is 62 lines on Clang (mixed callback error types), with the reason in the first of them. The totals include the compile command that the build tool echoes (one line) and any context lines before the first error ("In function …", "In file included from …"). The first-error counts do not include the command, but on Clang they include the trailing "1 error generated." when the first error is the only one. So the counts compare within one compiler's column; they are not absolute. Two harness rules decide what "in the first error" means (§9.1):
- On GCC a failed consteval call can explain itself in a nested error after "in 'constexpr' expansion of" lines, which the harness counts as part of the first error. For a rejected literal, GCC reports the call to the non-constexpr `literal_violates_invariant("…")` and quotes the source line that makes it: in `Tag::reject()` for the refined tolerances (which GCC also names), in the consteval constructor itself for `bracket`, `x_tol`, `width_tol`, `floored_width`, `max_iterations` and `evaluation_budget`. Either way the quoted line holds the reason as a string literal, and that quote counts.
- A deletion reason must be in the compiler's own message. Quoted source lines are left out of the match, so a reason cannot pass merely because a "declared here" note echoes the `NXX_DELETE("…")` line.

**What the spike changed in the design** (each item is in the section cited):
- Floating-point contraction is off inside the headers on Clang, clang-cl and em++ (§5.3, §6.1).
- `better` is a strict weak order, so the static and run-time chains merge to the same best estimate in any fold order (§6.7).
- The driver's success path is `detail::succeed` with `finish(p, sol) -> optional<failure>`, which avoids a GCC 16 false `-Wmaybe-uninitialized` (§6.7).
- Brent's `tol1` and stop reasons, so that `criterion` means the width criterion holds (§6.8, §9.3).
- The pole check takes its reference from the finite initial samples and rejects a non-finite fx, and its known limit is documented (§7.2).
- `expand` saturates at ±max and stalls when both ends are there; its crossing of 0 is a known limit for phase 3 (§7.2).
- secant falls back to x0 − h when x0 + h overflows (§7.2).
- The facades check "F invocable on the scalar" (`callable_v`) and delete the rest with a reason; `with_derivative` and `with_projection` exist only where they apply (§6.6).
- `any_solver`'s reasoned deleted constructor applies only to a solver that is callable with F but returns the wrong result type, so `any_solver` is not a conversion target for every type. Its concept checks invocability before copyability, which avoids a recursive constraint on Clang (§6.10).
- `evaluation_budget` is a class, checked in the source type like `max_iterations` (§6.2).
- A braced list or C array is a bracket only with two ends of a real type; `{x}` had been taken as `{x, 0}` (§6.6).
- MSVC: arrays are left out of the call operators' catch-alls and the `.on` catch-alls take a forwarding reference, so `f, {lo, hi}`, `f, arr` and `.on(arr)` all work on cl and a pointer still gets its reason (§6.6); `deriv::relative` has no `std::optional` member in its consteval constructor (§7.1); F is decayed before `const F&` is formed (§6.6).
- A braced `{lo, hi}` with a function of the wrong signature gets the facade's reason, as a `bracket<T>` or a pair does (§6.6).
- `NXX_DELETE` starts on the declarator's own line, where cl and the floor compilers point; the Clang header brackets save and restore the floating-point setting on every target that supports it (§5.3).
- `expand`'s public `rebuild(options)` accepts only `never{}`, so a searcher's stop criterion cannot be configured through it either (§7.2).
- `min_iterations` is a guard, accepted by a solver only under `&&` with a convergence test, on every path that sets a stop criterion: the constructors, `with_stop` and the public `rebuild(options)` (§6.8).
- `root_estimate` has a constructor that requires x and f(x); the uncertainty of an estimate with neither an enclosure nor a step is inf (§7.2).
- `best_x` is constrained to results whose solution has an x (§6.3).
- The compile-fail harness matches deletion reasons only in the compiler's own message, and a self-test case (`harness_quoted_reason`, which the harness must reject) checks that. The MSVC P2564 probe must fail with C7595 (§9.1).

**Open question found by the spike's audit:** GCC contracts `a * b + c` by default in C++ (ISO mode included) and has no pragma to stop it for a region, so on FMA targets a GCC build may differ from other platforms in the last ulp (§5.3). The presets are unaffected (baseline x86-64 has no FMA). Phase 1 decides between documenting `-ffp-contract=off` for consumers who need bit-identical results and adding it to the GCC interface flags.

**Minor review findings left open** (none reports a false success; each is for the phase that owns the code):
- `clamp_to` is a public aggregate, so reversed or NaN bounds are accepted. The review's probes failed rather than succeeded (Newton from 0.5 with `clamp_to{1.0, 0.0}` on x² − 2 ends with `zero_derivative`). A validated constructor belongs with the options in phase 3.
- Non-finite inputs get different codes: a NaN bracket end is `invalid_input`, a NaN guess `non_finite_input`, a non-finite x in `deriv::diff` `invalid_input`. Phase 1 settles the codes.
- The consteval literal forms `x_tol{abs, rel}` and `width_tol{abs, rel}` take plain `T`, so swapped arguments compile; `make(abs_tolerance, rel_tolerance)` has the role types (§6.2).
- Hosted CI does not upload `compile_time_report.txt` or `compile_fail_report.txt`; the numbers above were measured locally.

**Example.** `examples/quick_tour.cpp` shows the library as it is now: the one-call `solve`, choosing solvers and criteria, run-time configuration through `make()`, open methods with analytic and numeric derivatives, failures as values (no sign change, a pole, an exhausted budget, a fallible callback's own error), `first_of`, `then` and a run-time `any_solver` chain, derivatives, `steps_view` and an observer. It is built with the strict warning flags and runs as a smoke test on every preset except `integration`, which builds the consumer scenarios instead of the examples.

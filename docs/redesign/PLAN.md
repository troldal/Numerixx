# Numerixx 2: redesign plan

- **Date:** 2026-09-27
- **Status:** Approved on 2026-09-28, with every default in section 10 (and DESIGN §12) accepted. Phase 0 and the de-risking spike are done. The spike met its 11 exit criteria locally on all 12 presets (DESIGN §10.2 status, Appendix D), passed hosted CI on PR #3 (ci.yml, and the nightly dispatched on the spike branch on 2026-10-02, DESIGN §12.15), and was merged into master on 2026-10-02 (PR #3, 53d3384) with the nightly floor fixes. Phase 1 is the current phase. The spike's code is kept, and phases 1–3 continue from it (DESIGN §10.3). On 2026-10-04 the user approved the phase-1 core design (DESIGN §12.20): one name per quantity and `nxx::best`, one error code for every non-finite input, role-typed mixed tolerances, checked callback results, `nxx::better_than` as the order of failure payloads, and reasoned deletions in the combinators. It is written into DESIGN (mainly §6, §7.1, §7.2 and §10.3) and is not built yet. On 2026-10-06 the user accepted every recommendation of a simplicity review of that design (DESIGN §12.21): the same guarantees with fewer deleted declarations, combinator states and traits, built in three PRs. Phase 1 is re-estimated at 5–7.5 developer-days [est] (6.5–9.5 on 2026-10-04) and ends with the tag `v2.0.0-alpha.1`.
- **Companion documents:**
  - [`DESIGN.md`](DESIGN.md): the detailed design reference, covering every decision, the code sketches, per-module algorithm tables, CMake, the test strategy and the full roadmap. Section numbers there are stable; "§n" below refers to them.
  - [`prototype/`](prototype/): a throwaway feasibility prototype. It compiles and runs on nine configurations: GCC 16 and Clang 22 + libc++, each with and without `-fno-exceptions`; em++ 6.0.8 with `-fexceptions`, `-fno-exceptions` and `-fwasm-exceptions`; MSVC 19.51; and clang-cl 22.

This plan is the result of four steps:
1. A verified analysis of the current code (master and the unmerged `dev-reorg` branch) and of FXT.
2. Four independent designs: functional core, types first, sceptical pragmatist, and build/migration.
3. A synthesis of those designs.
4. Three adversarial reviews of the synthesis: a prototype-based C++ feasibility review, a numerical-analysis review, and an API-ergonomics review.

---

## 1. Summary

Numerixx 2 is a **general-purpose numerical library in the spirit of GSL, with a smaller scope**: header-only, C++23, MIT-licensed, and portable to WebAssembly. No single application domain drives the design. Every module is rewritten, and the library is built from four ideas:

1. **Solvers are immutable values.** A configured solver such as `roots::brent{}` is a small copyable object. Its algorithm is a pure `init`/`step` pair.
2. **One bounded iteration driver.** It replaces the five hand-written loops. Every exit reports the best estimate, iteration and evaluation counts, and a reason.
3. **One error model.** Every solver result is `std::expected<solution<Est>, failure<Est, UE>>`, and the library never throws. One-shot derivatives, smart constructors and linear solves return a plain `std::expected<T, E>`. A failure carries its cause, the best estimate so far, the counters, and the user's own callback error.
4. **Validated inputs.** Brackets, tolerances and budgets are refined types. Invalid literals are compile errors. Run-time values are validated once, by `make()` or in-band by the solver (an input error code, DESIGN §6.3), and a validated value cannot become invalid afterwards.

On top of these ideas:
- **Stop criteria are typed by what they can soundly judge.** `x_tol` compares successive iterates and does not compile on a bracketing solver, which uses `width_tol` instead.
- **Solvers compose.** `first_of` tries the next solver if one fails. `then` runs search, then solve, then polish. Chains are built at compile time by default; an opt-in type-erased `any_solver` also allows chains assembled at run time, for example from configuration.
- **Several operators return functions:** `derivative_of`, `integral_of`, `inverse_of`, interpolants and `minimizer_of`.

**Dependencies:** vcpkg is deleted. CPM fetches everything that is needed. The library itself needs **no Boost**, **no Blaze**, **no LAPACK/BLAS**, gcem, tl-expected or OpenMP. Linear algebra uses **Eigen 5.0.1**, fetched via CPM, behind a thin facade that returns `std::expected`. Eigen has no BLAS/LAPACK dependency and was verified under Emscripten.

**Scope:**
- v2.0 covers the eight existing areas: derivatives, 1-D roots, 1-D minimisation, polynomials, N-D roots, quadrature, interpolation, and a small linear-algebra facade.
- v2.1 adds N-D minimisation and nonlinear least squares.
- v2.2 adds ODE solvers.
- The core abstractions are designed so these later families reuse the driver, criteria, result types and combinators instead of adding parallel machinery; each adds only its own vocabulary. This is a sketch, not yet compiled (section 2).

**Effort:** about **52.5–77.5 focused developer-days** for v2.0, or 48.5–73.5 without the optional extra root solvers. Rough figures (±50 %) for later work: the planned v2.1/v2.2 families add 28–42 days, and the later families (Chebyshev, series acceleration, linear least squares) another 8–12, so 36–54 in all.

---

## 2. Scope and positioning

Details are in §1.1 and §10.5.

**Where Numerixx sits:**
- **GSL.** GSL is the model for breadth, but it is C, GPL-3.0-or-later, essentially `double`-only, and built around mutable workspace objects and a global error handler.
- **Boost.Math.** Its tools (roots, minima, quadrature) are header-only and generic over the scalar type, but they report errors through policies and exceptions, and they have no composable solver model.
- **Eigen.** Eigen (MPL-2.0) covers linear algebra, and its unsupported module has hybrid Powell and Levenberg–Marquardt solvers. These are MINPACK ports under the Minpack licence, which Numerixx uses only as test oracles.
- **Numerixx's niche:**
  - value semantics and composable solvers;
  - `std::expected` results that carry the best estimate and counters;
  - generic scalar types;
  - an MIT licence;
  - builds that run under Emscripten.

**Scope, using GSL's chapters as the map:**

| Status | Areas |
|---|---|
| **In v2.0** | numerical differentiation; 1-D root finding and bracket search; 1-D minimisation; polynomials (evaluation, algebra, roots); multidimensional root finding; 1-D quadrature (adaptive Gauss–Kronrod, tanh-sinh family, Romberg, Gauss–Legendre); interpolation (linear, cubic splines, monotone, rational); a small dense linear-algebra facade over Eigen |
| **Planned after v2.0** (new downstream modules; linear least squares extends `fit`) | v2.1: `multimin` (Nelder–Mead, BFGS/L-BFGS, nonlinear conjugate gradient) and `fit` (Levenberg–Marquardt nonlinear least squares). v2.2: `ode` (Dormand–Prince RK45 with dense output, then a stiff Rosenbrock/BDF method). Later: `chebyshev` (Chebyshev approximation), `series` (Richardson, Wynn ε, Levin u), and linear least-squares fitting in `fit` |
| **Candidates** (not scheduled; can be added without changing the core) | QUADPACK's weighted and oscillatory rules (QAWC, QAWS, QAWO, QAWF), QAGS extrapolation, QAGP break points and CQUAD; fixed Gauss rules other than Gauss–Legendre (Gauss–Laguerre, Gauss–Hermite, Gauss–Jacobi, …); 2-D interpolation (bilinear, bicubic); B-splines |
| **Out of scope** (use the standard library, Eigen or specialised libraries) | special functions (Boost.Math, or C++17 `<cmath>` where the standard library implements it), random numbers, quasi-random sequences and distributions (`<random>`), statistics, histograms and N-tuples, FFT and filtering, sparse matrices, eigensystems beyond Eigen, Monte Carlo integration, simulated annealing, wavelets, Hankel transforms, physical constants |
| **Not needed** (the standard library or Eigen covers it) | vectors and matrices, BLAS, permutations, combinations and multisets, sorting, complex numbers, elementary maths functions, IEEE utilities |

**Why the architecture extends** [sketch]:
- A minimiser is a solver over an N-D state with an extremum estimate.
- Nonlinear least squares reuses the system machinery with a residual vector.
- An ODE integrator is a solver whose step advances t with an error-controlled step size. Reaching t_end is a stop criterion, and dense output is a function-returning API: the solution is a callable.
- The driver, the criteria algebra, the result and failure types, the combinators and `steps_view` (which gives trajectories) apply unchanged. For example, `first_of(nonstiff, stiff)` is a meaningful chain for an integration.
- Each family is a new module downstream in the DAG, so no v2.0 module gains a dependency.

**FLAG, licensing.**
- Numerixx is MIT, so no GPL code may be ported, paraphrased or copied: nothing under the GPL, LGPL or AGPL, such as GSL or MPSolve (decided on 2026-09-30). Code may come only from Numerixx 1.x and the Boost.Math code DESIGN names; Numerical Recipes listings and code without a licence are excluded too (decided on 2026-10-02).
- Algorithms are implemented from the literature, with references cited in each header.
- The test oracles are Boost.Math, Eigen and high-precision reference tables. GPL projects (such as GSL) are not oracles, although values published in their documentation may serve as reference facts.
- Boost-derived code (Brent, TOMS748) keeps its BSL-1.0 notice.

---

## 3. What the analysis found

These findings are why incremental repair is not worth it.

- **Mutable CRTP state machines everywhere.**
  - There are five hierarchies (bracketing, polishing, searching, integration, multiroots), each with public setters and hidden flags.
  - There are five near-identical driver loops that share an off-by-one `maxiter` bug.
  - `multisolve` reports success when it runs out of iterations, and prints to `std::cout`.
- **The two "function-returning" APIs are broken.** `derivativeOf` and `integralOf` default-construct the function type and discard their argument, so only captureless lambdas work.
- **Incoherent error model.**
  - There are five unrelated error types.
  - There are 21 `throw`s inside APIs that return `expected`.
  - A Boost stack trace is captured on every error value, even for errors that are never thrown.
  - Errors from different solvers have different types, so solvers cannot be chained.
- **Weak termination.** Every test is an absolute residual check. There is no x or width tolerance, and sign tests are products that can underflow. Romberg and Simpson falsely converge on sin²(8πx).
- **Illegal states are accepted:**
  - unordered brackets;
  - `eps ≤ 0`, `maxiter ≤ 0`, and a `bool` passed as `maxiter`;
  - an empty system that "converges";
  - quadratics that return NaN roots as success;
  - a zero polynomial that cannot be told apart from a constant.
- **Confirmed correctness bugs:**
  - polynomial `operator-` when deg(lhs) < deg(rhs), and the test asserts the wrong answer;
  - `hessian()` is mathematically wrong;
  - the derivative step ignores the sign of x;
  - the default second and mixed derivatives are wrong by O(1).
- **Build rot.**
  - The vcpkg baseline is broken.
  - Boost is linked into every module only for a concept trait.
  - gcem, nlohmann-json, boost-math and openblas are dead dependencies.
  - The lowercase `integrate/`/`interpolate/` include paths do not match the directories on case-sensitive filesystems.
  - The `poly` and `roots` targets link each other in a cycle.
  - Only 2 of 4 test files compile.
- **The current code is on `dev-reorg`, not master.** It is 42 unmerged commits (2025-01-14) with C++23, `std::expected` in some modules, a real optimize module and a larger interpolate module. It also vendors 38.8 MB of Blaze and uses `std::stacktrace`, which is unavailable on Emscripten, and it adds new regressions. It is the richest source to port algorithms from, but not a base to build on.
- **Scalar genericity is claimed but never tested.** No test instantiates multiprecision or complex types, and several modules are hard-wired to `double`.

---

## 4. Your goals: what the plan does, and where I flag them

| Your goal | Delivered? | How, and the flag |
|---|---|---|
| **Functional style, possibly with FXT** | Yes, with nuance | Values in, values out: solvers, problems, criteria, results and errors are copyable values, and steps are pure functions. **FLAG:** purity holds at function boundaries only. Inside a step, Horner's rule or a tridiagonal solve, the code is an ordinary loop on locals, which is the only way to get competitive numerics in C++. **FLAG, FXT:** the library's own code uses only the `std::expected` members, and FXT is the optional user-facing pipe vocabulary (`<numerixx/pipes.hpp>`, target `numerixx::pipes`). The reasons: FXT has none of the combinators numerics needs yet (first success, bounded iteration, Kleisli composition, refined types), and today it does not compile under `-fno-exceptions` (`throw 0;` in two concept headers; a 2+2-line fix, FXT-1). Once the combinators are proven in Numerixx, their generic parts move to FXT (§8). |
| **Chain solvers: if one fails, the next tries** | Yes | `nxx::first_of(s1, s2, s3)`. Solvers with different inputs are unified by currying: `.on(input)` turns every solver into a callable `f → result` with the same result type. `then` (staging), `warm_fallback` (restart from the previous best) and `with_evaluation_budget` (one budget for a whole chain) complete the set. By default, chains are templates built at compile time: zero overhead, usable in `constexpr`. For chains assembled at run time, the opt-in `nxx::any_solver` type-erases a curried solver for a fixed callable type, and `first_of` accepts a range of them. **FLAG:** `any_solver` is built on `std::function`, because `std::move_only_function` is missing from libc++ 22 and hence from Emscripten. It may therefore allocate when a solver is wrapped or copied (never on a call with an existing `F`), and it costs an indirect call per alternative and per evaluation. **FLAG:** chaining only pays off when the methods fail in different ways. |
| **Functions as return values** | Yes | `deriv::derivative_of(f)`, `integrate::integral_of(f)`, `integrate::antiderivative(f, a, rule)`, `roots::inverse_of(f, bracket)`, interpolants, `optimize::minimizer_of(family, bracket, solver)`, `poly::derivative(p)`, and later ODE dense output. Each returns a named, copy-assignable callable. There are two conventions (§6.4). Scalar-valued ones such as `derivative_of` return `expected<T, fault>` and plug straight into solvers as callbacks. Estimate-valued ones such as `integral_of` or `inverse_of` return a full `result` (value, error estimate, counters) and become callbacks through `fn::value_of` [sketch]. **FLAG:** these recompute on every call (no memoisation), and they can fail. |
| **No exceptions; use optional/expected** | Yes | The library never throws and builds with `-fno-exceptions`. `expected` is for failures with a reason; `optional` is only where absence is the answer (for example `polynomial::degree()` of the zero polynomial). **FLAG:** exceptions thrown by user callbacks (for example from third-party code) propagate untouched. The library is exception-neutral with conditional `noexcept`; a blanket `noexcept` would call `std::terminate` on such a throw. `bad_alloc` from owning types (polynomials, interpolants, Eigen dynamic storage, `any_solver`) is not converted. Programmer errors (violated preconditions of the low-level stepping protocol) are assertions, not `expected`. |
| **Immutable objects** | Yes, by interface | No mutators, private members, and `with_*()` builders that return modified copies. **FLAG:** not `const` data members. They delete copy and move assignment, which breaks `std::expected`, containers and the driver itself. There are no exceptions to the rule. N-D solver states are values over Eigen vectors: fixed-size when the dimension is known at compile time, otherwise heap-allocated and copied per step. For non-trivial functions that cost is negligible next to function evaluations. |
| **Illegal states unrepresentable** | As far as C++ allows | **Tier A, compile time:** `nxx::bracket{2.0, 1.0}`, `nxx::max_iterations m = 0` or `= true`, and a negative tolerance literal do not compile. Neither do bisection given a guess, Newton without a derivative source, a width criterion on an open method, or `x_tol` on a bracketing solver (its guarantee would be false there). **Tier B, construction time:** `make() → expected` validates run-time values once. **Tier C, run time:** NaN mid-iteration, a singular Jacobian, stalls, poles and exhausted budgets are states no type can exclude, so they come back as errors. **FLAG:** in real applications almost every bracket and tolerance is a run-time value, so in practice the guarantee is "validated once at the boundary, never invalid afterwards". **FLAG:** MSVC does not implement P2564 (consteval escalation), so functions that forward tolerances must take the refined type, not a raw `double`. |
| **Fetch Boost via CPM** | Yes in form, but the library needs no Boost | **FLAG:** Boost is used today for exactly three things: one concept trait (`IsFloat`), a stack trace and a 2-D table. None of them justifies the dependency. Stack traces have no Emscripten backend and would be paid on every failed attempt in a fallback chain. Boost Math is not used at all. Plan: `IsFloat` becomes an open `scalar_traits`, the stack trace is dropped, and the table becomes two rows of an array. CPM fetches the **standalone** boostorg `config` + `multiprecision` (+ `math` for test oracles) 1.92.0 repositories, about 5 MB instead of 108 MB, and **only** for an optional `numerixx::multiprecision` adapter and the tests. They are never fetched under the CPM package name `Boost`, so that they cannot clash with a parent project that fetches full Boost. |
| **No vcpkg** | Yes | `vcpkg.json` and every vcpkg assumption are deleted. CPM 0.43.2 is committed with a pinned hash, and it reuses a parent project's CPM if one exists. |
| **Replace Blaze with a BLAS/LAPACK-free, Emscripten-capable linalg** | Yes | **Eigen 5.0.1**, fetched via CPM, is the backend of `numerixx::linalg` and `numerixx::multiroots`, and later of `multimin`, `fit` and the stiff `ode` solver. It is header-only with no BLAS/LAPACK dependency, and it was verified under em++ 6.0.8 (including `-fno-exceptions`) and with Boost multiprecision scalars. A thin `nxx::linalg` facade (LU, QR and Cholesky solves) returns `std::expected`. It checks dimensions and finiteness, and `lu_solve` also checks the condition estimate (rcond against n·ε), because Eigen's partial-pivot LU never reports singularity by itself. The facade returns concrete types only. **FLAG:** Eigen's expression templates dangle under `auto`, which clashes with an `auto`-heavy functional style; the concrete return types contain that. **FLAG:** Eigen costs compile time, measured at +2.5–6.4 s per translation unit once LU and N-D Newton are instantiated. That cost stays inside linalg/multiroots translation units, and users who need only the scalar modules set `NUMERIXX_WITH_LINALG=OFF` and never download Eigen. The spline tridiagonal solves stay hand-written: they are O(n) algorithms, not a library need. |

---

## 5. Architecture at a glance

### 5.1 Modules and dependencies (a DAG, no cycles)

```
core ─┬─► roots
      ├─► optimize
      ├─► poly
      ├─► integrate
      ├─► interpolate            (own O(n) tridiagonal solves)
      ├─► deriv ─────────┐
      ├─► linalg ────────┴──► multiroots   (linalg = facade over Eigen 5.0.1)
      └─► pipes (+ FXT)
adapter (leaf): multiprecision ─► core (+ standalone Boost; Eigen glue when linalg is enabled)
later:  multimin ─► optimize, multiroots;  fit ─► multiroots (+ poly once linear least squares lands);
        ode ─► multiroots;  chebyshev ─► poly;  series ─► core
```

- **Layout.** There is one include root, `<numerixx/…>`, and one thin INTERFACE target per module, `numerixx::<module>`: namespace = header = target.
- **Dependencies.** `core` depends only on the standard library. Eigen is fetched only when `linalg`/`multiroots` are enabled (`NUMERIXX_WITH_LINALG`, default ON). The umbrella header `<numerixx/numerixx.hpp>` leaves out linalg and multiroots, so Eigen's compile cost appears only where you include them explicitly. A layering test enforces the DAG.
- **The old cycles disappear:**
  - `roots` recognises derivatives structurally, so it no longer depends on `deriv` or `poly`;
  - `poly` does its own Newton polish;
  - `multiroots` no longer calls the polynomial solver to find a parabola vertex.

### 5.2 Core vocabulary (details in §6)

| Concept | Shape |
|---|---|
| Scalars | open, user-specialisable `nxx::scalar_traits<T>` (default: any inexact, non-integer type with `numeric_limits`; a user type needs a specialised `numeric_limits` too); every default tolerance is an expression in `T`, so it is attainable in `float`, `long double` and multiprecision |
| Refined inputs | `tolerance<T>`, `abs_tolerance<T>`, `rel_tolerance<T>`, `max_iterations`, `evaluation_budget`, `bracket<T>`, `sign_bracket<T>`, `interval<T>`; `consteval` literal constructors plus `make() → expected` |
| Solver protocol | `prepare(f, input)`, `init(p)`, `step(p, s)` (pure), `view(s)`, `estimate(s)`, `best(s)`, `intrinsic(s)`, `finish(p, r)` |
| Driver | `nxx::iterate(alg, problem[, observer])`: the only loop in the library |
| Stop criteria | values combined with `\|\|`/`&&`, **typed by what they can judge**: `x_tol`/`step_tol` for iterates, `width_tol`/`floored_width` for enclosures, `f_tol`, `max_evaluations`, `custom{λ}`; the iteration budget is separate and mandatory |
| Results | `result<Est, UE> = std::expected<solution<Est>, failure<Est, UE>>`; errors are cheap to copy and hold no heap memory of their own (a dynamically sized system's best estimate owns its vectors) |
| Errors | `errc` in ranges: input (1–31), numerical (32–63), callback (64+) |
| Combinators | `first_of`, `first_of_with(policy, …)`, `then`, `warm_fallback`, `with_evaluation_budget`: named, assignable class templates, usable in `constexpr`; plus the opt-in type-erased `any_solver` and `first_of(range)` for chains built at run time |
| Stepping and observation | `steps_view`, a lazy range of iteration states for manual stepping and tracing (composes with `std::views::take`/`transform`); `.with_observer(fn)` for logging; `.with_projection(p)` maps each proposed iterate into a box or domain before evaluation |

### 5.3 What it looks like

The spellings below are the design's (§6.11). The prototype compiled and ran the same chain and pipeline on all nine configurations, including inside a `static_assert`. It used older spellings: `r::secant{nxx::default_step{}, 5}`, `r::expand_out` and `r::bisection{nxx::x_tol{1e-4}}`. The design rejects the first and the last, and names the searcher `r::expand`. The spike compiles and runs this snippet as written (checked with GCC 16), and its tests cover each part on every preset (spike exit criterion 8).

```cpp
#include <numerixx/roots.hpp>
namespace r = nxx::roots;

constexpr auto f  = [](double x) { return x * x - 2.0; };
constexpr auto df = [](double x) { return 2.0 * x; };

// Three different solvers with three different inputs; if one fails, the next one tries.
constexpr auto chain = nxx::first_of(
    r::newton{}.with_derivative(df).on(0.0),        // fails: f'(0) = 0 -> errc::zero_derivative
    r::secant{}.with_budget(5).on(0.0),             // fails within its 5-iteration budget
    r::bisection{}.on(nxx::bracket{0.0, 2.0}));     // succeeds; the literal bracket is checked at compile time
static_assert(chain(f).has_value());                // the whole chain can run at compile time

auto res = chain(f);            // res->x == 1.41421356..., res->used counts the failed attempts too
auto x   = nxx::best_x(res);    // optional<double>: the solution, or the best estimate of a failure

// Staging: grow a window until f changes sign, get a coarse enclosure, polish with Newton.
constexpr auto pipeline = nxx::then(r::expand{}.on(nxx::bracket{2.0, 2.5}),
                                    r::bisection{nxx::width_tol{1e-4}},
                                    r::newton{}.with_derivative(df));

// Run-time inputs are first-class and validated in-band.
auto r1 = r::brent{}(f, {lo, hi});
auto r2 = r::brent{nxx::width_tol{1e-12}}.with_budget(60).on(nxx::bracket<double>::make(lo, hi))(f);
```

Everyday calls use the one-call facade, which is a thin composition of the core. The full list is in §6.13, and the normative canonical calls are in §6.14:

```cpp
nxx::roots::solve(f, {lo, hi});                           // bracketing default: Brent provisionally, fixed by corpus counts in the roots phase
nxx::roots::solve(f, df, x0);                             // expand + safeguarded Newton (rtsafe)
nxx::optimize::maximize(f, {0.0, 3.0}, nxx::optimize::golden{});    // -> extremum{x, fx}, never a bare double
nxx::deriv::central(f, x);                                // optimal relative step -> expected<double, fault>
nxx::integrate::quad(f, a, b);                            // adaptive Gauss-Kronrod 7-15; a > b is legal
nxx::interpolate::make_cubic_spline(xs, ys);              // expected<cubic_spline<double>, errc>
nxx::multiroots::solve(F, std::array{1.0, 0.5});          // damped Newton with an FD Jacobian until dogleg + Broyden pass the phase-5 corpus
```

Misuse is rejected at compile time with a reason. On GCC ≥ 15 and Clang ≥ 19 the reason is designed to appear in the first error; spike exit criterion 9 and the compile-fail tests check this. On the floor compilers and MSVC, the error names the file and line of the deleted declaration, and the reason starts on that line (`NXX_DELETE` is written on the declarator's own line; MSVC prints the location, not the source):

```cpp
r::newton{}(f, 1.0);                  // "newton needs a derivative: .with_derivative(df), .with_derivative(deriv::numeric{}), ..."
r::bisection{}(f, 1.0);               // "bracketing solvers need a bracket: pass {lo, hi}, nxx::bracket<T>::make(a, b), ..."
r::bisection{nxx::x_tol{1e-9}};       // "x_tol and step_tol compare successive iterates; bracketing methods converge on the enclosure: use width_tol{abs[, rel]} or floored_width{}"
nxx::bracket{2.0, 1.0};               // invalid literal
```

---

## 6. Module by module

The details are in §7. Each module's port includes a regression test for every confirmed bug it must not reintroduce.

| Module | Keep / fix | Add | Drop |
|---|---|---|---|
| **deriv** | Stencils become integer data. Duplicates are removed. **Per-stencil optimal relative step**, which fixes the O(1)-wrong second and mixed derivatives. The step is sign-aware: `factor·max(\|x\|, typical)`. | `diff_with_error`, Ridders extrapolation with an error estimate, noise-aware steps for functions computed with limited precision or by inner iterative solvers, `mixed`, `derivative_of`, and the `numeric{}` derivative policy for Newton | template-template `diff<ALGO>`, gcem |
| **roots** | Bisection with an overflow-safe midpoint, and representation-space bisection for extreme magnitudes. Illinois (Anderson–Björck variant), Ridders, Newton, **derivative-free** secant. | **Brent** (provisional default), **rtsafe** (safeguarded Newton). A pole check so that `tan` on [1, 2] is not "a root". Cycle and divergence detection for open methods. One `expand`/`scan`/`subdivide` search family whose output *is* a bracketing input. `inverse_of`, `steps_view`. Later and optional (phase 8): **TOMS748**, **ITP**, safeguarded Halley, and Steffensen (low priority; deleted if the corpus shows no benefit). | `fsolve`/`fdfsolve`/`search` drivers, CRTP bases, the `ResultProxy`/`StopToken` machinery, complex root solving (moved to poly) |
| **optimize** (1-D) | Golden section and the Brent minimiser (from dev-reorg), both ending on Brent's width test; `bracket_minimum` | `extremum{x, fx}`, a `maximizing(solver)` adaptor, `minimizer_of` | 1-D gradient descent, FD-of-FD Newton |
| **poly** | Canonical polynomial ("empty = zero, else leading ≠ 0"), Horner, algebra with `operator-` fixed, total `divmod`, stable quadratic and cubic that return sum types (never NaN), exact-zero trimming | deterministic **Aberth–Ehrlich** with polishing of every root, `antiderivative`, `compose`, `eval_with_derivatives` | Laguerre with `random_device`, the dependency on `roots` |
| **linalg** *(new)* | — | a thin facade over Eigen 5.0.1: vector and matrix aliases (fixed-size when the dimension is known at compile time), `lu_solve`/`qr_solve`/`cholesky_solve` returning `expected` (dimension and finiteness checks; `lu_solve` also checks rcond, `cholesky_solve` positive definiteness, and `qr_solve` reports the numerical rank), with `std::array`/`std::vector` accepted at the boundary | Blaze, LAPACK, OpenMP |
| **multiroots** | Damped Newton, now with scaling, a full-step convergence test, Armijo backtracking, distinct `singular`/`line_search_failed`/`stalled`/`local_minimum` errors, box projection, and the best iterate on every exit | `gradient_of`, `jacobian_of` (finite differences) and a correct `hessian_of` in `multiroots/derivatives.hpp` (they return Eigen types); Broyden, Powell dogleg (MINPACK rules); `system_of`; `solve` switches to dogleg + Broyden once they pass the phase-5 corpus | `MultiFunction(Array)`, steepest descent as a solver, "success at maxiter", `std::cout` |
| **integrate** | Trapezoid, Romberg and Simpson as one sample-reusing family with `min_iterations`, `(1 << level)` bounded by type, and a > b legal | **adaptive Gauss–Kronrod 7-15** with QUADPACK error rules (default), tanh-sinh / exp-sinh / sinh-sinh, `gauss_legendre` (plain function, no error estimate), `quad`, `integral_of` (fixed), `antiderivative` | `boost::multi_array`, CRTP bases |
| **interpolate** | Linear, natural/clamped splines, the monotone Hermite from dev-reorg (renamed PCHIP and fixed at the last knot) | not-a-knot and periodic splines, true Steffen, Floater–Hormann; knots validated once; coefficients computed eagerly (no `mutable` caches, so no data races); the out-of-range policy is a type | `makepoly` (Vandermonde + LAPACK), plain equispaced Lagrange |
| **func** (folded into core) | — | a few function adaptors in `nxx::fn` (`extend_linearly`, `counted`, `catching`) | the module and its unused `Function` wrapper |

**Genericity** (§3.5):
- `float`, `double` and `long double` are tested everywhere.
- Multiprecision (`cpp_bin_float_50`) works through the open scalar trait and gets an optional CI leg. It is used with Eigen through Boost's `eigen.hpp` glue.
- Complex numbers are supported only in `poly`.

**Testing against published suites** (§9.2): Alefeld–Potra–Shi for bracketing root finders, Moré–Garbow–Hillstrom for nonlinear systems (later also least squares and minimisation), and QUADPACK/Piessens and Bailey–Borwein-style test integrals (endpoint singularities, interior peaks, oscillatory integrands, infinite ranges). Reference values are generated at high precision; Boost.Math and Eigen serve as oracles.

---

## 7. Build and dependencies

The details and the full CMake sketch are in §4.

| Dependency | How | Needed by | Default |
|---|---|---|---|
| CPM 0.43.2 | committed `get_cpm.cmake`, hash-pinned; reuses a parent's CPM | build | always |
| FXT | pinned commit + SHA256, `DOWNLOAD_ONLY`, own `fxt::fxt` shim; local override `-DCPM_FXT_SOURCE=<path>` | `numerixx::pipes` only | `NUMERIXX_WITH_FXT=ON` |
| Eigen 5.0.1 | CPM, `DOWNLOAD_ONLY`, own target (skipped if a parent provides `Eigen3::Eigen`) | `numerixx::linalg`, `numerixx::multiroots` | `NUMERIXX_WITH_LINALG=ON` (turn off for scalar-only use) |
| Boost.Config + Multiprecision 1.92 (standalone) | CPM, non-`Boost` package names | `numerixx::multiprecision` adapter | OFF |
| Boost.Math 1.92 (standalone) | CPM | oracle tests only | OFF |
| doctest 2.5.3 | CPM, hash-pinned | tests | top-level only |
| google/benchmark 1.9.5 | CPM, hash-pinned | benchmarks | OFF |

- **Deleted:** `vcpkg.json`, gcem, tl-expected, Blaze, LAPACK, OpenBLAS, OpenMP, nlohmann-json, fmt, hwinfo, Boost.Stacktrace, Boost.MultiArray, the 172-file vendored Google Benchmark copy, and the `.idea` toolchain paths.
- **Good-subproject rules:**
  - tests, examples, benchmarks and docs are OFF when Numerixx is not the top-level project;
  - no global flags;
  - a parent's CPM, FXT and Eigen targets are reused;
  - CI tests a CPM parent and a FetchContent parent in both declaration orders, plus install + `find_package`.
- **Presets** for MSVC, clang-cl, GCC, Clang + libc++, ASan, no-exceptions, multiprecision, and Emscripten. Emscripten gets four presets: `emscripten` (`-fwasm-exceptions`), `emscripten-jsexcept` (JavaScript-based `-fexceptions`), `emscripten-noexcept`, and `emscripten-pthread` (`-fwasm-exceptions -pthread`).
- **CI:** GitHub Actions runs Windows (MSVC, clang-cl), Linux (GCC 16, Clang 22 + libc++), Emscripten tests under node, the consumer-build job, and a nightly floor-compiler job.
- **Windows MAX_PATH pitfalls** were reproduced and need workarounds:
  - keep a short `CPM_SOURCE_CACHE` (for example `C:\cpm`);
  - keep `EM_CACHE` short;
  - set `CMAKE_POLICY_DEFAULT_CMP0168=NEW`. A plain `cmake_policy()` does not reach CPM.

---

## 8. Roadmap

The details are in §10. Sizes are focused developer-days for one developer.

| # | Phase | Scope | Size |
|---|---|---|---|
| 0 | Skeleton | Tag `v1.0.0` (master, 5de1e07) and `v1.1.0-legacy` (dev-reorg tip, 8528e94) so that existing users can pin the old API; new CMake, presets and every CI leg (including multiprecision); delete the old tree | 2.5–3.5 |
| S | **De-risking spike** | Hosted CI green on all legs; the FXT-1 fix pinned; chains with fallible callbacks under clang-cl in CMake builds; CPM deduplication with a parent project in both declaration orders; the umbrella-header compile-time guard and a recorded linalg TU time; your decision on §12 items 1–8; criterion soundness; the canonical calls with run-time inputs; readable compile-fail diagnostics; regularity; derivative composition (11 exit criteria in §10.2; the §12 decisions were made on 2026-09-28) | 3–4 |
| 1 | Core vocabulary | scalar traits and maths helpers, refined types, error and result types, evaluation, `pipes`; plus the core design approved on 2026-10-04 and revised on 2026-10-06 (DESIGN §12.20, §12.21, §10.3) | 2.5–3.5, re-estimated on 2026-10-04 to 6.5–9.5 and on 2026-10-06 to 5–7.5 [est] |
| 2 | deriv | stencils, steps, `diff`, `diff_with_error`, `ridders`, `mixed`, `derivative_of`, the `numeric` policy | 3–5 |
| 3 | Driver + 1-D roots | criteria, driver, combinators, `any_solver` run-time chains, `steps_view`; bisection, Brent, Illinois, Ridders, rtsafe, secant, Newton; expand/scan/subdivide; `solve`, `inverse_of`; Alefeld–Potra–Shi suite | 10–14 |
| 4 | optimize (1-D) | golden, Brent-min, bracket_minimum, `maximize`, `minimizer_of` | 2.5–4.5 |
| 5 | linalg + multiroots | Eigen facade (`lu_solve`, `qr_solve`, `cholesky_solve`); gradient, FD Jacobian, true Hessian; damped Newton with scaling and projection; Broyden, dogleg; Moré–Garbow–Hillstrom suite | 8–15 |
| 6 | poly | polynomial, closed forms, Aberth | 4–6 |
| 7 | integrate + interpolate | G7K15, tanh-sinh family, Romberg; splines, PCHIP, Steffen, Floater–Hormann; QUADPACK-style test integrals | 8–11 |
| 8 | Optional roots | TOMS748, ITP, Halley, Steffensen; re-run the default-solver corpus | 4 |
| 9 | Multiprecision, docs, release | multiprecision adapter, docs, examples, benchmarks → `v2.0.0` | 5–7 |

- **Total for v2.0:** 52.5–77.5 days, or 48.5–73.5 without phase 8. §10.3 shows the arithmetic.
- **After the spike** (decided on 2026-09-29): phases 1–3 continue from the spike's code, with unchanged scope and acceptance criteria, except that phase 1 gained the core design approved on 2026-10-04 (DESIGN §12.20) and revised on 2026-10-06 (DESIGN §12.21), phase 2 gained the deriv items decided with it, not re-estimated, and phase 3 builds the items that design leaves to it, size unchanged (DESIGN §10.3). About 5–7.5 [est, re-estimated on 2026-10-06; 6.5–9.5 on 2026-10-04; 0.5–1 before], 2–3.5 and 6–9 days of them are left, and 44.5–67.5 days in all for phases 1–9 (40.5–63.5 without phase 8) [est; 46–69.5 and 42–65.5 on 2026-10-04]. DESIGN §10.3 lists what is done and what is left. The v2.0 total above is the 2026-09-28 plan and leaves out the work added to phase 1.
- **Order:** deriv comes right after the core vocabulary and before the driver, because it is small, needs only core, and is used by the numeric-derivative Newton in phase 3 and the FD Jacobians in phase 5. Phases 0 → S → 1 → 2 → 3 are sequential. After that, phases 4–7 depend only on core and the driver, so a second developer can run 6 and 7 in parallel with 4 and 5. Phase 8 may follow `v2.0.0`.
- **Start fresh on this branch.** Port algorithm bodies mostly from dev-reorg, with provenance noted in each commit. Never merge dev-reorg: it would put 38.8 MB of Blaze into history for good. The `prototype/` headers are the starting point for the core; its in-house LU (`nxx/linalg.hpp`) and zero-heap choices are not carried over (§10.1).
- **Migration from 1.x** (§10.4): `MIGRATION.md` maps the old API to the new one. For example:
  - `fsolve<Bisection>` becomes `roots::bisection`;
  - `fdfsolve<Newton>(f, df, x0)` becomes `newton{}.with_derivative(df)(f, x0)`;
  - `.result(fn)` becomes `transform`;
  - `fminimize` becomes `optimize::minimize` returning `extremum{x, fx}`;
  - `diff<ALGO>(f, x)` becomes `deriv::diff(f, x, stencil)`.

  There is no compatibility shim, because it would reproduce the bugs the redesign removes. Results change where 1.x was wrong: default second and mixed derivatives, success at maxiter, NaN roots.

**Beyond v2.0** (§10.5; rough, ±50 %):

| Release | Module (depends on) | Content | Size |
|---|---|---|---|
| v2.1 | `multimin` (optimize, multiroots) | Nelder–Mead, BFGS/L-BFGS with a strong-Wolfe line search, nonlinear CG; new `extremum_nd`, gradient-norm criterion | 8–12 |
| v2.1 | `fit` (multiroots) | Levenberg–Marquardt (trust-region form) on `qr_solve`; `least_squares_estimate` with covariance | 6–9 |
| v2.2 | `ode` (multiroots) | Dormand–Prince RK45 with dense output (the solution as a callable), then a stiff Rosenbrock/BDF method | 14–21 |
| later | `chebyshev` (poly), `series` (core), linear least squares in `fit` (adds fit → poly) | Chebyshev series; Richardson / Wynn ε / Levin u; basis-function fitting | 8–12 |

---

## 9. FXT work

The FXT items in priority order (details in §8):

1. **FXT-1 (do first; it is tiny):** replace `throw 0;` with `std::unreachable();` in `concepts/IsExpected.hpp:87-88` and `concepts/IsOptional.hpp:76-77`, and guard the throwing utilities (`attempt`, `failure`, `lazy`, formatting, enums) with `#if __cpp_exceptions`. The prototype showed that this 2+2-line change is necessary and sufficient for the pipes it exercised (`transform`, `and_then`, `value_or`, `match`, `tap`) under `-fno-exceptions` on GCC, Clang, em++ and clang-cl `/EHs-c-`. It matters only for the no-exceptions build mode with pipes. The patch is in `prototype/fxt-1.patch`. **Status:** the probe fix is troldal/FXT#1, which Numerixx pins; the guards are still to do.
2. **FXT-2:** CMake hygiene:
   - reuse a parent's CPM instead of downloading an unhashed one;
   - fetch tl-expected/tl-optional only when their options are ON;
   - add `cxx_std_23`, install/export and version tags.
3. **FXT-3:** support `std::expected<void, E>`, or document `fxt::unit`.
4. **FXT-4..7** come after the roots phase has settled their semantics:
   - generic `first_of_with`/`first_of`;
   - `collect_errors`;
   - `kleisli`;
   - `first_success` over a range.
5. **FXT-8 (optional) and FXT-9 (any time):** a lazy unfold view and `value_or_else`; self-contained headers, include guards, the `demos/LogicalOr.cpp` fix, and documentation of `operator||`.

---

## 10. Decisions

**Decided on 2026-09-28: every default below was accepted**, and so were the other items of §12 (19 in all, numbered in brackets here). These are the ones that change the most. The written decision on §12 items 1–8 was a spike exit criterion (§10.2, criterion 6), which this settles.

0. **Starting point** (decision D1, already settled in DESIGN; listed for visibility). A fresh tree on this branch, porting from dev-reorg and tagging it as legacy, rather than merging dev-reorg. **Default: fresh tree.**
1. **v2.0 scope** [§12.7]. The eight current modules in v2.0; `multimin` and `fit` in v2.1; `ode` in v2.2. The candidates in section 2 stay unscheduled, and the out-of-scope areas stay out. **Default: accept.**
2. **Linear algebra** [§12.2]. Eigen 5.0.1 behind the `nxx::linalg` facade, or hand-written kernels. Eigen costs +2.5–6.4 s per linalg/multiroots TU that uses it, and N-D solves cannot be `constexpr`. Hand-written kernels are constexpr and allocation-free, but they are more code to write and validate, and the v2.1 families would need QR too. **Default: Eigen behind the facade.**
3. **Boost** [§12.1]. None in the library; standalone Boost via CPM only for an optional multiprecision adapter and the test oracles. **Default: accept.**
4. **Reach of "illegal states unrepresentable"** [§12.3, §12.4]. Compile-time rejection for literals and for mismatched solvers, inputs and criteria (`x_tol` on a bracketing solver is a compile error with a reason); run-time inputs validated once in-band. **Default: accept.**
5. **API break and naming** [§12.5]. snake_case names; namespace = header = target for every module (`nxx::<module>`, `<numerixx/<module>.hpp>`, `numerixx::<module>`); builders (`secant{}.with_budget(5)`) instead of positional arguments; `integrate::quad`; no compatibility shim. **Default: accept.**
6. **Scale defaults** [§12.6]. Derivative steps are relative to |x| (1 at x = 0). Stopping tolerances keep an absolute floor at scale 1, so roots at 0 terminate. Both can be overridden through `typical`/`scale`. **Default: accept.**
7. **Features within the v2.0 modules** [§12.8]. `steps_view` as the only public manual-stepping API; defer `first_of_all`; drop `retry`, arclength continuation and `roots::sensitivity`; TOMS748, ITP, Halley and Steffensen in the optional phase 8; algorithm ids are an open enum. **Default: accept.**
8. **Run-time chains** [§12.19]. Offer `any_solver` and `first_of(range)` as an opt-in header, prototyped on all nine configurations. An empty chain is accepted and fails in-band with `invalid_input` when called. **Default: yes, opt-in header.**
9. **Multiprecision and complex** [§12.13, §12.14]. Keep multiprecision as an optional adapter and CI leg, and support complex numbers only in `poly`. **Default: keep both as described.**
10. **Where the combinators live** [§12.12]. Keep `first_of`/`then` in Numerixx until phase 3 has proven them, then upstream the generic parts to FXT. **Default: after phase 3.**
11. **Compiler floor** [§12.15]. GCC 14, Clang 19 + libc++, MSVC with `/std:c++latest` (19.51 tested), clang-cl, em++ ≥ 6.0.8. Deducing `this` sets GCC 14 (FXT uses it too, although its README states GCC 13 and Clang 17). Clang was raised from 18 to 19 on 2026-10-01, because Clang 18 rejects the refined literals. The nightly floor-compiler job checks the floor; all five of its legs passed on the spike branch on 2026-10-02. **Default: accept.**

---

## 11. Top risks

The full table is in §11 of the design.

| Risk | Mitigation |
|---|---|
| Scope creep towards "all of GSL" | the scope table in section 2 is the contract; each post-v2.0 family is its own module with its own corpus, and waits for `v2.0.0`; phase 8 may follow `v2.0.0`; candidates stay unscheduled; "out of scope" areas point to the standard library, Eigen or Boost.Math |
| GPL contamination (GSL, MPSolve and other GPL projects) | no GPL code is ported, paraphrased or copied; algorithms come from the literature, with references cited per header; oracles are Boost.Math, Eigen and reference tables (values published in the documentation of GPL projects may serve as reference facts) |
| FXT-1 delayed | only `numerixx::pipes` includes FXT; pin a patched fork commit until it lands |
| Learning curve (`.on`, `then`, criteria typed by view) | the one-call facade and the canonical calls are the front door; compile-fail tests check that diagnostics carry reasons |
| Template-heavy code on MSVC and clang-cl (diagnostics, P2564, mangling) | no clang-cl mangling issue found with these shapes; CI legs from day one; per-compiler deletion macro; forwarding rule for refined types |
| Eigen compile time (+2.5–6.4 s per TU that instantiates N-D solvers) and expression-template dangling under `auto` | the layering test confines Eigen to `linalg`/`multiroots` (and the multiprecision-linalg adapter); the umbrella header includes none of them; `NUMERIXX_WITH_LINALG=OFF` for scalar-only users; the facade returns concrete types only, and library code never binds `auto` to an Eigen expression; per-TU compile time is recorded per release (the 2 s CI guard covers only the umbrella header) |
| The facade misreports singularity (Eigen's partial-pivot LU never flags it, and rcond is not scale-invariant) | rcond compared against n·ε, plus a finiteness check of the solution → `errc::singular`; equilibrate before the rcond test if the corpus demands it; `qr_solve` for rank-deficient and least-squares systems; a corpus with singular and 1e10-column-scaled Jacobians; Eigen's `HybridNonLinearSolver` as oracle |

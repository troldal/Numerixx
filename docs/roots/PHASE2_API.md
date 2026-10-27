# Phase 2 API proposal: open methods and derivatives (Newton, Halley, Secant)

Status: approved by the user on 2026-10-27, P1 to P7 as recommended, with one change: the entry point is renamed
`solve<M>(...)` (it was `root<M>(...)` in Phase 1 and in the brief), so that the function and its result type `Root<T>`
are not told apart only by case. One review round (numerics-reviewer, simplicity-critic, usage-sketcher) is folded in
("Decisions from the reviews"). Architecture: `docs/roots/BRIEF.md` (approved 2026-10-25), sections 3.1 to
3.4; Phase 1: `PHASE1_API.md`.

## Summary

```cpp
namespace rt = nxx::roots;

// fn returns the value and its derivatives together: a std::pair, std::tuple or std::array of exactly T.
auto const fdf = [](double x) { return std::pair{ x * x - 2.0, 2.0 * x }; };
auto const a = rt::solve<rt::Newton>(fdf, 1.0);                    // from a guess: 7 evaluations
auto const b = rt::solve<rt::Newton>(fdf, rt::Bracket{ 1.0, 2.0 }); // safeguarded by the bracket: 8
auto const c = rt::solve<rt::Halley>([](double x) { return std::tuple{ x * x - 2.0, 2.0 * x, 2.0 }; }, 1.0); // 5
auto const d = rt::solve<rt::Brent>(fdf, rt::Bracket{ 1.0, 2.0 });  // a bracketing method uses the value only

// Or fn returns the value only, and each derivative is a callable of its own.
auto const f = [](double x) { return std::exp(x) - 2.0; };
auto const e = rt::solve<rt::Newton>(f, [](double x) { return std::exp(x); }, 0.0);
auto const g = rt::solve<rt::Secant>(f, 0.0);                     // no derivative: 9 evaluations

auto const h = rt::solve<rt::Newton>([](double x) { return std::pair{ (x * x - 2.0) * x + 2.0, 3.0 * x * x - 2.0 }; },
                                    0.0); // error::not_converged: Newton cycles between 0 and 1
```

- Three methods: `Newton`, `Halley` (from a guess, or safeguarded by a bracket) and `Secant` (from a guess).
- Derivatives in either form of BRIEF.md 3.4 (user decision): from fn as a tuple-like, or as separate callables.
- **Every result has a verified sign change** (P1, the main decision of this phase): an open method iterates from its
  guess until a sign change appears between two iterates, then finishes inside that bracket. `abs_error` means the
  same for every method: a bound on the distance to a sign change of fn as evaluated. This supersedes BRIEF.md 3.5
  ("open: the last step, an estimate") and 3.7's stop test for the open methods.
- **The safeguarded methods are Brent's loop** (Phase 1's `detail::brent`) with a Newton or Halley step in place of the
  interpolation; the open methods hand over to that loop, and Secant to Brent itself.
- A new error, `not_converged`.

## Public API (`roots.hpp`)

```cpp
namespace nxx::roots {

    enum class error {
        non_finite,     ///< A bracket end, the guess, a value of fn or a tolerance is not finite, or a tolerance is
                        ///< negative. Derivatives are not checked here (see the methods).
        no_sign_change,
        not_converged,  ///< An open method found no root from its guess: no sign change within its iteration cap, or
                        ///< no usable step (f' = 0 or not finite, Secant's two values equal, a step that would
                        ///< overflow).
    };

    struct Newton {}; ///< Newton's method: f and f'. From a guess, or safeguarded by a bracket.
    struct Halley {}; ///< Halley's method: f, f' and f''. From a guess, or safeguarded by a bracket.
    struct Secant {}; ///< The secant method: f only. From a guess.

    // 1. A bracket: Phase 1's overload, now also for Newton and Halley, and for fn returning a tuple-like.
    template<typename M, real T, typename F>
        requires detail::bracketing<M> && detail::gives<F, T, detail::values_needed<M>>
    [[nodiscard]] constexpr std::expected<Root<T>, error> solve(F &&fn, Bracket<T> const &bracket,
                                                               Criteria<std::type_identity_t<T>> const &criteria = {});

    // 2. A guess: Newton, Halley, Secant.
    template<typename M, real T, typename F>
        requires detail::open<M> && detail::gives<F, T, detail::values_needed<M>>
    [[nodiscard]] constexpr std::expected<Root<T>, error> solve(F &&fn, T const &guess,
                                                               Criteria<std::type_identity_t<T>> const &criteria = {});

    // 3 to 6. Separate derivatives: Newton with dfn, Halley with dfn and d2fn, each on a bracket or from a guess.
    // One-line forwarders to 1 and 2 through detail::Joined.
    template<typename M, real T, typename F, typename D>
        requires std::same_as<M, Newton> && nxx::detail::callable<F, T> && nxx::detail::callable<D, T>
    [[nodiscard]] constexpr std::expected<Root<T>, error> solve(F &&fn, D &&dfn, Bracket<T> const &bracket,
                                                               Criteria<std::type_identity_t<T>> const &criteria = {});
    template<typename M, real T, typename F, typename D>
        requires std::same_as<M, Newton> && nxx::detail::callable<F, T> && nxx::detail::callable<D, T>
    [[nodiscard]] constexpr std::expected<Root<T>, error> solve(F &&fn, D &&dfn, T const &guess,
                                                               Criteria<std::type_identity_t<T>> const &criteria = {});
    template<typename M, real T, typename F, typename D, typename D2>
        requires std::same_as<M, Halley> && nxx::detail::callable<F, T> && nxx::detail::callable<D, T>
                 && nxx::detail::callable<D2, T>
    [[nodiscard]] constexpr std::expected<Root<T>, error> solve(F &&fn, D &&dfn, D2 &&d2fn, Bracket<T> const &bracket,
                                                               Criteria<std::type_identity_t<T>> const &criteria = {});
    template<typename M, real T, typename F, typename D, typename D2>
        requires std::same_as<M, Halley> && nxx::detail::callable<F, T> && nxx::detail::callable<D, T>
                 && nxx::detail::callable<D2, T>
    [[nodiscard]] constexpr std::expected<Root<T>, error> solve(F &&fn, D &&dfn, D2 &&d2fn, T const &guess,
                                                               Criteria<std::type_identity_t<T>> const &criteria = {});

} // namespace nxx::roots
```

- **`detail::gives<F, T, K>`**: fn gives at least K values. Either fn's result converts to T (K = 1, Phase 1's
  `nxx::detail::callable`), or fn returns by value a `std::pair`, `std::tuple` or `std::array` with at least K
  elements, **every one of them exactly T**: `std::same_as<std::tuple_element_t<I, std::remove_cvref_t<R>>, T>`, which
  also rejects references (`std::tie` of fn's locals would dangle) and expression types. `values_needed` is 1 for
  `Bisection`, `Brent` and `Secant`, 2 for `Newton`, 3 for `Halley`. A tuple-like result is never read as a value.
- **`detail::bracketing`** gains `Newton` and `Halley`; **`detail::open`** is `Newton`, `Halley`, `Secant`.
- **fn's result is converted once**, where it is evaluated: `detail::evaluate<K>(fn, x)` returns `std::array<T, K>`
  (K = 1 for the value-only path), and the bracket's ends carry it, so the derivatives at the ends survive validation
  and the hand-over without a second evaluation. No method sees a pair or a tuple.
- **The separate callables** are value-only (`nxx::detail::callable`); a mixed form is not accepted. They are joined
  into one fn by `detail::Joined<F, Ds...>`, one variadic struct holding pointers (as `nxx::deriv`'s `Mirrored`, for
  MSVC's constant evaluation). Phase 3 lets them return a `std::expected` (`derivative_of`).
- **Overload resolution is unambiguous:** a callable is never a `Bracket` nor an `nxx::real`; a braced criteria list
  cannot deduce the T of a guess; a tuple-returning fn fails `callable` in overloads 3 to 6.
- **Misuse reads as an unsatisfied constraint** (`detail::gives<F, double, 2>`, `detail::open<Bisection>`), not a
  custom message (P7); the documentation of `F` names `derivative_of`, `Secant` and `Brent`.

## Behaviour

### A. Validation

As Phase 1 for the tolerances, a bracket and zero first. Only fn's value is checked: not finite at an end, at 0 or at
an iterate is `non_finite` (as in Phase 1). A guess not finite is `non_finite`. A derivative that is 0 or not finite
is never an error by itself: it makes the derivative step unusable (B: a bisection step; C: `not_converged`), so
f'(0) = inf at an end or at zero first does not fail a solve that Brent would finish (sqrt x - 0.5 on [0, 1];
cbrt x - 0.5 on [-1, 1]).

### B. The safeguarded methods: Brent's loop with a derivative step

Phase 1's `detail::brent`, unchanged in its bookkeeping (the swap that keeps |fn(b)| <= |fn(c)|, the reset of d and e
when c changes, the stop test on `distance_bound(c, b)`, the guard that bisects while the bracket is wider than
max / 16, the minimum step tol, `abs_error` = the bound on |c - b|), with the step chosen by the method:
- `Brent`: the inverse quadratic or linear interpolation, as now.
- `Newton`: `s = -f(b) / f'(b)`; `Halley`: Halley's step (C), from b, with the derivatives at b.
- **Acceptance is Brent's test written for a step s:** s points towards c, `|s| < 3 |m| / 2 - tol / 2` (well inside the
  bracket, so b + s cannot round onto c) and `|s| < |e| / 2` (less than half the step before last); and, as in Brent,
  a bisection step is forced when `|e| < tol` or `|f(a)| <= |f(b)|`. An unusable step (A) is a bisection step.
- **Termination:** Brent's argument holds unchanged: nested brackets, every step at least tol, bisection at least every
  2 log2(|e| / tol) steps; at most (k + 1)^2 - 2 evaluations for k bisection steps. No cap.
- **One loop for three methods** (and the hand-over is a call to it): about 70 lines that are not written twice, and
  one termination argument.
- **Source:** Brent 1973, ch. 4, with the textbook Newton and Halley steps. Boost.Math's `newton_raphson_iterate` and
  `halley_iterate` (BSL-1.0) also safeguard by a bracket and are a licensed reference, not copied. Not from Numerical
  Recipes' `rtsafe`, whose step is from the latest point; this one is from the better end, as Brent's.

### C. The open methods (`Newton`, `Halley`, `Secant` from a guess)

- **The step**, from the latest iterate:
  - `Newton`: `-f / f'`.
  - `Halley`: `-t / (1 - u)` with `t = f / f'` and `u = t (f'' / f') / 2`, the form without f'^2. Where u cannot be
    formed, `1 - u` is 0, or the Halley step points against Newton's (u > 1, which happens far from a root: x^3 - 2x + 2
    at 0.8 has u = 342), it takes Newton's step instead, as Boost.Math's Halley does (P5).
  - `Secant`: `(x - xp) (f / (fp - f))` through the latest two points: f and fp have one sign before the hand-over, so
    `fp - f` cannot overflow, and the ratio first keeps a large f from overflowing the product (Secant on
    1e300 (x - 1) from 1e7 converges in 4 evaluations; `f (x - xp)` first would overflow).
- **The minimum step** is tol(x), in the direction of the method's step taken from the signs (Newton, Halley: of f and
  f'; Secant: of f and the secant's slope), not from the sign of a step that may have underflowed to 0.
- **The hand-over (P1):** as soon as fn's values at two consecutive iterates have opposite signs, the method continues
  as the bracket loop of B on those two iterates (Secant as `Brent`), with their values already computed, without zero
  first (they are adjacent iterates, not the caller's bracket). A converging method hands over when an iterate crosses
  the root: early for x^2 - 2 from 1 (1.5 is past the root), late for x^2 - 2 from 10 (monotone from above, it crosses
  only at its minimum step). The cost is about one evaluation more than the textbook stop (7 against 6 for Newton on
  x^2 - 2 from 1).
- **`not_converged`:** no sign change within the cap, or no usable step (A), or an x + s that would overflow.
- **The cap is 4 x the bits of T** (`digits` times the smallest b with 2^b >= `radix`, an integer: 212 for double, 672
  for `cpp_bin_float_50`; P3). What it does not reach, documented with the remedy (a bracket, or `abs_tol`):
  - a multiple root at 0 from a guess with the default `abs_tol`: Newton on x^3 from 1 never changes sign and halves
    by 2/3 until x^3 underflows (612 steps); a bracket around 0 costs 3 evaluations (zero first);
  - a guess far from the root: Newton on x^2 - 2 from 1e100 halves its iterate about 330 times;
  - multiple roots beyond multiplicity about 5 (Newton), 3 (Secant), 11 (Halley).
- **Secant's second point** is `x0 (1 + eta)` with `eta = root_epsilon<4, T>()` (1.2e-4 in double; SciPy uses 1e-4),
  `eta` itself at x0 = 0, and `x0 (1 - eta)` where `x0 (1 + eta)` would overflow (P4). Relative, so no scale is assumed
  (the first proposal's `x0 + eta max(|x0|, 1)` left ln(x / 1e-12) from 1e-13 for 1.2e-4, seven decades away).

### D. The result

`Root<T>{value, abs_error}` as Phase 1: for every method a bound on the distance to a sign change of fn as evaluated,
at most 2 tol(value), or tol(value) where fn is exactly 0. What the documentation says besides:
- **An open method returns a root, not the root nearest the guess:** the first sign change its iterates cross (a cubic
  equation of state from a vapour-like volume may land on the liquid root, or across a pole).
- A pole is still returned as a root (`PHASE1_API.md` E): check |fn(value)| where fn may have poles; an open method
  can hand over to a bracket around a pole.
- From a guess that is already a root to within tol, the method takes one minimum step past it and stops: 2
  evaluations (x^2 - 2 from fl(sqrt 2)).
- The conditioning caveat (fn's rounding over |f'|) is as Phase 1's: on Wilkinson's W10 near 7 every method's bound
  is about 1e-15, and the computed fn's sign changes lie up to 2e-10 from 7.

### E. Overflow and constant evaluation

Every quotient, product and sum of the steps is checked before it is formed (GCC and Clang reject an overflow in
constant evaluation), by local helpers:
- a / b is formed only if `|a| / max < |b| / 2` (`nxx::deriv`'s reviewed form, conservative by a factor 2: a step
  above max / 2 is not formed; Newton on x - 1 from 1e308 is `not_converged`, a bracket solves it);
- a * b only if `|a| <= 1` or `|b| <= max / 2 / |a|` (Halley's u, Secant's product);
- x + s only if x and s have opposite signs, or `|x|, |s| <= max / 2`, or `|x| > max / 2` and `|s| <= max - |x|`
  (exact there, by Sterbenz's lemma).
In B no sum is checked: an accepted step lies inside the bracket, and the guard keeps every product of the
interpolation below max (Phase 1). `nxx::deriv`'s exact `sum_fits` and `product_fits` are not moved: these decisions
may be conservative.

## Prototype results

`.scratch/roots2/proto2.py`: Python, double, the behaviour of B to E on Phase 1's Brent (`proto3.py`); references from
mpmath. Evaluations of fn (one call returns every value the method needs). Every root is within `abs_error` of the
reference, except W10's (D) and the pole's (a pole).

| Case | Start | Newton | Halley | Secant | Brent |
|---|---|---|---|---|---|
| x^2 - 2 | x0 = 1 / [1, 2] | 7 / 8 | 5 / 6 | 9 | 9 |
| x^2 - 2 | x0 = 10 / 1e-3 / fl(sqrt 2) | 9 / 17 / 2 | 7 / 12 / 2 | 13 / 25 / 4 | |
| x^3 - 2x - 5 | x0 = 2 / 0 / 10 / [2, 3] | 6 / 21 / 11 / 7 | 5 / 12 / 8 / 6 | 8 / 15 / 15 | 8 |
| (x - 1)^3 | x0 = 2 / [0, 3] | 89 / 91 | 53 / 54 | 127 | 156 |
| x^20 - 1 | x0 = 0.5 / 2 / [0, 5] | 22 / 20 / 11 | 11 / 11 / 11 | 44 / 28 | 19 |
| x^3 - 2x + 2 | x0 = 0 (cycle) / -2 / [-2, -1] | `not_converged` / 7 / 8 | `not_converged` / 5 / 6 | 28 / 9 | 9 |
| 1/x - 1 | x0 = 0.5 / 3 / [0.5, 3] | 7 / `not_converged` / 9 | 2 / 3 / 6 | 10 / `not_converged` | 11 |
| x^3 (root at 0) | x0 = 1 / [-1, 2] | `not_converged` / 3 | `not_converged` / 3 | `not_converged` | 3 |
| 1e300 (x - 1) | x0 = 1e7 | 2 | 2 | 4 | |
| cos x - x | x0 = 0 / 1 / 3 / [0, 1] | 6 / 5 / 7 / 6 | 5 / 4 / 6 / 5 | 10 / 7 / 10 | 9 |
| e^x - 2 | x0 = 0 / 5 / [0, 1] | 7 / 10 / 7 | 4 / 6 / 5 | 9 / 15 | 8 |
| x e^x - 1 | x0 = 0 / 2 / [0, 1] | 8 / 9 / 8 | 5 / 5 / 6 | 10 / 12 | 9 |
| atan x (hand-over across 0) | x0 = 1 / 1.5 / [-1, 3] | 6 / 6 / 3 | 5 / 5 / 3 | 9 / 10 | 3 |
| Kepler, e = 0.9, M = 1 | x0 = 1 / 0 / pi / [0, pi] | 7 / 9 / 7 / 8 | 5 / 7 / 5 / 6 | 10 / 13 / 9 | 10 |
| ln(x / 1e-12) | x0 = 1e-11 / 1e-13 / [1e-20, 1] | `non_finite` / 8 / 20 | `non_finite` / 6 / 18 | `non_finite` / 11 | 45 |
| sqrt x - 0.5, f'(0) = inf | [0, 1] | 4 | 4 | | 4 |
| W10 (Wilkinson) near 7 | x0 = 7.2 / 6.8 / [6.7, 7.6] | 21 / 22 / 28 | 19 / 16 / 24 | 23 / 23 | 23 |

`non_finite` for ln(x / 1e-12) from 1e-11: every method's first step leaves fn's domain (x < 0), where fn is NaN.

## Tests (`tests/test_roots_open.cpp`, test cases `roots:...`)

- References as Phase 1 (closed forms, published constants, mpmath; 70 digits for the multiprecision types).
- **Accuracy:** `|value - ref| <= abs_error + c`, c derived per case as in Phase 1. For W10, whose c from a Horner
  error bound (about 7e-8) is too loose to catch a regression, the sign change itself is checked: fn at
  `value - abs_error` and `value + abs_error` has opposite signs or is 0.
- **Both derivative forms give the same result, bit for bit**, for the same fn; a tuple fn on `Bisection` and `Brent`
  gives the same result as the value-only fn (P6).
- **Counts:** exact for fn of `+ - * /` in double against `proto2.py`, assuming no FMA contraction in fn or the steps
  (true on x64 and wasm today; noted with the counts). The cycle's count follows from the cap: 1 + 212 evaluations.
- **Failures:** the cycle, the divergence (1/x - 1 from 3), x^3 from 1, Newton from 1e308, f' = 0 at the guess,
  Secant's flat start, a guess not finite, fn not finite at an iterate (ln from 1e-11); each the documented error.
- **Edges:** f' infinite at a bracket end and at 0 (sqrt, cbrt): bisection steps; a guess exactly at a root; the early
  and the late hand-over (x^2 - 2 from 1 and from 10); a hand-over across 0 (atan x from 1); a hand-over onto a pole;
  the Secant overflow case (1e300 (x - 1)).
- **Constant evaluation:** `static_assert`s for each overflow check of E (Newton from 1e308; Secant on 1e300 (x - 1);
  Halley with |f'| near max), and results of the three methods for float, double and long double.
- **Constraints**, as requires-expressions: `solve<Newton>` on a value-only fn, `solve<Bisection>` from a guess,
  `solve<Secant>` on a bracket, `solve<Halley>` on a pair, `solve<Newton>(f, dfn, d2fn, x0)`, a tuple of references, a
  tuple of expression templates, a tuple fn with a separate dfn: each rejected.
- **Types:** float, long double, `cpp_bin_float_50` with and without expression templates, `cpp_dec_float_50`: a few
  accuracy cases per method.
- The formatter writes `not_converged`.

## Documentation

- `roots.hpp`: the three methods, both derivative forms, the hand-over and its cost, "a root, not the nearest root",
  `not_converged` and the cap's limits with their remedies.
- `docs/roots/design.md`: this phase's decisions and the prototype's table; that P1 supersedes BRIEF.md 3.5 and 3.7 for
  the open methods. README: one Newton example.
- `demo/roots_demo.cpp`: a section on the open methods and on the two derivative forms.
- `docs/plan.md`: Phase 2's checkboxes.

## Decisions by the user (2026-10-27)

P1 to P7 approved as recommended. The entry point renamed from `root` to `solve` for every method, Phase 1's included
(breaking for Phase 1's callers; the library is unreleased). The options as they were put:

- **P1. The hand-over** (C): every open method finishes inside the first sign change its iterates cross, so every
  result has a bound. Recommended. Alternative: the textbook stop, |step| <= tol, with `abs_error` the last step as an
  estimate (BRIEF.md 3.5): one evaluation cheaper, but it accepts false roots (Secant on x^20 - 1 from 0.5 stopped at
  a point where fn was -1, after a step of 1e-84) and never stops on an ill-conditioned fn (Newton on W10 reached the
  cap); a stagnation rule was tried and is unsound (it stopped 3.7 from the root of x^3 - 2x - 5).
- **P2. One `not_converged`** for every failure of an open method. Recommended: the caller's action is the same,
  another guess or a bracket; the usage sketch would like f' = 0 told apart from the cap for diagnosis, which a later
  value can do (the enumeration is documented as extensible). Alternative: `zero_derivative` now.
- **P3. The cap: 4 x the bits of T, internal**, with its limits documented (C). Recommended. Alternatives: a larger
  multiple; `Criteria::max_iterations` (non-breaking later, appended). BRIEF.md 3.7 had a cap from the exponent range,
  which in `cpp_bin_float_50` would be about 2^31 iterations.
- **P4. Secant's second point** `x0 (1 + root_epsilon<4, T>())`, relative. Recommended. Alternative: a two-point start
  (deferred, as its own aggregate; see Deferred).
- **P5. Halley's step** in the form `-t / (1 - u)`, falling back on Newton's where it cannot be formed or points
  against it. Recommended. The fallback costs one case of the study: from 0 on x^3 - 2x + 2, raw Halley escapes
  Newton's cycle and converges in 9 evaluations, the fallback follows the cycle to `not_converged`.
- **P6. Brent and Bisection accept a tuple fn** and use its first element (BRIEF.md 3.4). Recommended.
- **P7. Constraint failures rather than custom messages** for misuse. Recommended; a `static_assert` catch-all would
  make the call look valid to a requires-expression.

## Decisions from the reviews

From the numerics review:
- A and B contradicted each other on a derivative that is not finite (and A's rule failed sqrt x - 0.5 on [0, 1]):
  only fn's value is an error; a derivative that is 0 or not finite makes the step unusable (A).
- The overflow checks were incomplete (Halley's 2 f' and t f'', Secant's f (x - xp)) and the sum rule both too strict
  and too loose: every check is stated, Secant's step reordered, the sum rule exact above max / 2 (E).
- The cap's claims missed a multiple root at 0, far guesses and Secant's and Halley's multiplicity limits: documented
  (C), and the change from BRIEF.md 3.7 called out (P3).
- Halley's step can point uphill (u > 1): Newton's step instead, as Boost.Math does (C, P5).
- B omitted Brent's forced bisection when |e| < tol and the reset of e when c changes, and an accepted step could round
  onto c: B is Brent's loop itself, with Brent's acceptance margin.
- The prototype's Secant hand-over ran zero first, its Newton hand-over did not: no zero first on a hand-over (C).
- The minimum step's direction came from the sign of a step that may underflow: from the signs of f and f' (C).
- Secant's second point assumed a scale of 1: relative (C, P4).
- The claim that every case ends in a bound or `not_converged` was false (ln from 1e-11 is `non_finite`), and the cost
  sentence contradicted its own example: both corrected, and the table now reports every case run.
- Tests: W10 on [6.5, 7.5] measured nothing (its first midpoint is exactly 7, a zero of the Horner evaluation): an
  asymmetric bracket, and the sign change checked directly; tests added for the overflow checks, f' infinite, a guess
  at a root, the hand-overs, the cycle's count.

From the simplicity review:
- **One bracket loop** for Brent, safeguarded Newton and safeguarded Halley, the step chosen by the method; the
  hand-over is a call to it.
- fn's result converted once (`evaluate<K>` returning `std::array`); one step function per method shared by B and C;
  one open loop for Newton and Halley; `Joined` variadic; overloads 3 to 6 as forwarders.
- The tuple rule rejected `std::tie` of locals only by accident of wording: every element must be exactly T, the
  element type itself (`tuple_element_t` of the cv-unqualified result).
- An overflowing x + s is `not_converged` (every value of fn was finite), and the enumerator's documentation says so.
- No x + s check in B; no named product predicate beyond the two uses; the deferred "safeguarded Secant" removed (it is
  Brent).
- Kept: the six overloads, `gives`, `open`, P1, P2, P3, the local overflow checks, P7.

From the usage sketch (7 snippets, 7 misuses):
- An open method returns a root, not the nearest one (a cubic equation of state can land on the other phase's root):
  documented (D).
- Polishing a Brent result with Newton costs 2 evaluations from a point already within tol: documented (D).
- `not_converged` cannot be diagnosed further: kept as one value (P2), with the causes in its documentation.
- A multiprecision tuple fn needs `std::pair<Big, Big>{...}`; the generic constraint message does not name the
  offending element: documented with an example.

## Deferred

- Derivative callables that return `std::expected`, and `nxx::deriv::derivative_of`: Phase 3.
- Mixed derivative forms (a tuple fn with a separate d2fn; a dfn returning f' and f''): new overloads, non-breaking.
- User tuple-like types: only by an explicit opt-in trait, never detected through `std::tuple_size`, which would change
  the meaning of a user type that converts to T and has a tuple protocol.
- New starts (a two-point start for Secant; a guess inside a bracket; fn's values supplied with a bracket): each as its
  own aggregate type, never as extra positional T parameters, which would capture a braced criteria list such as
  `{1.0}` (an identity conversion beats aggregate initialisation) and change the meaning of existing calls.
- `Criteria::max_iterations`; a separate `zero_derivative` error; moving the overflow predicates into `math.hpp` when a
  third module needs them; skipping the derivatives where the step will be a bisection (Phase 3's costly numerical
  derivatives).

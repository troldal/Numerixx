# `nxx::quadrature`: architecture brief

Status: approved by the user on 2026-10-30, D1 to D14 as recommended (draft 2). Not a phase API: it fixes the
architecture that the phase specifications (`PHASE<n>_API.md`) will refine. Draft 1 went through one review round
(`usage-sketcher`, `simplicity-critic`, `numerics-reviewer`); what it changed is in section 9, the decisions in
section 8. The module is `nxx::quadrature` in `quadrature.hpp` (D9); the roadmap's working name was `nxx::integrate`.

## 1. Scope

The definite integral of a real function of one real variable, f: T -> T for any `nxx::real` T, over a finite,
half-infinite or infinite interval, with an estimate of the error. Out of scope here: integrals in several dimensions,
oscillatory integrands of a known frequency (Filon, Levin, Ooura-Mori), Cauchy principal values, weight functions
(Gauss-Jacobi, -Laguerre, -Hermite), integration of tabulated data (the interpolation module), and complex-valued fn.

## 2. Prior work

### 2.1 Numerixx 1.x (survey of `v1.1.0-legacy`)

`numerixx/integrate/` (about 1,050 lines) and `demo/DemoIntegrate.cpp`; no tests.

- **Entry point.** `integrate<Solver>(fn, bounds, tol = eps, maxIter = 25)`, `Solver` a template-template
  parameter on a CRTP base; `Trapezoid`, `Simpson` and `Romberg`, all global (the whole grid halves), none adaptive.
  `integralOf<Algo>(fn)` returns a functor that ignores fn and throws.
- **Result.** `std::expected<T, Error>` with the bare value: no error estimate on success.
- **Termination.** `|I_k - I_{k-1}| < tol`, absolute, with `tol = epsilon` by default and no minimum level count. So
  sin^2 x on [0, 2 pi] returns 0 (the samples at 0, pi and 2 pi all vanish, T0 = T1 = 0), and any integral of
  magnitude above about 1 runs all 25 levels (2^25 evaluations) and fails.
- **Defects.** Validation throws (`lower >= upper`); nodes computed from `double` parameters (user types truncated;
  the multiprecision demo is commented out); Trapezoid and Simpson store every node (hundreds of MB at level 25);
  Simpson re-evaluates every point at every level; uninitialised reads at `maxIter <= 0`; unqualified `abs`;
  `boost::multi_array` for two rows of the Romberg table. Finite intervals only; an integrable singularity at an end
  fails ("initial estimate is not finite").
- **Worth keeping:** the `std::expected` return; the trapezoid halving that reuses the previous sum; the Romberg
  recurrence, which is correct as written. None of the code is carried over unchanged.

### 2.2 Boost.Math 1.92 (`boost/math/quadrature/`, BSL-1.0, code may be adapted with its notice)

- **Methods:** `gauss<Real, N>` (N = 7, 10, 15, 20, 25, 30; no error estimate), `gauss_kronrod<Real, N>` (N = 15 ...
  61; adaptive by recursive depth-first bisection, error `max(|K - G|, 2 eps |K|)`), `tanh_sinh`, `exp_sinh`,
  `sinh_sinh`, `trapezoidal`. Default tolerance sqrt(eps) everywhere; the documentation advises against asking for
  eps.
- **Nodes:** 36-digit `constexpr std::array` tables cast to the smallest built-in type that holds them, 115-digit
  strings for types up to 110 digits, and above that nodes computed on demand and cached in function-local statics.
  The double-exponential rules hard-code 8 levels for the built-in types and compute the rest with `tanh`, `sinh`,
  `cosh`, `exp`, `log` and `asinh`, cached lazily behind a `std::mutex` in a `shared_ptr` the copies share.
- **tanh-sinh's endpoint handling is the part to learn from:** the abscissas near +-1 are stored as their complements
  (the distance to the end, formed directly, never as 1 - x); fn may take `(x, xc)`, the distance to the nearer end;
  a per-end minimum distance in ulps keeps fn from ever being called at an end; points where fn is not finite are
  dropped from the tails and negligible tail terms pruned; the two extreme terms are added to the error estimate
  (`endpoint_error`), and growth at a narrowed end is an error; at least 4 levels are run; the test is
  `|I_k - I_{k-1}| <= tol L1`, against the integral of |f|, not |I|.
- **Errors:** policies, which throw by default; Gauss, Gauss-Kronrod and trapezoidal propagate NaN silently; on
  non-convergence the best estimate is returned with no signal.
- **Not usable as is:** the mutex (Emscripten), the throwing policies, the statics and vectors (constexpr), the
  string tables (construction from a string is not in `nxx::real`), the transcendental functions (`math.hpp` has
  none). Usable: the endpoint logic of tanh-sinh, restructured.
- **A defect, confirmed by reading:** the adaptive Gauss-Kronrod error is computed on [-1, 1] and never scaled to the
  subinterval, while the tolerance it is compared with and the L1 norm are (`gauss_kronrod.hpp:1201, 1209, 1222,
  1224`). Its (-inf, inf) map forms 1 - t^2 rather than (1 - t)(1 + t) (1243-1244).

### 2.3 The literature, and what it settles

The verified sources are listed in section 7. What they settle for the design:
- **Global adaptivity beats local.** Bisecting the subinterval with the largest error estimate (a priority queue)
  controls the total error better than giving each subinterval a share of the tolerance (Malcolm and Simpson 1975;
  QUADPACK's QAG; Shampine's `quadgk`, 2008). Boost's adaptive Gauss-Kronrod is local.
- **No error estimate is a bound.** Every estimate from finitely many samples is zero on an infinite-dimensional
  set of integrands, a narrow peak between the nodes among them (Gonnet 2012, section 6; Lyness 1983). What can be
  claimed: each panel integrates polynomials up to its degree exactly (3.6), up to the rounding of the nodes and
  weights and of the sum; the estimate is never below the rounding floor; and it is a heuristic that overestimates
  the error on the published test batteries in nearly all cases.
- **|K - G| is the robust estimator.** It estimates the Gauss rule's error, not the Kronrod result's, so on smooth fn
  it overstates the error of the returned value by orders of magnitude; that costs evaluations, not correctness.
  QUADPACK's `resasc min(1, (200 |K - G| / resasc)^1.5)` tracks Kronrod's own error on analytic fn, but the power
  1.5 holds only there: on fn with relative noise delta, |K - G| is about delta L1 and the scaled form claims
  (200 delta)^1.5 of it, 2.8e-15 for delta = 1e-12, and 3e-21 for a double-precision fn called from a
  `cpp_bin_float_50` integral (numerics review). Neither estimator sees a peak missed by every node (K = G).
- **Tolerances.** A relative test alone never ends on an integral equal to 0: the test is `err <= max(abs_tol,
  rel_tol |I|)` with a floor proportional to eps times the integral of |f| (rounding), as QUADPACK and Gander and
  Gautschi (2000) do.
- **Double-exponential rules** (Takahasi and Mori 1974; Mori and Sugihara 2001) are the best general method at high
  precision and for endpoint singularities (Bailey, Jeyabalan and Li 2005: best of three on 15 integrals at 400
  digits). Bailey's error estimate is `max(d1^2 / d2, 2 d1, d3, d4)` in decimal digits: d1 and d2 the last two level
  differences (d1^2 / d2 extrapolates the observed rate), d3 eps times the largest term (rounding), d4 the largest
  term at the ends of the sum (truncation). Boost uses the plain difference plus its endpoint term.
- **Romberg** needs more evaluations than Gauss-Kronrod for the same accuracy and is uncompetitive at high precision;
  its uses are teaching and equally spaced data. The trapezoid halving underneath it is the machinery the
  double-exponential rules use too.
- **QUADPACK's Fortran** carries no licence statement on netlib, so it is excluded twice over (our rule admits only
  Numerixx 1.x and Boost.Math); the QUADPACK book and its formulas are used freely.

## 3. The architecture

### 3.1 Entry point

```cpp
namespace nxx::quadrature {
    template<typename M = GaussKronrod<15>, real T, typename F>
    [[nodiscard]] constexpr std::expected<Integral<T>, Failure<T>>
    integrate(F &&fn, T const &a, T const &b, Criteria<std::type_identity_t<T>> const &criteria = {});

    // The fixed rule: a number, no estimate, no criteria (D4).
    template<typename M, real T, typename F> requires /* M is GaussLegendre<N> */
    [[nodiscard]] constexpr std::expected<T, error> integrate(F &&fn, T const &a, T const &b);
}
```

- **One function name for every method and every kind of interval**, as `solve` is for roots: the method is an empty
  tag, the set closed (as in `nxx::deriv` and `nxx::roots`). An infinite end (where T has one) selects exp-sinh or
  sinh-sinh under `TanhSinh`; the Gauss methods and Romberg reject it (3.5).
- **The ends as two positional arguments, a then b,** oriented: b < a gives the negated integral, a == b gives 0
  without calling fn (after the NaN check; inf == inf is 0 too). Every integration library writes it this way
  (QUADPACK, Boost, SciPy, MATLAB), and unlike a bracket the order carries meaning, so an aggregate `Interval{a, b}`
  would add a type without removing an ambiguity.
- **Later arguments are never bare spans or T:** break points, when they come, are a named aggregate (`Points{...}`)
  or a separate function. A `std::span<T const>` overload would make `integrate(fn, a, b, {1e-12, 1e-10})` ambiguous
  with `Criteria` once C++26 builds spans from braced lists (P2447); roots' rule against new positional T applies.
- **`[[nodiscard]]` and the default method** were missing from draft 1 (usage sketch).

### 3.2 The methods

| Method | Interval | Rule | Error estimate | Primary source |
|---|---|---|---|---|
| `GaussLegendre<N>` | finite | fixed N-point Gauss rule, any N in a range fixed by Phase 1 | none (D4) | textbook; nodes as in 3.4 |
| `GaussKronrod<15>` | finite | G7-K15 panels, globally adaptive; N = 15 only (D7) | sum of per-panel \|K - G\| (3.6) | Kronrod 1965; Laurie 1997; Piessens et al. 1983; Malcolm and Simpson 1975; Shampine 2008 |
| `TanhSinh` | finite, half-infinite, infinite | double exponential: tanh-sinh on [a, b], exp-sinh on a half line, sinh-sinh on the line | relative level difference plus the end term (3.6) | Takahasi and Mori 1974; Tanaka et al. 2009; Bailey et al. 2005 |
| `Romberg` | finite | trapezoid halving with Richardson extrapolation | from the table, after a minimum number of levels | Romberg 1955; Bauer, Rutishauser and Stiefel 1963 |

- **`GaussKronrod<15>` is the default:** smooth fn to full double precision in a few hundred evaluations, interior
  kinks and mild endpoint singularities by subdivision; Shampine's `quadgk` and QUADPACK make the same choice. The
  template parameter is kept with one value so that K21 ... K61 can be added later without breaking the callers who
  wrote `GaussKronrod<15>` (a plain `GaussKronrod` turned into a template breaks them; the naming rule forbids
  `GaussKronrod61`). It is a hook against a breaking change, not an extension point: a second value comes only on
  evidence that `TanhSinh` does not already serve high precision (D7).
- **`GaussLegendre<N>` is the fixed rule**, no adaptivity, so its result is a smooth function of parameters inside
  fn: the rule to use inside a root finder, an optimiser or a nested integral, where an adaptive result jumps when
  the subdivision changes (Lyness 1983). N is a range (as `nxx::deriv`'s stencil orders), not Boost's list of six:
  the generator works for any N, so a list is an arbitrary restriction.
- **`TanhSinh` is the one double-exponential method** and picks tanh-sinh, exp-sinh or sinh-sinh from the ends:
  there is no choice for the caller to make, and it saves two tags and a "wrong interval" error. The three node
  formulas exist either way. The name is tanh-sinh's because that is how the family is known; `DoubleExponential` is
  the alternative (D9).
- **`Romberg`** is a roadmap requirement (teaching, and the counterpart of Ridders in `nxx::deriv`: Richardson
  extrapolation in both); the literature gives it no numerical case beside Gauss-Kronrod, and its documentation says
  so. It evaluates fn at the ends (a closed rule), so an endpoint singularity is `non_finite`; it runs a minimum of 4
  levels (1.x's sin^2 failure). It needs no nodes, so it comes in Phase 1 as the first method with an estimate.

### 3.3 The result and the failure

```cpp
template<real T>
struct Integral
{
    T value;     // the integral
    T abs_error; // an estimate of |value - exact|, not a bound, never below the rounding floor (3.6)
    friend constexpr bool operator==(Integral const &, Integral const &) = default;
};

template<real T>
struct Failure
{
    error reason;                     // why there is no result
    std::optional<Integral<T>> best;  // for not_converged: the estimate reached, with its error estimate
    friend constexpr bool operator==(Failure const &, Failure const &) = default;
};
```

- **`abs_error` is a plain T,** present for every method that returns an `Integral` (the roadmap: every result
  carries an error estimate). Draft 1's `std::optional`, empty for the fixed rule, made every adaptive caller write
  `*r->abs_error` to serve one method, whose callers then got undefined behaviour from the same expression and
  criteria that were silently ignored; the method is known at compile time, so the answer moved there: the fixed rule
  returns a plain T (3.1). D4.
- **An estimate, never a bound**, and the documentation says so per method; unlike `nxx::roots::Root`, whose
  `abs_error` is a bound (2.3). It is never below the rounding floor, so an integral that is exactly 0 reports the
  floor, not 0 (as roots reports the tolerance at an exact zero of fn).
- **The best estimate on failure (D3).** Every integration library returns it (QUADPACK with `ier`, SciPy and MATLAB
  with a warning), and the cases are real: a noisy fn (a nested integral or an iterative solve inside fn) has an
  error floor that bisection cannot remove; a tight request reaches the cap; the caller wants to know whether it
  missed by a factor of 2 or 1e6. Roots' rule (failure returns only the code) was made where the last iterate of a
  failed open method is worthless; here it is an integral with an error estimate. Draft 1 argued that a richer error
  type breaks composition with `and_then`, but no composition across modules works without `transform_error` anyway
  (`nxx::roots::error` and this module's differ), and within the module nothing chains. Changing the error type later
  breaks every caller, so it is decided now. `best` is empty for `non_finite`. `Failure` is formattable through its
  reason.
- **No conversion to T**, as `Root`. Read by member name, so that members can be added: the number of evaluations
  (to compare methods) and the integral of |f| (L1 / |I| is the condition number of the sum; it is computed for the
  floor anyway) are wanted by the usage sketch; both are non-breaking to add, so they wait for the phase that needs
  them (D12).
- **Its own types**, not `nxx::deriv::Estimate` (the modules stay independent; `Estimate` converts to T).

### 3.4 Nodes and weights for any T

A spike (2026-10-30, scratchpad, not kept) settled what can be done from `+ - * /` alone, on clang-cl, cl, GCC 16 and
Clang 22, for double, GCC's 80-bit long double, `cpp_bin_float_50` and `cpp_dec_float_50` (expression templates on).

- **One generator for every T, arithmetic only.** Gauss-Legendre: a start, then Newton on P_n through the three-term
  recurrence, weights `2 / ((1 - x)(1 + x) P_n'(x)^2)`. The spike started from Tricomi's estimate (pi and cos from
  series in T); the zeros of P_n interlace those of P_(n-1), so Newton safeguarded by bisection between them needs no
  pi and no cos (simplicity review). Phase 1 measures the interlacing route's constexpr steps (an estimate is about
  100k for G30, against 14k) and takes it if it fits. Kronrod: the Stieltjes polynomial E_{n+1} as a Legendre series
  whose coefficients come from a triangular system (the orthogonality conditions, with Adams' closed form for the
  integral of a product of three Legendre polynomials; no linear solve), its zeros by safeguarded Newton between the
  interlacing Gauss nodes, the weights from the interpolatory formulas. Verified for G7 ... G30 and K15 ... K61
  against Boost's 115-digit tables.
- **Accuracy.** Nodes within 1.14 eps of the correctly rounded value for every T. Weights near +-1 lose up to 330 eps
  relative (at a node d ln w / dx = -2x / (1 - x^2), so a node's rounding is amplified; T arithmetic cannot avoid it),
  but the absolute error of an outer weight is about 2.7 times the node's, whatever n. Measured in the integral, on
  f = 1, e^x, Runge and positive end-weighted polynomials of the highest exact degree: the plain tables err by up to
  4 eps times the integral of |f| (24 eps absolute on K61's degree-90 polynomial), correctly rounded tables by up to
  1.15 (K21, double), a factor of about 3.5; on f = 1 and e^x the two differ by up to 8 eps, on Runge by 1. Against
  adversarial signs the bound is the sum of |w_hat - w| times max |f|, about 19 eps for G30 by the numerics review's
  estimate; Phase 1 computes that sum per rule and per T against the 115-digit tables, and the exact-polynomial tests
  take their tolerances from it. For tanh-sinh the plain double table's 506 eps in the far tail is a consistent shift
  in t (the complement and the weight come from the same rounded exp), which the trapezoid sum tolerates: no test
  integral moved by more than 8 eps in a single term. Computing in a double-word type and rounding gives correctly
  rounded tables but does not fit the default constexpr limits beyond G10 and K21. Not used.
- **Compile time, within the default limits** (Clang, clang-cl and cl 1,048,576 steps, GCC 2^25 operations): G30 in
  14k steps, K61 in 78k, the 769-node tanh-sinh table for double (h = 2^-7, t <= 6) in 573k on Clang and 719k on cl
  (69 % of its limit), 0.3 to 0.6 s of compile time. C++23's `static constexpr` in a constexpr function (P2647) and a
  non-literal static on a path not taken at compile time (P2242) both work on cl and clang-cl.
- **Run time, for types that are not literal** (Boost's): all twelve Gauss and Kronrod rules in 15 to 23 ms in
  `cpp_bin_float_50`, about 1.5 ms per rule.

The design that follows:
- **float and double:** each rule a `static constexpr` table, computed at compile time by the generator, once per
  program; no run-time cost, no data in the header.
- **long double:** as double where it is double (MSVC); otherwise the tanh-sinh table needs t up to about 8.9 (where
  `e^u cosh u` and `cosh^2 u` overflow, which is where the table must end, not at a fixed t), so the compile-time
  budget is measured in Phase 3, and where it does not fit, long double is computed per call as below. Emscripten's
  long double is binary128, on a CI toolchain: measured there too.
- **Every other T:** the same generator at the start of each `integrate` call (once per call, not per panel). The
  usage sketch puts the cost in context: a double integral by iterated calls in `cpp_bin_float_50` (about 150 inner
  calls) spends about 225 ms on nodes, several times its evaluations of a cheap fn. A cache across calls (a
  function-local static initialised once by the language, never mutated; P2242) removes it. It changes no API, so it
  is non-breaking either way; D10 asks whether it comes with Phase 2 or waits.
- **The alternative, Boost's tables as data** (36-digit literals for the built-in types, BSL notice): correctly
  rounded and no compile-time cost, but a second code path beside the generator, which other T still need, and a few
  thousand lines of constants, for a factor of about 3.5 on end-weighted integrands. Not recommended.
- Open for the phases: radix 10 (`cpp_dec_float`; its guard digits make it more accurate than its eps says, and the
  spike's double-word sums for pi and ln 2 assume radix 2); exp's overflow and subnormal ranges; the exp-sinh and
  sinh-sinh tables' step counts (Phase 4); in double, t <= 6 cuts the tanh-sinh tail at x = 1e-275, not at the
  smallest normal number.

### 3.5 Infinite intervals

- An end equal to +-infinity (`std::numeric_limits<T>::has_infinity` and the value infinite) is accepted by
  `TanhSinh` only: exp-sinh on a half line, sinh-sinh on the line, each with its own node formula. The Gauss methods
  and Romberg return `non_finite` for it, as for a NaN end, and their documentation points to `TanhSinh`. The roadmap
  gives infinite intervals to exp-sinh and sinh-sinh; draft 1's Phase 3 (Gauss rules through a map) had no use beyond
  them (simplicity review), and the cost is that the default method rejects `integrate(fn, 0.0, inf)`. D13.
- If the Gauss maps come later (deferred, non-breaking: an error becomes a success), the numerics review fixed the
  map: QUADPACK's x = a + (1 - t) / t puts the infinite end at t = 0, where T resolves down to denorm_min, while
  t / (1 - t) and t / (1 - t^2) put it at t = 1, where the spacing is eps / 2 and x stays below about 2 / eps (x^-1.5
  on [1, inf) then loses 3e-8 relative, the default tolerance); the Jacobian is formed as (f / t) / t, so an f that has
  underflowed to 0 gives 0, not 0 times inf; (-inf, inf) folds f(x) + f(-x) onto one half line, as QAGI does.
- The exp-sinh shift to [a, inf) helps the precision near a only together with fn taking the distance to the end
  (deferred, 3.8); until then fn sees a + t.
- What no rule handles and the documentation says: fn that does not decay (sin x on [0, inf) has no integral), slow
  algebraic decay (1 / x^1.01), oscillation on an infinite interval, and `max()` or 1e10 written for infinity (every
  node misses the mass, and the result converges to about 0).

### 3.6 Error estimates and termination

```cpp
template<real T>
struct Criteria
{
    T abs_tol = T(0);                 // for integrals near 0
    T rel_tol = /* sqrt(eps) of T */; // root_epsilon<2, T>()
    friend constexpr bool operator==(Criteria const &, Criteria const &) = default;
};
```

- **The test:** stop when the summed estimate is at most `max(abs_tol, rel_tol |I|, floor)`. As SciPy's `quad`,
  MATLAB's `integral` and QUADPACK (an absolute and a relative tolerance, either sufficient); roots' `abs_tol +
  rel_tol |x|` is the same up to a factor 2, and the max form lets abs_tol = 0 mean "relative only".
- **The floor**, per panel: about 50 eps times L1_p, the integral of |f| over the panel by the same rule (h times the
  sum of w_K,i |f_i| on the same evaluations, QUADPACK's `resabs`). That is the scale of the sum's own rounding, so it
  is right for cancelling integrands and for I = 0, where Boost's `2 eps |K|` is 0. The constant is not cosmetic: the
  plain-T weights make |K - G| on f = 1 up to 9 eps on I = 2 in `cpp_bin_float_50` (the spike's output), a fixed
  fraction of L1 that bisection cannot reduce, so a floor of "a few eps" never ends on a constant at a tolerance near
  eps. QUADPACK's 50 covers it with room; Phase 2 derives the constant from the measured sum of |w_K - w_G| errors and
  the n_K eps of the sum, and forms K - G with difference weights d_i = w_K,i - w_G,i whose sum is forced to 0 through
  the centre weight, so that a constant gives K - G = 0 up to the summation.
- **The defaults:** `rel_tol = sqrt(eps)` (1.5e-8 in double, 1e-25 in `cpp_bin_float_50`; as Boost and SciPy),
  `abs_tol = 0`. Not full precision, as roots defaults to: a quadrature cannot certify it, and |K - G| makes the
  result on smooth fn far more accurate than asked. A caller who chose a 50-digit type for its digits sets rel_tol;
  the documentation's multiprecision example shows it. D5.
- **Gauss-Kronrod's estimator is |K - G| per panel** (2.3), summed. Draft 1 proposed measuring it against QUADPACK's
  scaled form in Phase 2; the noise counterexample (2.3) and the simplicity review's count (a second sum and a power
  1.5, so a sqrt, for constants tuned in double) settle it without the study. The study, if ever, is deferred, with
  the numerics review's list of what it must include. D6.
- **The adaptive loop** (global): an initial mesh of a few panels (Shampine's `quadgk` starts with 10, so that one
  panel cannot accept K = G on a missed peak; Phase 2 fixes the count), then the panel with the largest estimate is
  bisected, at h = b/2 - a/2 and c = a/2 + b/2 so that [-max, max] does not overflow. Running sums of the value and
  the error are Neumaier-compensated (plain running updates drift by about n eps L1 over n panels, above the floor
  near 1,000 panels), and the result is a fresh compensated sum.
- **A panel is final** when a node of either child would not lie strictly inside the child's ends (the nodes then
  round onto the ends: on [1 - 2^-50, 1] K15's outer node 1 - 3.8e-18 rounds to 1), or when its estimate is at its
  floor. Final panels leave the heap; the loop ends when the heap is empty or the test passes; if the final panels'
  errors alone exceed the tolerance, the result is `not_converged` with the best estimate. Draft 1's "nodes are
  interior, so an end is never evaluated" was false in floating point: x^-0.99 on [0, 1] bisects into the subnormals
  until a node rounds to 0 and fn(0) = inf (numerics review).
- **Stagnation:** if k bisections in a row do not reduce the total estimate, the result is `not_converged` with the
  best estimate (a noisy fn otherwise runs to the cap; QUADPACK counts roundoff the same way). Phase 2 fixes k.
- **The cap** comes from T's exponent range and digits, not from its bits alone: x^-0.9 on [0, 1] needs about 200
  bisections in double at the default tolerance and about 765 at `cpp_bin_float_50`'s; roots' 4 times the bits (212
  for double) would fail the first. x^-0.9 and x^-0.99 are tests.
- **Overflow:** a value, an L1 or an estimate that overflows while every value of fn is finite (1e308 sin x on
  [0, 2 pi]) is `non_finite`, not a floor of inf that passes every test.
- **`TanhSinh`:** levels of halved step until the relative difference of successive levels meets the test, with at
  least 4 levels, plus the end-term test (Bailey's d4): whenever the sum is cut, at the end of the table or at a
  dropped non-finite value, the last term kept, times h, must be negligible against the tolerance, or the level is
  not converged. Without it the estimate fails silently: x^-0.99 on (0, 1) in double ends the table at x = 1e-275,
  missing a mass of 0.178 (1.8e-3 relative), while successive levels differ by only 2.2e-5 relative, and draft 1's
  squared difference (4.8e-10) would have reported success with an estimate six orders too small. The end term is
  1.12 there, so the test refuses it. Bailey's full formula (on digits) is an alternative to the plain difference;
  Phase 3 picks one, the end-term test is not optional.
- **`Romberg`:** the difference of the last two diagonal entries, after at least 4 levels.
- **The exactness degree** of the K(2n+1) panel is 3n + 1 for even n and 3n + 2 for odd n (K15: n = 7, degree 23);
  G(n) is exact to degree 2n - 1.
- Validation as in roots: a negative or non-finite tolerance is `non_finite`; tolerances below the floor are raised
  to it.

### 3.7 Errors

```cpp
enum class error {
    non_finite,    // an end is NaN (or infinite, for a method that rejects it), a tolerance is negative or not
                   // finite, fn returned a value that is not finite where the method cannot drop it, or a sum overflowed
    not_converged, // the tolerance was not met: the cap, stagnation, or the final panels' errors exceed it
};
```

- One enumeration for the module, formattable, documented as extensible, as in the other modules; carried in
  `Failure<T>::reason` (3.3), returned alone by the fixed rule.
- **fn not finite at a node, Gauss-Kronrod:** every K(2n+1) rule and tanh-sinh evaluate the centre of the interval,
  so log |x| or 1 / sqrt |x| on [-1, 1] hits its singularity on the first panel. Two choices (D14): `non_finite`, and
  the documentation says to split at the singularity; or bisect a splittable panel whose interior node gave a
  non-finite value, so that the point becomes an end of both children, whose nodes never reach it until they become
  final (then `non_finite`). The second costs one rule in the loop and turns a common failure into a result; a
  non-integrable singularity (1/x on [-1, 1]) still fails, by `not_converged` as its panels' errors grow.
- **fn not finite, `TanhSinh`:** a non-finite value in a tail is dropped and the tail cut there, as Boost does,
  because that is the endpoint singularity the method exists for; the end-term test (3.6) is what makes it sound: 1/x
  on (0, 1) has end terms of about pi cosh t (634 at t = 6), so it fails the test rather than passing on finiteness.
  A non-finite value elsewhere is `non_finite`. Phase 3 defines "tail".
- **`GaussLegendre`:** a non-finite value is `non_finite`; there is nothing to split.

### 3.8 fn

- `nxx::detail::callable<F, T>` (`math.hpp`), called through `std::invoke` with a `T const` lvalue, its result
  converted with `static_cast<T>`: as in the other modules.
- **The distance to the end, `(x, xc)`, is deferred** (simplicity review). One method would use it, and the caller
  has a workaround: substitute u = b - x so that the singularity moves to 0 (and split at the midpoint for
  singularities at both ends). `TanhSinh`'s documentation shows the reflection. When it comes, the one-argument form
  wins where both are viable (a generic variadic lambda accepts both), so adding it breaks nothing.
- **Interior singularities and break points:** not searched for; split at a known singularity or kink and add the
  parts. Break points as an argument (QUADPACK's QAGP, SciPy's `points`) are deferred, as a named aggregate (3.1).
- **fn returning an expression template** (`[](mp x) { return x * x; }` with Boost's expression templates on) returns
  an object that refers to the destroyed parameter: dangling, and not caught. The same holds in `nxx::deriv` and
  `nxx::roots` today, whose tuple results already require exactly T. A rule for value results (require the result to
  be exactly T when T is a class type) belongs in `docs/design.md` for every module; it is raised with the user
  separately, not decided here.

### 3.9 Inside

As in the other modules: `integrate` validates and dispatches once with `if constexpr` to `detail::gauss_legendre`,
`detail::gauss_kronrod`, `detail::tanh_sinh`, `detail::romberg`, each with its own loop. Gauss-Kronrod keeps its
panels in a `std::vector` used as a binary heap (`std::push_heap`; transient allocation, so constexpr), reserved from
the cap. No step interface, no shared mutable state. Everything constexpr; nothing throws; fn's exceptions pass
through.

### 3.10 What `math.hpp` gains

Nothing. `TanhSinh` needs exp (Phase 3); the interlacing route (3.4) avoids pi and cos for the Gauss rules, and if it
does not fit the constexpr limits, the spike's pi and cos come in Phase 1. The spike's versions, arithmetic only: pi
by Machin's formula and ln 2 as 2 atanh(1/3), summed in double-word arithmetic, correctly rounded; exp by Cody-Waite
reduction (k found by comparisons, since T to int is not an operation `nxx::real` has), a Taylor core on |r| <= 0.35
and 2^k by squaring: within 1.12 ulp on [-700, 700] in double and long double and 1.05 ulp in `cpp_bin_float_50`
(whose own exp reaches 885 ulp near |x| = 600). They stay in `quadrature.hpp`'s `detail` until a second module needs
one, by `docs/design.md`'s rule; then they move to `math.hpp` as public helpers. D11.

### 3.11 Shared rules this module follows

User types (`+ - * /`, comparisons, `nxx::math`, `T(1)`, no `auto` on arithmetic), constexpr throughout, nothing
throws, `std::expected`, callbacks named `fn`, qualified internal calls, Emscripten. The header includes `math.hpp`
only.

## 4. Usage

```cpp
namespace ig = nxx::quadrature;

// The default: adaptive Gauss-Kronrod. Cooling (T2 < T1) gives a negative duty without a branch.
auto const q = ig::integrate(cp, T1, T2).transform([n](ig::Integral<double> const &I) { return n * I.value; });

// Tolerances, by name (a braced {1e-10} would set abs_tol).
auto const s = ig::integrate(fn, 0.0, 2.0 * pi, {.abs_tol = 1e-12, .rel_tol = 1e-10});

// The best estimate when the tolerance was not met.
auto const r = ig::integrate(noisy, 0.0, 1.0, {.rel_tol = 1e-12});
if (!r && r.error().best) report(r.error().best->value, r.error().best->abs_error);

// An infinite interval: TanhSinh (the Gauss methods reject an infinite end).
double const inf = std::numeric_limits<double>::infinity();
auto const g = ig::integrate<ig::TanhSinh>([](double x) { return std::exp(-x * x); }, -inf, inf); // sqrt(pi)

// An endpoint singularity at 0; at b, reflect: u = b - x.
auto const t = ig::integrate<ig::TanhSinh>([](double x) { return std::log(x) / std::sqrt(x); }, 0.0, 1.0); // -4

// A fixed rule inside Newton: smooth in T2, a plain number, and the derivative is cp(T2) exactly.
auto const gap = [&](double T2) {
    auto const I = ig::integrate<ig::GaussLegendre<20>>(cp, T1, T2);
    return std::pair{I.value_or(nan) - heat, cp(T2)};
};

// Multiprecision: ask for the digits the type was chosen for.
auto const m = ig::integrate<ig::TanhSinh>(fn_mp, mp(0), mp(1), {.rel_tol = mp("1e-45")});

// Compile time, for fn built from + - * / (std::exp is not constexpr in C++23).
constexpr auto c = ig::integrate([](double x) { return 1 / (1 + x * x); }, 0.0, 1.0); // pi / 4
```

- A nested integral needs the inner tolerance tighter than the outer one, or the outer sees the inner's noise and
  reports `not_converged`; the documentation says so.
- `using namespace` for two modules makes `error` and `Criteria` ambiguous: a compile error, never a silent choice;
  the documentation uses the alias.

## 5. Phases

1. **The core, the fixed rule and Romberg:** `Integral`, `Failure`, `error`, `Criteria`, `integrate<M>`;
   `GaussLegendre<N>` with the generator and its compile-time tables (3.4); `Romberg`, the first method with an
   estimate, which needs no nodes. It fixes the API before the methods that use it most.
2. **Adaptive Gauss-Kronrod**, the default: the Kronrod generator, the heap, the floor, final panels, stagnation, the
   cap, the initial mesh, compensated sums; the cache across calls if D10 says so.
3. **exp (internal) and `TanhSinh` on finite intervals**, with the endpoint handling (2.2) and the end-term test.
4. **`TanhSinh` on half lines and the line** (exp-sinh, sinh-sinh).

Draft 1 had six phases; Gauss rules on infinite intervals left the plan (D13), Romberg moved into Phase 1, and the
distance-to-end fn was deferred. Each phase costs an API document, three reviews, five toolchains, CI and a Codex
review, so the decisions that would break callers later (the entry point, the result, the failure, the criteria)
come first.

## 6. Deferred

Non-breaking to add later unless marked:
- K21 ... K61 (`GaussKronrod<N>` keeps its parameter for them, D7); Gauss-Lobatto, Clenshaw-Curtis, Patterson's
  nested rules.
- Gauss rules on infinite intervals, with QUADPACK's map (3.5).
- fn taking the distance to the end, `(x, xc)` (3.8).
- Break points, as a named aggregate or a separate function, never a bare span (3.1).
- The result's extra members: the number of evaluations, the integral of |f| (D12).
- A cache of nodes across calls (if D10 defers it); compile-time tables for long double where they do not fit.
- QUADPACK's scaled estimator and Bailey's full formula, behind a measurement with the numerics review's protocol:
  false positives against the returned value, several tolerances, noisy and lower-precision fn, x^alpha, kinks and
  jumps, randomised Lyness-Kaganove peaks, the distribution of log10(estimate / error), the same floor and mesh for
  both, double and both 50-digit types.
- QUADPACK's epsilon-algorithm extrapolation (QAGS; Wynn 1956, de Doncker 1978) for endpoint singularities under
  Gauss-Kronrod; `TanhSinh` serves them now.
- `Criteria` members for the caps.
- Weight functions, oscillatory integrands, Cauchy principal values, several dimensions, complex-valued fn.
- A partial form for composition (`integral_of<M>(fn)`), as `derivative_of`, when a use appears.

## 7. Sources (verified by the literature survey, 2026-10-30)

- G. H. Golub, J. H. Welsch, Calculation of Gauss quadrature rules, Math. Comp. 23 (1969) 221-230.
- N. Hale, A. Townsend, Fast and accurate computation of Gauss-Legendre and Gauss-Jacobi quadrature nodes and
  weights, SIAM J. Sci. Comput. 35(2) (2013) A652-A674.
- F. Johansson, M. Mezzarobba, Fast and rigorous arbitrary-precision computation of Gauss-Legendre quadrature nodes
  and weights, SIAM J. Sci. Comput. 40(6) (2018) C726-C747 (reference nodes for tests).
- A. S. Kronrod, Nodes and Weights of Quadrature Formulas: Sixteen-Place Tables, Consultants Bureau, 1965.
- D. P. Laurie, Calculation of Gauss-Kronrod quadrature rules, Math. Comp. 66 (1997) 1133-1145.
- R. Piessens, E. de Doncker-Kapenga, C. W. Ueberhuber, D. K. Kahaner, QUADPACK, Springer, 1983.
- M. A. Malcolm, R. B. Simpson, Local versus global strategies for adaptive quadrature, ACM TOMS 1(2) (1975) 129-146.
- W. Gander, W. Gautschi, Adaptive quadrature - revisited, BIT 40(1) (2000) 84-101.
- P. Gonnet, A review of error estimation in adaptive quadrature, ACM Comput. Surv. 44(4) (2012), article 22.
- L. F. Shampine, Vectorized adaptive quadrature in MATLAB, J. Comput. Appl. Math. 211(2) (2008) 131-140.
- J. N. Lyness, When not to use an automatic quadrature routine, SIAM Review 25(1) (1983) 63-87.
- J. N. Lyness, J. J. Kaganove, Comments on the nature of automatic quadrature routines, ACM TOMS 2(1) (1976) 65-81.
- D. K. Kahaner, Comparison of numerical quadrature formulas, in J. R. Rice (ed.), Mathematical Software, Academic
  Press, 1971, 229-259.
- W. Romberg, Vereinfachte numerische Integration, Det Kongelige Norske Videnskabers Selskabs Forhandlinger 28(7)
  (1955) 30-36; F. L. Bauer, H. Rutishauser, E. Stiefel, New aspects in numerical quadrature, Proc. Symp. Appl. Math.
  15, AMS (1963) 199-218.
- H. Takahasi, M. Mori, Double exponential formulas for numerical integration, Publ. RIMS Kyoto Univ. 9(3) (1974)
  721-741; M. Mori, M. Sugihara, The double-exponential transformation in numerical analysis, J. Comput. Appl. Math.
  127 (2001) 287-296.
- K. Tanaka, M. Sugihara, K. Murota, M. Mori, Function classes for double exponential integration formulas, Numer.
  Math. 111 (2009) 631-655 (exp-sinh and sinh-sinh, as Boost cites them; to verify in Phase 4).
- D. H. Bailey, K. Jeyabalan, X. S. Li, A comparison of three high-precision quadrature schemes, Experimental
  Mathematics 14(3) (2005) 317-329 (the error estimate read from the 2003 preprint; the journal version to check in
  Phase 3).
- J. Vanherck, B. Soree, W. Magnus, Tanh-sinh quadrature for single and multiple integration using floating-point
  arithmetic, arXiv:2007.15057 (2020).
- Test batteries with closed forms: Kahaner's 21 (most closed forms; 8, 12, 13, 17, 18 need high-precision values),
  Bailey et al.'s 15, the QUADPACK book's parameterised families, Lyness and Kaganove's peaks, Gonnet's families;
  Boost.Math's tests as a list of cases (their closed forms re-derived, not their values copied).

## 8. Decisions

- **D1. One entry point `integrate<M>(fn, a, b[, criteria])`, methods as a closed set of empty tags, `GaussKronrod<15>`
  the default.** Recommended: yes.
- **D2. The ends as two positional, oriented arguments** (b < a negates, a == b is 0). Recommended: yes.
  Alternative: an aggregate `Interval{a, b}` as roots' `Bracket`.
- **D3. Failure carries the best estimate:** `std::expected<Integral<T>, Failure<T>>`, `Failure{reason, best}`, best
  present for `not_converged`. Recommended: yes (changed from draft 1; breaking to change later). Alternative: the
  code only, as roots.
- **D4. `Integral<T>{value, abs_error}` with a plain T `abs_error`; the fixed rule returns `std::expected<T, error>`
  and takes no criteria.** Recommended: yes (changed from draft 1). Alternatives: draft 1's `std::optional`; or the
  fixed rule estimates |G_N - G_(N-1)| for N - 1 more evaluations and returns an `Integral` like the others.
- **D5. Criteria `{abs_tol = 0, rel_tol = sqrt(eps)}`, the test `err <= max(abs_tol, rel_tol |I|, floor)`, the
  floor about 50 eps L1.** Recommended: yes. Alternative: roots' full-precision defaults.
- **D6. Gauss-Kronrod's estimator |K - G|, fixed now; the scaled form deferred behind a measurement.** Recommended:
  yes (changed from draft 1, which measured both in Phase 2).
- **D7. Methods: `GaussLegendre<N>` (N a range), `GaussKronrod<15>` (one value; the parameter kept so that K21 ...
  K61 need not break callers), `TanhSinh` (three kinds of interval), `Romberg`; four phases (section 5).**
  Recommended: yes. Alternative: a plain `GaussKronrod` (the two-uses rule read strictly; K21 ... K61 later as a new
  name).
- **D9. Names.** The namespace: `nxx::integrate` holding a function `integrate` stutters
  (`nxx::integrate::integrate`), and `using namespace nxx; integrate(fn, a, b)` finds the namespace, not the function,
  an ill-formed call with a confusing message (usage sketch); no other module names its namespace after its function
  (`nxx::roots::solve`, `nxx::deriv::derivative`). Recommended: `nxx::quadrature` (header `quadrature.hpp`, test
  prefix `quadrature:*`; Boost's term), the function `integrate`. Alternatives: `nxx::integrals` (as `nxx::roots`);
  keep `nxx::integrate` and name the function `quad`. The rest: `Integral`, `Failure`, `Criteria`, `GaussLegendre`,
  `GaussKronrod`, `TanhSinh` (or `DoubleExponential`), `Romberg`.
- **D10. Nodes and weights from one arithmetic-only generator, `static constexpr` tables for float and double; other
  T computed per call.** The cache across calls for other T: recommended in Phase 2, because nested integrals in
  multiprecision are a main use and the measured cost dominates them; alternative: defer it until a phase measures it.
- **D11. exp (and pi and cos, if needed) internal to the header** until a second module needs them. Recommended: yes.
- **D12. The result's extra members (evaluations, the integral of |f|) deferred,** non-breaking. Recommended: yes.
  Alternative: both now, as the usage sketch asks.
- **D13. Infinite ends for `TanhSinh` only; the Gauss methods return `non_finite`.** Recommended: yes (changed from
  draft 1). Alternative: Gauss-Kronrod through QUADPACK's map as well, so that the default method accepts them (a
  phase of its own).
- **D14. A non-finite value at an interior node of a splittable Gauss-Kronrod panel bisects it** (the point becomes
  an end). Recommended: yes. Alternative: `non_finite`, with the documentation saying to split by hand.

D8 of draft 1 (fn with the distance to the end) is deferred (3.8).

## 9. What the reviews changed (draft 1 to draft 2)

From the numerics review (all adopted):
- The floor: about 50 eps L1 per panel, L1 by the same rule, and difference weights; "a few eps" never ends on f = 1
  with plain-T weights (3.6).
- Nodes round onto the ends: final panels, which leave the heap (3.6).
- Bailey's estimate was misdescribed; the end-term test added to `TanhSinh`, without which x^-0.99 passes with an
  estimate six orders too small (2.3, 3.6).
- QUADPACK's scaled estimator fails on noisy fn (2.3, D6).
- Stagnation detection, a cap from the exponent range (3.6).
- The infinite-interval map, for when the Gauss maps come (3.5).
- The centre node hits a midpoint singularity (3.7, D14).
- Compensated sums, overflow as `non_finite`, the midpoint formed without overflow (3.6).
- The spike's numbers corrected: a factor of about 3.5, not 10, and up to 8 eps on f = 1 and e^x, not 3; the
  adversarial bound; long double's table ends where it overflows; Emscripten's binary128 (3.4).
- An initial mesh; Romberg's ends and minimum levels; the exactness degrees; abs_error never below the floor (3.2,
  3.3, 3.6).

From the usage sketch:
- The default method and `[[nodiscard]]` were missing from the declaration (3.1).
- The namespace finds itself under `using namespace nxx` (D9).
- The composition argument for draft 1's D3 does not hold; the best estimate is wanted (D3).
- `*r->abs_error` on every adaptive call and undefined behaviour on the fixed rule: the fixed rule returns a plain T
  (D4).
- Nested integrals in multiprecision spend most of their time on nodes (D10).
- fn returning an expression template dangles, in every module (3.8; raised separately).
- Documentation: the inner tolerance of a nested integral, designated initializers for criteria, `max()` for
  infinity, the default tolerance of a 50-digit type (3.5, 3.6, section 4).

From the simplicity review:
- Gauss rules on infinite intervals dropped from the plan (D13); the distance-to-end fn deferred (3.8).
- The estimator study replaced by |K - G| (D6); `GaussKronrod` reduced to N = 15 (D7); `GaussLegendre`'s N a range.
- The interlacing route for the Legendre zeros, which avoids pi and cos (3.4, 3.10).
- Six phases to four, with Romberg, a roadmap requirement, in Phase 1 (section 5).
- The roadmap's "every result carries an error estimate" against draft 1's empty optional (D4).
- Break points never as a bare span (3.1).
- Reproducing Boost's suspected defects cut (2.2: the scaling defect is confirmed by reading).
- Not adopted: none of its items was rejected outright; its "keep" list (D1, D2, the generator, internal helpers, the
  heap, `TanhSinh` as one tag) stands. Its keeping of draft 1's D3 is overruled by the usage sketch's case for the
  best estimate (D3).

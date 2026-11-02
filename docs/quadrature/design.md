# Design: quadrature (`nxx::quadrature`)

Decisions in force for `nxx::quadrature`. Each one says why; change it here before changing the code. What holds for
every module (the library, real types, the general test rules) is in `docs/design.md`; the architecture is in
`BRIEF.md`, the phases in `docs/plan.md`, and each phase's specification and evidence in `PHASE<n>_API.md`.

## Sources

- Gauss-Legendre: the rule, its weights `2 / ((1 - x^2) P_N'(x)^2)` and its truncation error
  `(b - a)^(2N+1) (N!)^4 / ((2N + 1) ((2N)!)^3) f^(2N)(xi)` are textbook. The zeros are isolated by a Sturm sequence
  (the three-term recurrence of an orthogonal family); no published listing was used. References for the tests:
  Boost.Math 1.92's 115-digit tables (BSL-1.0, values only, with its notice in `tests/quadrature_reference.hpp`) and
  Golub and Welsch (1969), computed with mpmath (below).
- Romberg, W. (1955), Vereinfachte numerische Integration; the error bound and the positivity of the extrapolation's
  weights from Bauer, F. L., Rutishauser, H. and Stiefel, E. (1963), New aspects in numerical quadrature, as Stoer and
  Bulirsch state the bound.
- Neumaier's compensated sum (1974); the rounding floor's constant 50 from QUADPACK (Piessens et al. 1983, the book);
  the running error analysis of the tests' error model after Higham (2002), Accuracy and Stability of Numerical
  Algorithms, ch. 3.
- Gauss-Kronrod (Phase 2): Kronrod, A. S. (1965), Nodes and Weights of Quadrature Formulas; Laurie, D. P. (1997),
  Calculation of Gauss-Kronrod quadrature rules, Math. Comp. 66 (the Stieltjes polynomial and the interpolatory
  weights); Adams' closed form for the integral of three Legendre polynomials; the global adaptive loop after Malcolm
  and Simpson (1975), Piessens et al. (1983, QUADPACK's QAG) and Shampine (2008, `quadgk`'s mesh of 10); Kahaner's
  21 problems (Kahaner 1971) as Gander and Gautschi (2000) list them, for the battery. Boost.Math's 115-digit G7-K15
  table (BSL-1.0, values only) as the tests' reference.
- Not used: GSL, QUADPACK's Fortran (no licence on netlib), Numerical Recipes. Numerixx 1.x's `integrate/` was
  surveyed (`BRIEF.md` 2.1); none of its code is kept.

## The interface (`quadrature.hpp`)

- **A method is a type** (Phase 3a, `PHASE3A_API.md`, approved 2026-11-02; amends `BRIEF.md` D1 and 3.9): a
  `fixed_rule` (`apply(fn, a, b)`), an `iterative_method` (`start(fn, a, b)`, `advance(fn, state)`, `settle(state)`)
  or an `estimating_method` (`estimate(fn, a, b, criteria)`), each a const member checked by a public concept
  `<M, T, F = T (*)(T const &)>` (F by default what `integrate` passes for a plain function: a pointer to it). `GaussLegendre<N>` is a fixed rule, `Romberg` and `GaussKronrod<15>` iterative methods;
  a user's type works as they do. Adding a method edits nothing: no list of methods, no dispatch on them. The
  contract (what a method is given and must return) is written once, in the concepts' comments.
- **Six forms of `integrate`:** the named forms of Phases 1 and 2, unchanged, `integrate<M>(fn, a, b)` for a fixed
  rule and `integrate<M = GaussKronrod<15>>(fn, a, b, criteria = {})` for a method with an estimate, now
  `integrate(M{}, ...)`; the method as a value, first, `integrate(m, fn, a, b[, criteria])`, for a method with
  parameters (a braced list is never deduced as a method; with M named, `integrate<M>({...}, fn, a, b)` takes one as
  the method); and with an observer, `integrate<M>(observer(o), fn, a, b[,
  criteria])` and `integrate(m, observer(o), fn, a, b[, criteria])`. Criteria for a fixed rule do not compile.
  `GaussKronrod<N>` requires N == 15 on the type, so `GaussKronrod<21>` fails where it is named (Phase 2, Q11);
  `GaussLegendre<N>`'s range is a constraint on its member, so `GaussLegendre<0>` is a type that `integrate` rejects.
- **Methods compute, `integrate` decides** (`PHASE3A_API.md` B, C): an iterative method reports a `Step<T>{current,
  attainable, ready, exhausted, trusted}` after its start and after each step, and `integrate` runs the loop: each step
  checked (a current result sound, attainable finite and not negative), then tested against `BRIEF.md` D5's request,
  `abs_error <= max(abs_tol, rel_tol |value|, attainable)` once ready, settled before any decision, advanced while not
  exhausted. The method keeps its own limits (Romberg ready from level 4 and attainable F; Gauss-Kronrod attainable 2F,
  exhausted by an empty heap, the cap or stagnation, untrusted after an undone descent, settled on sums formed
  afresh), so that no result changed: the 2,438 results of the protocol's battery (`docs/open-closed/quadrature-spike`)
  are those of Phases 1 and 2 apart from Q4 below.
- **The ends in order** (`PHASE3A_API.md` Q4): every method receives a < b, and `integrate` negates for b < a, so
  `integrate(fn, b, a)` is exactly `-integrate(fn, a, b)` for every method. It changed Romberg's `abs_error` on
  reversed ends where its plain L1 sum, summed in the other order, moved the floor (10 of the battery's float and
  double results), and Gauss-Legendre's and Romberg's values in `cpp_dec_float_50`'s guard digits (13 results); no
  status. It also gives Romberg's exactly 0 result on reversed ends the sign -0, as Gauss-Kronrod's already had (the
  independent review; -0 == 0).
- **Checked results** (Q5): in every build, a value that is not finite, an estimate that is not finite or negative,
  or a best with `non_finite`, is returned as `non_finite`; the library's methods never trip it.
- **Observing** (Q14): `observer(o)` wraps o, so that the overloads are told apart by type, never by calling an
  argument on trial (a bare callable broke a call with a generic lambda as fn, `bare_observer.cpp`). `o(step)` after
  the start and after each step, for the ends in order; the observed method shares its state, so the run is the same
  computation. Counting evaluations is the caller's (in fn).
- **No raw loops** in the code Phase 3a wrote or restructured (the user's request): Romberg's level adds its midpoints
  and extrapolates with `std::ranges::all_of` over `std::views::iota` (stopping at the first term that cannot be formed),
  Gauss-Kronrod's mesh likewise, its fresh sums with `std::ranges::for_each` over the panels not removed, its commit with
  `std::ranges::for_each`; the descent's work list and `integrate`'s loop are each the whole of their function. Phases
  1 and 2's other loops are unchanged.
- Both ends deduce T; criteria as `Criteria<std::type_identity_t<T>>`. fn as in every module (`nxx::detail::callable`).
- `Integral<T>{value, abs_error}`, `Failure<T>{reason, best}` (best present for `not_converged`),
  `Criteria<T>{abs_tol = 0, rel_tol = sqrt(eps)}`, `error{non_finite, not_converged}`; `error` and `Failure` format as
  the enumerator's name (Q7), through one `detail::name`.

## Behaviour

- **Validation, in order** (A): the tolerances (a method with an estimate; negative or not finite: `non_finite`, then
  `rel_tol` clamped to at most 1), the ends (NaN or infinite: `non_finite`), a == b (0, or `{0, 0}`, without calling
  fn). b < a negates: since Phase 3a in `integrate`, which passes the method a < b (it was the signed scale's).
- **The interval** (`detail::Interval`, built once by `detail::interval_of` and used by both methods): the nodes
  `c + hw u`, c = a / 2 + b / 2, hw = b / 2 - a / 2 (no overflow on [-max, max]), each clamped into
  [min(a, b), max(a, b)] by comparisons (`Interval::point`): where the halves round (on [D, 3D], D = denorm_min,
  c = hw = 2D) the map left the interval, reaching 0 and 4D (independent review, finding 2). A point may coincide with
  an end, as on any interval a few ulps wide. Romberg evaluates its ends as fn(a) and fn(b).
- **Whether b - a overflows** is decided before it is formed, since constant evaluation rejects the overflow: ends of
  opposite signs with `|a| / 2 + |b| / 2 > max / 2`. In radix 2 the halves are exact and the rounded half-sum exceeds
  max / 2 exactly when the rounded |a| + |b| exceeds max (finding 3).
- **The scale of each term** (B), formed before the sum, as `term(w, y)` = hw w y, w the weight relative to the
  half-width (a Gauss weight; Romberg's 2^(1-k), 1 for T_0's terms). The order of the products is chosen so that an
  intermediate is subnormal or overflows only if the term itself is:
  - b - a finite, |b - a| >= 1: `((b - a)(w / 2)) y`, the scale first. It lies between w / 2 and max, so a tiny y on a
    wide interval (1e-310 on [0, 1e10]: 1e-300) is not made subnormal first, which would cost about 1000 eps.
  - b - a finite, |b - a| < 1: `(b - a)((w / 2) y)`, the value first, so that a subnormal scale is never formed
    ([0, 3 denorm_min]: `(b - a) / 2^k` rounds 1.5 denorm_min to 2 and 3/8 denorm_min to 0 by level 3).
  - b - a overflows, so |hw| > max / 2: `(hw y) w` for |y| <= 1, where `hw y` cannot overflow and is normal unless
    y = 0; `hw (w y)` for |y| > 1, where |w y| >= w is normal; never `hw w`, which overflows for w > 1 (w_1 = 2 at
    N = 1). The first version formed `hw (w y)` for every y, which made the constant denorm_min on [-max, max]
    integrate to 0 (finding 1).

  The specification writes the scales as `s_i = (b - a)(w_i / 2)` and `h_k = (b - a) / 2^k`; the orders are
  algebraically equal and round differently in the normal range (by an ulp for a Gauss weight; Romberg's powers of 2
  are exact in radix 2, so its terms agree there). A clarification of B, adopted in the review round and completed by
  the independent review (findings 1 and 6). w / 2 and 2^(1-k) are exact in radix 2; in radix 10 they round, which
  Romberg's 14 eps L1 allowance covers.
- **Sums** (C): value sums Neumaier-compensated (`detail::Sum<T>`), Romberg's L1 sums plain. Finiteness is checked
  once per sum: GaussLegendre's result; Romberg's L1 per level (a non-finite value of fn makes it non-finite) and the
  estimate (R(k-1, k-1) is finite, so an overflow of T_k or of the extrapolation, values near max of alternating signs,
  makes it non-finite; `BRIEF.md` 3.6). A test shows both overflows of the extrapolation caught there: in the estimate
  at level 1, and in R(2, 2) at level 2.
- **Constant evaluation** (finding 3): an overflow, or arithmetic on a NaN, is not a constant expression, and GCC and
  Clang reject them at different operations, so a check after the arithmetic cannot report them. Under
  `if consteval` only, each value of fn is checked for finiteness, and each product of a term (`small_product`), each
  partial sum, L1 and Romberg's extrapolated values are kept within max / 2 (`large`) before the operation that could
  overflow, else `non_finite`. Conservative: a term or sum above max / 2 is `non_finite` in constant evaluation where
  run time gives the value (0.8 max on [0, 1]). At run time nothing changes: one check per sum.
- **No default constructor** is assumed: every array element is constructed from T(0) (`detail::zeros`), as
  `nxx::real` does not require T() (finding 4).
- **The floor** is `50 max(eps L1, denorm_min)`, L1 Romberg's trapezoid sum of |fn| at the level reached (for x^2 on
  [0, 3] at level 4, 9 + 9 / 512, not 9); abs_error = max(estimate, floor).
- **Samples that all vanish** give the value 0 with abs_error the floor, 50 denorm_min in double, as fn = 0 on
  [-max, max] does (the user's decision, recorded in `PHASE1_API.md`).

### GaussLegendre<N> (D)

- **The generator** (`detail::gauss_legendre_rule<N, T>`, + - * / in T, no pi, no cos): the positive zeros of P_N
  from the largest. Each is isolated by bisection on the Sturm count `detail::zeros_above` in [0, upper], capped at
  `digits log2(radix) + 8` steps; upper is 1 for x_1, then the previous isolation's lower end, whose count is exactly
  k - 1, so that x_k <= upper < x_(k-1) by counting alone (the previous computed zero only where the cap ended that
  isolation). Then `detail::bracketed_zero`: Newton on P_N safeguarded by bisection, the bracket's orientation from the
  zero's index (P_N > 0 above x_1, sign (-1)^k on (x_(k+1), x_k)); the first iterate, the midpoint, checked strictly
  inside the bracket; a Newton step taken only where |P| < |P'| (hi - lo), so that it cannot overflow and P' != 0
  (both needed in constant evaluation), and only where it lands strictly inside the bracket or rounds back to x
  (convergence: a step of 0); a stop when the step is at most 2 eps |x|, or when it stops shrinking once below
  sqrt(eps) |x|; the result asserted strictly inside its isolation bracket. Rejecting the step that rounds back to x
  bisected away from a converged zero (G7's x_2 one ulp off).
- **The stagnation stop never fires.** Instrumented over every N from 1 to 64 in float, double, long double,
  `cpp_bin_float_50` and `cpp_dec_float_50`, no Newton step ever failed to shrink, and with the rule removed every table
  is the same bit for bit and every test passes. It stays, as `PHASE1_API.md` D requires it, a guard for types whose
  rounding stalls Newton above 2 eps |x|; no test can need it.
- **Tables:** `detail::gauss_legendre_table<N, T>`, an `inline constexpr` variable template, used for float, double
  and long double (`detail::tabled`, named); the generator per call for every other T. `detail::Rule<T, N>` stores the
  non-negative half, ascending, 0 first for odd N.

### Romberg (E)

- T_0 and the levels as in E; two rows of the Richardson table; the estimate |R(k, k) - R(k-1, k-1)|; a stop at
  k >= 4 when estimate <= max(abs_tol, rel_tol |R(k, k)|, floor); the cap at 16 levels (65,537 evaluations), then
  `not_converged` with best = {R(16, 16), max(estimate, floor)}; no stagnation rule (Q3, Q4). The estimate overstates
  the error of R(k, k) by 30 to 10^4 times on smooth fn in double (the prototype).

## The tests' error model of the generator

`tests/test_quadrature_gauss.cpp` holds every computed node and weight to an a-posteriori first-order bound:

- `running_bound(n, x, u)`: the rounding errors of the computed P_n(x) and P_n'(x) by a running error analysis
  (Higham): the local error of each step of the recurrences from its own operations under fl(a op b) = (a op b)(1 + d),
  |d| <= u = eps / 2, each carried to degree n by its exact linear propagation (the homogeneous recurrence for an error
  in p or p', and its coupling c_k = (2k + 1) / (k + 1) from p into p'), computed by running the linearised
  recurrences from each degree. Local errors enter in absolute value, propagation with its signs.
- The node: |dx| <= (|residual| + error_p) / |P_N'(x)|, the residual the generator's own P_N at its result, in T.
- The weight, relative: 2 |x| / (1 - x^2) |dx| (d ln w / dx at a zero, by the Legendre equation), 6u for the formula,
  2 error_dp / |P_N'|.
- Measured against the 13 referenced rules, the largest ratio of actual error to bound: nodes 0.27 (float), 0.26
  (double), 0.24 (80-bit long double), 0.44 (`cpp_bin_float_50`, with or without expression templates), weights 0.15,
  0.11, 0.13, 0.09; `cpp_dec_float_50` 1e-11 (its guard digits). The node bound is 0.5 to 2 eps against a measured 0.003
  to 0.33 eps (N = 64: 0.51 eps at the inner node against 0.016, where the first model of this phase allowed 80 eps).
  The weight bound is looser, because the running bound on P' sums its local errors in absolute value (100 to 5,800 u
  relative against a few u actual): its sum over N = 64 is 484 eps against a measured 9.3. So the tests also pin the
  measured quality per referenced rule, max |dx_i| and sum |dw_i| within 1.5 times the values below, labelled as
  regression baselines (CLAUDE.md's exception for a bound on a quality, beside the correctness check against the
  references).
- The per-node rounding bound of D, `|hw| sum_i (|dw_i| |f_i| + w_i |f'_i| (|hw| |dx_i| + eps (|c| + 1.5 |hw x_i|)))`
  (the roundings of c, hw, hw x_i and c + hw x_i), plus `(2 eps + rho + 2 N eps^2) sum |t_i|` for the terms' three
  operations, fn's own rounding rho and the compensated sum: with the actual dx and dw where a reference exists
  (N = 1, 2, 3, 5, 7, 10, 15, 20, 25, 30, 40, 50, 64), the model's elsewhere.

## Evidence (Phase 1)

### The Romberg prototype

`scratchpad/romberg/proto3.py` (session scratchpad, not kept): the specification in Python floats (IEEE double),
mirroring the header's operation order (`detail::Interval::term`, the Neumaier sum, the level and Richardson updates,
the checks), so that the counts check a transliteration, as in `nxx::roots`; references from mpmath at 60 digits. The
C++ reproduces its evaluation counts exactly for fn of `+ - * /` and sqrt (`quadrature: romberg: evaluation counts in
double`); the rows with libm functions are checked for accuracy, not counts. Neither the scale order of the review
round nor the fixes of the independent review (the clamp, the order where b - a overflows) changed an earlier row.

| Case | Criteria | Status | Level | Evaluations | abs_error | Actual error |
|---|---|---|---|---|---|---|
| e^x on [0, 1] | default (1.5e-8) | ok | 4 | 17 | 3.35e-10 | 3.32e-14 |
| e^x on [0, 1] | rel 1e-14 | ok | 6 | 65 | 1.91e-14 | 7.75e-17 |
| e^x on [1, 0] | default | ok, 1 - e | 4 | 17 | 3.35e-10 | 3.32e-14 |
| 1 / (1 + x^2) on [0, 1] | default | ok | 5 | 33 | 2.91e-9 | 1.21e-11 |
| 1 / (1 + x^2) on [0, 1] | rel 1e-14; rel 0 | ok | 8 | 257 | 8.72e-15 | 3.06e-17 |
| 1 / (1 + x^2) on [1, 0] | default | ok | 5 | 33 | 2.91e-9 | 1.21e-11 |
| 1 / (1 + x^2) on [0, 1] | rel 1e300 (clamped to 1); abs 1e300 | ok | 4 | 17 | 1.72e-6 | 2.9e-9 |
| sin on [0, pi] | default | ok | 5 | 33 | 5.41e-9 | 1.32e-12 |
| sin^2 on [0, 2 pi] | default | ok | 7 | 129 | 5.66e-9 | 1.38e-12 |
| sin^2 on [0, 2 pi] | rel 1e-14 | ok | 9 | 513 | 3.49e-14 | 1.22e-16 |
| x^2 on [0, 3] and [3, 0] | default; rel 1e-14 | ok (the floor) | 4 | 17 | 1.00e-13 | 0 |
| x^3 on [-1, 1] (= 0) | default; abs 0 | ok (the floor) | 4 | 17 | 5.64e-15 | 0 |
| \|x - 1/5\| on [0, 1] | default | ok | 14 | 16,385 | 2.93e-9 | 5.87e-10 |
| \|x - 1/5\| on [0, 1] | rel 1e-14 | not_converged | 16 | 65,537 | 1.83e-10 | 3.67e-11 |
| e^(-1e4 (x - 1/2)^2) on [0, 1] | default | ok | 12 | 4,097 | 1.21e-12 | 2.92e-16 |
| cos(1000 x) on [0, 1] | default | ok | 14 | 16,385 | 3.41e-14 | 3.54e-18 |
| 1e-300 e^x on [0, 1] | default | ok | 4 | 17 | 3.35e-310 | 3.34e-314 |
| sqrt x on [0, 1] | default | ok | 16 | 65,537 | 7.47e-9 | 4.09e-9 |
| sqrt x on [0, 1] | rel 1e-14 | not_converged | 16 | 65,537 | 7.47e-9 | 4.09e-9 |
| 1 / sqrt x on [0, 1] | default | non_finite (fn(0)) | 0 | 2 | | |
| 1e308 on [0, 1e-10] | default | ok | 4 | 17 | 1.11e284 | 3.9e281 (4e-17 relative) |
| 1e-310 on [0, 1e10] | default | ok | 4 | 17 | 1.11e-314 | 7.0e-317 (7e-17 relative) |
| 0 on [-max, max] | default | ok | 4 | 17 | 2.47e-322 (50 denorm_min) | 0 |
| 1e-300 on [-max, max] | default | ok | 4 | 17 | 4.0e-6 | 1.1e-8 (3e-17 relative) |
| 1e300 on [0, 3 denorm_min] | default | ok | 4 | 17 | 1.65e-37 | 7.8e-40 |
| denorm_min on [-max, max] | default | ok | 4 | 17 | 1.97e-29 | 0 |
| 0 on [D, 3D], 1e300 outside (D = denorm_min) | default | ok, 0 | 4 | 17 | 2.47e-322 | 0 |
| 1e-320 e^x on [0, 1] | default | ok (the denorm_min floor) | 4 | 17 | 2.47e-322 | 1.1e-323 |
| 1e308 sin x on [0, 2 pi] | default | non_finite (L1 overflows) | 2 | 5 | | |
| -0.9 max at 0, 1; 0.9 max at 1/2 | default | non_finite (the estimate overflows) | 1 | 3 | | |
| -0.3 max at 0, 1; -0.7 max at 1/2; max at 1/4, 3/4 | default | non_finite (R(2, 2) overflows) | 2 | 5 | | |
| all-zero samples on [0, 1] | default | ok, 0 | 4 | 17 | 2.47e-322 | (a positive integral) |
| NaN at 1/2 | default | non_finite | 1 | 3 | | |
| infinite at 1 | default | non_finite | 0 | 2 | | |

### Constant evaluation of N = 64

The smallest passing limit (bisection to 2 %) for `detail::gauss_legendre_rule<64, T>()` alone, against the defaults
(Clang, em++ and cl 1,048,576 steps, GCC 2^25 = 33,554,432 operations):

| Compiler | double | long double |
|---|---|---|
| GCC 16 (long double 80-bit) | 1.98M operations (5.9 %) | 2.03M operations (6.1 %) |
| em++ (long double binary128) | 197k steps (18.8 %) | 221k steps (21.0 %) |
| clang-cl 22 | 197k steps (18.8 %) | the same type as double |
| cl | 186k steps (17.7 %) | the same type as double |

Binary128 is the dearest, just above a fifth of em++'s limit; every other case is below a fifth, as `PHASE1_API.md` D
expected.

### Run-time cost of the generator

GCC 16, -O2, per call of `detail::gauss_legendre_rule<N, T>()`: `cpp_bin_float_50` 0.11 ms at N = 7, 1.0 ms at N = 20,
10.6 ms at N = 64; with expression templates 0.98 ms at N = 20; `cpp_dec_float_50` 2.2 ms at N = 20 and 28 ms at
N = 64. A whole `integrate` call with a cheap fn costs the same to within 2 %.

### Accuracy of the computed rules

Against mpmath's Golub-Welsch values at every N from 1 to 64 (the script below), the maxima over N: nodes within 0.50
eps absolute and 1.86 eps relative; weights up to 765 eps relative (double, an outer weight: the amplification
2x / (1 - x^2)); the sum over the N nodes of |dw| at most 11.0 eps (float), 10.5 (double), 10.6 (long double,
80-bit), 11.5 (binary128) and 10.0 (`cpp_bin_float_50`). `cpp_dec_float_50` errs by less than the references' 60
digits resolve. The largest sum of |dw| / eps in each band of N:

| N | float | double | long double (80-bit) | long double (binary128) | cpp_bin_float_50 |
|---|---|---|---|---|---|
| 1 to 8 | 3.7 | 3.6 | 5.7 | 4.4 | 3.9 |
| 9 to 16 | 7.0 | 4.9 | 6.5 | 5.2 | 5.5 |
| 17 to 24 | 5.7 | 6.3 | 8.3 | 7.7 | 6.9 |
| 25 to 32 | 7.5 | 8.3 | 8.4 | 7.9 | 7.6 |
| 33 to 40 | 8.7 | 8.4 | 8.4 | 8.6 | 8.8 |
| 41 to 48 | 10.4 | 8.7 | 9.6 | 8.9 | 8.4 |
| 49 to 56 | 9.9 | 7.9 | 9.6 | 8.7 | 10.0 |
| 57 to 64 | 11.0 | 10.5 | 10.6 | 11.5 | 9.3 |

The regression baselines of the tests, at the referenced N (the maxima over the five types):

| N | 1 | 2 | 3 | 5 | 7 | 10 | 15 | 20 | 25 | 30 | 40 | 50 | 64 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| max \|dx\| / eps | 0 | 0.244 | 0.275 | 0.375 | 0.436 | 0.316 | 0.404 | 0.319 | 0.331 | 0.280 | 0.350 | 0.428 | 0.333 |
| sum \|dw\| / eps | 0 | 0 | 2.44 | 3.68 | 3.93 | 5.50 | 5.09 | 5.43 | 6.76 | 8.38 | 8.63 | 7.79 | 10.37 |

### The references for the generator

Boost.Math's 115-digit tables and mpmath's Golub-Welsch values agree to 5e-61 at N = 7, 10, 15, 20, 25, 30, and the
mpmath values give sum w = 2 and the integral of x^(2N-2) to within 2.6e-60 at every N. The script, run as
`python -I golub_welsch.py 1 2 3 5 40 50 64` for the tests' values and with every N from 1 to 64 for the measurements:

```python
# Gauss-Legendre nodes and weights by Golub and Welsch (1969): the eigenvalues of the symmetric Jacobi matrix of the
# Legendre recurrence (zero diagonal, off-diagonal k / sqrt(4 k^2 - 1)) are the nodes, and 2 v_0^2 the weights, v_0
# the first component of each normalised eigenvector. mpmath's eigsy at 80 digits, printed to 60.
# Usage: python -I golub_welsch.py N... > out ; prints "N i x w" for the non-negative nodes, ascending.
import sys
import mpmath as mp
mp.mp.dps = 80
def rule(n):
    J = mp.zeros(n, n)
    for k in range(1, n):
        b = k / mp.sqrt(4 * mp.mpf(k) ** 2 - 1)
        J[k - 1, k] = b; J[k, k - 1] = b
    E, Q = mp.eigsy(J)
    pts = sorted((E[i], 2 * Q[0, i] ** 2) for i in range(n))
    half = [(x, w) for x, w in pts if x > -mp.mpf(10) ** -70]
    if n % 2 == 1: half[0] = (mp.mpf(0), half[0][1])
    return half
if __name__ == '__main__':
    for a in sys.argv[1:]:
        n = int(a)
        for i, (x, w) in enumerate(rule(n)):
            print(n, i, mp.nstr(x, 60, min_fixed=-1, max_fixed=-1), mp.nstr(w, 60, min_fixed=-1, max_fixed=-1))
```

## Gauss-Kronrod (Phase 2, `PHASE2_API.md`)

- **The rule** (`detail::gauss_kronrod_rule<T>`, + - * / in T): GaussLegendre<7>'s nodes and weights, and the 8
  Kronrod nodes, the zeros of the Stieltjes polynomial E_8 = sum_k a_k P_k (`detail::stieltjes`: Adams' closed form,
  a triangular system with no linear solve; the conditions with even j vanish, condition j fixes a_(7-j)), each found
  by Phase 1's safeguarded Newton (`detail::bracketed_zero`, now generic in the polynomial) between two Gauss nodes or
  between the largest and 1, oriented by the index. Weights from the interpolatory formulas, 2 / (8 P_7 E_8') and
  w_G + 2 / (8 P_7' E_8). `detail::Kronrod<T>` stores the non-negative half, Gauss nodes at the even indices.
- **Tables** for float, double and long double (`detail::gauss_kronrod_table`); for every other T a **node cache**,
  for both rules (`with_gauss_legendre`, `with_gauss_kronrod`, through one `with_generated<generate>`: the generator
  in constant evaluation, else `cached<generate>()`), one function-local static per generator, so per rule and T
  whatever fn is (the first version's static lay in a function templated on fn's closure type, so each lambda
  regenerated), initialised once by the language and never changed. Safe because no `nxx::real` type changes its
  precision at run time (Boost's variable-precision types do not specialise `numeric_limits`; documented for user
  types).
- **A panel** (`detail::sample`, returning `std::expected<Panel<T>, fault>`): K15 and G7 on the same 15
  evaluations, two Neumaier sums, D = |K - G| (no difference weights, Q10), L1_p the plain sum of |term|, floor_p =
  `rounding_floor(L1_p)` = 50 max(eps L1_p, denorm_min) (shared with Romberg), e_p = max(D, floor_p), the excess
  e_p - floor_p stored with the panel. The order of evaluation: the centre, then c + d and c - d for the nodes
  ascending. One finiteness check per panel, on L1; a value or difference that overflows with L1 finite is
  `non_finite` (the prototype's check, kept; an implementation note of `PHASE2_API.md`).
- **The loop** (`detail::gauss_kronrod`): integrated over [min(a, b), max(a, b)] and negated where b < a (exactly);
  the mesh of 10 equal panels (empty ones skipped); a panel is final at creation if D <= floor or it cannot be split
  (`detail::splittable`: the midpoint, and both children's outer nodes, strictly inside); the heap by D, ties in
  creation order (a heap of indices into the panel store), growing on demand; three running Neumaier sums K, F, X
  with E = F + X; the stop test E <= max(abs_tol, rel_tol |K|, 2F); then the heap empty, the cap (64 b bisections
  and descent splits, b = digits times the bits of a digit rounded up, 12,800 in `cpp_dec_float_50`) and stagnation
  (32 bisections without a new minimum of X, the minimum starting at the completed mesh's X), each `not_converged`
  with best {K, max(E, F)}; the result {K, max(E, F)} (E can sit a hair below F through a Neumaier residue, so the
  max keeps abs_error at the floor or above exactly).
- **The panel store and the fresh sums** (independent review of 2098416, finding 1): removing a parent's
  contribution from a compensated sum is not exact (a small contribution absorbed beside much larger ones is lost
  when they are subtracted), so every committed panel is kept in commit order with a removed flag, final ones too.
  The running sums decide when to test; every exit (ok, the cap, stagnation, the heap empty) recomputes K, F and X
  afresh from the current panels, in commit order, and rechecks: where the fresh test fails where the running one
  passed, the running sums become the fresh ones and the loop continues. After a bisection whose removed contribution
  exceeds what remains of a running sum, the running sums are resynchronised, so that stagnation sees current values.
  The result and best come from the fresh sums. The store holds the mesh plus two panels per bisection (plus a
  descent's pieces): 6,794 for the comb at the cap in double (380 KB at 56 bytes; the heap alone held at most 1,138),
  about 5,650 for x^-0.875 in `cpp_bin_float_50` at its default (2.3 MB at 400 bytes). Time: x^-0.875 in
  `cpp_bin_float_50` 236 ms against 225 ms before, the comb at the cap 1.3 against 1.5 ms (noise), a kink at rel_tol
  0 4.1 against 3.3 microseconds. Resynchronisation fires often (179 times on x^-0.875 in double) and moved no count.
- **The compensated value in constant evaluation** (finding 2): `detail::compensated` refuses a sum above max / 2
  with a nonzero correction, whose final addition could overflow, for every compensated sum the methods form.
- **Values that are not finite** (D): a panel whose L1 is not finite is split, iteratively, at most 2 levels below
  the panel that started the descent, each split counted against the cap; a descent left with such a piece is undone
  (the parent stays, final, with its own estimate), and the run is tainted: it can no longer end ok. A mesh panel with
  no finite estimate is `non_finite`.
- **Constant evaluation**: a value or product `checked_term` refuses, or a sum above max / 2 (`detail::fits`), is
  `non_finite`, with no descent; so a `static_assert` and its run-time twin agree wherever both succeed. The tests run
  the bisection loop at compile time (a kink, 7 bisections) and check the run-time result equal to it bit for bit.
- **A gap:** no test has a descent that reaches the cap (a descent split when `work` is already at the cap), a path
  the prototype did not meet either.
- **The Gauss-Legendre tables are unchanged** by the generic `bracketed_zero` and the shared recurrence step: every
  rule, N = 1 to 64 in float, double and long double, is identical bit for bit to 930bb49's, on GCC (80-bit long
  double), clang-cl and em++ (binary128).
- **The stagnation stop of the Newton iteration** (Phase 1) is shared by the Kronrod zeros; it does not fire there
  either.

### The tests' error model of the K15 rule

As Phase 1's for Gauss-Legendre (`tests/test_quadrature_kronrod.cpp`, `kronrod_model`): a Kronrod node errs by at
most (|E_8 computed| + error_e) / |E_8'|, error_e the coefficients' errors against the exact rationals times |P_k|,
the computed P_k's running bounds times |a_k|, and the series' own operations; its weight c / (P_7 E_8') by the
node's error times |w'| plus the relative errors of P_7 and E_8'; at a Gauss node, Phase 1's model plus the same for
the term c / (P_7' E_8). Against Boost's 115-digit table the largest ratios of actual error to bound are, nodes and
weights: float 0.27 and 0.077, double 0.10 and 0.060, long double (80-bit) 0.17 and 0.061, `cpp_bin_float_50` (with
or without expression templates) 0.22 and 0.12, `cpp_dec_float_50` 4e-23 and 9e-24 (its guard digits). The sum of
|dw_K| over the 15 points is 3.7 (float), 2.9 (double), 3.1 (long double) and 6.1 (`cpp_bin_float_50`) eps; the
largest |dx| 0.44, 0.17, 0.33 and 0.44 eps.

The smooth-case bound (`kronrod_bound`): every panel lies in a mesh panel at most (b - a) / 10 wide, and K15, exact
to degree 23 with positive weights, errs on a panel of width L by at most 2 L E_23(f) <= 2 L (L / 2)^24
max |f^(24)| / (2^23 24!) (Chebyshev); plus the weights' errors sum |dw| max |f| (b - a) / 2, the nodes' and the
arguments' rounding max |f'| (b - a) ((b - a) / 20 max |dx| + eps (max(|a|, |b|) + 1.5 (b - a) / 20)), and
(4 eps + rho) L1 for the terms, fn's own rounding and the sums.

## Evidence (Phase 2)

### The prototype, with the C++ table

`scratchpad/gk-proto/gk3.py`, and `gk4.py` with the independent review's corrections (session scratchpad, not
kept): the loop of `PHASE2_API.md` C to E in Python floats,
mirroring the header operation for operation (the order of evaluation, `Interval::term`, the Neumaier sums, the heap
order, E = F + X, the running K as the value), with the C++ generator's double table loaded (`dump_k15.cpp`). Its
results equal the C++'s bit for bit, values and estimates, on every case the tests count:

| Case | Criteria | Status | Evaluations |
|---|---|---|---|
| 1 / (1 + x^2) on [0, 1] | default; 1e-14 | ok (the mesh alone) | 150 |
| \|x - 0.37\| on [0, 1] | default | ok | 360 |
| \|x - 0.37\| on [0, 1] | rel 0 | ok (the 2F test) | 660 |
| 1 / sqrt\|x\| on [-1, 3] | default | ok (one descent split) | 2,550 |
| sqrt(x - 0.5) on [0, 1] | default | non_finite (the bounded descent) | 105 |
| sqrt x on [0, 1] | default | ok | 360 |
| x^-0.75 on [0, 1] | default | ok | 2,790 |
| x^-0.875 on [0, 1] | default; 1e-12 | ok | 5,520; 10,710 |
| 1 / x^2 and 1 / (x - 0.5)^2 on [0, 1] | default | not_converged (stagnation) | 1,110 |
| 1 / (1e-4 + (x - 0.42)^2), NaN on \|x - 0.48753\| < 7.86e-5 | 1e-10 | not_converged (the taint; ok without it) | 630 |
| the same | 1e-12 | not_converged (final panels, with or without the taint) | 960 |
| \|127 x - floor(127 x) - 1/2\| on [0, 1] | rel 0 | not_converged (the cap, 3,392 bisections) | 101,910 |
| x^-(127/128) on [0, 1], 7 nested sqrt | default | not_converged (the heap empty: final panels) | 92,130 |
| max on \|x - 2.5\| < 0.25, 0 elsewhere, on [0, 50] | default | non_finite (L1 finite, G7's term overflows) | 15 |
| 1e100, 1 and 1e-100 at two nodes and elsewhere (review, finding 1) | rel 0 | ok, 1e-100 | 210 |
| 1e30, 1e10 and 1 at the same points | default | ok, 1 | 210 |

`counts3.py` holds every case, and reruns the holed ones without the taint rule. The counted powers are x^-0.75,
x^-0.875 and x^-(127/128), written with sqrt, in place of the spec's x^-0.9, which needs pow (a libm function, whose
last bits could move a decision). The cap is reached by the comb above, not by x^-0.9, which at rel_tol 0 in double
ends ok at the floor (17,460 evaluations). The cap's justification in `cpp_bin_float_50` holds: at its default,
x^-0.9 (by pow) takes 3,550 bisections of its 10,752 (the prototype at 168 bits: 3,472) and ends ok, x^-0.875 (by
sqrt) 2,820; the tests run the second.

### The battery (the prototype with its own double rule, `final2.py`; the tables `PHASE2_API.md` F summarises)

710 cases at four tolerances (Kahaner's 21, 200 Lyness-Kaganove peaks, 350 rational and Gaussian peaks, 21 powers,
50 kinks, 50 jumps, edge cases, noisy fn; references from closed forms or mpmath at 60 digits). False positives (the
estimate below the error) and false successes (ok, the error above both the tolerance and the estimate), and the
percentiles of log10(estimate / error), excluding the Gaussian peaks of width 1e-4 (narrower than the points'
spacing, every sample 0):

| rel_tol | ok | FP | FS | p5 | p25 | p50 | p75 | p95 |
|---|---|---|---|---|---|---|---|---|
| 1e-6 | 606 | 31 (5.1 %) | 11 (1.8 %) | -0.03 | 3.00 | 6.03 | 6.65 | 7.84 |
| 1.5e-8 (default) | 606 | 35 (5.8 %) | 12 (2.0 %) | -0.04 | 2.76 | 5.58 | 6.33 | 7.16 |
| 1e-12 | 603 | 31 (5.1 %) | 14 (2.3 %) | -0.00 | 1.33 | 2.43 | 2.88 | 3.69 |
| 0 | 591 | 28 (4.7 %) | 28 (4.7 %) | 0.02 | 0.74 | 1.50 | 1.85 | 2.65 |

The false successes are kinks and steps in a panel's unsampled outer 0.43 % (the worst: a step at 0.599763, estimate
4.4e-15, error 2.4e-4), the peaks of width 1e-3 at rel_tol 0, and powers beyond x^-0.65. Kahaner's 21 (status,
evaluations, estimate / actual error):

| Problem | 1.5e-8 | 1e-12 | rel 0 |
|---|---|---|---|
| K01 e^x [0,1] | ok 150, 1.9e-14 / 5.2e-16 | ok 150, 1.9e-14 / 5.2e-16 | ok 150, 1.9e-14 / 5.2e-16 |
| K02 step 0.3 [0,1] | ok 150, 7.8e-15 / 1.7e-16 | ok 150, 7.8e-15 / 1.7e-16 | ok 150, 7.8e-15 / 1.7e-16 |
| K03 sqrt x [0,1] | ok 360, 5.1e-09 / 2.9e-10 | ok 630, 4.5e-13 / 2.5e-14 | ok 810, 1.3e-14 / 3.0e-16 |
| K04 23/25 cosh - cos [-1,1] | ok 150, 6.0e-15 / 1.4e-16 | ok 150, 6.0e-15 / 1.4e-16 | ok 150, 6.0e-15 / 1.4e-16 |
| K05 1/(x^4+x^2+0.9) [-1,1] | ok 150, 1.8e-14 / 3.6e-16 | ok 150, 1.8e-14 / 3.6e-16 | ok 150, 1.8e-14 / 3.6e-16 |
| K06 x^1.5 [0,1] | ok 180, 2.0e-09 / 9.6e-12 | ok 330, 3.5e-13 / 1.5e-15 | ok 420, 6.4e-15 / 7.8e-17 |
| K07 1/sqrt x [0,1] | ok 1350, 2.1e-08 / 1.4e-08 | ok 2370, 1.8e-12 / 3.0e-13 | ok 3150, 4.2e-14 / 2.7e-15 |
| K08 1/(1+x^4) [0,1] | ok 150, 9.6e-15 / 2.0e-16 | ok 150, 9.6e-15 / 2.0e-16 | ok 150, 9.6e-15 / 2.0e-16 |
| K09 2/(2+sin 10pi x) [0,1] | ok 300, 2.6e-09 / 5.9e-16 | ok 900, 1.1e-12 / 3.7e-16 | ok 1350, 1.4e-14 / 1.5e-16 |
| K10 1/(1+x) [0,1] | ok 150, 7.7e-15 / 2.0e-16 | ok 150, 7.7e-15 / 2.0e-16 | ok 150, 7.7e-15 / 2.0e-16 |
| K11 1/(1+e^x) [0,1] | ok 150, 4.2e-15 / 1.1e-16 | ok 150, 4.2e-15 / 1.1e-16 | ok 150, 4.2e-15 / 1.1e-16 |
| K12 x/(e^x-1) [0,1] | ok 150, 8.6e-15 / 1.8e-16 | ok 150, 8.6e-15 / 1.8e-16 | ok 150, 8.6e-15 / 1.8e-16 |
| K13 sin(100 pi x)/(pi x) [0.1,1] | ok 2220, 1.0e-10 / 3.5e-17 | ok 4560, 8.6e-15 / 1.2e-16 | ok 4560, 8.6e-15 / 1.2e-16 |
| K14 sqrt50 e^(-50 pi x^2) [0,10] | ok 240, 3.3e-10 / 2.1e-16 | ok 330, 3.1e-13 / 1.0e-16 | ok 390, 1.0e-14 / 1.0e-16 |
| K15 25 e^(-25x) [0,10] | ok 210, 2.9e-09 / 4.4e-16 | ok 300, 4.3e-13 / 4.4e-16 | ok 420, 1.2e-14 / 2.2e-16 |
| K16 50/(pi(2500x^2+1)) [0,10] | ok 330, 2.1e-11 / 1.3e-16 | ok 510, 1.3e-13 / 1.3e-16 | ok 600, 9.6e-15 / 1.3e-16 |
| K17 50 sinc^2(50 pi x) [0.01,1] | ok 1110, 1.2e-09 / 4.4e-18 | ok 2340, 1.0e-13 / 4.4e-18 | ok 4140, 2.4e-15 / 4.4e-18 |
| K18 cos(cos x + ...) [0,pi] | ok 150, 1.6e-10 / 2.6e-16 | ok 330, 8.0e-13 / 1.5e-16 | ok 450, 1.9e-14 / 3.7e-16 |
| K19 log x [0,1] | ok 630, 1.5e-08 / 2.6e-09 | ok 1080, 8.6e-13 / 7.9e-14 | ok 1470, 2.1e-14 / 3.3e-16 |
| K20 1/(1.005+x^2) [-1,1] | ok 150, 1.7e-14 / 4.1e-16 | ok 150, 1.7e-14 / 4.1e-16 | ok 150, 1.7e-14 / 4.1e-16 |
| K21 sum sech(20^i (x - i/5)) [0,1] | ok 870, 1.5e-09 / 2.0e-17 | ok 1290, 1.3e-13 / 3.6e-17 | ok 1770, 3.0e-15 / 7.7e-18 |

### The parameter studies (prototype; the decisions in `PHASE2_API.md`)

- The mesh (Q1): false successes on 100 Gaussian peaks of width 1e-3 at the default: 36 with 1 panel, 2 with 2, 1 with
  4, none with 10; Kahaner's 21st missed with 1 or 2.
- The stop test (Q4): at rel_tol 0 the brief's test (E <= F) ended 32 of 50 kinks and 42 of 50 jumps
  `not_converged`; with 2F, 50 and 38 end ok with 4 times fewer evaluations.
- Stagnation (Q2), against the brief's rule on E and draft 1's reset on a tenfold rise: no false trip over the battery
  at four tolerances and five combs of 40 peaks; 1/x^2 and 1/(x - 0.5)^2 stop at 1,110 evaluations in double and 168
  bits; e^x (1 + 1e-10 h(x)) at 1,620 with an honest best, where it ran to the cap without a rule.
- The descent (Q8): unbounded, sqrt(x - 0.5) walked 1,063 levels (101,775 evaluations) in double; bounded at 2,
  `non_finite` in 105. The undo against carrying the ancestor's estimate: equal on the battery; on a NaN region the
  carried estimate was 3.4 times below the error, the undone one honest.
- Difference weights (Q10): no status changed and one count by 30 evaluations over the battery; |D| / (eps L1_p) on
  narrow panels of 1, e^x and 1 / (1 + x^2) at most 1.28, 1.62 and 1.52, 30 times below the floor's 50.

### Constant evaluation

The smallest passing limit (bisection to 1 %) for `integrate` in a `constexpr` variable, cheap fn: the mesh alone
(1 / (1 + x^2) on [0, 1], 150 evaluations) and a kink at 0.37 at rel_tol 0 (660 evaluations, 17 bisections); the
steps per bisection (30 evaluations) from the difference, and the bisections that fit the default limits after the
mesh:

| Compiler (default limit) | The mesh | Per bisection | Bisections that fit |
|---|---|---|---|
| GCC 16 (33,554,432 operations) | 365,649 | 76,478 | about 430 |
| Clang 22 (1,048,576 steps) | 44,136 | 9,064 | about 110 |
| em++ (1,048,576 steps) | 44,205 | 9,000 | about 110 |
| clang-cl (1,048,576 steps) | 33,915 | 6,965 | about 145 |
| cl (1,048,576 steps) | 25,329 | 5,155 | about 195 |

So a compile-time integral of a cheap fn reaches about 110 bisections (3,450 evaluations) on every compiler; a fn
dearer than a few operations proportionally fewer. The K15 table itself is computed separately, once: 12,130 steps on
em++ in binary128 long double (9,261 in double), 55,429 operations on GCC in 80-bit long double, 1.2 % and 0.2 % of
the limits.

### The node cache

`cpp_bin_float_50`, GCC 16 -O2: the generator costs 0.32 ms for G7-K15 and 1.0 ms for GaussLegendre<20>, which every
call paid before the cache. With it, and one cache per rule and T, a call of `integrate` on 1 / (1 + x^2) over [0, 1]
with a lambda of its own (50 distinct closure types, the cache built once beforehand by another) costs 0.435 ms (the
default method, 420 evaluations: 1.7 times faster than with the generator per call) and 0.012 ms (GaussLegendre<20>:
81 times faster). The first version's cache, per closure type, saved nothing on each lambda's first call.

## Tests

- References from closed forms at 70 digits (mpmath at 80), never Numerixx; Boost's and mpmath's nodes and weights
  (`tests/quadrature_reference.hpp`).
- GaussLegendre: every node and weight within the error model, and the per-rule regression baselines (13 values of N;
  float, double, long double, `cpp_bin_float_50` with and without expression templates, `cpp_dec_float_50`); at every N
  from 1 to 64 in float, double and long double, and at a sample of N in the 50-digit types, the weights sum to 2
  within sum |dw| (a real check where no reference exists) and x^(2N-2), x^(2N-1) on [0, 1] lie within the rounding
  bound; e^x and 1 / (1 + x^2) within the truncation bound (|f^(2N)| <= e, and <= (2N)! for 1 / (1 + x^2)) plus the
  rounding bound, in every type.
- Romberg: the bound of Bauer, Rutishauser and Stiefel at the level reached (from the evaluations) plus
  `(14 eps + 1.97 rho) L1`; where the truncation is negligible (e^x at 1e-14, 1/(1 + x^2) at 1e-14, the polynomials)
  that allowance is about 100 times the measured errors, a worst-case bound on rounding. The honesty check
  `|value - ref| <= abs_error` in the same loop, and alone on cos 1000 x, the peak, the kink and sqrt x; the counts
  against the prototype; the floors of x^2 and x^3 exactly (every sample and sum is exact in double).
- `static_assert` for N = 1, 2, 7, 64 (and 7 in float and long double) and Romberg's x^2 on [0, 3], with tolerances
  that a run-time test checks against the rounding bound; the failures in constant expressions.
- Gauss-Kronrod (`tests/test_quadrature_kronrod.cpp`): the rule against Boost's table within the error model, in
  every type; smooth cases within `kronrod_bound` in double and both 50-digit types, with the honesty check beside;
  Kahaner's 21 at the default, honest, the non-smooth ones within rel_tol |ref|; the counts above; the cap, stagnation
  (the hash noise, 1/x^2), the heap empty (x^-(127/128)), the 2F test, D14 on [-1, 3], the bounded descent (also in
  `cpp_bin_float_50`), the taint, the check of G7's overflow, each with its status, evaluations and best estimate in
  one table; the cap's justification in `cpp_bin_float_50`; the substitution x = u^10 for x^-0.9 cos x; static_asserts
  within their bounds, the bisection loop among them.
- The edge cases of both methods with an estimate as one template (`quadrature: estimating methods: edges`):
  orientation, a == b, NaN and infinite ends, tolerances, [-max, max], subnormal widths, the domain-checking callback,
  the scales 1e-300 and 1e300, x^3 with abs_error the floor, overflow, and the constant-evaluation cases.

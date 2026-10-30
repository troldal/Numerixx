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
- Not used: GSL, QUADPACK's Fortran (no licence on netlib), Numerical Recipes. Numerixx 1.x's `integrate/` was
  surveyed (`BRIEF.md` 2.1); none of its code is kept.

## The interface (`quadrature.hpp`)

- `integrate<M>(fn, a, b)` for the fixed rule `GaussLegendre<N>`, 1 <= N <= 64, returning `std::expected<T, error>`;
  `integrate<M>(fn, a, b, criteria = {})` for `Romberg`, returning `std::expected<Integral<T>, Failure<T>>`. Each
  overload is constrained to its methods (`detail::fixed_rule<M>`, `std::same_as<M, Romberg>`), so criteria for the
  fixed rule do not compile (`PHASE1_API.md`).
- Both ends deduce T; criteria as `Criteria<std::type_identity_t<T>>`. fn as in every module (`nxx::detail::callable`).
- `Integral<T>{value, abs_error}`, `Failure<T>{reason, best}` (best present for `not_converged`),
  `Criteria<T>{abs_tol = 0, rel_tol = sqrt(eps)}`, `error{non_finite, not_converged}`; `error` and `Failure` format as
  the enumerator's name (Q7), through one `detail::name`.
- No default method until `GaussKronrod<15>` exists (Phase 2).

## Behaviour

- **Validation, in order** (A): the tolerances (Romberg; negative or not finite: `non_finite`, then `rel_tol` clamped
  to at most 1), the ends (NaN or infinite: `non_finite`), a == b (0, or `{0, 0}`, without calling fn). b < a negates
  through the signed scale.
- **The interval** (`detail::Interval`, built once by `detail::interval_of` and used by both methods): the nodes
  `c + hw u`, c = a / 2 + b / 2, hw = b / 2 - a / 2 (no overflow on [-max, max]); Romberg evaluates its ends as fn(a)
  and fn(b).
- **The scale of each term** (B), formed before the sum, as `term(w, y)` = hw w y, w the weight relative to the
  half-width (a Gauss weight; Romberg's 2^(1-k), 1 for T_0's terms). The order of the products is chosen so that an
  intermediate is subnormal or overflows only if the term itself is:
  - b - a finite, |b - a| >= 1: `((b - a)(w / 2)) y`, the scale first. It lies between w / 2 and max, so a tiny y on a
    wide interval (1e-310 on [0, 1e10]: 1e-300) is not made subnormal first, which would cost about 1000 eps.
  - b - a finite, |b - a| < 1: `(b - a)((w / 2) y)`, the value first, so that a subnormal scale is never formed
    ([0, 3 denorm_min]: `(b - a) / 2^k` rounds 1.5 denorm_min to 2 and 3/8 denorm_min to 0 by level 3).
  - b - a overflows: `hw (w y)`, the value first, since `hw w` overflows for w > 1 (w_1 = 2 at N = 1; fn = 1e-300 on
    [-max, max] would give infinity, fn = 0 NaN). The one case left: a value of fn below about 1e-292 on a near-max
    interval makes `w y` subnormal and loses digits, though the term is normal.

  The specification writes the scales as `s_i = (b - a)(w_i / 2)` and `h_k = (b - a) / 2^k`; in the normal range both
  orders round alike. A clarification of B, recorded for the review (the review round adopted this rule). w / 2 and
  2^(1-k) are exact in radix 2; in radix 10 they round, which Romberg's 14 eps L1 allowance covers.
- **Sums** (C): value sums Neumaier-compensated (`detail::Sum<T>`), Romberg's L1 sums plain. Finiteness is checked
  once per sum: GaussLegendre's result; Romberg's L1 per level (a non-finite value of fn makes it non-finite) and the
  estimate (R(k-1, k-1) is finite, so an overflow of T_k or of the extrapolation, values near max of alternating signs,
  makes it non-finite; `BRIEF.md` 3.6). A test shows both overflows of the extrapolation caught there: in the estimate
  at level 1, and in R(2, 2) at level 2.
- **The floor** is `50 max(eps L1, denorm_min)`, L1 Romberg's trapezoid sum of |fn| at the level reached (for x^2 on
  [0, 3] at level 4, 9 + 9 / 512, not 9); abs_error = max(estimate, floor).
- **Samples that all vanish** give the value 0 with abs_error the floor, 50 denorm_min in double, as fn = 0 on
  [-max, max] does (the user's decision, recorded in `PHASE1_API.md`).

### GaussLegendre<N> (D)

- **The generator** (`detail::gauss_legendre_rule<N, T>`, + - * / in T, no pi, no cos): the positive zeros of P_N
  from the largest. Each is isolated by bisection on the Sturm count `detail::zeros_above` in [0, upper], capped at
  `digits log2(radix) + 8` steps; upper is 1 for x_1, then the previous isolation's lower end, whose count is exactly
  k - 1, so that x_k <= upper < x_(k-1) by counting alone (the previous computed zero only where the cap ended that
  isolation). Then `detail::legendre_zero`: Newton on P_N safeguarded by bisection, the bracket's orientation from the
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
double`); the rows with libm functions are checked for accuracy, not counts. The scale order of the review round
changed no row.

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

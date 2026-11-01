# Phase 3 API proposal: TanhSinh on finite intervals

Status: approved by the user on 2026-11-01, Q1 to Q10 as recommended (draft 2). One review round (numerics-reviewer, simplicity-critic, usage-sketcher),
with the prototype and the spike rerun on its questions, is folded in ("Decisions from the reviews"). Architecture:
`docs/quadrature/BRIEF.md` (approved 2026-10-30); Phase 2: `docs/quadrature/PHASE2_API.md` (done, 9cfb0f3). The
evidence is a Python prototype on Phase 2's battery and a C++ spike of exp and the node tables (Prototype, below).

## Summary

```cpp
#include <numerixx/quadrature.hpp>
namespace qd = nxx::quadrature;

// Endpoint singularities: the double-exponential rule samples ever closer to the ends, never at them.
auto const r = qd::integrate<qd::TanhSinh>([](double x) { return std::log(x); }, 0.0, 1.0);            // -1
auto const s = qd::integrate<qd::TanhSinh>([](double x) { return 1 / std::sqrt(x); }, 0.0, 1.0);       // 2

// Smooth fn at full precision in a multiprecision type: a few hundred evaluations.
auto const m = qd::integrate<qd::TanhSinh>([](mp const &x) -> mp { return fn(x); }, mp(0), mp(1),
                                           {.rel_tol = mp("1e-45")});

// Kinks, steps and narrow interior peaks: the default, GaussKronrod<15>, which subdivides. A fallback:
auto const k = qd::integrate<qd::TanhSinh>(g, 0.0, 1.0).or_else([&](auto const &) { return qd::integrate(g, 0.0, 1.0); });
```

- `TanhSinh`, an estimating method for finite intervals (BRIEF.md 3.2, D13; half lines and the line are Phase 4).
- An internal exp in T from `+ - * /`, on the domain the tables need (BRIEF.md 3.10, D11).
- Node tables in pieces per level: built at compile time for float, double and long double, cached at run time for
  other T.
- No 31-point Kronrod rule: `TanhSinh` serves high precision (closes PHASE2_API.md Q7).
- The default method stays `GaussKronrod<15>`. Everything constexpr; nothing throws.

## Public API (additions to `quadrature.hpp`)

```cpp
namespace nxx::quadrature {

    /// The tanh-sinh rule (Takahasi and Mori): the trapezoidal rule after x = tanh(sinh t), the step halved level by
    /// level. For endpoint singularities, smooth fn and full precision; weak on kinks, steps and narrow peaks inside
    /// the interval, where GaussKronrod<15> subdivides. Finite intervals (Phase 4: infinite ones).
    struct TanhSinh
    {
    };
}
```

- **`detail::estimating<M>`** admits `TanhSinh`, so the existing estimating overload,
  `std::expected<Integral<T>, Failure<T>> integrate<M>(fn, a, b, criteria)`, takes it; no new overload. `TanhSinh`
  stays a plain struct: a later option goes into `Criteria` or a new tag, so that `qd::TanhSinh` never breaks.
- **Comments** (usage sketch):
  - `integrate` and `TanhSinh`: choosing a method. An endpoint singularity or full precision: `TanhSinh`; kinks,
    steps, peaks or an unknown fn: the default, which fails less often, because `TanhSinh` fails honestly
    (`not_converged`) on interior features; fn too expensive for 150 evaluations: `GaussLegendre<N>`. The
    `or_else` fallback above (tested).
  - `TanhSinh`: the minimum cost (5 levels, about 70 to 140 evaluations); a singularity at an end other than 0 and
    the substitution u = b - x (F); constant evaluation needs fn built from `+ - * /` and costs about 275 constexpr
    steps per node on Clang for the tables, in every translation unit that uses it.
  - `error::not_converged`: TanhSinh's cause (the cap); `best` absent means the estimate is unbounded (the terms
    grow towards an end, as for 1/x: a divergent or hopeless integral) or the interval has no point of T strictly
    inside; then no value is an estimate: try another method or a substitution.
  - `error::non_finite`: replace "(an infinite interval is for a later method, TanhSinh)" by Phase 4's wording; for
    TanhSinh, fn not finite at a point that is not in a tail (D), with examples: 1/sqrt x on [0, 1] ok, 1/x
    `not_converged`, sqrt(x - 0.5) `non_finite`. An infinite end stays `non_finite` for every method until Phase 4.
- Nothing else in the public API changes.

## Behaviour

### A. Validation and the map

- As Phase 2 (PHASE2_API.md A): tolerances, then the ends (NaN or infinite: `non_finite`), a == b is `{0, 0}` without
  calling fn, b < a negates; Phase 1's `detail::Interval` with its scale order and clamp; in constant evaluation a
  value or product the Phase 1 checks (`checked_term`, `fits`, `compensated`) refuse is `non_finite`.
- **The points:** with m = a/2 + b/2 and hw = b/2 - a/2, the centre is m, the right points b - hw comp and the left
  points a + hw comp, comp the node's complement (B), in Phase 1's scale order. Each side's points are monotone in
  t, so **a point is evaluated only if a < p < b, and the first point of a side that is not ends that side** (a
  rounding cut, E). fn is never called at an end (prototype: no end evaluated in about 180,000 runs; the nearest
  points to 0 were 5.5e-318 for x^-0.99 and 6.2e-176 for log x). Boost's form (m +- hw x for the inner points) gave
  the same results on the battery, but on [0, 3 denorm_min] an inner point rounds onto b and needs a special case.
- **The centre not strictly inside** (an interval with no point of T strictly between its ends, b = nextafter(a)):
  then no point is (numerics review), so nothing can be sampled: `not_converged` with no best. Phase 2's methods
  return a value there because their clamped points are the ends; TanhSinh must not call fn at an end.

### B. Nodes, exp and the tables

- **The node at t** (c = 1, Q1): with e = exp(t), s = (e - 1/e)/2 and ch = (e + 1/e)/2 (both by their series below
  t = 1/2, where e - 1/e cancels), and q = exp(-2s): the complement comp = 1 - x = 2q / (1 + q), formed directly,
  never as 1 - x, and the weight w = ch 4q / (1 + q)^2 = dx/dt. Nothing overflows; q only underflows. Stored per
  node: (comp, w). Correction to BRIEF.md 3.4, which assumed c = pi/2 and so pi in T: x = tanh(c sinh t) is exact
  for any c > 0, and for c <= pi/2 the nearest singularity in t stays at Im t = pi/2, so the convergence rate is
  unchanged and only t_max grows, by ln(pi/2) (numerics review). Measured: c = 1 costs 15 to 30 % fewer evaluations
  on the smooth and endpoint cases in double (11 to 28 % at 168 bits), 5 to 12 % more on endpoint singularities, with
  no false success where pi/2 has two. No pi is needed.
- **Levels** (Q2): level 0 has t = 0, 1, 2, ... (step h = 1); level k >= 1 adds the odd multiples of 2^-k (step
  h = 2^-k), so each level reuses every earlier evaluation and adds as many again. At most 10 levels (h = 2^-9) for
  every type: in double 9 levels lose nothing on the battery, but in binary128 at rel_tol 1e-32 they lose 39 of 340
  cases, and 168 bits needs the 10th level too.
- **Where a table ends** (Q3): at the last node with q >= min(), T's smallest normal number (one comparison; at most 4
  nodes fewer than a cut on comp and w, in float). t_max: float 4.47, double 6.56, 80-bit long double and binary128
  9.34, `cpp_bin_float_50` 21.1. Nodes with t >= 0 to h = 2^-9: float 2,289, double 3,361, long double 4,781,
  `cpp_bin_float_50` 10,814. Not to `denorm_min`: in the subnormal tail the stored values lose their bits (1e6 eps
  relative in float, 5e13 in double) for no change in any result. Not shorter: in 168 bits a table cut at t = 9 turns
  x^-0.99 at rel_tol 1e-45 into `not_converged` after 7,030 evaluations, against ok in 468.
- **exp in T** (`detail::exp`, internal, D11), on the tables' domain only (simplicity review): arguments from
  ln min() up to t_max (at most about 21 in `cpp_bin_float_50`), results normal. ln 2 as a double-word pair from
  2 atanh(1/3); x = k ln 2 + r, with k found by a greedy comparison against a table of powers of 2 (T has no
  conversion to an integer), c1 = ln 2 cut so that k c1 is exact; an expm1 Horner core on r; one guard, x <= x_min
  (ln min() moved a few eps inward) gives 0, so that q >= min() and the search for t_max stops at q == 0. No NaN
  test, no overflow guard, no subnormal results (Phase 4's exp-sinh will need a wider exp). Error bound (spike,
  derived): the final rounding of 1 + p (0.5 ulp), Horner (0.73), the rounding of r (0.25): 1.48 ulp; measured on
  15,000 to 19,000 points per type at most 1.09 (float), 1.02 (double), 1.05 (80-bit), 0.99 (binary128), 1.05
  (`cpp_bin_float_50`). In `cpp_dec_float_50` 0.009 ulp, which holds because of its guard digits; a decimal type
  that rounds to exactly its digits is untested (documented). For comparison, Emscripten's `std::exp(long double)`
  overflows above 709, and Boost's exp in `cpp_bin_float_50` errs by 1,568 ulp at x = 1096. About 140 constexpr
  steps per call (double, Clang).
- **Accuracy of the tables** (spike, against mpmath at 220 digits): each value carries q = e^(-2s), and the error of
  s moves q by about 2s times it, so the outer nodes err by up to 474 eps relative in double (7,845 in 80-bit long
  double, 1.1e9 at t = 21 in `cpp_bin_float_50`), the same for complement and weight: a shift of the node in t,
  dt <= ds, at which comp and w are both evaluated, so the sum's error is of order eps L1 (numerics review). The
  a-priori bound, in eps (exp's error taken as 1.5 ulp): for t >= 1/2, 2.5 + 4 sinh t + 3.5 e^-t for comp and
  6.2 + 4 sinh t + 3.5 e^-t for w; for t < 1/2, 2.5 + (N + 0.5) sinh t and 6.2 + (N + 0.5) sinh t, N the series'
  additions (5 in float, 8 in double, 9 in 80-bit, 14 in binary128, 18 in `cpp_bin_float_50`). Every measured error
  is at most 0.49 of it.
- **Tables in pieces** (Q3): the unit is a piece of one level of at most 1,024 nodes, one initializer each, at
  compile time for float, double and long double (as Phases 1 and 2) and cached at run time for other T (Phase 2's
  `cached<generate>()` through `with_generated`, a function-local static initialised once by the language, never
  mutated), keyed per level, piece and T; one shape for both paths (simplicity review). The constexpr budget is per
  initializer, and a whole table does not fit: the spike's double table to 2^-9 took 94 % of Clang's and
  Emscripten's default 1,048,576 steps with c = pi/2, and long double's failed at 1.69M (Clang) and 1.79M
  (Emscripten). In pieces the largest uses at most: float 23 %, double 27 %, 80-bit long double 39 % (Clang),
  binary128 42 % (Emscripten), 8 to 14 % of GCC's limit. **Compile time:** every piece is instantiated wherever
  `integrate<TanhSinh>` is instantiated for a T with a table, whether a call reaches it or not (simplicity review):
  0.5 to 1.2 s for double, up to 1.5 s for binary128 long double, per translation unit (documented). Compile-time
  and run-time tables are bit-identical. **Run time, other T:** a call that converges at level 4 builds only the
  pieces up to level 4; about 14 us per node in `cpp_bin_float_50` (36 us in `cpp_dec_float_50`); the whole table,
  154 ms (343 ms), once per program.

### C. The sums

- **Terms:** a kept point gives v = h term(w, y), y = fn(p), term in Phase 1's scale order (hw w y, `Interval::term`)
  and multiplied by h before it is summed: the level's sum of new terms grows like V / h, so h applied after the sum
  would overflow (numerics review: 0.8 max on [0, 1] was `non_finite` at level 1, 2e307 on [0, 1] by level 3; now
  ok in 121 evaluations with honest estimates). Level 0 includes the centre, with w(0) = 1.
- **The level's value and L1:** V_0 = sum of v, V_k = V_(k-1) / 2 + sum of the new v; L1 likewise on |v|; Neumaier
  sums (Phase 1's `Sum`). The floor F_k = 50 max(eps L1_k, denorm_min) plus the restored slivers (E). A value, an L1
  or an estimate that overflows while every value of fn is finite is `non_finite`.
- **Level 0:** fn(m) first; not finite is `non_finite`. Then each side outward, t = 1, 2, ...: a cut where the point
  is not strictly inside (A) or fn is not finite in a tail (D).
- **Pruning** (Q4): after level 0, the outer terms of each side with |v| < eps L1 are dropped and the side cut there;
  at level k, a new node beyond a side's last kept node with |v| below the same threshold cuts the side. The
  prototype: 2.8 to 3.4 times fewer evaluations in double, no status changed (smooth fn need t up to about 3, not
  6.6).

### D. Values that are not finite (BRIEF.md 3.7: "Phase 3 defines tail")

- **A tail** is a point with comp < eps that lies beyond every kept point of its side. A non-finite value there is
  dropped and the side cut (the endpoint singularity the method exists for), in constant evaluation too; a non-finite
  value anywhere else is `non_finite`, with no best. Near an end of magnitude above about 1, most such points round
  onto the end and are never evaluated (not all: on [-1, 1] points with comp in [eps/4, eps) near 1 are), so a tail
  cut happens mostly at an end at or near 0.
- Against the looser rule (any point beyond the kept ones), which returns ok on sqrt(x - 0.5) over [0, 1] at rel_tol
  1e-4 after 1,631 evaluations: with comp < eps that case is `non_finite` in 5; sqrt(x - 0.01) likewise (6). A NaN
  region of width 0.1 is `non_finite` (15); 1/sqrt x with NaN below 1e-300 ends ok (131). A NaN region of width 0.02
  at 0.3 between the nodes is missed (ok, 120 evaluations): no sampling sees it (documented, as for Gauss-Kronrod).
- The end term (E) is what makes a cut sound: 1/x on (0, 1) is refused by it, not passed on finiteness.

### E. Termination and the estimate

- **The end term of a side, by how the side was cut:**
  - **By pruning:** 0 (simplicity review: the dropped term is below eps L1 = F/50).
  - **By a non-finite value** (D): from the last two kept terms of the side, Te (outermost) and Tp: 0 if Te = 0;
    unbounded if Te >= Tp; else h Te / (1 - r), r = Te / Tp, the geometric extrapolation of the terms beyond the cut
    (an upper bound for log-concave tails such as x^-alpha's, numerics review). A side with no kept node: h |w(0)
    y(0)|. Correction to BRIEF.md 3.6, which used h Te (Bailey's d4): x^-0.99 on (0, 1) at rel_tol 1e-4 then returns
    ok with an error of 0.077, 8 times the tolerance; near a singular end neighbouring terms decay with r close to 1,
    and the tail beyond the cut is Te r / (1 - r), 18 times Te at h = 2^-7.
  - **By rounding** (A: the next point rounds onto the end; numerics review, finding 1, with the prototype's
    correction): the sliver between the last kept point and the end, narrower than an ulp of the end, is part of the
    interval. Its share of the rule is restored: V gains term(h W, y_e), W the sum of the table's weights beyond the
    last kept node up to the current level (from the table, no fn calls; precomputed suffix sums or a loop, chosen by
    measurement), y_e the last kept value; and its size |term(h W, y_e)| joins the floor. The end term is then the
    excess of a power law through the last two kept points of distinct positions, at distances d_e < d_p from the
    end (formed as end - p, exact there) with values y_e, y_p: with R_e = d_e |y_e|, R_p = d_p |y_p|, rho = d_e /
    d_p and q = R_e / R_p, R_e max(0, (1 - rho) q / (rho (1 - q)) - 1), unbounded if q >= 1 or there are no two such
    points. It is 0 for a constant and an upper bound for x^a. Without this, fn = 1 on [1, 1 + 1e-9] ended
    `not_converged` (the sliver's mass, 2.2e-7 relative, reported as an error that never shrinks), as did 1 on
    [1000, 1001] at rel_tol 0 and 1e300 on [0, 3 denorm_min]; now ok in 89, 109 and 27 evaluations, exact. The
    reviewer's first form (the rectangle (end - p_e) y_e into V and the floor, no end term) gave false successes:
    (1 - x)^-0.9 ok after 120 evaluations with an error of 0.24 against an estimate of 0.030, since for x^a the mass
    beyond the last point is 1 / (1 + a) rectangles. With the power-law excess: `not_converged`, best 0.23 +- 0.34
    (honest; draft 1's estimate was 18 times below the error), and (1 - x)^-0.5 ok at the default with an error of
    1.0e-8 against an estimate of 2.5e-8 (draft 1's was below the error).
- **Unbounded** is a flag, not inf: GCC rejects a computed overflow to inf in a constant expression, and other
  compilers reject inf - inf and inf times 0 (spike), and `nxx::real` does not promise an infinity. Every product
  that could overflow (h Te / (1 - r), with 1 - r down to an ulp) is checked first, as Phase 1's `large` does.
- **The estimate:** E_k = max(|V_k - V_(k-1)|, F_k) + the two end terms. The plain level difference (Q5): successive
  levels roughly square the error, so the difference bounds the previous level's. Rejected: d1^2 / d2 (15 to 40 %
  fewer evaluations, 2 % of estimates below the error) and Bailey's formula on digits (needs log; 16 %).
- **The stop test** `E_k <= max(abs_tol, rel_tol |V_k|, 2 F_k)`, from level 4 on (Q2), after Romberg's level loop and
  minimum: Phase 2's 2F carries over (with F alone the test cannot pass, since E >= F plus an end term: 27 of 358
  cases end ok at rel_tol 0); level 4 rather than Boost's 3 costs 2 to 9 % and removes the remaining estimates below
  the error (3 at 1e-6).
- **No stagnation rule** (Q6; simplicity review): a rule on the end terms (above the target and not halving, two
  levels in a row) changed no status on the battery; without it, integrals that cannot converge run to the cap:
  1/x on (0, 1) 5,164 evaluations instead of 326 in double (13,157 instead of 825 at 168 bits), 1/x^2 3,022 instead
  of 198, with the same status and no best. Boost's rule (the error grows) wrongly ends 20 to 72 interior cases.
- **The cap:** 10 levels; then `not_converged`, with best unless the estimate is unbounded.
- **The result:** `{V_k, E_k}`; E >= F, so abs_error is never below the floor. **The best on failure** is the last
  level's {V, E}, absent when E is unbounded (Q7; the terms grow towards a cut, as for 1/x, 1/x^2, x^-0.999, so the
  value is no estimate), and absent when the centre is not strictly inside (A). `Integral::abs_error` stays finite,
  as in Phases 1 and 2.

### F. What the estimate is, and the documented limits

Measured by the prototype's final configuration on Phase 2's battery, split because tanh-sinh is one global rule: S,
the 360 smooth and endpoint cases (Kahaner's 21, the powers, edge and noisy cases, peaks of width 0.1 and 1e-2), and
I, the 250 interior features (kinks, steps, peaks of width 1e-3). The tables go into `docs/quadrature/design.md`,
only the numbers the decisions cite (simplicity review).
- **On S** (ok, estimates below the error, false successes; evaluations against Gauss-Kronrod's):

  | rel_tol | ok | below | false successes | median log10(estimate / error) | evaluations | Gauss-Kronrod |
  |---|---|---|---|---|---|---|
  | 1e-6 | 355 | 0 | 0 | 6.41 | 71,918 | 186,315 |
  | 1.5e-8 (default) | 355 | 0 | 0 | 5.60 | 76,402 | 206,175 |
  | 1e-12 | 353 | 0 | 0 | 1.76 | 126,486 | 288,285 |
  | 0 | 352 | 0 | 0 | 1.76 | 134,016 | 338,685 |

  Kahaner's 21: 19 ok at every tolerance, none with its estimate below its error; the step (K02) and the sech sum
  (K21) `not_converged` with honest estimates.
- **Endpoint singularities:** in double 3 to 9 times fewer evaluations than Gauss-Kronrod at the median; log x on
  [0, 1] in 112 evaluations with an error of 1.1e-16 (Gauss-Kronrod: 630, 2.6e-9); x^-0.9 at the default in 154,
  error 2.2e-15 (Gauss-Kronrod: 6,870, a false success). Smooth fn: comparable (median ratio 0.8 to 1.5).
- **High precision** (168 bits, 24 cases): at rel_tol 1e-45 28,937 evaluations in all, against 1,727,382 for G15-K31
  and 2,929,350 for G7-K15; within 1.05 times K31 on smooth fn, 20 to 580 times cheaper on endpoint singularities;
  at the default 17,349 against 684,604 and 533,640. So no K31 (Q8).
- **Interior features** (documented: use the default): ok on 61, 21, 11 and 10 of the 250 at the four tolerances,
  with 5 and 3 false successes at 1e-6 and the default (Gauss-Kronrod: 238 to 250 ok). The trapezoidal sum converges
  like h for a step and h^2 for a kink, so the levels cannot reach the tolerance; the failures are honest.
- **A singularity at an end that is not 0** is limited by the spacing of T near that end: (1 - x)^-0.5 ends ok at the
  default but `not_converged` at 1e-12, (1 - x)^-0.9 `not_converged` with an honest best. Documented with the
  substitution u = b - x, which moves it to 0, as a tested example.
- **Phase 2's edge cases** (the edge template), TanhSinh's status: 0.25 on [-max, max] ok; 1 on [-max, max]
  `non_finite` (the value overflows); 0 on [-max, max] {0, 50 denorm_min}; 1e300 on [0, 3 denorm_min] ok; e^x on
  [0, 1000 denorm_min] ok; 1e300 x / denorm_min on [denorm_min, 64 denorm_min] ok; reversed ends ok; 1e308 sin x
  `non_finite`; an odd fn ok with abs_error the floor; 0.8 max on [0, 1] ok; b = nextafter(a) `not_converged`, no
  best.
- fn's noise; the minimum cost of 5 levels (about 70 to 140 evaluations); the compile time per translation unit (B).

## Prototype

`scratchpad/ts-proto/` (session scratchpad, not kept; its tables go into `design.md`): the loop of A to E in double
and in mpmath at 24, 64, 113 and 168 bits, generic over the arithmetic, nodes from exp alone, Neumaier sums, Phase 1's
scale order, fn's exceptions as inf or NaN; Phase 2's battery, with GK15 from Phase 2's prototype and GK31 from
Boost's 115-digit table (checked to 4.2e-117). Its exp is mpmath's or Python's, so for the counted tests the phase
loads the C++ tables into the prototype first, as Phase 2 did. The spike, `scratchpad/spike-ts/`: exp, the tables'
budgets on the five toolchains, run-time generation, accuracy, and constant evaluation with inf.

## Tests (`tests/test_quadrature_tanh_sinh.cpp`, test cases `quadrature:tanh_sinh:*`)

- **exp** against mpmath at 60 digits on its domain (B) in float, double, long double and both 50-digit types,
  within the derived 1.5 ulp; 0 at and below x_min; `static_assert` at both ends of the domain.
- **The tables** against mpmath at 60 digits: comp and w within the a-priori bound of B, per node, in double and
  `cpp_bin_float_50`; the node counts per level.
- **Accuracy on smooth and endpoint cases** against closed forms or mpmath: `|value - ref| <= abs_error` and within
  the tolerance, on e^x, Kahaner's 1, 3, 6, 7, 19, x^-0.9, log x, sqrt(x) log x, in double, `cpp_bin_float_50` with
  and without expression templates and `cpp_dec_float_50` (at rel_tol 1e-45 too).
- **Counts** exactly in double, for fn of `+ - * /` and sqrt, against the prototype with the C++ tables loaded: e^x,
  1/sqrt x, x^-0.9 by sqrt, a peak; a handful, not the battery.
- **Termination:** the cap (hash noise at rel_tol 0; 1/x and 1/x^2: `not_converged`, no best), the geometric end term
  (x^-0.99 at 1e-4: `not_converged`, best present), the rounding cut (1 on [1, 1 + 1e-9] ok; (1 - x)^-0.9
  `not_converged` with best), the 2F test (rel_tol 0 ends ok on e^x).
- **Values not finite** (D): sqrt(x - 0.5) and a NaN region `non_finite`; 1/sqrt x with NaN below 1e-300 ok.
- **Edges:** Phase 2's edge template for the estimating methods run for TanhSinh too, with the statuses of F, the
  domain-checking callback failing if called at or outside an end; b = nextafter(a).
- **Constant evaluation:** `static_assert` on polynomial and rational cases with run-time twins; 1/x on [0, 1] in
  double, on GCC and Clang: `not_converged`, not a compile error; a NaN fn is `non_finite` there.
- **The substitution example** for (1 - x)^-0.9, and the `or_else` fallback.

## Documentation

- `quadrature.hpp`: `TanhSinh`'s Doxygen comment and the comments of the API section (above).
- `docs/quadrature/design.md`: the decisions in force, the sources (Takahasi and Mori 1974; Mori and Sugihara 2001;
  Bailey, Jeyabalan and Li 2005; Cody and Waite 1980 for exp's reduction), the numbers the decisions cite.
- `docs/sphinx/quadrature/reference.md`: the new declaration.
- `demo/quadrature_demo.cpp`: the strong endpoint singularity also by `TanhSinh`, beside the substitution.
- `docs/plan.md`: Phase 3's checkboxes; PHASE2_API.md's Q7 closed (no K31).

## Decisions for the user

- **Q1. The constant c = 1, not pi/2.** Recommended: yes (no pi in T; the convergence rate is the same; 15 to 30 %
  fewer evaluations on smooth and endpoint cases, 5 to 12 % more on endpoint singularities). Alternative: pi/2, the
  literature's choice (pi by Machin's formula costs 2.5 % of a table).
- **Q2. Levels: h = 1 down to 2^-9, at most 10 levels for every type; the first exit after level 4.** Recommended:
  yes. Alternative: 9 levels for float and double (no loss there; about half the compile time), 10 for the rest.
- **Q3. Tables to q >= min(), in pieces of at most 1,024 nodes per level: compile-time for float, double and long
  double, cached at run time per piece for other T.** Recommended: yes. Alternative: tables cut where comp < eps^2
  (4 times smaller and faster in `cpp_bin_float_50`; x^-0.99 at full precision fails).
- **Q4. Pruning of negligible outer terms.** Recommended: yes (2.8 to 3.4 times fewer evaluations, no status changed).
- **Q5. The estimate: the plain level difference, never below the floor, plus an end term per side by how it was
  cut** (0 after pruning; geometric after a non-finite value; the restored sliver and a power-law excess after
  rounding). Recommended: yes (corrects the brief's single end term; no false success on the battery).
- **Q6. No stagnation rule; the cap of 10 levels bounds every run.** Recommended: yes (no status changed on the
  battery). Alternative: the end-term rule, which makes failing runs on non-integrable ends 10 to 16 times cheaper
  (1/x: 326 evaluations instead of 5,164 in double; Gauss-Kronrod stops it at 1,110).
- **Q7. No best when the estimate is unbounded** (or nothing can be sampled). Recommended: yes (`abs_error` stays
  finite). Alternative: a best with abs_error = inf, which needs T to have an infinity.
- **Q8. No K31.** Recommended: yes.
- **Q9. The default stays `GaussKronrod<15>`;** the comments say when to use which, with the `or_else` fallback.
  Recommended: yes.
- **Q10. exp on the tables' domain only** (no overflow guard, no subnormal results, no NaN test), widened in Phase 4.
  Recommended: yes (bit-identical to the full version on the domain; less code to test).

## Decisions from the reviews

From the numerics review (the prototype and the spike rerun on each):
- Adopted: terms multiplied by h before summing, V_k = V_(k-1)/2 + the new terms (C; draft 1 overflowed on 0.8 max
  over [0, 1]); a zero end term is 0 and a side with no kept node has an end term (E; draft 1 never ended on fn = 0);
  the rounding cut's sliver (E), in the prototype's form after the reviewer's first form gave false successes;
  unbounded as a flag and the overflow checks before multiplying (E; the spike: GCC rejects a computed inf); the test
  tolerances derived (exp 1.5 ulp, the tables' per-node bound) instead of measured maxima plus 1 %, and the figures
  for c = 1; the tail rule in constant evaluation (D); the corrected wording on points near an end (D); the second
  case of no best (A, E); the guard digits behind the radix-10 figure (B).
- Confirmed: the node formula and c = 1, t_max and the counts, the level structure, the difference as an estimate,
  the geometric end term as an upper bound on log-concave tails, pruning's bound, the min() cut, the consistent shift
  in t, exp's guard, the per-piece budget, monotone points, the 2F test, the cap, the tail rule.
From the simplicity review:
- Adopted: the claim that only reached levels are compiled corrected (B); no stagnation rule (E, Q6); exp on the
  tables' domain only (B, Q10); a pruned side's end term 0 and h Te / (1 - r) (E; no result changed); one
  (level, piece) shape for both paths (B); the nextafter case stated as the centre not strictly inside (A);
  design.md limited to the cited numbers; Phase 1 and 2 helpers reused by name.
- Not adopted: 9 levels for the tabled types (Q2's alternative; binary128 needs 10), unless you choose it.
From the usage sketch:
- Adopted: the method choice and the `or_else` fallback in the comments, tested; what an absent best means; the
  `non_finite` comment corrected for Phase 4 and given the tail examples; the substitution and the constant
  evaluation cost in `TanhSinh`'s comment; the minimum cost.
- Not adopted: a combined integrator (two uses not shown); fn returning `std::expected` (raised before, scheduled
  after Milestone 4).

## Deferred

- Infinite intervals: Phase 4 (exp-sinh, sinh-sinh), with a wider exp.
- fn taking the distance to the end, `(x, xc)` (BRIEF.md 3.8): the substitution covers it.
- exp in `math.hpp` as a public helper, when a second module needs it (D11).
- A stagnation rule (Q6), Bailey's estimate, d1^2 / d2 (E: measured, not adopted).

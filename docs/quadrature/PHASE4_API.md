# Phase 4 API proposal: TanhSinh on half lines and the line

Status: approved by the user on 2026-11-03, Q1 to Q8 as recommended (draft 3). Draft 1's review round (numerics,
simplicity, usage) and the numerics review of draft 2 are folded in ("Decisions from the reviews"); every figure below
is the prototype's on the design approved here, rerun after each review.

It builds on:
- **`BRIEF.md`** 3.2 and 3.5 (approved): `TanhSinh` picks tanh-sinh, exp-sinh or sinh-sinh from the ends; the Gauss
  methods and Romberg reject an infinite end (D13).
- **`PHASE3A_API.md`** Q10 (approved): infinite ends stay outside the method contract until this phase makes them
  opt-in per method.
- **`PHASE3_API.md`** and **`PHASE3B_API.md`** (approved, implemented): tanh-sinh on [a, b], its rules, state and
  prototype. All of it stands; this phase adds two maps on the same machinery.

## Summary

```cpp
#include <numerixx/quadrature.hpp>
namespace qd = nxx::quadrature;
double const inf = std::numeric_limits<double>::infinity();

auto const g = qd::integrate<qd::TanhSinh>([](double x) { return std::exp(-x * x); }, -inf, inf);    // sqrt(pi)
auto const e = qd::integrate<qd::TanhSinh>([](double x) { return x * x * std::exp(-x); }, 0.0, inf); // 2
auto const p = qd::integrate<qd::TanhSinh>([](double x) { return std::pow(x, -1.5); }, 1.0, inf);    // 2
auto const l = qd::integrate<qd::TanhSinh>([](double x) { return std::exp(x); }, -inf, 0.0);         // 1
auto const o = qd::integrate<qd::TanhSinh>(qd::observer(print), fn, 0.0, inf);                      // observed too

auto const d = qd::integrate(fn, 0.0, inf); // the default, GaussKronrod<15>: non_finite, as today
```

- **A method opts in to infinite ends** with a static member `infinite_ends = true` (Q1). `TanhSinh` does; for every
  other method an infinite end stays `non_finite`.
- **`TanhSinh`'s map follows from the ends:** tanh-sinh on [a, b] (Phase 3, unchanged), exp-sinh on [a, inf) and
  (-inf, b], sinh-sinh on (-inf, inf).
- No new form of `integrate`, no new tag, no change to `Step`, `Integral`, `Failure` or `Criteria`.

## Public API (changes to `quadrature.hpp`)

```cpp
namespace nxx::quadrature {

    struct TanhSinh
    {
        /// It takes infinite ends: exp-sinh on a half line, sinh-sinh on the line.
        static constexpr bool infinite_ends = true;

        // start, advance, settle: as Phase 3b.
    };
}
```

- **The opt-in** is `M::infinite_ends`, a `static constexpr bool`, read by the estimating forms of `integrate` (a fixed
  rule takes none). Absent or false: an infinite end is `non_finite`, as today. `detail::Observed` forwards it, so an
  observed `TanhSinh` takes infinite ends; a user's combinator around a method must forward it too, and the
  concepts' comment says so.
- **The contract for a method that opts in** (the concepts' comment): it may receive a = -inf, b = +inf or both, always
  a < b after `integrate`'s ordering. Then b - a and a/2 + b/2 are not finite: the method must not form them, nor call
  fn at an infinite end.
- **Comments:**
  - `error::non_finite`: an end is NaN, or infinite for a method that does not take infinite ends (every method but
    `TanhSinh`); for `TanhSinh`, a value not finite where D's tail rule does not apply.
  - `TanhSinh`: the three maps; F's limits; the cost (e^-x on [0, inf) in 100 evaluations, at most 13 at level 0 in
    double, 43 in `cpp_bin_float_50`; level 0 calls fn at |x| up to about 1e131 in double, later levels up to about
    1e305 on a side that is not pruned).
  - `integrate`: an infinite end is for `TanhSinh` only. The `or_else` fallback is for finite intervals (the default
    rejects an infinite end, so it would return TanhSinh's failure unchanged); on a half line, split at a finite c:
    the default on [a, c], `TanhSinh` on [c, inf).
  - `GaussKronrod`, `Romberg`, `GaussLegendre`: an infinite end is `non_finite`; use `TanhSinh`.

## Behaviour

### A. What integrate does with an infinite end

As `PHASE3A_API.md` A, with step 2 changed: NaN is `non_finite`; an infinite end is `non_finite` unless
`M::infinite_ends`; a == b (finite, or the same infinity) is `{0, 0}` without calling fn; then the ends are put in
order, (inf, a) as (a, inf) and (inf, -inf) as (-inf, inf), the result negated.

### B. The maps and their table

- **The map:** [a, inf), x = a + E; (-inf, b], x = b - E; (-inf, inf), x = +-sinh u; with u = c sinh t, E = exp(u)
  and c = 3/2 (Q2), exact in every radix; tanh-sinh keeps c = 1.
- **One table for both maps,** per t >= 0, from exp on its present domain (no widening): q = exp(-u), then E = 1/q;
  sinh t and cosh t by their series below t = 1/2, as Phase 3. Its columns: q, E and c cosh t, and for the 167 nodes
  where u < 1/2 (t <= 167/512, in every type) the line's x by sinh's series.
  - exp-sinh: the node toward infinity E, weight (c cosh t) E; the node toward the finite end q, weight (c cosh t) q.
    t = 0 is one node, x = 1, weight c.
  - sinh-sinh: x = (E - q) / 2 where u >= 1/2 (bit for bit the value of a table of its own), else the series column;
    weight ((c cosh t) E + (c cosh t) q) / 2.
  - Accuracy against mpmath: below t = 1/2, x within 2.1 eps and the weights within 1.9; beyond, the error of u
    (477 eps at t = 7.07 in double, 8,010 in binary128 at c = 1) shifts each node consistently in t, the same for x
    and its weight, as Phase 3. Without the series, x near 0 errs by up to 182 eps.
- **Where it ends, per side:** toward infinity at the last node with c cosh t <= (max / 2) q, decided before E is formed
  (so 1/q cannot overflow for any T, and E and its weight stay at most max / 2); toward the finite end at the last
  node with q at least min(). Nodes with t > 0, to 2^-9:

  | | toward infinity: t_max, nodes | toward the end: t_max, nodes |
  |---|---|---|
  | double | 6.842, 3,503 | 6.850, 3,507 |
  | binary128 | 9.623, 4,927 | 9.625, 4,928 |
  | 168 bits | 21.41, 10,961 | 21.41, 10,961 |

  A shorter table changes results (at c = 1 in double, the side toward infinity cut at t = 5 turned 1/(1 + x^2)^0.6
  on the line not_converged and x^-1.01 lost its best).
- **Levels, pieces, budget:** as Phase 3 (levels 0 to 9, pieces of at most 1,024 nodes, compile-time for float,
  double and long double, cached for other T). The spike measured separate tables at c = 1: their largest piece took
  at most 45 % of a toolchain's default constexpr budget (binary128, em++), 31 % in double on Clang, and both together
  about 1.1 s (GCC) and 0.9 s (Clang) more per translation unit than tanh-sinh alone. The one table's budget and
  compile time are measured in step 1, before step 2 builds on them.
- **A latent fix of Phase 3b** in the same step: `exp_constants` forms its powers of 2 only where they fit (the spike's
  checked type showed 2^128 formed where max_exponent is 100; no result used it), and the tables' capacity bound uses
  a ceiling for the bits of a radix-10 digit.

### C. The sums

- **The term** is h w y: no half-width, the weights carry the scale. Phase 3b's shortcut past `checked_term`
  (`safe_value`, which assumes weights at most 1) does not apply; each product is checked within max / 2 in constant
  evaluation (`small_product`).
- **The centre:** a + 1 for exp-sinh (b - 1 reflected), 0 for sinh-sinh, weight c. A centre not strictly inside (a + 1
  rounding onto a, from |a| about 2 / eps) leaves no point to sample: `not_converged` without a best, as Phase 3.
- **The order:** each level evaluates the side toward infinity first (on (-inf, b] the side toward -inf, on the line
  the positive side), then the other; level 0 the centre before both. Otherwise Phase 3 and `PHASE3B_API.md` D: the
  sweeps, the level-k threshold eps L1_(k-1) / 2, the raw sum, the floor with the slivers.

### D. How each side ends, and its end term

- **A side toward a finite end:** Phase 3's rules, with the distance q in place of hw comp: a rounding cut where
  a + q rounds onto a (the sliver term(h W, y_e) and the power-law excess), the table's end (geometric), a tail (not
  finite at a distance below eps, beyond the kept points: geometric), pruning (0).
- **A side toward infinity** (no rounding cut, sliver or W; its weights reach max / 2): pruning (0); the table's end and
  a tail (the geometric term h Te / (1 - r)); and two rules of its own, settled by the prototype.
  - **The tail rule (Q3):** a value not finite at a node beyond every kept point is a tail (the side cut there) where
    E (sinh u on the line) exceeds 1/eps, or where the side's outermost kept |term| is already below the pruning
    threshold; anywhere else it is `non_finite`.
  - **The underflow rule (Q4):** a value below min() is no evidence of a negligible term where the weight reaches max / 2
    (without the rule 1e-15 x^-1.02 on [1, inf) ended ok with an error of 8.9e-7 at the default, and 1e-250 x^-1.01 ok
    with an error of 19 % at 1e-3). A node beyond every kept point with |y| < min(), after a kept value of at least
    min(), cuts the side there and is not kept: as pruning if w min() is below the threshold, else as the table's end
    (the geometric term over the last two kept values).
  - **Level 0, in order:** sweep each side outward, storing each node's term (a point not strictly inside ends the
    side; a value not finite is a tail where E > 1/eps, is pending toward infinity, else `non_finite`); form the
    threshold eps L1 from the centre and every stored term; decide each pending value (a tail if the side's last
    stored |term| is below it, else `non_finite`); apply the underflow rule to the stored nodes in increasing t; then
    prune from the outside in (Phase 3).
  - **Level k:** for each new node in increasing t before the side's cut: a value not finite is decided by the tail
    rule (against the threshold eps L1_(k-1) / 2 in term scale), then the underflow rule, then pruning, as listed.
- **Soundness:** for x^-p, p > 1, the terms in t are about c cosh t exp(-(p - 1) u), log-concave wherever
  (p - 1) c sinh t cosh^2 t > 1, which holds at every cut measured for p - 1 above 6e-4, so the geometric term bounds
  the mass beyond. At a prune point the decay rate in t is about ln(1/eps) for a dominant power-law tail, so the mass
  dropped is about 7 eps L1 at h = 2^-9, inside the floor; for a small slowly decaying part (e^-x + 1e-12 (1 + x)^-2)
  it can reach about 60 eps L1, as `PHASE3B_API.md` Q5's caveat.

### E. Termination

Unchanged: E_k = max(|V_k - V_(k-1)|, F_k) + the end terms, the test against max(abs_tol, rel_tol |V_k|, 2 F_k)
from level 4, 10 levels, no best where the estimate is unbounded.

### F. What the estimate is, and the documented limits

On S, 26 smooth, endpoint and decay cases on half lines and the line (e^-x, x^k e^-x, e^-x^2, 1/(1 + x^2),
x^-0.5 e^-x, log(x) e^-x, e^-x cos x, x^-1.5, log(x) / x^2, sech x, 1/(1 + x^4), ...): every case ok at 1e-6, the
default, 1e-12 (1e-30 and 1e-45 at 113 and 168 bits) and 0, in double, 113 and 168 bits, with no estimate below its
error and no false success; 3,031 to 4,689 evaluations in all in double, 20,218 at 168 bits and rel_tol 0. In double
at the default: e^-x on [0, inf) in 100 evaluations, e^-x^2 on the line in 99, x^-1.5 on [1, inf) in 141, x^2 e^-x in
156. Of nine further cases (scales, far ends, off-centre peaks, an interior kink), the only flagged one is the peak at
x0 = 1000 below.

The limits, in `TanhSinh`'s comment, measured in double at c = 3/2 unless noted:
- **The maps centre at 0 (the line) or a (a half line) at scale 1.**
  - A peak far from 0 on the line is missed. e^-(x - x0)^2 is ok up to x0 = 39 (22 at 168 bits, whose default
    tolerance is tighter), then an honest not_converged; from about 180 the outcome alternates between that, a best
    below its error and a false success (ok, value 0); from about 1,600 every run fails, in double as a false success
    (ok, 0 +- 2.5e-322). The same at c = 1. Split at the peak, or substitute x = x0 + u: each half of e^-(x - 1e3)^2
    split at 1e3 is ok in 85.
  - A half line far from 0 is cut by rounding near a, and abs_error grows with ulp(a): a = 1e10, ok +- 1.1e-6 (error
    4.3e-8); a = 1e15, ok +- 0.094 (error 0.0047), far above the request; from |a| about 1e16 no point is inside.
    Substitute x = a + u with fn written in u.
  - A narrow feature on the line, e^-(x / 1e-3)^2, is an honest not_converged in 33 evaluations; substituted, ok in 99.
    A wide one costs evaluations, not accuracy: e^(-x / L) on [0, inf) in 81 to 868 for L from 1e-2 to 1e6.
- **A value not finite nearer than the tail rule allows** is `non_finite`: x^k e^-x is ok up to k = 50 in double and
  16 in float, `non_finite` beyond (x^k overflows before e^-x underflows); write fn so that it stays finite where its
  factors do not, as x^k e^-x = exp(k log x - x). e^x e^-2x is ok.
- **Values below min() toward infinity** are not used: 1e-300 x^-1.5 on [1, inf) is an honest not_converged (error
  2.8e-3), 1e-250 x^-1.01 likewise. Scale fn.
- **Slow decay** (x^-1.01 on [1, inf)): an honest not_converged in double (best 99.91, error 0.09); ok at 113 and 168
  bits.
- **Oscillation and no integral** (sin x / x, sin x, 1/x, 1 on [0, inf)): not_converged; their bests, where any, are
  not reliable (sin x / x has one below its error).
- **A singular end other than 0** is limited by the spacing near it, as Phase 3: (x - 1)^-0.5 e^-x on [1, inf) is ok at
  1e-6 and an honest not_converged at the default; with u = x - 1 written in fn, ok in 111.
- **A large number written for infinity** integrates that finite interval (e^-x on [0, 1e10]: ok, exact; on [0, 1e300]:
  an honest not_converged): the user's integral, not the infinite one (correction to `BRIEF.md` 3.5, which said it
  converges to about 0).

## Prototype and spike

- **The prototype** (`es.py`, session scratchpad; kept beside `ts.py`, Q6): Phase 3's machinery generic over double and
  mpmath at 113 and 168 bits, every rule above as an option and this design as its default; its finite path equals
  `ts.py` bit for bit. Its tables come from Python's exp: the counted C++ tests load a dump of the C++ table, as `ts.py`
  does.
- **The spike** (session scratchpad): the tables in a copy of the header, accuracy, budgets on the five toolchains,
  compile and run times. Its widened exp is not needed (B).

## The work, in two steps

Each step committed, pushed and checked on the five toolchains; the user approves the next.
1. **The table** (internal), with Phase 3b's latent fix: tests against mpmath within the derived bounds, the counts per
   level, compile-time and run-time tables bit-identical (built on the heap, beyond Emscripten's stack); the budgets
   and the compile time in `design.md`.
2. **The infinite ends:** the opt-in in `integrate` and `Observed`, `TanhSinh`'s maps and sides, the comments, the
   tests below, `design.md`, the reference page, the demo.

## Tests

- **The opt-in:** every other method returns `non_finite` for an infinite end (the edge template); an observed
  `TanhSinh` takes them; a user's method with `infinite_ends` receives them ordered; a == b infinite is {0, 0} without
  calling fn; reversed infinite ends negate exactly.
- **Accuracy** against closed forms with the honesty check: S in double; two or three cases per map in
  `cpp_bin_float_50` (with and without expression templates) and `cpp_dec_float_50` at 1e-25 and 1e-45.
- **Counts** exactly in double against the prototype with the C++ table loaded, in total and per level, for a handful
  of cases per map.
- **The sides:** toward infinity, pruning (e^-x), the table's end (x^-1.01: honest best), the tail at level 0 (x^k e^-x
  ok; sqrt(10 - x) e^-x `non_finite`) and at level k (the prototype's cases: x^-2 sqrt((x - X1)(x - X2)) / x ok,
  e^-x sqrt((x - 20.49)(x - 24.59)) `non_finite`), the underflow rule at level 0 and k (1e-250 x^-1.01, 1e-15 x^-1.02:
  honest); toward the finite end, a table end (x^-0.99 e^-x on [0, inf) at 1e-4: honest best), a singular end other
  than 0 and its substitution; 1/x on [0, inf) and [1, inf), 1 on [0, inf): not_converged, no best.
- **The limits' remedies:** the split and x = x0 + u for the peak, x = a + u for the far half line, each ok.
- **Constant evaluation:** 1/(1 + x^2) on the line and 1/(1 + x)^2 on [0, inf) with run-time twins; a case to the cap
  whose budget is measured first; a class type's cost per node recorded, as Phase 3b.

## Decisions for the user

- **Q1. The opt-in: a static member `infinite_ends`,** forwarded by `Observed`. Recommended: yes, as Phase 3a's Q10
  decided (opt-in per method).
  - **The tension:** today `TanhSinh` is its only use, so the two-uses rule would hard-code `TanhSinh` in `integrate`
    (simplicity review); the member is non-breaking to add later, and a later finer opt-in (half lines but not the
    line) would break code that wrote the member now. Q10 chose the open protocol; the member is its smallest form.
  - Alternative: `TanhSinh` named in `integrate`, the member deferred.
- **Q2. The constant c of the two new maps: 3/2.** Recommended: yes. Exact in every radix, below pi/2, so no pi in T;
  tanh-sinh keeps c = 1.
  - **Evaluations:** on S, 12 % fewer than c = 1 in double (15,446 against 17,514 over the four tolerances), 5 % at
    113 bits, 8 % at 168.
  - **Honesty:** on the 35 cases, fewer false successes than c = 1 in every type: the interior kink e^-|x - 1| at 1e-6
    is a false success at c = 1 in all three types, honest at 3/2; at 113 bits the peak at x0 = 1000 is a false success
    at c = 1, not at 3/2. The off-centre limit is the same at both (F).
  - **Against it:** x^k e^-x stays ok to k = 70 in double at c = 1 (50 at 3/2), though to only 8 in float (16 at 3/2);
    e^x e^-2x is `non_finite` at c = 1 in double; sin x / x gets a best below its error at either c.
  - Alternatives: c = 1, one constant for all three maps; pi/2, the literature's, which needs pi in T and measures as
    3/2.
- **Q3. The tail rule toward infinity** (D). Recommended: yes. Alternatives: E > 1/eps alone (x^k e^-x `non_finite`
  from k = 27 in double at c = 1); any point beyond the kept ones (ok, false, on sqrt(10 - x) e^-x and two others).
- **Q4. The underflow rule toward infinity** (D), with the geometric term at the table's end and a tail. Recommended:
  yes. Alternatives refused: pruning on the term (the false successes in D); the rule for y = 0 only, keeping nonzero
  subnormal values (1e-15 x^-1.02 at c = 1: a best below its error); an end term of 0 or h Te at the table's end (false
  successes on x^-1.5 and x^-1.2 with NaN beyond 1e20, bests for 1/x and 1 on [0, inf)).
- **Q5. One table, from exp's present domain** (B). Recommended: yes: no widened exp, one generator, three columns and
  the line's 167 series nodes; no status or count changes against separate tables, values within 1 or 2 ulp.
  Alternative: exp widened to ln max() and a table per map (draft 1).
- **Q6. The prototype kept,** `es.py` beside `ts.py`, reading a dump of the C++ table for the counted cases.
  Recommended: yes.
- **Q7. Two steps** (the table, the infinite ends). Recommended: yes.
- **Q8. The limits documented, not engineered:** no centre or scale on `TanhSinh`; the substitution in its comment.
  Recommended: yes. Alternative (usage sketch): `TanhSinh{.centre = x0, .scale = L}` now, since a density centred at
  1,000 is not exotic and its false success is silent; non-breaking to add later (a defaulted member keeps
  `integrate<TanhSinh>` and `TanhSinh{}` working).

## Decisions from the reviews

From the numerics review of draft 1 (each reproduced on the prototype):
- Adopted: `Observed` forwards the opt-in (verified: an observed run was `non_finite`); the underflow rule (reproduced:
  1e-15 x^-1.02 ok with an error of 8.9e-7 at the default, 1e-250 x^-1.01 ok with 19 % at 1e-3; refined in the next
  round, D); the wide tail rule (x^27 e^-x was `non_finite` in double at c = 1); no `safe_value` shortcut; the line's
  narrow features and the far half line's worst case in F; the side order on (-inf, b]; the tests it listed
  (expression templates, a finite-side table end, a singular end other than 0, 1/x on [0, inf), a measured
  constant-evaluation cap); the contract's b - a and a/2 + b/2; the radix-10 ceiling; the corrected figures (level 5
  costs 74 % at 1e-6; 12 levels change e^-x on [0, 1e300]; the cuts' log-concavity condition).
- Not adopted: unbounded end terms where the last two kept values differ in sign: it removed some of oscillation's
  dishonest bests and not others.
From the simplicity review of draft 1:
- Adopted: no widened exp and one table (measured: no status or count changes); two steps; c without pi; the opt-in
  read by the estimating forms only; the limits merged in the comment; the tests trimmed; the document shortened; the
  tension in Q1 stated.
- Not adopted: `TanhSinh` hard-coded (Q1's alternative, for the user).
From the usage sketch of draft 1:
- Adopted: `Observed` forwards the opt-in; the fallback's limit and the split in `integrate`'s comment; the
  substitutions for a far half line and a peak; finite factors for x^k e^-x; the contract on b - a.
- Not adopted: a public trait to `static_assert` the opt-in, and a separate error for an infinite end (the comment of
  `non_finite` says it; a trait waits for a second use); a sum combinator for split integrals (deferred); a centre and
  scale now (Q8).
From the numerics review of draft 2 (each rerun on the prototype):
- Adopted: Q2 rerun with both values of c on every measure, which turned the recommendation from 1 to 3/2 (draft 2
  rested on the off-centre peak sampled too coarsely; finely sampled, its limit is the same at both); the table's end
  at max / 2 decided before 1/q is formed (the prototype had kept one node more, at max); the underflow rule written
  exactly as the prototype implements it, and refined so that it acts only after a kept normal value (at c = 3/2 the
  first form had turned off-centre peaks into false successes); the cost of the underflow rule as a limit (F), and its
  variant for y = 0 only measured and refused (Q4); the table's columns; the level-k tail cases for the tests; the
  figures requoted from runs of this design (x^k e^-x to k = 50 in double, not 60; the weights within 1.9 eps; fn's
  largest arguments at level 0).

## Deferred

- Gauss-Kronrod and the other methods on infinite intervals through a map (`BRIEF.md` D13).
- A centre and scale for the maps (Q8); a sum of integrals over split intervals.
- fn of the distance to the end (`BRIEF.md` 3.8).
- Oscillatory integrals on half lines (Ooura and Mori's double-exponential formula for Fourier integrals).
- exp widened to ln max(), when a module needs it.

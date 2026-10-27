# Phase 6 API proposal: the deferred methods (`Ridders`, `RegulaFalsi`, `Illinois`, `Steffensen`)

Status: approved by the user on 2026-10-27, P1 to P5 as recommended. One review round (numerics-reviewer, simplicity-critic, usage-sketcher) is
folded in ("Decisions from the reviews"). Architecture: `docs/roots/BRIEF.md` (approved 2026-10-25), section 3.2:
these four were deferred, not dropped, as teaching examples (user decision); the user asked for all four. Phases 4 and
2 hold the loops they extend: `PHASE4_API.md` (`enclose`, AndersonBjorck's loop and guard) and `PHASE2_API.md` (the
open methods and the hand-over).

## Summary

```cpp
namespace rt = nxx::roots;

auto const f = [](double x) { return (x - 1.0) * (x - 1.0) * (x - 1.0); };    // a triple root
auto const a = rt::solve<rt::Brent>(f, rt::Bracket{ 0.0, 3.0 });       // 156 evaluations
auto const b = rt::solve<rt::Ridders>(f, rt::Bracket{ 0.0, 3.0 });     //  86
auto const c = rt::solve<rt::Illinois>(f, rt::Bracket{ 0.0, 3.0 });    // 148
auto const d = rt::solve<rt::RegulaFalsi>(f, rt::Bracket{ 0.0, 3.0 }); // 145 (guarded: see B)

auto const g = [](double x) { return x * x - 2.0; };
auto const s = rt::solve<rt::Steffensen>(g, 1.0);  // sqrt 2, 10 evaluations (Secant 9, Newton 7)
auto const t = rt::solve<rt::Steffensen>([](double x) { return 1000.0 * (x * x - 2.0); }, 1.0); // -sqrt 2, 27
```

- Four tags, no new entry point: `Ridders`, `RegulaFalsi` and `Illinois` on a `Bracket`, `Steffensen` from a guess,
  with Phase 1's validation, zero first, tolerances, stop test, result and errors, and Phase 2's hand-over.
- **`Ridders`** (1979): per iteration the midpoint, then the zero of an exponential fit through the three values, both
  points through Phase 4's `enclose`; at most 2 evaluations per halving. Its square root is a new `detail::sqrt`, from
  + - * / (P5).
- **`RegulaFalsi`** and **`Illinois`**: AndersonBjorck's loop with the retained end's factor fixed at 1 and at 1/2,
  with its guard, a bisection after three steps that have not halved the bracket (P2): guarded, so RegulaFalsi does
  not show the textbook stall, whose cost its counts show instead.
- **`Steffensen`**: z = x + f(x), and the secant through x and z as the step; 2 evaluations per iteration; Brent at the
  first sign change. **Scale-dependent**: f is added to x, so the scale of f decides the cost and even which root is
  found (1000 (x^2 - 2) from 1 finds -sqrt 2, x^2 - 2 finds sqrt 2) (P3).
- **What they cost** (Prototype results): on 27 of Alefeld, Potra and Shi's test cases, Brent 309 evaluations,
  AndersonBjorck 327, Illinois 382, Ridders 416, RegulaFalsi 645; Ridders is the cheapest of all on odd multiple and
  flat roots ((x - 1)^3: 86 against Brent's 156). Steffensen costs as much as Secant or more from the same guess, far
  more from a far guess (x^2 - 2 from 10: 39 against 13). They are documented as the teaching examples BRIEF.md 3.2
  calls them; Brent stays the default.

## Public API (`roots.hpp`)

```cpp
namespace nxx::roots {

    struct Ridders {};     ///< Ridders' method (1979): the midpoint and an exponential fit. A bracket.
    struct RegulaFalsi {}; ///< False position, guarded by bisection (the published method's stall is cut off).
    struct Illinois {};    ///< False position with the retained end's value halved, guarded. A bracket.
    struct Steffensen {};  ///< Steffensen's method: z = x + f(x), the secant through x and z. A guess; scale-dependent.

} // namespace nxx::roots
```

- `detail::bracketing` gains `Ridders`, `RegulaFalsi` and `Illinois`; `detail::open` gains `Steffensen`;
  `values_needed` is 1 for all four.
- `detail::anderson_bjorck` becomes `detail::false_position<M>`, the factor chosen by M (three uses: AndersonBjorck,
  Illinois, RegulaFalsi); AndersonBjorck's results are unchanged bit for bit (m = 1 and 1/2 are exact in every radix,
  and nothing in the loop depends on m < 1: g_a keeps its sign, so the step's ratio stays in (0, 1)).
- Secant's step (the chord's zero from the latest point, its minimum step and direction, `PHASE2_API.md` C) becomes a
  helper shared by Secant and Steffensen; Secant's results are unchanged bit for bit.
- `Ridders` and `nxx::deriv::Ridders` share a name in different namespaces (P4); the documentation of each says it is
  not the other, and the existing mentions of the derivative method in `roots.hpp` and `docs/roots/design.md` are
  qualified as `nxx::deriv::Ridders`.
- `detail::sqrt(x)` (P5): x scaled into [1/radix^2, 1] by even powers of the radix, squared up (O(log^2) of the
  exponent: about 500 operations in `cpp_bin_float_50`, where single factors would take 10^9); Newton's iteration from
  1 in the correction form y + (s/y - y)/2, as `root_epsilon`; unscaled exactly. 0, +inf and NaN return themselves; a
  negative x is a precondition (asserted) and returns itself, so that the loops end. **Faithful** (within one ulp):
  the textbook form (y + s/y)/2 is 3 ulps off in radix 10 (y + s/y falls in the decade above), the correction form 1
  (exhaustive in decimal at 4, 5 and 6 digits; 300,000 random doubles none beyond one ulp of the correctly rounded
  `math.sqrt`; a quarter one ulp off, sqrt 2 among them). Correctly rounded would need the exact residual x - y^2, which
  Dekker's splitting gives in radix 2 with known precision, but not for every `nxx::real` (radix 10, or a type whose
  working precision is not `digits`).

## Behaviour

### A. `Ridders` (P1)

Ridders (1979): for a bracket [a, b] with values of opposite signs and its midpoint m, the exponential fit's zero is
x4 = m + (m - a) sign(f(a) - f(b)) f(m) / sqrt(f(m)^2 - f(a) f(b)). In exact arithmetic it lies between m and the end e
whose value has the other sign than f(m), so in the half the midpoint step keeps.
1. Each iteration: `enclose` at the midpoint (the bracket halves). The new bracket has m at one end and e at the other;
   the end the step dropped is d. Then x4 from that state (f(e) and f(d) are the old ends' values), through `enclose`
   (Phase 4's placement, stop test, exact zero), which keeps the computed point inside even where rounding puts it a
   few ulps past m.
2. **Formed without cancellation:** x4 = e + w (m - e), w = -g_e g_d / (root (root + |g_m|)), with g = f / s,
   s = max(|f(e)|, |f(d)|, |f(m)|) (no square can overflow; root = sqrt(g_m^2 - g_e g_d) is in [|g_m|, sqrt 2]), the
   quotient formed by `quotient` (where root underflows, it fails and the point is the midpoint). The textbook form,
   1 - |r| with r = f(m) / root, rounds r to 1 where f(e) f(d) is below eps of f(m)^2: on x - 1e-300 over [0, 1] it
   placed x4 on the end, and the search bisected through 1922 evaluations; with w, 4.
3. **Worst case:** each iteration halves the bracket at a cost of at most 2 evaluations: twice bisection's count.

### B. `RegulaFalsi` and `Illinois` (P2)

Phase 4's AndersonBjorck loop (`PHASE4_API.md` D) with the retained end's factor m = 1 (`RegulaFalsi`) or m = 1/2
(`Illinois`, Dowell and Jarratt 1971), and the same guard and placement. Without the guard, RegulaFalsi on x^20 - 1
over [0, 5] ran past 200,000 evaluations (the linear rate there needs about 1e13 steps); with it, 45. On smooth fn the
guard changes little (x^2 - 2: 20 with it, 23 without; cos x - x: 16 and 15). The documentation says it is guarded,
and that the guard is why it can be cheaper than Illinois on a triple root (145 against 148).

### C. `Steffensen` (P3)

From a guess x with f(x):
1. **z = x + f(x)**, but at least tol(x) away from x in f(x)'s direction (x + f(x) rounds to x where |f(x)| is below
   the resolution of x, and the chord would have no width); x + f(x) that would overflow: `not_converged`.
2. fn at z: not finite, `non_finite` (Phase 2's rule for every point: z leaves fn's domain easily, since the units
   are mixed); 0, the root; of the other sign than f(x), the hand-over.
3. **The step:** the secant's through (x, f(x)) and (z, f(z)), -(f(x) / (f(z) - f(x))) (z - x), which is Steffensen's
   -f(x)^2 / (f(z) - f(x)) where z - x = f(x). z - x is the difference as formed (exact by Sterbenz where |z - x| <=
   |x|), else the intended offset (where z - x could overflow). At least tol(x) long, in the direction of decreasing
   |f| along the chord, from the signs of f(x), f(z) - f(x) and z - x (Phase 2's rule). f(z) = f(x), or a quotient or
   product that cannot be formed: `not_converged`.
4. fn at the new point: as in 2, with the sign change against f(x) handed over; else it becomes x.
5. **The hand-over:** Brent's loop on the two points of the first sign change, x and z or x and the new point (as
   Secant, `PHASE2_API.md` C). The cap is Phase 2's (212 iterations in double), so at most 2 x 212 + 1 evaluations
   before the hand-over; Brent's loop after it has its own bound, with no cap (a jump at 0 from 1: 1359 evaluations
   in all; independent review).
6. **Scale-dependent:** f's scale sets z's distance from x, so it changes the cost (1000 (x^2 - 2) from 1: 27
   evaluations, x^2 - 2: 10) and which root is found: z = 1 - 1000 = -999 brackets -sqrt 2 with x = 1; from 0.5 on
   x^2 - 2 the first step lands at -1.83, and the root found is -sqrt 2 too, where Secant and Newton find sqrt 2.
   Documented with these examples; a root, not the nearest root, as for every open method.

## Prototype results

`.scratch/roots6/proto6.py` (Python, double; on Phase 4's prototype revision 2; `detail::sqrt` ported step for step),
references from mpmath. Evaluations of fn, the ends and zero first included. Every bracketing result is within its
`abs_error` of the reference (4 eps |root| of slack for the libm rows); Steffensen's roots are listed.

| Case | Bracket | Brent | AndersonBjorck | Ridders | Illinois | RegulaFalsi |
|---|---|---|---|---|---|---|
| x^2 - 2 | [1, 2] | 9 | 9 | 12 | 11 | 20 |
| x^3 - 2x - 5 | [2, 3] | 8 | 9 | 14 | 10 | 23 |
| x^20 - 1 | [0, 5] | 19 | 30 | 18 | 30 | 45 |
| (x - 1)^3 | [0, 3] | 156 | 152 | 86 | 148 | 145 |
| x - 1e-300 | [0, 1] | 3 | 4 | 4 | 4 | 4 |
| x - 1 | [-max, max] / [0, max] | 8 / 7 | 13 / 12 | 9 / 8 | 13 / 12 | 13 / 12 |
| x^3 - 2x + 2 | [-2, -1] | 9 | 10 | 14 | 12 | 20 |
| 1/x - 1 | [0.5, 3] | 11 | 5 | 12 | 12 | 25 |
| W10 near 7 | [6.7, 7.6] | 23 | 18 | 18 | 24 | 14 |
| a jump at 0 | [-1, 2] | 1076 | 51 | 778 | 51 | 1076 |
| 1/(x - 1/3), a pole | [0, 1] | 55 | 95 | 79 | 95 | 197 |
| sqrt x - 0.5 | [0, 1] | 4 | 9 | 12 | 12 | 22 |
| cos x - x | [0, 1] | 9 | 8 | 12 | 9 | 16 |
| e^x - 2 | [0, 1] | 8 | 8 | 10 | 10 | 18 |
| x e^x - 1 | [0, 1] | 9 | 9 | 10 | 12 | 23 |
| Kepler, e = 0.9, M = 1 | [0, pi] | 10 | 10 | 14 | 12 | 24 |
| ln(x / 1e-12) | [1e-20, 1] | 45 | 18 | 46 | 23 | 79 |
| atan(x - 1)^3 | [0, 3] | 155 | 153 | 86 | 149 | 147 |

Alefeld, Potra and Shi's problems (`PHASE4_API.md`'s set, powers as products):

| Problem | Brent | AndersonBjorck | Ridders | Illinois | RegulaFalsi |
|---|---|---|---|---|---|
| #1, sin x - x/2 | 9 | 9 | 14 | 10 | 24 |
| #2, n = 1 / 5 / 10 | 16 / 13 / 12 | 23 / 17 / 18 | 18 / 18 / 16 | 19 / 14 / 17 | 26 / 17 / 21 |
| #4, x^n - 1, n = 4 / 8 / 12 | 15 / 15 / 17 | 19 / 22 / 23 | 16 / 14 / 16 | 21 / 23 / 24 | 34 / 37 / 37 |
| #5, sin x - 0.5 | 10 | 10 | 16 | 11 | 18 |
| #6, n = 5 / 20 / 100 | 11 / 13 / 15 | 10 / 13 / 16 | 12 / 12 / 14 | 12 / 13 / 16 | 24 / 26 / 29 |
| #7, n = 5 / 10 / 20 | 10 / 9 / 9 | 9 / 8 / 8 | 18 / 18 / 18 | 10 / 9 / 9 | 17 / 14 / 12 |
| #8, n = 5 / 10 / 20 | 9 / 11 / 13 | 10 / 11 / 13 | 14 / 16 / 14 | 10 / 13 / 12 | 24 / 26 / 22 |
| #9, n = 5 / 15 / 20 | 8 / 7 / 7 | 10 / 8 / 8 | 12 / 14 / 14 | 9 / 9 / 9 | 11 / 8 / 8 |
| #10, n = 5 / 10 / 20 | 8 / 10 / 13 | 9 / 10 / 14 | 12 / 16 / 18 | 10 / 19 / 23 | 25 / 30 / 33 |
| #11, n = 5 / 15 / 20 | 12 / 13 / 14 | 7 / 6 / 6 | 16 / 16 / 18 | 15 / 17 / 16 | 35 / 32 / 33 |
| #12, n = 3 | 10 | 10 | 16 | 12 | 22 |
| Total (27 cases) | 309 | 327 | 416 | 382 | 645 |

Steffensen from a guess, against Secant and Newton (Phase 2's), evaluations and the root found:

| Case | x0 | Steffensen | Secant | Newton |
|---|---|---|---|---|
| x^2 - 2 | 1 / 10 / 0.5 | 10 / 39 / 11 (-sqrt 2) | 9 / 13 / 11 | 7 / 9 / 8 |
| x^3 - 2x - 5 | 2 / 3 | 9 / 45 | 8 / 11 | 6 / 8 |
| 1000 (x^2 - 2) | 1 | 27 (-sqrt 2) | 9 | 7 |
| (x^2 - 2) / 1000 | 1 | 9 | 9 | 7 |
| 1/x - 1 | 0.5 / 3 | 10 / `not_converged` (12) | 10 / `not_converged` (11) | 7 / `not_converged` (10) |
| cos x - x | 0 / 1 | 9 / 7 | 10 / 7 | 6 / 5 |
| e^x - 2 | 0 / 1 | 10 / 13 | 9 / 8 | 7 / 6 |
| Kepler, e = 0.9, M = 1 | 1 | 10 | 10 | 7 |
| e^(20 - 5000 / T) - 1 | 350 | 91 | 16 | 11 |

Where Steffensen's first iteration already brackets a root (x^2 - 2 from 1, cos x - x from 0, Kepler), the count is
mostly Brent's; its own convergence shows from far guesses.

## Tests (`tests/test_roots_classic.cpp`, test cases `roots: ...`)

- **Counts, exactly**, for fn of `+ - * /` and sqrt in double: the rows of the first table less the libm rows, and
  problems #2, #4, #7, #8, #9 and #11, for Ridders, Illinois and RegulaFalsi (comments on rows mark Ridders'
  cancellation, x - 1e-300, and the guard, x^20 - 1); Steffensen's `+ - * /` rows with the root's sign. Phase 2's
  Secant and Phase 4's AndersonBjorck counts unchanged (the shared code).
- **Accuracy:** Phase 4's conditioning terms on the same rows; the libm rows within `abs_error` and 2 ulps.
- **Steffensen's edges:** f(x) below the resolution of x (z's minimum step); x + f(x) overflowing (`not_converged`);
  fn not finite at z (`non_finite`); the minimum step with f(x) < 0.
- **`detail::sqrt`:** exact on perfect squares and even powers of the radix (`static_assert`); a short list in float
  and double within one ulp of `std::sqrt` (correctly rounded by IEEE 754); 0, a subnormal, max; one tiny and one huge
  value in `cpp_bin_float_50`, and values in `cpp_dec_float_50` against 70-digit references.
- **Constant evaluation:** Ridders and Steffensen in float, double and long double; RegulaFalsi and Illinois in double
  (their loop's arithmetic is Phase 4's).
- **Types:** Ridders and Steffensen in `cpp_bin_float_50` with and without expression templates and
  `cpp_dec_float_50`; RegulaFalsi and Illinois once each with expression templates.
- **Constraints:** each new tag rejected at the wrong start (one line each).

## Documentation

- `roots.hpp`: each tag with its source, its cost against Brent and what it teaches (Ridders: an exponential fit,
  quadratic per iteration at 2 evaluations, the cheapest on odd multiple roots, which need a sign change; RegulaFalsi:
  the stalled end, cut off by the guard; Illinois: the simplest cure; Steffensen: Newton without a derivative, at the
  price of scale dependence, with the -sqrt 2 example); one line added to the guide to choosing in `solve`.
- `docs/roots/design.md`: the sources, A to C, the tables; BRIEF.md 3.2's rows for the four corrected to the measured
  costs; the mentions of the derivative method qualified.
- `demo/roots_demo.cpp`: the new tags in section 2's comparison, and Steffensen's scale in section 10's.
- `docs/plan.md`: Phase 6's checkboxes.

## Changes during implementation

- No change to the API or the behaviour. Every count and root of the tables for `+ - * /` and sqrt is reproduced
  exactly (`tests/test_roots_classic.cpp`), and Phase 2's Secant and Phase 4's AndersonBjorck counts are unchanged
  (the shared `secant_move` and `false_position<M>`).
- The prototype's Steffensen minimum-step direction was aligned with Secant's rule (the one the C++ shares) before its
  counts became test data; no count changed.
- The test of the minimum step with f(x) < 0 is not a separate case: the rule is Secant's `secant_move`, now shared.
- `demo/roots_demo.cpp`: Ridders, RegulaFalsi and Illinois on (x - 1)^3 in section 2; Steffensen on x^2 - 2 and
  1000 (x^2 - 2), with Secant's result, in section 10; section 11's label names `nxx::deriv`'s Ridders.

## Independent review (OpenAI Codex, 840e04b)

After the implementation (840e04b, CI green), a read-only review by OpenAI Codex. A cold read by Kimi K3.1 Pro through
the Cline CLI did not complete (a dropped connection, then the ClinePass usage limit; its first run had also read the
Codex reproductions in `.scratch`, so it was not independent); its partial notes found Codex's third finding
independently and nothing else. The user approved all three fixes.
1. **`sum_fits` rejected sums that fit**: it required |s| <= max / 2 wherever |x| <= max / 2, so Steffensen on max - x
   from 0 gave up after one evaluation, though its first point z = 0 + max is the root. Now exact on the magnitudes
   (the smaller at most max less the larger, exact by Sterbenz). It only admits sums the old test wrongly rejected;
   every earlier count is unchanged.
2. **`detail::sqrt` overflowed in a type too narrow to hold radix^2** (a Boost type whose max is 3.98): sqrt 2 came out
   1.9e-6. Now such a type is not scaled, and Newton starts at max(1, s), above the root.
3. **"At most 2 x 212 + 1 evaluations" for Steffensen** holds before the hand-over only; Brent's loop after it has no
   cap (1359 evaluations on a jump at 0 from 1). Reworded.
Also added: `detail::sqrt` exact on perfect squares and even powers of 10 in `cpp_dec_float_50` (a test gap Kimi's notes
named). Codex checked and found correct, among others: Secant and AndersonBjorck unchanged bit for bit (every evaluated
point compared with the parent commit on five cases), the square root faithful (all positive `float16` values, a million
random doubles against exact comparisons, decimal sweeps at 1 to 4 digits), Ridders' bookkeeping and Steffensen's
directions. Noted for the library, not fixed: `nxx::real` does not require round-to-nearest, which the faithful square
root and the Sterbenz arguments assume.

## Decisions for the user

- **P1. Ridders as two `enclose` steps** (the midpoint, then the exponential fit's zero), in the cancellation-free
  form. Recommended: Ridders' algorithm with Phase 4's placement and stop test; the textbook form loses roots near an
  end (1922 evaluations on x - 1e-300). Alternative: the textbook formula, documented.
- **P2. RegulaFalsi and Illinois guarded**, as AndersonBjorck. Recommended: a library method must end in bounded time;
  the published regula falsi does not (past 200,000 evaluations on x^20 - 1). Its stall shows in its counts against
  Illinois and AndersonBjorck (645, 382, 327 on the paper's set). Alternative: unguarded, documented as possibly
  unbounded. Final for behaviour: an unguarded variant would be a new tag.
- **P3. Steffensen as the classic step** (z = x + f(x)) with the minimum steps and Secant's hand-over at z. Recommended:
  the method the name promises, its scale dependence (cost and root) documented with numbers. Without the hand-over at
  z, 1000 (x^2 - 2) from 1 is `not_converged` after 425 evaluations instead of -sqrt 2 after 27. Alternative: a
  scale-invariant variant (z = x + h, h from x), which is a secant method under another name.
- **P4. The name `Ridders`** in `nxx::roots`, as `nxx::deriv::Ridders` in `nxx::deriv`. Recommended: each is Ridders'
  method in its module; only a `using namespace` of both modules makes the unqualified name ambiguous (a compile
  error, not a silent choice). Alternative: `RiddersMethod`, or `ExponentialFit`.
- **P5. The square root as `detail::sqrt` in `roots.hpp`.** Recommended: Ridders is its one use; a public
  `nxx::math::sqrt` would be a contract (its domain, faithful rounding, its assert), costly to withdraw, and adding it
  later breaks nothing. Alternative: public `nxx::math::sqrt` in `math.hpp` now (the first proposal).

## Decisions from the reviews

From the numerics review:
- Steffensen returned -sqrt 2 in two rows (x^2 - 2 from 0.5, 1000 (x^2 - 2) from 1), unchecked: the table lists the
  root, C states that scale changes the root found, and the tests assert its sign.
- The square root in the form (y + s/y)/2 is not faithful in radix 10 (3 ulps, checked exhaustively at 4 to 6
  digits): now the correction form (1 ulp), as `root_epsilon`; ported to the prototype, every Ridders count regenerated
  (the paper's total 416, from 418).
- The square root looped forever on +inf or a negative x without asserts: now 0, +inf, NaN and a negative x return
  themselves.
- Steffensen's minimum step pointed the wrong way for f(x) < 0 (now from the signs of f(x), f(z) - f(x) and z - x);
  the slope used the intended offset, not z - x as formed (now z - x where exact); fn not finite at z is `non_finite`.
- "Costs more than Secant" was contradicted by the table: "as much as Secant or more". Ridders' x4 lies in the kept half
  only in exact arithmetic (the placement keeps the computed one inside). The per-problem table added. The reason
  given against correct rounding (an FMA) was wrong: Dekker's splitting, but not for every type.

From the simplicity review:
- The square root is `detail` (P5): one use.
- Ridders' x4 is computed from the state after the midpoint step (no capture; f(e) and f(d) are its ends' and d's
  values), its w by `quotient` (which also covers a root that underflows, where w was 0/0).
- Secant's step shared with Steffensen; the duplicate tests became row comments; constant-evaluation and type tests
  cut for RegulaFalsi and Illinois (their arithmetic is Phase 4's); the sqrt sweeps stay as evidence, not tests; the
  test file named by its content; one line in the guide to choosing; the demo extends sections 2 and 10.
- The derivative method's mentions in `roots.hpp` and `design.md` qualified as `nxx::deriv::Ridders`.

From the usage sketch (6 snippets, 6 misuses):
- Steffensen's silent other root (the sketch's first finding) confirmed and documented (above).
- RegulaFalsi's one-line documentation says it is guarded.
- Ridders' advantage holds for odd multiplicity (an even root has no sign change): said.
- Steffensen's worst case is in evaluations (2 x 212 + 1 before the hand-over), not iterations.

## Sources

- Ridders, C. J. F. (1979). A new algorithm for computing a single root of a real continuous function. *IEEE
  Transactions on Circuits and Systems* 26(11), 979-980. Not read (the publisher's paywall); the formula
  x4 = x3 + (x3 - x1) sign(f1 - f2) f3 / sqrt(f3^2 - f1 f2), x3 the midpoint, as the R package pracma documents it
  (`ridders`, citing the paper), matched by other independent statements (searched 2026-10-27). Not from Numerical
  Recipes' `zriddr`.
- Dowell, M. and Jarratt, P. (1971). A modified regula falsi method for computing the root of an equation. *BIT* 11,
  168-174: the Illinois factor 1/2 (as the numerics review confirms; not read).
- Steffensen, J. F. (1933). Remarks on iteration. *Skandinavisk Aktuarietidskrift* 16, 64-72: the step as textbooks
  state it (not read).

## Deferred

- Public `nxx::math::sqrt` (P5), when a second module needs a square root; non-breaking.
- Public concepts for the methods (`nxx::roots::bracketing_method<M>`, `open_method<M>`), so that a caller's generic
  code can constrain on them (the usage sketch); non-breaking.
- An evaluation count in `Root` (the teaching comparisons count fn by a wrapper); non-breaking (a member added).
- Pegasus and King's variants of false position (a fourth use of `false_position<M>`'s factor); unguarded regula falsi
  or a scale-invariant Steffensen, only as new tags.
- Steffensen's hand-over on the narrower (z, next) where z lies between x and the next iterate.
- A correctly rounded square root (Dekker's splitting in radix 2).

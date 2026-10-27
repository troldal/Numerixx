# Phase 4 API proposal: more bracketing methods (`Toms748`, `AndersonBjorck`)

Status: approved by the user on 2026-10-27, P1 and P2 as recommended. One review round (numerics-reviewer, simplicity-critic, usage-sketcher) is
folded in ("Decisions from the reviews"). Architecture: `docs/roots/BRIEF.md` (approved 2026-10-25), sections 3.2 and
3.8; Phase 1: `PHASE1_API.md` (validation, zero first, the stop test, `distance_bound`), which these methods share.

## Summary

```cpp
namespace rt = nxx::roots;

auto const kepler = [](double x) { return x - 0.9 * std::sin(x) - 1.0; };     // e = 0.9, M = 1
auto const a = rt::solve<rt::Brent>(kepler, rt::Bracket{ 0.0, std::numbers::pi });          // 10 evaluations
auto const b = rt::solve<rt::Toms748>(kepler, rt::Bracket{ 0.0, std::numbers::pi });        //  8
auto const c = rt::solve<rt::AndersonBjorck>(kepler, rt::Bracket{ 0.0, std::numbers::pi }); // 10

auto const ln = [](double x) { return std::log(x / 1e-12); };
auto const d = rt::solve<rt::Brent>(ln, rt::Bracket{ 1e-20, 1.0 });          // 45
auto const e = rt::solve<rt::Toms748>(ln, rt::Bracket{ 1e-20, 1.0 });        // 53
auto const f = rt::solve<rt::AndersonBjorck>(ln, rt::Bracket{ 1e-20, 1.0 }); // 18
```

- Two tags, no new API: `solve<Toms748>` and `solve<AndersonBjorck>` on a `Bracket`, with Phase 1's validation, zero
  first, tolerances, stop test, result and errors. A tuple fn is accepted and its first element used, as for Brent.
  The result is a point (the end with the smaller |fn|) and a bound, not the final bracket (Boost.Math's
  `toms748_solve` returns the bracket).
- **`Toms748`**: Alefeld, Potra and Shi's Algorithm 4.2 (1995) as printed (P1): inverse cubic interpolation with
  quadratic and secant fallbacks, a double-length secant step, and a bisection whenever an iteration has not halved
  the bracket.
- **`AndersonBjorck`**: regula falsi with Anderson and Bjorck's rescaling of the retained end (1973), and a bisection
  after three steps that have not halved the bracket (P2).
- **Neither is cheaper than Brent in general** (Prototype results): on 27 of the paper's own test cases at full
  precision Brent takes 309 evaluations, Toms748 325, AndersonBjorck 327 (the paper reports 4% fewer for Algorithm 4.2
  than Brent over its 154 cases, Table II). What they offer: a worst case of about 4 evaluations per halving, where
  Brent's is quadratic in the number of halvings; and AndersonBjorck is far cheaper where fn is strongly curved near
  one end (ln(x / 1e-12): 18 against 45; 1/x - 1: 5 against 11) and on a walk to a jump at 0 (51 against 1076).
  Brent stays the recommended default (BRIEF.md 3.2).

## Public API (`roots.hpp`)

```cpp
namespace nxx::roots {

    struct Toms748 {};        ///< Alefeld, Potra and Shi's Algorithm 748 (4.2). A bracket only; no derivative form.
    struct AndersonBjorck {}; ///< Regula falsi with Anderson and Bjorck's rescaling, guarded. A bracket only.

} // namespace nxx::roots
```

- `detail::bracketing` gains both: overload 1 (a bracket) accepts them; the guess and derivative overloads reject them
  by the existing constraints (`open`, `same_as<M, Newton>`).
- `solve` dispatches to `detail::toms748` and `detail::anderson_bjorck`, separate loops (BRIEF.md 3.8) on
  `Start<T, 1>` (no K parameter: one value), sharing one new helper, `detail::enclose` (B), and Phase 1's `prepare`,
  `distance_bound`, `midpoint`, `quotient`, `product` and `sum_fits`. Brent's loop is not touched (Newton and Halley
  depend on it bit for bit).
- The documentation of the tags and of `solve` gets a short guide to choosing (P3 below is not a decision; BRIEF.md
  decided it): Brent by default; Toms748 for its worst case of about 4 evaluations per halving; AndersonBjorck where fn
  is strongly curved near one end of the bracket, or for a walk to 0; for a value-only fn, never a derivative method.

## Behaviour

### A. Shared with Phase 1

Validation, zero first, `Criteria`, tol(x) = abs_tol + rel_tol |x|, `error` and `Root` are Phase 1's; the ends may
come in either order. The result is the end with the smaller |fn| (Brent's b) with `abs_error` the bound on the final
bracket's width, at most 2 tol(value), or tol(value) where fn is exactly 0. While the bracket is wider than max / 16,
every step is a bisection step, and no other step is formed: every difference of points is then at most max / 4.

### B. One step: `detail::enclose` (the paper's `bracket`)

Both methods take every point through one helper, the paper's `bracket` subroutine as modified for its tolerance
(Alefeld, Potra and Shi 1995, p. 341), on a state passed in and returned: the ends a < b with fn's values, and the end
it drops, d.
1. The point: the midpoint while the bracket is wider than max / 16; otherwise the method's point (computed only
   then).
2. **Placement:** delta = 0.7 tol(u), u the end with the smaller |f| (the paper: lambda = 0.7, tole = 2 |u| macheps +
   tol, Brent's), at least `denorm_min` (in `cpp_bin_float` and `cpp_dec_float`, `denorm_min()` is `min()`, and
   0.7 min flushes to 0), and max / 4 where tol(u) exceeds max / 2. Where the bracket is at most 4 delta wide
   (`width / 4 <= delta`, which cannot overflow) the point is the midpoint; otherwise it is kept at least 2 delta inside
   both ends. **Where the result is not strictly inside, the midpoint** (2 delta rounds away at a far end of much
   larger magnitude than u: fl(2.5 - 2 delta) = 2.5 for u = 0.3), and where the method has no point (a check failed).
   A point at least 2 delta inside either ends the search or removes at least 1.4 tol (Brent's minimum step is tol).
3. Evaluate; fn not finite is `non_finite`, fn exactly 0 is the root. Keep the half with the sign change; the dropped
   end becomes d.
4. The stop test of A on the new bracket.
A midpoint not strictly inside would mean the ends are adjacent numbers, whose distance is below 2 tol(u): the stop
test has fired before. It is asserted, and a release build returns the bracket with its bound (Phase 2's guard).

### C. `Toms748` (P1)

Algorithm 4.2 as printed (pp. 334-335), checked against the scan the first author publishes:
1. A secant step (4.2.1-4.2.2).
2. Each iteration n = 2, 3, ...: the inverse cubic point (`ipzero`, pp. 333-334) where fn's four values at a, b, d, e
   are distinct (not at n = 2) and the point is inside, otherwise `Newton-Quadratic` (p. 330) with 2 steps; e := d;
   `enclose` (4.2.3-4.2.4). Again with 3 steps as the fallback (4.2.5-4.2.6). A double-length secant step from the end
   with the smaller |f|, the midpoint where it would move more than half the bracket; e := d; `enclose`
   (4.2.7-4.2.10). If the bracket is not below mu = 1/2 of its width at the start of the iteration: e := d and a
   bisection step (4.2.11).
3. **Arithmetic:** every quotient, product and sum or difference that involves fn's values or the interpolations'
   terms is checked before it is formed (`quotient`, `product`, `sum_fits`): fb - fa, fd - fb, fe - fd, fe - fb,
   fd - fa, fe - fa; `Newton-Quadratic`'s B, A, fa + (B + A (c - b))(c - a), B + A ((c - a) + (c - b)) (not 2c - a -
   b, which overflows above max / 2) and c - step; `ipzero`'s nine terms, each a distance times a ratio of values,
   the ratio first, and their sums. A check that fails makes the point unusable, and the paper's fallbacks apply:
   inverse cubic, then quadratic, then secant, then the midpoint (B). `Newton-Quadratic` leaves for the secant as
   soon as an iterate leaves (a, b) (in exact arithmetic its iterates never do), which also bounds c - a and c - b.
   The secant's and the double-length secant's ratios are below 1 and 1/2 in magnitude (fa and fb have opposite
   signs), so their products with b - a cannot overflow.
4. **Worst case:** each iteration costs at most 4 evaluations and at least halves the bracket (the paper, p. 335), so
   with k the number of halvings bisection needs down to the final tolerance, at most 4 k + 1 evaluations after the
   ends; Brent's bound is (k + 1)^2 - 2.

### D. `AndersonBjorck` (P2)

Regula falsi on the bracket, keeping the latest point b and the retained end a, whose value is replaced by a scaled
value g_a:
1. The point is the zero of the line through (a, g_a) and (b, f_b), `b - (f_b / (f_b - g_a)) (b - a)`: the ratio is in
   (0, 1); f_b - g_a, of opposite signs, is checked; `enclose` (B).
2. If f_c has the sign of f_b, a is retained and g_a is multiplied by m = 1 - f_c / f_b, or by 1/2 where m is not
   positive (Illinois' factor; Anderson and Bjorck 1973, as stated by King 1983): formed as `1 - f_c / f_b` only where
   |f_c| < |f_b| (then 0 < f_c / f_b < 1, nothing overflows), else 1/2, which is the same rule. Otherwise b becomes
   the retained end, with its true value. c becomes b. g_a keeps its sign and only shrinks; where it underflows to 0
   the point is a, which B moves inside.
3. **The guard:** a bisection step after three consecutive steps that have not brought the bracket below half its
   width at the last halving; a guard bisection counts as a halving (its rounding cannot then cost a fifth
   evaluation), and its point updates a, b and g_a as any other. Without the guard the published method took
   65,109,847 evaluations on x^20 - 1 over [0, 5]: the latest point creeps up from 0 (5e-14, 1e-13, 2e-13, ...,
   0.007 after 40 evaluations) while the end at 5 is retained. At most 4 evaluations per halving with it.
4. Where the guard does not act (simple roots of smooth fn) the counts are the published method's; it acts on the
   multiple root (x - 1)^3 (152 against 126 unguarded) and on atan(x - 1)^3 (153 against 122), and makes x^20 - 1
   finite (30).

### E. Poles, jumps and multiprecision

- A pole inside the bracket is returned as a root (Phase 1); a point that lands exactly on it, where fn is infinite,
  is `non_finite` (Toms748 on 1/(x - 1/3) over [0, 1] lands on fl(1/3)). Check |fn(value)| where fn may have poles.
- A sign change near 0 that fn does not evaluate to 0 is approached down to T's smallest numbers (Phase 1's walk to
  0): about 1000 evaluations in double for Brent and Toms748, about 2e9 halvings in `cpp_bin_float_50`; set
  `abs_tol` there. AndersonBjorck's walk is geometric in the exponent (51 evaluations in double): its rescaled end
  makes each secant point a factor 2^-k closer to 0. The double counts do not carry over to other types.
- Constant evaluation of a walk of about 1000 steps may need a compiler's step limit raised; the tests use `abs_tol`.

## Prototype results

`.scratch/roots4/proto4.py` (Python, double, revision 2; every check of C and D with the C++'s semantics, on the
Phase 2 prototype's helpers). References from mpmath at 50 digits. Evaluations of fn, the ends and zero first
included. Every root is within its `abs_error` of the reference, with 4 eps |root| of slack for the libm rows only.

Cross-check: in Boost.Math's arrangement and under its placement and stop rule, the prototype takes Boost.Math
1.92.0's own counts (`toms748_solve`, BSL-1.0, built and run) on 7 of 8 cases, and one fewer on the eighth.

| Case | Bracket | Brent | Toms748 | AndersonBjorck |
|---|---|---|---|---|
| x^2 - 2 | [1, 2] | 9 | 9 | 9 |
| x^3 - 2x - 5 | [2, 3] | 8 | 9 | 9 |
| x^20 - 1 | [0, 5] | 19 | 21 | 30 |
| (x - 1)^3 | [0, 3] | 156 | 149 | 152 |
| x - 1e-300 | [0, 1] | 3 | 3 | 4 |
| x - 1 | [-max, max] / [0, max] | 8 / 7 | 8 / 7 | 13 / 12 |
| x^3 - 2x + 2 | [-2, -1] | 9 | 9 | 10 |
| 1/x - 1 | [0.5, 3] | 11 | 6 | 5 |
| W10 near 7 | [6.7, 7.6] | 23 | 18 | 18 |
| a jump at 0 | [-1, 2] | 1076 | 981 | 51 |
| 1/(x - 1/3), a pole | [0, 1] | 55 | `non_finite` (84) | 95 |
| sqrt x - 0.5 | [0, 1] | 4 | 6 | 9 |
| cos x - x | [0, 1] | 9 | 7 | 8 |
| e^x - 2 | [0, 1] | 8 | 8 | 8 |
| x e^x - 1 | [0, 1] | 9 | 9 | 9 |
| x e^-x - 0.01 | [0, 1] | 8 | 7 | 7 |
| Kepler, e = 0.9, M = 1 | [0, pi] | 10 | 8 | 10 |
| ln(x / 1e-12) | [1e-20, 1] | 45 | 53 | 18 |
| atan(x - 1)^3 | [0, 3] | 155 | 154 | 153 |

Alefeld, Potra and Shi's test problems (Table I, with its parameters; #3, #13, #14 and #15 contain 0 and zero first
settles or splits them), powers written as products:

| Problem | Brent | Toms748 | AndersonBjorck | Bisection |
|---|---|---|---|---|
| #1, sin x - x/2 | 9 | 10 | 9 | 52 |
| #2, n = 1 / 5 / 10 | 16 / 13 / 12 | 16 / 13 / 13 | 23 / 17 / 18 | 52 / 51 / 50 |
| #4, x^n - 1 on [0, 5], n = 4 / 8 / 12 | 15 / 15 / 17 | 13 / 16 / 18 | 19 / 22 / 23 | 55 |
| #5, sin x - 0.5 | 10 | 10 | 10 | 54 |
| #6, n = 5 / 20 / 100 | 11 / 13 / 15 | 10 / 11 / 14 | 10 / 13 / 16 | 55 / 57 / 60 |
| #7, n = 5 / 10 / 20 | 10 / 9 / 9 | 9 / 7 / 8 | 9 / 8 / 8 | 57 / 59 / 61 |
| #8, n = 5 / 10 / 20 | 9 / 11 / 13 | 11 / 11 / 12 | 10 / 11 / 13 | 54 / 55 / 55 |
| #9, n = 5 / 15 / 20 | 8 / 7 / 7 | 11 / 9 / 9 | 10 / 8 / 8 | 61 / 68 / 69 |
| #10, n = 5 / 10 / 20 | 8 / 10 / 13 | 8 / 11 / 13 | 9 / 10 / 14 | 53 |
| #11, n = 5 / 15 / 20 | 12 / 13 / 14 | 19 / 20 / 18 | 7 / 6 / 6 | 55 / 56 / 57 |
| #12, n = 3 | 10 | 5 | 10 | 58 |
| Total (27 cases) | 309 | 325 | 327 | 1520 |

- Boost.Math's arrangement of Algorithm 4.2 (P1's alternative) totals 321 on the 27 cases; on the first table it
  takes 19 on x^20 - 1, 151 on (x - 1)^3, 12 on 1/x - 1, 20 on W10, 9 on cos x - x, 85 on the pole.
- AndersonBjorck's guard after two steps: 19 on x^2 - 2 and 22 on x^3 - 2x - 5 (against 9); after four: 203 on
  (x - 1)^3 and 206 on atan(x - 1)^3 (against 152 and 153).

## Tests (`tests/test_roots_bracketing.cpp`, test cases `roots:...`)

- **Counts:** the prototype's, exactly, for fn of `+ - * /` and sqrt in double: the first table less the libm rows
  (each row commented with the rule it exercises: the wide bracket, the walk to 0, the guard on x^20 - 1, the pole),
  and problems #2, #4, #7, #8, #9 and #11; assuming no FMA contraction (Phase 2's assumption). The pole's outcome
  for Toms748 (`non_finite` after 84) is asserted as such.
- **Accuracy:** `|value - ref| <= abs_error + c` on the same rows, and on the libm rows with Phase 1's 2-ulp budget.
- **Edges:** fn exactly 0 at an interior point; a reversed bracket (the same `Root`); a tuple fn (the same result as
  the value fn, bit for bit); B's two fallbacks, each in a unit test of `enclose`'s placement with the inputs that
  reach it (delta flushing to 0 in `cpp_bin_float_50`; 2 delta rounding away at a far end).
- **Constant evaluation:** one `static_assert` per method in float, double and long double; the wide brackets in
  double.
- **Types:** float, long double, `cpp_bin_float_50` with and without expression templates (`T(7) / T(10)` and the
  interpolations' intermediates), one case per method in `cpp_dec_float_50` (radix 10).

## Documentation

- `roots.hpp`: the two tags (sources, cost against Brent in one line, the worst case, "a bracket only"); the guide to
  choosing in `solve`'s documentation; the file `@brief` and overload 1's `@tparam M`, which list the methods.
- `docs/roots/design.md`: the sources, B to D, the tables and the Boost.Math cross-check. `docs/roots/BRIEF.md` 3.2:
  Toms748's row corrected to Algorithm 4.2 (at most 4 evaluations per iteration, asymptotically 3, index 1.6686; the
  row had 4.1's figures). README: no change.
- `demo/roots_demo.cpp`: Toms748 and AndersonBjorck in section 2 (the comparison of methods).
- `docs/plan.md`: Phase 4's checkboxes.

## Changes during implementation

- The tests are in their own file, `tests/test_roots_toms_anderson.cpp` (as Phases 2 and 3), not in
  `test_roots_bracketing.cpp`. Every count of both tables for `+ - * /` and sqrt is reproduced exactly.
- The pole's `non_finite` for Toms748 is tested at run time only: fn itself divides by zero at fl(1/3), which is not a
  constant expression.
- Toms748's steps go through a small function, `detail::advance`, rather than a lambda capturing fn: MSVC rejects a
  lambda's capture of a function reference in constant evaluation (as at `Joined`).
- `demo/roots_demo.cpp` section 2 adds cos x - x and ln(x / 1e-12) for the two methods (7 / 8 and 53 / 18
  evaluations, Brent 9 and 45).

## Decisions for the user

- **P1. Toms748 is Algorithm 4.2 as printed.** Recommended: the method the name and the citation promise. Its first
  iteration has one interpolation step fewer than Boost.Math's arrangement, which takes a quadratic step before the
  loop and falls back from a failed inverse cubic to 3 Newton steps where 4.2.3 prescribes 2. Alternative: Boost.Math's
  arrangement (321 against 325 on the paper's problems; mixed on the first table), which would let the C++ be compared
  with Boost.Math's counts under Boost's rules as well.
- **P2. AndersonBjorck with a bisection guard after three steps** that have not halved the bracket. Recommended: the
  published method has no worst-case bound (65 million evaluations on x^20 - 1); three steps change nothing on simple
  roots of smooth fn. Alternatives: two steps (twice the evaluations on smooth fn) or four (a third more on multiple
  roots); or the method as published, documented as unguarded.

Decided by the brief or by Phase 1, and stated rather than put to the user: the placement of B (the paper's rule, with
our tolerance); the result as a point and a bound (`Root`'s contract); Brent as the default (BRIEF.md 3.2).

## Decisions from the reviews

From the numerics review:
- B's delta underflowed to 0 in `cpp_bin_float_50` (`denorm_min()` is `min()` there), and 2 delta could round away at a
  far end of larger magnitude; either left a point on an end, where Phase 2's guard would have stopped with the whole
  bracket. Now: delta at least `denorm_min`, and the midpoint wherever the placed point is not strictly inside (B).
- The overflow plan was incomplete: these methods, unlike Brent, form differences of fn's values of opposite signs;
  also 2c - a - b, 4 delta and m's quotient. Every checked operation is now listed (C.3, D), and no step is formed from
  a wide bracket (B.1).
- The prototype followed Boost.Math's arrangement, not Algorithm 4.2 as printed: now the paper's (P1), every count
  regenerated, Boost's kept as the cross-check and the alternative.
- Exact counts for problems #2 and #7 used Python's `**` (libm's pow): the problems' powers are now products, in the
  prototype and the tests alike; #4, #8, #9 and #11 join the exact counts.
- B's rationale was wrong (a point closer to an end can end the search); the worst case is stated with k defined; a
  guard bisection counts as a halving; the problems are Table I's with its parameters, and the excluded ones named;
  the stop test is on p. 340; BRIEF.md 3.2's row had Algorithm 4.1's figures; the accuracy slack applies to libm rows
  only.

From the simplicity review:
- One helper for the step both methods take (`enclose`), the loops separate; no K parameter on the new loops.
- The decisions cut to P1 and P2; the "Algorithm 4.1 variant" removed from the deferred list.
- Tests cut to what the new loops do: the guards are rows of the count table; zero at an end or at 0, and the
  constraints, are Phase 1's and 2's code, already tested.
- The cost figures in one place (`docs/roots/design.md`), a line in each tag's documentation.

From the usage sketch (6 snippets, 6 misuses):
- How to choose was scattered: a short guide in `solve`'s documentation, and the summary now shows a case each method
  wins and the totals.
- `Root` is a point and a bound, not a bracket, which a Boost.Math user may expect (and a structured binding of a
  `Root` gives value and abs_error): stated in the tags' documentation.
- The multiprecision walk to 0, and constant evaluation of long walks: stated (E).
- A reversed bracket: tested.

## Independent reviews (OpenAI Codex and Kimi K3.1 Pro, 8951e7c)

After the implementation (8951e7c, CI run 40511081791), two read-only reviews with the same brief: OpenAI Codex, and,
at the user's request, Kimi K3.1 Pro through the Cline CLI (plan mode) as a cold read, without Codex's findings.
Neither found a defect in `roots.hpp`: both checked Algorithm 4.2 against the paper, the arithmetic checks, the
placement and its fallbacks, termination, the stop test and `abs_error`, AndersonBjorck's rescaling and guard, the user
types and the constraints, and Codex ran extreme inputs in constant evaluation on GCC and Clang and 400 random cases
per method. Their findings were in the tests, and both are fixed:
1. **The constant-evaluation bounds for float and double had no derivation** (both): now the contract's bound,
   2 tol(value) + c with the reference's rounding, 8 eps for every type, derived in a comment.
2. **The no-subnormal test accepted an invalid bound** (Codex): it checked `abs_error <= 4 min` only; it now checks
   the value (0) and `abs_error <= 2 tol(0) = 2 min`.

## Sources

- Alefeld, G. E., Potra, F. A. and Shi, Y. (1995). Algorithm 748: Enclosing zeros of continuous functions. *ACM
  Transactions on Mathematical Software* 21(3), 327-344. Read from the scan on the first author's page
  (na.math.kit.edu): `bracket` and `Newton-Quadratic` p. 330, `ipzero` pp. 333-334, Algorithms 4.1 and 4.2
  pp. 334-335, the stop test p. 340, the modified `bracket` and Table I p. 341, Tables II and VI pp. 342-343.
- Anderson, N. and Bjorck, A. (1973). A new high order method of regula falsi type for computing a root of an
  equation. *BIT* 13, 253-264. Not read (behind the publisher's paywall), nor were King (1983) and Galdino (2011):
  the rescaling is as the Wikipedia article "Regula falsi" states it, citing them. The user may confirm it against the
  paper or Dahlquist and Bjorck's textbook.
- Boost.Math 1.92.0, `boost/math/tools/toms748_solve.hpp` (BSL-1.0): built and run as a cross-check; no code taken.

## Deferred

- `RegulaFalsi` and `Illinois` (Phase 6): D's step with m = 1 or m = 1/2 always; the factor becomes a parameter of the
  loop then, when it has a second use.
- Returning the final bracket in `Root` (BRIEF.md's deferred extra members).

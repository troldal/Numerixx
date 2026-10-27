# Phase 5 API proposal: searching for a bracket (`search`)

Status: approved by the user on 2026-10-27, P1 to P4 as recommended, with one change: the function is named
`search(fn, x0)`, not `expand(fn, x0)` (BRIEF.md's name), so that another algorithm that does not expand can share
the name later (a tag, `search<M>(fn, x0)`, with the outward search as the default: a defaulted leading template
parameter, non-breaking for every `search(fn, x0)` call). One review round (numerics-reviewer, simplicity-critic, usage-sketcher) is
folded in ("Decisions from the reviews"). Architecture: `docs/roots/BRIEF.md` (approved 2026-10-25), sections 3.1
and 3.9 and decision D7 (searching as one concrete function, `expand(fn, x0)` there); Phase 1: `PHASE1_API.md` (`Bracket`, `error`).

## Summary

```cpp
namespace rt = nxx::roots;

// A bracket from a guess, then a bracketing method on it.
auto const f = [](double t) { return std::exp(20.0 - 5000.0 / t) - 1.0; };   // a vapour-pressure-like curve
auto const r = rt::search(f, 300.0)                                          // 14 evaluations: [231.7, 254.4]
                   .and_then([&](rt::Bracket<double> b) { return rt::solve<rt::Brent>(f, b); }); // 8 more: 250

// Outside fn's domain, return NaN: the search then approaches the domain's edge instead of crossing it.
auto const g = [](double v) {    // van der Waals, CO2-like, 280 K, 4 MPa; a pole at the co-volume b = 4.27e-5
    if (!(v > 4.27e-5)) return std::numeric_limits<double>::quiet_NaN();
    return 8.314 * 280.0 / (v - 4.27e-5) - 0.364 / (v * v) - 4e6;
};
auto const v = rt::search(g, 1e-4); // [3.6e-4, 4.9e-4] in 29 evaluations; without the NaN, a bracket around the pole
```

- **`search(fn, x0)`**: evaluates fn at x0, then at x0 - h and x0 + h for growing h, until fn's sign differs from its
  sign at x0; returns the narrowest bracket that shows it, its ends in order.
- **The rules** (P1, P2): a first step |x0| / 50 (1/50 at 0), growing by 3/2 per step; at a non-finite value of fn,
  that side approaches fn's domain edge by bisection, at most `digits` probes (53 in double), instead of stopping
  (which changes BRIEF.md 3.9's rule); no bound on the search region (P3: return NaN outside fn's domain).
- **Errors:** `non_finite` (x0 or fn(x0) not finite) and `no_sign_change` (none found: P4). Composes with `solve`
  through `and_then`; a `no_sign_change` from the composition is always search's, since its bracket always has one.
- **What it cannot find** (documented): a root of even multiplicity ((x - 1)^2); two roots within one step (x^2 - 2
  from 100: the points 23.1 and -15.3 straddle both); a root beyond a region where fn is not finite on the same side;
  a sign change beyond the cap (x - 1 from 1e-300). **A pole brackets like a root** (tan x; van der Waals' b without
  the NaN), and `solve` then returns the pole: check |fn(value)|. The bracket need not hold the root nearest x0, and
  where several roots lie in it (a cubic equation of state inside its two-phase region), `solve` returns one of them.

## Public API (`roots.hpp`)

```cpp
namespace nxx::roots {

    /**
     * A bracket around a sign change of fn, searched outward from x0: fn at x0, then at x0 - h and x0 + h for
     * h = |x0| / 50 (1/50 at x0 = 0) growing by 3/2.
     */
    template<real T, typename F>
        requires detail::gives<F, T, 1>
    [[nodiscard]] constexpr std::expected<Bracket<T>, error> search(F &&fn, T const &x0);

} // namespace nxx::roots
```

- T comes from x0 (`search(f, 0)` is T = int, rejected: write 0.0). A tuple fn is accepted and its first element used
  (`gives<F, T, 1>`, as the bracketing methods), so `search(fdf, x0)` composes with `solve<Newton>(fdf, b)`.
- The result's ends are in order, lower <= upper; equal (both x) where fn is exactly 0 at x.
- fn's values at the ends are not returned: `solve` evaluates them again (2 evaluations). The return type stays
  `Bracket<T>`; a saving would come through a composed entry point (`bracket<M>(fn, x0)`, deferred), not by changing
  it.
- `error::no_sign_change`'s documentation is extended to `search`: no sign change between the points evaluated,
  within the cap or before both sides ended.

## Behaviour

1. **Validation:** x0 finite, else `non_finite`. fn(x0) not finite: `non_finite` (the sign at x0 is the reference).
   fn(x0) exactly 0: the bracket {x0, x0}, which `solve` returns at once as a root.
2. **Steps:** h_0 = |x0| / T(50), or T(1) / T(50) where x0 = 0, and at least `denorm_min` (|x0| / 50 underflows to 0
   for |x0| below 25 denorm_min). Each step tries x0 - h, then x0 + h, checking the sign after each evaluation. A side
   whose point would overflow (`sum_fits(x0, -+h)` fails) ends. h grows by 3/2 after a step while h <= max / (3/2);
   past that no larger step exists, and the search ends.
3. **A sign change** (fn at the new point of the other sign than fn(x0)) returns [the side's last point, the new
   point], in order: fn has x0's sign at the former, since every earlier finite point of that side had it. fn exactly
   0 at the new point returns {x, x}.
4. **Domain edges (P2):** a non-finite value of fn makes that side bisect (Phase 1's `midpoint`) between its last
   finite point and the failing one, checking the sign at each finite midpoint, instead of growing; the side ends
   after `digits` probes (53 in double, 168 in `cpp_bin_float_50`), or where the midpoint is no longer strictly between
   the two. Each probe is that side's evaluation for the step.
5. **The cap:** at most `iteration_cap<T>()` steps (4 times the bits of T: 212 in double, 672 in `cpp_bin_float_50`),
   reaching h_0 1.5^211 = 1.4e37 h_0; so at most 2 cap + 1 evaluations (425 in double, 1345 in `cpp_bin_float_50`),
   then `no_sign_change`. A search to max in `cpp_bin_float_50`, whose exponents reach 2^31, would take about 3.7e9
   steps.
6. **Constant evaluation:** constexpr; x0 +- h is formed only where `sum_fits` admits it, h grows only below
   max / (3/2), the bisection is `midpoint`'s. Constants are `T(3) / T(2)` and `T(50)`; h_0 = |x0| / T(50) (not
   |x0| (1/50), whose rounding differs).

### Why these rules

- **First step relative to x0 (P1):** x0 sets the scale, so a quantity in its own units searches at its own scale:
  van der Waals from v0 = 1e-3 m^3/mol finds the root at 4.4e-4 in 20 evaluations, where max(|x0|, 1) / 50 = 0.02
  (`nxx::deriv`'s scale convention) steps to v < 0 at once and finds nothing in 425. The cost: from a tiny x0 a root
  at scale 1 lies beyond the cap (x - 1 from 1e-300), documented: start at 0 or at the expected scale.
- **Growth 3/2 (P1):** small enough to keep a root and a nearby pole apart where 2 does not (1/x - 1 from 3: with
  growth 2 the points 1.08 and -0.84 straddle the root at 1 and the pole at 0, no sign change; found with 3/2 and
  sqrt 2), large enough that a far root costs little (x - 1e6 from 1: 91 evaluations; sqrt 2: 107; 2: 55). Exact in
  every radix, where sqrt 2 needs a square root, which `nxx::real` does not offer.
- **Domain edges approached, a bounded number of probes (P2):** stopping a side at a non-finite value (BRIEF.md 3.9;
  MATLAB's `fzero` ends its search with exit flag -3 when it meets a NaN or Inf value) never finds a root near a domain edge: ln(x / 1e-12)
  from 1 steps from 0.23 to -0.15 (NaN) and stops; probing finds [8.7e-13, 3.7e-12] in 96 evaluations, sqrt x - 1e-3
  from 1 likewise in 60. The cost, where the edge lies on the side away from the root: one probe per step until the
  limit (ln x - 30 from 1: 150 evaluations; stopping: 97; unlimited probing: 171).
- **No search bound (P3):** where fn is finite outside its domain, the caller returns NaN there, and the search
  approaches the edge: van der Waals with NaN at v <= b finds the root from 1e-4 and 5e-5 in 29 and 33 evaluations
  (a bound at b would take 25 and 23; without either, the search brackets the pole at b). The library gives the same
  advice for numerical derivatives (`roots.hpp`, Newton's documentation).

## Prototype results

`.scratch/roots5/proto5d.py` (Python, double), with Phase 1's Brent on the bracket found. Evaluations of fn:
`search` + `solve<Brent>` (which evaluates the ends again).

| Case | x0 | search + Brent | Bracket | Root |
|---|---|---|---|---|
| x^2 - 2 | 1 / 0 | 19 + 8 / 24 + 8 | [1.342, 1.513] / [-1.730, -1.153] | 1.414 / -1.414 |
| x^3 - 2x - 5 | 0 | 27 + 9 | [1.730, 2.595] | 2.095 |
| 1/x - 1 | 3 | 20 + 9 | [0.693, 1.462] | 1 |
| x - 1e6 | 1 | 91 + 3 | [7.46e5, 1.12e6] | 1e6 |
| x - 1e-300 | 1 | 22 + 4 | [-0.153, 0.231] | 1e-300 |
| x + 1e-300 | denorm_min | 266 + 3 | [-1.08e-300, -7.2e-301] | -1e-300 |
| x - 1 | max / 2 | 22 + 5 | [-1.38e307, 2.08e307] | 1 |
| sqrt x - 1e-3 | 1 | 60 + 4 | [2.97e-7, 1.03e-6] | 1e-6 |
| van der Waals | 1e-3 | 20 + 8 | [2.31e-4, 4.87e-4] | 4.37e-4 |
| van der Waals, NaN at v <= b | 1e-4 / 5e-5 | 29 + 8 / 33 + 8 | [3.59e-4, 4.89e-4] / [3.42e-4, 4.88e-4] | 4.37e-4 |
| van der Waals, finite below b | 1e-4 / 5e-5 | 20 + 67 / 12 + 61 | around the pole at b | b / `non_finite` (Brent lands on b) |
| cos x - x | 0 | 21 + 7 | [0.513, 0.769] | 0.739 |
| Kepler, e = 0.9, M = 1 | M | 23 + 7 | [1.769, 2.153] | 1.862 |
| e^x - 1e10 | 0 | 39 + 14 | [19.7, 29.6] | 23.03 |
| atan x - 1.5 | 0 | 37 + 8 | [13.1, 19.7] | 14.10 |
| e^(20 - 5000 / T) - 1 | 300 | 14 + 8 | [231.7, 254.4] | 250 |
| ln(x / 1e-12) | 1 / 1e-13 | 96 + 8 / 35 + 7 | [8.7e-13, 3.7e-12] / [9.8e-13, 1.4e-12] | 1e-12 |
| ln x - 30 (edge away from the root) | 1 | 150 + 7 | [8.25e12, 1.24e13] | 1.07e13 |
| tan x | 1.4 | 13 + 61 | [1.542, 1.613] | pi / 2, the pole |
| x^2 - 2 (two roots in a step) | 100 | 425, `no_sign_change` | | |
| x^2 + 1 (no root); (x - 1)^2 (even) | 1; 0 | 425, `no_sign_change` | | |
| x - 1 (beyond the cap) | 1e-300 | 425, `no_sign_change` | | |

Alternatives on the same cases: growth sqrt 2 finds the same roots at up to 16 evaluations more (2 fewer on
ln(x / 1e-12) from 1); growth 2 is cheaper but misses 1/x - 1 from 3; a first step max(|x0|, 1) / 50 finds x - 1 from
1e-300 but not van der Waals from 1e-3; stopping a side at a non-finite value misses ln(x / 1e-12) and sqrt x - 1e-3
from 1.

## Tests (`tests/test_roots_search.cpp`, test cases `roots:search:...`)

- **Counts and brackets of `search`, exactly**, for the `+ - * /` and sqrt rows (x^2 - 2, the cubic, 1/x - 1,
  x - 1e6, x -+ 1e-300, x - 1 from max / 2, sqrt, van der Waals with and without the NaN), in double, the test
  functions written so that they cannot contract to FMA (as Phases 2 to 4 assume); one `no_sign_change` case at 425.
- **ln(x / 1e-12) from 1** (a domain edge through libm's -inf and NaN): the bracket contains 1e-12 and fn's signs at
  its ends differ.
- **Composition:** `search(...).and_then(solve<Brent>)` gives the root, and `no_sign_change` passes through.
- **Edges:** fn(x0) = 0 ({x0, x0}, and `solve` on it); fn(x0) and x0 not finite; x0 = denorm_min; x0 = +-max and
  max / 2 (no overflow); a tuple fn.
- **Constant evaluation:** a `static_assert` per type (float, double, long double), and the domain edge in double with
  an explicit `quiet_NaN()` (sqrt or log of a negative number is not a constant expression).
- **Types:** `cpp_bin_float_50` with expression templates and `cpp_dec_float_50`: a bracket holding the root; the cap
  in `cpp_bin_float_50` (x^2 + 1: 1345 evaluations).

## Documentation

- `roots.hpp`: `search`, with the rules, the composition, what it cannot find (poles first, with the |fn(value)| check
  and the NaN idiom for a domain), `no_sign_change` extended; the file comment.
- `docs/roots/design.md`: P1 to P4, the table; that P2 changes BRIEF.md 3.9's rule. README: one line.
- `demo/roots_demo.cpp`: a section on searching (a bracket from a guess, a domain edge, a pole, a miss).
- `docs/plan.md`: Phase 5's checkboxes.

## Changes during implementation

- The probe limit is the bits of T (`iteration_cap<T>() / 4`), not `std::numeric_limits<T>::digits`: they agree in
  binary (53 in double, 168 in `cpp_bin_float_50`), but for a decimal type `digits` counts decimal digits, which would
  give it about a quarter of the halvings.
- From x0 near max the search reaches down to about max / 3 (no step exceeds max / (3/2)): x - 1 from max is
  `no_sign_change`, x - max / 2 from max is found. Documented with `search` and tested.
- The tests use the test cases `roots:search:*` in `tests/test_roots_search.cpp`; every count and bracket end of the
  `+ - * /` and sqrt rows is reproduced exactly, the ends written as hexadecimal literals.
- `demo/roots_demo.cpp` section 12: a bracket and its root, a domain edge, van der Waals with and without the NaN, and
  a miss.

## Decisions for the user

- **P1. A first step |x0| / 50 (1/50 at 0, at least denorm_min) and growth 3/2.** Recommended: x0 sets the scale; 3/2
  is exact in every type and keeps a root and a pole apart where 2 does not. Alternatives: max(|x0|, 1) / 50 (assumes
  a scale of 1; misses van der Waals from 1e-3); sqrt 2 (needs a square root); 2 (cheaper, misses 1/x - 1 from 3).
- **P2. A side approaches a domain edge by bisection, at most `digits` probes**, rather than stopping. Recommended:
  finds roots near an edge (logarithms, square roots, volumes above a co-volume); the limit bounds the cost where the
  edge is on the far side (150 against 97 on ln x - 30). Alternatives: BRIEF.md 3.9's rule, the side stops; unlimited
  probing (171).
- **P3. No search bound; outside fn's domain, fn returns NaN.** Recommended: the NaN form finds what a bound finds at a
  few evaluations more (29 and 33 against 25 and 23), states an open edge (v > b) exactly where a bound needs a fudge
  (b (1 + 1e-9)), and leaves `Bracket`'s meaning (an interval with a sign change) alone. Alternative: an optional
  `within` region (proposed first), with infinite ends meaning unbounded and its own type rather than `Bracket`; it
  also gives a one-sided search (x^2 - 2 from 1 above 1: 10 evaluations against 19). Deferred, non-breaking: a
  trailing defaulted parameter.
- **P4. No sign change found is `no_sign_change`.** Recommended: the same meaning as for a bracket, and the caller's
  remedy is the same; MATLAB's `fzero` reports the same case separately (exit flag -6). Alternative: `not_converged`
  (BRIEF.md's first assignment), whose meaning is an open method's failure.

## Decisions from the reviews

From the numerics review:
- The guard of x0 + h against the distance to an end used `distance_bound`, an upper bound: h within an ulp or two of
  it could form a point beyond the end (an overflow at max). With no `within` (P3) the only ends are +-max, and the
  check is `sum_fits` (2).
- h_0 = |x0| / 50 underflows to 0 for |x0| below 25 denorm_min, and the search evaluated x0 until the cap: now at
  least denorm_min (x + 1e-300 from denorm_min: 266 evaluations).
- P2's "no cost elsewhere" was false: an edge on the far side costs a probe per step (171 against 97 on ln x - 30
  from 1); now limited to `digits` probes (150), and the cost stated.
- Two examples used growth 2's points (36 and -28; 0.36 and -0.28): corrected to 3/2's. The reach is 1.5^211 h_0; a
  search to max in `cpp_bin_float_50` would take 3.7e9 steps; `no_sign_change` also when both sides ended; Brent, not
  search, lands on van der Waals' pole; h_0 is |x0| / T(50), so that the exact counts hold; the test functions avoid
  FMA contraction; the constant-evaluation domain edge returns NaN explicitly.
- MATLAB: the first step and sqrt 2 are not verified (MathWorks' code is proprietary); its public documentation
  lists exit flag -3 for a NaN or Inf met while searching and -6 for no sign change, cited now (checked against the
  page).

From the simplicity review:
- `within` cut and deferred (P3): its domain use is met by the NaN form, its one-sided use only saves evaluations; it
  would have reused `Bracket` for a region without a sign change, needed a fudge at an open edge, a precondition and
  its release-build behaviour, and six tests.
- The overflow rule is `sum_fits`, not a new distance comparison.
- The deferred `solve<M>(fn, x0)` for a bracketing M is struck: `solve<Newton>(fn, x0)` already means Newton from a
  guess; only `bracket<M>(fn, x0)` remains.
- Tests: `search`'s counts only (Brent's are Phases 1 and 4's); most libm rows to the demo; one cap case per type.

From the usage sketch (8 snippets, 7 misuses):
- x0 outside `within` arose from data (v0 = RT / P below b at high pressure), and infinite ends were the natural
  spelling of a half-line: both gone with `within` (P3).
- The result's order is now guaranteed (lower <= upper).
- Poles lead the documentation of what it cannot find, with the |fn(value)| check and the NaN idiom; several roots in
  one bracket (a cubic equation of state) documented.
- `search(f, 0)` deduces int: documented.

## Sources

- MATLAB's `fzero` documentation (MathWorks, public): from a scalar start it searches intervals around x0 that grow
  until fn's sign changes; a NaN or Inf value met while searching ends it with exit flag -3, and no sign change is
  exit flag -6 (both quoted from the page, 2026-10-27). Its first
  step and growth factor (|x0| / 50 and sqrt 2 as commonly described, BRIEF.md 3.9) are not verified: MathWorks' code
  is proprietary and not read. The rules here rest on the measurements above.
- R's `uniroot` (`extendInt`: "yes", "upX", "downX", extending a given interval) and Boost.Math's
  `bracket_and_solve_root` (multiplicative growth in a direction the caller states), from their public documentation.
  Not used as code.
- Not used: Numerical Recipes' `zbrac` (BRIEF.md: from these descriptions and first principles only).

## Deferred

- A search bound (`within`) or a direction (one-sided search): a trailing defaulted parameter, non-breaking; its own
  type rather than `Bracket` (P3).
- `bracket<M>(fn, x0)`, composing `search` and `solve<M>` and reusing the ends' values (BRIEF.md 3.1).
- A first step or growth chosen by the caller.
- Every root in an interval (BRIEF.md 3.9: a grid scan, later Chebyshev proxy root-finding).

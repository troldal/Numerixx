# Phase 2 API proposal: adaptive Gauss-Kronrod, the default method

Status: approved by the user on 2026-10-31, Q1 to Q12 as recommended (draft 2). One review round (numerics-reviewer,
simplicity-critic, usage-sketcher), with the prototype rerun on its questions, is folded in ("Decisions from the
reviews"). Architecture: `docs/quadrature/BRIEF.md` (approved 2026-10-30); Phase 1: `docs/quadrature/PHASE1_API.md`
(done, 930bb49).

## Summary

```cpp
#include <numerixx/quadrature.hpp>
namespace qd = nxx::quadrature;

auto const f = [](double x) { return std::exp(x); };

// The default method: adaptive G7-K15, 10 panels to start (150 evaluations), more where the estimate is largest.
auto const r = qd::integrate(f, 0.0, 1.0);                                  // e - 1; abs_error an estimate
auto const s = qd::integrate(f, 0.0, 1.0, {.rel_tol = 1e-14});
auto const t = qd::integrate<qd::GaussKronrod<15>>(f, 0.0, 1.0);            // the same, named

// A kink, a mild endpoint singularity: subdivision finds them.
auto const k = qd::integrate([](double x) { return std::abs(x - 0.3); }, 0.0, 1.0);
auto const u = qd::integrate([](double x) { return std::log(x); }, 0.0, 1.0);  // -1

// Failure with the best estimate: fn noisier than the request (stagnation), or the cap.
if (!r && r.error().best) report(r.error().best->value, r.error().best->abs_error);
```

- `GaussKronrod<15>`, the default method of the estimating overload (BRIEF.md D1, D7).
- The adaptive loop of BRIEF.md 3.6, with the parameters a prototype fixed and five corrections to the brief's
  design that the prototype and the reviews found (B, D, E).
- A node cache for types without a compile-time table, for every rule (Phase 1's Gauss-Legendre included).
- The README example and a demo program, which Phase 1 deferred until the default method existed.
- Everything constexpr; nothing throws.

## Public API (additions to `quadrature.hpp`)

```cpp
namespace nxx::quadrature {

    /// The G7-K15 Gauss-Kronrod pair, globally adaptive: the default method. Only N = 15 exists (Q7, Q11).
    template<int N>
        requires(N == 15)
    struct GaussKronrod
    {
    };

    /// A method with an estimate: GaussKronrod<15> (the default) or Romberg.
    template<typename M = GaussKronrod<15>, real T, typename F>
        requires detail::estimating<M> && nxx::detail::callable<F, T>
    [[nodiscard]] constexpr std::expected<Integral<T>, Failure<T>>
    integrate(F &&fn, T const &a, T const &b, Criteria<std::type_identity_t<T>> const &criteria = {});
}
```

- **The second overload gains the default method;** its constraint becomes `detail::estimating<M>` (Romberg or
  `GaussKronrod<15>`), the concept Phase 1 deferred until a second method existed. A default template argument breaks
  no existing call.
- **`GaussKronrod` is constrained on its N** (Q11), so `GaussKronrod<21>` fails where it is named, with the
  constraint in the message, instead of an overload set that does not match. Removing the constraint when K21 comes
  breaks nobody.
- **Comments:** `error::not_converged` gains Gauss-Kronrod's three causes (the cap, stagnation, final panels whose
  errors exceed the tolerance), each with `Failure::best`; `error::non_finite` says that an infinite end means
  `TanhSinh` (Phase 4) and that fn not finite on a region no short bisection avoids is the other cause; `integrate`'s
  `@return` says that a result is returned when the request is met or the estimate is at the floor, so a caller who
  needs the request reads `abs_error`; `Criteria`'s comment says that Gauss-Kronrod accepts at twice the floor
  (Romberg at the floor), and that `Criteria` and `error` collide with `nxx::roots`' under two using-directives, with
  different default tolerances.
- Nothing else in the public API changes. `Integral`'s extra members stay deferred (BRIEF.md D12).

## Behaviour

### A. Validation and the map

As Phase 1 (PHASE1_API.md A, B, and its independent review's corrections): tolerances, then the ends (NaN or
infinite: `non_finite`), a == b is `{0, 0}` without calling fn, b < a negates; `detail::Interval` with its scale
order and the clamp of every point into [min(a, b), max(a, b)]. Equal adjacent mesh points (an interval too narrow
for 10 panels) give empty panels, which are skipped. In constant evaluation, a value or product the Phase 1 checks
refuse (`checked_term`) is `non_finite`, as in Phase 1: the descent of D is a run-time path only, so a
`static_assert` and its run-time twin agree wherever both succeed.

### B. The rule

- **G7-K15** on each panel: Phase 1's 7 Gauss-Legendre nodes and weights, and 8 Kronrod nodes, the zeros of the
  Stieltjes polynomial E_8, interlaced with them. E_8 as a Legendre series whose coefficients come from the triangular
  orthogonality system (Adams' closed form for the integral of three Legendre polynomials; with n = 7 odd, only odd j
  give conditions, and condition j fixes coefficient n - j); its zeros by Phase 1's safeguarded Newton, each between
  two Gauss nodes or between the largest and 1, the bracket's orientation from the node's index. Kronrod weights from
  the interpolatory formulas: 2 / (8 P_7 E_8') at a new node, and w_G + 2 / (8 P_7' E_8) at a Gauss node.
- **The estimate's difference** `D = K - G` from two Neumaier sums over the same evaluations (Q10). Correction
  (simplicity review, prototype R4): BRIEF.md 3.6 proposed difference weights with their sum forced to 0 through the
  centre; on the whole battery they changed no status and one count by 30 evaluations, and the largest
  |D| / (eps L1_p) on narrow panels of f = 1, e^x and 1 / (1 + x^2) is 1.28, 1.62 and 1.52 with K - G (0.47 to 0.68
  forced), 30 times below the floor's 50 in double, 11 times in `cpp_bin_float_50` with Phase 1's plain-T weights.
- **Accuracy of the generated rule** (spike, against Boost.Math's 115-digit table): nodes within 0.44 eps absolute in
  every type; Kronrod weights within 0.89 eps absolute; sum |dw_K| at most 6.1 eps.
- **Cost:** at most 12,079 constexpr steps for the table (em++, binary128 long double), 1.2 % of the tightest default
  limit. Tables for float, double and long double, as in Phase 1. The phase measures the steps per bisection, so
  that the documentation can say how far a compile-time integral reaches.
- **The node cache** (Q6, BRIEF.md D10): for every T without a table, each rule (`GaussLegendre<N>` and K15) is a
  function-local `static const`, built on first use on the path not taken in constant evaluation (`if consteval
  { build } else { static ... }`). Initialised once, thread-safely, by the language; never mutated. The spike:
  compiles without warnings on all five toolchains; e^x on [0, 1] with one panel in `cpp_bin_float_50` 0.50 ms a call
  without, 0.095 ms with. Safe because no `nxx::real` type can change its precision at run time: Boost.Multiprecision
  has no header-only variable-precision type, and its MPFR and GMP ones do not specialise `numeric_limits`. A user
  type that declares a constant `digits` and changes its precision at run time would see stale nodes: documented.

### C. Panels, the heap and the sums

- **The initial mesh: 10 equal panels** (Q1, as Shampine's `quadgk`), so 150 evaluations before any bisection.
  The prototype: false successes on 100 Gaussian peaks of width 1e-3 at the default tolerance were 36 with 1 panel,
  2 with 2, 1 with 4, none with 10; Kahaner's 21st problem was missed with 1 or 2. No mesh cures a feature narrower
  than the nodes' spacing (documented).
- **A panel's estimate** e_p = max(|D_p|, floor_p), floor_p = 50 max(eps L1_p, denorm_min), L1_p = |hw| sum
  w_K,i |f_i| (Phase 1's floor).
- **A panel is final when it is created** if |D_p| <= floor_p, or if it cannot be split: an outer node of either
  child would not lie strictly inside that child (which also covers a midpoint that is not strictly inside). A final
  panel never enters the heap; it adds to the running sums. Finality at the floor cannot stop early: splitting
  cannot bring a panel's estimate below its floor, since L1 adds up (numerics review).
- **The heap** (`std::vector`, `std::push_heap`, transient allocation, so constexpr) holds the active panels ordered
  by |D_p| (the same order as e_p, since only panels above their floor enter it), grows on demand (no reserve: a
  reserve from the cap would allocate about 1.7 MB per call in `cpp_bin_float_50`), and the largest is bisected.
- **Three running sums, Neumaier-compensated** (simplicity review): the value K, the floor F, and the excess
  X = sum (e_p - floor_p) over all panels, active and final, each updated by subtracting a parent's contribution
  and adding its children's; the estimate is E = F + X. Panels at their floor add exactly 0 to X, and the same bits
  are subtracted that were added, so X does not drift. The result's value is the compensated running K.

### D. Values that are not finite (BRIEF.md D14, corrected twice)

- One finiteness check per panel, on L1_p.
- **A splittable panel whose L1_p is not finite is bisected,** iteratively, **at most 2 levels below the panel that
  started the descent** (Q8), each split counted against the cap. A singularity at a panel's centre node becomes an
  end of both children: log |x| and 1 / sqrt |x| on [-1, 3], where the centre of the mesh panel [p2, p3] is exactly
  0 in double, take one split (1,320 and 2,550 evaluations).
- **A descent that ends in a piece still not finite is undone:** the panel it bisected becomes final with its own
  last finite estimate (its parent's, for a child created by the descent), and **the run can no longer end ok**,
  only `not_converged` with best (Q8). A panel of the initial mesh with no finite estimate gives `non_finite`.
- Corrections: draft 1's descent was unbounded (on sqrt(x - 0.5) over [0, 1] it walked to 0: 1,063 levels and
  101,775 evaluations in double, about 2^31 levels in `cpp_bin_float_50`; bounded: `non_finite` in 105
  evaluations); the brief's D14 ended `non_finite` with no best where fn overflows near a singular end (x^-0.99 on
  [0, 1]); now `not_converged` with best 99.924 +- 1.4e-3 (the true error is 0.076: the end panel is beyond the
  documented x^-0.65 limit). The rule that a run with an undone piece cannot end ok changed no status in any case
  the prototype ran (the battery, 60 overflowing peaks, 40 NaN regions, 1,200 runs with NaN micro-regions), but a
  dropped piece with a small share of E could otherwise let the rest converge; it costs nothing.

### E. Termination

- **The stop test** `E <= max(abs_tol, rel_tol |K|, 2 F)` (Q4), that is, X <= max(abs_tol, rel_tol |K|, 2 F) - F.
  Correction (prototype): the brief's test with F needs X = 0, every panel at its floor, which a kink or a jump never
  reaches (its |D_p| / L1_p does not fall with the width). At rel_tol = 0 the brief's test ended 32 of 50 kinks and
  42 of 50 jumps `not_converged`; with 2 F, 50 and 38 of them end ok with 4 times fewer evaluations. On a smooth fn
  the 2 F test cannot accept an error above max(request, 2 F) (numerics review).
- **The heap empty** (every panel final): ok if the test passes, else `not_converged` with best.
- **The cap** (Q3): 64 b bisections and descent splits together, b = digits log2(radix) (3,392 in double, 1,536 in
  float, about 10,750 in `cpp_bin_float_50`), then `not_converged` with best. A change to BRIEF.md 3.6, which said
  "from T's exponent range and digits": `cpp_bin_float`'s exponent range is about 2^31, so a cap from it does not
  bound anything. The binding case is x^-0.9 on [0, 1] at `cpp_bin_float_50`'s default tolerance, 3,472 bisections.
- **Stagnation** (Q2): count the bisections since X last reached a new minimum; at 32, `not_converged` with best.
  Correction: the brief's rule watched E, which a peak found late raises and which barely moves at the floor; draft 1
  watched X with a reset when X grew tenfold, which never tripped on a divergent integral (1/x^2 on [0, 1] in 168 bits
  ran to the cap: 322,710 evaluations). Without the reset (prototype R3, against two variants): no false trip over the
  battery at four tolerances and five combs of 40 peaks; 1/x^2 and 1/(x - 0.5)^2 stop at 1,140 evaluations in double
  and 168 bits; e^x (1 + 1e-10 h(x)), h a hash of x's bits (noise that never resolves), stops at 1,620 with an honest
  best (1.2e-11 against an actual 1.7e-13), where it ran to the cap without a rule.
- **The result:** `{K, E}`; E >= F, so abs_error is never below the floor.

### F. What the estimate is, and the documented limits

Measured by the prototype on 710 cases at four tolerances (Kahaner's 21, 200 Lyness-Kaganove peaks, 350 rational and
Gaussian peaks, 21 powers, 50 kinks, 50 jumps, edge cases, noisy fn; references from closed forms or mpmath at 60
digits; the tables go into `docs/quadrature/design.md`): Kahaner's 21 all converge at the default with no estimate
below its error; false successes (ok, the error above both the tolerance and the estimate), excluding Gaussian peaks
of width 1e-4, are 1.8 % at 1e-6, 2.0 % at the default, 2.3 % at 1e-12, 4.7 % at 0; the median estimate is 10^5.6
times the error at the default (|K - G| estimates G7's error, not K15's).

Documented in the header and on the site:
- a feature narrower than the nodes' spacing is missed (most Gaussian peaks of width 1e-4 on [0, 1]: every sample
  is 0);
- each panel's outer 0.43 % at each end, 0.85 % around each boundary between panels, is sampled by neither rule, so a
  kink or a step there is invisible (a step at 0.599763 lies between the outer node 0.599573 of [0.5, 0.6] and its
  end: estimate 4.4e-15, error 2.4e-4);
- **an endpoint singularity stronger than about x^-0.65:** |K - G| understates the end panel's error by a ratio that
  bisection does not change (1.34 at alpha = -0.7, 4.92 at -0.9, 53.7 at -0.99). Until `TanhSinh` (Phase 3), the
  documentation gives the substitution that removes it, x = a + u^k for (x - a)^alpha, with k (1 + alpha) >= 1, and
  u = b - x for the other end, with a tested example (Q9);
- fn's noise (as for Romberg), and a nested integral's inner tolerance tighter than the outer one;
- full precision in a multiprecision type costs: at 1e-45 in 50 digits K15 needs 2,250 (e^x) to 88,440 (a peak)
  evaluations, 7 to 45 times K31's, because |K - G| converges at G7's order; at the default (1e-25) it is within 1 to
  5 times K31 (Q7);
- the minimum of 150 evaluations: for an expensive fn that must cost less, `GaussLegendre<N>` (no estimate) or
  `Romberg` (17 at least).

## Prototype

`scratchpad/gk-proto/` (session scratchpad, not kept; its tables go into `design.md`): the loop of C to E in double
mirroring Phase 1's operation order, a 168-bit mpmath version, rules from mpmath checked against Boost's 115-digit
table to 4e-117. For the counted tests the phase loads the C++ generator's double table into the prototype first
(numerics review: the prototype's own double rule rounds differently, and an ulp moves the floor decisions at
rel_tol = 0).

## Tests (`tests/test_quadrature_kronrod.cpp`, test cases `quadrature:kronrod:*`)

- **The rule** against Boost.Math's 115-digit G7-K15 table (BSL notice), in float, double, long double and both
  50-digit types, within an a-posteriori bound as in Phase 1.
- **Accuracy on smooth cases:** `|value - ref| <=` the derived rounding bound (about 4 eps L1 plus the weights'
  6.1 eps max |f| |b - a| / 2) where truncation is below it, with the honesty check `|value - ref| <= abs_error`
  beside it; Kahaner's 21 at the default against their closed forms (mpmath for 8, 12, 13, 17, 18) within rel_tol
  |ref|, the non-smooth ones only. In double, `cpp_bin_float_50` with and without expression templates, and
  `cpp_dec_float_50`.
- **Counts** exactly in double for fn of `+ - * /` and sqrt, against the prototype with the C++ table loaded: the
  mesh's 150, a kink, the D14 cases, x^-0.9 at two tolerances; a handful, not the battery.
- **Termination:** the cap (x^-0.9 at a tolerance it cannot meet), stagnation (the hash noise and 1/x^2: status,
  evaluations, best), final panels (x^-0.99: `not_converged`, best present; its estimate is not asserted honest),
  the 2 F test (a kink at rel_tol = 0 ends ok), D14 on [-1, 3] (log |x|, 1 / sqrt |x|), the bounded descent
  (sqrt(x - 0.5) on [0, 1]: `non_finite` within 105 evaluations; in `cpp_bin_float_50` too), the taint (a case with an
  undone piece never ends ok).
- **Edges:** Phase 1's edge tests run for both estimating methods as one template (orientation, a == b, NaN and
  infinite ends, tolerances, [-max, max], subnormal widths, the domain-checking callback, 1e-300 and 1e300 scales,
  an odd fn with abs_error the floor, overflow).
- **Constant evaluation:** `static_assert` on a few polynomial and rational cases within the step limits (the
  phase reports steps per bisection), with run-time twins; a NaN fn is `non_finite` there.
- **The cache:** no test beyond the multiprecision results themselves (a deterministic generator gives the same
  nodes cached or not).

## Documentation

- `quadrature.hpp`: `GaussKronrod`'s Doxygen comment (the cost of 150 evaluations at least, what the estimate is and
  is not, the limits of F briefly), the default method and the comments of the API section.
- `docs/quadrature/design.md`: the decisions in force, the sources (Kronrod 1965; Laurie 1997; Piessens et al. 1983;
  Malcolm and Simpson 1975; Shampine 2008; Gander and Gautschi 2000 for the battery), the prototype's tables.
- `docs/sphinx/quadrature/reference.md`: the new declarations.
- `README.md`: a three-line example of the default method.
- `demo/quadrature_demo.cpp` (`numerixx_quadrature_demo`), four cases (Q5): a smooth integral with the default, a
  strong endpoint singularity with the substitution, a failure and its best estimate, a multiprecision integral with
  its `rel_tol`. The fixed rule inside a root finder, nested integrals, Romberg and compile time go to the site's
  chapter at the end of the milestone, whose examples are tested.
- `docs/plan.md`: Phase 2's checkboxes.

## Decisions for the user

- **Q1. The initial mesh: 10 panels (150 evaluations at least).** Recommended: yes. Alternative: 4 panels (60
  evaluations; peaks beside a mesh point missed by up to 300 times, about 1 % more misses).
- **Q2. Stagnation: 32 bisections since the excess over the floor last reached a new minimum.** Recommended: yes
  (corrected twice: the brief watched the total estimate; draft 1 reset on a tenfold rise and missed divergence).
- **Q3. The cap: 64 bisections and descent splits per bit of T** (3,392 in double). Recommended: yes (a change to
  the brief's "from the exponent range"). Alternative: 32 per bit.
- **Q4. The stop test with twice the floor.** Recommended: yes (corrected from the brief).
- **Q5. The README example and a four-case demo now.** Recommended: yes.
- **Q6. The node cache for types without a table, for every rule.** Recommended: yes.
- **Q7. K15 only; K31 waits for Phase 3's `TanhSinh` measurements in multiprecision.** Recommended: yes.
- **Q8. Non-finite values: a descent of at most 2 levels; a piece no descent avoids is undone and the run can no
  longer end ok.** Recommended: yes (corrected twice). Alternative: carry the ancestor's estimate and drop the piece
  with value 0 (equal on the battery; on a NaN region its error was 3.4 times its estimate, against an honest one with
  the undo).
- **Q9. Endpoint singularities stronger than about x^-0.65:** documented, with the substitution x = a + u^k as a
  tested example, until `TanhSinh`. Recommended: yes.
- **Q10. D = K - G from two sums; no difference weights.** Recommended: yes (a change to the brief's 3.6 prose; no
  result changed on the battery).
- **Q11. `GaussKronrod<N>` constrained to N == 15 on the type,** for a clear diagnostic. Recommended: yes.
- **Q12. One `not_converged` for the cap, stagnation and final panels** (the usage sketch asked to tell stagnation,
  noise, from the others). Recommended: yes, one value now, documented with how to read `best`; `error` is documented
  as extensible, and a value of its own for stagnation can come when a caller needs to act on it.

## Decisions from the reviews

From the numerics review:
- Adopted: the descent bounded and counted against the cap (D); a run with an undone piece cannot end ok (D); the
  stagnation rule that trips on divergence (E, via the prototype's comparison); the D14 tests moved to [-1, 3], where
  the descent really happens; the counted tests against the prototype with the C++ table; the end gap 0.43 %, not
  0.85 % (F); constant evaluation: a refused value is `non_finite`, no descent (A); derived bounds, not rel_tol |ref|,
  on smooth cases, and the loop tested in both 50-digit types; no honesty assertion on x^-0.99; no reserve; empty
  mesh panels skipped; the heap-empty exit; `Criteria`'s 2 F; the wording on Emscripten's threads.
- Confirmed: the Kronrod weights and E_8's triangular system, the brackets and orientation, finality at the floor,
  the heap order, the Neumaier sums' adequacy, X's robustness, the 2 F test, the cap arithmetic, the API.
From the simplicity review:
- Adopted: D = K - G, no difference weights (B, Q10); three running sums with E derived (C); the compensated running
  K as the result (draft 1 also asked for a fresh sum of panels it had forgotten); one unsplittability predicate; the
  study material shortened here and moved to `design.md`; Phase 1's edge tests as one template; no tautological
  default-method or cache tests; a four-case demo; the cap's change to the brief marked.
- Not adopted: carrying the ancestor's estimate instead of the undo (Q8: less honest on a NaN region); dropping the
  smooth-case bounds for rel_tol |ref| (the numerics review: 10^5.6 too loose to catch a regression).
From the usage sketch:
- Adopted: the substitution recipe for strong endpoint singularities, with a tested example (F, Q9); `GaussKronrod`
  constrained on its N (Q11); the comments on `non_finite`'s infinite end, on "ok at the floor", on `Criteria`'s and
  `error`'s collisions; the Phase 3 / Phase 4 wording (finite intervals in Phase 3, infinite in 4).
- Not adopted now: `Failure<T>` for the fixed rule (BRIEF.md D4, approved); a value of its own for stagnation (Q12);
  combining `Integral`s of pieces (with break points, deferred); a nested-integral guide (the site's chapter).

## Implementation notes

Two deviations of the implementation, accepted by the implementation's review round: a panel whose L1 is finite but
whose 15-point value or difference K - G overflows (the 7-point weights exceed the 15-point ones at the Gauss nodes)
makes the run `non_finite`, a check beside D's one per panel on L1_p (tested: max on |x - 2.5| < 0.25 over [0, 50]);
and the cap's b counts the bits of a digit rounded up, as Phase 1's isolation cap does, so in radix 10 it is 64 x 4
bits per digit (`cpp_dec_float_50`: 12,800 against E's about 10,630).

## Deferred

- K21 ... K61 (Q7); QUADPACK's scaled estimator and the epsilon algorithm (BRIEF.md 6).
- `Integral`'s extra members (evaluations, L1); a value of `error` for stagnation (Q12).
- Break points, and summing the `Integral`s of pieces.
- Everything in BRIEF.md section 6.

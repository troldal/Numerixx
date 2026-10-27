# Design: roots (`nxx::roots`)

Decisions in force for `nxx::roots`. Each one says why; change it here before changing the code. What holds for every
module (the library, real types, the general test rules) is in `docs/design.md`; the architecture is in `BRIEF.md`,
the phases in `docs/plan.md`, and each phase's specification and evidence in `PHASE<n>_API.md`.

## Sources

- Brent, R. P. (1973). *Algorithms for Minimization without Derivatives*, Prentice-Hall (Dover reprint 2002), ch. 4:
  the procedure `zero` (section 6), its tolerance (eq. 2.9), the error bound of its result (eq. 2.18) and its
  worst case (eq. 3.4). Checked against the scan Brent publishes on his website (Phase 1). After Dekker (1969).
- Bisection: textbook.
- Alefeld, G. E., Potra, F. A. and Shi, Y. (1995). Algorithm 748: Enclosing zeros of continuous functions. *ACM TOMS*
  21(3), 327-344: Algorithm 4.2, its subroutines and its modified `bracket` (Phase 4), read from the scan on the first
  author's page. Boost.Math's `toms748_solve` (BSL-1.0) built and run as a cross-check, no code taken.
- Anderson, N. and Bjorck, A. (1973), *BIT* 13, 253-264: the rescaling, as King (1983) and the Wikipedia article
  "Regula falsi" state it; the paper itself was not read (Phase 4).
- Not used: Numerical Recipes (`zbrent`, `rtsafe`, `zbrac`, `zbrak`), GSL, and the Forsythe, Malcolm and Moler
  `zeroin` listing; Numerixx 1.x's root finders were surveyed (`BRIEF.md` section 2) but none of their code is kept.

## The interface (`roots.hpp`)

- The entry point is `solve<M>(...)` (user decision, 2026-10-27; `root<M>(...)` in Phase 1 and the brief), so that it
  is not told apart from its result type `Root<T>` only by case.
- One entry point, `solve<M>(fn, Bracket[, Criteria])`, for every bracketing method M (`Bisection`, `Brent`); a closed
  set of empty tags, `detail::bracketing<M>`, as in `nxx::deriv`. The open methods (Phase 2) take a guess as their
  start, by an overload; `BRIEF.md` 3.1.
- `Bracket<T>{lower, upper}`: the ends in either order (swapped, not an error: an interval is unambiguous).
- `Root<T>{value, abs_error}`, its own type rather than `nxx::deriv::Estimate<T>`, without a conversion to T; read by
  member name, so that members can be added (`BRIEF.md` 3.5).
- `Criteria<T>{abs_tol, rel_tol}`, taken as `Criteria<std::type_identity_t<T>>` so that T comes from the bracket and a
  braced list of doubles converts.
- `error`: `non_finite`, `no_sign_change`, `not_converged` (Phase 2: an open method found no root from its guess);
  extensible, formattable.
- fn is constrained by `nxx::detail::callable` (`math.hpp`), shared with `nxx::deriv`, or, where it gives
  derivatives, by `detail::gives<F, T, K>`: a `std::pair`, `std::tuple` or `std::array` returned by value, every
  element exactly T (a reference or an expression type could dangle), at least K of them. A tuple-like result is never
  read as a value; a bracketing method on a tuple fn uses its first element.
- Six overloads (Phase 2): a bracket or a guess, each with fn alone, with fn and dfn (Newton), or with fn, dfn and d2fn
  (Halley). The separate callables are joined by `detail::Joined` into one fn returning an array, so each method has
  one implementation, and each callable is evaluated once per point.

## Behaviour

- **Validation, in order:** tolerances (negative or not finite: `non_finite`), then clamped to [denorm_min, max / 4]
  and [2 eps, 1]; the ends (not finite: `non_finite`), put in order; fn at the ends (not finite: `non_finite`; 0: that
  end is the root; the same sign: `no_sign_change`).
- **Zero first:** a bracket containing 0 in its interior is split at 0 after evaluating fn(0), so a root at 0 costs
  one evaluation and the rest sees one sign; skipped where fn(0) is not finite (removable singularities: sin(x)/x,
  x/(e^x - 1)). User decision P2.
- **The stop test:** the bracket at most 2 tol(x) wide, tol(x) = abs_tol + rel_tol |x| (saturated at max once
  rel_tol |x| exceeds max / 2, so that the sum cannot overflow), or fn exactly 0. The defaults, also the minimum, are
  `abs_tol = denorm_min` and `rel_tol = 2 eps`: full precision at every magnitude, no scale assumed (MATLAB and SciPy
  assume 1, and a root at 1e-12 loses its digits there). User decision P1. Brent's own tolerance is 2 macheps |b| + t
  with macheps half of epsilon for rounded arithmetic, so his relative term is eps |b|; ours is twice it, and he
  advises against less.
- **The walk to 0** (documented, user decision P6): a sign change near 0 that fn does not evaluate to exactly 0 (a
  root fn misses at 0, a tiny root, a jump at 0) is approached down to T's smallest numbers, about 1100 halvings in
  double and about 2e9 in `cpp_bin_float_50`. The remedy is `abs_tol`; a fast descent toward 0 is deferred.
- **Distances:** every abs_error, stop test and guard uses `detail::distance_bound(u, v)`, an upper bound on
  |u - v| formed without overflow: exact where the difference is (ends of one sign whose computed difference is below
  the smaller magnitude, Sterbenz's lemma; in radix 2 also its boundary case; an end at 0), which is every stop with
  the default tolerances; otherwise the computed distance raised by one unit, `d + max(d eps, denorm_min)`;
  denorm_min where ends that differ give a difference of 0 (a type without subnormals flushes it); max where the
  distance may exceed what T holds. Added after the independent review (below).
- **Bisection:** midpoints `a + (b - a) / 2` (one sign) or `(a + b) / 2` (opposite signs), neither overflowing;
  `abs_error` the larger of the bounds on the distances to the ends, a true bound even where half the width rounds.
- **Brent:** the procedure `zero`, with the caller's tolerance, the stop test on a bound on |c - b| against 2 tol(b),
  and a bisection step whenever that bound, the bound on |b - a| or |e| exceeds max / 16, so that no product of the
  interpolation overflows (3 m q reaches 12 |c - b|; an infinite one would accept a step rather than reject it). The
  quotients cannot overflow: |fb| < |fa| <= |fc| in the interpolation. `abs_error` is the bound on |c - b|.
- **The result:** `abs_error` bounds the distance to a sign change of fn as evaluated and is at most 2 tol(value),
  saturating at max where the distance exceeds what T holds (only with tolerances near max);
  where fn is exactly 0 it is tol(value), not 0 (user decision P4; e^x - 2 is exactly 0 at a point 2.3e-17 from ln 2).
  The conditioning of the root, fn's rounding over |f'|, is the caller's.
- **Poles** are returned as roots, as GSL, Boost and SciPy do, and documented (user decision P3). The rule first
  proposed (|fn| at the final points above |fn| at the given ends) was unsound both ways; a local-growth rule is
  deferred, with the three counterexamples of `PHASE1_API.md` E as its tests.
- **No iteration cap** for the bracketing methods: Bisection ends at T's resolution, Brent within (k + 1)^2 - 2
  evaluations; a cap without an error to report it would return the last iterate silently.
- Everything is constexpr; nothing throws, and fn's exceptions pass through.

## Open methods and derivatives (Phase 2)

- **Every result has a verified sign change** (user decision P1, superseding `BRIEF.md` 3.5 and 3.7 for the open
  methods): an open method iterates from its guess until fn's values at two consecutive iterates have opposite signs,
  then continues as the bracket loop on them (Newton and Halley safeguarded, Secant as Brent), without zero first. The
  textbook stop, |step| <= tol, accepted a false root in the study (Secant on x^20 - 1 from 0.5, a step of 1e-84 where
  fn was -1) and never stopped on Wilkinson's W10; a stagnation rule stopped 3.7 from a root. The hand-over costs
  about one evaluation over the textbook stop.
- **One bracket loop:** Phase 1's Brent with the step chosen by the method: the interpolation, or Newton's or
  Halley's step from b, accepted by Brent's test written for a step (towards c, |s| < 3 |m| / 2 - tol / 2,
  |s| < |e| / 2). Brent's bookkeeping, guard and termination argument hold for all three; Phase 1's results are
  unchanged bit for bit.
- **Derivatives:** only fn's value is checked; a derivative that is 0 or not finite makes the step unusable (on a
  bracket a bisection step, from a guess `not_converged`), so f'(0) = inf at an end or at zero first does not fail a
  solve.
- **Halley's step** `-t / (1 - u)`, t = f / f', u = t (f'' / f') / 2, without f'^2; Newton's step where u cannot be
  formed, 1 - u is 0, or the step points against Newton's (Boost.Math does the same). It costs one case: from 0 on
  x^3 - 2x + 2, Halley follows Newton's cycle to `not_converged` (user decision P5).
- **Secant's step** `(x - xp) (f / (fp - f))`, the ratio first (fp - f cannot overflow before the hand-over); its
  second point `x0 (1 + eps^(1/4))`, relative (P4), or x0 + tol(x0) where that rounds back to x0. Every point it
  evaluates is checked, the last before the cap included.
- **The minimum step** tol(x), in the direction of the method's step taken from the signs (of f and f'; of f, fp - f
  and x - xp), not from a step that may have underflowed.
- **Overflow:** every quotient (directly where |b| >= 1, else only if `|a| < (max / 2) |b|`), product (`|a| <= 1`
  or `|b| <= (max / 2) / |a|`) and x + s (opposite signs, both at most max / 2, or exact above max / 2) is checked
  before it is formed; a step that fails is unusable. The quotient's one-line test `|a| / max < |b| / 2` rejected
  ordinary quotients whose operands underflow (independent review, `PHASE2_API.md`).
- **Progress:** the bracket loop stops where its next point would not land strictly inside the bracket (it is then at
  the tolerance's resolution); `distance_bound` measures ends of opposite signs above max / 4 by halves rather than
  saturating. Together they closed an endless loop under huge tolerances (independent review).
- **The cap** of the open methods is 4 times the bits of T (212 for double, 672 for `cpp_bin_float_50`; P3). Beyond
  it: a multiple root at 0 with the default abs_tol, a guess very far from the root, multiplicities above about 5
  (Newton), 3 (Secant), 11 (Halley). A bracket solves them. One `not_converged` for every failure (P2).
- **An open method returns a root,** the first sign change its iterates cross, not necessarily the one nearest the
  guess; documented.

## Numerical derivatives (Phase 3)

- **A separate derivative may return a std::expected** whose value converts to T (user decision P3): any
  `std::expected<U, E>`, so `nxx::roots` names no type of `nxx::deriv`; `nxx::deriv::derivative_of`'s callable and a
  lambda around `derivative(fn, x, h)` are the uses. fn itself still returns a value or a tuple (deferred).
- **A failed derivative counts as one that is not finite** (P2, superseding `BRIEF.md` 3.4's `non_finite`): on a
  bracket a bisection step, from a guess `not_converged`, for Halley's f'' Newton's step. Stored as 0, which every
  `nxx::real` has and which has the same effect in both steps (tested against a derivative returning infinity, bit for
  bit). Final for behaviour: an error later would turn roots found today into errors.
- **The cost is documented, not prevented:** with `Central<2>` a point costs 3 calls of fn for Newton, 6 for Halley
  (f'' evaluates fn(x) again; sharing it is deferred, P6). Newton then costs 1.5 to 2.3 times Secant from the same
  guess and 1.7 to 3.5 times Brent on the same bracket (6.1 on ln(x / 1e-12) over [1e-20, 1]); Halley 1.5 to 4.2 times
  Secant; with `nxx::deriv::Ridders` 6 to 11 times Secant.
- **The bound never depends on the derivative; the iterates can.** `Central<2>`'s error is absolute (h^2 f''' / 6,
  h^2 = 3.7e-11 near |x| <= 1, plus rounding). Where |f'| is far above it a run takes as many points as the analytic one (they differ in their last bits); where
  comparable, the outcome changes either way: (x - 1)^3 from 2 and x^2 - 2 from 1e-12 become `not_converged`,
  x^3 - 1e-3 from 0 converges. A bracket always converges.
- **The default step assumes fn varies on the scale max(|x|, 1):** ln(x / 1e-12) from 1e-13 has a failed f' at the
  guess (`not_converged`), and on [1e-20, 1] every point below 6e-6 has one (pure bisection, 276 calls).

## More bracketing methods (Phase 4)

- **`Toms748` is Algorithm 4.2 as printed** (user decision P1): a secant step; then per iteration the inverse cubic
  (`ipzero`) or `Newton-Quadratic` with 2 steps, the same with 3 steps, a double-length secant step, and a bisection
  where the iteration has not halved the bracket. Boost.Math's arrangement (a quadratic step before the loop, a
  fallback of 3 steps in 4.2.3) differs; it was the prototype's first form and is the cross-check.
- **`AndersonBjorck` is guarded** (P2): a bisection after three steps that have not halved the bracket since the last
  halving, a guard bisection counting as one. The published method took 65,109,847 evaluations on x^20 - 1.
- **One step for both, `detail::enclose`** (the paper's modified `bracket`): delta = 0.7 tol(u), at least denorm_min
  (Boost.Multiprecision's denorm_min() is min(), and 0.7 min flushes to 0); the midpoint where the bracket is at most
  4 delta wide, where the method has no point, and where the point kept 2 delta inside is not strictly inside (2 delta
  rounds away at an end of much larger magnitude than u). No point is formed while the bracket is wider than max / 16.
- **Arithmetic:** unlike Brent's, these steps form differences of fn's values of opposite signs; every quotient,
  product and sum on fn's values and the interpolations' terms is checked, and a failed check makes the point
  unusable, so the paper's fallbacks apply (inverse cubic, quadratic, secant, midpoint).
- **The cost** (`PHASE4_API.md`): neither is cheaper than Brent in general (309, 325 and 327 evaluations on 27 of the
  paper's cases at full precision; the paper reports 4% fewer for 4.2 over its own set); both bound the worst case at
  about 4 evaluations per halving. AndersonBjorck is far cheaper where fn is strongly curved near one end and on a walk
  to a sign change at 0 (51 evaluations against Brent's 1076). Brent stays the default.

## Searching (Phase 5)

- **`search(fn, x0)`** (the brief's `expand`, renamed by the user so that a search that does not expand can share the
  name; a second algorithm would be a defaulted leading tag, `search<M>(fn, x0)`): fn at x0, then at x0 - h and
  x0 + h, h = |x0| / 50 (1/50 at 0, at least denorm_min) growing by 3/2, until the sign differs from fn(x0); the result
  is the narrowest bracket that shows it, its ends in order, or {x, x} where fn is exactly 0 at x (user decision P1).
- **Domain edges:** a non-finite value makes that side bisect towards its last finite point, at most as many probes as
  T has bits (P2; BRIEF.md 3.9 had the side stop, which misses roots near an edge).
- **No search bound** (P3): where fn is finite outside its domain, the caller returns NaN there; a bound (`within`)
  is deferred. **No sign change found is `no_sign_change`** (P4).
- **Limits:** at most 4 times the bits of T steps (2 cap + 1 evaluations); h grows only where 1.5 h fits, so no step
  passes max (from max the search reaches down to 0.23 max); x0 +- h formed only where `sum_fits` admits it; the first
  step |x0| where denorm_min() is 0. What no outward search
  finds, and the documentation says: an even root, two roots within one step, a root beyond a non-finite region on
  the same side, a root beyond the cap; a pole brackets like a root.

## The teaching methods (Phase 6)

- **`Ridders`** (Ridders 1979; formula as pracma documents it): per iteration `enclose` at the midpoint, then at the
  exponential fit's zero computed from the state that step leaves, in the cancellation-free form x4 = e + w (m - e),
  w = -g_e g_d / (root (root + |g_m|)) (user decision P1; the textbook form bisected through 1922 evaluations on
  x - 1e-300). Its square root is `detail::sqrt` (P5): faithful, from + - * /, Newton in the correction form (the
  textbook form is 3 ulps off in radix 10). Not `nxx::deriv::Ridders` (P4).
- **`RegulaFalsi`** and **`Illinois`**: `false_position<M>`, AndersonBjorck's loop with the factor 1 or 1/2, guarded
  (P2): unguarded, regula falsi ran past 200,000 evaluations on x^20 - 1.
- **`Steffensen`**: z = x + f(x) (at least tol(x) from x), the secant step through x and z shared with Secant
  (`secant_move`), Brent at the first sign change between x and z or the next iterate (P3). f's scale changes the cost
  and the root found (1000 (x^2 - 2) from 1: -sqrt 2).
- The cost (`PHASE6_API.md`): on 27 of the paper's cases Brent 309, AndersonBjorck 327, Illinois 382, Ridders 416,
  RegulaFalsi 645; Ridders 86 on (x - 1)^3 against Brent's 156. Teaching examples; Brent stays the default.

## Evidence (Phase 1)

An independent Python implementation of the specification in double (`.scratch/roots1/proto3.py`, references from
mpmath at 60 digits) gave these counts, which the C++ reproduces exactly for fn of `+ - * /`:

| Case | Bracket | Bisection | Brent |
|---|---|---|---|
| x^2 - 2 | [1, 2] | 52 | 9 |
| x^3 - 2x - 5 (Wallis) | [2, 3] | 51 | 8 |
| x^20 - 1 | [0, 5] | 55 | 19 |
| (x - 1)^3 | [0, 3] | 54 | 156 |
| x - 1e-300 | [0, 1] | 1049 | 3 |
| x - 1 | [-max, max] | 1077 | 8 |
| x - 1 | [0, max] | 1076 | 7 |
| x^3 (root at 0, zero first) | [-1, 2] | 3 | 3 |
| a jump at 0 | [-1, 2] | 1076 | 1076 |

With libm fns the prototype's results (cos x - x 9, e^x - 2 8, x e^x - 1 9, ln(x / 1e-12) 45 evaluations for Brent)
were all within `abs_error` of the references; the C++ tests check accuracy there, not counts.

## Evidence (Phase 2)

The independent Python implementation `.scratch/roots2/proto2.py` (double, on `.scratch/roots1/proto3.py`'s Brent)
gives these evaluation counts, which the C++ reproduces exactly for fn of `+ - * /` (Newton / Halley / Secant):

| Case | From a guess | On a bracket (Newton / Halley) |
|---|---|---|
| x^2 - 2 | x0 = 1: 7 / 5 / 9; 10: 9 / 7 / 13; 1e-3: 17 / 12 / 25; fl(sqrt 2): 2 / 2 / 4 | [1, 2]: 8 / 6 |
| x^3 - 2x - 5 | x0 = 2: 6 / 5 / 8; 0: 21 / 12 / 15; 10: 11 / 8 / 15 | [2, 3]: 7 / 6 |
| (x - 1)^3 | x0 = 2: 89 / 53 / 127 | [0, 3]: 91 / 54 |
| x^20 - 1 | x0 = 0.5: 22 / 11 / 44; 2: 20 / 11 / 28 | [0, 5]: 11 / 11 |
| x^3 - 2x + 2 | x0 = 0: `not_converged` (213) / `not_converged` / 28; -2: 7 / 5 / 9 | [-2, -1]: 8 / 6 |
| 1/x - 1 | x0 = 0.5: 7 / 2 / 10; 3: `not_converged` (10) / 3 / `not_converged` | [0.5, 3]: 9 / 6 |
| 1e300 (x - 1) | x0 = 1e7: 2 / 2 / 4 | |
| W10 (Wilkinson) | x0 = 7.2: 21 / 19 / 23; 6.8: 22 / 16 / 23 | [6.7, 7.6]: 28 / 24 (Brent 23) |

## Evidence (Phase 3)

`.scratch/roots3/proto3.py` (Phase 2's prototype with `nxx::deriv`'s `Central<2>` f' and f'' reproduced step for step)
gives these calls of f, which the C++ reproduces exactly for fn of `+ - * /` and sqrt (Newton / Halley, numerical;
the analytic counts in brackets):

| Case | From a guess | On a bracket |
|---|---|---|
| x^2 - 2 | x0 = 1: 21 / 30 (7 / 5); 10: 27 / 42 (9 / 7); 1e-12: `not_converged` 3 / 6 (47 / 31) | [1, 2]: 24 / 36 (8 / 6) |
| x^3 - 2x - 5 | x0 = 2: 18 / 30 (6 / 5); 10: 33 / 48 (11 / 8) | [2, 3]: 21 / 36 (7 / 6) |
| (x - 1)^3 | x0 = 2: `not_converged` 639 / 1278 (89 / 53) | [0, 3]: 549 / 1008 (91 / 54) |
| x^3 - 1e-3 | x0 = 0: 105 / 204 (`not_converged`) | [0, 1]: 33 / 54 (11 / 9) |
| 1/x - 1 | x0 = 0.5: 21 / 18 (7 / 2) | [0.5, 3]: 27 / 36 (9 / 6) |
| x^20 - 1 | x0 = 0.5: 66 / 66 (22 / 11); 2: 60 / 66 (20 / 11) | [0, 5]: 33 / 66 (11 / 11) |
| sqrt x - 0.5 | x0 = 1: 15 / 30 (4 / 5) | [0, 1]: 12 / 24 (4 / 4) |
| ln(x / 1e-12) | x0 = 1e-13: `not_converged` 3 / 6 (8 / 6) | [1e-20, 1]: 276 / 552 (20 / 18) |

`nxx::deriv::Ridders` (the library's, measured, not a reference): Newton from a guess / on a bracket 53 / 58 calls on x^2 - 2,
69 / 78 on cos x - x, 91 / 95 on e^x - 2, 107 / 112 on Kepler.

## Evidence (Phase 4)

`.scratch/roots4/proto4.py` (revision 2: the C++'s checks for every operation) gives these counts, which the C++
reproduces exactly for fn of `+ - * /` and sqrt (Brent / Toms748 / AndersonBjorck):

| Case | Bracket | Counts |
|---|---|---|
| x^2 - 2 | [1, 2] | 9 / 9 / 9 |
| x^3 - 2x - 5 | [2, 3] | 8 / 9 / 9 |
| x^20 - 1 | [0, 5] | 19 / 21 / 30 |
| (x - 1)^3 | [0, 3] | 156 / 149 / 152 |
| x - 1 | [-max, max] | 8 / 8 / 13 |
| 1/x - 1 | [0.5, 3] | 11 / 6 / 5 |
| W10 near 7 | [6.7, 7.6] | 23 / 18 / 18 |
| a jump at 0 | [-1, 2] | 1076 / 981 / 51 |
| 1/(x - 1/3) | [0, 1] | 55 / `non_finite` (84) / 95 |
| sqrt x - 0.5 | [0, 1] | 4 / 6 / 9 |
| Alefeld, Potra and Shi's problems #1, #2, #4 to #12 (27 cases) | | 309 / 325 / 327 in all |

Under Boost.Math's arrangement, placement and stop rule, the prototype takes Boost.Math 1.92.0's own counts on 7 of 8
cases (one fewer on the eighth).

## Evidence (Phase 5)

`.scratch/roots5/proto5d.py` gives these evaluations and brackets, which the C++ reproduces exactly (to the last bit
of each end) for fn of `+ - * /` and sqrt: x^2 - 2 from 1 (19) and 0 (24); x^3 - 2x - 5 from 0 (27); 1/x - 1 from 3
(20; growth 2 misses it); x - 1e6 from 1 (91); x - 1e-300 from 1 (22); x + 1e-300 from denorm_min (266); x - 1 from
max / 2 (22); sqrt x - 1e-3 from 1 (60, through a domain edge); van der Waals from 1e-3 (20), from 1e-4 and 5e-5
around its pole (20, 12) and, with NaN outside v > b, to the root (29, 33); `no_sign_change` after 425 for x^2 - 2
from 100, x^2 + 1, (x - 1)^2 and x - 1 from 1e-300.

## Evidence (Phase 6)

`.scratch/roots6/proto6.py` (with `detail::sqrt` ported step for step) gives the counts of `PHASE6_API.md`'s tables,
which the C++ reproduces exactly for fn of `+ - * /` and sqrt: Ridders / Illinois / RegulaFalsi on x^2 - 2 12 / 11 /
20, (x - 1)^3 86 / 148 / 145, x^20 - 1 18 / 30 / 45, the jump at 0 778 / 51 / 1076, and problems #2, #4, #7, #8, #9,
#11; Steffensen 10 (sqrt 2) from 1, 39 from 10, 11 (-sqrt 2) from 0.5, 27 (-sqrt 2) on 1000 (x^2 - 2) from 1.

## Tests

- References from closed forms and published constants, at 70 digits for the multiprecision types (more digits than
  the type: a 50-digit reference errs by up to ten epsilons of `cpp_bin_float_50`), cross-checked with mpmath.
- Accuracy: `|value - ref| <= abs_error + c` with c, the conditioning term, derived per case from fn's error budget
  (its own rounding for `+ - * /`, 2 ulps for a libm function) over |f'(root)|; and `abs_error <= 2 tol(value)`.
- Counts exactly for fn of `+ - * /` in double, against the prototype; none where a libm or a multiprecision type's
  rounding is involved.
- The edge cases in double; a few accuracy cases and a `static_assert` per other type (no `static_assert` for the
  Boost types, which are not literal types).
- Numerical derivatives (`tests/test_roots_numerical.cpp`, `roots:numerical:*`): counts are calls of the user's f.
  A derivative returning an expected against the plain one, a failed derivative against one returning infinity, and
  a derivative that fails only sometimes, chosen so that a step that would be accepted is refused, with a check that
  the run differs from the analytic one (so the test is not vacuous).

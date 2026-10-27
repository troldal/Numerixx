# Design: roots (`nxx::roots`)

Decisions in force for `nxx::roots`. Each one says why; change it here before changing the code. What holds for every
module (the library, real types, the general test rules) is in `docs/design.md`; the architecture is in `BRIEF.md`,
the phases in `docs/plan.md`, and each phase's specification and evidence in `PHASE<n>_API.md`.

## Sources

- Brent, R. P. (1973). *Algorithms for Minimization without Derivatives*, Prentice-Hall (Dover reprint 2002), ch. 4:
  the procedure `zero` (section 6), its tolerance (eq. 2.9), the error bound of its result (eq. 2.18) and its
  worst case (eq. 3.4). Checked against the scan Brent publishes on his website (Phase 1). After Dekker (1969).
- Bisection: textbook.
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
  second point `x0 (1 + eps^(1/4))`, relative (P4).
- **The minimum step** tol(x), in the direction of the method's step taken from the signs (of f and f'; of f, fp - f
  and x - xp), not from a step that may have underflowed.
- **Overflow:** every quotient (`|a| / max < |b| / 2`), product (`|a| <= 1` or `|b| <= (max / 2) / |a|`) and x + s
  (opposite signs, both at most max / 2, or exact above max / 2) is checked before it is formed; a step that fails is
  unusable.
- **The cap** of the open methods is 4 times the bits of T (212 for double, 672 for `cpp_bin_float_50`; P3). Beyond
  it: a multiple root at 0 with the default abs_tol, a guess very far from the root, multiplicities above about 5
  (Newton), 3 (Secant), 11 (Halley). A bracket solves them. One `not_converged` for every failure (P2).
- **An open method returns a root,** the first sign change its iterates cross, not necessarily the one nearest the
  guess; documented.

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

## Tests

- References from closed forms and published constants, at 70 digits for the multiprecision types (more digits than
  the type: a 50-digit reference errs by up to ten epsilons of `cpp_bin_float_50`), cross-checked with mpmath.
- Accuracy: `|value - ref| <= abs_error + c` with c, the conditioning term, derived per case from fn's error budget
  (its own rounding for `+ - * /`, 2 ulps for a libm function) over |f'(root)|; and `abs_error <= 2 tol(value)`.
- Counts exactly for fn of `+ - * /` in double, against the prototype; none where a libm or a multiprecision type's
  rounding is involved.
- The edge cases in double; a few accuracy cases and a `static_assert` per other type (no `static_assert` for the
  Boost types, which are not literal types).

# Phase 1 API proposal: the bracketing core, Bisection and Brent

Status: approved by the user on 2026-10-25 (P1 to P5 as recommended; P6: the walk to 0 documented, the fast descent
deferred). One review round (numerics-reviewer, simplicity-critic, usage-sketcher) is folded in ("Decisions from the
reviews"). Architecture: `docs/roots/BRIEF.md` (approved 2026-10-25).

## Summary

```cpp
#include <numerixx/roots.hpp>
namespace rt = nxx::roots;

auto const f = [](double x) { return x * x - 2.0; };
auto const r = rt::root<rt::Brent>(f, rt::Bracket{1.0, 2.0});      // sqrt(2) to about 1e-16, 9 evaluations
if (r) { double const x = r->value; double const bound = r->abs_error; }
auto const s = rt::root<rt::Bisection>(f, rt::Bracket{2.0, 1.0}); // either order; 52 evaluations
auto const t = rt::root<rt::Brent>(f, rt::Bracket{1.0, 2.0}, {.rel_tol = 1e-8}); // a looser tolerance
auto const u = rt::root<rt::Brent>(f, rt::Bracket{2.0, 3.0});      // error::no_sign_change

static_assert(rt::root<rt::Bisection>([](double x) { return x - 0.375; }, rt::Bracket{0.0, 1.0})->value == 0.375);
```

- A new header, `include/numerixx/roots.hpp`, with the namespace `nxx::roots`. It includes `math.hpp` only.
- Two methods, `Bisection` and `Brent`, one entry point, `root<M>(fn, Bracket[, criteria])`.
- Everything is constexpr, and nothing throws.

## Public API (`roots.hpp`)

```cpp
namespace nxx::roots {

    /// Why root() returned no root. More values will be added (Phase 2: not_converged), and a listed cause may move to
    /// a new value: a switch should have a default. Formattable: std::format("{}", e) writes the enumerator's name.
    enum class error {
        non_finite,     ///< A bracket end, a value of fn or a tolerance is not finite, or a tolerance is negative.
        no_sign_change, ///< fn has the same sign at both ends of the bracket, and is 0 at neither.
    };

    /// An interval [lower, upper] in which fn changes sign. The ends may be given in either order.
    template<real T>
    struct Bracket
    {
        T lower;
        T upper;
        friend constexpr bool operator==(Bracket const &, Bracket const &) = default;
    };

    /// A root, and a bound on its distance to a sign change of fn as evaluated (at most 2 tol(value); see F). Read it
    /// by member name: members may be added (BRIEF.md 3.5).
    template<real T>
    struct Root
    {
        T value;
        T abs_error;
        friend constexpr bool operator==(Root const &, Root const &) = default;
    };

    /// When to stop: when the bracket is at most 2 tol(x) wide, tol(x) = abs_tol + rel_tol |x|. The defaults ask for
    /// the full precision of T at every magnitude, and are also the minimum: smaller values are raised to them.
    template<real T>
    struct Criteria
    {
        T abs_tol = T(std::numeric_limits<T>::denorm_min());
        T rel_tol = T(2) * T(std::numeric_limits<T>::epsilon());
        friend constexpr bool operator==(Criteria const &, Criteria const &) = default;
    };

    struct Bisection {}; ///< Halves the bracket: one bit per evaluation, guaranteed.
    struct Brent {};     ///< Brent's method: inverse quadratic and linear interpolation, safeguarded by bisection.

    template<typename M, real T, typename F>
        requires detail::bracketing<M> && nxx::detail::callable<F, T>
    [[nodiscard]] constexpr std::expected<Root<T>, error> root(F &&fn, Bracket<T> const &bracket,
                                                               Criteria<std::type_identity_t<T>> const &criteria = {});

} // namespace nxx::roots

template<> struct std::formatter<nxx::roots::error> : std::formatter<std::string_view> { /* the enumerator's name */ };
```

- **`nxx::detail::callable<F, T>`** is `nxx::deriv`'s concept (fn invocable with `T const &`, result convertible to
  T), moved from `deriv.hpp` to `math.hpp` for both modules (the shared-core rule: the second module needs it). Every
  constraint spells it `nxx::detail::callable`, since inside `nxx::deriv` or `nxx::roots` the name `detail::` finds
  the module's own namespace first. Wrapped in `@cond INTERNAL`, as `deriv.hpp`'s details are. Nothing public changes.
- **`detail::bracketing<M>`** is true for `Bisection` and `Brent`; a closed set, as `nxx::deriv`'s methods.
- **No `not_converged`, and no iteration cap, yet.** Both methods end by construction (C, D). A cap could only
  return the last iterate silently (1.x's defect) without an error value to report it; Phase 2 adds both with the
  open methods. This amends BRIEF.md 3.6 and 3.7, which planned an internal cap for every method.
- **No `pole`.** The rule proposed for it was unsound (E); a pole is returned as a root, as GSL, Boost and SciPy do,
  and documented.

## Behaviour

### A. Validation at the boundary, in this order

1. **Tolerances.** Negative or not finite: `non_finite` (as `nxx::deriv` reports a step that is not positive and
   finite; an assert would vanish under `NDEBUG`, and a NaN tolerance would then never stop Brent). Then each is
   clamped: `abs_tol` to [`denorm_min`, `max / 4`], `rel_tol` to [`2 eps`, 1]. The lower clamps make Brent's minimum
   step move b (a step of `tol >= 2 eps |b|` is at least two units of T's spacing in radix 2 and 10; `rel_tol = eps`
   would also do, `eps / 2` would not); the upper ones keep `tol(x)` and `2 tol(x)` finite, which constant
   evaluation needs. `denorm_min() > 0` is asserted (every standard and Boost type has it).
2. **Ends.** Not finite: `non_finite`. Put in order (`lower > upper` is swapped; BRIEF.md 3.3).
3. **fn at both ends.** Not finite: `non_finite`. Exactly 0: that end is the root. Otherwise the same sign at both
   ends: `no_sign_change`. Equal ends follow from these rules.
4. **Zero first.** If the bracket contains 0 in its interior, fn(0) is evaluated next (documented on `root`): 0 is the
   root if fn(0) is 0; if fn(0) is finite, the bracket becomes [lower, 0] or [0, upper]; if it is not finite, zero
   first is skipped and the bracket is kept as given.

Why zero first (P2):
- **A root at 0 costs one evaluation.** With the default `abs_tol`, the full precision of T at every magnitude, a root
  at 0 is otherwise approached down to T's smallest numbers: in the first prototype Brent on x^3 over [-1, 2] took 1090
  evaluations (the triple root makes it linear) until x^3 underflowed, and in `cpp_bin_float_50` (exponent range
  about 2^31) the same walk takes about 2e9.
- **The rest sees a bracket of one sign** in the usual case, where differences of points cannot overflow.
- **Skipped where fn is not finite at 0:** removable singularities are common (sin(x)/x, x/(e^x - 1), (e^x - 1)/x),
  and an evaluation the caller did not ask for must not fail their solve. The bracket then keeps both signs, which the
  midpoint and Brent's guard (C, D) handle.
- It costs one evaluation for a bracket around 0 whose root is not 0; in a symmetric bracket it also halves it.

### B. The stop test, and the walk to 0

`tol(x) = abs_tol + rel_tol |x|`, with the clamped tolerances. A method stops when its bracket is at most `2 tol(x)`
wide (Bisection: half the bracket, `max(m - a, b - m)`, at most `tol(m)`; Brent: `|c - b| <= 2 tol(b)`, with `c - b`
exact for points this close), or when fn is exactly 0 at a point.

The defaults assume no scale (MATLAB's `fzero` assumes 1, with an absolute floor of eps; SciPy's `brentq` an absolute
2e-12), so a root at 1e-12, a mole fraction of a trace component say, comes out to about 2 eps relative, not 2e-12
absolute. The cost of that choice: a sign change near 0 that fn does not evaluate to exactly 0 is approached down to
T's smallest numbers, about `log2(width / denorm_min)` halvings. That covers a root that fn misses at 0, a genuine
tiny root (1049 evaluations for 1e-300 in double with Bisection; Brent needs 3) and a jump at 0, such as a sign
function (1076 in double for both methods, since Brent bisects at a jump). In double this is at most about 1100
evaluations; in `cpp_bin_float_50` about 2e9 and in `cpp_dec_float_50` about 2e8, which does not end in practice. The
documentation says so and gives the remedy, `abs_tol` at the scale where 0 is close enough. A fast descent toward 0
(trial points `b radix^-k`, k doubling, then bisection; about 31 more evaluations at most) would bound it by the
exponent's digits rather than its range; it is deferred (user decision P6).

### C. Bisection

- The midpoint: `a + (b - a) / 2` for ends of one sign, `(a + b) / 2` for ends of opposite signs (only when zero first
  was skipped); neither overflows. Exact in binary when a and b are within a factor of 2, else rounded, harmlessly.
- It stops by the test of B, when the midpoint is no longer strictly between the ends (T's resolution), or when fn is
  0 at a midpoint. A value of fn not finite inside the bracket: `non_finite`.
- **Result:** the midpoint of the final bracket, `abs_error = max(m - a, b - m)` (both differences exact, by
  Sterbenz's lemma; half the width can round to 0 among subnormals, `2x - 5 denorm_min` on [0, 3 denorm_min] being an
  example); fn 0 at a midpoint: that point.
- **Cost:** 2 evaluations, one more for zero first, then one per halving: about `log2(width / (2 tol))`, so 52 for
  [1, 2] in double.
- **Termination:** every step moves an end strictly inward, in any radix.

### D. Brent

- Brent's procedure `zero` (Brent 1973, *Algorithms for Minimization without Derivatives*, ch. 4; after Dekker 1969),
  implemented from the book: the bracket [b, c] with |fn(b)| <= |fn(c)|, a step by inverse quadratic interpolation
  through three points or by the secant through two, accepted only if `2p < 3mq - |tol q|` and `p < |e q / 2|` (well
  inside the bracket, and shrinking faster than the step before last), else a bisection step; every step at least
  `tol(b)`. Not from Numerical Recipes (`zbrent`, whose acceptance test differs), GSL, or the Forsythe, Malcolm and
  Moler `zeroin` listing. To check against the book while implementing: whether Brent's `macheps` is epsilon or half
  of it (our default `rel_tol = 2 eps` stands either way), the error bound he states for the result, and the
  worst-case count (about `(k + 1)^2 - 2` for k bisection steps).
- **No overflow is formed.** The quotients cannot overflow: in the interpolation |s|, |r| < 1 and |q| <= 1 (the
  swap and the reset order the values). The products can: `3 m q` reaches 12 |c - b| and `2p` 4 (|c - b| + |b - a|),
  so `[0, max]` overflows `3 m q` at the first step (and an infinite `3 m q` would accept a step, not reject it). The
  guard: when |c - b|, |b - a| or |e| exceeds `max / 16` (decided by comparisons before any difference is formed), the
  step is a bisection step. Below about 1e307 in double Brent is unchanged.
- **Result:** b, with `abs_error = |c - b|` (fn(b) and fn(c) keep opposite signs, so a sign change lies between them;
  the difference is exact); fn(b) exactly 0: b.
- **Cost:** 2 (+1) evaluations, then one per step: 8 or 9 on the smooth cases below, 156 on the triple root of
  (x - 1)^3 (Bisection 54), at most about the square of Bisection's count.
- **Termination:** the brackets are nested, every step is at least `tol(b)`, and an accepted step stays strictly
  inside; the guard only adds bisection steps.

### E. Poles: returned as roots, documented

A pole inside the bracket changes sign like a root, and both methods converge to it. The proposal's rule (report a
pole when |fn| at the final points exceeds |fn| at the given ends) was unsound both ways (both reviews):
- a false pole: (x - 1/3) exp(-1000 (x - 1/3)^2) on [0, 1], whose end values are 2e-49 while the final points' are
  at least 5.6e-17, and a resonance, (x^2 - 2) / ((x^2 - 2)^2 + 1e-24) on [1, 2], whose final values exceed 1e8 at a
  well-conditioned root;
- a missed pole: 1/(x - 1/3) + 1e20 (x - 1/3)^3 on [0, 1], whose end values exceed anything near the pole.

So Phase 1 returns the pole's location, as GSL, Boost and SciPy do, and the documentation says to check |fn(value)|
when fn may have poles. A sound rule judges growth locally (|fn| at the final points against |fn| when the bracket was
about 2^10 times wider: a root's values shrink as the bracket closes, a pole's grow); it is deferred, with these three
functions as its tests. Reporting a pole later turns a result into an error: the documentation says now that it may.

### F. The result's error

- `abs_error` bounds the distance from `value` to a sign change of fn as evaluated, and is at most `2 tol(value)`.
- The true root may differ by fn's own rounding divided by |f'| (the root's conditioning), which no method can see.
- Where fn evaluates to exactly 0, `abs_error` is the requested tolerance `tol(value)`, not a bound (P4): 0 would
  claim an exact root, and in the prototype e^x - 2 evaluates to exactly 0 at a point 2.3e-17 from ln 2. An exact or
  underflowed zero is only as good as fn's evaluation: `(x + 1e10) - (1e10 + 1/3)` is 0 on an interval about 2e-6
  wide; `(x - 1e-110)^3` underflows to 0 at x = 0, 1e-110 from its root. The documentation shows both.

## Prototype results

`.scratch/roots1/proto3.py`: a throwaway Python implementation of A to F in double, with mpmath references at 60
digits. Every result is within `abs_error` of the reference, and `abs_error <= 2 tol(value)`.

| Case | Bracket | Bisection: evaluations | Brent: evaluations | Brent: actual error |
|---|---|---|---|---|
| x^2 - 2 | [1, 2] | 52 | 9 | 1.3e-16 |
| cos x - x | [0, 1] | 53 | 9 | 8.0e-17 |
| e^x - 2 | [0, 1] | 51 | 8 | 2.3e-17 (fn exactly 0) |
| x e^x - 1 | [0, 1] | 53 | 9 | 3.3e-17 (fn exactly 0) |
| x^3 (root at 0) | [-1, 2] | 3 | 3 | 0 |
| x - 1e-300 | [0, 1] | 1049 | 3 | 0 |
| ln(x / 1e-12) | [1e-20, 1] | 92 | 45 | 0 |
| (x - 1)^3 | [0, 3] | 54 | 156 | 1.1e-16 |
| x^20 - 1 | [0, 5] | 55 | 19 | 0 |
| x - 1 | [-max, max] | 1077 | 8 | 0 |
| x - 1 | [0, max] | 1076 | 7 | 0 |
| sin(x)/x - 0.5 (NaN at 0: zero first skipped) | [-1, 3] | 55 | 9 | 7.5e-18 |
| x/(e^x - 1) - 0.5 (NaN at 0) | [-1, 3] | 54 | 12 | 1.5e-16 |
| 2x - 5 denorm_min | [0, 3 denorm_min] | 3 | 3 | 0 |
| (x^2 - 2) / ((x^2 - 2)^2 + 1e-24) | [1, 2] | 52 | 58 | 1.3e-16 |
| a jump at 0 (sign function) | [-1, 2] | 1076 | 1076 | |
| 1/(x - 1/3) (a pole, returned) | [0, 1] | 54 | 55 | |

## Tests (`tests/test_roots_bracketing.cpp`, test cases `roots:...`)

- **References** from closed forms (sqrt 2, ln 2, 1e-300 as `fl(1e-300)`, 1, 0) and published constants (the Dottie
  number, OEIS A003957; the omega constant W(1), OEIS A030178), at 50 digits, cross-checked with mpmath; never from
  Numerixx.
- **Accuracy:** `|value - ref| <= abs_error + c`, with c, the conditioning term, derived per case: fn's error budget
  (an exact `+ - * /` fn: its rounding; a libm fn: a stated budget of 2 ulps) times the spacing of fn's values near the
  root, over |f'(root)|. And `abs_error <= 2 tol(value)` always.
- **Cost:** evaluations counted by a counting fn, asserted exactly only where they cannot depend on a libm: fn built
  from `+ - * /`, in double; Bisection against the derived count (C), Brent against the prototype (an independent
  implementation, as `nxx::deriv`'s were; its values in a comment). For libm fns a range below the derived count; for
  the other types no counts beyond Bisection's derivation.
- **Types:** the edge cases (A to F) in double; in float, long double, `cpp_bin_float_50` with and without expression
  templates and `cpp_dec_float_50` (radix 10), a few accuracy cases with fn built from `+ - * /` and one
  `static_assert` each.
- **Edges:** swapped and equal ends; no sign change; a NaN or infinite end, end value and interior value; negative,
  NaN, below-floor and huge tolerances; a root at 0 by zero first; fn not finite at 0 (sin(x)/x); the subnormal case
  of C; [-max, max] and [0, max]; a jump; a flat root; the three functions of E (the pole returned, the resonance's
  root found).
- **Constant evaluation:** `static_assert` on cases with exact dyadic answers, and the same cases at run time
  (exact results, so FMA contraction at run time cannot change them).
- **fn's exceptions pass through** unchanged.

## Documentation

- `roots.hpp`: Doxygen comments for every public name: the cost, the evaluation at 0, the walk to 0 and its remedy
  (B), poles (E), the conditioning caveat and the exact-zero cases (F).
- `docs/roots/design.md`: the decisions of this phase, with the prototype's table (a link into `.scratch/` would lead
  nowhere), and the primary sources.
- `README.md`: a three-line `nxx::roots` example.
- `docs/plan.md`: Phase 1's checkboxes.
- No demo yet: it would repeat the README; it comes when the module has more methods.

## Decisions by the user (2026-10-25)

P1 to P5 approved as recommended; P6: document the walk to 0, defer the fast descent. The options as they were put:

- **P1. Defaults `abs_tol = denorm_min`, `rel_tol = 2 eps`**, also the minimum (full precision at every magnitude, no
  scale assumed). Recommended. Alternative: an absolute floor on a scale of 1, as MATLAB (`eps max(|x|, 1)`) and SciPy
  (2e-12) do: cheaper near 0, but a root at 1e-12 loses its digits. Changing it later changes every caller's results,
  so it is decided now.
- **P2. Zero first** (A4), skipped where fn(0) is not finite. Recommended: it is what makes P1 affordable for a root
  at 0. Decided now, for the same reason.
- **P3. No `pole` in Phase 1:** a pole is returned as a root, documented; a local-growth rule deferred. Recommended.
  Alternative: design and test the local rule in this phase.
- **P4. `abs_error = tol(value)` where fn is exactly 0** (F). Recommended. Decided now.
- **P5. `callable` moves to `math.hpp`** for both modules. Recommended. Alternative: a copy in `roots.hpp`.
- **P6. The walk to 0** (B): document it now and defer a fast descent toward 0 (recommended: Phase 1 stays the
  textbook methods), or add the fast descent in this phase, which bounds multiprecision costs too.

## Decisions from the reviews

From the numerics review:
- The pole rule was unsound both ways (three counterexamples, E): dropped (P3).
- The overflow risk is in Brent's products, not its quotients, and an infinite product can accept a step; the
  prototype's own [0, max] case overflowed: the guard of D.
- Zero first turned a removable singularity at 0 into an error (sin(x)/x): skipped where fn(0) is not finite (A4).
- Zero first does not stop the walk to 0 for a jump at 0 or a tiny root: the walk is documented with its multiprecision
  cost, and a fast descent offered (P6). B's earlier claim that only roots fn misses at 0 pay is corrected.
- Bisection's `abs_error` as half the width could be 0 among subnormals: `max(m - a, b - m)` (C).
- Brent stops on `|c - b| <= 2 tol` rather than on a rounded half-difference, which in radix 10 could exceed the bound.
- Huge tolerances could overflow `tol`: clamped (A1).
- The conditioning term of the tests is derived per case with a libm budget, not `eps |ref|`; counts are exact only
  where no libm is involved, and Brent's count has an independent model only for double.
- The floors' rationale corrected (A1); `abs_error <= 2 tol(value)` stated; the exact-zero caveats documented (F).
- To check against Brent's book while implementing: `macheps`, the stated error bound, the worst-case count (D).

From the simplicity review:
- `pole` dropped (the same finding); the overflow guard kept as a few comparisons, without margin tests.
- No demo in Phase 1; the type matrix limited (Tests); the README example kept to three lines.
- The prototype's table goes into `design.md`, since `.scratch/` is not in the repository.
- `Criteria` gets `operator==`, as the other aggregates have.
- Kept as proposed: zero first, the floors, `non_finite` for bad tolerances, P4, moving `callable`, no cap.
- Not adopted: dropping Brent's count tests for being from a prototype. The test rule allows independent prototypes,
  and `nxx::deriv`'s tests rest on them; the counts are limited to libm-free cases instead.

From the usage sketch (5 snippets, 7 misuses):
- `non_finite` for a negative tolerance reads oddly: kept (the `nxx::deriv` precedent), documented on the enumerator.
- `abs_error` is easily read as the error of the true root: F's caveats go into the documentation and the README.
- The defaults surprise callers used to an absolute tolerance, both ways (a trace root needs nothing; a root near 0
  that fn misses needs `abs_tol`): the README shows both.
- `auto [x, e] = *r` compiles and would break when a member is added: documented (read by member name), as approved
  in the brief.
- Asked for: the evaluation count, a bracket search, fn(lower) and fn(upper) supplied: deferred (BRIEF.md 6).

## Deferred

- An iteration cap and `not_converged`: Phase 2, with the open methods.
- The fast descent toward 0 (P6): trial points `b radix^-k`, k doubling, then bisection, bounding the walk to 0 by
  the exponent's digits rather than its range; with Roots.jl's bisection of the representation as the model for
  hardware floats. Non-breaking: results stay within the same tolerance, only the points and counts change.
- `error::pole` with a local-growth rule (E). Behaviour change: a result becomes an error; documented now.
- Infinite end values accepted for their sign (fn = log x on [0, 1]): today `non_finite`; a behaviour change, an error
  becoming a result; `non_finite`'s documentation says it may happen.
- The overflow guard's cases in constant evaluation (a guard that never fires there is not tested).
- fn(lower) and fn(upper) supplied with the bracket; the residual or the evaluation count in `Root` (BRIEF.md 6).
- A demo program.

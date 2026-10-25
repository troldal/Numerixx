# Phase 5 API proposal: one-sided Ridders, and the float warning from f''

Status: approved by the user on 2026-10-24.

The user has decided:
- one-sided Ridders as studied (`ONESIDED-RIDDERS.md`, prototype in `.scratch/onesided`):
  - a `Forward<1>` base of order N;
  - the rounding floor carried through the tableau;
  - fast descent after an untrusted result;
  - an initial step of max(|x|, 1) / 32 for every N;
- the float warning for every derivative of order 2 to 4 and every method, the existing ones included, since float
  is not a good type for derivatives above the first;
- the fast restart for the existing `Ridders` too (item C).

One review round (simplicity-critic, numerics-reviewer) is folded in ("Decisions from the reviews").

## Summary

```cpp
using namespace nxx::deriv;
auto const f = [](double x) { return x < 1.0 ? std::nan("") : std::log(x); }; // undefined below 1
auto const g = [](double x) { return x > 2.0 ? std::nan("") : std::log(3.0 - x); }; // undefined above 2

auto const a = derivative<Ridders>(f, 1.0);           // non_finite: every central sample reaches below 1
auto const b = derivative<ForwardRidders>(f, 1.0);    // 1 to about 1e-13 relative, with an error estimate
auto const c = derivative<ForwardRidders, 2>(f, 1.0); // -1 to about 1e-10
auto const d = derivative<BackwardRidders>(g, 2.0);   // -1: samples at 2 and below only
auto const e = derivative<Central<2>, 2>([](float x) { return x * x; }, 1.0f); // now a warning: f'' in float
```

- There are two new methods. They are accepted wherever `Ridders` is: `derivative<M, N>(fn, x[, h])` and
  `default_step<M, N>(x)`.
- Existing double and multiprecision calls are unchanged, bit for bit, except those that restart after an untrusted
  result (item C).
- This is breaking for float: f'', f''' and f'''' in a coarse type now warn, by any method.

## Public API (`deriv.hpp`)

```cpp
namespace nxx::deriv {

    /// Ridders' extrapolation of the Forward<1> rule of order N: fn is evaluated at x and above only. With an error
    /// estimate, like Ridders. For fn undefined or different below x: x at a domain edge, the lower end of a
    /// correlation's range, a one-sided derivative at a kink.
    struct ForwardRidders {};

    /// ForwardRidders mirrored: fn is evaluated at x and below only.
    struct BackwardRidders {};

    // derivative<M, N>(fn, x, h), derivative<M, N>(fn, x) and default_step<M, N>(x) accept both, unchanged
    // otherwise; detail::method<M> becomes stencil<M> || one of the three Ridders methods.

} // namespace nxx::deriv
```

- **Two tags, not `Ridders<Side>`.**
  - Making `Ridders` a template would break every existing `derivative<Ridders>` call.
  - The two tags mirror `Forward<P>` and `Backward<P>`.
  - The names follow the naming rule: PascalCase, no digits.
- **Default initial step:** `default_step<ForwardRidders, N>(x)` and `default_step<BackwardRidders, N>(x)` are
  max(|x|, 1) / 32 for every N. This is one constant, not a per-N array, because the rules reach N h.
  - **N = 1 to 3:** in the sweep (double), s/32 and s/64 are the starts with no gross errors. s/32 is the more
    accurate of the two, by 1.4 to 2 times. s/16 has one gross error for N = 2 and 3.
  - **N = 4:** s/32 has the fewest gross errors (5, against 9 at s/16 and 34 at s/64). Its worst case is 1.8 times
    the estimate, against 12 at s/16.
- **`exact_weights`:** checked against the method's base table (below). For the one-sided methods that is the
  `Forward<1>` table of order N, whose largest integers are 1, 2, 3 and 6 for N = 1 to 4. Any type that holds
  `Central<2>`'s table holds it.
- **Cost:** fn(x) once, then N evaluations per level, with the same cap of 32 levels: at most 33, 65, 97 and 129
  evaluations for N = 1 to 4.

### A: the algorithm

**One driver on the base stencil.** `detail::base<M>` maps `Ridders` to `Central<2>` and both one-sided methods to
`Forward<1>`. It replaces the `std::conditional_t` in `largest_weight_integer` and selects the driver,
`detail::ridders<S, N>` (`ridders<N>` now). Most of what differs between the bases follows from S's weight table:

| Follows from the table of S | `Central<2>` | `Forward<1>` |
|---|---|---|
| Reach for the overflow check, q = max(2K, denominator) (`stencil_derivative`'s formula, factored out) | 2, 2, 4, 4 | 2N |
| fn(x) needed (the table has node 0) | for even N | always, once per call, kept across restarts |
| Resolution check: consecutive points of the nodes and 0 distinct | x - 2h ... x + 2h, the current checks | x ... x + N h |

Two things are chosen by S directly:

| Chosen by S | `Central<2>` | `Forward<1>` |
|---|---|---|
| Column divisor | 4^j - 1 (an even error expansion) | 2^j - 1 (every power of h) |
| Rounding floor of entry (i, j) | the floor of level i (as now) | carried through the tableau |

**The rest is shared**, with one change for all methods. A restart, after a non-finite level or after an untrusted
result, now starts at half the last step evaluated (item C); today the untrusted restart starts at h_start / 2. The
shared rules are:
- the stall of two levels;
- the untrusted test (factors 4 and 64);
- exact agreement;
- the diagonal check;
- the cap of 32 levels;
- the fallback to the last untrusted result.

**`BackwardRidders` is the `Forward<1>` driver run on the mirror.**
- It calls `ridders<Forward<1>, N>` on `t -> fn(-t)` at -x, then negates the value for odd N. The estimate is
  unchanged.
- Negation is exact, so the points fn receives are fl(x - k h), and the mirror holds by construction.
- A `Backward<1>` table would not do this: `weighted_sum` adds in node order, so that table would start its sum from
  f(x - N h), not f(x). Floating-point addition is not associative, so the results would differ.

**The levels**, with the library's conventions and its helpers split, not duplicated:
- **N = 1:** the arithmetic of `difference` is split from its evaluation. `difference` itself and the forward level
  call it, the forward level with f(x) cached. This uses the actual distance and the same guards.
- **N = 2:** the arithmetic of `second_difference` becomes a function of three points x0 < x1 < x2 and their values,
  keeping its expression ((f2 - f1) / h_hi - (f1 - f0) / h_lo) / ((h_hi + h_lo) / 2) and its guards.
  - The central level calls it with the cached value in the middle; the forward level with the cached value at x0.
  - `second_floor` is already the Lagrange-weight floor for three points.
  - The result is exact for quadratics when the offsets are exact (Sterbenz: 0 < 2h <= |x|).
- **N = 3, 4:** the current body, templated on S.
  - The nominal step and S's table.
  - The floor is sum |w_k| eps |f_k| / h^N, plus the rounding of the points, eps/2 (|x| + K h) |f'| sum |w_k| / h^N.
    K is the table's reach: 2, as now, or N.
  - The one branch is f': (f(x + h) - f(x - h)) / 2h for central, (f(x + h) - f(x)) / h for forward.

**The propagated floor** (one-sided only):
- floor(i, 0) is the level's floor;
- floor(i, j) = (2^j floor(i, j-1) + floor(i-1, j-1)) / (2^j - 1). The Richardson coefficients alternate in sign, so
  this is exactly sum |c_k| floor_k;
- it saturates at every step: `saturating_quotient(saturating_sum(saturating_product(...), ...), ...)`, and max stays
  max ("unknown");
- the arrays are built with `detail::filled`, for a T without a default constructor;
- the estimate of an entry is max(change, floor(i, j)), and the untrusted test uses the best entry's floor(i, j);
- the evidence: with the level's floor, f' was covered in 97.5% of the cases, with a worst case of 8.9 times the
  estimate and 11 of 800 poles missed. With the propagated floor, 99.4%, worst 2.6, and no pole missed;
- it is computed for every base and selected once by S. For central this is wasted arithmetic, about one tableau's
  worth, which only matters for multiprecision.

**Fast descent** (all methods):
- A one-sided base never sees the NaN beyond a domain edge. It gets only untrusted results, and halving the initial
  step once per tableau ran out of levels. sqrt and log at 1e-6 and 1e-8 were 80 to 100% off, with the error above
  the estimate.
- No step is skipped: the steps tried are still h0 / 2^k, so the resolution check and the level cap work as now.
- **The trade-off.** A restart throws away the previous tableau's last levels, so a case that the slow descent
  already handled can lose some accuracy. Covered cases lost up to 8 times: f'' of log at 0.01, one-sided, 8.5e-11
  against 7.1e-10; f' of 1/(x - 0.01) at 0, central, 3.8e-15 against 2.1e-14.
- **Over the 15 edge cases per order** (one-sided, double), fast descent did better overall:

  | Rule | Uncovered | Geometric mean error, N = 1 / 4 | Evaluations, N = 1 / 4 |
  |---|---|---|---|
  | slow descent (h_start / 2), as now | 4 to 7 | 1.6e-9 / 1.5e-3 | 341 / 1427 |
  | half the last step (proposed) | 0 | 2.0e-13 / 3.4e-6 | 290 / 1059 |
  | the last step | 0 to 2 | 2.4e-12 / 1.6e-5 | 303 / 1175 |
  | twice the last step | 2 to 4 | 5.1e-11 / 7.2e-5 | 320 / 1307 |

  - Restarts that overlap one or two levels are worse on every measure.
  - The random study is identical under all four rules.
- **Its reach**, from the default step:
  - with fn's edge at a distance x, the one-sided method stays accurate down to about x = 1e-9, which is
    2^-30 max(|x|, 1). At that distance f' is 6e-9 off and f'''' 1e-3;
  - closer than about 1e-10 to 1e-12, depending on N, the levels run out and the result is wrong;
  - in every such case of the study (sqrt and log, 1e-12 to 1e-14, N = 1 to 4), the result is flagged: abs_error is
    at least |value| / 4, mostly above |value|;
  - a step on fn's own scale works closer.

### B: the float warning from f''

- `detail::noisy<N, T>` becomes `precision_warning && N >= 2 && coarse<T>`, for every method. Today it is N == 4.
- The message: "derivatives above the first in a type with 7 or fewer significant digits keep few correct digits:
  use double, or silence this warning at the call or with NUMERIXX_NO_PRECISION_WARNING".
- The evidence, as median relative errors in float on the random study (the C++ harness, x drawn as in the
  prototypes):

  | Method | f' | f'' | f''' | f'''' |
  |---|---|---|---|---|
  | `Central<2>`, default step | 2.9e-5 | 4.4e-4 | 5.1e-3 | 1.5e-2 |
  | `Ridders` | 1.7e-6 | 1.1e-5 | 1.3e-4 | 9.6e-4 |
  | `ForwardRidders` | 9.5e-6 | 7.9e-4 | 2.3e-2 (worst 27 times its estimate) | 1.8 |

  - The header's comment on `noisy` quotes 2e-3 for Ridders' f''''. That figure came from the Python prototype's
    draws; the comment is rewritten with these numbers.
- **Breaking:** float f'' and f''' calls now warn, and with warnings as errors they fail to build. In the repository,
  each deliberate float call gets the silencing pragma at the call, as the float f'''' tests have now:
  - `tests/test_orders.cpp`: the exact float f'' cases (lines 235, 236, 419, 420), and Ridders' float f'' and f'''
    accuracy (564 to 567);
  - `tests/test_deriv.cpp`: `agrees_for_methods<3>(float)` becomes `<1>`; the default-step agreement in float for
    N >= 2 moves into a silenced block;
  - `tests/precision_warning.cpp`: `third_float` moves to the cases that must warn;
  - `docs/deriv-book/tools/generate_data.cpp`: its float data for N = 2 and 3. This is book tooling, outside the
    default build, and is handled with the book's update.
  - The demo and the book's examples have no float call above f' that is not already silenced.

### C: the same fast descent for `Ridders` (decided by the user)

The study found that central `Ridders` fails on poles close to x, where fn is finite on both sides. Measured on
1/(t - a) at 0, with |a| = 1e-3 or 1e-4:
- **f''' and f'''' at 1e-3 and 1e-4, and f'' at 1e-4:** the relative error is 1 to 11, with an estimate that does not
  cover it, after every evaluation.
- **f'' at 1e-3:** 2e-2 off, covered by its estimate of 0.57.
- **f' at 1e-4:** 6% off, with a large estimate.
- **The documented edge limit** (2^-34 for f', 2^-31 for f'' to f'''') holds today for non-finite points only.
- **With fast descent:** these cases are 1.5e-15 to 3.8e-9 off and covered, in 32 to 88 evaluations. The limit then
  holds for finite poles too.

Effect on everything else (from the prototype, which reproduces the library bit for bit):
- **Random study (2,700 cases, 800 poles):** unchanged for double N = 1 to 4 and float N = 1 to 3. Float N = 4 has
  more gross cases (214 against 88). That is noise, and now warned of (item B).
- **Existing test cases:** results change only where an untrusted restart happens.

  | Case | Now | With fast descent |
  |---|---|---|
  | f', 1/(x - 0.01) at 0 | 3.8e-15, 30 evaluations | 2.1e-14, 26 |
  | f', sin(64 pi x) at 0 | exact, 24 | exact, 20 |
  | f'' to f'''', 1/(x - 0.01) at 0 | 55, 108, 109 evaluations | the same values, with 31, 60, 61 |
  | f', Heaviside at 0 (level cap) | untrusted, 64 | untrusted (another value), 64 |
  | log near 0, the other restarts | | unchanged |

  The references in the first and last rows come from the Python prototype, which gets the same change.

## Tests

As in Phases 3 and 4, references come from closed forms and from an independent Python prototype. The C++ prototype
in `.scratch/onesided` uses the library's helpers, so it is not independent. The Python prototype
(`.scratch/phase4/riddersn.py`) is extended with the forward base, the propagated floor and fast descent. It must
reproduce the C++ results for double exactly, as Phase 4's did.

- **Polynomials** (`static_assert` and at run time, binary types), with dyadic x and a power-of-two h:
  - x^N gives f^(N) exactly at the first two levels, which agree exactly. The estimate is the propagated floor
    2F(1, 0) + F(0, 0).
  - x^(N+1) is exact after one column.
  - Float is exact at x = 3 and h = 1/32 for N = 4.
- **Accuracy and coverage** against closed forms: exp, sin, log and the Runge function for N = 1 to 4 in double.
  - Per case, from the prototype, as in Phase 4: the error is within 6 times the estimate and within 4 times the
    prototype's error, the estimate is within twice the prototype's, and the evaluations are within the prototype's
    plus 4.
  - Float: accuracy for N = 1 only, and coverage within 6 times the estimate for N = 2 (worst case 4.0). The higher
    orders in float are the noise that item B warns about.
- **Edge cases:**
  - fn NaN below x = 1: `Ridders` gives non_finite, and the one-sided method gives the derivative;
  - sqrt and log at 1e-4 and 1e-8, through fast descent;
  - a pole 1e-4 behind x;
  - a kink at x, where the one-sided method returns the one-sided derivative;
  - beyond the reach, log at 1e-12: the result is flagged (4 abs_error >= |value|). This is not a coverage test.
- **The side:** a fn that records its arguments. Every argument is >= x for `ForwardRidders` and <= x for
  `BackwardRidders`, for N = 1 to 4, including after restarts. This is the contract the tags promise. The NaN-below
  case does not prove it, since a sample below x would only cause a restart.
- **The mirror:** one run-time check that `BackwardRidders` equals `ForwardRidders` on f(-t) at -x. It holds by
  construction.
- **Failures and cost:**
  - a non-finite fn(x) gives non_finite after 1 evaluation;
  - fn finite at x only exhausts the levels for N = 1 to 4, with exactly 33, 65, 97 and 129 evaluations;
  - x near max, and h below the resolution of x, cause no evaluation beyond fn(x).
- **Multiprecision** (`cpp_bin_float_50` with and without expression templates, `cpp_dec_float_50`): coverage
  against closed forms only, without a regression baseline.
- **The warning:**
  - a float f'' must warn, one CTest compile case per overload;
  - a float f' must not;
  - `static_assert` on `noisy<N, T>` for N = 1 to 4 in float and double.
  - `noisy` does not depend on the method, so the methods are not multiplied out.
- **Item C:** the close poles above for `Ridders`, N = 1 to 4, and the changed references.
- **Central bit identity:** the prototype's harness compares 21,600 central cases against the refactored header (the
  same values, estimates and counts, outside the cases that restart after an untrusted result). This is run once, and
  the result recorded here.

## Documentation

- The header's Doxygen for the two tags, `default_step`, the overloads and the warning.
- `docs/deriv/design.md`:
  - the one-sided method and how it differs from `Ridders`;
  - the warning from f'';
  - fast descent, with the edge limit now valid for finite poles, and the one-sided reach;
  - the passing of a step on fn's scale near an edge (for example `derivative<ForwardRidders, 4>(fn, x, x / 8)`).
- `docs/deriv/DERIV_TUTORIAL.md`, the README and the demo: one example each of `ForwardRidders` at a domain edge.
- `docs/plan.md`: Phase 5 with its checkboxes.
- Then the independent Codex review, with its adopted findings recorded here.
- The book is a separate step, after the library change is committed:
  - a section in the Ridders chapter;
  - Part II's reference;
  - the warning's new scope;
  - the version references.

Caveats to document:
- the one-sided method is less accurate than central: by 16, 260, 1300 and 3000 times for N = 1 to 4;
- 2 of 400 poles ahead of x and inside the initial span are uncovered for N = 4, by up to 4.9 times;
- the reach above;
- the float accuracy table above.

## Decisions from the reviews

From simplicity-critic:
- **The driver takes the base stencil, not the method.** The reach, the need for fn(x) and the resolution check
  follow from the table, so only the column divisor and the floor are chosen by S.
- **`BackwardRidders` runs the forward driver on the mirror**, not a third variant.
- **Helpers are split, not duplicated:** the arithmetic of `difference` and `second_difference`, and one N >= 3
  level.
- **One restart rule** for every restart and every method.
- **A smaller test plan:** warning cases per overload, not per method; float accuracy for N = 1 only; no
  multiprecision regression baseline; the cost test merged with the exhausted-levels test; a new test that fn is
  sampled on one side only.

From numerics-reviewer:
- **The mirror** cannot use a `Backward<1>` table and stay bit for bit (counterexample: N = 3 with values 1, 1, 0,
  2^54, where the order of addition decides a tie). That is solved by construction (above).
- **The one-sided reach is now stated:** measured, flagged beyond, and tested once.
- **The numbers are for the configuration that ships:** fast descent with the last fallback. They are identical in
  double to the variant with the best fallback; in float N = 2 there are 12 gross cases, not 5.
- **The cost of fast descent in individual cases is stated**, with the alternatives measured.
- **Float tests** are limited to the orders where they test something.
- **The propagated floor saturates** at every step, and its arrays need no default constructor.
- **The N = 2 level** keeps `second_difference`'s division by (h_hi + h_lo) / 2. A division by x2 - x0 after
  multiplying by 2 could overflow first (0, 0 and 1.2e308 at h = 1).
- **Corrections:** item C's description of the failures, the rationale for s/32, and the example for
  `BackwardRidders`.
- **A harness reference had the wrong sign** (central f'', 1/(x - 0.01)). Corrected, the row's error is 1.3e-12, as
  the existing test's reference.

## Changes during implementation

- **Central bit identity, measured.** The study's 28,000 central cases (nine functions at 300 points and 800 poles,
  N = 1 to 4, double and float) were dumped from the old and the new header. 5,233 changed: every one of the 800 poles
  for N = 2 to 4 and 110 for N = 1, in both types, and 213 float f'''' cases, all of them after an untrusted result.
  No double result on the nine functions changed. The C++ prototype with only the restart changed reproduces the new
  header on all 28,000 cases, so the restart is the only change.
- **The Python prototype** (`.scratch/onesided/ridders5.py`, numpy, from the textbook rules) reproduces the new
  header for double bit for bit on 28,000 cases, both bases: values, estimates and evaluations. Its first run differed
  only for x^7 - x, where it used `t ** 7` and the C++ repeated multiplication; numpy's exp, sin, log, atan and sqrt
  round as the C++ libraries here. Every new test reference comes from it.
- **fn finite at x only** costs 33, 65, 33 and 33 evaluations, not 33 to 129: a level stops at its first non-finite
  value, except the f'' level, which evaluates both of its points first. The cap of 1 + 32 N is reached with every
  evaluation made, and is tested on log at 1e-12, beyond the reach.
- **The float polynomial case** is x^4 at 1 with h = 1/8, not at 3 with h = 1/32 as in the review: there the node
  3 + 1/32 gives 97^4 / 2^20, which needs 27 bits.
- **The warning's tests** silence deliberate float calls at the call (`tests/silence.hpp`); a call inside a test
  template cannot be silenced from the line that instantiates it. The float f'' and f''' accuracy tests of `Ridders`
  call the method directly in a silenced block, and share only the checks.
- **The mirror is a struct, not a lambda** (`detail::Mirrored`, holding a pointer to fn): MSVC 19.51 rejects, in
  constant evaluation, a lambda that captures a function reference by reference ("a read of an uninitialized
  symbol"), which a static_assert of `BackwardRidders` on a plain function hit; GCC, Clang and clang-cl accept it.
  Reproduced in isolation (`.scratch/onesided/msvc/repro.cpp`).
- **Callables:** both methods are tested with a member-function pointer, a callable that rejects rvalues, a type
  without a default constructor and one whose `numeric_limits` is written member by member, and with the
  `exact_weights` constraint on a 2-digit type (f''' accepted, f'''' rejected, as for `Ridders`).

## Independent review (OpenAI Codex, f20696f)

No Phase 5 defect. It confirmed the column factors, the propagated floor (indexing, saturation, no stale entries
after a restart), the forward levels and their floors, the reach and resolution checks, the mirror (signed zero, the
sign of odd orders, the unchanged estimate), the restart rule, the callables and the warning's overload selection.

One inherited limitation, reproduced: a final difference quotient that overflows is formed before it is checked, so
GCC rejects the constant evaluation, where at run time the level counts as non-finite and the method restarts below
it. `derivative<ForwardRidders>` of (max t)(64 t) at 0 returns 0 at run time but does not compile as a constant
expression. It predates Phase 5 and is documented (the header, `design.md`). Left as it is (user decision): closing
it needs a conservative check that would turn stencil results between about max / 2 and max into non_finite.

## Not adopted

- **The propagated floor for `Ridders`:** it raised coverage slightly (f' from 99.0% to 99.3%; the worst f'' case
  from 3.4 to 2.5 times), but it changes nearly every central estimate, every Phase 3 and 4 reference, and the
  book's figures.
- **Returning the untrusted result with the smallest estimate** when the levels run out: it helped only float f''
  (12 gross cases fell to 5), which now warns.
- **Divided differences on the actual nodes for N = 3, 4:** 1.5 to 2 times more accurate in double, but worse
  coverage in float, and a departure from the central base's convention.
- **A `Forward<2>` base:** no more accurate than `Forward<1>` in the first study, at a higher cost.
- **Restarts that overlap one or two levels:** worse than half the last step on every measure (above).

## Deferred

None of these is a breaking change, so none needs a hook now.
- the propagated floor for `Ridders`;
- a `Forward<2>` base;
- divided differences for N = 3, 4;
- the best-estimate fallback;
- a `Backward<1>` table instead of the mirror, should a type with asymmetric rounding ever need it.

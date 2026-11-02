# Phase 3b API proposal: TanhSinh through the method protocol

Status: approved by the user on 2026-11-02, Q1 to Q5 as recommended (draft 2). Draft 1's review round (simplicity, usage, numerics) is folded in
("Decisions from the reviews"), with the prototype rebuilt on its questions (D); the numerics review rerun on draft 2
is folded in too.

The decisions it builds on:
- **`PHASE3_API.md`** (approved 2026-11-01, Q1 to Q10 as recommended): TanhSinh's numerics, A to F.
- **`PHASE3A_API.md`** (approved 2026-11-02, Q1 to Q15 as recommended), section G: TanhSinh is an iterative method
  of the protocol, and `integrate` owns the loop, the stop test and the order of the ends.

This document adds the declaration, the steps, the state, the points `PHASE3_API.md` leaves implicit (settled by the
prototype), the order of the work, and the protocol's tests. It reopens none of `PHASE3_API.md`'s numerics.

## Summary

```cpp
#include <numerixx/quadrature.hpp>
namespace qd = nxx::quadrature;

// Endpoint singularities: the double-exponential rule samples ever closer to the ends, never at them.
auto const r = qd::integrate<qd::TanhSinh>([](double x) { return std::log(x); }, 0.0, 1.0);      // -1
auto const s = qd::integrate<qd::TanhSinh>([](double x) { return 1 / std::sqrt(x); }, 0.0, 1.0); // 2

// High precision in a multiprecision type: a few hundred evaluations.
auto const m = qd::integrate<qd::TanhSinh>([](mp const &x) -> mp { return fn(x); }, mp(0), mp(1),
                                           {.rel_tol = mp("1e-45")});

// Each level reported to an observer: PHASE3A_API.md's Printer.
auto const o = qd::integrate<qd::TanhSinh>(qd::observer(printer), fn, 0.0, 1.0);
```

- `TanhSinh` is an iterative method: `start` is level 0, `advance` the next level, `settle` nothing. The existing
  forms of `integrate` take it, named or as a value, observed or not. No form of `integrate` is added or changed.
- `detail::exp` and the node tables in pieces, as `PHASE3_API.md` B approved.
- The default stays `GaussKronrod<15>`. Everything constexpr; nothing throws.
- The fallback from TanhSinh to the default is Q4.

## Public API (additions to `quadrature.hpp`)

```cpp
namespace nxx::quadrature {

    /// The tanh-sinh rule (Takahasi and Mori 1974): the trapezoidal rule after x = tanh(sinh t), the step halved
    /// level by level. (The full comment: "Comments", below.)
    struct TanhSinh
    {
        /// The protocol's start: level 0.
        template<real T, typename F>
        constexpr std::expected<detail::TanhSinhState<T>, Failure<T>> start(F &fn, T const &a, T const &b) const;

        /// The protocol's step: the next level.
        template<real T, typename F>
        constexpr std::expected<void, Failure<T>> advance(F &fn, detail::TanhSinhState<T> &s) const;

        /// The protocol's settle: nothing to refine.
        template<real T>
        constexpr std::expected<void, Failure<T>> settle(detail::TanhSinhState<T> &s) const;
    };
}
```

- **The failures:** `start` and `advance` return `error::non_finite`, with no best, where:
  - fn is not finite at the centre, or at a point that is not in a tail (`PHASE3_API.md` D);
  - a value, an L1 sum or the estimate overflows while fn is finite;
  - in constant evaluation, Phase 1's checks (`checked_term`, `fits`, `compensated`) refuse a term of a finite
    value, or a sum (the tail rule holds there too).

  The method never returns `not_converged` itself; the loop does (B).
- **The state** is `detail::TanhSinhState<T>` (C): independent of fn's type, of fixed size, no allocation.

### Comments

`PHASE3_API.md`'s "Comments" stand, plus:
- **`integrate`, choosing a method:** TanhSinh fails *usually* honestly on interior features (`PHASE3_API.md` F: 5 and
  3 false successes in 250 at 1e-6 and the default). Cost: TanhSinh's minimum, 5 levels, is about 70 to 140
  evaluations against the default's 150; on smooth fn the two are comparable (median ratio 0.8 to 1.5), with an
  endpoint singularity TanhSinh is 3 to 9 times cheaper; below about 100 evaluations, `GaussLegendre<N>`. The
  default is not TanhSinh because of interior features. The fallback: Q4.
- **`TanhSinh`:**
  - its steps: report i is level i; level 0 has no current; a later report without one means the estimate is
    unbounded; `ready` from report 4, `exhausted` only at report 9; the evaluations roughly double per level (e^x on
    [0, 1]: 10, 18, 33, 62, 121 in all after each level). A failure in `start` (`non_finite` at level 0, as for sqrt(x - 0.5)) gives no
    report: record fn's arguments to find it. A run to the cap (1/x: 5,159 evaluations) can be stopped by an
    observer that throws.
  - a singularity at an end other than 0: substitute u = b - x (or u = x - a) and write fn in u, as fn(b - u)
    rounds b - u and loses what the substitution gains;
  - a success's `abs_error` can exceed the request where a rounding cut's sliver dominates the floor: 1 on
    [1, 1 + 1e-9] ends ok with `abs_error` 2.0e-7 relative: the two slivers, each under half an ulp of its end, so up
    to an ulp over the width;
  - compile time: every piece of the tables is instantiated wherever `integrate<TanhSinh>` is, 0.5 to 1.5 s per
    translation unit; instantiate it in one;
  - constant evaluation: fn rational or built from one's own constexpr functions (`std::sqrt` and the like are not
    constexpr on Clang and MSVC); about 250 constexpr steps per evaluation on Clang with a tabled T; for a T without
    a table, the generation of each level reached as well (C);
  - high precision: the first call in a T without a table builds the pieces it reaches, about 14 us per node in
    `cpp_bin_float_50` (36 us in `cpp_dec_float_50`), 154 ms (343 ms) for all; a type that changes its precision at
    run time sees its first tables.
- **`error::not_converged`:** an absent best can also mean convergence too slow for T's range (x^-0.999 on [0, 1]):
  substitute x = u^k with fn rewritten in u (k u^(k(1 - alpha) - 1); fn(u^k) underflows and gains nothing), or use a
  type with a wider exponent range (80-bit long double gives a best, `cpp_bin_float_50` converges); the default is no
  better (its best 509 +- 0.91 for 1,000).
- **`GaussKronrod`:** "substitute ... until TanhSinh" becomes "substitute ..., or use TanhSinh".
- **Comments that would go stale,** updated: the file's @brief and its paragraph on the iterative methods,
  `iterative_method`'s first sentence, `integrate`'s `@tparam M` lists, `Step::attainable`, `Step::ready`,
  `Step::exhausted` and `Criteria` (the floor that can exceed the request).

## Behaviour

### A. What TanhSinh receives

`integrate` validates the criteria and the ends, returns {0, 0} for a == b without calling fn, and passes a < b,
negating the result for b < a (`PHASE3A_API.md` A). This replaces only that part of `PHASE3_API.md` A; Phase 1's
`Interval` with its scale order, the constant-evaluation refusals, the points and the centre not strictly inside
stand.

### B. The steps

| | `start`: level 0 | `advance`: level k = 1 to 9 |
|---|---|---|
| fn evaluated at | the centre, then per side t = 1, 2, ... up to and including the first node cut by a non-finite value; a node not strictly inside is not evaluated (at most 13 in double) | per side, the odd multiples of 2^-k below its outermost kept node, and the next one if it is in the table, before the side's cut and strictly inside |
| `current` | none: level 0 has no estimate | {V_k, E_k}, or none where E_k is unbounded |
| `attainable` | 2 F_0 | 2 F_k |
| `ready` | no | k >= 4 |
| `exhausted` | only where no point is strictly inside: then fn is never called | k = 9 |

`trusted` is always yes. V_k, F_k and E_k are `PHASE3_API.md` C and E's, with D3 below; F_k includes the slivers.
The loop then reproduces `PHASE3_API.md` E's termination: the same test at the same levels, so the same results.

### C. The state

`detail::TanhSinhState<T>` holds no table and no reference to fn:
- **The level:** k, h = 2^-k, the raw sum R_k and L1_k as plain T (each level's Neumaier `Sum` is local, as in
  `RombergState`), and V_(k-1), the previous reported value.
- **Per side:**
  - the cut: the first node not sampled, or none;
  - its kind (below);
  - the outermost kept point, which starts as the centre: its node, point, value and |term|;
  - the kept point with the largest t below it (Tp of the geometric end term);
  - the kept point with the largest t farther from the end (the power law's second point);
  - W, the table's weights beyond the outermost kept node.
- **The `Step`.**

**The kinds of a side and how they change:** the table's end (uncut), rounding, pruning, a tail. Table end becomes
any kind; rounding becomes rounding, pruning or a tail; pruning and a tail become each other. Pruning and a tail never
become rounding: their cut node was strictly inside, so every node before it is. The prototype shows table end, then
rounding, then a tail on x^-0.99 over [0, 1e-40].

**Each step needs only its own level's nodes:**
- At most one new node lies beyond a side's outermost kept node, for every kind (by induction on the levels).
- W_k = W_(k-1) plus level k's weights beyond the outermost kept node, added one by one in increasing t after the
  side's sweep, at every level while the side is uncut or cut by rounding (a pruned or tail side never needs it
  again). No node of an earlier level lies between the old and the new outermost node, so this is the suffix sum of
  `PHASE3_API.md` E, term for term: the prototype finds it equal bit for bit wherever it is used in the counted cases.
- Each step reads its level's pieces once, for both sides and W:
  - float, double and long double: from the compile-time table;
  - another T at run time: from the run-time cache;
  - another T in constant evaluation: generated as a temporary.

So `PHASE3A_API.md` G's question (copy the tables into the state, or regenerate them) does not arise (Q1).

**Constant evaluation for a T without a table** generates each level reached in the evaluation that calls
`integrate`, at about 275 constexpr steps per node of a double-like table (Clang; `PHASE3_API.md` B), on top of the
run, within one top-level budget:

| levels | nodes | steps for the tables |
|---|---|---|
| 0 to 4 | 106 | 29k |
| 0 to 5 | 211 | 58k |
| 0 to 7 | 841 | 231k |
| 0 to 9 | 3,361 | 924k |

A run past level 7 likely exceeds Clang's default of 1,048,576; a class type costs more per node. Documented, not
engineered; steps 2 and 3 measure a constexpr wrapper of double at level 4 and to the cap.

### D. What `PHASE3_API.md` leaves implicit, settled by the prototype

The prototype (Q3) is rebuilt in double from `PHASE3_API.md` A to E and this section, through the protocol: `start`,
`advance` and `PHASE3A_API.md` B's loop. It reads Phase 3's dump of the C++ double table (bit for bit the spike's).
On the 38 counted cases of Phase 3 it reproduces Phase 3's prototype bit for bit (status, value, `abs_error`,
evaluations), and that prototype matched the pre-3a C++ bit for bit.

1. **A side that reaches the table's end uncut gets the geometric end term,** as after a non-finite value.
   - **Why it arises:** x^-0.99 in double overflows only below max^(-1/0.99) = 4.28e-312, a subnormal, so fn is
     finite at every node of a table that ends at min(). The side near 0 is never cut, and the mass below the table
     is 100 min^0.01 = 0.0838.
   - **Soundness:** for x^-alpha the terms are log-concave beyond t_e >= 0.88 wherever they decrease there, so
     h Te / (1 - r) bounds the mass beyond the table; where they do not decrease, the term is unbounded.
   - **x^-0.99 on [0, 1] at rel_tol 1e-4:**

     | end term at the table's end | result | error |
     |---|---|---|
     | geometric (D1) | `not_converged`, best 99.9164 +- 0.0993 | 0.0836 |
     | 0, as after pruning | ok +- 0.0084 (325 evaluations) | 0.0751: false success |
     | h Te (Bailey) | ok +- 0.0069 | 0.0818: false success |

     x^-0.97 at the default: ok +- 4.3e-8 with an error of 9.4e-9 (geometric); with 0, ok +- 5.7e-9, a false success.
     x^-(127/128) (a counted case, 7 square roots): best 127.496 +- 0.625, error 0.504.
   - **`PHASE3_API.md` E's 0.077** came from a table into the subnormals: with it, the same run cuts that side at
     a tail, best 99.9231 +- 0.0911, error 0.0769. Moving the table's end to min() made it a table-end case; D1 keeps
     that move from changing the status (the best and the count change slightly: 99.9164 +- 0.0993 in 5,159
     evaluations against 5,168).
2. **Pruning,** as Phase 3's prototype and the pre-3a C++ do:
   - **Level 0:** after both sweeps, the threshold is eps times the L1 of the whole sweep (the centre and every
     finite term of both sides, h = 1). Each side drops its run of outermost terms below it, from the outside in, and
     is cut at the innermost dropped node (kind pruning, whatever the kind before). A small term with a larger one
     beyond it is kept. The reading "cut at the first small term going outward" is a false success on
     (x - x2)/sqrt x over [0, 1], x2 the left point at t = 2 (its term exactly 0): ok +- 4.1e-9, error 2.5e-5. D2's
     reading gives the exact value in 129 evaluations.
   - **V_0 and L1_0** are the sums of the kept terms only; F_0 follows from L1_0.
   - **Level k:** a new node beyond the side's outermost kept node is pruned where |h_k term| < eps L1_(k-1) / 2,
     eps times the running L1 at the start of the level. It is known before the level starts, the same for both
     sides and independent of the order of evaluation (Q5).
   - **The order:** level 0 evaluates the centre, the right side outward, then the left side. Level k sweeps the
     right side, then the left, each from the centre outward, so a side's new node beyond comes last. Each level's
     new terms go into a local Neumaier sum in that order, and R_k = R_(k-1)/2 + its value; L1 likewise.
3. **The raw sum and the slivers:** R_k = R_(k-1)/2 + the new terms holds no sliver. The reported value is
   V_k = R_k + the slivers term(h_k W, y_e), and |V_k - V_(k-1)| uses the previous reported value. Draft 1, which
   kept the slivers in the running sum, counts each again at the next level: 1 on [1, 1 + 1e-9] then errs by 6.8e-8
   relative (covered by its floor), and 1 on [1000, 1001] at rel_tol 0 by 5.3e-14; with D3 both are exact.

### E. Observing

`integrate<qd::TanhSinh>(qd::observer(o), fn, a, b[, criteria])` reports level 0 and then each level, as B gives them.
Level 0 has no current: `PHASE3A_API.md`'s Printer shows a dash. The evaluations are counted by the caller. Numbers
specific to the method (the level, each side's cut) stay deferred (`PHASE3A_API.md` Q14).

## The work, in three steps (Q2)

Each step is committed and pushed, and checked on the five toolchains; the user approves the next.
1. **`detail::exp`** (internal): `PHASE3_API.md` B and Q10, with its tests.
2. **The tables in pieces** (internal): `PHASE3_API.md` B and Q3, with its tests; the largest piece's share of each
   toolchain's constexpr budget and the cost per node in constant evaluation for a T without a table (C) recorded in
   `design.md`. The prototype's table is then dumped from these tables.
3. **`TanhSinh`:** the struct and the comments above, the tests below, `design.md`, the reference page and the demo.

Then the phase's independent review.

## Tests

**`tests/test_quadrature_tanh_sinh.cpp`** (test cases `quadrature:tanh_sinh:*`): `PHASE3_API.md`'s list, plus:
- **The steps,** observed in double:
  - level 0 has no current and is not ready; `ready` from level 4;
  - a run to the cap reports 10 steps, and only the last is exhausted;
  - with no point strictly inside: one report, exhausted, and fn never called.
- **Counts,** exactly in double per level, through the observer, against the prototype: the count cases of
  `PHASE3_API.md`, which include D1's x^-(127/128).
- **D's cases:** (x - x2)/sqrt x ok and exact (D2); x^-0.99 on [0, 1e-40] at 1e-4 `not_converged` with an honest
  best; 1 on [1, 1 + 1e-9] (the reference (1.0 + 1e-9) - 1.0) and on [1000, 1001] exact (D3). The tests assert the
  results; the kinds of the cuts are not reported, so they stay the prototype's evidence.
- **The substitution:** (1 - x)^-0.9 as approved, and 1/sqrt(1 - x^2) on [-1, 1], split at 0, each half as the
  integral of 1/sqrt(u (2 - u)) over [0, 1]: pi.
- **The fallback pattern of Q4.**
- **Constant evaluation:** 1/x on [0, 1] to the cap on all five toolchains.

**`tests/test_quadrature_protocol.cpp`:** TanhSinh joins the lines Romberg is on: `iterative_method` (one
`static_assert`, double); observed equals unobserved, in constant evaluation too; the last report equals the result
for a < b; named and as a value; reversed ends exactly negated.

## Documentation

- `PHASE3_API.md`'s list: the comments, `design.md` (the sources, the numbers the decisions cite, D's evidence), the
  reference page (a `doxygenstruct` for `TanhSinh`), the demo.
- `docs/plan.md`: Phase 3b's checkboxes, one per step.
- The prototype in `docs/quadrature/tanh-sinh-prototype/` (Q3).

## Decisions for the user

- **Q1. The state holds no table** (C): each step reads its level's pieces once, W is kept incrementally.
  Recommended: yes. Equal bit for bit to the suffix sum; nothing copied.
  - Alternative: copy the pieces into the state (up to 3,361 node pairs per call in double).
  - Alternative: regenerate earlier levels at each step (exp per node again).
- **Q2. Three steps** (exp, the tables, the method), each checked on CI before the next. Recommended: yes.
  Alternative: one step.
- **Q3. The prototype kept,** in double only, in `docs/quadrature/tanh-sinh-prototype/`: `ts.py` (about 350 lines
  without the review's variants, Python floats, no mpmath) and the table dump from step 2. Recommended: yes. The counts are the only check of the
  cut and pruning logic, and D rests on it.
  - Phase 3's prototype survives in an old session's scratchpad, contrary to draft 1, but it covers five arithmetics
    and the battery, and is not written to the protocol.
  - Alternative: kept in a scratchpad again.
- **Q4. The `or_else` fallback** (revisits `PHASE3_API.md` Q9's). The approved one, any failure retried with the
  default, turns TanhSinh's honest failures on endpoint singularities into the default's false successes.
  On 19 cases (TanhSinh by the prototype, the default by today's C++, double):

  | pattern | false successes | dishonest bests returned | interior features solved |
  |---|---|---|---|
  | P0, approved: any failure to the default | 2 | 4 | 8 of 8 |
  | P1, the usage sketch's: `not_converged` only; both fail: the smaller estimate | 2 | 4 | 8 of 8 |
  | P2, the code below | 1 | 0 | 8 of 8 |

  - P2 calls the default only where TanhSinh is `not_converged` with a best. It returns the default's success only
    if it agrees with that best within the two estimates, and otherwise TanhSinh's failure.
  - The false successes are x^-0.99 at 1e-4 (the default ok +- 0.0099, error 0.53; P2 rejects it) and (1 - x)^-0.9
    at 1e-2 (the default ok +- 0.092, error 0.46).
  - The dishonest bests are the default's: on (1 - x)^-0.9 at the default and at 1e-3, x^-0.999 and 1/x.
  - P2's limit (numerics review): the guard rejects a default contradicted by TanhSinh's best, not an error both
    share. An accepted value is bounded only by 2 e_t + e_k. At an end other than 0 both methods lose the same mass
    inside an ulp of the end (0.25 for (1 - x)^-0.9), so the two agree and both are wrong; the remedy there is the
    substitution, which the comment names.

  Recommended: P2, as a tested comment that states its limit:

  ```cpp
  auto const t = qd::integrate<qd::TanhSinh>(g, a, b);
  auto const r = t.or_else([&](qd::Failure<double> const &f) -> decltype(t) {
      if (f.reason != qd::error::not_converged || !f.best) return t; // non_finite, or no estimate at all
      auto const k = qd::integrate(g, a, b);                          // interior features: the default
      if (k && std::abs(k->value - f.best->value) <= f.best->abs_error + k->abs_error) return k;
      return t;
  });
  ```

  - Alternative: no fallback in the comments, only the choice of method, saying that the default is no remedy
    for TanhSinh's failure at an endpoint.
  - Alternative: P0 or P1 (above).
- **Q5. The level-k pruning threshold** eps L1_(k-1) / 2 (D2), as Phase 3's prototype and the pre-3a C++, on which
  `PHASE3_API.md` F's figures were measured. Recommended: yes.
  - The numerical argument (numerics review, for `design.md`): with a log-concave tail of decay rate lambda, the
    mass dropped beyond a pruned node is below eps L1 / (2 h lambda). Where the weights' decay causes the drop,
    lambda is about 31 at h = 2^-9, so at most about 8 eps L1 per side, inside the floor of 50 eps L1. For a fn about
    1e-13 L1 near an end, pruned at t of 1 to 1.6 (lambda 2.3 to 4), levels 8 and 9 drop 56 to 111 eps L1 per side,
    above the floor. That is harmless, since a run reaches level 8 only with E far above F, and the alternative
    doubles every figure.
  - Alternative (the simplicity review's): eps L1_(k-1). On 48 cases it changes no status (no estimate falls below
    its error on the 10 with a reference), and it saves up to 63 evaluations a case (1/x: 5,096 instead of 5,159),
    0.4 % in all; the count tests would then differ from Phase 3's recorded ones.

## Decisions from the reviews

From the simplicity review:
- Adopted: "Comments" lists only what is new, the member comments shortened, the `state_t` sentence dropped; draft
  1's Q4 and Q5 dropped as decisions (D1 a clarification, the constant-evaluation cost documented); no step adding the
  printing example's output; the tests merged into `test_quadrature_protocol.cpp`'s lines and one `static_assert`;
  the prototype in double only, reading the C++ table (Q3); the slivers' share of the floor dropped from the state,
  and the raw sum (D3); the pruning threshold named (D2); the trims (A, B, the steps' tests, "No Raw Loops",
  Deferred).
- Not adopted: eps L1_(k-1) as the level-k threshold (Q5).

From the usage sketch:
- Adopted: the fallback reproduced, and a question (Q4), with a test of the chosen pattern; "usually honestly"; the
  absent best of x^-0.999; both substitutions, the warning on fn(b - u) and the arcsine example (as a test); a
  success's `abs_error` above the request, confirmed (1.5e-8 requested, 2.0e-7 returned on [1, 1 + 1e-9]), in the
  comments of `Criteria`, `Step::attainable` and `TanhSinh`; the observing comments, the `Printer`, the advice for a
  failure in `start` and the throwing observer; the cost guidance, the compile time, "rational fn or your own
  constexpr functions", the steps per evaluation, the constexpr 1/x test on all five toolchains; "high precision"
  rather than "full precision", `-> mp` kept, the first call's cost; the comments that would go stale; the state
  documented as independent of fn's type.
- Corrected: "for the fewest evaluations on a smooth fn, TanhSinh": the battery's 2.7 times comes from the endpoint
  cases; on smooth fn the median ratio is 0.8 to 1.5.
- Deferred: a `SplitAt<Inner>` combinator in the guide, to the module's documentation pages.

From the numerics review (each settled or confirmed by the prototype, D):
- Adopted: the raw sum (D3: draft 1's reading errs by 6.8e-8 relative on [1, 1 + 1e-9], the review estimated
  2.2e-7); D1's evidence replaced (fn finite down to 4.28e-312, the mass 0.0838, the old table a tail cut); the four
  kinds and their changes, W kept incrementally (C; table end, rounding, tail observed on [0, 1e-40]);
  the four pruning details pinned (D2, the false success of the other level-0 reading reproduced); the
  constant-evaluation table (C); A's scope; the centre as each side's innermost kept point (C); the evaluation rule
  per level (B), and "the last report equals the result" for a < b.
- Confirmed: B's fields and the loop's equivalence; level 0's 13 evaluations at most; 3,361 nodes (7, 7, 13, 26, 53,
  105, 210, 420, 840, 1,680); at most one new node beyond; the incremental W (bit for bit, not only up to rounding);
  D2's halving.
- Draft 2's rerun, adopted (each reproduced or derived):
  - Q4: P2's guard misses an error both methods share. Reproduced: (1 - x)^-0.9 at 1e-2 is a false success under
    every pattern. The table now has 19 cases, describes P2 by its code and states its limit.
  - The comparison with Phase 3's prototype hid an old best wherever the new one had none. Made strict, it is
    still 38 of 38.
  - W kept only while a side is uncut or rounding: the same results, and one loop fewer per level in constant
    evaluation.
  - Q5's figures and its numerical argument; D1 keeps the status, not the best.
  - The x^-0.999 advice (fn rewritten in u; a wider exponent range).
  - The wording: cumulative counts, a term of a finite value, the two slivers up to an ulp, the test's reference
    (1.0 + 1e-9) - 1.0, and the kinds of the cuts as prototype evidence only.
  - A prototype line labelled Bailey that ran the geometric term was dropped.
- Noted, not raised against 3b: a node where fn is exactly 0 can turn a side's cut into pruning and drop its end term
  (`PHASE3_API.md` C's approved rule).

## Deferred

As `PHASE3_API.md` (infinite intervals, fn of the distance to the end, exp as a public helper, a stagnation rule) and
`PHASE3A_API.md` (numbers specific to the method in the reports).

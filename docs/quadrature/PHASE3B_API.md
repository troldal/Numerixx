# Phase 3b API proposal: TanhSinh through the method protocol

Status: draft 1 (2026-11-02), for review.

The decisions it builds on:
- **`PHASE3_API.md`** (approved 2026-11-01, Q1 to Q10 as recommended): TanhSinh's numerics. These are the node
  formula with c = 1, the levels, the tables in pieces, exp in T, the sums, pruning, the tails, the end terms, the
  estimate and the cap.
- **`PHASE3A_API.md`** (approved 2026-11-02, Q1 to Q15 as recommended), section G: TanhSinh becomes an iterative
  method of the protocol. `integrate` owns the loop and the stop test, and the ends arrive in order.

This document is the specification the implementation follows. It adds:
- the public declaration (as code);
- TanhSinh's steps and what each reports;
- its state, which answers the question G left open (copy the tables into the state or regenerate them);
- two points `PHASE3_API.md` leaves implicit;
- the order of the work;
- the tests and documents with the protocol's additions.

It reopens none of `PHASE3_API.md`'s numerics or of its decisions Q1 to Q10.

## Summary

```cpp
#include <numerixx/quadrature.hpp>
namespace qd = nxx::quadrature;

// Endpoint singularities: the double-exponential rule samples ever closer to the ends, never at them.
auto const r = qd::integrate<qd::TanhSinh>([](double x) { return std::log(x); }, 0.0, 1.0);      // -1
auto const s = qd::integrate<qd::TanhSinh>([](double x) { return 1 / std::sqrt(x); }, 0.0, 1.0); // 2

// Smooth fn at full precision in a multiprecision type: a few hundred evaluations.
auto const m = qd::integrate<qd::TanhSinh>(fn, mp(0), mp(1), {.rel_tol = mp("1e-45")});

// Each level reported to an observer (Phase 3a).
auto const o = qd::integrate<qd::TanhSinh>(qd::observer(print), fn, 0.0, 1.0);

// Kinks, steps and narrow interior peaks: the default, GaussKronrod<15>, which subdivides. A fallback:
auto const k = qd::integrate<qd::TanhSinh>(g, 0.0, 1.0).or_else([&](auto const &) { return qd::integrate(g, 0.0, 1.0); });
```

- `TanhSinh` is an iterative method: `start` is level 0, `advance` the next level and `settle` nothing. The existing
  forms of `integrate` take it, named or as a value, observed or not. No form of `integrate` is added or changed.
- `detail::exp` and the node tables in pieces, as `PHASE3_API.md` B approved.
- The default stays `GaussKronrod<15>`. Everything constexpr; nothing throws.

## Public API (additions to `quadrature.hpp`)

```cpp
namespace nxx::quadrature {

    /// The tanh-sinh rule (Takahasi and Mori 1974): the trapezoidal rule after x = tanh(sinh t), the step halved
    /// level by level; an iterative method, one step per level. (The full comment: "Comments", below.)
    struct TanhSinh
    {
        /// Level 0: the centre, then each side outward at t = 1, 2, ... until it is cut. No current result yet.
        template<real T, typename F>
        constexpr std::expected<detail::TanhSinhState<T>, Failure<T>> start(F &fn, T const &a, T const &b) const;

        /// The next level k: the odd multiples of 2^-k inside each side's kept extent, and at most one node beyond.
        template<real T, typename F>
        constexpr std::expected<void, Failure<T>> advance(F &fn, detail::TanhSinhState<T> &s) const;

        /// Nothing to refine.
        template<real T>
        constexpr std::expected<void, Failure<T>> settle(detail::TanhSinhState<T> &s) const;
    };
}
```

- **The members** are documented as Romberg's are: each says which part of the protocol it is.
- **The failures:**
  - `start` and `advance` return `error::non_finite`, with no best, in three cases: fn is not finite at the centre;
    fn is not finite at a point that is not in a tail (`PHASE3_API.md` D); or a value, an L1 sum or the estimate
    overflows while fn is finite.
  - The method never returns `not_converged` itself; the loop does (B).
- **The state** is `detail::TanhSinhState<T>` (C), reachable as `state_t<TanhSinh, T, F>` for a combinator.
- **Comments.** As `PHASE3_API.md` "Comments" (the usage sketch), with three of Phase 3a's comments updated:
  - `integrate`: choosing a method.
    - An endpoint singularity, or full precision: `TanhSinh`.
    - Kinks, steps, peaks or an unknown fn: the default, which fails less often. `TanhSinh` fails honestly on
      interior features, with `not_converged`.
    - fn too expensive for 150 evaluations: `GaussLegendre<N>`.
    - The `or_else` fallback above (tested).
  - `TanhSinh`:
    - the minimum cost: 5 levels, about 70 to 140 evaluations;
    - a singularity at an end other than 0, and the substitution u = b - x;
    - constant evaluation: fn must be built from `+ - * /`; the cost of the tables per translation unit; and, for a
      T without a table, the cost of generating them in the same evaluation (Q4).
  - `error::not_converged`: TanhSinh's cause, the cap. An absent best means one of two things. Either the estimate
    is unbounded (the terms grow towards an end, as for 1/x), or the interval has no point of T strictly inside.
  - `error::non_finite`: the parenthesis "an infinite interval is for a later method, TanhSinh" becomes Phase 4's
    wording. The comment gains TanhSinh's tail examples:
    - 1/sqrt x on [0, 1]: ok;
    - 1/x: `not_converged`;
    - sqrt(x - 0.5): `non_finite`.
  - `GaussKronrod`: "substitute ... until TanhSinh" becomes "substitute ..., or use TanhSinh".

## Behaviour

### A. What TanhSinh receives

Phase 3a's `integrate` validates the criteria and the ends before TanhSinh is called. It returns {0, 0} for a == b
without calling fn, and puts the ends in order, negating the result for b < a. TanhSinh therefore receives
a < b, both finite. This replaces `PHASE3_API.md` A's first bullet; its points (the centre m = a/2 + b/2 and the
complement form, each side's points monotone in t, fn never at an end) stand.

### B. The steps

| | `start`: level 0 | `advance`: level k = 1 to 9 |
|---|---|---|
| fn evaluated at | the centre, then each side at t = 1, 2, ... until the side is cut (at most 2 floor(t_max) + 1: 13 in double) | the odd multiples of 2^-k inside each side's kept extent, and at most one node beyond it per side (C) |
| `current` | none: level 0 has no estimate | {V_k, E_k}, or none where E_k is unbounded |
| `attainable` | 2 F_0 | 2 F_k |
| `ready` | no | k >= 4 |
| `exhausted` | only where no point is strictly inside: then fn is never called | k = 9 (10 levels) |
| `trusted` | yes | yes |

V_k, F_k and E_k are `PHASE3_API.md` C and E's: the level's value, its floor (with the restored slivers) and its
estimate. The loop (`PHASE3A_API.md` B) then reproduces `PHASE3_API.md` E's termination:
- **From level 4 on,** it tests `E_k <= max(abs_tol, rel_tol |V_k|, 2 F_k)`. On a pass it settles (nothing to do)
  and returns {V_k, E_k}.
- **At level 9,** with the test failed, the result is `not_converged`. Its best is {V_9, E_9}, absent where E_9 is
  unbounded.
- **With no point strictly inside,** `start`'s step is exhausted and has no current, so the result is
  `not_converged` without a best.

The test is the same one, at the same levels, so the results are those `PHASE3_API.md` specified.

### C. The state

`detail::TanhSinhState<T>` holds no table and no reference to fn; it is a value of fixed size, with no allocation:
- **The level:** its number k, h = 2^-k, V_k and L1_k (Phase 1's Neumaier `Sum`), V_(k-1), and the slivers' share
  of the floor.
- **Per side:**
  - where the side is cut: the index of its last kept node in level k's grid;
  - how it was cut: by pruning, by a non-finite value or the table's end (D), or by rounding;
  - for the geometric end term: the terms of its outermost two kept nodes, without h;
  - for the power-law end term: the last kept point, and the nearest kept point at a different position, each as its
    distance to the end and its value;
  - W: the weights beyond the last kept node, for the sliver restored after a rounding cut.
- **The `Step`.**

**Each step needs only its own level's nodes,** for three reasons:
- **The new nodes are level k's.**
- **At most one new node lies beyond a side's last kept node.** At level k - 1 the side was cut at the node after
  its last kept node, t_last + 2^-(k-1). The only node of level k between the two is t_last + 2^-k. Nothing beyond
  the cut is evaluated:
  - after a cut by rounding, the points beyond also round onto the end, since points are monotone in t;
  - after a cut by pruning or a tail, `PHASE3_API.md` C and D cut the side at the first such node.
- **W is kept incrementally:** W_k = W_(k-1) + the weights of level k's nodes beyond the last kept node.
  - The grid of the earlier levels has no node strictly between t_last(k-1) and t_last(k), which is either
    t_last(k-1) or t_last(k-1) + 2^-k. So W_k sums exactly the weights `PHASE3_API.md` E names: every level's nodes
    beyond the last kept node, up to the current level.
  - Only the order of summation differs from a suffix sum. The difference is a rounding of a term that already
    enters the floor.
  - A side's cut can change from rounding to pruning or to a tail at a later level, never back. So W is needed only
    while the cut stays a rounding cut.

**What G asked, copy or regenerate, therefore does not arise.** Each level's pieces are read once, by the step that
reaches that level, through the same pattern as Phases 1 and 2:
- for float, double and long double, from the compile-time table;
- for another T at run time, from the run-time cache;
- for another T in constant evaluation, generated as a temporary.

**Constant evaluation with a T that has no table** (a user's literal type) generates each level in the evaluation
that calls `integrate`. At about 275 constexpr steps per node (Clang; `PHASE3_API.md` B):
- a run that stops at level 4 or 5 generates some 200 nodes of a double-like table, about 60,000 steps;
- a run to the cap generates all 3,361, about 920,000 steps, near Clang's default limit of 1,048,576.

This is documented, not engineered (Q4). The tabled types have no such cost: their tables are constants of their
own.

### D. Two points `PHASE3_API.md` leaves implicit

1. **A side that reaches the table's end uncut gets the geometric end term**, as a side cut by a non-finite value
   does (Q5).
   - **When it happens:** fn is finite and not negligible at the last node. x^-0.99 on [0, 1] in double is the
     example: fn is finite down to the smallest normal number, min, and the mass below it is 100 min^0.01 ≈ 0.084.
   - **Why it is implicit:** `PHASE3_API.md` E names three ways a side is cut: pruning, a non-finite value and
     rounding.
   - **The evidence:** its own figures for x^-0.99 at rel_tol 1e-4 fit only this reading. With Bailey's h Te at the
     cut, they report an error of 0.077, close to the mass beyond the table; with the geometric term, the result
     is `not_converged` with a best.
   - **To confirm:** the rebuilt prototype (Q3).
2. **Level 0's pruned terms leave V_0, L1_0 and F_0.** The outer terms of each side below eps L1_0 are dropped, and
   the level's sums are those of the kept terms. Level 1 halves V_0 and adds its new terms. The rebuilt prototype
   does the same, so that the count tests compare the same computation.

### E. Observing

`integrate<qd::TanhSinh>(qd::observer(o), fn, a, b[, criteria])` reports level 0 and then each level, as B's table
gives them.
- **Level 0** has no current. Phase 3a's printing example shows it as a dash.
- **The evaluations per level** are counted by the caller, as in that example.
- **Numbers specific to the method** (the level, the cut of each side) stay deferred (`PHASE3A_API.md` Q14).

No sample output here: it comes from the implementation. Step 3 runs Phase 3a's printing example on TanhSinh and adds
its output to this document.

## The work, in three steps (Q2)

Each step is committed and pushed, and checked on the five toolchains; the user approves the next.
1. **`detail::exp`** (internal; `PHASE3_API.md` B and Q10).
   - Tests against mpmath within 1.5 ulp in float, double, long double and both 50-digit types.
   - 0 at and below x_min.
   - `static_assert`s at both ends of the domain.
2. **The tables in pieces** (internal; `PHASE3_API.md` B and Q3).
   - Each node against mpmath within the a-priori bound, and the node counts per level.
   - The compile-time and run-time tables bit-identical.
   - The largest piece's share of each toolchain's constexpr budget recorded in `design.md`.
3. **`TanhSinh`.**
   - The public struct and the comments of the API section.
   - The tests below, against closed forms, mpmath and the prototype's counts.
   - `design.md`, the reference page and the demo.
   - The output of the printing example (E).

Then the phase's independent review.

Following Phase 3a, each side's outward sweep and each level's sum are written with standard algorithms where they
fit ("No Raw Loops"). One example is `std::ranges::all_of` over `std::views::iota` for a sweep that stops at a cut.

## Tests (`tests/test_quadrature_tanh_sinh.cpp`, test cases `quadrature:tanh_sinh:*`)

`PHASE3_API.md`'s list stands, with the protocol's additions:
- **The concept:** `static_assert(qd::iterative_method<qd::TanhSinh, T>)` for float, double, long double and both
  50-digit types. In `tests/test_quadrature_protocol.cpp`, TanhSinh joins the `observable` and the "named and as a
  value" checks.
- **The steps**, observed in double:
  - level 0 has no current and is not ready;
  - `ready` from level 4;
  - a run to the cap reports 10 steps, and only the last is exhausted;
  - the evaluations counted per level match B's table;
  - the last report equals the result, or the best;
  - with no point strictly inside: one report, exhausted, and fn never called.
- **Observed equals unobserved,** in both forms, in constant evaluation too.
- **Reversed ends** give exactly the negated result (Phase 3a's Q4, for every method).
- **The table's end** (D1): x^-0.99 on [0, 1] in double at rel_tol 1e-4 is `not_converged` with a best that covers
  the error. This is already in `PHASE3_API.md`'s termination tests.
- **Counts** exactly in double, against the rebuilt prototype with the C++ tables loaded (Q3).

## Documentation

- `PHASE3_API.md`'s list: `TanhSinh`'s comment; `design.md` with the sources and the numbers the decisions cite;
  the reference page (a `doxygenstruct` for `TanhSinh`); the demo.
- `docs/plan.md`: Phase 3b's checkboxes, one per step.
- This document: the decisions taken, the prototype's confirmation of D, and the example's output.

## Decisions for the user

- **Q1. The state holds no table** (C): each step reads its own level's pieces once, and W is kept incrementally.
  Recommended: yes.
  - It answers G's question with neither of its two options.
  - Results are the same as with a suffix sum, up to the rounding of a term inside the floor.
  - **Alternative:** copy the pieces into the state, as Gauss-Kronrod copies its rule. That is up to 3,361 node
    pairs in double, copied per call.
  - **Alternative:** regenerate the earlier levels at each step: exp per node again.
- **Q2. Three steps** (exp, the tables, the method), each checked on CI before the next. Recommended: yes.
  - The first two are internal and can be reviewed on their own.
  - **Alternative:** one step.
- **Q3. The prototype rebuilt, and kept.** It is rebuilt in Python with mpmath from `PHASE3_API.md` A to E and this
  document, and kept in `docs/quadrature/tanh-sinh-prototype/`, as Phase 3a's spike is kept. Recommended: yes.
  - **Why:** Phase 3's prototype lived in a session scratchpad and is gone. The count tests need it, and D needs
    its confirmation.
  - **Its tables:** loaded from the C++ ones, as Phases 2 and 3 did.
  - **Alternative:** rebuilt in a scratchpad again.
  - **Alternative:** no count tests, closed forms only.
- **Q4. The constant-evaluation cost for a T without a table is documented, not engineered** (C). Recommended: yes.
  - **Alternative:** in constant evaluation, generate only the nodes a step evaluates. That is cheaper there, but a
    second path through the code.
- **Q5. A side that reaches the table's end gets the geometric end term** (D1), subject to the prototype's
  confirmation. Recommended: yes.
  - **Alternative:** an end term of 0, as after pruning. The mass beyond the table (0.084 for x^-0.99 in double)
    would then enter no estimate.

## Deferred

As `PHASE3_API.md`:
- infinite intervals (Phase 4, with a wider exp);
- fn of the distance to the end;
- exp as a public helper;
- a stagnation rule.

As `PHASE3A_API.md`: numbers specific to the method in the reports.

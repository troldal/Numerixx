# Phase 3a API proposal: the method protocol

Status: draft 3, for review (numerics-reviewer, simplicity-critic, usage-sketcher), then the user. History:
- **Draft 1:** one member per method.
- **Draft 2:** after the first review round ("Decisions from the reviews", below).
- **Draft 3:** the user's request (2026-11-02) that the methods only compute while `integrate` decides when to
  stop (single responsibility), and that the numbers of each iteration can be reported (debugging).

Draft 3 adds the iterative protocol, `Step`, `Progress` and `observed()` (sections B, C and E, and Q11 to Q14), and
keeps every decision of draft 2. The earlier decisions it builds on:
- `docs/open-closed/ANALYSIS.md` (2026-11-02): the methods open to users (D1), quadrature first (D4).
- `docs/quadrature/BRIEF.md` (approved 2026-10-30), amended in two places: D1 ("methods as a closed set of empty
  tags") and 3.9 ("No step interface").
- Phase 3 (`PHASE3_API.md`, TanhSinh, approved 2026-11-01, not yet implemented) becomes Phase 3b. It adds `TanhSinh`
  through this protocol (G).

The evidence is a spike on a patched copy of `quadrature.hpp` (Spike, below).

## Summary

```cpp
#include <numerixx/quadrature.hpp>
namespace qd = nxx::quadrature;

// Unchanged: every call that compiles today compiles, with the same results except on reversed ends (Q4).
auto const r = qd::integrate(fn, 0.0, 1.0);
auto const s = qd::integrate<qd::Romberg>(fn, 0.0, 1.0, {.rel_tol = 1e-12});
auto const g = qd::integrate<qd::GaussLegendre<20>>(fn, 0.0, 1.0);

// New: each iteration reported, for any iterative method, the library's or a user's.
auto const p = qd::integrate(qd::observed(qd::GaussKronrod<15>{}, [](qd::Progress<double> const &p) {
                                 std::println("{:3} {:5} {:.15f} ± {:.3g}", p.iteration, p.evaluations,
                                              p.step.current->value, p.step.current->abs_error);
                             }),
                             fn, 0.0, 1.0, {.rel_tol = 1e-10});

// New: a method of one's own. An iterative method computes and reports; integrate() decides when to stop.
struct Trapezoid
{
    template<nxx::real T> struct State { /* ... */ qd::Step<T> step; };
    template<nxx::real T, typename F>
    constexpr std::expected<State<T>, qd::Failure<T>> start(F &fn, T const &a, T const &b) const;  // a < b
    template<nxx::real T, typename F>
    constexpr std::expected<void, qd::Failure<T>> advance(F &fn, State<T> &s) const;           // one more level
};
auto const t = qd::integrate<Trapezoid>(fn, 0.0, 1.0, {.rel_tol = 1e-10});

// A fixed rule, or a method that does not step (a recursion, a combinator), is one member, as in draft 2.
auto const n = qd::integrate(Midpoint{.panels = 1000}, fn, 0.0, 1.0);
auto const u = qd::integrate(SquareAtLowerEnd<qd::GaussKronrod<15>>{}, fn, 0.0, 1.0, {.rel_tol = 1e-12});
```

- **Three public concepts:**
  - `fixed_rule`: `apply(fn, a, b)`.
  - `iterative_method`: `start(fn, a, b)` and `advance(fn, state)`, optionally `settle(state)`.
  - `estimating_method`: `estimate(fn, a, b, criteria)`, the whole computation in one member.
- **An iterative method computes; `integrate` decides.** The method reports a `Step` after its start and after each
  step: its current result, the least error it can certify, and whether it is ready, exhausted and trusted.
  `integrate` owns the loop and the test against the criteria. The method never sees the criteria.
- **`Romberg`** steps by level and **`GaussKronrod<15>`** by bisection. Their arithmetic is unchanged, and so are
  their results: the battery's 2,438 results are identical to draft 2's.
- **`observed(m, observer)`** reports each iteration, `Progress{iteration, evaluations, step}`, of any iterative
  method. It is an adaptor, not a parameter of `integrate`, and without it nothing is counted or called.
- **From draft 2, unchanged:**
  - `integrate` validates, puts the ends in order (b < a negates), and checks the result's finiteness.
  - It takes a method named, `integrate<M>(fn, ...)` as today, or as a value, first: `integrate(m, fn, ...)`.
  - The closed lists and the `if constexpr` dispatch are removed.
  - No ADL. Everything constexpr. Nothing throws.

## Public API (changes to `quadrature.hpp`)

```cpp
namespace nxx::quadrature {

    // ... error, Integral, Failure, Criteria as today; then:

    /// Where an iterative method stands after its start or a step.
    template<real T>
    struct Step
    {
        std::optional<Integral<T>> current; ///< the value and estimate now; empty where it has no usable estimate
        T attainable;                       ///< the least error it can certify: Romberg F, Gauss-Kronrod 2F
        bool ready;                         ///< it has done its minimum work (Romberg: level 4)
        bool exhausted;                     ///< it cannot go on: its cap, stagnation, every panel final
        bool trusted;                       ///< a result that meets the request may be returned as a success
    };

    // ... detail; then the concepts:

    template<typename M, typename T, typename F = T(T const &)>
    concept fixed_rule = requires(M const &m, F &fn, T const &a, T const &b) {
        { m.apply(fn, a, b) } -> std::same_as<std::expected<T, error>>;
    };

    template<typename M, typename T, typename F = T(T const &)>
    concept estimating_method = requires(M const &m, F &fn, T const &a, T const &b, Criteria<T> const &criteria) {
        { m.estimate(fn, a, b, criteria) } -> std::same_as<std::expected<Integral<T>, Failure<T>>>;
    };

    /// start returns std::expected<S, Failure<T>> for a state S of the method's own, whose member step is current.
    template<typename M, typename T, typename F = T(T const &)>
    concept iterative_method = requires(M const &m, F &fn, T const &a, T const &b) {
        { m.start(fn, a, b) } -> /* std::expected<S, Failure<T>> */;
    } && requires(M const &m, F &fn, /* S */ &s) {
        { m.advance(fn, s) } -> std::same_as<std::expected<void, Failure<T>>>;
        { std::as_const(s).step } -> std::same_as<Step<T> const &>;
    };

    /// What an observer sees after the start (iteration 0) and after each step.
    template<real T>
    struct Progress
    {
        int iteration;
        std::size_t evaluations; ///< of fn, so far
        Step<T> step;
    };

    /// An iterative method whose progress is reported; itself an iterative method.
    template<typename M, typename Observer>
    struct Observed { M method; Observer observer; /* start, advance, settle: M's, counted and reported */ };

    template<typename M, typename Observer>
    constexpr Observed<M, std::decay_t<Observer>> observed(M const &m, Observer &&observer);

    template<int N> struct GaussLegendre { /* apply, as in draft 2 */ };
    struct Romberg { /* start: T_0 and level 1; advance: the next level */ };
    template<int N> requires(N == 15) struct GaussKronrod { /* start: the mesh; advance: a bisection; settle */ };

    // The four forms of integrate, as in draft 2; the estimating two take an iterative or an estimating method.
    template<typename M, real T, typename F>
    requires nxx::detail::callable<F, T> && fixed_rule<M, T, std::remove_reference_t<F>>
    [[nodiscard]] constexpr std::expected<T, error> integrate(M const &m, F &&fn, T const &a, T const &b);

    template<typename M, real T, typename F>
    requires nxx::detail::callable<F, T>
             && (iterative_method<M, T, std::remove_reference_t<F>> || estimating_method<M, T, std::remove_reference_t<F>>)
    [[nodiscard]] constexpr std::expected<Integral<T>, Failure<T>>
        integrate(M const &m, F &&fn, T const &a, T const &b, Criteria<std::type_identity_t<T>> const &criteria = {});

    template<typename M, real T, typename F> /* + std::default_initializable<M> */
    [[nodiscard]] constexpr std::expected<T, error> integrate(F &&fn, T const &a, T const &b);

    template<typename M = GaussKronrod<15>, real T, typename F> /* + std::default_initializable<M> */
    [[nodiscard]] constexpr std::expected<Integral<T>, Failure<T>>
        integrate(F &&fn, T const &a, T const &b, Criteria<std::type_identity_t<T>> const &criteria = {});
}
```

- **The state is the method's own type,** opaque to `integrate` except for its member `step`. It is returned by value
  from `start`, so it is movable, and it lives for one call of `integrate`. The library's states are in `detail`
  (`RombergState<T>`, `KronrodState<T>`).
- **The methods** move after `detail` and define their members inline, as in draft 2. Gauss-Kronrod's state copies
  its rule (24 values of T), because for a type without a table, in constant evaluation, the rule is a temporary.
- **A type with both an iterative and an estimating form** is driven: the iterative form wins, with no ambiguity. A
  type that is both a fixed rule and an estimating method is ambiguous when named, as in draft 2.
- **The concepts** take `<M, T, F>`, F defaulted, as in draft 2, so that `static_assert(qd::iterative_method<MyMethod,
  double>)` beside a method is the usual check.

## Behaviour

### A. What `integrate` does

In order, as in draft 2:
1. The tolerances, for a method with an estimate (negative or not finite: `non_finite`); `rel_tol` clamped to 1.
2. The ends (NaN or infinite: `non_finite`). a == b: 0, or `{0, 0}`, without calling fn.
3. The ends in order (Q4). The method receives `(min, max)`, and for b < a `integrate` negates the value, and a
   best's value.
4. The method:
   - a fixed rule: `apply`;
   - an iterative method: the loop of B;
   - an estimating method: `estimate`, with the validated criteria (Q6).
5. The result checked in every build (Q5): a value must be finite; an `Integral`, or a best, must have a finite value
   and a finite, non-negative `abs_error`; a best may come only with `not_converged`. Anything else is returned as
   `non_finite` without a best.

### B. The loop (`integrate` on an iterative method)

```
s = start(fn, a, b)                                   (a Failure ends the run)
loop:
    passed = s.step.ready and s.step.current and abs_error <= max(abs_tol, rel_tol |value|, s.step.attainable)
    if passed or s.step.exhausted:
        settle(s), if the method has it              (a Failure ends the run)
        if the settled step passes: success if trusted, else not_converged with it as the best
        if it did not pass before settling: not_converged, with s.step.current as the best (or none)
        otherwise go on from the settled step        (the running sums had drifted: Gauss-Kronrod)
    else:
        advance(fn, s)                                (a Failure ends the run)
```

- **The test against the request** is `BRIEF.md` D5's, the same for every method: `abs_error <= max(abs_tol,
  rel_tol |value|, attainable)`. It is equivalent to each method's test today. Romberg tests its estimate E against
  `max(abs_tol, rel_tol |value|, F)`, and Gauss-Kronrod against `max(..., 2F)`. The step's `abs_error` is max(E, F),
  and F <= attainable <= the right-hand side, so `max(E, F) <= W` exactly when `E <= W`. The comparisons are exact,
  so the decisions are the same bit for bit.
- **Termination is the method's:** the loop advances only while the step is not exhausted. So the contract (D)
  requires every iterative method to report `exhausted` within finitely many steps. The library's do: a cap, plus
  for Gauss-Kronrod stagnation and an empty heap.
- **`settle`** is Gauss-Kronrod's check against sums formed afresh, which Phase 2 makes before every decision to stop
  (independent review, 2098416). A method without it is decided on its step as it stands.

### C. What stays in the method, and why

The loop owns the user's request; the method owns its own limits, which no driver can know:
- **`ready`, the minimum work:** Romberg's level 4 (fewer samples of sin^2 on [0, 2 pi] all vanish); TanhSinh's
  level 4.
- **`attainable`, the least error it can certify:** Romberg's floor F; Gauss-Kronrod's and TanhSinh's 2F (a kink or
  a jump never brings every panel to its floor).
- **`exhausted`:** the caps (Romberg 16 levels; Gauss-Kronrod 64 bisections per bit of T; TanhSinh 10 levels),
  Gauss-Kronrod's stagnation (32 bisections without a new minimum of its excess), and its empty heap.
- **`trusted`:** Gauss-Kronrod's undone descent, a piece where fn is not finite that no short bisection avoided. Such
  a run cannot end as a success.
- **`current` empty:** TanhSinh's unbounded end term, or an interval with no point strictly inside. The run then
  ends `not_converged` without a best, as `PHASE3_API.md` Q7 decided.

One common rule for these would change when the methods stop. A common minimum or floor, for example, would move
Romberg's or Gauss-Kronrod's results.

### D. The contract

Written once, in the concepts' Doxygen comments, as in draft 2.

- **Given to every method:**
  - a < b, both finite, possibly extreme: b - a may overflow, the ends may be subnormal or adjacent, or of opposite
    signs so that b - a rounds.
  - fn as an lvalue: call it with `std::invoke` and a T, and convert the result with `T(...)`.
  - An estimating method also gets the criteria: finite, non-negative, `rel_tol <= 1`, `abs_tol` unclamped, possibly
    both 0. An iterative method never sees them.
- **Required of every method:**
  - fn is evaluated only at points of [a, b]; clamp a point formed as a + h t.
  - A finite result (checked, A.5), or an error: `non_finite` without a best, or `not_converged` with or without one.
  - Its parameters are its own preconditions (Q9).
  - fn's exceptions, and the method's own, pass through.
- **Required of an iterative method:**
  - `start` and every `advance` leave `step` describing the state.
  - It reports `exhausted` within finitely many steps.
  - `advance` is never called on an exhausted state, and nothing is called after a Failure.
  - `current`, where present, satisfies A.5.
  - `settle` may only make `step` more exact; it must not take a step.
- **Required of an estimating method:** it terminates for every criteria it can be given, both tolerances 0
  included.
- **Required in meaning, not checkable:**
  - `abs_error` estimates |value - exact| and covers the rounding of the method's own sum.
  - `attainable` is an error the method can actually certify.
  - For an estimating method, success means its own stop test was met.
- **The library's methods only** (Q8): the rounding floor (for a != b), `abs_tol = 0` working, a best with every
  `not_converged` of Romberg and Gauss-Kronrod, fn never evaluated at an end (GaussLegendre, TanhSinh), and
  `non_finite` in constant evaluation where an overflow or a NaN would otherwise be a compile error.
- **Infinite ends** stay outside this contract (Q10); Phase 4 makes them opt-in per method.

### E. Progress (`observed`)

- `observed(m, observer)` wraps an iterative method. Its `start` and `advance` call m's with fn wrapped in a counter,
  then call `observer(Progress<T>{iteration, evaluations, step})`: iteration 0 after the start, then one per step.
  Its `settle` forwards to m's and reports nothing; the result `integrate` returns shows the settled values.
- It is itself an iterative method, so it works with every form of `integrate`, and an observed user method needs
  nothing more from its author.
- **Cost:** without `observed`, nothing is counted or called. With it, a counter per evaluation and one call per
  step. In constant evaluation it works with an observer that is itself constexpr (a lambda that records nothing, or
  one that writes to a constexpr-usable container).
- **The evaluations** come from the counter, not from the method, so every method's count is right without its
  author's help. This also covers `BRIEF.md` D12's deferred "evaluations" for those who want them, without a new
  member of `Integral`.
- **Not in the protocol:** the observer cannot stop the run (Q14), and method-specific numbers (Gauss-Kronrod's panel
  count, Romberg's tableau) are not reported; `Step` holds what every iterative method has.
- An estimating method (one member) and a fixed rule have no steps to report. `observed` requires an iterative
  method.

### F. The library's methods

- **`GaussLegendre<N>::apply`:** unchanged from draft 2.
- **`Romberg`:** `start` forms T_0 from the ends and level 1, and `advance` forms the next level. The level is
  Phase 1's loop body, line for line, with the stop test removed. Its step:
  - current: {R(k, k), max(estimate, F)};
  - attainable F;
  - ready from level 4;
  - exhausted at level 16;
  - always trusted.
- **`GaussKronrod<15>`:** `start` builds the initial mesh, `advance` bisects the panel with the largest estimate, and
  `settle` forms the sums afresh. Phase 2's loop is kept line for line, with its locals as members of
  `detail::KronrodState<T>` and the stop test moved to the loop. Its step:
  - current: {K, max(E, F)}, negated for reversed ends, which `integrate` no longer passes;
  - attainable 2F;
  - always ready;
  - exhausted by an empty heap, the cap or stagnation;
  - trusted unless a descent was undone.

### G. Phase 3b (TanhSinh) through the protocol

`PHASE3_API.md` stays approved, except:
- **Its "Public API" bullet:** `TanhSinh` is an iterative method (`start`: level 0 with its outward sweep; `advance`:
  the next level), so the existing estimating forms of `integrate` take it.
- **Its A:** "b < a negates" is done by `integrate`; TanhSinh receives a < b. The centre not strictly inside gives a
  start whose step is exhausted with no current: `not_converged` without a best.
- **Its E:** the stop test `E_k <= max(abs_tol, rel_tol |V_k|, 2 F_k)` from level 4 moves to the loop as ready
  (level >= 4) and attainable 2F_k. The cap of 10 levels is `exhausted`. An unbounded estimate is a step without
  current, so the best is absent (its Q7).
- **Progress:** an observed TanhSinh reports each level.

Its behaviour otherwise, the decisions Q1 to Q10 and the tests are unchanged. `docs/plan.md` splits Phase 3 into 3a
and 3b.

### H. What a method of one's own cannot use yet

As in draft 2: `Interval`, `Sum`, `rounding_floor`, the overflow checks for constant evaluation and the node tables
stay in `detail`; the guide page gives the floor's recipe (50 eps L1 from the method's own samples), which an
iterative method reports as its attainable error.

## Spike

`docs/open-closed/quadrature-spike/` (README there):
- `quadrature.patch.py ... normalize` writes draft 2 into a copy of `quadrature.hpp`.
- `steps.py` then writes draft 3: `Step`, the stepped Romberg and Gauss-Kronrod, the loop, `iterative_method`,
  `Progress` and `observed`.

Built with GCC 13.3 and Boost 1.83, and the programs also with Clang 18, which needs `-U__cpp_concepts
-D__cpp_concepts=202002L` for libstdc++ 13's `<expected>`. The container's network refuses Boost 1.92, and the
project's five toolchains are not available here; the phase itself builds and tests on all five.

- **The existing quadrature tests** (`test_quadrature_gauss.cpp`, `_romberg.cpp`, `_kronrod.cpp`), built unchanged
  against draft 3, pass exactly as against today's header and drafts 1 and 2. That is 25 of 26 test cases. The one
  failure is the same in every build: the `GaussLegendre<2>` long double node baseline, which was recorded on the
  project's toolchains and does not hold for GCC 13's 80-bit long double. The tests' `static_assert`s, Gauss-Kronrod
  and Romberg in constant evaluation among them, now run through the loop of B.
- **Bit for bit** (`dump.cpp`): the battery covers 8 integrands on 6 intervals (three reversed); 4 Gauss-Legendre
  rules; Gauss-Kronrod at 3 tolerances; Romberg at 3; float, double, long double, `cpp_bin_float_50` and
  `cpp_dec_float_50`; edges and three constant-expression results: 2,438 results.
  - Draft 1, and draft 2 without ordering the ends: identical to today.
  - Draft 2 with the ends ordered (Q4): 23 differ, all on reversed ends, with no status among them. Romberg's
    `abs_error` alone in 9 float results and 1 double result (the floor's L1 sum); in `cpp_dec_float_50`, 7
    Gauss-Legendre values and 6 Romberg results, in the guard digits (about the 70th significant digit).
    Gauss-Kronrod: none.
  - Draft 3: identical to draft 2 with the ends ordered, all 2,438 results, so the stepped methods and the loop of B
    change nothing.
- **Orientation today** (`orient.cpp`): 3,000 intervals (positive, negative, across 0, up to 2e6 wide), each
  integrated forwards and reversed, the whole result compared. Gauss-Legendre<64> and Gauss-Kronrod<15> give the
  exact negation every time. Romberg gives the exact negation in value and status every time, and a different
  `abs_error` in 374 of 6,000 runs.
- **Progress** (`progress.cpp`):
  - **Same results observed:** `observed()` gives the same result as the unobserved method, in `static_assert`s
    (Gauss-Kronrod on x^3 - x over [-1, 2]; Romberg reversed at 1e-10) and at run time.
  - **Gauss-Kronrod on a peak** of width 1e-3 at rel_tol 1e-10: 15 reports (the mesh, then 14 bisections), with
    evaluations 150 to 570. The 570 equals both a count kept by fn and the unobserved run's.
  - **Romberg on e^x at 1e-12:** a report per level, ready from level 4, the result at level 5 (33 evaluations,
    1.7182818284590453).
  - **A user's iterative method** (a trapezoidal rule halving its step, about 40 lines): it gets the stop test and
    the reports without writing either. `static_assert(qd::iterative_method<Trapezoid, double>)`; x^3 - x over
    [-1, 2] in a constant expression.
- **Methods of one's own** (`user_methods.cpp`, unchanged from draft 2, passing against draft 3): the fixed rule
  `Midpoint`, the combinator `SquareAtLowerEnd<Inner>` over the iterative library methods, the spy, the checked
  results, and the concepts' rejections. On 1/sqrt x and cos x / sqrt x over [0, 1] at 1e-12, the substitution takes
  150 evaluations against Gauss-Kronrod's 2,370 and 2,400.

## Tests (`tests/test_quadrature_protocol.cpp`, test cases `quadrature:protocol:*`)

As in draft 2:
- the existing suite unchanged;
- both forms agreeing;
- reversal exact for every method;
- the spy;
- the midpoint rule's closed forms and its domain check;
- the combinator's derived tolerances, in double and `cpp_bin_float_50` with and without expression templates;
- the checked result;
- the concepts.

New:
- **The loop's branches** with a scripted iterative method, whose steps come from a list:
  - not ready: it goes on;
  - passing but untrusted: `not_converged` with the best;
  - exhausted and failing: `not_converged` with the best, and without one where current is empty;
  - passing before settle and failing after: it goes on from the settled step;
  - a Failure from `start`, `advance` or `settle`: returned as it is;
  - `advance` never called on an exhausted step.
- **`observed` changes nothing:** `integrate(observed(M{}, quiet), ...) == integrate<M>(...)` for Romberg and
  Gauss-Kronrod, as `static_assert`s in double with run-time twins; reversed too.
- **The reports:** iterations 0, 1, 2, ... in order; evaluations non-decreasing and, at the end, equal to a count kept
  by fn; the last report's step equal to the result for Romberg; for Gauss-Kronrod, equal unless settle moved it.
- **A user's iterative method** (the trapezoidal rule) against closed forms: x^3 - x over [-1, 2] = 9/4, its error
  (b - a) h^2 f'' / 12 in closed form at each level; e^x over [0, 1] within its estimate and the tolerance.

## Documentation

As in draft 2, and:
- `quadrature.hpp`: the comments of `Step`, `Progress`, `iterative_method` and `observed`. The contract of D goes in
  the concepts. The file comment says that a method computes and `integrate` decides, and shows `observed` in its
  example.
- `docs/quadrature/design.md`: "The interface" records the loop and what stays in the method (C); BRIEF.md 3.9
  marked as amended.
- `docs/sphinx/quadrature/reference.md`: `Step`, `Progress`, the three concepts and `observed`.
- The guide page on writing a method shows the three forms, with a user's iterative method as its main example, and
  `observed` for debugging.

## Decisions for the user

From draft 2 (unchanged unless noted):
- **Q1. Member names:**
  - `apply` (fixed rule) and `estimate` (estimating method);
  - draft 3 adds `start`, `advance` and `settle` (iterative method).

  Recommended: yes.
- **Q2. Concepts:**
  - `fixed_rule`, `estimating_method` and, new, `iterative_method`;
  - each `<M, T, F = T(T const &)>`, no copy requirement.

  Recommended: yes.
- **Q3. A method value goes first,** `integrate(m, fn, a, b[, criteria])`, beside today's named forms, which keep
  their signatures. Recommended: yes.
- **Q4. `integrate` puts the ends in order and negates.** Recommended: yes.
  - **Cost:** 23 of the battery's 2,438 results change, all on reversed ends, no status. That covers Romberg's
    `abs_error` in float and double, and values beyond `cpp_dec_float_50`'s 50 digits.
  - **Alternative:** the method orients itself, bit for bit with today.
- **Q5. The result's finiteness checked in every build,** a broken result returned as `non_finite`. Recommended: yes.
- **Q6. The validated criteria passed as `Criteria<T>`** to an estimating method. Recommended: yes.
- **Q7. The protocol as Phase 3a, before TanhSinh (3b).** Recommended: yes.
- **Q8. The statements that hold only for the library's methods say so.** Recommended: yes.
- **Q9. A method's parameters are its own preconditions;** no new enumerator now. Recommended: yes.
- **Q10. Finite ends only;** infinite ends opt-in per method in Phase 4. Recommended: yes.

New in draft 3:
- **Q11. `integrate` owns the loop and the test against the request; the method keeps its own limits** (`ready`,
  `attainable`, `exhausted`, `trusted`, `settle`; C). Recommended: yes. The results are unchanged, and the test is
  `BRIEF.md` D5's.
  - **Alternative, a common rule for the limits:** one floor multiple and one minimum. That changes when Romberg and
    Gauss-Kronrod stop.
  - **Alternative, draft 2:** each method owns its loop, with no reports.
- **Q12. `Step` as above:** `current` optional, `attainable`, `ready`, `exhausted`, `trusted`. Recommended: yes.
  - **`current` is optional** so that TanhSinh's "no best when the estimate is unbounded" (its Q7) needs no change
    later.
  - **The names are open to change.** `attainable` could be `floor`, but the floor is F while Gauss-Kronrod's
    attainable is 2F; `trusted` could be `acceptable`.
- **Q13. The one-member `estimate` form stays beside the iterative one,** for methods that do not step: a
  recursion, a combinator over a library method, or a port of an existing routine. Such methods are not observed.
  Recommended: yes.
  - **Alternative:** the iterative form only. That is one concept fewer, but a recursive method must become a work
    list, and a combinator must forward start and advance to its inner method.
- **Q14. Progress through `observed(m, observer)`,** an adaptor, reporting `Progress{iteration, evaluations, step}`
  after the start and after each step; the observer cannot stop the run. Recommended: yes.
  - **Alternative:** an observer parameter of `integrate` (a fifth argument, in every form).
  - **Deferred:** a stopping observer, or a public stepping interface (`start`, then `advance` in the caller's own
    loop), which the same protocol allows later without a change.

## Decisions from the reviews

Draft 2's round (on draft 1):

From the numerics review:
- **Adopted:**
  - fn evaluated only in [a, b] required of every method, and the spike's methods clamp and refuse an overflowing
    width. Draft 1's midpoint rule called fn at +inf on [-max, max], and its substitution at b + 2 on
    [-1, 2^53 + 2].
  - The orientation of methods tied to an end, answered by Q4.
  - Termination for every criteria, and parameters as the method's preconditions.
  - The edge inputs and constant evaluation in the contract.
  - Q8 extended, and the a == b inconsistency of the floor's wording fixed.
  - The tests' tolerances derived, a spy, an expression-template instantiation, a reversed combinator, the domain
    check.
  - The spike's presentation corrected.
  - `orient.cpp` compares whole results, which found Romberg's `abs_error`.
- **Confirmed:** the checked part for all four methods, the edges and constant evaluation; the members' arguments;
  the clamped criteria; the spike's values; Q4's derivation.

From the simplicity review:
- **Adopted:**
  - the methods moved after `detail` with inline members;
  - no `std::copy_constructible`;
  - the contract written once, in the header;
  - "a method is of one kind";
  - the tests trimmed;
  - "an lvalue";
  - Q2's reason restated.
- **Not adopted:**
  - Q5's "no check", in favour of the release check (a NaN was returned as a success under `NDEBUG`);
  - dropping the second integrand of the substitution.

From the usage sketch:
- **Adopted:**
  - the method value first (Q3);
  - the ends ordered (Q4);
  - the release check (Q5);
  - parameters (Q9);
  - infinite ends recorded (Q10);
  - `<M, T, F>` with F defaulted;
  - no copy requirement;
  - the concepts' comments on diagnostics;
  - the floor's recipe and the guide's style.
- **Not adopted:**
  - a precedence between the kinds (a sentence instead; in draft 3 an iterative form is preferred over an
    estimating one);
  - publishing `rounding_floor` now.

Draft 3's round: pending.

## Deferred

- The building blocks (`Interval`, `Sum`, `rounding_floor`, the overflow checks, the tables) as public API (H).
- A stopping observer, a public stepping interface, and method-specific progress (Q14, E).
- Infinite ends as an opt-in (Q10): Phase 4.
- `Adaptive<Pair>`, the globally adaptive loop over any embedded pair (`ANALYSIS.md` 6.3).
- The same protocol for `nxx::roots` and `nxx::deriv` (`ANALYSIS.md` 10). In roots the tolerance enters Brent's and
  TOMS 748's steps (`docs/roots/BRIEF.md` 3.8), so there the loop cannot own the test alone.

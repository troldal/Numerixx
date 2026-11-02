# Phase 3a API proposal: the method protocol

Status: approved by the user on 2026-11-02, Q1 to Q15 as recommended (draft 5). History:
- **Draft 1:** one member per method.
- **Draft 2:** after the first review round.
- **Draft 3:** at the user's request (2026-11-02), methods only compute, `integrate` decides when to stop (single
  responsibility), and each iteration's numbers can be reported (debugging).
- **Draft 4:** after draft 3's review round (numerics-reviewer, simplicity-critic, usage-sketcher; "Decisions from
  the reviews", below). The observer receives the `Step` itself; `settle` is required; the loop checks every step;
  `Step` has defaults.
- **Draft 5:** the user's syntax for observing (2026-11-02):
  - the forms `integrate<M>(qd::observer(o), fn, a, b[, criteria])` and `integrate(m, qd::observer(o), fn, a,
    b[, criteria])`;
  - the observer wrapped by `qd::observer`, not passed bare;
  - draft 4's `observed()` moved to `detail` as the mechanism.

  No new review round: the change is the user's choice of syntax, and the spike checks it (Q14).

The decisions it builds on:
- `docs/open-closed/ANALYSIS.md` (2026-11-02): the methods open to users (D1), quadrature first (D4).
- `docs/quadrature/BRIEF.md` (approved 2026-10-30), amended in two places: D1 ("methods as a closed set of empty
  tags") and 3.9 ("No step interface").
- Phase 3 (`PHASE3_API.md`, TanhSinh, approved 2026-11-01, not yet implemented) becomes Phase 3b and adds
  `TanhSinh` through this protocol (G).

The evidence is a spike on a patched copy of `quadrature.hpp` (Spike, below).

## Summary

```cpp
#include <numerixx/quadrature.hpp>
namespace qd = nxx::quadrature;

// Unchanged: every call that compiles today compiles, with the same results except on reversed ends (Q4).
auto const r = qd::integrate(fn, 0.0, 1.0);
auto const s = qd::integrate<qd::Romberg>(fn, 0.0, 1.0, {.rel_tol = 1e-12});
auto const g = qd::integrate<qd::GaussLegendre<20>>(fn, 0.0, 1.0);

// New: each step reported, for any iterative method, the library's or a user's.
auto const print = [](qd::Step<double> const &s) {
    if (s.current) std::println("{:.15f} ± {:.3g}", s.current->value, s.current->abs_error);
};
auto const p = qd::integrate<qd::Romberg>(qd::observer(print), fn, 0.0, 1.0, {.rel_tol = 1e-12});

// New: a method of one's own. An iterative method computes and reports; integrate() decides when to stop.
struct Trapezoid
{
    template<nxx::real T> struct State { /* ... */ qd::Step<T> step; };
    template<nxx::real T, typename F>
    constexpr std::expected<State<T>, qd::Failure<T>> start(F &fn, T const &a, T const &b) const;  // a < b
    template<nxx::real T, typename F>
    constexpr std::expected<void, qd::Failure<T>> advance(F &fn, State<T> &s) const;           // one more level
    template<nxx::real T>
    constexpr std::expected<void, qd::Failure<T>> settle(State<T> &) const { return {}; }      // nothing to refine
};
auto const t = qd::integrate<Trapezoid>(fn, 0.0, 1.0, {.rel_tol = 1e-10});
auto const v = qd::integrate(Trapezoid{}, qd::observer(print), fn, 0.0, 1.0);  // observed too, as a method value

// A fixed rule, or a method that does not step (a recursion, a combinator), is one member, as in draft 2.
auto const n = qd::integrate(Midpoint{.panels = 1000}, fn, 0.0, 1.0);
auto const u = qd::integrate(SquareAtLowerEnd<qd::GaussKronrod<15>>{}, fn, 0.0, 1.0, {.rel_tol = 1e-12});
```

- **Three public concepts:**
  - `fixed_rule`: `apply(fn, a, b)`.
  - `iterative_method`: `start(fn, a, b)`, `advance(fn, state)` and `settle(state)`.
  - `estimating_method`: `estimate(fn, a, b, criteria)`, the whole computation in one member.
- **An iterative method computes; `integrate` decides.** The method reports a `Step` after its start and after each
  step: its current result, the least error it can certify, and whether it is ready, exhausted and trusted.
  `integrate` owns the loop and the test against the criteria. The method never sees the criteria.
- **`Romberg`** steps by level and **`GaussKronrod<15>`** by bisection. Their arithmetic is unchanged, and so are
  their results: the battery's 2,438 results are identical to draft 2's.
- **`qd::observer(o)`** marks o as an observer: `o(step)` is called after the start and after each step of any
  iterative method, named (`integrate<M>(qd::observer(o), fn, ...)`) or given as a value (`integrate(m,
  qd::observer(o), fn, ...)`). It shares the method's state, so an observed run computes exactly what an unobserved
  one does.
- **From draft 2, unchanged:**
  - `integrate` validates, puts the ends in order (b < a negates), and checks the result's finiteness.
  - It takes a method named, `integrate<M>(fn, ...)` as today, or as a value, first: `integrate(m, fn, ...)`.
  - The closed lists and the `if constexpr` dispatch are removed.
  - No ADL. Everything constexpr. Nothing throws.

## Public API (changes to `quadrature.hpp`)

```cpp
namespace nxx::quadrature {

    // ... error, Integral, Failure, Criteria as today; then:

    /// Where an iterative method stands after its start or a step; the defaults are the neutral values.
    template<real T>
    struct Step
    {
        std::optional<Integral<T>> current{}; ///< the value and estimate now; empty where it has no usable estimate
        T attainable = T(0);   ///< the least error it can certify (Romberg F, Gauss-Kronrod 2F); finite, not negative
        bool ready = true;     ///< it has done its minimum work (Romberg: level 4)
        bool exhausted = false; ///< it cannot go on: its cap, stagnation, every panel final
        bool trusted = true;   ///< a result that meets the request may be returned as a success
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

    /// The state of an iterative method: what start returns, inside its std::expected.
    template<typename M, typename T, typename F = T(T const &)>
    using state_t = typename decltype(std::declval<M const &>().start(
        std::declval<F &>(), std::declval<T const &>(), std::declval<T const &>()))::value_type;

    template<typename M, typename T, typename F = T(T const &)>
    concept iterative_method = requires(M const &m, F &fn, T const &a, T const &b) {
        { m.start(fn, a, b) } -> /* std::expected<S, Failure<T>> */;
    } && requires(M const &m, F &fn, state_t<M, T, F> &s) {
        { m.advance(fn, s) } -> std::same_as<std::expected<void, Failure<T>>>;
        { m.settle(s) } -> std::same_as<std::expected<void, Failure<T>>>;
        { std::as_const(s).step } -> std::same_as<Step<T> const &>;
    };

    /// o, marked as an observer: o(Step<T> const &), returning void, is called after the start and after each step.
    template<typename O>
    struct Observer { /* the callable: by reference for an lvalue, by value for an rvalue */ };

    template<typename O>
    constexpr Observer<O> observer(O &&o);

    template<int N> struct GaussLegendre { /* apply, as in draft 2 */ };
    struct Romberg { /* start: T_0 and level 1; advance: the next level; settle: nothing */ };
    template<int N> requires(N == 15) struct GaussKronrod { /* start: the mesh; advance: a bisection; settle */ };

    // The four forms of integrate, as in draft 2; the estimating two take an iterative or an estimating method.
    // Then two forms with an observer, for an iterative method.
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

    template<typename M = GaussKronrod<15>, typename O, real T, typename F>
    requires std::default_initializable<M> && nxx::detail::callable<F, T>
             && iterative_method<M, T, std::remove_reference_t<F>>
    [[nodiscard]] constexpr std::expected<Integral<T>, Failure<T>> integrate(Observer<O> o, F &&fn, T const &a,
        T const &b, Criteria<std::type_identity_t<T>> const &criteria = {});

    template<typename M, typename O, real T, typename F>
    requires nxx::detail::callable<F, T> && iterative_method<M, T, std::remove_reference_t<F>>
    [[nodiscard]] constexpr std::expected<Integral<T>, Failure<T>> integrate(M const &m, Observer<O> o, F &&fn,
        T const &a, T const &b, Criteria<std::type_identity_t<T>> const &criteria = {});
}
```

- **The state is the method's own type,** opaque to `integrate` except for its member `step`. It is returned by value
  from `start`, so it is movable, and it lives for one call of `integrate`. `state_t` names it for a combinator
  that keeps an inner method's state. The library's states are in `detail` (`RombergState<T>`, `KronrodState<T>`).
- **`Step`'s defaults are the neutral values:** a method sets only the limits it has, with designated initialisers,
  `Step<T>{.current = ..., .attainable = F, .ready = k >= 4}`.
- **The observed forms forward** to the value form of `integrate` with `detail::observed(m, o)`, an iterative
  method that wraps m and shares its state (draft 4's `observed()`, now internal). They are told from the other
  forms by the type `Observer<O>`, never by calling an argument on trial (Q14).
- **The methods** move after `detail` and define their members inline, as in draft 2. Gauss-Kronrod's state copies
  its rule (24 values of T), because for a type without a table, in constant evaluation, the rule is a temporary.
- **A type of two kinds** is outside the contract: `integrate` promises nothing about which form it uses.
- **Calling a method's members directly** (`Romberg{}.start(fn, a, b)`, then `advance`) compiles: they are public,
  because `integrate` calls them. It is outside the contract, though, which assumes `integrate`'s validation and
  ordering of the ends and its use of `settle` before any decision. A public stop test for one's own loop is
  deferred.

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
    check s.step: current, where present, as A.5; attainable finite and not negative; else non_finite, no best
    passed = s.step.ready and s.step.current and abs_error <= max(abs_tol, rel_tol |value|, s.step.attainable)
    if passed or s.step.exhausted:
        settle(s), then check s.step again          (a Failure ends the run)
        if the settled step passes: success if trusted, else not_converged with it as the best
        if it did not pass before settling: not_converged, with s.step.current as the best (or none)
        otherwise go on from the settled step        (the running sums had drifted: Gauss-Kronrod)
    else:
        advance(fn, s)                                (a Failure ends the run)
```

- **The test against the request** is `BRIEF.md` D5's, the same for every method: `abs_error <= max(abs_tol,
  rel_tol |value|, attainable)`. It is equivalent to each method's test today:
  - Romberg tests its estimate E against `max(abs_tol, rel_tol |value|, F)`, and Gauss-Kronrod against
    `max(..., 2F)`. The step's `abs_error` is max(E, F), and W, the right-hand side, is formed by `larger`, which
    returns one of its arguments without rounding, so F <= attainable <= W exactly and `max(E, F) <= W` exactly when
    `E <= W`.
  - 2F is the same expression in both versions; round(2F) >= F in any radix. So the decisions are the same bit for
    bit (numerics review: no path differs).
- **Each step is checked before it is tested.** The check is A.5's, applied to `current`, plus `attainable` finite
  and not negative. A method's broken step is thus `non_finite` before any arithmetic on it: in constant evaluation a
  NaN in `rel_tol |value|` would otherwise be a compile error, and an infinite `attainable` would accept anything.
  The library's steps always pass the check, so no result changes.
- **Termination is the method's:** the loop advances only while the step is not exhausted, and settles at most twice
  per advance. So it ends exactly when the method reports `exhausted` (D).
- **`settle`** is Gauss-Kronrod's check against sums formed afresh, which Phase 2 makes before every decision to stop
  (independent review, 2098416). The other methods' `settle` does nothing. It is required, not detected (Q15): a
  `settle` with a near-miss signature, or a combinator that forgets to forward it, would otherwise be skipped
  silently. Without it, 22 of the battery's 2,438 results change (simplicity review).

### C. What stays in the method, and why

The loop owns the user's request; the method owns its own limits, which no driver can know:
- **`ready`, the minimum work:** Romberg's level 4 (fewer samples of sin^2 on [0, 2 pi] all vanish); TanhSinh's
  level 4. Kept as a field, not folded into `start`, so that the steps before it are reported: the first levels are
  what a user debugging Romberg wants to see.
- **`attainable`, the least error it can certify:** Romberg's floor F; Gauss-Kronrod's and TanhSinh's 2F (a kink or
  a jump never brings every panel to its floor).
- **`exhausted`:** the caps (Romberg 16 levels; Gauss-Kronrod 64 bisections per bit of T; TanhSinh 10 levels),
  Gauss-Kronrod's stagnation (32 bisections without a new minimum of its excess), and its empty heap.
- **`trusted`:** Gauss-Kronrod's undone descent, a piece where fn is not finite that no short bisection avoided. Such
  a run cannot end as a success. It keeps Phase 2's behaviour of bisecting on to a better best before failing;
  exhausting at once would change those runs' best and their evaluation count.
- **`current` empty:** TanhSinh's level 0, which has no estimate yet, and its unbounded end term or an interval
  with no point strictly inside. The run then ends `not_converged` without a best, as `PHASE3_API.md` Q7 decided.

One common rule for these would change when the methods stop. A common minimum or floor, for example, would move
Romberg's or Gauss-Kronrod's results.

### D. The contract

Written once, in the concepts' Doxygen comments, as in draft 2.

- **Given to every method:**
  - a < b, both finite, possibly extreme: b - a may overflow, the ends may be subnormal or adjacent, or of opposite
    signs so that b - a rounds.
  - fn as an lvalue, valid for the call only: call it with `std::invoke` and a T, and convert the result with
    `T(...)`.
  - An estimating method also gets the criteria: finite, non-negative, `rel_tol <= 1`, `abs_tol` unclamped, possibly
    both 0. An iterative method never sees them.
- **Given to an iterative method:**
  - `advance` is never called on an exhausted step.
  - Nothing is called after a Failure.
  - `settle` may be called twice between steps.
- **Required of every method:**
  - fn is evaluated only at points of [a, b]; clamp a point formed as a + h t, and form a midpoint as a/2 + b/2.
  - A finite result (checked, A.5), or an error: `non_finite` without a best, or `not_converged` with or without one.
  - Its parameters are its own preconditions (Q9).
  - fn's exceptions, and the method's own, pass through.
- **Required of an iterative method:**
  - Its state keeps no reference, pointer or copy of fn: fn is passed to every call, and need not be the same object
    each time.
  - `start` and every `advance` leave `step` describing the state; the step passes the check of B.
  - It reports `exhausted` within finitely many steps. A method that never does loops forever at run time, and hits
    the compiler's limit in constant evaluation, as an estimating method that never returns would.
  - `settle` refines `step` without taking a step. It leaves `ready` and `exhausted` as they are, and calling it twice
    gives the same step as calling it once.
- **Required of an estimating method:** it terminates for every criteria it can be given, both tolerances 0
  included.
- **Required in meaning, not checkable:**
  - `abs_error` estimates |value - exact| and covers the rounding of the method's own sum.
  - `attainable` is an error the method can actually certify. A huge one would let it stop anywhere, as an estimating
    method may declare success at any `abs_error`.
  - For an estimating method, success means its own stop test was met. That is the one exception to "methods
    compute, integrate decides".
- **The library's methods only** (Q8): the rounding floor (for a != b), `abs_tol = 0` working, a best with every
  `not_converged` of Romberg and Gauss-Kronrod, fn never evaluated at an end (GaussLegendre, TanhSinh), and
  `non_finite` in constant evaluation where an overflow or a NaN would otherwise be a compile error.
- **Infinite ends** stay outside this contract (Q10); Phase 4 makes them opt-in per method.

### E. Observing (`observer`)

- **The forms:** `integrate<M>(qd::observer(o), fn, a, b[, criteria])`, with M named (`GaussKronrod<15>` by
  default), and `integrate(m, qd::observer(o), fn, a, b[, criteria])`, with a method value.
- **What it calls:** after the method's `start` and after each `advance`, `o(step)` with the method's `Step<T>`, for
  the ends in order (a < b), so a reversed integral reports values of the opposite sign to the result. It reports
  nothing after `settle`; the result `integrate` returns shows the settled values.
- **It shares the method's state** and passes fn on unchanged. An observed run is therefore the same instantiation
  of the method and computes exactly what an unobserved one does.
- **The observer:**
  - It must be callable as `o(Step<T> const &)` and return `void`. A `static_assert` names both requirements; a
    result is reserved for a stopping observer (deferred).
  - An lvalue observer is kept by reference, so what a reusable recorder records is there afterwards. An rvalue is
    kept by value.
  - It is called as a non-const callable, so a mutable lambda can count its calls.
  - Its exceptions pass through. That is also the way out of a user's method that never reports `exhausted`.
- **Counting is the caller's:** the evaluations are counted in fn (two lines; the spike's progress program does it),
  and the steps in the observer. The library's methods have exact counts there, which the tests assert: Romberg's
  report i at 2^(i+1) + 1 evaluations, Gauss-Kronrod's at 150 + 30 i without descents.
- **Where it works:** in both forms, in constant evaluation with a constexpr observer, and on any iterative method,
  a user's included. On a fixed rule or an estimating method it is rejected softly: no overload is viable, with no
  error inside the header.
- **Not in the protocol:** the observer cannot stop the run (Q14), and method-specific numbers (Gauss-Kronrod's
  panels, Romberg's tableau) are not reported.

### Example: printing each step

The state at each iteration, as the observer receives it: every field of the `Step`, with the step's number and the
evaluations so far, both counted by the caller. The program is the spike's `print_steps.cpp`, built against the
draft-5 header; the output below is its run (GCC 13.3, x86-64).

```cpp
#include <numerixx/quadrature.hpp>

#include <cmath>
#include <cstdio>

namespace qd = nxx::quadrature;

// Prints one Step: the iteration, the evaluations of fn so far, and every field of the step.
struct Printer
{
    std::size_t const *evaluations; // counted by fn
    int iteration = 0;

    void operator()(qd::Step<double> const &s)
    {
        std::printf("%4d %6zu  ", iteration++, *evaluations);
        if (s.current)
            std::printf("%22.17f %10.3g", s.current->value, s.current->abs_error);
        else
            std::printf("%22s %10s", "-", "-");
        std::printf(" %10.3g  %-5s %-9s %s\n", s.attainable, s.ready ? "yes" : "no", s.exhausted ? "yes" : "no",
            s.trusted ? "yes" : "no");
    }
};

void header(char const *title)
{
    std::printf("%s\n%4s %6s  %22s %10s %10s  %-5s %-9s %s\n", title, "step", "evals", "value", "abs_error",
        "attainable", "ready", "exhausted", "trusted");
}

template<typename R>
void footer(R const &r)
{
    if (r) std::printf("result: %.17f +- %.3g\n\n", r->value, r->abs_error);
    else std::printf("no result: reason %d\n\n", static_cast<int>(r.error().reason));
}

int main()
{
    std::size_t evaluations = 0;

    // Romberg on e^x over [0, 1] at rel_tol 1e-12: one step per level.
    auto const e = [&evaluations](double x) {
        ++evaluations;
        return std::exp(x);
    };
    header("Romberg, e^x on [0, 1], rel_tol 1e-12");
    Printer romberg{ &evaluations };
    footer(qd::integrate<qd::Romberg>(qd::observer(romberg), e, 0.0, 1.0, { .rel_tol = 1e-12 }));

    // Gauss-Kronrod on a peak of width 1e-3 at 0.3 over [0, 1] at rel_tol 1e-10: the mesh, then one step per bisection.
    evaluations = 0;
    auto const peak = [&evaluations](double x) {
        ++evaluations;
        double const d = x - 0.3;
        return 1.0 / (d * d + 1e-6);
    };
    header("GaussKronrod<15>, 1 / ((x - 0.3)^2 + 1e-6) on [0, 1], rel_tol 1e-10");
    Printer kronrod{ &evaluations };
    footer(qd::integrate<qd::GaussKronrod<15>>(qd::observer(kronrod), peak, 0.0, 1.0, { .rel_tol = 1e-10 }));
}
```

Output:

```
Romberg, e^x on [0, 1], rel_tol 1e-12
step  evals                   value  abs_error attainable  ready exhausted trusted
   0      3     1.71886115187659283       0.14   1.95e-14  no    no        yes
   1      5     1.71828268792475725   0.000578   1.92e-14  no    no        yes
   2      9     1.71828182879453029   8.59e-07   1.91e-14  no    no        yes
   3     17     1.71828182845907840   3.35e-10   1.91e-14  yes   no        yes
   4     33     1.71828182845904531   3.31e-14   1.91e-14  yes   no        yes
result: 1.71828182845904531 +- 3.31e-14

GaussKronrod<15>, 1 / ((x - 0.3)^2 + 1e-6) on [0, 1], rel_tol 1e-10
step  evals                   value  abs_error attainable  ready exhausted trusted
   0    150  3205.52671587585973612   1.22e+03   7.12e-11  yes   no        yes
   1    180  3162.35187766923809249        679   7.02e-11  yes   no        yes
   2    210  3119.17703946259962322        140   6.93e-11  yes   no        yes
   3    240  3128.71314052564457597        140   6.95e-11  yes   no        yes
   4    270  3127.97134058085566721       75.2   6.95e-11  yes   no        yes
   5    300  3137.50744164390016522       75.5   6.97e-11  yes   no        yes
   6    330  3136.76564169911125646       10.9   6.97e-11  yes   no        yes
   7    360  3136.79841664353580200       6.37   6.97e-11  yes   no        yes
   8    390  3136.83119158795989279        1.9   6.97e-11  yes   no        yes
   9    420  3136.83097700711550715      0.969   6.97e-11  yes   no        yes
  10    450  3136.83076242627112151     0.0421   6.97e-11  yes   no        yes
  11    480  3136.83076228576510402     0.0211   6.97e-11  yes   no        yes
  12    510  3136.83076214525908654   0.000117   6.97e-11  yes   no        yes
  13    540  3136.83076214526909098   5.86e-05   6.97e-11  yes   no        yes
  14    570  3136.83076214527955017   1.35e-07   6.97e-11  yes   no        yes
result: 3136.83076214527955017 +- 1.35e-07
```

How to read it:
- **Romberg** reports one step per level, starting at level 1 (3 evaluations); level k costs 2^k + 1. `ready` turns
  on at level 4 (step 3), the minimum before any decision (C). Step 3 does not meet the request (3.35e-10 against
  1e-12 times 1.72); step 4 does (3.31e-14), so the loop settles (nothing to do for Romberg) and returns it.
  `attainable` is the floor F, about 1.9e-14, far below the request, so it plays no part here. The result,
  1.71828182845904531, is 7.4e-17 from e - 1 = 1.71828182845904523536 (mpmath), well inside its estimate.
- **Gauss-Kronrod** reports the initial mesh of 10 panels as step 0 (150 evaluations), then one bisection per step,
  30 evaluations each. The estimate does not fall steadily: some bisections cut it by 23 to 430 times (steps 10,
  12 and 14), others halve it or leave it almost unchanged (steps 3 and 5), depending on which panel holds the
  largest estimate when it is split. `attainable` is twice the floor, about 7e-11. At step 14 the estimate,
  1.35e-7, meets the request, 1e-10 times 3137 = 3.1e-7: the loop settles (the sums formed afresh, here the same,
  as the result line shows) and returns it. The result, 3136.83076214527955, is 2.2e-11 from
  1000 (atan 700 + atan 300) = 3136.83076214530129340 (mpmath), within its estimate, which overstates the error by a
  factor of about 6,000, as Gauss-Kronrod's estimate usually does.
- **`exhausted` and `trusted`** stay no and yes throughout: neither run reaches a cap or stagnates, and no descent
  was undone. A run that runs out of steps (its cap, or stagnation) shows `exhausted` turning on in its last step.
- **The values are for the ends in order.** For b < a, `integrate` negates the result, not the reports (E).

### F. The library's methods

- **`GaussLegendre<N>::apply`:** unchanged from draft 2.
- **`Romberg`:** `start` forms T_0 from the ends and level 1; `advance` forms the next level; `settle` does nothing.
  The level is Phase 1's loop body, line for line, with the stop test removed. Its step:
  - current: {R(k, k), max(estimate, F)};
  - attainable F;
  - ready from level 4;
  - exhausted at level 16;
  - always trusted.
- **`GaussKronrod<15>`:** `start` builds the initial mesh, `advance` bisects the panel with the largest estimate, and
  `settle` forms the sums afresh. Phase 2's loop is kept line for line, with its locals as members of
  `detail::KronrodState<T>` and the stop test moved to the loop. Its orientation code is removed: `integrate` passes
  a < b, so it was dead. Its step:
  - current: {K, max(E, F)};
  - attainable 2F;
  - always ready;
  - exhausted by an empty heap, the cap or stagnation;
  - trusted unless a descent was undone.

### G. Phase 3b (TanhSinh) through the protocol

`PHASE3_API.md` stays approved, except:
- **Its "Public API" bullet:** `TanhSinh` is an iterative method, so the existing estimating forms of `integrate`
  take it:
  - `start`: level 0 with its outward sweep, reporting no current (level 0 has no estimate);
  - `advance`: the next level;
  - `settle`: nothing;
  - always trusted.
- **Its A:** "b < a negates" is done by `integrate`; TanhSinh receives a < b. The centre not strictly inside gives a
  start whose step is exhausted with no current: `not_converged` without a best.
- **Its B:** the state must own what it uses across `advance`. For a type without a table, in constant evaluation,
  the generated pieces are temporaries, so the state copies them, as Gauss-Kronrod copies its rule, or regenerates
  them. The rounding cut's suffix sums need every earlier level's table at each step, so regenerating repeats that
  cost; the phase measures which is cheaper.
- **Its E:** the stop test `E_k <= max(abs_tol, rel_tol |V_k|, 2 F_k)` from level 4 moves to the loop as ready
  (level >= 4) and attainable 2F_k. The cap of 10 levels is `exhausted`. An unbounded estimate is a step without
  current, so the best is absent (its Q7).
- **Observing:** `integrate<TanhSinh>(qd::observer(o), ...)` reports each level.

Its behaviour otherwise, the decisions Q1 to Q10 and the tests are unchanged. `docs/plan.md` splits Phase 3 into 3a
and 3b.

### H. What a method of one's own cannot use yet

As in draft 2: `Interval`, `Sum`, `rounding_floor`, the overflow checks for constant evaluation and the node tables
stay in `detail`; the guide page gives the floor's recipe (50 eps L1 from the method's own samples), which an
iterative method reports as its attainable error.

## Spike

`docs/open-closed/quadrature-spike/` (README there):
- `quadrature.patch.py ... normalize` writes draft 2 into a copy of `quadrature.hpp`.
- `steps.py` then writes drafts 3 to 5: `Step`, the stepped Romberg and Gauss-Kronrod, the loop, `iterative_method`
  and the observer. Its last two sections hold drafts 4's and 5's changes.

Built with GCC 13.3 and Boost 1.83, and the programs also with Clang 18, which needs `-U__cpp_concepts
-D__cpp_concepts=202002L` for libstdc++ 13's `<expected>`. The container's network refuses Boost 1.92, and the
project's five toolchains are not available here; the phase itself builds and tests on all five.

- **The existing quadrature tests** (`test_quadrature_gauss.cpp`, `_romberg.cpp`, `_kronrod.cpp`), built unchanged
  against draft 4, pass exactly as against today's header and every earlier draft. That is 25 of 26 test cases. The
  one failure is the same in every build: the `GaussLegendre<2>` long double node baseline, which was recorded on
  the project's toolchains and does not hold for GCC 13's 80-bit long double. The tests' `static_assert`s, Gauss-Kronrod
  and Romberg in constant evaluation among them, now run through the loop of B.
- **Bit for bit** (`dump.cpp`): the battery covers 8 integrands on 6 intervals (three reversed); 4 Gauss-Legendre
  rules; Gauss-Kronrod at 3 tolerances; Romberg at 3; float, double, long double, `cpp_bin_float_50` and
  `cpp_dec_float_50`; edges and three constant-expression results: 2,438 results.
  - Draft 1, and draft 2 without ordering the ends: identical to today.
  - Draft 2 with the ends ordered (Q4): 23 differ, all on reversed ends, with no status among them. Romberg's
    `abs_error` alone in 9 float results and 1 double result (the floor's L1 sum); in `cpp_dec_float_50`, 7
    Gauss-Legendre values and 6 Romberg results, in the guard digits (about the 70th significant digit).
    Gauss-Kronrod: none.
  - Draft 3: identical to draft 2 with the ends ordered.
  - Draft 4: identical to draft 2 with the ends ordered: the per-step checks, the required `settle` and the
    observer change no result.
  - Draft 5: identical to draft 2 with the ends ordered: the observed forms change no result.
- **Orientation today** (`orient.cpp`): 3,000 intervals (positive, negative, across 0, up to 2e6 wide), each
  integrated forwards and reversed, the whole result compared. Gauss-Legendre<64> and Gauss-Kronrod<15> give the
  exact negation every time. Romberg gives the exact negation in value and status every time, and a different
  `abs_error` in 374 of 6,000 runs.
- **The loop and observing** (`progress.cpp`; the `static_assert`s run in constant evaluation on GCC and Clang):
  - **Observed equals unobserved:** as `static_assert`s, in both forms, for Gauss-Kronrod (named by default),
    Romberg (reversed, named and as a value) and the user's trapezoidal rule. An observer on `GaussLegendre<7>` is
    rejected softly in both forms.
  - **Today's calls still resolve** with a generic lambda as fn (`[](auto const &x) { return x * x; }`), named, as
    a value and by default: no overload calls an argument on trial.
  - **Every branch of the loop,** with a scripted method whose steps come from a list:
    - not ready: it goes on;
    - passing but untrusted: `not_converged` with the best;
    - exhausted: `not_converged` with the best, and without one where current is empty, with `advance` never
      called there;
    - passing before settle and failing after: it goes on to the next step;
    - a NaN current, an infinite or a negative attainable: `non_finite`, without any arithmetic on the NaN.
  - **A user's iterative method** (a trapezoidal rule halving its step, about 45 lines, which meets D): x^3 - x over
    [-1, 2] gives exactly 9/4 + 27/2^26 at rel_tol 1e-6, its closed form at level 12 (the error is 3h^2/4 exactly,
    and every sum is exact in double), as a `static_assert`. On e^x at 1e-10 it stops at level 16, as derived.
  - **At run time:**
    - Gauss-Kronrod on a peak of width 1e-3 at rel_tol 1e-10: 15 steps, the evaluations exactly 150 + 30 i, 570 in
      all, and the result equal to the unobserved run's.
    - Romberg on e^x at 1e-12: a step per level through a mutable lambda, ready from level 4, the result at level 5.
    - A reusable recorder passed as an lvalue: it keeps its 16 steps afterwards.
    - An observer returning `bool`: the `static_assert` message, "observer must return void (a result is reserved
      for stopping the run)".
- **Printing each step** (`print_steps.cpp`): the example above, run as shown. Draft 5's syntax prints exactly what
  draft 4's `observed()` printed.
- **The bare observer, rejected** (`bare_observer.cpp`): built as is, a call that compiles today
  (`integrate<qd::Romberg>(fn, 0.0, 1.0, {.rel_tol = 1e-12})` with a generic lambda fn) compiles and runs. With
  `-DBARE`, which adds an overload for a bare callable told apart by whether it can be called with a `Step<T>`, the
  same call fails inside the user's lambda: "no match for 'operator*' (operand types are 'const
  nxx::quadrature::Step<double>' and 'const nxx::quadrature::Step<double>')" (Q14).
- **Methods of one's own** (`user_methods.cpp`, unchanged from draft 2, passing against draft 5): the fixed rule
  `Midpoint`, the combinator `SquareAtLowerEnd<Inner>` over the iterative library methods, the spy, the checked
  results, and the concepts' rejections. On 1/sqrt x and cos x / sqrt x over [0, 1] at 1e-12, the substitution takes
  150 evaluations against Gauss-Kronrod's 2,370 and 2,400.

## Tests (`tests/test_quadrature_protocol.cpp`, test cases `quadrature:protocol:*`)

Each in constant evaluation (`static_assert`) and at run time, unless it needs a library function or a type that is
not constexpr.

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
- **The loop's branches** with the scripted method:
  - every branch the spike shows;
  - passing, then failing after settle, while exhausted. This is Phase 2's double settle: settle is called twice,
    and the best is the settled step's.
  - exhausted before ready, with a best;
  - ready with an empty current while not exhausted (it goes on);
  - a Failure from `start`, `advance` or `settle`, returned as it is;
  - the scripted method run with an observer too, the branch where settle moves the step included, since that is
    where the observed method forwards `settle`.
- **An observer changes nothing:** `integrate<M>(qd::observer(quiet), ...) == integrate<M>(...)` and `integrate(M{},
  qd::observer(quiet), ...) == integrate(M{}, ...)` for Romberg and Gauss-Kronrod, reversed too.
- **Today's calls with a generic lambda as fn** resolve as before, in every form (no overload calls an argument on
  trial).
- **The reports:**
  - Romberg's report i at 2^(i+1) + 1 evaluations of a counting fn, its last step equal to the result.
  - Gauss-Kronrod's at 150 + 30 i on a case without descents, and its last step equal to the result on a mesh-only
    case, where the running and fresh sums add the same panels in the same order.
  - An lvalue recorder keeps its steps; a mutable lambda observer compiles and counts.
- **A user's iterative method** (the trapezoidal rule): x^3 - x over [-1, 2] exactly 9/4 + 27/2^26 at rel_tol 1e-6;
  e^x over [0, 1] at 1e-10 within its estimate (T_h - (e - 1) = (e - 1)(h^2/12 - h^4/720 + ...), plus the plain
  sum's rounding, about 1.2e-11 at level 16).
- Not tested (compile errors): an observer of the wrong signature or result, and a state without `step`. Their
  diagnostics are checked by hand in the spike.

## Documentation

As in draft 2, and:
- `quadrature.hpp`: the comments of `Step`, `state_t`, `iterative_method`, `Observer` and `observer()`, and of the two
  observed forms of `integrate`. The contract of D goes in the concepts, with the usual causes of a failed concept:
  a missing `const`, the exact return types, `settle`, and the `step` member. The file comment says that a method
  computes and `integrate` decides, and shows `qd::observer` in its example.
- `docs/quadrature/design.md`: "The interface" records the loop and what stays in the method (C); BRIEF.md 3.9
  marked as amended.
- `docs/sphinx/quadrature/reference.md`: `Step`, `state_t`, the three concepts, `Observer` and `observer`.
- The guide page on writing a method:
  - the three forms, with a user's iterative method as its main example;
  - `qd::observer` for debugging: a table, a recorder, counting evaluations in fn, slicing fn's recorded points
    between steps to see where Gauss-Kronrod bisects;
  - choosing at compile time whether a method can be observed: `if constexpr (qd::iterative_method<M, double>)`;
  - an observable combinator that forwards `start`, `advance` and `settle` to an inner method, keeping
    `state_t<Inner, T, G>` and rebuilding its function G in each call;
  - throwing from the observer as the way out of a method that never ends.

## Decisions for the user

From draft 2 (unchanged unless noted):
- **Q1. Member names:**
  - `apply` (fixed rule) and `estimate` (estimating method);
  - draft 3 adds `start`, `advance` and `settle` (iterative method).

  Recommended: yes.
- **Q2. Concepts:**
  - `fixed_rule`, `estimating_method` and `iterative_method`, each `<M, T, F = T(T const &)>`, with no copy
    requirement;
  - draft 4 adds the alias `state_t<M, T, F>`.

  Recommended: yes.
- **Q3. A method value goes first,** `integrate(m, fn, a, b[, criteria])`, beside today's named forms, which keep
  their signatures. Recommended: yes. Draft 5 adds the two observed forms (Q14), six forms in all.
- **Q4. `integrate` puts the ends in order and negates.** Recommended: yes.
  - **Cost:** 23 of the battery's 2,438 results change, all on reversed ends, no status. That covers Romberg's
    `abs_error` in float and double, and values beyond `cpp_dec_float_50`'s 50 digits.
  - **Alternative:** the method orients itself, bit for bit with today.
- **Q5. The result's finiteness checked in every build,** a broken result returned as `non_finite`; draft 4 checks
  each step of the loop too (B). Recommended: yes.
- **Q6. The validated criteria passed as `Criteria<T>`** to an estimating method. Recommended: yes.
- **Q7. The protocol as Phase 3a, before TanhSinh (3b).** Recommended: yes.
- **Q8. The statements that hold only for the library's methods say so.** Recommended: yes.
- **Q9. A method's parameters are its own preconditions;** no new enumerator now. Recommended: yes.
- **Q10. Finite ends only;** infinite ends opt-in per method in Phase 4. Recommended: yes.

From draft 3, as revised:
- **Q11. `integrate` owns the loop and the test against the request; the method keeps its own limits** (`ready`,
  `attainable`, `exhausted`, `trusted`, `settle`; C). Recommended: yes. The results are unchanged, and the test is
  `BRIEF.md` D5's.
  - **Alternative, a common rule for the limits:** one floor multiple and one minimum. That changes when Romberg and
    Gauss-Kronrod stop.
- **Q12. `Step`:** `current` optional, `attainable`, `ready`, `exhausted`, `trusted`, with the neutral values as
  defaults. Recommended: yes.
  - **`current` is optional** for TanhSinh's Q7.
  - **The defaults:** a method that forgets `trusted` or `ready` is then not stuck at `not_converged` (usage sketch).
  - **The names are open to change:** `attainable` could be `floor`, but Gauss-Kronrod's is 2F; `trusted` could be
    `acceptable`.
- **Q13. The one-member `estimate` form stays beside the iterative one,** for methods that do not step: a
  recursion, a combinator over a library method, a port of an existing routine. Such methods are not observed.
  Recommended: yes. An iterative combinator takes about three times the code of the 15-line `SquareAtLowerEnd`.
- **Q14. Observing:** `integrate<M>(qd::observer(o), fn, a, b[, criteria])` and `integrate(m, qd::observer(o), fn,
  a, b[, criteria])`, the user's syntax and choice of both forms (draft 5). The observer receives the `Step`, after
  the start and after each step, and cannot stop the run. Recommended: yes.
  - **The wrapper** (`qd::observer(o)`, of type `Observer<O>`), the user's decision after the comparison below:
    - The overloads are told apart by the argument's type. A bare callable can only be told apart by trying to call
      it with a `Step<T>`, which for a generic lambda instantiates its body: with the constraints in their natural
      order, a call that compiles today failed inside the user's integrand (`bare_observer.cpp`). Ordering the
      constraints so that fn and the method are checked first avoids it for today's overloads, but leaves
      correctness to that order.
    - It names the argument's role: two callables side by side are told apart by position only.
    - A swapped observer and fn fail to compile either way; the bare form's first error is inside the integrand.
    - A stopping observer, if it comes, can be a second wrapper rather than a meaning guessed from a return type.
    - Cost: about a dozen characters per call, one small public type and one function.
  - **The Step** (draft 4, changed from draft 3's `Progress{iteration, evaluations, step}`): draft 3 counted
    evaluations by wrapping fn, which made an observed run a different instantiation, brought a `detail` type into
    users' argument-dependent lookup, and left a method that kept fn in its state with a dangling pointer. The
    counts are two lines of the caller's code (E).
  - **Not chosen:**
    - draft 4's public adaptor `integrate(observed(m, o), fn, ...)`, now internal;
    - a bare callable (above);
    - a trailing observer after a defaulted `criteria`;
    - draft 3's `Progress`.
  - **Deferred:** a stopping observer; a public stop test, `meets(step, criteria)`, for an observer that shows the
    distance to the request or for one's own loop over `start` and `advance`; method-specific numbers.

New in draft 4:
- **Q15. `settle` required of every iterative method,** `return {};` where there is nothing to refine.
  Recommended: yes.
  - **Why:** detected, as in draft 3, a `settle` whose signature is a near miss (not `const`, say), or a combinator
    that forgets to forward it, is skipped without a diagnostic. Gauss-Kronrod without `settle` changes 22 of the
    battery's 2,438 results.
  - **Alternative:** detected, with the trap documented.

## Decisions from the reviews

Draft 4's round (on draft 3):

From the numerics review:
- **Adopted:**
  - the state keeps no reference, pointer or copy of fn (D), and draft 3's counting wrapper, which made such a state
    dangle under `observed`, is gone (Q14);
  - every step checked in the loop, `current` as A.5 and `attainable` finite and not negative, before any arithmetic
    (B);
  - the tests: the trapezoid's exact closed form, 3h^2/4, and exact equality at level 12; exact evaluation counts;
    named Gauss-Kronrod cases instead of "equal unless settle moved it"; the missing branches; `static_assert`
    twins;
  - the trapezoid's midpoint as a/2 + b/2, so it stays in [a, b] near max (D);
  - no type in `detail` passed to users' methods (the wrapper is gone);
  - the reports for the ends in order (E);
  - G: level 0 without current, always trusted, `settle` nothing, the tables owned by the state;
  - D: settle idempotent and leaving `ready` and `exhausted`; the loop's promises under "Given"; the consequence of
    never exhausting;
  - E and Summary: a mutable observer, the `current` check in the example, "the value forms".
- **Confirmed:**
  - the stepped Romberg and Gauss-Kronrod are Phases 1 and 2 on every path, including the double settle and
    constant evaluation;
  - the request test's equivalence, including E < F and F = 0;
  - termination;
  - G preserves every decision of `PHASE3_API.md`, and nothing in it uses the tolerance inside a step;
  - the observer changes no result.

From the simplicity review:
- **Adopted:**
  - the observer receives the `Step`, `Progress` dropped (Q14);
  - `settle` made required (Q15) rather than detected;
  - Gauss-Kronrod's dead orientation code removed;
  - `observed`'s type moved to `detail`, `observed()` kept (in draft 5 `observed()` itself moved to `detail`, behind
    `qd::observer`);
  - "the value forms";
  - calling members directly stated as outside the contract, and the deferred item restated as a public stop test;
  - a type of two kinds promised nothing;
  - the example's `current` check;
  - the scripted method through `observed`;
  - the battery's result filled in;
  - "the one exception" for an estimating method.
- **Not adopted:**
  - `ready` folded into `start`, which would hide Romberg's first three levels from the reports (C);
  - dropping the trapezoid's closed-form check, kept as one exact `static_assert` at its stopping level, which the
    numerics review derived.

From the usage sketch:
- **Adopted:**
  - an lvalue observer kept by reference, called as non-const, with a `static_assert` naming its signature and
    requiring `void` (E);
  - fn valid for the call only (D);
  - `Step`'s defaults (Q12);
  - `observed()` on a method that does not step rejected softly;
  - a hang's way out documented (E);
  - `state_t` public (Q2);
  - `attainable` checked (B);
  - the guide's patterns (Documentation).
- **Not adopted:**
  - a public `meets(step, criteria)` now (deferred, Q14);
  - accepting a `step()` member function beside the data member;
  - a library-wide safety cap.

Draft 2's round (on draft 1):

From the numerics review:
- **Adopted:**
  - fn evaluated only in [a, b] required of every method, and the spike's methods clamp and refuse an overflowing
    width;
  - the orientation of methods tied to an end (Q4);
  - termination for every criteria, and parameters as preconditions;
  - the edge inputs and constant evaluation in the contract;
  - Q8 extended;
  - the tests' tolerances derived, a spy, an expression-template instantiation, a reversed combinator, the domain
    check;
  - `orient.cpp` comparing whole results.
- **Confirmed:** the checked part for all four methods; the members' arguments; the clamped criteria; Q4's
  derivation.

From the simplicity review:
- **Adopted:**
  - the methods after `detail` with inline members;
  - no `std::copy_constructible`;
  - the contract once, in the header;
  - "a method is of one kind";
  - the tests trimmed;
  - "an lvalue".
- **Not adopted:**
  - Q5's "no check";
  - dropping the second integrand of the substitution.

From the usage sketch:
- **Adopted:**
  - the method value first (Q3);
  - the ends ordered (Q4);
  - the release check (Q5);
  - parameters (Q9);
  - infinite ends (Q10);
  - `<M, T, F>` with F defaulted;
  - no copy requirement;
  - the concepts' diagnostics;
  - the floor's recipe.
- **Not adopted:** a precedence between the kinds; publishing `rounding_floor` now.

## Deferred

- The building blocks (`Interval`, `Sum`, `rounding_floor`, the overflow checks, the tables) as public API (H).
- A stopping observer, a public stop test (`meets`), and method-specific numbers in the reports (Q14, E).
- Infinite ends as an opt-in (Q10): Phase 4.
- `Adaptive<Pair>`, the globally adaptive loop over any embedded pair (`ANALYSIS.md` 6.3).
- The same protocol for `nxx::roots` and `nxx::deriv` (`ANALYSIS.md` 10). In roots the tolerance enters Brent's and
  TOMS 748's steps (`docs/roots/BRIEF.md` 3.8), so there the loop cannot own the test alone.

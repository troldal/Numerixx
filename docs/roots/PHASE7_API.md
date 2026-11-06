# Phase 7 API proposal: the method protocol

Status: approved by the user on 2026-11-06, Q1 to Q19 as recommended (draft 2): the one-call concepts deferred (Q14),
F's fix to today's header in this phase (Q16), and the state-by-reference clause in quadrature's `iterative_method`
too (Q17). Draft 2 followed draft 1's review round (numerics-reviewer, simplicity-critic, usage-sketcher; "Decisions
from the reviews"). History:
- **Draft 1:** the stepped and the one-call protocol, `Step` with an `exhausted` flag, an open method required to keep
  a sign change once found.
- **Draft 2:** the review round folded in. `Step` is five numbers, and a method that cannot go on returns
  `not_converged` from `advance`. An open method may lose a sign change it found. Every bracketed step must lie
  strictly inside the last, so every method ends. Bisection's stop test is never later than today's. The contract
  states the two placement rules a convergent method needs. The concepts reject a state taken by value. A latent bug in
  today's header (a bound above 2 tol in types without subnormals) is fixed (F). The questions the reviews raised for
  the user are added (Q14 to Q19).

The decisions it builds on:
- `docs/open-closed/ANALYSIS.md` (2026-11-01): the methods open to users (D1), quadrature first (D4); its section 10
  step 2 is this phase.
- `docs/quadrature/PHASE3A_API.md` (approved 2026-11-02): the quadrature protocol, the model this phase follows. A
  method computes and reports; the entry point owns the loop and the stop test; `observer(o)` watches the steps.
- `docs/roots/BRIEF.md` (approved 2026-10-25), amended in two places: 3.1 ("methods as a closed set of empty tags")
  and 3.8 (no customization).
- The user's choices of 2026-11-05, on the assessment of the roots code (U1 to U4, below).

The evidence is a spike in `docs/open-closed/roots-spike/` (Spike, below).

## The user's choices (2026-11-05)

These four were asked after an assessment of `roots.hpp` and answered as recommended. They are recorded with the
reasons given then.

- **U1. Both protocols, designed together and built in steps.** The stepped protocol decides names and types that
  the one-call protocol also uses: what a method receives (`Tolerance`, `Start`), what `solve` checks, and what a
  `Step` holds. Publishing the one-call protocol alone could force a breaking change later. The one-call protocol
  stays as the analogue of quadrature's `estimating_method`, for a method that does not step. (Revised by Q14, on the
  simplicity review: designed here, shipped later.)
- **U2. `solve` verifies the sign change itself.** "Every result is a verified sign change with a bound"
  (`design.md`, Phase 2 P1) is the module's central promise. With the check it holds for any method that reports its
  bracket truthfully, not only for the library's carefully written ones. It costs two sign comparisons and one
  `distance_bound` per step and no evaluation of fn, and it keeps `distance_bound` and its edge cases inside `solve`.
  Quadrature cannot do this: its error estimate has to come from the method.
- **U3. The method value first,** `solve(m, fn, bracket)`, as `integrate(m, fn, a, b)`. A user who has learned one
  module knows the other. The alternative in `ANALYSIS.md`, the method last after a defaulted `criteria`, makes a
  caller write `solve(fn, b, {}, Pegasus{})`.
- **U4. The observer moves to the shared core,** `nxx::observer` and `nxx::Observer` in `math.hpp`, with
  using-declarations in `nxx::quadrature` and `nxx::roots`. It does not depend on the module: it wraps a callable and
  passes it the module's step. `docs/design.md`'s rule, "a shared core grows only when a second module needs
  something", applies, and `nxx::optimize` would be the third. Since `qd::observer` and `rt::observer` are then the
  same entity, using-directives for both namespaces do not make it ambiguous, unlike `Criteria` and `error`.

## Summary

```cpp
#include <numerixx/roots.hpp>
namespace rt = nxx::roots;

auto const f = [](double x) { return x * x - 2.0; };

// Unchanged: every call that compiles today compiles, with the same results (but Bisection's, Q3, and the fix of F).
auto const r = rt::solve<rt::Brent>(f, rt::Bracket{ 1.0, 2.0 });
auto const s = rt::solve<rt::Secant>(f, 1.0);

// New: a method as a value, first.
auto const t = rt::solve(rt::Toms748{}, f, rt::Bracket{ 1.0, 2.0 }, { .rel_tol = 1e-12 });

// New: each step reported, for any stepped method, the library's or a user's.
auto const print = [](rt::Step<double> const &s) {
    std::println("[{}, {}] point {}", s.lower, s.upper, s.point);
};
auto const o = rt::solve<rt::Brent>(rt::observer(print), f, rt::Bracket{ 1.0, 2.0 });

// New: a method of one's own. It reports its bracket after each step; solve() checks it and decides when to stop.
struct Trisection
{
    template<nxx::real T> struct State { T a, fa, b, fb; rt::Step<T> step; };
    template<nxx::real T, typename F>
    constexpr std::expected<State<T>, rt::error> start(F &fn, rt::Start<T> const &s, rt::Tolerance<T> const &) const;
    template<nxx::real T, typename F>
    constexpr std::expected<void, rt::error> advance(F &fn, State<T> &s, rt::Tolerance<T> const &) const;
};
static_assert(rt::bracketing_method<Trisection, double>); // names the member at fault, where solve() would not
auto const u = rt::solve<Trisection>(f, rt::Bracket{ 1.0, 2.0 });
```

- **Two public concepts:**
  - `bracketing_method`: `start(fn, start, tolerance)` and `advance(fn, state, tolerance)` on a bracket, each leaving
    a `Step` in `state.step`;
  - `open_method`: the same from a guess, `start(fn, guess, tolerance)`.

  The one-call protocol (a method that does the whole computation in one member) is designed here and deferred (Q14).
- **A stepped method computes; `solve` decides.** After the start and after each step the method reports the two
  points it holds, fn's values there, and its point. `solve` checks the report, returns the point once it meets the
  request, and forms `abs_error` itself. A method that cannot go on returns `not_converged`.
- **Every method ends.** Each bracketed step must lie strictly inside the last, and an open method's steps without a
  sign change count against the iteration cap, so no method, a user's included, can loop for ever.
- **The library's eleven methods become stepped methods.** The closed lists `detail::bracketing` and `detail::open`,
  the table `values_needed` and the `if constexpr` dispatch in `solve` are removed. Their arithmetic is unchanged.
- **`rt::observer(o)`** marks o as an observer: `o(step)` after the start and after each step of any stepped method.
- **From quadrature, unchanged in spirit:** validation stays in `solve`; a method is named, `solve<M>(fn, ...)` as
  today, or given as a value, first; no ADL; everything constexpr; nothing throws.

## Public API (changes to `roots.hpp` and `math.hpp`)

```cpp
namespace nxx {  // math.hpp, moved from nxx::quadrature (U4)
    template<typename O> struct Observer { /* the callable: by reference for an lvalue, by value for an rvalue */ };
    template<typename O> constexpr Observer<O> observer(O &&o);
}

namespace nxx::roots {

    enum class error { non_finite, no_sign_change, not_converged, invalid_step };  // invalid_step new (Q5)

    // Bracket, Root, Criteria: unchanged.

    /// The tolerances after validation, tol(x) = absolute + relative |x|; promoted from detail, unchanged.
    template<real T> struct Tolerance { T absolute; T relative; constexpr T at(T const &x) const; };

    /// A validated bracket: lower < upper, fn's values finite, nonzero, of opposite signs; K values at each end.
    template<real T, std::size_t K = 1>
    struct Start { T lower; std::array<T, K> at_lower; T upper; std::array<T, K> at_upper; };   // names: Q15

    /// What a stepped method reports after its start and after each step.
    template<real T>
    struct Step
    {
        T lower; T f_lower;       // the two points it holds, in order, and fn's values there
        T upper; T f_upper;
        T point;                  // its estimate: lower, upper or a point between them
    };

    /// How many values a method needs from fn at a point: a static integral M::values if declared (Newton 2,
    /// Halley 3), else 1.                                                                             name: Q15
    template<typename M> inline constexpr std::size_t values_of = /* ... */;

    template<typename M, typename T, typename F = /* T (*)(T const &), or std::array<T, K> for K values */>
    using bracketing_state_t = /* what m.start(fn, Start<T, values_of<M>>, Tolerance<T>) returns, inside expected */;
    template<typename M, typename T, typename F = /* as above */>
    using open_state_t = /* what m.start(fn, T guess, Tolerance<T>) returns, inside expected */;

    template<typename M, typename T, typename F = /* as above */>
    concept bracketing_method = requires(M const &m, F &fn, Start<T, values_of<M>> const &s, Tolerance<T> const &t) {
        { m.start(fn, s, t) } -> /* std::expected<S, error> */;
    } && requires(M const &m, F &fn, bracketing_state_t<M, T, F> &st, Tolerance<T> const &t) {
        { m.advance(fn, st, t) } -> std::same_as<std::expected<void, error>>;
        { std::as_const(st).step } -> std::same_as<Step<T> const &>;
        requires !requires { m.advance(fn, std::move(st), t); };   // the state by reference, not by value
    };
    template<typename M, typename T, typename F = /* ... */> concept open_method = /* the same from a T guess */;

    using nxx::Observer;
    using nxx::observer;

    // The eleven methods: types with const member templates; Newton and Halley declare `values` and work from both
    // starts. Bisection, Brent, Toms748, AndersonBjorck, Ridders, RegulaFalsi, Illinois, Newton, Halley, Secant,
    // Steffensen.

    // solve: five lines, each on a Bracket or from a guess (12 overloads, 6 of them new; Q8).
    solve<M>(fn, bracket | guess[, criteria]);                          // today
    solve<M>(fn, dfn[, d2fn], bracket | guess[, criteria]);             // today, keyed on values_of<M>
    solve(m, fn, bracket | guess[, criteria]);                          // new
    solve<M>(observer(o), fn, bracket | guess[, criteria]);             // new
    solve(m, observer(o), fn, bracket | guess[, criteria]);             // new
}
```

## Behaviour

### A. What `solve` does

On a bracket, in order, as today:
1. The tolerances validated and clamped (`non_finite` for a negative or non-finite one).
2. The bracket validated, put in order, fn evaluated at its ends and, for a bracket around 0, at 0 (zero first).
   `non_finite`, `no_sign_change`, or a root found here, as today.
3. The method: `start`, then the loop of B, inside the caller's bracket.

From a guess: the tolerances, then the guess (`non_finite` if it is not finite; today each method checked it), then
the method, as above but with no caller's bracket.

### B. The loop

```
s = start(...)                                           (an error ends the run)
loop:
    check s.step: all five numbers finite (else non_finite); lower <= point <= upper; inside the caller's bracket
        (on a bracket)                                   (else invalid_step)
    point an end where fn is exactly 0:                  return {point, tol(point)}
    bracketed (fn nonzero at both ends, opposite signs):
        lower < upper, and strictly inside the last bracketed step, if the last step was one   (else invalid_step)
        larger = max(distance_bound(point, lower), distance_bound(upper, point))
        if distance_bound(lower, upper) <= 2 tol(point) or larger <= tol(point):   return {point, larger}
    not bracketed:
        on a bracket:                                    invalid_step (the sign change was lost)
        from a guess: the step counts against the cap;   at the cap (4 times the bits of T): not_converged
    advance(...)                                         (an error ends the run; not_converged where it cannot go on)
```

- **The stop test is `Criteria`'s documented rule,** "the bracket at most 2 tol(x) wide", for every method, or the
  point within tol of both ends (Q3).
  - Brent, Toms748, AndersonBjorck, Ridders, RegulaFalsi, Illinois and Newton and Halley on a bracket test exactly
    the width today, at the end with the smaller |f|, which is their point. The larger distance from an end is the
    width, and `distance_bound` is symmetric, so their stop and their `abs_error` are the same bit for bit; the second
    clause implies the first for them.
  - Bisection tests today the larger distance from its midpoint against tol(midpoint), which is the second clause.
    The first can only stop it earlier: where the midpoint rounds, the half-width exceeds half the width, and the
    width meets 2 tol a halving sooner. Each such result is within 2 tol, `Root`'s bound. Without the second clause,
    the unit `distance_bound` adds to a width that Sterbenz's lemma does not cover could also stop it later than today,
    and in one case turn a returned pole into `non_finite` (numerics review, finding 2).
- **Every method ends** (numerics review, finding 3). A bracketed step after a bracketed step must lie strictly inside
  it: neither end outside, not both the same. A sequence of strictly nested brackets of a finite set of numbers is
  finite. From a guess, steps without a sign change count against the cap, so a method cannot alternate for ever
  either. The library's methods meet the rule: each evaluates a point strictly inside its bracket and keeps one end.
  Draft 1 left termination to the method, as quadrature does; roots can check it, because its steps are brackets.
- **No step without progress returns a root.** A method that cannot place a point strictly inside its bracket, or
  cannot form a step from a guess, returns `not_converged` from `advance`. The library's "no point strictly inside"
  exits returned that bracket as a root today. In IEEE types they are not reached once the stop test has failed:
  tol(b) >= 2 eps |b| is at least two units of T's spacing at b (denorm_min at 0), so a failed test leaves more than
  four units between the ends, and Brent's next point (under 3/4 of the width by interpolation, half by bisection,
  tol(b) toward c) and Bisection's midpoint lie strictly inside. In types without subnormals they were reached (F);
  with F's fix the methods make progress there too, and any case left returns `not_converged` rather than a root
  whose bound exceeds 2 tol.
- **The cap counts only the steps without a sign change** from a guess, as today: Newton's, Halley's and Steffensen's
  loops of 4 times the bits of T iterations, and Secant's, whose start evaluates two points. After the library's
  hand-over there is no cap, as on a bracket (`design.md`, "No iteration cap"). The cap now also bounds a user's open
  method, and an open method may lose a sign change it found (simplicity review, finding 3): every success is still a
  bracketed step that met the test.

### C. What stays in the method, and why

- **The tolerance inside the iteration.** Unlike quadrature's methods, which never see the criteria, a root finder's
  steps depend on the request: Brent's minimum step is tol(b); Toms748 and AndersonBjorck keep their points 0.7 tol(u)
  inside the ends (`place`); the open methods step at least tol(x). Without it Brent's termination argument fails, and
  a one-sided method stalls at the stop test (D). So `start` and `advance` receive the validated `Tolerance`. The stop
  test is still `solve`'s alone.
- **Its point.** The end with the smaller |f| (Brent, Toms748 and its family; Brent's b after its swap) or the
  midpoint (Bisection). The choice is the method's, ties included.
- **The hand-over of the library's open methods.** Newton and Halley continue in Brent's loop with their own step,
  Secant and Steffensen in Brent's. The method's state holds the finisher's state and forwards to it. A user's open
  method may do the same with a library method's public members, or not hand over at all.
- **When it cannot go on.** No point strictly inside, or no usable step: f' = 0 or not finite, a step that would
  overflow, Secant's or Steffensen's values equal. It returns `not_converged`.

### D. The contract

Written once, in the concepts' Doxygen comments, and shown on the guide page.

- **Given to every method:**
  - the validated `Tolerance`;
  - fn as an lvalue, valid for the call only: evaluate it as the library does (`std::invoke`, the values converted
    to T at once), and keep no reference, pointer or copy of it in the state;
  - on a bracket, a `Start`: lower < upper, finite, fn's value finite, nonzero and of opposite signs at the ends, the
    caller's bracket after zero first; K = `values_of<M>` values at each end, the derivatives unchecked;
  - from a guess, a finite guess.
- **Given to a stepped method:** `advance` is called only on a step that does not meet the request; nothing is called
  after an error.
- **Required of every method:** fn evaluated only inside the caller's bracket, on a bracket; a non-finite value of fn
  returned as `non_finite`; fn's exceptions, and the method's own, pass through.
- **Required of a stepped method:**
  - the state taken by reference: `start` and every `advance` leave `step` describing it, with finite numbers and
    lower <= point <= upper;
  - on a bracket, every step a bracket inside the caller's and strictly inside the last: lower < upper, fn's values at
    the ends nonzero and of opposite signs, or point an end where fn is exactly 0;
  - f_lower and f_upper are fn's values at lower and upper: its true values, not a method's scaled ones (Pegasus's).
    `solve` cannot check this; it is what makes the result a verified sign change;
  - `not_converged` from `advance` where it cannot go on.
- **Required for convergence, not checked** (usage review, finding 1): on a bracket, keep each new point about
  tol(point) inside the ends, or a method that moves one end stalls at the stop test, the other end never moving
  (textbook Pegasus took 33 evaluations on the cubic, kept tol inside 10); from a guess, steps of at least tol(x), or
  an iterate converging from one side never brackets (a Chebyshev step: `not_converged` after 213 evaluations, with
  the minimum step 5). The library's methods do both.
- **The arithmetic assumed:** IEEE round-to-nearest with subnormals, or a type that flushes below min() as Boost's do.
  Flush-to-zero modes (`-ffast-math`) are outside the derivations (numerics review, finding 8).

### E. Observing (`observer`)

As quadrature's (`PHASE3A_API.md` E), with the roots `Step`:
- the forms `solve<M>(rt::observer(o), fn, bracket | guess[, criteria])` and `solve(m, rt::observer(o), fn, ...)`;
- `o(step)` after the method's start and after each successful `advance`, for the ends in order; an observed run is
  the same instantiation of the method and computes exactly what an unobserved one does;
- an lvalue observer kept by reference, an rvalue by value; called as a non-const callable; its exceptions pass
  through; it must return void (a `static_assert` says so), a result being reserved for a stopping observer
  (deferred);
- from a guess it sees the iterates (lower and upper the last two, equal after the start; point the latest) and,
  after the hand-over, the bracket.

### F. The library's methods, and a fix to today's header

Each `detail` loop of Phases 1 to 6 becomes a state and one step, its arithmetic and order of evaluations unchanged,
its stop test moved to the loop:

| Method | State | One step | Point |
|---|---|---|---|
| `Bisection` | the bracket and fn's values | fn at the midpoint, the half with the sign change kept | the midpoint |
| `Brent`, `Newton`, `Halley` (bracket) | Brent's a, b, c, d, e with K values each | Brent's iteration, the next b evaluated, then the top of his loop (c reset, b and c swapped) | b |
| `Toms748` | the bracket subroutine's state, the stage of Algorithm 4.2, the bracket at the iteration's start | one evaluation: the secant step, the two interpolation steps, the double-length secant step or the bisection | the end with the smaller \|f\| |
| `AndersonBjorck`, `RegulaFalsi`, `Illinois` | the retained and latest points, the scaled value, the guard's counters | one evaluation | as Toms748 |
| `Ridders` | the bracket subroutine's state, and the midpoint just evaluated | the midpoint, or the fit's point | as Toms748 |
| `Newton`, `Halley` (guess) | the iterate with K values, then Brent's state | the derivative step, or Brent's after the hand-over | the latest iterate, then b |
| `Secant` | the last two iterates, then Brent's state | the secant step, or Brent's | the same |
| `Steffensen` | the iterate, then Brent's state | one iteration (up to two evaluations), or Brent's | the same |

- Toms748 and Ridders evaluate several times per iteration and tested after each evaluation; one step is one
  evaluation, so their loops become small state machines over the stage. The numerics review checked each stage
  against today's loop.
- The internal strategy switches by type identity (`std::same_as<M, Brent>` inside `brent<M>`, the factor chain in
  `false_position<M>`) become internal enumerations (`Safeguard`, `Scaling`). AndersonBjorck, RegulaFalsi and Illinois
  derive from one internal template, Newton and Halley from another; their public types and names are unchanged.

**The fix (numerics review, finding 1; Q16).** In a type without subnormals (`cpp_bin_float`, `cpp_dec_float`, whose
denorm_min() is min() and which flush results below min() to 0), today's header returns a bound above 2 tol, against
`Root`'s documented bound, on a bracket of one sign whose ends lie less than 2 min apart: on [min, 2.5 min] around 1.25
min at the default criteria, every bracketing method returns min with abs_error 2.5 min against 2 tol = 2 min (and
Toms748's family reaches an assert in Debug builds). (b - a) / 2 flushes to 0, so the midpoint is an end, and Brent's m
is 0, so his minimum step goes the wrong way. Two changes:
- `midpoint`, where a + (b - a) / 2 is not strictly inside, tries (a + b) / 2: both ends are then tiny, so the sum fits
  (a guard keeps it from being formed near max, where constant evaluation would reject the overflow). Between adjacent
  IEEE numbers no point exists, so nothing changes there.
- Brent's minimum step takes its direction from b < c rather than m > 0: the same wherever m is not 0.

Bit for bit wherever today's code makes progress; the case becomes a test for the nine bracketing methods in both
Boost types, et_on and et_off.

### G. The call forms (Q8)

Five lines on the reference page, each on a bracket or from a guess: 12 overloads, today's 6 and 6 new.
- `solve<M>(fn, bracket | guess[, criteria])`: today's;
- `solve<M>(fn, dfn[, d2fn], bracket | guess[, criteria])`: today's four, their constraint changed from
  `std::same_as<M, Newton>` (or `Halley`) to `values_of<M> == 2` (or 3), so a user's method that declares `values`
  gets them;
- `solve(m, fn, bracket | guess[, criteria])`: a method given as a value, for a method with parameters;
- `solve<M>(rt::observer(o), fn, ...)` and `solve(m, rt::observer(o), fn, ...)`.

The separate-derivative forms stay named only: with a method value or an observer, fn is one callable returning the
values (the roots pages already show the pair-returning lambda around `derivative_of`). The usage review found the gap
a trap: `solve(m, f, df, bracket)` fails with "'const lambda' is not derived from 'const Bracket<T>'". Adding the four
value-first separate-derivative forms (without an observer) closes it, 16 overloads in all; see Q8.

A named form needs a default-constructible method (`std::default_initializable<M>`); a method with required
parameters is given as a value.

## Spike

`docs/open-closed/roots-spike/` (README there): `assemble.py` writes the proposal's headers from today's and three
chunks. GCC 16.2 (the Codespace; the project's other toolchains are not available there). Draft 2's evidence, on the
spike as revised after the review round:
- **The existing tests, unchanged, against the spike's headers:** every roots test file, the quadrature protocol,
  examples and Romberg files (the observer moved), and `derivative_of`'s: 52 test cases and 1,374 assertions pass, including the `static_assert`s that pin evaluation counts in
  constant evaluation. One test names
  `detail::Tolerance`; the spike keeps an alias for it, and the phase moves the test to the public name.
- **The manual's programs** (`docs/sphinx/examples/roots_*.cpp`, `basics_*.cpp`, 12), which compare the output their
  pages quote, evaluation counts included, and assert their pages' claims: all pass.
- **Bit for bit** (`dump.cpp`): 17 functions (smooth, multiple, tiny and huge roots, a pole, a step, brackets around 0,
  brackets up to +-1e308, removable singularities), their brackets and guesses, five tolerances (the defaults,
  rel_tol 1e-8, abs_tol 1e-6 alone, both huge, rel_tol 1/4 alone), every method, and the separate-derivative forms;
  in float, double, long double, `cpp_bin_float_50` and `cpp_dec_float_50`: 10,930 results, each with its count of
  evaluations. 
  - **10,923 are identical** to today's, value, `abs_error`, error and count: every result of Brent, Toms748,
    AndersonBjorck, Ridders, RegulaFalsi, Illinois, Newton and Halley (both starts), Secant and Steffensen; the
    separate-derivative forms; F's fix changes none of them.
  - **7 differ, all Bisection** (Q3), each one halving earlier, where the final midpoint rounds; none stops later than
    today, the second clause of the stop test keeping today's stop wherever it came first:
    | Case | Type | Today: n, abs_error | Spike: n, abs_error |
    |---|---|---|---|
    | x^2 - 2 on [0, 10], defaults and rel_tol 1e-8 | float | 27, 2.38e-7 | 26, 3.58e-7 |
    | x^2 - 2 on [0, 10], defaults | double | 56, 4.44e-16 | 55, 6.66e-16 |
    | x^2 - 2 on [0, 10], defaults | long double | 67, 2.17e-19 | 66, 3.25e-19 |
    | e^x - 2 on [-700, 700], defaults | long double | 75, 1.08e-19 | 74, 1.63e-19 |
    | x^2 - 2 on [0, 10], defaults | `cpp_bin_float_50` | 171, 5.35e-51 | 170, 1.60e-50 |
    | e^x - 2 on [-5, 5], defaults | `cpp_bin_float_50` | 172, 7.41e-51 | 171, 8.02e-51 |

    Each new `abs_error` is within 2 tol(value), `Root`'s documented bound (in double, tol(1.414) = 6.3e-16, so 2 tol =
    1.26e-15 > 6.66e-16): the bracket is at most 2 tol wide, which is the rule `Criteria` states. The last halving
    today buys a bound tighter than requested.
- **The numerics review's cases** (`numerics_findings.cpp`): on [min, 2.5 min] around 1.25 min at the default criteria, with F's
  fix, every bracketing method returns a value within 2 tol = 2 min of the root, in `cpp_bin_float_50` (Bisection 1.75
  min, abs_error 1 min; Brent, Toms748, AndersonBjorck and Ridders min, abs_error min) and in `cpp_dec_float_50` (abs_error
  0.75 min, Brent's 2 min), in a Debug build as with NDEBUG: no assert is reached. Today's header returns abs_error
  2.5 min for all five; draft 1's spike returned `not_converged`. Bisection's two cases return today's results: x - 2.4
  on [1, 3] at rel_tol 0.5 {2, 1}, and 1/x on [-1, 1] at abs_tol 1 the pole at 0 with abs_error 1.
- **The loop and observing** (`protocol.cpp`, all in constant evaluation and again at run time), passing:
  - every branch of B with a scripted method whose steps come from a list: the stop at the first step that meets
    the request, `abs_error` at an end and at an interior point, an exact zero at the point, fn 0 at an end that is not
    the point, a NaN or an infinity (`non_finite`), a point outside its points, a step outside the caller's bracket, a
    lost sign change on a bracket, opposite signs at one point, the same bracket twice, a bracket not nested in the last
    (each `invalid_step`), one end kept and the other moved (nested), an error from `advance`; from a guess, the cap
    (in double 212 steps without a sign change after the start: a bracketed 211th or 212th step succeeds, a 213th is too
    late), a sign change found, lost and found again (allowed), and two bracketed steps not nested (`invalid_step`);
  - observed equals unobserved for all eleven methods, named and as a value, on a bracket and from a guess; Bisection
    reports once after its start and once per evaluation;
  - the concepts accept the library's methods for the right starts, all eleven passing the state-by-reference clause,
    and reject `int`, Brent from a guess, Secant on a bracket and a method taking its state by value; `values_of` is 1
    for a method with a data member or a static table called `values`; today's calls with a generic lambda as fn
    resolve, named and as a value.
- **Methods of one's own** (`sketches/`, the usage review, on draft 1): 13 programs, among them Pegasus and a
  trisection on a bracket, a Chebyshev step from a guess with Brent as its finisher, one-call methods, a parameter given
  as a value, misuse and its diagnostics, constant evaluation, `cpp_bin_float_50` (Decisions from the reviews).

## Tests (`tests/test_roots_protocol.cpp`, test cases `roots:protocol:*`)

Each in constant evaluation (`static_assert`) and at run time where the type allows:
- the existing suite unchanged (the one `detail::Tolerance` moved to `rt::Tolerance`);
- every method named and as a value giving the same result;
- the loop's branches with the scripted method, as `protocol.cpp`, with the numerics review's additions: cap - 1, cap
  and cap + 1 steps without a sign change, then a bracketed one; a bracketed open step with lower == upper;
- the maximum evaluation counts from the cap: Newton and Halley 1 + cap, Secant 2 + cap, Steffensen 1 + 2 cap;
- Bisection's stop test: the cases of the numerics review's finding 2 (x - 2.4 on [1, 3] at rel_tol 0.5 returns {2, 1};
  1/x on [-1, 1] at abs_tol 1 returns the pole at 0 with abs_error 1, as documented), an odd-unit bracket, a saturated
  tolerance; and the battery's changed results, each within 2 tol;
- F's fix: [min, 2.5 min] around 1.25 min and around 2 min (x / min - 2), the nine bracketing methods in
  `cpp_bin_float_50` (et_on and et_off) and `cpp_dec_float_50`, each against its closed-form root with abs_error at
  most 2 tol;
- `distance_bound` symmetric, as `static_assert`s in the five types (pairs of one sign, of opposite signs above max / 4,
  with 0, near min); brackets with -0.0 ends;
- observing: the counts of reports, an lvalue recorder, a mutable lambda;
- user methods: Pegasus (kept tol inside), a trisection, a user method with `values = 2` through the
  separate-derivative form, each checked against a closed form;
- the concepts: the library's methods for the right starts; `int`, an empty tag, a non-const `advance`, a state taken
  by value rejected; `values_of` on a data member called `values`; today's calls with a generic lambda as fn.

## Documentation

- `roots.hpp`: the file comment (a method computes, `solve` decides; methods of one's own); `Tolerance`, `Start`,
  `Step`, `values_of`, the state aliases and the concepts, with the contract of D in the concepts; each method's
  members; the new forms of `solve`; `invalid_step`, documented as "the method broke its contract", so that a caller
  of someone else's method knows whom to blame; `not_converged` no longer "an open method" only.
- `math.hpp`: `Observer` and `observer`, their comments moved from `quadrature.hpp`.
- `docs/roots/design.md`: "The interface" records the protocol, the loop, the stop test, the progress rule and the
  arithmetic assumed; `BRIEF.md` 3.1 and 3.8 marked as amended.
- `docs/open-closed/ANALYSIS.md`: section 11 records this phase.
- The roots chapter of the manual: the reference page gains the new entities; a guide page, "A method of your own",
  built and run as a test, with:
  - a stepped method (Pegasus) and a trisection with a parameter, given as a value;
  - the two rules for convergence (D) with the numbers that show why;
  - the recipes a correct step needs (usage review, finding 5): the overflow-free midpoint, the sign test without a
    product (`(fa < 0) != (fb < 0)`, where fa fb underflows), the secant point with the quotient first, fn's true
    values in the step;
  - `static_assert(rt::bracketing_method<M, double>)` beside the method, whose diagnostic names the member at fault
    where `solve`'s does not (usage review, finding 3);
  - a finisher held in the state (`std::optional<rt::bracketing_state_t<rt::Brent, T>>`), as one way to end an open
    method;
  - `rt::observer` for debugging.

  The quadrature pages that name `qd::observer` stay correct (it is the same function).

## Implementation steps

Draft 1 split the work by start; draft 2 adds the observer as its own first step (simplicity review, finding 8), since
it changes `quadrature.hpp` and is checked on its own by quadrature's tests:
0. **The observer in `math.hpp`,** with using-declarations in `nxx::quadrature`; and `iterative_method`'s
   state-by-reference clause (Q17). Quadrature's tests unchanged.
1. **On a bracket:** the public types, the stepped concepts, the loop and its checks; F's fix; the nine bracketing
   methods stepped; the value-first and observed forms on a bracket; the separate-derivative forms on `values_of`;
   `tests/test_roots_protocol.cpp` for all of it. The guess path keeps its dispatch until step 2 (the only step
   where two mechanisms coexist).
2. **From a guess:** the four open methods stepped with the hand-over in their state; the cap in the loop; the forms
   from a guess; the documentation, including the guide page.

Each step is built and tested on the five toolchains, locally and in CI, with the docs build clean, before the next.

## Decisions for the user

From draft 1, as revised:
- **Q1. Names.** The members `start`, `advance` (stepped); the concepts `bracketing_method`, `open_method`; the
  aliases `bracketing_state_t`, `open_state_t`. Recommended: yes. `start` and `advance` are quadrature's; Newton
  overloads `start` on its argument (a `Start` or a guess), as `solve` itself is overloaded. (The one-call names: Q14.)
- **Q2. Public types.** `Tolerance` promoted unchanged; `Start<T, K>` without today's internal `root` member (a root
  found by validation never reaches a method); `Step<T>` with five numbers; `values_of`. Recommended: yes. (Their
  names: Q15.)
- **Q3. The stop test:** the bracket at most 2 tol(point) wide, `Criteria`'s rule, or the point within tol(point) of
  both ends; `abs_error` formed by `solve`. Recommended: yes.
  - **Cost:** the same 7 of the battery's 10,930 results as in draft 1, all Bisection on wide brackets, each one
    evaluation fewer with a larger `abs_error` still within 2 tol (Spike); none later than today, and the numerics
    review's two cases return today's results. No test and no page of the manual changes.
  - **Alternative:** a test that depends on where the point lies (the width at an end, the larger distance inside),
    bit for bit with today. It would ask a user's method with an interior point (a trisection) for more than `Criteria`
    states.
- **Q4. A method that cannot go on returns `not_converged` from `advance`;** `Step` has no flag (simplicity review,
  finding 2). Recommended: yes. Draft 1's `exhausted` duplicated the error channel: in roots a failure carries no best,
  unlike quadrature's `Failure`, and the observer saw the unchanged step again (numerics review, finding 6).
- **Q5. A broken step or result is a new error, `invalid_step`,** documented as "the method broke its contract".
  Recommended: yes. The library's methods never return it. Quadrature returns `non_finite` for its broken steps, all of
  whose checks are of finiteness; here the checks are of position, sign and progress, and `non_finite` would mislead.
  Non-finite numbers in a step stay `non_finite`. `error`'s documentation already announces new values.
- **Q6. The library's open methods hand over inside their state; the cap is in the loop,** counting the steps without
  a sign change. Recommended: yes. The finisher stays each method's choice; a parameter for it (`Secant{.finish =
  Toms748{}}`) is deferred.
- **Q7. The guess checked by `solve`,** not by each method. Recommended: yes; no result changes.
- **Q8. Twelve overloads** (G), the separate-derivative forms named only, the reference and the guide saying plainly
  that a method value or an observer takes one callable returning the values. Recommended: yes, as the simplicity
  review. Alternative: the four value-first separate-derivative forms too (16 overloads), which closes the trap the
  usage review found (`solve(m, f, df, bracket)`).
- **Q9. The observer in `math.hpp`** (U4). Recommended: yes.
- **Q11. The strategies inside the library become internal enumerations and base templates** (F). Recommended: yes;
  publishing them as policies (`FalsePosition<Scaling>`) stays deferred (`ANALYSIS.md` 6.4).
- **Q12. The helpers stay in `detail`:** `midpoint`, `distance_bound`, `same_sign`, the overflow-safe quotient and
  product, `place`, `evaluate`. With `solve` forming the bound (U2), a user's method must keep its points inside, keep
  them tol from the ends, and evaluate fn. The guide page gives the recipes (Documentation). Recommended: yes for this
  phase. The usage review proposes two public helpers instead: a point kept tol inside a step's ends (protocol
  knowledge, it argues, not arithmetic), and `values<K>(fn, x)` for a K-valued fn; both deferred, to add when a user
  asks.
- **Q13. Four implementation steps** (Implementation steps). Recommended: yes.

New in draft 2:
- **Q14. The one-call concepts:** in this phase, as U1 decided, or deferred (simplicity review, finding 1).
  Recommended: defer, with the design kept in this document.
  - **For deferring:** they are the one path where a result is the method's word rather than a verified sign change,
    which U2 was chosen to guarantee; U1's purpose, no later break, is met by fixing `Tolerance`, `Start` and `Root`
    now, and adding them later only widens `solve`'s constraint; no method in view needs them (quadrature's had its
    combinator); the usage review's one-call sketches (a hand-written bisection, `search` then a finisher) were as
    easy stepped or as plain functions.
  - **If kept:** the member `estimate` rather than `solve` (which collides with the free function of another
    contract); results checked against `Root`'s bound (numerics review, finding 4); a type with both kinds of members
    outside the contract; step 3.
- **Q15. Two names.** `values` (the member a method declares) or `values_needed`; `Start` or another name (`Ends`,
  `Bracketed`) with members matching `Step`'s (simplicity review, findings 4 and 5). Recommended: keep `values` and
  `Start`. `values_of` now counts only a static integral `values`, so a table of that name is not misread (numerics
  review, finding 7); `Start`'s members are arrays of K values where `Step`'s are single values, so different names
  keep them apart; `start(fn, s, t)` reads as quadrature's `start(fn, a, b)`.
- **Q16. F's fix to today's header** (numerics review, finding 1), in this phase. Recommended: yes: today's header
  breaks `Root`'s documented bound in types without subnormals, and the protocol's loop would otherwise turn those
  cases into `not_converged`.
- **Q17. The concepts reject a state taken by value** (usage review, finding 2), in roots; and the same clause in
  quadrature's `iterative_method`, which has the same hole (a change to `quadrature.hpp`, step 0). Recommended: yes for
  both.
- **Q18. Every bracketed step strictly inside the last** (numerics review, finding 3): termination for every method,
  checked by `solve`. Recommended: yes.
- **Q19. An open method may lose a sign change it found** (simplicity review, finding 3): every success is still a
  verified sign change; a user's open method need not embed a finisher. Recommended: yes.

Draft 1's Q10 (the one-call concepts beside the stepped ones) is now Q14.

## Decisions from the reviews

Draft 1's round: `numerics-reviewer`, `simplicity-critic` and `usage-sketcher` (the last two run as general-purpose
agents with the role in their brief, since the Codespace defines only `numerics-reviewer`). Each finding was recorded as
it arrived, with its assessment, and the outcome is noted after it: adopted in draft 2, or a question for the user.
Draft 2's spike was verified again after the changes (Spike).

**simplicity-critic:**
1. **Defer the one-call concepts** (`bracketing_solver`, `open_solver`). They are the one path where a result is the
   method's word rather than a verified sign change, which U2 was chosen to guarantee. U1's purpose (no later break) is
   met by fixing `Tolerance`, `Start` and `Root` now; adding the concepts later only widens `solve`'s constraint.
   Quadrature's `estimating_method` had a client (`SquareAtLowerEnd`); roots names none. If kept: the member `estimate`
   (not `solve`, which collides with the free function of another contract), a type with both kinds of members outside
   the contract, one naming scheme. *Assessment: a question for the user, since it revisits U1 (both protocols); the
   argument about U2 is sound. Outcome: Q14, recommended defer.*
2. **Drop `Step::exhausted`;** `advance` returns `error::not_converged` where it cannot go on. In roots a failure has no
   best to carry, unlike quadrature's `Failure`. The spike already sets the flag inside `advance` and returns success
   with the step unchanged; the loop re-tests a step it has already tested. Removing it deletes a public member, Q4 and
   two clauses of the contract. *Outcome: adopted (Q4); the numerics review's finding 6 (a repeated report) goes with it.*
3. **Do not require an open method to keep its sign change once found.** `solve` returns only on a bracketed step at
   most 2 tol wide, so every success is still verified. The rule forces every user's open method to embed a bracketing
   finisher, and a plain Newton-like iteration fails with `invalid_step` the first time a bracketed pair is followed
   by one that is not. Unbracketed steps keep counting against the cap, so only a method that stays bracketed without
   narrowing can loop forever, as on a bracket. *Outcome: adopted (Q19), with the numerics review's progress rule between consecutive bracketed
   steps (Q18), so that termination does not depend on keeping the sign change.*
4. **`values` is a generic member name** (a method with a table called `values` would be misread inside `values_of`);
   `values_needed`, today's internal name, is clearer. Or keep the hook in `detail` and publish it later.
   *Assessment: keep it public, since a user's Halley variant (Schröder, Chebyshev) is the protocol's
   purpose. Outcome: the name is Q15; the misreading is fixed by the numerics review's finding 7.*
5. **`Start` reads badly in `m.start(fn, start, t)`,** and its members differ from `Step`'s (`at_lower` against
   `f_lower`); every one-value method writes `s.at_lower[0]`. A name that is not the member's (`Ends<T, K>`,
   `Bracketed<T, K>`), with members matching `Step`. *Outcome: Q15, recommended keep.*
6. **The state aliases public only if users need them:** with finding 3 adopted, only a user's combinator or the
   deferred finisher parameter would. *Outcome: kept public: a user's open method may still hold a library finisher (the usage review did).*
7. **Q8's count is inconsistent:** the spike has 12 overloads, 6 of them new, and the alternative would be 24, not
   18. Present them as the reference page would, five lines each taking a bracket or a guess. *Outcome: adopted (G).*
8. **Move the observer first, as its own step,** since it touches quadrature only and is checked in isolation; say
   whether step 1 includes the one-call concepts; avoid two dispatch mechanisms coexisting between steps.
   *Outcome: adopted (Implementation steps, step 0; the one-call concepts step 3, if kept).*
9. **Keep `invalid_step`,** documented as "the method broke its contract", so that a caller of someone else's method
   knows whom to blame; its causes shrink with findings 1 and 3. *Outcome: adopted (Q5, Documentation).*

**usage-sketcher:** 13 programs in `docs/open-closed/roots-spike/sketches/`, each compiled with the
tests' warnings and run, none using `detail`, every result checked against a closed form (sqrt 2, the root of x^3 - 2x -
5, cos x = x, 1/sqrt(ln 2), a 50-digit root in `cpp_bin_float_50`). Today's calls, the value form, observing (lvalue and
rvalue, generic, the `static_assert` messages), the run-time checks (every broken method it wrote was caught), constant
evaluation (observed runs included), Boost types and one-call combinators work as proposed. Its findings:
1. **The contract leaves out the two rules a convergent method needs.** D says only "the bracket narrows to 2 tol
   within finitely many steps".
   - On a bracket, a new point must stay about tol inside the ends, or a one-sided method stalls at the stop test:
     textbook Pegasus took 33 evaluations on the cubic against AndersonBjorck's 9, its latest point at the root while
     the other end stayed put; kept tol(point) inside, 10. A user's safeguarded Newton (`values = 2`) took 23, then
     7; a section at 1/100 of the bracket returned `not_converged`.
   - From a guess, a step must be at least tol(x), or an iterate converging from one side never brackets: a
     Chebyshev step (`values = 3`) returned `not_converged` after 213 evaluations on the cubic from 2; with the minimum
     step, 5, as Halley.

   Proposal: state both rules in D and the concepts' comments, show them on the guide page, and reconsider Q12 for
   placement: a public helper keeping a point tol inside a step's ends is protocol knowledge, not arithmetic.
   *Outcome: the rules in D and on the guide page; the helper deferred (Q12).*
2. **A state taken by value in `advance` (a forgotten `&`) satisfies the concepts and loops forever** on a copy that
   never changes (`s8b_by_value.cpp`; reproduced: a million steps on an unchanged bracket). Proposal: the concepts
   also require that `advance` not accept an rvalue state, `!requires { m.advance(fn, std::move(st), t); }`, which
   rejects a state taken by value or by const reference; Brent, Toms748, Bisection, Ridders and Illinois pass it.
   *Outcome: adopted in roots, all eleven methods passing (Spike); quadrature's `iterative_method` has the
   same form and the same hole, and the same clause there is Q17.*
3. **Mistakes found through `solve` never name the member at fault:** a non-const `advance`, a missing `step` or a
   wrong `start` give 96 lines of "no matching function" ending in "no operand of the disjunction is satisfied";
   `static_assert(rt::bracketing_method<M, double>)` beside the method names the invalid expression in 25. The guide
   page should recommend it (as quadrature's `user_methods.cpp` does). Today's misuse diagnostics also grew from 6
   candidates to 12. *Outcome: adopted (Documentation, Summary).*
4. **Q8's gap is a trap:** `solve(m, f, df, bracket)`, written by a user who learned "the method first", fails with
   "'const lambda' is not derived from 'const Bracket<T>'". Proposal: the four value-first separate-derivative forms,
   or say plainly that a value or an observer needs one callable returning the values. *Outcome: Q8.*
5. **The arithmetic helpers belong on the guide page at least.** A correct step needed four of `detail`'s pieces: an
   underflow-free sign test (`fa * fx < 0` on 1e-200 times the cubic gave `invalid_step`), a quotient-first secant
   (the textbook form took 1,001 evaluations on [-1e308, 1e308], quotient first 6), the overflow-free midpoint, and
   `Step` holding fn's true values rather than a method's scaled ones. *Outcome: the recipes on the guide page
   (Documentation); publishing the helpers stays Q12.*
6. **Evaluating a K-valued fn is fiddly** (a pair, tuple or array of at least K elements; structured bindings break
   on a longer tuple). Proposal: a public `values<K>(fn, x)` (today's `detail::evaluate`), or document `std::get`.
   *Outcome: deferred, with Q12.*
7. **The open method's hand-over works with the public API alone** (`std::optional<rt::bracketing_state_t<rt::Brent,
   T>>`, `rt::Brent{}.start(...)`, `.advance(...)`, a 3-valued fn passed to the 1-valued Brent), but needs explaining:
   the iterates ordered in `Step`, the finisher's report copied after each step, `bracketing_state_t`'s fn type.
   *Outcome: with Q19 the hand-over is optional; the page shows it as one way (Documentation).*
8. **Smaller:** `exhausted` without a new evaluation is the same as returning `not_converged` from `advance` (as
   simplicity finding 2); "inside the caller's bracket" is checked only after fn was called outside it; a
   non-narrowing method can be stopped from the observer only by throwing or exiting (an argument for the deferred
   stopping observer); the `Tolerance` to `Criteria` round trip in a combinator could be a conversion.

**numerics-reviewer** (findings 1 and 2 reproduced with `numerics_findings.cpp`, against today's header
and the spike's, GCC 16.2, with and without NDEBUG):
1. **B's derivation is wrong in types without subnormals; the "unreachable" exit is reachable, and today's header
   breaks `Root`'s bound there.** `cpp_bin_float` and `cpp_dec_float` flush results below min() to 0 and have
   denorm_min() == min(). With ends of one sign whose difference d lies in (min, 2 min) and is not below the smaller
   end, `distance_bound` returns d + min > 2 min, while the default tol is min (rel |x| flushes to 0): the stop test
   fails, yet (b - a) / 2 flushes to 0, so `midpoint` returns an end; Brent's m is 0, its step rejected, and b - tol is
   outside the bracket; `place` returns the midpoint, an end. Reproduced on [min, 2.5 min] with fn(x) = x / min - 1.25
   (root 1.25 min), default criteria, in both `cpp_bin_float_50` and `cpp_dec_float_50`:
   - **today**, Bisection, Brent, Toms748, AndersonBjorck and Ridders all return value min with abs_error 2.5 min,
     above 2 tol(min) = 2 min, `Root`'s documented bound; in a Debug build, Toms748 and its family reach the assert in
     `enclose` (reported, not run here);
   - **the spike** returns `not_converged` for all five.

   Proposed fixes, bit for bit wherever today's code makes progress: `midpoint` falls back on (a + b) / 2 where a +
   (b - a) / 2 is not strictly inside (both ends then tiny, so the sum cannot overflow; no such point exists between
   adjacent IEEE numbers, so nothing changes there); Brent's minimum step takes its direction from `b < c` rather than
   `m > 0` (the same wherever m is not 0). Then rewrite B's derivation, relax `narrow`'s assert, and test the case for
   the nine bracketing methods in both Boost types. *Outcome: adopted in the spike (F), verified (Spike); a fix to today's header,
   Q16.*
2. **Bisection's new stop test differs where the midpoint is exact too, and can turn a root into an error.** The raised
   width (`distance_bound` adds a unit where Sterbenz does not apply) can exceed 2 tol where today's exact halves meet
   tol. Reproduced: x - 2.4 on [1, 3] at rel_tol 0.5, today {2, 1} with no interior evaluation, the spike {2.5, 0.5}
   after one; 1/x on [-1, 1] at abs_tol 1 (zero first skipped, fn(0) infinite), today {0, 1}, the pole returned as
   documented, the spike `non_finite` after evaluating fn(0). Where the midpoint rounds, the spike stops a halving
   earlier (the battery's 7). Under saturation (tol above max / 2) the width test always passes. Proposed: stop where
   the width is at most 2 tol *or* the larger distance from the point to an end is at most tol, never later than
   today. *Outcome: adopted (B, Q3); for the end-point methods the second clause implies the first, so they are
   unchanged; both cases return today's results (Spike).*
3. **No progress check:** a stepped method that reports the same bracketed step for ever passes every check. Proposed:
   each bracketed step after a bracketed step strictly nested in it (neither end outside, not both equal), else
   `invalid_step`. Every library method keeps one end and evaluates strictly inside, and strictly nested brackets of a
   finite set of numbers end: termination for any method, for two comparisons. *Outcome: adopted (B, Q18), between consecutive bracketed steps, the cap bounding the rest (simplicity
   finding 3).*
4. **The one-call check does not enforce `Root`'s contract:** an `abs_error` above 2 tol(value) passes, and a
   negative one is reported as `non_finite`. Proposed: both `invalid_step`. *Outcome: adopted in the spike, for the case
   that Q14 keeps the one-call concepts.*
5. **Documentation that becomes false:** `error::not_converged` says "an open method", and bracketing methods can now
   return it; B leans on `narrow`'s assert. *Outcome: adopted (B, Documentation).*
6. **The observer sees a repeated step** when `advance` only sets `exhausted`. *Outcome: gone with Q4.*
7. **`values_of` picks up any member named `values`** (`std::array<double, 3> values;` is then a hard error).
   Proposed: require a static integral constant. *Outcome: adopted (Spike: a data member and a static table
   called `values` count 1).*
8. **Two assumptions to pin:** `distance_bound` symmetric in every type (by case analysis under round-to-nearest; to
   be asserted, `cpp_dec_float` included), and IEEE round-to-nearest with subnormals, or types that flush like Boost's
   (FTZ/DAZ under -ffast-math breaks the arguments again). *Outcome: adopted (Tests; D, "The arithmetic assumed").*

Checked and found correct: Bisection's points and order; Brent's start, normalisation, steps, minimum step, zero and
hand-over; Toms748's stages (`first`, a0 and b0, the shifts, the halving test before the reset); Ridders' stages; false
position's report from `narrow`'s ordered bracket, ties to the upper end; the open methods' caps (Newton, Halley and
Steffensen `cap` advances, Secant 2 + cap evaluations) and hand-over order; the guess checked after the tolerances;
the bound formed by `run`; the rules for user types. Its additions to the test plan: finding 1's case for the nine
bracketing methods in both Boost types (et_on and et_off); finding 2's cases; `distance_bound`'s symmetry; the
scripted cases of a repeated step, a non-nested bracket, fn 0 at an end that is not the point, cap - 1, cap and cap + 1
unbracketed steps; the one-call results; the maximum counts from the cap (Newton 1 + cap, Secant 2 + cap, Steffensen
1 + 2 cap); brackets with -0.0 ends; `values_of` on a non-static `values`.

Right as proposed, in the simplicity review's view: `solve` owning the stop test and `abs_error`; the guess checked by `solve`; the
tolerance passed to the methods; the hand-over inside the library's open methods; the cap in the loop; `Tolerance`
promoted rather than `Criteria` passed; `Step`'s point; the internal strategies and helpers; the observer limited to
stepped methods; the separate-derivative forms keyed on the count of values; the concepts shaped like quadrature's.

## Deferred

- The safe-arithmetic helpers as public API (`ANALYSIS.md` 6.6; Q12).
- Policy templates (`FalsePosition<Scaling>`, `Safeguarded<Step>`; `ANALYSIS.md` 6.4).
- A finisher parameter for the open methods (Q6).
- A stopping observer, and a public stop test (`meets(step, tolerance)`), as in quadrature. The usage review found
  that a user's method that never narrows could be stopped from an observer only by throwing; Q18's progress rule now
  ends it with `invalid_step`.
- A public helper keeping a point tol inside a step's ends, and `values<K>(fn, x)` for a K-valued fn (Q12).
- A conversion from `Tolerance` back to `Criteria`, for a combinator that calls `solve` (usage review, finding 8).
- The one-call concepts (Q14): `bracketing_solver` and `open_solver`, a member `estimate(fn, start | guess,
  tolerance)` returning `std::expected<Root<T>, error>`, the result checked by `solve` for finiteness, position and
  `Root`'s bound (the numerics review's finding 4, in the spike), a type with both kinds of members outside the
  contract. Adding them later widens `solve`'s constraints and breaks nothing.
- The separate-derivative forms with a method value or an observer (Q8).

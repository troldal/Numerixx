# Phase 3 API proposal: Ridders' extrapolation

Status: approved by the user on 2026-10-21 (uniform API, restart, `Estimate<T>` implicitly convertible to T, the
untrusted-result restart and the exact-agreement stop) and implemented. The design went through two review rounds
(simplicity-critic, usage-sketcher, numerics-reviewer, a compile spike) and five prototypes.

## Summary

Ridders' method is one more method of `derivative`, and every method returns the same type: a value with an optional
error estimate, implicitly convertible to the value. Ridders fills in the estimate; the stencils leave it empty, since
they cannot estimate their error without extra evaluations.

```cpp
using namespace nxx::deriv;
auto const f = [](double x) { return std::exp(x); };

auto const a = derivative(f, 1.0);                // Central<2>; a->abs_error is empty
auto const b = derivative<Central<4>>(f, 1.0);    // the same type
auto const c = derivative<Ridders>(f, 1.0);       // c->value is e to about 1e-15, *c->abs_error about 6e-14
double const v = *c;                              // the conversion to T
```

Every call is greppable as `derivative<`; the method is the first template argument, as for the stencils.

## Public API (`deriv.hpp`)

```cpp
namespace nxx::deriv {

    /// A derivative estimate: the value, and an estimate of its absolute error when the method provides one.
    template<real T>
    struct Estimate
    {
        T value;
        std::optional<T> abs_error; // an estimate of |value - exact|, not a guaranteed bound; empty for stencils

        /// Implicit, so an Estimate can be used where a T is expected: double d = *derivative(fn, x).
        constexpr operator T() const { return value; }

        /// Compares both members; without it, == would compare only the values through the conversion.
        friend constexpr bool operator==(Estimate const &, Estimate const &) = default;
    };

    /// Ridders' extrapolation (tag): central differences at h, h/2, h/4, ..., extrapolated until rounding stops the
    /// improvement, with an error estimate.
    struct Ridders
    {};

    // Central<P>, Forward<P>, Backward<P> as in Phase 2.

    /// The derivative of fn at x by method M with step h: for a stencil the step, for Ridders the initial step.
    /// Precondition: h > 0 and h finite.
    template<typename M = Central<2>, real T, typename F>
    requires detail::method<M> && detail::callable<F, T>
    [[nodiscard]] constexpr std::expected<Estimate<T>, error> derivative(F &&fn, T const &x,
                                                                        std::type_identity_t<T> const &h);

    /// The same with the method's default step: eps^(1/(P+1)) * max(|x|, 1) for a stencil of order P,
    /// max(|x|, 1) / 32 for Ridders.
    template<typename M = Central<2>, real T, typename F>
    requires detail::method<M> && detail::callable<F, T>
    [[nodiscard]] constexpr std::expected<Estimate<T>, error> derivative(F &&fn, T const &x);

} // namespace nxx::deriv
```

- `detail::method<M>` is closed: `detail::stencil<M> || std::same_as<M, Ridders>`, not a trait users could
  specialise. `detail::callable<F, T>` is the constraint `derivative` already spells out
  (`std::invocable<F &, T const &>`, result convertible to T), written once.
- `error` keeps its one value, `non_finite`. Its documentation says that it may gain values; until then a `switch` on
  it should have a `default`.
- `Estimate` is a plain aggregate (designated initialisers work, confirmed on all four toolchains) and stays
  `{value, abs_error}`: adding a field later would break structured bindings and comparisons.
- No new requirement on T: the Ridders tableau is filled from `T(0)`, so T need not be default-constructible.

### Implementation structure

The two existing overloads stay the only public ones, each with one `if constexpr (std::same_as<M, Ridders>)`:
- the step overload calls `detail::ridders(fn, x, h)`, or the current stencil body, moved to a detail function, whose
  result is wrapped as `Estimate{d, std::nullopt}`;
- the default-step overload picks `max(|x|, 1) / 32` or `default_step<order>`, sharing the `max(|x|, 1)` scale.

The two-node central difference is pulled into one detail helper returning the value and the rounding floor `r` (the
floor needs fn's two values, which the current branch discards). It has two uses: the two-node stencils and every
Ridders level, which thereby share the overflow guard and the finiteness checks.

### Why `std::optional` for the error

- **An estimate from every method** would mislead: a stencil can only bound its truncation error with extra
  evaluations (Boost.Math spends one or two per call), and the rounding part alone is smallest exactly when
  truncation is largest.
- **A result object that computes the error on request** would have to keep fn and x, so it could dangle, and would
  evaluate fn behind the user's back.
- **Different return types per method** (the second revision) are the surprise the uniform API removes.

An empty optional says exactly what is known. Code that relies on an estimate uses Ridders; generic code checks
`e.abs_error.has_value()` first.

### Implicit conversion to T

Confirmed on clang-cl, `cl`, GCC 16 and Emscripten (spike): with T = double, `double d = *r`, `r == 6.0`,
`r.transform([](double df) { ... })`, arithmetic, passing `*r` to a function taking T, structured bindings and
constant evaluation all work, and `*r == T(6)` is not ambiguous with the defaulted `==`. To document:
- `auto d = *r` deduces `Estimate`, not T.
- Operators that are function templates ignore user conversions: with Boost.Multiprecision, `*r + y`, `*r - y` and
  `abs(*r - y)` do not compile; write `T(*r)` or `r->value`. `std::format` likewise needs `r->value`.
- The conversion chains on: `if (*r)` tests the value against 0, and `int n = *r` narrows silently. Restricting the
  conversion to T exactly (a deleted conversion template) was tried and rejected: Clang then finds `r == 6.0` and
  `*r * 2.0` ambiguous.
- With an expression-template T, a `transform` lambda must return T explicitly (`-> T`), or it returns an expression
  that refers to its destroyed parameter; this holds for any `std::expected<T>`, not only for `Estimate`.
- An empty `abs_error` compares less than any value: test `has_value()` before `*e.abs_error < tol`.

## Algorithm (Ridders)

From Ridders (1982), Richardson extrapolation and its a posteriori error estimate (Burden et al., *Numerical
Analysis*; Sidi, *Practical Extrapolation Methods*, 2003).

1. **Levels.** Level i (counted from the last restart) is the central difference
   `D_i = (f(x + h_i) - f(x - h_i)) / ((x + h_i) - (x - h_i))` with `h_i = h_start / 2^i`, with its rounding floor
   `r_i = (eps |f(x + h_i)| + eps |f(x - h_i)|) / ((x + h_i) - (x - h_i))` (eps applied term by term, so the sum cannot
   overflow, and saturated at max where the division would: a large fn over a tiny step). 2 evaluations per level.
2. **Tableau**, in Burden's difference form: `A[i][0] = D_i` and
   `A[i][j] = A[i][j-1] + (A[i][j-1] - A[i-1][j-1]) / (4^j - 1)` for 1 <= j <= i. Column j removes the h^(2j) term of the
   central difference's even expansion. The difference form cannot overflow below max (the product form did, for
   float at j = 7 with f' near 2^115) and is exact when the two parents are equal. Only the previous row is kept.
3. **Error estimate** of an extrapolated entry (j >= 1; `A[i][0]` has none): `max(|A[i][j] - A[i-1][j-1]|, r_i)`, the
   change from its diagonal parent (the textbook successive-approximation estimate), never less than the rounding
   floor. The best entry is the one with the smallest estimate since the last restart.
4. **Stopping.** The tableau stops
   - when an entry equals its diagonal parent exactly: it is then the best, with the floor as its estimate (x^3 at 2
     stops after 6 evaluations; a derivative that is exactly 0 after 4 to 6);
   - after two consecutive levels with no entry of a smaller estimate;
   - when h_i reaches the resolution of x (`x + h_i == x`, `x - h_i == x` or `h_i == 0`, checked before evaluating);
   - after 32 attempted levels in total, restarts included (at most 64 evaluations).
5. **Restart below a failing step.** A level fails if fn is not finite at x +- h_i, if `D_i` or a tableau entry is not
   finite (finite values whose difference overflows), or if h_i fails the overflow check for x near max. The tableau,
   the best entry and the stall counter are then discarded, and it restarts at `h_start = h_i / 2`: fn misbehaves
   within h_i of x, so every earlier level's interval contained the bad point.
6. **Result.** The best entry, with `abs_error = max(its estimate, |best - A[last][last]|)`: it must also agree with
   the most recent diagonal entry, which catches a best entry that only looked good. A result needs at least two
   successful levels since the last restart.
7. **Restart on an untrusted result.** If `4 abs_error >= |value|` (at most two correct bits) and
   `abs_error > 64 r_best` (more than rounding can explain), fn misbehaves within h_start without producing a non-finite
   value (a pole between the samples, aliasing): restart at `h_start / 2`, within the same cap. When the cap runs out,
   the last result is returned with its large estimate; so is a result whose tableau ended at the resolution of x,
   where a restart could not go lower. A derivative that is 0, whose estimate is near the floor, is not affected. The
   test is written with divisions (`abs_error >= |value| / 4`, `abs_error / 64 > r_best`), which cannot overflow when
   the floor saturates.
8. **Inputs.** x not finite: `non_finite`, before any evaluation. If no level succeeds within the cap: `non_finite`.

Arithmetic: the factors 4^j are built by repeated multiplication in T (exact powers of two in binary; exact up to
about 4^83 in the 50-digit decimal type; `4^j - 1` rounds in float for j >= 13, by at most eps / 2, which does not
matter); only `+ - * /`; constexpr.

## The constants, and why

| Constant | Value | Why |
|---|---|---|
| Step ratio | 2 | Each column then divides by 4^j - 1. A power of the radix in binary, so the steps h / 2^i are exact. In the prototype as accurate as the ratio 1.4 that is often used. |
| Initial step | `max(|x|, 1) / 32` | s/8 and s/32 were equally accurate and s/32 needed 2 to 4 fewer evaluations; it also tolerates faster variation and closer domain edges. s/2 was fooled by sin(20x); s/128 lost a digit. |
| Stall | two levels without improvement | One is not enough: a Runge-function case then stopped early with an error of 3.7e-5 (worst 2.6e-12 with two levels, over 3,200 random cases). |
| Rounding floor | `r_i`, factor 1 | Without it the estimate was below the actual error in 58% of the random cases (it is the minimum of noisy differences); with it and the final check, in 0.7 to 1.0%. A larger factor only makes the estimate more conservative. |
| Untrusted result | `4 abs_error >= |value|`, `abs_error > 64 r_best` | With the factor 1 instead of 4, 11 of 400 poles 1/(x - a), a in (1/128, 1/32), stayed uncovered (estimate 0.85 to 0.98 of a wrong-signed value); with 4 or 16, none. The factor 64 over the floor: results identical for 16 to 256; no random case and no zero derivative restarted. |
| Level cap | 32 | A guard: the prototype stops after 5 to 11 levels for double and float, 10 to 19 for 50 and 100 digits, and used up to 18 near a domain edge. |

These differ from the Numerical Recipes routine (`dfridr`), whose code and constants we may not use: its ratio is 1.4,
it stops on a growth factor, it uses the product form, its estimate is the larger of two differences without a
rounding floor or a final check, and it neither restarts nor stops on exact agreement nor has a default step.

## Prototype results

Python: numpy `float32` and `float64`, `decimal` for 50 and 100 digits; default initial step; the algorithm above.

| Case | Relative error | Relative estimate | Evaluations |
|---|---|---|---|
| exp at 1 | 9.8e-16 | 6.0e-14 | 12 |
| sin at 0.5 | 8.7e-15 | 3.1e-14 | 12 |
| x^3 at 2 | 0 | 9.5e-15 | 6 |
| Runge function at 0.3 | 5.6e-15 | 4.9e-14 | 16 |
| log at 0.2 | 1.6e-14 | 7.3e-14 | 16 |
| log at 1e-4 (domain edge, restarts) | 1.6e-14 | 2.8e-13 | 36 |
| log at 0.01 (domain edge, restarts) | 5.0e-14 | 1.8e-13 | 24 |
| 1/(x - 1/128) at 0 (pole on a sample, restarts) | 5.3e-15 | 5.7e-14 | 26 |
| 1/(x - 0.01) at 0 (pole between samples, untrusted restart) | 3.8e-15 | 3.9e-14 | 30 |
| sin(64 pi x) at 0 (aliasing, untrusted restart) | 0 | 2.2e-16 | 24 |
| float: exp at 1 | 7.2e-6 | 1.9e-5 | 10 |
| float: sin at 0.5 | 4.4e-6 | 3.3e-5 | 10 |
| float: Runge function at 0.3 | 3.3e-7 | 5.7e-6 | 10 |
| float: log at 0.2 | 1.1e-6 | 9.8e-6 | 12 |
| exp at 1, 50 digits | 5.1e-47 | 4.1e-46 | 20 |
| exp at 1, 100 digits | 2.0e-95 | 1.3e-94 | 30 |

Random cases (nine functions, x in [0.2, 3], 300 each): the estimate covered the actual error in 99.3% (double) and
99.0% (float); the worst actual error was 3.6 (double) and 3.1 (float) times the estimate; the median estimate was
6.9 and 6.8 times the actual error; at most 22 and 18 evaluations. All 400 poles 1/(x - a), a in (1/128, 1/32), are
covered in both types. For comparison, `Central<6>` reaches about 1e-14 for exp at 1 in 6 evaluations, without an
estimate.

## Behaviour to document

- **The estimate is an estimate.** It assumes fn is accurate to about eps; a noisier fn (an iterative solver, say)
  makes it optimistic. It is 0 only when eps |f| underflows or f is 0 at both points.
- **fn must be smooth over h.** Ridders restarts when it notices that it is not (non-finite values, or a result with
  at most two correct bits), but fn can still vary in a way its samples do not show. A large `abs_error` relative to
  `value` means: pass a smaller h.
- **Domain edges.** With the default step Ridders needs the nearest bad point of fn to be at least about
  2^-34 max(|x|, 1) away (the cap of 32 halvings, less the levels a result needs); closer, or exactly at an edge, use
  `derivative<Forward<P>>` or `derivative<Backward<P>>`.
- **Cost.** Ridders: 2 evaluations per level, typically 10 to 16 for double, 20 to 38 for 50 or 100 digits, more
  after restarts, at most 64. A stencil: its number of nodes.

## Breaking changes

All pre-release:
- `derivative` returns `std::expected<Estimate<T>, error>` for every method. With the conversion, `double d = *r`,
  comparisons with a T and `transform` with a T parameter still compile; `auto`, `std::format` and multiprecision
  arithmetic need `->value` or `T(*r)`. The tests, the demo (its `show` uses `std::println`), the README and the
  tutorial change where needed.
- The first template parameter of `derivative` is a method rather than a stencil; every stencil is still accepted.

## Tests

- **Polynomials** (`static_assert` and run time, binary types only: decimal division need not round correctly): with
  x and h chosen per type so that every level is exact (x = 1/2, h = 1/4 for float; x = 2, h = 1/16 for double), x^2
  and x^3 give the exact value with `abs_error` equal to the floor.
- **Accuracy** against closed forms (exp, sin, the Runge function, log) for double and float, and for the restart
  cases: the actual error is within 6 times the estimate (the prototype's worst over 5,400 cases was 3.6) and within
  4 times the prototype's actual error for that case (4 eps where that was 0); the estimate is at most twice the
  prototype's; the number of evaluations at most the prototype's plus 4 (libm implementations round fn's values
  differently, which can move a stall by a level). A linear fn (finite values whose difference first overflows) has
  no truncation error, so its estimate is the floor, eps |f'|, within a factor 4.
- **Multiprecision** (`cpp_bin_float_50` with and without expression templates, `cpp_dec_float_50`): coverage as
  above; the estimate bounded per type from the C++ build's first run (the decimal prototype does not model
  `cpp_dec_float`'s guard digits), with the margin stated.
- **Restarts**: log near 0, a pole on a sample, a pole between samples (non-dyadic a), aliasing; fn never finite gives
  `non_finite` within 64 evaluations; finite values whose difference overflows (x = 1e300, fn(t) = 5e9 (t - 1e300)).
- **Failures**: non-finite x; x near `max`; a tiny h that reaches the resolution of x (ends the tableau, no assert).
- **Cost**: evaluations counted, at most 64.
- **Uniform type**: every stencil returns an `Estimate` with an empty `abs_error`; the Phase 1 and 2 tests keep their
  cases, adapted where the conversion does not apply.
- **Conversion**: the cases of the spike, as `static_assert` for double.
- **Callables**: the user-defined real type, the member-function pointer, and a T without a default constructor.

## Decisions from the reviews and the user

- **Uniform API** (user): `derivative<Ridders>`, one return type for every method, implicitly convertible to T.
- **Restart below a failing step at any level** (user, generalised from level 0): a later failure has the same cause.
- **Restart on an untrusted result, stop on exact agreement** (prototypes 4 and 5): poles between samples (46 of 400
  uncovered before, 0 after), aliasing (now exact after one restart), zero derivatives (up to 64 evaluations before, 4
  to 6 after).
- **Difference form, rounding floor, final check, resolution check before each level, non-finite entries as failures,
  floor that cannot overflow, reset on restart, results only from extrapolated entries** (numerics-reviewer), each
  with its counterexample.
- **Estimate formula**: the single textbook difference rather than the larger of two; identical results.
- **Two-level stall** (simplicity-critic asked for evidence): the one-level stall's 3.7e-5 failure.
- **Closed `detail::method`, one dispatch point per overload, a shared two-node helper, `Estimate` fixed, defaulted
  `==` kept** (simplicity-critic).
- **Not adopted** (usage-sketcher): `relative_error()`, `error_or_inf()` or `within(tol)` members, a
  `has_error_estimate<M>` trait, a formatter, a field with the levels used; none has a current use, and a field would
  break `Estimate`.

## Deferred

- **A tolerance** that stops Ridders early once the estimate is below it. A run-time value cannot travel through the
  tag, and making `Ridders` a template would break `derivative<Ridders>`: a new tag type or a trailing overload
  argument, when there is a use.
- **Variants of Ridders** (one-sided, another ratio or cap, another base stencil): a new class template, with
  `Ridders` becoming an alias of its default, is non-breaking; templating `Ridders` itself would be breaking.
- **Error estimates for stencils** from extra evaluations: only as a new tag. Filling the optional for the existing
  stencils would change their cost without the caller asking (breaking in behaviour).
- **More `error` values** (domain, overflow): breaking for an exhaustive `switch`; the documentation reserves them.
- **Second derivatives** (Phase 4): Ridders extends to them with the second central difference as the base; keeping
  `method` closed and M first leaves the shape open.
- **Moving `Estimate<T>` to the shared core** when a second module returns one, with an alias in `nxx::deriv`.

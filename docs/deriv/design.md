# Design: derivatives (`nxx::deriv`)

Decisions in force for `nxx::deriv`. Each one says why; change it here before changing the code. What holds for every
module (the library, real types, the general test rules) is in `docs/design.md`; the phases are in `docs/plan.md`.

## Derivatives (`deriv.hpp`)

- One entry point, `derivative<M, N>(fn, x[, h])`, for every method M (the stencils, `Ridders`, `PHASE3_API.md`, and
  its one-sided forms `ForwardRidders` and `BackwardRidders`, `PHASE5_API.md`) and every derivative order N from 1 to
  4 (`PHASE4_API.md`); N defaults to 1, so `derivative<M>` is unchanged. The method concept `detail::method` is closed
  (`stencil<M> || extrapolation<M>`), not a trait users could specialise. Each public overload forwards to one detail
  function, which dispatches once with `if constexpr`.
- The order is an integer or the scoped enum `order` (`order::second`, user decision), so the parameter is `auto N`:
  a scoped enum does not convert to `int`. `detail::valid_order<N>` accepts any integer type except `bool`, and
  `order`, from 1 to 4; it is converted to `int` once, below the public overloads. The cap is 4 (user decision): the
  textbook tables stop there, and the accuracy is already about 1e-4 for f'''' with `Central<2>` in double.
- A derivative above the first in a type whose epsilon is above 1e-7 (float; tested with integers,
  `radix^(digits-1) < 10^7`) warns at compile time, for every method (user decisions: f'''' in Phase 4, f'' and f'''
  added in Phase 5, since float is not a good type for derivatives above the first). The median relative errors in
  float are 4e-4 to 2e-2 for f'' to f'''' with `Central<2>`, and 1e-5 to 1e-3 with `Ridders`. The mechanism is a twin
  of each public overload, additionally constrained by `detail::noisy` and marked `[[deprecated]]`, which subsumption
  selects; a
  deprecated helper inside the body was rejected, because its warning is reported inside the header, so a pragma at
  the call cannot silence it, and GCC and MSVC suppress it for a system-header include. `NUMERIXX_NO_PRECISION_WARNING`
  turns it off; it must be defined alike in every translation unit.
- Every method returns `std::expected<Estimate<T>, error>`, `Estimate<T> = {value, std::optional<T> abs_error}`: one
  return type for every method (user decision). The stencils leave `abs_error` empty, because a stencil can only bound
  its truncation error with extra evaluations, and the rounding part alone would mislead. `Estimate` converts
  implicitly to T (user decision); it stays an aggregate. Restricting the conversion to T exactly (a deleted
  conversion template) would stop `if (*r)` and `int n = *r`, but makes `r == 6.0` ambiguous on Clang, so it was not
  done. The defaulted `==` compares both members; without it, `==` would compare only the values, through the
  conversion.
- Stencils are the class templates `Central<P>` (even P, 2 to 10), `Forward<P>` and `Backward<P>` (1 to 6), the first
  template argument of `derivative`, defaulting to `Central<2>` (`PHASE2_API.md`). P is the accuracy order, the O(h^P)
  of the error, not the number of points: for the N-th derivative `Central<P>` has P + N - 1 nodes (without 0 for odd
  N), the one-sided rules P + N. They are empty tags; their tables are internal (`detail::make_table<N>`), because no
  current use needs user-written stencils. The caps are the same for every N and are design decisions, not measured
  limits: by the model's exponent alone, `Central<12>` would cut the error of f' by a factor of about 1.7 in double and
  about 5 in a 50-digit type, and the one-sided rounding factor `S = sum |w| / den`, 27.7 at P = 6 for f' (4019 for f''''), grows about as
  fast as the order pays (`Forward<7>` would lower the best error for f' by about 14 %, `Forward<8>` not at all).
- The weights of every rule come from one consteval generator, the Lagrange formula
  `w_k = N! [t^N] prod_{j != k} (t - n_j) / (n_k - n_j)` (Burden et al.), in exact long long arithmetic over one common
  denominator; each family only supplies its nodes. It replaced the Phase 2 closed forms, which it reproduces bit for
  bit (checked once when it was introduced): the one-sided rules of higher orders have no short closed form, so one
  generator builds every family and order. Every value must fit 24 bits, exact in float (the largest numerator is
  6,222,216), or constant evaluation fails.
- The one error so far is `non_finite`: x is not finite, the step is not positive and finite, fn returned an
  infinite or NaN value, a difference, term, sum or quotient would overflow, a point `x + k h` would overflow, or h
  is below the resolution of x. Its documentation reserves more values, and says that a listed cause may move to a
  new one. `std::formatter<error>` writes the enumerator's name; it is in `deriv.hpp`, which therefore includes
  `<format>` (measured on a minimal file: 0.22 to 0.56 s on GCC 16, 0.32 to 0.61 s on MSVC; free where the caller
  prints anyway). Splitting `non_finite` by cause is deferred until a use needs to tell them apart. A step that is not
  positive and finite was first a precondition (`assert`); it is now validated, because
  `derivative` is the public boundary and a step is often computed from data: under `NDEBUG` the assert vanished,
  and then `Forward<1>` with a negative step returned a backward difference (`Backward<1>` a forward one), a step
  of -inf made the two-node rules call fn at infinite points, and in a constant expression a bad step was a
  compile error, not a value. One comparison per call.
- Failures are detected before the arithmetic that would produce them: with `K = max |node|` and `Q = max(2K, den)`,
  the input check `h <= max / 2Q` and `|x| <= max - Q h` keeps every point at least about `K h` below max and the
  divisor finite (for `Central<2>`, Q = 2: `h <= max/4`, `|x| <= max - 2h`). The margin matters: with Q = K, rounding
  could tie a point up to infinity one ulp below max. A zero distance is rejected (with more nodes, `x + h == x` or
  `x - h == x`: at a power of the radix the spacing differs on the two sides), and each of fn's values is checked
  before it is combined. So fn is only called at finite points, and failures are error values in constant
  expressions too (GCC rejects floating-point overflow and NaN in constant evaluation, Clang NaN and division by zero).
  Differences, terms and partial sums are checked before they are formed (`sum_fits`, which decides a + b by
  a / radix + b / radix, exact; and |f_k| <= max / |numerator| before each term), so neither the NaN of inf - inf nor
  an overflow before the final quotient is ever formed; the second difference also checks its one-sided quotients,
  conservatively (|a| / max < b / 2). A term c f_k is decided exactly (`product_fits`: below fl(max / |c|) it fits,
  above it does not, and the tie by the scaled product, exact like `sum_fits`). The divisions by `den * h` and h are
  ordered so that none underflows spuriously (h first when `den * h > 1 > h`). The rounding of a value of fn is
  `eps |f|`, but at least `denorm_min` for a nonzero subnormal value, where rounding is absolute. The exception is an overflow of the final quotient (or of the repeated divisions
  by h for N >= 2), which is checked afterwards: guarding it exactly needs overflow tests accurate to the last unit, so
  in constant evaluation on GCC it is a compile error rather than an error value. Overflow is detected, not avoided:
  halving the operands would make it impossible, but halving is not exact in decimal types and a step near max has no
  numerical meaning.
- `fn` must be invocable with `T const &` and return something convertible to `T`; the result is converted with
  `static_cast<T>` before any arithmetic, so expression templates do not escape. Every call is
  `std::invoke(fn, point)` with `point` a `T const` lvalue, exactly what the constraint promises: member-function
  pointers work, and so do callables that reject rvalues.
- The default-step overload rejects a non-finite x before computing the step, so that it is an error value in
  constant expressions too: the step would be NaN or infinite, and Clang rejects arithmetic that produces a NaN at
  compile time.
- The step parameter is `std::type_identity_t<T>`, so `T` is deduced from `x` alone and a `double` step converts to a
  multiprecision `T`.
- Two-node rules (`Central<2>`, `Forward<1>`, `Backward<1>`, f' only) divide by the distance between the points
  actually evaluated, `x_hi - x_lo`, rather than by `(hi - lo) h`, so `Central<2>` is the Phase 1 central difference
  bit for bit. `Central<2>` for f'' likewise uses the actual offsets, `((f+ - f0) / h+ - (f0 - f-) / h-) /
  ((h+ + h-) / 2)`, exact for quadratics whatever the rounding of the points (the offsets are exact for h <= |x|, by
  Sterbenz's lemma); with the nominal step the Ridders prototype's worst case for float was 8.7 against 1.4 times the
  estimate. Rules with more nodes divide the weighted sum once by `den * h`, then N - 1 more times by h, never by h^N
  (h^4 underflows below 1.2e-77): dividing by the outer-node distance corrects only the outer-node rounding, which for
  more than two weights enlarges the error (4x for `Central<4>`); the divisions keep polynomial results exact for
  dyadic x and h, and cannot overflow spuriously. Neighbouring points must differ (at h near 0.6 ulp(x), x + h and
  x + 2h coincide), and a rule whose integer weights T cannot hold exactly (radix^digits not above the largest value)
  is rejected by the constraint `detail::exact_weights` of `derivative`, so that a requires-expression sees it; an
  extrapolation method is checked against the table of its base stencil (`detail::base<M>`).
- The default step is `eps^(1/(p+N)) * max(|x|, 1)` for a rule of order p and the N-th derivative. It balances
  truncation error (O(h^p)) against rounding error (O(eps / h^N)) when fn varies on a scale of `max(|x|, 1)`, giving an
  error of about `eps^(p/(p+N))`. The constant `(N S / (p C))^(1/(p+N))` is left out. For f' it lies between 1.4 and 2.0
  and costs at most a factor of about 1.6 (1.63 for `Central<10>`); for higher orders it costs more by the model (3.5
  for `Central<2>` and f'', 11 for f''''), but on test functions that vary faster than `max(|x|, 1)` half the default
  step was better for every order, f' included, so one convention is kept for all N. Callers pass a step when fn varies
  faster; the sensitivity grows with N.
- `default_step<M, N>(x)` returns the step that `derivative<M, N>(fn, x)` uses (for an extrapolation method, the
  initial step), so a
  caller can repeat a rule at 2h and h/2 or apply it by hand without rebuilding the formula. The overload without a
  step calls it, so the two calls give the same result for every x: a non-finite x gives a non-finite step, and both
  return `non_finite`.
- The two-point difference is one helper, `detail::difference`, returning the value and the values and spacing for the
  rounding floor `(eps |f_lo| + eps |f_hi|) / (x_hi - x_lo)`; the two-node stencils and every f' Ridders level use it,
  as `detail::second_difference` serves `Central<2>` for f'' and the f'' Ridders levels. Their arithmetic is split from
  the evaluation (`difference_of`, `second_difference_of`), so that the forward levels, which have fn(x) already, use
  the same checks.

## Ridders (`derivative<Ridders>`, `ForwardRidders`, `BackwardRidders`)

The specification, with the evidence for each rule and constant, is `PHASE3_API.md` ("Algorithm" and "The
constants"); the code is implemented from Ridders (1982) and Richardson extrapolation, not from Numerical Recipes.
In short (for every order N, `PHASE4_API.md`; the one-sided forms and the restart rule, `PHASE5_API.md`):
- The base is the `Central<2>` rule of order N (an even error expansion for every N): the difference and the second
  difference above for N = 1 and 2, the nominal-step rules for N = 3 and 4, whose floor is
  `eps sum |w_k| |f_k| / (den h^N)` plus the rounding of the points, `eps/2 (|x| + 2h) |f'| sum |w_k| / (den h^N)`
  with f' from the values at x +- h, each saturated at every division (without the second part, f''' of exp(10x) on
  [1, 3] was covered in 96.3% of the prototype's cases, with it 99.0%). fn(x), part of the base for even N, is evaluated
  once and kept across restarts; if it is not finite, the result is `non_finite` after one evaluation.
- Rules at `h_start / 2^i`, in a tableau in Burden's difference form (it cannot overflow below max and
  is exact for equal parents), keeping only the previous row; the rows are filled from `T(0)`, so T need not be
  default-constructible.
- Each extrapolated entry's estimate is its change from its diagonal parent, floored at the rounding floor; the result
  is the best entry, with `abs_error = max(its estimate, |best - latest diagonal|)`.
- Stops on exact agreement with the diagonal parent, after two levels without improvement, at the resolution of x
  (checked after the overflow check, with q = 2 for N <= 2 and 4 for N >= 3, so x +- h and x +- 2h are only formed when
  finite; for N >= 3 the inner and outer points must also differ), or after 32 attempted levels (at most 64, 65, 128
  and 129 evaluations for N = 1 to 4).
- Restarts below the failing step when a level fails (fn not finite, a non-finite entry, the overflow check), and,
  when the result is untrusted (`4 abs_error >= |value|` and `abs_error > 64 floor`), below the smallest step the
  tableau tried; the restart discards the tableau, the best entry and the stall counter. Until Phase 5 the untrusted
  restart started at half the previous initial step: it re-evaluated the same large steps, and the levels ran out
  before the step reached fn's scale, so poles 1e-3 and 1e-4 from x (finite on both sides) left f'' to f'''' off by
  1 to 11 times with an estimate that did not cover it. No step is skipped (the steps are still `h0 / 2^k`); a case
  the slower descent handled can lose up to 8 times in accuracy, and restarts overlapping one or two levels did worse
  on every measure.
- The default initial step is `max(|x|, 1) / 32` for f' and `max(|x|, 1) / 4` for N >= 2: the base's rounding grows
  as 1/h^N, so a larger start pays (for f'' 1.6 times as accurate as s/8 and 7 times as s/32).
- The nearest point where fn is not finite must be at least about `2^-34 max(|x|, 1)` from x for f' and `2^-31` for
  f'' to f'''' (the 32 levels, less those a result needs); with the restart rule above this holds for poles too.

One driver, `detail::ridders<S, N>`, serves all three methods, on the base stencil S = `detail::base<M>`: `Central<2>`
for `Ridders`, `Forward<1>` for the one-sided methods (user decision, from the study in `PHASE5_API.md`). The reach
for the overflow check (`detail::reach`, shared with the stencils), whether fn(x) is needed, and the points that must
differ follow from S's table. Two things are chosen by S:
- The column factor: 4^j for the even central expansion, 2^j for the forward rule, whose error has every power of h.
- The floor of an entry: centrally the newest level's, as tuned in Phase 3; one-sided it is carried through the
  tableau, `(r floor(i, j-1) + floor(i-1, j-1)) / (r - 1)`, saturating (`detail::propagated_floor`), exactly
  `sum |c_k| floor_k` since the coefficients alternate. Columns divided by 2^j - 1 amplify the newest level's rounding
  about 3.5 times rather than 1.4: with the level's floor the prototype's forward f' was covered in 97.5% of the cases
  (worst 8.9 times the estimate), with the propagated one in 99.4% (worst 2.6). For `Ridders` it would change nearly
  every estimate for a small gain, so it was not adopted there.

The one-sided methods:
- The forward levels follow the central conventions: the actual distance for f', the second difference on the actual
  offsets of x, x + h, x + 2h for f'', and the nominal step with the node-rounding term, `eps/2 (|x| + N h) |f'|`, for
  f''' and f''''. fn(x) is evaluated once per call. Cost: at most 33, 65, 97 and 129 evaluations.
- `BackwardRidders` runs the forward driver on `fn(-t)` at `-x` and negates odd orders: negation is exact, so fn
  receives exactly `x - k h`, and the mirror holds by construction. A `Backward<1>` table would sum from `x - N h` and
  round differently.
- The initial step is `max(|x|, 1) / 32` for every N (`detail::one_sided_divisor`): the rules reach N h, and s/32 was
  the more accurate of the two starts without gross errors for f' to f''', and the best for f''''.
- Their reach from the default step: fn must vary on a scale of at least about `2^-30 max(|x|, 1)` (sqrt or log at
  1e-9 from 0). Closer, the levels run out and the result is wrong, but flagged by an estimate of at least |value| / 4.
- They are less accurate than `Ridders` (medians 2e-13, 1e-10, 4e-8 and 2e-6 for f' to f'''' in double, against 1e-14
  to 7e-10), and for fn defined on one side only: where both sides exist, `Ridders` is better, also near a NaN edge,
  whose NaN guides its restarts.

## Tests

- Ridders: the actual error must be within 6 times the estimate, and the estimate within twice the independent
  Python prototype's value for the case (for the multiprecision types, which the prototype does not model, twice the
  first C++ run's: a regression baseline, as the test rule in CLAUDE.md allows, labelled with its commit). The Phase 5
  prototype (`.scratch/onesided/ridders5.py`) reproduces the C++ results of both bases for double bit for bit
  (28,000 cases). The one-sided methods' multiprecision tests check coverage and a trusted result, without a
  baseline. Polynomials with dyadic x and h are exact, with the rounding floor of the predicted level as the
  estimate (binary types only).
- Stencil tables are checked by `static_assert` against the moment conditions of their order and a reduced
  denominator, with a moment function in the test file, independent of the generator (a Taylor check against the
  interpolation formula), and against published rows (Fornberg 1988, Chapra ch. 23). Run-time tolerances are `(C + 2S) eps^(p/(p+1))`:
  truncation plus rounding of fn's values and of the points, a bound that the measured errors stay 4 to 16 times
  below for double; `(C + 2S) eps^(p/(p+N))` for the N-th derivative.
- The low-precision warning is tested by CTest cases that build a float f'' with warnings as errors, one per overload:
  without silencing each must fail with the warning's text, and silenced at the call or by the macro it must build.
  Tests that make deliberate float calls above f' silence them at the call (`tests/silence.hpp`), never inside a
  template, where GCC and MSVC report the warning.

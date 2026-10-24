# Design

Decisions in force. Each one says why; change it here before changing the code.

## Library

- Header-only: the CMake target `numerixx::numerixx` is an INTERFACE library with no dependencies. The tests fetch
  doctest and Boost.Multiprecision with CPM.
- One header per module (`deriv.hpp`), plus `math.hpp` for what the modules share. A shared core grows only when a
  second module needs something.

## Real types (`math.hpp`)

- `nxx::real<T>` accepts a type that specialises `std::numeric_limits` as non-integer, is constructible from `int`, is
  copyable and totally ordered, and has `+ - * /` and unary minus whose results convert to `T`. This admits the
  built-in floating types and Boost.Multiprecision types, with or without expression templates. It also requires the
  `numeric_limits` members the code reads: `radix` (at least 2) and `digits` (at least 1) as constant expressions,
  `epsilon()` and `max()` convertible to `T`, and `denorm_min()` from which `T` can be constructed, so that a type
  lacking one fails at the concept instead of inside the header. `denorm_min()` need only construct `T` because the
  code reads it as `T(denorm_min())`; a type that inherits `numeric_limits<double>` and constructs from `double`
  only explicitly is accepted.
- The type must also be bounded and inexact (`numeric_limits` `is_bounded`, `!is_exact`): `isfinite` needs a finite
  `max()` and `root_epsilon` a positive `epsilon()`. This rejects rationals such as `cpp_rational`, whose epsilon is 0
  (the scaling loop never ended). The concept is not split into arithmetic and default-step parts: exact types have no
  current use, and even the explicit-step derivative needs `isfinite`.
- The helpers in `nxx::math` (`abs`, `isfinite`, `root_epsilon`) use only comparisons and `+ - * /`. Why:
  - `<cmath>` is not constexpr for most functions in C++23 and does not accept user types.
  - ccmath (v0.3.0) failed on em++ and on llvm-mingw's libc++, and supports only the built-in floating types.
  - The same operations give the same results at compile time and run time, on every compiler.
- `isfinite(x)` is `x == x && |x| <= max`: NaN is the only value unequal to itself.
- `root_epsilon<N, T>()` scales epsilon into (radix^-N, 1] by powers of radix^N, which is exact, then runs Newton's
  iteration for `y^N = a` from 1 down to the root and stops when the iterate stops decreasing. One function serves
  every default step (N = p + n for a stencil of order p and the n-th derivative); it replaced `cbrt_epsilon`. Any N >= 2 works: each
  scaling step multiplies by the radix one factor at a time and stops once the value passes 1, so radix^N, which
  overflows for large N (2^128 for float), is never formed.

## Derivatives (`deriv.hpp`)

- One entry point, `derivative<M, N>(fn, x[, h])`, for every method M (the stencils and `Ridders`, `PHASE3_API.md`)
  and every derivative order N from 1 to 4 (`PHASE4_API.md`); N defaults to 1, so `derivative<M>` is unchanged. The
  method concept `detail::method` is closed (`stencil<M> || same_as<M, Ridders>`), not a trait users could
  specialise. Each public overload forwards to one detail function, which dispatches once with `if constexpr`.
- The order is an integer or the scoped enum `order` (`order::second`, user decision), so the parameter is `auto N`:
  a scoped enum does not convert to `int`. `detail::valid_order<N>` accepts any integer type except `bool`, and
  `order`, from 1 to 4; it is converted to `int` once, below the public overloads. The cap is 4 (user decision): the
  textbook tables stop there, and the accuracy is already about 1e-4 for f'''' with `Central<2>` in double.
- A 4th derivative in a type whose epsilon is above 1e-7 (float; tested with integers, `radix^(digits-1) < 10^7`)
  warns at compile time (user decision): its result is mostly rounding noise. The mechanism is a twin of each public
  overload, additionally constrained by `detail::noisy` and marked `[[deprecated]]`, which subsumption selects; a
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
  limits: by the model's exponent alone, `Central<12>` would gain about 0.2 digit for f' in double and 0.7 in a 50-digit
  type, and the one-sided rounding factor `S = sum |w| / den`, 27.7 at P = 6 for f' (4019 for f''''), grows about as
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
- Numerixx never throws and never catches: an exception thrown by fn or by T's operations passes through unchanged.
  Catching would break builds without exceptions, such as Emscripten's default.
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
  is rejected by the constraint `detail::exact_weights` of `derivative`, so that a requires-expression sees it;
  Ridders is checked against the `Central<2>` table it uses.
- The default step is `eps^(1/(p+N)) * max(|x|, 1)` for a rule of order p and the N-th derivative. It balances
  truncation error (O(h^p)) against rounding error (O(eps / h^N)) when fn varies on a scale of `max(|x|, 1)`, giving an
  error of about `eps^(p/(p+N))`. The constant `(N S / (p C))^(1/(p+N))` is left out. For f' it lies between 1.4 and 2.0
  and costs at most a factor of about 1.6 (1.63 for `Central<10>`); for higher orders it costs more by the model (3.5
  for `Central<2>` and f'', 11 for f''''), but on test functions that vary faster than `max(|x|, 1)` half the default
  step was better for every order, f' included, so one convention is kept for all N. Callers pass a step when fn varies
  faster; the sensitivity grows with N.
- `default_step<M, N>(x)` returns the step that `derivative<M, N>(fn, x)` uses (for Ridders, the initial step), so a
  caller can repeat a rule at 2h and h/2 or apply it by hand without rebuilding the formula. The overload without a
  step calls it, so the two calls give the same result for every x: a non-finite x gives a non-finite step, and both
  return `non_finite`.
- The two-point difference is one helper, `detail::difference`, returning the value and the values and spacing for the
  rounding floor `(eps |f_lo| + eps |f_hi|) / (x_hi - x_lo)`; the two-node stencils and every f' Ridders level use it,
  as `detail::second_difference` serves `Central<2>` for f'' and the f'' Ridders levels.

## Ridders (`derivative<Ridders>`)

The specification, with the evidence for each rule and constant, is `PHASE3_API.md` ("Algorithm" and "The
constants"); the code is implemented from Ridders (1982) and Richardson extrapolation, not from Numerical Recipes.
In short (for every order N, `PHASE4_API.md`):
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
- Restarts below the failing step when a level fails (fn not finite, a non-finite entry, the overflow check), and
  below the starting step when the result is untrusted (`4 abs_error >= |value|` and `abs_error > 64 floor`); the
  restart discards the tableau, the best entry and the stall counter.
- The default initial step is `max(|x|, 1) / 32` for f' and `max(|x|, 1) / 4` for N >= 2: the base's rounding grows
  as 1/h^N, so a larger start pays (for f'' 1.6 times as accurate as s/8 and 7 times as s/32).

## Tests

- doctest. Every constexpr function has `static_assert` tests next to run-time tests of the same cases.
- Ridders: the actual error must be within 6 times the estimate, and the estimate within twice the independent
  Python prototype's value for the case (for the multiprecision types, which the prototype does not model, twice the
  first C++ run's). Polynomials with dyadic x and h are exact, with the rounding floor of the predicted level as the
  estimate (binary types only).
- Stencil tables are checked by `static_assert` against the moment conditions of their order and a reduced
  denominator, with a moment function in the test file, independent of the generator (a Taylor check against the
  interpolation formula), and against published rows (Fornberg 1988, Chapra ch. 23). Run-time tolerances are `(C + 2S) eps^(p/(p+1))`:
  truncation plus rounding of fn's values and of the points, a bound that the measured errors stay 4 to 16 times
  below for double; `(C + 2S) eps^(p/(p+N))` for the N-th derivative.
- The low-precision warning is tested by CTest cases that build a float f'''' with warnings as errors: without
  silencing it must fail with the warning's text, and silenced at the call or by the macro it must build.
- Multiprecision cases use standalone Boost.Multiprecision 1.92 with Boost.Config (header-only): the full Boost CMake
  tarball builds Boost.Container, which fails on em++. Standalone mode has no `cbrt`, so references use `pow`.

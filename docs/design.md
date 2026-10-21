# Design

Decisions in force. Each one says why; change it here before changing the code.

## Library

- Header-only: the CMake target `numerixx::numerixx` is an INTERFACE library with no dependencies. The tests fetch
  doctest and Boost.Multiprecision with CPM.
- One header per module (`deriv.hpp`), plus `math.hpp` for what the modules share. A shared core grows only when a
  second module needs something.

## Real types (`math.hpp`)

- `nxx::real<T>` accepts a type that specialises `std::numeric_limits` as non-integer, is constructible from `int`, is
  totally ordered, and has `+ - * /` and unary minus whose results convert to `T`. This admits the built-in floating
  types and Boost.Multiprecision types, with or without expression templates.
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
  every default step (N = p + 1 for a stencil of order p); it replaced `cbrt_epsilon`. Any N >= 2 works: each
  scaling step multiplies by the radix one factor at a time and stops once the value passes 1, so radix^N, which
  overflows for large N (2^128 for float), is never formed.

## Derivatives (`deriv.hpp`)

- Stencils are the class templates `Central<P>` (even P, 2 to 10), `Forward<P>` and `Backward<P>` (1 to 6), the first
  template argument of `derivative`, defaulting to `Central<2>` (`PHASE2_API.md`). They are empty tags; their tables
  are internal (`detail::make_table`), because no current use needs user-written stencils. The caps are numerical:
  beyond `Central<10>` the gain is under half a digit even for 50-digit types, and the one-sided rounding factor
  `S = sum |w| / den` is 27.7 at P = 6.
- The weights come from closed forms, computed at compile time in exact integer arithmetic and put over one common
  denominator, so they are exact in every real type (the largest numerator is 2100). Fornberg's general algorithm is
  not needed for three fixed families.
- Results are `std::expected<T, nxx::deriv::error>`. The one error so far is `non_finite`: fn returned an infinite or
  NaN value, the sum or quotient overflowed, a point `x + k h` would overflow, or h is below the resolution of x. A
  non-positive step is a precondition violation (`assert`).
- Failures are detected before the arithmetic that would produce them: with `K = max |node|` and `Q = max(2K, den)`,
  the input check `h <= max / 2Q` and `|x| <= max - Q h` keeps every point at least about `K h` below max and the
  divisor finite (for `Central<2>`, Q = 2: `h <= max/4`, `|x| <= max - 2h`). The margin matters: with Q = K, rounding
  could tie a point up to infinity one ulp below max. A zero distance is rejected (with more nodes, `x + h == x` or
  `x - h == x`: at a power of the radix the spacing differs on the two sides), and each of fn's values is checked
  before it is combined. So fn is only called at finite points, and failures are error values in constant
  expressions too (GCC rejects floating-point overflow and NaN in constant evaluation, Clang NaN and division by zero).
  The exception is an overflow of the final difference, sum or quotient, which is checked afterwards: guarding it exactly
  needs overflow tests accurate to the last unit, so in constant evaluation on GCC it is a compile error rather than an
  error value. Overflow is detected, not avoided: halving the operands would make it impossible, but halving is not
  exact in decimal types and a step near max has no numerical meaning.
- Numerixx never throws and never catches: an exception thrown by fn or by T's operations passes through unchanged.
  Catching would break builds without exceptions, such as Emscripten's default.
- `fn` must be invocable with `T const &` and return something convertible to `T`; the result is converted with
  `static_cast<T>` before any arithmetic, so expression templates do not escape. Every call is
  `std::invoke(fn, point)` with `point` a `T const` lvalue, exactly what the constraint promises: member-function
  pointers work, and so do callables that reject rvalues.
- The default-step overload rejects a non-finite x before computing the step: a NaN x would give a NaN step and trip
  the `h > 0` precondition, which stays an `assert` for caller-supplied steps.
- The step parameter is `std::type_identity_t<T>`, so `T` is deduced from `x` alone and a `double` step converts to a
  multiprecision `T`.
- Two-node rules (`Central<2>`, `Forward<1>`, `Backward<1>`) divide by the distance between the points actually
  evaluated, `x_hi - x_lo`, rather than by `(hi - lo) h`, so `Central<2>` is the Phase 1 central difference bit for
  bit. Rules with more nodes divide the weighted sum once by `den * h`: dividing by the outer-node distance corrects
  only the outer-node rounding, which for more than two weights enlarges the error (4x for `Central<4>`); one division
  keeps polynomial results exact for dyadic x and h.
- The default step is `eps^(1/(p+1)) * max(|x|, 1)` for a rule of order p. It balances truncation error (O(h^p))
  against rounding error (O(eps / h)) when fn varies on a scale of `max(|x|, 1)`, giving an error of about
  `eps^(p/(p+1))`. The order-dependent constant `(S / (p C))^(1/(p+1))` is left out; it lies between 1.4 and 2.0 and
  costs at most a factor 1.64. Callers pass a step when fn varies faster.

## Tests

- doctest. Every constexpr function has `static_assert` tests next to run-time tests of the same cases.
- Stencil tables are checked by `static_assert` against the moment conditions of their order, with a moment function
  in the test file, independent of the closed forms in the library. Run-time tolerances are `(C + 2S) eps^(p/(p+1))`:
  truncation plus rounding of fn's values and of the points, a bound that the measured errors stay 4 to 16 times
  below for double.
- Multiprecision cases use standalone Boost.Multiprecision 1.92 with Boost.Config (header-only): the full Boost CMake
  tarball builds Boost.Container, which fails on em++. Standalone mode has no `cbrt`, so references use `pow`.

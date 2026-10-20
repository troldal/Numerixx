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
- The helpers in `nxx::math` (`abs`, `isfinite`, `cbrt_epsilon`) use only comparisons and `+ - * /`. Why:
  - `<cmath>` is not constexpr for most functions in C++23 and does not accept user types.
  - ccmath (v0.3.0) failed on em++ and on llvm-mingw's libc++, and supports only the built-in floating types.
  - The same operations give the same results at compile time and run time, on every compiler.
- `isfinite(x)` is `x == x && |x| <= max`: NaN is the only value unequal to itself.
- `cbrt_epsilon<T>()` scales epsilon into (radix^-3, 1] by powers of radix^3, which is exact, then runs Newton's
  iteration from 1 down to the root and stops when the iterate stops decreasing.

## Derivatives (`deriv.hpp`)

- Results are `std::expected<T, nxx::deriv::error>`. The one error so far is `non_finite`: fn returned an infinite or
  NaN value, or the difference overflowed. A non-positive step is a precondition violation (`assert`).
- `fn` must be invocable with `T const &` and return something convertible to `T`; the result is converted with
  `static_cast<T>` before any arithmetic, so expression templates do not escape.
- The step parameter is `std::type_identity_t<T>`, so `T` is deduced from `x` alone and a `double` step converts to a
  multiprecision `T`.
- The central difference divides by `(x + h) - (x - h)`, the distance between the points actually evaluated, rather
  than by `2h`.
- The default step is `eps^(1/3) * max(|x|, 1)`. It balances truncation error (O(h^2)) against rounding error
  (O(eps / h)) when fn varies on a scale of `max(|x|, 1)`, giving an error of about `eps^(2/3)` times fn''' on that
  scale. Callers pass a step when fn varies faster.

## Tests

- doctest. Every constexpr function has `static_assert` tests next to run-time tests of the same cases.
- Multiprecision cases use standalone Boost.Multiprecision 1.92 with Boost.Config (header-only): the full Boost CMake
  tarball builds Boost.Container, which fails on em++. Standalone mode has no `cbrt`, so references use `pow`.

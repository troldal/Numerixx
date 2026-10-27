# Design

Decisions in force for every module. Each one says why; change it here before changing the code. Each module's own
decisions are in `docs/<module>/design.md` (`docs/deriv/design.md`).

## Library

- Header-only: the CMake target `numerixx::numerixx` is an INTERFACE library with no dependencies. The tests fetch
  doctest and Boost.Multiprecision with CPM.
- One header per module (`deriv.hpp`, `roots.hpp`), plus `math.hpp` for what the modules share. A shared core grows
  only when a second module needs something.
- Numerixx never throws and never catches: an exception thrown by fn or by T's operations passes through unchanged.
  Catching would break builds without exceptions, such as Emscripten's default.
- fn is constrained by `nxx::detail::callable<F, T>` in `math.hpp` (invocable with a `T const` lvalue, its result
  convertible to T), shared by the modules since `nxx::roots` was the second to need it; spelled with `nxx::` in the
  modules, whose own `detail` namespaces would otherwise hide it.

## Real types (`math.hpp`)

- `nxx::real<T>` accepts a type that specialises `std::numeric_limits` as non-integer, is constructible from `int`, is
  copyable and totally ordered, and has `+ - * /` and unary minus whose results convert to `T`. This admits the
  built-in floating types and Boost.Multiprecision types, with or without expression templates. It also requires the
  `numeric_limits` members the code reads: `radix` (at least 2) and `digits` (at least 1) as constant expressions,
  `epsilon()` and `max()` convertible to `T`, and `denorm_min()` from which `T` can be constructed, so that a type
  lacking one fails at the concept instead of inside the header. `denorm_min()` need only construct `T` because the
  code reads it as `T(denorm_min())`; a type that inherits `numeric_limits<double>` and constructs from `double`
  only explicitly is accepted.
- `nxx::real` first requires an object type (`std::is_object_v<T>`): overload resolution may try a function type
  as a deduced T (a derivative passed as a free function to `nxx::roots::solve`), and the standard library's
  `numeric_limits` of a function type is ill-formed, a hard error rather than a failed constraint.
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

## Tests

- doctest, one executable (`numerixx_tests`); a module's test cases carry its prefixes (the table in `CLAUDE.md`).
  Every constexpr function has `static_assert` tests next to run-time tests of the same cases.
- Multiprecision cases use standalone Boost.Multiprecision 1.92 with Boost.Config (header-only): the full Boost CMake
  tarball builds Boost.Container, which fails on em++. Standalone mode has no `cbrt`, so references use `pow`.

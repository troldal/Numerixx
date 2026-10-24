# How `nxx::deriv` works: a guided tour

This tutorial walks through the current implementation of numerical derivatives in Numerixx 2: two headers,
`include/numerixx/math.hpp` and `include/numerixx/deriv.hpp`. It starts with how to call the library, then explains the
mathematics behind the central difference and its default step, goes through the code line by line, including the C++
details that let the same code run on `double`, on Boost.Multiprecision types and at compile time, and ends with the
higher-order and one-sided stencils (section 8), Ridders' extrapolation (section 9) and derivatives of order 2 to 4
(section 10).

## 1. Using it

```cpp
#include <numerixx/deriv.hpp>

auto const f = [](double x) { return std::exp(x); };

auto const d = nxx::deriv::derivative(f, 1.0);         // default step
auto const e = nxx::deriv::derivative(f, 1.0, 1e-4);   // explicit step h

if (d) std::println("{}", d->value);                   // 2.7182818284569006 (exact: 2.7182818284590451)

using nxx::deriv::Central, nxx::deriv::Forward, nxx::deriv::Ridders;
auto const g = nxx::deriv::derivative<Central<4>>(f, 1.0);   // higher order: 4 evaluations
auto const k = nxx::deriv::derivative<Forward<2>>(f, 1.0);   // one-sided: never evaluates f below x
auto const r = nxx::deriv::derivative<Ridders>(f, 1.0);      // adapts its step; r->abs_error estimates the error

using nxx::deriv::order;
auto const d2 = nxx::deriv::derivative<Central<2>, 2>(f, 1.0);          // f'': (f(x+h) - 2f(x) + f(x-h)) / h^2
auto const d4 = nxx::deriv::derivative<Ridders, order::fourth>(f, 1.0); // f'''', with an error estimate
```

There are two overloads, with and without a step. The first optional template argument chooses the method: a stencil
(section 8) or Ridders (section 9); by default it is `Central<2>`, the central difference of sections 2 to 7. The
second chooses the derivative order, 1 (the default) to 4, as an integer or an `order` (section 10). Every
method returns `std::expected<Estimate<T>, nxx::deriv::error>`: either an `Estimate`, or `error::non_finite` if no
finite estimate exists (section 5 lists the cases). An `Estimate` holds `value` and `abs_error`, an
`std::optional<T>` that only Ridders fills in, and converts implicitly to T, so `double slope = *d;` works. Numerixx
itself never throws. The demo program (`demo/deriv_demo.cpp`) runs every example in this tutorial.

## 2. The formula: the central difference

### Where it comes from

Expand f around x by Taylor's theorem, one step h to each side:

```
f(x + h) = f(x) + h f'(x) + h^2/2 f''(x) + h^3/6 f'''(x) + ...
f(x - h) = f(x) - h f'(x) + h^2/2 f''(x) - h^3/6 f'''(x) + ...
```

Subtract the second line from the first. The even terms cancel, and dividing by 2h leaves

```
(f(x + h) - f(x - h)) / 2h = f'(x) + h^2/6 f'''(x) + O(h^4)
```

So the *central difference* on the left approximates f'(x) with a *truncation error* of about `h^2/6 * f'''(x)`. It
is "second order": halving h divides the error by four. Two evaluations of f are all it costs.

A quadratic has f''' = 0, so the central difference is exact for it at any h. For x^3, f''' = 6, so the error is
exactly h^2. The tests use both facts.

### Why h cannot just be made tiny

The formula subtracts two nearly equal numbers. Each value of f carries a rounding error of roughly `eps * |f|`,
where `eps` is machine epsilon (2.2e-16 for `double`). The subtraction keeps those errors, and the division by 2h
magnifies them. The *rounding error* of the estimate is about `eps * |f| / h`, and it grows as h shrinks.

The total error is the sum of a term that grows with h and a term that shrinks with it:

```
E(h) ~ h^2/6 * |f'''|  +  eps * |f| / h
```

The demo prints this trade-off for f = sin at x = 1:

```
h = 1e-1   error 9.0e-04      truncation: h^2/6 * cos(1) = 9.0e-04, as predicted
h = 1e-3   error 9.0e-08
h = 1e-5   error 1.2e-11      near the minimum
h = 1e-8   error 2.9e-09      rounding takes over
h = 1e-12  error 3.0e-05
```

### The best step

Setting the derivative of E(h) to zero gives `h = (3 eps |f| / |f'''|)^(1/3)`, and an error at that step of order
`eps^(2/3)`. In other words, the central difference can give about two thirds of the type's digits: roughly 10 to
11 of the 16 digits of `double`.

The library does not know f''' or f, so it assumes that f varies on a scale of `max(|x|, 1)`: near x = 0 the step is
absolute, and for large |x| it is relative to x. Dropping the constant 3^(1/3) gives the default step

```
h = eps^(1/3) * max(|x|, 1)          (6.1e-6 for double at x = 1)
```

The demo's "default h" line shows the result: an error of 2.9e-12 for sin at 1, as good as the best entry of the
table. When f varies much faster than the scale `max(|x|, 1)` (for example `exp(100 x)`), the assumption fails and
the caller should pass a step.

## 3. The real-number concept (`math.hpp`)

Everything is templated on a type T that must satisfy `nxx::real`:

```cpp
template<typename T>
concept real = std::numeric_limits<T>::is_specialized && !std::numeric_limits<T>::is_integer
               && std::numeric_limits<T>::is_bounded && !std::numeric_limits<T>::is_exact
               && std::constructible_from<T, int> && std::totally_ordered<T> && requires(T const a, T const b) {
                      { a + b } -> std::convertible_to<T>;
                      ...
                      { -a } -> std::convertible_to<T>;
                  };
```

Read it as a checklist:

- **`numeric_limits` is specialised and not integer.** The algorithm needs `epsilon()` and `max()`, and an integer
  type would make the difference quotient meaningless. `derivative(f, 1)` with an `int` argument does not compile.
- **Bounded and inexact**, like a floating-point type. `isfinite` compares with `max()`, and `root_epsilon` needs a
  positive `epsilon()`. An exact rational type such as Boost's `cpp_rational` has epsilon 0 and no meaningful
  `max()`, so it is rejected at compile time.
- **Constructible from `int`.** Constants are written `T(1)`, `T(2)`, which works for every real type. A `double`
  literal such as `1.0` would silently lose precision in a 50-digit type, so the code never uses one.
- **Totally ordered, with `+ - * /` and unary minus.** That is all the algorithm uses.
- **"Converts to T" rather than "is T".** Boost.Multiprecision with expression templates returns an expression object
  from `a + b`, which converts to T but is not T. The concept accepts that, and section 5 explains how the code copes
  with it.

`float`, `double`, `long double`, `cpp_bin_float_50` (with and without expression templates) and `cpp_dec_float_50`
all qualify; the tests run on each.

## 4. Three small helpers (`nxx::math`)

The derivative needs |x|, a finiteness test and the cube root of epsilon. The standard library has `std::abs`,
`std::isfinite` and `std::cbrt`, but in C++23 most of `<cmath>` is not constexpr, and it does not accept user types.
So Numerixx writes its own, using only comparisons and `+ - * /`. As a bonus, they give the same answer at compile time
and at run time on every compiler.

### `abs`

```cpp
template<real T>
constexpr T abs(T const &x)
{
    if (x < T(0)) return T(-x);
    return x;
}
```

`T(-x)` turns a possible expression template back into a T.

### `isfinite`

```cpp
template<real T>
constexpr bool isfinite(T const &x)
{ return x == x && math::abs(x) <= (std::numeric_limits<T>::max)(); }
```

- `x == x` is false only for NaN, the one value that compares unequal to itself.
- `|x| <= max` is false for an infinity.
- `(std::numeric_limits<T>::max)()` is in parentheses so that a `max` macro (from `<windows.h>`) cannot expand there.

### `root_epsilon`

The default step of the central difference needs `eps^(1/3)`; a stencil of order p needs `eps^(1/(p+1))` (section 8).
`root_epsilon<N, T>()` computes `eps^(1/N)` in two stages. This section follows the case N = 3.

**Stage 1: scale epsilon into (radix^-3, 1].** Epsilon is tiny (2^-52 for `double`, about 1e-49 for a 50-digit type).
The code multiplies it by radix^3 (8 for binary types, 1000 for decimal ones) until one more multiplication would pass
1, and counts the multiplications in `scale`:

```cpp
T a = std::numeric_limits<T>::epsilon();
T scale = T(1);
while (a * radix_n <= T(1)) {      // radix_n = radix^N
    a = a * radix_n;
    scale = scale * radix;
}
```

Multiplying by a power of the radix only changes the exponent, so this is exact. For `double`, 17 multiplications by
8 = 2^3 turn 2^-52 into a = 2^-1 = 0.5, with scale = 2^17. Since `eps = a / scale^3`, the answer is
`cbrt(eps) = cbrt(a) / scale`.

**Stage 2: Newton's iteration for cbrt(a).** To solve `y^3 = a`, Newton's method iterates

```
y <- (2y + a / y^2) / 3          in general, for y^N = a:  y <- ((N-1) y + a / y^(N-1)) / N
```

Starting at y = 1, which is at or above the root because a <= 1, every iterate stays above the root (by the
inequality of arithmetic and geometric means) and decreases. In exact arithmetic it would decrease forever towards the
root; in floating point it reaches the root to within rounding and then stops decreasing. So the loop runs until the
next iterate is not smaller:

```cpp
T y = T(1);
while (true) {
    T y_n1 = T(1);
    for (int i = 1; i < N; ++i) y_n1 = y_n1 * y;          // y^(N-1)
    T const next = (T(N - 1) * y + a / y_n1) / T(N);
    if (!(next < y)) break;
    y = next;
}
return y / scale;
```

Stage 1 is what makes stage 2 fast. With a in (radix^-N, 1], the root lies in (1/radix, 1], close to the start, so
Newton converges in a handful of steps. Started from 1 against eps itself, it would need dozens of slow steps to come
down to 6e-6.

For `double` the result is `cbrt(0.5) / 2^17 = 0.7937 / 131072 = 6.055e-6`, which the tests check against
`std::cbrt(epsilon)` to a relative error of at most 2 eps; the other roots are checked against `pow(epsilon, 1/N)`.

## 5. The derivative (`deriv.hpp`)

`derivative` handles every method and order (sections 8 to 10 cover the others): it forwards to one detail
function, where one `if constexpr` sends Ridders to `detail::ridders` and the stencils to
`detail::stencil_derivative`, whose result is wrapped as `Estimate{d, std::nullopt}`. For the default `Central<2>` and
the first derivative, the compile-time checks and branches reduce to this (shown as one function returning the value;
in the code the two-point part is the helper `detail::difference`, which Ridders shares):

```cpp
template<typename M = Central<2>, auto N = 1, real T, typename F>
requires detail::method<M> && detail::valid_order<N> && detail::callable<F, T>
[[nodiscard]] constexpr std::expected<T, error> derivative_value(F &&fn, T const &x, std::type_identity_t<T> const &h)
{
    assert(h > T(0) && math::isfinite(h));
    // With |x| + 2h <= max, the points and their distance are finite, so nothing overflows before fn is called.
    T const largest = (std::numeric_limits<T>::max)();
    if (!(h <= largest / T(4) && math::abs(x) <= largest - T(2) * h)) return std::unexpected(error::non_finite);
    T const xm = x - h;
    T const xp = x + h;
    T const dx = xp - xm;
    if (dx == T(0)) return std::unexpected(error::non_finite); // h is below the resolution of x
    T const fm = static_cast<T>(std::invoke(fn, xm));
    T const fp = static_cast<T>(std::invoke(fn, xp));
    if (!math::isfinite(fm) || !math::isfinite(fp)) return std::unexpected(error::non_finite);
    if (!detail::sum_fits(fp, T(T(0) - fm))) return std::unexpected(error::non_finite); // fp - fm would overflow
    T const d = (fp - fm) / dx;
    if (!math::isfinite(d)) return std::unexpected(error::non_finite); // the quotient overflowed
    return d;
}
```

Each line has a reason.

### The signature

- **`typename M = Central<2>`** and **`detail::method<M>`**: the method is a type, the first template argument, so
  `derivative(f, x)` keeps working while `derivative<Central<4>>(f, x)` picks another rule. The concept accepts only
  `Central`, `Forward` and `Backward` in their allowed orders, and `Ridders`; anything else is rejected by the
  constraint, not by a hard error.
- **`auto N = 1`**: the derivative order, second (section 10). It is `auto` so that both an integer and the scoped
  enum `order` are accepted. Since the second template argument is the order, T can no longer be given explicitly
  as `derivative<M, double>(...)`; it is always deduced from x.
- **`requires std::invocable<F &, T const &>`**: fn must be callable with a `T`. The `F &` is there because the body
  calls `fn(...)` on a named reference, which is an lvalue. The result must convert to T, so fn may return an
  expression template or another type that converts.
- **`std::type_identity_t<T> const &h`**: T is deduced from x only. Without this, `derivative(f, bin50(1), 0.5)` would
  deduce T as `bin50` from x and as `double` from h, and fail. `type_identity_t` makes h a "non-deduced context",
  so the `double` 0.5 simply converts to `bin50`.
- **`[[nodiscard]]`**: ignoring the result of a derivative is always a bug.
- **`constexpr`**: if fn is constexpr, the whole computation can run at compile time (section 6).

### The precondition

`assert(h > T(0))`. A zero or negative step is a programming error, not a run-time condition, so it is checked with
`assert` rather than reported through `std::expected`.

### Dividing by the points actually evaluated

The obvious formula divides by `2 * h`. The code divides by `dx = xp - xm` instead. The reason is rounding: `x + h` is
rounded to the nearest representable number, so the points where f is actually evaluated are generally not exactly
2h apart. At x = 1 with h = 1e-5, each point can be off by up to half an ulp of 1, about 1.1e-16, which is a relative
error of about 1e-11 in the spacing, as large as the error the method is trying to achieve. The difference
`xp - xm` is the true distance between the two evaluated points (and computing it is exact in the usual case where the
points are close), so the quotient is the exact slope of the secant through the two points that were evaluated.

### Expression templates

Boost.Multiprecision can return *expression templates*: `a + b` returns a lightweight object that records the
operation and is evaluated only when assigned to a number. That avoids temporaries, but it means `auto d = fp - fm;`
would store an unevaluated expression holding references, which can dangle. Hence two rules in the code:

- every intermediate is declared `T const`, never `auto`, so it is evaluated on the spot;
- fn's result goes through `static_cast<T>` before any arithmetic, since fn may itself return an expression (for
  example `mp::exp(x)`).

### The error checks

Every failure becomes `error::non_finite`, and each one is checked *before* the arithmetic that would produce it:

1. **The input check**, `h <= max/4` and `|x| <= max - 2h`, guarantees that `x + h`, `x - h` and their distance are
   finite. Without it, `derivative(x / 2, 0.0, max)` evaluated f at the finite points `+-max`, but their distance
   overflowed to infinity and the result was a "successful" `0.5 * max / inf = 0`. Neither side of the check can
   overflow: `2h <= max/2`, and `max - 2h` only subtracts.
2. **A zero distance** means h is below the resolution of x (`x + h` and `x - h` both round to x, for example h = 1e-300
   at x = 1). Dividing by it would give 0/0.
3. **fn's values** are checked before they are subtracted: if fn returned NaN (for example `log` at a negative point)
   or an infinity, `inf - inf` is never formed.
4. **The difference** is checked before it is formed: both values can be finite while their difference overflows
   (`max - (-max)`). `detail::sum_fits(a, b)` decides whether a + b is finite without forming it: dividing by the
   radix is exact, so a + b exceeds max exactly when a / radix + b / radix exceeds max / radix, and that sum cannot
   overflow.
5. **The result** is checked last: the quotient can overflow when the distance is tiny.

Checking first rather than last matters for two reasons. fn is only ever called at finite points. And in a constant
expression, compilers reject arithmetic that produces an infinity or a NaN: GCC rejects floating-point overflow and NaN,
Clang rejects NaN and division by zero. Because the checks come first, those failures are ordinary error values even in
a `static_assert`. The one exception is check 5: an overflow of the final quotient is detected after the fact, so in
constant evaluation on GCC it is a compile error rather than an error value.

### Exceptions

Numerixx never throws and never catches. If fn, or an operation of T, throws, the exception passes through
`derivative` unchanged. Catching it would make the library unusable in builds without exceptions, such as
Emscripten's default.

### The default-step overload

```cpp
template<int Exponent, real T>
constexpr T default_step(T const &x)                     // in nxx::deriv::detail
{ return math::root_epsilon<Exponent, T>() * detail::scale(x); }      // scale(x) = max(|x|, 1)

template<typename M, int N, real T, typename F>          // in nxx::deriv::detail, behind the public overload
constexpr std::expected<Estimate<T>, error> derivative_default(F &fn, T const &x)
{
    if (!math::isfinite(x)) return std::unexpected(error::non_finite);
    if constexpr (std::same_as<M, Ridders>) {
        constexpr int divisor = ridders_divisors[N - 1];  // 32 for f', 4 for f'' to f''''
        return detail::derivative_with_step<M, N>(fn, x, T(detail::scale(x) / T(divisor)));
    } else {
        constexpr int exponent = detail::make_table<N>(M{}).order + N;
        return detail::derivative_with_step<M, N>(fn, x, detail::default_step<exponent>(x));
    }
}
```

For `Central<2>` (order 2) and f' this is section 2's `eps^(1/3) * max(|x|, 1)`; for Ridders the initial step is
`max(|x|, 1) / 32` (section 9). Every internal call is written `detail::...`, qualified: an unqualified call would
also look in the namespaces of T and F (argument-dependent lookup), where a user type could provide its own function
of the same name and silently take over.

## 6. At compile time

Since everything is constexpr and uses only `+ - * /`, a derivative of a constexpr function is a constant expression:

```cpp
constexpr auto cube = [](double x) { return x * x * x; };
static_assert(nxx::deriv::derivative(cube, 2.0, 0.25) == 12.0625);   // 3 * 2^2 + h^2, exactly
```

The comparison is exact: with h = 1/4 every intermediate value is a short binary fraction, and the error of the
central difference for x^3 is exactly h^2 = 1/16 (section 2). `std::exp` and friends are not constexpr, so this works
only for fn built from arithmetic.

## 7. How it is tested

`tests/test_deriv.cpp` and `tests/test_stencils.cpp` follow three rules from the project:

- **Exact cases as `static_assert`.** Quadratics and cubics with power-of-two steps, where the result is known exactly
  (section 6), are checked at compile time; the same cases also run at run time.
- **References from closed forms.** `exp' = exp`, `sin' = cos`, `(x^3)' = 3x^2`, never values computed by Numerixx.
- **Every tolerance explained.** Each check carries a comment with its error estimate, for example "exp varies faster
  than the scale |x| = 20: h = 1.2e-4 and the relative error is h^2 / 6 = 2.4e-9", so a tolerance is never just a number
  that happened to pass.

The same tests run on `float`, `double` and the three Boost.Multiprecision types, on five toolchains (clang-cl, MSVC,
GCC, Clang and Emscripten). For the stencils, a `static_assert` also checks every table against the moment conditions
of section 8, with a moment function written in the test file, independent of the library's formulas.

## 8. Stencils

The central difference is the simplest of a family of *stencils* (`PHASE2_API.md`). This section explains how they
work, then how the code implements them.

A stencil is a fixed recipe for f'(x) from a few values of f near x: a list of offsets (nodes) and weights,

```
f'(x) ~= [ w_1 f(x + n_1 h) + w_2 f(x + n_2 h) + ... ] / (den * h)
```

The Phase 1 central difference is the simplest one: nodes {-1, +1}, weights {-1, +1}, den 2.

### Where the weights come from

Expand each sample by Taylor's theorem around x:

```
f(x + n h) = f + n h f' + (n h)^2/2! f'' + (n h)^3/3! f''' + ...
```

Multiply by its weight and add them all. Grouped by derivative, each term collects a *moment* of the stencil:

```
sum_k w_k f(x + n_k h) = M_0 f + M_1 h f' + M_2 h^2/2! f'' + M_3 h^3/3! f''' + ...

where  M_j = sum_k w_k * n_k^j
```

For the left side to equal `den * h * f'`, the weights must satisfy:
- **M_0 = 0:** the f term drops out.
- **M_1 = den:** f' gets the right coefficient.
- **M_2 = M_3 = ... = M_p = 0:** as many higher terms as possible cancel.

The first term that does not cancel is the truncation error:

```
error = M_(p+1) / ((p+1)! * den) * h^p * f^(p+1)(x)
```

p is the accuracy order: halving h divides the error by 2^p. N points give N unknown weights, enough to satisfy N
equations, from M_0 up to M_(N-1).

### Example: the four-point central stencil

Take nodes -2, -1, 1, 2. For a first derivative the weights are antisymmetric, `w(-n) = -w(n)`, so call them a, b,
-b, -a. Antisymmetry makes every even moment (M_0, M_2, M_4, ...) zero automatically, which leaves the odd ones:

```
M_1 = -2a - b - b - 2a   = -4a - 2b   = den
M_3 = -8a - b - b - 8a   = -16a - 2b  = 0     ->  b = -8a
                                               ->  den = 12a; take a = 1
```

The weights are {1, -8, 8, -1} with den 12:

```
f'(x) ~= (f(x-2h) - 8 f(x-h) + 8 f(x+h) - f(x+2h)) / 12h
```

M_4 is zero by symmetry, and the first surviving moment is M_5 = -64 + 16 = -48. The error is
`-48 / (5! * 12) * h^4 f^(5) = -h^4/30 * f^(5)`: order 4.

This is why central stencils pay off. Symmetry cancels the even terms for free, so **each extra pair of points gains
two orders**: two points give order 2, four give order 4, six give order 6.

### Example: a forward stencil

Take nodes 0, 1, 2, never to the left of x, with weights a, b, c:

```
M_0 = a + b + c  = 0
M_1 = b + 2c     = den
M_2 = b + 4c     = 0     ->  b = -4c,  den = -2c; take c = -1
```

The weights are {-3, 4, -1} with den 2:

```
f'(x) ~= (-3 f(x) + 4 f(x+h) - f(x+2h)) / 2h
```

The error is `M_3 / (3! * 2) * h^2 f''' = -h^2/3 * f'''`: order 2. With no symmetry to exploit, **one-sided stencils
gain one order per extra point**. That is why the forward rule needs three points to match the two-point central
difference, and why its error constant is twice as large (1/3 against 1/6).

### The price: rounding

Each sample carries a rounding error of about `eps * |f|`. The stencil multiplies those errors by its weights and
divides by `den * h`, so the rounding error is about

```
S * eps * |f| / h,      S = sum |w_k| / den
```

| Stencil  | Order | S    | Best error, about                    |
|----------|-------|------|--------------------------------------|
| central2 | 2     | 1    | eps^(2/3)                            |
| central4 | 4     | 1.5  | eps^(4/5)                            |
| central6 | 6     | 1.83 | eps^(6/7)                            |
| forward2 | 2     | 4    | eps^(2/3), about 3x central2's error |
| forward4 | 4     | 10.7 | eps^(4/5), about 7x central4's error |

Balancing the two errors, as section 2 did for p = 2, gives the best step `h ~ eps^(1/(p+1)) * scale` and a best error
of order `eps^(p/(p+1))`. Higher order allows a larger step and loses fewer digits, but the gains shrink: for `double`,
order 6 already loses only about 2 of the 16 digits. One-sided stencils have large weights and amplify rounding much
more, so they are for points where f cannot be evaluated on both sides.

Higher order also assumes that f is smooth across the whole stencil, which grows with the number of points. A kink or a
fast feature nearby can make central6 worse than central2. Ridders' method addresses this in Phase 3, with an error
estimate.

### Two consequences

- **Exact on polynomials.** A stencil of order p is exact for polynomials up to degree p, since every term it does not
  cancel involves f^(p+1) or higher. With power-of-two steps the result is even bit-exact in binary floating point,
  provided the code divides once, by `den * h`; multiplying by a rounded `1/12` would break that.
- **Richardson extrapolation in disguise.** Combining two central differences to cancel their h^2 terms,
  `(4 D(h) - D(2h)) / 3`, expands to exactly the four-point formula above. Fixed stencils and fixed-level Richardson
  extrapolation are the same thing; Ridders adds adaptivity and the error estimate.

The same moment equations give second derivatives: require M_2 = 2 * den and the other low moments zero. On nodes
{-1, 0, 1} that gives {1, -2, 1}, the classic `(f(x+h) - 2 f(x) + f(x-h)) / h^2`. So Phase 4 can reuse Phase 2's table
type and evaluator.

### In the code

The stencils are three class templates, `Central<P>` (even P from 2 to 10), `Forward<P>` and `Backward<P>` (P from 1 to
6). They are empty: the type only names the rule. The caps are design decisions: beyond `Central<10>` the model's error
would fall by a factor of only about 1.7 for f' in double and about 5 in a 50-digit type, and the one-sided rounding
factor S is already 27.7 at P = 6.

**The tables** are computed at compile time by `detail::make_table<N>`, for the N-th derivative (section 10). Each
family only supplies its nodes, and one generator turns nodes into weights by the Lagrange formula: the weight of node
k is `N! [t^N] prod_{j != k} (t - n_j) / (n_k - n_j)`, the N-th derivative at 0 of the polynomial through the
nodes, which is the solution of the moment conditions above. For `Central<4>` and node -1: `(t + 2)(t - 1)(t - 2) =
t^3 - t^2 - 4t + 4` has the t coefficient -4, and `(-1 + 2)(-1 - 1)(-1 - 2) = 6`, so the weight is -4/6 = -2/3, the
-8/12 derived above. The arithmetic is exact (`long long`, `std::gcd`, `std::lcm`), and the weights are put over one
common denominator. (Phase 2 used a closed form per family; the generator reproduces those tables bit for bit, and
also covers the higher orders, whose one-sided rules have no short closed form.)

Every value must fit 24 bits, so that it is exact in every real type, even `float`; the largest numerator is 2100 for
f' (`Central<10>`) and 6,222,216 for f''''.

**The evaluator** has two branches, chosen at compile time:

- Two nodes (`Central<2>`, `Forward<1>`, `Backward<1>`): `(f_hi - f_lo) / (x_hi - x_lo)`, dividing by the distance
  between the points actually evaluated, as in section 5. `Central<2>` is the Phase 1 code, bit for bit.
- More nodes: `sum_k numerator_k * f(x + n_k h)`, divided once by `denominator * h` (and N - 1 more times by h for the
  N-th derivative). Dividing by the outer-node distance would correct the rounding of the outer points only, and with
  more than two weights it makes the error larger (4x for `Central<4>`). One division keeps the results exact for
  polynomials when x and h are short binary fractions, which the tests use. Before each term is formed, `|f_k| <= max
  / |numerator_k|` is checked, and the partial sum with `sum_fits`: a constant 1e308 under `Central<4>` used to form
  -8e308 and then inf - inf, a compile error in a constant expression.

**The input check** of section 5 generalises with `K = max |node|` and `Q = max(2K, denominator)`: `h <= max / (2Q)`
and `|x| <= max - Q h`. Every point lies within K h of x, so at least about K h below max, and the divisor is at most
`Q h <= max / 2`. The margin is needed: with no margin (`Q = K`), `max - Q h` can round up by half an ulp, and a point
one ulp below max can then round up to infinity. For `Central<2>`, Q = 2: exactly the check of section 5. With more
nodes, the zero distance becomes "h does not move x on one of its sides" (`x + h == x` or `x - h == x`); both sides
are needed because at a power of the radix the spacing differs on the two sides of x. Each value of fn is checked
before it is added to the sum.

**The default step** is `root_epsilon<P + 1>() * max(|x|, 1)` for f', as derived above, without the stencil's
constant `(S / (p C))^(1/(p+1))`. Those constants lie between 1.4 and 2.0, and leaving one out costs at most a factor
of about 1.6 in error (1.63 for `Central<10>`). For the N-th derivative it is `root_epsilon<P + N>()` (section 10).

## 9. Ridders' extrapolation

A stencil needs a good step, and cannot say how accurate its result is. Ridders' method (`derivative<Ridders>`)
solves both: it computes central differences at a sequence of steps and extrapolates them, until rounding stops the
improvement, and it estimates its own error on the way. The full specification, with the evidence for every rule, is
`PHASE3_API.md`.

### Extrapolating the central difference

Section 2 showed that the central difference has an error expansion in even powers of h:

```
D(h) = f'(x) + c_2 h^2 + c_4 h^4 + c_6 h^6 + ...
```

So D(h) and D(h/2) share the same c_2, and `(4 D(h/2) - D(h)) / 3` cancels the h^2 term. That is Richardson
extrapolation (the four-point stencil of section 8 is exactly this). Repeating it gives a tableau: `A[i][0] = D(h/2^i)`
and, in Burden's difference form,

```
A[i][j] = A[i][j-1] + (A[i][j-1] - A[i-1][j-1]) / (4^j - 1)
```

where column j has cancelled the terms up to h^(2j). Each new level i halves the step and adds one row; only the
previous row is kept.

### Knowing when to stop

Moving down the tableau lowers the truncation error, but each halving of h doubles the rounding error of D. So each
extrapolated entry gets an error estimate, its change from its diagonal parent `|A[i][j] - A[i-1][j-1]|`, but never
less than the rounding floor `eps (|f(x+h)| + |f(x-h)|) / 2h` of the difference it came from. Without the floor, the
estimate (the smallest of many noisy differences) was too small in more than half of the prototype's cases. The
method keeps the entry with the smallest estimate and stops after two levels that improve nothing, when an entry
equals its parent exactly (nothing left to learn), or at the resolution of x. The result's `abs_error` is also at
least the distance from the best entry to the latest diagonal entry: a best entry that only looked good is caught
there.

### Restarting below a bad step

Two things make a step bad. If fn is not finite at x +- h (a domain edge, a pole on a sample), every earlier level's
interval contained the bad point, so the tableau is discarded and restarted at h/2. And if the result has at most two
correct bits while its estimate is far above the rounding floor, fn misbehaves within the step without producing an
infinity (a pole between the samples, or aliasing: sin(64 pi x) vanishes at every sample of h = 1/32 and 1/64), so
the method restarts below the smallest step it tried. (Until Phase 5 it restarted at half its starting step, which
re-evaluated the same large steps: with a pole 1e-4 from x the levels ran out first.) At most 32 levels are attempted
in all, 64 evaluations of fn.

### What it costs and gives

For double, typically 10 to 16 evaluations and an error near 1e-15 relative, with an estimate that covered the
actual error in 99% of the prototype's 2,700 random cases and was about 7 times conservative at the median. Compare
`Central<6>`: 6 evaluations, about 1e-14, no estimate. The estimate assumes fn is accurate to about eps; a noisier fn
makes it optimistic.

## 10. Higher derivatives

`derivative<M, N>` estimates the N-th derivative, for N from 1 to 4, with every method. N is an integer or an
`order`, a scoped enum (`order::second` is 2):

```cpp
auto const a = derivative<Central<2>, 2>(f, x);              // (f(x+h) - 2f(x) + f(x-h)) / h^2
auto const b = derivative<Central<2>, order::third>(f, x);   // (f(x+2h) - 2f(x+h) + 2f(x-h) - f(x-2h)) / 2h^3
auto const c = derivative<Backward<2>, 4>(f, x);             // Chapra's (3f(x) - 14f(x-h) + ... - 2f(x-5h)) / h^4
auto const d = derivative<Ridders, 4>(f, x);                 // f'''' with an error estimate
```

### The rules

The weights come from the same generator as for f' (section 8), with N! [t^N] instead of the t coefficient, and the
moment conditions become M_N = N! denominator, M_j = 0 for every other j up to N + P - 1. P stays the accuracy order,
the O(h^P) of the error, so it is not the number of points: `Central<P>` uses the P + N - 1 integers around 0 (without
0 for odd N, since an odd derivative gives the centre no weight), the one-sided rules P + N. For the third derivative,
`Central<2>` is Chapra's four-point O(h^2) formula and `Central<4>` his six-point O(h^4) one.

### The step, and why accuracy falls with the order

The rounding error of fn's values is divided by h^N, so the error model becomes `C h^P |f^(P+N)| + S eps |f| / h^N`,
smallest near `h = eps^(1/(P+N))`, with an error of about `eps^(P/(P+N))`. So the default step is
`root_epsilon<P + N>() * max(|x|, 1)`, and the accuracy falls quickly with N: for `Central<2>` in double about 4e-11
for f', 6e-8 for f'', 2e-6 for f''' and 1e-4 for f''''. In float, f'''' is mostly rounding noise, which is why it
warns at compile time (below). The step's sensitivity grows with N too: for f'''', misjudging fn's scale by a factor 2
costs up to 16 in error.

The code divides by h N times instead of once by h^N: h^4 underflows for h below 1.2e-77 in double, and the repeated
divisions cannot overflow spuriously (for h >= 1 each quotient is at most the one before, otherwise the last is the
largest).

### f'' from the actual offsets

As the two-node f' rules divide by the distance between the points actually evaluated, `Central<2>` for f'' uses the
actual offsets `h+ = (x + h) - x` and `h- = x - (x - h)`:

```
f'' ~= ((f(x+h) - f(x)) / h+ - (f(x) - f(x-h)) / h-) / ((h+ + h-) / 2)
```

This is exact for quadratics whatever the rounding of the points; the plain formula picks up an error of about
eps |x f'| / h^2 from that rounding, as much as its own error at the default step. In Ridders' prototype, the plain
formula let the estimate fall 8.7 times below the actual error for float, against 1.4.

### Ridders

Ridders' method extends unchanged: its base is the `Central<2>` rule of order N, whose error expansion has even powers
of h for every N, so the same tableau applies. The rounding floor of a level is the same sum `eps sum |w_k| |f_k|`
over `denominator * h^N`, saturated at every division. fn(x), part of the base for even N, is evaluated once. The
default initial step is larger, `max(|x|, 1) / 4` instead of `/ 32`, since the rounding grows as 1/h^N: for f'' it is
7 times as accurate as `/ 32`. For double the method typically reaches 4e-13, 3e-11 and 7e-10 relative for f'' to
f'''', in 13 to 45 evaluations, at most 129.

### One side only: `ForwardRidders` and `BackwardRidders`

Ridders' method samples both sides of x, so it fails where fn is undefined on one side: at a domain edge, or at the
end of a correlation's range. `ForwardRidders` extrapolates the `Forward<1>` rule of order N instead, which samples x
and above only, and gives an error estimate as `Ridders` does; `BackwardRidders` samples x and below.

```cpp
auto const f = [](double x) { return x < 1.0 ? std::nan("") : std::log(x); }; // undefined below 1
auto const a = derivative<Ridders>(f, 1.0);           // non_finite: its samples reach below 1
auto const b = derivative<ForwardRidders>(f, 1.0);    // 1 to about 1e-13, with an error estimate
auto const c = derivative<BackwardRidders>(g, 2.0);   // g defined up to 2 only
```

The forward rule's error has every power of h, `D(h) = f^(N) + c_1 h + c_2 h^2 + ...`, so column j of the tableau
divides by `2^j - 1` and gains one order per level instead of two. That makes it less accurate than `Ridders`: about
2e-13, 1e-10, 4e-8 and 2e-6 relative for f' to f'''' in double, against 1e-14 to 7e-10, still 10 to 25 times better
than the best one-sided stencil. Those columns also amplify the newest level's rounding about 3.5 times, against 1.4
for the central ones, so the rounding floor of each entry is carried through the tableau rather than taken from the
newest level; with it the estimate covered the error in 99.4% of the prototype's random cases for f'.

A one-sided method never sees the NaN beyond an edge, so near one it learns only from untrusted results, and
restarting below the smallest step tried is what lets it reach the scale of fn. From the default step
`max(|x|, 1) / 32` it handles sqrt or log down to about 1e-9 from 0; closer, its result is wrong but flagged (an
estimate of at least |value| / 4): pass a step on fn's own scale, such as `x / 8`. `BackwardRidders` is
`ForwardRidders` applied to fn(-t) at -x, which gives exactly the points x - k h. Where fn is defined on both sides,
use `Ridders`, also near an edge.

### The low-precision warning

A derivative above the first in a type whose epsilon is above 1e-7 (float) gives a compiler warning, by every method:
the median relative error of f'' in float is 4e-4 with `Central<2>` and 1e-5 with `Ridders`, and f'''' is 2e-2 and
1e-3. (Phase 4 warned for f'''' only.) Standard C++ cannot
emit a warning from a template instantiation (`#warning` acts at preprocessing), so the warning comes from
`[[deprecated]]`: each public overload has a twin with the extra
constraint `detail::noisy<N, T>`, which overload resolution prefers for such calls because it is more constrained. The
warning is reported at the caller's line, so a pragma there silences it, and `NUMERIXX_NO_PRECISION_WARNING` turns it
off everywhere.

## 11. Limits, and what comes next

- **Fixed rules only.** The methods are the three stencil families and the three Ridders methods; a user cannot
  supply their own
  weights, and a method cannot be chosen at run time except by a `switch` over `derivative<M>` calls.
- **Orders up to 4, one at a time.** Several orders from shared evaluations (for Halley's method, say) and mixed
  partial derivatives are deferred.
- **Only the Ridders methods estimate their error.** The stencils leave `abs_error` empty.
- **The step assumes f varies on the scale `max(|x|, 1)`.** For a stencil and a faster-varying function, pass a step,
  or use Ridders.

`DERIV_TECHNIQUES.md` surveys the alternatives (complex step, automatic differentiation, contour integrals) and why
they are deferred.

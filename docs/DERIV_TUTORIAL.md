# How `nxx::deriv` works: a guided tour

This tutorial walks through the current implementation of numerical derivatives in Numerixx 2: two headers,
`include/numerixx/math.hpp` and `include/numerixx/deriv.hpp`. It starts with how to call the library, then explains the
mathematics behind the central difference and its default step, goes through the code line by line, including the C++
details that let the same code run on `double`, on Boost.Multiprecision types and at compile time, and ends with the
higher-order and one-sided stencils (section 8).

## 1. Using it

```cpp
#include <numerixx/deriv.hpp>

auto const f = [](double x) { return std::exp(x); };

auto const d = nxx::deriv::derivative(f, 1.0);         // default step
auto const e = nxx::deriv::derivative(f, 1.0, 1e-4);   // explicit step h

if (d) std::println("{}", *d);                         // 2.7182818284569006 (exact: 2.7182818284590451)

using nxx::deriv::Central, nxx::deriv::Forward;
auto const g = nxx::deriv::derivative<Central<4>>(f, 1.0);   // higher order: 4 evaluations
auto const k = nxx::deriv::derivative<Forward<2>>(f, 1.0);   // one-sided: never evaluates f below x
```

There are two overloads, with and without a step. The optional template argument chooses the stencil (section 8); by
default it is `Central<2>`, the central difference of sections 2 to 7. Both return `std::expected<T, nxx::deriv::error>`: either the
estimate, or `error::non_finite` if no finite estimate exists (section 5 lists the cases). Numerixx itself never
throws. The demo program
(`demo/deriv_demo.cpp`) runs every example in this tutorial.

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

`derivative` handles every stencil (section 8 covers the general case). For the default `Central<2>`, its compile-time
checks and branches reduce to this:

```cpp
template<typename S = Central<2>, real T, typename F>
requires detail::stencil<S> && std::invocable<F &, T const &>
         && std::convertible_to<std::invoke_result_t<F &, T const &>, T>
[[nodiscard]] constexpr std::expected<T, error> derivative(F &&fn, T const &x, std::type_identity_t<T> const &h)
{
    assert(h > T(0));
    // With |x| + 2h <= max, the points and their distance are finite, so nothing overflows before fn is called.
    T const largest = (std::numeric_limits<T>::max)();
    if (!(h <= largest / T(4) && math::abs(x) <= largest - T(2) * h)) return std::unexpected(error::non_finite);
    T const xm = x - h;
    T const xp = x + h;
    T const dx = xp - xm;
    if (dx == T(0)) return std::unexpected(error::non_finite); // h is below the resolution of x
    T const fm = static_cast<T>(fn(xm));
    T const fp = static_cast<T>(fn(xp));
    if (!math::isfinite(fm) || !math::isfinite(fp)) return std::unexpected(error::non_finite);
    T const d = (fp - fm) / dx;
    if (!math::isfinite(d)) return std::unexpected(error::non_finite); // the difference or the quotient overflowed
    return d;
}
```

Each line has a reason.

### The signature

- **`typename S = Central<2>`** and **`detail::stencil<S>`**: the stencil is a type, the first template argument, so
  `derivative(f, x)` keeps working while `derivative<Central<4>>(f, x)` picks another rule. The concept accepts only
  `Central`, `Forward` and `Backward` in their allowed orders; anything else is rejected by the constraint, not by a
  hard error.
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
4. **The result** is checked last: both values can be finite while their difference overflows (`max - (-max)`), or the
   quotient can overflow when the distance is tiny.

Checking first rather than last matters for two reasons. fn is only ever called at finite points. And in a constant
expression, compilers reject arithmetic that produces an infinity or a NaN: GCC rejects floating-point overflow and NaN,
Clang rejects NaN and division by zero. Because the checks come first, those failures are ordinary error values even in
a `static_assert`. The one exception is check 4: an overflow of the final difference or quotient is detected after the
fact, so in constant evaluation on GCC it is a compile error rather than an error value.

### Exceptions

Numerixx never throws and never catches. If fn, or an operation of T, throws, the exception passes through
`derivative` unchanged. Catching it would make the library unusable in builds without exceptions, such as
Emscripten's default.

### The default-step overload

```cpp
template<int Order, real T>
constexpr T default_step(T const &x)                     // in nxx::deriv::detail
{
    T const scale = math::abs(x) < T(1) ? T(1) : math::abs(x);
    return math::root_epsilon<Order + 1, T>() * scale;
}

template<typename S = Central<2>, real T, typename F>
requires /* as above */
[[nodiscard]] constexpr std::expected<T, error> derivative(F &&fn, T const &x)
{ return nxx::deriv::derivative<S>(fn, x, detail::default_step<detail::make_table(S{}).order>(x)); }
```

For `Central<2>` (order 2) this is section 2's `eps^(1/3) * max(|x|, 1)`. The call is written `nxx::deriv::derivative`, fully qualified: an
unqualified call would also look in the namespaces of T and F (argument-dependent lookup), where a user type could
provide its own `derivative` and silently take over.

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

The stencils are three class templates, `Central<P>` (even P from 2 to 10), `Forward<P>` and `Backward<P>` (P from 1
to 6). They are empty: the type only names the rule. The caps are numerical: beyond `Central<10>` the gain is under
half a digit even for 50-digit types, and the one-sided rounding factor S is already 27.7 at P = 6.

**The tables** are computed at compile time by `detail::make_table`, from closed forms in exact integer arithmetic
(`long long`, `std::gcd`, `std::lcm`), and put over one common denominator:

- Central, P = 2m: the coefficient of `f(x + k h) - f(x - k h)` is `(-1)^(k+1) (m!)^2 / (k (m-k)! (m+k)!)`. For m = 2
  this gives 2/3 and -1/12, the {1, -8, 8, -1}/12 derived above.
- Forward, P: `w_0 = -(1 + 1/2 + ... + 1/P)` and `w_k = (-1)^(k+1) C(P, k) / k`. For P = 2: -3/2, 2, -1/2.
- Backward: Forward mirrored.

The largest numerator is 2100 (`Central<10>`), exact in every real type, even `float`.

**The evaluator** has two branches, chosen at compile time:

- Two nodes (`Central<2>`, `Forward<1>`, `Backward<1>`): `(f_hi - f_lo) / (x_hi - x_lo)`, dividing by the distance
  between the points actually evaluated, as in section 5. `Central<2>` is the Phase 1 code, bit for bit.
- More nodes: `sum_k numerator_k * f(x + n_k h)`, divided once by `denominator * h`. Dividing by the outer-node distance
  would correct the rounding of the outer points only, and with more than two weights it makes the error larger (4x
  for `Central<4>`). One division keeps the results exact for polynomials when x and h are short binary fractions,
  which the tests use.

**The input check** of section 5 generalises with `K = max |node|` and `Q = max(2K, denominator)`: `h <= max / (2Q)`
and `|x| <= max - Q h`. Every point lies within K h of x, so at least about K h below max, and the divisor is at most
`Q h <= max / 2`. The margin is needed: with no margin (`Q = K`), `max - Q h` can round up by half an ulp, and a point
one ulp below max can then round up to infinity. For `Central<2>`, Q = 2: exactly the check of section 5. With more
nodes, the zero distance becomes "h does not move x on one of its sides" (`x + h == x` or `x - h == x`); both sides
are needed because at a power of the radix the spacing differs on the two sides of x. Each value of fn is checked
before it is added to the sum.

**The default step** is `root_epsilon<P + 1>() * max(|x|, 1)`, as derived above, without the stencil's constant
`(S / (p C))^(1/(p+1))`. Those constants lie between 1.4 and 2.0, and leaving one out costs at most a factor 1.64 in
error.

## 9. Limits, and what comes next

- **Fixed rules only.** The stencils are the three families; a user cannot supply their own weights, and a stencil
  cannot be chosen at run time except by a `switch` over `derivative<S>` calls.
- **No error estimate.** The caller cannot tell how accurate a result is. Ridders' extrapolation (Phase 3) provides
  one and is robust to a poor initial step.
- **The step assumes f varies on the scale `max(|x|, 1)`.** For faster-varying functions, pass a step.

`DERIV_TECHNIQUES.md` surveys the alternatives (complex step, automatic differentiation, contour integrals) and why
they are deferred.

# A walkthrough of `include/numerixx/deriv.hpp`

## How to read this document

This is a walkthrough of the implementation: what each part of the code does, how the parts fit together, and why
they are written the way they are. The mathematics behind the methods (finite differences, truncation and rounding,
the choice of step, Richardson extrapolation) is in the book, and this document takes it as known.

The header is 1951 lines (as of d89b50a). About two thirds of it is doc comments, and most of those record *why*
something is the way it is: a bug that was found, a measurement from a prototype, a compiler that rejected the
obvious form. This document explains the code in the order you need to understand it rather than the order of the
file.

Two short sections come first, A and B, on two constraints that shape almost every line: the code must work in
constant expressions, and with user-defined real types. Then sections 0 to 14 go through every declaration. Each
section shows the exact code of every function it explains, copied from the header with the doc comments removed (this
text explains them; the `//` comments in the code are kept), so you don't need the header open. The line
numbers still name where it is. Inside the explanations, some fragments are repeated in a shortened form; they are
marked as such. Short notes marked **C++ note** explain a language feature where the code uses one that isn't
everyday C++.

The numbers in the examples come from running the library (GCC 16, double). The tables that show the inside of
Ridders' method, which the library doesn't expose, come from the Python prototype `.scratch/onesided/ridders5.py`,
which gives the same results bit for bit. I checked every example result against the library.

---

# Before the code

## A. Constant evaluation, and why the code checks before it computes

Every function in the header is `constexpr`, and some are `consteval`. That means a derivative can be computed while
the program is compiled, as in

```cpp
static_assert(derivative<Central<2>>([](double x) { return x * x; }, 3.0, 0.5)->value == 6.0);
```

The test suite has many such `static_assert`s. They prove that the arithmetic is exact where it should be, and they
run on every compiler.

**The cost is a strict rule.** In a constant expression, an operation that produces a NaN is not allowed, and on GCC
an overflow to infinity is not allowed either. Such an operation doesn't give a NaN or inf value that you could test
afterwards; it makes the expression non-constant, which is a **compile error**. At run time, `x / 0.0` gives `inf`
and you can test it; in a `static_assert` it fails the build.

The library promises that every failure is an *error value* (`error::non_finite`), at run time and at compile time
alike. Hence the habit you will see everywhere below: **check whether an operation would overflow, and only then
do it**. Section 6 contains the helpers that answer "would `a + b` overflow?" without computing `a + b`. Without
constant evaluation, half of that code would be unnecessary.

## B. User-defined real types and expression templates

The algorithms accept any type `T` that satisfies the concept `nxx::real` (`math.hpp:38`): `float`, `double`,
`long double`, or a user type such as Boost.Multiprecision's `cpp_bin_float_50` (50 decimal digits). This shapes the
code in three ways.

**The concept** (math.hpp, lines 29-52; its `//` comment explains each requirement):

```cpp
// Constructible from int, because constants are written T(1), T(2): a double literal would silently lose precision
// in a wider type. The operators' results need only convert to T, because Boost.Multiprecision with expression
// templates returns an expression object from a + b, not a T. The numeric_limits members are those the headers
// read: radix and digits in constant expressions (the exact-weights and float-warning constraints of derivative),
// radix at least 2 because root_epsilon and the overflow guards scale by it; epsilon() and max() are
// copy-initialised into T, while denorm_min() is read as T(denorm_min()), so a type that inherits
// numeric_limits<double> and constructs from double only explicitly is accepted. An object type first: for a
// function type, which overload resolution may try as a deduced T (a derivative passed as a free function), the
// standard library's numeric_limits is ill-formed, a hard error rather than a failed constraint.
template<typename T>
concept real = std::is_object_v<T> && std::numeric_limits<T>::is_specialized && !std::numeric_limits<T>::is_integer
               && std::numeric_limits<T>::is_bounded && !std::numeric_limits<T>::is_exact
               && std::constructible_from<T, int> && std::copyable<T> && std::totally_ordered<T>
               && requires(T const a, T const b) {
                      requires std::numeric_limits<T>::radix >= 2 && std::numeric_limits<T>::digits >= 1;
                      { std::numeric_limits<T>::epsilon() } -> std::convertible_to<T>;
                      { (std::numeric_limits<T>::max)() } -> std::convertible_to<T>;
                      requires std::constructible_from<T, decltype(std::numeric_limits<T>::denorm_min())>;
                      { a + b } -> std::convertible_to<T>;
                      { a - b } -> std::convertible_to<T>;
                      { a * b } -> std::convertible_to<T>;
                      { a / b } -> std::convertible_to<T>;
                      { -a } -> std::convertible_to<T>;
                  };
```

**1. Constants are written `T(1)`, `T(2)`, never `1.0`.** A `double` literal has 53 bits. In a 50-digit type,
`T(0.1)` is the 50-digit value of the double nearest 0.1, not 0.1. Small integers are exact in every type, so the
code only ever writes integers.

**2. Expression templates, and why the code never writes `auto` for a number.** Boost.Multiprecision, by default,
doesn't compute `a + b` when you write it. It returns a small object that *remembers* "a plus b", and the sum is
computed later, when the object is assigned to a number. This lets a chain like `a + b * c` run without temporaries.
The catch: the object holds *references* to `a` and `b`. Write `auto s = a + b;` and `s` is that object, not a
number. If `a` was a temporary, `s` now refers to a destroyed value. So the code writes `T const s = T(a + b);`, or
wraps an expression in `T(...)`, to force evaluation on the spot. You'll see `T(...)` wrapped around many expressions
that look as if they didn't need it; that is why.

**3. Calls into the library are qualified (`detail::fits(...)`, `math::abs(...)`).** C++ looks up an unqualified
function name in the namespaces of its arguments too: this is argument-dependent lookup, ADL. If a user's type `T`
lives in a namespace that also has a function `fits`, an unqualified `fits(x, h, q)` might call the user's function.
Qualifying the name prevents that.

---

# The code

## 0. The map

| Lines     | What                                                                        | Section |
|-----------|-----------------------------------------------------------------------------|---------|
| 1-38      | file comment, includes                                                      |         |
| 50-109    | public types: `error`, `Estimate<T>`, `order`                               | 1       |
| 112-219   | compile-time gates: valid orders, exactness, the float warning              | 3       |
| 222-330   | public method tags                                                          | 1       |
| 332-389   | which types are methods; `base<M>`                                          | 3       |
| 391-558   | weight tables: `Table`, `lagrange_table`, `make_table`                      | 4       |
| 560-629   | default stencil step; `exact_weights`                                       | 5, 3    |
| 631-836   | overflow-safe arithmetic                                                    | 6       |
| 838-1143  | building blocks: differences, rounding floors, weighted sums                | 7       |
| 1145-1219 | `stencil_derivative`: the fixed-step methods                                | 8       |
| 1221-1575 | Ridders: constants, `Level`, `ridders_level`, `propagated_floor`, `ridders` | 9       |
| 1577-1640 | `Mirrored`, `derivative_with_step`                                          | 10, 2   |
| 1644-1705 | `default_step`, `derivative_default`                                        | 5, 2    |
| 1707-1824 | the four public `derivative` overloads                                      | 2       |
| 1826-1911 | `derivative_of`                                                             | 11      |
| 1915-1951 | `std::formatter<error>`                                                     | 12      |

From `math.hpp` the header uses:

- the concept `nxx::real` and `nxx::detail::callable`;
- `math::abs` and `math::isfinite`, written with comparisons only so that they work for any `real` and in constant
  expressions;
- `math::root_epsilon`.

**The includes** (deriv.hpp, lines 22-38):

```cpp
#pragma once

#include <numerixx/math.hpp>

#include <array>
#include <concepts>
#include <cstddef>
#include <expected>
#include <format>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <string_view>
#include <type_traits>
#include <utility>
```

**The two helpers** (math.hpp, lines 83-98):

```cpp
template<real T>
constexpr T abs(T const &x)
{
    if (x < T(0)) return T(-x); // T(...) evaluates a possible expression template
    return x;
}

// A NaN is the only value that compares unequal to itself. max is parenthesised so that a max macro, as defined
// by <windows.h>, cannot expand there.
template<real T>
constexpr bool isfinite(T const &x) { return x == x && math::abs(x) <= (std::numeric_limits<T>::max)(); } // NOLINT
```

**The call path.** A call goes through three layers.

1. **The public `derivative` overloads** check, at compile time, that the call makes sense.
2. **Two internal functions** check the run-time inputs and choose the algorithm.
3. **One of two engines** does the work.

```
derivative<M, N>(fn, x [, h])                     public, 4 overloads (2 plain + 2 deprecated twins)
 |-- without h: detail::derivative_default        checks x, then h = default_step<M, N>(x)
 `-- detail::derivative_with_step                 checks h, then chooses by M:
      |-- a stencil (Central/Forward/Backward<P>) -> detail::stencil_derivative<M, N>
      |-- Ridders, ForwardRidders                 -> detail::ridders<base<M>, N>
      `-- BackwardRidders                         -> detail::ridders<Forward<1>, N> on fn(-t) at -x,
                                                     the value negated for odd N
```

Everything below `derivative` lives in `namespace detail`. It isn't part of the public API, and the Doxygen comments
there are hidden from the user documentation by `@cond INTERNAL`. The header opens and closes `namespace detail`
four times, between the public declarations; those namespace lines, and `namespace nxx::deriv` around everything,
are the only code lines this document leaves out.

---

## 1. The public vocabulary (lines 43-109, 222-330)

### `error` (50-58)

```cpp
enum class error {
    non_finite
};
```

There is only one enumerator. Every failure maps to it:

- `x` isn't finite;
- the step isn't positive and finite;
- `fn` returned inf or NaN;
- something would overflow;
- the step is so small that `x + h` rounds back to `x`.

A single code is enough for now because the caller's response is the same in every case: this point and step don't
give a derivative. The doc comment says more values may be added later. That is why it tells callers to test for
success (`if (r)`) rather than compare with `error::non_finite`, and to give a `switch` on it a `default`. Section 12
explains how it prints.

### `Estimate<T>` (78-95)

```cpp
template<real T>
struct Estimate
{
    T value;

    std::optional<T> abs_error;

    constexpr operator T() const { return value; }

    friend constexpr bool operator==(Estimate const &, Estimate const &) = default;
};
```

This is what a successful call returns: the derivative, and an estimate of its absolute error when the method can
give one. The stencils can't, so `abs_error` is a `std::optional`, empty for them.

**It is a plain aggregate.** There are no constructors, so the code builds it with braces:
`Estimate<T>{ d, std::nullopt }`.

**The implicit conversion to `T`.** This lets you write `double s = *r;` and use the estimate where a number is
expected, as most callers want just the value. The doc comment lists where the conversion surprises:

- `if (*r)` *compiles*, but it tests the value against 0 through the conversion, not whether there is a result.
  Clang with `-Wfloat-conversion` rejects it.
- `auto v = *r;` gives you an `Estimate`, not a `double`.
- `int n = *r;` silently truncates.

**The defaulted `==`.** Without it, `a == b` on two Estimates would still compile: both would convert to `T` and
only the values would be compared, ignoring `abs_error`. The defaulted operator compares both members, which is what
equality should mean.

### `order` (104-109)

```cpp
enum class order : int {
    first = 1,
    second = 2,
    third = 3,
    fourth = 4
};
```

It exists only for readability: `derivative<Ridders, order::second>` reads better than `derivative<Ridders, 2>`, and
the two mean the same. It is a *scoped* enum, so `first` and `second` don't leak into `nxx::deriv`, and an `order`
doesn't convert to `int` by accident. That last property is why the public templates take their order as `auto N`
rather than `int N` (section 2).

### The method tags (222-330)

```cpp
template<int P>
struct Central
{
    static_assert(detail::valid_central(P), "Central<P> needs an even order P from 2 to 10");
};

template<int P>
struct Forward
{
    static_assert(detail::valid_one_sided(P), "Forward<P> needs an order P from 1 to 6");
};

template<int P>
struct Backward
{
    static_assert(detail::valid_one_sided(P), "Backward<P> needs an order P from 1 to 6");
};

struct Ridders
{
};

struct ForwardRidders
{
};

struct BackwardRidders
{
};
```

**C++ note: tag types.** These six types have no members and no code. They are *tags*: their only job is to be
passed as a template argument, `derivative<Ridders>(...)`, so that the compiler can choose an algorithm with
`if constexpr` or overloading. At run time they don't exist. The algorithms themselves are free functions in
`detail`.

**What the `static_assert` does, and doesn't do.** The assertion in `Central<P>` fires only when the *body* of
`Central<3>` is instantiated, for example by writing `Central<3>{}`. Merely naming `Central<3>` as a template
argument, as in `derivative<Central<3>>(fn, x)`, doesn't instantiate its body. That call is rejected earlier and more
quietly, by the constraint `method<M>` (section 3), so `derivative` with a bad order is "no matching function", not a
hard error. The test `tests/test_stencils.cpp:109-112` relies on exactly this. The `static_assert` is a second line of
defence for someone who creates an invalid tag object directly.

**The doc comments are the user-facing specification.** The comments on these tags (accuracy, cost, default step,
limits, examples) are what users read in the documentation. This document explains the code that makes those statements
true.

---

## 2. Following a call from the top (lines 1707-1824, 1617-1640, 1695-1703)

Take `derivative<Central<2>>(fn, 1.0)` and follow it down.

### The four public overloads

| Lines | Signature                                                          | Extra constraint         |
|-------|--------------------------------------------------------------------|--------------------------|
| 1742  | `derivative(F &&fn, T const &x, std::type_identity_t<T> const &h)` | none                     |
| 1773  | `derivative(F &&fn, T const &x)`                                   | none                     |
| 1796  | same as 1742, marked `[[deprecated("...")]]`                       | `&& detail::noisy<N, T>` |
| 1817  | same as 1773, marked `[[deprecated("...")]]`                       | `&& detail::noisy<N, T>` |

**The code** (deriv.hpp, lines 1740-1824; the two `//` lines belong to the first overload):

```cpp
// h is a non-deduced context (type_identity_t), so T comes from x alone and a double step converts to a
// multiprecision T instead of making deduction fail. N is auto, because the scoped order does not convert to int.
template<typename M = Central<2>, auto N = 1, real T, typename F>
requires detail::method<M> && detail::valid_order<N> && nxx::detail::callable<F, T>
         && detail::exact_weights<M, N, T>
[[nodiscard]] constexpr std::expected<Estimate<T>, error>
    derivative(F &&fn, T const &x, std::type_identity_t<T> const &h)
{
    // Qualified, so that argument-dependent lookup cannot pick a function from the namespace of T or F.
    return detail::derivative_with_step<M, static_cast<int>(N)>(fn, x, h);
}

template<typename M = Central<2>, auto N = 1, real T, typename F>
requires detail::method<M> && detail::valid_order<N> && nxx::detail::callable<F, T>
         && detail::exact_weights<M, N, T>
[[nodiscard]] constexpr std::expected<Estimate<T>, error> derivative(F &&fn, T const &x)
{ return detail::derivative_default<M, static_cast<int>(N)>(fn, x); }

template<typename M = Central<2>, auto N = 1, real T, typename F>
requires detail::method<M> && detail::valid_order<N> && nxx::detail::callable<F, T>
         && detail::exact_weights<M, N, T> && detail::noisy<N, T>
[[deprecated(
    "derivatives above the first in a type with 7 or fewer significant digits keep few correct digits: "
    "use double, or silence this warning at the call or with NUMERIXX_NO_PRECISION_WARNING")]] [[nodiscard]]
constexpr std::expected<Estimate<T>, error> derivative(F &&fn, T const &x, std::type_identity_t<T> const &h)
{ return detail::derivative_with_step<M, static_cast<int>(N)>(fn, x, h); }

template<typename M = Central<2>, auto N = 1, real T, typename F>
requires detail::method<M> && detail::valid_order<N> && nxx::detail::callable<F, T>
         && detail::exact_weights<M, N, T> && detail::noisy<N, T>
[[deprecated(
    "derivatives above the first in a type with 7 or fewer significant digits keep few correct digits: "
    "use double, or silence this warning at the call or with NUMERIXX_NO_PRECISION_WARNING")]] [[nodiscard]]
constexpr std::expected<Estimate<T>, error> derivative(F &&fn, T const &x)
{ return detail::derivative_default<M, static_cast<int>(N)>(fn, x); }
```

Each one is a one-line body that forwards to `detail`. All four carry the same core constraint:

```cpp
template<typename M = Central<2>, auto N = 1, real T, typename F>
requires detail::method<M> && detail::valid_order<N> && nxx::detail::callable<F, T>
         && detail::exact_weights<M, N, T>
```

Read it as: M is one of the six methods, N is an order from 1 to 4, `fn` can be called with a `T` and returns
something convertible to `T`, and `T` can hold the method's integer weights exactly. If any of these fails, the
overload isn't a candidate, and the compiler reports that no overload matches.

### Why there are "deprecated twins"

The user decided that a derivative above the first in `float` (or any type as coarse) should produce a compiler
warning, because the result has few correct digits. The question was how to produce a warning that:

1. appears at the **user's** line, not deep inside the header;
2. can be silenced at that line with the usual pragma;
3. leaves the result unchanged.

A `static_assert` is an error, not a warning. A `#warning` can't depend on template arguments. But
`[[deprecated("message")]]` on a function makes every compiler warn **wherever that function is called**, with the
message, at the call site.

So for every plain overload there is a twin with the same signature, the same body, the `[[deprecated]]` attribute,
and one extra constraint, `noisy<N, T>`, which is true only for "N >= 2 in a coarse type". The compiler picks the
twin only in exactly those cases, and the next note explains how.

**C++ note: constraint subsumption.** When two overloads both match and are otherwise equally good, C++20 picks the
*more constrained* one. One constraint is more constrained than another if it implies it, in the formal sense the
standard calls *subsumption*. The twin's constraint is `A && B && C && D && noisy`; the plain one is `A && B && C &&
D`, written with the very same concepts. So the twin's constraint implies the plain one, and when both are
satisfied, the twin wins. When `noisy` is false, the twin isn't viable, and the plain overload is the only
candidate. No ambiguity is possible either way.

The twins' doc comments say how to silence the warning: `#pragma` with `-Wdeprecated-declarations` (GCC, Clang) or
C4996 (MSVC) around the call, or define `NUMERIXX_NO_PRECISION_WARNING` to switch it off everywhere.

### The details of the signature

**`auto N = 1` instead of `int N = 1`.** A non-type template parameter declared `auto` accepts a value of any type.
That lets `N` be `2`, `std::size_t(2)` (handy in generic code over index sequences) or `order::second`, which, being
a scoped enum, wouldn't convert to an `int` parameter. `valid_order<N>` (section 3) then checks that it is one of the
allowed types and values, and the body converts it with `static_cast<int>(N)` once, so everything below works with a
plain `int`.

**`std::type_identity_t<T> const &h`.** Without it, `T` would be deduced from both `x` and `h`. Then
`derivative(fn, cpp_bin_float_50(1), 1e-3)` would deduce `T` as `cpp_bin_float_50` from `x` and as `double` from
`h`, a conflict, so the call fails.

**C++ note: non-deduced contexts.** `std::type_identity_t<T>` is just `T`, but written through a type alias that the
compiler can't "see through" when it deduces template arguments. A parameter of that form doesn't take part in
deduction. So `T` comes from `x` alone, and `h` is then converted to that `T`. A `double` step works for any `T`.

**`F &&fn`, but `fn` is passed on as a plain lvalue.**

**C++ note: forwarding references.** `F &&` in a deduced context binds to anything: a temporary lambda, a named
function object, a function. It would usually be passed on with `std::forward<F>(fn)`. Here it is deliberately
passed as `fn`, a named lvalue, because the algorithms call `fn` many times. Moving from it on the first call would
break the second.

For the same reason the check `nxx::detail::callable<F, T>` (`math.hpp:68`) asks whether `F &` (an lvalue) can be
invoked with a `T const &` (a constant lvalue). It checks the exact way the library will call it. The algorithms
always call `fn` with a *named* `T const` variable, never a temporary. So a user callable that only accepts lvalues
works, and so does a member function pointer such as `&MyType::square`, through `std::invoke`.

**The concept** (math.hpp, lines 68-70):

```cpp
template<typename F, typename T>
concept callable =
    std::invocable<F &, T const &> && std::convertible_to<std::invoke_result_t<F &, T const &>, T>;
```

**The order of the constraints matters.** `exact_weights<M, N, T>` builds the weight table of method M for order N,
which is only possible for a valid method and order. C++ checks a `&&` of constraints left to right and stops at the
first false one, like `&&` at run time. Placing `exact_weights` after `method<M>` and `valid_order<N>` means the
table is only built when it can be. Otherwise an invalid M would produce a hard error inside the table generator
instead of a clean "no matching overload". The comment at lines 621-622 records this.

**`[[nodiscard]]`.** It makes the compiler warn if the result is ignored. A computed derivative that nobody reads is
almost certainly a bug.

**The body is one line.** For instance:

```cpp
return detail::derivative_with_step<M, static_cast<int>(N)>(fn, x, h);
```

It is qualified (`detail::...`) because of ADL (section B).

### `derivative_default` (1695-1703): the overload without a step

```cpp
template<typename M, int N, real T, typename F>
constexpr std::expected<Estimate<T>, error> derivative_default(F &fn, T const &x)
{
    // Checked first, so that a non-finite x is an error value in constant expressions too: its step would be
    // NaN or infinite, and Clang rejects arithmetic that produces a NaN at compile time.
    if (!math::isfinite(x)) return std::unexpected(error::non_finite);
    // The public default_step, so that derivative(fn, x) and derivative(fn, x, default_step(x)) agree.
    return detail::derivative_with_step<M, N>(fn, x, nxx::deriv::default_step<M, N>(x));
}
```

This function computes the default step and passes it on. Two details matter.

- **It checks `x` first.** The default step is computed from `|x|`, and for a NaN `x` it would be a NaN. Clang
  refuses NaN-producing arithmetic in constant evaluation (section A), so without this check
  `derivative(fn, NaN)` in a `static_assert` would be a compile error rather than an error value.
- **It calls the public `default_step`.** It doesn't compute the step privately. So `derivative(fn, x)` and
  `derivative(fn, x, default_step(x))` are the same computation by construction, a promise the documentation makes.
  A user can call `default_step` to see the step, or to try twice and half that step.

**C++ note: `std::expected` and `std::unexpected`.** `std::expected<V, E>` holds either a `V` (success) or an `E`
(failure). `return std::unexpected(error::non_finite);` creates the failure state. On the caller's side, `if (r)`
tests for success, `r->value` and `*r` read the success value, and `r.error()` reads the error. `.transform(f)`
applies `f` to the success value and passes a failure through unchanged; the library uses it in a few places.

### `derivative_with_step` (1617-1640): validate the step, then choose the engine

**The code** (deriv.hpp, lines 1617-1640):

```cpp
template<typename M, int N, real T, typename F>
constexpr std::expected<Estimate<T>, error> derivative_with_step(F &fn, T const &x, T const &h)
{
    // The step is the caller's input, often computed: validated here, before fn is called, the same way in
    // debug and release builds and in constant expressions (comparisons only, so a NaN is no arithmetic).
    if (!(h > T(0)) || !math::isfinite(h)) return std::unexpected(error::non_finite);
    if constexpr (std::same_as<M, BackwardRidders>) {
        if (!math::isfinite(x)) return std::unexpected(error::non_finite);
        // ForwardRidders on fn(-t) at -x. Negation is exact and rounding symmetric, so fn receives exactly the
        // points x - k h, and the result is that of a backward base by construction; a Backward<1> table would
        // sum its terms from x - N h and round differently. The N-th derivative changes sign for odd N.
        Mirrored<F, T> mirrored{ std::addressof(fn) };
        return detail::ridders<Forward<1>, N>(mirrored, T(-x), h).transform([](Estimate<T> e) {
            if constexpr (N % 2 == 1) e.value = T(-e.value);
            return e;
        });
    } else if constexpr (detail::extrapolation<M>) {
        if (!math::isfinite(x)) return std::unexpected(error::non_finite);
        return detail::ridders<detail::base<M>, N>(fn, x, h);
    } else {
        return detail::stencil_derivative<M, N>(fn, x, h).transform(
            [](T const &d) { return Estimate<T>{ d, std::nullopt }; });
    }
}
```

**First, the step is checked:**

```cpp
if (!(h > T(0)) || !math::isfinite(h)) return std::unexpected(error::non_finite);
```

It is written as `!(h > 0)` rather than `h <= 0` on purpose. Every comparison with a NaN is false, so `h <= 0` is
false for a NaN `h` and would let it through; `!(h > 0)` is true for a NaN and rejects it. This is a run-time check,
not an `assert`. A step is a user input, often computed, so a bad one is an expected failure rather than a
programming error, and inputs are validated at the boundary.

**Then `if constexpr` chooses the engine** (shortened from the code above):

```cpp
if constexpr (std::same_as<M, BackwardRidders>) {
    // check x; run the forward Ridders engine on the mirrored function (section 10)
} else if constexpr (detail::extrapolation<M>) {
    // check x; return detail::ridders<detail::base<M>, N>(fn, x, h);
} else {
    return detail::stencil_derivative<M, N>(fn, x, h).transform(
        [](T const &d) { return Estimate<T>{ d, std::nullopt }; });
}
```

**C++ note: `if constexpr`.** The condition is evaluated at compile time, and the branches not taken are discarded:
not compiled for this instantiation. So `derivative<Central<2>>` contains no Ridders code at all, and no run-time
test of the method happens. Each method and order is its own separately compiled function.

**The stencil branch converts the result.** `stencil_derivative` returns `std::expected<T, error>`, a plain number.
`.transform` wraps a success into an `Estimate` with no error estimate and passes a failure through unchanged.

**Why the stencils don't check `x` here.** Their engine starts with `fits()`, which rejects a non-finite `x`
anyway: `abs(NaN) <= ...` is false. The Ridders engine evaluates `fn(x)` before it calls `fits`, so it needs the
check first.

---

## 3. Compile-time gates (lines 112-219, 332-389, 597-629)

These pieces decide, while the call is compiled, whether it is allowed and whether it warns. They are all
`constexpr` functions, variable templates or concepts, and none of them runs at run time.

### The allowed orders (123-148)

```cpp
constexpr bool valid_central(int p) { return p >= 2 && p <= 10 && p % 2 == 0; }

constexpr bool valid_one_sided(int p) { return p >= 1 && p <= 6; }
```

**Why central orders are even.** A central rule's error has even powers of `h` only, so its order is always
even.

**Why the caps.** They are design decisions backed by numbers, as the comments explain. A higher order needs more
points and larger weights. The rounding error grows with the sum of the absolute weights, and beyond these caps the
gain in accuracy is small. For example, `Forward<7>` would improve the best f' error by only about 14%, and
`Forward<8>` not at all.

`valid_order<N>` (145-148) is the concept for the derivative order:

```cpp
template<auto N>
concept valid_order =
    (std::integral<decltype(N)> && !std::same_as<decltype(N), bool> && N >= 1 && N <= 4)
    || (std::same_as<decltype(N), order> && static_cast<int>(N) >= 1 && static_cast<int>(N) <= 4);
```

Any integer type except `bool` (`derivative<M, true>` would be silly), or the enum `order`, with a value from 1 to 4.
An integer is compared *in its own type*. Had it been cast to `int` first, a huge `std::size_t` value could wrap
around into 1..4.

### Which types are methods (339-389)

```cpp
template<typename S>
inline constexpr bool is_stencil = false;

template<int P>
inline constexpr bool is_stencil<Central<P>> = valid_central(P);

template<int P>
inline constexpr bool is_stencil<Forward<P>> = valid_one_sided(P);

template<int P>
inline constexpr bool is_stencil<Backward<P>> = valid_one_sided(P);

template<typename S>
concept stencil = is_stencil<S>;

template<typename M>
concept extrapolation =
    std::same_as<M, Ridders> || std::same_as<M, ForwardRidders> || std::same_as<M, BackwardRidders>;

template<typename M>
concept method = stencil<M> || extrapolation<M>;
```

**C++ note: variable templates with partial specialisation.** `is_stencil<S>` is a `bool` constant that depends on a
type. The general version says `false`. The specialisations say: "if S has the form `Central<P>` for some P, the
value is `valid_central(P)`". This is the standard way to ask "is S an instance of the template `Central`, and with
which P?". Matching the form `Central<P>` only looks at the type's *name*; it doesn't instantiate `Central<P>`'s
body, so its `static_assert` doesn't fire. A bad P simply gives `false`.

**`method` is closed on purpose.** It lists exactly six kinds of types. A user can't create a new method by writing a
new tag type: `derivative` would reject it. This was a deliberate design decision (`docs/design.md`). An open "method
protocol" was discussed and deferred until a second algorithm in the library needs it.

**`base<M>` (387-389)** maps each method to the stencil whose weight table it evaluates:

```cpp
template<typename M>
using base = std::
    conditional_t<std::same_as<M, Ridders>, Central<2>, std::conditional_t<extrapolation<M>, Forward<1>, M>>;
```

`std::conditional_t<cond, A, B>` is a compile-time "if": the type A if `cond` is true, else B. So:

- `Ridders` maps to `Central<2>`, because Ridders extrapolates the simplest central rule;
- `ForwardRidders` and `BackwardRidders` map to `Forward<1>` (`BackwardRidders` too, because it runs the forward
  engine on a mirrored function, section 10);
- a stencil maps to itself.

Two places use it: the exact-weights check (just below), and the dispatch in `derivative_with_step`.

### Can `T` hold the weights exactly? (160-169, 607-629)

**The code** (deriv.hpp, lines 160-169, 607-629):

```cpp
template<typename T>
constexpr bool holds_integers(long long v)
{
    long long power = 1;
    for (int i = 0; i < std::numeric_limits<T>::digits; ++i) {
        power *= std::numeric_limits<T>::radix;
        if (power > v) return true;
    }
    return false;
}

template<typename M, int N>
inline constexpr long long largest_weight_integer = [] {
    constexpr auto table = detail::make_table<N>(detail::base<M>{});
    long long v = table.denominator;
    for (int const c : table.numerators) v = c < -v ? -c : (c > v ? c : v);
    return v;
}();

template<typename M, auto N, typename T>
concept exact_weights = detail::holds_integers<T>(detail::largest_weight_integer<M, static_cast<int>(N)>);
```

Each rule's weights are stored as integers over a common denominator (section 4). If `T` can't represent those
integers exactly, the weights are rounded, the cancellation of the error terms is no longer exact, and the rule
silently loses its order. Every type with at least 24 significant bits, float and above, holds every table, but a
16-bit "half" type does not. Three pieces handle this:

- **`holds_integers<T>(v)` (160)** decides whether every integer from 0 to `v` is exact in `T`. A floating-point type
  with `digits` significant digits in base `radix` represents every integer below `radix^digits` exactly. The
  function multiplies `power` by the radix, `digits` times, and returns `true` as soon as `power > v`. It works in
  `long long` and stops early, so it never overflows; for `double` it would otherwise compute 2^53. Integer
  arithmetic is used because `std::numeric_limits<T>::epsilon()` isn't guaranteed to be `constexpr` for a user
  type, while `radix` and `digits` are.
- **`largest_weight_integer<M, N>` (607-613)** is the largest integer in method M's table for order N: the largest
  numerator in absolute value, or the denominator, whichever is larger.

  ```cpp
  template<typename M, int N>
  inline constexpr long long largest_weight_integer = [] {
      constexpr auto table = detail::make_table<N>(detail::base<M>{});
      long long v = table.denominator;
      for (int const c : table.numerators) v = c < -v ? -c : (c > v ? c : v);
      return v;
  }();
  ```

  **C++ note: an immediately invoked lambda.** `[] { ... }()` defines a lambda and calls it at once. It is a way to
  initialise a constant with a computation that needs a loop. Since the variable is `constexpr`, the call happens at
  compile time.

  Over all rules the largest value is 6,222,216 (`Central<10>`, f''''), under 2^24 = 16,777,216.
- **`exact_weights<M, N, T>` (628-629)** is the concept that ties them together:
  `holds_integers<T>(largest_weight_integer<M, N>)`. It is one of `derivative`'s constraints. For `float` and wider
  it always holds. For a narrower type it removes the high-order rules from overload resolution, so they don't
  compile.

### The float warning (171-217)

**The code** (deriv.hpp, lines 180-201, 216-217):

```cpp
        template<typename T>
        constexpr bool coarse()
        {
            long long power = 1;
            for (int i = 1; i < std::numeric_limits<T>::digits; ++i) {
                power *= std::numeric_limits<T>::radix;
                if (power >= 10'000'000) return false;
            }
            return true;
        }

#ifdef NUMERIXX_NO_PRECISION_WARNING
        inline constexpr bool precision_warning = false;
#else
        inline constexpr bool precision_warning = true;
#endif

template<auto N, typename T>
concept noisy = precision_warning && static_cast<int>(N) >= 2 && detail::coarse<T>();
```

**`coarse<T>()` (180-189)** is true when the type has 7 or fewer significant decimal digits, that is, when epsilon is
above 1e-7. Epsilon is `radix^(1 - digits)`, so the test is `radix^(digits - 1) < 10^7`, computed with integers for
the same reason as above. For `float`, `digits` is 24 and `2^23 = 8,388,608 < 10^7`, so it is coarse. For `double`,
`2^52` is far above, so it isn't. It also catches decimal32 and bfloat16. The comment notes that `digits10` was not
used, because it doesn't count binary and decimal types alike.

**`precision_warning` (191-201)** is a `bool` constant: `true` normally, `false` if the macro
`NUMERIXX_NO_PRECISION_WARNING` is defined before the header is included. The comment warns that the macro must be
defined the same way in every source file of a program, best on the compiler command line. Otherwise one file would
call the deprecated twin and another the plain overload for the same arguments: two different functions for the same
call, which violates the One Definition Rule.

**`noisy<N, T>` (216-217)** is the condition the twins add: `precision_warning && N >= 2 && coarse<T>()`. The comment
above it records why: the measured median errors in float, for example 4e-4 relative for `Central<2>` f'' and 1e-5
for Ridders, and noise for `ForwardRidders` f''''.

---

## 4. Weight tables: where the formulas come from (lines 391-558)

Section B derived two formulas by hand. The library needs dozens: five central orders, six forward and six backward
orders, each for four derivative orders. Rather than type them in, which invites typos, it **generates** each table
at compile time from its list of points.

### What a table holds (400-414)

```cpp
template<std::size_t Size>
struct Table
{
    std::array<int, Size> nodes;

    std::array<int, Size> numerators;

    int denominator;

    int order;
};
```

The meaning is

```
f^(N)(x)  ~=  sum_k numerators[k] * fn(x + nodes[k] * h)  /  (denominator * h^N)
```

For example, `Central<2>` for f' has the nodes {-1, 1}, the numerators {-1, 1}, the denominator 2 and the order 2:
`(-fn(x - h) + fn(x + h)) / (2h)`.

**Why integers over one denominator rather than fractions or doubles?**

- **Integers are exact** in every real type: the weight -2/3 can't be stored exactly in binary, but -8 and 12 can.
- **One denominator means one division at the end** instead of one per weight.

### `factorial` and `table_value_too_large` (422-433)

**The code** (deriv.hpp, lines 422-433):

```cpp
constexpr long long factorial(std::size_t n)
{
    long long f = 1;
    for (std::size_t i = 2; i <= n; ++i) f *= static_cast<long long>(i);
    return f;
}

inline void table_value_too_large() {}
```

`factorial(n)` is the ordinary loop; the tables need at most `4!`.

`table_value_too_large()` is an empty, ordinary (not `constexpr`) function. The generator calls it if a value
doesn't fit 24 bits.

**C++ note: the non-constexpr call trick.** During constant evaluation, calling a function that isn't `constexpr` is
not allowed. Since the generator is `consteval`, it is *always* evaluated at compile time, so reaching that call is
a compile error that points at the call. That is how the generator reports "a weight would be inexact in float" with
nothing at run time. No current table reaches it; it guards against a future change.

### `lagrange_table(nodes, n, p)` (462-506): computing the weights

**The code** (deriv.hpp, lines 462-506):

```cpp
template<std::size_t Size>
consteval Table<Size> lagrange_table(std::array<int, Size> const &nodes, int n, int p)
{
    std::array<long long, Size> num{};
    std::array<long long, Size> den{};
    long long common = 1;
    for (std::size_t k = 0; k < Size; ++k) {
        // poly = prod_{j != k} (t - n_j), coefficients from t^0 up, built one factor at a time;
        // q = prod_{j != k} (n_k - n_j), nonzero for distinct nodes.
        std::array<long long, Size> poly{};
        poly[0] = 1;
        std::size_t degree = 0;
        long long q = 1;
        for (std::size_t j = 0; j < Size; ++j) {
            if (j == k) continue;
            ++degree;
            for (std::size_t i = degree; i > 0; --i) poly[i] = poly[i - 1] - nodes[j] * poly[i];
            poly[0] = -nodes[j] * poly[0];
            q *= nodes[k] - nodes[j];
        }
        long long w = detail::factorial(static_cast<std::size_t>(n)) * poly[static_cast<std::size_t>(n)];
        long long const g = std::gcd(w, q);
        w /= g;
        q /= g;
        if (q < 0) {
            w = -w;
            q = -q;
        }
        num[k] = w;
        den[k] = q;
        common = std::lcm(common, q);
    }
    constexpr long long limit = 1LL << 24;
    if (common >= limit) detail::table_value_too_large();
    Table<Size> t{};
    for (std::size_t k = 0; k < Size; ++k) {
        long long const v = num[k] * (common / den[k]); // exact: den[k] divides common
        if (v >= limit || v <= -limit) detail::table_value_too_large();
        t.nodes[k] = nodes[k];
        t.numerators[k] = static_cast<int>(v);
    }
    t.denominator = static_cast<int>(common);
    t.order = p;
    return t;
}
```

**What it computes.** The weight of node k is the n-th derivative at 0 of the Lagrange polynomial of that node (the book
derives it):

```
w_k = n! * (coefficient of t^n in prod_{j != k} (t - n_j))  /  prod_{j != k} (n_k - n_j)
```

The code computes this exactly, in integers: the numerator polynomial, the denominator product, then the reduced
fraction.

**The code, for each node `k`:**

1. **Build the numerator polynomial.** It starts from `poly = 1` (`poly[0] = 1`, degree 0) and multiplies by
   `(t - n_j)` for every `j != k`. The coefficients are stored from `t^0` upwards. Multiplying a polynomial by
   `(t - a)` shifts every coefficient up one place (the `t` part) and subtracts `a` times the old coefficient (the
   `-a` part):

   ```cpp
   for (std::size_t i = degree; i > 0; --i) poly[i] = poly[i - 1] - nodes[j] * poly[i];
   poly[0] = -nodes[j] * poly[0];
   ```

   The loop runs from the top down so that each `poly[i - 1]` is still the old value when it is read.
2. **Build the denominator.** `q` multiplies up `prod (n_k - n_j)` in the same loop. It is never 0, because the nodes
   are distinct.
3. **Form and reduce the weight.** `w = n! * poly[n]`. The fraction `w / q` is reduced by their greatest common
   divisor, and the sign is moved to the numerator so that the denominator is positive.

After all nodes, the weights are put over a common denominator. `common` is the least common multiple of the reduced
denominators, and each numerator becomes `num[k] * (common / den[k])`. This is exact, because `den[k]` divides
`common`. Each result is checked against 2^24 (see above), and the table is filled in.

**Worked example, as in the comment at 450-452: `Central<4>` for f', nodes -2, -1, 1, 2, n = 1.** Take node `k`
= -1, so the other nodes are -2, 1 and 2.

- `poly` starts as `1`.
- Times `(t + 2)`: `t + 2`.
- Times `(t - 1)`: `t^2 + t - 2`.
- Times `(t - 2)`: `t^3 - t^2 - 4t + 4`. The coefficient of `t^1` is -4.
- `q = (-1 - (-2)) * (-1 - 1) * (-1 - 2) = 1 * (-2) * (-3) = 6`.
- `w = 1! * -4 = -4`, and `-4 / 6` reduces to `-2 / 3`.

Doing the same for the other nodes gives the weights 1/12, -2/3, 2/3, -1/12. Their lcm is 12, so the numerators are
{1, -8, 8, -1} and the denominator is 12. That is the textbook formula
`(f(x - 2h) - 8 f(x - h) + 8 f(x + h) - f(x + 2h)) / (12 h)`.

**The comment's note on sizes.** The largest intermediate value is `12!`, for `Central<10>` and f'''', which fits
easily in `long long`.

**C++ note: `consteval`.** A `consteval` function can only run at compile time; calling it with run-time arguments is
an error. That fits here: tables are always fixed by template arguments, and the run-time code only ever reads the
finished integers.

### `make_table<N>(Tag{})` (516-558): which nodes each family uses

**The code** (deriv.hpp, lines 516-558):

```cpp
template<int N, int P>
consteval Table<std::size_t{ P + N - 1 }> make_table(Central<P>)
{
    constexpr int count = P + N - 1;
    constexpr int m = count / 2; // count is even for odd N (P is even), odd for even N
    std::array<int, std::size_t{ count }> nodes{};
    std::size_t i = 0;
    for (int k = -m; k <= m; ++k)
        if (k != 0 || N % 2 == 0) nodes[i++] = k;
    return detail::lagrange_table(nodes, N, P);
}

template<int N, int P>
consteval Table<std::size_t{ P + N }> make_table(Forward<P>)
{
    std::array<int, std::size_t{ P + N }> nodes{};
    for (std::size_t i = 0; i < nodes.size(); ++i) nodes[i] = static_cast<int>(i);
    return detail::lagrange_table(nodes, N, P);
}

template<int N, int P>
consteval Table<std::size_t{ P + N }> make_table(Backward<P>)
{
    std::array<int, std::size_t{ P + N }> nodes{};
    for (std::size_t i = 0; i < nodes.size(); ++i) nodes[i] = static_cast<int>(i) - (P + N - 1);
    return detail::lagrange_table(nodes, N, P);
}
```

Three overloads choose the nodes for each family, then call `lagrange_table`.

**C++ note: tag dispatch.** A caller writes `make_table<N>(S{})`: it creates an empty object of the tag type and
passes it. Overload resolution picks the overload whose parameter type matches (`Central<P>`, `Forward<P>` or
`Backward<P>`) and deduces `P` from it. The tag object itself is never used. This is how generic code such as
`make_table<N>(detail::base<M>{})` gets the right table for any method.

**`Central<P>` (516-526)** uses `P + N - 1` nodes, symmetric around 0. For odd N, node 0 is left out (an odd
derivative gives it weight 0 anyway); for even N it is included. The code walks `k` from `-m` to `m` with
`m = count / 2` and skips 0 for odd N:

| Rule                | count | nodes        | numerators / denominator |
|---------------------|-------|--------------|--------------------------|
| `Central<2>`, f'    | 2     | -1, 1        | -1, 1 / 2                |
| `Central<2>`, f''   | 3     | -1, 0, 1     | 1, -2, 1 / 1             |
| `Central<2>`, f'''  | 4     | -2, -1, 1, 2 | -1, 2, -2, 1 / 2         |
| `Central<2>`, f'''' | 5     | -2 ... 2     | 1, -4, 6, -4, 1 / 1      |
| `Central<4>`, f'    | 4     | -2, -1, 1, 2 | 1, -8, 8, -1 / 12        |

**`Forward<P>` (535-541)** uses the nodes `0, 1, ..., P + N - 1`: `P + N` points, one more than a central rule (symmetry
saves a central rule one point), and one of
them is `x` itself. `Forward<1>` gives the classic forward differences: {-1, 1}, {1, -2, 1}, {-1, 3, -3, 1},
{1, -4, 6, -4, 1}, each over 1.

**`Backward<P>` (552-558)** uses the nodes `-(P + N - 1), ..., 0`. No special mirroring code is needed: the same
generator, given the negative nodes, produces the right weights automatically. As the comment notes, they come out
as `(-1)^N` times the forward weights in reverse.

---

## 5. Default steps (lines 560-595, 1245-1265, 1644-1680; `math.hpp:111`)

**The code** (deriv.hpp, lines 569-571, 593-595 and 1666-1680; the divisors, lines 1245-1265, are shown in section 9):

```cpp
template<real T>
constexpr T scale(T const &x)
{ return math::abs(x) < T(1) ? T(1) : math::abs(x); }

template<int Exponent, real T>
constexpr T stencil_step(T const &x)
{ return math::root_epsilon<Exponent, T>() * detail::scale(x); }

template<typename M = Central<2>, auto N = 1, real T>
requires detail::method<M> && detail::valid_order<N>
[[nodiscard]] constexpr T default_step(T const &x)
{
    constexpr int n = static_cast<int>(N);
    if constexpr (std::same_as<M, Ridders>) {
        constexpr int divisor = detail::ridders_divisors[static_cast<std::size_t>(n - 1)];
        return T(detail::scale(x) / T(divisor));
    } else if constexpr (detail::extrapolation<M>) {
        return T(detail::scale(x) / T(detail::one_sided_divisor));
    } else {
        constexpr int exponent = detail::make_table<n>(M{}).order + n;
        return detail::stencil_step<exponent>(x);
    }
}
```

These functions compute the default steps; the book explains how they were chosen.

**`scale(x)` (569-571) returns max (|x|, 1).** The default step is proportional to it, which makes the step relative to
`x` for large `|x|` (at `x = 1e6` a step of 6e-6 would be below the spacing of the doubles near `x`) and absolute
near 0 (at `x = 1e-20` a relative step would be absurdly small). It encodes the assumption that `fn` varies on a
scale of about `max(|x|, 1)`. When that is false, the user should pass a step.

**`stencil_step<P+N>(x)` (593-595) returns eps^ (1/ (P+N)) times `scale(x)`,** the step that balances truncation
against rounding when `fn` varies on that scale, without the constant factors; the comment at 579-586 records what
dropping them costs. The exponent `P + N` is a template parameter, because `root_epsilon` takes it as one. In double
at `x = 1` this gives 6.06e-6 for `Central<2>` f' (eps^ (1/3)), 2^-13 for its f'' (eps^ (1/4)), 7.4e-4 for
`Central<4>` f' (eps^ (1/5)) and 1.49e-8 for `Forward<1>` f' (eps^ (1/2)).

**`math::root_epsilon<N, T>()` (`math.hpp:111-153`) computes `eps^(1/N)`.** It can't call `std::pow`, which isn't
`constexpr` and doesn't exist for every `real` type, so it uses + - * / only, in two steps.

1. **Scale epsilon into (radix^-N, 1].** It multiplies eps by `radix^N` as long as the result stays at most 1, and
   counts the factors in `scale`. Multiplying by a power of the radix only changes the exponent, so this is exact.
   For double and N = 3, `2^-52` becomes `scaled = 0.5` after 17 multiplications by 8, with `scale = 2^17`. Then
   `eps^(1/3) = 0.5^(1/3) / 2^17`.
2. **Newton's iteration for the cube root (in general the N-th root) of `scaled`,** starting from 1. Each step is
   `y <- y + (scaled / y^(N-1) - y) / N`. It starts above the root and decreases, and stops when it no longer
   decreases: by then it has reached the root to within rounding. The scaling keeps the root between `1/radix` and
   1, so only a few iterations are needed.

The comments prove that the result is within an ulp (a unit in the last place), and explain why the correction is
written as `y + (...)/N` rather than the textbook `((N - 1) y + scaled / y^(N-1)) / N`: that form rounds a larger
sum and could stop 1.8 ulps above the root.

**The code** (math.hpp, lines 111-153):

```cpp
template<int N, real T>
constexpr T root_epsilon()
{
    static_assert(N >= 2);
    T const radix = static_cast<T>(std::numeric_limits<T>::radix);

    // Multiplying by a power of the radix changes only the exponent, so scaled and scale stay exact, and eps =
    // scaled / scale^N gives eps^(1/N) = scaled^(1/N) / scale. For double and N = 3, 2^-52 becomes scaled = 0.5
    // with scale = 2^17. Each step tries scaled * radix^N one factor of radix at a time and stops as soon as
    // the value passes 1, so no intermediate exceeds radix: forming radix^N itself would overflow for large N
    // (2^128 for float), which constant evaluation rejects. For such N no step is taken, and scaled stays eps.
    T scaled = std::numeric_limits<T>::epsilon();
    T scale = T(1);
    while (true) {
        T trial = scaled;
        for (int i = 0; i < N && trial <= T(1); ++i) trial = trial * radix;
        if (!(trial <= T(1))) break;
        scaled = trial;
        scale = scale * radix;
    }

    // Newton's iteration for y^N = scaled: y <- y + (scaled / y^(N-1) - y) / N. From y = 1, at or above the
    // root since scaled <= 1, every iterate stays above the root (by the inequality of arithmetic and geometric
    // means) and decreases; in floating point it stops decreasing once it reaches the root to within rounding.
    // The scaling put the root in (1/radix, 1], so this takes a few steps; started from 1 against eps itself,
    // it would take dozens. y^(N-1) cannot underflow: if eps was scaled, then radix^N <= 1/eps and y > 1/radix,
    // so y^(N-1) > radix^(1-N) >= radix eps; otherwise scaled = eps and y^(N-1) >= eps^((N-1)/N) >= eps.
    //
    // The update is written as a correction to y. Near the root scaled / y^(N-1) - y is exact (the operands are
    // within a factor 2), and the rounding errors of scaled / y^(N-1), up to about N/2 ulps from the N - 1
    // products, are divided by N, so the result is within an ulp: at most 0.58 ulps for float and double, N = 2
    // ... 5000, against 60-digit roots. The textbook form ((N - 1) y + scaled / y^(N-1)) / N rounds a sum of
    // size N y, worth up to an ulp of y, so its stopping test could fire 1.8 ulps above the root.
    T y = T(1);
    while (true) {
        T y_n1 = T(1);
        for (int i = 1; i < N; ++i) y_n1 = y_n1 * y;
        T const next = T(y + T(scaled / y_n1 - y) / T(N));
        if (!(next < y)) break;
        y = next;
    }
    return y / scale;
}
```

**Ridders' initial step (1255).** `ridders_divisors = {32, 4, 4, 4}`, so the initial step is `scale(x) / 32` for f'
and `scale(x) / 4` for f'' to f''''. This is only the *starting* step, since the method halves it level by level.
The comment records the prototype measurements behind these numbers. For f', `s/2` was fooled by `sin(20x)` and
`s/128` lost a digit. For higher orders a larger start pays, because the rounding error grows like `1/h^N`. The
divisors are powers of two, so `x + k h` and the halvings are as exact as they can be.

**The one-sided initial step (1265).** `one_sided_divisor = 32`, for every order, tuned in Phase 5 (the comment has
the numbers).

**The public `default_step<M, N>(x)` (1666-1680)** chooses one of these with `if constexpr`. For a stencil it reads
the accuracy order from the table, `make_table<n>(M{}).order`, so the exponent `P + N` comes from the same place as
the rule's weights and can't disagree with them.

---

## 6. Overflow-safe arithmetic (lines 631-836)

Section E explained why: in constant evaluation an overflow or a NaN is a compile error, so the code must find out
whether an operation would overflow *without performing it*. This section is a small toolkit for that. Each helper
is short, but each one exists because the obvious version failed somewhere.

**Two words used throughout.** `max` is the largest finite value of `T`, `std::numeric_limits<T>::max()`, about
1.8e308 for double. **Saturating** means "if the true result would exceed `max`, return `max` instead". The
saturating helpers are used only for *rounding floors* (error estimates), where `max` reads as "huge or unknown". The
*values* of derivatives are never saturated: they get `non_finite` instead.

### `fits(x, h, q)` (647-652): can the points and the divisor be formed?

```cpp
template<real T>
constexpr bool fits(T const &x, T const &h, int q)
{
    T const largest = (std::numeric_limits<T>::max)();
    return h <= largest / T(2 * q) && math::abs(x) <= largest - T(q) * h;
}
```

A rule evaluates `fn` at `x + k h` for nodes up to `K` in magnitude, and divides by about `denominator * h`. Both must
be finite. With `q = max(2K, denominator)` (computed by `reach()`, section 7):

- **`h <= max / (2q)`** keeps `q h`, and with it every `k h` and the divisor, at most `max / 2`;
- **`|x| <= max - q h`** keeps every point `x + k h` at least `K h` below `max`.

The comment explains why the margin is `2K` and not just `K`. `max - q h` is itself rounded, and with no margin a
point could round *up* to infinity. The check can't overflow itself: it only divides `max`, multiplies `q` by an `h`
already known to be small, and subtracts from `max`. As a bonus, it rejects a NaN or infinite `x`.

### `sum_fits(a, b)` (666-672): is `a + b` finite?

```cpp
template<real T>
constexpr bool sum_fits(T const &a, T const &b)
{
    T const radix = static_cast<T>(std::numeric_limits<T>::radix);
    T const largest = (std::numeric_limits<T>::max)();
    return math::abs(T(a / radix + b / radix)) <= T(largest / radix);
}
```

Two finite numbers can have a sum beyond `max`, for example `1e308 + 1e308`, and computing it would overflow, which
is exactly what must be avoided. The trick: divide both by the radix (2 for binary types) first. Dividing by 2 only
lowers the exponent, so it is exact, and now `a/2 + b/2` can't overflow, since each half is at most `max/2`. Then
`a + b > max` exactly when `a/2 + b/2 > max/2`. Subtraction is the same test with `-b`.

### `distance(a, b)` (682-687): `|a - b|`, saturated

**The code** (deriv.hpp, lines 682-687):

```cpp
template<real T>
constexpr T distance(T const &a, T const &b)
{
    if (!detail::sum_fits(a, T(T(0) - b))) return (std::numeric_limits<T>::max)();
    return math::abs(T(a - b));
}
```

`|a - b|` can overflow when `a` and `b` have opposite signs (`1e308 - (-1e308)`). The function checks with `sum_fits`
and returns `max` if it would overflow. Ridders uses it to measure the change between tableau entries.

### `saturating_quotient(a, b)` (698-704) and `quotient_fits(a, b)` (717-722): division

**The code** (deriv.hpp, lines 698-704, 717-722):

```cpp
template<real T>
constexpr T saturating_quotient(T const &a, T const &b)
{
    T const largest = (std::numeric_limits<T>::max)();
    if (T(a / largest) >= T(b / T(2))) return largest;
    return a / b;
}

template<real T>
constexpr bool quotient_fits(T const &a, T const &b)
{
    T const largest = (std::numeric_limits<T>::max)();
    return T(math::abs(a) / largest) < T(b / T(2));
}
```

A large value divided by a tiny step can overflow (`1e300 / 1e-300`). Both helpers decide by comparing `a / max`
with `b / 2`. Both of those divisions only make numbers smaller, so neither can overflow.

- `saturating_quotient(a, b)` returns `max` when `a / b` would be too large, otherwise `a / b`. It is used for
  floors.
- `quotient_fits(a, b)` just answers yes or no, conservatively: it also says no for quotients between `max / 2` and
  `max`. It is used for the intermediate quotients of the second difference.

### `saturating_sum(a, b)` (732-737)

**The code** (deriv.hpp, lines 732-737):

```cpp
template<real T>
constexpr T saturating_sum(T const &a, T const &b)
{
    if (!detail::sum_fits(a, b)) return (std::numeric_limits<T>::max)();
    return T(a + b);
}
```

`a + b` for non-negative `a` and `b`, saturated, using `sum_fits`.

### `product_fits(a, b)` (754-763) and `saturating_product(a, b)` (773-778): multiplication

**The code** (deriv.hpp, lines 754-763, 773-778):

```cpp
template<real T>
constexpr bool product_fits(T const &a, T const &b)
{
    T const radix = static_cast<T>(std::numeric_limits<T>::radix);
    T const largest = (std::numeric_limits<T>::max)();
    T const q = largest / b;
    if (a < q) return true;
    if (a > q) return false;
    return T(T(q / radix) * b) <= T(largest / radix);
}

template<real T>
constexpr T saturating_product(T const &a, T const &b)
{
    if (b > T(1) && !detail::product_fits(a, b)) return (std::numeric_limits<T>::max)();
    return T(a * b);
}
```

Is `a * b <= max`, for `b >= 1`? The natural test is to compare `a` with `q = max / b`. But `q` is rounded, and the
comment records two earlier versions that were each wrong in one direction:

- `a >= q` (reject) wrongly rejected `max * 1`;
- `a > q` (reject) wrongly accepted `6 * fl(max / 6)`, which rounds up to infinity in double.

The current version reasons as follows:

- If `a < q`, the product fits, even if `q` was rounded up: the value just below `q` is at most the exact `max / b`.
- If `a > q`, it doesn't fit.
- Only `a == q` is undecided. Then it computes `(q / radix) * b` and compares it with `max / radix`. The value is
  about `max / radix`, so it can't overflow, and the comparison is exact up to the product's own rounding.

`saturating_product` uses this to return `max` or the product.

### `value_rounding(f)` (791-798): how wrong is one value of `fn`?

```cpp
template<real T>
constexpr T value_rounding(T const &f)
{
    T const e = std::numeric_limits<T>::epsilon() * math::abs(f);
    T const tiny = T(std::numeric_limits<T>::denorm_min());
    if (f != T(0) && e < tiny) return tiny;
    return e;
}
```

This is the basic unit of every rounding floor: a value `f` of `fn` is assumed correct to about `eps |f|`. There is
one refinement. Below the normal range, around 1e-308 in double, numbers are "subnormal" and their spacing no longer
shrinks with `|f|`: it stays at `denorm_min`, about 5e-324. So the error of a nonzero value is at least
`denorm_min`. The comment records the case that prompted this. Without it, a `fn` with subnormal values got a floor
of exactly 0, and Ridders reported a result that was 192 `denorm_min` off with an error estimate of 0.

### `divide_by_powers<N, Saturating>(a, dh, h)` (819-836): dividing by `h^N` safely

**The code** (deriv.hpp, lines 819-836):

```cpp
template<int N, bool Saturating, real T>
constexpr T divide_by_powers(T const &a, T const &dh, T const &h)
{
    auto const divide = [](T const &p, T const &q) -> T {
        if constexpr (Saturating)
            return detail::saturating_quotient(p, q);
        else
            return T(p / q);
    };
    T v = a;
    if (dh > T(1) && h < T(1)) {
        for (int i = 1; i < N; ++i) v = divide(v, h);
        return divide(v, dh);
    }
    v = divide(v, dh);
    for (int i = 1; i < N; ++i) v = divide(v, h);
    return v;
}
```

Every rule ends by dividing a weighted sum by `denominator * h^N`. The obvious code, `a / (denominator * pow(h, N))`,
has a problem: `h^N` can underflow to 0. For instance, `h^4` is 0 in double for `h` below about 1e-77, while the
quotient itself is perfectly representable. So the function never forms `h^N`. With `dh = denominator * h` it
computes `a / dh`, then divides by `h` another `N - 1` times.

**The order of the divisions matters too.** If `dh > 1` and `h < 1`, dividing by `dh` first makes the value smaller
before the divisions by `h` make it larger again, and that first quotient can underflow to 0. The comment records a
real case: `Central<10>` for f'''' with a tiny step returned 0 for a true result of 2.6e-308. So in that case the
divisions by `h` come first and `dh` last. The intermediates then exceed the result by at most the factor `dh`,
which can only overflow spuriously for a result within that factor of `max`, and that is reported as `non_finite`
rather than as a silent wrong 0.

The template flag `Saturating` chooses the division: `saturating_quotient` for floors, plain `/` for values. A value
is checked with `isfinite` afterwards.

---

## 7. Building blocks: differences and their rounding (lines 838-1143)

These functions compute the parts that the two engines (sections 8 and 9) are made of. Each one evaluates `fn`, or
takes values already evaluated, checks everything before computing, and returns either a result or `non_finite`.
Several of them also return what Ridders needs to estimate the rounding error.

### First differences (845-927)

**The code** (deriv.hpp, lines 845-859, 871-876, 889-900, 917-927):

```cpp
template<real T>
struct Difference
{
    T value;

    T f_lo;

    T f_hi;

    T dx;
};

template<real T>
constexpr T rounding_floor(Difference<T> const &d)
{
    return detail::saturating_quotient(
        T(detail::value_rounding(d.f_lo) + detail::value_rounding(d.f_hi)), d.dx);
}

template<real T>
constexpr std::expected<Difference<T>, error> difference_of(T const &f_lo, T const &f_hi, T const &dx)
{
    // Checked before subtracting, so neither inf - inf nor an overflowing difference is formed: constant
    // evaluation rejects arithmetic that produces a NaN (and on GCC an overflow), so the failure would be a
    // compile error instead of an error value.
    if (!math::isfinite(f_lo) || !math::isfinite(f_hi)) return std::unexpected(error::non_finite);
    if (!detail::sum_fits(f_hi, T(T(0) - f_lo))) return std::unexpected(error::non_finite);
    T const d = (f_hi - f_lo) / dx;
    if (!math::isfinite(d)) return std::unexpected(error::non_finite); // the quotient overflowed
    return Difference<T>{ d, f_lo, f_hi, dx };
}

template<real T, typename F>
constexpr std::expected<Difference<T>, error> difference(F &fn, T const &x_lo, T const &x_hi)
{
    T const dx = x_hi - x_lo;
    if (dx == T(0)) return std::unexpected(error::non_finite); // h is below the resolution of x
    // fn may return an expression template, and every intermediate is a T const, never auto: an unevaluated
    // expression would hold references to its operands, which can dangle.
    T const f_lo = static_cast<T>(std::invoke(fn, x_lo));
    T const f_hi = static_cast<T>(std::invoke(fn, x_hi));
    return detail::difference_of(f_lo, f_hi, dx);
}
```

**`Difference<T>` (845-859)** holds the result of `(f_hi - f_lo) / dx` and the three inputs it came from. The
stencils only use `.value`; Ridders also needs `f_lo`, `f_hi` and `dx` for its rounding floor.

**`difference_of(f_lo, f_hi, dx)` (889-900)** does the arithmetic, in a carefully checked order:

1. both values must be finite: `fn` may have returned inf or NaN;
2. `f_hi - f_lo` must fit (`sum_fits`): two huge values of opposite sign;
3. only now form `(f_hi - f_lo) / dx`;
4. the quotient must be finite: a tiny `dx`.

The comment at 892-894 explains why the order matters: if the check came after the subtraction, `inf - inf` would
already have produced a NaN, which in constant evaluation is a compile error rather than an error value. It is a
separate function from `difference` because the forward f' level of Ridders already has `fn(x)` and needs only this
arithmetic.

**`difference(fn, x_lo, x_hi)` (917-927)** evaluates `fn` at both points and calls `difference_of`. Two choices here
are deliberate.

- **It divides by the distance between the points as actually evaluated, `x_hi - x_lo`, not by the nominal `2h`.**
  The computer can't represent most values of `x + h`; it stores the nearest double. So the point `fn` actually
  receives is slightly off, and the true spacing of the two points differs from `2h`. The comment gives the size: at
  `x = 1`, `h = 1e-5`, the spacing can differ from `2h` by about 1e-11 relative, which is as large as the error the
  rule is aiming for. Dividing by the true spacing makes this error disappear for a two-point rule: the
  rule is then the exact slope of the line through the two points actually evaluated.
- **The values are stored as `T const f_lo = static_cast<T>(std::invoke(fn, x_lo));`.**
    - `std::invoke` calls `fn` whether it is a function, a lambda or a member function pointer.
  - `static_cast<T>` converts the result, which may be an expression template (section B), into a real `T` at once.
    - The comment at 922-923 makes the `auto` point: an unevaluated expression holding references can dangle.

If the two points round to the same number (`dx == 0`), the step is below the resolution of `x`, and the function
returns `non_finite` before calling `fn`.

**`rounding_floor(d)` (871-876)** estimates how much of `d.value` is rounding noise: each value is off by about
`eps |f|` (`value_rounding`), and the difference divides that by `dx`, so

```
floor = (eps |f_lo| + eps |f_hi|) / dx
```

with saturation in case `dx` is tiny. Ridders uses this "floor" so that it never claims an error smaller than
rounding can explain.

### Second differences (934-1040)

**The code** (deriv.hpp, lines 934-951, 973-993, 1009-1018, 1033-1040):

```cpp
template<real T>
struct SecondDifference
{
    T value;

    T f_lo;

    T f_hi;

    T h_lo;

    T h_hi;
};

template<real T>
constexpr std::expected<SecondDifference<T>, error>
    second_difference_of(T const &f_lo, T const &f_mid, T const &f_hi, T const &h_lo, T const &h_hi)
{
    // Every difference is checked before it is formed, as in difference(): two first differences that
    // overflow with the same sign would otherwise give inf - inf.
    if (!math::isfinite(f_lo) || !math::isfinite(f_mid) || !math::isfinite(f_hi))
        return std::unexpected(error::non_finite);
    if (!detail::sum_fits(f_hi, T(T(0) - f_mid)) || !detail::sum_fits(f_mid, T(T(0) - f_lo)))
        return std::unexpected(error::non_finite);
    T const up = T(f_hi - f_mid);
    T const down = T(f_mid - f_lo);
    if (!detail::quotient_fits(up, h_hi) || !detail::quotient_fits(down, h_lo))
        return std::unexpected(error::non_finite);
    T const d_hi = up / h_hi;
    T const d_lo = down / h_lo;
    if (!detail::sum_fits(d_hi, T(T(0) - d_lo))) return std::unexpected(error::non_finite);
    T const d = T(d_hi - d_lo) / T(T(h_hi + h_lo) / T(2));
    if (!math::isfinite(d)) return std::unexpected(error::non_finite); // the quotient overflowed
    return SecondDifference<T>{ d, f_lo, f_hi, h_lo, h_hi };
}

template<real T, typename F>
constexpr std::expected<SecondDifference<T>, error>
    second_difference(F &fn, T const &x, T const &f0, T const &x_lo, T const &x_hi)
{
    T const h_lo = x - x_lo;
    T const h_hi = x_hi - x;
    T const f_lo = static_cast<T>(std::invoke(fn, x_lo));
    T const f_hi = static_cast<T>(std::invoke(fn, x_hi));
    return detail::second_difference_of(f_lo, f0, f_hi, h_lo, h_hi);
}

template<real T>
constexpr T second_floor(SecondDifference<T> const &d, T const &f_mid)
{
    T const e0 = detail::value_rounding(f_mid);
    T const hi = detail::saturating_quotient(T(detail::value_rounding(d.f_hi) + e0), d.h_hi);
    T const lo = detail::saturating_quotient(T(e0 + detail::value_rounding(d.f_lo)), d.h_lo);
    return detail::saturating_quotient(detail::saturating_sum(hi, lo), T(T(d.h_hi + d.h_lo) / T(2)));
}
```

**`second_difference_of(f_lo, f_mid, f_hi, h_lo, h_hi)` (973-993)** computes f'' from three points that need not be
evenly spaced:

```
( (f_hi - f_mid) / h_hi  -  (f_mid - f_lo) / h_lo )  /  ( (h_hi + h_lo) / 2 )
```

**What the formula means.** The two inner quotients are the slopes on the right and on the left. Their difference,
divided by the distance between the midpoints of the two intervals, is the rate of change of the slope, which is f''.
With `h_lo = h_hi = h` it becomes the textbook `(f_hi - 2 f_mid + f_lo) / h^2`.

**Why allow uneven offsets?** For the same reason as above: the evaluated points `x - h` and `x + h` are rounded, so
the real offsets `h_lo` and `h_hi` usually differ slightly from `h` and from each other. Using the real offsets makes
the formula exact for quadratics, whatever the rounding of the points.

**Why the offsets are exact.** By Sterbenz's lemma, the difference of two floating-point numbers within a factor of 2
of each other is exact. For `x_hi = fl(x + h)` with `0 < h <= |x|`, `x_hi` and `x` are that close, so
`h_hi = x_hi - x` has no rounding. The comment states the same conditions.

**The checks before every step.** Two first differences that each overflow with the same sign would give
`inf - inf`. So the function checks:

1. all three values finite;
2. both differences fit (`sum_fits`);
3. both inner quotients fit (`quotient_fits`);
4. their difference fits.

Only then does it divide by the half-sum of the offsets, and it checks that result too. `h^2` is never formed: in
double it underflows for `h` below about 1.5e-154.

**`second_difference(fn, x, f0, x_lo, x_hi)` (1009-1018)** computes the real offsets `x - x_lo` and `x_hi - x`,
evaluates the two outer points, and takes `f0 = fn(x)` from the caller, who has it already. The comment to the
previous function records that in the prototype, using the nominal step instead of the real offsets made the worst
case for float 8.7 times the error estimate, against 1.4.

**`second_floor(d, f_mid)` (1033-1040)** propagates each value's rounding through the same formula, so it is the
rounding floor of a second difference:

```
( (eps|f_hi| + eps|f_mid|) / h_hi  +  (eps|f_mid| + eps|f_lo|) / h_lo )  /  ( (h_hi + h_lo) / 2 )
```

It saturates at every step, because each step can overflow on its own. The comment gives one example for each case:
`fn = 1e300` with `h = 1e-300` overflows the inner division, while `fn = 1` with `h = 1e-170` overflows only the last.

### Weighted sums (1047-1113): the general case

**The code** (deriv.hpp, lines 1047-1061, 1084-1113):

```cpp
template<real T>
struct WeightedSum
{
    T value;

    T magnitude;

    T f_minus;

    T f_plus;
};

template<bool Magnitude, std::size_t Size, real T, typename F>
constexpr std::expected<WeightedSum<T>, error>
    weighted_sum(F &fn, Table<Size> const &t, T const &x, T const &h, std::optional<T> const &f0)
{
    T sum = T(0);
    T magnitude = T(0);
    T f_minus = T(0);
    T f_plus = T(0);
    for (std::size_t k = 0; k < Size; ++k) {
        int const node = t.nodes[k];
        int const c = t.numerators[k];
        T fk = T(0);
        if (node == 0 && f0) {
            fk = *f0;
        } else {
            T const xk = node == 0 ? x : T(x + T(node) * h); // a T const lvalue, as for two nodes
            fk = static_cast<T>(std::invoke(fn, xk));
        }
        if (!math::isfinite(fk)) return std::unexpected(error::non_finite);
        if (node == -1) f_minus = fk;
        if (node == 1) f_plus = fk;
        T const size = T(c < 0 ? -c : c);
        if (c != 0 && !detail::product_fits(math::abs(fk), size)) return std::unexpected(error::non_finite);
        T const term = T(T(c) * fk);
        if (!detail::sum_fits(sum, term)) return std::unexpected(error::non_finite);
        sum = T(sum + term);
        if constexpr (Magnitude) magnitude = T(magnitude + size * detail::value_rounding(fk));
    }
    return WeightedSum<T>{ sum, magnitude, f_minus, f_plus };
}
```

Rules with more than three points use their table directly.

**`WeightedSum<T>` (1047-1061)** holds four things:

- the sum itself;
- its "magnitude", the sum of each term's rounding, for Ridders' floor;
- the values at nodes -1 and +1, if the rule has them. Ridders uses these to estimate `f'` (section 9).

**`weighted_sum<Magnitude>(fn, table, x, h, f0)` (1084-1113)** loops over the nodes. For each:

1. **Get the value `fk`.** If the node is 0 and the caller passed `f0` (an `std::optional`), use it: Ridders
   evaluates `fn(x)` once and reuses it for every level. Otherwise form the point `x + node * h` as a named `T const`
   (`x` itself for node 0) and evaluate `fn` there.
2. **Check `fk`.** It must be finite.
3. **Remember the neighbours.** Record `fk` if the node is -1 or +1.
4. **Add the term.** Check that `|fk| * |c|` fits (`product_fits`), form the term `c * fk`, check that the partial sum
   plus the term fits (`sum_fits`), and add.
5. **Accumulate the rounding, if requested.** If `Magnitude` is true, add `|c| * eps|fk|` to the magnitude.

The doc comment records the case that led to this level of checking. A constant function `fn = 1e308` under
`Central<4>` (numerators 1, -8, 8, -1) produced the term `-8e308`, which is `-inf`, and then `+8e308`, which is `+inf`,
and `-inf + inf` is NaN, a compile error in a `static_assert`. With the checks, it is `non_finite`.

`Magnitude` is a template parameter, not a run-time flag, so the stencils, which pass `false`, don't pay for the
extra arithmetic: `if constexpr` removes it.

### `largest_node` and `reach` (1122-1143)

**The code** (deriv.hpp, lines 1122-1143):

```cpp
template<std::size_t Size>
constexpr int largest_node(Table<Size> const &t)
{
    int const lo = t.nodes[0];
    int const hi = t.nodes[Size - 1];
    return -lo > hi ? -lo : hi;
}

template<std::size_t Size>
constexpr int reach(Table<Size> const &t)
{
    int const k = 2 * detail::largest_node(t);
    return k > t.denominator ? k : t.denominator;
}
```

- **`largest_node(t)`** returns K, the largest node in absolute value. Every point the rule uses lies within `K h` of
  `x`.
- **`reach(t)`** returns `q = max(2K, denominator)`, the value `fits()` needs.

For `Central<2>` that is 2, 2, 4, 4 for f' to f''''; for `Forward<1>` it is `2N`.

---

## 8. `stencil_derivative`: the fixed-step methods (lines 1145-1219)

**The code** (deriv.hpp, lines 1162-1219), explained below:

```cpp
template<typename S, int N, real T, typename F>
constexpr std::expected<T, error> stencil_derivative(F &fn, T const &x, T const &h)
{
    constexpr auto table = detail::make_table<N>(S{});
    constexpr std::size_t size = table.nodes.size();
    constexpr int lo = table.nodes[0];
    constexpr int hi = table.nodes[size - 1];
    if (!detail::fits(x, h, detail::reach(table))) return std::unexpected(error::non_finite);

    if constexpr (size == 2) {
        // With two nodes (f' only), the weights are -1 and +1 over hi - lo, so the difference is the whole
        // formula.
        auto const point = [&x, &h](int k) -> T {
            if (k == 0) return x;
            return T(x + T(k) * h);
        };
        auto const d = detail::difference(fn, point(lo), point(hi));
        if (!d) return std::unexpected(d.error());
        return d->value;
    } else {
        // h is below the resolution of x. Both sides are checked: at a power of the radix the spacing differs
        // on the two sides of x, so h can move x one way and not the other (at x = -1, h = 0x1.8p-54 moves x
        // up but not down), which would leave Backward and Central rules with coinciding points and a wrong
        // result.
        if (T(x + h) == x || T(x - h) == x) return std::unexpected(error::non_finite);
        if constexpr (N == 2 && std::same_as<S, Central<2>>) {
            // The three-point f'' from the actual offsets, as the two-node f' rules use the actual distance.
            T const f0 = static_cast<T>(std::invoke(fn, x));
            if (!math::isfinite(f0)) return std::unexpected(error::non_finite);
            auto const d = detail::second_difference(fn, x, f0, T(x - h), T(x + h));
            if (!d) return std::unexpected(d.error());
            return d->value;
        } else {
            // With more nodes, dividing by the outer-node distance would cancel the outer-node rounding only
            // for two-node rules and enlarge it for the others (4x for Central<4>), so the divisor is
            // denominator * h, then h N - 1 more times, never h^N (h^4 underflows below 1.2e-77 in double):
            // with dyadic x and a power-of-two h the result is then exact for polynomials of degree up to
            // N + order - 1. The divisions are ordered so that they neither underflow nor overflow
            // spuriously (divide_by_powers). Node rounding adds about
            // S / 2 * eps * |x f'| / h^N, the same order as fn's rounding. Neighbouring points must differ:
            // at h near 0.6 ulp(x), x + h and x + 2h round to the same value, and the rule would return a
            // finite but wrong value.
            for (std::size_t k = 1; k < size; ++k) {
                int const a = table.nodes[k - 1];
                int const b = table.nodes[k];
                T const x_a = a == 0 ? x : T(x + T(a) * h);
                T const x_b = b == 0 ? x : T(x + T(b) * h);
                if (!(x_a < x_b)) return std::unexpected(error::non_finite);
            }
            std::optional<T> const none;
            auto const s = detail::weighted_sum<false>(fn, table, x, h, none);
            if (!s) return std::unexpected(s.error());
            T const d = detail::divide_by_powers<N, false>(s->value, T(T(table.denominator) * h), h);
            if (!math::isfinite(d)) return std::unexpected(error::non_finite); // the quotient overflowed
            return d;
        }
    }
}
```

This is the engine for `Central<P>`, `Forward<P>` and `Backward<P>`. It applies one rule at one step.

The comment at 1148-1151 summarises why the weighted sum, divided by `denominator * h^N`, is the derivative (the
moment conditions; see the book).

The code begins:

```cpp
constexpr auto table = detail::make_table<N>(S{});            // built at compile time
if (!detail::fits(x, h, detail::reach(table))) return std::unexpected(error::non_finite);
```

The table is a `constexpr` local, so it is computed once by the compiler, and the run-time code just reads its
integers. `fits` rules out overflow of the points and, as noted, a non-finite `x`.

Then come three branches, chosen at compile time by the size of the table and the rule.

**Branch 1: two nodes (`size == 2`).** This covers f' by `Central<2>`, `Forward<1>` and `Backward<1>`. Two points
always have weights -1 and +1 over the distance between them, so the rule is just a difference, and `difference()`
does everything: the true spacing, the checks, and `non_finite` if the two points coincide. The small lambda
`point(k)` returns `x` for node 0 (so `Forward<1>` evaluates at `x` itself, not at `x + 0 * h`) and `x + k h`
otherwise.

**Branch 2: first, a check for steps below the resolution of x.**

```cpp
if (T(x + h) == x || T(x - h) == x) return std::unexpected(error::non_finite);
```

If `h` is smaller than about half the spacing of the doubles near `x`, then `x + h` rounds back to `x`, and the rule
would combine values at coinciding points. The comment explains why both sides are tested. Just below a power of two
the spacing of the doubles is different on the two sides of `x` (it halves below the power). So a step can move `x`
in one direction and not the other: at `x = -1`, a step of `0x1.8p-54` moves `x` up but not down. Testing only one
side would let a `Backward` or `Central` rule run with two points that coincide.

**Branch 2a: `Central<2>` for f''** uses `second_difference` with the real offsets, as in section 7. It evaluates
`fn(x)` itself first.

**Branch 2b: every other rule.**

1. **Neighbouring points must differ.** The loop compares each pair of consecutive points. At a step around 0.6 ulp
   of `x`, `x + h` and `x + 2h` can round to the same number. The check above passes in that case, but the rule
   would combine two identical points and return a finite, wrong value.
2. **The weighted sum:** `weighted_sum<false>`, with no rounding magnitude and no cached `fn(x)`.
3. **Divide** by `denominator * h`, then by `h` another `N - 1` times (`divide_by_powers`), and check that the result
   is finite.

**Why the nominal step `h` here, when the two-point rules use the true spacing?** The comment at 1195-1203 explains
it. With more than two points, the rounding errors of the different points don't share one distance, so dividing by
the outer-node distance doesn't cancel them; it *enlarges* them, four times for `Central<4>`. With the nominal `h`,
a power-of-two step and an `x` that is exactly representable, every point is exact, and the rule is exact for
polynomials of low enough degree. The tests rely on that: their `static_assert`s compare such results with exact
values. The rounding of the points adds about `S / 2 * eps * |x f'| / h^N`, the same order as `fn`'s own rounding.

---

## 9. Ridders' extrapolation (lines 1221-1575)

This is the most complex part of the file. A short summary of the tableau the code builds comes first, then the
pieces, then the main function from top to bottom, then three worked examples.

### What the code implements

The engine builds Richardson's tableau. Row `i` is a **level**: the base rule `D` at the step `h_start / 2^i`, and
the entries extrapolated from it:

```
A[0][0]
A[1][0]  A[1][1]
A[2][0]  A[2][1]  A[2][2]
A[3][0]  A[3][1]  A[3][2]  A[3][3]

A[i][0] = D(h_start / 2^i)
A[i][j] = A[i][j-1] + (A[i][j-1] - A[i-1][j-1]) / (r - 1),     r = 4^j (central) or 2^j (one-sided)
```

Each entry needs only its left neighbour, in the same row, and the entry diagonally up-left, its **diagonal parent**,
in the previous row. That is why the code keeps just two rows. The last entry of a row, `A[i][i]`, is its **diagonal
entry**, the most extrapolated one. The factor is `4^j` for the central base, whose error has even powers
of `h` only, and `2^j` for the forward base, whose error has every power.

The method estimates the error of every entry, keeps the entry with the smallest estimate, and stops when further
levels only add rounding. Each halving of `h` multiplies the base rule's rounding by `2^N`.

### The constants (1221-1265)

**The code** (deriv.hpp, lines 1229-1234, 1245-1265):

```cpp
template<std::size_t N, real T>
constexpr std::array<T, N> filled(T const &v)
{
    return [&v]<std::size_t... I>(std::index_sequence<I...>) { return std::array<T, N>{ ((void)I, v)... }; }(
               std::make_index_sequence<N>{});
}

inline constexpr int ridders_levels = 32;

inline constexpr std::array<int, 4> ridders_divisors{ 32, 4, 4, 4 };

inline constexpr int one_sided_divisor = 32;
```

- **`filled<N>(v)` (1229-1234)** returns a `std::array<T, N>` with every element equal to `v`.

  **C++ note: why not `std::array<T, 32> a{};`?** That needs `T` to be default-constructible, which `nxx::real` does
  not require. `filled` builds the array by repeating `v` through a parameter pack: `((void)I, v)...` expands to
  `v, v, v, ...`, once per index in `std::make_index_sequence<N>`. The `(void)I,` part only uses `I` to make the
  expansion happen.
- **`ridders_levels = 32` (1245)** is the total number of levels the method may compute, *including restarts*. It
  caps the cost, but it is only a guard: the comment notes that a typical f' result takes 5 to 11 levels. With 2
  evaluations per level for central f', that caps the cost at 64. The other limits are in the comment.
- **`ridders_divisors`, `one_sided_divisor`** are the initial steps from section 5.

### `Level<T>` and `ridders_level<S, N>` (1272-1361): one entry of column 0

**The code** (deriv.hpp, lines 1272-1280, 1310-1361):

```cpp
template<real T>
struct Level
{
    T value;

    T floor;
};

template<typename S, int N, real T, typename F>
constexpr std::expected<Level<T>, error> ridders_level(F &fn, T const &x, T const &f0, T const &h)
{
    constexpr bool central = std::same_as<S, Central<2>>;
    if constexpr (N == 1 && central) {
        auto const d = detail::difference(fn, T(x - h), T(x + h));
        if (!d) return std::unexpected(d.error());
        return Level<T>{ d->value, detail::rounding_floor(*d) };
    } else if constexpr (N == 1) {
        T const x_1 = T(x + h);
        T const dx = T(x_1 - x);
        T const f_1 = static_cast<T>(std::invoke(fn, x_1));
        auto const d = detail::difference_of(f0, f_1, dx);
        if (!d) return std::unexpected(d.error());
        return Level<T>{ d->value, detail::rounding_floor(*d) };
    } else if constexpr (N == 2 && central) {
        auto const d = detail::second_difference(fn, x, f0, T(x - h), T(x + h));
        if (!d) return std::unexpected(d.error());
        return Level<T>{ d->value, detail::second_floor(*d, f0) };
    } else if constexpr (N == 2) {
        T const x_1 = T(x + h);
        T const x_2 = T(x + T(2) * h);
        T const f_1 = static_cast<T>(std::invoke(fn, x_1));
        T const f_2 = static_cast<T>(std::invoke(fn, x_2));
        auto const d = detail::second_difference_of(f0, f_1, f_2, T(x_1 - x), T(x_2 - x_1));
        if (!d) return std::unexpected(d.error());
        return Level<T>{ d->value, detail::second_floor(*d, f_1) };
    } else {
        constexpr auto table = detail::make_table<N>(S{});
        constexpr int reach = detail::largest_node(table);
        std::optional<T> const cached = f0;
        auto const s = detail::weighted_sum<true>(fn, table, x, h, cached);
        if (!s) return std::unexpected(s.error());
        constexpr int weights = [](auto const &t) {
            int w = 0;
            for (int const c : t.numerators) w += c < 0 ? -c : c;
            return w;
        }(table);
        T const eps = std::numeric_limits<T>::epsilon();
        T const slope = central
                            ? detail::saturating_quotient(detail::distance(s->f_plus, s->f_minus), T(T(2) * h))
                            : detail::saturating_quotient(detail::distance(s->f_plus, f0), h);
        T const spread = T(T(eps / T(2)) * T(math::abs(x) + T(reach) * h)) * T(weights);
        T const divisor = T(table.denominator) * h;
        T const value = detail::divide_by_powers<N, false>(s->value, divisor, h);
        T const rounding = detail::divide_by_powers<N, true>(s->magnitude, divisor, h);
        T const points =
            detail::divide_by_powers<N, true>(detail::saturating_product(spread, slope), divisor, h);
        if (!math::isfinite(value)) return std::unexpected(error::non_finite); // the quotient overflowed
        return Level<T>{ value, detail::saturating_sum(rounding, points) };
    }
}
```

A `Level` is the base rule's estimate at one step (`value`) and its rounding floor (`floor`). `ridders_level`
computes it. The base `S` is `Central<2>` or `Forward<1>`, and there are five cases, all chosen at compile time:

| Case             | Rule                                                                                | Rounding floor               |
|------------------|-------------------------------------------------------------------------------------|------------------------------|
| N = 1, central   | `difference(fn, x - h, x + h)`                                                      | `rounding_floor`             |
| N = 1, forward   | `difference_of(f0, fn(x + h), (x + h) - x)`, reusing `fn(x)`                        | `rounding_floor`             |
| N = 2, central   | `second_difference(fn, x, f0, x - h, x + h)`                                        | `second_floor(d, f0)`        |
| N = 2, forward   | `second_difference_of(f0, fn(x + h), fn(x + 2h), ...)`: the middle point is `x + h` | `second_floor(d, fn(x + h))` |
| N = 3 or 4, both | `weighted_sum<true>` over S's table, with the nominal step                          | two parts, below             |

For f' and f'' it reuses the building blocks of section 7, with the true spacing and offsets. The comment records
what the nominal step cost in the prototype: the worst case was 22 times the error estimate instead of 3.6 for
central f', and 8.7 instead of 1.4 for float f''.

**The floor for f''' and f''''** has two parts.

- **`rounding`: `fn`'s own rounding.** `sum_k |w_k| eps |f_k|` (the magnitude from `weighted_sum<true>`), divided by
  `denominator * h^N`.
- **`points`: the rounding of the points.** These rules use the nominal step, so each point `x + k h` is off by up to
  half an ulp: about `eps/2 * (|x| + K h)`, with K the largest node. That moves `fn`'s value there by about that much
  times the slope `f'`. Summed with the weights, it is

  ```
  eps/2 * (|x| + K h) * |f'| * sum_k |w_k|   /  (denominator * h^N)
  ```

  The slope `f'` is estimated from values the sum already has: `(f(x+h) - f(x-h)) / 2h` centrally, `(f(x+h) - f(x))
  / h` forward. That is what `f_minus` and `f_plus` in `WeightedSum` are for.

Every step saturates, and the two parts are added with `saturating_sum`. The comment gives the measured effect of
the second part: the prototype's central f''' of `exp(10x)` on [1, 3] was covered by the estimate (estimate at least
the actual error) in 96.3% of the cases without it, and 99.0% with it.

**A naming detail.** Inside this branch, `constexpr int reach = detail::largest_node(table);` (1339) is K, not the
`q` returned by the function `detail::reach()`. The local name shadows the function's name, harmlessly, since the
function is always called with its `detail::` qualification. `weights` is the sum of absolute numerators, computed by
a small `constexpr` lambda.

### `propagated_floor(f_new, f_old, r)` (1376-1383): the floor of an extrapolated entry

An entry is `A[i][j] = (r a - b) / (r - 1)`, where `a = A[i][j-1]` and `b = A[i-1][j-1]` (the difference form,
rearranged). If `a` carries a rounding error of up to `fa` and `b` up to `fb`, then `A[i][j]` carries up to
`(r fa + fb) / (r - 1)`. The function computes that, saturating each step:

```cpp
template<real T>
constexpr T propagated_floor(T const &f_new, T const &f_old, T const &factor)
{
    T const largest = (std::numeric_limits<T>::max)();
    T const sum = detail::saturating_sum(detail::saturating_product(f_new, factor), f_old);
    if (sum == largest) return largest;
    return T(sum / T(factor - T(1)));
}
```

The middle line matters. If the sum saturated to `max`, dividing it by `r - 1` would give a large but finite number
that *looks* like a real floor. Keeping `max` keeps the meaning "unknown".

The comment notes why this is exact rather than just a bound: unrolled, each entry is a combination of levels whose
coefficients alternate in sign, so the formula gives exactly `sum |coefficient_k| * floor_k`.

**Who uses it.** Only the one-sided methods (see "Find the best entry" below). Their columns divide by `2^j - 1`,
which is small (1, 3, 7, ...), so rounding is amplified about 3.5 times by the time it reaches a useful entry. In the
central tableau, dividing by `4^j - 1` (3, 15, 63, ...) amplifies it only about 1.4 times.

### `ridders<S, N>(fn, x, h)` (1420-1575), from top to bottom

**The code** (deriv.hpp, lines 1420-1575). The whole function comes first; the steps below go through it piece by piece,
and their fragments are shortened, with checks summarised as comments:

```cpp
template<typename S, int N, real T, typename F>
constexpr std::expected<Estimate<T>, error> ridders(F &fn, T const &x, T const &h)
{
    constexpr auto table = detail::make_table<N>(S{});
    constexpr std::size_t size = table.nodes.size();
    // The points of a level reach x +- K h (K = 1 or 2 centrally, N forward).
    constexpr int q = detail::reach(table);
    constexpr bool one_sided = table.nodes[0] == 0;
    constexpr int ratio = one_sided ? 2 : 4;
    // The nodes and 0 are the consecutive integers lo ... hi for both bases: these points must differ.
    constexpr int lo = table.nodes[0] < 0 ? table.nodes[0] : 0;
    constexpr int hi = table.nodes[size - 1] > 0 ? table.nodes[size - 1] : 0;
    constexpr bool needs_f0 = [](auto const &t) {
        for (int const k : t.nodes)
            if (k == 0) return true;
        return false;
    }(table);
    // fn(x) is part of every level where the table has node 0 (even N centrally, always forward): evaluated
    // once, and kept across restarts. If it is not finite, no step avoids it.
    T const f0 = needs_f0 ? static_cast<T>(std::invoke(fn, x)) : T(0);
    if (!math::isfinite(f0)) return std::unexpected(error::non_finite);
    auto const point = [&x](T const &step, int k) -> T {
        if (k == 0) return x;
        return T(x + T(k) * step);
    };

    std::array<T, ridders_levels> prev = detail::filled<ridders_levels>(T(0));
    std::array<T, ridders_levels> row = prev;
    // The floors of the entries, carried through the tableau; used by the one-sided methods only.
    std::array<T, ridders_levels> floor_prev = prev;
    std::array<T, ridders_levels> floor_row = prev;
    std::optional<Estimate<T>> fallback; // the last untrusted result, if the levels run out
    T h_start = h;
    int budget = ridders_levels;

    while (budget > 0) {
        // One tableau, from h_start. The best entry and the stall counter belong to it.
        bool have_best = false;
        bool restart = false;
        bool resolution = false;
        T best = T(0);
        T best_error = T(0);
        T best_floor = T(0);
        T diagonal = T(0);
        int stall = 0;
        T h_i = h_start;
        for (std::size_t i = 0; budget > 0; ++i, h_i = T(h_i / T(2))) {
            // Checked before evaluating, and after the overflow check, so that the points are only formed when
            // they are finite: below the resolution of x the points coincide, and a tiny h could underflow to
            // 0. Neighbouring points must differ too (at h near 0.6 ulp, x + h and x + 2h round to the same
            // value).
            bool const fit = detail::fits(x, h_i, q);
            bool coincide = false;
            if (fit)
                for (int k = lo; k < hi && !coincide; ++k) coincide = point(h_i, k) == point(h_i, k + 1);
            if (!(h_i > T(0)) || coincide) {
                resolution = true;
                break;
            }
            --budget;
            std::expected<Level<T>, error> const level =
                fit ? detail::ridders_level<S, N>(fn, x, f0, h_i)
                    : std::expected<Level<T>, error>(std::unexpect, error::non_finite);
            bool finite = level.has_value();
            if (finite) {
                row[0] = level->value;
                floor_row[0] = level->floor;
                T factor = T(1);
                // Both the difference and the sum are checked before they are formed: finite entries of
                // opposite sign can differ by more than max.
                for (std::size_t j = 1; j <= i && finite; ++j) {
                    factor = T(factor * T(ratio));
                    finite = detail::sum_fits(row[j - 1], T(T(0) - prev[j - 1]));
                    if (!finite) break;
                    T const correction = T((row[j - 1] - prev[j - 1]) / (factor - T(1)));
                    finite = detail::sum_fits(row[j - 1], correction);
                    if (!finite) break;
                    row[j] = T(row[j - 1] + correction);
                    floor_row[j] = detail::propagated_floor(floor_row[j - 1], floor_prev[j - 1], factor);
                }
            }
            // fn is not finite within reach of x (or an entry would exceed max): every earlier level's
            // interval contained the bad point, so the whole tableau goes, and a new one starts below the
            // failing step.
            if (!finite) {
                restart = true;
                h_start = T(h_i / T(2));
                break;
            }

            // The estimate of an extrapolated entry is its change from its diagonal parent (the textbook
            // successive-approximation estimate), but never less than the rounding floor: the minimum of noisy
            // differences is biased low. Without the floor, the estimate was below the actual error in 58% of
            // the prototype's random cases. An entry equal to its parent ends the tableau: nothing is left to
            // learn (a derivative that is exactly 0 then costs 4 to 6 evaluations instead of 64).
            bool improved = false;
            bool exact = false;
            for (std::size_t j = 1; j <= i; ++j) {
                T const entry_floor = one_sided ? floor_row[j] : level->floor;
                T const change = detail::distance(row[j], prev[j - 1]);
                if (change == T(0)) {
                    best = row[j];
                    best_error = entry_floor;
                    best_floor = entry_floor;
                    have_best = exact = true;
                    break;
                }
                T const estimate = change < entry_floor ? entry_floor : change;
                if (!have_best || estimate < best_error) {
                    best = row[j];
                    best_error = estimate;
                    best_floor = entry_floor;
                    have_best = improved = true;
                }
            }
            diagonal = row[i];
            std::swap(prev, row);
            std::swap(floor_prev, floor_row);
            if (exact) break;
            // Each halving multiplies the rounding error of D_i by 2^N, so once two levels improve nothing,
            // more levels only add rounding. One level is not enough: in the prototype a Runge-function case
            // then stopped early with an error of 3.7e-5.
            if (i > 0) {
                stall = improved ? 0 : stall + 1;
                if (stall == 2) break;
            }
        }
        if (restart) continue;
        if (!have_best) break; // the resolution of x or the level cap came before two successful levels

        // The best entry must also agree with the most recent diagonal entry, which catches one that only
        // looked good (aliasing). This cut the prototype's uncovered cases from 3.9% to under 1%.
        T const disagreement = detail::distance(best, diagonal);
        T const abs_error = best_error < disagreement ? disagreement : best_error;
        Estimate<T> const result{ best, abs_error };

        // At most two correct bits, and more error than rounding explains: fn misbehaves within h_start without
        // producing a non-finite value (a pole near x, aliasing, or, for a one-sided base, fn varying on a
        // scale far below h_start), so start again. A derivative that is 0 has an estimate near the floor and
        // is not affected. In the prototype, the factor 1 instead of 4 left 11 of 400 poles 1/(x - a)
        // uncovered for f' (32 for f''), 4 and 16 none; 16 to 256 for the floor gave identical results.
        // Written with divisions, which cannot overflow: the floor, and so the estimate, saturates at max.
        bool const untrusted = abs_error >= T(math::abs(best) / T(4)) && T(abs_error / T(64)) > best_floor;
        if (!untrusted || resolution) return result;
        fallback = result;
        // The new tableau starts below the smallest step this one tried, h_i (the loop ended on a stall or
        // exact agreement, before h_i was halved), as after a non-finite level: restarting at h_start / 2
        // re-evaluated the same large steps, and the levels ran out before the step reached fn's scale. With
        // it, central Ridders missed poles 1e-3 and 1e-4 from x (f'' to f'''' off by 1 to 11, uncovered) and
        // the one-sided methods failed at sqrt and log at 1e-6 and 1e-8. A case the slower descent handled can
        // lose up to 8 times in accuracy; restarts overlapping one or two levels did worse on every measure.
        h_start = T(h_i / T(2));
    }
    if (fallback) return *fallback;
    return std::unexpected(error::non_finite);
}
```

The function runs one tableau after another (restarts), each one level after another, until it can return a
result it trusts or runs out of levels.

#### Setup (1423-1453)

Facts fixed at compile time, from S's table:

- **`q = reach(table)`**, for `fits`.
- **`one_sided = (table.nodes[0] == 0)`.** The forward table starts at node 0; the central one at -1 or -2.
- **`ratio`:** 2 for one-sided, 4 for central. It is the base of `r = ratio^j`.
- **`lo` and `hi`:** the lowest and highest node, with 0 included. All the points from `x + lo h` to `x + hi h`
  must be distinct at every level; that is checked below.
- **`needs_f0`:** whether the table has node 0. That is even N centrally, and always for forward.

The value at `x` comes first:

```cpp
T const f0 = needs_f0 ? static_cast<T>(std::invoke(fn, x)) : T(0);
if (!math::isfinite(f0)) return std::unexpected(error::non_finite);
```

If node 0 is used, `fn(x)` is evaluated **once** here and reused by every level and every restart: it doesn't depend
on `h`. If it isn't finite, no step can help, so the method fails at once.

Then the working storage:

- **`prev`, `row`:** the previous and current rows of the tableau. Each new entry needs only its left neighbour (same
  row) and its diagonal parent (previous row), so two rows suffice, swapped after each level.
- **`floor_prev`, `floor_row`:** the same for the propagated floors (used by the one-sided methods).
- **`fallback`:** the last result that wasn't trusted. It is returned if the levels run out.
- **`h_start`:** the initial step of the current tableau.
- **`budget`:** the levels left, starting at 32.

#### The outer loop: one tableau per pass (1455-1572)

It runs while budget remains. At the start of each pass it resets the per-tableau state:

- `have_best`, `best`, `best_error`, `best_floor`: the best entry so far and its estimate;
- `diagonal`: the latest diagonal entry;
- `stall`: how many levels in a row brought no improvement;
- the flags `restart` and `resolution`.

#### The inner loop: one level per pass (1466-1546)

```cpp
for (std::size_t i = 0; budget > 0; ++i, h_i = T(h_i / T(2))) {
```

Level `i` uses the step `h_i = h_start / 2^i`. The halving happens at the *end* of each pass, which turns out to
matter for restarts.

**Step 1: can this step be used at all? (1471-1478)**

```cpp
bool const fit = detail::fits(x, h_i, q);
bool coincide = false;
if (fit)
    for (int k = lo; k < hi && !coincide; ++k) coincide = point(h_i, k) == point(h_i, k + 1);
if (!(h_i > T(0)) || coincide) { resolution = true; break; }
```

The step must be positive: after many halvings it could underflow to 0. And no two neighbouring points may
coincide. When either fails, we've reached the **resolution** of `x`: the step can't get any smaller and remain
meaningful. The flag `resolution` records that, and the tableau ends. The points are only formed when `fit` holds,
so forming them can't overflow.

**Step 2: compute the level (1479-1482)**

```cpp
--budget;
level = fit ? detail::ridders_level<S, N>(fn, x, f0, h_i) : (non_finite);
```

The level costs one unit of budget. A step too large to fit (only possible for `x` near `max`, or a huge `h`)
counts as a non-finite level without calling `fn`. It is then handled like any other non-finite level, by
restarting with a smaller step.

**Step 3: extend the row (1483-1500)**

```cpp
row[0] = level->value;  floor_row[0] = level->floor;
T factor = T(1);
for (j = 1; j <= i; ++j) {
    factor = factor * ratio;                              // r = ratio^j
    // check that row[j-1] - prev[j-1] fits
    T const correction = (row[j-1] - prev[j-1]) / (factor - 1);
    // check that row[j-1] + correction fits
    row[j] = row[j-1] + correction;
    floor_row[j] = propagated_floor(floor_row[j-1], floor_prev[j-1], factor);
}
```

This is exactly the tableau formula from "What the code implements". Level `i` has `i + 1` entries, `row[0]` to
`row[i]`. Both the
difference and the sum are checked first: entries of opposite sign can be far apart.

**Step 4: a non-finite level restarts the tableau (1504-1508)**

```cpp
if (!finite) { restart = true; h_start = T(h_i / T(2)); break; }
```

Typically `fn` returned NaN or inf at one of the points, as when the points of `log` reach below 0. Why discard the
*whole* tableau, not just this level? The levels go from large steps to small ones, and every earlier level spanned
a *wider* interval around `x` than this one. If `fn` is bad somewhere within this level's interval, it was bad
within theirs too: their values came from a function that isn't smooth there, and an extrapolation built on them
means nothing. So a fresh tableau starts below the failing step, at `h_i / 2`.

**Step 5: find the best entry (1515-1534)**

```cpp
for (j = 1; j <= i; ++j) {
    T const entry_floor = one_sided ? floor_row[j] : level->floor;
    T const change = detail::distance(row[j], prev[j - 1]);
    if (change == T(0)) { best = row[j]; best_error = best_floor = entry_floor; exact = true; break; }
    T const estimate = change < entry_floor ? entry_floor : change;
    if (!have_best || estimate < best_error) { best = row[j]; best_error = estimate; best_floor = entry_floor;
                                               improved = true; }
}
```

For each new extrapolated entry, the method estimates its error and remembers the entry with the smallest estimate
seen in this tableau.

- **The estimate is the change from the diagonal parent**, `|A[i][j] - A[i-1][j-1]|`. It is the textbook
  a-posteriori estimate, and pessimistic: the parent is the less accurate of the two.
- **The estimate is never below the floor.** The rounding noise in an entry is about its floor, so no estimate
  smaller than that can be trusted. And the method picks the *smallest* of many noisy estimates, which biases the
  choice towards estimates that are too low by chance. The comment gives the measurement: without the floor, the
  estimate was below the actual error in 58% of the prototype's random cases.
- **Which floor.** `entry_floor` is the propagated floor for the one-sided methods, and the newest level's own floor
  for central `Ridders`, as tuned in Phase 3.
- **Exact agreement ends the tableau.** If an entry equals its parent exactly (`change == 0`), there is nothing more
  to learn. This happens for polynomials, and for a derivative that is exactly 0. The comment notes that such a
  derivative then costs 4 to 6 evaluations instead of 64.

**Step 6: bookkeeping (1535-1538)**

```cpp
diagonal = row[i];
std::swap(prev, row);  std::swap(floor_prev, floor_row);
if (exact) break;
```

`row[i]`, the last entry of the row, is the newest diagonal entry, the most extrapolated one. It is kept for the
check after the tableau. The swap makes the current row the "previous" one for the next level; the old previous row
becomes scratch space.

**Step 7: the stall rule (1542-1545)**

```cpp
if (i > 0) { stall = improved ? 0 : stall + 1; if (stall == 2) break; }
```

If two levels in a row bring no better estimate, rounding has started to win, and more levels will only make it
worse: each halving multiplies the rounding of `D` by `2^N`. So the tableau ends. Why two and not one? The comment
records a Runge-function case that stopped too early, with an error of 3.7e-5, when one level was allowed.

#### After a tableau (1547-1572)

```cpp
if (restart) continue;          // a non-finite level: start the next tableau (h_start was set in step 4)
if (!have_best) break;          // not even two levels: give up
```

`have_best` needs at least two levels, because the first level has no extrapolated entry. It can be false when the
resolution of `x` was reached almost at once, or the budget ran out.

**The diagonal check (1552-1554):**

```cpp
T const disagreement = detail::distance(best, diagonal);
T const abs_error = best_error < disagreement ? disagreement : best_error;
```

The best entry must also agree with the most recent diagonal entry. An entry can look good by accident: two
neighbouring entries happen to agree while both are wrong. This is aliasing, which is typical for oscillating
functions. Requiring agreement with the latest, most extrapolated entry catches most such cases; the comment says it
cut the prototype's uncovered cases from 3.9% to under 1%.

**The untrusted test (1562):**

```cpp
bool const untrusted = abs_error >= T(math::abs(best) / T(4)) && T(abs_error / T(64)) > best_floor;
```

A result is **untrusted** when both of these hold:

1. **`abs_error >= |best| / 4`.** The estimated error is at least a quarter of the value, so the result has at most
   two correct bits by its own estimate.
2. **`abs_error / 64 > best_floor`.** The estimated error is more than 64 times what rounding explains. So it isn't
   rounding that went wrong; `fn` must be doing something within the tableau's steps that the method can't follow.
   That could be a pole close to `x`, aliasing, or, for a one-sided method, `fn` varying on a scale far below
   `h_start`.

Condition 2 protects a derivative that is truly 0 or tiny. Its estimate is about its floor, so it fails condition
2 and is accepted even though condition 1 holds. The tests are written with divisions rather than multiplications
(`abs_error / 64`, not `64 * best_floor`) because the floor may be `max` and a multiplication could overflow. The
comment records how the factors 4 and 64 were chosen.

**Return, or restart (1563-1571):**

```cpp
if (!untrusted || resolution) return result;
fallback = result;
h_start = T(h_i / T(2));
```

A trusted result is returned. An untrusted one is also returned if the tableau ended at the resolution of `x`: no
smaller step is possible, so there is nothing to retry.

Otherwise the result is kept as a `fallback`, and a new tableau starts **below the smallest step this tableau
tried**. Here the detail from the top of the inner loop matters: the loop ended with `break` (on a stall or exact
agreement), and `break` skips the `h_i` halving in the `for` header. So `h_i` still holds the smallest step that was
actually evaluated, and `h_i / 2` is the next one down.

This "fast descent" was added in Phase 5, for all three Ridders methods. Before it, the restart was at
`h_start / 2`, which re-evaluated almost the same large steps again and ran out of budget before reaching the scale
of a pole close to `x`. The comment lists what that broke (central Ridders missed poles 1e-3 and 1e-4 from `x`; the
one-sided methods failed for `sqrt` and `log` at 1e-6 and 1e-8) and what was tried instead.

**The end (1573-1574):** when the budget is spent, the last untrusted result is better than nothing, so it is
returned if there is one. Otherwise the call fails with `non_finite`.

### Worked example 1: a smooth function, f' of exp at x = 1 from h = 1/32

The table shows the error `A[i][j] - e` of each entry, from the prototype trace; the library's result is identical.

**`Ridders` (central, r = 4^j):**

```
i  h      floor    A[i][0]   A[i][1]   A[i][2]   A[i][3]   A[i][4]   A[i][5]   best so far
0  2^-5   1.9e-14  +4.4e-4                                                   -
1  2^-6   3.9e-14  +1.1e-4   -5.4e-9                                         A[1][1], est 4.4e-4
2  2^-7   7.7e-14  +2.8e-5   -3.4e-10  +2.8e-14                              A[2][2], est 5.4e-9
3  2^-8   1.5e-13  +6.9e-6   -2.1e-11  -2.2e-15  -2.7e-15                    A[3][3], est 1.5e-13 (the floor)
4  2^-9   3.1e-13  +1.7e-6   -1.3e-12  -1.7e-14  -1.8e-14  -1.8e-14          no improvement (stall 1)
5  2^-10  6.2e-13  +4.3e-7   -2.4e-13  -1.6e-13  -1.7e-13  -1.7e-13  -1.7e-13 no improvement (stall 2): stop
```

Read it like this.

- **Column 0 is the plain central difference.** Its error falls by 4 with each halving of `h`, as `h^2` predicts.
  Even at the smallest step it is only 4e-7.
- **Each column to the right removes one more error term.** At level 1, combining two errors of 4.4e-4 and 1.1e-4 (4
  times smaller) cancels almost all of it: 5.4e-9. By `A[3][3]` the error is 2.7e-15, about the accuracy of the
  value of `e` itself.
- **The estimates are pessimistic early on.** At level 1 the estimate is 4.4e-4: the change from `A[0][0]`, whose
  error it really measures. The actual error of `A[1][1]` is 5.4e-9. Each new level's estimate is roughly the
  previous best entry's actual error.
- **The floor doubles at each level**, because for f' rounding grows like `1/h`. At level 3 the change,
  `|A[3][3] - A[2][2]| = 3e-14`, is below the floor of 1.5e-13, so the estimate becomes the floor. The method
  refuses to claim more accuracy than rounding allows.
- **Then rounding wins.** At levels 4 and 5 the entries get worse (1.8e-14, then 1.7e-13), and no estimate improves.
  After two such levels the tableau stops.
- **The diagonal check:** `|A[3][3] - A[5][5]|` is about 1.6e-13, slightly above the best estimate, so
  `abs_error = 1.6e-13`.
- **The untrusted test:** 1.6e-13 is far below `|best| / 4`, so the result is trusted.

The result is 2.7182818284590424 (actual error 3e-15) with `abs_error` 1.6e-13. It cost 12 evaluations: 6 levels
of 2.

**`ForwardRidders`, same function and step (r = 2^j):**

```
i  h      A[i][0]   ...   best entry's error   estimate
0  2^-5   +4.3e-2
...
5  2^-10  +1.3e-3   ...   A[5][5]: +1.9e-13     6.5e-12 (the propagated floor)
6, 7      no improvement: stop
```

The first column is much worse: first order, 4.3e-2 at the first step, against 4.4e-4 centrally. Each column
removes only one power of `h`, so the best entry, `A[5][5]` with an error of 1.9e-13, is reached after 6 levels
instead of 4.

- **The estimate is the propagated floor.** At that entry it is 6.5e-12, five times the level's own floor of 1.2e-12,
  because columns divided by `2^j - 1` amplify rounding, as described under `propagated_floor`.
- **The result** is 2.7182818284592369, with `abs_error` 6.5e-12. It costs 9 evaluations: `fn(1)` once, then 8
  levels of 1.

### Worked example 2: a non-finite value, f' of log at x = 0.01 (the default step)

This is the example in the doc comment of `Ridders`. The default initial step is `max(0.01, 1) / 32 = 0.03125`, so
the first level evaluates `log` at `0.01 - 0.03125 = -0.02125`, where it is NaN. What happens, level by level:

```
tableau 1:  h = 1/32    x - h = -0.02125   log is NaN  -> non-finite level: restart at h/2
tableau 2:  h = 1/64    x - h = -0.005625  log is NaN  -> non-finite level: restart at h/2
tableau 3:  h = 1/128   x - h = 0.0021875  fine; D = 134.2
            h = 1/256                            D = 105.6
            ...                                  (the column-0 errors fall by 4 per level)
            h = 2^-16 (1.5e-5)                   D = 100.00008
            two levels without improvement: stop
```

Each failed level still cost two evaluations: `difference()` evaluates both points before it checks them. The
result is 100.00000000000504 (exact 100, relative error 5e-14), with `abs_error` 1.8e-11, for 24 evaluations: 2
failed levels and 10 good ones, 2 evaluations each.

Note that the first column is far from 100 at first (134 at h = 1/128, where the step is almost as large as `x`), yet
the extrapolated entries converge.

### Worked example 3: an untrusted result, a pole close to x

`fn(t) = 1 / (t - 1.001)` at `x = 1`. The pole is at a distance of 0.001, far less than the default step of 1/32.
The exact derivative is `-1 / (1 - 1.001)^2`, about -1,000,000.

```
tableau 1:  h = 1/32     D = +1,025       (the points straddle the pole: x + h = 1.031 > 1.001)
            h = 1/64     D = +4,113
            h = 1/128    D = +16,657
            h = 1/256    D = +70,132      two levels without improvement: stop
            best = +5,142, abs_error = 88,400
            untrusted: 88,400 >= |5,142| / 4 = 1,286, and 88,400 / 64 = 1,381 > floor (9e-13)
            -> fallback = this result; restart at h_i / 2 = 1/512
tableau 2:  h = 1/512    D = +355,278     (x + h = 1.00195, still beyond the pole)
            h = 1/1024   D = -21,586,298  (now below the pole, but very close to it)
            h = 1/2048   D = -1,313,057
            ...
            h = 2^-20    D = -1,000,000.9
            stop; best = -1,000,000.0000002214, abs_error = 5.8e-8: trusted
```

Every value of `fn` was finite: the points never hit the pole exactly. So no non-finite level occurred, and without
the untrusted test the method would have returned `+5,142`, wrong in sign and size, from the first tableau. The
untrusted test spotted it: an error estimate larger than the value, and nothing to do with rounding. The restart
went straight to the next smaller step, 1/512, rather than back to 1/64.

The final result has a relative error of 1.2e-15 and cost 32 evaluations: 4 levels in the first tableau and 12 in
the second, 2 evaluations each. Because the restart continues the halving where the first tableau stopped, the
steps form one unbroken sequence from 1/32 down to 2^-20. Only the tableau around them was rebuilt.

---

## 10. `BackwardRidders` and `Mirrored` (lines 1577-1603, 1623-1632)

`BackwardRidders` must evaluate `fn` only at `x` and below. The trick: reflect the function. Let `g(t) = fn(-t)`.
Then the points above `-x` for `g` are the points below `x` for `fn`, and by the chain rule

```
g^(N)(-x) = (-1)^N * fn^(N)(x)
```

So the code runs the forward engine on `g` at `-x`, and flips the sign for odd N (shortened from the first branch of
`derivative_with_step`, shown in full in section 2):

```cpp
Mirrored<F, T> mirrored{ std::addressof(fn) };              // g(t) = fn(-t)
return detail::ridders<Forward<1>, N>(mirrored, T(-x), h)
       .transform([](Estimate<T> e) { if constexpr (N % 2 == 1) e.value = T(-e.value); return e; });
```

`abs_error` is a magnitude and stays as it is.

**Why mirror instead of writing a backward engine?** The comment at 1625-1627 gives the reasons.

- **The points are exactly right.** The forward engine evaluates `g` at `fl(-x + k h)`, so `fn` receives
  `-fl(-x + k h)`. Negation is exact, and IEEE rounding is symmetric (`fl(-a) = -fl(a)`), so that equals
  `fl(x - k h)`: precisely the point a backward rule would use.
- **All the code and tuning are shared.** The forward engine, its floors and its stopping rules are reused as they
  are, so the two methods have the same accuracy, cost and limits by construction.
- **A `Backward<1>` table would differ.** Its nodes would be -N to 0, so it would sum its terms starting from
  `x - N h`, which rounds differently. And the engine's one-sided logic tests "the first node is 0".

**`Mirrored<F, T>` (1587-1603)** is the function `g`:

```cpp
template<typename F, real T>
struct Mirrored
{
    F *fn;

    constexpr T operator()(T const &t) const
    {
        T const point = -t;
        return static_cast<T>(std::invoke(*fn, point));
    }
};
```

It holds a pointer to the user's `fn` and, when called with `t`, calls `fn(-t)`, with `-t` stored as a named
`T const` lvalue as everywhere else.

**Why a struct with a pointer, and not the obvious lambda `[&fn](T const &t) { return fn(-t); }`?** When the user
passes a free function, `F` is a *function type*, and the lambda captures a reference to a function. MSVC 19.51
rejects that in constant evaluation with "a read of an uninitialized symbol", while GCC and Clang accept it. It was
reproduced in isolation in Phase 5. With a pointer, `F *` is a function pointer for a function and an object pointer
for a lambda, and it works on all compilers. `std::addressof(fn)` is used instead of `&fn` in case a user's callable
type overloads `operator&`.

---

## 11. `derivative_of`: a derivative as a callable (lines 1826-1911)

**The code** (deriv.hpp, lines 1840-1870, 1908-1911):

```cpp
template<typename M, int N, typename G>
struct DerivativeOf
{
    G fn;

    template<typename H>
    constexpr DerivativeOf(std::in_place_t, H &&f) : fn(std::forward<H>(f))
    {}

    template<typename Self, real T>
    requires nxx::detail::callable<decltype((std::declval<std::remove_reference_t<Self> &>().fn)), T>
             && detail::exact_weights<M, N, T>
    constexpr std::expected<Estimate<T>, error> operator()(this Self &&self, T const &x)
    { return detail::derivative_default<M, N>(self.fn, x); }
};

template<typename M = Central<2>, auto N = 1, typename F>
requires detail::method<M> && detail::valid_order<N> && std::constructible_from<std::decay_t<F>, F>
[[nodiscard]] constexpr auto derivative_of(F &&fn)
{ return detail::DerivativeOf<M, static_cast<int>(N), std::decay_t<F>>{ std::in_place, std::forward<F>(fn) }; }
```

Some algorithms take a function as an argument. Newton's method for roots needs `f'`, for instance. `derivative_of`
packages "the derivative of `fn`" as such a function:

```cpp
auto d = derivative_of<Ridders, 2>(fn);   // d(x) gives the same result as derivative<Ridders, 2>(fn, x)
auto r = d(1.0);                          // a std::expected<Estimate<double>, error>
```

**`derivative_of(F &&fn)` (1908-1911)** returns a `DerivativeOf<M, N, std::decay_t<F>>` holding a *copy* of `fn`.

- **`std::decay_t<F>`** is "the type you'd get by passing `fn` by value": a reference is removed, and a function type
  becomes a function pointer. So a lambda is copied (or moved, if passed as a temporary), and a free function is
  stored as a pointer.
- **Why a copy.** It means `d` stays valid when the original `fn` goes out of scope. The doc comment explains how to
  share instead of copying: pass `std::ref(fn)`, and keep `fn` alive.
- **The constraint.** It is only `std::constructible_from<std::decay_t<F>, F>`. At this point the real type `T` isn't
  known (it is decided at each call `d(x)`), so the callable and exactness checks can't be made yet.

**`DerivativeOf<M, N, G>` (1840-1870)** has one member, `G fn`, and two functions.

**The constructor:**

```cpp
template<typename H>
constexpr DerivativeOf(std::in_place_t, H &&f) : fn(std::forward<H>(f)) {}
```

- **Why a constructor at all.** Without one, `DerivativeOf` would be an aggregate and `fn` would be *copy-initialised*,
  which refuses a type whose copy or move constructor is `explicit`. The constraint on
  `derivative_of` promises that any `G` constructible from `F` works, so the member must be *direct-initialised*,
  `fn(std::forward<H>(f))`. The Codex review of Phase 3 found this gap.
- **The `std::in_place` tag.** A constructor template taking just `H &&` would also match when someone copies a
  non-const `DerivativeOf`, and it would win over the copy constructor. Requiring the tag as the first argument rules
  that out.

**The call operator:**

```cpp
template<typename Self, real T>
requires nxx::detail::callable<decltype((std::declval<std::remove_reference_t<Self> &>().fn)), T>
         && detail::exact_weights<M, N, T>
constexpr std::expected<Estimate<T>, error> operator()(this Self &&self, T const &x)
{ return detail::derivative_default<M, N>(self.fn, x); }
```

**C++ note: the explicit object parameter ("deducing this", C++23).** `this Self &&self` makes the object the call
was made on into an ordinary, deduced parameter. If `d` is a non-const object, `Self` is `DerivativeOf &` and
`self.fn` is non-const; if `d` is const, `self.fn` is const. One function covers both, and `fn` is called with the
same constness as `d`. That matters for a `mutable` lambda, which can only be called when non-const: it works
through a non-const `d` and is rejected through a const one, the same rule as `std::bind_front`.

**The constraint.** It spells out `derivative`'s conditions directly: `fn` callable with `T`, and exact weights. It
doesn't just ask "would `derivative<M, N>(self.fn, x)` compile?". Asking that would run overload resolution on
`derivative`, which for float f'' would select the deprecated twin and warn.

**The body.** It calls `detail::derivative_default` directly, which is why **calls through `derivative_of` never
warn**. The doc comment explains why that is intended: the call `d(x)` usually happens inside some algorithm in a
header, where the user couldn't put a pragma, so a warning there could not be silenced.

---

## 12. `std::formatter<nxx::deriv::error>` (lines 1915-1951)

**The code** (deriv.hpp, lines 1916-1951):

```cpp
namespace std {

    template<>
    struct formatter<nxx::deriv::error> : formatter<string_view>
    {
        template<typename FormatContext>
        auto format(nxx::deriv::error e, FormatContext &ctx) const
        {
            // No default: a new enumerator makes the compilers warn here (-Wswitch) until its name is added.
            string_view name = "unknown";
            switch (e) {
            case nxx::deriv::error::non_finite:
                name = "non_finite";
                break;
            }
            return formatter<string_view>::format(name, ctx);
        }
    };

} // namespace std
```

This lets you print an error: `std::println("no estimate: {}", r.error())` prints `no estimate: non_finite`.

- **It inherits from `std::formatter<std::string_view>`**, so it gets the string's parsing of format options for
  free: `{:>12}` right-aligns the name in 12 characters.
- **Its `format` maps the enumerator to its name** with a `switch` and hands the name to the string formatter.
- **The `switch` has no `default:` case, on purpose.** If someone adds an enumerator to `error` and forgets this
  function, compilers warn (`-Wswitch`) that the new value isn't handled, and with warnings as errors the build
  stops. A `default:` would silence that. The initial `"unknown"` only covers a value outside the enumeration,
  produced by a cast.

---

## 13. Habits you will see everywhere, and why

| Habit                                                | Why                                                                                                                                       |
|------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------|
| `T(1)`, `T(2)`, never `1.0`                          | a `double` literal would silently reduce a 50-digit `T` to 53 bits (section B)                                                            |
| `T const v = T(a + b)`, never `auto`                 | Boost.Multiprecision returns an expression template from `a + b`; `auto` would store an unevaluated object holding references (section B) |
| `static_cast<T>(std::invoke(fn, p))`                 | `std::invoke` handles member pointers like `&Real::square`; the cast evaluates an expression-template result at once                      |
| points stored as named `T const` before calling `fn` | `fn` always receives a `T const` lvalue, the exact form `callable` checks, so a callable that rejects temporaries works                   |
| `(std::numeric_limits<T>::max)()`                    | the parentheses stop a `max` macro (from `<windows.h>`) from expanding                                                                    |
| `detail::...`, `math::...` written out               | ADL can't pick a same-named function from the namespace of a user's `T` or `F` (section B)                                                |
| check, *then* compute                                | in constant evaluation an overflow or NaN is a compile error, not a value (section A)                                                     |
| `!(h > T(0))` instead of `h <= T(0)`                 | every comparison with NaN is false, so this form rejects a NaN as well                                                                    |
| `if constexpr` dispatch, no virtual functions        | each method and order is a separate function, with no run-time choice                                                                     |
| `constexpr auto table = make_table<N>(S{})`          | each table is built by the compiler; the run-time code only reads integers                                                                |
| saturate floors, never values                        | a floor of `max` reads as "unknown"; a value that would overflow is `non_finite`                                                          |

---

## 14. Check yourself

Each answer is hidden under its question; CLion's Markdown preview shows them as expandable sections.

<details><summary>1. Why does <code>derivative&lt;Central&lt;2&gt;, 2&gt;(fn, 1.0f)</code> warn, but <code>derivative_of&lt;Central&lt;2&gt;, 2&gt;(fn)(1.0f)</code> doesn't?</summary>

The first call matches both the plain overload and its deprecated twin. The twin's constraint adds `noisy<2, float>`
to the plain one, so it subsumes it and wins, and calling a `[[deprecated]]` function warns at the call.
`DerivativeOf::operator()` calls `detail::derivative_default` directly, so it never goes through overload resolution
of `derivative`.
</details>

<details><summary>2. How many points does <code>Central&lt;4&gt;</code> use for f''? For f'''?</summary>

P + N - 1. That is 5 for f'' (nodes -2 ... 2, with 0) and 6 for f''' (-3 ... 3, without 0).
</details>

<details><summary>3. Why does the two-point rule divide by <code>x_hi - x_lo</code>, while <code>Central&lt;4&gt;</code> divides by <code>12 h</code>?</summary>

For two points, dividing by their real distance makes the rule the exact slope through the points actually
evaluated, so the rounding of the points cancels. With more points, the points' rounding errors don't share a single
distance, and dividing by the outer distance would amplify them (4 times for `Central<4>`). With the nominal step
and a power-of-two `h`, the rule stays exact for polynomials.
</details>

<details><summary>4. In <code>ridders</code>, a level is non-finite. Why throw away the whole tableau and not just that level?</summary>

The bad point is within this level's interval around `x`, and every earlier level used a larger step, so its
interval contained the bad point too. Their values came from a function that isn't smooth there, and the
extrapolation built on them is meaningless. The next tableau starts at `h_i / 2`.
</details>

<details><summary>5. Why does the one-sided engine carry floors through the tableau (<code>propagated_floor</code>), while the central one uses the newest level's floor?</summary>

One-sided columns divide by 2^j - 1 (1, 3, 7, ...), which amplifies the newest level's rounding about 3.5 times. The
central ones divide by 4^j - 1 (3, 15, 63, ...), about 1.4 times. With the level's own floor, forward f' was covered
in 97.5% of the prototype's cases (worst 8.9 times the estimate); with the propagated floor, in 99.4% (worst 2.6).
Central was tuned with the level floor in Phase 3 and keeps it.
</details>

<details><summary>6. Where is fn(x) evaluated for <code>Ridders</code> f'? For <code>Ridders</code> f''?</summary>

For f', nowhere: `Central<2>`'s f' table has no node 0, so `needs_f0` is false. For f'', once, at the top of
`ridders`, and it is reused by every level and every restart.
</details>

<details><summary>7. Why is <code>BackwardRidders</code> not simply <code>ridders&lt;Backward&lt;1&gt;, N&gt;</code>?</summary>

The engine's one-sided logic tests that the first node is 0 (for `one_sided`, the ratio 2, the propagated floor),
and estimates the slope forward from `f(x)` and `f(x + h)`. A backward table has neither. The mirror reuses the
forward engine unchanged and still passes `fn` exactly the points `x - k h`. A `Backward<1>` table would also sum
its terms in a different order and round differently.
</details>

<details><summary>8. A half-precision type has 11 significant bits. What stops <code>derivative&lt;Central&lt;10&gt;, 4&gt;</code> from compiling for it, and why isn't it a hard error?</summary>

`exact_weights<Central<10>, 4, half>`: `largest_weight_integer` is 6,222,216, and 11 bits hold integers only up to

2048. The constraint fails, so no overload is viable: "no matching function". It isn't a hard error because the
      constraint is checked after `method` and `valid_order`, so the table is only built for a method and order that
      exist.

</details>

<details><summary>9. In worked example 1, why is the final estimate 1.6e-13 when the best entry's own estimate was 1.5e-13?</summary>

The diagonal check: `|A[3][3] - A[5][5]|` = |-2.7e-15 - (-1.7e-13)|, about 1.6e-13, is larger, and `abs_error` is the
larger of the two.
</details>

<details><summary>10. In worked example 3, why did the second tableau start at 1/512 and not at 1/64?</summary>

The first tableau ended on a stall at h_i = 1/256: the `break` skips the halving in the loop header, so `h_i` is the
smallest step tried. The restart is at `h_i / 2 = 1/512`. Restarting at `h_start / 2 = 1/64` would re-evaluate steps
that are known to straddle the pole and spend the budget before reaching the pole's scale.
</details>

<details><summary>11. Why is a derivative that is exactly 0 not flagged as untrusted, although its error estimate is far more than a quarter of |0|?</summary>

The untrusted test also needs `abs_error / 64 > best_floor`. For a zero derivative the estimate is about the
rounding floor, so the second condition fails and the result is accepted. Usually the tableau has already ended on
exact agreement (`change == 0`).
</details>

<details><summary>12. Why does <code>sum_fits</code> divide both operands by the radix?</summary>

Dividing by the radix only lowers the exponent, so it is exact, and two halves can't overflow when added. So
`a/2 + b/2 <= max/2` decides `a + b <= max` without ever forming `a + b`, which in constant evaluation would be a
compile error if it overflowed.
</details>

---

## Glossary

| Term                     | Meaning here                                                                                                |
|--------------------------|-------------------------------------------------------------------------------------------------------------|
| **node**                 | a point of a rule in units of the step: node `k` is `x + k h`                                               |
| **stencil**              | a rule with fixed nodes and weights, applied at one step                                                    |
| **rounding floor**       | the library's estimate of the rounding error of a value or an entry                                         |
| **level**                | one evaluation of the base rule at one step in Ridders' method: column 0 of a tableau row                   |
| **tableau, entry**       | Richardson's triangle of extrapolated values; an entry is one value in it                                   |
| **diagonal entry**       | the last entry of a row, the most extrapolated one                                                          |
| **stall**                | a level that doesn't improve the best error estimate; two in a row end a tableau                            |
| **resolution of x**      | the smallest step for which the rule's points are still distinct numbers                                    |
| **untrusted**            | a Ridders result with at most two correct bits by its own estimate, and an error far above rounding         |
| **restart**              | a new tableau with a smaller initial step, after a non-finite level or an untrusted result                  |
| **saturate**             | return `max` instead of overflowing; used for floors only                                                   |
| **constant evaluation**  | running code at compile time, in a `constexpr` or `static_assert` context                                   |
| **expression template**  | an object returned by `a + b` in Boost.Multiprecision that computes the sum only when converted to a number |

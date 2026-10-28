# A walkthrough of `include/numerixx/roots.hpp`

## How to read this document

This is a walkthrough of the implementation: what each part of the code does, how the parts fit together, and why
they are written the way they are. It is the companion of `DERIV-WALKTHROUGH.md` and follows its plan. There is no
book on `nxx::roots`, so this document also says, briefly, what each method computes. The decisions are recorded in
`docs/roots/design.md`, and the measurements behind them in `docs/roots/PHASE1_API.md` to `PHASE6_API.md`.

The header is 2340 lines (as of d89b50a). More than 40% of it is doc comments, and many of them record *why* something
is the way it is: a bug a review found, a count from a prototype, a compiler that rejected the obvious form. This
document explains the code in the order you need to understand it rather than the order of the file.

Three short sections come first, A to C. A and B are the two constraints that shape almost every line: the code must
work in constant expressions, and with user-defined real types. `DERIV-WALKTHROUGH.md` explains both at length; here
they are summarised, with what is particular to root finding. C says what a result of this header means, which every
algorithm below is built to guarantee. Then sections 0 to 17 go through every declaration.

Each section shows the exact code of every function it explains, copied from the header with the doc comments
removed (this text explains them; the `//` comments in the code are kept), so you don't need the header open. The line
numbers still name where it is. Inside the explanations, some fragments are repeated in a shortened form; they are
marked as such. Short notes marked **C++ note** explain a language feature where the code uses one that isn't
everyday C++, and that `DERIV-WALKTHROUGH.md` doesn't already explain.

The numbers in the examples come from running the library (GCC 16, double), with `fn` wrapped so that it records
every point it is called at. I checked every example result against the library.

---

# Before the code

## A. Constant evaluation: check, then compute

Every function in the header is `constexpr`, so a root can be found while the program is compiled. The tests do this
for many cases, for example (`tests/test_roots_bracketing.cpp`):

```cpp
static_assert(rt::solve<rt::Brent>(square_minus_4, rt::Bracket{ 0.0, 4.0 })->value == 2.0);
```

The rule this imposes is the one from `DERIV-WALKTHROUGH.md` A: in a constant expression an operation that overflows
or produces a NaN is not a value you can test afterwards but a **compile error** (GCC rejects both, Clang the NaN).
The library promises an error value instead, at compile time as at run time. So every operation that *could*
overflow is checked before it is done.

A root finder has more such operations than a derivative, because its points are anywhere in the range of `T`:

- the distance between two bracket ends of opposite signs, `b - a`, can exceed `max` (a bracket `[-max, max]` is
  allowed);
- an interpolation multiplies and divides values of `fn`, which can be as large as `max`;
- an open method adds a step to `x`, and the step can be huge where `f'` is tiny;
- `search` grows its step geometrically.

Section 5 holds the helpers that answer "can this be formed?" without forming it. Where the answer is no, each method
has a fallback, usually a bisection step, so that a check never ends a search early.

## B. User-defined real types

The algorithms accept any `T` that satisfies `nxx::real` (`math.hpp`, lines 29-52), such as Boost.Multiprecision's
`cpp_bin_float_50` and `cpp_dec_float_50`. `DERIV-WALKTHROUGH.md` B shows the concept and the three habits it
imposes: constants are written `T(1)`, never `1.0`; a number is never declared `auto` (Boost's expression templates);
calls into the library are qualified (`detail::...`, `math::...`) against argument-dependent lookup.

Four properties of real types matter more here than in `nxx::deriv`.

1. **There is no square root.** `nxx::real` asks for `+ - * /` and comparisons, and nothing more. Ridders' method
   needs a square root, so the header has its own, `detail::sqrt` (section 12).
2. **`denorm_min()` need not be a subnormal.** Boost's binary types have no subnormals, and their `denorm_min()` is
   `min()`. A test type in `tests/test_roots_search.cpp` even has `denorm_min() == 0`. The smallest tolerance is
   `denorm_min`, so the code guards against it being too small to matter (sections 9 and 14).
3. **The radix may be 10.** In `cpp_dec_float_50` halving is not exact (1/2 is exact, but `x / 2` can round), and
   `b - a` for nearby ends is exact by Sterbenz's lemma in every radix, but its boundary case only in radix 2. The
   comments say where an argument holds only in radix 2 (section 5).
4. **The exponent range may be enormous.** `cpp_bin_float_50` reaches about 10^ (+-6.5e8). A bracket that must be
   halved down to the smallest positive number then takes about 2e9 steps instead of about 1100 in double (the
   comment on `Criteria`, section 1).

## C. What a result means

Every method returns `std::expected<Root<T>, error>`. A `Root` is a point `value` and a bound `abs_error`, and the
guarantee behind it is the same for every method:

> `fn`, as evaluated, changes sign within `abs_error` of `value`, and `abs_error` is at most `2 tol(value)`.

**A bracket** is an interval `[a, b]` with `fn(a)` and `fn(b)` of opposite signs. If `fn` is continuous it has a root
inside (the intermediate value theorem); if it isn't, there is still a sign change inside, at a jump or a pole. Every
bracketing method keeps such an interval and narrows it. Each new point `c` is evaluated, and the half with the sign
change is kept: `[a, c]` if `fn(c)` has `fn(b)`'s sign, `[c, b]` otherwise. The methods differ only in how they
choose `c`.

**The stop test.** `tol(x) = abs_tol + rel_tol |x|`, with the caller's tolerances (section 4). A method stops when
the bracket is at most `2 tol(u)` wide, where `u` is the end with the smaller `|fn|`, and returns `u` with the bracket's
width as `abs_error`. It also stops when `fn` is exactly 0 at a point, and then returns that point with `tol(point)`
as `abs_error` (section 1 says why not 0).

**Open methods** start from a guess, not a bracket. They iterate until two consecutive iterates have values of
opposite signs, and from then on they *are* a bracketing method on those two points (the **hand-over**, section 13).
So their results are verified sign changes too, with the same guarantee. The textbook stop, "the step was smaller than
the tolerance", isn't used: in the Phase 2 study it accepted a point where `fn` was -1.

**What no method can see:** `fn`'s own rounding. The root of the exact function may lie `|rounding of fn| / |f'|`
away from the sign change of `fn` as evaluated. The bound is about the computed `fn`, which is all a method has.

---

# The code

## 0. The map

| Lines     | What                                                                                      | Section |
|-----------|-------------------------------------------------------------------------------------------|---------|
| 1-49      | file comment, includes                                                                    |         |
| 53-160    | public types: `error`, `Bracket<T>`, `Root<T>`, `Criteria<T>`                             | 1       |
| 162-293   | public method tags                                                                        | 1       |
| 298-441   | which tags are methods; what `fn` may return                                              | 2, 3    |
| 442-543   | `largest`, `iteration_cap`, `Tolerance`, `twice`, `checked_tolerance`                     | 4       |
| 545-705   | overflow-safe arithmetic: `distance_bound`, `midpoint`, `quotient`, `product`, `sum_fits` | 5       |
| 707-809   | `evaluate`, `derivative_value`, `Joined`                                                  | 3       |
| 811-885   | `Start`, `prepare`: validating a bracket, zero first                                      | 6       |
| 887-918   | `bisection`                                                                               | 7       |
| 920-1098  | `newton_step`, `halley_step`, `derivative_step`, `brent`                                  | 8       |
| 1100-1124 | `difference`, `sum`                                                                       | 5       |
| 1126-1290 | `stop_at`, `Enclosure`, `Enclosed`, `wide`, `place`, `enclose`                            | 9       |
| 1292-1501 | Toms748: `secant_point`, `newton_quadratic`, `inverse_cubic`, ..., `toms748`              | 10      |
| 1503-1575 | `false_position<M>`: AndersonBjorck, Illinois, RegulaFalsi                                | 11      |
| 1577-1699 | `sqrt`, `ridders_point`, `ridders`                                                        | 12      |
| 1701-1895 | the open methods: `hand_over`, `open_derivative`, `open_secant`, `open_steffensen`        | 13      |
| 1900-2017 | the two main `solve` overloads: a bracket, a guess                                        | 2       |
| 2019-2146 | the four `solve` overloads with separate derivatives                                      | 3       |
| 2148-2294 | `SearchSide`, `outward`, `search`                                                         | 14      |
| 2298-2340 | `std::formatter<error>`                                                                   | 15      |

From `math.hpp` the header uses:

- the concept `nxx::real` and `nxx::detail::callable`;
- `math::abs` and `math::isfinite`, written with comparisons only (shown in `DERIV-WALKTHROUGH.md` section 0);
- `math::root_epsilon<4, T>()`, the fourth root of epsilon, for Secant's second point.

**The includes** (roots.hpp, lines 31-48):

```cpp
#pragma once

#include <numerixx/math.hpp>

#include <array>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <expected>
#include <format>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
```

`<cassert>` is new compared with `deriv.hpp`: the header asserts three preconditions (sections 4, 9 and 12). `<tuple>`
is for the derivatives `fn` may return with its value (section 3).

**The call path.** Every call goes through the same three layers as in `nxx::deriv`: a public `solve` checks the call
at compile time; `detail` functions validate the run-time inputs; one engine per method does the work.

```
solve<M>(fn, bracket [, criteria])                    public; M bracketing
 |-- checked_tolerance(criteria)                      validate and clamp the tolerances
 |-- prepare<K>(fn, bracket, tolerance)               validate the ends, evaluate fn there, zero first
 `-- by M:
      |-- Bisection                       -> bisection
      |-- Brent, Newton, Halley           -> brent<M>           (Brent's procedure; M chooses the step)
      |-- Toms748                         -> toms748            (these three take every point
      |-- AndersonBjorck, Illinois,       -> false_position<M>   through enclose: section 9)
      |   RegulaFalsi
      `-- Ridders                         -> ridders

solve<M>(fn, guess [, criteria])                      public; M open
 |-- checked_tolerance(criteria)
 `-- by M:
      |-- Newton, Halley                  -> open_derivative<M> --\
      |-- Secant                          -> open_secant        ---> at the first sign change:
      `-- Steffensen                      -> open_steffensen    --/   hand_over -> brent<M or Brent>

solve<M>(fn, dfn [, d2fn], bracket or guess [, criteria])   Newton, Halley with separate derivatives
 `-- joins fn and its derivatives into one callable (Joined), then calls one of the two above

search(fn, x0)                                        public
 `-- outward: a bracket, for solve
```

Everything below `solve` lives in `namespace detail`, hidden from the user documentation by `@cond INTERNAL`. The
header opens and closes `namespace detail` twice; those lines, and `namespace nxx::roots` around everything, are the
only code lines this document leaves out.

---

## 1. The public vocabulary (lines 53-293)

### `error` (53-81)

```cpp
enum class error {
    non_finite,

    no_sign_change,

    not_converged
};
```

Three causes, each with its own response from the caller:

- **`non_finite`:** an input or a value of `fn` is inf or NaN, or a tolerance is negative. Fix the input, or the
  function's domain.
- **`no_sign_change`:** `fn` has the same sign at both ends of the bracket (and is 0 at neither), or `search` found no
  sign change. Find another bracket.
- **`not_converged`:** an open method found no root from its guess. Use a bracket.

As in `nxx::deriv`, the doc comment says that values will be added and a cause may move to a new value, so a caller
should test for success, and a `switch` on it should have a `default`. Note what the first enumerator's comment says
about derivatives: they are *not* checked, because a derivative that is 0 or not finite only makes a derivative step
unusable, and the method then takes another step (section 8).

### `Bracket<T>` (83-103)

```cpp
template<real T>
struct Bracket
{
    T lower;

    T upper;

    friend constexpr bool operator==(Bracket const &, Bracket const &) = default;
};
```

A plain aggregate, so `Bracket{1.0, 2.0}` deduces `Bracket<double>` by class template argument deduction. The ends may
come in either order: an interval is unambiguous, so `prepare` swaps them rather than reporting an error (section 6).

### `Root<T>` (105-129)

```cpp
template<real T>
struct Root
{
    T value;

    T abs_error;

    friend constexpr bool operator==(Root const &, Root const &) = default;
};
```

The result of every method. The doc comment carries the guarantee from section C, and two details.

- **At an exact zero, `abs_error` is `tol(value)`, not 0.** A point where `fn` evaluates to exactly 0 is only as good
  as `fn`'s evaluation. The comment's example, `(x + 1e10) - (1e10 + 1.0 / 3)`, is 0 on an interval about 2e-6 wide,
  because adding `1e10` loses the low bits of `x`. Brent on `[0, 1]` hits such a point after 3 evaluations and
  returns 0.33333396911621094, 6e-7 from 1/3, with `abs_error` 1.5e-16. Neither 0 nor 1.5e-16 is the true error; the
  bound is about the sign change of `fn` as evaluated, and the tolerance is the honest statement of what was asked.
- **No conversion to `T`, unlike `nxx::deriv::Estimate`.** The members are read by name, `r->value`, so members can be
  added later (a count of evaluations is one deferred candidate) without breaking callers. A structured binding
  `auto [v, e] = *r;` would break, which is why the comment advises against it.

### `Criteria<T>` (131-160)

```cpp
template<real T>
struct Criteria
{
    T abs_tol = T(std::numeric_limits<T>::denorm_min());

    T rel_tol = T(T(2) * T(std::numeric_limits<T>::epsilon()));

    friend constexpr bool operator==(Criteria const &, Criteria const &) = default;
};
```

**C++ note: default member initialisers in an aggregate.** `Criteria` has no constructor, yet its members have
defaults. Since C++14 an aggregate may have default member initialisers, and since C++20 it can be initialised with
designated initialisers: `{.abs_tol = 1e-12}` sets `abs_tol` and leaves `rel_tol` at its default. The default
argument `criteria = {}` of `solve` uses both defaults.

**The defaults are also the minimum.** `abs_tol = denorm_min` and `rel_tol = 2 eps` ask for the full precision of `T`
at every magnitude, and assume no scale. MATLAB's `fzero` and SciPy assume a scale of 1, so a root at 1e-12 loses its
digits there; here it comes out to about 2 eps relative. Smaller tolerances are raised to these, larger ones capped at
`max / 4` and 1 (section 4).

**Their cost** is documented in the comment, and is worth seeing once. A sign change near 0 that `fn` doesn't hit
exactly is approached down to `T`'s smallest numbers:

```
solve<Bisection>(x - 1e-200, [0, 1])         717 evaluations
solve<Brent>(sign(x) jump at 0, [-1, 2])    1076 evaluations, value -9.9e-324
```

Both are correct, and both are halvings of the bracket all the way down from 1 to the root's magnitude. The comment
gives the remedy: set `abs_tol` to the distance from 0 that is close enough.

### The method tags (162-293)

```cpp
struct Bisection
{
};

struct Brent
{
};

struct Toms748
{
};

struct AndersonBjorck
{
};

struct Ridders
{
};

struct RegulaFalsi
{
};

struct Illinois
{
};

struct Newton
{
};

struct Halley
{
};

struct Secant
{
};

struct Steffensen
{
};
```

Eleven empty tags, as in `nxx::deriv` (`DERIV-WALKTHROUGH.md` section 1 explains tag types). What each method does,
in one line each; the sections below explain the code.

| Tag              | Start            | Each new point                                                                      | Section |
|------------------|------------------|-------------------------------------------------------------------------------------|---------|
| `Bisection`      | bracket          | the midpoint                                                                        | 7       |
| `Brent`          | bracket          | inverse quadratic or linear interpolation, or the midpoint (Brent 1973)             | 8       |
| `Toms748`        | bracket          | inverse cubic, Newton on a quadratic, a double secant step, the midpoint (APS 1995) | 10      |
| `AndersonBjorck` | bracket          | the secant through the latest point and the retained end, its value rescaled        | 11      |
| `Ridders`        | bracket          | the midpoint, then the zero of an exponential fit (Ridders 1979)                    | 12      |
| `RegulaFalsi`    | bracket          | the secant through the bracket's ends                                               | 11      |
| `Illinois`       | bracket          | as `RegulaFalsi`, the retained end's value halved                                   | 11      |
| `Newton`         | bracket or guess | `x - f / f'`                                                                        | 8, 13   |
| `Halley`         | bracket or guess | `x - t / (1 - t f'' / (2 f'))`, `t = f / f'`                                        | 8, 13   |
| `Secant`         | guess            | the secant through the latest two iterates                                          | 13      |
| `Steffensen`     | guess            | Newton with the slope `(f(x + f(x)) - f(x)) / f(x)`                                 | 13      |

**The doc comments are the user-facing specification.** Each states the cost on test cases, the worst case and when
to use the method. Two names need care. `Ridders` is not `nxx::deriv::Ridders`, the derivative method of the same
author; inside the header the latter is always written out in full. And `Toms748` is "Algorithm 748" of ACM TOMS, by
Alefeld, Potra and Shi; the comment notes that it returns a point and a bound, not the final bracket as Boost.Math's
`toms748_solve` does.

**Teaching examples.** `Ridders`, `RegulaFalsi`, `Illinois` and `Steffensen` are in the library because they are in
every textbook, and the comments say plainly that another method is better: Brent on a bracket, Secant from a guess.

---

## 2. Following a call from the top (lines 298-323, 1900-2017)

Take `solve<Brent>(fn, Bracket{1.0, 2.0})` and follow it down.

### The two main overloads

**The code** (roots.hpp, lines 1947-1971 and 2004-2017; the two `//` lines above the first belong to it):

```cpp
// criteria is a non-deduced context (type_identity_t), so T comes from the bracket alone and a braced list of
// doubles converts to a multiprecision T.
template<typename M, real T, typename F>
requires detail::bracketing<M> && detail::gives<F, T, detail::values_needed<M>>
[[nodiscard]] constexpr std::expected<Root<T>, error>
    solve(F &&fn, Bracket<T> const &bracket, Criteria<std::type_identity_t<T>> const &criteria = {})
{
    constexpr std::size_t K = detail::values_needed<M>;
    // Qualified, so that argument-dependent lookup cannot pick a function from the namespace of T or F.
    std::expected<detail::Tolerance<T>, error> const tolerance = detail::checked_tolerance(criteria);
    if (!tolerance) return std::unexpected(tolerance.error());
    std::expected<detail::Start<T, K>, error> const start = detail::prepare<K>(fn, bracket, *tolerance);
    if (!start) return std::unexpected(start.error());
    if (start->root) return *start->root;
    if constexpr (std::same_as<M, Bisection>)
        return detail::bisection(fn, *start, *tolerance);
    else if constexpr (std::same_as<M, Toms748>)
        return detail::toms748(fn, *start, *tolerance);
    else if constexpr (std::same_as<M, AndersonBjorck> || std::same_as<M, Illinois> || std::same_as<M, RegulaFalsi>)
        return detail::false_position<M>(fn, *start, *tolerance);
    else if constexpr (std::same_as<M, Ridders>)
        return detail::ridders(fn, *start, *tolerance);
    else
        return detail::brent<M>(fn, *start, *tolerance);
}

template<typename M, real T, typename F>
requires detail::open<M> && detail::gives<F, T, detail::values_needed<M>>
[[nodiscard]] constexpr std::expected<Root<T>, error>
    solve(F &&fn, T const &guess, Criteria<std::type_identity_t<T>> const &criteria = {})
{
    std::expected<detail::Tolerance<T>, error> const tolerance = detail::checked_tolerance(criteria);
    if (!tolerance) return std::unexpected(tolerance.error());
    if constexpr (std::same_as<M, Secant>)
        return detail::open_secant(fn, guess, *tolerance);
    else if constexpr (std::same_as<M, Steffensen>)
        return detail::open_steffensen(fn, guess, *tolerance);
    else
        return detail::open_derivative<M>(fn, guess, *tolerance);
}
```

**The constraints.** Read the first as: M is a bracketing method, and `fn` gives at least as many values at a point as
M needs. The second: M is an open method, with the same condition on `fn`. The concepts are just lists (lines
298-323):

```cpp
template<typename M>
concept bracketing = std::same_as<M, Bisection> || std::same_as<M, Brent> || std::same_as<M, Toms748>
                     || std::same_as<M, AndersonBjorck> || std::same_as<M, Ridders>
                     || std::same_as<M, RegulaFalsi> || std::same_as<M, Illinois> || std::same_as<M, Newton>
                     || std::same_as<M, Halley>;

template<typename M>
concept open = std::same_as<M, Newton> || std::same_as<M, Halley> || std::same_as<M, Secant>
               || std::same_as<M, Steffensen>;

template<typename M>
inline constexpr std::size_t values_needed = std::same_as<M, Newton>   ? 2
                                             : std::same_as<M, Halley> ? 3
                                                                       : 1;
```

**The two overloads can't both match.** The bracket overload takes a `Bracket<T> const &`, and `T` is deduced from it;
the guess overload deduces `T` from a `T const &`. A `Bracket<double>` isn't a `real`, so `solve<Newton>(fdf,
Bracket{1.0, 2.0})` deduces `T = double` in the first overload and fails `real<Bracket<double>>` in the second. A
guess `1.0` can't deduce a `Bracket<T>`. Newton and Halley are in both lists, and the start chooses.

**`values_needed<M>`** is how many values `fn` must give at a point: 1 (the value) for most methods, 2 for Newton (with
`f'`), 3 for Halley (with `f''`). `detail::gives<F, T, K>` (section 3) checks that `fn` gives at least K.

**`Criteria<std::type_identity_t<T>>`.** `T` is a non-deduced context there (`DERIV-WALKTHROUGH.md` section 2 has
the C++ note), so `T` comes from the bracket or the guess alone, and `{.rel_tol = 1e-10}` converts to a
`Criteria<cpp_bin_float_50>`. Without it, the braced list couldn't be deduced at all.

**`F &&fn`, passed on as the lvalue `fn`.** As in `nxx::deriv`: `fn` is called many times, so it is never forwarded
(moved from), and the engines take it as `F &`.

**The body of the bracket overload** does three things in order:

1. `checked_tolerance` validates and clamps the tolerances (section 4). A bad tolerance is `non_finite` before `fn` is
   ever called.
2. `prepare<K>` validates the bracket and evaluates `fn` at its ends, and at 0 if the bracket contains it (section
   6). It may already find the root, at an end or at 0; then `start->root` holds it.
3. `if constexpr` chooses the engine. Every engine receives the validated `Start` and the `Tolerance`.

**The body of the guess overload** has no `prepare`: each open method validates its guess itself, because the first
evaluation is part of its iteration.

**C++ note: `std::unexpected(x.error())`.** Two lines like

```cpp
if (!tolerance) return std::unexpected(tolerance.error());
```

pass a failure on. The types differ (`expected<Tolerance<T>, error>` and `expected<Root<T>, error>`), so the error is
taken out and wrapped again. `.and_then` would do the same, but a chain of lambdas reads worse here, where each step
needs several earlier results.

**The last `else` is `brent<M>` for three methods.** Brent, Newton and Halley on a bracket are one loop with three
ways of choosing a step (section 8). `if constexpr` discards the branches not taken, so `solve<Bisection>` contains no
Brent code at all.

---

## 3. What `fn` may return (lines 325-441, 707-809, 2019-2146)

Newton needs `f` and `f'` at a point, Halley also `f''`. They can come from `fn` itself, returned together, or from
separate callables. Both forms reach the same engine: the separate form is turned into the first.

### The tuple form: `tuple_valued` (325-383)

```cpp
template<typename R>
inline constexpr bool is_standard_tuple = false;

template<typename A, typename B>
inline constexpr bool is_standard_tuple<std::pair<A, B>> = true;

template<typename... Ts>
inline constexpr bool is_standard_tuple<std::tuple<Ts...>> = true;

template<typename U, std::size_t N>
inline constexpr bool is_standard_tuple<std::array<U, N>> = true;

template<typename R, typename T, std::size_t... I>
consteval bool elements_are(std::index_sequence<I...>)
{ return (std::same_as<std::tuple_element_t<I, R>, T> && ...); }

template<typename F, typename T>
concept tuple_valued =
    std::invocable<F &, T const &> && !std::is_reference_v<std::invoke_result_t<F &, T const &>>
    && !std::is_volatile_v<std::invoke_result_t<F &, T const &>>
    && is_standard_tuple<std::remove_cv_t<std::invoke_result_t<F &, T const &>>>
    && detail::elements_are<std::remove_cv_t<std::invoke_result_t<F &, T const &>>, T>(
        std::make_index_sequence<std::tuple_size_v<std::remove_cv_t<std::invoke_result_t<F &, T const &>>>>{});
```

`fn` may return a `std::pair`, a `std::tuple` or a `std::array`, **by value, with every element exactly `T`**:

```cpp
auto const fdf = [](double x) { return std::pair{x * x - 2.0, 2.0 * x}; };   // (f, f')
```

**Why exactly `T`, not convertible to `T`.** The comment explains: an element of reference type (`std::tie(f, df)`
of `fn`'s locals) refers to locals that die when `fn` returns, and an element of an expression type (a Boost
expression template, `std::pair{x * x - 2, 2 * x}` without the `T(...)`) holds references to the same locals. Reading
it after `fn` returned is too late to convert it. Requiring `T` itself rules both out at compile time. For the same
reason the result itself must not be a reference or volatile.

**How the check works.** `is_standard_tuple<R>` is a variable template, false in general and true for the three
standard tuple-likes, by partial specialisation (the C++ note in `DERIV-WALKTHROUGH.md` section 3). `elements_are`
then checks every element type:

**C++ note: `std::index_sequence` and a fold expression.** `std::make_index_sequence<3>{}` is an empty object of type
`std::index_sequence<0, 1, 2>`. Passing it to `elements_are(std::index_sequence<I...>)` deduces the pack `I...` as
`0, 1, 2`, and the body expands `(std::same_as<std::tuple_element_t<I, R>, T> && ...)` into
`same_as<element 0, T> && same_as<element 1, T> && same_as<element 2, T>`. It is the standard way to loop over the
elements of a tuple at compile time. `consteval` (`DERIV-WALKTHROUGH.md` section 4) makes sure it never runs at run
time.

### The value form, and `gives` (385-441)

```cpp
template<typename F, typename T>
concept value_valued = nxx::detail::callable<F, T>
                       && !is_standard_tuple<std::remove_cvref_t<std::invoke_result_t<F &, T const &>>>;

template<typename R>
inline constexpr bool is_expected = false;

template<typename U, typename E>
inline constexpr bool is_expected<std::expected<U, E>> = true;

// is_expected comes first, so that *r is formed only for an expected.
template<typename D, typename T>
concept derivative_valued =
    value_valued<D, T>
    || (std::invocable<D &, T const &> && is_expected<std::remove_cvref_t<std::invoke_result_t<D &, T const &>>>
        && std::convertible_to<decltype(*std::declval<std::invoke_result_t<D &, T const &>>()), T>);

template<typename F, typename T, std::size_t K>
concept gives =
    (K == 1 && value_valued<F, T>)
    || (tuple_valued<F, T> && std::tuple_size_v<std::remove_cv_t<std::invoke_result_t<F &, T const &>>> >= K);
```

- **`value_valued<F, T>`:** `fn` returns something convertible to `T` (the shared `callable`), and it isn't a
  tuple-like. The exclusion matters for a type `T` that can be constructed from a pair: a tuple result is never read
  as a single value.
- **`derivative_valued<D, T>`:** what a *separate* derivative may return: a value, or a `std::expected` whose value
  converts to `T`. That is what `nxx::deriv::derivative_of(fn)` returns (`DERIV-WALKTHROUGH.md` section 11), so a
  numerical derivative can be passed directly. Any `std::expected<U, E>` is accepted, so `nxx::roots` names no type
  of `nxx::deriv`.
- **`gives<F, T, K>`:** for K = 1, a value or a tuple (a bracketing method on a tuple `fn` uses its first element);
  for K = 2 or 3, a tuple of at least K elements.

**Why `*std::declval<...>()`.** The constraint checks the very expression the body converts (`derivative_value`,
below): `*r` on the result as returned, with its reference and constness. The comment records why: a conversion the
constraint admits must be the one the body performs, and the Phase 3 review found a gap between the two.

**The order inside `derivative_valued` matters.** `is_expected<...>` comes before the `convertible_to` that forms
`*r`. A `&&` of constraints stops at the first false one, so `*r` is only formed for an `expected`; for any other type
it could be a hard error.

### Reading the values: `evaluate` (707-740)

```cpp
template<std::size_t K, real T, typename F, std::size_t... I>
constexpr std::array<T, K> evaluate(F &fn, T const &x, std::index_sequence<I...>)
{
    if constexpr (tuple_valued<F, T>) {
        std::invoke_result_t<F &, T const &> const values = std::invoke(fn, x);
        return std::array<T, K>{ std::get<I>(values)... };
    } else {
        // A value: expression templates are converted to T before any arithmetic.
        return std::array<T, K>{ static_cast<T>(std::invoke(fn, x)) };
    }
}

template<std::size_t K, real T, typename F>
constexpr std::array<T, K> evaluate(F &fn, T const &x)
{ return detail::evaluate<K>(fn, x, std::make_index_sequence<K>{}); }
```

Every engine calls `fn` through `detail::evaluate<K>(fn, x)`, which returns a `std::array<T, K>`: index 0 is the
value, 1 is `f'`, 2 is `f''`. A tuple result is stored as returned (exactly `T`, so nothing dangles) and its first K
elements copied out; a value is converted with `static_cast<T>` at once, which also evaluates an expression template.
The engines then work with arrays of the same shape for every form of `fn`: `vb[0]` is always the value.

### Separate derivatives: `derivative_value` and `Joined` (742-809)

```cpp
template<real T, typename D>
constexpr T derivative_value(D &dfn, T const &x)
{
    if constexpr (value_valued<D, T>) {
        return static_cast<T>(std::invoke(dfn, x));
    } else {
        // Called once, and the result kept as returned: a copy would require a copyable expected, and a const
        // copy would refuse a value that converts only as an rvalue (independent review of Phase 3).
        decltype(auto) r = std::invoke(dfn, x);
        if (!r) return T(0);
        return static_cast<T>(*std::forward<decltype(r)>(r));
    }
}
```

A derivative that returns a value is converted; one that returns a `std::expected` is tested, and a failure becomes
**0**. Why 0, and not an error? The comment explains: a failed derivative is treated as one that isn't finite, and 0
has the same effect in both steps. `f' = 0` makes Newton's step unusable, as `f' = inf` would; `f'' = 0` makes
Halley's step exactly Newton's, as `f'' = inf` would. And every `real` has a 0, while a NaN or an infinity need not
exist. The method then takes a bisection step on a bracket, or returns `not_converged` from a guess.

**C++ note: `decltype(auto)` and forwarding a result.** `decltype(auto) r = std::invoke(dfn, x);` declares `r` with
exactly the type the call returns: a value if it returns by value, a reference if it returns a reference. Then
`*std::forward<decltype(r)>(r)` dereferences `r` as an rvalue if it was returned by value, so a value type that
converts to `T` only when moved (an rvalue-only conversion) still works, and as an lvalue if it was a reference. The
comment gives the reasons against a copy, `auto const r = ...`: it would require a copyable `expected`, and a const
copy would refuse a value that converts only as an rvalue (Phase 3 review).

```cpp
template<real T, typename F, typename... Ds>
struct Joined
{
    F *fn;

    std::tuple<Ds *...> derivatives;

    template<std::size_t... I>
    constexpr std::array<T, 1 + sizeof...(Ds)> values(T const &x, std::index_sequence<I...>) const
    {
        return std::array<T, 1 + sizeof...(Ds)>{ static_cast<T>(std::invoke(*fn, x)),
            detail::derivative_value(*std::get<I>(derivatives), x)... };
    }

    constexpr std::array<T, 1 + sizeof...(Ds)> operator()(T const &x) const
    { return values(x, std::index_sequence_for<Ds...>{}); }
};
```

`Joined` turns `fn` plus one or two derivative callables into one callable that returns `std::array<T, 2>` or
`std::array<T, 3>`, the tuple form above. Each callable is called once per point. The pack `Ds...` holds the
derivatives' types, and `values` expands `derivative_value(*std::get<I>(derivatives), x)...` over them with an index
sequence, as `elements_are` does.

**Why pointers, not references or a lambda.** `fn` may be a function, so `F` may be a function type, and a lambda
capturing a function reference fails on MSVC in constant evaluation. The same problem, and the same fix, as
`nxx::deriv`'s `Mirrored` (`DERIV-WALKTHROUGH.md` section 10).

### The four overloads with separate derivatives (2019-2146)

**The code** (roots.hpp, lines 2063-2071, 2085-2093, 2109-2121 and 2137-2146: Newton on a bracket and from a guess,
then Halley on a bracket and from a guess):

```cpp
template<typename M, real T, typename F, typename D>
requires std::same_as<M, Newton> && detail::value_valued<F, T> && detail::derivative_valued<D, T>
[[nodiscard]] constexpr std::expected<Root<T>, error>
    solve(F &&fn, D &&dfn, Bracket<T> const &bracket, Criteria<std::type_identity_t<T>> const &criteria = {})
{
    detail::Joined<T, std::remove_reference_t<F>, std::remove_reference_t<D>> const joined{ std::addressof(fn),
        { std::addressof(dfn) } };
    return nxx::roots::solve<M>(joined, bracket, criteria);
}

template<typename M, real T, typename F, typename D>
requires std::same_as<M, Newton> && detail::value_valued<F, T> && detail::derivative_valued<D, T>
[[nodiscard]] constexpr std::expected<Root<T>, error>
    solve(F &&fn, D &&dfn, T const &guess, Criteria<std::type_identity_t<T>> const &criteria = {})
{
    detail::Joined<T, std::remove_reference_t<F>, std::remove_reference_t<D>> const joined{ std::addressof(fn),
        { std::addressof(dfn) } };
    return nxx::roots::solve<M>(joined, guess, criteria);
}

template<typename M, real T, typename F, typename D, typename D2>
requires std::same_as<M, Halley> && detail::value_valued<F, T> && detail::derivative_valued<D, T>
         && detail::derivative_valued<D2, T>
[[nodiscard]] constexpr std::expected<Root<T>, error> solve(F &&fn,
    D &&dfn,
    D2 &&d2fn,
    Bracket<T> const &bracket,
    Criteria<std::type_identity_t<T>> const &criteria = {})
{
    detail::Joined<T, std::remove_reference_t<F>, std::remove_reference_t<D>, std::remove_reference_t<D2>> const
        joined{ std::addressof(fn), { std::addressof(dfn), std::addressof(d2fn) } };
    return nxx::roots::solve<M>(joined, bracket, criteria);
}

template<typename M, real T, typename F, typename D, typename D2>
requires std::same_as<M, Halley> && detail::value_valued<F, T> && detail::derivative_valued<D, T>
         && detail::derivative_valued<D2, T>
[[nodiscard]] constexpr std::expected<Root<T>, error>
    solve(F &&fn, D &&dfn, D2 &&d2fn, T const &guess, Criteria<std::type_identity_t<T>> const &criteria = {})
{
    detail::Joined<T, std::remove_reference_t<F>, std::remove_reference_t<D>, std::remove_reference_t<D2>> const
        joined{ std::addressof(fn), { std::addressof(dfn), std::addressof(d2fn) } };
    return nxx::roots::solve<M>(joined, guess, criteria);
}
```

Each wraps the callables in a `Joined` and calls the main `solve` again. The constraint
`std::same_as<M, Newton> && value_valued<F, T> && derivative_valued<D, T>` limits the first two to Newton with a plain
`fn` (the last two to Halley, with `D2` checked too); a tuple-returning `fn` uses the main overloads.

So these two calls run the same code:

```cpp
rt::solve<rt::Newton>([](double x) { return std::pair{x * x - 2.0, 2.0 * x}; }, 1.0);
rt::solve<rt::Newton>([](double x) { return x * x - 2.0; }, [](double x) { return 2.0 * x; }, 1.0);
```

With a numerical derivative, `solve<Newton>(f, nxx::deriv::derivative_of(f), 1.0)` takes the same 7 points as the
analytic one on `x^2 - 2`, but each costs 3 calls of `f` (`Central<2>`: `f(x - h)`, `f(x + h)`, and `f(x)` itself), 21
in all. The doc comment of the first overload compares such costs with Secant and Brent, which need no derivative.

**C++ note: why these overloads don't collide with the main ones.** `solve<Newton>(fn, dfn, Bracket{1.0, 2.0})`
has three arguments, like `solve<Newton>(fdf, bracket, criteria)`. Tried against the main overload, its second
argument, the lambda `dfn`, would have to deduce a `Bracket<T>`, which fails. In the other direction, the separate
overload would take the `Bracket` as `D` and then fail to deduce `Bracket<T>` from a `Criteria`. Deduction failure
removes an overload quietly, so exactly one remains.

---

## 4. Tolerances (lines 442-543)

### `largest` and `iteration_cap` (442-472)

```cpp
template<real T>
constexpr T largest()
{
    // Parenthesised, so that a max macro, as defined by <windows.h>, cannot expand.
    return T((std::numeric_limits<T>::max)());
}

template<real T>
consteval int iteration_cap()
{
    // In 64 bits: an admitted radix may exceed 2^30, where an int shift would overflow.
    long long bits_per_digit = 0;
    long long reach = 1;
    while (reach < std::numeric_limits<T>::radix) {
        reach *= 2;
        ++bits_per_digit;
    }
    long long const cap = 4LL * std::numeric_limits<T>::digits * bits_per_digit;
    return cap > (std::numeric_limits<int>::max)() ? (std::numeric_limits<int>::max)() : static_cast<int>(cap);
}
```

`largest<T>()` is `max`, written once with the parentheses that stop a `<windows.h>` `max` macro.

`iteration_cap<T>()` is the open methods' cap, and `search`'s: 4 times the bits of `T`. The number of bits is
`digits` times the bits per digit, the smallest `b` with `2^b >= radix`: 1 for binary, 4 for decimal. That is 212 for
double, 672 for `cpp_bin_float_50`. Why 4 times: Newton doubles the correct bits per step near a simple root, so it
needs far fewer; the cap is reached only by a cycle, a divergence, or a multiple root (each step then gains a fixed
fraction of a bit), and the doc comment of `Newton` lists those. The arithmetic is in `long long` and clamped,
because an admitted radix could be large enough to overflow an `int`.

**There is no cap on the bracketing methods.** Bisection ends at `T`'s resolution, Brent within `(k + 1)^2 - 2`
evaluations, Toms748 and the false-position methods within about 4 per halving. A cap without an error to report it
would return the last point silently (`design.md`, Behaviour).

### `Tolerance` and `twice` (474-512)

```cpp
template<real T>
struct Tolerance
{
    T absolute;

    T relative;

    constexpr T at(T const &x) const
    {
        // scaled is at most |x|, since relative <= 1. Above max / 2 it saturates: max - absolute would round,
        // and could admit a sum that overflows; below it the sum is at most 3/4 max, since absolute <= max / 4.
        T const scaled = T(relative * math::abs(x));
        if (scaled > T(detail::largest<T>() / T(2))) return detail::largest<T>();
        return T(absolute + scaled);
    }
};

template<real T>
constexpr T twice(T const &t)
{
    if (t > T(detail::largest<T>() / T(2))) return detail::largest<T>();
    return T(T(2) * t);
}
```

`Tolerance` is the validated form of `Criteria`, and `at(x)` is `tol(x) = absolute + relative |x|`, **saturated at
`max`**. The comment gives the argument that it can't overflow: `relative <= 1`, so `scaled <= |x| <= max`; if it
exceeds `max / 2` the result is `max`; otherwise `absolute + scaled <= max / 4 + max / 2`. Why not test `absolute +
scaled > max` directly? Because that sum is what can't be formed. And why not `scaled > max - absolute`? Because
`max - absolute` rounds back to `max` for any small `absolute`, and would let an overflowing sum through.

`twice(t)` is `2 t`, saturated the same way. The stop test compares a width with `twice(tol)`.

### `checked_tolerance` (514-543)

```cpp
template<real T>
constexpr std::expected<Tolerance<T>, error> checked_tolerance(Criteria<T> const &criteria)
{
    // An assert would vanish under NDEBUG, and a NaN tolerance would then never stop Brent: validated here.
    if (!math::isfinite(criteria.abs_tol) || !math::isfinite(criteria.rel_tol) || criteria.abs_tol < T(0)
        || criteria.rel_tol < T(0))
        return std::unexpected(error::non_finite);
    T const tiny = T(std::numeric_limits<T>::denorm_min());
    assert(tiny > T(0) && "the real type's denorm_min() must be positive");
    T const absolute_cap = T(detail::largest<T>() / T(4));
    T const relative_floor = T(T(2) * T(std::numeric_limits<T>::epsilon()));
    T absolute = criteria.abs_tol;
    if (absolute < tiny) absolute = tiny;
    if (absolute > absolute_cap) absolute = absolute_cap;
    T relative = criteria.rel_tol;
    if (relative < relative_floor) relative = relative_floor;
    if (relative > T(1)) relative = T(1);
    return Tolerance<T>{ absolute, relative };
}
```

**Validated with a run-time check, not an `assert`.** A NaN tolerance would make every comparison `width <= 2 tol`
false, and Brent would never stop. An `assert` disappears under `NDEBUG`. A tolerance is the caller's input, so it is
validated at the boundary, as in `nxx::deriv`.

**Then clamped.** The lower clamps make every method end: Brent's minimum step is `tol(b)`, and with `rel_tol >= 2
eps` it is at least two units of `T`'s spacing at `b` in radix 2 and 10, so a step of `tol` always moves `b`. Brent
himself requires at least `eps |b|`. The upper clamps keep `tol` and `2 tol` meaningful.

**The `assert`** is a different matter: `denorm_min() > 0` is a property of the type, not of the input, so a type
that violates it is a programming error, not an error value. The lower clamp relies on it: at 0, `tol` is just
`abs_tol`, and the argument that a minimum step moves the point needs it positive. (`search`, which uses no
tolerance, accepts such a type; section 14.)

**C++ note: `assert` in a `constexpr` function.** `assert` expands to a call of a non-`constexpr` function when the
condition is false. In constant evaluation that call makes the expression non-constant, so a failed `assert` is a
compile error there, and a passing one costs nothing. At run time it is the usual `assert`.

---

## 5. Overflow-safe arithmetic (lines 545-705, 1100-1124)

These helpers decide whether a difference, sum, product or quotient can be formed before forming it. They are the
reason section A matters: without constant evaluation, most of them could be replaced by "compute, then test with
`isfinite`".

### `same_sign` (545-553)

```cpp
template<real T>
constexpr bool same_sign(T const &a, T const &b)
{ return (a < T(0)) == (b < T(0)); }
```

For two nonzero values. Written with comparisons rather than `a * b > 0`, which could overflow or underflow to 0.

### `distance_bound(u, v)` (555-617): an upper bound on `|u - v|`

```cpp
template<real T>
constexpr T distance_bound(T const &u, T const &v)
{
    if (u == v) return T(0);
    T const max = detail::largest<T>();
    T const tiny = T(std::numeric_limits<T>::denorm_min());
    T const au = math::abs(u);
    T const av = math::abs(v);
    if (au == T(0)) return av; // a distance from 0 is exact
    if (av == T(0)) return au;
    T distance = T(0);
    T units = T(1);
    if ((u < T(0)) == (v < T(0))) {
        // Of one sign, |u - v| <= max(|u|, |v|): finite.
        distance = math::abs(T(u - v));
        if (distance == T(0)) return tiny;
        T const smaller = au < av ? au : av;
        if (distance < smaller) return distance;
        // The boundary case of Sterbenz's lemma, the larger end exactly twice the smaller, is exact too. It is
        // recognised in radix 2 only, where doubling is exact; in radix 10 it is raised, conservatively.
        if constexpr (std::numeric_limits<T>::radix == 2) {
            T const larger = au < av ? av : au;
            if (distance == smaller && (smaller > T(max / T(2)) || larger <= T(T(2) * smaller)))
                return distance;
        }
    } else {
        // Of opposite signs the distance is |u| + |v|: finite (at most max / 2) when both are at most max / 4;
        // above that, by halves. Saturating there instead left a bracket around 0 under huge tolerances wider
        // than its true width forever (independent review of Phase 2).
        T const quarter = T(max / T(4));
        if (au <= quarter && av <= quarter) {
            distance = T(au + av);
        } else {
            T const half_sum = T(T(au / T(2)) + T(av / T(2)));
            if (half_sum > T(max / T(2))) return max;
            distance = T(T(2) * half_sum);
            units = T(2);
        }
    }
    T raise = T(T(distance * T(std::numeric_limits<T>::epsilon())) * units);
    if (raise < tiny) raise = tiny;
    if (raise >= T(max - distance)) return max; // max - distance is exact above max / 2 (Sterbenz)
    return T(distance + raise);
}
```

This is the most important helper: every `abs_error`, every stop test and every guard uses it. It returns an upper
bound on `|u - v|` that is **never below** the true distance, formed without overflow. Its cases:

1. **Equal points, or one point 0:** the distance is exact (`0`, or the other's magnitude).
2. **One sign** (0 counts as positive): `|u - v| <= max(|u|, |v|)`, so the difference can be formed. It is exact
   when the computed difference is below the smaller magnitude. That is Sterbenz's lemma: if `v/2 <= u <= 2v`, then
   `u - v` is exact in floating point. This is the case at every stop with the default tolerances, where the ends are
   a few units apart, so `abs_error` there is exactly the bracket's width. In radix 2 the boundary case, the larger
   end exactly twice the smaller, is exact too; in radix 10 it isn't recognised, which only costs a slightly larger
   bound.
3. **A difference of 0 between ends that differ:** a type without subnormals (Boost's `cpp_bin_float`) flushes a
   difference below its smallest number to 0. The distance is then at most `denorm_min`, which for such a type is
   `min`.
4. **Opposite signs:** the distance is `|u| + |v|`. It is formed directly when both are at most `max / 4`, and by
   halves otherwise. If the halves sum to more than `max / 2`, the distance exceeds `max`, and `max` is returned.
5. **Rounding:** where the difference may have rounded (cases 2 and 4 when not exact), the result is raised by
   `max(d eps, denorm_min)` (twice that when it was computed by halves), at least one unit of `T`'s spacing at `d`.
   The raised value is formed only if it can't pass `max`: `max - distance` can always be formed, and it is exact
   above `max / 2` (Sterbenz), the only place where the test can fail.

The comment on case 4 records a review finding. An earlier form returned `max` for any opposite-sign distance above
`max / 4`. Under huge tolerances a bracket around 0 then stayed "wider than it is" forever, and Brent looped. Measuring
by halves fixed it.

### `bounded_difference(u, v, limit)` (619-635)

```cpp
template<real T>
constexpr T bounded_difference(T const &u, T const &v, T const &limit)
{
    // A bound at most limit (at most max / 16) means the true distance is too: u - v is finite.
    if (detail::distance_bound(u, v) > limit) return detail::largest<T>();
    return T(u - v);
}
```

`u - v`, or `max` where `|u - v|` exceeds `limit`. Brent uses it for its initial step `b - a`, with `limit = max / 16`.
Returning `max` there makes Brent take bisection steps (its guard tests `|e| > max / 16`) until the bracket is
narrow enough for the difference to be formed (section 8).

### `midpoint(a, b)` (637-648)

```cpp
template<real T>
constexpr T midpoint(T const &a, T const &b)
{
    if ((a < T(0)) != (b < T(0))) return T(T(a + b) / T(2));
    return T(a + T(T(b - a) / T(2)));
}
```

The textbook `(a + b) / 2` overflows for ends of one sign near `max`; `a + (b - a) / 2` overflows for ends of
opposite signs far apart. Each form is safe where the other isn't: ends of opposite signs have `|a + b| <=
max(|a|, |b|)`, ends of one sign have `|b - a| <= max(|a|, |b|)`. So the function picks by the signs.

### `quotient(a, b)` and `product(a, b)` (650-683)

```cpp
template<real T>
constexpr std::optional<T> quotient(T const &a, T const &b)
{
    if (b == T(0)) return std::nullopt;
    T const ab = math::abs(b);
    if (ab >= T(1)) return T(a / b);
    if (!(math::abs(a) < T(T(detail::largest<T>() / T(2)) * ab))) return std::nullopt;
    return T(a / b);
}

template<real T>
constexpr std::optional<T> product(T const &a, T const &b)
{
    T const aa = math::abs(a);
    if (aa > T(1) && math::abs(b) > T(T(detail::largest<T>() / T(2)) / aa)) return std::nullopt;
    return T(a * b);
}
```

Both return `std::optional<T>`: the result where it can be formed, nothing otherwise.

- **`quotient`:** where `|b| >= 1`, `|a / b| <= |a|`, so the division is safe. Otherwise it requires `|a| < (max /
    2) |b|`, which guarantees `|a / b| < max / 2`. `(max / 2) |b|` itself can't overflow for `|b| < 1`, and if it
  underflows to 0 the test fails, conservatively. The comment records a review finding: an earlier one-line test,
  `|a| / max < |b| / 2`, divided `|a|` by `max` and underflowed to 0 for tiny `a`, rejecting quotients like
  `denorm_min / denorm_min` that are perfectly representable.
- **`product`:** where `|a| <= 1`, `|a b| <= |b|`. Otherwise it requires `|b| <= (max / 2) / |a|`.

Both keep their results below about `max / 2`, so one more addition of a like quantity can't overflow either.

**C++ note: `std::optional` as "maybe a number".** The interpolation formulas of sections 10 to 12 chain many such
operations. Each step returns `std::optional<T>`, and the next step takes it and returns nothing if its input was
nothing. A failed check anywhere makes the whole formula return nothing, and the method then uses its fallback point.
It is a hand-written form of the "maybe" monad; `std::optional::and_then` would do the same with lambdas.

### `sum_fits(x, s)` (685-705): can `x + s` be formed?

```cpp
template<real T>
constexpr bool sum_fits(T const &x, T const &s)
{
    if ((x < T(0)) != (s < T(0))) return true;
    T const half = T(detail::largest<T>() / T(2));
    T const ax = math::abs(x);
    T const as = math::abs(s);
    T const larger = ax < as ? as : ax;
    T const smaller = ax < as ? ax : as;
    if (larger <= half) return true;
    return smaller <= T(detail::largest<T>() - larger);
}
```

The cases:

- **Opposite signs:** `|x + s| <= max(|x|, |s|)`, always safe.
- **One sign, both at most `max / 2`:** the sum is at most `max`.
- **One sign, the larger above `max / 2`:** then `max - larger` is exact (Sterbenz again: `larger` and `max` are
  within a factor 2), and the sum fits exactly when `smaller <= max - larger`.

The comment records the latest review finding (Phase 6). The first form tested `x`'s magnitude rather than the larger
one: it required `|s| <= max / 2` wherever `|x| <= max / 2`, so it rejected `0 + max`, which fits. Steffensen on
`f(x) = max - x` from 0 then gave up, although its very first point, `z = 0 + max`, is the root. The test is now
symmetric in its arguments, and exact.

### `difference` and `sum` (1100-1124)

```cpp
template<real T>
constexpr std::optional<T> difference(T const &x, T const &y)
{
    if (!detail::sum_fits(x, T(-y))) return std::nullopt;
    return T(x - y);
}

template<real T>
constexpr std::optional<T> sum(std::optional<T> const &x, std::optional<T> const &y)
{
    if (!x || !y || !detail::sum_fits(*x, *y)) return std::nullopt;
    return T(*x + *y);
}
```

Two small wrappers on `sum_fits`, in the `std::optional` style, for the interpolations: `x - y` checked, and the sum of
two maybe-values. Unlike Brent's, the interpolations of Toms748 and Ridders form differences of `fn`'s values of
*opposite* signs, which can exceed `max`; these check them.

---

## 6. Starting from a bracket: `prepare` (lines 811-885)

### `Start<T, K>` (811-833)

```cpp
template<real T, std::size_t K>
struct Start
{
    std::optional<Root<T>> root;

    T lower;

    std::array<T, K> at_lower;

    T upper;

    std::array<T, K> at_upper;
};
```

What every bracketing engine starts from: the bracket in order, `fn`'s K values at each end, and possibly a root
already found. The values come with the ends so that no engine evaluates an end again.

### `prepare<K>(fn, bracket, tolerance)` (835-885)

```cpp
template<std::size_t K, real T, typename F>
constexpr std::expected<Start<T, K>, error>
    prepare(F &fn, Bracket<T> const &bracket, Tolerance<T> const &tolerance)
{
    if (!math::isfinite(bracket.lower) || !math::isfinite(bracket.upper))
        return std::unexpected(error::non_finite);
    T lower = bracket.lower;
    T upper = bracket.upper;
    if (upper < lower) {
        lower = bracket.upper;
        upper = bracket.lower;
    }
    std::array<T, K> at_lower = detail::evaluate<K>(fn, lower);
    std::array<T, K> at_upper = detail::evaluate<K>(fn, upper);
    if (!math::isfinite(at_lower[0]) || !math::isfinite(at_upper[0])) return std::unexpected(error::non_finite);
    if (at_lower[0] == T(0))
        return Start<T, K>{ Root<T>{ lower, tolerance.at(lower) }, lower, at_lower, upper, at_upper };
    if (at_upper[0] == T(0))
        return Start<T, K>{ Root<T>{ upper, tolerance.at(upper) }, lower, at_lower, upper, at_upper };
    if (detail::same_sign(at_lower[0], at_upper[0])) return std::unexpected(error::no_sign_change);
    if (lower < T(0) && T(0) < upper) {
        T const zero = T(0);
        std::array<T, K> const at_zero = detail::evaluate<K>(fn, zero);
        if (at_zero[0] == T(0))
            return Start<T, K>{ Root<T>{ zero, tolerance.at(zero) }, lower, at_lower, upper, at_upper };
        if (math::isfinite(at_zero[0])) {
            if (detail::same_sign(at_zero[0], at_lower[0])) {
                lower = zero;
                at_lower = at_zero;
            } else {
                upper = zero;
                at_upper = at_zero;
            }
        }
    }
    return Start<T, K>{ std::nullopt, lower, at_lower, upper, at_upper };
}
```

In order:

1. **The ends must be finite**, or `non_finite`.
2. **Put them in order.** A reversed bracket is swapped, not rejected.
3. **Evaluate `fn` at both ends.** A value that isn't finite is `non_finite`. Only the value is checked, never the
   derivatives: an infinite `f'` at an end only makes a derivative step unusable there, and shouldn't fail the
   solve.
4. **An exact 0 at an end is the root**, returned with `tol(end)` as its bound.
5. **The same sign at both ends is `no_sign_change`.**
6. **Zero first.** If 0 lies strictly inside, `fn(0)` is evaluated. An exact 0 there is the root. Otherwise 0 replaces
   the end whose value has `fn(0)`'s sign, so the bracket is halved at 0, and the rest of the method works on a bracket
   of one sign.

**Why zero first.** Roots at 0 are common (odd functions, `x f(x)`). With the default `abs_tol = denorm_min`, a method
that doesn't hit 0 exactly would walk down to the smallest numbers, about 1100 halvings in double. Evaluating 0 first
costs one evaluation and finds such a root at once:

```
solve<Brent>(sin, [-1, 2])           eval 1: x = -1    f = -0.84
                                     eval 2: x = 2     f = +0.91
                                     eval 3: x = 0     f = 0         -> value 0, abs_error 4.9e-324 (tol(0))
```

When the root isn't at 0, the halving at 0 isn't wasted either. For `x - 0.25` on `[-1, 2]`, `fn(0) = -0.25` has
`fn(-1)`'s sign, so the bracket becomes `[0, 2]`, and Brent's first step, the secant through `(0, -0.25)` and
`(2, 1.75)`, lands on 0.25 exactly: 4 evaluations.

**Why it is skipped where `fn(0)` isn't finite.** Removable singularities at 0 are common: `x / (e^x - 1)` is NaN at
0 when evaluated naively, though its limit is 1. The caller didn't ask for `fn(0)`, so its failure must not fail the
solve. The bracket is then left as it is.

**What zero first doesn't help.** A sign change at 0 that `fn` doesn't evaluate as 0, such as a jump, is still
approached all the way down (the 1076 evaluations of section 1): `fn(0)` has one of the two signs, and the bracket
becomes `[-1, 0]` or `[0, 2]` with the sign change at its end.

---

## 7. Bisection (lines 887-918)

```cpp
template<std::size_t K, real T, typename F>
constexpr std::expected<Root<T>, error> bisection(F &fn, Start<T, K> const &s, Tolerance<T> const &tolerance)
{
    T a = s.lower;
    T fa = s.at_lower[0];
    T b = s.upper;
    while (true) {
        T const m = detail::midpoint(a, b);
        // The larger of the distances to the ends is a true bound even where half the width rounds, as among
        // the subnormals or with ends far apart. A midpoint not strictly inside means T's resolution.
        T const left = detail::distance_bound(m, a);
        T const right = detail::distance_bound(b, m);
        T const half = right < left ? left : right;
        if (!(a < m && m < b) || half <= tolerance.at(m)) return Root<T>{ m, half };
        T const fm = detail::evaluate<K>(fn, m)[0];
        if (!math::isfinite(fm)) return std::unexpected(error::non_finite);
        if (fm == T(0)) return Root<T>{ m, tolerance.at(m) };
        if (detail::same_sign(fm, fa)) {
            a = m;
            fa = fm;
        } else {
            b = m;
        }
    }
}
```

The reference method: each step halves the bracket, one bit per evaluation, whatever `fn` is.

- **The midpoint** is `detail::midpoint` (section 5), which can't overflow.
- **The bound** is the larger of the bounds on the distances from `m` to the ends. In exact arithmetic both are half
  the width. In floating point the halving can round, and with ends far apart (or among the subnormals, where the
  midpoint rounds) the two distances differ; the larger is a true bound on the distance from `m` to the sign change.
- **The stop test** is `half <= tol(m)`, which is the common test of section C (bracket at most `2 tol`) measured from
  the midpoint.
- **The midpoint is returned without being evaluated.** The sign change is somewhere in `[a, b]`, so it is within
  `half` of `m`, whatever `fn(m)` is. Evaluating it would cost one more call and add nothing to the guarantee.
- **A midpoint not strictly inside** means the ends are adjacent numbers of `T`: the bracket can't be split, and that
  is `T`'s resolution. The loop ends there even if the tolerance asks for more.
- **Only `fa` is kept.** The new point replaces the end whose value has its sign, and to know which, only one end's
  sign is needed.

On `x^2 - 2` over `[1, 2]`: 2 evaluations at the ends, then 50 midpoints, until the half-width `2^-51 = 4.4e-16` is
below `tol(m) = 2 eps * 1.414 = 6.3e-16`. The result is 1.4142135623730954, 1 ulp above the closest double to the
root, with `abs_error` 4.4e-16.

---

## 8. Brent's procedure, and Newton and Halley on a bracket (lines 920-1098)

`brent<M>` is Brent's procedure `zero` (Brent 1973, chapter 4), implemented from the book with three changes listed in
its doc comment. It serves three methods: for `M = Brent` each step is an interpolation; for `M = Newton` or `Halley`
it is the derivative step. Everything else, the bookkeeping, the safeguards and so the termination argument, is
Brent's for all three.

### The derivative steps (920-973)

```cpp
template<real T, std::size_t K>
constexpr std::optional<T> newton_step(std::array<T, K> const &v)
{
    static_assert(K >= 2);
    if (v[1] == T(0) || !math::isfinite(v[1])) return std::nullopt;
    std::optional<T> const ratio = detail::quotient(v[0], v[1]);
    if (!ratio) return std::nullopt;
    return T(-*ratio);
}

template<real T>
constexpr std::optional<T> halley_step(std::array<T, 3> const &v)
{
    std::optional<T> const newton = detail::newton_step(v);
    if (!newton) return std::nullopt;
    if (!math::isfinite(v[2])) return newton;
    std::optional<T> const curvature = detail::quotient(v[2], v[1]); // f'' / f'
    if (!curvature) return newton;
    std::optional<T> const u = detail::product(T(-*newton), T(*curvature / T(2))); // t f'' / (2 f')
    if (!u) return newton;
    T const denominator = T(T(1) - *u);
    if (denominator == T(0)) return newton;
    std::optional<T> const step = detail::quotient(*newton, denominator);
    if (!step || (*step > T(0)) != (*newton > T(0))) return newton;
    return step;
}

template<typename M, real T, std::size_t K>
constexpr std::optional<T> derivative_step(std::array<T, K> const &v)
{
    if constexpr (std::same_as<M, Halley>)
        return detail::halley_step(v);
    else
        return detail::newton_step(v);
}
```

Each returns `std::optional<T>`: the step, or nothing where it can't be formed.

- **Newton's step** is `-f / f'`. Nothing if `f'` is 0 or not finite, or the quotient isn't representable (`quotient`,
  section 5). `static_assert(K >= 2)` documents that it needs `f'`.
- **Halley's step** is `-t / (1 - u)` with `t = f / f'` and `u = t f'' / (2 f')`. The textbook writes it as
  `-2 f f' / (2 f'^2 - f f'')`; the form here avoids `f'^2`, which overflows sooner. Each intermediate is checked. Where
  any check fails, or `1 - u` is 0, or the step would point *against* Newton's (`u > 1`, which happens far from a
  root), it falls back on Newton's step, as Boost.Math's Halley does.
- **`derivative_step<M>`** chooses between them at compile time.

### Brent's variables

Brent's procedure keeps three points, and the code keeps his names:

| Name | Meaning                                                                            |
|------|------------------------------------------------------------------------------------|
| `b`  | the best estimate: the end of the bracket where `fn` is smaller in magnitude       |
| `c`  | the other end: `fn(b)` and `fn(c)` have opposite signs, so `[b, c]` is the bracket |
| `a`  | the previous `b` (it may equal `c`)                                                |
| `d`  | the step just taken                                                                |
| `e`  | the step before it                                                                 |

`va`, `vb`, `vc` are the arrays of values at those points (`va[0]` is `fn(a)`; for Newton `va[1]` is `f'(a)`).

### `brent<M>(fn, s, tolerance)` (975-1098)

```cpp
template<typename M, std::size_t K, real T, typename F>
constexpr std::expected<Root<T>, error> brent(F &fn, Start<T, K> const &s, Tolerance<T> const &tolerance)
{
    T const big = T(detail::largest<T>() / T(16));

    T a = s.lower;
    std::array<T, K> va = s.at_lower;
    T b = s.upper;
    std::array<T, K> vb = s.at_upper;
    T c = a;
    std::array<T, K> vc = va;
    T d = detail::bounded_difference(b, a, big);
    T e = d;
    while (true) {
        // [b, c] is the bracket, with fn(b) and fn(c) of opposite signs; a is the previous b.
        if (detail::same_sign(vb[0], vc[0])) {
            c = a;
            vc = va;
            d = detail::bounded_difference(b, a, big);
            e = d;
        }
        if (math::abs(vc[0]) < math::abs(vb[0])) {
            a = b;
            b = c;
            c = a;
            va = vb;
            vb = vc;
            vc = va;
        }
        T const fa = va[0];
        T const fb = vb[0];
        T const fc = vc[0];
        T const tol = tolerance.at(b);
        T const width = detail::distance_bound(c, b);
        if (width <= detail::twice(tol)) return Root<T>{ b, width };

        bool const wide = width > big || detail::distance_bound(b, a) > big || math::abs(e) > big;
        // Half the bracket, towards c: exact Brent below big, the overflow-free midpoint above it.
        T const m = wide ? T(detail::midpoint(b, c) - b) : T(T(c - b) / T(2));
        if (wide || math::abs(e) < tol || math::abs(fa) <= math::abs(fb)) {
            d = m;
            e = m;
        } else if constexpr (std::same_as<M, Brent>) {
            T const ratio = T(fb / fa);
            T p = T(0);
            T q = T(0);
            if (a == c) {
                // Linear interpolation (the secant).
                p = T(T(2) * m * ratio);
                q = T(T(1) - ratio);
            } else {
                // Inverse quadratic interpolation.
                T const qa = T(fa / fc);
                T const r = T(fb / fc);
                p = T(ratio * T(T(T(2) * m * qa * T(qa - r)) - T(T(b - a) * T(r - T(1)))));
                q = T(T(qa - T(1)) * T(r - T(1)) * T(ratio - T(1)));
            }
            if (p > T(0))
                q = T(-q);
            else
                p = T(-p);
            T const before_last = e; // Brent's s := e
            e = d;
            // Accept the interpolated point only if it falls well inside the bracket and the steps shrink. A
            // NaN or q = 0 fails both tests, and falls back on bisection.
            if (T(T(2) * p) < T(T(T(3) * m * q) - math::abs(T(tol * q)))
                && p < math::abs(T(before_last * q / T(2)))) {
                d = T(p / q);
            } else {
                d = m;
                e = m;
            }
        } else {
            std::optional<T> const step = detail::derivative_step<M>(vb);
            T const before_last = e;
            e = d;
            // Brent's acceptance for a step: towards c, well inside the bracket (so that b + step cannot round
            // onto c), and shorter than half the step before last.
            T const inside = T(T(T(3) * math::abs(m)) / T(2) - T(tol / T(2)));
            if (step && (*step > T(0)) == (m > T(0)) && math::abs(*step) < inside
                && math::abs(*step) < T(math::abs(before_last) / T(2))) {
                d = *step;
            } else {
                d = m;
                e = m;
            }
        }
        T next = b;
        if (math::abs(d) > tol)
            next = T(b + d);
        else if (m > T(0))
            next = T(b + tol);
        else
            next = T(b - tol);
        // A point not strictly inside the bracket means it is already at the tolerance's resolution: no step
        // can make progress, and evaluating there again would repeat the iteration forever.
        if (!((b < next && next < c) || (c < next && next < b))) return Root<T>{ b, width };
        a = b;
        va = vb;
        b = next;
        vb = detail::evaluate<K>(fn, b);
        if (!math::isfinite(vb[0])) return std::unexpected(error::non_finite);
        if (vb[0] == T(0)) return Root<T>{ b, tolerance.at(b) };
    }
}
```

**Setup.** `a` and `c` start at the lower end, `b` at the upper. `d` and `e` start as `b - a`, through
`bounded_difference`, which gives `max` instead of a difference that can't be formed.

**Each pass of the loop:**

1. **Restore the bracket (1009-1014).** If the last step left `fn(b)` and `fn(c)` with one sign, the sign change is now
   between `b` and the previous `b`, which is `a`. So `c = a`, and the step memory is reset.
2. **Make `b` the best point (1015-1022).** If `|fn(c)| < |fn(b)|`, `b` and `c` swap, and `a` becomes the old `b`.
3. **The stop test (1023-1028).** `width = distance_bound(c, b)`, and if it is at most `2 tol(b)` the result is `b` with
   `width` as its bound. Brent's book tests `|c - b| / 2 <= tol`; the doc comment explains the change: the bound is
   exact at every stop with the default tolerances, while the half may round in radix 10.
4. **The guard (1030-1032).** `big = max / 16`. While the bracket, the last step, or the step before is wider than
   `big`, the step is a bisection, and `m` (half the bracket, towards `c`) is formed by the overflow-free midpoint.
   Below `big`, `m = (c - b) / 2` exactly as Brent writes it. The doc comment gives the reason for 16: the
   interpolation forms `3 m q`, which reaches 12 `|c - b|`, and `2 p`, which reaches 4 (`|c - b| + |b - a|`). With
   every distance below `max / 16`, nothing overflows.
5. **Choose the step (1033-1080).** A bisection if the guard is on, if the step before last was already below `tol`
   (`|e| < tol`: interpolation has stopped helping), or if the last step didn't reduce `|f|` (`|fa| <= |fb|`).
   Otherwise:
    - **`M = Brent`:** with only two distinct points (`a == c`) the secant, otherwise inverse quadratic interpolation
      through the three points. Both are written as a fraction `p / q`, with `p >= 0` after the sign change, and the
      division is done only after the acceptance test. The ratios `fb / fa`, `fa / fc`, `fb / fc` can't overflow: in
      this branch `|fb| < |fa| <= |fc|`.
    - **`M = Newton, Halley`:** the derivative step from `b`.
6. **Accept the step, or bisect.** Brent's two conditions, for the interpolation written on `p` and `q`, for a
   derivative step written on the step `s` itself:
    - **Well inside the bracket:** `2p < 3 m q - |tol q|`, that is `|s| < 3 |m| / 2 - tol / 2`. The step goes at most
      three quarters of the way to `c`, less half a tolerance, so `b + s` can't round onto `c`. A derivative step must
      also point towards `c`; for the interpolation the first test implies it (with `p >= 0`, it forces `q` to have
      `m`'s sign).
    - **Shrinking:** `|s| < |e| / 2`, half the step *before last*. This is what bounds Brent's worst case: if steps
      stop shrinking, every second step is a bisection.
      A rejected step, a NaN or a `q` of 0 (both tests fail for those) becomes a bisection: `d = e = m`.
7. **Take the step, at least `tol` (1081-1087).** A step shorter than `tol` is replaced by `tol` towards `c`. Near the
   root, interpolation steps become tiny, and steps of `tol` make sure the bracket still shrinks: the last step lands
   just on the other side of the root, and the bracket is closed.
8. **Progress (1088-1090).** If the new point isn't strictly between `b` and `c`, the bracket is already at the
   tolerance's resolution (a huge tolerance, or ends that are adjacent numbers). Evaluating there would repeat the
   same iteration forever, so the result is returned as it stands. This came from an independent review (Phase 2).
9. **Evaluate (1091-1096).** `a` takes the old `b`, and `fn` is evaluated at the new `b`. Not finite is `non_finite`;
   exactly 0 is the root.

**C++ note: `else if constexpr` after a run-time `if`.** Lines 1033-1080 read

```cpp
if (wide || ...) { /* bisection */ }
else if constexpr (std::same_as<M, Brent>) { /* interpolation */ }
else { /* derivative step */ }
```

The first `if` is a run-time test; the `else` branch holds an `if constexpr`. That is legal: `else if constexpr (...)`
is just an `else` whose statement is an `if constexpr`. For `M = Brent` the derivative branch is discarded, and for
Newton the interpolation branch is; the run-time bisection test stays in both.

### Worked example 1: Brent on `x^2 - 2` over `[1, 2]`

The library's points, with the kind of step that produced each (from an instrumented copy of the header):

```
eval  x                       f(x)        step
1     1                       -1          prepare: an end
2     2                       +2          prepare: an end
3     1.3333333333333333      -0.22       secant (a == c: only two points so far)
4     1.4190476190476189      +1.4e-2     inverse quadratic
5     1.4140715109573241      -4.0e-4     secant (c was reset to a in step 1)
6     1.4142133199955025      -6.9e-7     secant
7     1.4142135623731364      +1.2e-13    inverse quadratic
8     1.4142135623730949      -4.4e-16    secant
9     1.4142135623730956      +1.3e-15    the minimum step: b + tol(b)
-> value 1.4142135623730949, abs_error 6.7e-16
```

- **Superlinear convergence.** The error goes 8e-2, 5e-3, 1.4e-4, 2.4e-7, 4e-14: the number of correct digits grows
  by a factor of about 1.6 to 2 per step.
- **The minimum step closes the bracket.** After eval 8, `b = 1.4142135623730949` is the closest double to the root,
  but the other end `c` is still eval 7's point, 4e-14 away. The interpolation now proposes a step of about 1.6e-16,
  shorter than `tol(b) = 6.3e-16`, so the step is `tol` towards `c`. Eval 9 lands 3 units of `T` above, with the other
  sign. Now `[b, c]` is 6.7e-16 wide, below `2 tol = 1.26e-15`: the stop test returns `b`, and the bound is exactly the
  width, by Sterbenz (section 5).

### Worked example 2: Newton on the same bracket

```
eval  x                       f(x)        step
1     1                       -1          prepare
2     2                       +2          prepare
3     1.5                     +0.25       bisection: Newton's step from b = 1 rejected
4     1.4166666666666667      +6.9e-3     Newton
5     1.4142156862745099      +6.0e-6     Newton
6     1.4142135623746899      +4.5e-12    Newton
7     1.4142135623730951      +4.4e-16    Newton
8     1.4142135623730945      -1.6e-15    the minimum step: b - tol(b)
-> value 1.4142135623730951, abs_error 6.7e-16
```

After the swap, `b = 1` (`|f| = 1 < 2`). Newton's step from 1 is `+0.5`, but the shrinking test requires it to be
below half of `e = b - a = 1`, and `0.5 < 0.5` is false, so the step is a bisection, which happens to reach the same
point 1.5. From then on every Newton step is accepted, and the quadratic convergence shows in the errors: 2.5e-3,
2.1e-6, 1.6e-12.

---

## 9. One step of a bracket: `enclose` (lines 1126-1290)

Toms748, the false-position family and Ridders share one function for "evaluate a point and keep the half with the
sign change": `enclose`, which is the paper's subroutine `bracket` (Alefeld, Potra and Shi 1995) with a tolerance of
our own. Each method only proposes points.

### `stop_at` (1126-1144)

```cpp
template<real T>
constexpr std::optional<Root<T>>
    stop_at(T const &a, T const &fa, T const &b, T const &fb, Tolerance<T> const &tolerance)
{
    T const u = math::abs(fa) < math::abs(fb) ? a : b;
    T const width = detail::distance_bound(a, b);
    if (width <= detail::twice(tolerance.at(u))) return Root<T>{ u, width };
    return std::nullopt;
}
```

The common stop test of section C, on a bracket `[a, b]` with its values: `u` is the end with the smaller `|f|`, and
if the bracket is at most `2 tol(u)` wide, `u` is the root.

### `Enclosure` and `Enclosed` (1146-1197)

```cpp
template<real T>
struct Enclosure
{
    T a;

    T fa;

    T b;

    T fb;

    T d;

    T fd;

    T e;

    T fe;
};

template<real T>
struct Enclosed
{
    std::optional<Root<T>> root;

    Enclosure<T> next;

    T c;

    T fc;
};
```

`Enclosure` is a method's state: the bracket `a < b` with its values, and the ends dropped by the last two steps, `d`
and then `e`. Toms748 interpolates through all four points; the others ignore `d` and `e`, or, in Ridders' case, use
`d` once (section 12). `Enclosed` is the outcome of one step: the new state, the point evaluated and its value, and a
root if the step ended the search.

### `wide` (1199-1207)

```cpp
template<real T>
constexpr bool wide(Enclosure<T> const &s)
{ return detail::distance_bound(s.a, s.b) > T(detail::largest<T>() / T(16)); }
```

Brent's guard of section 8, for these methods: while the bracket is wider than `max / 16`, every step is a bisection
and no other point is even formed. Every difference of points is then at most `max / 4`.

### `place(s, proposed, tolerance)` (1209-1245): where a point is evaluated

```cpp
template<real T>
constexpr T place(Enclosure<T> const &s, std::optional<T> const &proposed, Tolerance<T> const &tolerance)
{
    T const u = math::abs(s.fa) < math::abs(s.fb) ? s.a : s.b;
    T const t = tolerance.at(u);
    T delta = T(detail::largest<T>() / T(4));
    if (t <= T(detail::largest<T>() / T(2))) {
        delta = T(T(T(7) / T(10)) * t);
        T const tiny = T(std::numeric_limits<T>::denorm_min());
        if (delta < tiny) delta = tiny;
    }
    T const middle = detail::midpoint(s.a, s.b);
    if (!proposed || T(detail::distance_bound(s.a, s.b) / T(4)) <= delta) return middle;
    T c = *proposed;
    T const margin = T(T(2) * delta);
    if (c <= T(s.a + margin))
        c = T(s.a + margin);
    else if (c >= T(s.b - margin))
        c = T(s.b - margin);
    if (s.a < c && c < s.b) return c;
    return middle;
}
```

A method's proposed point isn't evaluated as it is. Two dangers: a point very close to an end gains almost nothing,
and a point whose rounding puts it on or outside an end breaks the bracket. The paper's rule, implemented here with
`delta = 0.7 tol(u)`:

- **A margin of `2 delta`.** A point closer than `2 delta` to an end is moved to `2 delta` inside. Then either the
  root lies in the sliver of width `2 delta = 1.4 tol` next to that end, and the next stop test ends the search, or
  the step removes at least `1.4 tol` of the bracket. Either way it makes progress.
- **The midpoint** where there is no proposed point, or where the bracket is at most `4 delta` wide (the two margins
  would overlap), or where the moved point still isn't strictly inside. The last case is real: at an end of much
  larger magnitude than `u`, adding `2 delta` can round back to the end itself.
- **`delta` is at least `denorm_min`.** In a type whose `denorm_min()` is `min()` (Boost's), `0.7 min` flushes to 0,
  and a margin of 0 would allow a point on an end. And where `tol(u)` exceeds `max / 2`, `delta` is `max / 4`, which
  can be formed.

The margin is visible at the end of every trace below: the last point is the closest double to the root moved
`2 delta`, about 4 units of `T`, inside.

### `enclose(fn, s, proposed, shift, tolerance)` (1247-1290)

```cpp
template<real T, typename F>
constexpr std::expected<Enclosed<T>, error> enclose(F &fn,
    Enclosure<T> const &s,
    std::optional<T> const &proposed,
    bool shift,
    Tolerance<T> const &tolerance)
{
    T const c = detail::wide(s) ? detail::midpoint(s.a, s.b) : detail::place(s, proposed, tolerance);
    // Unreachable: a midpoint not strictly inside means adjacent ends, which met the stop test before.
    assert(s.a < c && c < s.b);
    if (!(s.a < c && c < s.b))
        return Enclosed<T>{ Root<T>{ c, detail::distance_bound(s.a, s.b) }, s, c, s.fa };
    Enclosure<T> next = s;
    if (shift) {
        next.e = s.d;
        next.fe = s.fd;
    }
    T const fc = detail::evaluate<1>(fn, c)[0];
    if (!math::isfinite(fc)) return std::unexpected(error::non_finite);
    if (fc == T(0)) return Enclosed<T>{ Root<T>{ c, tolerance.at(c) }, next, c, fc };
    if (detail::same_sign(s.fa, fc)) {
        next.d = s.a;
        next.fd = s.fa;
        next.a = c;
        next.fa = fc;
    } else {
        next.d = s.b;
        next.fd = s.fb;
        next.b = c;
        next.fb = fc;
    }
    return Enclosed<T>{ detail::stop_at(next.a, next.fa, next.b, next.fb, tolerance), next, c, fc };
}
```

1. **The point:** the midpoint while the bracket is wide, otherwise `place`'s point.
2. **An `assert` that it is strictly inside.** The comment says why it can't fail: a midpoint not strictly inside
   means adjacent ends, and those met the stop test on the previous step. In a release build the `if` after it still
   returns safely, with the bracket as the bound.
3. **`shift`:** whether `e` takes `d`'s value before `d` records the newly dropped end. Toms748 needs it at some steps
   and not others, to keep the paper's assignments (section 10); the others pass `false`.
4. **Evaluate, then keep the half with the sign change.** The end replaced by `c` becomes the new `d`.
5. **Stop** on an exact zero, or on `stop_at`.

---

## 10. Toms748 (lines 1292-1501)

Alefeld, Potra and Shi's Algorithm 4.2, as printed in the paper (a user decision of Phase 4; Boost.Math's
`toms748_solve` arranges the first steps differently, and was used as a cross-check only). Each iteration takes up to
four points:

1. an **interpolation** with `k = 2`: the inverse cubic through `a, b, d, e` if it can be formed and lands inside,
   otherwise **Newton-Quadratic** with 2 Newton steps;
2. the same with `k = 3`;
3. a **double-length secant step** from the better end;
4. a **bisection**, but only if the iteration hasn't halved the bracket.

So each iteration at least halves the bracket, at a cost of at most 4 evaluations: the worst case is about 4
evaluations per halving, where Brent's grows with the square of the halvings.

### The interpolations (1292-1426)

```cpp
template<real T>
constexpr std::optional<T> secant_point(T const &a, T const &fa, T const &b, T const &fb)
{
    std::optional<T> const spread = detail::difference(fb, fa);
    if (!spread) return std::nullopt;
    std::optional<T> const r = detail::quotient(fa, *spread);
    if (!r) return std::nullopt;
    return T(a - T(*r * T(b - a))); // |r| < 1, and b - a at most max / 16
}
```

**`secant_point`:** the zero of the line through `(a, fa)` and `(b, fb)`, written as `a - r (b - a)` with
`r = fa / (fb - fa)`. With values of opposite signs, `|r| < 1`. The difference `fb - fa` of values of opposite signs
can overflow, so it goes through `difference`.

```cpp
template<real T>
constexpr std::optional<T> newton_quadratic(Enclosure<T> const &s, int k)
{
    std::optional<T> const secant = detail::secant_point(s.a, s.fa, s.b, s.fb);
    std::optional<T> const rise = detail::difference(s.fb, s.fa);
    std::optional<T> const rise_d = detail::difference(s.fd, s.fb);
    if (!rise || !rise_d) return secant;
    std::optional<T> const b1 = detail::quotient(*rise, T(s.b - s.a)); // f[a, b]
    std::optional<T> const bd = detail::quotient(*rise_d, T(s.d - s.b)); // f[b, d]
    if (!b1 || !bd) return secant;
    std::optional<T> const curve = detail::difference(*bd, *b1);
    if (!curve) return secant;
    std::optional<T> const a2 = detail::quotient(*curve, T(s.d - s.a)); // f[a, b, d]
    if (!a2 || *a2 == T(0)) return secant;
    T c = detail::same_sign(*a2, s.fa) ? s.a : s.b;
    for (int i = 0; i < k; ++i) {
        // P(c) = fa + (B + A (c - b)) (c - a), P'(c) = B + A ((c - a) + (c - b)): 2c - a - b would overflow
        // above max / 2.
        std::optional<T> const ab = detail::product(*a2, T(c - s.b));
        std::optional<T> const slope = ab ? detail::sum(b1, ab) : std::nullopt;
        std::optional<T> const term = slope ? detail::product(*slope, T(c - s.a)) : std::nullopt;
        std::optional<T> const value = term ? detail::sum(std::optional<T>(s.fa), term) : std::nullopt;
        std::optional<T> const along = detail::product(*a2, T(T(c - s.a) + T(c - s.b)));
        std::optional<T> const derivative = along ? detail::sum(b1, along) : std::nullopt;
        std::optional<T> const step =
            value && derivative ? detail::quotient(*value, *derivative) : std::nullopt;
        if (!step || !detail::sum_fits(c, T(-*step))) return secant;
        T const next = T(c - *step);
        if (!(s.a < next && next < s.b)) return secant;
        c = next;
    }
    return c;
}
```

**`newton_quadratic(s, k)`:** the quadratic through `(a, fa)`, `(b, fb)`, `(d, fd)` in Newton's form,

```
P(x) = fa + f[a, b] (x - a) + f[a, b, d] (x - a)(x - b)
```

with the divided differences `b1 = f[a, b]` and `a2 = f[a, b, d]`. Its zero isn't computed by the quadratic formula,
which needs a square root and cancels badly, but by `k` Newton steps on `P`. They start from the end where `P` is
convex towards the root (`a` if `a2` has `fa`'s sign, else `b`), where Newton's method converges monotonically. Every
quantity is checked, and anything that can't be formed returns the secant point. So does an iterate that leaves
`(a, b)`, which in exact arithmetic can't happen. `P'` is written as `b1 + a2 ((c - a) + (c - b))` rather than with
`2c - a - b`, which would overflow above `max / 2`.

```cpp
template<real T>
constexpr std::optional<T> inverse_cubic(Enclosure<T> const &s)
{
    auto const term = [](std::optional<T> const &dx, T const &f, std::optional<T> const &df) -> std::optional<T> {
        if (!dx || !df) return std::nullopt;
        std::optional<T> const r = detail::quotient(f, *df);
        if (!r) return std::nullopt;
        return detail::product(*dx, *r);
    };
    auto const minus = [](std::optional<T> const &x, std::optional<T> const &y) -> std::optional<T> {
        if (!x || !y) return std::nullopt;
        return detail::difference(*x, *y);
    };
    std::optional<T> const q11 = term(T(s.d - s.e), s.fd, detail::difference(s.fe, s.fd));
    std::optional<T> const q21 = term(T(s.b - s.d), s.fb, detail::difference(s.fd, s.fb));
    std::optional<T> const q31 = term(T(s.a - s.b), s.fa, detail::difference(s.fb, s.fa));
    std::optional<T> const d21 = term(T(s.b - s.d), s.fd, detail::difference(s.fd, s.fb));
    std::optional<T> const d31 = term(T(s.a - s.b), s.fb, detail::difference(s.fb, s.fa));
    std::optional<T> const q22 = term(minus(d21, q11), s.fb, detail::difference(s.fe, s.fb));
    std::optional<T> const q32 = term(minus(d31, q21), s.fa, detail::difference(s.fd, s.fa));
    std::optional<T> const d32 = term(minus(d31, q21), s.fd, detail::difference(s.fd, s.fa));
    std::optional<T> const q33 = term(minus(d32, q22), s.fa, detail::difference(s.fe, s.fa));
    std::optional<T> const c = detail::sum(std::optional<T>(s.a), detail::sum(detail::sum(q31, q32), q33));
    if (!c || !(s.a < *c && *c < s.b)) return std::nullopt;
    return c;
}
```

**`inverse_cubic(s)`:** the paper's `ipzero`. Instead of fitting `f` as a function of `x` and solving for its zero,
it fits `x` as a function of `f` through the four points `(fa, a)`, `(fb, b)`, `(fd, d)`, `(fe, e)` and evaluates that
cubic at `f = 0`, by the Aitken-Neville scheme. No equation is solved; it is one evaluation of an interpolating
polynomial. Each term is "a distance times a ratio of values", and the lambda `term` forms the ratio first and the
product second, both checked. A point that isn't strictly inside `(a, b)` is rejected.

```cpp
template<real T>
constexpr std::optional<T> interpolate(Enclosure<T> const &s, int k, bool quadratic)
{
    bool const distinct = s.fa != s.fb && s.fa != s.fd && s.fa != s.fe && s.fb != s.fd && s.fb != s.fe
                          && s.fd != s.fe;
    if (!quadratic && distinct) {
        std::optional<T> const c = detail::inverse_cubic(s);
        if (c) return c;
    }
    return detail::newton_quadratic(s, k);
}
```

**`interpolate(s, k, quadratic)`:** the inverse cubic needs four *distinct* values (otherwise a division by 0 in
the scheme), and the first iteration doesn't have `e` yet (`quadratic = true`). Otherwise it falls back on
Newton-Quadratic.

```cpp
template<real T>
constexpr std::optional<T> double_secant(Enclosure<T> const &s)
{
    bool const at_a = math::abs(s.fa) < math::abs(s.fb);
    T const u = at_a ? s.a : s.b;
    T const fu = at_a ? s.fa : s.fb;
    std::optional<T> const spread = detail::difference(s.fb, s.fa);
    if (!spread) return std::nullopt;
    std::optional<T> const r = detail::quotient(fu, *spread); // |r| <= 1/2: fa and fb have opposite signs
    if (!r) return std::nullopt;
    T const width = T(s.b - s.a);
    T const c = T(u - T(T(T(2) * *r) * width));
    if (math::abs(T(c - u)) > T(width / T(2))) return std::nullopt;
    return c;
}
```

**`double_secant(s)`:** step 3 of the iteration. From the end `u` with the smaller `|f|`, take *twice* the secant step.
If the secant is accurate, this overshoots the root on purpose, so the root ends up near the end of the new bracket
that will be retained, and the far end moves. That is what makes the bracket shrink from both sides, where a plain
secant method converges from one side only. If the step would move more than half the bracket, it returns nothing,
and `enclose` takes the midpoint.

### `advance` and `toms748` (1428-1501)

```cpp
template<real T, typename F, typename P>
constexpr std::expected<Enclosed<T>, error>
    advance(F &fn, Enclosure<T> const &s, P const &point, bool shift, Tolerance<T> const &tolerance)
{
    std::optional<T> const proposed = detail::wide(s) ? std::optional<T>() : std::optional<T>(point(s));
    return detail::enclose(fn, s, proposed, shift, tolerance);
}
```

A step whose point is computed only when the bracket isn't wide, so that no interpolation runs on a bracket where its
arithmetic could overflow. `point` is a lambda that computes the proposal from the state.

**Why a function and not a lambda inside `toms748`?** A lambda capturing `fn` by reference fails on MSVC in constant
evaluation when `fn` is a function, the same problem as `Joined` (section 3). A function that takes `fn` as a
parameter works everywhere.

```cpp
template<real T, typename F>
constexpr std::expected<Root<T>, error> toms748(F &fn, Start<T, 1> const &start, Tolerance<T> const &tolerance)
{
    T const lower = start.lower;
    T const upper = start.upper;
    T const f_lower = start.at_lower[0];
    T const f_upper = start.at_upper[0];
    if (std::optional<Root<T>> const r = detail::stop_at(lower, f_lower, upper, f_upper, tolerance)) return *r;
    // d and e are read only once steps have set them (the first iteration uses Newton-Quadratic, without e).
    Enclosure<T> s{ lower, f_lower, upper, f_upper, upper, f_upper, upper, f_upper };
    std::expected<Enclosed<T>, error> r = detail::advance(
        fn, s, [](Enclosure<T> const &e) { return detail::secant_point(e.a, e.fa, e.b, e.fb); }, false,
        tolerance); // 4.2.1-4.2.2
    for (bool first = true;; first = false) {
        if (!r) return std::unexpected(r.error());
        if (r->root) return *r->root;
        s = r->next;
        T const a0 = s.a;
        T const b0 = s.b;
        // 4.2.3-4.2.4: e := d before the bracket.
        r = detail::advance(
            fn, s, [first](Enclosure<T> const &e) { return detail::interpolate(e, 2, first); }, true, tolerance);
        if (!r) return std::unexpected(r.error());
        if (r->root) return *r->root;
        // 4.2.5-4.2.6.
        r = detail::advance(
            fn, r->next, [](Enclosure<T> const &e) { return detail::interpolate(e, 3, false); }, false,
            tolerance);
        if (!r) return std::unexpected(r.error());
        if (r->root) return *r->root;
        // 4.2.7-4.2.10, e := d-bar (4.1.9).
        r = detail::advance(
            fn, r->next, [](Enclosure<T> const &e) { return detail::double_secant(e); }, true, tolerance);
        if (!r) return std::unexpected(r.error());
        if (r->root) return *r->root;
        // 4.2.11: a bisection where the iteration has not halved the bracket, e := d-hat.
        if (detail::distance_bound(r->next.a, r->next.b) < T(detail::distance_bound(a0, b0) / T(2))) continue;
        r = detail::advance(
            fn, r->next, [](Enclosure<T> const &e) { return std::optional<T>(detail::midpoint(e.a, e.b)); },
            true, tolerance);
    }
}
```

The loop follows the paper's step numbers, written in the comments:

- **4.2.1-4.2.2:** one secant step on the initial bracket. Before it, `d` and `e` are set to `upper` as placeholders;
  the first interpolation uses Newton-Quadratic only, which doesn't read `e`.
- **Each iteration:** remember the bracket `[a0, b0]`, take the two interpolation steps, then the double secant, and
  if the bracket isn't below half of `b0 - a0`, a bisection.
- **`shift`:** the interpolation with `k = 2`, the double secant and the bisection pass `true`. The comments name the
  paper's assignments they implement (`e := d` at 4.2.4, `e := d-bar` and `e := d-hat`), which keep `e` the dropped
  end the next inverse cubic needs.

**C++ note: an `if` with a declaration.** `if (std::optional<Root<T>> const r = detail::stop_at(...)) return *r;`
declares `r` inside the condition. The condition is `r` converted to `bool`, so it is true if `r` holds a value, and
`r` exists only in the `if` statement. C++17 also allows `if (init; condition)`; this older form is enough when the
variable itself is the condition.

**C++ note: `for (bool first = true;; first = false)`.** A loop without a condition, run until a `return`, whose
variable is `true` on the first pass and `false` afterwards. The lambda `[first](...)` captures it by value, so the
first interpolation knows it may use only the quadratic.

### Worked example: Toms748 on `x^2 - 2` over `[1, 2]`

```
eval  x                       f(x)        step
1     1                       -1          prepare
2     2                       +2          prepare
3     1.3333333333333333      -0.22       4.2.1: secant
                                          iteration 1, bracket [1.333, 2]
4     1.4166666666666665      +6.9e-3     Newton-Quadratic, k = 2 (no e yet)
5     1.4142241834420135      +3.0e-5     inverse cubic
6     1.4142023159964394      -3.2e-5     double secant: overshoots, the root is now bracketed tightly
                                          bracket [1.41420, 1.41422], below half: no bisection
                                          iteration 2
7     1.4142135623731005      +1.6e-14    inverse cubic
8     1.4142135623730949      -4.4e-16    inverse cubic
9     1.4142135623730958      +2.2e-15    double secant, moved 2 delta inside a by place
-> value 1.4142135623730949, abs_error 8.9e-16
```

Eval 6 shows the double secant at work: from `u = 1.41422` (eval 5) twice the secant step lands *below* the root, so
both ends of the bracket are now within 1.1e-5 of it. On `x^20 - 1` over `[0, 5]`, where the start is poor, the
bisection step of 4.2.11 runs twice (21 evaluations; Brent takes 19).

---

## 11. False position: AndersonBjorck, Illinois, RegulaFalsi (lines 1503-1575)

**Regula falsi** ("false position") takes the zero of the secant through the bracket's ends, and keeps the half with
the sign change. On a convex function, one end never moves: every new point lands on the same side of the root,
and the bracket shrinks from one side only, linearly and often very slowly. The two cures differ in one factor:
when the same end is **retained** again, its value is scaled down, so that the next secant is pulled towards it.

| Method           | Factor on the retained end's value, when retained again |
|------------------|---------------------------------------------------------|
| `RegulaFalsi`    | 1 (none)                                                |
| `Illinois`       | 1/2                                                     |
| `AndersonBjorck` | `m = 1 - f(c) / f(b)`, or 1/2 where `m` isn't positive  |

```cpp
template<typename M, real T, typename F>
constexpr std::expected<Root<T>, error>
    false_position(F &fn, Start<T, 1> const &start, Tolerance<T> const &tolerance)
{
    T retained = start.lower;
    T f_retained = start.at_lower[0];
    T scaled = f_retained;
    T latest = start.upper;
    T f_latest = start.at_upper[0];
    if (std::optional<Root<T>> const r = detail::stop_at(retained, f_retained, latest, f_latest, tolerance))
        return *r;
    T reference = detail::distance_bound(retained, latest);
    int unhalved = 0;
    while (true) {
        bool const ordered = retained < latest;
        Enclosure<T> const s{ ordered ? retained : latest, ordered ? f_retained : f_latest,
            ordered ? latest : retained, ordered ? f_latest : f_retained, retained, f_retained, retained,
            f_retained };
        bool const bisect = unhalved >= 3;
        std::optional<T> proposed;
        if (!detail::wide(s)) {
            if (bisect) {
                proposed = detail::midpoint(s.a, s.b);
            } else if (std::optional<T> const spread = detail::difference(f_latest, scaled)) {
                // f_latest and scaled have opposite signs, so 0 < r < 1.
                if (std::optional<T> const r = detail::quotient(f_latest, *spread))
                    proposed = T(latest - T(*r * T(latest - retained)));
            }
        }
        std::expected<Enclosed<T>, error> const step = detail::enclose(fn, s, proposed, false, tolerance);
        if (!step) return std::unexpected(step.error());
        if (step->root) return *step->root;
        if (!detail::same_sign(step->fc, f_latest)) {
            retained = latest;
            f_retained = f_latest;
            scaled = f_latest;
        } else if constexpr (std::same_as<M, AndersonBjorck>) {
            T const m = math::abs(step->fc) < math::abs(f_latest) ? T(T(1) - T(step->fc / f_latest))
                                                                  : T(T(1) / T(2));
            scaled = T(scaled * m);
        } else if constexpr (std::same_as<M, Illinois>) {
            scaled = T(scaled * T(T(1) / T(2)));
        } // RegulaFalsi: the factor 1
        latest = step->c;
        f_latest = step->fc;
        T const width = detail::distance_bound(retained < latest ? retained : latest,
            retained < latest ? latest : retained);
        if (bisect || width <= T(reference / T(2))) {
            reference = width;
            unhalved = 0;
        } else {
            ++unhalved;
        }
    }
}
```

**The state** isn't an `Enclosure` with `a < b` but the method's own view: `retained` (the end that stays) with its
true value `f_retained` and its scaled value `scaled`, and `latest`, the newest point. Each pass builds the `Enclosure`
that `enclose` needs, in order (`d` and `e` are filled with the retained end and not used).

**The point** is the secant's zero through `(latest, f_latest)` and `(retained, scaled)`:
`latest - r (latest - retained)` with `r = f_latest / (f_latest - scaled)`. The two values have opposite signs, so
`0 < r < 1` and the point is inside. Both the difference and the quotient are checked; if either fails, `proposed`
stays empty and `enclose` takes the midpoint.

**After the step (1553-1565):**

- If `f(c)` has the other sign than `f_latest`, the sign change is between `latest` and `c`: `latest` becomes the
  retained end, with its true value, and nothing is scaled.
- Otherwise the old retained end is retained again, and `scaled` is multiplied by the method's factor. Anderson and
  Bjorck's `m = 1 - f(c) / f(b)` is formed only where `|f(c)| < |f(b)|`, so the quotient is below 1 and nothing
  overflows. In that case `m` is in `(0, 1)`. Otherwise `m` would be at most 0, and 1/2 is used.
- `c` becomes `latest`.

The comment adds an observation: nothing depends on the factor being below 1. `scaled` keeps its sign, so `r` stays in
`(0, 1)` whatever the factor.

**The guard (1566-1573).** The published methods have no guarantee: on `x^20 - 1` over `[0, 5]` Anderson-Bjorck as
published took 65 million evaluations, regula falsi over 200,000. So a bisection is forced after three consecutive
steps that haven't brought the bracket below half its width at the last halving (`reference`). A guard bisection
counts as a halving. The worst case is then about 4 evaluations per halving, as for Toms748.

**C++ note: `if constexpr` with no `else`.** The chain `else if constexpr (AndersonBjorck) ... else if constexpr
(Illinois) ...` has no final branch: for `RegulaFalsi` nothing happens, which is the factor 1, as the comment says.

### Worked example 1: regula falsi stalls, and the guard

`RegulaFalsi` on `x^2 - 2` over `[1, 2]`, 20 evaluations:

```
eval  x                       f(x)        step
3     1.3333333333333335      -0.22       secant      (2 retained)
4     1.4000000000000001      -4.0e-2     secant      (2 retained)
5     1.411764705882353       -6.9e-3     secant      (2 retained: three steps without halving)
6     1.7058823529411766      +0.91       guard bisection of [1.4118, 2]
7     1.4139844617092119      -6.5e-4     secant      (1.7059 retained)
8     1.4141921443070291      -6.1e-5     secant
9     1.4142115601832654      -5.7e-6     secant
10    1.560046956562221       +0.43       guard bisection
...   (the same pattern at 11-14 and 15-18)
19    1.4142135623730949      -4.4e-16    secant
20    1.4142135623730958      +2.2e-15    secant, moved 2 delta inside by place
```

Every secant point lands below the root, on the convex side: the left end converges, linearly (each step gains a
factor of about 6 to 10), while the right end stays where it is. Only the guard moves it. Illinois takes 11
evaluations on the same problem, Anderson-Bjorck 9, the same as Brent.

### Worked example 2: where Anderson-Bjorck wins

`ln(x / 1e-12)` on `[1e-20, 1]` is strongly curved near the left end: the root is at 1e-12, and on a linear scale
almost all of the bracket has large positive values. Brent takes 45 evaluations, Anderson-Bjorck 18. Its points go
0.4, 9e-3, 3.4e-5, 4.1e-8, 3.2e-11, 5.1e-14: each secant step moves several decades towards the root, because the
rescaling keeps pulling the secant towards the retained end at 1e-20. Then it converges superlinearly: evals 15 to 17
have errors 7.5e-18, 2.4e-20 and 9e-26, and eval 18 is the double nearest 1e-12, where `fn` is exactly 0.

---

## 12. Ridders' method (lines 1577-1699)

**What it computes** (Ridders 1979). Given a bracket `[x1, x2]` and its midpoint `x3`, find the factor `e^Q` for which
`f(x) e^(Q x)` is a straight line through the three points, and take that line's zero. The result is

```
x4 = x3 + (x3 - x1) sign(f1 - f2) f3 / sqrt(f3^2 - f1 f2)
```

It needs a square root, which `nxx::real` doesn't provide, and in this form it loses roots near an end. The header
solves the first with `detail::sqrt` and the second with a rearranged formula.

### `sqrt(x)` (1577-1633)

```cpp
template<real T>
constexpr T sqrt(T const &x)
{
    assert(!(x < T(0)));
    if (!(x > T(0)) || !math::isfinite(x)) return x;
    T const r = static_cast<T>(std::numeric_limits<T>::radix);
    T const top = detail::largest<T>();
    T s = x;
    T scale = T(1);
    // A type too narrow to hold radix^2 (a Boost type whose max is 3.98) is not scaled, and Newton starts at
    // max(1, s), above the root (independent review of Phase 6: radix^2 overflowed, and sqrt 2 came out 1.9e-6).
    bool const scalable = r <= T(top / r);
    T const r2 = scalable ? T(r * r) : r;
    // p = q^2 is squared up while it fits and keeps s on its side of the interval; s is scaled by p, scale by q.
    while (scalable && s > T(1)) {
        T p = r2;
        T q = r;
        while (p <= T(top / p) && T(s / p) > p) {
            p = T(p * p);
            q = T(q * q);
        }
        s = T(s / p);
        scale = T(scale * q);
    }
    while (scalable && s < T(T(1) / r2)) {
        T p = r2;
        T q = r;
        while (p <= T(top / p) && T(s * p) < T(T(1) / p)) {
            p = T(p * p);
            q = T(q * q);
        }
        s = T(s * p);
        scale = T(scale / q);
    }
    T y = s > T(1) ? s : T(1);
    while (true) {
        T const next = T(y + T(T(T(s / y) - y) / T(2)));
        if (!(next < y)) break;
        y = next;
    }
    return T(y * scale);
}
```

A square root from `+ - * /` only, the same at compile time and run time, **faithful**: within one unit in the last
place of the true root. (Correctly rounded would need the exact residual `x - y^2`, which Dekker's splitting gives in
radix 2 but not for every `real`.)

1. **0, inf and NaN are returned as they are**, and so is a negative `x`, which is a precondition, asserted. Returning
   rather than looping is what lets the loops below assume a positive finite `s`.
2. **Scale `x` into `[1/radix^2, 1]`** by even powers of the radix: `s = x / radix^(2k)` and `scale = radix^k`, both
   exact, since multiplying by a power of the radix changes only the exponent. Then `sqrt(x) = sqrt(s) * scale`.
   Dividing by `radix^2` once per step would take as many steps as the exponent is large, about 10^9 in
   `cpp_bin_float_50`. So the inner loops square `p = q^2` up (radix^2, radix^4, radix^8, ...) while it fits and
   doesn't overshoot, and the number of steps is logarithmic in the exponent.
3. **Newton's iteration for `y^2 = s`** from `y = 1`, written as a correction: `y + (s / y - y) / 2`. Started above
   the root, the iterates decrease monotonically (by the inequality of the means), so the loop ends when an iterate
   stops decreasing: rounding has reached the root. The comment explains the form. The textbook `(y + s / y) / 2`
   is up to 3 ulps off in radix 10, where `y + s / y` lands in the next decade and rounds coarsely; near the root
   `s / y - y` is small and exact, and the correction form stays within one ulp (checked exhaustively in decimal at 4
   to 6 digits).
4. **Unscale**, exactly.

**The narrow type.** The Phase 6 review found a Boost type with 8 bits whose largest value is 3.98, where `radix^2 =
4` overflows. `scalable` is false for such a type: it isn't scaled at all, and Newton starts at `max(s, 1)`, which is
still above the root. Before the fix, `sqrt(2)` in that type came out 1.9e-6.

**An example, in double.** `x = 2`: `p = 4` doesn't fit twice (`2 / 4 > 4` is false), so `s = 0.5`, `scale = 2`.
Newton from 1: 0.75, 0.70833, 0.7071078, 0.70710678119, ..., and `* 2` gives 1.4142135623730949. That is one ulp below
the correctly rounded 1.4142135623730951 (`std::sqrt`): faithful, as promised, not correctly rounded.

### `ridders_point` (1635-1667)

```cpp
template<real T>
constexpr std::optional<T> ridders_point(Enclosure<T> const &s, T const &m, T const &fm)
{
    bool const m_is_a = s.a == m;
    T const e = m_is_a ? s.b : s.a;
    T const fe = m_is_a ? s.fb : s.fa;
    T scale = math::abs(fe);
    if (math::abs(s.fd) > scale) scale = math::abs(s.fd);
    if (math::abs(fm) > scale) scale = math::abs(fm);
    T const ge = T(fe / scale);
    T const gd = T(s.fd / scale);
    T const gm = T(fm / scale);
    T const root = detail::sqrt(T(T(gm * gm) - T(ge * gd)));
    std::optional<T> const w = detail::quotient(T(-T(ge * gd)), T(root * T(root + math::abs(gm))));
    if (!w) return std::nullopt;
    return T(e + T(*w * T(m - e)));
}
```

It is called with the state *after* `enclose` evaluated the midpoint `m`. In that state the bracket has `m` at one end
and the old end `e` with the other sign at the other, and the end the midpoint step dropped is `d`. So the old bracket
was `{e, d}`, and the formula becomes

```
x4 = e + w (m - e),   w = 1 - |f(m)| / root = -g_e g_d / (root (root + |g_m|))
```

with `root = sqrt(g_m^2 - g_e g_d)`, and `g = f / s` for `s = max(|f(e)|, |f(d)|, |f(m)|)`.

- **Why `e + w (m - e)`.** The textbook form moves from `m` towards `e` by `|m - e| |f(m)| / root`, so `x4` is `e` plus
  the rest, `(1 - |f(m)| / root) (m - e)`.
- **Why the second form of `w`.** `1 - |f(m)| / root` cancels: where `f(e) f(d)` is below `eps` of `f(m)^2`, the
  quotient rounds to 1 and `w` to 0, which puts the point exactly on `e`. Multiplying by `(root + |g_m|) / (root +
  |g_m|)` gives `(root^2 - g_m^2) / (root (root + |g_m|))`, and `root^2 - g_m^2 = -g_e g_d`, which is positive (`e` and
  `d` have values of opposite signs) and computed without cancellation. The comment gives the measured effect: on
  `x - 1e-300` over `[0, 1]` the textbook form bisected through 1922 evaluations; this form takes 4.
- **Why scale by `s`.** `f(m)^2` overflows for `|f(m)|` above `sqrt(max)`. Divided by the largest of the three
  magnitudes, every `g` is at most 1, so the squares are at most 1 and `root` is in `[|g_m|, sqrt 2]`.
- **Nothing** where `w` can't be formed: `root` underflows, and the quotient check fails. `enclose` then takes the
  midpoint.

### `ridders` (1669-1699)

```cpp
template<real T, typename F>
constexpr std::expected<Root<T>, error> ridders(F &fn, Start<T, 1> const &start, Tolerance<T> const &tolerance)
{
    T const lower = start.lower;
    T const upper = start.upper;
    T const f_lower = start.at_lower[0];
    T const f_upper = start.at_upper[0];
    if (std::optional<Root<T>> const r = detail::stop_at(lower, f_lower, upper, f_upper, tolerance)) return *r;
    Enclosure<T> s{ lower, f_lower, upper, f_upper, upper, f_upper, upper, f_upper };
    while (true) {
        // The midpoint: enclose's point where none is proposed.
        std::expected<Enclosed<T>, error> const half = detail::enclose(fn, s, std::optional<T>(), false, tolerance);
        if (!half) return std::unexpected(half.error());
        if (half->root) return *half->root;
        std::optional<T> const proposed =
            detail::wide(half->next) ? std::optional<T>() : detail::ridders_point(half->next, half->c, half->fc);
        std::expected<Enclosed<T>, error> const fit =
            detail::enclose(fn, half->next, proposed, false, tolerance);
        if (!fit) return std::unexpected(fit.error());
        if (fit->root) return *fit->root;
        s = fit->next;
    }
}
```

Each iteration is two `enclose` steps: the midpoint (no proposed point), then Ridders' point from the state the
midpoint step left. So every iteration at least halves the bracket, for at most 2 evaluations: at most twice
Bisection's count, guaranteed, which is the method's selling point. And near a simple root the fitted point converges
quadratically per iteration.

### Worked example: Ridders on `x^2 - 2` over `[1, 2]`, and the cancellation

```
eval  x                       f(x)        step
3     1.5                     +0.25       midpoint
4     1.4129611720221511      -3.5e-3     Ridders' point
5     1.4564805860110757      +0.12       midpoint
6     1.4142132868484725      -7.8e-7     Ridders' point
7     1.4353469364297742      +6.0e-2     midpoint
8     1.4142135623579404      -4.3e-11    Ridders' point
9     1.4247802493938573      +3.0e-2     midpoint
10    1.4142135623730949      -4.4e-16    Ridders' point
11    1.4194969058834761      +1.5e-2     midpoint
12    1.4142135623730958      +2.2e-15    Ridders' point, moved 2 delta inside by place
-> value 1.4142135623730949, abs_error 8.9e-16
```

The fitted points converge fast (errors 1.3e-3, 2.8e-7, 1.5e-11, 1e-16), but every one lands below the root, so the
upper end only moves through the midpoints, which halve the distance each time. Half of the 12 evaluations are spent
on the upper end. That is why Ridders costs more than Brent in general (416 against 309 evaluations on 27 of the
paper's cases), and why it is a teaching example.

On `x - 1e-300` over `[0, 1]`: the midpoint 0.5, then Ridders' point is exactly 1e-300, where `fn` is 0. 4
evaluations, against the textbook formula's 1922.

---

## 13. The open methods (lines 1701-1895)

### `hand_over` (1701-1723)

```cpp
template<typename M, std::size_t K, real T, typename F>
constexpr std::expected<Root<T>, error> hand_over(F &fn,
    T const &x,
    std::array<T, K> const &vx,
    T const &y,
    std::array<T, K> const &vy,
    Tolerance<T> const &tolerance)
{
    if (x < y) return detail::brent<M>(fn, Start<T, K>{ std::nullopt, x, vx, y, vy }, tolerance);
    return detail::brent<M>(fn, Start<T, K>{ std::nullopt, y, vy, x, vx }, tolerance);
}
```

When an open method finds two consecutive iterates whose values have opposite signs, it hands them to the bracket loop
of section 8: `brent<M>` for Newton and Halley (which then keep taking their derivative step, safeguarded), and
`brent<Brent>` for Secant and Steffensen. The two points are put in order and wrapped in a `Start`, with the values
already known. There is no zero first here: they are adjacent iterates, not a bracket the caller chose.

### `open_derivative<M>` (1725-1763): Newton and Halley from a guess

```cpp
template<typename M, real T, typename F>
constexpr std::expected<Root<T>, error> open_derivative(F &fn, T const &guess, Tolerance<T> const &tolerance)
{
    constexpr std::size_t K = values_needed<M>;
    if (!math::isfinite(guess)) return std::unexpected(error::non_finite);
    T x = guess;
    std::array<T, K> vx = detail::evaluate<K>(fn, x);
    if (!math::isfinite(vx[0])) return std::unexpected(error::non_finite);
    if (vx[0] == T(0)) return Root<T>{ x, tolerance.at(x) };
    for (int i = 0; i < detail::iteration_cap<T>(); ++i) {
        std::optional<T> const step = detail::derivative_step<M>(vx);
        if (!step) return std::unexpected(error::not_converged);
        T const tol = tolerance.at(x);
        bool const upward = !detail::same_sign(vx[0], vx[1]); // Newton's direction, -f / f'
        T const move = math::abs(*step) > tol ? *step : (upward ? tol : T(-tol));
        if (!detail::sum_fits(x, move)) return std::unexpected(error::not_converged);
        T const next = T(x + move);
        std::array<T, K> const vn = detail::evaluate<K>(fn, next);
        if (!math::isfinite(vn[0])) return std::unexpected(error::non_finite);
        if (vn[0] == T(0)) return Root<T>{ next, tolerance.at(next) };
        if (!detail::same_sign(vn[0], vx[0])) return detail::hand_over<M>(fn, x, vx, next, vn, tolerance);
        x = next;
        vx = vn;
    }
    return std::unexpected(error::not_converged);
}
```

1. **The guess must be finite**, and `fn` there finite; an exact 0 is the root.
2. **Each iteration:** the step (nothing: `not_converged`), at least `tol(x)` long, and `x + move` checked with
   `sum_fits`. Evaluate. Not finite: `non_finite`. Exactly 0: the root. Opposite signs to `fn(x)`: **hand over**.
3. **The cap** of section 4 ends a run that never crosses a sign change: `not_converged`.

**The minimum step's direction comes from the signs.** Newton's step `-f / f'` points up when `f` and `f'` have
opposite signs. The step itself may have underflowed to 0, which has no direction, so `upward` is computed from the
signs, not from the step.

**Why the minimum step at all.** Near a root the steps become tiny, and without it an iterate could land on the same
side every time. Steps of at least `tol` guarantee that the iterates cross the root once they are within `tol` of it,
which produces the sign change, and the hand-over then meets the stop test at once.

### Worked example: Newton and Halley from 1, on `x^2 - 2`

```
Newton (f, f')                                  Halley (f, f', f'')
eval  x                    f(x)                 eval  x                    f(x)
1     1                    -1                   1     1                    -1
2     1.5                  +0.25  sign change   2     1.3999999999999999   -4.0e-2
      -> hand over: brent<Newton> on [1, 1.5]   3     1.4142131979695431   -1.0e-6
3     1.4166666666666667   +6.9e-3              4     1.4142135623730949   -4.4e-16
4     1.4142156862745099   +6.0e-6              5     1.4142135623730956   +1.3e-15  minimum step:
5     1.4142135623746899   +4.5e-12                   sign change, hand over: [b, c] already
6     1.4142135623730951   +4.4e-16                   meets the stop test
7     1.4142135623730945   -1.6e-15  minimum step
```

Newton's first step overshoots, and the sign change appears at once. From there it is `brent<Newton>`: every Newton
step is accepted. Halley's steps come from below and never cross the root; its fourth point is the closest double, and
only the minimum step of `tol` produces the sign change. The doc comment of `Newton` calls this "about one
evaluation more than the textbook stop": the price of a verified sign change.

### `secant_second_point` (1765-1782)

```cpp
template<real T>
constexpr T secant_second_point(T const &x0, Tolerance<T> const &tolerance)
{
    T const eta = math::root_epsilon<4, T>();
    if (x0 == T(0)) return eta;
    if (math::abs(x0) > T(detail::largest<T>() / T(T(1) + eta))) return T(x0 * T(T(1) - eta));
    T const second = T(x0 * T(T(1) + eta));
    if (second == x0) return T(x0 + tolerance.at(x0)); // tiny x0: tol(x0) is about denorm_min
    return second;
}
```

Secant needs two points to start. The second is `x0 (1 + eta)` with `eta = eps^(1/4)`, about 1.2e-4 in double: a
relative perturbation, so no scale is assumed (a fixed `x0 + 1e-4` would be huge for `x0 = 1e-10` and lost for `x0 =
1e10`). At `x0 = 0` it is `eta` itself. Near `max` it goes down instead, so that it can be formed. Among the
subnormals, where `x0 (1 + eta)` rounds back to `x0`, it is `x0 + tol(x0)`.

`eps^(1/4)` is `math::root_epsilon<4, T>()`, close to SciPy's 1e-4 in double and tied to `T`'s own precision in
other types. Why relative: the first proposal, `x0 + eta max(|x0|, 1)`, left `ln(x / 1e-12)` from 1e-13 for
1.2e-4, seven decades away (`PHASE2_API.md` C).

### `secant_move` (1784-1810): the step shared by Secant and Steffensen

```cpp
template<real T>
constexpr std::optional<T>
    secant_move(T const &x, T const &f, T const &dx, T const &difference, Tolerance<T> const &tolerance)
{
    if (difference == T(0)) return std::nullopt;
    std::optional<T> const ratio = detail::quotient(f, difference);
    if (!ratio) return std::nullopt;
    std::optional<T> const step = detail::product(dx, *ratio);
    if (!step) return std::nullopt;
    T const tol = tolerance.at(x);
    bool const upward = ((dx > T(0)) == (f > T(0))) == (difference > T(0));
    T const move = math::abs(*step) > tol ? *step : (upward ? tol : T(-tol));
    if (!detail::sum_fits(x, move)) return std::nullopt;
    return move;
}
```

The secant step from `x` through `(x, f)` and `(xp, fp)`:

```
step = (x - xp) * f / (fp - f)
```

- **`fp - f` can't overflow**: before the hand-over the two values have one sign, so their difference is at most the
  larger magnitude. The caller forms it.
- **The ratio first.** `f / (fp - f)` and then the product with `dx`, both checked. The comment's example: with
  `1e300 (x - 1)` from `1e7`, forming `dx * f` first would overflow.
- **The minimum step**, `tol(x)` long, with its direction from the signs of `dx`, `f` and the difference: the step's
  sign is the product of the three signs, and `upward` computes it without the step, which may have underflowed.
- **`x + move` must fit**, or nothing.

### `open_secant` (1812-1851)

```cpp
template<real T, typename F>
constexpr std::expected<Root<T>, error> open_secant(F &fn, T const &guess, Tolerance<T> const &tolerance)
{
    if (!math::isfinite(guess)) return std::unexpected(error::non_finite);
    T xp = guess;
    std::array<T, 1> vp = detail::evaluate<1>(fn, xp);
    if (!math::isfinite(vp[0])) return std::unexpected(error::non_finite);
    if (vp[0] == T(0)) return Root<T>{ xp, tolerance.at(xp) };
    T x = detail::secant_second_point(xp, tolerance);
    T dx = T(x - xp); // small: x0 eta, or eta
    std::array<T, 1> vx = detail::evaluate<1>(fn, x);
    for (int i = 0;; ++i) {
        // The checks come first, so that the last point evaluated before the cap is checked too.
        if (!math::isfinite(vx[0])) return std::unexpected(error::non_finite);
        if (vx[0] == T(0)) return Root<T>{ x, tolerance.at(x) };
        if (!detail::same_sign(vx[0], vp[0])) return detail::hand_over<Brent>(fn, xp, vp, x, vx, tolerance);
        if (i == detail::iteration_cap<T>()) return std::unexpected(error::not_converged);
        std::optional<T> const step = detail::secant_move(x, vx[0], dx, T(vp[0] - vx[0]), tolerance);
        if (!step) return std::unexpected(error::not_converged);
        T const move = *step;
        xp = x;
        vp = vx;
        x = T(x + move);
        dx = move;
        vx = detail::evaluate<1>(fn, x);
    }
}
```

**The checks come first in the loop.** Every point evaluated is checked before the next step is computed, including
the last one before the cap: `i == cap` is tested after the checks. Otherwise a sign change at the last allowed point
would be reported as `not_converged`.

**The bookkeeping.** `xp, vp` is the previous point, `x, vx` the latest, `dx = x - xp` the last move. Each step shifts
them along, and `dx` becomes the move itself, rather than `x - xp` recomputed after rounding.

**Worked example: Secant from 1 on `x^2 - 2`:** the second point is `1.0001220703125` (`1 + 2^-13`); the secant through
the two lands at `1.49997`, where `fn = +0.25`: a sign change after 3 evaluations. Brent finishes on `[1.000122,
1.49997]` in 6 more: 9 in all.

### `open_steffensen` (1853-1895)

```cpp
template<real T, typename F>
constexpr std::expected<Root<T>, error> open_steffensen(F &fn, T const &guess, Tolerance<T> const &tolerance)
{
    if (!math::isfinite(guess)) return std::unexpected(error::non_finite);
    T x = guess;
    std::array<T, 1> vx = detail::evaluate<1>(fn, x);
    if (!math::isfinite(vx[0])) return std::unexpected(error::non_finite);
    if (vx[0] == T(0)) return Root<T>{ x, tolerance.at(x) };
    for (int i = 0; i < detail::iteration_cap<T>(); ++i) {
        T const tol = tolerance.at(x);
        T const offset = math::abs(vx[0]) > tol ? vx[0] : (vx[0] > T(0) ? tol : T(-tol));
        if (!detail::sum_fits(x, offset)) return std::unexpected(error::not_converged);
        T const z = T(x + offset);
        std::array<T, 1> const vz = detail::evaluate<1>(fn, z);
        if (!math::isfinite(vz[0])) return std::unexpected(error::non_finite);
        if (vz[0] == T(0)) return Root<T>{ z, tolerance.at(z) };
        if (!detail::same_sign(vz[0], vx[0])) return detail::hand_over<Brent>(fn, x, vx, z, vz, tolerance);
        T const width = math::abs(offset) <= math::abs(x) ? T(z - x) : offset;
        std::optional<T> const move = detail::secant_move(x, vx[0], T(-width), T(vz[0] - vx[0]), tolerance);
        if (!move) return std::unexpected(error::not_converged);
        T const next = T(x + *move);
        std::array<T, 1> const vn = detail::evaluate<1>(fn, next);
        if (!math::isfinite(vn[0])) return std::unexpected(error::non_finite);
        if (vn[0] == T(0)) return Root<T>{ next, tolerance.at(next) };
        if (!detail::same_sign(vn[0], vx[0])) return detail::hand_over<Brent>(fn, x, vx, next, vn, tolerance);
        x = next;
        vx = vn;
    }
    return std::unexpected(error::not_converged);
}
```

**What it computes.** Newton's step with the slope estimated from `fn` alone: `z = x + f(x)`, and the slope
`(f(z) - f(x)) / (z - x)`. Since `z - x = f(x)`, the step is `-f(x)^2 / (f(z) - f(x))`. It converges quadratically
near a simple root without a derivative, at 2 evaluations per iteration.

**The code** reuses the secant step: it is the secant through `(x, f(x))` and `(z, f(z))`, so `secant_move` with
`dx = x - z = -width` and the difference `f(z) - f(x)`.

- **`offset`** is `f(x)`, but at least `tol(x)` in `f(x)`'s direction, so that `z` differs from `x`.
- **`width`** is `z - x` as formed where `|offset| <= |x|`, where the comment's argument (Sterbenz) makes it exact:
  the true distance between the points actually evaluated, as Secant uses. Otherwise the intended `offset` is used,
  because there `z - x` could overflow.
- **Two hand-overs:** between `x` and `z`, and between `x` and the next iterate.

**Scale dependence.** `f(x)` is added to `x`, a length, so the method depends on `fn`'s units. The doc comment's
example shows it:

```
Steffensen on x^2 - 2 from 1:          z = 1 + (-1) = 0, f(0) = -2: one sign
                                       secant through (1, -1), (0, -2): next = 2, f = +2: sign change
                                       -> hand over on [1, 2]: Brent's trace from its eval 3 on; 10 in all

Steffensen on 1000 (x^2 - 2) from 1:   z = 1 + (-1000) = -999, f = +9.98e8: sign change at once
                                       -> hand over on [-999, 1], which brackets -sqrt 2, not +sqrt 2
                                       Brent: 25 more evaluations; value -1.4142135623730951
```

On `[-999, 1]` Brent spends 18 evaluations alternating: a secant step lands near 1, where `|fn|` hasn't decreased
(`|fa| <= |fb|`), so the next step is a bisection: -499, -249, -124, ..., -3.1. Only when the bracket is a few units
wide does the interpolation take over (evals 21 to 27).

---

## 14. `search`: finding a bracket (lines 2148-2294)

`search(fn, x0)` looks outward from a guess for a sign change, and returns a bracket for a bracketing method to
finish.

### `SearchSide` (2151-2173)

```cpp
template<real T>
struct SearchSide
{
    int direction;

    bool alive;

    T last;

    std::optional<T> edge;

    int probes;
};
```

The search runs on two sides, below and above `x0`, alternately. Each side has a direction, a flag, the last point
where `fn` was finite (and so had `fn(x0)`'s sign), and, once it has met a point where `fn` wasn't finite, that point
as `edge`, with a count of the probes towards it.

### `outward(fn, x0)` (2175-2249)

```cpp
template<real T, typename F>
constexpr std::expected<Bracket<T>, error> outward(F &fn, T const &x0)
{
    if (!math::isfinite(x0)) return std::unexpected(error::non_finite);
    T const f0 = detail::evaluate<1>(fn, x0)[0];
    if (!math::isfinite(f0)) return std::unexpected(error::non_finite);
    if (f0 == T(0)) return Bracket<T>{ x0, x0 };
    // |x0| / 50 rather than |x0| (1/50), whose rounding differs; at least denorm_min, since |x0| / 50
    // underflows to 0 for |x0| below 25 denorm_min.
    T h = x0 == T(0) ? T(T(1) / T(50)) : T(math::abs(x0) / T(50));
    if (h == T(0)) h = T(std::numeric_limits<T>::denorm_min());
    // A type whose denorm_min() is 0, where the standard says min() (independent reviews of Phase 5).
    if (h == T(0)) h = math::abs(x0);
    T const growth = T(T(3) / T(2));
    T const half_max = T(detail::largest<T>() / T(2));
    // Probes towards a domain edge: the bits of T, enough to approach an edge at 0 by as many halvings.
    int const probe_limit = detail::iteration_cap<T>() / 4;
    std::array<SearchSide<T>, 2> sides{ SearchSide<T>{ -1, true, x0, std::nullopt, 0 },
        SearchSide<T>{ 1, true, x0, std::nullopt, 0 } };
    for (int step = 0; step < detail::iteration_cap<T>(); ++step) {
        for (SearchSide<T> &side : sides) {
            if (!side.alive) continue;
            T x = x0;
            if (side.edge) {
                if (side.probes >= probe_limit) {
                    side.alive = false;
                    continue;
                }
                x = detail::midpoint(side.last, *side.edge);
                bool const between =
                    side.last < *side.edge ? side.last < x && x < *side.edge : *side.edge < x && x < side.last;
                if (!between) {
                    side.alive = false;
                    continue;
                }
                ++side.probes;
            } else {
                T const s = side.direction < 0 ? T(-h) : h;
                if (!detail::sum_fits(x0, s)) {
                    side.alive = false;
                    continue;
                }
                x = T(x0 + s);
            }
            T const fx = detail::evaluate<1>(fn, x)[0];
            if (!math::isfinite(fx)) {
                side.edge = x;
                continue;
            }
            if (fx == T(0)) return Bracket<T>{ x, x };
            if (!detail::same_sign(fx, f0)) {
                if (side.last < x) return Bracket<T>{ side.last, x };
                return Bracket<T>{ x, side.last };
            }
            side.last = x;
        }
        if (!sides[0].alive && !sides[1].alive) return std::unexpected(error::no_sign_change);
        // h grows by 3/2 where that fits: h (3/2) below max / 2, and above it h + h / 2, checked (the same
        // value in binary, where h / 2 is exact). A guard h <= max / (3/2) let 1.5 h overflow where the quotient
        // rounds up, as in double and std::float16_t (independent reviews of Phase 5). Else no larger step.
        if (h <= half_max)
            h = T(h * growth);
        else if (detail::sum_fits(h, T(h / T(2))))
            h = T(h + T(h / T(2)));
        else
            return std::unexpected(error::no_sign_change);
    }
    return std::unexpected(error::no_sign_change);
}
```

1. **`x0` and `fn(x0)` must be finite.** An exact 0 at `x0` is returned as the bracket `{x0, x0}`, which `solve`
   returns at once as a root (`prepare` sees `fn` exactly 0 at an end).
2. **The first step** is `h = |x0| / 50`, or 1/50 at `x0 = 0`: the search starts at `x0`'s own scale, so a
   temperature of 300 is searched in steps of 6, a pressure of 1e5 in steps of 2000. Written `|x0| / 50` rather than
   `|x0| * (1/50)`, whose rounding differs. If it underflows to 0, `denorm_min`; and if *that* is 0 (a type whose
   `denorm_min()` is 0; the test type `Flushed`), `|x0|`.
3. **Each step** evaluates `x0 - h` and then `x0 + h`, for each side still alive, with `x0 + s` checked by `sum_fits`.
   A side whose next point can't be formed ends.
    - **Not finite:** the point becomes the side's `edge`, and the side switches to probing.
    - **Exactly 0:** `{x, x}`.
    - **The other sign:** the bracket is `{last, x}` in order: the narrowest bracket that shows the sign change.
    - **The same sign:** `last = x`, and the search goes on.
4. **Probing towards an edge.** A side that has met a non-finite value bisects between its last finite point and the
   edge, at most as many probes as `T` has bits. A root can lie just inside a domain's edge (`log x + 3` has its root
   at 0.05, near the edge at 0), and a side that simply stopped at the first NaN would miss it.
5. **Growing the step.** `h` grows by 3/2 per step. Above `max / 2` the product `1.5 h` could overflow, so there it is
   `h + h / 2`, checked; if that doesn't fit, no larger step exists and the search ends. The comment records the
   review finding: an earlier guard `h <= max / 1.5` let `1.5 h` overflow where the quotient `max / 1.5` rounds up, as
   in double and `std::float16_t`.
6. **The cap** is `iteration_cap<T>()` steps, so at most `2 cap + 1` evaluations (425 in double); then, or when both
   sides have ended, `no_sign_change`.

### `search` (2254-2294)

```cpp
template<real T, typename F>
requires detail::gives<F, T, 1>
[[nodiscard]] constexpr std::expected<Bracket<T>, error> search(F &&fn, T const &x0)
{ return detail::outward(fn, x0); }
```

The public function only adds the constraint: `gives<F, T, 1>`, so a tuple `fn` works too (its first element). `T` is
deduced from `x0`; the doc comment warns that `search(fn, 0)` deduces `int`, which isn't a `real`.

**What it can't find**, as the doc comment lists: a root of even multiplicity (`fn` doesn't change sign), two roots
within one step (`x^2 - 2` from 100: both lie between two points of one side, and `fn` is positive at both), a root
beyond a region where `fn` isn't finite, and anything beyond the cap. A pole changes sign like a root, so `tan x` from
1.4 brackets the pole at `pi / 2`.

### Worked example 1: `exp(20 - 5000 / t) - 1` from 300

The root is at `t = 250` (`20 = 5000 / t`).

```
eval  x            f(x)       side, h
1     300          +27.0      x0
2     294          +19.0      below, h = 6
3     306          +37.9      above
4     291          +15.7      below, h = 9
5     309          +44.6      above
6     286.5        +11.8      below, h = 13.5
7     313.5        +56.5      above
8     279.75       +7.4       below, h = 20.25
9     320.25       +79.4      above
10    269.625      +3.3       below, h = 30.375
11    330.375      +128.8     above
12    254.4375     +0.42      below, h = 45.5625
13    345.5625     +251.4     above
14    231.65625    -0.79      below, h = 68.34375: the other sign
-> [231.65625, 254.4375], then Brent: 250 after 8 more evaluations
```

### Worked example 2: a domain edge, `log x + 3` from 1

The steps grow from 0.02 until the lower side reaches -0.153, where `log` is NaN (eval 22). That side then probes:
the midpoint of its last finite point 0.231 and the edge -0.153 is 0.0389, where `log x + 3 = -0.25`, the other sign.
The bracket `[0.0389, 0.2311]` holds the root `e^-3 = 0.0498`; 24 evaluations.

---

## 15. `std::formatter<nxx::roots::error>` (lines 2298-2340)

```cpp
template<>
struct formatter<nxx::roots::error> : formatter<string_view>
{
    template<typename FormatContext>
    auto format(nxx::roots::error e, FormatContext &ctx) const
    {
        // No default: a new enumerator makes the compilers warn here (-Wswitch) until its name is added.
        string_view name = "unknown";
        switch (e) {
        case nxx::roots::error::non_finite:
            name = "non_finite";
            break;
        case nxx::roots::error::no_sign_change:
            name = "no_sign_change";
            break;
        case nxx::roots::error::not_converged:
            name = "not_converged";
            break;
        }
        return formatter<string_view>::format(name, ctx);
    }
};
```

The same as `nxx::deriv`'s (`DERIV-WALKTHROUGH.md` section 12), with three names: `std::println("no root: {}",
r.error())` prints `no root: no_sign_change`. The `switch` has no `default`, so that a new enumerator makes the
compilers warn here until its name is added.

---

## 16. Habits you will see everywhere, and why

`DERIV-WALKTHROUGH.md` section 13 lists the habits the two headers share (`T(1)`, never `auto` on a number,
`std::invoke`, qualified calls, check then compute, `!(h > 0)`). These are particular to `roots.hpp`:

| Habit                                                                     | Why                                                                                                                      |
|---------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------|
| every distance through `distance_bound`                                   | an upper bound, never below the true distance, formed without overflow; exact at every default stop (section 5)          |
| every quotient, product and sum through `quotient`, `product`, `sum_fits` | `max` is a legitimate value of `fn` and of the ends; a failed check selects a fallback point, never an error (section 5) |
| `std::optional<T>` for a proposed point                                   | "no point" is a normal outcome; the caller falls back on the secant, the midpoint or a bisection                         |
| bisection whenever something is wider than `max / 16`                     | no interpolation is formed where its products could overflow (sections 8, 9)                                             |
| a step of at least `tol` towards the other end                            | near the root the steps become tiny; the minimum step closes the bracket (sections 8, 13)                                |
| directions from signs, not from a step                                    | a tiny step can underflow to 0, which has no direction (section 13)                                                      |
| a point must be strictly inside the bracket                               | otherwise `T`'s resolution is reached, and evaluating again would loop (sections 7, 8, 9)                                |
| `T(T(1) / T(2))`, `T(T(7) / T(10))`                                       | fractions are formed from integers in `T` itself, never from a `double` literal                                          |
| derivatives never checked, only the value                                 | a bad derivative only makes its step unusable; the method takes another (sections 3, 6, 8)                               |
| pointers to `fn`, or `fn` as a parameter, never a lambda capturing it     | MSVC rejects a captured function reference in constant evaluation (`Joined`, `advance`)                                  |

---

## 17. Check yourself

Each answer is hidden under its question; CLion's Markdown preview shows them as expandable sections.

<details><summary>1. Why does <code>solve&lt;Bisection&gt;</code> on <code>[1, 2]</code> evaluate <code>fn</code> 52 times, and why isn't the value it returns one of the points it evaluated?</summary>

2 ends, then 50 midpoints: after 50 halvings the half-width 2^-51 = 4.4e-16 is below tol (m) = 6.3e-16. The returned
point is the next midpoint, which is within half the width of the sign change whatever `fn` is there, so evaluating it
would add nothing to the guarantee.
</details>

<details><summary>2. Brent on <code>sin</code> over <code>[-1, 2]</code> returns 0 after 3 evaluations with <code>abs_error</code> 4.9e-324. Where do the 3 evaluations and the bound come from?</summary>

`prepare` evaluates the two ends, then, because 0 is strictly inside, `fn(0)`, which is exactly 0. An exact zero is
returned with `tol(0) = abs_tol = denorm_min = 4.9e-324` as its bound.
</details>

<details><summary>3. In Brent's trace on <code>x^2 - 2</code>, eval 8 is already the closest double to the root. Why is there a 9th evaluation?</summary>

The bracket `[b, c]` is still 4e-14 wide, because `c` is eval 7's point. The stop test needs it at most 2 tol (b) =
1.26e-15. The interpolation proposes a step shorter than tol, so the minimum step, tol towards `c`, is taken; it lands
on the other side of the root, and the bracket is then 6.7e-16 wide.
</details>

<details><summary>4. Why does <code>distance_bound</code> exist, instead of <code>abs(b - a)</code>?</summary>

`b - a` overflows for ends of opposite signs far apart, which is a compile error in constant evaluation; and where it
doesn't overflow it may round below the true distance, which would make `abs_error` a false bound. `distance_bound` is
never below the true distance, never overflows, and is exact in the common case (Sterbenz).
</details>

<details><summary>5. A false-position step's proposed point is 1e-17 above the lower end of a bracket around 1.4. Where is <code>fn</code> evaluated?</summary>

At `a + 2 delta`, with `delta = 0.7 tol(u)`, about 8.8e-16 inside (`place`). Either the root lies in that sliver and
the next stop test ends the search, or the step removes at least 1.4 tol of the bracket.
</details>

<details><summary>6. Why does <code>RegulaFalsi</code> on <code>x^2 - 2</code> take 20 evaluations where <code>AndersonBjorck</code> takes 9?</summary>

`x^2 - 2` is convex, so every secant point lands below the root and the upper end is retained again and again. Regula
falsi then converges linearly from one side, and only the guard bisections (after three steps without halving) move
the upper end. Anderson-Bjorck scales the retained end's value down each time it is retained again, which pulls the
secant past the root and moves the other end.
</details>

<details><summary>7. Why is Ridders' point computed as <code>e + w (m - e)</code> with <code>w = -g_e g_d / (root (root + |g_m|))</code>?</summary>

The textbook form needs `1 - |f(m)| / root`, which rounds to 0 when `f(e) f(d)` is below eps of `f(m)^2`, and then
puts the point exactly on the end `e`. Multiplying numerator and denominator by `root + |g_m|` turns the difference
into `-g_e g_d`, a product of values of opposite signs, with no cancellation. The scaling by the largest magnitude keeps
the squares from overflowing.
</details>

<details><summary>8. Newton from a guess returns a bracket-verified root. Where does the bracket come from?</summary>

The iterates: as soon as two consecutive iterates have values of opposite signs, they are a bracket, and `hand_over`
continues with `brent<Newton>` on it. A step of at least tol makes sure the iterates cross the root once they are
within tol of it.
</details>

<details><summary>9. Why does Steffensen on <code>1000 (x^2 - 2)</code> from 1 return <code>-sqrt 2</code>?</summary>

Its first point is `z = x + f(x) = 1 - 1000 = -999`, where `fn` is positive. That is a sign change between 1 and -999,
so it hands over at once to Brent on `[-999, 1]`, which contains only the negative root.
</details>

<details><summary>10. <code>search(log x + 3, 1)</code> meets a NaN at -0.153. Why doesn't that side just stop?</summary>

The root, at 0.0498, lies between the side's last finite point 0.231 and the edge. A side that stopped would miss every
root near a domain's edge. It bisects towards the edge instead, and its first probe, 0.0389, already has the other
sign.
</details>

<details><summary>11. Why does <code>detail::sqrt</code> use <code>y + (s / y - y) / 2</code> rather than <code>(y + s / y) / 2</code>?</summary>

In radix 10, `y + s / y` is about `2 y` and can fall in the decade above `y`, where it rounds to one digit fewer; the
textbook form is then up to 3 ulps off. Near the root `s / y - y` is small and exact, so the correction form adds a
small, accurately computed correction to `y` and stays within one ulp.
</details>

<details><summary>12. Why is a failed derivative from <code>derivative_of</code> stored as 0, and not as NaN?</summary>

0 has the same effect as a non-finite derivative in both steps (Newton's step can't be formed with f' = 0; Halley's
step becomes Newton's with f'' = 0), and every `real` has a 0, while a NaN or an infinity need not exist in a user's
type.
</details>

---

## Glossary

| Term                 | Meaning here                                                                                                 |
|----------------------|--------------------------------------------------------------------------------------------------------------|
| **bracket**          | an interval at whose ends `fn` has values of opposite signs, so it contains a sign change                    |
| **sign change**      | a point where `fn`, as evaluated, changes sign: a root, a jump, or a pole                                    |
| **tol(x)**           | `abs_tol + rel_tol * abs(x)`, saturated at `max`; a method stops when the bracket is at most `2 tol(u)` wide |
| **u**                | the end of the bracket with the smaller absolute value of `fn`; the point returned                           |
| **stop test**        | the bracket at most `2 tol(u)` wide, or `fn` exactly 0                                                       |
| **zero first**       | evaluating `fn(0)` when the bracket contains 0, and splitting the bracket there                              |
| **retained end**     | in false position, the end that stays when the new point lands on the other end's side                       |
| **minimum step**     | a step of `tol` towards the other end, taken where the method's own step is shorter                          |
| **hand-over**        | an open method continuing as a bracketing method once two consecutive iterates have opposite signs           |
| **delta, margin**    | `0.7 tol(u)`, and `2 delta`: how far inside the bracket `place` keeps a point                                |
| **wide**             | a bracket (or step) above `max / 16`, where only bisection steps are taken                                   |
| **Sterbenz's lemma** | if `y / 2 <= x <= 2 y`, then `x - y` is exact in floating point                                              |
| **faithful**         | within one unit in the last place of the true value                                                          |
| **saturate**         | return `max` instead of overflowing; used for tolerances and distances                                       |
| **edge**             | in `search`, the nearest point where `fn` was not finite, towards which a side probes                        |

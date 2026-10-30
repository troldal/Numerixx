# Working with Numerixx

Each module of Numerixx solves its own kind of problem, but they all work the same way. You pass a function, written as
a lambda or any other callable and called `fn` throughout this documentation, together with the numbers that start the
computation. The type of those numbers sets the type of the whole computation. You get back a result or the reason there
is none, never an exception. And you can run the same call at compile time. This page explains these conventions once,
so the module chapters need not repeat them. Most of its examples use `nxx::roots`, shortened to `rt` by
`namespace rt = nxx::roots;`; the same rules hold for `nxx::deriv`, shortened to `dv`, whose examples come at the end of
[Results](#results).

## Real types

Numerixx works in any *real type*: a type that behaves like a floating-point number. The built-in ones are `float`,
`double` and `long double`; others come from libraries such as Boost.Multiprecision, whose numbers can carry 50 digits
or more. You choose the type by what you pass. A computation runs in the type of its starting numbers, the ends of a
bracket or the point of a derivative, and `fn` takes a number of that type and returns one. Write the starting numbers
in the type you want: `Bracket{1.0f, 2.0f}` is a `float` bracket, while `Bracket{1.0, 2.0}` is a `double` one. (A
*bracket* is an interval at whose ends `fn` has opposite signs, in which a root is sought.)

`long double` differs between platforms. With MSVC and clang-cl it is the same as `double`; with GCC and Clang on x86-64
it is the x87's extended precision, with a 64-bit significand; with Emscripten it is quadruple precision, with 113 bits.
Take tolerances and expected accuracies from `std::numeric_limits<long double>`, not from a fixed number of digits.

### Multiprecision types

A Boost.Multiprecision number works like a `double`, with three differences you have to handle in `fn`. Here is a solve
of $\cos x = x$ in `cpp_dec_float_50`, a type with 50 decimal digits. The program includes the type's header beside the
library's:

```{literalinclude} examples/basics_types.cpp
:start-after: // <includes>
:end-before: // </includes>
```

and solves

```{literalinclude} examples/basics_types.cpp
:start-after: // <multiprecision>
:end-before: // </multiprecision>
:dedent:
```

A Boost number does not print through `std::println` as a `double` does: its member `str(n)` gives its first `n` digits
as a string. The program prints the root to all 50 digits:

```{literalinclude} examples/basics_types.cpp
:language: text
:start-after: output = R"(
:end-before: )";
```

The three differences are these.

- **Return the type, not the expression.** Some Boost types, `cpp_dec_float_50` among them, use *expression templates*:
  `cos(x) - x` is not a number but an object that records the operation and evaluates it later, and it may refer to
  `fn`'s locals or to a parameter taken by value, which are gone before Numerixx reads it. Converting it to the type,
  `dec50(cos(x) - x)`, makes `fn` return a number. This is a rule, checked when the call compiles: for a type that is a
  class, `fn` must return the type itself (or a reference to it), or a built-in number such as a `double`, and a call
  whose `fn` returns an expression does not compile. Types without expression templates, such as `cpp_bin_float_50`,
  return the type from `cos(x) - x` already, so the conversion there, `bin50(cos(x) - x)`, is a harmless habit. For
  the built-in types the rule is looser: a `fn` for `double` may return anything that converts to `double`, such as a
  `float`.
- **Call mathematical functions unqualified.** `std::cos` is only for the built-in types, and `std::cos(x)` does not
  compile for a Boost number. Write `using std::cos;` and then `cos(x)`, with a using-declaration for each function you
  call (`using std::sin;` for `sin`): the call finds Boost's function for a Boost type and the standard one for a
  built-in type, so the same `fn` works for both.
- **Write constants that are exact in the type.** An integer such as `2` converts exactly. A decimal literal such as
  `0.1` is a `double`, already rounded to 16 digits before it reaches your type; write `dec50("0.1")` or `dec50(1) / 10`
  instead.

### A type of your own

A type of your own works if it satisfies the concept {cpp:concept}`nxx::real`: it specialises `std::numeric_limits` as a
bounded, inexact, non-integer type, is constructible from `int`, is totally ordered, and has `+ - * /`. The [core
reference](core.md) lists every requirement. When a type does not satisfy it, the call does not compile, and the
compiler's message names the constraint that failed. Look for `nxx::real` in it. The most common case is a start written
as integers: `Bracket{1, 2}` asks for a bracket of `int`s, and `derivative(fn, 2)` for a derivative at an `int` point,
and `int` is not a real type.

## Results

Every algorithm returns a [`std::expected`](https://en.cppreference.com/w/cpp/utility/expected): either the result, or
an error saying why there is none. Test it before you read it. Here the first bracket holds the root of $x^2 - 2$ and
the second does not:

The program includes both modules' headers:

```{literalinclude} examples/basics_results.cpp
:start-after: // <includes>
:end-before: // </includes>
```

```{literalinclude} examples/basics_results.cpp
:start-after: // <test>
:end-before: // </test>
:dedent:
```

```{literalinclude} examples/basics_results.cpp
:language: text
:start-after: test_output = R"(
:end-before: )";
```

After the test, read the result with `->` or `*`. An error formats as its name, so you can print or log it directly, as
the second line shows. Each module has its own error enumeration, listed in its reference.

Numerixx itself never throws, and never catches. Exceptions can still reach you from outside the library. An exception
thrown by `fn`, or by the operations of a real type of your own, passes through unchanged. And `value()`, the standard
library's checked access, throws `std::bad_expected_access` when there is no result, which is why you test first.

```{note}
More errors may be added in a later version, and a cause may move from one error to another. Test for success rather
than for a particular error, and give a `switch` on the error a `default`.
```

### Composing steps

A computation often has several steps, any of which can fail. `std::expected` chains them without a test after each
step:

- `and_then` runs a next step that can fail, only if the previous one succeeded;
- `transform` runs a next step that cannot fail, such as reading a member;
- `value_or` gives a plain value at the end, with a fallback for any failure;
- `or_else` runs a fallback step, only if the previous one failed.

This finds a bracket from a guess, solves in it and keeps the root's value, with NaN if any step fails:

```{literalinclude} examples/basics_results.cpp
:start-after: // <compose>
:end-before: // </compose>
:dedent:
```

and this tries Newton's method from a guess, then Brent's method in a bracket if Newton's method fails, as it does here
because the derivative is 0 at the guess:

```{literalinclude} examples/basics_results.cpp
:start-after: // <or-else>
:end-before: // </or-else>
:dedent:
```

```{literalinclude} examples/basics_results.cpp
:language: text
:start-after: compose_output = R"(
:end-before: )";
```

### What a result holds

Each module's result has a `value` and an `abs_error`, but `abs_error` means something different in each, so read the
module's description of it.

- **A root** ({cpp:struct}`nxx::roots::Root`): `abs_error` is a *bound*, always present, on the distance from `value` to
  a point where `fn`, as computed, changes sign. [Roots of functions](roots/index.md#what-you-get-back) explains what it
  covers.
- **A derivative** ({cpp:struct}`nxx::deriv::Estimate`): `abs_error` is an *estimate* of the error, not a guarantee, and
  a `std::optional`: Ridders' methods fill it, the plain finite differences leave it empty, so test it before reading it
  with `*`. An `Estimate` converts to its value, so it can be used where a number is expected. That has a trap: for a
  result `d`, `if (*d)` compiles, and tests whether the derivative is 0. Test `d` itself.

Here are both kinds of derivative of $1/x$ at 2, which is $-1/4$:

```{literalinclude} examples/basics_results.cpp
:start-after: // <estimate>
:end-before: // </estimate>
:dedent:
```

```{literalinclude} examples/basics_results.cpp
:language: text
:start-after: estimate_output = R"(
:end-before: )";
```

Read the members by name, not with a structured binding: a later version may add members.

## Compile time

Every algorithm is `constexpr`, so with a `constexpr` function it runs in a constant expression, and its result can be
checked by `static_assert` or used as a constant:

```{literalinclude} examples/basics_types.cpp
:start-after: // <constexpr>
:end-before: // </constexpr>
:dedent:
```

Two things limit this. First, `fn` must be `constexpr` itself. Arithmetic and comparisons are, but the mathematical
functions of `<cmath>` mostly are not: C++23 makes a few `constexpr`, such as `std::abs` and `std::floor`, but not
`std::sin`, `std::exp` or `std::sqrt` (GCC accepts some of them as an extension, which other compilers reject). Second,
the type must be a *literal type*, usable in constant expressions: the built-in types are, while the
Boost.Multiprecision types on this site, `cpp_bin_float_50` and `cpp_dec_float_50`, are not.

A constant gives exactly what the same call gives at run time, as long as the compiler does not fuse a multiplication
and an addition into one instruction (a fused multiply-add) at run time. Compilers for x86-64 and WebAssembly do not by
default, but GCC and Clang do as soon as you target hardware that has the instruction (`-march=native`, `-mfma`);
`-ffp-contract=off` turns that off. Fused arithmetic also changes the last digits of results, and the counts of
evaluations, that this site quotes.

```{note}
A compiler limits the work in one constant expression. If it reports that a solve exceeds the limit, raise it:
`-fconstexpr-steps` for Clang, `-fconstexpr-ops-limit` for GCC, `/constexpr:steps` for MSVC.
```

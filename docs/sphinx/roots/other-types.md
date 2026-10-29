# Solve in other types

Every method works in any type that satisfies {cpp:concept}`nxx::real`: `float`, `long double`, Boost.Multiprecision
numbers, or a type of your own. The type comes from the start, the bracket's ends or the guess, so write them in the
type you want.

## `float` and `long double`

```{literalinclude} ../examples/roots_other_types.cpp
:start-after: // <float>
:end-before: // </float>
:dedent:
```

Write the literals in the type as well (`2.0f`, `2.0L`): `Bracket{1.0, 2.0}` would be a `double` bracket. The result is
accurate to the precision of the type. That precision differs by platform for `long double`: it is the same as `double`
with MSVC and clang-cl, the x87's 80-bit extended precision (a 64-bit significand) with GCC and Clang on x86-64, and
quadruple precision (a 113-bit significand) with Emscripten. Take tolerances and expectations from
`std::numeric_limits<long double>`, not from a fixed number of digits.

## Boost.Multiprecision

```{literalinclude} ../examples/roots_other_types.cpp
:start-after: // <multiprecision>
:end-before: // </multiprecision>
:dedent:
```

The root comes out to the 50 digits of the type. Three rules make a multiprecision fn correct:

- Return the type itself, `bin50(x * x - 2)`, not the arithmetic expression. Some Boost types, such as
  `cpp_dec_float_50`, use expression templates (`cpp_bin_float_50` does not): `x * x - 2` is then an object that may
  refer to fn's locals or to a parameter taken by value, which are gone once fn returns. For a pair or tuple of
  derivatives, name the element types, `std::pair<T, T>`:

  ```{literalinclude} ../examples/roots_other_types.cpp
  :start-after: // <expression-templates>
  :end-before: // </expression-templates>
  :dedent:
  ```

- Call mathematical functions unqualified, after `using std::sin;`, so that the type's own `sin` is found.
- Set `abs_tol` where a sign change may lie near 0: with the default, the method approaches it through the type's
  whole exponent range ([Set tolerances](tolerances.md)). The open methods' iteration cap also grows with the
  type's precision.

## At compile time

Every method is `constexpr`. In a constant expression, fn must be `constexpr` too: `+ - * /` and comparisons, but not,
portably, `std::sin`, `std::exp`, `std::log` or `std::sqrt`, which C++23 does not make `constexpr` (GCC accepts some as
an extension). The type must be a literal type, which excludes the Boost types:

```{literalinclude} ../examples/roots_other_types.cpp
:start-after: // <constexpr>
:end-before: // </constexpr>
:dedent:
```

The constant gives exactly what the same call gives at run time, on builds that do not contract `x * x - 2.0` into a
fused multiply-add (the default on x86-64 and WebAssembly).

The program needs:

```{literalinclude} ../examples/roots_other_types.cpp
:start-after: // <includes>
:end-before: // </includes>
```

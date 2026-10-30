# Solve in other types

A `double` holds about 16 significant digits, which suits most problems. Sometimes you want a different type: `float` to
save memory or time where 7 digits are enough, or a multiprecision type when 16 digits are not, say to check a `double`
computation against a more precise one. Every method works in any *real type*, a type that behaves like a floating-point
number. [Working with Numerixx](../basics.md#real-types) explains what a real type is and the rules for using one; this
page covers what is particular to solving.

## The type comes from the start

A solve runs in the type of its start: the ends of the *bracket* (an interval at whose ends `fn` has opposite signs), or
the guess. Write them in the type you want, and write `fn` to take that type:

```{literalinclude} ../examples/roots_other_types.cpp
:start-after: // <float>
:end-before: // </float>
:dedent:
```

Write the literals in the type as well, `2.0f` and `2.0L`. `Bracket{1.0, 2.0}` is a `double` bracket whatever `fn`
takes, and a `Bracket{1, 2}` of `int`s does not compile. A `float` function with a `double` bracket compiles (MSVC warns
of the conversion), but it rounds every point to `float`: its root is then only as accurate as a `float`, about
$3 \times 10^{-8}$ from $\sqrt 2$ here, while `abs_error` still reports a `double`'s bound. Each result is accurate to
the precision of its type: about 7 digits for `float`; for `long double`, whatever the platform gives it ([Real
types](../basics.md#real-types)).

## Boost.Multiprecision

A Boost.Multiprecision type works the same way. Here is $\sqrt 2$ in `cpp_bin_float_50`, a binary type with 50 decimal
digits:

```{literalinclude} ../examples/roots_other_types.cpp
:start-after: // <multiprecision>
:end-before: // </multiprecision>
:dedent:
```

```{literalinclude} ../examples/roots_other_types.cpp
:language: text
:start-after: output = R"(
:end-before: )";
```

The root comes out to all 50 digits. `fn` returns `bin50(x * x - 2)`, the type itself rather than the arithmetic
expression, as every `fn` in a Boost type should ([Multiprecision types](../basics.md#multiprecision-types) explains
why).

The same holds for the derivatives that Newton's and Halley's methods take. Name the element types of the pair or tuple,
`std::pair<T, T>`, so that each element is converted to the type before `fn` returns. `cpp_dec_float_50`, a decimal
type, uses expression templates, where it matters most:

```{literalinclude} ../examples/roots_other_types.cpp
:start-after: // <expression-templates>
:end-before: // </expression-templates>
:dedent:
```

Without the element types, `std::pair{x * x - 2, 2 * x}` would hold two expressions, and the call would not compile: the
elements must be exactly the type.

Two things change in a multiprecision type, because of its precision and its range:

- **Set `abs_tol` when a root may lie near zero.** With the default tolerances, a method approaches a sign change near
  zero down to the type's smallest numbers. A multiprecision type's exponent range is enormous, so that walk can take
  about $2 \times 10^9$ halvings in `cpp_bin_float_50`, and the solve seems never to end. [Set
  tolerances](tolerances.md#a-root-near-zero-set-abs_tol) shows how.
- **The open methods try longer.** From a guess, Newton's, Halley's, the secant and Steffensen's methods give up with
  `not_converged` after a number of steps that grows with the type's precision: 4 times its number of bits, 212 for a
  `double`. A solve that fails in a multiprecision type therefore takes longer to fail.

A solve in a multiprecision type also cannot run at compile time, since those types are not literal types ([Compile
time](../basics.md#compile-time)).

The program needs:

```{literalinclude} ../examples/roots_other_types.cpp
:start-after: // <includes>
:end-before: // </includes>
```

# High precision

A `double` carries about 16 significant digits. When you need more, an integral to 40 digits to check another
computation or a constant to be tabulated, compute in a multiprecision type such as Boost.Multiprecision's
`cpp_bin_float_50`, which carries 50. Every method of `nxx::quadrature` works in such a type, but three things change:
the default request, the way `fn` is written, and which method is cheap. The examples write `qd` for
`nxx::quadrature` and `bin50` for `cpp_bin_float_50`; [Working with Numerixx](../basics.md#multiprecision-types)
explains the type itself.

## Ask for the digits you want

The default criteria ask for a relative error of $\sqrt{\varepsilon}$, where $\varepsilon$ is the type's machine
epsilon: half its digits, about 25 in `cpp_bin_float_50`. The methods do not certify more than you ask for, so a
caller who chose 50 digits for their digits must say so. The program includes the type's header beside the
library's:

```{literalinclude} ../examples/quadrature_high_precision.cpp
:start-after: // <includes>
:end-before: // </includes>
```

Here the default method, adaptive Gauss-Kronrod, integrates $1/(1 + x^2)$ over $[0, 1]$, which is $\pi/4$, with the
default request and with a request of 40 digits, and `TanhSinh` does the same with 40:

```{literalinclude} ../examples/quadrature_high_precision.cpp
:start-after: // <request>
:end-before: // </request>
:dedent:
```

```{literalinclude} ../examples/quadrature_high_precision.cpp
:language: text
:start-after: request_output = R"(
:end-before: )";
```

With the default request the default method stops after 420 evaluations, with an estimate of $1.8 \times 10^{-26}$.
Its value is in fact accurate to more than 40 digits, but nothing in the result tells you that: the estimate is what you
can rely on. Asked for 40 digits, it meets the request, at 5,430 evaluations. It cuts the interval into *panels*,
pieces it integrates one at a time, and its estimate comes from the 7-point Gauss rule inside each panel, which gains
digits slowly as the panels shrink, so full precision in a 50-digit type is costly.
`TanhSinh` gains digits much faster as it halves its step, and meets the same request with 325 evaluations, about a
sixteenth. For many digits, choose `TanhSinh`.

## Write `fn` in the type

`fn` takes a `bin50` and returns one, and every number in it must carry the type's 50 digits:

```{literalinclude} ../examples/quadrature_high_precision.cpp
:start-after: // <pendulum>
:end-before: // </pendulum>
:dedent:
```

This is the pendulum's integral from the [tutorial](tutorial.md), at an amplitude of $90^\circ$. It follows the
habits that keep `fn` correct in any multiprecision type, with or without *expression templates* (objects that
record an operation and evaluate it later), which [Working with Numerixx](../basics.md#multiprecision-types)
explains: constants written in the type, so `bin50(1) / 2` rather than `0.5`, and `bin50("1e-40")` for a tolerance,
since a literal such as `1e-40` is a `double` rounded to 16 digits; mathematical functions called unqualified after
`using std::sqrt;`; a named type for every intermediate, `bin50 const w`, never `auto`; and a result converted to the
type, `return bin50(...)`. `cpp_bin_float_50` has no expression templates, so the last two cost nothing here, and the
same `fn` then works for a type that has them.

The program prints

```{literalinclude} ../examples/quadrature_high_precision.cpp
:language: text
:start-after: pendulum_output = R"(
:end-before: )";
```

To check it, compute $K = \pi / (2\,\mathrm{AGM}(1, k'))$ in the same type. $\pi$ itself must come in the type too:
`std::numbers::pi` is a `double`, so take $4 \arctan 1$ with Boost's `atan`:

```{literalinclude} ../examples/quadrature_high_precision.cpp
:start-after: // <reference>
:end-before: // </reference>
:dedent:
```

The two agree to within $10^{-48}$, so all 40 printed digits are right.

## The first call

For `float`, `double` and `long double` the methods' nodes and weights are tables computed when the program is
compiled. A type such as `cpp_bin_float_50` has no such table: the first call in that type computes the nodes it needs,
to the type's precision, from arithmetic alone, and keeps them for every later call. A program that integrates many
times in the type pays that cost once. The nodes are computed for the precision the type has at the first call, so a
type whose precision can be changed at run time, such as Boost's `mpfr_float`, keeps the nodes of its first precision.

```{note}
The cost of each evaluation of `fn` grows with the precision too: a `bin50` operation is many times slower than a
`double` one. Count the evaluations, as on [Choose a method](choose-a-method.md), before you plan a large computation
at high precision.
```

# Higher derivatives

The curvature of a function, the rate of change of a slope, a Taylor expansion to a few terms: each needs a derivative
above the first. `nxx::deriv` computes the second, third and fourth derivatives with every method. They cost more
evaluations of `fn`, and, more importantly, they are far less accurate than the first derivative. This page shows how to
ask for them, how accurate they are, and why `float` is not enough.

## Asking for the N-th derivative

The examples write `dv` for `nxx::deriv`. The order is the second template argument of `derivative`, after the method:
an integer from 1 to 4, or a name from the enumeration {cpp:enum}`nxx::deriv::order` (`order::second` is 2), whichever
reads better. Here are $f''$ and $f''''$ of $1/x$ at 2, which are $2/x^3 = 0.25$ and $24/x^5 = 0.75$:

```{literalinclude} ../examples/deriv_higher.cpp
:start-after: // <orders>
:end-before: // </orders>
:dedent:
```

```{literalinclude} ../examples/deriv_higher.cpp
:language: text
:start-after: orders_output = R"(
:end-before: )";
```

`Central<2>` for $f''$ is the classic second difference, $(f(x + h) - 2 f(x) + f(x - h)) / h^2$, at three points. The
default step changes with the order, so you need not choose a different one.

## How accuracy falls with the order

A derivative of order $N$ divides the rounding error of `fn`'s values by $h^N$ rather than by $h$, so each order loses
digits. For the central difference, the balanced step is $\varepsilon^{1/(2 + N)}$ and the error about
$\varepsilon^{2/(2 + N)}$: in `double`, about 10 digits for $f'$, 8 for $f''$, 5 or 6 for $f'''$ and 4 or 5 for $f''''$
(the constants the default step leaves out cost up to a digit at the higher orders). The relative errors of `Central<2>`
and Ridders' method on $1/x$ at 2 show the fall:

```{literalinclude} ../examples/deriv_higher.cpp
:language: text
:start-after: table_output = R"(
:end-before: )";
```

Ridders' method keeps far more digits at every order, with an error estimate, but it costs more: typically 13 to 23
evaluations for $f''$ and 25 to 45 for $f'''$ and $f''''$ ({cpp:struct}`nxx::deriv::Ridders`), against 3, 4 and 5 for
`Central<2>`. Above the first derivative, it is the method to use unless evaluations are expensive.

## Not in `float`

A `float` holds about 7 significant digits, and a second derivative loses about half of them. So a derivative above the
first in `float`, or in any type with 7 or fewer digits, gives a compiler warning, by every method (though not through
`derivative_of`, [Pass a derivative as a function](derivative-of.md#its-copy-of-fn)):

```{literalinclude} ../examples/deriv_higher.cpp
:start-after: // <float>
:end-before: // </float>
:dedent:
```

```{literalinclude} ../examples/deriv_higher.cpp
:language: text
:start-after: float_output = R"(
:end-before: )";
```

The result, $0.2500806$ for $0.25$, keeps only three or four digits. Use `double`. If you have checked that the accuracy
is enough for your purpose, silence the warning at the call, with your compiler's pragma for deprecation warnings
(`-Wdeprecated-declarations` for GCC and Clang, C4996 for MSVC), or everywhere, by defining
`NUMERIXX_NO_PRECISION_WARNING`. Define it on the compiler's command line, so that every file sees it alike, or before
the header, as this page's one-file program does:

```{literalinclude} ../examples/deriv_higher.cpp
:start-after: // <silence>
:end-before: // </silence>
```

The program needs:

```{literalinclude} ../examples/deriv_higher.cpp
:start-after: // <includes>
:end-before: // </includes>
```

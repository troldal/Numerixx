# Higher derivatives

The curvature of a function, the rate of change of a slope, a Taylor expansion to a few terms: each needs a derivative
above the first. `nxx::deriv` computes the second, third and fourth derivatives with every method. They are far less
accurate than the first derivative, and not usable in `float`. This section shows how to ask for them and what
accuracy to expect.

## Asking for the N-th derivative

The derivative's order $N$ is the second template argument of `derivative`, after the method: an integer from 1 to 4,
or a name from the enumeration {cpp:enum}`nxx::deriv::order`. Here are $f''$ and $f''''$ of $1/x$ at 2, which are
$2/x^3 = 0.25$ and $24/x^5 = 0.75$:

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

`Central<2>` for $f''$ is the classic second difference, $(f(x + h) - 2 f(x) + f(x - h)) / h^2$, at three points.

## How accuracy falls with the order

A derivative of order $N$ divides the rounding error of `fn`'s values by $h^N$ rather than by $h$. The default step
grows with $N$ to balance that, but each order still loses digits. In a `double`, the central difference keeps about 10
digits for $f'$, 8 for $f''$, 5 or 6 for $f'''$ and 4 or 5 for $f''''$. The relative errors of `Central<2>` and
Ridders' method on $1/x$ at 2 show the fall:

```{literalinclude} ../examples/deriv_higher.cpp
:language: text
:start-after: table_output = R"(
:end-before: )";
```

Ridders' method keeps far more digits at every order, with an error estimate, for more evaluations
([Costs](costs.md)). Above the first derivative, it is the method to use unless evaluations are expensive.

## Not in `float`

A `float` holds about 7 significant digits, and a second derivative loses about half of them. So a derivative above the
first in `float`, or in any type with 7 or fewer digits, gives a compiler warning, by every method:

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

The result, $0.2500806$ for $0.25$, keeps only three or four digits. Use `double`. If you have checked that the
accuracy is enough, define `NUMERIXX_NO_PRECISION_WARNING` on the compiler's command line, or before the header in
every file:

```{literalinclude} ../examples/deriv_higher.cpp
:start-after: // <silence>
:end-before: // </silence>
```

To silence one call instead, wrap it in your compiler's pragma for deprecation warnings. `derivative_of` gives no
warning, though its accuracy is just as poor.

The program needs:

```{literalinclude} ../examples/deriv_higher.cpp
:start-after: // <includes>
:end-before: // </includes>
```

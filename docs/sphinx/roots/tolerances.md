# Set tolerances

A method stops when the bracket around the root is narrow enough. A *bracket* is an interval at whose ends `fn` has
opposite signs, and "narrow enough" is set by two *tolerances* added together: a relative one, `rel_tol`, which asks for
a number of significant digits, and an absolute one, `abs_tol`, which asks for a distance. At the current estimate $x$
of the root they give

$$\text{tol}(x) = \text{abs\_tol} + \text{rel\_tol}\,|x|$$

and the method stops once the bracket is at most $2\,\text{tol}(x)$ wide, or when `fn` is exactly 0. Whichever of the
two terms is larger decides. Every result keeps its side of this: its `abs_error` is at most $2\,\text{tol}$ at the root
it returns.

The defaults ask for all but the last digit or so of the type, whatever the size of the root, so you rarely need to
change them. This page is about the cases where you do: a root near zero, and a root you need only to a few digits.

You pass tolerances as a {cpp:struct}`nxx::roots::Criteria`, usually as a braced list after the start:
`solve<Brent>(fn, bracket, {.abs_tol = 1e-12})`. The list converts to the type of the bracket, so the same `double`
values work for a multiprecision type.

## The defaults

The default `abs_tol` is the smallest positive number of the type, and the default `rel_tol` is twice its machine
epsilon (the gap between 1 and the next number of the type). The absolute term is then negligible for any root that is
not extremely close to zero, and the relative term asks for the full precision of the type. The defaults assume nothing
about the scale of your problem: a root at $10^{-12}$ keeps all its digits.

```{note}
The defaults are also the minimum. A smaller tolerance is raised to the default, and gives the default result. A
negative tolerance, or one that is not finite, gives {cpp:enumerator}`nxx::roots::error::non_finite`.
```

## A root near zero: set `abs_tol`

The exception is a root near zero. With the defaults, a method must first close in on the root's own scale before the
relative tolerance can stop it. For the root of $x - 10^{-300}$ in $[0, 1]$, bisection halves its way down through a
thousand powers of two: 1049 evaluations, almost all of them spent getting near zero rather than on digits. If a
distance of $10^{-12}$ from the root is close enough for you, say so:

```{literalinclude} ../examples/roots_tolerances.cpp
:start-after: // <near-zero>
:end-before: // </near-zero>
:dedent:
```

```{literalinclude} ../examples/roots_tolerances.cpp
:language: text
:start-after: near_zero_output = R"(
:end-before: )";
```

With `abs_tol = 1e-12` the same search takes 41 evaluations. The result is within $2 \times 10^{-12}$ of the root, which
for a root at $10^{-300}$ says only that it is near zero. That is the trade: you get an answer at the scale you named,
not at the root's own scale.

How far a method walks depends on `fn`. Brent's method reaches the root of $x - 10^{-300}$ in 3 evaluations, because its
interpolation lands on a straight line's root at once. But a function that jumps at 0, or that misses a root at 0 by
rounding, gives the interpolation nothing to work with, and Brent's method walks down to the type's smallest numbers as
bisection does: 1076 evaluations on a jump at 0 in $[-1, 2]$, where {cpp:struct}`nxx::roots::AndersonBjorck` takes 51.

In a multiprecision type the walk is far longer, because the type's exponent range is far wider: about $2 \times 10^9$
halvings in `cpp_bin_float_50`, long enough that a solve seems never to end. There, always set `abs_tol` when a root may
lie near zero. From a guess, an open method on a multiple root at 0, such as that of $x^3$, gives
{cpp:enumerator}`nxx::roots::error::not_converged` with the default `abs_tol`; setting `abs_tol` fixes that too.

```{note}
A root at exactly 0, where your function returns exactly 0, is found at once: a bracket that contains 0 inside it is
split there first.
```

## `abs_tol` acts at every magnitude

The absolute term is added everywhere, not only near zero. With `abs_tol = 1e-12`, a root at $10^{-9}$ is guaranteed
only to within about $2 \times 10^{-12}$, which is about three significant digits:

```{literalinclude} ../examples/roots_tolerances.cpp
:start-after: // <every-magnitude>
:end-before: // </every-magnitude>
:dedent:
```

```{literalinclude} ../examples/roots_tolerances.cpp
:language: text
:start-after: every_magnitude_output = R"(
:end-before: )";
```

A method may do better than its guarantee: Brent's method finds this straight line's root to within a unit or two in the
last place. But `abs_error`, $10^{-12}$ here, is what you can rely on. Choose `abs_tol` from the scale of the roots you
expect, not as a general accuracy setting.

## Fewer digits for less work: `rel_tol`

A larger `rel_tol` asks for fewer significant digits, so the method stops earlier. The saving is large for bisection,
which gains one binary digit per evaluation, and small for the interpolating methods, which gain digits quickly once
they are near the root:

```{literalinclude} ../examples/roots_tolerances.cpp
:start-after: // <fewer-digits>
:end-before: // </fewer-digits>
:dedent:
```

```{literalinclude} ../examples/roots_tolerances.cpp
:language: text
:start-after: fewer_digits_output = R"(
:end-before: )";
```

With `rel_tol = 1e-6`, bisection on $x^2 - 2$ over $[1, 2]$ takes 21 evaluations instead of 52, and the result is within
a relative $2 \times 10^{-6}$ (that is, $2\,\text{rel\_tol}$) of $\sqrt 2$.

The program needs:

```{literalinclude} ../examples/roots_tolerances.cpp
:start-after: // <includes>
:end-before: // </includes>
```

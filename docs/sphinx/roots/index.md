# Roots of functions

Many problems end in an equation that has to be solved numerically: the position of a planet in its orbit at a given
time, the volume a real gas fills at a given pressure, the temperature at which a liquid boils. Each can be written as

$$f(x) = 0$$

for some function $f$. A solution is called a *root* of $f$, and a numerical method finds one by evaluating $f$ at a
sequence of points that close in on it. `nxx::roots` does this for a function of one real variable. You give it the
function, as a lambda or any other callable, and a place to start. It gives back the root and a bound on its error.
[What you get back](#what-you-get-back) says exactly what that bound covers.

## Two ways to start

The surest start is a *bracket*: an interval $[a, b]$ at whose ends $f$ has opposite signs. If $f$ is continuous on
$[a, b]$, it must cross zero somewhere inside, so the bracket holds a root. A *bracketing method* keeps such an interval
at every step and shrinks it, so it always ends at a sign change of $f$. The simplest, bisection, halves the interval
with every evaluation. Brent's method, the one to start with, fits a curve through the last few points and jumps to
where the curve crosses zero. It bisects instead whenever the jump would leave the bracket or the bracket is not
shrinking fast enough. To pin the root of $x^2 - 2$ in $[1, 2]$ to a few units in the last place of a `double`,
bisection needs 52 evaluations and Brent's method 9 ([Costs](costs.md) compares every method).

Here is Brent's method on that problem. The function is a lambda, the bracket is a `Bracket`, and the method is named as
a template argument:

```{literalinclude} ../examples/roots_overview.cpp
:start-after: // <brent>
:end-before: // </brent>
:dedent:
```

It prints

```{literalinclude} ../examples/roots_overview.cpp
:language: text
:start-after: brent_output = R"(
:end-before: )";
```

The root is $\sqrt 2$ to within one unit in the last place of a `double`, and `abs_error` is the bound on its error.

If the bracket holds several roots, you get one of them. If $f$ has a pole inside, as $\tan x$ does at $\pi/2$, the sign
change may be the pole, and the method returns that instead ([Handle failures](handle-failures.md) shows how to tell).

Often you know only roughly where the root is. Then start from a *guess* with an *open method*. Newton's method replaces
$f$ by its tangent line at the current point and steps to where that line crosses zero; it needs the derivative $f'$,
either yours or one that `nxx::deriv` computes. Near a root where $f$ crosses zero at an angle this is very fast; far
from one it can wander. Numerixx's open methods therefore watch for the first step across which $f$ changes sign, and
finish inside that bracket, so a root they return carries the same guarantee as a bracketing method's. Unlike a
bracketing method they can fail, when their steps never cross a sign change, and the sign change they cross first need
not be the root nearest your guess.

When you have neither a bracket nor a good guess, `search` looks outward from a rough guess until $f$ changes sign and
returns that bracket, which you then pass to `solve` ([Find a bracket](find-a-bracket.md)). It cannot find a root where
$f$ touches zero without crossing it, as $x^2$ does at 0, unless `fn` happens to be exactly 0 at a point it tries.

Every method is called the same way, `solve<Method>(fn, start)`, where `fn` is your function and the start is a
`Bracket` or a guess. Every method works in any real type: `float`, `double`, or a multiprecision type ([Working with
Numerixx](../basics.md)).

## What you get back

A solve returns either a root or the reason it did not find one, in a `std::expected`. The root has two parts: `value`,
the root itself, and `abs_error`, a bound on its distance to a point where $f$, *as your program computes it*, changes
sign. With the default tolerances the bound is a few units in the last place: for a `double`, a relative error of at
most about $10^{-15}$. Where `fn` is exactly 0 at the root, `abs_error` is not 0 but the tolerance there: a point where
`fn` is 0 is only as good as `fn`'s evaluation. At a root at 0 the default tolerance is the type's smallest positive
number.

The bound is about $f$ as computed. Rounding in $f$ moves its sign change away from the exact root, by about $f$'s
rounding error divided by its slope, and no method can see that. For an $f$ computed in a few steps and crossing zero
steeply, the shift is smaller than the bound: about $2 \times 10^{-16}$ for Kepler's equation in the
[tutorial](tutorial.md). Where the terms of $f$ cancel, or its slope at the root is small, the shift can be far larger.
Take $(x - 1)^3$, written out and written as a product:

```{literalinclude} ../examples/roots_overview.cpp
:start-after: // <cubic>
:end-before: // </cubic>
:dedent:
```

```{literalinclude} ../examples/roots_overview.cpp
:language: text
:start-after: cubic_output = R"(
:end-before: )";
```

Written out as $x^3 - 3x^2 + 3x - 1$, the function is a sum of terms near 1 that cancel to almost nothing, and the
rounding in that sum is larger than $(x - 1)^3$ itself near the root. As computed, it comes out as 0 or with the wrong
sign, over and over, anywhere within about $7 \times 10^{-6}$ of 1. Brent's method lands on one of those points,
$4 \times 10^{-6}$ from 1, where the computed value is exactly 0, and returns it with `abs_error` $4 \times 10^{-16}$,
the tolerance there: a statement about $f$ as computed, not about the root. Written as a product, the same function is
computed with the right sign everywhere, and the root comes out within its bound of 1. When a root matters to many
digits, write $f$ so that it does not cancel near the root.

## Where to go next

The [tutorial](tutorial.md) solves Kepler's equation step by step. [Choose a method](choose-a-method.md) says when a
method other than Brent's does better, and the other guide pages each solve one task. The [reference](reference.md)
documents every name.

```{toctree}
:maxdepth: 1

tutorial
```

```{toctree}
:caption: Guide
:maxdepth: 1

choose-a-method
tolerances
find-a-bracket
restricted-domain
derivatives
other-types
handle-failures
```

```{toctree}
:caption: Reference
:maxdepth: 1

costs
reference
```

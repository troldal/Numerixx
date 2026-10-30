# Choose a method

`nxx::roots` has eleven methods, and Brent's method is the one to start with: it needs only a bracket and the values of
`fn`, it is fast on smooth functions, and its worst case is bounded. This page is for the cases where another method
does better. Which one depends on three things: what you start from (a bracket or a guess), whether you can compute the
derivative, and how much one evaluation of `fn` costs.

A *bracket* is an interval at whose ends `fn` has opposite signs; an *evaluation* is one call of `fn`. Every method is
called the same way, `solve<Method>(fn, start)`, so trying another method that starts from the same kind of thing is a
one-word change. Newton's and Halley's methods also need derivatives ([Give derivatives](derivatives.md)).

## From what you have to a method

| You have | Use | Why |
|---|---|---|
| A bracket, and `fn`'s values only | {cpp:struct}`nxx::roots::Brent` | The recommended method: fast on a smooth `fn`; at worst about the square of bisection's count |
| A bracket, and a worst case to bound | {cpp:struct}`nxx::roots::Toms748` | About Brent's cost, with at most about 4 evaluations per halving of the bracket |
| A bracket, `fn` strongly curved near one end (a logarithm, $1/x$), or a sign change at 0 | {cpp:struct}`nxx::roots::AndersonBjorck` | Far fewer evaluations than Brent's method there |
| A bracket, and a root of odd multiplicity, such as that of $(x - 1)^3$ | {cpp:struct}`nxx::roots::Bisection` | The interpolating methods slow down there; of them {cpp:struct}`nxx::roots::Ridders` is the cheapest |
| A bracket, a rough `fn`, and a count you must predict | {cpp:struct}`nxx::roots::Bisection` | One bit per evaluation, whatever `fn` is |
| A guess, and $f'$ (and $f''$) cheap to compute | {cpp:struct}`nxx::roots::Newton` (or {cpp:struct}`nxx::roots::Halley`) | Fewest evaluations from a good guess |
| A guess, and `fn`'s values only | {cpp:struct}`nxx::roots::Secant` | No derivative needed |
| A guess far from the root, or none | {cpp:func}`nxx::roots::search`, then a bracketing method | See [Find a bracket](find-a-bracket.md) |

The bracketing methods differ in how they choose the next point. Bisection takes the midpoint, so it gains one binary
digit per evaluation, whatever `fn` looks like. The others, Brent's method among them, *interpolate*: they fit a line or
a curve through the points so far and try where it crosses zero, which gains digits much faster once the bracket is
small and `fn` looks like a smooth curve. Brent's method and Toms748 fall back on bisection when interpolating does not
pay, which bounds their worst case: about the square of bisection's count for Brent's method, and about 4 evaluations
per halving of the bracket for Toms748.

Newton's and Halley's methods also take a bracket instead of a guess. The bracket then guards their steps as in Brent's
method: a step that would leave it, or that does not shrink it fast enough, is replaced by bisection. If `fn` gives only
values, prefer Brent's method or the secant method to Newton's method with a numerical derivative, which costs several
calls of `fn` per point ([Give derivatives](derivatives.md)).

```{note}
{cpp:struct}`nxx::roots::Ridders`, {cpp:struct}`nxx::roots::RegulaFalsi`, {cpp:struct}`nxx::roots::Illinois` and
{cpp:struct}`nxx::roots::Steffensen` are included as teaching examples. Each is a classic method, and each is beaten
by one of the methods above on the functions of the [costs table](costs.md), with one exception: on $x^{20} - 1$,
Ridders' method (18 evaluations) beats every other method that needs only values (Brent's takes 19), and loses only
to Newton's and Halley's methods, which need derivatives.
```

## Multiple roots

A root $r$ has *multiplicity* $m$ when $f$ behaves near it like $c\,(x - r)^m$ with $c \ne 0$. At a root of odd
multiplicity, such as that of $(x - 1)^3$, $f$ still changes sign, but it is so flat that the interpolating methods lose
their speed: Brent's method needs 156 evaluations on it, where bisection needs 54 ([costs](costs.md)). Use bisection
there.

A root of even multiplicity, such as that of $(x - 1)^2$, has no sign change at all. No bracket shows it, and the
methods find it only where `fn` happens to evaluate to exactly 0.

## Measure the cost on your own function

The table is a guide; your own `fn` decides. The methods do not report how many evaluations they made, but you can count
them: wrap `fn` in a lambda that increments a counter and then calls it.

```{literalinclude} ../examples/roots_choose_a_method.cpp
:start-after: // <count>
:end-before: // </count>
```

Then run each candidate method on your `fn` and bracket. Here the function is $x^3 - 2x - 5$, on $[2, 3]$:

```{literalinclude} ../examples/roots_choose_a_method.cpp
:start-after: // <compare>
:end-before: // </compare>
:dedent:
```

It prints

```{literalinclude} ../examples/roots_choose_a_method.cpp
:language: text
:start-after: output = R"(
:end-before: )";
```

Bisection pays about one evaluation per binary digit of the root; the interpolating methods need about a sixth as many
here. The [costs table](costs.md) gives every method's count on five such functions.

The program needs:

```{literalinclude} ../examples/roots_choose_a_method.cpp
:start-after: // <includes>
:end-before: // </includes>
```

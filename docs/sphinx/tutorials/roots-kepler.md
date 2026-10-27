# Solving Kepler's equation

This tutorial solves an equation that has no closed-form solution, step by step: Kepler's equation, which gives a
planet's position in its orbit. You will write the function, find a bracket for its root, solve it, read the result
and its error bound, handle a failure, and solve it again from a guess with Newton's method. It continues where
[Getting started](getting-started.md) ends.

## The equation

For an orbit of eccentricity e, the mean anomaly M (which grows uniformly with time) and the eccentric anomaly E (which
gives the position) satisfy

```text
E - e sin E = M
```

Given e and M, we want E. Write the equation as a function whose root is E, with e and M captured by the lambda:

```{literalinclude} ../examples/roots_kepler.cpp
:start-after: // <function>
:end-before: // </function>
:dedent:
```

The program needs these headers:

```{literalinclude} ../examples/roots_kepler.cpp
:start-after: // <includes>
:end-before: // </includes>
```

## A bracket, and the root in it

A bracketing method needs an interval in which the function changes sign. The equation gives one: E - M = e sin E,
and |sin E| <= 1, so the root lies within e of M, in [M - e, M + e]. The function's slope, 1 - e cos E, is at least
1 - e > 0, so it increases everywhere and there is exactly one root. Solve with Brent's method:

```{literalinclude} ../examples/roots_kepler.cpp
:start-after: // <bracket>
:end-before: // </bracket>
:dedent:
```

{cpp:func}`nxx::roots::solve` returns a `std::expected`: test it before reading it. The result, a
{cpp:struct}`nxx::roots::Root`, has two members, read by name: `value`, the root, and `abs_error`, a bound on its
distance to the point where the function, as computed, changes sign. With the default tolerances that bound is about
the precision of a double at the root.

## Checking the result

E for e = 0.5 and M = 1 is 1.4987011335178483140... (computed to 100 digits). Print the actual error beside
`abs_error`, and substitute E back into the equation:

```{literalinclude} ../examples/roots_kepler.cpp
:start-after: // <check>
:end-before: // </check>
:dedent:
```

The actual error is within `abs_error`, give or take the function's own rounding: `std::sin` is accurate to a few
units in the last place, which moves the computed sign change by about 2e-16 here, and no method can see beyond that.
The residual `kepler(E)` is close to 0.

## When there is no sign change

Give a bracket that does not contain the root, [2, 3], where the function is positive at both ends:

```{literalinclude} ../examples/roots_kepler.cpp
:start-after: // <failure>
:end-before: // </failure>
:dedent:
```

The result holds an error instead of a root, {cpp:enumerator}`nxx::roots::error::no_sign_change`, which formats as its
name. Numerixx never throws; every failure comes back this way. [Handle failures](../how-to/roots-handle-failures.md)
shows how to recover.

## From a guess, with Newton's method

Newton's method needs the derivative, 1 - e cos E, and starts from a guess instead of a bracket; M is a good one.
Return the value and the derivative together, as a `std::pair` of two doubles (a literal `1` instead of `1.0` would
make the second an `int`, which `solve` rejects):

```{literalinclude} ../examples/roots_kepler.cpp
:start-after: // <newton>
:end-before: // </newton>
:dedent:
```

Newton's method iterates from the guess until two iterates have opposite signs, then finishes inside that bracket, so
its result carries the same kind of bound as Brent's. It gives the same root to the same precision.

## Next

- [Choose a method](../how-to/roots-choose-a-method.md) for the other methods, and what each costs.
- [Give derivatives](../how-to/roots-derivatives.md) for the other ways to pass f', including a numerical one.
- [Find a bracket](../how-to/roots-find-a-bracket.md) when no bracket is known.

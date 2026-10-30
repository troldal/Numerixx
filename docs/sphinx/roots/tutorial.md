# Solving Kepler's equation

A planet moves around the sun on an ellipse, fast near the sun and slowly far from it. Its position at a given time
follows from an equation Kepler wrote down four centuries ago, and that equation has no solution in closed form: it has
to be solved numerically. This tutorial solves it step by step. You will write the equation as a function, find an
interval that must hold the solution, solve it with Brent's method, check the result, see what happens when the interval
is wrong, and solve it again from a guess with Newton's method. It continues where [Getting
started](../getting-started.md) ends; [Roots of functions](index.md) gives the ideas behind each step.

## The equation

Two angles describe where the planet is. The *mean anomaly* $M$ grows uniformly with time, so it is what you know for a
given date. The *eccentric anomaly* $E$ gives the position on the ellipse, so it is what you want. For an orbit of
eccentricity $e$ (0 for a circle, close to 1 for a long thin ellipse) they are related by Kepler's equation:

$$E - e \sin E = M$$

Given $e$ and $M$, you want $E$. To use a root finder, move everything to one side, so that $E$ is a root of

$$f(E) = E - e \sin E - M$$

In code, $f$ becomes a lambda, `kepler`. It takes $E$ and captures $e$ and $M$, so the solver sees a function of one
variable:

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

A *bracket* is an interval at whose ends the function has opposite signs, so that it crosses zero somewhere inside.
Kepler's equation gives one directly. It says $E - M = e \sin E$, and $|\sin E| \le 1$, so $E$ lies within $e$ of $M$:
in $[M - e, M + e]$. There is exactly one root there, because the slope of $f$, $1 - e \cos E$, is at least $1 - e > 0$:
$f$ increases everywhere, and crosses zero once.

Pass the function and the bracket to {cpp:func}`nxx::roots::solve`, with the method as a template argument. Brent's
method is the one to start with:

```{literalinclude} ../examples/roots_kepler.cpp
:start-after: // <bracket>
:end-before: // </bracket>
:dedent:
```

`solve` returns a `std::expected`: either a root, or the reason there is none. The `if (r)` tests which, and only then
does the code read the result. The result, a {cpp:struct}`nxx::roots::Root`, has two members. `value` is the root.
`abs_error` is a bound on its distance to the point where `kepler`, as the program computes it, changes sign. With the
default tolerances that bound is a few units in the last place of `value`: here at most about $1.3 \times 10^{-15}$. The
program prints $E \approx 1.49870113351785$, to 17 digits.

## Checking the result

For $e = 0.5$ and $M = 1$, $E$ is 1.4987011335178483140... (computed to 100 digits by an independent program). A
`double` holds only about 16 of those digits, and the root's error lies in the last of them, so a reference rounded to a
`double` would hide it. The program therefore keeps the reference as the nearest `double` and the rest:

```{literalinclude} ../examples/roots_kepler.cpp
:start-after: // <reference>
:end-before: // </reference>
:dedent:
```

Print the actual error beside `abs_error`, and put $E$ back into the equation:

```{literalinclude} ../examples/roots_kepler.cpp
:start-after: // <check>
:end-before: // </check>
:dedent:
```

The actual error is within `abs_error`, give or take about $2 \times 10^{-16}$. That small extra comes from the function
itself. `std::sin` is accurate to a few units in the last place, so the computed $f$ is slightly off, which moves its
sign change by about its rounding error divided by its slope. No method can see beyond that, which is why `abs_error`
describes the function as computed. The residual, `kepler(E)`, is close to 0. Its size follows from the distance to the
root: $f$'s slope is at most $1 + e = 1.5$, so $|f(E)|$ is at most 1.5 times (`abs_error` plus $2 \times 10^{-16}$),
plus $f$'s own rounding of about $2 \times 10^{-16}$: a few times $10^{-15}$ at most.

The exact digits printed for the error and the residual can differ between platforms, because each standard library
computes `std::sin` in its own way; the bounds above hold on all of them.

## When there is no sign change

What if the bracket is wrong? Give it $[2, 3]$, which lies above the root, so that $f$ is positive at both ends:

```{literalinclude} ../examples/roots_kepler.cpp
:start-after: // <failure>
:end-before: // </failure>
:dedent:
```

The result now holds an error instead of a root, {cpp:enumerator}`nxx::roots::error::no_sign_change`, and the program
prints `no root: no_sign_change`: an error formats as its name. Nothing is thrown; every failure comes back this way.
[Handle failures](handle-failures.md) shows how to recover from one.

## From a guess, with Newton's method

Sometimes you have no bracket, only a guess. For Kepler's equation, $M$ itself is a good guess, since $E$ lies within
$e$ of it. Newton's method starts from a guess, but it needs the derivative,

$$f'(E) = 1 - e \cos E$$

Return the value and the derivative together from one lambda, as a `std::pair` of two doubles:

```{literalinclude} ../examples/roots_kepler.cpp
:start-after: // <newton>
:end-before: // </newton>
:dedent:
```

Both elements must be `double`. A bare literal, as in `std::pair{f, 1}`, would make the second element an `int`, and
`solve` would not compile.

Newton's method steps from the guess along the tangent of $f$. As soon as $f$ has opposite signs at two consecutive
steps, the root is bracketed, and the method finishes inside that bracket. That is why its result carries the same kind
of bound as Brent's: you get the same root, to the same precision.

## Next

- [Choose a method](choose-a-method.md): the other methods, and when each does better than Brent's.
- [Give derivatives](derivatives.md): other ways to pass $f'$, including a numerical one.
- [Find a bracket](find-a-bracket.md): what to do when you have neither a bracket nor a good guess.

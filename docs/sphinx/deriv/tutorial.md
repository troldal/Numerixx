# The slope of a computed function

Numerical derivatives earn their keep on functions that have no formula you can differentiate. This tutorial takes one
such function and differentiates it step by step: the position of a planet as a function of time, which is itself
computed by solving Kepler's equation. You will write the function, take its derivative with the default method, check
the result against a closed form, get a more accurate result with an error estimate from Ridders' method, and take the
second derivative. It builds on the [roots tutorial](../roots/tutorial.md), which solves Kepler's equation; you do not
need to have read it.

## The function

A planet's *mean anomaly* $M$ grows uniformly with time; its *eccentric anomaly* $E$ gives its position on the ellipse.
For an orbit of eccentricity $e$ they are related by Kepler's equation,

$$E - e \sin E = M$$

which has no solution in closed form. For each $M$, $E(M)$ is the root of $E - e \sin E - M$, which lies in
$[M - e, M + e]$. So $E(M)$ is a function you can evaluate, by a root solve, but not write down. How fast does the
eccentric anomaly advance with time? That is the derivative $dE/dM$.

In code, $E(M)$ is a lambda that solves Kepler's equation with `nxx::roots`. The examples write `rt` for `nxx::roots`
and `dv` for `nxx::deriv`:

```{literalinclude} ../examples/deriv_kepler.cpp
:start-after: // <function>
:end-before: // </function>
:dedent:
```

If the solve fails, the lambda returns NaN, which the derivative does not take for a number: a simple method reports an
error, Ridders' method tries smaller steps. The program needs both headers:

```{literalinclude} ../examples/deriv_kepler.cpp
:start-after: // <includes>
:end-before: // </includes>
```

## The derivative

Pass the function and the point to {cpp:func}`nxx::deriv::derivative`. Without a method, it uses the *central
difference*, $(f(M + h) - f(M - h)) / 2h$, with a default step $h$ of about $6 \times 10^{-6}$ here ([Numerical
derivatives](index.md#the-idea-and-its-catch) explains the choice):

```{literalinclude} ../examples/deriv_kepler.cpp
:start-after: // <central>
:end-before: // </central>
:dedent:
```

Like every algorithm in Numerixx, `derivative` returns a `std::expected`, so the code tests it before reading `value`.
At $M = 1$ it prints $dE/dM$ to 15 decimals; the first ten are 1.0373620219.

## Checking it

This function has a closed-form derivative after all, which makes it a good test. Differentiate Kepler's equation with
respect to $M$: $E' - e \cos E \, E' = 1$, so

$$\frac{dE}{dM} = \frac{1}{1 - e \cos E}$$

Compare:

```{literalinclude} ../examples/deriv_kepler.cpp
:start-after: // <exact>
:end-before: // </exact>
:dedent:
```

The relative error is far below $10^{-9}$, about what the central difference gives on a smooth function. That holds
although each value of $E(M)$ is itself only accurate to a few units in the last place, since it comes from a root
solve; those errors are divided by the step and add up to no more than about $3 \times 10^{-10}$. The exact digits
depend on the platform's `std::sin`, so the program checks the error against that bound rather than a printed value.

## A better result, with an error estimate

The central difference cannot tell you how accurate it is. *Ridders' method* can: it computes the central difference at
a sequence of shrinking steps and extrapolates the results to a step of zero, estimating its error on the way. It is
more accurate too. Ask for it by name:

```{literalinclude} ../examples/deriv_kepler.cpp
:start-after: // <ridders>
:end-before: // </ridders>
:dedent:
```

The result now has an `abs_error`, Ridders' estimate of its error, a `std::optional` that this method always fills. The
estimate is below $10^{-11}$, and the actual error is smaller still: more than ten times smaller than the central
difference's, for about a dozen evaluations of $E(M)$ instead of two.

The estimate assumes that each value of $E(M)$ is accurate to about the type's epsilon, and a root solve gives a few
times that. Here it still covers the actual error, on every platform the program is tested on; for a function computed
less accurately, check the digits as [Handle failures](handle-failures.md#a-result-that-is-not-accurate-enough)
describes.

## The second derivative

The derivative's order is the second template argument. Differentiating the closed form once more gives
$d^2E/dM^2 = -e \sin E \, (dE/dM)^3$, to check against:

```{literalinclude} ../examples/deriv_kepler.cpp
:start-after: // <second>
:end-before: // </second>
:dedent:
```

The second derivative is harder: rounding is now divided by $h^2$. Ridders' estimate is below $10^{-10}$, and again
covers the actual error.

## Next

- [Choose a method](choose-a-method.md): what the other methods cost and give.
- [Higher derivatives](higher-derivatives.md): orders 2 to 4, and how fast their accuracy falls.
- [Pass a derivative as a function](derivative-of.md): give $dE/dM$ to an algorithm that takes a function.

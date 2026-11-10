# The slope of a computed function

Numerical derivatives are for functions that have no formula you can differentiate. This tutorial differentiates one
such function, the position of a planet in its orbit, which is itself computed by solving Kepler's equation.

## The function

The *mean anomaly* $M$ measures time along the orbit, and the *eccentric anomaly* $E$ the planet's position; for an
orbit of eccentricity $e$, Kepler's equation relates them ([the roots tutorial](../roots/tutorial.md) says more):

$$E - e \sin E = M$$

For each $M$, $E(M)$ is the root of $E - e \sin E - M$, which lies in $[M - e, M + e]$. You can evaluate it, but not
write it down. Its derivative $dE/dM$ says how fast the position advances with time.

In code, $E(M)$ is a lambda that solves Kepler's equation with `nxx::roots`. The examples write `rt` for `nxx::roots`
and `dv` for `nxx::deriv`:

```{literalinclude} ../examples/deriv_kepler.cpp
:start-after: // <function>
:end-before: // </function>
:dedent:
```

If the solve fails, the lambda returns NaN, which no method takes for a number ([Handle failures](handle-failures.md)).
The program needs both headers:

```{literalinclude} ../examples/deriv_kepler.cpp
:start-after: // <includes>
:end-before: // </includes>
```

## The derivative

Pass the function and the point to {cpp:func}`nxx::deriv::derivative`. Without a method, it uses a *stencil*, a fixed
formula with a step $h$: the *central difference*, $(f(M + h) - f(M - h)) / 2h$. Its default step is about
$6 \times 10^{-6}$ here ([Numerical derivatives](index.md#the-idea-and-its-catch) explains the choice):

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

The relative error is below $10^{-9}$. Each $E(M)$ is off from the root solve by a few units in its last digit, which
adds at most about $3 \times 10^{-10}$. The digits depend on the platform's `std::sin`, so the program checks that bound
rather than printed output.

## A better result, with an error estimate

The central difference cannot tell you how accurate it is. *Ridders' method* can: it extrapolates the central difference
from a sequence of shrinking steps to a step of zero, and estimates its error on the way. Ask for it by name:

```{literalinclude} ../examples/deriv_kepler.cpp
:start-after: // <ridders>
:end-before: // </ridders>
:dedent:
```

The result now has an `abs_error`, the method's estimate of its error. The estimate is below $10^{-11}$, and covers the
actual error. That error is more than ten times smaller than the central difference's, for about a dozen evaluations of
$E(M)$ instead of two.

The estimate assumes that each value of $E(M)$ is accurate to about the type's epsilon. A root solve gives a few times
that, and the estimate still covers the error here. For a function computed less accurately, check the digits as
[Handle failures](handle-failures.md#a-result-that-is-not-accurate-enough) describes.

## The second derivative

The derivative's order is the second template argument. Differentiating the closed form once more gives
$d^2E/dM^2 = -e \sin E \, (dE/dM)^3$, to check against:

```{literalinclude} ../examples/deriv_kepler.cpp
:start-after: // <second>
:end-before: // </second>
:dedent:
```

Ridders' estimate is below $10^{-10}$, and again covers the actual error. A second derivative keeps fewer digits than
the first ([Higher derivatives](higher-derivatives.md)).

## Next

- [Choose a method](choose-a-method.md): what the other methods cost and give.
- [Higher derivatives](higher-derivatives.md): orders 2 to 4, and how fast their accuracy falls.
- [Pass a derivative as a function](derivative-of.md): give $dE/dM$ to an algorithm that takes a function.

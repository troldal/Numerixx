# Integrals of functions

Many quantities are integrals that cannot be written in closed form: the period of a pendulum swinging through a large
angle, the probability that a normally distributed measurement exceeds a limit, the heat needed to warm a substance
whose heat capacity depends on temperature. Sometimes the integrand has no antiderivative you can write down; sometimes
it exists only as a program, a function you can call but not differentiate or integrate by hand. Either way the integral

$$I = \int_a^b f(x)\,dx$$

has to be computed from values of $f$. `nxx::quadrature` does this for a real function of one variable. You give it the
function, as a lambda or any other callable, and the ends of the interval; it gives back the integral and an estimate of
its error.

## A weighted sum, and an error estimate

A *quadrature rule* replaces the integral by a weighted sum of values of $f$ at chosen points, the *nodes*:

$$\int_a^b f(x)\,dx \approx \sum_{i=1}^{N} w_i\, f(x_i)$$

The trapezoidal rule, with equally spaced nodes, is the simplest. A *Gauss rule* chooses its nodes and weights so that
the sum is exact for every polynomial of degree up to $2N - 1$, which makes it very accurate for a smooth $f$ with few
points. A rule alone does not tell you how accurate its sum is. To estimate that, compute two sums of different accuracy
and take their difference: the *Gauss-Kronrod* pair does it with 15 nodes, 7 of which give a Gauss rule of their own, so
the second sum costs nothing extra.

The default method, adaptive Gauss-Kronrod, starts by cutting the interval into 10 equal *panels* and applying the pair
on each, 150 evaluations of $f$. Then, while the sum of the panels' estimates is too large, it cuts the panel with the
largest estimate in half. The work goes where $f$ is hard: a peak, a kink, a singularity at an end. Here it integrates
$4 / (1 + x^2)$ over $[0, 1]$, which is $\pi$. The examples write `qd` for `nxx::quadrature`:

```{literalinclude} ../examples/quadrature_overview.cpp
:start-after: // <first>
:end-before: // </first>
:dedent:
```

```{literalinclude} ../examples/quadrature_overview.cpp
:language: text
:start-after: first_output = R"(
:end-before: )";
```

The ten panels were enough: the estimate met the request at once, and no panel was cut. The value is $\pi$ to the last
digit of a `double`, while `abs_error` says $3.5 \times 10^{-14}$. That is typical. The estimate is the 7-point rule's
error, not the 15-point result's, so it usually overstates the error by orders of magnitude; and it never falls below a
*rounding floor*, about $10^{-14}$ relative in a `double`, under which the rounding of the values and of the sum would
dominate any estimate.

## Four methods

| Method | When to use it |
|---|---|
| `GaussKronrod<15>`, the default | Most integrals on a finite interval: it finds the hard places by itself |
| `TanhSinh` | A singularity at an end, such as $\log x$ or $1 / \sqrt x$ at 0; many digits in a multiprecision type; an infinite interval, which no other method takes |
| `GaussLegendre<N>` | A fixed rule of $N$ points, without an error estimate: an integral inside a root solve, an optimiser or another integral, where the result must change smoothly with the parameters of $f$ |
| `Romberg` | Teaching and comparison: trapezoidal sums with halved steps, extrapolated |

Name the method as a template argument, `integrate<Method>(fn, a, b)`. On the same integral:

```{literalinclude} ../examples/quadrature_overview.cpp
:start-after: // <methods>
:end-before: // </methods>
:dedent:
```

```{literalinclude} ../examples/quadrature_overview.cpp
:language: text
:start-after: methods_output = R"(
:end-before: )";
```

The 10-point Gauss rule is accurate to 13 digits from 10 evaluations, but it returns a plain number: nothing tells you
so. `TanhSinh` is as accurate as the default here, at a similar cost; it pays off when $f$ is singular at an end, where
it can need ten times fewer evaluations ([Choose a method](choose-a-method.md)). Romberg's method stops after 33
evaluations, because its estimate, $1.2 \times 10^{-8}$, meets the request; its actual error is more than a hundred
times smaller. Every estimate here is larger than the actual error, as it should be, but none is a guaranteed bound
([What you get back](#what-you-get-back)).

## The call

Every method is called the same way, `integrate<Method>(fn, a, b, criteria)`. `fn` takes an `x` and returns $f(x)$;
the ends `a` and `b` set the type the integral is computed in, so write them as `0.0` and `1.0` for a `double`, and
they must have the same type. The ends may come in either order: `integrate(fn, b, a)` is the negated integral, and
`integrate(fn, a, a)` is 0 without calling `fn`. An infinite end, `std::numeric_limits<double>::infinity()`, works with
`TanhSinh` alone ([Infinite intervals](infinite-intervals.md)). The criteria are optional: by default the methods ask
for a relative error of about $1.5 \times 10^{-8}$ in a `double`, half its digits, not all of them
([Set tolerances](tolerances.md)). `integrate(m, fn, a, b)` takes the method as a value instead, for a method with
parameters, such as one of your own ([A method of your own](method-of-your-own.md)).

## What you get back

A call with a method that estimates its error returns a `std::expected`: either an
{cpp:struct}`nxx::quadrature::Integral`, whose members are `value` and `abs_error`, or a
{cpp:struct}`nxx::quadrature::Failure`, which says why there is no integral. `abs_error` estimates the distance from
`value` to the exact integral. It is an estimate, not a bound: no method that samples $f$ at finitely many points can
bound its error, since a narrow peak between the points goes unseen. Compare it with the accuracy you need before you
rely on `value` ([Working with Numerixx](../basics.md#what-a-result-holds) says how to read a result).

A failure has a `reason`, {cpp:enumerator}`nxx::quadrature::error::non_finite` or
{cpp:enumerator}`nxx::quadrature::error::not_converged`, and, where the method ran out of work, its `best` estimate.
The integral of $1/x$ over $[0, 1]$ is infinite:

```{literalinclude} ../examples/quadrature_overview.cpp
:start-after: // <failure>
:end-before: // </failure>
:dedent:
```

```{literalinclude} ../examples/quadrature_overview.cpp
:language: text
:start-after: failure_output = R"(
:end-before: )";
```

The default method kept cutting panels near 0 until its estimate stopped improving, and gave up. Its `best` holds a
finite number, but for an integral that does not exist, that number means nothing. [Handle failures](handle-failures.md)
shows what each reason means and what to try next. `GaussLegendre<N>` returns a `std::expected` that holds a plain
number or an `error`, since it has no estimate to fail on.

## What it does not do

`nxx::quadrature` integrates a function you can call at any point of the interval. It does not integrate tabulated data,
complex or vector-valued functions, or functions of several variables, though you can nest one integral inside another
([Choose a method](choose-a-method.md#an-integral-inside-another-algorithm)). It has no method for principal values or
for oscillating integrands on an infinite interval, such as $\sin x / x$ on $[0, \infty)$.

## Where to go next

The [tutorial](tutorial.md) computes the period of a pendulum, step by step. [Choose a method](choose-a-method.md) says
which method suits which integrand, and the other guide pages each solve one task. [Costs](costs.md) compares the
methods on a few functions, and the [reference](reference.md) documents every name.

```{toctree}
:maxdepth: 1

tutorial
```

```{toctree}
:caption: Guide
:maxdepth: 1

choose-a-method
tolerances
rough-integrands
infinite-intervals
high-precision
handle-failures
method-of-your-own
```

```{toctree}
:caption: Reference
:maxdepth: 1

costs
reference
```

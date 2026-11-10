# Numerical derivatives

Sometimes you need the derivative of a function you can only evaluate: a simulation, interpolated data, or the root of
an equation. `nxx::deriv` estimates its first to fourth derivative at a point $x$ from its values near $x$. You pass the
function, as any callable, and the point; you get back the derivative and, from some methods, an estimate of its error.

## The idea, and its catch

The *central difference* evaluates $f$ one *step* $h$ to each side
of $x$:

$$f'(x) \approx \frac{f(x + h) - f(x - h)}{2h}$$

Its *truncation error* is the error of the formula itself. It falls fourfold each time you halve the step, so a small
step seems best. The catch is rounding. Each value of $f$ is rounded to about $\varepsilon |f|$, where $\varepsilon$ is
the type's *machine epsilon*, $2.2 \times 10^{-16}$ for a `double`. The formula subtracts two nearly equal values and
divides by $2h$. So its *rounding error* is about $\varepsilon |f| / h$, and it grows as the step shrinks.

```{eval-rst}
.. plot:: plots/deriv_trade_off.py
   :alt: On logarithmic axes, with smaller steps to the right, the relative error of the central difference of the
         square root at 2 falls in a straight line from 3.5e-2 at a step of 1 to 4.1e-12 at 5e-6, then rises raggedly
         to 2.5e-3 at 1e-13. A dashed vertical line marks the default step, 1.2e-5, at the bottom of the V. A
         horizontal line marks the error of Ridders' method, 5.8e-14, below the whole curve.
   :caption: The relative error of the central difference of :math:`\sqrt{x}` at 2 against its step :math:`h`, with
             the default step and the error of Ridders' method.
```

The figure shows both errors on $\sqrt x$ at 2. Down to about $h = 10^{-5}$, truncation rules, and the error falls in a
straight line. Below that, rounding takes over: at $h = 10^{-13}$ fewer than four digits are left. The best step
balances the two, and keeps about two thirds of the type's digits.

You need not find that step yourself. The *default step* is $h = \varepsilon^{1/3} \max(|x|, 1)$, the balance for a
function that varies on the *scale* $\max(|x|, 1)$. Here it is $1.2 \times 10^{-5}$, the dashed line. The examples
write `dv` for `nxx::deriv`:

```{literalinclude} ../examples/deriv_overview.cpp
:start-after: // <default>
:end-before: // </default>
:dedent:
```

```{literalinclude} ../examples/deriv_overview.cpp
:language: text
:start-after: default_output = R"(
:end-before: )";
```

That is about as good as the best step in the figure. A function that varies much faster needs a step of its own ([Choose a
step](step.md)).

## Two kinds of method

A *stencil* is a fixed formula of this kind: a few points around $x$, a weight for each, and a step. The central
difference is `Central<2>`, the default. `Central<4>` to `Central<10>` use more points; on a smooth function
`Central<4>` and `Central<6>` are more accurate. `Forward<P>` and `Backward<P>` use points on one side of $x$ only. A
stencil costs a fixed number of evaluations, but cannot tell you how accurate it is.

*Ridders' method* computes the central difference at a sequence of shrinking steps, and extrapolates the results to a
step of zero. It finds a good step by itself and estimates its error, typically for 10 to 16 evaluations of `fn`.
`ForwardRidders` and `BackwardRidders` do the same from one side of $x$; with `Ridders`, they are *the Ridders
methods*:

```{literalinclude} ../examples/deriv_overview.cpp
:start-after: // <ridders>
:end-before: // </ridders>
:dedent:
```

```{literalinclude} ../examples/deriv_overview.cpp
:language: text
:start-after: ridders_output = R"(
:end-before: )";
```

Its error, the horizontal line in the figure, is about a hundred times smaller than the central difference's. Its
estimate, `abs_error`, is a few times the actual error. [Choose a method](choose-a-method.md) compares the methods.

Every method is called as `derivative<Method, N>(fn, x)`, where `N` is the derivative's order, 1 (the default) to 4.
`derivative(fn, x, h)` takes a step of your own, and `derivative_of<Method, N>(fn)` turns `fn` into its derivative
([Pass a derivative as a function](derivative-of.md)).

## What you get back

A call returns a `std::expected`. On success it holds an {cpp:struct}`nxx::deriv::Estimate` with two members: `value`,
the derivative, and `abs_error`, a `std::optional` estimate of its error, not a bound. The Ridders methods fill it; the
stencils leave it empty. On failure it holds {cpp:enumerator}`nxx::deriv::error::non_finite`, for instance when `fn`
returned NaN where a method needed a value ([Handle failures](handle-failures.md)). Only values that are not finite are
detected: finite but wrong values give a wrong derivative, without an error. [Working with
Numerixx](../basics.md#what-a-result-holds) says how to read a result.

## Tutorial

The tutorial differentiates a function that is itself computed by a root solve.

```{toctree}
:maxdepth: 1

tutorial
```

## Guide

Each section solves one task. Start with [Choose a method](choose-a-method.md).

```{toctree}
:maxdepth: 1

choose-a-method
step
one-sided
higher-derivatives
derivative-of
handle-failures
```

## Reference

The methods' costs side by side, the API, and the book on `nxx::deriv`, which explains the theory in full.

```{toctree}
:maxdepth: 1

costs
reference
book
```

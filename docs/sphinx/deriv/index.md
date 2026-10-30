# Numerical derivatives

Sometimes you need the derivative of a function you can only evaluate: a simulation, a table with interpolation, a model
too long to differentiate by hand, or a function that is itself the result of a numerical method, such as the root of an
equation. `nxx::deriv` estimates the first to fourth derivative of such a function at a point $x$ from its values at a
few points near $x$. You give it the function, as a lambda or any other callable, and the point; it gives back the
derivative, and, if you ask for Ridders' method, an estimate of its error.

## The idea, and its catch

The derivative is the limit of a difference quotient. The *central difference* takes one step $h$ to each side of $x$:

$$f'(x) \approx \frac{f(x + h) - f(x - h)}{2h}$$

Taylor's theorem shows that its *truncation error*, the error of the formula itself, is about
$\tfrac{h^2}{6} |f'''(x)|$: halve the step and the error falls by four. So a small step seems best. The catch is
rounding. Each value of $f$ carries a rounding error of about $\varepsilon |f|$, where $\varepsilon$ is the type's
machine epsilon ($2.2 \times 10^{-16}$ for a `double`). The formula subtracts two nearly equal values, which keeps those
errors, and divides by $2h$, which magnifies them. The *rounding error* of the result is about $\varepsilon |f| / h$,
and it grows as the step shrinks.

Here is the central difference of $\sqrt x$ at 2, for steps from $10^{-1}$ down to $10^{-11}$:

```{literalinclude} ../examples/deriv_overview.cpp
:start-after: // <steps>
:end-before: // </steps>
:dedent:
```

```{literalinclude} ../examples/deriv_overview.cpp
:language: text
:start-after: steps_output = R"(
:end-before: )";
```

Down to about $h = 10^{-5}$ the error is the truncation error, falling a hundredfold with each tenfold smaller step.
Below it, rounding takes over, and the error grows again: at $h = 10^{-11}$ only five digits are left. The best step
balances the two, and gives about two thirds of the type's digits, 10 or 11 of a `double`'s 16.

You rarely have to find that step yourself. Without one, the central difference takes
$h = \varepsilon^{1/3} \max(|x|, 1)$ for a first derivative, the balance for a function that varies on the scale
$\max(|x|, 1)$ (each method and order has its own default). The examples write `dv` for `nxx::deriv`:

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

That is as good as the best step in the table. A function that varies much faster than that scale needs a step of its
own ([Choose a step](step.md)).

## Two kinds of method

A *stencil* is a fixed formula of this kind: a few points around $x$ and a weight for each. The central difference is
the simplest, `Central<2>`, with two points and an error that falls as $h^2$. `Central<4>`, `Central<6>`, `Central<8>`
and `Central<10>` use more points, and their error falls as $h^4$ to $h^{10}$; `Forward<P>` and `Backward<P>`, for $P$
from 1 to 6, use points on one side of $x$ only, for a function that is defined only there. A stencil costs a known
number of evaluations of $f$, but it cannot tell you how accurate its result is.

*Ridders' method* computes the central difference at a sequence of shrinking steps and extrapolates the results to a
step of zero, stopping when rounding stops the improvement. It finds a good step by itself, typically keeps 13 to 15
digits of a `double`, and estimates its own error, at the cost of typically 10 to 16 evaluations of $f$ for a first
derivative ({cpp:struct}`nxx::deriv::Ridders`). `ForwardRidders` and `BackwardRidders` do the same from one side of $x$;
together they are *Ridders' methods*:

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

Its error is about a hundred times smaller than the central difference's, and its estimate, `abs_error`, is a few times
larger than the actual error. [Choose a method](choose-a-method.md) compares the methods.

Every method is called the same way, `derivative<Method, N>(fn, x)`, where `fn` is your function, `x` the point and `N`
the order of the derivative, 1 (the default) to 4. `derivative(fn, x, h)` takes a step of your own.
`derivative_of<Method, N>(fn)` turns `fn` into its derivative, a function you can pass on, for example to Newton's
method in `nxx::roots` ([Pass a derivative as a function](derivative-of.md)).

## What you get back

A call returns a `std::expected`: either an {cpp:struct}`nxx::deriv::Estimate`, or the error
{cpp:enumerator}`nxx::deriv::error::non_finite` when no finite estimate exists, for instance when `fn` returns NaN at a
stencil's point, or at every step Ridders' method tries ([Handle failures](handle-failures.md)). Only values that are
not finite are detected: a function that gives finite but wrong numbers gives a wrong derivative. The estimate has two
members: `value`, the derivative, and `abs_error`, a `std::optional`. Ridders' methods fill it with an estimate of the
error; the stencils leave it empty, since they cannot estimate their error without more evaluations. The estimate is not
a guaranteed bound: it assumes that `fn` is accurate to about the type's epsilon, and it can be too small for a function
that is noisier than that. [Working with Numerixx](../basics.md#what-a-result-holds) says how to read it.

## Where to go next

The [tutorial](tutorial.md) differentiates a function that is itself computed by a root solve, step by step. [Choose a
method](choose-a-method.md) says which method suits which function, and the other guide pages each solve one task. The
book on `nxx::deriv` explains the theory in full ([the book](book.md)), and the [reference](reference.md) documents
every name.

```{toctree}
:maxdepth: 1

tutorial
```

```{toctree}
:caption: Guide
:maxdepth: 1

choose-a-method
step
one-sided
higher-derivatives
derivative-of
handle-failures
```

```{toctree}
:caption: Reference
:maxdepth: 1

accuracy
reference
book
```

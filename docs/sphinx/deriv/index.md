# Numerical derivatives

Sometimes you need the derivative of a function you can only evaluate: a simulation, a table with interpolation, a model
too long to differentiate by hand. `nxx::deriv` estimates derivatives of order 1 to 4 of such a function from its values
at a few points near $x$. It offers *finite differences*, fixed formulas that cost a known number of evaluations, and
*Ridders' method*, which repeats a finite difference with shrinking steps and extrapolates to a zero step, returning an
estimate of its own error as well.

A derivative is computed the same way as a root, `derivative<Method>(fn, x)`, and returns a `std::expected` that holds
the estimate and, for Ridders' method, an estimate of its error: a `std::optional`, empty for the finite differences
([Working with Numerixx](../basics.md#what-a-result-holds) says what they mean). `derivative_of(fn)` turns `fn` into its
derivative as a function, which [Give derivatives](../roots/derivatives.md) uses to give Newton's method a numerical
$f'$.

The ideas, the choice of method and step, and the accuracy you can expect are explained at length in [the book](book.md)
on `nxx::deriv`. This chapter of the documentation will gain an overview, a tutorial and guide pages like those of
[Roots of functions](../roots/index.md).

```{toctree}
:maxdepth: 1

book
reference
```

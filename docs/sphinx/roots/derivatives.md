# Give derivatives

Newton's method replaces $f$ by its tangent line at the current point, and steps to where that line crosses zero. To
draw the tangent it needs the derivative $f'$. Halley's method fits a curve that also matches the curvature, so it needs
the second derivative $f''$ as well. In return, both close in on a root in very few steps once they are near it. This
page shows the three ways to give them the derivatives: from `fn` itself, as separate functions, or computed
numerically.

The examples use $f(x) = x^3 - 2x - 5$, with $f'(x) = 3x^2 - 2$ and $f''(x) = 6x$.

## From `fn` itself

Return the value and the derivatives together, as a `std::pair` (the value and $f'$) or a `std::tuple` (the value, $f'$
and $f''$); a `std::array` works too. One call then gives everything at a point, which suits a function whose value and
derivatives share work, such as a common subexpression or a table lookup:

```{literalinclude} ../examples/roots_derivatives.cpp
:start-after: // <pair>
:end-before: // </pair>
:dedent:
```

For Halley's method, return a tuple with $f''$ as well:

```{literalinclude} ../examples/roots_derivatives.cpp
:start-after: // <tuple>
:end-before: // </tuple>
:dedent:
```

Newton's method takes 7 evaluations in the bracket $[2, 3]$ (an interval at whose ends $f$ has opposite signs) and 6
from the guess 2, and Halley's method 5 from 2; one evaluation gives the value and the derivatives.

Every element of the pair or tuple must be exactly the real type. A bare literal, as in `std::pair{f, 1}`, makes the
second element an `int`, and the call does not compile. The rule is strict on purpose: in a multiprecision type an
element that only converts to the type could be an expression that refers to `fn`'s locals, gone once `fn` returns
([Working with Numerixx](../basics.md#multiprecision-types)).

## As separate functions

If you already have $f'$ as a function of its own, pass it after `fn`:

```{literalinclude} ../examples/roots_derivatives.cpp
:start-after: // <separate>
:end-before: // </separate>
:dedent:
```

Each function is then called once per point, and the result is the same as with the pair. For Halley's method, pass
$f''$ after $f'$.

## Numerically

When $f'$ has no formula you can write, {cpp:func}`nxx::deriv::derivative_of` turns `fn` into its numerical derivative,
a function you pass in the place of $f'$. It needs `numerixx/deriv.hpp`:

```{literalinclude} ../examples/roots_derivatives.cpp
:start-after: // <numerical>
:end-before: // </numerical>
:dedent:
```

This costs calls of `fn`. The default numerical derivative is a *central difference*: it evaluates `fn` a small step to
each side of the point, so each point of Newton's method costs 3 calls of `fn`, the value and the two sides. On
$x^3 - 2x - 5$, Newton's method with `derivative_of` makes 21 calls in $[2, 3]$ and 18 from the guess 2, where the
secant method, which needs no derivative at all, makes 8 from 2. For a function of values only, Brent's method or the
secant method is usually cheaper; the numerical derivative is a convenience.

For Halley's method, give $f''$ as `nxx::deriv::derivative_of<nxx::deriv::Central<2>, 2>(fn)` after $f'$. A point then
costs 6 calls of `fn`: the value, 2 for $f'$, and 3 for $f''$, which evaluates `fn` at the point again. That makes 30
from the guess 2 here.

```{note}
`derivative_of` keeps its own copy of `fn`. A counter that `fn` captures by reference, as in [Choose a
method](choose-a-method.md#measure-the-cost-on-your-own-function), is shared by the copy and counts every call. State
held inside a `mutable` lambda is not: pass `std::ref(fn)` to share it.
```

Whichever way you give the derivatives, the bound on the result does not depend on them: `abs_error` bounds the distance
to a sign change of `fn`, which the method verifies with `fn`'s values alone. A poor derivative can cost evaluations, or
make a solve from a guess fail, but it cannot make a root it returns less accurate. Near the edge of `fn`'s domain a
numerical derivative can fail; see [Solve on a restricted
domain](restricted-domain.md#a-numerical-derivative-near-the-edge).

The program needs:

```{literalinclude} ../examples/roots_derivatives.cpp
:start-after: // <includes>
:end-before: // </includes>
```

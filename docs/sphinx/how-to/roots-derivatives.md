# Give derivatives

{cpp:struct}`nxx::roots::Newton` needs f', and {cpp:struct}`nxx::roots::Halley` f' and f''. There are three ways to
give them.

## From fn itself

Return the value and the derivatives together, as a `std::pair` (value, f') or a `std::tuple` (value, f', f''); a
`std::array` works too. Every element must be exactly the real type: `std::pair{f, 1}` makes the second an `int`.
One call then gives everything at a point, which suits a fn whose value and derivatives share work:

```{literalinclude} ../examples/roots_derivatives.cpp
:start-after: // <pair>
:end-before: // </pair>
:dedent:
```

```{literalinclude} ../examples/roots_derivatives.cpp
:start-after: // <tuple>
:end-before: // </tuple>
:dedent:
```

On x^3 - 2x - 5, Newton takes 7 evaluations in the bracket [2, 3] and 6 from the guess 2, and Halley 5 from 2, where
one evaluation gives the value and the derivatives.

## As separate callables

Pass the derivative after fn. Each callable is evaluated once per point, and the result is the same as with the pair:

```{literalinclude} ../examples/roots_derivatives.cpp
:start-after: // <separate>
:end-before: // </separate>
:dedent:
```

## Numerically

When f' has no formula, {cpp:func}`nxx::deriv::derivative_of` turns fn into its numerical derivative (include
`numerixx/deriv.hpp`):

```{literalinclude} ../examples/roots_derivatives.cpp
:start-after: // <numerical>
:end-before: // </numerical>
:dedent:
```

This costs calls of fn: with the default central difference, each point evaluates fn there and at two nodes beside it, 3
calls in all. For Halley, give f'' as `nxx::deriv::derivative_of<nxx::deriv::Central<2>, 2>(f)` after f'; a point then
costs 6 calls (30 from the guess 2 on x^3 - 2x - 5). On x^3 - 2x - 5, Newton with `derivative_of` makes 21 calls in [2,
3] and 18 from 2, where {cpp:struct}`nxx::roots::Secant` makes 8 from 2 without any derivative. For a fn of value only,
Brent or Secant is usually cheaper; the numerical derivative is a convenience. The bound on the result never depends on
the derivative. Near the edge of fn's domain the derivative can fail: see [Solve on a restricted
domain](roots-restricted-domain.md).

The program needs:

```{literalinclude} ../examples/roots_derivatives.cpp
:start-after: // <includes>
:end-before: // </includes>
```

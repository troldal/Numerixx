# Solve on a restricted domain

Many functions are defined only on part of the line: a logarithm for x > 0, a physical model above some limit. Numerixx
treats a value that is not finite (NaN or infinite) as a point outside fn's domain. What a method does with such a
point depends on the method.

## Return NaN outside the domain

If fn's formula is finite outside its domain, where it means nothing, make it return NaN there. Otherwise a method
takes the meaningless value at its word. Van der Waals' equation of state is finite below the co-volume b, where no
gas exists:

```{literalinclude} ../examples/roots_restricted_domain.cpp
:start-after: // <nan>
:end-before: // </nan>
:dedent:
```

From 1e-4, {cpp:func}`nxx::roots::search` meets NaN on its lower side, approaches the edge of the domain instead of
stopping there, and finds the root at v = 4.3666e-4 m^3/mol after 29 evaluations.

## Prefer a bracket

A bracketing method never leaves its bracket: once the bracket lies inside the domain, every point it evaluates does.
An open method steps where its formula points and can leave the domain; the first value that is not finite ends the
solve with {cpp:enumerator}`nxx::roots::error::non_finite`:

```{literalinclude} ../examples/roots_restricted_domain.cpp
:start-after: // <leave>
:end-before: // </leave>
:dedent:
```

Newton's first step from 20 on ln x - 1 lands at -19.9, where the logarithm is NaN, so it fails; `search` and Brent from
the same guess find e. When an open method returns `non_finite`, search for a bracket instead.

## A numerical derivative near the edge

{cpp:func}`nxx::deriv::derivative_of` evaluates fn at points around x, a step away on each side. Near the edge of the
domain those points can fall outside it, and the derivative fails there. From a guess, Newton then ends with
{cpp:enumerator}`nxx::roots::error::not_converged`; in a bracket it takes a bisection step instead and still converges:

```{literalinclude} ../examples/roots_restricted_domain.cpp
:start-after: // <stencil>
:end-before: // </stencil>
:dedent:
```

The program needs:

```{literalinclude} ../examples/roots_restricted_domain.cpp
:start-after: // <includes>
:end-before: // </includes>
```

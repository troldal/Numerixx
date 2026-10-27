# Handle failures

Numerixx never throws. A solve that finds no root returns a `std::expected` holding an
{cpp:enum}`nxx::roots::error` instead of a {cpp:struct}`nxx::roots::Root`. Test for success, then decide what to try
next.

## Recover by trying something else

Each error suggests a next step:

- {cpp:enumerator}`nxx::roots::error::no_sign_change`: the bracket shows no sign change. Widen it, or search from a
  guess ([Find a bracket](roots-find-a-bracket.md)). From `search` itself, try another guess, or check that fn changes
  sign at all: an even root never does.
- {cpp:enumerator}`nxx::roots::error::not_converged`: an open method found no root from its guess. Give it a bracket.
- {cpp:enumerator}`nxx::roots::error::non_finite`: a value was not finite, often a point outside fn's domain ([Solve
  on a restricted domain](roots-restricted-domain.md)), or a start or tolerance was not finite, or a tolerance was
  negative.

`or_else` runs its fallback only on failure. Here Newton from 0 cycles on x^3 - 2x + 2 and gives `not_converged`;
the fallback searches for a bracket and solves in it:

```{literalinclude} ../examples/roots_handle_failures.cpp
:start-after: // <recover>
:end-before: // </recover>
:dedent:
```

## React to the error

More errors may be added, and a cause may move from one error to another in a later version. Test for success rather
than for a particular error, and give a `switch` on the error a `default`. An error formats as its name, so it can be
printed or logged directly. To reduce a result to a number, use `transform` and `value_or`:

```{literalinclude} ../examples/roots_handle_failures.cpp
:start-after: // <react>
:end-before: // </react>
:dedent:
```

Read a result with `->` or `*` after testing it. Its `value()` member throws `std::bad_expected_access` on a failure;
that is the standard library's exception, not Numerixx's. An exception thrown by fn itself passes through unchanged.

## Check for a pole

A pole, where fn jumps between infinities of opposite signs, changes sign like a root, and a bracketing method may
return it as one. The size of |fn(value)| alone does not tell them apart: a well-conditioned root can have a large
|fn| nearby, and a pole a small one. How |fn| changes near the result does: it grows away from a root and falls away
from a pole. Compare it a little to each side:

```{literalinclude} ../examples/roots_handle_failures.cpp
:start-after: // <pole-test>
:end-before: // </pole-test>
```

and apply it to a result:

```{literalinclude} ../examples/roots_handle_failures.cpp
:start-after: // <pole>
:end-before: // </pole>
:dedent:
```

A point that lands exactly on a pole, where fn is infinite, gives `non_finite` instead; which happens depends on the
method.

The program needs:

```{literalinclude} ../examples/roots_handle_failures.cpp
:start-after: // <includes>
:end-before: // </includes>
```

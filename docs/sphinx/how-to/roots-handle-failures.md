# Handle failures

Numerixx never throws. A solve that finds no root returns a `std::expected` holding an
{cpp:enum}`nxx::roots::error` instead of a {cpp:struct}`nxx::roots::Root`. Test for success, then decide what to try
next.

## Recover by trying something else

Each error suggests a next step:

- {cpp:enumerator}`nxx::roots::error::no_sign_change`: the bracket shows no sign change. Widen it, or search from a
  guess ([Find a bracket](roots-find-a-bracket.md)).
- {cpp:enumerator}`nxx::roots::error::not_converged`: an open method found no root from its guess. Give it a bracket.
- {cpp:enumerator}`nxx::roots::error::non_finite`: a value was not finite, often a point outside fn's domain ([Solve
  on a restricted domain](roots-restricted-domain.md)), or a start or tolerance was not finite.

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

A pole, where fn jumps from minus to plus infinity, changes sign like a root, and a bracketing method may return it as
one. Where fn may have poles, check |fn(value)|: small at a root, large at a pole.

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

# Handle failures

A solve can fail: the *bracket*, the interval at whose ends `fn` should have opposite signs, may hold no sign change;
Newton's method may wander off; `fn` may return NaN. A solve can also succeed and still not give you what you wanted,
when the sign change it found is a pole. This page shows how to tell what went wrong and what to try next, and how to
check a result for a pole.

A failed solve returns a `std::expected` that holds an {cpp:enum}`nxx::roots::error` instead of a
{cpp:struct}`nxx::roots::Root`. Nothing is thrown. [Working with Numerixx](../basics.md#results) explains how to test a
`std::expected` and chain steps on it; this page is about what the errors of `nxx::roots` mean.

## What each error means

Each error points to a next step:

- {cpp:enumerator}`nxx::roots::error::no_sign_change`: `fn` has the same sign at both ends of the bracket. Widen the
  bracket, or search for one from a guess ([Find a bracket](find-a-bracket.md)). When `search` itself returns it, try
  another guess, or check that `fn` changes sign at all: a root where `fn` touches zero without crossing it, as $x^2$
  does at 0, never does.
- {cpp:enumerator}`nxx::roots::error::not_converged`: an *open method*, one that starts from a guess rather than a
  bracket, found no root. Its steps cycled or diverged; or it had no usable step, because the derivative was 0, a
  numerical derivative failed near the edge of `fn`'s domain ([Solve on a restricted
  domain](restricted-domain.md#a-numerical-derivative-near-the-edge)), or the secant method met two equal values; or it
  ran out of steps, as on a multiple root at 0 with the default tolerances ([Set
  tolerances](tolerances.md#a-root-near-zero-set-abs_tol)). Give it a bracket instead.
- {cpp:enumerator}`nxx::roots::error::non_finite`: a value was not finite. Most often `fn` was evaluated outside its
  domain ([Solve on a restricted domain](restricted-domain.md)); it can also be a start or a tolerance that is not
  finite, or a tolerance that is negative.

## Recover by trying something else

`or_else` runs a fallback only when the step before it failed. Here Newton's method from the guess 0 cycles on
$x^3 - 2x + 2$: its steps go from 0 to 1 and back to 0, for ever, so it gives `not_converged`. The fallback searches for
a bracket from the same guess and solves in it with Brent's method:

```{literalinclude} ../examples/roots_handle_failures.cpp
:start-after: // <recover>
:end-before: // </recover>
:dedent:
```

`fdf` returns a pair, the value and the derivative, for Newton's method. `search` and Brent's method use only its first
element, the value, so the same function serves all three. The fallback receives the error, and could choose its step by
it; this one ignores it. To end with a number, report any error and fall back on NaN:

```{literalinclude} ../examples/roots_handle_failures.cpp
:start-after: // <react>
:end-before: // </react>
:dedent:
```

```{literalinclude} ../examples/roots_handle_failures.cpp
:language: text
:start-after: react_output = R"(
:end-before: )";
```

The recovery worked, so the program prints only the root. The `switch` has a `default`, because a later version may add
errors.

## Check for a pole

A *pole* is a point where $f$ grows without bound, as $1/x$ does at 0. Where $f$ jumps between large values of opposite
signs, it changes sign there, just as at a root, and a bracketing method may return the pole as a root. The bound
`abs_error` is still true: the result is close to a sign change of `fn`. It is the wrong kind of sign change.

The size of $|f|$ at the result alone does not tell the two apart: a well-conditioned root can have a large $|f|$
nearby, and a pole a small one. How $|f|$ changes near the result usually does. Away from a root $|f|$ grows; away from
a pole it falls. So compare $|f|$ at the result with $|f|$ a little to each side:

```{literalinclude} ../examples/roots_handle_failures.cpp
:start-after: // <pole-test>
:end-before: // </pole-test>
```

and apply the test to a result. Here $1/(x - 1/3)$ changes sign at its pole, and Brent's method returns it:

```{literalinclude} ../examples/roots_handle_failures.cpp
:start-after: // <pole>
:end-before: // </pole>
:dedent:
```

```{literalinclude} ../examples/roots_handle_failures.cpp
:language: text
:start-after: pole_output = R"(
:end-before: )";
```

The test is a heuristic. It works where `fn`'s rounding is small next to the change in $|f|$ over that distance. Where
`fn` loses most of its digits to cancellation near the result, as the written-out cubic on [Roots of
functions](index.md#what-you-get-back) does, the three values are mostly rounding, and the test's verdict is not to be
trusted.

```{note}
A point that lands exactly on a pole, where `fn` is infinite, gives `non_finite` instead of a result; which happens
depends on the method and the bracket.
```

The program needs:

```{literalinclude} ../examples/roots_handle_failures.cpp
:start-after: // <includes>
:end-before: // </includes>
```

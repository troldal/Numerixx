# Handle failures

A derivative can fail outright, when `fn` returns NaN where a method needs a value. It can also return a number that is
not accurate enough; for the Ridders methods, the error estimate tells you. This section shows how to recognise each
case, and what to try next.

## When there is no derivative

A failed call returns a `std::expected` holding an {cpp:enum}`nxx::deriv::error` instead of an
{cpp:struct}`nxx::deriv::Estimate`. [Working with Numerixx](../basics.md#results) shows how to test the result, and why
to test for success rather than for a particular error. There is one,
{cpp:enumerator}`nxx::deriv::error::non_finite`: no finite estimate exists. Its causes:

- `fn` returned NaN or an infinity where the method needed a value. A stencil fails as soon as one of its points gives
  one. The Ridders methods first restart with smaller steps, and fail only when no step they try avoids the bad points.
  The usual cause is a point outside `fn`'s domain: use a one-sided method there ([Near the edge of a
  domain](one-sided.md)).
- `x` is not finite, or a step you gave is not positive and finite.
- The step is below the resolution of `x`: the method's points round to the same number, as they do for $h = 10^{-300}$
  at 1.
- An intermediate value would overflow.

Here the central difference's points reach below 0, where $\sqrt x$ is NaN:

```{literalinclude} ../examples/deriv_handle_failures.cpp
:start-after: // <failure>
:end-before: // </failure>
:dedent:
```

```{literalinclude} ../examples/deriv_handle_failures.cpp
:language: text
:start-after: failure_output = R"(
:end-before: )";
```

A wrong but finite value is not detected at all ([A formula that goes
on](one-sided.md#a-formula-that-goes-on)).

## A result that is not accurate enough

Before you use a result of the Ridders methods, compare its estimate, `abs_error`, with the accuracy you need: here, a
relative error of $10^{-8}$. The test also checks that there is an estimate, since the stencils leave it empty:

```{literalinclude} ../examples/deriv_handle_failures.cpp
:start-after: // <trust>
:end-before: // </trust>
:dedent:
```

```{literalinclude} ../examples/deriv_handle_failures.cpp
:language: text
:start-after: trust_output = R"(
:end-before: )";
```

The derivative of $\sqrt x$ at $10^{-11}$ is $1.581 \times 10^5$. The function changes on the scale of $10^{-11}$, far
below the step `ForwardRidders` starts from, $1/32$. Its result is 7% off, and its estimate says that it is nowhere
near $10^{-8}$. With a step on the function's own scale, $x/8$, the result is good to about 12 digits.

Closer still, at $10^{-12}$, the result is 70% off. Its estimate is then at least a quarter of the value, which is how
the Ridders methods flag a result they cannot vouch for. The result at $10^{-11}$ did not raise that flag, so test
against the accuracy you need.

```{note}
Where the derivative is near 0, compare `abs_error` with a scale of your own rather than with the value. And the
estimate assumes that `fn` is accurate to about the type's epsilon. A less accurate `fn`, such as a solve to a loose
tolerance or interpolated data, can make it too small: compare the results at two steps, or of two methods.
```

The program needs:

```{literalinclude} ../examples/deriv_handle_failures.cpp
:start-after: // <includes>
:end-before: // </includes>
```

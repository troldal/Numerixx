# Handle failures

A derivative can fail outright, when `fn` returns NaN where a method needs a value, or a step is unusable. It can also
return a number that is not accurate enough, and for Ridders' methods the error estimate is how you tell. This page
shows how to recognise each case, and what to try next.

The examples write `dv` for `nxx::deriv`.

## When there is no derivative

A failed call returns a `std::expected` holding an {cpp:enum}`nxx::deriv::error` instead of an
{cpp:struct}`nxx::deriv::Estimate`. Nothing is thrown. [Working with Numerixx](../basics.md#results) explains how to
test the result; this page is about what the errors mean. Today there is one,
{cpp:enumerator}`nxx::deriv::error::non_finite`: no finite estimate exists. Its causes:

- `fn` returned NaN or an infinity where the method needed a value. A stencil fails as soon as one of its points gives
  one. Ridders' methods first restart with smaller steps, and fail only when no step they try avoids the bad points. The
  usual cause is a point outside `fn`'s domain: use a one-sided method there ([Near the edge of a
  domain](one-sided.md)).
- `x` is not finite, or a step you gave is not positive and finite.
- The step is below the resolution of `x`: the method's points round to the same number, as they do for $h = 10^{-300}$
  at 1.
- An intermediate value would overflow: a point $x + kh$, a difference of values, or the result.

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

The error formats as its name. More errors may be added in a later version, so test for success rather than for
`non_finite`.

```{warning}
Only values that are not finite are detected. A function that returns a finite but wrong number, such as a formula
used outside its range, or that has a kink or a jump among the points, gives a wrong derivative without an error, and
possibly with a small estimate. [Near the edge of a domain](one-sided.md#a-formula-that-goes-on) shows a case.
```

## A result that is not accurate enough

Ridders' methods estimate their error, in `abs_error`. Before you use a result, compare that estimate with the accuracy
you need: here, a relative error of $10^{-8}$. The test also checks that there is an estimate at all, since the stencils
leave it empty:

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
below the steps `ForwardRidders` starts from ($1/32$), and it cannot come down far enough: its result is 7% off, and its
estimate, 11% of the value, says that it is nowhere near $10^{-8}$. With a step on the function's own scale, $x/8$, the
result is good to about 12 digits.

Closer still, at $10^{-12}$, the method gives up: its result is 70% off, and its estimate is at least a quarter of the
value, which is how Ridders' methods flag a result their own estimate cannot vouch for. Do not rely on that flag alone:
the result at $10^{-11}$, 7% off, did not raise it. Test against the accuracy you need.

```{note}
A relative test means nothing where the derivative is 0 or near it: there, compare `abs_error` with a scale of your own.
And the estimate assumes that `fn` is accurate to about the type's epsilon. A function computed with less accuracy, such
as one that solves an equation to a loose tolerance or reads interpolated data, can make it too small. For such a
function, compare the results at two steps, or with two methods, before trusting the digits.
```

The program needs:

```{literalinclude} ../examples/deriv_handle_failures.cpp
:start-after: // <includes>
:end-before: // </includes>
```

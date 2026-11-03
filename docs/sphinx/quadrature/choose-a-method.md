# Choose a method

The default method, adaptive Gauss-Kronrod, integrates most functions on a finite interval well, and it finds the
hard places by itself. This page says when another method does better, how to measure the cost on your own function, and
how to use an integral inside another algorithm, such as a root solve or a second integral. The examples write `qd` for
`nxx::quadrature`.

## From what you have to a method

| You have | Use |
|---|---|
| A function on a finite interval, smooth or with peaks, kinks or steps you have not located | `GaussKronrod<15>`, the default: `integrate(fn, a, b)` |
| A kink, a step or a peak at a point you know | The default on each side of the point ([Singularities, kinks and steps](rough-integrands.md)) |
| A singularity at an end, such as $\log x$ or $1/\sqrt x$ at 0 | `TanhSinh`, with the singular end moved to 0 if it is elsewhere ([Singularities, kinks and steps](rough-integrands.md)) |
| An infinite interval | `TanhSinh`, the one method that takes an infinite end ([Infinite intervals](infinite-intervals.md)) |
| A multiprecision type, and many digits | `TanhSinh` ([High precision](high-precision.md)) |
| An integral inside another algorithm, whose result must change smoothly with a parameter of `fn` | `GaussLegendre<N>`, or the default with a tight tolerance ([below](#an-integral-inside-another-algorithm)) |
| A function so costly that 150 evaluations are too many | `GaussLegendre<N>`, which costs exactly $N$, with no error estimate |
| A wish to watch the steps of an extrapolation | `Romberg` |

A *singularity* here is a point where the integrand grows without bound while its integral stays finite, as
$1/\sqrt x$ does at 0. A *kink* is a jump in the slope, a *step* a jump in the value.

## Measure the cost on your own function

Each method's cost is the number of times it calls `fn`, and the way to measure it is to count the calls in `fn`
itself, through a counter that the lambda captures by reference. Here the default and `TanhSinh` integrate $1/\sqrt x$
over $[0, 1]$, whose integral is 2, with the default request, a relative error of about $1.5 \times 10^{-8}$:

```{literalinclude} ../examples/quadrature_choose_a_method.cpp
:start-after: // <cost>
:end-before: // </cost>
:dedent:
```

```{literalinclude} ../examples/quadrature_choose_a_method.cpp
:language: text
:start-after: cost_output = R"(
:end-before: )";
```

Both meet the request. The default cuts the interval into *panels*, pieces it integrates one at a time, and here it cuts
the panel next to 0 in half, again and again, 30 evaluations a cut; `TanhSinh` places its points ever closer to 0 from
the start, and gets the integral to the last digit with ten times fewer evaluations. On a smooth function the two cost
about the same ([Costs](costs.md)). Output quoted on this site assumes that the compiler does not fuse a multiplication
and an addition into one instruction ([Working with Numerixx](../basics.md#compile-time)); with fused arithmetic the
last digits and the counts can differ.

`integrate<TanhSinh>` also has a cost at compile time: its node tables are built in every translation unit that uses it,
which makes those files noticeably slower to compile. If several files integrate with it, put the calls in one of them.

## An integral inside another algorithm

An integral is often one step of a larger computation: a root solve for the upper limit at which an integral reaches a
given value, an optimiser that tunes a parameter inside `fn`, or a double integral computed as an integral of
integrals. Three things then matter.

**The outer algorithm sees only what the integral returns.** Suppose you want the $x$ at which
$F(x) = \int_0^x \sqrt{t}\,dt$ equals $2/3$; the answer is 1. Newton's method in `nxx::roots` needs $F$ and its
derivative, and the derivative of an integral with respect to its upper limit is the integrand itself,
$F'(x) = \sqrt x$, exactly. So `fdf` returns both, as a pair:

```{literalinclude} ../examples/quadrature_choose_a_method.cpp
:start-after: // <limit>
:end-before: // </limit>
:dedent:
```

```{literalinclude} ../examples/quadrature_choose_a_method.cpp
:language: text
:start-after: limit_output = R"(
:end-before: )";
```

With the default request the integral is off by $2.9 \times 10^{-10}$, and so is the root, whose `abs_error` of
$4.4 \times 10^{-16}$ says nothing about it: the root solver bounds its distance to a sign change of $F$ *as computed*,
and cannot see the integral's error. Its own default asks for about two units in the last place, far more than the
integral's default delivers. Either tighten the integral's request to match, as the second line does, or loosen the
root's ([Set tolerances](../roots/tolerances.md) in the roots chapter), and read the result as no more accurate than
the integral. Each solve here took six integrals.

**A function must return a number.** `fn`, here `F`, cannot return a failure, so it maps one to NaN. `nxx::roots`
treats a NaN as a point outside the domain. `nxx::quadrature` reports
{cpp:enumerator}`nxx::quadrature::error::non_finite` where inner failures cover a stretch of the outer interval, but the
default method steps around a single failed point and can return a number that silently leaves it out ([Handle
failures](handle-failures.md)). If an inner failure must not pass, count the failures in `fn` and check the count after
the outer call.

**The result should change smoothly with the parameter.** An adaptive method chooses its points from `fn`, so as a
parameter of `fn` changes, the panels it cuts can change too, and its result moves in small jumps rather than smoothly
(Lyness 1983). The default is usually so far inside its request that the jumps lie near the rounding level, but nothing
guarantees that, and an algorithm that takes differences of the result at small steps, such as a secant method, a
numerical derivative or an optimiser, can stall on them. `GaussLegendre<N>` evaluates `fn` at the same $N$ points
whatever the parameter, so its result is as smooth as `fn`. It has no error estimate: compare two values of $N$, or
compare it once with the default, to choose $N$.

A double integral is an integral whose integrand is itself an integral. Give the inner one a tighter tolerance than
the outer one: the outer method sees the inner errors as noise in its integrand, and its estimate assumes an integrand
accurate to about the precision of the type. Here is $\int_0^1 \int_0^1 dx\,dy / (1 + x + y) = 3 \ln 3 - 4 \ln 2$:

```{literalinclude} ../examples/quadrature_choose_a_method.cpp
:start-after: // <nested>
:end-before: // </nested>
:dedent:
```

```{literalinclude} ../examples/quadrature_choose_a_method.cpp
:language: text
:start-after: nested_output = R"(
:end-before: )";
```

The cost multiplies: 150 inner integrals of 150 evaluations each. Count it before you nest a third level.

```{note}
This page uses `nxx::roots` beside `nxx::quadrature`, through the aliases `rt` and `qd`. Both modules have a `Criteria`
and an `error`, so do not bring both in with `using namespace` ([Working with
Numerixx](../basics.md#what-a-result-holds)).
```

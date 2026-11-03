# Handle failures

An integration can end without an integral: the input was invalid, the function returned something that is not a
number, or the method ran out of work before its estimate met your request. This page shows how to read a failure,
when its best estimate is worth keeping, how to tell a divergent integral, how to recover, and how to watch a method's
steps. The examples write `qd` for `nxx::quadrature`.

## What each error means

A method with an error estimate returns a `std::expected` whose error side is a {cpp:struct}`nxx::quadrature::Failure`.
Its `reason` is one of two values, and it formats as the value's name:

```{literalinclude} ../examples/quadrature_handle_failures.cpp
:start-after: // <reasons>
:end-before: // </reasons>
:dedent:
```

```{literalinclude} ../examples/quadrature_handle_failures.cpp
:language: text
:start-after: reasons_output = R"(
:end-before: )";
```

{cpp:enumerator}`nxx::quadrature::error::non_finite` means that a number was not finite where the method cannot accept
one. The input can cause it: an end that is NaN, an infinite end for any method but `TanhSinh`, or a tolerance that is
negative or not finite. So can the function: $\sqrt{x - 0.5}$ is NaN on half of $[0, 1]$, and Romberg's method evaluates
$f$ at the ends, where $1/x$ is infinite. So does a sum that overflows. A different tolerance will not help; change the
interval, the method, or the way $f$ is written.

Not every value that is not finite fails. The default method works on *panels*, equal pieces of the interval, and cuts a
panel in half, a *bisection*, where its estimate is largest. A value that is not finite at one of its points makes it
cut that panel and sample the two halves afresh, at points that usually miss the bad one, and go on. A region where $f$
is NaN cannot be cut away, and fails:

```{literalinclude} ../examples/quadrature_handle_failures.cpp
:start-after: // <isolated>
:end-before: // </isolated>
:dedent:
```

```{literalinclude} ../examples/quadrature_handle_failures.cpp
:language: text
:start-after: isolated_output = R"(
:end-before: )";
```

One extra bisection, 30 evaluations, and the result does not show the bad point at all; a point no node lands on is
never seen. So a success does not prove that $f$ is finite everywhere. `TanhSinh` accepts a value that is not finite
only very close to an end, where a singularity sits ([Singularities, kinks and steps](rough-integrands.md)).

{cpp:enumerator}`nxx::quadrature::error::not_converged` means that the method ran out of work before its estimate met
your request. Each method has a cap: 16 *levels* for Romberg and 10 for `TanhSinh`, where a level halves the spacing
of the points, and for the default a number of bisections that grows with the type's precision. The default also
gives up when 32 bisections in a row fail to improve its estimate, as on $1/x$. The integral may need another method,
a substitution, a split or a request your function can support; or it may not exist.

## A failure's best

A method that runs out of work can still return the last estimate it reached, as the failure's `best`, an
{cpp:struct}`nxx::quadrature::Integral`. Two failures show the range. The first integrates $1/(1 + x)$ with a relative
noise of $10^{-10}$ in its values (a hash of the bits of $x$, as an iterative solver stopped at that tolerance would
leave), asked for $10^{-12}$; the second $1/x$ over $[0, 1]$, which is infinite:

```{literalinclude} ../examples/quadrature_handle_failures.cpp
:start-after: // <best>
:end-before: // </best>
:dedent:
```

```{literalinclude} ../examples/quadrature_handle_failures.cpp
:language: text
:start-after: best_output = R"(
:end-before: )";
```

The noisy run's best happens to be within its `abs_error` of $\ln 2$. But `abs_error` cannot see noise in $f$, so what
you can count on is `abs_error` plus the noise's bound, $10^{-10}$ times the integral of $|f|$: about
$7 \times 10^{-11}$ in all. If that is close enough, the best serves. The second is a finite number for an integral
that does not exist. Treat a best as a lead, not a result: use it only where you know the integral exists.

## A divergent integral

No error says "this integral is infinite", and the signs are easy to misread. The default's best on $1/x$ is the same
at every request, but that is no sign of convergence: on a panel $[0, h]$ both of its rules give the same sum for every
$h$, so cutting panels near 0 changes nothing. Worse, asked for little, the default reports success:

```{literalinclude} ../examples/quadrature_handle_failures.cpp
:start-after: // <diverge>
:end-before: // </diverge>
:dedent:
```

```{literalinclude} ../examples/quadrature_handle_failures.cpp
:language: text
:start-after: diverge_output = R"(
:end-before: )";
```

At a relative tolerance of 0.2 its first 150 evaluations give $9.33 \pm 1.8$, which meets the request. The real signs
are two. `TanhSinh` returns no best at all on $1/x$, since its bound on what lies beyond its last points grows without
limit. And the integral changes when you move the singular end: over $[10^{-k}, 1]$ it grows like $k \ln 10$, without
settling, where a convergent integral would settle on a value.

## Recover

- **Split** the interval at a kink, a step or a peak whose position you know, and add the pieces
  ([Singularities, kinks and steps](rough-integrands.md)).
- **Substitute** a variable that removes a singularity, or moves it to 0 ([the tutorial](tutorial.md) does both).
- **Ask for less**, when $f$ itself is accurate to fewer digits than its type: the noisy function above succeeds at
  $10^{-10}$.
- **Another method.** `TanhSinh` suits a singularity at an end but not a kink inside, where the default cuts panels.
  Try `TanhSinh` first and fall back to the default, accepting the default's success only where it agrees with
  `TanhSinh`'s best within the two estimates, since at a strong singularity the default can succeed falsely.

Here is that fallback on $1/\sqrt x + |x - 1/3|$, singular at 0 and with a kink at $1/3$, whose integral is
$2 + 5/18$:

```{literalinclude} ../examples/quadrature_handle_failures.cpp
:start-after: // <fallback>
:end-before: // </fallback>
:dedent:
```

```{literalinclude} ../examples/quadrature_handle_failures.cpp
:language: text
:start-after: fallback_output = R"(
:end-before: )";
```

`TanhSinh` gave up at the kink, with a best; the default's result agreed with it, and is within its `abs_error` of the
exact value. `or_else` runs the fallback only on a failure ([Working with Numerixx](../basics.md#composing-steps)).

## Watch the steps

The default, `TanhSinh` and Romberg's method work in steps, bisections or levels. To see why a run failed, pass an
*observer*, a callable that receives each {cpp:struct}`nxx::quadrature::Step`, wrapped by
{cpp:func}`nxx::quadrature::observer`, before `fn`. A step holds the `current` estimate, if any, and flags: `ready` once
the method has done the minimum work for its estimate to be trusted, `exhausted` when it cannot go on. Count the steps
in the observer and the evaluations in `fn`:

```{literalinclude} ../examples/quadrature_handle_failures.cpp
:start-after: // <observe>
:end-before: // </observe>
:dedent:
```

```{literalinclude} ../examples/quadrature_handle_failures.cpp
:language: text
:start-after: observe_output = R"(
:end-before: )";
```

At no level did `TanhSinh` have an estimate, since its bound beyond its last point near 0 was unbounded; after the
tenth it was exhausted, without a best. On an integral that converges, the estimate falls level by level instead. The
observer changes nothing in the run. An exception thrown by `fn` or by an observer passes through `integrate`
unchanged ([Working with Numerixx](../basics.md#results)), so throwing from an observer can stop a run that takes too
long.

The output on this page is the same on every toolchain, as long as the compiler does not fuse multiplications and
additions ([Working with Numerixx](../basics.md#compile-time)).

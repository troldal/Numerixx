# Infinite intervals

Some integrals run to infinity: the probability in the tail of a distribution, a Laplace transform, the energy radiated
over all wavelengths. You cannot place a finite number of points on an infinite interval and hope to cover it, so a
method must first map the interval onto one it can handle. `TanhSinh` does this, and it is the only method in
`nxx::quadrature` that takes an infinite end. The examples write `qd` for `nxx::quadrature`, and the output they quote
uses only arithmetic and comparisons, so it is the same on every platform built without fused multiply-add ([Working
with Numerixx](../basics.md#compile-time)).

## An infinite end

Write an infinite end as `std::numeric_limits<double>::infinity()`, or its negative. On the whole line, $1 / (1 + x^2)$
integrates to $\pi$, and on $[0, \infty)$ to $\pi / 2$:

```{literalinclude} ../examples/quadrature_infinite_intervals.cpp
:start-after: // <line>
:end-before: // </line>
:dedent:
```

```{literalinclude} ../examples/quadrature_infinite_intervals.cpp
:language: text
:start-after: line_output = R"(
:end-before: )";
```

On a half line, `TanhSinh` substitutes $x = a + e^{u}$ with $u = \tfrac32 \sinh t$ (the *exp-sinh* map), and on the
whole line $x = \sinh u$ (the *sinh-sinh* map). Equal steps in $t$ then become points that reach out towards infinity
faster than exponentially, and the decay of $f$ there makes the far terms negligible. The method still never evaluates
$f$ at infinity itself.

Every other method returns {cpp:enumerator}`nxx::quadrature::error::non_finite` for an infinite end. The call compiles
and fails when it runs:

```{literalinclude} ../examples/quadrature_infinite_intervals.cpp
:start-after: // <others>
:end-before: // </others>
:dedent:
```

```{literalinclude} ../examples/quadrature_infinite_intervals.cpp
:language: text
:start-after: others_output = R"(
:end-before: )";
```

## A finite part and an infinite tail

`TanhSinh` is weak on a kink or a step inside the interval ([Singularities, kinks and steps](rough-integrands.md)), and
an infinite interval is no exception. If the rough spot is at a known place, split there: give the finite piece to the
default method and the infinite tail to `TanhSinh`. Here the integrand doubles at $x = 1$:

```{literalinclude} ../examples/quadrature_infinite_intervals.cpp
:start-after: // <split>
:end-before: // </split>
:dedent:
```

```{literalinclude} ../examples/quadrature_infinite_intervals.cpp
:language: text
:start-after: split_output = R"(
:end-before: )";
```

The split gives $3\pi/4$ to the last digit or two.

## Where the maps look

Both maps are centred: the line's at 0, a half line's at its finite end, and both at a scale of about 1. They place
most of their points within a few units of the centre and spread the rest ever more thinly outwards. An integrand
whose mass lies far from the centre is seen only by a few widely spaced points. A peak at $x_0 = 10$ on the line still
works, at about 2,000 evaluations; at $x_0 = 40$ it no longer does.

```{eval-rst}
.. plot:: plots/quadrature_line.py
   :alt: Vertical lines at the points TanhSinh places on the line in its first five levels, dense around 0 and ever
         sparser outwards, and the peak 1 / (1 + (x - x0)^2) at x0 = 0, 10 and 40. The peak at 0 is covered by many
         lines, the peak at 10 by one, and the peak at 40 falls between two lines about 11 apart.
   :caption: The points ``TanhSinh`` places on :math:`(-\infty, \infty)` in its first five levels, which every run
             takes, and the peak :math:`1 / (1 + (x - x_0)^2)` at three centres.
```

The figure shows where the runs start. Each *level* of `TanhSinh` halves the spacing of its points, and every run
takes the first five: here the 133 points of the run on $1 / (1 + x^2)$, 61 of them within the figure's range, where
the map puts them in the same places whatever the integrand. Within the half-width of the peak at 0, where it is at
least half its height, there are 17, and the run on it ends there. At 10 there is one, the points about 2 apart; four
more levels bring them 0.13 apart, with 15 within the half-width, and the run succeeds after 2,069 evaluations. At 40
there is none: the two nearest are about 11 apart, five times the peak's width. The tenth and last level still leaves
them 0.36 apart, with 5 within the half-width, and the run ends `not_converged` after 4,119 evaluations.

Split at the peak, or substitute $u = x - x_0$ and write the integrand in $u$:

```{literalinclude} ../examples/quadrature_infinite_intervals.cpp
:start-after: // <far_peak>
:end-before: // </far_peak>
:dedent:
```

```{literalinclude} ../examples/quadrature_infinite_intervals.cpp
:language: text
:start-after: far_peak_output = R"(
:end-before: )";
```

A half line far from 0 has the opposite problem. Its points near the finite end $a$ are $a + e^{-u}$, which rounds to
the `double` numbers near $a$, spaced about $a \times 2 \times 10^{-16}$ apart; for $a = 10^{17}$ even the map's centre,
$a + 1$, rounds to $a$, and the method stops before it evaluates $f$ at all. The result degrades with $a$:

```{literalinclude} ../examples/quadrature_infinite_intervals.cpp
:start-after: // <far_end>
:end-before: // </far_end>
:dedent:
```

```{literalinclude} ../examples/quadrature_infinite_intervals.cpp
:language: text
:start-after: far_end_output = R"(
:end-before: )";
```

Both are reported as successes, but with an `abs_error` far above the default request. Compare `abs_error` with what you
need before you use the value, and substitute $x = a + u$, with the integrand written in $u$, to move the end to 0.

## Writing the integrand for large x

`TanhSinh` calls the integrand at very large $x$, beyond $10^{100}$ in a `double`. A value that is not finite there is
accepted as a negligible tail only far out, or once the terms are already negligible; nearer, it is an error. A product
of a growing and a decaying factor can produce such a value:

```{literalinclude} ../examples/quadrature_infinite_intervals.cpp
:start-after: // <overflow>
:end-before: // </overflow>
:dedent:
```

$e^{-x}$ underflows to exactly 0 from $x \approx 745$, and $x^{51}$ overflows to infinity from about
$1.1 \times 10^{6}$; beyond that the product is $\infty \times 0$, which is NaN. The method's first points jump from
$x \approx 230$, where the integrand is still $2 \times 10^{20}$, straight to $x \approx 3.4 \times 10^{6}$, so it meets
the NaN before it has seen the terms become negligible, and the first form returns `non_finite`. With $x^{50}$ the same
NaN is accepted, because the integrand before it is a hundred times smaller; do not read that as a rule about the
power, since it depends on where the points fall. The second form computes the same function through logarithms,
where nothing overflows, and gives $51!$ to within its estimate. Write such an integrand so that it stays finite
wherever its value does.

A value that underflows to 0 is harmless, as long as the terms are negligible by then: $e^{-x^2/2}$ over
$[3, \infty)$, the tail of a normal distribution, is 0 in a `double` beyond $x \approx 38.6$, and integrates to within
its estimate. What fails is an integrand tiny everywhere, such as $10^{-300} x^{-1.5}$ over $[1, \infty)$: values below
the smallest normal `double` are not trusted towards infinity, and the result is `not_converged`. Scale such an
integrand first.

## What does not work

Three kinds of integrand fail on an infinite interval, each with
{cpp:enumerator}`nxx::quadrature::error::not_converged`:

- **Slow decay.** $x^{-1.01}$ over $[1, \infty)$ integrates to 100, but its tail decays so slowly that a `double` runs
  out of range before the terms become negligible: the best estimate is about 99.9, with an error estimate that says
  so. More digits do not help, a wider exponent range does: the tail beyond the largest `double` still holds about
  0.08 of the integral, while in `cpp_bin_float_50`, whose exponents reach far further, the same integral succeeds.
- **A large number written for infinity.** $[0, 10^{300}]$ is a finite interval to `TanhSinh`, mapped by the
  tanh-sinh rule onto its own scale. $1 / (1 + x^2)$ has its mass within a few units of 0, a relative width of about
  $10^{-299}$ at one end of that interval, which the points cannot resolve: `not_converged`. Written as `infinity()`,
  the same integral succeeds. A moderate number such as $10^{10}$ works, but integrates exactly what it
  says: $1 / (1 + x^2)$ over $[0, 10^{10}]$ is $\pi/2 - 10^{-10}$.
- **Oscillation.** $\sin x / x$ over $[0, \infty)$ converges only because its positive and negative lobes cancel, and no
  method here sums lobes. Integrate it lobe by lobe with the default method and sum the series, or use a method
  made for oscillating integrands.

[Handle failures](handle-failures.md) says what a failure's best estimate is worth.

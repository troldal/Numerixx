# The period of a pendulum

A pendulum's period is almost independent of how far it swings: that is what made pendulum clocks possible. Almost,
because the textbook formula $T_0 = 2\pi\sqrt{L/g}$ holds only for small swings. At larger amplitudes the period grows,
and its exact value is an integral that has no closed form. This tutorial computes it step by step. You will write the
integrand as a function, integrate it with the default method, check the result against an independent value, see how
an end of the interval that is not exactly a `double` costs accuracy, rewrite the integral with a singularity at an
end, integrate that with the tanh-sinh method, and see what happens when the integral does not exist. [Integrals of
functions](index.md) explains the ideas behind each step.

## The integral

A pendulum of length $L$ released from rest at an angle $\theta_0$ from the vertical swings with the period

$$T = 4\sqrt{L/g}\; K(k), \qquad k = \sin(\theta_0 / 2)$$

where $K$ is the *complete elliptic integral of the first kind*:

$$K(k) = \int_0^{\pi/2} \frac{d\varphi}{\sqrt{1 - k^2 \sin^2\varphi}}$$

For small swings $k$ is near 0, $K$ is near $\pi/2$, and $T$ is near $T_0$. The ratio $T / T_0 = 2K / \pi$ says how much
longer the true period is.

It pays to write the integrand with the *complementary modulus* $k' = \cos(\theta_0 / 2)$, so that $k^2 + k'^2 = 1$:

$$K = \int_0^{\pi/2} \frac{d\varphi}{\sqrt{k'^2 + k^2 \cos^2\varphi}}$$

The two forms are equal, but they round differently. Near $\theta_0 = 180^\circ$, $1 - k^2 \sin^2\varphi$ subtracts two
nearly equal numbers, while $k'^2 + k^2\cos^2\varphi$ adds two positive ones. In code, $K$ becomes a function of $k'$
that builds the integrand as a lambda and integrates it over $[0, \pi/2]$ with the default method. The examples
write `qd` for `nxx::quadrature`:

```{literalinclude} ../examples/quadrature_pendulum.cpp
:start-after: // <smooth>
:end-before: // </smooth>
:dedent:
```

The program needs these headers:

```{literalinclude} ../examples/quadrature_pendulum.cpp
:start-after: // <includes>
:end-before: // </includes>
```

## The integral at 90 degrees

For a swing of $90^\circ$, $k' = \cos 45^\circ = \sqrt{1/2}$:

```{literalinclude} ../examples/quadrature_pendulum.cpp
:start-after: // <first>
:end-before: // </first>
:dedent:
```

{cpp:func}`nxx::quadrature::integrate` returns a `std::expected`: either an {cpp:struct}`nxx::quadrature::Integral`,
or the reason there is none. The `if (K)` tests which, and only then does the code read the result. `value` is the
integral; `abs_error` estimates its distance from the exact value. The program prints $K \approx 1.854074677301372$
with an `abs_error` of about $2 \times 10^{-14}$. The default request is a relative error of $1.5 \times 10^{-8}$; the
method met it with its first 150 evaluations, and its estimate fell to the *rounding floor*, the level below which the
rounding of the values and of the sum would dominate any estimate.

## Checking the result

An estimate is not a proof, so check the result against a value computed another way. For $K$ there is a classic one:
$K = \pi / (2\,\mathrm{AGM}(1, k'))$, where the *arithmetic-geometric mean* $\mathrm{AGM}(a, b)$ is the common limit of
replacing $a$ and $b$ by their arithmetic mean $(a + b) / 2$ and their geometric mean $\sqrt{ab}$, over and over. The
two sequences meet within a handful of steps, to within a few rounding errors:

```{literalinclude} ../examples/quadrature_pendulum.cpp
:start-after: // <reference>
:end-before: // </reference>
:dedent:
```

The loop stops when $a$ and $b$ agree to a few units in the last place rather than when they are equal, since rounding
can leave them alternating between two neighbouring numbers forever. Its error is below $2 \times 10^{-15}$ relative
for every $k'$ in this tutorial. Compare:

```{literalinclude} ../examples/quadrature_pendulum.cpp
:start-after: // <check>
:end-before: // </check>
:dedent:
```

The difference is a few units in the last place at most, far below `abs_error`. The estimate overstates the error, as
it usually does.

Now the period for three amplitudes, written exactly: $k' = \cos 30^\circ$, $\cos 45^\circ$ and $\cos 60^\circ$:

```{literalinclude} ../examples/quadrature_pendulum.cpp
:start-after: // <periods>
:end-before: // </periods>
:dedent:
```

```{literalinclude} ../examples/quadrature_pendulum.cpp
:language: text
:start-after: periods_output = R"(
:end-before: )";
```

At $60^\circ$ the pendulum is 7% slower than the small-angle formula says, at $120^\circ$ 37%.

## Close to the top

As $\theta_0$ approaches $180^\circ$, $k'$ approaches 0, and $K$ grows like $\ln(4/k')$: a pendulum released almost
upside down hangs near the top for a long time. The integrand then has a sharp peak at $\varphi = \pi/2$, where
$\cos\varphi$ is near 0 and the integrand near $1/k'$. The default method copes with the peak, by cutting more
panels around it, and still meets the request, but it loses digits it would otherwise keep: at $k' = 10^{-6}$ the
result is off by more than $10^{-11}$, at $k' = 10^{-8}$ by more than $10^{-9}$.

The cause is that the peak sits where the `double` numbers are coarse. $\pi/2$ is not a `double`:
`std::numbers::pi / 2.0` lies about $6 \times 10^{-17}$ below it, so the interval misses a sliver in which the integrand
is about $1/k'$, worth about $6 \times 10^{-17} / k'$ of the integral: $6 \times 10^{-11}$ at $k' = 10^{-6}$, most of
the error there. And every point the method places near the end is rounded to a neighbouring `double`, by up to
$10^{-16}$, where the integrand changes fast; that costs about as much again, in either direction. Both grow like
$1/k'$. They also explain something stranger. Asked for a relative error of $10^{-14}$ at $k' = 10^{-4}$, the method
reports an `abs_error` smaller than its distance from $K$. Its estimate is honest about the integral over the interval
it was given, and that interval is not $[0, \pi/2]$.

The remedy is to put the peak at an end that is exact. Measuring the angle from the top, $\psi = \pi/2 - \varphi$, turns
$\cos\varphi$ into $\sin\psi$ and moves the peak to $\psi = 0$; then the result is within $10^{-13}$ of $K$ at both
$k'$. There is another way, which shows what the second method is for. Substitute $t = \sin\varphi$, which gives

$$K = \int_0^1 \frac{dt}{\sqrt{(1 - t^2)(1 - k^2 t^2)}}$$

This integrand has a *singularity* at $t = 1$: it grows without bound there, like $1/\sqrt{1 - t}$, while its integral
stays finite. A singularity at an end is what the tanh-sinh method, `TanhSinh`, is made for: its points crowd towards
the ends of the interval without ever reaching them. But they can only crowd as close to 1 as the `double` numbers
there allow, $10^{-16}$ apart, and written in $t$ the result is about $10^{-8}$ off at $60^\circ$: the same lesson as
$\pi/2$, from the other side. Move the singularity to 0, where the numbers are far denser, by substituting
$u = 1 - t$, and write $1 - t^2 = u(2 - u)$ in $u$ itself, so that nothing is computed as a difference near the
singular end. The integrand then needs no $\pi$ and no cosine at all:

```{literalinclude} ../examples/quadrature_pendulum.cpp
:start-after: // <algebraic>
:end-before: // </algebraic>
:dedent:
```

The function counts its evaluations through a reference, so that the caller can see the cost. Close to the top:

```{literalinclude} ../examples/quadrature_pendulum.cpp
:start-after: // <near>
:end-before: // </near>
:dedent:
```

```{literalinclude} ../examples/quadrature_pendulum.cpp
:language: text
:start-after: near_output = R"(
:end-before: )";
```

Every value agrees with the arithmetic-geometric mean to within a few units in the last place, at a few hundred
evaluations. The default method on this form would also succeed, but at four to eleven times the cost: it handles the
singularity by cutting panels, while `TanhSinh` was built for it. [Singularities, kinks and steps](rough-integrands.md)
says more about such integrands. These values and counts, like every output this manual quotes, are the same on every
platform whose compiler does not fuse a multiplication and an addition into one instruction ([Working with
Numerixx](../basics.md#compile-time)).

## At the top

At $\theta_0 = 180^\circ$ exactly, $k' = 0$: the pendulum is balanced upside down and never comes back. The period is
infinite, and so is the integral, since the integrand becomes $1 / (u(2 - u))$, which grows too fast at 0:

```{literalinclude} ../examples/quadrature_pendulum.cpp
:start-after: // <divergent>
:end-before: // </divergent>
:dedent:
```

```{literalinclude} ../examples/quadrature_pendulum.cpp
:language: text
:start-after: divergent_output = R"(
:end-before: )";
```

`TanhSinh` ran through all its levels, 5,159 evaluations, and returned
{cpp:enumerator}`nxx::quadrature::error::not_converged`. It also returned no `best` estimate: its bound on what lies
beyond its last points grew without limit, so it had nothing it could call an estimate. There is no separate error for
"this integral diverges", and other methods can return a best, or even a success at a loose request, for an integral
that does not exist: [Handle failures](handle-failures.md) shows how to tell.

## Next

[Choose a method](choose-a-method.md) says which method suits which integrand, [Set tolerances](tolerances.md) how to
ask for more digits, and [Singularities, kinks and steps](rough-integrands.md) how to treat integrands like this one.

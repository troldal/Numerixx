# Singularities, kinks and steps

A smooth integrand is easy: a Gauss rule fits it with a polynomial, and a handful of points per panel give every digit.
Real integrands are often rougher. They have a *singularity*, a point where they grow without bound while their
integral stays finite, such as $1/\sqrt{x}$ at 0; a *kink*, where the slope jumps, such as $|x - c|$; a *step*, where
the value jumps; or a narrow peak. Where the rough spot sits decides what to do: at an end, inside at a place you know,
or somewhere you do not know. The examples write `qd` for `nxx::quadrature`, and the output they quote uses only
arithmetic, square roots and comparisons, so it is the same on every platform built without fused multiply-add
([Working with Numerixx](../basics.md#compile-time)).

## A singularity at an end

The default method, `GaussKronrod<15>`, handles a singularity at an end by cutting the panel next to it in half,
again and again, until the panels' estimates are small enough. `TanhSinh` is built for it instead: its points crowd
towards the ends faster than exponentially, without ever reaching them, so it never evaluates $f$ at the singularity
and needs no cutting. On $1/\sqrt{x}$ over $[0, 1]$, whose integral is 2:

```{literalinclude} ../examples/quadrature_rough_integrands.cpp
:start-after: // <endpoint>
:end-before: // </endpoint>
:dedent:
```

```{literalinclude} ../examples/quadrature_rough_integrands.cpp
:language: text
:start-after: endpoint_output = R"(
:end-before: )";
```

Both meet the default request, a relative error of $1.5 \times 10^{-8}$, but `TanhSinh` gives every digit at a tenth
of the cost. On $\log x$ over $[0, 1]$ the pattern repeats: about 110 evaluations by `TanhSinh` against about 600.

```{eval-rst}
.. plot:: plots/quadrature_tanh_sinh_map.py
   :alt: Left: x against t from -4 to 4, an S-shaped curve from 0 to 1. Dots at steps of 1/2 in t, each carried
         across to the vertical x axis by a line, crowd towards x = 0 and x = 1. Right: on a logarithmic scale, the
         integrand in t rises from below 1e-10 at t = -4 to a maximum of about 0.87 near t = -0.75, and falls to about
         1e-22 at t = 4.
   :caption: The substitution :math:`x = (1 + \tanh(\sinh t)) / 2` that ``TanhSinh`` makes on :math:`[0, 1]` (left),
             and the integrand in :math:`t`, :math:`f(x(t))\,dx/dt` for :math:`f(x) = 1/\sqrt{x}` (right, on a
             logarithmic scale). The dots are equal steps of :math:`1/2` in :math:`t`.
```

The figure shows why `TanhSinh` is so much cheaper here. It substitutes $x = (1 + \tanh(\sinh t)) / 2$, which maps
the whole line of $t$ onto $(0, 1)$, and adds up the integrand in $t$, $f(x(t))\,dx/dt$, at equal steps: the
trapezoidal rule. On the left, equal steps in $t$ land ever closer to the ends. The point at $t = -3$ is already
$x \approx 2 \times 10^{-9}$, too close to 0 to see here, and each further step multiplies the distance to 0 by a
factor that itself shrinks. On the right is that integrand. Where $1/\sqrt x$ grows without bound, $dx/dt$ falls
faster, so the product stays finite and is below $10^{-10}$ at both $t = -4$ and $t = 4$. On the singular side it
falls more slowly, since $1/\sqrt x$ takes back half of the decay of $dx/dt$, but still faster than exponentially.
The trapezoidal rule converges very fast on a smooth function that dies off like this. So `TanhSinh` halves its step
in $t$ everywhere at once, reusing every earlier point, until two successive sums agree as closely as you asked.

The default method has a limit. For a singularity stronger than about $(x - a)^{-0.65}$ its estimate falls short
of its error, by a factor that further cutting does not change. On $x^{-3/4}$ over $[0, 1]$, whose integral is 4, it
reports success with an error above both its `abs_error` and the request; `TanhSinh` is accurate to within its
estimate. Use `TanhSinh` for such an integrand, or substitute $x = a + u^k$ with $k(1 + \alpha) \ge 1$ for a
singularity like $(x - a)^\alpha$, which makes the integrand in $u$ bounded.

## At an end other than 0

`TanhSinh`'s points can come only as close to an end as the `double` numbers there allow. Near 0 those numbers are
dense, down to $10^{-308}$; near 1 they are $1.1 \times 10^{-16}$ apart. A singularity at 1 is therefore integrated
only up to a sliver short of it, and for $1/\sqrt{1 - x}$ that sliver holds up to about
$2\sqrt{1.1 \times 10^{-16}}$, or $2 \times 10^{-8}$, of the integral; the method's correction for what lies beyond its
last point recovers part of it, and the result here is $10^{-8}$ off. Substitute $u = 1 - x$ and write the integrand
in $u$ itself, so that the singular end is 0:

```{literalinclude} ../examples/quadrature_rough_integrands.cpp
:start-after: // <other_end>
:end-before: // </other_end>
:dedent:
```

```{literalinclude} ../examples/quadrature_rough_integrands.cpp
:language: text
:start-after: other_end_output = R"(
:end-before: )";
```

Writing `fn(1.0 - u)` would not help: computing $1 - u$ rounds it back onto the same coarse numbers near 1. The formula
itself must be in $u$, as the [tutorial](tutorial.md) did for the pendulum.

## Kinks, steps and peaks inside

Inside the interval, the default method finds a kink, a step or a peak by itself, since the panel around it keeps
the largest estimate and is cut until it is small. It pays for it in evaluations. `TanhSinh` is weak there: it has no
panels to cut, and on a kink it usually runs through all its levels and returns `not_converged`. If you know where the
rough spot is, split the interval there, so that it becomes an end of both pieces, and add the two results and their
`abs_error`s. For $|x - 1/3|$:

```{literalinclude} ../examples/quadrature_rough_integrands.cpp
:start-after: // <kink>
:end-before: // </kink>
:dedent:
```

```{literalinclude} ../examples/quadrature_rough_integrands.cpp
:language: text
:start-after: kink_output = R"(
:end-before: )";
```

At once, the method needed 330 evaluations and gave 8 digits; split, each piece is a straight line, which the first
150 evaluations integrate exactly. Split at the same `double` that the integrand uses, here `c`: a split at a
different rounding of $1/3$ would leave the kink a hair inside one piece.

The kink sits at $1/3$ for a reason. The default method's first step cuts $[a, b]$ into 10 equal panels, and a kink
exactly on a panel boundary, $a + k(b - a)/10$, costs nothing extra: $|x^2 - 1|$ on $[0, 2]$, with its kink at 1, takes
the first 150 evaluations and is exact. Elsewhere an unsplit kink can cost more than evaluations at a tight request:
on $|x^2 - 2|$ over $[0, 2]$ at a relative tolerance of $10^{-13}$, the method reports success with an `abs_error` of
$6.1 \times 10^{-14}$ while its actual error is $7.1 \times 10^{-14}$. Split at the kink, the error is within a few
units in the last place (0 here). The kink at
$\sqrt 2$ is not a `double`, so split at the nearest one, `std::sqrt(2.0)`: the error this leaves is bounded by the
jump in the slope of $f$ times the distance to the true kink, here below $10^{-15}$.

A singularity inside the interval is the hardest case, and splitting matters most. On $1/\sqrt{|x - 1/3|}$ the
default method reports success with an error above its estimate and above the request. Split at $1/3$, each piece
has the singularity at an end other than 0, so write each in its distance $u$ from $1/3$, as above:

```{literalinclude} ../examples/quadrature_rough_integrands.cpp
:start-after: // <inside>
:end-before: // </inside>
:dedent:
```

```{literalinclude} ../examples/quadrature_rough_integrands.cpp
:language: text
:start-after: inside_output = R"(
:end-before: )";
```

## What no estimate sees

Every method learns about $f$ only at its points. A feature narrower than the spacing of the points, falling between
them, is invisible: here a box of width $2 \times 10^{-4}$ that none of the first 150 points hits.

```{literalinclude} ../examples/quadrature_rough_integrands.cpp
:start-after: // <narrow>
:end-before: // </narrow>
:dedent:
```

```{literalinclude} ../examples/quadrature_rough_integrands.cpp
:language: text
:start-after: narrow_output = R"(
:end-before: )";
```

Every value was 0, so the method returned 0 with its smallest possible estimate, and called it a success; `TanhSinh`
does the same. No method that samples can rule this out. If your integrand has a narrow feature at a place you know,
split there; if you do not know where it is, no quadrature method will find it for you.

## Next

[Infinite intervals](infinite-intervals.md) treats ends at infinity; [Handle failures](handle-failures.md) says what to
do when a method gives up.

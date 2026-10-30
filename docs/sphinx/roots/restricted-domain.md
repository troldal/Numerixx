# Solve on a restricted domain

Many functions are defined only on part of the line: a logarithm for $x > 0$, a square root for $x \ge 0$, a physical
model only where its quantities make sense. Outside that *domain*, a method that evaluates `fn` gets either NaN or a
number that means nothing. This page shows how to keep a solve inside the domain, and what each method does when it
steps outside.

Numerixx reads a value of `fn` that is not finite, NaN or infinite, as a sign that the point lies outside the domain.
`search` then steps back towards the last point where `fn` was finite; `solve` stops with an error. (One point is
exempt: `solve` first tries 0 in a bracket around 0, and skips it where `fn` is not finite there, as $\sin x / x$ is
not.)

## Return NaN outside the domain

If `fn`'s formula gives finite numbers outside its domain, where they mean nothing, make `fn` return NaN there.
Otherwise a method takes those numbers at their word, and may return a root, or a pole, that is no solution.

Van der Waals' equation of state is an example. It gives the pressure of a real gas from its *molar volume* $v$, the
volume one mole of it fills:

$$p(v) = \frac{RT}{v - b} - \frac{a}{v^2}$$

Here $b$, the *co-volume*, is the volume the molecules themselves take up, so no gas has $v \le b$. The formula does not
know that: below $b$ it gives finite, meaningless pressures, and at $b$ it has a pole, where it jumps from $-\infty$ to
$+\infty$. To find the volume of carbon dioxide at 280 K and 4 MPa, solve $p(v) - 4\,\text{MPa} = 0$, with NaN below
$b$:

```{literalinclude} ../examples/roots_restricted_domain.cpp
:start-after: // <nan>
:end-before: // </nan>
:dedent:
```

```{literalinclude} ../examples/roots_restricted_domain.cpp
:language: text
:start-after: nan_output = R"(
:end-before: )";
```

The search starts from $10^{-4}$. On its lower side it meets NaN, and probes towards the edge of the domain, where it
finds no sign change. On its upper side it finds the bracket of the root, after 29 evaluations in all, and Brent's
method finds the root in it. Without the NaN, the same code brackets the pole at $b$, where the formula changes sign,
and returns $v = b$: a volume at which the gas cannot exist.

```{note}
At 280 K carbon dioxide is below its critical temperature, about 304 K by these constants, and there van der Waals'
equation can have three roots at a given pressure: a liquid's volume, a vapour's, and an unphysical one between them.
A search returns the first sign change it meets, which need not be the one you want. At 4 MPa the equation has a
single root.
```

## Prefer a bracket

The two kinds of method behave differently near the edge of a domain. A bracketing method keeps the root in a *bracket*,
an interval at whose ends `fn` has opposite signs, and never leaves it: once the bracket lies inside the domain, with
finite values of `fn` at its ends, every point at which the method evaluates `fn` does too (a numerical derivative's
points are another matter, [below](#a-numerical-derivative-near-the-edge)). An open method steps wherever its formula
points, and can step outside. The first value that is not finite ends the solve with
{cpp:enumerator}`nxx::roots::error::non_finite`. Here is Newton's method on $\ln x - 1$, whose root is $e$, from the
guess 20:

```{literalinclude} ../examples/roots_restricted_domain.cpp
:start-after: // <leave>
:end-before: // </leave>
:dedent:
```

```{literalinclude} ../examples/roots_restricted_domain.cpp
:language: text
:start-after: leave_output = R"(
:end-before: )";
```

Far from the root the logarithm is flat, so Newton's first step overshoots: from 20 it lands at $-19.9$, where the
logarithm is NaN. A search and Brent's method from the same guess find $e$. When an open method returns `non_finite`,
search for a bracket instead.

## A numerical derivative near the edge

{cpp:func}`nxx::deriv::derivative_of` turns `fn` into its numerical derivative, which evaluates `fn` at points a small
step away from $x$ on each side. Near the edge of the domain those points can fall outside it, even when $x$ itself is
inside, and the derivative then fails. From a guess, Newton's method ends with
{cpp:enumerator}`nxx::roots::error::not_converged`. In a bracket it takes a bisection step instead, and still converges.
Here the root is $10^{-12}$ and the numerical derivative's step, about $6 \times 10^{-6}$, reaches below 0:

```{literalinclude} ../examples/roots_restricted_domain.cpp
:start-after: // <stencil>
:end-before: // </stencil>
:dedent:
```

```{literalinclude} ../examples/roots_restricted_domain.cpp
:language: text
:start-after: stencil_output = R"(
:end-before: )";
```

From a guess, a one-sided difference avoids the problem. `nxx::deriv::Forward<2>` evaluates `fn` only at $x$ and above
it, so near a domain that starts at 0 its points stay inside:

```{literalinclude} ../examples/roots_restricted_domain.cpp
:start-after: // <one-sided>
:end-before: // </one-sided>
:dedent:
```

```{literalinclude} ../examples/roots_restricted_domain.cpp
:language: text
:start-after: one_sided_output = R"(
:end-before: )";
```

For a domain that ends at an upper edge, `nxx::deriv::Backward<2>` evaluates `fn` only at $x$ and below it.

The program needs:

```{literalinclude} ../examples/roots_restricted_domain.cpp
:start-after: // <includes>
:end-before: // </includes>
```

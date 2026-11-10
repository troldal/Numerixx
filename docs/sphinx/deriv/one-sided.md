# Near the edge of a domain

Many functions are defined only on part of the line, their *domain*. At its *edge*, a method that evaluates `fn` on both
sides of $x$ reaches outside. One-sided methods do not.

## Return NaN outside the domain

A method can only tell that a point lies outside the domain if `fn` returns NaN (or an infinity) there. Here
$1/(1 + x)$ is defined for $x \ge 0$ only; its derivative at the edge, 0, is $-1$:

```{literalinclude} ../examples/deriv_one_sided.cpp
:start-after: // <edge>
:end-before: // </edge>
:dedent:
```

```{literalinclude} ../examples/deriv_one_sided.cpp
:language: text
:start-after: edge_output = R"(
:end-before: )";
```

The central difference and Ridders' method get NaN below 0 and report {cpp:enumerator}`nxx::deriv::error::non_finite`.
A *one-sided* method evaluates `fn` at $x$ and on one side only: the stencils (fixed formulas) `Forward<P>` and
`Backward<P>`, and the Ridders methods `ForwardRidders` and `BackwardRidders`. `Forward<2>` keeps about 10
digits here, `ForwardRidders` about 13, with an error estimate, `abs_error`.

## Near the edge

A little inside, Ridders' method works, at a price (the central difference still fails at $10^{-6}$). The program
compares it with `ForwardRidders` at three distances, with the relative error, estimate and evaluations:

```{literalinclude} ../examples/deriv_one_sided.cpp
:language: text
:start-after: distances_output = R"(
:end-before: )";
```

Far from the edge, Ridders' method is the more accurate. It restarts below every step that reaches past the edge, more
often the nearer $x$: at $10^{-10}$ its error is $2 \times 10^{-10}$, its estimate $8 \times 10^{-6}$, for all 64
evaluations it may make; `ForwardRidders` stays below $10^{-12}$ at 10. Yet on $\sqrt{1 - x}$ near 1, which changes
fast near its edge, Ridders' method is the more accurate one.

So try both. Keep the smaller estimate's result if the two differ by less than the sum of their
estimates; if not, trust neither until you have checked at another initial step ([Choose a step](step.md)). A formula
often cancels near an edge (subtracts nearly equal numbers and loses digits), as $1 - x^2$ does near 1: write it as
$(1 - x)(1 + x)$, or its noise makes an estimate too small.

Where `fn` changes fast near its edge, avoid the one-sided stencils: their step does not follow `fn`'s scale, which
shrinks with the distance to the edge, and nothing warns you. That holds for Newton's method
too: on $\sqrt x - 10^{-6}$ from $2 \times 10^{-12}$, `Forward<2>` gives 525 for 353553, and Newton leaves the
domain. Start it in a bracket, where a bad step falls back to bisection (the [roots
chapter](../roots/restricted-domain.md#a-numerical-derivative-near-the-edge)).

```{note}
From their default step, the one-sided Ridders methods lose accuracy where `fn` changes on a scale below about
$10^{-9} \max(|x|, 1)$ ({cpp:struct}`nxx::deriv::ForwardRidders`), as $\sqrt x$ does near 0 ([Handle
failures](handle-failures.md#a-result-that-is-not-accurate-enough)). And an estimate is not a proof: on a polynomial,
NaN below 0, with zeros at the points Ridders' method tries after its restarts, it returns 0 with an estimate of 0 for
a derivative of 1.
```

## A formula that goes on

A formula does not return NaN outside its domain by itself; here $1/(1 + x)$ goes on below 0 with slope $-3$:

```{literalinclude} ../examples/deriv_one_sided.cpp
:start-after: // <kink>
:end-before: // </kink>
:dedent:
```

```{literalinclude} ../examples/deriv_one_sided.cpp
:language: text
:start-after: kink_output = R"(
:end-before: )";
```

At the kink, Ridders' method returns $-2$, the mean of the two slopes, with an estimate of only $1.2 \times 10^{-8}$:
nothing tells it the kink is there. `ForwardRidders` looks only above 0 and gets $-1$. So return NaN outside the
domain, or use a one-sided method.

The program needs:

```{literalinclude} ../examples/deriv_one_sided.cpp
:start-after: // <includes>
:end-before: // </includes>
```

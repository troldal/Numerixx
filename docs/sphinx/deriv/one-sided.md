# Near the edge of a domain

Many functions are defined only on part of the line, such as $\sqrt x$ or a model valid above some limit. The part where
a function is defined is its *domain*. A central method evaluates `fn` on both sides of $x$, so at the *edge* of the
domain one of its points falls outside. One-sided methods avoid that.

## Return NaN outside the domain

The methods can only tell that a point lies outside the domain if `fn` returns NaN (or an infinity) there. The example
is $1/(1 + x)$, defined here for $x \ge 0$ only. Its derivative at the edge, 0, is $-1$:

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

The central difference gets NaN at $-h$ and reports {cpp:enumerator}`nxx::deriv::error::non_finite`; so does Ridders'
method. A *one-sided* method evaluates `fn` only at $x$ and on one side of it: `Forward<P>` and `ForwardRidders` above,
`Backward<P>` and `BackwardRidders` below. The stencil `Forward<2>` keeps about 10 digits here, and `ForwardRidders`
about 13, with an error estimate.

## Near the edge

A little inside the domain, the central methods work, but are harder to trust. The program also compares Ridders'
method and `ForwardRidders` at three distances from the edge: the relative error, the relative estimate, and the
evaluations of `fn`:

```{literalinclude} ../examples/deriv_one_sided.cpp
:language: text
:start-after: distances_output = R"(
:end-before: )";
```

Far from the edge, Ridders' method is the more accurate. Closer, it restarts below each step that met NaN; its estimate
grows, and at $10^{-10}$ it makes all the 64 evaluations it may. `ForwardRidders` keeps its accuracy and its 10
evaluations. On a function that changes fast near its edge, such as $\sqrt{1 - x}$ near 1, Ridders' method is the more
accurate one.

So near an edge, try both. Keep the one with the smaller estimate only if the two agree within their estimates. An
estimate is not a proof: on a function whose values happen to vanish at the points a method tries, it can even be 0
for a wrong result.

Avoid the one-sided stencils there: their step does not shrink near the edge, and they have no estimate to warn you. A
root solve meets the same problem ([A numerical derivative near the
edge](../roots/restricted-domain.md#a-numerical-derivative-near-the-edge)).

## A formula that goes on

A formula used outside its domain does not return NaN by itself. Here $1/(1 + x)$ goes on below 0 as another formula,
with slope $-3$ instead of $-1$:

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

Ridders' method sees a kink at 0, and returns $-2$, the mean of the two slopes, with a small estimate. `ForwardRidders`
looks only above 0, and gets $-1$. So return NaN outside the domain, or use a one-sided method.

```{note}
From their default step, the one-sided Ridders methods need `fn` to vary on a scale of at least about
$10^{-9} \max(|x|, 1)$; closer to a singularity, such as $\sqrt x$ near 0, see [Handle
failures](handle-failures.md#a-result-that-is-not-accurate-enough).
```

The program needs:

```{literalinclude} ../examples/deriv_one_sided.cpp
:start-after: // <includes>
:end-before: // </includes>
```

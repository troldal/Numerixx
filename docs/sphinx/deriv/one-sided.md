# Near the edge of a domain

Many functions are defined only on part of the line: a square root or a logarithm for $x > 0$, a correlation between the
ends of its range, a physical model above some limit. The part where a function is defined is its *domain*. A central
method evaluates `fn` on both sides of $x$, so at the edge of the domain one of its points falls outside it. This page
shows what happens there, and the methods that look only one way. The examples write `dv` for `nxx::deriv`.

## Return NaN outside the domain

The methods can only tell that a point lies outside the domain if `fn` says so, by returning NaN (or an infinity). Make
it do that. A formula that goes on giving numbers outside its domain gives a wrong derivative, with no error
([below](#a-formula-that-goes-on)).

The example is $1/(1 + x)$, defined here for $x \ge 0$ only, with NaN below 0, and its derivative at the edge, 0, which
is $-1$:

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

The central difference evaluates `fn` at $-h$, gets NaN, and reports {cpp:enumerator}`nxx::deriv::error::non_finite`;
Ridders' method, also central, does the same, since every step it tries reaches below 0. A *one-sided* method evaluates
`fn` only at $x$ and on one side of it. `Forward<P>` and `ForwardRidders` look at $x$ and above, for a domain that
starts at $x$; `Backward<P>` and `BackwardRidders` at $x$ and below, for a domain that ends there.

`Forward<2>` is the one-sided stencil of order 2, with three points, $x$, $x + h$ and $x + 2h$; it gives about 10 digits
here. `ForwardRidders` extrapolates a one-sided difference as Ridders' method extrapolates the central one, and gives
about 13 digits with an error estimate.

## Near the edge

At a small distance from the edge, the central methods can work. The central difference with its default step,
$6 \times 10^{-6}$, still reaches below 0 from $10^{-6}$, but Ridders' method meets the NaN, restarts with smaller
steps, and finds the derivative:

```{literalinclude} ../examples/deriv_one_sided.cpp
:start-after: // <near>
:end-before: // </near>
:dedent:
```

```{literalinclude} ../examples/deriv_one_sided.cpp
:language: text
:start-after: near_output = R"(
:end-before: )";
```

Its estimate, about $10^{-9}$, is larger than usual: to stay inside the domain its steps must stay below $10^{-6}$,
where rounding, which grows as $1/h$, is larger. A one-sided method can take larger steps away from the edge, and on
this function `ForwardRidders` is more accurate at $10^{-6}$; on a function that changes fast near its edge, such as
$\sqrt{1 - x}$ near 1, Ridders' method can be the more accurate one. So close to an edge, try both. If they agree to
within their estimates, keep the one with the smaller estimate; if they do not, trust neither until you have checked,
for example at another step. An estimate is not a proof: on a function whose values happen to vanish at the points a
method tries, it can even be 0 for a wrong result.

A one-sided *stencil* with its default step is riskier near an edge where `fn` changes fast: its step,
$\varepsilon^{1/(P+1)} \max(|x|, 1)$, does not shrink with the distance to the edge, and it has no estimate to warn you.
Prefer the Ridders methods there.

## A formula that goes on

A formula used outside its domain does not return NaN by itself. Here the function is $1/(1 + x)$ above 0, as before,
but below 0 it goes on as another formula, with slope $-3$ instead of $-1$:

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

Ridders' method now sees a function with a kink at 0, and returns $-2$, the mean of the two slopes, with an estimate of
about $10^{-8}$ that gives no hint of the error. Nothing told it where the domain ends. `ForwardRidders` looks only
above 0, and gets $-1$. Where `fn` is defined on one side only, either return NaN on the other side, or use a one-sided
method; a kink or a jump among a method's points gives a wrong result without an error.

```{note}
A one-sided method from its default step, $\max(|x|, 1)/32$, needs `fn` to vary on a scale of at least about
$10^{-9} \max(|x|, 1)$. Closer to a singularity, such as $\sqrt x$ near 0, its result is wrong, and its estimate says
so; [Handle failures](handle-failures.md#a-result-that-is-not-accurate-enough) shows how to tell, and the remedy: a step
on `fn`'s own scale.
```

The program needs:

```{literalinclude} ../examples/deriv_one_sided.cpp
:start-after: // <includes>
:end-before: // </includes>
```

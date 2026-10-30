# Choose a step

A stencil evaluates `fn` at points a step $h$ apart. Too large a step, and the formula's truncation error, the error of
the formula itself, dominates; too small, and rounding does ([Numerical derivatives](index.md#the-idea-and-its-catch)
shows the trade-off). The default step balances the two under one assumption: that `fn` varies on the scale
$\max(|x|, 1)$. Near $x = 0$ that scale is 1; for a large $|x|$ it is $|x|$. For the central difference and a first
derivative the default is $\varepsilon^{1/3} \approx 6 \times 10^{-6}$ times that scale. This page is about the
functions that break the assumption, and what to do then. The examples write `dv` for `nxx::deriv`.

## When the default step is too large

A function that changes much faster than $\max(|x|, 1)$, such as a narrow peak, a steep front, or $\exp(100 x)$, needs a
smaller step. Take a peak of width $a = 10^{-4}$, and its derivative at $2 \times 10^{-4}$:

```{literalinclude} ../examples/deriv_step.cpp
:start-after: // <scale>
:end-before: // </scale>
:dedent:
```

With the default step:

```{literalinclude} ../examples/deriv_step.cpp
:start-after: // <default>
:end-before: // </default>
:dedent:
```

The step, $6 \times 10^{-6}$, is about a sixteenth of the peak's width, so the truncation error is large: the result
keeps only three or four digits (the output is at the end of the page). The default step is $\varepsilon^{1/3}$ times
the assumed scale; for this function, use $\varepsilon^{1/3}$ times its own scale, $a$.
{cpp:func}`nxx::math::root_epsilon` gives $\varepsilon^{1/3}$:

```{literalinclude} ../examples/deriv_step.cpp
:start-after: // <scaled>
:end-before: // </scaled>
:dedent:
```

The third argument of `derivative` is the step. With it the result keeps about 10 digits, as the central difference does
on a function that suits its default. For a stencil of order $P$ the balanced step is $\varepsilon^{1/(P+1)}$ times the
scale ($\varepsilon^{1/5}$ for `Central<4>`), and for the $N$-th derivative $\varepsilon^{1/(P+N)}$ times it;
{cpp:func}`nxx::deriv::default_step` gives the default for any method and order.

## Or let Ridders' method find the step

You may not know the scale of `fn`. Ridders' method starts from a large step, $\max(|x|, 1)/32$ for a first derivative,
and halves it until rounding stops the improvement, so it finds a suitable step by itself:

```{literalinclude} ../examples/deriv_step.cpp
:start-after: // <ridders>
:end-before: // </ridders>
:dedent:
```

Here it keeps 13 digits, with an estimate that covers its error, at the cost of 36 evaluations: more than its usual 10
to 16, because it had to come down from a step of $1/32$ to the peak's scale. If you know the scale, give Ridders'
method a step: it starts from that step instead, and needs fewer evaluations, 14 here:

```{literalinclude} ../examples/deriv_step.cpp
:start-after: // <ridders-step>
:end-before: // </ridders-step>
:dedent:
```

The results:

```{literalinclude} ../examples/deriv_step.cpp
:language: text
:start-after: output = R"(
:end-before: )";
```

```{note}
A step that is too small for $x$, so that the method's points round to the same number (at 1, $h = 10^{-300}$ does),
gives the error `non_finite`; so does a step that is
not positive and finite. Below the balanced step, rounding grows as $1/h$ long before that.
```

The program needs:

```{literalinclude} ../examples/deriv_step.cpp
:start-after: // <includes>
:end-before: // </includes>
```

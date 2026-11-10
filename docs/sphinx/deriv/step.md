# Choose a step

A stencil, a fixed formula such as the central difference, evaluates `fn` at points a step $h$ from $x$. The default
step balances truncation against rounding error ([Numerical derivatives](index.md#the-idea-and-its-catch)) under one
assumption: that `fn` varies on the *scale* $\max(|x|, 1)$, the distance over which `fn` changes noticeably. For the
central difference and a first derivative it is $\varepsilon^{1/3} \approx 6 \times 10^{-6}$ times that scale, where
$\varepsilon$ is the type's machine epsilon. Some functions break the assumption.

## When the default step is too large

A function that changes much faster than $\max(|x|, 1)$, such as a narrow peak or a steep front, needs a smaller step.
Take a peak of width $a = 10^{-4}$, and its derivative at $2 \times 10^{-4}$:

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

```{literalinclude} ../examples/deriv_step.cpp
:language: text
:start-after: default_output = R"(
:end-before: )";
```

The result keeps only three or four digits. The figure shows why: the central difference's error on this peak and on
$\sqrt x$ at 2, against the step.

```{eval-rst}
.. plot:: plots/deriv_scale.py
   :alt: Two V-shaped curves of one shape on logarithmic axes, smaller steps to the right. For the square root at 2 the
         lowest error, 4.1e-12, is at a step of 5e-6, and the default step, 1.2e-5, is at the bottom of the V. The
         peak's error stays at 1 or above down to a step of 2e-4, then falls to its lowest, 1.2e-11, at 1e-9, about four
         decades further right. The peak's default step, 6.1e-6, is circled on its falling side, at an error of 3.4e-4.
   :caption: The relative error of the central difference against its step on :math:`\sqrt{x}` at 2 and on the peak
             :math:`x / (x^2 + a^2)`, :math:`a = 10^{-4}`, at :math:`2a`, each with its default step circled.
```

The two curves have one shape, but the peak's lies about four decades to the right, at smaller steps: the best step
depends on the function's scale, not on $x$ alone. The default step sits at the bottom of the square root's curve. On
the peak it is far too large, high on the side where truncation error dominates. For this function, use
$\varepsilon^{1/3}$ times its own scale, $a$, as the step. {cpp:func}`nxx::math::root_epsilon` gives
$\varepsilon^{1/3}$, and the third argument of `derivative` is the step:

```{literalinclude} ../examples/deriv_step.cpp
:start-after: // <scaled>
:end-before: // </scaled>
:dedent:
```

```{literalinclude} ../examples/deriv_step.cpp
:language: text
:start-after: scaled_output = R"(
:end-before: )";
```

Now the result keeps about 10 digits, at the bottom of the peak's curve. For another method or order,
{cpp:func}`nxx::deriv::default_step`, called as `default_step<Method, N>(1.0)` (in another real type, `T(1)`), gives the
default step for a function of scale 1; multiply it by your function's scale.

## Or let Ridders' method find the step

You may not know the scale of `fn`. Ridders' method starts from a large step, $\max(|x|, 1)/32$ for a first derivative,
and halves it until rounding stops the improvement:

```{literalinclude} ../examples/deriv_step.cpp
:start-after: // <ridders>
:end-before: // </ridders>
:dedent:
```

```{literalinclude} ../examples/deriv_step.cpp
:language: text
:start-after: ridders_output = R"(
:end-before: )";
```

It keeps 13 digits, with an estimate that covers its error, at the cost of 36 evaluations: it had to come down from
$1/32$ to the peak's scale. If you know the scale, pass a step, and the method starts from it instead, here with 14
evaluations:

```{literalinclude} ../examples/deriv_step.cpp
:start-after: // <ridders-step>
:end-before: // </ridders-step>
:dedent:
```

```{literalinclude} ../examples/deriv_step.cpp
:language: text
:start-after: ridders_step_output = R"(
:end-before: )";
```

```{note}
A step too small for the points to differ, or not positive and finite, gives `non_finite` ([Handle
failures](handle-failures.md)).
```

The program needs:

```{literalinclude} ../examples/deriv_step.cpp
:start-after: // <includes>
:end-before: // </includes>
```

# nxx::deriv: evaluations and accuracy per method

How many times each method calls `fn` for a first derivative with its default step, in `double`, and the relative error
of the result. The functions have closed-form derivatives:

- $1/x$ at 2 and $\sqrt x$ at 2, smooth on the scale $\max(|x|, 1)$ that the default step assumes;
- a narrow peak, $x / (x^2 + a^2)$ with $a = 10^{-4}$, at $2 \times 10^{-4}$, which varies about ten thousand times
  faster ([Choose a step](step.md)).

Each cell gives the evaluations, then the relative error:

```{literalinclude} ../examples/deriv_costs.cpp
:language: text
:start-after: output = R"(
:end-before: )";
```

On the peak the higher-order stencils lose every digit; Ridders' method keeps 13 at 36 evaluations.

The table is printed by `docs/sphinx/examples/deriv_costs.cpp`, which checks that every estimate of the Ridders
methods covers its error. The numbers hold for builds that do not fuse a multiplication and an addition ([Working with
Numerixx](../basics.md#compile-time)).

The Ridders methods' typical and maximum evaluations for each derivative order are in the reference
({cpp:struct}`nxx::deriv::Ridders`, {cpp:struct}`nxx::deriv::ForwardRidders`).

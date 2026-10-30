# nxx::deriv: evaluations and accuracy per method

How many times each method calls `fn` for a first derivative with its default step, in `double`, and the relative error
of the result, on three functions with closed-form derivatives:

- $1/x$ at 2, smooth on the scale the default step assumes;
- $\sqrt x$ at 2, likewise;
- a narrow peak, $x / (x^2 + a^2)$ with $a = 10^{-4}$, at $2 \times 10^{-4}$, which varies about ten thousand times
  faster than that scale ([Choose a step](step.md)).

Each cell gives the evaluations, then the relative error:

```{literalinclude} ../examples/deriv_accuracy.cpp
:language: text
:start-after: output = R"(
:end-before: )";
```

`docs/sphinx/examples/deriv_accuracy.cpp` prints this table and compares it with the one above, and checks that every
estimate of the Ridders methods covers its actual error. The functions use only arithmetic and `std::sqrt`, which is
correctly rounded, so the table is the same on every platform built without fused multiply-add.

Three things stand out:

- **Higher order is not always better.** `Central<4>` and `Central<6>` gain one to three digits over `Central<2>` on the
  smooth functions, but `Central<10>` gives some of that back, and on $1/x$ all of it and more.
- **On the peak, the higher-order stencils fail completely.** Their default steps, $\varepsilon^{1/(P+1)}$, are larger
  than the peak's width, and their results have no correct digit. `Central<2>` keeps three; Ridders' method, which finds
  its own step, keeps nearly all, for 36 evaluations.
- **The one-sided methods are less accurate** than the central ones, and cost one evaluation more for the stencils.

For derivatives of order 2 to 4, [Higher derivatives](higher-derivatives.md) gives the errors of `Central<2>` and
Ridders' method on $1/x$.

Ridders' method typically makes 10 to 16 evaluations for a first derivative, and at most 64; 13 to 23 for a second
derivative, and 25 to 45 for a third or fourth, at most 65, 128 and 129. The one-sided Ridders methods typically make 10
for a first derivative, at most 33 ({cpp:struct}`nxx::deriv::Ridders`, {cpp:struct}`nxx::deriv::ForwardRidders`).

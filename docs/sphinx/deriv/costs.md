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

`docs/sphinx/examples/deriv_costs.cpp` prints this table, and checks that every estimate of the Ridders methods covers
its error. It holds for builds that do not fuse a multiplication and an addition.

Ridders' method typically makes 10 to 16 evaluations for $f'$ (at most 64), 13 to 23 for $f''$ (at most 65), and 25
to 45 for $f'''$ and $f''''$ (at most 128 and 129). The one-sided Ridders methods typically make 10, 17, 22 and 29 for
$f'$ to $f''''$, at most 33, 65, 97 and 129 ({cpp:struct}`nxx::deriv::Ridders`,
{cpp:struct}`nxx::deriv::ForwardRidders`).

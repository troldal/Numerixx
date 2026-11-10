# Choose a method

The default, the central difference `Central<2>`, costs two evaluations of `fn` and gives about 10 correct digits of a
`double`. That holds where `fn` varies on the *scale* the default step assumes: the distance over which `fn` changes
noticeably, taken to be about $\max(|x|, 1)$. When you need more digits, an error estimate, or a function defined on
one side of $x$ only, choose another method.

A *stencil* is a fixed formula: a few points around $x$, a weight for each, and a step $h$. Its *truncation error* is
the error of the formula itself, apart from rounding. Its *accuracy order* $P$ says how fast that error falls with the
step: as $h^P$. Trying another method is a one-word change, `derivative<Method>(fn, x)`.

## From what you need to a method

| You need | Use | Why |
|---|---|---|
| A few correct digits, at the fewest evaluations | {cpp:struct}`nxx::deriv::Central`, `Central<2>` (the default) | 2 evaluations for $f'$; about 10 digits of a `double` |
| More digits at a fixed cost, `fn` smooth on the scale $\max(\lvert x \rvert, 1)$ | `Central<4>` or `Central<6>` | 4 or 6 evaluations; one to three digits more |
| The most digits, or an error estimate | {cpp:struct}`nxx::deriv::Ridders` | Typically 13 to 15 digits ([Costs](costs.md)), with `abs_error`; typically 10 to 16 evaluations for $f'$ |
| $x$ at the edge of `fn`'s domain | {cpp:struct}`nxx::deriv::ForwardRidders`, {cpp:struct}`nxx::deriv::BackwardRidders`; at a fixed cost {cpp:struct}`nxx::deriv::Forward`, {cpp:struct}`nxx::deriv::Backward` | Points on one side only ([Near the edge of a domain](one-sided.md)) |
| $x$ near the edge | `Ridders`, and `ForwardRidders` or `BackwardRidders` for the side the domain lies on | Try both, and compare ([Near the edge](one-sided.md#near-the-edge)) |
| `fn` varying faster than $\max(\lvert x \rvert, 1)$ | `Ridders`, or a stencil with a step of your own | [Choose a step](step.md) |

A higher accuracy order is not always better. Its points spread wider and its default step is larger, so it can be less
accurate even on a function that varies slowly (`Central<10>` on $1/x$, [Costs](costs.md)), and on one that changes
fast it can lose every digit.

## Measure on your own function

The table is a guide; your own `fn` decides. The methods do not report how many evaluations they made, so wrap `fn` in
a lambda that counts its calls:

```{literalinclude} ../examples/deriv_choose_a_method.cpp
:start-after: // <count>
:end-before: // </count>
```

Here the function is $1/(1 + x^2)$ at 0.5, whose derivative is $-2x/(1 + x^2)^2 = -0.64$:

```{literalinclude} ../examples/deriv_choose_a_method.cpp
:start-after: // <compare>
:end-before: // </compare>
:dedent:
```

It prints

```{literalinclude} ../examples/deriv_choose_a_method.cpp
:language: text
:start-after: view output = R"(
:end-before: )";
```

`Central<4>` and `Central<6>` gain about two digits over `Central<2>`, and Ridders' method a little more, with an error
estimate, for 14 evaluations.

## More digits than a `double`

A `double` limits every method to its 16 digits. In a multiprecision type the same calls give more. Here is the
derivative of $\sqrt x$ at 2, $1/(2\sqrt 2) = 0.35355339059327376220042218105242451964\ldots$, by Ridders' method in
`cpp_bin_float_50`, a Boost.Multiprecision type with 50 decimal digits ([Working with
Numerixx](../basics.md#multiprecision-types) gives the rules for `fn` in such a type):

```{literalinclude} ../examples/deriv_choose_a_method.cpp
:start-after: // <more-digits>
:end-before: // </more-digits>
:dedent:
```

```{literalinclude} ../examples/deriv_choose_a_method.cpp
:language: text
:start-after: more_digits_output = R"(
:end-before: )";
```

The estimate is near $10^{-46}$, and all 40 printed digits are correct. A stencil gains less: its default step scales
with the type too, and keeps about two thirds of the type's digits.

The program needs:

```{literalinclude} ../examples/deriv_choose_a_method.cpp
:start-after: // <includes>
:end-before: // </includes>
```

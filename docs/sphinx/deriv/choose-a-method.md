# Choose a method

The default, the central difference `Central<2>`, costs two evaluations of `fn` and gives about 10 or 11 correct digits
of a `double` on a function that varies on the scale $\max(|x|, 1)$. That is often enough. This page is for the cases
where it is not: when you need more digits, an error estimate, a function defined on one side of $x$ only, or the fewest
evaluations of an expensive `fn`.

A *stencil* is a fixed formula: a few points around $x$, a weight for each, and a step $h$. Its *truncation error* is
the error of the formula itself, apart from rounding, and its *order* $P$ says how fast that error falls with the step:
as $h^P$. Every method is called the same way, `derivative<Method>(fn, x)`, so trying another is a one-word change. The
examples write `dv` for `nxx::deriv`.

## From what you need to a method

| You need | Use | Why |
|---|---|---|
| A few correct digits, at the fewest evaluations | {cpp:struct}`nxx::deriv::Central`, `Central<2>` (the default) | 2 evaluations for $f'$; about 10 digits of a `double` |
| More digits at a fixed cost, `fn` smooth on the scale $\max(|x|, 1)$ | `Central<4>` or `Central<6>` | 4 or 6 evaluations; one to three digits more |
| The most digits, or an error estimate | {cpp:struct}`nxx::deriv::Ridders` | Typically 13 to 15 digits, with `abs_error`; typically 10 to 16 evaluations for $f'$ |
| `fn` defined on one side of $x$ only, $x$ at the edge | {cpp:struct}`nxx::deriv::ForwardRidders` or {cpp:struct}`nxx::deriv::BackwardRidders` | An error estimate, from points on one side; see [Near the edge of a domain](one-sided.md) |
| $x$ near the edge of `fn`'s domain | `Ridders` and the one-sided Ridders method | Try both, and compare; see [Near the edge](one-sided.md#near-the-edge) |
| The same at a fixed cost | {cpp:struct}`nxx::deriv::Forward` or {cpp:struct}`nxx::deriv::Backward` | `Forward<2>`: 3 evaluations; typically fewer digits than `Central<2>` |
| `fn` varying faster than $\max(|x|, 1)$ | `Ridders`, or a stencil with a step of your own | See [Choose a step](step.md) |

The higher-order stencils are not always better. They assume that `fn` is smooth across all their points, which spread
wider with the order, and their default steps are larger: on a function that changes faster than the scale
$\max(|x|, 1)$, they can be far worse than `Central<2>` ([Accuracy](accuracy.md) shows a case where they lose every
digit). Ridders' method finds a suitable step by itself, which makes it the safest choice when you know little about
`fn`.

## Measure on your own function

The table is a guide; your own `fn` decides. The methods do not report how many evaluations they made, but you can count
them: wrap `fn` in a lambda that counts its calls, and compare the methods' errors on a point where you know the
derivative.

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

Each higher-order stencil gains about two digits over `Central<2>`, though `Central<6>` is no better than `Central<4>`
here, and Ridders' method a little more, with an error estimate besides, for 14 evaluations. Where `fn` is cheap, that
is a good trade; where one evaluation takes a simulation run, a stencil's fixed cost may matter more. The [accuracy
table](accuracy.md) compares the main methods on three functions.

## More digits than a `double`

A `double` limits every method to its 16 digits. In a multiprecision type the same calls give more: here is the
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

The estimate is near $10^{-46}$, and all 40 printed digits are correct. A stencil gains less: the central difference's
default step scales with the type too, and keeps about two thirds of its digits.

The program needs:

```{literalinclude} ../examples/deriv_choose_a_method.cpp
:start-after: // <includes>
:end-before: // </includes>
```

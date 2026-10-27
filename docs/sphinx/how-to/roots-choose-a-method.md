# Choose a method

Every method is called the same way, `solve<Method>(fn, start)`, so trying another method of the same kind (from a
bracket, or from a guess) is a one-word change; Newton and Halley also need the derivatives ([Give
derivatives](roots-derivatives.md)). Decide by what you have and what fn costs.

| You have | Use | Why |
|---|---|---|
| A bracket, fn of value only | {cpp:struct}`nxx::roots::Brent` | The default: superlinear on a smooth fn; at worst about the square of Bisection's count |
| A bracket, and a worst case to bound | {cpp:struct}`nxx::roots::Toms748` | About Brent's cost, with at most about 4 evaluations per halving of the bracket |
| A bracket, fn strongly curved near one end (a logarithm, 1/x), or a sign change at 0 | {cpp:struct}`nxx::roots::AndersonBjorck` | Far fewer evaluations than Brent there |
| A bracket, a root of odd multiplicity, such as (x - 1)^3 | {cpp:struct}`nxx::roots::Bisection`, or {cpp:struct}`nxx::roots::Halley` with f' and f'' | The interpolating methods slow down there; of them {cpp:struct}`nxx::roots::Ridders` is the cheapest |
| A bracket, a rough fn, a count you must predict | {cpp:struct}`nxx::roots::Bisection` | One bit per evaluation, whatever fn is |
| A guess, and f' (and f'') cheap to compute | {cpp:struct}`nxx::roots::Newton` (or {cpp:struct}`nxx::roots::Halley`) | Fewest evaluations from a good guess |
| A guess, fn of value only | {cpp:struct}`nxx::roots::Secant` | No derivative needed |
| A guess far from the root, or none | {cpp:func}`nxx::roots::search`, then a bracketing method | See [Find a bracket](roots-find-a-bracket.md) |

Newton and Halley also take a bracket, and are then safeguarded by it. For a fn of value only, prefer Brent or Secant
to Newton with a numerical derivative, which costs several calls of fn per point ([Give
derivatives](roots-derivatives.md)). {cpp:struct}`nxx::roots::Ridders`, {cpp:struct}`nxx::roots::RegulaFalsi`,
{cpp:struct}`nxx::roots::Illinois` and {cpp:struct}`nxx::roots::Steffensen` are teaching examples.

A root of even multiplicity, such as that of (x - 1)^2, has no sign change: no bracket shows it, and the methods find
it only where fn happens to evaluate to exactly 0.

## Measure the cost on your own fn

The methods report no count. Wrap fn in a lambda that counts its calls:

```{literalinclude} ../examples/roots_choose_a_method.cpp
:start-after: // <count>
:end-before: // </count>
```

and compare the methods on your fn and bracket:

```{literalinclude} ../examples/roots_choose_a_method.cpp
:start-after: // <compare>
:end-before: // </compare>
:dedent:
```

On x^3 - 2x - 5 over [2, 3], Bisection takes 51 evaluations, Brent 8, Toms748 and AndersonBjorck 9 each. The
[costs table](../reference/roots-costs.md) gives every method's count on five such functions.

The program needs:

```{literalinclude} ../examples/roots_choose_a_method.cpp
:start-after: // <includes>
:end-before: // </includes>
```

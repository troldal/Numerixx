# Choose a method

Every method is called the same way, `solve<Method>(fn, start)`, so trying another is a one-word change. Decide by
what you have and what fn costs.

| You have | Use | Why |
|---|---|---|
| A bracket, fn of value only | {cpp:struct}`nxx::roots::Brent` | The default: superlinear on a smooth fn, never much slower than bisection |
| A bracket, and a worst case to bound | {cpp:struct}`nxx::roots::Toms748` | About Brent's cost, with at most about 4 evaluations per halving of the bracket |
| A bracket, fn strongly curved near one end (a logarithm, 1/x), or a sign change at 0 | {cpp:struct}`nxx::roots::AndersonBjorck` | Far fewer evaluations than Brent there |
| A bracket, a root of odd multiplicity or a flat fn | {cpp:struct}`nxx::roots::Ridders` | Cheapest on (x - 1)^3 and its kind |
| A bracket, a rough fn, a count you must predict | {cpp:struct}`nxx::roots::Bisection` | One bit per evaluation, whatever fn is |
| A guess, and f' (and f'') cheap to compute | {cpp:struct}`nxx::roots::Newton` (or {cpp:struct}`nxx::roots::Halley`) | Fewest evaluations from a good guess |
| A guess, fn of value only | {cpp:struct}`nxx::roots::Secant` | No derivative needed |
| A guess far from the root, or none | {cpp:func}`nxx::roots::search`, then a bracketing method | See [Find a bracket](roots-find-a-bracket.md) |

Newton and Halley also take a bracket, and are then safeguarded by it. For a fn of value only, prefer Brent or Secant
to Newton with a numerical derivative, which costs several calls of fn per point ([Give
derivatives](roots-derivatives.md)). {cpp:struct}`nxx::roots::RegulaFalsi`, {cpp:struct}`nxx::roots::Illinois` and
{cpp:struct}`nxx::roots::Steffensen` are teaching examples.

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

# nxx::roots: evaluations per method

The number of calls of fn each method makes to reach full precision with the default tolerances, in `double`. One call
gives fn's value, together with f' for Newton and f' and f'' for Halley. "Newton, numerical" is Newton with f' from
{cpp:func}`nxx::deriv::derivative_of`, counting every call of fn, its nodes included.

The counts come from independent implementations of the methods (`docs/roots/design.md`, Evidence) and are asserted by
`docs/sphinx/examples/roots_costs.cpp`, which also checks each result against the closed-form root. They hold for fn
written as that program writes it, on x86-64 and WebAssembly builds without FMA contraction; another form of the same
fn, or contracted arithmetic, can change them by a few. With a library function such as `std::sin` in fn, counts also
differ between platforms.

## From a bracket

| fn | Bracket | Bisection | Brent | Toms748 | AndersonBjorck | Ridders | Illinois | RegulaFalsi | Newton | Halley | Newton, numerical |
|---|---|---|---|---|---|---|---|---|---|---|---|
| x^2 - 2 | [1, 2] | 52 | 9 | 9 | 9 | 12 | 11 | 20 | 8 | 6 | 24 |
| x^3 - 2x - 5 | [2, 3] | 51 | 8 | 9 | 9 | 14 | 10 | 23 | 7 | 6 | 21 |
| x^20 - 1 | [0, 5] | 55 | 19 | 21 | 30 | 18 | 30 | 45 | 11 | 11 | 33 |
| (x - 1)^3 | [0, 3] | 54 | 156 | 149 | 152 | 86 | 148 | 145 | 91 | 54 | 549 |
| 1/x - 1 | [0.5, 3] | 54 | 11 | 6 | 5 | 12 | 12 | 25 | 9 | 6 | 27 |

## From a guess

| fn | Guess | Newton | Halley | Secant | Steffensen | Newton, numerical |
|---|---|---|---|---|---|---|
| x^2 - 2 | 1 | 7 | 5 | 9 | 10 | 21 |

Steffensen's count, and which root it finds, depend on fn's scale: on 1000 (x^2 - 2) from 1 it finds -sqrt(2).

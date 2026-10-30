# nxx::roots: evaluations per method

How many times each method calls `fn` to reach full precision with the default tolerances, in `double`, on five test
functions. When one evaluation of `fn` is expensive, this is the cost that matters. [Choose a
method](choose-a-method.md) turns these counts into advice, and shows how to measure them on your own function.

One call gives `fn`'s value, together with $f'$ for Newton's method and $f'$ and $f''$ for Halley's. "Newton, numerical"
is Newton's method with $f'$ from {cpp:func}`nxx::deriv::derivative_of`, counting every call of `fn`, those the
numerical derivative makes included.

The counts come from independent implementations of the methods (`docs/roots/design.md`, Evidence), and
`docs/sphinx/examples/roots_costs.cpp` asserts each of them and checks each result against the closed-form root. They
hold for `fn` written as that program writes it, on x86-64 and WebAssembly builds that do not fuse a multiplication and
an addition into one instruction. Another form of the same function, or fused arithmetic, can change them by a few. With
a library function such as `std::sin` in `fn`, counts also differ between platforms.

## From a bracket

| $f$ | Bracket | Bisection | Brent | Toms748 | AndersonBjorck | Ridders | Illinois | RegulaFalsi | Newton | Halley | Newton, numerical |
|---|---|---|---|---|---|---|---|---|---|---|---|
| $x^2 - 2$ | [1, 2] | 52 | 9 | 9 | 9 | 12 | 11 | 20 | 8 | 6 | 24 |
| $x^3 - 2x - 5$ | [2, 3] | 51 | 8 | 9 | 9 | 14 | 10 | 23 | 7 | 6 | 21 |
| $x^{20} - 1$ | [0, 5] | 55 | 19 | 21 | 30 | 18 | 30 | 45 | 11 | 11 | 33 |
| $(x - 1)^3$ | [0, 3] | 54 | 156 | 149 | 152 | 86 | 148 | 145 | 91 | 54 | 549 |
| $1/x - 1$ | [0.5, 3] | 54 | 11 | 6 | 5 | 12 | 12 | 25 | 9 | 6 | 27 |

## From a guess

| $f$ | Guess | Newton | Halley | Secant | Steffensen | Newton, numerical |
|---|---|---|---|---|---|---|
| $x^2 - 2$ | 1 | 7 | 5 | 9 | 10 | 21 |

Steffensen's count, and which root it finds, depend on $f$'s scale: on $1000\,(x^2 - 2)$ from 1 it finds $-\sqrt 2$.

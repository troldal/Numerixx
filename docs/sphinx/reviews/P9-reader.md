# Phase 9, stage 5: docs-reader on the deriv chapter (commit d00357e)

`docs-reader` (Opus, high), 2026-11-10: learning `nxx::deriv` (and `nxx::roots` for task 5) from the manual alone.
Saved by the main session as the agent returned it. Programs in the session scratchpad (`p9-review/docs-reader/`:
`t1.cpp` to `t5.cpp`, `t6_snippet.cpp`, `extract.py`).

**Summary.** The pages alone got me to correct, compiling code for all five tasks, and every program compiled first
time under GCC 16 with `-Wall -Wextra -Werror`. Every result matches its independent value. Two problems:
- **Task 3:** `BackwardRidders`' estimate (1.63e-10) is smaller than its actual error (2.7e-10). The root cause is
  that `sqrt(1 - x*x)` is noisy near 0.999, and no page tells the reader that a formula with cancellation counts as a
  "less accurate `fn`".
- **A snippet that does not compile:** derivative-of.md line 25, copied as written, fails with `'derivative' was not
  declared in this scope`.

What I read: `index.md`, `getting-started.md`, `basics.md`, `core.md`, every page in `deriv/`, and from `roots/` the
pages `index.md`, `derivatives.md`, `restricted-domain.md` and `choose-a-method.md`. Of the examples only the regions
the `literalinclude` markers show, extracted by a script that applies each directive's markers. Also the built
reference pages `deriv/reference.html`, `roots/reference.html` and `core.html`. Reference values from mpmath at 40
digits at the exact double inputs, and from each closed form in double.

## Task 1: d/dx x^x at 1.5

- Docs enough: yes (getting-started.md lines 20-24, deriv/index.md lines 63-73). Guesses: none.
- Error 2.75e-14 against the closed form (2.71e-14 against mpmath); estimate 1.39e-13, so it covers the error, at
  about 5 times the actual error, as deriv/index.md line 76 says.

## Task 2: second derivative of ln(1 + x^2) at 0.3, and a 1e-8 test

- Docs enough: yes (higher-derivatives.md lines 10-11 and 41-43; `good_enough` from handle-failures.md lines 42-49).
- Relative error 3.0e-12; estimate 6.75e-11 (4.4e-11 relative); good to 1e-8 by the estimate, and truly. Control:
  `Central<2>` has relative error 1.24e-8, no estimate, `good_enough` false, correctly.
- Guesses: whether to compare the estimate against `|value|` or the exact value (the page uses the value;
  handle-failures.md line 66 adds a caveat near 0). The default step for f'' in Ridders: step.md line 70 gives
  `max(|x|, 1)/32` for a first derivative and is silent for higher orders; only the reference says `/4`.

## Task 3: d/dx sqrt(1 - x^2) at 0.999, NaN beyond 1

- Docs enough for choosing the method: yes (choose-a-method.md line 21; one-sided.md lines 42-46, whose line 44 even
  predicts that Ridders' method can be the more accurate one on such a function). fn as the task gives it:
  `std::sqrt(1.0 - t*t)`, NaN beyond 1.

| Method | Error | Estimate | Calls |
|---|---|---|---|
| Ridders | 4.9e-12 | 1.580e-10 | 30 |
| BackwardRidders | **2.7e-10** | **1.628e-10** (does not cover) | 19 |
| Central<2> | 1.0e-4 | none | 2 |
| Backward<2> | 2.0e-4 | none | 3 |

- Root cause of the miss, proven: against mpmath, `sqrt(1 - x*x)` in double has relative errors up to 1.58e-14 near
  0.999, about 71 epsilon, from the cancellation in `1 - x*x`. Written as `sqrt((1 - x)(1 + x))`, Ridders gives
  4.5e-13 (estimate 1.3e-12) and BackwardRidders 5.3e-13 (estimate 5.7e-11): both estimates cover. The reference says
  Ridders "is optimistic for a noisier fn".
- The page that would have prevented it: handle-failures.md lines 66-68, the note "A less accurate `fn`, such as a
  solve to a loose tolerance or interpolated data, can make it too small", names only computed or tabulated
  functions. One clause saying a closed formula that cancels (such as 1 - x^2 near 1) is also a less accurate `fn`
  would have warned me; one-sided.md would be a good place too, since functions near a domain edge often cancel so.
- A guess that decides the outcome: one-sided.md line 46, "Keep the one with the smaller estimate only if the two
  agree within their estimates." The results differ by 2.71e-10: within the sum of the estimates (3.2e-10), not within
  either alone (1.6e-10). Under the "sum" reading I keep Ridders, which is right; under the "each" reading they
  disagree, and the page does not say what to do then.

## Task 4: d/dx exp(-(x/1e-3)^2) at 5e-4

- Docs enough: yes; step.md, "When the default step is too large" (lines 9-66), told me what to do (lines 49-51:
  eps^(1/3) times the function's own scale a, with `nxx::math::root_epsilon` and `derivative`'s third argument), and
  "Or let Ridders' method find the step" (lines 68-99) the alternative.
- Default step 6.06e-6: relative error 3.1e-5. Own step `root_epsilon<3, double>() * a` = 6.06e-9: 4.1e-11. Ridders
  from its default start: 7.0e-15 (estimate 3.3e-14); from a/8: 6.4e-15 (estimate 5.7e-14).
- Guesses: the template order of `root_epsilon`, copied from the example; core.html confirms `template<int N, real T>`.

## Task 5: x e^x = 2 by Newton with a numerical derivative, against Brent

- Docs enough: yes (roots/derivatives.md lines 54-61 show `rt::solve<rt::Newton>(f, nxx::deriv::derivative_of(f),
  rt::Bracket{...})`; the guess form appears only in restricted-domain.md and the reference).
- Newton from 1.0 with the default `Central<2>`: 0.85260550201372542 +- 4.4e-16; with `derivative_of<dv::Ridders>` the
  same, +- 5.6e-16. Brent on [0, 1]: 0.85260550201372554 +- 3.8e-16. |Newton - Brent| = 1.1e-16, within the sum of
  the bounds. Against Lambert W(2) from mpmath: Newton 1 ulp off, Brent 0.
- Guess: whether the guess form takes a separate derivative; derivatives.md shows only the bracket form.

## Snippet that does not compile

derivative-of.md line 25: "For a step of your own, write the lambda yourself: `[&](double x) { return derivative(fn,
x, h); }`." Copied into user code with `namespace dv = nxx::deriv;`, GCC says `error: 'derivative' was not declared
in this scope; did you mean 'nxx::deriv::derivative'?` (no argument-dependent lookup for a lambda and a `double`). It
should read `dv::derivative(fn, x, h)`. The reference text for `derivative_of` has the same unqualified form, although
the roots reference qualifies it correctly.

## Sentences I had to read twice, and other confusions

1. deriv/index.md line 59, "It finds a good step by itself and estimates its error, typically for 10 to 16
   evaluations of `fn`": "typically for" is awkward.
2. one-sided.md line 46: what "agree within their estimates" means, and what to do when the two do not agree.
3. one-sided.md line 50, "Avoid the one-sided stencils there: their step does not shrink near the edge": why would a
   step that stays inside the domain need to shrink? The reason, that the function's scale shrinks to the distance
   from the edge, is left implicit. Task 3 bears it out: `Backward<2>` has an error of 2e-4.
4. one-sided.md lines 31-32, "the central difference ... still fails at 1e-6": only on the second pass did I see that
   1e-6 is a distance from the edge, not a step.
5. basics.md lines 164-165, "`if (*d)` compiles, and tests whether the derivative is 0": it tests whether it is
   nonzero. Clear in meaning, inverted in wording.
6. choose-a-method.md lines 20-21: the rows "at the edge" and "near the edge" look nearly alike; "the one-sided
   Ridders method" does not say which one (chosen by the side the domain is on, which only the linked page explains).
7. step.md lines 5-6, "For the central difference it is eps^(1/3) ~ 6e-6 times that scale": holds for f' only; higher
   orders use eps^(1/(P+N)), which the chapter gives only as "The default step grows with N" (higher-derivatives.md
   line 31).
8. Cosmetic, in the source only: choose-a-method.md line 5 runs past the wrap width.

## Programs' outputs

```
== t1
value      2.5820042746129221
closed     2.5820042746129497  (mpmath 2.5820042746129492)
abs_error  1.39e-13
actual     2.75e-14  (vs mpmath 2.71e-14)
estimate covers error: true
== t2
value       1.5318575877405061
closed      1.5318575877451392
abs_error   6.75e-11  (relative 4.41e-11)
actual rel  3.02e-12
good to 1e-8 by the estimate: true; truly: true
Central<2>: 1.531857568770647, rel 1.24e-08, good_enough false
== t3
closed: -22.343905770087080 (mpmath -22.343905770087083)
fn = sqrt(1 - x*x):
  Ridders          -22.343905770082227  err 4.9e-12  est 1.580e-10  covers true  calls 30
  BackwardRidders  -22.343905770353427  err 2.7e-10  est 1.628e-10  covers false  calls 19
  Central<2>       -22.344008289025375  err 1.0e-04  no estimate     calls 2
  Backward<2>      -22.343703040947599  err 2.0e-04  no estimate     calls 3
fn = sqrt((1 - x)(1 + x)):
  Ridders          -22.343905770087535  err 4.5e-13  est 1.301e-12  covers true  calls 30
  BackwardRidders  -22.343905770087609  err 5.3e-13  est 5.667e-11  covers true  calls 18
== t4
closed            -778.800783071405
default step      h = 6.06e-06
  Central<2>      -778.776985530175  rel 3.1e-05
own step          h = 6.06e-09
  Central<2>      -778.800783039462  rel 4.1e-11
Ridders, default  -778.8007830714  rel 7.0e-15  est rel 3.3e-14
Ridders, a / 8    -778.8007830714  rel 6.4e-15  est rel 5.7e-14
== t5
Newton (Central<2>) 0.85260550201372542 +- 4.4e-16
Newton (Ridders)    0.85260550201372542 +- 5.6e-16
Brent [0, 1]        0.85260550201372554 +- 3.8e-16
|Newton - Brent| = 1.1e-16, within the bounds' sum: true
|Newton - W(2)|  = 1.1e-16, |Brent - W(2)| = 0.0e+00
== t6_snippet (derivative-of.md line 25, verbatim)
t6_snippet.cpp:13:43: error: 'derivative' was not declared in this scope; did you mean 'nxx::deriv::derivative'?
```

Measured noise of `sqrt(1 - x*x)` on 0.999 +- 2e-4: at most 1.58e-14 relative, about 71 epsilon.

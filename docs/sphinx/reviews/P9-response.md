# Phase 9, stage 6: the author's response to the decision letter

`docs-author`, 2026-11-10, on `P9-decision.md` (9a04d32), with the user's decisions U1 to U4 as recommended. Each item
of the letter, with what was done. U1 and U2 (`include/numerixx/deriv.hpp`, `docs/deriv/design.md`,
`tests/test_onesided.cpp`) are the main session's; the pages that touch that behaviour (`handle-failures.md`, the note on
`one-sided.md`) agree with the reworded comment and needed no change beyond the letter's items.

## For the user

| Id | Done |
|---|---|
| U1, U2 | The main session's. `handle-failures.md` already states the gap (7% off, estimate 11%, the flag not raised); unchanged on that point |
| U3 | `handle-failures.md`'s note lists "a formula that cancels, such as $1 - x^2$ near 1". `one-sided.md`'s advice gains the sentence, worded "or its rounding can make an estimate optimistic" rather than "both estimates": run, Ridders' estimate covers its error on $\sqrt{1 - x^2}$ at 0.999 (error 4.9e-12, estimate 1.6e-10) and only `BackwardRidders`' misses (error 2.7e-10, estimate 1.6e-10). `deriv_one_sided.cpp` asserts the miss (return 26) and the cure on $\sqrt{(1 - x)(1 + x)}$, both methods' estimates covering (27); reference $-x / \sqrt{(1 - x)(1 + x)}$ |
| U4 (a), (b) | `basics.md`: the overflow sentence moved after the fused multiply-add sentences it had split; "the Ridders methods" and "the stencils (fixed formulas such as `Central<2>`)"; `if (*d)` "tests whether the derivative is nonzero" (`Estimate` converts to its `value`, a `double`, then to `bool`) |
| U4 (c), (d) | Not this phase (the letter) |

## Merged findings

| Id | Done |
|---|---|
| G1 | `step.md`: "about four decades to the right, at smaller steps"; "the side where truncation error dominates". The plan's text is the editor's |
| G2 | The lambda is in `deriv_derivative_of.cpp` as its own block, `<own-step>`, rather than inside `<callable>`, so that the page explains it in its own paragraph and `<callable>`'s quoted output is unchanged: `double const h = 1.0 / 8.0;` and `[&](double x) { return dv::derivative<dv::Ridders>(fn, x, h); }`. `main` asserts its result at 2 within its estimate of $1 / (2\sqrt 2)$ (return 6) |
| G3 | `index.md` and `choose-a-method.md` as the letter words them. `deriv_costs.cpp` keeps each cell's numbers and asserts: `Central<4>` and `Central<6>` below `Central<2>` on $1/x$ and $\sqrt x$ (3); `Central<10>` above `Central<2>` on $1/x$ (4); `Central<4>`, `<6>`, `<10>` at least 0.5 off on the peak (5) |
| G4 | `step.md`: "`default_step`, called as `default_step<Method, N>(1.0)`, gives the default step for a function of scale 1; multiply it by your function's scale". `deriv_step.cpp` asserts `default_step<Central<2>>(1.0) * a == root_epsilon<3, double>() * a` (11) |
| G5 | `one-sided.md`, adapted: "Where `fn` changes fast near its edge, as $\sqrt x$ does near 0, avoid the one-sided stencils: their step does not follow `fn`'s scale, which shrinks to the distance from the edge, and they have no estimate to warn you. Where the derivative only steers Newton's method, `Forward<2>` is enough", with the roots link |
| G6 | `deriv_one_sided.cpp` now uses $\sqrt{1 - x}$ (NaN above 1), $f' = -1 / (2\sqrt{1 - x})$: Ridders' error below `BackwardRidders`' at 0.999 (13) and `Backward<2>` more than 10% off at $1 - 10^{-7}$ (14), both passing; the page's sentence stays, as "is the more accurate one" |
| G7 | Alt text as the letter; "two or three digits are left". `deriv_trade_off.cpp` asserts the last error between $10^{-3}$ and $10^{-2}$ (11) |
| G8 | `tutorial.md` as the letter |
| G9 | `deriv_kepler.cpp` asserts `default_step<Central<2>>(M)` between $6 \times 10^{-6}$ and $6.1 \times 10^{-6}$ (8) |
| G10 | The note moved after "Near the edge", worded as the letter, linking `nxx::deriv::ForwardRidders` and Handle failures; it now also holds the "not a proof" sentence (C9) |
| G11 | `one-sided.md`: "At the kink at 0, Ridders' method returns $-2$ ... with an estimate of only $1.2 \times 10^{-8}$: nothing tells it the kink is there" (the estimate is in the quoted output) |
| G12 | `one-sided.md`: the sum-of-estimates rule and what to do when they differ. `deriv_one_sided.cpp` asserts it on the $10^{-10}$ row (25): the two differ by less than the sum, the smaller estimate is `ForwardRidders`', and its error is the smaller |

## The critic

| Id | Done |
|---|---|
| C1 | G1; "truncation error dominates" on `step.md` and `index.md` |
| C2, C3 | G5, G2 |
| C4 | `handle-failures.md`: the stencil and the Ridders methods defined in the first cause. `one-sided.md`: "the stencils (fixed formulas) `Forward<P>` and `Backward<P>`, and the one-sided Ridders methods". The scale as "the distance over which `fn` changes noticeably" on `index.md`, `step.md` and `choose-a-method.md` |
| C5 | `choose-a-method.md`: reflowed; the fewest-evaluations clause cut; "Trying another method is a one-word change"; ", where it is defined" dropped |
| C6 | The preambles cut on `higher-derivatives.md`, `handle-failures.md`, `tutorial.md` and `step.md`. On `choose-a-method.md` the "This section is for ..." sentence becomes the page's job in one sentence ("When you need more digits, an error estimate, or a function defined on one side of $x$ only, choose another method"), since the rules ask the first paragraph to state it |
| C7 | `step.md`: the opening links the overview for the trade-off; "a step $h$ from $x$"; the note as the letter; G4 |
| C8 | `tutorial.md`: one sentence defines $M$ and $E$ and links the roots tutorial; the equation and the bracket stay; "a `std::optional` that Ridders' method always fills" cut |
| C9 | The domain becomes a clause ("their *domain*") rather than a sentence, so the term stays defined; "not a proof" in a note with the asserted case in a clause; the $10^{-10}$ numbers in the text, asserted (23: Ridders' error between $10^{-10}$ and $10^{-9}$, the one-sided method below $10^{-12}$ at 10; the estimate and 64 as before, 21); G11; G10 |
| C10 | The link to Choose a step; "today" cut; the finite-value sentence stays (rejected) |
| C11 | G7 |
| C12 | "Only a non-finite value of `fn` is detected; a finite but wrong value gives a wrong derivative, without an error"; "This manual calls the three *the Ridders methods*" as its own sentence, then "Here is Ridders' method on $\sqrt x$ at 2:"; G3. The scaffolding is U4(c) |
| C13 | "`Central<2>` for $f''$ uses three points, $x$ and $x \pm h$"; the two paragraphs merged into the letter's sentence |
| C14 | The last paragraph of `derivative-of.md` split into two |
| C15 | `costs.md`: "On the peak the higher-order stencils lose every digit; Ridders' method keeps 13 at 36 evaluations" (G3's check 5, and 6: Ridders' peak error below $10^{-12}$ at 36 evaluations); "The numbers hold for builds ..." |

## The docs-reviewer

| Id | Done |
|---|---|
| D6 | "It restarts below every step that reaches past the edge, the more often the closer $x$ is", with C9's numbers. `deriv_one_sided.cpp` counts the NaN evaluations of each Ridders run and asserts that there are some at every distance and more at each smaller one (24) |
| D8 | "Call such a copy through a non-`const` variable, as `copied` is" |
| D9 | The pragma's warnings named, as in the header |
| D10 | "`Central<2>` keeps about ..." |
| D11 | The cell links Costs |
| D15 | "the function, `fn`, as any callable" |
| D16 | "passes `derivative_of(fn)` to Newton's method" |
| D17 | "can outlive" |
| D18 | Every changed page reflowed to 120 columns (tables excepted) |

## The numerics-reviewer

| Id | Done |
|---|---|
| N11 | Measured first: the central difference in `cpp_bin_float_50` at 2 is 5.7e-35 relative. The bound becomes $10^{-32}$, the lower $10^{-40}$ unchanged, with the model and the measurement in the comment |

## The reader

| Id | Done |
|---|---|
| R1 | "at a cost of typically 10 to 16 evaluations" |
| R4 | "at $10^{-6}$ from the edge" |
| R6 | "`Ridders`, and `ForwardRidders` or `BackwardRidders` for the side the domain lies on" |
| R7 | "For the central difference and a first derivative it is ..." |
| task 2 | `higher-derivatives.md`: "Ridders' method starts from $\max(\lvert x \rvert, 1) / 4$ for these orders instead of $/32$". `deriv_higher.cpp` asserts `default_step<Ridders, 2>(2.0) == 0.5` and `default_step<Ridders, 1>(2.0) == 2 / 32` (10) |

## Budgets

Prose words by the editor's count (inline code and formulas one word each, headings counted, fenced blocks and notes
not), and with the notes counted:

| Page | Budget | Editor's count | With notes |
|---|---|---|---|
| index | 700 | 651 | 651 |
| tutorial | 700 | 505 | 505 |
| choose-a-method | 550 | 478 | 478 |
| step | 450 | 353 | 371 |
| one-sided | 500 | 446 | 499 |
| higher-derivatives | 400 | 337 | 337 |
| derivative-of | 300 | 269 | 269 |
| handle-failures | 500 | 388 | 459 |
| costs | 200 | 147 | 147 |
| book | 60 | 27 | 27 |

## Second round

On the "Second round" section of `P9-decision.md` (`P9-docs-reviewer-2.md`, `P9-numerics-2.md`). D2-12 and the
`include/`, `tests/` and `design.md` items are the main session's.

| Id | Done |
|---|---|
| D2-1 | `one-sided.md`: "if `fn` returns NaN (or an infinity) there" |
| D2-2 | "So return NaN outside the domain, or use a one-sided method." restored after "and gets $-1$" |
| D2-3 | `deriv_derivative_of.cpp`: the lambda captures by value, `[fn, h]`; its assertion unchanged. `derivative-of.md`, after the block: "For Ridders' method `h` is the initial step. The lambda holds copies of `fn` and `h`, so, like `derivative_of`, it can outlive them; with `[&]` it could not." |
| D2-4 | "...trust neither until you have checked at another initial step, the third argument of `derivative` ([Choose a step](step.md))" |
| D2-5 | `tutorial.md`: "off from the root solve by a few units in its last digit, which adds at most about $3 \times 10^{-10}$" |
| D2-6 | "cancels near an edge (subtracts nearly equal numbers and loses digits) ... or its noise makes an estimate too small"; the program's comment quotes the new wording |
| D2-7 | The first paragraph ends "One-sided methods evaluate `fn` on one side only."; "with an error estimate, `abs_error`"; the table's lead-in names "the relative error, the estimate relative to the derivative, and the evaluations" |
| D2-8 | The note: "on a polynomial, NaN below 0, whose zeros lie at the points Ridders' method tries after its restarts" |
| D2-9 | "`Forward<2>` or `Backward<2>` is enough" |
| N2-7 | `step.md`: "`default_step<Method, N>(1.0)` (in another real type, `T(1)`)"; the assertion stands |
| N2-8 | The note: "lose accuracy from their default step where `fn` changes on a scale below ..." |
| N2-9 | `deriv_one_sided.cpp`'s comment: up to about 125 eps relative in $1 - x^2$, about 60 eps after the square root; no assertion changed |
| N2-10 | `deriv_kepler.cpp`'s two comments say "below 1e-9" |

Budgets, by the editor's count (notes not counted; with them in parentheses), for the pages that changed: one-sided 500
of 500 (564), derivative-of 297 of 300, tutorial 510 of 700, step 358 of 450 (376).

## Codex

On `P9-codex.md` (e70ed13), the user's decisions of 2026-11-10: all three findings and the roots sentence adopted.

| Item | Done |
|---|---|
| 1 (major) | `one-sided.md`: the sentence that `Forward<2>` or `Backward<2>` "is enough" for Newton's method is replaced by "That holds for Newton's method too: on $\sqrt x - 10^{-6}$ from $2 \times 10^{-12}$, `Forward<2>` gives 525 for 353553, and Newton leaves the domain. Start it in a bracket, where a bad step falls back to bisection (the roots chapter)", following the sentence on one-sided stencils where `fn` changes fast near its edge. `deriv_one_sided.cpp` asserts: the slope at $2 \times 10^{-12}$ prints as 525 and the exact one as 353553 (28); Newton's first step lands below 0 (29); `solve<Newton>` with `derivative_of<Forward<2>>` from the guess returns `non_finite` (30); in the bracket $[0, 10^{-11}]$ it converges to $10^{-12}$ within `abs_error` plus $8 \times 10^{-28}$, justified in the comment (31). To stay within 500 words, the domain's example, the central difference's step and the roots link's text were shortened |
| 2 (minor) | `choose-a-method.md`: "A stencil gains less, since its default step scales with the type too: the central difference keeps about two thirds of the type's digits" (the program's check 7 is on `Central<2>`) |
| 3 (minor) | `deriv_kepler.cpp` checks against $E'(1) = 1.0373620218936458704990245805$ and $E''(1) = -0.556713032668587779463634931485$ (mpmath at 60 digits, as the comment says), each allowing half an ulp of its rounding to a double; the central difference's bound, Ridders' coverage, the gain and $f''$'s coverage are now against these. The page's closed forms stay, and a check (9) confirms that they agree with the constants within $10^{-14}$, justified by the solve's bound on $E$ |
| roots sentence | `roots/restricted-domain.md`: "From a guess, a one-sided difference can avoid the problem", and after the `Backward<2>` sentence: "It is not a cure. Its step does not shrink near the edge, so where `fn` changes fast there its slope can be far off, and Newton's step can leave the domain: on $\sqrt x - 10^{-6}$ from $2 \times 10^{-12}$ the solve fails with `non_finite`. A bracket safeguards Newton's method: a step that leaves it is replaced by bisection, and the same solve in $[0, 10^{-11}]$ converges." `roots_restricted_domain.cpp` asserts both (20, 21). Not "where `fn` varies slowly near the edge, as in its $\ln(x / 10^{-12})$ example": $\ln$ is steep at 0 too, and the forward slope there is also far off; Newton survives because its overshoot points away from the edge. So the page says what is true of both cases instead |

Words (headings and notes counted; inline code and formulas one word each): `one-sided.md` 500 of 500,
`choose-a-method.md` 481 of 550; `roots/restricted-domain.md`, outside this phase's budgets, 708 before and 771 after.

# Phase 9, stage 5: the editor's decisions on the reviews of the deriv chapter (d00357e)

`docs-editor`, 2026-11-10, on `P9-critic.md` (C), `P9-docs-reviewer.md` (D), `P9-numerics.md` (N) and
`P9-reader.md` (R). Duplicates are merged under one decision; each adopted finding says what the author changes and
what a program asserts. The reviewers agree on the chapter's shape and numbers; what they found is one false sentence
about figure 2, one untested line of code, one header comment the manual's own program contradicts, and wording.

**The user's decisions (2026-11-10): U1 to U4 all as recommended.** U1: the `ForwardRidders` comment and
`docs/deriv/design.md` reworded, and a test at 1e-11 added to `test_onesided.cpp` (closes the book's U-17). U2:
`derivative_of`'s comment qualified as `nxx::deriv::derivative<M, N>(fn, x, h)`. U3: the cancelling-formula clause
and sentence, asserted by `deriv_one_sided.cpp`. U4: the three one-line edits to `basics.md` in this phase; the
scaffolding and the roots note wait for the roots phase.

## For the user

**U1. The `ForwardRidders` comment, `deriv.hpp:306-307` (N2), with `docs/deriv/design.md:190-191`, which says the
same.** The comment promises that closer than $2^{-30} \max(|x|, 1)$ a wrong result is "flagged by an estimate of at
least |value| / 4". Evidence against it: `deriv_handle_failures.cpp` asserts $\sqrt x$ at $10^{-11}$ is 6 to 8% off
with an estimate below $|value| / 4$ (the estimate, 11%, does cover the error); the book's review recorded the same gap
on $\log x$ at $2^{-37}$ (error 0.224, estimate 0.193: neither covered nor flagged, by 0.4%) as D-18 and left it open
as U-17 for the library (`docs/deriv-book/PLAN.md` 12.3, 12.8); `tests/test_onesided.cpp:289-299` tests the flag at
$10^{-12}$ only. Recommendation: change both sentences to "closer, the levels run out and the result is wrong:
usually flagged by an estimate of at least |value| / 4, otherwise within a few times its estimate (sqrt at 1e-11
from 0: 7% off, estimate 11%); pass a step on fn's own scale", and add a test in `test_onesided.cpp` at $10^{-11}$
asserting the error above 5%, $4\,$`abs_error` below $|value|$ and the error within twice the estimate, so that the
new sentence is pinned as the old one was. This closes U-17. `handle-failures.md` already states the gap correctly
and needs no change.

**U2. The `derivative_of` comment, `deriv.hpp:1893` (D2, R).** The lambda it suggests, `[&](double x) { return
derivative<M, N>(fn, x, h); }`, does not compile in user code: `derivative` is unqualified and a `double` brings no
namespace for ADL. Recommendation: `nxx::deriv::derivative<M, N>(fn, x, h)`. No test: the manual's program will hold
the compiled lambda (A2 below), and the comment is prose.

**U3. A formula that cancels is a noisy `fn` (R, task 3).** The reader's $\sqrt{1 - x^2}$ at 0.999 carries a relative
error of up to 71 $\varepsilon$ from $1 - x^2$, and `BackwardRidders`' estimate (1.63e-10) missed its error (2.7e-10);
written as $\sqrt{(1 - x)(1 + x)}$ both methods' estimates cover. No page says that a closed formula can be a "less
accurate `fn`". Recommendation, two places: the note in `handle-failures.md` lists "a solve to a loose tolerance or
interpolated data" and gains "or a formula that cancels, such as $1 - x^2$ near 1"; `one-sided.md`'s advice gains one
sentence, since functions cancel near their edges: "A formula often cancels near an edge, as $1 - x^2$ does near 1:
write it as $(1 - x)(1 + x)$, or its noise makes both estimates optimistic." `deriv_one_sided.cpp` asserts both: on
$\sqrt{1 - x^2}$ at 0.999 `BackwardRidders`' estimate is below its error, and on $\sqrt{(1 - x)(1 + x)}$ both
estimates cover (references: $-x / \sqrt{1 - x^2}$). About 40 words; `one-sided.md` stays within 500 with C9's cuts.

**U4. Outside this chapter.** (a) `basics.md:230-236` (D): the overflow sentence was inserted between a sentence on
fused multiply-add and its continuation, so "do not by default" now reads as being about `nxx::quadrature`; (b)
`basics.md:162-165` (D, R5): "Ridders' methods" and "the plain finite differences" for the Ridders methods and the
stencils, and "`if (*d)` ... tests whether the derivative is 0" is inverted (it is true when the derivative is not 0).
Recommendation: three one-line edits by the main session now, rather than waiting for the shared pages' phase; they
are wording, and the reader met (b). (c) The Tutorial, Guide and Reference headings with a sentence each over the
toctrees (C12): about 45 words on each of the three overviews, and the PDF needs the headings; leave until the roots
and quadrature phases and decide for all three together. (d) `roots/derivatives.md` shows `derivative_of` with Newton's
method in a bracket only; the reader had to guess the guess form (R, task 5): a note for the roots phase.

## Merged findings (three or more reviews, or two with the same fix)

| Id | Finding | Decision | Author |
|---|---|---|---|
| G1 | `step.md:47` says the peak's curve lies "to the left"; the axis puts smaller steps to the right (C1, D1, N1) | adopt | "but the peak's lies about four decades to the right, at smaller steps"; the plan's figure 2 text corrected the same way (editor) |
| G2 | `derivative-of.md:25`: a code line in prose, unqualified, without `<Method, N>`, which does not compile (C3, D2, R) | adopt | The lambda goes into `deriv_derivative_of.cpp`'s `<callable>` block as `[&](double x) { return dv::derivative<dv::Ridders>(fn, x, h); }` with a step `h` of 1/8, shown through `literalinclude`; `main` asserts its result at 2 within its estimate of $1 / (2\sqrt 2)$ |
| G3 | `index.md:54-55`: "`Central<4>` and `Central<6>` are more accurate" is unasserted, implies 8 and 10 are not, and "smooth" alone is false for the peak (D4, N4, N3, C12d) | adopt | `index.md`: "`Central<4>`, `<6>`, `<8>` and `<10>` use more points; on a function smooth on the scale $\max(\lvert x \rvert, 1)$, `Central<4>` and `Central<6>` are the more accurate, and the highest orders can be worse ([Costs](costs.md))". `choose-a-method.md:24-25`: "...so it can be less accurate even on a function that varies slowly (`Central<10>` on $1/x$, [Costs](costs.md)), and on one that changes fast it can lose every digit". `deriv_costs.cpp` asserts, from its cells: on $1/x$ and $\sqrt x$ `Central<4>` and `Central<6>` below `Central<2>`; `Central<10>` above `Central<2>` on $1/x$; `Central<4>`, `<6>`, `<10>` at or above 0.5 on the peak |
| G4 | `step.md:65-66`: "`default_step` ... which you can scale the same way" gives a step $\max(\lvert x \rvert, 1)$ times too large for $\lvert x \rvert > 1$ (D3, N7, C7d) | adopt | "`default_step<Method, N>(1.0)` gives the default step for a function of scale 1; multiply it by your function's scale." `deriv_step.cpp` asserts `default_step<Central<2>>(1.0) * a` equals the page's `root_epsilon<3, double>() * a` (an API fact, not a numerical reference) |
| G5 | `one-sided.md:50`: "Avoid the one-sided stencils there: their step does not shrink near the edge" reads as general advice; `Forward<2>` works at the edge and at $10^{-10}$ inside it, and the linked roots section recommends it for Newton's method; the reason, that `fn`'s scale shrinks to the distance from the edge, is implicit (C2, N6, D7, R3) | adapt | "Where `fn` changes fast near its edge, as $\sqrt x$ does near 0, a one-sided stencil's step does not follow `fn`'s scale, which shrinks to the distance from the edge, and it has no estimate to warn you; use the Ridders methods there. Where the derivative only steers Newton's method, `Forward<2>` is enough ([A numerical derivative near the edge](../roots/restricted-domain.md#a-numerical-derivative-near-the-edge))." The `Backward<2>` case at $1 - 10^{-7}$ stays asserted (G6) |
| G6 | `one-sided.md:43` names $\sqrt{1 - x}$, the program asserts $\sqrt{1 - x} + x$ (D12, N5); D ran $\sqrt{1 - x}$ itself at 0.999: Ridders 2.3e-14 against `BackwardRidders` 1.9e-13 | adopt | `deriv_one_sided.cpp` asserts on $\sqrt{1 - x}$ (NaN above 1): Ridders' error below `BackwardRidders`' at 0.999, and `Backward<2>` more than 10% off at $1 - 10^{-7}$, rechecking the bound on that function; the page's sentence stays |
| G7 | Figure 1's alt text says the straight line reaches 4.1e-12 at 5e-6, where the point is already rounding-dominated; "fewer than four digits" at 2.5e-3 is loose (C11, N10) | adopt | Alt: "falls in a straight line to about 5e-12 at 1e-5, reaches its lowest, 4.1e-12, at 5e-6, then rises raggedly to 2.5e-3 at 1e-13"; text: "two or three digits are left". `deriv_trade_off.cpp` asserts the last error between $10^{-3}$ and $10^{-2}$ |
| G8 | `tutorial.md:67-68`: "far below $10^{-9}$" is asserted only as below; the $3 \times 10^{-10}$ is the solve's error, not rounding (D13, N8) | adopt | "The relative error is below $10^{-9}$. Each $E(M)$ is a few ulps off from the root solve, which adds at most about $3 \times 10^{-10}$." No program change |
| G9 | `tutorial.md:40`: "about $6 \times 10^{-6}$" is not asserted by `deriv_kepler.cpp` (D14, N9b) | adopt | `deriv_kepler.cpp` asserts `default_step<Central<2>>(M)` between $6 \times 10^{-6}$ and $6.1 \times 10^{-6}$ |
| G10 | The $10^{-9} \max(\lvert x \rvert, 1)$ note on `one-sided.md` is asserted nowhere, sits after the kink section, and is about near-edge behaviour (N9, C9e) | adopt | Move it after "Near the edge", as: "The one-sided Ridders methods lose accuracy where `fn` changes on a scale below about $10^{-9} \max(\lvert x \rvert, 1)$ ({cpp:struct}`nxx::deriv::ForwardRidders`), as $\sqrt x$ does near 0 ([Handle failures](handle-failures.md#a-result-that-is-not-accurate-enough))"; the number is the header's $2^{-30}$, linked |
| G11 | `one-sided.md:71`: "Ridders' method sees a kink" contradicts the header (a kink at $x$ is invisible to it); the small estimate's meaning is unsaid (D5, C9d) | adopt | "At the kink at 0, Ridders' method returns $-2$, the mean of the two slopes, with an estimate of only $1.2 \times 10^{-8}$: nothing tells it the kink is there." The estimate is in the quoted output |
| G12 | `one-sided.md:46`: "agree within their estimates" is undefined (the reader's sum reading gave the right answer; the each-estimate reading would not), and the page does not say what to do when they disagree (R2, task 3) | adopt | "Keep the one with the smaller estimate if the two differ by less than the sum of their estimates; if they differ by more, trust neither until you have checked at another step." `deriv_one_sided.cpp` asserts the rule on the near-edge table's $10^{-10}$ row (the two agree, the smaller estimate is the one-sided method's and its error is the smaller) |

## The critic (C)

| Id | Finding | Decision | Author |
|---|---|---|---|
| C1 | figure 2 direction; "truncation rules" idiom | adopt | G1; "the side where truncation error dominates" on `step.md` and `index.md` |
| C2 | one-sided stencils advice | adapt | G5 |
| C3 | the lambda in prose | adopt | G2 |
| C4 | "stencil" and "the Ridders methods" used undefined on `one-sided.md` and `handle-failures.md`; "scale" defined only by its formula | adopt | `handle-failures.md:14-15`: "A stencil, a fixed formula such as `Central<2>`, fails ... The Ridders methods (`Ridders`, `ForwardRidders`, `BackwardRidders`) first restart..."; the same clauses at first use on `one-sided.md`; on `index.md`, `step.md` and `choose-a-method.md` the scale is "the distance over which `fn` changes noticeably, assumed to be about $\max(\lvert x \rvert, 1)$" |
| C5 | `choose-a-method.md:3-11`: the 150-column line; "the fewest evaluations of an expensive `fn`" as a reason the default is not enough; the repeated "called the same way"; ", where it is defined" | adopt | Reflow; cut the clause; "Trying another method is a one-word change"; drop ", where it is defined" (the editor's stage-4 addition) |
| C6 | Preambles that restate the headings (`higher-derivatives.md:5-6`, `handle-failures.md:4-5`, `tutorial.md:4-6`, `choose-a-method.md:5`, `step.md:6-7`) | adopt | Cut each; the opening sentences already state the page's job, as the rule requires |
| C7 | `step.md` repeats the overview (lines 3-5) and `handle-failures.md` (the note); "a step $h$ apart" is wrong for a central stencil; the `default_step` sentence | adopt | "The default step balances truncation against rounding error ([Numerical derivatives](index.md#the-idea-and-its-catch)) under one assumption: ..."; "at points a step $h$ from $x$"; the note becomes "A step too small for the points to differ, or not positive and finite, gives `non_finite` ([Handle failures](handle-failures.md))"; G4 |
| C8 | `tutorial.md:11-18` repeats the roots tutorial's anomalies; line 82 repeats the overview | adapt | One sentence defines $M$ and $E$ (readers land anywhere), the equation and the bracket stay, the rest links to the roots tutorial; "a `std::optional` that Ridders' method always fills" cut |
| C9 | `one-sided.md`: the domain definition; the "not a proof" sentence; name the error in the table's text; the kink's estimate; the note's place; "above" and "below" | adopt, except the last | Cut "The part where a function is defined is its *domain*"; "An estimate is not a proof ..." becomes a `{note}` naming the asserted case in a clause; the table's text names the error, estimate and evaluations at $10^{-10}$ (asserted: error between $10^{-10}$ and $10^{-9}$, estimate $3 \times 10^{-6}$ to $3 \times 10^{-5}$, 64; the one-sided method below $10^{-12}$ at 10); G11; G10. Rejected: "above $x$" and "below $x$" compare values, and the manual's rule is about position in the document |
| C10 | `handle-failures.md`: link "a step on the function's own scale" to Choose a step; cut "today"; "A wrong but finite value is not detected at all" repeats the overview | adopt, adopt, reject | The link; "today" cut. The sentence stays: the page on failures must say what is not a failure, and it carries the link |
| C11 | figure 1's alt text and "fewer than four digits" | adopt | G7 |
| C12 | `index.md`: "Only values that are not finite are detected" ambiguous; "they are *the Ridders methods*" inside the code lead-in; the scaffolding; link Costs for the orders | adopt, adapt, user, adopt | "Only a non-finite value of `fn` is detected; a finite but wrong value gives a wrong derivative, without an error." The definition becomes its own sentence before the lead-in: "This manual calls the three *the Ridders methods*." U4(c). G3 |
| C13 | `higher-derivatives.md`: the second-difference formula; merge lines 41-43 | adapt, adopt | "`Central<2>` for $f''$ uses three points, $x$ and $x \pm h$" (the formula is theory the reader does not act on); "Use Ridders' method above the first derivative unless evaluations are expensive: it keeps far more digits and gives an estimate ({cpp:struct}`nxx::deriv::Ridders` gives the counts)" |
| C14 | `derivative-of.md:48-50`: two ideas in 45 words | adopt | Split as proposed |
| C15 | `costs.md`: the peak's `1e+00` cells unexplained though `choose-a-method.md` sends the reader here; "It holds" | adopt | "On the peak the higher-order stencils lose every digit; Ridders' method keeps 13 at 36 evaluations." Asserted in G3's additions, plus Ridders' peak error below $10^{-12}$ and 36 evaluations; "The numbers hold for builds..." |

## The docs-reviewer (D)

| Id | Finding | Decision | Author |
|---|---|---|---|
| D1-D5, D7, D12-D14 | merged | | G1, G2, G4, G3, G11, G5, G6, G8, G9 |
| D6 | "Closer, it restarts" implies no restart at $10^{-2}$; Ridders restarts at every distance (2, 15 and 29 evaluations at NaN points) | adopt | "It restarts below every step that reaches past the edge, and the closer $x$ is, the more often: its error and estimate grow, and at $10^{-10}$ it makes all 64 evaluations it may." Joined with C9's numbers |
| D8 | A callable copied from a non-const `fn` cannot be called through a `const` result; the page writes `auto copied` without saying why | adopt | One clause after the state example: "call such a copy through a non-`const` variable, as `copied` is" |
| D9 | The pragma's warning is not named | adopt | "`-Wdeprecated-declarations` for GCC and Clang, C4996 for MSVC" (the header's names) |
| D10 | "the central difference keeps about 10 digits for $f'$, 8 for $f''$": the term is `Central<2>` for $f'$ only | adopt | "`Central<2>` keeps about..." |
| D11 | "Typically 13 to 15 digits" unlinked | adopt | Link the cell to [Costs](costs.md), whose rows pin it |
| D15 | `fn` first used without introduction on `index.md` | adopt | Line 5: "the function, `fn`, as any callable" |
| D16 | `derivative-of.md:26` "passes such a function" follows the own-step lambda | adopt | "passes `derivative_of(fn)` to Newton's method" |
| D17 | "the result outlives `fn`" | adopt | "can outlive" |
| D18 | Lines over 120 columns (`choose-a-method.md:5`, `index.md:48`) | adopt | Reflow |
| out of scope | `basics.md:230-236`, `basics.md:162-163` | user | U4 |

## The numerics-reviewer (N)

| Id | Finding | Decision | Author |
|---|---|---|---|
| N1, N3-N10 | merged | | G1, G3, G3, G6, G5, G4, G8, G10 and G9, G7 |
| N2 | The `ForwardRidders` header comment | user | U1 |
| N11 | Three loose assertions: `cpp_bin_float_50` error below $10^{-30}$ for "two thirds of the digits" (model $3 \times 10^{-34}$); Kepler's 8 to 20 calls for "about a dozen"; 5 to 7 digits for "5 or 6" | adopt, reject, reject | `deriv_choose_a_method.cpp`: the bound becomes $10^{-32}$ (the author checks the measured error first, and keeps the model in the comment). The Kepler range is deliberately wide for the platforms' `std::sin`, as its comment says; the digits are pinned by the printed output |

## The reader (R)

| Id | Finding | Decision | Author |
|---|---|---|---|
| task 3 | A cancelling formula is a noisy `fn`; the "agree" rule | user; adopt | U3; G12 |
| snippet | the lambda | adopt | G2 |
| R1 | "typically for 10 to 16 evaluations" is awkward | adopt | "at a cost of typically 10 to 16 evaluations" |
| R2, R3 | "agree within their estimates"; why a step must shrink | adopt | G12; G5 |
| R4 | "still fails at $10^{-6}$" read as a step | adopt | "at $10^{-6}$ from the edge" |
| R5 | `basics.md`'s `if (*d)` sentence inverted | user | U4(b) |
| R6 | The "at the edge" and "near the edge" rows look alike; which one-sided method | adopt | The near-edge row: "`Ridders` and `ForwardRidders` or `BackwardRidders`, for the side the domain lies on" |
| R7 | `step.md:5-6`'s $\varepsilon^{1/3}$ holds for $f'$ only | adopt | "For the central difference and a first derivative it is..." |
| task 2 | Ridders' starting step for $f''$ to $f''''$ ($\max(\lvert x \rvert, 1) / 4$) is only in the reference | adopt | `higher-derivatives.md`, after "The default step grows with $N$": "Ridders' method starts from $\max(\lvert x \rvert, 1) / 4$ for these orders instead of $/32$." `deriv_higher.cpp` asserts `default_step<Ridders, 2>(2.0)` is 0.5 |
| task 5 | The guess form of Newton's method with `derivative_of` is shown nowhere in roots | user | U4(d) |
| R8 | the long line | adopt | D18 |

## Counts

Items after merging: 50. Adopted 34 (G1 to G4 and G6 to G12; C4 to C7, C10, C12's wording and link, C14, C15; D6,
D8 to D11, D15 to D18; N11's multiprecision bound; R1, R4, R6, R7 and the reader's task-2 gap). Adapted 5 (G5, C8,
C9's note, C12's definition, C13's formula). Rejected 3 (C9's "above" and "below", C10's sentence on finite values,
N11's two ranges). For the user 8, in four groups: U1 (N2), U2 (the header part of D2 and the reader's snippet), U3
(the reader's task 3), U4 (the two `basics.md` items, R5, the scaffolding of C12, the reader's task-5 gap).

## The author's changes in brief

Pages: `index.md` (G3, C4, C12, D15, R1, the idiom); `tutorial.md` (G8, C6, C8); `choose-a-method.md` (G3, C4, C5,
D11, R6, D18); `step.md` (G1, G4, C4, C6, C7, R7); `one-sided.md` (G5, G6, G10, G11, G12, C4, C9, D6, R4, and U3 if
approved); `higher-derivatives.md` (C6, C13, D9, D10, task 2); `derivative-of.md` (G2, C14, D8, D16, D17);
`handle-failures.md` (C4, C6, C10, and U3 if approved); `costs.md` (C15). Programs: `deriv_derivative_of.cpp` (G2),
`deriv_costs.cpp` (G3, C15), `deriv_step.cpp` (G4), `deriv_one_sided.cpp` (G6, G12, C9, and U3), `deriv_trade_off.cpp`
(G7), `deriv_kepler.cpp` (G9), `deriv_higher.cpp` (task 2), `deriv_choose_a_method.cpp` (N11). Every page stays within
its budget: the cuts (C5 to C10, C13) outweigh the additions (G5, G12, C4, U3). The editor corrects the plan's figure 2
text (G1). The second, targeted round: `docs-reviewer` and `numerics-reviewer` on `one-sided.md` and its program,
which carry most of the changes, and on the header comments if the user approves U1 and U2.

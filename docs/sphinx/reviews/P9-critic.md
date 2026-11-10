# Phase 9, stage 5: docs-critic on the deriv chapter (commit d00357e)

`docs-critic` (Sonnet, high), 2026-11-10. Saved by the main session as the agent returned it.

Method: I read `docs/sphinx/CLAUDE.md`, `PHASE9_PLAN.md`, `reviews/P9-plan-critic.md`, all ten pages, the two plot
scripts, the programs' tables and outputs, `roots/index.md`, `roots/costs.md` and `quadrature/index.md` and
`quadrature/costs.md`. I could not run a word counter, so the counts are hand counts of prose. They sit within about
5% of the editor's, except `handle-failures`, where I count more.

**Verdict.** The chapter is in good shape. Every page is within budget and the plan's cuts were made. What is left is
two wrong or misleading statements, one rule violation, a few terms used before they are defined, and about 300 words
of preamble and repetition. The pages overall come to about 3,650 words against a planned 4,360, so the budgets are no
longer the issue.

## Word counts against budget

| Page | Budget | Count (mine, editor's) | Verdict |
|---|---|---|---|
| index | 700 | ~620 (622) | within |
| tutorial | 700 | ~590 (573) | within |
| choose-a-method | 550 | ~490 (486) | within |
| step | 450 | ~365 (358) | within |
| one-sided | 500 | ~470 (442) | within, but the least slack |
| higher-derivatives | 400 | ~345 (352) | within |
| derivative-of | 300 | ~265 (261) | within |
| handle-failures | 500 | ~420 (397) | within |
| costs | 200 | ~135 | within |
| book | 60 | 27 | within |

## Findings, most valuable first

**1. `step.md` line 47: the text says the opposite of the figure.** "the peak's lies about four decades to the left"
is wrong. Smaller steps are to the right (`plots/deriv_scale.py` line 53), and the alt text itself says "four decades
further right". The plan's description has the same slip.
- Rewrite: "but the peak's lies at steps about four decades smaller".
- Also at line 49, "the side where truncation rules" becomes "the side where truncation error dominates". The same
  idiom is in `index.md` line 28.

**2. `one-sided.md` line 50: "Avoid the one-sided stencils there: their step does not shrink near the edge, and they
have no estimate to warn you."** It misleads.
- Counter-case: at 1e-10 inside the edge, `Forward<2>` evaluates only at x and above, so every point is inside the
  domain. It works, and the page itself shows `Forward<2>` working at the edge, with about 10 digits.
- A stencil that points at the edge fails loudly with `non_finite`, not silently.
- What is true is that a stencil has no estimate to take part in the "try both, keep the smaller estimate" rule.
- Rewrite: "A stencil has no estimate to compare, so use the Ridders methods near an edge." Ask `numerics-reviewer` to
  confirm.

**3. `derivative-of.md` line 25: untested code in prose.** The rule is that code appears only through
`literalinclude` from a tested program. No program contains `[&](double x) { return derivative(fn, x, h); }`.
- Fix, shortest: "For a step of your own, pass a lambda that calls `derivative(fn, x, h)`."
- The alternative is to add the lambda to `// <callable>` with an assertion.

**4. Terms used before they are defined on the page.**
- `one-sided.md` line 25 uses "the stencil `Forward<2>`" and line 50 "one-sided stencils" without defining stencil.
- `handle-failures.md` lines 14-15 and 44 use "A stencil fails..." and "the stencils leave it empty". Lines 15, 42 and
  62 use "the Ridders methods" and line 58 `ForwardRidders`, none defined on the page.
- Fix for `handle-failures` line 14: "A stencil, a fixed formula such as `Central<2>`, fails as soon as...". At line
  15: "The Ridders methods (`Ridders`, `ForwardRidders`, `BackwardRidders`) first restart...". Add a similar clause on
  `one-sided`.
- "Scale" is defined only by its formula (`index.md` line 32, `choose-a-method.md` line 4, `step.md` line 5). "A
  function that varies on the scale max(|x|, 1)" means nothing to a reader who lands there. Proposed clause: "...a
  function that changes noticeably only over distances of about max(|x|, 1), its *scale*".

**5. `choose-a-method.md` lines 3-11.**
- Line 5 is a reflow accident, about 150 characters against the 120 of the other lines.
- "or the fewest evaluations of an expensive `fn`" is listed as a reason the default is not enough. The default is
  the cheapest method, and the table's first row says so. Cut it.
- "Every method is called the same way, `derivative<Method>(fn, x)`, so trying another is a one-word change" repeats
  `index.md` line 78. Cut to "Trying another method is a one-word change."
- Table row 4, "x at the edge of `fn`'s domain, where it is defined": drop ", where it is defined".
- The plan's Stage 4 says this page defines "domain" at first use; it does not. It does not need to (finding 9).

**6. Preambles that restate headings (about 75 words).**
- `higher-derivatives.md` lines 5-6: "This section shows how to ask for them and what accuracy to expect."
- `handle-failures.md` lines 4-5: "This section shows how to recognise each case, and what to try next."
- `tutorial.md` lines 4-6: the roadmap sentence repeats the headings; keep the roots-tutorial sentence.
- `choose-a-method.md` line 5: "This section is for when that is not enough:" can open the list directly.
- `step.md` lines 6-7: "This section is about functions that break the assumption, and what to do then."
- Each page's first sentences already state its job. Delete these sentences or fold them into the opening.

**7. `step.md` repeats the overview and `handle-failures.md` (about 55 words).**
- Lines 3-5 restate `index.md` lines 13-16. Replace with "The default step balances truncation against rounding
  error ([Numerical derivatives](index.md#the-idea-and-its-catch)) under one assumption".
- Line 3, "points a step h apart", is wrong for a central stencil, whose points are 2h apart. Say "at points a step h
  from x".
- The `{note}` at lines 101-104, on `non_finite` for a step of 1e-300 or not positive, repeats `handle-failures.md`
  bullets 2 and 3. Replace with one sentence and a link, or cut it.
- Lines 65-66, "which you can scale the same way", is vague. Say: "`default_step<Method, N>(1.0)` gives the default
  step for a function of scale 1; multiply it by your function's scale."

**8. `tutorial.md`: repeats the roots tutorial and the overview (about 40 words).**
- Lines 11-18 explain mean and eccentric anomaly in about 95 words; the roots tutorial does the same. Keep the
  equation and the root bracket, and link to the roots tutorial.
- Line 82, "a `std::optional` that Ridders' method always fills", repeats `index.md` line 84.

**9. `one-sided.md` (about 45 words saved, 20 added).**
- Line 4, "The part where a function is defined is its *domain*": the reader knows calculus. Cut. Keep "edge".
- Lines 46-48, "An estimate is not a proof: ... it can even be 0 for a wrong result." An unillustrated edge case in the
  main line; `index.md` already says an estimate is "not a bound". Cut it, or make it a `{note}`.
- Lines 41-44, the text for the table, should name the error as well as the estimate: "Closer, it restarts below each
  step that met NaN, so its cost, error and estimate all grow: at 1e-10 it makes all 64 evaluations it may, with an
  error of 2e-10 and an estimate of 8e-6. `ForwardRidders` keeps 1e-13 and 10 evaluations."
- Lines 71-72: add that the small estimate does not warn you: "...returns -2, the mean of the two slopes, with an
  estimate of only 1.2e-8: nothing warns you."
- The `{note}` at lines 74-78 follows the kink section but is about near-edge and singular behaviour. Move it after
  "Near the edge" and simplify: "The one-sided Ridders methods lose accuracy where `fn` changes on a scale below about
  1e-9 max(|x|, 1), as sqrt(x) does near 0 (link to Handle failures)." The placement is only a suggestion.
- Line 25, "`Forward<P>` and `ForwardRidders` above, `Backward<P>` and `BackwardRidders` below": "above" and "below"
  collide with the manual's "never below" rule for position. Say "at larger x ... at smaller x".

**10. `handle-failures.md`.**
- Line 59, "With a step on the function's own scale, x/8": link it to [Choose a step](step.md).
- Line 12, "The enumeration has one value today": cut "today". Lines 37-38, "A wrong but finite value is not detected
  at all (link)", repeats `index.md` lines 87-88. Keep the link and drop the clause, or cut it.

**11. Figures.**
- Figure 1 (`deriv_trade_off`) earns its place.
- Alt and text fixes: "falls in a straight line ... to 4.1e-12 at 5e-6" overstates, since the table has 5.1e-12 at
  1e-5 and 4.09e-12 at 5e-6, and the last point is no longer on the truncation line. Say "to about 5e-12 at 1e-5, its
  lowest 4.1e-12 at 5e-6". "fewer than four digits are left" at 2.5e-3 should be "two or three digits left".
- Figure 2 (`deriv_scale`) earns its place as the one picture of "scale". Besides finding 1, "The two curves have one
  shape" hides the peak's plateau at error 1; harmless, since the plan dropped it deliberately.
- Alt text for both gives numbers that match the tables. The near-edge table earns its place (see finding 9).

**12. `index.md` small items (about 25 words).**
- "Only values that are not finite are detected" is ambiguous. Say "Only a non-finite value of `fn` is detected; a
  finite but wrong value gives a wrong derivative, without an error."
- Lines 60-61, "with `Ridders`, they are *the Ridders methods*:", puts a definition into the sentence that leads into a
  code block. Drop it and name the three methods in line 87, or move the definition.
- The Tutorial, Guide and Reference headings with their sentences run to about 45 words for toctrees. Low value;
  roots and quadrature have the same layout, so only worth doing if they change together.
- Lines 53-56: "`Central<4>` and `Central<6>` are more accurate" invites the question why not 8 and 10. Add a link to
  Costs.

**13. `higher-derivatives.md` (about 28 words).**
- Line 26, the classic second difference formula, is not needed to choose a step or read a result. Cut, or reduce to
  "uses three points".
- Lines 41-43: merge into "Use Ridders' method above the first derivative unless evaluations are expensive; it keeps
  far more digits and gives an estimate (counts: link)".

**14. `derivative-of.md` lines 48-50, the last sentence.** About 45 words and two ideas. Split: "A lambda that
captures its counter by reference needs no `std::ref`: the copy shares the counter. The captured object must still be
alive when you call the derivative."

**15. `costs.md`.** The cell "`1e+00`" for `Central<4>` on the peak is unexplained, and `choose-a-method.md` line 25
sends the reader here for "a case". Add: "On the peak, the higher-order stencils lose every digit (error 1); Ridders'
method keeps 13 at 36 evaluations." "It holds for builds that do not fuse..." should say "The numbers hold for
builds...".

## Voice against roots and quadrature

The deriv chapter is now terser, more second person and more colon-led than `roots/index.md` and
`quadrature/index.md`, which are still long; expected before their phases. To align then: "an estimate, not a bound"
(deriv) against "a bound on its error" (roots): the libraries really differ, so keep it. "The Ridders methods" and
"stencil" are chapter terms and need not cross modules.

## Keep as it is

- The overview's one idea, truncation against rounding, with figure 1 and the default-step formula.
- Figure 2, its caption and its alt text, once the "left" slip is fixed.
- The near-edge table and the "try both" advice.
- `choose-a-method.md`: the table, the counting lambda, the "higher order is not always better" reason, the
  multiprecision example.
- `step.md`: the code blocks and outputs, the 36-against-14 evaluations contrast, the scaled-step example.
- `higher-derivatives.md`: the table, the "divides rounding by h^N" reason and the `float` section.
- `derivative-of.md`: the callable example and the `std::ref` state example.
- `handle-failures.md`: the causes list, the 7% and 70% case with the "do not rely on the flag alone" advice, the
  near-0 and noisy `fn` note.
- `book.md` at 27 words; `costs.md` as a pure reference with the pointer to the header's counts.
- The tutorial's structure, the second derivative as call plus check, and the "Next" list.

Net: about 330 words cut and about 60 added, leaving the chapter near 3,400 prose words.

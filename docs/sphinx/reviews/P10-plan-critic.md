# Phase 10, stage 1: docs-critic on the plan

`docs-critic` (Sonnet, high), 2026-11-11, on `docs/sphinx/PHASE10_PLAN.md` as the editor wrote it. Saved by the main
session as the agent returned it; line numbers refer to that version of the plan and to the current pages. The critic
had no shell: its counts are hand estimates, about 5-10% above the plan's counter, used only to set targets.

**Verdict.** The shape is right: no merges, one job per page, the programs kept, and the Phase 9 lessons that cost
nothing learned (two figures instead of four, `rt` into the programs, costs at 200, a terms list, a "considered and not
proposed" list). Not the main one: the budgets are too loose, the counter is not auditable, one figure does not earn
its place, and the one merge it proposes contradicts itself. On the plan's budgets the chapter ends near 6,250 words;
proposed about 5,400.

## Findings, most valuable first

**1. Budgets larger than the job needs (plan lines 3-4, 32-45; saves about 850 words).** The plan cuts 20% (7,823 to
about 6,250). Phase 9 planned 25% and delivered about 37%; the deriv pages came in at 80-90% of their budgets. Overviews:
the rule is about 800 and deriv's delivered 622.

| Page | Now | Plan | Proposed |
|---|---|---|---|
| index | 981 | 750 | 650 |
| tutorial | 876 | 700 | 620 |
| choose-a-method | 734 | 600 | 520 |
| tolerances | 771 | 600 | 520 |
| find-a-bracket | 511 | 450 | 380 |
| restricted-domain | 772 | 600 | 520 |
| derivatives | 568 | 500 | 420 |
| other-types | 528 | 450 | 350 |
| handle-failures | 697 | 600 | 560 |
| method-of-your-own | 1065 | 800 | 650 |
| costs | 293 | 200 | 210 |
| total | 7,796 | 6,250 | about 5,400 |

`costs.md` goes up by ten: about 110 of its 200 words are table cells and the Steffensen footnote. Do not cut the
tables.

**2. The counter is not in the repository (plan lines 28-30).** "Counted by one script ... the editor's, not in the
repository" recreates the Phase 9 problem the counting rule exists to prevent. Commit the script (for example
`docs/sphinx/count_words.py`) before stage 3 and name it in `docs/sphinx/CLAUDE.md`.

**3. Drop figure 1, `roots_closing_in` (plan lines 139-152).** It does not replace words: the "52 and 9" sentence stays
as the figure's text, and `costs.md`, `choose-a-method.md` and `tolerances.md` already carry the fact. Brent's width
goes from 5e-3 to 4.8e-3 across evaluations 5 and 6, a flat step the reader will ask about, which pulls in Brent's
internals (theory cap). The optional third series (Q3) shows the plan is unsure of the figure's fact. Its program would
assert counts `roots_overview.cpp` and `roots_costs.cpp` already assert. If the user wants a picture on the overview,
figure 2 is the one.

**4. Keep figure 2, `roots_cancellation`, but make it replace words and trim its checks (plan lines 153-165).** The
scatter shows what a table cannot: noise filling a window around the root against a smooth cubic, and it carries the
chapter's second idea (`abs_error` bounds `fn` as computed). It must replace the paragraph at `index.md` lines 88-94
(about 130 words): a caption, two "what to see" sentences and the Brent result. Keep the `<cubic>` block and output.
Keep three assertions (more than 50 sign changes in the window written out, exactly one for the product, the farthest
wrong sign near 7e-6); `roots_overview.cpp` asserts the rest. Use 201 points, not 401. Alt text with the asserted
counts only.

**5. `method-of-your-own.md`: 650, not 800 (plan lines 43, 47-50, 92-98).** The observer is about 90 words of the 265 to
be cut. This page serves the smallest audience; quadrature's 924 is not a precedent. Where the words go: the opening
and "Two kinds of method" into one paragraph of about 55; the Pegasus history to the program's comment ("Dowell and
Jarratt (1972)" is in `PHASE6_API.md` and the spike files, not in `docs/roots/design.md`: check it against the
Citations rule, or drop the name); the "four habits" paragraph to one sentence; "What `solve` checks" to three
sentences and the link; the `inside = 0` paragraph from about 120 to 70; "An open method" from about 230 to 170 (keep
`values = 3` and the `min_step = 0` failure). Q2: do not split. Keep the observer's move.

**6. `restricted-domain.md`: the merge contradicts itself, so do not merge (plan lines 38-39, 75-80).** The plan merges
"under the second heading" and keeps the anchor `#a-numerical-derivative-near-the-edge`, linked from deriv's
`one-sided.md`, `handle-failures.md` line 23 and `derivatives.md` line 80; the two sections are different ideas (an open
method steps out; a stencil's points leave the domain). Keep three `##` headings and cut inside them. A larger option
for the user: drop the `<one-sided>` block and its output (about 120 words) for one sentence and a link to
`deriv/one-sided.md` ("a one-sided stencil works only where `fn` is tame near the edge"), keeping the sqrt(x)
counter-case; otherwise keep it and budget 560.

**7. `other-types.md`: drop the expression-templates block, 350 words (plan lines 41, 85-87).** Lines 49-61 repeat
`derivatives.md` lines 34-37 and `basics.md#multiprecision-types`: state it once, in `derivatives.md`, with one link.
Cut the compile-time sentence. Q4: drop the MSVC clause; keep "a `float` `fn` on a `double` bracket compiles but is
accurate only to `float`". Keep the multiprecision output.

**8. `tolerances.md`: replace the "`abs_tol` acts at every magnitude" section with one sentence (plan lines 37,
70-72).** "`abs_tol` applies at every magnitude, so choose it from the scale of the roots you expect, not as an accuracy
setting." Drop "Brent's does better than its guarantee": it invites reliance. If the user wants the demonstration,
keep the block and cut the prose after it.

**9. Cross-page repetition (about 150 words).** "A root where `fn` touches zero without crossing it" with x^2 at 0 on
four pages: define once in `choose-a-method.md` ("Multiple roots"), link from the others. Newton's tangent on three
pages: keep the overview's sentence. The `int` pitfall in full on `tutorial.md` and `derivatives.md`: keep it in
`derivatives.md`. "Every method is called the same way" on `index.md` and `choose-a-method.md`: keep it in the
overview. The stop test on `tolerances.md` and `method-of-your-own.md`: link. "Bracket" defined on nine pages: keep it
as a clause of at most 12 words.

**10. `choose-a-method.md`: 520.** Merge table rows 4 and 5 (both Bisection). Open straight with the three things the
choice depends on (no "this section is for..." preamble). The note on the teaching methods from about 85 to 35 words.

**11. `handle-failures.md`: keep the observer, sized (saves about 40).** Hold the new section to about 60 words and the
first six steps. Cheaper alternative if tight: assert the cycle without quoting output. Cut the preamble (about 110
words: `basics.md`'s).

**12. `index.md` to 650 and `tutorial.md` to 620.** Index: the opening from about 130 to 55; the open-method paragraph
from about 135 to 90 (no tangent explanation); `search` two sentences and a link; the pole and several-roots paragraph
to a clause and the link. Tutorial: the roadmap to one sentence; "Checking the result" from about 270 to 170 (cut the
residual-bound arithmetic and the platform sentences; keep the two-doubles reference).

**13. `find-a-bracket.md` to 380.** Cut "tedious"; "Search then solve" prose to 30 words; drop the "a fiftieth, growing
by half" detail with the second example; the touches-zero bullet per finding 9.

**14. `derivatives.md` to 420.** Finding 9's cuts save about 70; keep the Halley sentence; accept the guess form (P9 U4
(d)).

**15. Q1, the Tutorial / Guide / Reference headings (about 90 words over three overviews).** Keep the headings (the PDF
needs them); drop the sentences, since the toctree titles say the same; keep one only where it adds ("Start with
[Choose a method]"). Check whether a toctree `:caption:` gives the PDF the names without headings.

**16. The plan's own text.** "What the chapter does well" to the three lines of what goes; figure descriptions to
what to see, the data and the one fact (P9-plan-critic finding 13 again). Keep "Evidence gaps".

**17. Q5, "Ridders'" in two chapters.** The recommendation is fine.

## Theory cap

No theory over the cap beyond Pegasus's history, the tutorial's residual-bound arithmetic, and holding Brent's "square
of bisection's count" and Toms748's "4 per halving" in both `choose-a-method`'s table and its paragraph (the row is
enough).

## Keep as it is

The page structure; the overview's two ideas, the `<brent>` block, figure 2 with the `<cubic>` block; the terms list
("bound", never "estimate"; the teaching methods; Ridders' linked at first use; no "below"); `rt` into the includes;
the tutorial's structure; `choose-a-method`'s table, counting lambda and "Multiple roots"; `tolerances`' near-zero
section (1049 against 41), `rel_tol` (52 against 21) and the minimum-defaults note; `find-a-bracket`'s blocks and the
kelvin example; `restricted-domain`'s NaN example, the Newton-leaves block, the `stencil` block and the sqrt(x)
counter-case; `handle-failures`' error list, `or_else` and the pole test; `method-of-your-own`'s Pegasus and Chebyshev
blocks, `values = 3`, the two convergence rules, the 33-against-10 stall; `costs.md` as a pure reference; the
programs table; the reader's tasks; the "Considered and not proposed" list. Open questions: Q1 headings kept,
sentences shortened; Q2 no split, 650; Q3 no, and drop figure 1; Q4 yes; Q5 yes.

## Main session's note

Finding 6's "contradicts itself" does not hold: the plan's "the second heading" is the second of the two it names, "A
numerical derivative near the edge", whose anchor deriv links to, so the anchor would stay. Whether to merge remains a
judgment for the user.

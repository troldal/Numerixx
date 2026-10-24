# R2-06-ridders: response to review round 2

Author: `book-author-methods` (Fable 5.1), stage 7. File: `chapters/06-ridders.tex` (the editor's one-line change to
the comment of `figures/06-ridders-coverage.tex`, T4, is kept). Compiled with `tools/check_chapter.sh`
(`COMPILED OK`, 8 pages standalone, unchanged; Figure 6.3 alone on its page, 6.4 and 6.5 on the following pages);
`tools/check.py` reports no warning. Net change: +4 lines, as the letter's plan allows.

| Finding | Decision | What changed |
|---|---|---|
| T6 (evidence-01) | modified | 6.1: "Ridders' paper supplies the idea; the rules, the ratio 2, a rounding floor, a two-level stall, a final check and restarts, are the library's own, chosen on a prototype (Sidebar 6.2), and every statistic in this chapter is of them". The paper's "factor 1.4" and "stops when the error grows" are gone; Numerical Recipes is not cited |
| numerics-04 | accept | 6.3, the scope sentence: "which the method survives by restarting, or, where the samples stay finite, by taking more levels, Sidebar 6.2" |
| T2 (codex-01) | modified | 6.4: the paragraph is now "A jump or a kink inside the span". The large-estimate statement is restricted to a jump in f ("Where f itself jumps, it also returns a large estimate ... The estimate flags a jump; the stencil does not"), and one sentence adds the kink: "A kink, a jump in f' with f continuous, at x itself shows no symptom at all: in exact arithmetic the central difference of c\|t\| is zero at every step, so every rule at every step returns the mean of the two one-sided slopes plus the smooth part's own error, Ridders' tableau converges to that mean as it would for a smooth function, with the estimate at the floor, and nothing in the samples can show it; know the branch boundaries and differentiate the branch you mean (7.8)". No data set and no counterexample numbers, as decided |
| T10 (numerics-07a, evidence-08) | accept | 6.4: "(derived from the initial step and the level cap; for f' also measured here, Section 10.4)"; the ten-digit distance reads `2^{-\valRiddersEdgeTenDigitsExponent} s` instead of the typed 2^-30 |
| evidence-03 | accept | 6.5: "returns about three more digits than Central<2> (a median of `\valExpCentralTwoRelError` on e^x, about `\valHeadlineFirstDigits` digits, against Ridders' `\valExpRiddersFirstRelError`, about `\valHeadlineRiddersFirstDigits`), as many as Central<4> (...), none over Central<6> (...)"; 6.6: "about three more digits than Central<2> (`\valHeadlineRiddersFirstDigits` against `\valHeadlineFirstDigits`) and as many as Central<4>". 7.3 and 10.3 say the same (practice author) |
| T4 (numerics-05, evidence-05) | accept | Sidebar 6.2: "occurred in `\valRiddersZeroErrorFirst` cases for f', which lie on `\valRiddersZeroFunctionsFirst` functions, `\valRiddersZeroOnQuinticFirst` of them on x^5 and the rest values that happen to equal the rounded closed form, and in `\valRiddersZeroErrorSecond` for f'', all on x^5" (x^5 is counted among the seven, as the macro does) |
| numerics-03 | modified | Sidebar 6.2: "a smaller initial step usually restores the margin: rerun at a quarter of the default initial step, `\valRiddersMarginQuarterCovered` of the study's `\valRiddersMarginCases` uncovered cases over all four orders are covered, and at a sixty-fourth `\valRiddersMarginSixtyFourthCovered`" (the editor's DS9 extension: 11 of 14 and 14 of 14; the sidebar says "over all four orders" because the 14 include the one uncovered f''' case) |

Nothing else in the chapter was objected to in round 2.

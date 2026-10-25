# P5-07-using: the Phase 5 update (PLAN.md 12.8)

Author: `book-author-practice` (Opus 5.5). File: `chapters/07-using.tex`. `tools/check_chapter.sh`: `COMPILED OK`,
13 pages standalone (14 before: about 12 lines more of text, but the page breaks after Table 7.1 moved and the last
page filled); `tools/check.py`: 0 errors, 0 warnings. Examples `code/07-using-*.cpp` unchanged; all ten build and
pass (`ctest --preset gcc -R numerixx_book_07`), and none quotes a statement Phase 5 changed.

| Item of 12.8 | What changed |
|---|---|
| 7.1, the call | The method list names `ForwardRidders` and `BackwardRidders`; "for Ridders' method and its one-sided forms, the initial step" |
| 7.2, the result | `abs_error` "filled only by Ridders' method and its one-sided forms" |
| Table 7.1, "f defined on one side" | Adds "with an error estimate, `ForwardRidders` or `BackwardRidders`": `\valForwardRiddersEvalsMedianFirst` evaluations, at most `\valForwardRiddersMaxEvalsFirst`; `\valExpForwardRiddersFirstRelError`, covering in `\valForwardRiddersCoverageFirst` % of the study; evidence also `sec:ridders-onesided` |
| Table 7.1, float | "never float" now also in the second-derivative row, and the caption says "Above f', float warns" (the third/fourth row keeps "never float") |
| Table 7.1, kink row | "or with an estimate `ForwardRidders` or `BackwardRidders`"; "any method across the break returns a wrong value" |
| 7.8, kink cure | "above a break with `Forward<P>`, or `ForwardRidders` for an error estimate; below it with `Backward<P>` or `BackwardRidders` ... At the break itself the one-sided methods return that side's derivative, provided f(x), which they sample too, belongs to the branch" (header, the `derivative` comment on kinks; `test_onesided.cpp` "a kink at x") |
| 7.8, a poor Ridders value | Cure: "`ForwardRidders` or `BackwardRidders` at an edge" |
| 7.8, domain edge | The nearest bad point is "a NaN or a pole" (the edge limit now holds for poles); "Exactly at an edge, use `ForwardRidders` or `BackwardRidders` for a value with an error estimate (`sec:ridders-onesided`), `Forward<P>` or `Backward<P>` inside an iteration" |
| 7.8, float | Renamed "Derivatives above the first in float": `\valFloatSecondCentralRelError` for f'' and `\valFloatFourthCentralRelError` for f'''' with `Central<2>`; "Every such call, by every method, warns at compile time" |

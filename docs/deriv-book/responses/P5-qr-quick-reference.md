# P5-qr-quick-reference: the Phase 5 update (PLAN.md 12.8)

Author: `book-author-practice` (Opus 5.5). File: `chapters/qr-quick-reference.tex`. `tools/check_chapter.sh`:
`COMPILED OK`, 1 page (the right column keeps about a sixth of the page free); `tools/check.py`: 0 errors, 0
warnings. The one-sided row of the method summary is the editor's shared file (P5-8) and appears unchanged.

| Item of 12.8 | What changed |
|---|---|
| The result | `*d->abs_error` is "the error estimate of the three Ridders methods" |
| The warning | Failures: "A float derivative above the first warns at compile time, by every method" (was: a float fourth derivative) |
| The methods list | New item: "`ForwardRidders`, `BackwardRidders`: Ridders on one side of x only, for an edge, with an error estimate; a median of `\valForwardRiddersEvalsMedianFirst` evaluations for f', at most 33, 65, 97 and 129 (macros); about `\valExpForwardRiddersFirstRelError` for f'" |
| The edge bullet | Ridders' item: "Near a bad point (a NaN or a pole)"; the pitfall "Evaluable is not valid" adds "or `ForwardRidders` or `BackwardRidders` for an error estimate" (its wording of style decision 16 unchanged) |

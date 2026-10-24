# R2-07-using: response to review round 2

Author: `book-author-practice` (Opus 5.5), stage 7. Files: `chapters/07-using.tex`, `code/07-using-expected.cpp`
(E11). E11 builds with GCC 16 and Clang 22 under the book's `-Werror` flags and prints `true false` as before;
`ctest --preset gcc -R numerixx_book` passes 24 of 24. `tools/check_chapter.sh`: `COMPILED OK`, 14 pages standalone
with 5 lines on the last, as before the revision (658 printed lines against 655); `tools/check.py`: no warning.

| Finding | Decision | What changed |
|---|---|---|
| T1 (engineer-01, developer-01, numerics-01, library-01) | accept / modified | `ridders-check`: `accept(res, atol, rtol)` is `res && res->abs_error && *res->abs_error <= atol + rtol * std::abs(res->value)`; the calls pass `0.0, 1e-8`; E11's `main` also asserts that a stencil result fails the check. 7.2 introduces it as "compare its estimate with the error you can accept, in two tolerances" and follows it with the guide: rtol = 10^-n for n significant figures; atol, in the derivative's units, decides near a zero of f' (7.5); no check inside an iteration (7.4.1); the check says the estimate is small, not that the value is right where the estimate is blind, as at a kink (7.8) |
| T2 (codex-01) | modified | 7.5: "A jump inside the span does not survive this check ...; a kink at x does: every step returns the mean of the two slopes (6.4)". 7.8's kink pitfall: "a kink at x itself shows no symptom at all, since every rule returns the mean of the two slopes; Cure: know the branch boundaries, and differentiate the branch you mean". Table 7.1's last row: "only Ridders' estimate shows a jump, nothing shows a kink at x" |
| codex-04, cline-14 | accept | Figure 7.2's caption: "a plateau of J/h times the sum of the weights on one side of the break: J/(2h) for Central<2> (`\valKinkInsideCentralRelError`, h = `\valKinkCentralStep` K); 7J/(12h) for Central<4> while the break lies within h of x (`\valKinkInsideCentralFourRelError`, h = `\valKinkCentralFourStep` K) and -J/(12h) between h and 2h, the two levels the figure shows. So it falls as 1/h across the rules ... only Ridders' estimate (dashed) shows a jump" |
| evidence-04, cline-12 | accept | 7.4.2: "m x n calls each costing two model evaluations with Central<2>: `\valJacobianBalanceCallEvals` model evaluations here, against `\valJacobianBalanceHandEvals` by hand" |
| evidence-03 | accept | 7.3: Central<4>'s "about two more digits" replaced by the medians (`\valExpCentralTwoRelError`, `\valExpCentralFourRelError`); Ridders "about three more digits than Central<2> (14 against 11, by the macros), as many as Central<4>, none over Central<6>" |
| library-02 | accept | 7.9: "Let the lambda return the inner value, or a NaN when the inner call fails (`d ? d->value : std::numeric_limits<double>::quiet_NaN()`), so that the outer call reports non_finite" |
| library-04 | accept | 7.1: "a double step is diagnosed under -Wconversion unless float holds it exactly" |
| cline-13 | accept | 7.1's margin note: "a pointer to a const member function of T itself" |
| cline-16 | modified | 7.6: "the library converts a result at once, so an expression over the parameter alone is safe, but one over the lambda's own locals would refer to destroyed values (Sidebar 7.1); the explicit type is the habit that is always safe" |
| cline-17 | accept | 7.2: `r.value_or(0.0)` and `r->abs_error.value_or(0.0)` |
| cline-18 | modified | 7.5: "is the stencil's error to within Ridders' own error, `\valCpRiddersActual` J/(mol K) here, below its estimate; the stencil's actual error is ..." |
| cline-19 | accept | 7.1: "a step of another type converts to it" |
| Cline remark (Table 7.1, row 4) | accept | Evaluations "2, or Ridders' `\valRiddersEvalsMedianFirst`" |
| Cuts named in section 3 | done | Sidebar 7.1's second paragraph from nine lines to six; Table 7.1's caption by a line; 7.4.3's fixed-point paragraph from seven lines to four (the two remedies kept, the iteration-count remark dropped); the two trailing comment lines of `estimate-pitfalls` (the `.value()` trap stays in the prose); 7.9's mixed-partials and deferred-features paragraphs by a line each (the caching sentence shortened, the list condensed). Also: the `ridders-check` comment by a line, the guide and the "rule of higher order" paragraph tightened, to hold the chapter's length |

Not changed: the editor's `\needspace` commands are kept as they are (one before `check-by-steps` leaves a gap on the
standalone build, not necessarily in the book); the "Did it work" paragraph and the rest of the chapter are as in
round 1.

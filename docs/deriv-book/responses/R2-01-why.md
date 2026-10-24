# R2-01-why: responses of the fundamentals author

Chapter 1, `chapters/01-why.tex`, revised against `reviews/R2-decision.md` (section 3, row `01-why.tex`, +2 lines).
Wrapped PDF: 6 pages before and after; the additions fall in the slack of the last page.

| Finding | Decision | What changed |
|---|---|---|
| T1 (engineer-01, developer-01, numerics-01, library-01; codex-01 in part) | accept | 1.5, "Is it right?": the relative-only test is replaced by the two-tolerance check of the decision, `res && res->abs_error && *res->abs_error <= atol + rtol * std::abs(res->value)`, so that a stencil result (no estimate) fails it rather than dereferencing an empty `optional`. One-sentence guide: `rtol` = 10^-n for a value reported to n significant figures; `atol` is the error you can accept in the derivative's units, and it decides near a zero of f'; inside an iteration no check is needed. Then: "A passed check says that the estimate is small, not that the value is right where the estimate is blind, as at a kink" (T2's sentence) |
| engineer-03 | (editor) | Table 1.1 follows the two checks; the anchor and its comment are the editor's and were left as they are |

# R2-qr-quick-reference: response to review round 2

Author: `book-author-practice` (Opus 5.5), stage 7. File: `chapters/qr-quick-reference.tex`. `tools/check_chapter.sh`:
`COMPILED OK`, one page.

| Finding | Decision | What changed |
|---|---|---|
| T1 | accept | "Two checks": "For Ridders, also require `r->abs_error && *r->abs_error <= atol + rtol * std::abs(r->value)`: rtol = 10^-n for n figures, atol in the derivative's units decides near f' = 0; a stencil has no estimate to test. Inside an iteration no check is needed" |
| T2 | modified | The pitfall: "A kink or a range break within a step of x gives a wrong value from every method, without failing; Ridders' estimate flags a jump, nothing flags a kink at x: differentiate the branch you mean" |
| developer-03 | accept | "Your own step": "(Central<2>, f'; in general `root_epsilon<P + N, T>()`)" |
| Cuts named in section 3 | done | "Where to look" cut to two lines; the Ridders bullet loses its last clause ("(for f')") |

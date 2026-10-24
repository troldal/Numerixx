# R2-a1-floating-point: responses of the fundamentals author

Appendix A, `chapters/a1-floating-point.tex`, revised against `reviews/R2-decision.md` (section 3, row
`a1-floating-point.tex`, +2 lines). Wrapped PDF: 2 pages before and after.

| Finding | Decision | What changed |
|---|---|---|
| codex-06 | accept | A.4 states the two resolution rules: "If neither x + h nor x - h moves, the library reports `error::non_finite`. If only one moves, as can happen near a power of two, a two-node rule for f' returns a finite value that is mostly rounding, and every other rule reports `non_finite` (Sidebar 3.3, Section 8.2)." This matches `deriv.hpp`: `difference()` rejects only coinciding points; the multi-node branch checks both sides |
| cline-04 | accept | A.4: `fl(x + h)` "differs from x + h by up to half the spacing at the result: ulp(x)/2, or ulp(x) when x + h crosses a power of two (the two-sided spacing of Sidebar 3.3)" |

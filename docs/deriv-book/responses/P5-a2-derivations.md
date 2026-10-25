# P5-a2-derivations: the Phase 5 update (PLAN 12.8)

Author: `book-author-methods`, Phase 5 update. File: `chapters/a2-derivations.tex`. `tools/check_chapter.sh`:
`COMPILED OK`, 4 pages standalone (unchanged).

| Item | What changed |
|---|---|
| B.3, the forward tableau's divisor | New paragraph at the end of B.3: the `Forward<1>` rule of order N has an expansion in every power of h; the induction of B.3 with `c h_i^j` and the parent's `c 2^j h_i^j` gives the divisor 2^j - 1 and one order per column; unrolled, the coefficients of an entry alternate in sign, so the inherited rounding is at most `sum |gamma_k| r_k`, which is what the carried floor's recursion computes (referred to from 6.7) |
| PLAN 1.1's first cut (B.2'

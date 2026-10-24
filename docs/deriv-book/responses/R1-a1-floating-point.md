# R1-a1-floating-point: responses of the fundamentals author

Appendix A, `chapters/a1-floating-point.tex`, revised against `reviews/R1-decision.md` (section 3, row
`a1-floating-point.tex`). Wrapped PDF: 2 pages before and after (budget 0); A.3 and A.4 were tightened to make room.

| Finding | Decision | What changed |
|---|---|---|
| library-01 | accept | A.4: a step below ulp(x)/2 does not move x; if neither x + h nor x - h moves, the library reports `non_finite`; if only one moves, as near a power of two, the result is finite but mostly rounding (Sidebar 3.3); hence the precondition, a step well above ulp(x). The old "no difference can be formed, which the library reports as non_finite" for any step below ulp/2 is gone |
| numerics-17f | accept | A.5: Sterbenz's lemma is stated for two numbers of the same sign ("if a, b > 0 and a/2 <= b <= 2a, then fl(a - b) = a - b; likewise for a, b < 0"), cited as `\cite[Thm~4.3.1]{sterbenz1974}` per the reference table |
| developer-01 | accept | A.6 lists the `numeric_limits` members the module uses: `is_specialized`, `is_integer`, `is_bounded` and `is_exact`, which the concept `nxx::real` checks, and `epsilon()`, `max()`, `radix`, `digits` and `denorm_min()`, which the algorithms use (pointer to 8.7, where E23 is quoted). The gap between the two lists is U-5, for the author |

Also tightened, with no change of content: A.3's overflow sentence ("the library checks what could overflow, before it
is formed; the final quotient after").

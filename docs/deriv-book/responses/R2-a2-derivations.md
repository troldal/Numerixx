# R2-a2-derivations: response to review round 2

Author: `book-author-methods` (Fable 5.1), stage 7. File: `chapters/a2-derivations.tex`. Compiled with
`tools/check_chapter.sh` (`COMPILED OK`, 4 pages standalone, unchanged); `tools/check.py` reports no warning.
Net change: +1 line, as the letter's plan allows.

| Finding | Decision | What changed |
|---|---|---|
| T8 (numerics-06, codex-08) | accept | B.1: "since each xi_k tends to x with h and f^(N+P) is continuous, the remainder is M_(N+P) h^P f^(N+P)(x)/(N+P)! + o(h^P) for every rule, and the correction is O(h^(P+1)) when f has one more bounded derivative, which the error model assumes throughout". The appendix's opening assumption names that derivative: "f has as many continuous derivatives as the expansions use, and one more bounded derivative where the error model's O(h^(P+1)) is claimed (B.1)" |
| codex-02, cline-01 | accept | B.4: E(h*) = (1 + N/P) S eps \|f\| / h*^N (the truncation term is N/P times the rounding term at the minimum; 3R/2 for N = 1, P = 2, as Chapter 3 has it). The cost paragraph's fractions N/(P+N) and P/(P+N) of E(h*) were already consistent with this and are unchanged |
| cline-02 | accept | B.4: "the largest of the four rules the design record quotes" replaces "the largest of the library's rules" (Central<10> for f'''' would cost 18 by the same formula, which nothing quotes) |
| cline-03 | accept | B.4: "the 'at most a factor 1.64' of the library's documentation for f' is the cost at kappa rounded to 1.8 (Central<10>, where the unrounded kappa gives `\valStepCostCentralTenFirst`)" replaces "is a rounding of kappa = 1.8, where the formula gives 1.63" |

Nothing else in the appendix was objected to in round 2.

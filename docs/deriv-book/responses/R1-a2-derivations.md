# R1-a2-derivations: response to review round 1

Author: `book-author-methods` (Fable 5.1), stage 6. File: `chapters/a2-derivations.tex`. Compiled with
`tools/check_chapter.sh` (`COMPILED OK`, 4 pages standalone, as before); `tools/check.py` reports no warning in
the file.

| Finding | Decision | What changed |
|---|---|---|
| T8 (numerics-04) | accept | B.1 is restructured. The Taylor argument now ends with the asymptotic remainder `M_(N+P) h^P f^(N+P)(x)/(N+P)! + O(h^(P+1))`, stated to hold for every rule and to be what the error model uses. A new paragraph gives the single-point form: by the intermediate-value theorem where the products `w_k n_k^(N+P)` have one sign, which is the case for Eq. (5.1) (products 1, 0, 1) and not for Eq. (4.3) (products -8/3, 2/3, 2/3, -8/3), nor for any central rule for f' above order 2, nor for some one-sided rules; and for the first derivative by the interpolation error formula (Burden et al., Sect. 4.1), with x among the nodes (added with weight 0 where the rule omits it), differentiated once at the node x where the node product vanishes, which for Central<4> gives the `+ h^4 f^(5)(xi)/30` of Eq. (4.3) (the product over -2, -1, 1, 2 is 4, m = 5). The one-sided rules are covered by the same formula; for the higher derivatives with mixed signs the book uses only the asymptotic form. Chapter 4's sidenote at Eq. (4.3) now points here |
| T8 (numerics-01) | accept | B.1's reference to Eq. (4.3) is consistent with the corrected sign |
| library-17, evidence-21 | accept | B.4: "These agree with the design record's figures of 3.5, 11 and 15 (`docs/PHASE4_API.md`); the 1.64 that the library's documentation gives for Central<10> is a rounding of kappa = `\valStepConstantCentralTenFirst`, where the formula gives `\valStepCostCentralTenFirst`" replaces "These are the figures the header and the design record state, rounded". The half-step result is attributed to "the design record's prototype (`docs/PHASE4_API.md`)" |
| Citations (evidence reference table) | accept | `\cite[Sect.~1.1]{burden2016}` at Taylor's theorem, `\cite[Sect.~4.1]{burden2016}` at the interpolation error formula and the Lagrange formula, `\cite[Sect.~4.2]{burden2016}` at the tableau's induction, `\cite[ch.~23]{chapra2015}` at Chapra's row |

Not changed, and why: B.2, B.3 and the formulas of B.4 were checked by the numerics reviewer and reported correct;
only the citation forms changed there.

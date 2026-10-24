# R1-09-tables: response to review round 1

Author: `book-author-practice` (Opus 5.5), stage 6. File: `chapters/09-tables.tex`. `tools/check_chapter.sh`:
`COMPILED OK`, 4 pages standalone, as before. `tools/check.py`: no warning in the file.

| Finding | Decision | What changed |
|---|---|---|
| T10 (numerics-07, evidence-05) | modified | 9.2: "The cap at P = 10 is a design decision (docs/PHASE2_API.md): by the model's exponent alone, a Central<12> would gain about `\valCapGainCentralTwelveFirstDouble` digit for f' in double and `\valCapGainCentralTwelveFirstMultiprecision` in a 50-digit type (for f'''', `\valCapGainCentralTwelveFourthDouble` and `\valCapGainCentralTwelveFourthMultiprecision`), over a stencil on which f must be smooth that much wider, and the library stops there". 9.3: "The cap at P = 6 is a design decision (docs/PHASE2_API.md): beyond it the amplification grows about as fast as the higher order pays, by the model, so the library stops there" (qualitative, as decided; the "which is why" phrasing is gone) |
| library-20 | accept | 9.5: "The level cap and the initial-step divisors are read from the library; the ratio, the stall and the two factors of the untrusted-result test are literals of its code, transcribed by the generator". The introduction says the same in a parenthesis. Table 9.10's untrusted-result row is the generator's new text ("estimate at least |A|/4 and more than 64 times the floor") |
| evidence-25 (citation form) | accept | The introduction cites "the published tables of Fornberg\cite{fornberg1988} and Chapra.\cite[ch.~23]{chapra2015}" |

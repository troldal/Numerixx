# P5-04-stencils: the Phase 5 update (PLAN 12.8)

Author: `book-author-methods`, Phase 5 update. File: `chapters/04-stencils.tex`. `tools/check_chapter.sh`:
`COMPILED OK`, 6 pages standalone (unchanged).

| Item (PLAN 12.8) | What changed |
|---|---|
| 4.4, Ridders at an edge (lines ~263-267) | Kept the measured failure of `Ridders` at 500 K (`\valEdgeRiddersStatus`, `\valEdgeRiddersEvals`) and added: "Its one-sided form `ForwardRidders` samples at 500 K and above only, and returns the derivative there with an error estimate (6.7)" |
| 4.5, choosing P | "For a value to report with an error estimate, use Ridders, or at an edge its one-sided forms" |

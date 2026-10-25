# P5-06-ridders: the Phase 5 update (PLAN 12.8)

Author: `book-author-methods`, Phase 5 update. Files: `chapters/06-ridders.tex`, the new
`figures/06-ridders-onesided-reach.tex` (F19, DS23) and the new `code/06-ridders-onesided.cpp` (E24). Compiled with
`tools/check_chapter.sh` (`COMPILED OK`, 9 pages standalone against 8 at release; the last page is full, so any
further line in this chapter adds a page); `tools/check.py`: 0 errors, 0 warnings; `ctest --preset gcc -R
numerixx_book`: 25 of 25 pass.

| Item (PLAN 12.8) | What changed |
|---|---|
| New section (P5-2, F19, E24) | 6.7 "One-sided Ridders", `sec:ridders-onesided`, after 6.6 and before the sidebars (so its margin figure is far from the full-width F12 and F14, style decision 19). Uses: the one-sided situations (domain edge, end of a correlation's range, right- or left-hand slope at a kink); the `Forward<1>` base of order N, fn at x and above only, restarts included; `BackwardRidders` as `ForwardRidders` on f(-t) at -x; the expansion in every power of h, divisor 2^j - 1, one order per level (pointing to B.3); the floor carried through the tableau, `r_{i,j} = (2^j r_{i,j-1} + r_{i-1,j-1}) / (2^j - 1)`; everything else as Ridders; initial step `s / \valForwardRiddersDivisor` for every order; cost fn(x) once and N per level, at most `\valForwardRiddersMaxEvalsFirst` to `...Fourth`, typically `\valForwardRiddersEvalsMedianFirst` and `...Fourth`; the study's medians, accuracy factors, coverage and worst ratio (First and Fourth, the extremes) with the study's scope; when to prefer Ridders; E24 with `\valEdgeRiddersStatus`, `\valEdgeRiddersEvals` and the `\valEdgeForwardRidders...` macros; the reach on log x at 2^-k (ten digits, the errors at k = 30, covered to, flagged from, the silent case `\valForwardRiddersReachSilentExponentFirst` stated as measured, never "always flagged", D-18); the rule of thumb 2^-30 s or a step on f's scale |
| F19 | `figures/06-ridders-onesided-reach.tex`: relative error (solid) and estimate (dashed) for f' (ink) and f'''' (accent), central Ridders' f' error (muted, dotted), a line at `\valForwardRiddersReachCapExponentFirst` (all evaluations spent) and a circle at the silent case, both from macros. The caption says "further out, to the last k measured, the result is flagged ... At one distance for f' (the circle) it is neither covered nor flagged" |
| E24 | `code/06-ridders-onesided.cpp`, snippet `onesided-edge` (5 lines, the calls only, E6's wrapper above the markers and named in the prose); main asserts Ridders `non_finite`, each one-sided actual error within 6 times its estimate, relative errors at most 1e-10 and 1e-7. Output matches the macros (2.1e-13, estimate 2.9e-12 relative, 9 evaluations; 4.5e-9, 3.4e-7, 17; Ridders non_finite after 64) |
| Intro | Last sentence names the one-sided forms |
| 6.4, edge (12.8: "exactly at an edge ... the one-sided rules of 4.4 are the tool") | "exactly at an edge it does not, and its one-sided forms (6.7) or the one-sided stencils of 4.4 are the tool" |
| 6.4, poles (12.8: restart rule; edge limit now valid for poles) | "restarts below the smallest step it tried when it fires"; new sentence: the restart descends as after a non-finite value, so the edge limit holds for a pole too, with `\valPoleCovered` of `\valPoleCases` covered, `\valPoleRelError{Min,Max}`, `\valPoleEvals{Min,Max}`, and f'' at the two poles (`\valPoleAheadRiddersSecond{RelError,Evals}`, `\valPoleBehindRiddersSecond{RelError,Evals}`) |
| Sidebar 6.2 (restart "from half the initial step") | "restarts from half the smallest step the tableau tried, as after a failed level"; the design record now named as `docs/PHASE3_API.md` to `docs/PHASE5_API.md` |
| Sidebar 6.3 (restart; "about two dozen evaluations") | "restarts below the smallest step it tried" (no step named, D-20); "measured with a relative error of `\valAliasingRelError` and a relative estimate of `\valAliasingRelEstimate`, after `\valAliasingEvals` evaluations in all" |

Length: the needspace before Sidebar 6.3 went from 10 to 5 lines; a sentence in Sidebar 6.2 pointing to 6.7 was
drafted and cut for length.

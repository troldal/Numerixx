# P5-10-guarantees: the Phase 5 update (PLAN.md 12.8)

Author: `book-author-practice` (Opus 5.5). File: `chapters/10-guarantees.tex`. `tools/check_chapter.sh`:
`COMPILED OK`, 6 pages standalone (5 before, whose last page was nearly full); the text grew by about half a page,
and Table 10.4 now sits alone at the top of page 6. `tools/check.py`: 0 errors, 0 warnings. The editor's layout
edits of Table 10.1 (P5-7: `\tabcolsep` 2 pt, group gaps 6 pt) are kept; reviewed, readable, no overfull box.

| Item of 12.8 | What changed |
|---|---|
| 10.1, "six methods" | "seven methods ... (`BackwardRidders` is the mirror of `ForwardRidders`)" |
| Table 10.1, caption | "Every float derivative above the first warns; these were computed with the warning silenced" (the asterisks of DS15 are now on N = 2 to 4) |
| 10.1, text | A sentence on the `ForwardRidders` row: its columns gain one order rather than two, and its median error on the random study is `\valForwardRiddersAccuracyFactor{First..Fourth}` (14, 212, 1261, 3007) times Ridders'; the float columns "show why every float derivative above the first warns, whatever the method" |
| 10.2, the rules | The floor: "for the one-sided forms, the floor carried through the tableau"; the restart: "at half the smallest step it tried after a step at which f is not finite and after an untrusted result" |
| 10.2, the one-sided study | New paragraph: coverage `\valForwardRiddersCoverage{First..Fourth}`, worst ratios `\valForwardRiddersWorstRatio{...}`, median errors `\valForwardRiddersMedianRelError{...}` (measured; `sec:ridders-onesided`). No restart count (D-20) |
| Table 10.2, caption | The restarts are those after a non-finite value; one after an untrusted result "continues the halving of the step and cannot be told apart" (D-20) |
| 10.3, cost | `ForwardRidders` named with Ridders as measured on the random study, and "`BackwardRidders` costs what `ForwardRidders` costs on the mirrored function"; the caption likewise. The two new rows come from `table-cost.tex` |
| 10.4, Types | "A derivative above the first in a type with an epsilon above 1e-7 keeps few correct digits, by every method, and warns at compile time" |
| 10.4, Domain edges, poles | "The limit holds for a pole, where f is finite on both sides, too": the library's two test poles, `\valPoleCovered` of `\valPoleCases` covered, errors `\valPoleRelErrorMin` to `\valPoleRelErrorMax`, `\valPoleEvalsMin` to `\valPoleEvalsMax` evaluations |
| 10.4, Domain edges, at the edge | "Exactly at an edge, the one-sided methods apply: `Forward<P>` and `Backward<P>`, or with an error estimate `ForwardRidders` and `BackwardRidders`", with their reach from DS23: covered up to k = `\valForwardRiddersReachCoveredExponent{...}`, flagged from k = `\valForwardRiddersReachFlaggedExponent{...}` to `\valForwardRiddersReachLastK`, and the one case between (f' at k = 37, D-18), stated as measured, never as "always flagged" |
| 10.5, where fn is called | "`ForwardRidders` never below x and `BackwardRidders` never above, restarts included, at most `\valForwardRiddersMaxEvals{First..Fourth}` times" (`test_onesided.cpp`, "fn is evaluated on one side of x only, restarts included"; "failures and cost") |
| 10.5, Types | "the float derivatives above the first with the warning silenced at the call" (PHASE5_API.md, `tests/silence.hpp`) |
| 10.6, Codex | "Phases 4 and 5 also by another vendor's model (OpenAI Codex)"; the record is `PHASE3_API.md` to `PHASE5_API.md` and `design.md` |
| Table 10.4 | New row `test_onesided.cpp` (its seven test cases summarised); `test_orders.cpp` adds "poles close to x"; `precision_warning.cpp` "four builds: a float f'' warns through each overload, and a pragma or the macro silences it" (`tests/CMakeLists.txt`); the text says "four build tests" |

Note for the editor: the statement "flagged ... to 45" holds for every order because each
`\valForwardRiddersReachFlaggedThrough{First..Fourth}` equals `\valForwardRiddersReachLastK`; if a regeneration
changed that, the sentence would need the per-order macros.

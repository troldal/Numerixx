# P5-09-tables: the Phase 5 update (PLAN.md 12.8)

Author: `book-author-practice` (Opus 5.5). File: `chapters/09-tables.tex`. `tools/check_chapter.sh`: `COMPILED OK`,
4 pages standalone (as before); `tools/check.py`: 0 errors, 0 warnings. The editor's layout edits (P5-6: Table 9.10
in three columns, `[!htb]`) are kept; reviewed, and the table still fits on the page of 9.4's table.

| Item of 12.8 | What changed |
|---|---|
| 9.5, text | The table collects the constants "of Ridders' method and of its one-sided forms"; the design record is `PHASE3_API.md` to `PHASE5_API.md`; "the one-sided evaluations per level" are named among the values read from the library (the generator reads them from the `Forward<1>` table and `detail::one_sided_divisor`, DS14); "with the evaluations per level it gives the maximum evaluations" (no longer "two or four") |
| Table 9.10, caption | "The constants of Ridders' method and, in the last column, of its one-sided forms (`sec:ridders-onesided`; `BackwardRidders` is `ForwardRidders` mirrored). They differ in the column divisor, since the error of the forward difference has every power of h, not only the even ones, and in the rounding floor, which the one-sided forms carry through the tableau." |

Layout: a first draft explained the two differences in the text (five lines more), which pushed Table 9.10 to a page
of its own (5 pages); the explanation moved to the margin caption, which costs no text height.

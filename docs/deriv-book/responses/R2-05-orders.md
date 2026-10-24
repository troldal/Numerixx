# R2-05-orders: response to review round 2

Author: `book-author-methods` (Fable 5.1), stage 7. File: `chapters/05-orders.tex`. Compiled with
`tools/check_chapter.sh` (`COMPILED OK`, 5 pages standalone, unchanged; the chapter's last page had slack);
`tools/check.py` reports no warning. Net change: +3 lines, as the letter's plan allows.

| Finding | Decision | What changed |
|---|---|---|
| engineer-04 | accept | 5.3: "yet on test functions that vary faster than s, half the default step was better for every order (Sidebar 3.2), so the library keeps one convention ..."; the repository path `docs/PHASE4_API.md` is no longer in the main text (it stays in Sidebar 3.2 and in Sidebar 5.2) |
| codex-09 | accept | Sidebar 5.2: the Sterbenz argument now covers the whole interval: "for 0 < h <= \|x\|/2, x_+ and x lie within a factor 2 of each other and have the same sign, and so do x and x_-, so the subtractions round nothing; for \|x\|/2 < h <= \|x\| the subtraction x - h is itself exact by the same lemma, so x_- = x - h and the offset h_- is h, a representable number, with x_- = 0 at h = \|x\|, while x_+ stays within a factor 2 of x as before". The conclusion (exact offsets for every 0 < h <= \|x\|) is unchanged |
| codex-05 | accept | Sidebar 5.2: "every rule with more nodes divides by d h and then by h again, N - 1 times, in the order that avoids a spurious underflow or overflow (by h first when d h > 1 > h), for the same reason", the clause the header's comment states and the book had omitted |

Nothing else in the chapter was objected to in round 2.

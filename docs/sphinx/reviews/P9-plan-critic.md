# Phase 9, stage 1: docs-critic on the plan

`docs-critic` (Sonnet, high), 2026-11-09, on `docs/sphinx/PHASE9_PLAN.md` as the editor wrote it. Saved by the main
session as the agent returned it; line numbers refer to that version of the plan and to the current pages.

**Verdict.** The plan is sound in shape: one job per page, no merges of the guide pages, and it keeps the programs. It
is too generous in four places. It cuts only 8% of the words. Five budgets are raised above the current counts. Two
of the four figures are not needed. Seven pages carry book links.

## Findings, most valuable first

**1. Drop figure 4 (`deriv_orders`) (plan lines 133-142, 56, 160-161).** Saves one figure, one program, one script,
two validated palette colours, and the figure sentences on the page.
- The page already has a four-line table giving the digits per order (10, 8, 5-6, 4-5).
- The figure is figure 1's V drawn four times, shifted. The only new content is that the rounding branch is steeper
  for higher N, and the reader never acts on that because the default step handles it.
- It needs two further colours for a four-series plot, and the plan itself calls it the one most open to a cut.
- The page then keeps the table and loses the formula sentence, as the plan already intends.

**2. Replace figure 3 (`deriv_edge`) (plan lines 120-132, 52-55, 157-158) with a small table.** Saves one figure, one
script, about 10 plan lines, and the risk of misreading it.
- It is the heaviest figure: two panels, four series, a marker, and 11 distances. Its message is "near the edge, try
  both and keep the smaller estimate only when they agree". Three or four rows (d = 1e-2, 1e-6, 1e-10; error,
  estimate, evaluations for each method) carry the same message.
- The only fact words handle poorly is that Ridders' cost climbs from 12 to 64 evaluations near the edge. A column of
  evaluations in that table shows it.
- The figure plots one function's ranking (central Ridders better far from the edge, worse near it). The text must
  then hedge it with the sqrt(1 - x) counter-case, which is not drawn (lines 131-132). A figure that needs its
  opposite asserted elsewhere to avoid being read as a rule does not earn its place.
- Keep the `deriv_edge` d-sweep as a program asserting the table (or extend `deriv_one_sided.cpp`).

**3. Keep figures 1 and 2, trimmed.**
- Figure 1 (`deriv_trade_off`): keep. It is the chapter's one idea (truncation falling, rounding rising), and it
  replaces the 6-row table, its code loop and about 70 words. Trim it to the V, the default-step line, and Ridders'
  result as one horizontal line. Drop the "seven marked steps" and the evaluation-count machinery that reads steps
  from the points `fn` receives (lines 101-109): no sentence of the page needs "14 evaluations at those seven steps".
  If the author wants the point that Ridders extrapolates from steps all worse than its result, one sentence says it.
- Figure 2 (`deriv_scale`): keep, because "scale" is the one hard idea on `step.md`. It repeats figure 1's shape
  shifted four decades, so it must not carry extras: drop the "flat top near 1" and the 0.5 / "no correct digit"
  assertions (lines 116-119) unless the text mentions them, and ask for step-page budget savings in return.
- Net: two figures instead of four, two new programs instead of four.

**4. Open question 3 (the steps code block).** Drop it, keep only the default call and the Ridders call. The block is
a loop plus a seven-line output that the figure shows better.

**5. Budgets are larger than the job needs (plan lines 25-36).**
- The plan says "5,832 now, 5,350 planned" (line 12). The table's budgets sum to 5,478, because 5,350 leaves out
  `book.md`'s 128. Either way the cut is 6-8%.
- Five pages are budgeted at or above their current counts: `tutorial` 773 to 850, `step` 521 to 550,
  `higher-derivatives` 532 to 550, `derivative-of` 446 to 450, `costs` 327 to 350. A page that gains a figure should
  shrink, not grow.
- Proposed, about 4,350 in all:

| Page | Plan | Proposed |
|---|---|---|
| `index` | 800 | 700 |
| `tutorial` | 850 | 700 |
| `choose-a-method` | 600 | 550 |
| `step` | 550 | 450 |
| `one-sided` | 600 | 500 |
| `higher-derivatives` | 550 | 400 |
| `derivative-of` | 450 | 300 |
| `handle-failures` | 600 | 500 |
| `costs` | 350 | 200 |
| `book` | 128 | 60 |

- Where the savings come from:
  - `index`: "What you get back" (about 160 words) to about 60: two members, "only non-finite values are detected",
    and a link to `basics.md`; that rule is repeated in `handle-failures.md` and `one-sided.md`, so say it once here
    and link. The Tutorial / Guide / Reference scaffolding shrinks to the toctrees plus one line.
  - `tutorial`: the plan says "nothing new" yet raises the budget. Redefining Ridders' method is a clause, not three
    sentences; the root-solve rounding sentences in "Checking it" become one sentence plus the bound. Optional: the
    "second derivative" section repeats `higher-derivatives.md`.
  - `step`: the formulas eps^(1/(P+1)) and eps^(1/(P+N)) also appear in `costs.md` and `higher-derivatives.md`. Keep
    one: "`default_step` gives the default for any method". The "a sixteenth of the width" paragraphs and the
    "output follows the last example" aside go when the figure carries the explanation.
  - `higher-derivatives`: cut the eps^(2/(2+N)) sentence; the Ridders-cost paragraph, which `costs.md` repeats; the
    `order` enum's explanation to one clause; the warning-silencing paragraph (pragmas for three compilers) to the
    macro plus one sentence.
  - `handle-failures`: link to `basics.md` for "test for success, not for a particular error"; shorten the
    `{warning}` to a link to `one-sided.md#a-formula-that-goes-on`; the 7% and 70% cases can lose about 40 words.

**6. `derivative-of.md`: shrink to about 300 words.**
- The `{note}` on `float` (about 55 words) repeats `higher-derivatives.md`; keep one sentence there and drop the note.
- The Newton/Halley paragraph becomes one link to `roots/derivatives.md#numerically`.
- The `std::ref` section is the real content; keep the callable example and the state example only.
- Stronger option, offered as an alternative: drop the page, put the callable example in the overview and the copy
  semantics in `derivative_of`'s Doxygen comment (a change to `include/`, which needs the user's approval).
- `roots/derivatives.md` already has the copy-of-`fn` note. Settle now which page keeps it: the full text here, the
  roots note shortened to a link.

**7. `costs.md`: cut the "Three things stand out" bullets and the last paragraph (about 170 of 327 words).** Each
repeats another page (`choose-a-method.md`, `step.md`, the overview, `higher-derivatives.md`). `roots/costs.md` and
`quadrature/costs.md` carry the table with one line of commentary or none. Accept the rename to `costs.md`
(question 1), but do not let it grow the page.

**8. `one-sided.md`.** The table of finding 2 saves more than the plan's 175-word cut. Keep the advice ("try both,
keep the smaller only when they agree, an estimate is not a proof"). Fold the one-sided stencil paragraph and the
`{note}` on the 1e-9 scale into one sentence each, or move the second to the reference. The `{note}` here and the
warning in `handle-failures.md` say the same thing; keep it here.

**9. Book links (question 4): the overview and `book.md` only.** `roots/` and `quadrature/` have no book links; the
manual stands alone; the book's chapter names change in Milestone 5. Shrink `book.md` to the download sentence and
the PDF note (about 60 words).

**10. `choose-a-method.md`.** Merge the near-duplicate rows (one-sided; near the edge; the same at a fixed cost) into
"at the edge" and "near the edge". In "Measure on your own function", keep the counting-lambda technique, which no
other page shows, and cut the commentary to one sentence.

**11. "The examples write `dv` for `nxx::deriv`" on seven pages (about 85 words).** Move `namespace dv = nxx::deriv;`
into each program's "The program needs:" block, and name `dv` once in the overview and the tutorial.

**12. Terms.** A good set. "Accuracy order P" is needed only on `choose-a-method`; the overview can say "`Central<4>`
and higher use more points and are more accurate on smooth functions". Define "scale" once, in the overview; other
pages link.

**13. The plan's own text.** The figure descriptions are mostly lists of assertions that belong in the programs; keep
what to see, the data and the one fact each figure must show. Drop the paragraph on the tableau figure left to the
book.

## Theory cap

Within the cap apart from three leaks: `index.md` line 16 gives the constant h^2/6 |f'''| (the plan says constants go
to the book); the P/N step formulas on three pages (finding 5). `step.md`'s Ridders starting step max(|x|, 1)/32 is a
number the reader can need: keep it.

## What I would keep as it is

- The page structure: nine pages, no merges of the guide pages, `reference.md` generated; the one-job-per-page table.
- The overview's one idea (truncation against rounding) and the default-step formula eps^(1/3) max(|x|, 1).
- Figures 1 and 2, trimmed.
- The programs table: every kept program, its independent reference, the dyadic steps and `+ - * /` / `std::sqrt`
  for figure data.
- The rename `accuracy` to `costs`; the terms list; the reader's tasks; the Ridders-methods naming.
- `basics.md`'s "plain finite differences" (question 5): leave for the shared pages' phase.
- The `handle-failures` causes list, the 7% / 70% case and its "do not rely on the flag alone" advice, the
  `std::ref` state example, the `step.md` code blocks and outputs, `choose-a-method`'s multiprecision example, and
  the Ridders-on-peak 36-versus-14 evaluations.

## Net effect

Figures 4 to 2; new programs 4 to 2, plus one table assertion added to `deriv_one_sided`; words about 5,478 to about
4,350. Questions: 1 yes, 2 no (drop figure 4), 3 drop the block, 4 overview and `book.md` only, 5 leave, 6 fine.

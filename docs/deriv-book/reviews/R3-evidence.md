# R3-evidence: review round 3 (targeted, the Phase 5 update), evidence

Reviewer: `book-reviewer-evidence` (Claude Opus 5.5), 2026-10-25. Object: the uncommitted Phase 5 update on
`56d27a1` (`git diff 56d27a1 -- docs/deriv-book`) and the PDF `build/gcc/docs/deriv-book/pdf/main.pdf` (90 pages),
read against PLAN.md 4.3, 4.4, 5 and 12.8, `data/README.md` and `docs/PHASE5_API.md`.

Runs (gcc preset): `numerixx_book_generator` (rebuilt: no work; newer than `generate_data.cpp` and `deriv.hpp`, and
`include/` has no uncommitted change) twice into the scratchpad: exit 0, 41 files each, byte-identical to each other
and to `data/` (`README.md` is not generated). `tools/check.py`: 0 errors, 0 warnings. `ctest --preset gcc -R
numerixx_book`: 25 of 25. E24 run by hand: its output equals the macros. A probe of `Ridders` on poles at
$a = \pm 2^{-k}$ (scratch program, for R3-evidence-02). I did not run the full `numerixx_book` target, since it
rewrites `data/`; the generator check above covers it.

## Findings

### R3-evidence-01 (major): the version references predate Phase 5

- **Location:** `chapters/00-frontmatter.tex:50` ("as of commit `f8e98e849b4f`"); `code/getting-started.cmake:9`
  (`GIT_TAG f8e98e849b4f...`), quoted in 7.1.
- **Problem and evidence.** `f8e98e8` is the release revision of the book (24 October, 22:40).
  `git merge-base --is-ancestor f20696f f8e98e8` fails, so the stated commit has no `ForwardRidders` or
  `BackwardRidders`, has the old restart rule, and warns only for float $f''''$. A reader who follows 7.1 gets a
  library for which 6.7, E24, the new rows of Tables 8.2, 9.10, 10.1 and 10.3, and the warning text of 5.4 and 8.6
  are false, and the call in E24's snippet does not compile. PLAN 12.8 schedules "the version references at
  release". It is recorded here so that the release cannot pass without it.
- **Correction.** At release, pin both to the released commit (at or after `f20696f`, the Phase 5 code; `710bda5`
  records its review), in the same edit.

### R3-evidence-02 (major): "the edge limit holds for a pole" goes beyond its evidence, and half of it is untrue

- **Location:** `06-ridders.tex:196-201` (6.4); `10-guarantees.tex:177-181` (10.4, marked "measured");
  `qr-quick-reference.tex:48` ("Near a bad point (a NaN or a pole) it returns a value down to $2^{-34}$ ..."); also
  `07-using.tex:409-411`.
- **Problem.** The limit of 6.4 has two halves: a value is returned down to $2^{-34} s$ for $f'$, and "below that the
  result is `non_finite`". The evidence quoted (`\valPole...`, 8 of 8 covered) comes from the tests' two distances,
  $10^{-3}$ and $10^{-4}$. These are seven orders of magnitude inside the limit ($2^{-34} = 5.8 \cdot 10^{-11}$), so
  they show that such poles are handled; they do not show where the limit lies. 10.4 nevertheless labels the sentence
  "measured". The design record states the conclusion (`PHASE5_API.md:198-199`, "The limit then holds for finite
  poles too") without a measurement at the limit. The second half does not hold. Near a pole, $f$ is finite, so the
  method never fails: when the levels run out, it returns its untrusted fallback.
- **Evidence.** My probe (`Ridders`, default step, $1/(t - a)$ at 0, $a = \pm 2^{-k}$) gave these results:
  - $f'$ is covered down to $a = 2^{-34}$ (relative error \num{2.2e-2}, relative estimate 0.36, 64 evaluations),
    as for $\log x$.
  - At $2^{-36}$, $2^{-38}$ and $2^{-40}$ a *value* is returned, with relative error 1.00 and relative estimates of
    \num{4.9e-3} to \num{8.4e-2}. It is flagged only by the quarter test, because the value is near 0. It is not
    `non_finite`.
  - $f''$ and $f''''$ are covered at $2^{-30}$, and are wrong and flagged at $2^{-32}$.

  R3-codex-01 reports the same thing independently.
- **Correction.**
  - 6.4 and 10.4: "The restart descends as after a non-finite value (`docs/PHASE5_API.md`, item C), so a pole is
    resolved as close as a NaN edge. At the tests' two poles, \num{\valPoleCovered} of \num{\valPoleCases} were
    covered (measured). Closer than the limit, the method does not fail: it returns a wrong value whose estimate is at
    least a quarter of it." In 10.4, mark the limit for poles "derived", or add generated readings at a pole $2^{-34}$
    and $2^{-36}$ away.
  - The quick reference, the same way: "a usable value down to ...; closer, a NaN gives `non_finite` and a pole a
    flagged wrong value".

### R3-evidence-03 (minor): F19 states a claim that the plot does not show

- **Location:** `figures/06-ridders-onesided-reach.tex`; its caption, `06-ridders.tex:266-275`.
- **Problem.** The caption says that, further out, "the result is flagged by an estimate of at least a quarter of the
  value returned". The plot draws the estimate relative to the *exact* derivative, and beyond the reach the returned
  value tends to 0. The flagged region therefore looks like the opposite:
  - the $f''''$ estimate (accent, dashed) falls from 0.62 at $k = 33$ to \num{3.5e-11} at $k = 45$, while its error
    stays at 1.0;
  - the $f'$ estimate falls below the error from $k = 37$.

  The reader sees an estimate many orders of magnitude below the error and is told that it is a flag. The data set
  has the reading (`fr1_flagged`, `fr4_flagged`), but the figure does not use it. The file's comment concedes that
  only the caption says it.
- **Correction.** Mark the flagged $k$ from those columns, for example with ticks along the top edge in each
  curve's colour, or plot the estimate relative to the absolute value returned. At the least, add to the caption:
  "the dashed curves are relative to the exact derivative; the value returned tends to 0, so an estimate far below
  the error can still be a quarter of the value".

### R3-evidence-04 (minor): the AI statement does not cover the update

- **Location:** `00-frontmatter.tex:54-61`.
- **Problem.** The statement is true for the release (`56d27a1`). The update was drafted by the Claude authoring
  roles (`responses/P5-*.md`), integrated by the editor (PLAN 12.8), and reviewed in round 3 by the Claude reviewers
  and by Codex (`reviews/R3-codex.md`). The record shows no Cline review of it. So "OpenAI Codex and Cline with
  DeepSeek V4.1 Flash reviewed the result independently" fully applies only to the first release. The file's header
  comment requires the statement to be made current at release.
- **Correction.** Add one sentence at release, for example: "It was updated for the library's Phase 5 (the
  one-sided Ridders methods, the restart below the smallest step tried, the wider precision warning) in the same way,
  with an independent review by OpenAI Codex." Name Cline only if it reviews the update.

### R3-evidence-05 (minor): the all-powers Richardson extrapolation has no citation

- **Location:** `06-ridders.tex:284-286`; `a2-derivations.tex:143-150`.
- **Problem.** The divisor $2^j - 1$ for a base whose error has every power of $h$ is textbook Richardson
  extrapolation. CLAUDE.md says "Cite for every established method and result". B.3 cites Burden et al.,
  Section 4.2, for the central case two paragraphs earlier, but the one-sided paragraph cites nothing.
- **Correction.** Cite `burden2016` (Section 4.2). To my recollection, its presentation starts from an error in
  every power of $h$, but I could not confirm the passage online this round, so check it against the copy.
  Otherwise cite `sidi2003` (Chapter 1, the polynomial Richardson process). Both entries are already verified.

### R3-evidence-06 (minor): the staged `book.pdf` is the pre-update release

- **Location:** `docs/deriv-book/book.pdf` (staged as a new file).
- **Problem.** It has 88 pages, was created 2026-10-24 at 22:46, and its text contains "ForwardRidders" 0 times; the
  build has 90 pages. Committed with the update, the repository's PDF would contradict its sources.
- **Correction.** Replace it with the final build at release, or keep it out of the update's commit.

### R3-evidence-07 (editorial): "to $k = 45$" is the end of the study, not of the flagged run

- **Location:** `06-ridders.tex:272, 317`; `10-guarantees.tex:190`.
- **Problem.** `\valForwardRiddersReachFlaggedThroughMin` = 45 = `\valForwardRiddersReachLastK`: the results are
  flagged through the last $k$ measured. "From $k = 38$ to $k = 45$" reads as a bounded interval. The P5-06 response
  describes the caption as "to the last k measured"; the integrated caption lost that.
- **Correction.** Write "to the end of the study, $k = \num{\valForwardRiddersReachLastK}$", and let the generator
  fail if the flagged run ends earlier, so that the wording stays true.

### R3-evidence-08 (editorial): a worst ratio of "1.0" for an uncovered case

- **Location:** `10-guarantees.tex:113` (`\valForwardRiddersWorstRatioSecond`).
- **Problem.** $f''$ has one uncovered case (coverage 99.96 %), so its worst ratio lies between 1.00 and 1.05. It is
  printed as "1.0", which reads as covered beside 10.2's "a ratio below 1 is covered". The coverage macros already
  avoid the same trap (`coverage_percent`, R1-numerics-17e).
- **Correction.** Print a ratio above 1 that would round to 1.0 with one more digit.

### R3-evidence-09 (editorial): wording and stale comments

- `06-ridders.tex:305` "from about 16 K" reads as a temperature: "from an initial step of about 16 K".
- `09-tables.tex:153-156` lists the constants read from the library and those transcribed, but not the three new
  rows (column divisor, floor, restart), which the generator writes as text from the header's documentation. Name
  them among the transcribed ones.
- `code/06-ridders-onesided.cpp:6` and `PLAN.md:688` say "measured 0.073"; the program prints a ratio of 0.074.
- `code/05-orders-float.cpp:5` says its printed numbers are the book's through `\valFloatFourth...AtOne`, which 8.6
  no longer quotes (P5-13).

## References

No `\cite` was added or moved; `refs.bib` is unchanged since `56d27a1`.

| Entry | Status | Checked against |
|---|---|---|
| `burden2016` | Verified (R1, R2), unchanged. Section 4.2 still supports Sidebar 6.1 and B.3 as cited; proposed for the one-sided divisor (R3-evidence-05, to be confirmed in the copy) | R1-evidence table (library catalogues) |
| `sidi2003` | Verified (R1), unchanged; the alternative for R3-evidence-05 | R1-evidence table |
| `chapra2015`, `fornberg1988` | Verified (R1, R2), unchanged; 10.6's citations stand at their sentence | R1-evidence table |
| `ridders1982`, `richardson1911`, `richardson1927`, `higham2002`, `goldberg1991`, `sterbenz1974`, `ieee754-2019`, `lyness1967`, `squire1998`, `martins2003`, `griewank2008`, `more2011`, `more2012`, `dennis1996`, `nocedal2006`, `gill1983`, `chase1998`, `nist-webbook`, `shomate1954`, `colebrook1939`, `iso14882-2024` | Verified (R1, R2), unchanged; not cited in the changed text | R1- and R2-evidence tables |

No new material is reproduced, so no licence question arises.

## Checked and correct

- **E24.** `e24_onesided_edge` uses E6's wrapper (now the shared `edge_H`), the default step, a counting wrapper and
  the closed forms `shomate_cp` and `shomate_derivative(..., 2)`.
  - The macros: `\valEdgeForwardRidders{RelError, RelEstimate, Evals}` = 2.1e-13, 2.9e-12, 9; `...Second` = 4.5e-9,
    3.4e-7, 17; `...InitialStep` = 16 (15.625, "about"). E6's `\valEdgeRiddersStatus` and `\valEdgeRiddersEvals`
    (non-finite, 64) are reused (one computation, P5-3).
  - The example prints the same numbers and asserts PLAN section 5's list: `non_finite`; actual error at most 6
    times the estimate (measured 0.074 and 0.013); relative errors at most 1e-10 and 1e-7. It uses the default step
    and the same wrapper.
- **The one-sided study (DS9).** The cases are the same `study_points(f)` as Ridders' (a refactor; the DS9 CSVs are
  unchanged).
  - Median errors 1.5e-13 to 3.5e-6. Factors 14, 212, 1261, 3007, as medians over the medians 1.0e-14 to 1.2e-9 of
    the new `\valRiddersMedianRelError...`, which is consistent after rounding.
  - Coverage 99.96, 99.96, 99.96 and 99.9 %, with Min and Max 99.9 and 99.96. Worst ratios 1.3, 1.0, 1.6 and 1.7,
    with Max 1.7.
  - Evaluations, median (range): 10 (7--11) to 29 (21--37). They are the same in `table-cost.tex`.
- **DS14.** `\valForwardRiddersDivisor` = 32 (`detail::one_sided_divisor`). `\valForwardRiddersMaxEvals...` = 33,
  65, 97, 129 = 1 + 32 N, checked by the generator against a counted run on log at 1e-12. Tables 8.2 and 9.10 use the
  same macros.
- **DS23, recomputed from the CSV.**
  - Cap: k = 25 (f') and 26 (f'''').
  - Covered to 36, 34, 33, 32; flagged from 38, 35, 34, 33, through 45.
  - The one silent case, k = 37: error 0.224 and estimate 0.193, so 4 x 0.193 = 0.772 < 0.776.
  - At k = 30: 9.4e-9 and 0.0014.
  - Ten digits to 29: every k <= 29 is below 1e-10, by DS10's definition.
  - The central `rid1` column is identical to `ridders-edge.csv` for k = 2 to 40, is `nan` from 35 (the limit 34),
    and is more accurate than `fr1` at every k where it returns a value.
- **F19.** It reads the DS23 columns that its caption names. Its ticks 2^-10 to 2^-40 are about 1e-3 to 1e-12.
  "Keeps its accuracy for some halvings after the evaluations are spent" holds (f' at 3e-13 to 3.3e-12 for k = 25 to
  28; f'''' to 28), and so does "central Ridders is more accurate".
- **Sidebar 6.3.** `\valAliasing...` comes from `Ridders` with the default step 1/32 on sin(64 pi x) at 0, against
  64 pi (the exact derivative of the coded function): relative error 0, relative estimate 2.2e-16, 20 evaluations.
  That matches `PHASE5_API.md:209` (24 before Phase 5, 20 after).
- **Poles.** `\valPole...`: 8 cases, 8 covered, errors 1.6e-15 to 3.8e-9, 32 to 88 evaluations. -n!/a^(n+1) is
  right for 1/(t - a) at 0. Every case matches a test (`test_ridders.cpp:261-262`, `test_orders.cpp:580-585`). For
  what they cannot show, see R3-evidence-02.
- **DS6.** `\valFloat{First,Second,Third}CentralRelError` = 5.0e-6, 1.4e-4, 4.0e-4.
  `\valFloatSecondRiddersRelError` = 2.7e-5, with estimate 5.3e-4. `\valFloat{Second,Fourth}ForwardRiddersRelError`
  = 0.0022 and 5.7. All are identical to the float cells of Table 10.1 (the same eleven points, `float_default`), and
  every float call above f' is made under the pragma.
- **DS15.** The `ForwardRidders` row is computed by `matrix_row`, like the other rows.
  `\valExpForwardRiddersFirstRelError` = 1.7e-13 is its double f' cell and `frid_double` at N = 1 (T0, Table 7.1,
  the quick reference). The asterisks are on float N = 2 to 4. Figure 5.x's box encloses all nine float points for
  N = 2 to 4 (2.0e-5 to 0.014); its corner coordinates are typed layout values.
- **Typed numbers.** All are exact or attributed:
  - 500 K (the fit's range);
  - 1e-7 and "7 significant digits" (`deriv.hpp:1810`);
  - s/32 and s/4 (the header's constants);
  - 1e-3 and 1e-4 (the tests' poles);
  - k = 30 and 2^-30 s (`deriv.hpp:302`, "the documented reach");
  - 64 pi, and the factors 4 and 64.

  "Since the rules reach N h" is the design record's reason (`PHASE5_API.md:61`).
- **Claims checked against the tests and the record.**
  - Table 10.4: `test_onesided.cpp` (its test cases) and `precision_warning`'s four builds
    (`tests/CMakeLists.txt:44-64`).
  - 10.6, "Phases 4 and 5 also by ... Codex" (`PHASE5_API.md:340`).
  - The quoted warning is identical to `deriv.hpp:1810-1811`.
  - Table 10.2's restart wording rests on D-20 (identical DS9 data in double).
- **Unused macros.** `\valScaleRiddersDefaultEvals` (56 to 32) and the dropped `\valPoleAhead/Behind...` and
  `\valFloatFourth...AtOne` are no longer quoted anywhere.

## Verdict

The evidence pipeline holds. Every new number in the changed text is a `\val` macro computed as its sentence says.
Two clean runs of the generator reproduce `data/` byte for byte, the new example prints and asserts what the text and
the plan say, and the typed numbers are attributed constants. Two findings must be fixed before release:
- the version pin (R3-evidence-01), which the plan already schedules;
- the pole sentence of 6.4, 10.4 and the quick reference (R3-evidence-02). It extends a limit measured on NaN edges
  to poles on the strength of two distances far inside it, and implies a failure beyond the limit that does not
  happen.

F19's caption makes a claim that the plot visually contradicts (R3-evidence-03). The AI statement needs a line for
the update (R3-evidence-04). The rest is minor or editorial.

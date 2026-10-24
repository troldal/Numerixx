# R2-evidence: review round 2 (targeted), evidence

Reviewer: `book-reviewer-evidence` (Claude Opus 5.5), 2026-10-24. Object: the revision `5432abf`, the PDF
`build/gcc/docs/deriv-book/pdf/main.pdf` (88 pages), read against `reviews/R1-decision.md`, the responses and
PLAN.md sections 4.4, 9 and 12.

Runs (gcc preset): the generator twice into two scratch directories, byte-identical to each other and to the
committed `data/` (40 files); `numerixx_book`: no diff in `data/`, PDF up to date, `check.py` 0 errors, 0 warnings
(the new `signature` check included); `ctest -R numerixx_book`: 24 of 24. The library is unchanged since round 1.

## Round-1 findings

| Id | Status | Remark |
|---|---|---|
| 01 AI statement | Partly | Rewritten and mostly true to the record; the set-up is attributed to the wrong model (R2-02) |
| 02 histogram tail | Resolved | Caption, 6.5 and macros; the "restarts" label is gone |
| 03 calculator ratio | Resolved | `\valHandCalcTruncationHundredth` (4.5e-5, e 0.01^2/6) |
| 04 exact stops | Partly | Renamed and split; one clause quotes the function count for the wrong set (R2-05) |
| 05 cap rationale | Resolved | Design decision with `\valCapGain...` (0.22, 0.70; 0.56, 1.8); U-7 |
| 06 more2011 | Resolved | Cited at the procedure; `more2012` kept at the step |
| 07 complex step | Resolved | `squire1998`, `martins2003`, Lyness and Moler's priority |
| 08 Table 7.1 labels | Resolved | Model and measured separated in rows 1, 5, 6; evidence pointers right |
| 09 Ridders' cost and gain | Partly | "2.3 to 7" typed by decision; "about two more digits" now sits beside "14 against 11" (R2-03) |
| 10 macro misuse at the edge | Resolved | `\valRiddersEdgeEvalsInterior` (16, k = 2) |
| 11 factor of ten | Resolved | `\valModelRatio...` (6.3, 5.6, 5.8, 32, 37) |
| 12 study scope | Resolved | Per order, in the caption of Figure 6.3, 6.3 and 10.2 |
| 13 edge rule | Resolved | Both readings in 6.4, 7.8, 10.4 and the quick reference; two slips in 6.4 (R2-08) |
| 14 whiskers | Resolved | `\valDefaultStepSpreadAbove` 3.5 (Forward<3>; my round-1 3.1 missed it), `...Below` 53 |
| 15 "4 to 16" | Resolved | Moved to 10.6, attributed to the tests' comment |
| 16 Figure 7.2 | Resolved | Plateau about J/(2h c_p) with the steps; checked: 0.768/(2 x 0.0030 x 29.58) = 4.3, and 7J/(12h c_p) = 0.041 for Central<4> |
| 17 Ridders' paper | Partly | "In the form the library implements it" and Sidebar 6.1 fixed; a new sentence attributes unseen content to the paper (R2-01) |
| 18 edge cross-reference | Resolved | Ridders at 500 K generated: non-finite after 64 evaluations |
| 19 JANAF comparison | Resolved | `\valCpTableValue` cited at Part II, p. 1621; `\valCpFitVsTableRelDiff` 1.6e-4 |
| 20 phase reviews | Resolved | |
| 21 design-record figures | Resolved | `docs/PHASE4_API.md` named; B.4's 1.64 sentence corrected |
| 22 numbers in words | Resolved | The three slips corrected; model digits by macro |
| 23 typed caps | Resolved | |
| 24 Figure 3.5 | Resolved | |
| 25 bibliography | Partly | Titles print, notes gone, `ieee754-2019` corrected in the entry; two rendering problems remain (R2-06) |
| 26 licence | Open | With the author (U-9) |
| 27 `data/README.md` | Resolved | |
| 28 typed declarations | Resolved | Exception stated; `check.py` checks the `signature` lines against the headers |
| 29 Colebrook | Resolved | "in its common form" |

`refs.bib` as corrected: `sterbenz1974` (series, "Englewood Cliffs, NJ", ISBN; Thm 4.3.1 at the citations),
`ieee754-2019`, `more2011` and `martins2003` (both match Crossref), `colebrook1939` (note gone) and the removed
notes are all as the round-1 table asked. The remaining comments are accurate. Every citation sits at the sentence
it supports, with the right qualifier where one is given (Burden 1.1, 4.1, 4.2; Chapra ch. 23; Sterbenz Thm 4.3.1;
Higham Thm 2.5; Chase Part II p. 1621). The two unqualified citations at 3.3 (`nocedal2006`) and 3.4 (`dennis1996`)
are fine as they are.

## New findings

### R2-evidence-01 (major): content attributed to Ridders' paper that nobody has seen

- **Location:** `chapters/06-ridders.tex:38-39` (6.1).
- **Problem and evidence.** "Ridders' paper reduces the step by a factor 1.4 and stops when the error grows". The paper
  was not seen (paywall; `refs.bib` says so). The ratio 1.4 and the growth stop are the constants of the Numerical
  Recipes routine `dfridr`, as the design record says (`docs/PHASE3_API.md:167-169`: "These differ from the Numerical
  Recipes routine (`dfridr`) ...: its ratio is 1.4, it stops on a growth factor"). R1-numerics-11 took them for the
  paper's.
- **Suggested change.** "The routine usually taken from Ridders' paper, as in Numerical Recipes, reduces the step by
  1.4 and stops when the error grows; the library uses ..." (the book may cite Numerical Recipes, but not reproduce
  its code). Or drop the paper's constants: "the library's rules are its own: the ratio 2, ...".

### R2-evidence-02 (major): the AI statement gives the set-up to the wrong model

- **Location:** `chapters/00-frontmatter.tex:59-61`.
- **Problem and evidence.** "coordinated by a main session (Anthropic's Claude, Fable 5.1) that set up the book's
  folder, its layout, its build and checks and the skeleton of its data generator". The set-up is commit `7704191`
  (2026-10-23 09:27), co-authored by Claude Opus 5.5, as is `1ffb38f`. Every commit from the plan on (`e1965dd` to
  `5432abf`) is co-authored by Claude Fable 5.1. R1-decision section 4, item 1, gave the set-up to Opus 5.5; PLAN.md
  S6-1 changed it to Fable 5.1. The rest of the statement matches PLAN.md sections 9 and 12, the agent definitions
  and the record: roles and models, the stage-1 reviews, the stage-2 verification, round 1, and the reviews to come
  (round 2, Codex, Cline with DeepSeek V4.1 Flash 1017).
- **Suggested change.** "The book's folder, layout, build and checks and the skeleton of its data generator were set
  up in a session on Claude Opus 5.5; from the plan on, a main session on Claude Fable 5.1 coordinated the agents and
  reported to the author". Confirm with the author if the session switched models.

### R2-evidence-03 (minor): "about two more digits" against "14 against 11"

- **Location:** `chapters/06-ridders.tex:239-243` (6.5), `:254` (6.6); `chapters/07-using.tex:87-90` (7.3);
  `chapters/10-guarantees.tex:115-118`.
- **Problem and evidence.** The medians are Central<2> 5.8e-12 (11.2 digits), Central<4> 6.7e-14 (13.2), Ridders
  1.3e-14 (13.9) and Central<6> 8.0e-15 (14.1). With the headline macros now rounded, the same sentences say Ridders
  gives "about 14 digits against 11", which is three, and in the same breath "about two more digits than
  Central<2>" (6.5, 6.6). 7.3 also says Central<4> gives "about two more digits" and Ridders "as many as Central<4>",
  while quoting 14 for Ridders.
- **Suggested change.** "about three more digits than Central<2> (11 against 14) and well under one more than
  Central<4>" in 6.5 and 6.6, or quote the medians only and drop the digit counts there.

### R2-evidence-04 (minor): the Jacobian's call count is quoted for a different quantity

- **Location:** `chapters/07-using.tex:185-187`.
- **Problem and evidence.** "one `derivative` call per entry ..., at m x n calls, \valJacobianBalanceCallEvals here".
  m x n = 4 calls of `derivative`. The macro (8) counts model evaluations: two per call with Central<2>
  (`generate_data.cpp`, `e22_jacobian`, `JacobianBalanceCallEvals`).
- **Suggested change.** "at m x n calls, each costing two model evaluations with Central<2>:
  \valJacobianBalanceCallEvals model evaluations here, against \valJacobianBalanceHandEvals by hand".

### R2-evidence-05 (minor): "the rest ... spread over 7 functions" counts x^5 too

- **Location:** `chapters/06-ridders.tex:326-327` (Sidebar 6.2).
- **Problem and evidence.** "(3 on x^5, the rest values that happen to equal the rounded closed form, spread over
  \valRiddersZeroFunctionsFirst functions)". The macro (7) counts every function with a zero-error case, x^5 included
  (`generate_data.cpp`, the `seen` array in `ds9_ridders_coverage`). The other 10 cases fall on 6 functions (log 4,
  e^(-x^2) 2, Runge, 1/x, e^x and sin 1 each).
- **Suggested change.** "the 13 spread over 7 functions, 3 on x^5", or have the generator count the functions other
  than x^5.

### R2-evidence-06 (minor): two problems remain in the bibliography as printed

- **Location:** PDF p. 80; `preamble.tex` (`achemso`); `refs.bib` (`ieee754-2019`).
- **Problem and evidence.** (a) `achemso` turns every `\cite[...]{key}` into a separate numbered note. The list has 40
  entries for 23 works, 17 of them pointers ("Ref. 20, ch. 23" five times, "Ref. 8, Sect. 4.1" three times, ...), and
  a reader following a superscript lands on a note, not a work. (b) `@techreport` prints only "IEEE Standard for
  Floating-Point Arithmetic; 2019; DOI" (entry 34), without the standard's number or publisher. Tested in a scratch
  document with the book's `achemso` options: `@manual` drops the DOI, and `@book` with `series = {IEEE Std 754-2019
  (Revision of IEEE 754-2008)}`, `publisher` and `address` prints all of them.
- **Suggested change.** (a) Put the qualifier in the text next to an unqualified citation ("Burden et al.,\cite{burden2016}
  Section 4.1"), or accept the notes as the journal style's way (the author's choice; a preamble question).
  (b) Use `@book` with `series` as above.

### R2-evidence-07 (editorial): "every number ... produced by the library itself" has attributed exceptions

- **Location:** `chapters/00-frontmatter.tex:39-41` ("The evidence"), `:75-76`.
- **Problem.** A few computed figures are typed with an attribution, by decision: the design record's prototype
  figures (Sidebar 5.2: 8.7 and 1.4; B.4: 1.2 to 1.8, 3.5, 11, 15; Sidebar 6.2: 16 to 256), the tests' "4 to 16" (10.6)
  and the "2.3 to 7" of 6.5.
- **Suggested change.** "... produced by the library itself, apart from a few figures of the library's design record
  and tests, which the text attributes".

### R2-evidence-08 (editorial): two slips in the edge rule of 6.4

- **Location:** `chapters/06-ridders.tex:177-180`.
- **Problem.** "(derived from the initial step and the level cap, and measured here)" applies to 2^-34 for f' only;
  10.4 says rightly that 2^-31 is derived, not measured. And "2^{-30}" is typed, while 10.4 and the quick reference use
  `\valRiddersEdgeTenDigitsExponent`.
- **Suggested change.** "(derived from the initial step and the level cap; for f' also measured here)" and the macro.

### R2-evidence-09 (editorial): a stale comment

- **Location:** `figures/06-ridders-coverage.tex:4`: "An actual error of exactly 0 (an exact stop)".
- **Suggested change.** "(a bit-exact result)", as in `data/README.md`.

## New macros traced to the generator

Each recomputed from its code in `tools/generate_data.cpp` and, where possible, from the CSVs or by hand; each matches
the sentence that quotes it:
- `RiddersTail{Cases,MinEvals,Restarts,Functions}Second`: cases above the median count (17): 332 at 19 or more, 11
  restarts, e^(-x^2) 117, 1/x 69, Runge 67, log 46, sqrt 33 (recounted from `ridders-coverage-N2.csv`).
- `RiddersZeroError*`, `RiddersZeroOnQuinticFirst` (3), `RiddersZeroFunctionsFirst` (7, see R2-05).
- `RiddersUncoveredMain*First` and `Second` (e^(-x^2), 8 at x = 1.62 to 2.93 and 4 at x = 2.08 to 2.87; the other x^5),
  `RiddersCoverageThird` (99.96, from 2699 of 2700).
- `MisjudgeCost{CentralTwoFirst,CentralTwoFourth,CentralFourFourth}`: max of (N 2^-P + P 2^N) and (N 2^P + P 2^-N),
  over P + N: 1.7, 5.5, 8.0.
- `CapGainCentralTwelve*`: -log10(eps) (12/(12+N) - 10/(10+N)) with eps = 2^-52 and 2^-167: 0.22, 0.56, 0.70, 1.8.
- `ModelDigitsCentralTwo*Double`: 52 log10(2) times 2/(2+N): 10.4, 7.8, 6.3, 5.2.
- `ModelRatio*`: the DS15 medians against eps^(P/(P+N)): central rules below by at most 6.3 (Central<2>, f'), above by
  at most 5.6 (Central<6>, f''''); Forward<2> up to 5.8 (f''''); Forward<4> 32 and 37.
- `NoiseEstimate*` (DS22): `book::estimate_noise` at eight points (the evaluation count is the number of points, by
  construction); friction 4.3e-11 against the actual RMS deviation 4.9e-11; noisy e^x 6.2e-9 against
  eta/sqrt(3) = 5.8e-9; clean e^x 9.2e-17.
- `JacobianBalance*` (E22): 4 hand-applied model calls, 8 through `derivative` (see R2-04), both 2.4e-11.
- `ExpRiddersFirstRelError` (1.3e-14, the DS15 Ridders median for f'), `RiddersEdgeLimitRelEstimate` (0.10, the k = 34
  row), `RiddersEdgeEvalsInterior` (16, the k = 2 row), `EdgeRiddersStatus` and `EdgeRiddersEvals` (non-finite, 64).
- `CpTableValue` (32.697: a cited constant from Chase 1998, Part II, p. 1621, not a computation) and
  `CpFitVsTableRelDiff` (abs(32.69173 - 32.697)/32.697 = 1.6e-4), `HandCalcTruncationHundredth` (4.5e-5),
  `DefaultStepSpread{Above,Below}` (3.5, 53).
- Headline digits rounded (`lround` of -log10 of the median): 11, 5, 14 and 9. They are quoted by macro wherever they
  occur (1.5, T0, 3.6, 4 intro, 5.2, Figure 5.1, 6.5, 7.3, Table 7.1, 10.1, 10.3, the quick reference).

Revised captions: Figure 3.5 (floors at or below eps^(2/3)), Figure 4.2 (with the text's 3.5 and 53), Figure 6.3
(scope per order, zero-error count), Figure 6.5 (the tail, label removed) and Figure 7.2 (the plateau with the steps)
all say what their data show. No typed computed number slipped in beyond those noted in R2-07; the few new typed
figures are exact or arithmetic (293 / 20 = 15, the tick labels 2^-k, the toolchain versions as in `ci.yml`, the
60-digit references of `test_math.cpp`).

## Verdict

The revision resolves 23 of the 29 round-1 findings outright and 5 in part (01, 04, 09, 17, 25); finding 26 stays
with the author.
`data/` is current and reproducible, every new macro traced is computed as its sentence says, and the corrected
references are right. Two points need a fix before release. A sentence attributes Numerical Recipes' constants to
Ridders' unseen paper (R2-01), and the AI statement gives the set-up to the wrong model (R2-02). The other seven
findings are wording and rendering. Counts: 2 major, 4 minor, 3 editorial.

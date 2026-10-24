# R2-numerics: targeted review of the revision (round 2)

Reviewer: `book-reviewer-numerics` (Fable), 2026-10-24, on the revised build (commit `5432abf`, 88-page PDF).

Read: `reviews/R1-decision.md` (the decisions on R1-numerics-01 to 17 and themes T7 to T13), the responses that cite
them (`responses/R1-01-why.md`, `02-taylor`, `03-rounding`, `04-stencils`, `05-orders`, `06-ridders`, `09-tables`,
`10-guarantees`, `a1-floating-point`, `a2-derivations`, `qr-quick-reference`), the diff `aeec7e6..5432abf` of every
chapter source, of `tools/generate_data.cpp`, `data/values.tex`, `data/README.md`, `macros.tex`, the three Chapter 6
figure files, the new `code/noise.hpp` (`estimate_noise`), `code/colebrook.hpp`, `code/03-rounding-noise-estimate.cpp`
and `code/07-using-expected.cpp` (the `ridders-check` snippet), the current `chapters/07-using.tex`, PLAN.md section
4.4, and the PDF pages of 3.6 and of Chapter 6 (Figures 6.3 to 6.5, Sidebars 6.1 to 6.3).

Recomputed (scratchpad `r1num/`): every new macro of PLAN.md 4.4 against the generator's code and my round-1
numbers (`checks.py`, `tail.py`); two scratch programs against the library itself (`margin2.cpp`: Ridders on the
exact uncovered points of the study with smaller initial steps; `margin.cpp`: the noise recipe at spacings on and off
the function's scale). Everything else was re-read, not recomputed.

## Round-1 findings: status

| Id | Status | Note |
|---|---|---|
| 01 | resolved | Eq. (4.3) reads `+ (h^4/30) f^(5)(xi)`; the sentence after it gives the sign's meaning (the rule underestimates f' where f^(5) > 0); the sidenote points to B.1 for the single xi |
| 02 | resolved | Fig. 6.5's caption and 6.5 quote the generated tail (332 cases of 19 or more evaluations, five functions with counts, 11 restarts all at the most evaluations); the "restarts" label is gone from the figure (checked in the PDF) |
| 03 | resolved | "An actual error of exactly 0" in Fig. 6.3 and Sidebar 6.2, with the split 3 on x^5 and 7 functions for f', all 8 on x^5 for f''; "the study records values, not the stop reason". One wording nit: R2-numerics-05 |
| 04 | resolved | B.1 gives the asymptotic form for every rule, names the products (1, 0, 1 for Eq. 5.1; -8/3, 2/3, 2/3, -8/3 for Eq. 4.3; every central rule for f' above order 2; some one-sided rules), and derives the single xi of (4.3) from the interpolation error formula with the node product 4 and m = 5. Correct. One nit on the O(h^(P+1)): R2-numerics-06 |
| 05 | resolved | 10.1 quotes `\valModelRatio...`: 6.3 (central, N <= 3, below), 5.6 (central, f'''', above), 5.8 (Forward<2>), 32 and 37 (Forward<4>, f''' and f''''); recomputed from `accuracy-matrix.csv` (my 33 is their 32 before rounding) |
| 06 | resolved | `\valMisjudgeCost...` = 1.7, 5.5, 8.0 (recomputed: 1.67, 5.50, 8.03) in 3.5, Sidebar 3.2, Sidebar 5.1 and 5.3; the 16 is named as the rounding term alone |
| 07 | resolved (modified) | Caps as design decisions of `docs/PHASE2_API.md`; `\valCapGain...` = 0.22, 0.56, 0.70, 1.8 by the exponent alone (recomputed; with the constants S and C the model gives 0.20 and 0.68 for f', 2.4 for f'''' at 50 digits, so "by the model's exponent alone" is the honest tag); the one-sided sentence is qualitative and attributed |
| 08 | resolved | 4.3 and Sidebar 3.1: the outer nodes' rounding is quadrupled (1/12 to 4/12), the inner nodes' untouched |
| 09 | resolved | Fig. 6.3's caption, 6.3 and 10.2 scope the "no singularity" statement to f' and name the restarts for the higher orders. One wording nit: R2-numerics-04 |
| 10 | resolved | Sidebar 6.2's floor reads `sum_k |c_k| / (d h^N)` with the integer numerators |
| 11 | resolved | 6.1: ratio 1.4 and the stop on growth in the paper; ratio 2, floor, stall, final check and restarts in the library; "every statistic in this chapter is of the latter"; the introduction says "in the form the library implements it" |
| 12 | resolved (modified) | Sidebar 6.2 names the uncovered cases by macro: 8 of 9 are e^(-x^2) at x in [1.6, 2.9] (the other x^5), all 4 for f'' e^(-x^2) at [2.1, 2.9], the Runge function at 0.54 for f'''; recomputed from the CSVs (my round-1 "all nine" missed the x^5 case at 0.69, ratio 1.008; the generator is right). The remedy added with it is a new finding: R2-numerics-03 |
| 13 | resolved | 3.6: cancellation inside f, k eps |f|, treated as the noise level |
| 14 | resolved | Sidebar 6.1: 4^j exact, 4^j - 1 rounds from j = 13 in float and 27 in double |
| 15 | resolved | 6.3: the parent's error is 4^j times the left neighbour's, the entry's own smaller still |
| 16 | resolved | 6.4: 24 against `\valRiddersEdgeEvalsInterior` = 16 (log x at 1/4, no restart); four for the two failed levels, the rest a tableau from 1/128 (16 + 4 + 4, consistent with ten levels against eight) |
| 17 | resolved | (a) "first four digits"; (b) 12! reached, 12! 4! a bound; (c) "13 halvings"; (d) the decades left out of Fig. 4.1's caption; (e) 99.96 in Table 10.2 (`coverage_percent` prints two decimals only where one would read 100.0 with a case uncovered); (f) Sterbenz for a, b of one sign in A.5 and Sidebar 5.2; (g) "up to a known factor", the factor sqrt((2k)!/k!^2) named and divided out; (h) the model digits 10.4, 7.8, 6.3, 5.2 by macro, the headline digits rounded (11, 5, 14, 9); (i) "ten orders of magnitude above the 2.7e-13 it reports 50 K from the break"; (j) "2.3 to 7 times"; (k) the two readings of the edge rule in 6.4, 7.8, 10.4 and the quick reference, with `\valRiddersEdgeLimitRelEstimate` = 0.10 (the CSV's 0.1026 at k = 34) |

Nothing from round 1 is unresolved. The three nits (05, 06 of this round and the "survives by restarting") are
editorial and listed below so that they are not lost.

## New findings

### R2-numerics-01 (minor). The Ridders tolerance test is relative, and fails where the derivative is near zero

Location: `chapters/01-why.tex` lines 193-194 (1.5, "Is it right?"); `code/07-using-expected.cpp` lines 111-114
(`ridders-check`, quoted in 7.2); `chapters/qr-quick-reference.tex` ("For Ridders, also require ..."); the pitfall "A
Ridders result that is a value but poor" in 7.8.

Problem: the new check `*r->abs_error <= tol * std::abs(r->value)` is a relative test on the value. Where the
derivative is near zero it rejects good results: the book's own example of 7.2, cos x at 0, returns exactly 0 with
`abs_error` at the rounding floor, so the test is false for every `tol`; sin x near pi/2 (7.5) has a relative
error of 7.4e-7 at an absolute error of 2.7e-12, and 7.5 tells the reader to judge such a result by its absolute
error, "in units of f per unit x". The check of 1.5 and the quick reference therefore contradicts 7.5 at exactly
the points where the estimate matters most (a stationary point, a sensitivity that happens to vanish), and it is
presented as the general idiom ("use r->value only if ...").

Evidence: `abs_error` is absolute by construction (8.3); 7.5's `check-absolute` numbers; the floor is positive for
every nonzero f (Sidebar 6.2), so a zero derivative never passes a relative test.

Correction: test against an absolute tolerance built from the problem's scale, as 7.5 does:
`*r->abs_error <= tol * std::max(std::abs(r->value), scale)` with `scale` about |f| / max(|x|, 1) (or the
accuracy you need in the derivative's units), and say in 1.5 and the quick reference "relative to the value, or
to the scale |f|/max(|x|, 1) where the derivative may vanish (7.5)". One line in the snippet, one clause in each
of the three places.

### R2-numerics-02 (minor). The noise recipe's preconditions and failure modes are not stated

Location: `chapters/03-rounding.tex` 3.6, the paragraph "If the noise level is not known, it can be measured" and
the `noise-estimate` listing (PDF page 20); `code/noise.hpp` (the header comment says what the chapter does not).

Problem: the text states the precondition "points spaced far more closely than the scale on which f varies" and
calls the recipe "a simplified form of More and Wild's procedure", but not what the simplification drops nor what
the recipe does when its assumptions fail. Measured with the book's own `estimate_noise` (`margin.cpp`): on
sin(50 x) (1 + 1e-8 u) at x = 1, a spacing h = 0.1 (on the function's own scale) "settles" at k = 1 and reports a
noise of 2.2 |f|, a wrong value without any failure; h = 1e-2 and 1e-3 report no estimate (k = 0, sigma NaN);
h = 1e-4 and the default step 6.1e-6 report 4.7e-9 and 5.6e-9 against eta / sqrt(3) = 5.8e-9. On the clean
e^(-x^2) at x = 2.5, h = 0.5 gives 0.25 |f| at k = 1 and h = 0.01 no estimate (the clean differences keep shrinking
and no order settles within the seven available). So the recipe can (a) return a wrong level when the spacing is
on f's scale, which the full procedure guards against with its sign test and its spacing adjustment, and (b)
report nothing (k = 0) when the spacing is too small for the smooth part to fall below the noise within the eight
points, which the example checks (`if (e.k > 0)`) but the chapter does not mention. The three cases the book runs
are all in the recipe's working range, so its numbers are right; the reader who applies it elsewhere is not told
how to see that it is not.

Correction: two sentences after the listing: "The recipe reports `k = 0` and no value when no three orders agree,
which happens when the spacing is so small that the smooth part has not dropped below the noise within the eight
points; a spacing on the scale of f itself can instead make the orders agree by accident and report a wrong level
(More and Wild's full procedure tests the signs of the differences and adjusts the spacing). Use the default step
of Central<2>, as here, and check `k`." Optionally name sigma as a root mean square, as the generator's comment
does, since the chapter's eta is a half-width (the factor sqrt(3) is immaterial for the step).

### R2-numerics-03 (minor). "A smaller initial step restores the margin" is stated as a fact

Location: `chapters/06-ridders.tex` Sidebar 6.2 (PDF page 40), the sentence added for R1-numerics-12.

Problem: the remedy is given without evidence and holds usually, not always. Measured on the study's 13 uncovered
cases with their exact x (`margin2.cpp`, which reproduces the study's ratios at the default step to three digits):
with a quarter of the default initial step 10 of the 13 are covered (ratios 0.02 to 0.41; three return the same
entry as before, since an entry depends only on the levels below it, and keep their ratios 1.41, 1.16 and 1.36);
with a sixteenth, 12 of 13 (one at 1.09); with a sixty-fourth, all 13 (0.00 to 0.61). The ratio at a single point
is luck in the rounding, as 4.2 says of the stencils, so the sentence should be tagged as a rule of thumb or
measured.

Correction: "a smaller initial step usually restores the margin (on the study's 13 uncovered cases a quarter of
the default step covers 10, a sixty-fourth all)", or attribute it to the argument of B.3 (the change from the
parent is conservative while the expansion is asymptotic) with "as a rule".

### R2-numerics-04 (editorial). "Which the method survives by restarting"

Location: `chapters/06-ridders.tex` 6.3, the sentence on the study's scope (PDF page 37, top).

Problem: for the higher orders the initial step reaches a singularity at small x, "which the method survives by
restarting". Sidebar 6.2 says, correctly, that log x and sqrt(x) restart (NaN samples) while 1/x does not: its pole
lies between finite samples and the tableau takes more levels. Correction: "by restarting, or, where the samples
stay finite, by taking more levels (Sidebar 6.2)".

### R2-numerics-05 (editorial). "Spread over 7 functions"

Location: Sidebar 6.2, the zero-error sentence: "(3 on x^5, the rest values that happen to equal the rounded closed
form, spread over 7 functions)". `\valRiddersZeroFunctionsFirst` = 7 counts x^5; the rest lie on six other
functions. Correction: "spread over six other functions" (generate the count without x^5), or "13 cases on 7
functions, 3 of them on x^5".

### R2-numerics-06 (editorial). B.1's `+ O(h^(P+1))` needs one more derivative

Location: `chapters/a2-derivations.tex` B.1: "since each xi_k tends to x with h, the remainder is M h^P f^(N+P)(x)
/ (N+P)! + O(h^(P+1)) for every rule". With f^(N+P) continuous, f^(N+P)(xi_k) - f^(N+P)(x) is o(1), so the
remainder is o(h^P); O(h^(P+1)) needs f^(N+P+1) bounded, one derivative more than the expansion uses. The
appendix's preamble ("as many continuous derivatives as the expansions use") does not cover it. Correction: "o(h^P),
and O(h^(P+1)) when f has one more bounded derivative, as the error model assumes".

### R2-numerics-07 (editorial). Two small attributions

(a) `chapters/06-ridders.tex` 6.4: "about 2^-34 s away for f' and 2^-31 s for the higher orders ... (derived from
the initial step and the level cap, and measured here, Section 10.4)": "measured here" is true of 2^-34 only;
10.4 says the second is derived, not measured. Say so in 6.4 too.

(b) Two generated macros used by the text are not in PLAN.md 4.4 nor in `data/README.md`:
`\valRiddersEdgeLimitRelEstimate` (DS10; 6.4) and `\valExpRiddersFirstRelError` (DS15; 6.5, 7.3). Both are
computed correctly in `generate_data.cpp` (the relative estimate at the last x that returns a value, 0.10; Ridders'
median relative error for f' in double, 1.3e-14). An entry each in the two tables.

## Checked and correct (round 2)

The new and revised passages verified, so that they need not be re-read:

- Eq. (4.3) and its sentence; B.1's two arguments (the products, the interpolation error formula with x among the
  nodes, the node product 4 and m = 5 giving + h^4 f^(5)(xi)/30; the one-sided rules covered by the same formula;
  the asymptotic form for the rest); the Chapter 4 sidenote.
- 4.3 and Sidebar 3.1 on the outer-node distance (1/12 to 4/12, inner nodes untouched); Sidebar 3.1's "a fixed
  multiple of the nominal step".
- The generated constants of PLAN.md 4.4, each recomputed: `MisjudgeCost` 1.7, 5.5, 8.0 (the larger of the halved
  and doubled cases at the balance); `CapGainCentralTwelve` 0.22, 0.56, 0.70, 1.8 from -log10(eps) (12/(12+N) -
  10/(10+N)); `ModelDigitsCentralTwo` 10.4, 7.8, 6.3, 5.2; the headline digits 11, 5, 14, 9 as lround(-log10 of
  the medians 5.8e-12, 1.2e-5, 1.3e-14, 6.8e-10); `ModelRatio` 6.3, 5.6, 5.8, 32, 37; `DefaultStepSpread` 3.5
  (Forward<3>) and 53 (Forward<1>) from `default-step-accuracy.csv`; `HandCalcTruncationHundredth` 4.5e-5, and the
  remainder of that row's 2.2e-4 is the six-digit rounding of the two values (+3.5e-6 over 0.02 = 1.7e-4);
  `CpFitVsTableRelDiff` 1.6e-4; `EdgeRidders` non-finite after 64 evaluations (32 failed levels of two, the cap);
  the DS9 tail, zero-error and uncovered macros for N = 1 and 2 against my round-1 counts; `RiddersEdgeEvalsInterior`
  16 and `RiddersEdgeLimitRelEstimate` 0.10 against `ridders-edge.csv`; the coverage 99.96 for f'''.
- DS22 and `estimate_noise`: gamma_k = k!^2/(2k)! = 1/C(2k, k) by the recursion gamma_k = gamma_(k-1) k/(2(2k-1)),
  which removes the growth of the k-th differences of independent noise (variance C(2k, k) sigma^2); the in-place
  difference table; the acceptance "three consecutive orders within a factor 4", which is More and Wild's
  criterion; the three reported values (4.3e-11 against the 4.9e-11 root mean square of the bisection output's
  deviation, both below the half-width 1e-10 and near eta / sqrt(3) = 5.8e-11; 6.2e-9 against 5.8e-9 for the noisy
  e^x; 9.2e-17, about 0.4 eps, for the clean e^x). In the recipe's working range the numbers are right (R2-02 is
  about what the text says outside it).
- 3.2's "factor of only 21, not 100: truncation 4.5e-5, the rest rounding"; 3.3's "about log10(|f| / (w |f'|))
  digits" and the formation-enthalpy example (a few hundred kJ/mol against T0 cp = 33 kJ/mol: about one digit);
  3.5's kelvin against Celsius (293/20 = 15, about one digit on the rounding branch) and the 1.7 of the scale
  misjudged by two; 3.6's cautions (a residual tolerance is not eta; cancellation as k eps); 3.8's sixth question.
- Chapter 6: the 24 = 16 + 4 + 4 evaluations at x = 0.01; the two readings of the edge rule (a value to 2^-34 s
  with the estimate at a tenth of the value, ten digits to 2^-30 s); the tail of Fig. 6.5 (332, 19, 11, five
  functions with counts); the 2.3 to 7 cost and the digits against Central<2>, <4> and <6> (1.3e-14 against
  5.8e-12, 6.7e-14, 8.0e-15); Sidebar 6.1's bit counts; Sidebar 6.2's floor with c_k, the untrusted test "at |A|
  (factor 1)"; Figures 6.4 and 6.5 as rendered (regions labelled, no "restarts" label).
- 7.4.3: eta is the tolerance on 1/sqrt(lambda), which gives lambda a relative error of at most eta (a factor 2 on
  the bisection's half-width), so "the tolerance is on the output itself" holds; the fixed-point sentence now says
  the estimate did not cover the error (2.7e-11 against 1.4e-12) and gives the action. 7.9's mixed-partials rule
  of thumb: the inner Central<2> error of about eps^(2/3) acts as eta for the outer call, whose step is then
  eps^(2/9) max(|V|, 1) and whose error about eps^(4/9) = 1.1e-7 in double, at 2 x 2 = 4 evaluations: consistent
  with the model of 3.6. E22's steps and counts (2 and 8 model calls, errors about 1e-11).
- 10.1 (the ratios), 10.2 (the scope), 10.3 (the gain), 10.4 (the caps as decisions, the two readings), 10.5
  (the caps by macro, h > 0 for the one-sided guarantee); 9.2 and 9.3's cap sentences; 4.5; 5.2's digits; 5.3 and
  Sidebar 5.1's factors; A.5's Sterbenz with signs; A.4's two cases at the resolution of x.

## Verdict

Every round-1 finding is resolved, fifteen as decided and two in the modified form the decision letter chose; the
sign of Eq. (4.3), the proof gap in B.1, the two misreadings of the study's data and the loose factors are gone, and
the quantitative statements that replaced them are generated and, where I recomputed them, right. Three minor
points remain, all in material added at this round: the Ridders tolerance check is relative and contradicts 7.5
near a zero of the derivative (R2-01); the noise recipe's text does not say what the simplification drops or how
its failure shows (R2-02), although its three quoted values are correct and the recipe is sound in its working
range; and the "smaller initial step" remedy is a rule of thumb stated as a fact (R2-03, measured: usually, not
always). Four editorial slips follow. With those addressed, the numerics of the book withstand scrutiny.

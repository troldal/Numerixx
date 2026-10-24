# R2-decision: the editor's decisions on review round 2

Decisions of 2026-10-24: the author had instructed the main session to carry the process to the end without further checkpoints, and the main session took the recommended option on every point of section 5(a): (1) citation qualifiers moved into the prose, no change of citation style; (2) the acceptance check in two tolerances; (3) the kink wording; (4) the DS9 extension for the smaller-initial-step sentence; (5) the Ridders-paper constants dropped, no Numerical Recipes citation; (6) the AI statement with the set-up on Claude Opus 5.5 and the coordination from the plan on by Claude Fable 5.1, as git records; (7) the version statement and the CMake pin name the full hash of the final-fixes commit, a tag may replace it later. Points (8) the repository test rule and (9) a header documentation note on kinks change the library or its rules and are left to the author as U-14 and U-15.

Editor: `book-editor` (Fable 5.1), 2026-10-24. Object: the six round-2 reviews of the revised book at `5432abf`
(88 pages): `R2-engineer.md` (25 round-1 findings: 22 resolved, 3 partly; 4 new: 2 minor, 2 editorial),
`R2-developer.md` (18 resolved; 3 new: 1 minor, 2 editorial), `R2-numerics.md` (17 resolved; 7 new: 3 minor, 4
editorial), `R2-library.md` (21 resolved, 1 partly; 7 new: 2 minor, 5 editorial), `R2-evidence.md` (23 resolved, 5
partly, 1 with the author; 9 new: 2 major, 4 minor, 3 editorial) and `R2-codex.md` (OpenAI Codex, independent; 10
findings: 1 major, 8 minor, 1 editorial, with compiled counterexamples) and `R2-cline.md` (Cline with DeepSeek V4.1
Flash, independent, extracted from a timed-out run in several passes; the editor assigned the ids R2-cline-01 to -21
in order of appearance to the findings, skipping the working notes: 2 major, 8 minor, 10 editorial and one
cross-scope note on the header).

Every new finding is decided below by its id: **accept**, **modified** (accepted in a modified form, with the form
stated), **reject** (with the reason) or **defer**. Where two or more reviewers raise the same point it is decided
once, in section 1, and the per-finding rows point there. The partly resolved round-1 findings are closed by the
decisions on the new findings that carry their residue (engineer-03 by T9, engineer-19 and -25 by engineer-04 and T3,
library-16 by library-05, evidence-01, -04, -09, -17 and -25 by T5, T4, evidence-03, T6 and T7).

Counts: engineer 4 accept; developer 3 accept; numerics 5 accept, 2 modified; library 7 accept; evidence 7 accept,
1 modified, 1 for the author (T7); Codex 9 accept, 1 modified; Cline 18 accept, 2 modified, 1 reject (the missing
year of a database entry). Every Codex and Cline finding that needed reproduction was reproduced before deciding
(sections 2.6 and 2.7). Where Cline and Codex or a Claude reviewer found the same thing, the point is decided once
and both rows say so: the B.4 multiplier (codex-02, cline-01), Figure 7.2's plateau (codex-04, cline-14), the tests'
references (codex-07, cline-05), Table 9.4's caption (codex-10, cline-06), the Jacobian's call count (evidence-04,
cline-12) and the section qualifiers of three citations (T7, cline-09).

Nothing is revised yet (stage 7's final fixes, by the authors, per `responses/README.md`). Section 3 is the revision
plan by chapter file; section 4 the editor's own items; section 5 what needs the author; section 6 the page budget.

## 1. Recurring themes, decided once

### T1. The Ridders acceptance check is relative only (engineer-01, developer-01, numerics-01, library-01; codex-01 in part)

Four reviewers converge on the new check of 1.5, 7.2 (`ridders-check`) and the quick reference, `*r->abs_error <=
tol * std::abs(r->value)`: it rejects every good result near a zero of the derivative (cos at 0, the book's own
example three lines earlier; sin near pi/2, which 7.5 says to judge by its absolute error), gives no guide to `tol`,
and, generic as written, dereferences an empty `optional` on a stencil result (undefined behaviour; `true` in a
release build, an assertion under `_GLIBCXX_ASSERTIONS`, developer-01's probe). Codex adds that a result wrong by 50 %
can pass it (T2). **Decision: accept, once, in the two-tolerance form** (library-01's first option): the helper reads

    auto const accept = [](auto const &res, double atol, double rtol) {
        return res && res->abs_error
            && *res->abs_error <= atol + rtol * std::abs(res->value);
    };

A stencil result fails it, which is the honest answer: it has no estimate to test. The text of 1.5, 7.2 and the quick
reference says, in one sentence each: `rtol` is 10^-n for a value reported to n significant figures; `atol` is the
error you can accept in the derivative's units and is what decides near a zero of f' (7.5); inside an iteration no
check is needed, the iteration tolerates an approximate slope (7.4.1); and the check says that the estimate is small,
not that the value is right where the estimate is blind (T2). The scale form of numerics-01, `tol * max(|value|,
|f| / max(|x|, 1))`, is not used: it needs the size of f, which the helper does not have, and the standard
absolute-plus-relative form says the same thing with what the caller knows. The example keeps its two calls, with
`atol` 0 and `rtol` 1e-8, and prints `true false` as now (the editor confirms the output at integration).

### T2. What the checks cannot see: a kink (codex-01; with 6.4, 7.5, 7.8, Table 7.1, Figure 7.2, the quick reference)

Reproduced (section 2.6): f(t) = t + |t + 2^-40| at 0 has f'(0) = 2; Ridders returns 1.000000000067909 with
`abs_error` 2.15e-10, which passes the relative check at 1e-8; g(t) = t + |t| at 0 returns exactly 1 from `Central<2>`
at h = 0.125, 0.25 and 0.5 and from Ridders with an estimate of one epsilon, although g'(0) does not exist. The book
says in 6.4 that "a range break or a kink within h0 of x" gives a large estimate, in 7.5 that "a kink or a jump inside
the span does not survive" the 2h and h/2 check, and on the quick reference and in Table 7.1 that "only Ridders'
estimate flags it". All three are true of a jump in f (the measured nitrogen break, J = 0.77 J/mol, where the plateau
J/(2h) doubles as the step halves) and false of a kink, a jump in f' with f continuous: at x itself the central
difference of c |t| is zero at every step, so every rule returns the mean of the two one-sided slopes plus the smooth
part's own error, the tableau converges to that mean like a smooth function, and the estimate sits at the floor; the
same holds for a kink so close to x that no level resolves it (Codex's 2^-40, which 32 halvings from 1/32 never reach).
This is exact arithmetic, not a measurement, so no data set is added and the counterexample's numbers are not
quoted. **Decision: modified.** 6.4 restricts the large-estimate statement to a jump and adds the kink sentence (the
mean of the slopes, the estimate at the floor, "nothing in the samples can show it"); 7.5 says "a jump inside the span
does not survive this check ...; a kink at x does: every step returns the mean of the two slopes (6.4)"; 7.8's kink
pitfall gains "a kink at x itself shows no symptom at all: know the branch boundaries", which is its cure already
(`Forward<P>` above a break, `Backward<P>` below); Table 7.1's last row and Figure 7.2's caption say "only Ridders'
estimate shows a jump"; the quick reference's pitfall reads "A kink or a range break within a step of x gives a wrong
value from every method, without failing; Ridders' estimate flags a jump, nothing flags a kink at x: differentiate the
branch you mean". 1.5's bullet stays (a kink is listed among the silent failures; the pointer to 7.8 now lands on the
whole truth). Codex's "present these checks as heuristics requiring a smooth, sufficiently resolved function" is
taken as the one sentence of T1. A documentation note for the library is U-15 (section 5b).

### T3. The CMake pin (engineer-25, developer-02, library-07)

`code/getting-started.cmake` pins `GIT_TAG aeec7e6`, the stage-4 integration commit, as a short hash; the front
matter's "The version" still says TODO; CMake asks for a full hash, and a short one fails under `GIT_SHALLOW`; the
library needs CMake 3.25. **Decision: accept, once.** At release the editor pins what "The version" names, and both
name the same thing. Since a commit's hash is not known before the commit, the recommendation to the author (section
5a, item 7) is a tag created at release and pinned by name in both places; failing a tag, the full 40-character hash
of the commit before the release commit, which the headers equal. The snippet's comment becomes "a tag or a full
commit hash; CMake 3.25 or later", on its own line so that every line stays within 60 characters. Until release the
file carries the full hash of `aeec7e6` (`git rev-parse`), so that the shape readers copy is the right one.

### T4. The zero-error residue (evidence-05, numerics-05, evidence-09)

Sidebar 6.2's "(3 on x^5, the rest ... spread over 7 functions)" counts x^5 among the seven; the comment in
`figures/06-ridders-coverage.tex` still says "an exact stop". **Decision: accept, once.** The sentence becomes "the 13
lie on 7 functions, 3 of them on x^5, the rest values that happen to equal the rounded closed form", with the three
existing macros and no new one; the figure comment says "a bit-exact result", as `data/README.md` does (the editor
makes that one-line change in the methods author's figure file, with notice, as at stage 2).

### T5. The AI statement (evidence-02, evidence-07)

The statement gives the set-up to Fable 5.1. The main session verified in git that the set-up commits `7704191` and
`1ffb38f` are co-authored by Claude Opus 5.5 and every commit from the plan `e1965dd` on by Claude Fable 5.1, which
is what R1-decision section 4 said and PLAN.md S6-1 changed. **Decision: accept.** The editor rewrites the second
paragraph's opening: "The work ran as a staged process with AI models in named roles. The book's folder, its layout,
its build and checks and the skeleton of its data generator were set up in a session on Claude Opus 5.5; from the plan
on, a main session on Claude Fable 5.1 coordinated the agents and reported to the author. An editor (Claude Fable
5.1) wrote ..." Evidence-07: "produced by the library itself" gains "apart from a few figures of the library's design
record and tests, which the text attributes" in "The evidence" and in the statement's last paragraph. The author
approves the text, as every version.

### T6. Ridders' paper (evidence-01; R1's T15)

6.1 says "Ridders' paper reduces the step by a factor 1.4 and stops when the error grows". Nobody has seen the paper
(paywall, as `refs.bib` records); the ratio 1.4 and the stop on growth are the constants of the Numerical Recipes
routine `dfridr`, which the design record (`docs/PHASE3_API.md`) compares with; R1's T15 took them for the paper's.
**Decision: modified, the reviewer's second option.** The sentence attributes nothing to the paper beyond the idea:
"Ridders' paper supplies the idea; the rules, the ratio 2, a rounding floor, a two-level stall, a final check and
restarts, are the library's own, chosen on a prototype (Sidebar 6.2), and every statistic in this chapter is of
them". Numerical Recipes is not cited: the book may cite it, but a new entry would need verification and a reference
list that is already full, for a comparison that lives in the design record. For the author, section 5a, item 5.

### T7. The bibliography as printed (evidence-06; R1-evidence-25)

(a) `achemso` turns every `\cite[qualifier]{key}` into its own numbered entry: 40 entries for 23 works, 17 of them
pointers ("Ref. 20, ch. 23" five times), and a superscript lands on a pointer. (b) `@techreport` prints
`ieee754-2019` without its number or publisher. **Decision: (b) accept** (editor: `@book` with `series = {IEEE Std
754-2019 (Revision of IEEE 754-2008)}`, `publisher` and `address`, which the reviewer tested). **(a) is for the
author** (section 5a, item 1), with three options: (i) the qualifiers move into the prose beside an unqualified
citation ("Burden et al.\cite{burden2016} (Section 4.1)"), a mechanical change at 17 places that the editor makes at
integration, after which the list has 23 entries and every superscript lands on a work; (ii) a different citation
style, a preamble design change; (iii) accept the pointer entries as the journal style's way. The editor recommends
(i): no design change, no page change (the list stays one page), and the reader's expectation of a book is met.
Cline-09 belongs here: `refs.bib`'s comments ask for the section at the citing sentence for `dennis1996` (5.4 at
3.4, chapter 5 at 7.4.1) and `nocedal2006` (8.1 at 3.3), and the three citations carry none since the stage-6 close
(PLAN S6-7, to keep the list on one page). Under (i) the sections go into the prose at those three sentences; under
(ii) or (iii) the qualifiers return to the `\cite` and the list grows by three entries. Either way the entries'
comments and the citations agree again.

### T8. B.1's remainder needs one more derivative (numerics-06, codex-08)

Both reviewers: with f^(N+P) merely continuous, replacing f^(N+P)(xi_k) by f^(N+P)(x) leaves o(h^P), not O(h^(P+1));
Codex's f(t) = sgn(t) |t|^(7/2) is C^3 with a central first difference of h^(5/2) at 0. **Decision: accept, once.** B.1
says "o(h^P), and O(h^(P+1)) when f has one more bounded derivative, which the error model assumes throughout", and
the appendix's opening assumption names that derivative.

### T9. The noise recipe's last step, limits and source (engineer-02, numerics-02; closes engineer-03)

The engineer cannot finish the recipe (no line turns `e.sigma` into the step, nothing says what to do at k = 0, nor
where the header is and whether it may be copied); the numerics reviewer measured that the recipe reports a wrong
level without failure when the spacing is on f's own scale (sin 50x at h = 0.1: 2.2 |f|) and nothing (k = 0) when the
spacing is too small for the smooth part to fall below the noise within eight points, which the example checks and
the chapter does not say. **Decision: accept both.** The `noise-estimate-use` snippet gains the two lines `double
const eta = e.sigma / std::abs(f(Re));` and `double const h_eta = std::cbrt(eta) * Re;` after its `k` check (the
`if (e.k > 0)` moves inside the markers as the guard the text then explains); three sentences follow the listing: the
recipe reports k = 0 and no value when no three orders agree, which happens when the spacing is so small that the
smooth part has not dropped below the noise within the eight points, so retry with ten times the spacing or fall back
to the solver's tolerance; a spacing on the scale of f itself can make the orders agree by accident and report a
wrong level, which More and Wild's full procedure guards against with a sign test and a spacing adjustment, so use
the default step of `Central<2>`, as here, and check k; and the header is `docs/deriv-book/code/noise.hpp` of the
repository, under its MIT licence, to copy. (The numerics reviewer's optional note that sigma is a root mean square
while eta is a half-width is already in the text's "close to the root mean square eta / sqrt(3)".)

### T10. The edge rule's two slips (numerics-07a, evidence-08)

6.4's "(derived from the initial step and the level cap, and measured here)" is true of 2^-34 for f' only, and
"2^-30" is typed where 7.8, 10.4 and the quick reference use `\valRiddersEdgeTenDigitsExponent`. **Decision: accept,
once**: "(derived from the initial step and the level cap; for f' also measured here, Section 10.4)" and the macro.

## 2. Decisions by finding

### 2.1 R2-engineer (4: 4 accept)

| Id | Decision | Where and what |
|---|---|---|
| 01 | accept | T1 |
| 02 | accept | T9 |
| 03 | accept | Editor, at integration: Table 1.1 anchored before the "Did it work?" paragraph (or after the list) so that the three silent failures stay together on one page |
| 04 | accept | 5.3: "yet on test functions that vary faster than s, half the default step was better for every order (Sidebar 3.2)"; the repository path stays in the sidebar |

### 2.2 R2-developer (3: 3 accept)

| Id | Decision | Where and what |
|---|---|---|
| 01 | accept | T1: `res->abs_error` tested before it is read; a stencil result fails the check |
| 02 | accept | T3 |
| 03 | accept | Quick reference, "Your own step": "(`Central<2>`, f'; in general `root_epsilon<P + N, T>()`)" |

### 2.3 R2-numerics (7: 5 accept, 2 modified)

| Id | Decision | Where and what |
|---|---|---|
| 01 | modified | T1: the absolute-plus-relative form rather than the scale form |
| 02 | accept | T9 |
| 03 | modified | Sidebar 6.2: "a smaller initial step usually restores the margin", and the claim is generated rather than tagged: the generator reruns the 13 uncovered cases (9 for f', 4 for f'') at a quarter and a sixty-fourth of the default initial step and writes how many are then covered (`\valRiddersMarginQuarterCovered`, `\valRiddersMarginSixtyFourthCovered`; the reviewer's own run gives 10 and 13), so the sentence reads "on the study's 13 uncovered cases a quarter of the default step covers 10, a sixty-fourth all" from the data. For the author (5a, item 4): a small DS9 extension for one sentence |
| 04 | accept | 6.3: "which the method survives by restarting, or, where the samples stay finite, by taking more levels (Sidebar 6.2)" |
| 05 | accept | T4 |
| 06 | accept | T8 |
| 07 | accept | (a) T10; (b) editor: the two macros (`\valRiddersEdgeLimitRelEstimate`, `\valExpRiddersFirstRelError`) enter PLAN.md 4.4 and `data/README.md` |

### 2.4 R2-library (7: 7 accept)

| Id | Decision | Where and what |
|---|---|---|
| 01 | accept | T1 (its first option) |
| 02 | accept | 7.9, one sentence: "return the inner value, or a NaN when the inner call fails (`d ? d->value : std::numeric_limits<double>::quiet_NaN()`), so that the outer call reports `non_finite`" |
| 03 | accept | Editor: `check.py` joins the `\code` lines of a declaration up to the line ending in `;` or `};` and looks the normalised declaration up in the normalised header, so that a change confined to one overload is caught |
| 04 | accept | 7.1 and 8.2: "a double step is diagnosed under `-Wconversion` unless `float` holds it exactly" |
| 05 | accept | `code/03-rounding-step.cpp` `explicit-step`: the two-line `if` braced (closes R1-library-16) |
| 06 | accept | Editor, at integration: Table 8.2 anchored before the last paragraph of 8.5 (or `[t]`), so that it no longer splits a sentence of 8.6 |
| 07 | accept | T3, with "CMake 3.25 or later" in the comment |

### 2.5 R2-evidence (9: 7 accept, 1 modified, 1 for the author)

| Id | Decision | Where and what |
|---|---|---|
| 01 | modified | T6 (the second option, no Numerical Recipes citation) |
| 02 | accept | T5 |
| 03 | accept | 6.5, 6.6, 7.3 and 10.3 say the same thing: "about three more digits than `Central<2>` (14 against 11), as many as `Central<4>`, none over `Central<6>`"; 7.3 drops "about two more digits" for `Central<4>` in favour of the medians it quotes |
| 04 | accept | 7.4.2: "one `derivative` call per entry, m x n calls each costing two model evaluations with `Central<2>`: \valJacobianBalanceCallEvals model evaluations here, against \valJacobianBalanceHandEvals by hand" (the editor's close-out wording "m x n calls, 8 here" made the slip worse; it goes with it). Cline-12 is the same finding |
| 05 | accept | T4 |
| 06 | for the author (a); accept (b) | T7 |
| 07 | accept | T5 |
| 08 | accept | T10 |
| 09 | accept | T4 (editor) |

### 2.6 R2-codex (10: 9 accept, 1 modified; all reproduced)

Reproduction (scratchpad `r2codex/repro.cpp`, GCC 16 against `include/numerixx/deriv.hpp` at the working tree,
whose headers equal `5432abf`): codex-01 as quoted in T2 (value 1.000000000067909, `abs_error` 2.1500912161798169e-10,
passes the relative check at 1e-8; `Central<2>` returns exactly 1 at h = 0.125, 0.25, 0.5; Ridders 1 with
2.2e-16); codex-03: `Central<4>` on t^5 at 0 with h = 0.5 returns -0.25 = -4 h^4, the sign of M_5 = -4, where the
prose's +C h^P f^(5) predicts +0.25; codex-06: at x = -1, h = 0x1.8p-54 (x + h moves, x - h does not) `Central<2>` for
f' returns 1, `Central<4>` for f' and `Central<2>` for f'' return `non_finite`, as `deriv.hpp:1028` checks both sides
for every rule but the two-node f' rules; codex-05: `divide_by_powers` (`deriv.hpp:733`) divides by h first and by dh
last when dh > 1 > h, the branch the regression of `tests/test_orders.cpp:289` exists for (`Central<10>`, N = 4,
h = 19 2^-22: dividing denorm_min by 453600 h = 2.05 first underflows to zero), observable only at the subnormal
boundary, where the test places it; codex-02 by algebra: at the minimum of a h^P + b / h^N the truncation term is
N/P times the rounding term, so E(h*) = (1 + N/P) R, which is 3R/2 for N = 1, P = 2 as Chapter 3's special case has
it, not the (1 + P/N) R of B.4; codex-04 from the weights and `data/kink.csv`: `Central<4>`'s jump contribution is
J/h times the sum of the weights on one side of the break, 7/12 for 0 < delta < h and -1/12 for h < delta < 2h, which
with J = 0.77 J/mol, h = 0.37 K and c_p = 29.58 gives 0.041 and 0.0059, the CSV's 0.04092 at delta = 0.28 and 0.005843
at delta = 0.5; codex-07 from `tests/test_ridders.cpp:313` and `test_orders.cpp:632`, whose multiprecision estimate
ceilings are "twice the first C++ run"; codex-08 and -09 by the argument (B.1's continuity gives o(h^P); Sterbenz on x
and x_- needs h <= |x|/2, and x = 1, h = 0.75 gives x_- = 0.25); codex-10 by reading. Nothing failed to reproduce, so
nothing is rejected.

| Id | Decision | Where and what |
|---|---|---|
| 01 | modified | T2 (and T1's sentence on what the check says); U-15 |
| 02 | accept | B.4: "the two terms stand in the ratio N : P, so E(h*) = (1 + N/P) S eps |f| / h*^N"; the editor checks the sentence's two other uses of the ratio. Cline-01 is the same finding, with the generator (`generate_data.cpp:1499`) and Chapter 3's 3R/2 as witnesses |
| 03 | accept | 4.1, 2.3 and the glossary: "the leading error term is +- C h^P f^(N+P), its sign that of M_(N+P), which the tables give; C is its magnitude" (methods, fundamentals, editor) |
| 04 | accept | Figure 7.2's caption: "a plateau of J/h times the sum of the weights on one side of the break: J/(2h) for `Central<2>`, 7J/(12h) for `Central<4>` while the break lies within h of x and -J/(12h) between h and 2h, the two levels the figure shows; so it falls as 1/h across the rules, whatever their order" (the macros quoted stay; `\valKinkInsideCentralFourRelError` is 7J/(12 h c_p), as R1-evidence-16 checked). Cline-14 is the same finding, with the same 0.04095 and the ratio 105 against 122 |
| 05 | accept | Sidebar 5.2 and 8.5: "divides by dh and then by h, N - 1 times, in the order that avoids a spurious underflow or overflow (by h first when dh > 1 > h)". No library item: the header's own comment (`deriv.hpp:1041`) states the ordering; the book had omitted the clause |
| 06 | accept | A.4: "if only one moves, as can happen near a power of two, a two-node rule for f' returns a finite value that is mostly rounding, and every other rule reports `non_finite` (Sidebar 3.3, Section 8.2)"; U-13 |
| 07 | accept | 10.6: "Their reference values for the derivatives come from closed forms, ...; the ceilings on `abs_error` for the multiprecision types are twice the estimates of the first run, regression baselines that `docs/design.md` records as such, not independent references"; U-14. Cline-05 is the same finding |
| 08 | accept | T8 |
| 09 | accept | Sidebar 5.2: "for 0 < h <= |x| / 2 by Sterbenz's lemma on x and x_+-; for |x| / 2 < h <= |x| the subtraction x - h is itself exact by the same lemma, so x_- = x - h and the offset is h, a representable number; at h = |x|, x_- = 0". The conclusion stands; the argument covers the whole interval |
| 10 | accept | Table 9.4's caption: "every numerator fits in 24 bits, so every rule's integer weights are exact in `float`". Cline-06 is the same finding |

### 2.7 R2-cline (21: 18 accept, 2 modified, 1 reject; every claim reproduced)

Ids, in order of appearance (the pass over the appendices, Chapters 9 and 10, the back matter, `data/README.md` and
`refs.bib` first, then the pass over Chapters 7 and 8; the model's working notes between them carry no ids):
01 a2:152-153 the multiplier; 02 a2:178-179 "the largest of the library's rules"; 03 a2:180-181 "1.64 is a rounding
of kappa"; 04 a1:73 "up to ulp(x)/2"; 05 10:213-216 "never from Numerixx itself"; 06 09:78-79 "every rule is exact in
float"; 07 99:29 "spans K h on each side"; 08 `data/README.md:13-14` the coverage decimals; 09 `refs.bib:200, 212` the
sections at the sentence (its pass-2 item on 07:156 is the same point); 10 `refs.bib` `nist-webbook` without a year;
11 the cross-scope note on `deriv.hpp:385` "12! times 4!"; 12 07:186-187 "m x n calls"; 13 07:33 and 08:78 "a member
function of T"; 14 07:363-366 Figure 7.2's plateau; 15 08:270 "T is always deduced from x"; 16 07:286 "-> mp50
matters"; 17 07:71-72 `d.value_or`; 18 07:260-261 "to within it"; 19 07:29-30 "a wider step"; 20 08:248 the members
the concept does not check; 21 08:39-40 "reads `detail` only to generate the tables". A remark in the working notes
(Table 7.1's fourth row offers Ridders but counts 2 evaluations) is taken without an id.

Reproduction (scratchpad `r2codex/cline.py`, `cline.cpp`, `cline_et.cpp`, GCC 16, Boost.Multiprecision from the CPM
cache): cline-02 by the generator's own formula with the constants of Tables 9.1 and 9.4: the cost
(N kappa^-P + P kappa^N) / (P + N) is 15.135 for `Central<4>` and f'''' but 18.151 for `Central<10>` and f''''
(S = 44.012, kappa = 2.2452), so "the largest of the library's rules" is false (the generator computes the cost for
the five rules the book quotes only; `Central<6>` and `<8>` for f'''' give 7.81 and 9.77, below 15.1, so Cline's
16.88 and 17.72 for them do not reproduce, but its headline case does); cline-03: the cost for `Central<10>` and f'
is 1.634 at kappa = 1.7975 and 1.637 at kappa = 1.8, so the header's "at most a factor 1.64" (`deriv.hpp:521`, which
never names `Central<10>`) is the cost at the rounded kappa, not "a rounding of kappa"; cline-04: x = 1 - 2^-53,
h = 2^-52: x + h = 1 + 2^-53 is a tie and rounds to 1, an error of 2^-53 = ulp(x), twice the stated bound; cline-13:
`detail::callable<R (R::*)() const, R>` holds and `callable<R (R::*)(), R>` does not, since `fn` is invoked with a
`T const` lvalue; cline-15: `derivative<Central<2>, 1, float>(sq, 1.0)` compiles and returns 2 in `float`; cline-16:
a lambda returning `exp(x)` for `number<cpp_bin_float<50>, et_on>` without `-> mpet` compiles and gives e, because
the library converts every result at once (`deriv.hpp:804, 865, 963, 1031, 1200`); a lambda returning an expression
over its own local compiled and happened to give the right value too, which is undefined behaviour, not safety;
cline-01, -05, -06, -12 and -14 as in section 2.6; cline-07, -09, -10, -11, -17 to -21 by reading the lines named.
Cline-18 holds as the book says it (|d - r| differs from the stencil's actual error by Ridders' own error, 2.0e-12,
below its 2.7e-12 estimate), and only the two-digit printing hides it.

| Id | Decision | Where and what |
|---|---|---|
| 01 | accept | codex-02 |
| 02 | accept | B.4: "the largest of the four rules the design record quotes" in place of "the largest of the library's rules" (`Central<10>` for f'''' costs 18 by the same formula, which nothing quotes) |
| 03 | accept | B.4: "the 'at most a factor 1.64' of the library's documentation for f' is the cost at kappa = 1.8 rounded (`Central<10>`, where the exact kappa 1.7975 gives 1.63)". The sentence's slip was the R1 letter's own wording (library-17, evidence-21) |
| 04 | accept | A.4: "differs from x + h by up to half the spacing at the result: ulp(x)/2, or ulp(x) when x + h crosses a power of two (the two-sided spacing of Sidebar 3.3)" |
| 05 | accept | codex-07 |
| 06 | accept | codex-10 |
| 07 | accept | The notation table: "K: the largest |n_k| of a rule; the rule reaches at most K h from x" (editor); 4.3 has the same slip, "which spans K h on each side of x", and reads "which reaches K h from x on the sides it samples" (methods) |
| 08 | accept | `data/README.md`: the convention gains its exception, "coverage with one decimal, two where one would read 100.0 with a case uncovered" (editor) |
| 09 | accept | T7: the sections at the three sentences, in the prose under option (i) or in the `\cite` otherwise |
| 10 | reject | `nist-webbook` has no `year` because NIST's recommended citation of SRD 69 has none (a continuously updated database); the retrieval date in the note is its date, as the ACS style prints web resources. The entry stays |
| 11 | accept (library) | U-16: the header's "the largest intermediate is 12! times 4!" is a bound the generator never reaches (4.1 and R1-numerics-17b have the right figure, 12!); the book is unchanged |
| 12 | accept | evidence-04 |
| 13 | accept | 7.1's margin note and 8.2: "a pointer to a const member function of T itself" (a non-const one is not invocable on the `T const` lvalue the library passes; the library's own test uses `constexpr Real square() const`) |
| 14 | accept | codex-04 |
| 15 | accept | 8.8: "T is deduced from x, never from h, unless it is given as the third template argument" |
| 16 | modified | 7.6: "give a lambda an explicit `-> mp50`: the library converts a result at once, so an expression over the parameter alone is safe, but one over the lambda's own locals would refer to destroyed values (Sidebar 7.1); the explicit type is the habit that is always safe". The sidebar's explanation stands |
| 17 | accept | 7.2: `r.value_or(0.0)` and `r->abs_error.value_or(0.0)`, the listing's name |
| 18 | modified | 7.5: "is the stencil's error to within Ridders' own error, \valCpRiddersActual here, below its estimate" so that the claim can be checked from the numbers quoted (the two-digit macros hide a margin of 2.0e-12) |
| 19 | accept | 7.1: "a step of another type converts to it" |
| 20 | accept | 8.7: "the algorithms also read `radix`, `digits`, `epsilon()`, `max()` and `denorm_min()`, which the concept does not check", as Sidebar 7.1 and the listing already say |
| 21 | accept | 8.1: "the book names it only where the declarations of 8.2 do and to generate the tables of Chapter 9" |
| -- | accept | Table 7.1, row 4: evaluations "2, or Ridders' 14" (the remark in the working notes) |

## 3. Revision plan by chapter file (stage 7, final fixes)

Each author writes `responses/R2-<chapter-id>.md` as for round 1. The editor's items (section 4) that change macros
or snippets come first; the editor circulates the two new macro names before the methods author starts.

| File | Owner | Items | Lines |
|---|---|---|---|
| `01-why.tex` | fundamentals | 1.5: the two-tolerance check with the one-sentence guide (T1); Table 1.1's anchor is the editor's (engineer-03) | +2 |
| `02-taylor.tex` | fundamentals | 2.3: the sign of the leading term (codex-03) | 0 |
| `03-rounding.tex`, `03-rounding-noise-estimate.cpp` | fundamentals | 3.6: the two snippet lines, the k check inside the markers, the three sentences, the header's path and licence (T9); `explicit-step` braced (library-05). Paid for in the chapter: the "Ridders' method does not help here" paragraph and 3.7 lose a line each; the editor trims the noise header's comments (section 4) | +9, -8 |
| `04-stencils.tex` | methods | 4.1: the sign of the leading term (codex-03); 4.3: "reaches K h from x on the sides it samples" (cline-07) | 0 |
| `05-orders.tex` | methods | 5.3: the design-record citation to Sidebar 3.2 (engineer-04); Sidebar 5.2: the Sterbenz argument over the whole interval (codex-09), the conditional division order (codex-05) | +3 |
| `06-ridders.tex` | methods | 6.1: the paper against the library without the paper's constants (T6); 6.3: "or by taking more levels" (numerics-04); 6.4: the kink sentence (T2), the edge rule's slips (T10); 6.5 and 6.6: "about three more digits" (evidence-03); Sidebar 6.2: the zero-error wording (T4), the margin sentence from the new macros (numerics-03) | +4 |
| `07-using.tex`, `07-using-expected.cpp` | practice | 7.1: the `float` step wording (library-04), "a step of another type" (cline-19), "a const member function" in the margin note (cline-13); 7.2: the two-tolerance helper and its guide (T1), `r` for `d` (cline-17); Table 7.1: row 4's evaluations (the Cline remark), the last row (T2); 7.4.2: the call count (evidence-04, cline-12); 7.5: a jump, not a kink (T2), "to within Ridders' own error" (cline-18); 7.6: the `-> mp50` sentence (cline-16); 7.8: the kink's missing symptom, Figure 7.2's caption (T2, codex-04, cline-14); 7.9: the inner failure as NaN (library-02). Paid for in the chapter: Sidebar 7.1's second paragraph (-3), Table 7.1's caption (-1), 7.4.3's fixed-point paragraph (-2), the two trailing comment lines of `estimate-pitfalls` (-2), 7.9's mixed-partials paragraph (-1), 7.9's deferred-features paragraph (-1) | +11, -10 |
| `08-api.tex` | practice | 8.1: where the book names `detail` (cline-21); 8.2: the `float` step wording (library-04), "a const member function" (cline-13); 8.5: the conditional division order (codex-05); 8.7: `epsilon()` and `max()` in the list (cline-20); 8.8: T given third (cline-15); Table 8.2's anchor is the editor's (library-06) | +1 |
| `09-tables.tex` | practice | Table 9.4's caption (codex-10) | 0 |
| `10-guarantees.tex` | practice | 10.3: "about three more digits" (evidence-03); 10.6: the estimate ceilings (codex-07). Paid for: 10.7's list loses two lines so that Table 10.4 stays on its page | +2, -2 |
| `a1-floating-point.tex` | fundamentals | A.4: the two resolution rules (codex-06); the rounding of x + h at a power of two (cline-04) | +2 |
| `a2-derivations.tex` | methods | B.1: the remainder's regularity (T8); B.4: the multiplier 1 + N/P (codex-02, cline-01), "the largest of the four the design record quotes" (cline-02), the 1.64 sentence (cline-03) | +1 |
| `99-backmatter.tex` | editor | The notation table's K (cline-07); the glossary's C (codex-03) | 0 |
| `qr-quick-reference.tex` | practice | "Two checks": the two-tolerance form and "a stencil has no estimate to test" (T1); the kink pitfall (T2); "Your own step" (developer-03). Paid for: "Where to look" cut to two lines, the Ridders bullet loses its last clause; the page stays one page | +3, -3 |

## 4. The editor's items

1. **The AI statement** (T5) and "The evidence" (evidence-07); the author approves the text.
2. **`refs.bib`**: `ieee754-2019` as `@book` with series, publisher and address (T7b). If the author chooses T7's
   option (i), the 17 qualified citations move their qualifier into the prose at integration, in every chapter, as a
   citation-format change; the reference numbers change with the collapse to 23 entries.
3. **The generator and data**: the DS9 extension for numerics-03 (`\valRiddersMarginQuarterCovered`,
   `\valRiddersMarginSixtyFourthCovered`), if approved; PLAN.md 4.4 and `data/README.md` entries for them and for the
   two stage-6 macros (numerics-07b); the README's coverage convention with its exception (cline-08); the comment of
   `figures/06-ridders-coverage.tex` (T4); the comments of `code/noise.hpp` trimmed by five lines to pay for T9 (the
   header is the editor's).
4. **`tools/check.py`**: declarations compared whole (library-03).
5. **`code/getting-started.cmake`**: the full hash now, the tag or release hash at release, the comment on its own
   line (T3).
6. **Layout at integration**: Table 1.1 (engineer-03) and Table 8.2 (library-06) anchored; the `\needspace` guards
   re-checked after the chapters' changes; the quick reference verified as one page.
7. **The glossary**: the entry for C (codex-03).
8. **`.gitignore`**: `tools/__pycache__/` added (done with this letter), so that `check.py`'s cache never enters a
   run or a commit.
9. **PLAN.md**: 12.3 gains U-13 to U-15 (section 5b); 12.5 records the round-2 decisions when they are applied; 4.4
   and 5 gain the macro and snippet changes above.

## 5. For the author

### (a) Contested or scope-changing points, with the editor's recommendation

1. **The bibliography** (T7a): (i) qualifiers in the prose, 23 entries, no design change (recommended); (ii) a
   different citation style, a preamble change; (iii) the pointer entries as they are.
2. **The acceptance check in two tolerances** (T1): `atol + rtol * |value|` in 1.5, 7.2 and the quick reference, with
   the one-sentence guide. Recommend accept: it is the standard form, it is right near a zero derivative, and it is
   safe on a stencil result.
3. **What the checks cannot see** (T2): the quick reference and Table 7.1 lose "only Ridders' estimate flags it" as a
   general statement, in favour of "a jump; nothing flags a kink at x". Recommend accept: it is the arithmetic, and
   the cure the book already gives (know the branch boundaries) is the only one.
4. **A DS9 extension for one sentence** (numerics-03): the 13 uncovered cases rerun at two smaller initial steps so
   that "a smaller initial step usually restores the margin" quotes measured counts. Recommend accept (a few lines in
   the generator, seconds of run time); the alternative is the sentence tagged "as a rule", without numbers.
5. **Ridders' paper** (T6): the constants 1.4 and "stop on growth" are dropped from 6.1 rather than attributed to
   Numerical Recipes with a new citation. Recommend the drop; the comparison with `dfridr` stays in the design record.
6. **The AI statement** (T5): the set-up on Opus 5.5, the coordination from the plan on by Fable 5.1, as git records;
   the text needs approval.
7. **The pin at release** (T3): a tag created at release and named in both "The version" and the CMake snippet
   (recommended; the author names it), or the full hash of the commit before the release commit.
8. **The repository's test rule** (codex-07, U-14): the tests' multiprecision estimate ceilings come from the first
   run, an exception to "never from Numerixx itself" that `docs/design.md` documents but the rule in the repository's
   CLAUDE.md does not; the author decides whether the rule says so.
9. **A documentation note on kinks** (codex-01, U-15): recommend a one-line note in the header's documentation of the
   smoothness assumption; no code change is possible or proposed.

### (b) Library issues outside the book, to add to PLAN.md 12.3

- **U-13** (codex-06): the public documentation of `error::non_finite` says "h is below the resolution of x" for
  every rule; the two-node f' rules fail only when both points coincide, every other rule when either side does not
  move (`deriv.hpp:1028`, an internal comment). A documentation note, no code change.
- **U-14** (codex-07): `tests/test_ridders.cpp` and `tests/test_orders.cpp` bound the multiprecision types' `abs_error`
  by twice the first run's estimates, regression baselines rather than independent references; the tests' comments and
  `docs/design.md` say so, the repository's test rule does not.
- **U-15** (codex-01): a kink at x, or so close that no level resolves it, is invisible to every rule and to Ridders'
  estimate (the central difference returns the mean of the one-sided slopes at every step). A documentation note on
  the smoothness assumption of `derivative`.
- **U-16** (cline-11): the header's comment at `deriv.hpp:385`, "the largest intermediate is 12! times 4!
  (`Central<10>`, n = 4)", names a bound the generator never reaches; the largest intermediate is 12! (the book's
  4.1, verified in round 1 by numerics-17b and again by Cline). A one-word documentation fix.
- codex-05 adds no item: the header's comment at `deriv.hpp:1041` states the conditional order of the divisions; the
  book had omitted the clause.

## 6. The page budget

The build is at 88 pages, the limit. The additions of section 3 sum to about 39 lines, the cuts named there to about
24; the difference falls in chapters with slack at their last page (Chapters 1, 5, 6, 8, Appendices A and B, the
front matter), so the expected total is **88**. The chapters that are full carry their own cuts: Chapter 3 (T9: +9,
-8, the noise header's comments included), Chapter 7 (+11, -10), Chapter 10 (+2, -2) and the quick reference (+3, -3,
one page by rule). The reference list stays one page under every option of T7. If the stage-7 integration build is
89, the next cuts, in order, are the `estimate-pitfalls` listing of 7.2 reduced to its four code lines (-6), 3.8's
checklist items to one line each (-4), and then PLAN.md 1.1's order (B.2's worked example, Chapter 9's tables to
`\footnotesize`, Sidebar 7.1 to half its length).

# R1-decision: the editor's decisions on review round 1

Approved by the author on 2026-10-23: every decision below, the four additions (E21, E22, the complete first program with the CMake snippet, E23), the two checks, rounded headline digits, the pragma move, the three preamble changes, captions keeping their data paths, and the rewritten AI statement to be approved with the stage-6 build. Review round 2 adds a Cline review with DeepSeek V4.1 Flash 1017 beside the Codex review (the author's choice of model, 2026-10-24; the decision letter had named MiMo-V2.6-Pro).

Editor: `book-editor` (Fable 5.1), 2026-10-23. Object: the five round-1 reviews of the integrated book at `aeec7e6`
(83 pages): `R1-engineer.md` (25 findings), `R1-developer.md` (18), `R1-numerics.md` (17, of which 17 bundles eleven
slips), `R1-library.md` (22), `R1-evidence.md` (29, plus the reference table with two corrected entries and two
proposed additions). Every finding is decided below by its id: **accept**, **modified** (accepted in a modified form,
with the form stated), **reject** (with the reason) or **defer** (with why). Where two or more reviewers raise the
same point it is decided once, in section 1, and the per-finding rows point there.

Counts: engineer 13 accept, 12 modified, 0 reject, 0 defer; developer 14, 4, 0, 0; numerics 15, 2, 0, 0; library 20,
2, 0, 0; evidence 23, 5, 0, 1 (the licence question goes to the author). No finding is rejected outright; the parts
rejected inside a modified acceptance are named in the rows (the data-file paths in the captions, `return 1;` in the
listings, a generated accuracy for mixed partials, the `check.py` extension to integers).

Nothing is revised yet (stage 6, by the authors, per `responses/README.md`). Section 3 is the revision plan by chapter
file; section 4 the editor's own items; section 5 what needs the author; section 6 the page budget.

## 1. Recurring themes, decided once

### T1. "The one check" is not enough (engineer-01, developer-02, developer-15, library-02; also engineer-09)

Three reviewers converge: section 1.5 and the quick reference call `if (d)` "the one check", while the book itself
shows silent wrong values (the range break, the narrow feature, the solver output) and, for Ridders, a value that is a
fallback with a large `abs_error`. **Decision: accept.** Section 1.5 splits "The one check" into two short paragraphs:
*Did it work* (`if (d)`, with the `else` that reports the error, developer-15) and *Is it right* (the three silent
failures in one line each with their pointers, 3.5, 3.6, 7.8, and "if in doubt, repeat at 2h and h/2", 7.5). For
Ridders the book teaches a second test, `*r->abs_error <= tol * abs(r->value)` before the value is used where it
matters, as a compiled snippet in 7.2 (E11 extended) and a pitfall in 7.8 ("a Ridders result with a value can still be
poor; read its estimate"), and the quick reference carries both checks. The first listing of 1.4 gains the `else`.
The pitfalls listing of 7.2 moves its reads inside `if (r)` (library-02). One house idiom in Part I (engineer-09):
`d->value`, and for the estimate `*r->abs_error` with a margin gloss at its first use in 1.4 ("`abs_error` is a
`std::optional`; `*` reads it"); the implicit conversion is shown only in 7.2 and Chapter 8, and 7.5's
`digits(*d2, *d)` becomes `digits(d2->value, d->value)`.

### T2. The Jacobian is promised, not delivered (engineer-02, developer-14; engineer-10 adjacent)

Chapter 1 opens with "a Jacobian is the slope of a set of balances"; 7.4.2 computes a 1 x 2 gradient of a scalar and
sends the reader to the source file for the loop; 7.9 then defers Jacobians. **Decision: modified, for the author's
approval (scope: a new example).** 7.4.2's loop is renamed "parameter sensitivities" and printed (developer-14). A new
example E22 (`code/07-using-jacobian.cpp`, practice author) shows a genuine 2 x 2 Jacobian of two residuals (an energy
and a mass balance with closed-form partials as the answer key), in the two ways the engineer names, and the text
states plainly which the library supports: (a) one `derivative` call per entry, which the library supports today and
which costs m x n x (evaluations) calls of the whole model; (b) the central stencil applied by hand to the vector
residual at x_j +- h (two model calls per column, h = `root_epsilon<3, double>() * max(|x_j|, 1)`), which is the
formula of Chapter 2 without the library's checks and is what one writes when the model is expensive. The listing is
(b), a six-line loop, with (a) stated in one sentence and its cost; 7.9's "a loop of calls does the job meanwhile"
points to it, and "several outputs from shared evaluations" stays named as deferred. Cost: about 0.6 page. Mixed
partials and caching (engineer-10): modified, a short paragraph in 7.9 with the nested-call approach, its cost (the
inner call's evaluations times the outer's) and the model's rule of thumb for its accuracy (the inner error about
eps^(2/3) acts as the noise level of the outer call, so about eps^(4/9), roughly 10^-7 in double); the generated
accuracy value the reviewer asks for is **not** added (it needs a data set for one paragraph, and the statement is a
rule of thumb by construction); one sentence on caching an expensive call so that several outputs share evaluations.

### T3. The solver-noise recipe (engineer-03, engineer-08; numerics-13, numerics-17g, evidence-06)

The cure for a noisy f is clear, but the reader cannot obtain eta when the solver's tolerance is on a residual rather
than on the output, and the fallback (a difference table) is one paragraph without code. **Decision: modified, for the
author's approval (scope: a new example and data set).** Section 3.6 gains (i) one sentence that a residual tolerance
is not eta and that eta should be measured when the convergence test is not on the output; (ii) a compiled recipe E21
(`code/03-rounding-noise-estimate.cpp`, fundamentals author, with a data set DS22 from the editor): evaluate f at
eight points spaced on the scale of the noise's effect, form the difference table, and take the k-th differences' root
mean square divided by sqrt((2k)!)/k! as the estimate of eta |f| (More and Wild's ECnoise rule, cited as `more2011`,
evidence-06, which is the paper that holds the procedure; `more2012` stays at the step), run on the friction factor of
7.4.3 so that the text can say that it recovers about 10^-10 at a stated cost in evaluations, as generated values; and
(iii) "estimates eta |f| up to a known factor" (numerics-17g). Cancellation inside f (numerics-13): accept, one
sentence in 3.6 ("if f is computed by cancellation its error is eps times the largest intermediate, k eps |f|; treat
k eps as eta"). "The solver matters" (engineer-08): accept; the paragraph leads with the conclusion (a solver can be
smooth but wrong by its tolerance, and no estimate sees it), says that the estimate did not cover the error here, and
gives the action (converge the inner solver well below the derivative accuracy wanted, keep the iteration count fixed
across the stencil's points if the solver allows, or differentiate the converged equations implicitly). Cost: about
0.7 page in Chapter 3.

### T4. A complete first program and how to get the library (engineer-04, developer-11)

**Decision: modified, for the author's approval (one exception to the listing rule).** The first listing of 1.4
becomes a complete program (includes, `main`, the lambda, the call, `if (d) ... else ...`, print), quoted from E1 as
today; one sentence says that `book::` is the namespace of the book's own example header. Section 7.1 gains a short
"Getting the library" paragraph: there is no install rule and no `find_package` config, so `add_subdirectory`, CPM or
FetchContent is the way in; a CMake snippet (`FetchContent_Declare` with a pinned tag, `FetchContent_MakeAvailable`,
`target_link_libraries(app PRIVATE numerixx::numerixx)`), quoted with `\lstinputlisting` from a file
`code/getting-started.cmake` that the book's build does not configure; this is the one listing in the book that is
not a compiled example, and PLAN.md section 8 records it as such if the author agrees. The toolchains tested are
named with their versions from `.github/workflows/ci.yml` (GCC 16, Clang 22 with libc++, MSVC and clang-cl of Visual
Studio 2026, Emscripten at the pinned `EMSDK_VERSION`), as "the toolchains the tests run on", not as minima. Cost:
about 0.5 page.

### T5. The symbol s changes meaning in Chapter 3 (engineer-06)

**Decision: modified.** In 3.3 the derivation uses w for the scale on which f varies (h* about 3^(1/3) eps^(1/3) w);
3.4 introduces s = max(|x|, 1) as the library's guess for w; 3.5 then needs no new symbol. The later chapters' "the
default step assumes that f varies on the scale s" is the assumption, correctly worded, and stays; Sidebar 3.2, 4.3,
5.2, 5.3 and B.4 are checked for the same reading. The notation table's entry for w becomes "the scale on which f
varies; s is the library's guess for it".

### T6. Developer and test detail in the Part I text (engineer-07, engineer-19, engineer-20; library-08, evidence-15)

**Decision: modified.** 4.1's integer-numerator storage, the binary-fraction test and "Central<4> on x^4 returns
exactly 4" move into Sidebar 4.1; 4.3 keeps the practical statement (the third rounding term from the points) and the
test budget moves to 10.6, where it is stated as the tests state it: for a selection of rules and orders, two-node
rules included (library-08), and "the tests' comment reports 4 to 16 on its three points" (evidence-15: attributed,
not generated, because the book's DS7 medians give a different band and the sentence is about the tests). Sidebar 3.1
writes "a fixed multiple of h" instead of "dh" (d is defined in Chapter 4). "As the header warns" and "the header's
own example" become plain statements (engineer-19); "The example checks ..." sentences are deleted where they carry no
lesson (engineer-20: 2.3, 4.4, 5.6). The data-file path at the end of each caption **stays** (the rejected part of
engineer-19): it is the evidence trail the book promises in its front matter (PLAN.md 4.1), and a list of figures in
the back matter would cost a page for the same information.

### T7. The histogram tail and the "exact stops" (numerics-02, evidence-02; numerics-03, evidence-04)

Both reviewers recomputed the same two misreadings of DS9. **Decision: accept both.** The tail of the N = 2 histogram
(19 and 21 evaluations, 332 cases) is the slower convergence of e^(-x^2) beyond x = 2 and of 1/x, log x, sqrt(x) and
the Runge function at small x; the 11 restarts all take 21 evaluations and are a small part of the last bar. The
generator writes the tail's size, its share of restarts and the functions in it as macros (DS9 extension), the
caption and 6.5 quote them, and the "restarts" label in `figures/06-ridders-cost.tex` goes. The "13 and 8 exact
stops, all on x^5" become "cases with an actual error of exactly 0": for f'' all on x^5, for f' three on x^5 and the
rest bit-exact coincidences on six other functions; the macros are renamed `\valRiddersZeroError...` and the
generator adds the split (on x^5, number of functions). Exact stops are **not** claimed anywhere: the generator does
not observe the stop reason, and logging it would need a change to the library's interface. The caption of Figure 6.3
says "an actual error of exactly 0".

### T8. The sign in Eq. (4.3) and Appendix B.1 (numerics-01, numerics-04)

**Decision: accept both.** The remainder of Eq. (4.3) is + (h^4/30) f^(5)(xi): the chapter's own M_5 = -4, Table
9.1's sign column, the test on x^5 and Burden all say so. B.1's single-point argument needs the products w_k n_k^(N+P)
of one sign, which holds for Eq. (5.1) but not for the central rules for f' above order 2; B.1 states that the Taylor
argument gives the asymptotic form in every case, that the single-xi form of (4.3) follows from the interpolation
error formula (Burden et al., cited), and names the central rules above order 2 and some one-sided rules as the cases
where the hypothesis fails.

### T9. The cost of misjudging the scale (numerics-06, library-09; 3.5 against Sidebar 5.1 and 5.3)

The two chapters contradict each other and each reviewer proposes a different repair. **Decision: accept, decided once
as follows.** A step off by a factor 2 multiplies the rounding term by 2^N (halved) or the truncation term by 2^P
(doubled); at the model's balance the total error grows by at most (N 2^-P + P 2^N)/(P + N) or (N 2^P + P 2^-N)/(P +
N), whichever is larger: 1.7 for f' and 5.5 for f'''' with Central<2>, 8.0 for Central<4> and f''''. The header's 16
is the growth of the rounding term alone and is named as such. The factors are generated (`\valMisjudgeCost...`, DS13
extension) and quoted in Sidebar 5.1, 5.3 and 3.5, so the three places agree by construction.

### T10. The rationales of the order caps (numerics-07, evidence-05)

**Decision: modified.** The caps are presented as design decisions with the model's figures, tagged as such and
attributed to `docs/PHASE2_API.md`: beyond P = 10 the next central order would gain, by the model's exponent, about
0.2 digit in double and 0.7 digit in a 50-digit type for f' (generated from eps alone, `\valCapGain...`), and more for
the higher derivatives, so the library stops there; beyond P = 6 the one-sided rules' amplification grows about as
fast as the exponent pays, and the library stops there. The one-sided figures the numerics reviewer computed (order 7
gains a seventh of a decade, order 8 nothing) are **not** quoted: they need the weights of rules the library does not
have, and the generator will not grow a second Lagrange solver for two numbers; the sentence is qualitative and
attributed. The header's and `design.md`'s "half a digit" and "would cost more than the order gains" are reported to
the author as U-7.

### T11. "Never throws" and `.value()` (developer-03, developer-07, library-05)

**Decision: accept, once.** 1.3 and 7.2 say "the library itself throws nothing: a failure is a value; an exception
from `fn` or from `T` passes through, and `derivative` is not `noexcept`". `d.value()` (the `std::expected` accessor,
which throws `std::bad_expected_access` on failure, or terminates under `-fno-exceptions`) joins the traps of 7.2, 8.3
and the quick reference, against `d->value`, the estimate's member. Reported to the author as U-8 (a documentation
note for the header).

### T12. Statements loosened or mislabelled, raised twice

- The "factor of ten" of Table 10.1 (numerics-05, evidence-11): accept, once: "within a factor of about six for the
  central rules for f' to f''' and up to six times above the model for f''''; the one-sided rules of higher order lie
  above it by their amplification S, about 35 times for Forward<4> at N = 3 and 4". The ratios are generated (DS15
  extension) so that the sentence follows a regeneration.
- The study's scope "no singularity within h0" (numerics-09, evidence-12): accept, once: "for f'; for f'' to f'''' the
  initial step s/4 (and the span 2 h0 for f''' and f'''') reaches past zero for x below 0.25 (0.5), which are the
  restarts of Table 10.2". Caption of Figure 6.3, 6.3 and 10.2.
- The edge "rule" 2^-34 (numerics-17k, evidence-13): accept, once: both readings named wherever the rule is stated
  (6.4, 7.8, 10.4, quick reference): "a value is returned down to 2^-34 s (derived and measured), with the estimate
  risen to a tenth of the value; for ten digits keep the nearest bad point at least 2^-30 s away (measured for f')".
- Ridders' cost and gain (numerics-17j, evidence-09): accept, once: "two to seven times as many evaluations as a
  stencil for f'; about two more digits than Central<2>, as many as Central<4>, none over Central<6>; and an
  estimate", with the macros of Table 10.1. 6.6 likewise.
- The 24 against 14 evaluations at the edge (numerics-16, evidence-10): accept, once: the comparison is with log x
  without a restart (16 evaluations, DS10's k = 2 row, generated as `\valRiddersEdgeEvalsInterior`), four of the extra
  are the two failed levels and the rest a tableau that starts at 1/128, nearly the size of x, and needs ten levels.

### T13. Numbers in words, and "digits" (numerics-17h, evidence-22, evidence-23)

**Decision: modified.** Counts and levels stated as facts use macros (`\valTableauLevelsUsed`, the level and column of
the returned entry, the caps `\valRiddersMaxEvals...` in 8.2, 8.5, 10.5 and the quick reference). "Digits" are made
consistent: the headline macros become the rounded value of -log10 of the median rather than its floor (so "about 5
digits" for a fourth derivative at 1.2e-5, not 4; Ridders' "about 13" becomes 14 for 1.3e-14), and the model digits of
5.2 are generated to one decimal (10.4, 7.8, 6.3, 5.2) in place of "roughly eleven, eight, six and five". The three
slips (three and a half against 3.8 decades; eleven against ten orders; "below 4.3e-12" against "at most") are
corrected. Extending `check.py` to integers and short decimals is **deferred**: on this book it would flag every
"two evaluations" and "level 3" and bury the real warnings; the editor's audit script lists the words-for-numbers
instead, once, at stage 6.

### T14. The typed declarations of Chapter 8 (library-15, evidence-28)

**Decision: modified.** The `signature` blocks stay (quoting from `include/` would need markers in a public header,
which the library's rules forbid the book to add). `tools/check.py` gains a check that every code line of a
`signature` block occurs, whitespace-normalised, in `deriv.hpp` or `math.hpp`, so that a drift fails the build; PLAN.md
section 8 records the declarations as the one typed exception to the listing rule, checked by `check.py`, and the
front matter's "every code listing is quoted from a compiled program" says "every code listing but the declarations
of Chapter 8, which are transcribed from the header and checked against it at every build".

### T15. The implemented method against Ridders' paper (numerics-11, evidence-17)

**Decision: accept, once.** One sentence in 6.1: "Ridders' paper reduces the step by a factor 1.4 and stops when the
error grows; the library uses the ratio 2, a rounding floor, a two-level stall, a final check and restarts, chosen on
a prototype (Sidebar 6.2): the idea is Ridders', the rules are the library's", and the chapter's statistics are said
to be of the latter. Sidebar 6.1: "Eq. (6.1) generalises to the weighted combination ...; Burden et al. write the
equal correction form, which the library computes" (the uncited "textbook" goes).

## 2. Decisions by finding

### 2.1 R1-engineer (25: 13 accept, 12 modified)

| Id | Decision | Where and what |
|---|---|---|
| 01 | accept | T1. 1.5 "Did it work / Is it right"; Table 1.1's caption glosses w and eta in plain words |
| 02 | modified | T2, for the author. E22 with the hand-applied stencil printed and the library's way stated with its cost; 7.4.2 renamed |
| 03 | modified | T3, for the author. E21 and DS22; the residual-tolerance sentence |
| 04 | modified | T4, for the author. A complete first program; "Getting the library" in 7.1 with a quoted CMake file; `book::` explained; the CI toolchains named with versions, not as minima |
| 05 | modified | A sixth checklist question in 3.8 with the rule of thumb (expect to lose log10(|f| / (s |f'|)) digits; subtract a reference) and one sentence in 3.3 giving the order of magnitude for an enthalpy with a formation term; the generated water case is **not** added (it would bring new thermochemical data to verify for one number; the rule needs none) |
| 06 | modified | T5 |
| 07 | modified | T6 |
| 08 | accept | T3 |
| 09 | modified | T1: the house idiom `d->value`, `*r->abs_error` with a gloss; the conversion only in 7.2 and Chapter 8 |
| 10 | modified | T2: mixed partials and caching in 7.9, the model's rule of thumb, no generated value |
| 11 | accept | 3.5: two sentences with a units example (K against degrees C: s = 20 against 293, a step about 15 times smaller) and its consequence in digits, with "about" |
| 12 | accept | 6.4, 7.4.3 and 7.5 lead with the conclusion; "digits of |f|/s" in plain words |
| 13 | accept | 3.5 "for a smooth interpolant (a cubic spline) it is the node spacing" with a pointer to 7.9; 7.9 adds the spline advice (its own derivative; or keep the stencil within one interval, or accept the kink error of 7.8) |
| 14 | accept | "monadic operations" replaced by the plain description; the C++ terms of Chapter 8 (aggregate, designated initialiser, structured binding, lvalue, requires-expression) added to the glossary with one line each |
| 15 | accept | The quick reference gains the two lines of code for the step (`root_epsilon<3, double>() * w`; `std::cbrt(eta) * std::max(std::abs(x), 1.0)`); 3.5 gains a margin note at the first `root_epsilon` |
| 16 | modified | 7.6: the two rules in the text in one sentence, the why in the sidebar. Sidebar 3.3 stays after 3.4 (it is about x + h, which 3.4 introduces) and 3.1 points to it |
| 17 | modified | For the author (a layout move): 5.4 keeps the two sentences and the cure; the pragma listing is quoted in 8.6 instead, which already points back; E8 stays compiled |
| 18 | accept | Editor: cleveref's appendix naming (preamble), so the references read "Appendix A" |
| 19 | modified | T6: the header phrases go; the caption paths stay (rejected part, reason in T6) |
| 20 | accept | T6 |
| 21 | accept | A margin gloss at the first O(h^2) in 2.2 |
| 22 | accept | 6.4 reworded as suggested |
| 23 | accept | Editor: glossary entries for cancellation, amplification, rounding floor, moment conditions, Richardson extrapolation, aliasing, unit roundoff, untrusted result, each with its section |
| 24 | modified | Figure 6.4: x tick labels in x (10^-8, 10^-10) beside k, the three regions labelled, one vertical line dropped; the two panels stay (the evaluations and the error are the two things the figure shows) |
| 25 | accept | No action now; both markers go at release (the version statement is filled, the UNVERIFIED notes are removed after this round's verification, section 4) |

### 2.2 R1-developer (18: 14 accept, 4 modified)

| Id | Decision | Where and what |
|---|---|---|
| 01 | accept | Sidebar 7.1, 8.7 and A.6 list every `numeric_limits` member the module uses (`is_specialized`, `is_integer`, `is_bounded`, `is_exact`, `epsilon()`, `max()`, `radix`, `digits`, `denorm_min()`); a compiled minimal user type E23 (`code/07-using-usertype.cpp`, not inheriting `numeric_limits<double>`), of which the member list is quoted; U-5 to the author |
| 02 | accept | T1 |
| 03 | accept | T11; U-8 |
| 04 | accept | 8.2: the check is an `assert`, compiled out under `NDEBUG`, after which the behaviour is unspecified (a negative step gives `Forward<1>` a backward difference, library-03); validate a computed step where it is computed, or use the overload without a step; U-6 |
| 05 | accept | 7.7: GCC evaluates some `<cmath>` functions in constant expressions as an extension; test a compile-time derivative on every toolchain shipped |
| 06 | accept | 7.1 and 8.2: a wider step converts to T; for `float` write `1e-3f`, since the narrowing is diagnosed under `-Wconversion` |
| 07 | accept | T11 |
| 08 | modified | E12 returns `std::unexpected` for non-convergence too, with a small solver enumeration that wraps `nxx::deriv::error`; the listing grows by a few lines; the sentinel goes |
| 09 | accept | 7.2: `error` has no `std::formatter`; print it with `static_cast<int>` or a `switch` with `default`; two lines on diagnosis (evaluate `fn` at x +- h yourself; compare h with the spacing of x); U-11 |
| 10 | accept | The margin note and 8.2: a member function of T itself; for a member of your model, capture the object in a lambda |
| 11 | modified | T4 |
| 12 | accept | The `check-by-steps` snippet writes `std::max(std::abs(x), 1.0)` and the text says that the default step is `root_epsilon<P + N, T>() * max(|x|, 1)` with no accessor; U-12 |
| 13 | accept | 10.5 gains a guarantee: no global or static state, `fn` called on the calling thread one call at a time, concurrent calls as safe as `fn`; verified by inspection |
| 14 | modified | T2 |
| 15 | accept | T1 |
| 16 | modified | E8 differentiates e^x at 1 instead of sin at 0.5, so that the text can quote this call's own values, which exist as the `\valFloatFourth...AtOne` macros of DS6; no new data set |
| 17 | accept | 7.2 and 8.3: "GCC accepts `if (*r)` without a warning; `int n = *r` is a `-Wfloat-conversion` warning on both" (with library-22) |
| 18 | accept | 7.7 and 8.3: "`==` with a value compares `value` only, through the conversion, and is false on failure" (with library-10) |

### 2.3 R1-numerics (17: 15 accept, 2 modified)

| Id | Decision | Where and what |
|---|---|---|
| 01 | accept | T8 |
| 02 | accept | T7 |
| 03 | accept | T7 |
| 04 | accept | T8 |
| 05 | accept | T12 (with evidence-11) |
| 06 | accept | T9 (with library-09) |
| 07 | modified | T10 |
| 08 | accept | 4.3 and Sidebar 3.1: dividing by the outer-node distance would quadruple the outer nodes' point rounding (1/12 to 4/12) and leave the inner nodes' untouched; only the two-node rules remove it |
| 09 | accept | T12 (with evidence-12) |
| 10 | accept | Sidebar 6.2: sum_k |c_k| / (d h^N) |
| 11 | accept | T15 |
| 12 | modified | One sentence in Sidebar 6.2 (not the running text): every uncovered case is e^(-x^2) at x above 1.6 (the Runge function at 0.54 for f'''), where the derivatives grow like (2x)^k and the first levels are not yet asymptotic; a smaller initial step restores the margin. The function and the smallest x are generated (`\valRiddersUncoveredFunction...`, DS9 extension) |
| 13 | accept | T3 |
| 14 | accept | Sidebar 6.1: 4^j exact, 4^j - 1 rounds once 2j exceeds the significand (j >= 13 in float, 27 in double), by at most half an epsilon |
| 15 | accept | 6.3 reworded (the factor 4^j relates the parent to the entry's left neighbour) |
| 16 | accept | T12 (with evidence-10) |
| 17 | accept | All eleven: (a) "four digits"; (b) 12! as the largest intermediate reached, 12! 4! a bound; (c) "13 halvings"; (d) say which decades are meant (the model's 3.8, the measured 3.5), T13; (e) print 99.96 or "one case 0.1 % over"; (f) Sterbenz with a, b of one sign; (g) T3; (h) T13; (i) "ten orders"; (j) T12; (k) T12 |

### 2.4 R1-library (22: 20 accept, 2 modified)

| Id | Decision | Where and what |
|---|---|---|
| 01 | accept | Sidebar 3.3, A.4, 8.2 and 7.2: the two-node rules fail only when both points round back to x; near a power of two a step can move x on one side only, and the central difference then divides by the half-ulp distance and returns a value that is mostly rounding; hence the precondition |
| 02 | accept | T1 |
| 03 | accept | 8.2 (with developer-04) and 10.5: "Forward<P> never below x and Backward<P> never above, for a step h > 0 (the precondition)"; U-6 |
| 04 | accept | The comment in E14 and 7.6: `cpp_bin_float_50` has no expression templates (`et_off`); `-> mp50` matters for `number<..., et_on>` and is a safe habit; `mp50(*r)` is needed for every Boost type |
| 05 | accept | T11 |
| 06 | accept | 8.5 and 8.8: the weight check is a `static_assert` in the body, not a constraint; a requires-expression reports such a call as valid; U-10 |
| 07 | accept | 8.2: a value is also returned, with its large `abs_error`, when the step reaches the resolution of x; an `abs_error` of a quarter of the value or more marks a result the method does not trust |
| 08 | accept | T6: 4.3, Table 10.5 ("stencils of each family") and 10.5 (`root_epsilon` for selected N from 2 to 2000 and the multiprecision types) |
| 09 | accept | T9 |
| 10 | accept | 8.3: comparing with a T compares the value only; `abs_error` saturates at `max()` where the estimate would overflow and then carries no information (with developer-18) |
| 11 | accept | "One call of each form" |
| 12 | accept | 5.1: "the two forms agree in exact arithmetic whenever the points are exact" |
| 13 | accept | Table 4.1's Ridders row as suggested, or a pointer to Table 8.2 |
| 14 | accept | Sidebar 6.2: the untrusted test at |A| (factor 1) instead of |A|/4 |
| 15 | modified | T14 |
| 16 | modified | Multi-line `if`s in quoted snippets are braced and style decision 18 is narrowed to multi-line `if`s (single-line ones build clean on both compilers and read better); `return 1;` **stays** in the examples (rejected part): the examples are tests, which 1.4 and the front matter say, and a `return` from `main` is the honest form of "stop" |
| 17 | accept | B.4 (with evidence-21): "these agree with the design record's 3.5, 11 and 15; the header's 1.64 for Central<10> is a rounding of kappa = 1.8, where the formula gives 1.63" |
| 18 | accept | Table 8.1 placed before 8.2's heading, or the lead-in ends without a colon |
| 19 | accept | Sidebar 7.1: rejected by `is_exact` (and `is_bounded`), with the reason |
| 20 | accept | Editor: the generator writes the untrusted-result row as "estimate at least |A|/4 and more than 64 times the floor"; 9.5 says which constants are read from the library (the cap, the divisors) and which are literals of its code, transcribed |
| 21 | accept | Table 7.1: `derivative<Central<2>, 2>`, `derivative<Ridders, 2>` |
| 22 | accept | With developer-17 |

### 2.5 R1-evidence (29: 23 accept, 5 modified, 1 for the author)

| Id | Decision | Where and what |
|---|---|---|
| 01 | accept | Editor: the AI statement rewritten (section 4); the author approves the text |
| 02 | accept | T7 |
| 03 | accept | 3.2: "falls by a factor of only 21, not the 100 that truncation alone would give: at h = 0.01 truncation is about 4.5e-5 of the 2.2e-4, the rest is already rounding"; the truncation term generated (`\valHandCalcTruncationHundredth`, HT extension) |
| 04 | accept | T7 |
| 05 | accept | T10; U-7 |
| 06 | accept | T3: `more2011` added and cited at the procedure; `more2012` stays at the step, 7.9 and 10.7 |
| 07 | accept | Sidebar 1.1: `squire1998` after the formula; `martins2003` added for the automatic-differentiation remark; "the complex step is a special case of Lyness and Moler's earlier contour-integral methods, which also give higher derivatives" |
| 08 | accept | Table 7.1: "model eps^(2/3) = 3.7e-11, measured median 5.8e-12 (e^x)" with the existing macros in rows 1, 5 and 6; row 1's evidence points to Figure 4.2 and Table 10.1 |
| 09 | accept | T12 |
| 10 | accept | T12 |
| 11 | accept | T12 |
| 12 | accept | T12 |
| 13 | accept | T12 |
| 14 | accept | 4.2: "up to about three times above the median and up to tens of times below it", generated as the DS7 extremes' ratios |
| 15 | modified | T6: attributed to the tests' comment and its three points, not generated |
| 16 | accept | Figure 7.2's caption: the plateau J/(h c_p) with `\valKinkCentralStep` and `\valKinkCentralFourStep`, which exist and were unused; the 1/h statement is made across the two rules, which the data support |
| 17 | accept | T15 |
| 18 | accept | The second option: the generator runs Ridders at 500 K with the NaN wrapper (E6 extension: status and evaluations) and 4.4 quotes the outcome; the cross-reference to Chapter 6 says "near such an edge (log x close to 0)" |
| 19 | modified | The JANAF table value at 1000 K enters the generator as a cited constant (`chase1998`, Part II, p. 1621) and the relative difference of the fit is generated (`\valCpFitVsTableRelDiff`); 1.4 quotes it with the citation. (Modified: the reviewer offers either; the generated form is chosen so that the sentence follows the data) |
| 20 | accept | 10.6: "each phase was reviewed by specialised review agents, Phase 4 also by another vendor's model (OpenAI Codex)" |
| 21 | accept | `docs/PHASE4_API.md` named at Sidebar 5.2, Sidebar 3.2, 5.3 and B.4; the half-step result attributed to the design record's prototype; B.4's sentence per library-17 |
| 22 | modified | T13 |
| 23 | accept | T13 |
| 24 | accept | Figure 3.5's caption: "each floor lies at or below eps^(2/3) of its type (dashed), with the rounding scatter up to a decade or more below it" |
| 25 | modified | Editor: `ieee754-2019` as `@manual` with organisation and number; the section numbers move from the notes to the citing sentences (`\cite[Sect. 5.4]{dennis1996}`); the UNVERIFIED notes are removed per the reference table; `articletitle=true` is a preamble change **for the author** (section 5) |
| 26 | for the author | U-9: the NIST SRD redistribution terms for the 24 coefficients in the MIT-licensed `shomate.hpp`; the attribution in the header comment stays whatever the answer |
| 27 | accept | `data/README.md`: "the reference values of the data sets avoid `pow`; the model figures of DS13 use it, at two or three digits" |
| 28 | modified | T14 |
| 29 | accept | "the Colebrook equation, in its common form" |

The reference table is accepted in full: `sterbenz1974` (address, series, ISBN, Theorem 4.3.1 p. 138 at the citing
sentences), `ieee754-2019` (entry type and fields), `colebrook1939` verified (the note goes), `more2011` and
`martins2003` added as proposed; `burden2016` cited with sections (1.1, 4.1, 4.2), `chapra2015` with "ch. 23",
`higham2002` with "Thm 2.5" where the lemma is used; `shomate1954` and `iso14882-2024` stay uncited and unprinted.

## 3. Revision plan by chapter file (stage 6)

Each author writes `responses/R1-<chapter-id>.md` listing, per accepted finding, what changed. The editor's items
(section 4) come first where a chapter needs a new macro, since the chapters quote them; the editor regenerates
`data/` and circulates the new macro names before the authors start.

| File | Owner | Items (by finding id and theme) | Pages |
|---|---|---|---|
| `01-why.tex` | fundamentals | 1.3 "throws nothing itself" (T11); 1.4 complete program with `else`, `book::` explained, the `*r->abs_error` gloss, the house idiom (T1, T4, engineer-09); Sidebar 1.1 citations and history (evidence-07); the JANAF comparison as a generated value with `chase1998` (evidence-19); 1.5 "Did it work / Is it right", the Ridders check, the caption glosses for w and eta (T1); headline digits by the rounded macros (T13) | +0.5 |
| `02-taylor.tex` | fundamentals | The O(h^2) gloss (engineer-21); drop "the example checks" (engineer-20) | 0 |
| `03-rounding.tex` | fundamentals | T5 (w in 3.3, s as the guess in 3.4; Sidebar 3.2 wording); 3.2 the calculator's ratio and "four digits" (evidence-03, numerics-17a); Sidebar 3.1 "a fixed multiple of h" and the outer-distance correction (T6, numerics-08); Sidebar 3.3 the resolution asymmetry, and a pointer from 3.1 (library-01, engineer-16); 3.3 the large-|f| sentence (engineer-05); 3.5 the units example, the `root_epsilon` margin note, the smooth-interpolant sentence, the misjudgement factor by macro (engineer-11, 15, 13; T9); 3.6 the residual-tolerance sentence, E21 and its numbers, "up to a known factor", `more2011`, cancellation inside f (T3); 3.7 Figure 3.5's caption (evidence-24); 3.8 the sixth question (engineer-05); the prototype's half-step result attributed to `docs/PHASE4_API.md` (evidence-21) | +0.9 |
| `04-stencils.tex` | methods | Eq. (4.3) sign (T8); the developer detail into Sidebar 4.1 and the test budget to 10.6 (T6); 4.3 the outer-distance correction (numerics-08); Sidebar 4.1 "12!" (numerics-17b); Figure 4.1's decades and 4.2's whiskers (numerics-17d, evidence-14); Table 4.1's Ridders row (library-13); 4.4 Ridders at 500 K from the new macros and the Chapter 6 pointer (evidence-18); 4.4 the cap as a design decision (T10); "as the header" phrases (engineer-19) | -0.2 |
| `05-orders.tex` | methods | The misjudgement factors by macro in Sidebar 5.1 and 5.3 (T9); the model digits by macro (T13); 5.4 keeps two sentences, the listing goes to 8.6 if the author agrees, E8 on e^x at 1 (engineer-17, developer-16); "as the header warns" (engineer-19); 5.1 "agree in exact arithmetic" (library-12); Sidebar 5.2 names `docs/PHASE4_API.md` (evidence-21); 5.6 drop "the example checks" (engineer-20) | -0.4 |
| `06-ridders.tex` | methods | 6.1 the paper against the library (T15); 6.3 the 4^j sentence (numerics-15); Figure 6.3's caption and 6.3: zero-error cases, the scope for f' (T7, T12); 6.4 lead with the conclusion, the 24-against-16 comparison, the edge rule's two readings, "ten orders", the ambiguous sentence (engineer-12, T12, numerics-17i, engineer-22); 6.5 the tail of the histogram by macros, cost 2.3 to 7, the gain against Central<4> (T7, T12); 6.6 the gain (T12); Figure 6.4's readability and Figure 6.5's label (engineer-24, T7); Sidebar 6.1 Burden and 4^j - 1 (T15, numerics-14); Sidebar 6.2 the floor formula, the constants' attribution, the uncovered-case class, the zero-error cases, "13 halvings", 99.96 % (numerics-10, library-14, numerics-12, T7, numerics-17c, 17e) | +0.2 |
| `07-using.tex` | practice | 7.1 getting the library, the `float` step (T4, developer-06); 7.2 "throws nothing itself", `.value()`, printing `error`, the member pointer, `if (*r)` on GCC, the Ridders tolerance check as a snippet, E11's reads inside the test, "monadic" in plain words (T11, developer-09, 10, 17, T1, engineer-14); Table 7.1 model against measured, the code forms (evidence-08, library-21); 7.4.1 E12 with `std::unexpected` for non-convergence (developer-08); 7.4.2 renamed, the loop printed, E22 (T2); 7.4.3 the conclusion first and the action (T3); 7.5 the house idiom, `max(|x|, 1)`, plain words (T1, developer-12, engineer-12); 7.6 the two rules in the text, `et_off` (engineer-16, library-04); 7.7 the GCC extension, `==` with a value (developer-05, 18); 7.8 the Ridders pitfall, the edge rule's two readings, Figure 7.2's caption, "at most" (T1, T12, evidence-16, evidence-22); 7.9 mixed partials and caching, the spline advice, E22 as the Jacobian stopgap (engineer-10, 13, T2); Eq. (7.1) "in its common form" (evidence-29) | +1.6 |
| `08-api.tex` | practice | 8.2 the `assert` under `NDEBUG`, both points at the resolution, the untrusted result at the resolution, the `float` step, the member pointer, the caps by macro (developer-04, library-01, 07, developer-06, 10, evidence-23); 8.3 `==` with a T, the saturated `abs_error`, `.value()`, `int n = *r` (library-10, T11, library-22); 8.5 and 8.8 the `static_assert` against the constraints (library-06); 8.6 the pragma listing from 5.4 (engineer-17, if approved); 8.7 the `numeric_limits` members (developer-01); "one call of each form", Table 8.1's placement (library-11, 18) | +0.5 |
| `09-tables.tex` | practice | 9.3 the cap as a design decision (T10); 9.5 which constants are transcribed (library-20); Table 9.10's row text comes from the generator | 0 |
| `10-guarantees.tex` | practice | 10.1 the model-to-measured ratios per family by macro (T12); 10.2 the scope for f' (T12); 10.4 the caps as design decisions with the model's gains, the edge rule's two readings (T10, T12); 10.5 "for a step h > 0", thread safety, the caps by macro, `root_epsilon` for selected N (library-03, developer-13, evidence-23, library-08); 10.6 the test budget from 4.3, Table 10.5 "stencils of each family", the phase reviews (T6, library-08, evidence-20) | +0.3 |
| `a1-floating-point.tex` | fundamentals | A.4 both points at the resolution (library-01); A.5 Sterbenz with a, b of one sign (numerics-17f); A.6 the `numeric_limits` members the module uses (developer-01) | 0 |
| `a2-derivations.tex` | methods | B.1 the single-point remainder by the interpolation error formula, the cases named (T8); B.4 the cost figures' attribution (library-17, evidence-21) | +0.1 |
| `qr-quick-reference.tex` | practice | The two checks and the Ridders tolerance (T1); `.value()` among the pitfalls (T11); the two lines of code for the step (engineer-15); the caps by macro (evidence-23); the edge rule's two readings (T12); the rounded headline digits (T13) | 0 |

The quick-reference page must stay one page; if the additions overflow it, the "Where to look" paragraph is cut
first.

## 4. The editor's items (before the authors start, except where noted)

1. **The AI statement** (evidence-01): rewritten so that it is true to the record: the coordinating main session
   (Claude Opus 5.5), which set up the folder, the layout adaptation, the build, `check.py` and the generator's
   skeleton and coordinates the agents; the editor (Fable 5.1); the authors (fundamentals and practice Opus 5.5,
   methods Fable 5.1); the stage-1 reviews of the plan by the engineer reader and the numerics reviewer; the stage-2
   verification of the thread's data and the generator's first results by the evidence reviewer; round 1 by the five
   reviewers; round 2 and the Codex review to come; and the sentence on what was generated: "the numbers, the plotted
   data and the bodies of the generated tables were computed by running the library through the book's data
   generator, written by the editor; no model typed them. The code listings were written by the authors and are
   compiled and run as tests; the declarations of Chapter 8 are transcribed from the header and checked against it"
   (T14). The sentence on the references becomes "every reference was checked by the evidence reviewer against the
   work or a catalogue" once the notes are gone. The author approves the text.
2. **`refs.bib`**: the two corrected entries, the two additions, the notes removed, the section-at-citation form
   (evidence-25 and the reference table).
3. **Preamble and macros** (for the author's approval, section 5): cleveref's appendix naming (engineer-18);
   `articletitle=true` (evidence-25). Glossary entries (engineer-23, engineer-14); the notation table's w (T5); the
   front matter's listing sentence (T14).
4. **The generator and data sets** (each named so that the authors can quote them):
   - DS9 extension (T7, numerics-12): `\valRiddersTailCasesSecond`, `\valRiddersTailRestartsSecond`,
     `\valRiddersTailFunctionsSecond` (the function names), `\valRiddersZeroError{First..Fourth}` (renamed from
     `RiddersExactCases`), `\valRiddersZeroOnQuinticFirst`, `\valRiddersZeroFunctionsFirst`,
     `\valRiddersUncoveredFunction{First,Second,Third}` and `\valRiddersUncoveredMinX{First,Second}`; the
     "restarts" label removed from `figures/06-ridders-cost.tex` (the methods author's file; the editor removes it
     with notice, as for F2 at stage 2); the coverage of f''' printed to two decimals (numerics-17e).
   - DS10 extension (T12): `\valRiddersEdgeEvalsInterior` (the k = 2 row).
   - E6 extension (evidence-18): `\valEdgeRiddersStatus`, `\valEdgeRiddersEvals` (Ridders at 500 K with the NaN
     wrapper).
   - HT extension (evidence-03): `\valHandCalcTruncationHundredth`.
   - DS13 extension (T9, T10, T13): `\valMisjudgeCost{CentralTwoFirst,CentralTwoFourth,CentralFourFourth}`,
     `\valCapGainCentralTwelveFirst{Double,Multiprecision}`, `\valModelDigitsCentralTwo{First..Fourth}Double`;
     the headline macros rounded instead of floored (their values change: `HeadlineFourthDigits` 4 to 5,
     `HeadlineRiddersFirstDigits` 13 to 14; every sentence quoting them is re-read).
   - DS15 extension (T12): the model-to-measured ratios per family and order, or their extremes
     (`\valModelRatioCentralMax...`, `\valModelRatioForwardFourThird`, `...Fourth`).
   - DS7 extension (evidence-14): the extremes' ratios to the median (`\valSpreadAboveMax`, `\valSpreadBelowMax`).
   - DS1b extension (evidence-19): `\valCpTableValue` (32.697, `chase1998`) and `\valCpFitVsTableRelDiff`.
   - DS14 (library-20): the untrusted-result row's text.
   - DS22 and E21 (T3), E22 (T2), E23 (developer-01): if the author approves; the generator recomputes E21's
     recovered eta and its cost and E22's errors.
5. **`tools/check.py`** (T14): the `signature`-block check against the headers. The extension to integers is deferred
   (T13).
6. **PLAN.md**: section 8 gains the two recorded exceptions (the transcribed declarations; the CMake file) and the
   narrowed decision 18; 12.3 gains U-5 to U-12 (done with this letter); 4.3 and 5 gain the new sets and examples once
   approved; `data/README.md` per evidence-27 and the new macros.

## 5. For the author

### (a) Contested or scope-changing points, with the editor's recommendation

1. **E21, the measured-eta recipe** (T3): a new compiled example and data set, about 0.7 page in Chapter 3. Recommend
   accept: it closes the one situation of Table 1.1 the engineer could not carry out, and the procedure has a
   citation (`more2011`).
2. **E22, a real Jacobian** (T2): a new example, about 0.6 page in Chapter 7, which prints a six-line hand-applied
   stencil that is not a library call. Recommend accept: it answers the question Chapter 1 raises, says plainly what
   the library supports, and keeps "shared evaluations" deferred.
3. **A complete first program and the "Getting the library" CMake snippet** (T4): the snippet is the one listing not
   compiled by the book's build; it would be quoted from `code/getting-started.cmake` and recorded as an exception.
   Recommend accept (about 0.5 page). Alternative: state the three CMake lines in prose.
4. **E23, a minimal user type** (developer-01): a compiled example that tests the documented requirements, with a
   short quoted member list (about 0.3 page). Recommend accept: it is the only way to keep the requirement list true.
5. **The Ridders tolerance check in 1.5 and on the quick reference** (T1): changes the book's headline message from
   "one check" to two. Recommend accept.
6. **Headline digits rounded, not floored** (T13): "about 4 digits" becomes "about 5" for a fourth derivative at
   1.2e-5, and Ridders' "about 13" becomes 14. Recommend accept (the rounded figure is what "about" means).
7. **The pragma listing moved from 5.4 to 8.6** (engineer-17). Recommend accept.
8. **Preamble changes**: cleveref's appendix naming so that references read "Appendix A" (engineer-18);
   `articletitle=true` so that article titles print (evidence-25); and the stage-4 "Appendix A" heading fix (S4-12)
   still awaiting approval. Recommend accept all three.
9. **Style decision 18 narrowed to multi-line `if`s** (library-16). Recommend accept.
10. **Mixed partials** (engineer-10): a paragraph with the model's rule of thumb and no generated value. Recommend
    accept as such.
11. **The data-file paths stay in the captions** (engineer-19, rejected part). The author may overrule; the
    alternative is a list of figures in the back matter (about one page).
12. **The AI statement** (evidence-01): the rewritten text of section 4 needs the author's approval, as every version
    does.

### (b) Library issues outside the book, added to PLAN.md 12.3

- **U-5** (developer-01): `nxx::real` accepts a type that the code then rejects: `deriv.hpp` and `math.hpp` use
  `numeric_limits<T>::radix`, `digits` and `denorm_min()`, which the concept does not require; the tests do not see it
  because their minimal type inherits `numeric_limits<double>`. Either the concept requires them or the header stops
  using them.
- **U-6** (developer-04, library-03): the `h > 0` precondition is an `assert`, compiled out under `NDEBUG`; a release
  build then returns a backward difference for `Forward<1>` with a negative step (the multi-node rules and Ridders
  return `non_finite`).
- **U-7** (numerics-07, evidence-05): the header's and `design.md`'s rationales for the caps, "under half a digit
  beyond P = 10 even for 50-digit types" and "beyond 6 the amplification would cost more than the order gains", are
  not supported by the model they cite (0.7 digit at 50 digits for Central<12>; Forward<7> still gains by the model).
- **U-8** (developer-03): `std::expected::value()` throws `std::bad_expected_access` on failure (or terminates under
  `-fno-exceptions`) and converts silently to T; the header's documentation and `design.md`'s "nothing throws" could
  name it.
- **U-9** (evidence-26): the 24 Shomate coefficients in the MIT-licensed `code/shomate.hpp` are NIST Standard
  Reference Data ("(c) U.S. Secretary of Commerce ... All rights reserved"); the author confirms NIST's redistribution
  terms, and the attribution in the header comment stays.
- **U-10** (library-06): the weight-exactness check is a `static_assert` in `stencil_derivative`'s body, not a
  constraint, so a requires-expression reports a call with `std::float16_t` and `Central<10>` as valid and the call
  then fails inside the library.
- **U-11** (developer-09): `nxx::deriv::error` has no `std::formatter`, and `non_finite` covers six causes that a
  caller cannot tell apart.
- **U-12** (developer-12): there is no public accessor for the default step; a reader who repeats a stencil at 2h and
  h/2 must rebuild `root_epsilon<P + N, T>() * max(|x|, 1)` by hand.
- Already recorded: U-3 (the header's "the operator templates", D-14) and U-4 (the 1.64 rounding, D-16).

## 6. The page budget

The build is at 83 pages, inside the band of 74 to 88. The additions of section 3 sum to about +4.3 pages, the
removals and moves to about -0.6 (the developer detail into sidebars is a move, net 0; the pragma listing moves from
Chapter 5 to Chapter 8, net 0): expected total **about 87 pages**, inside the band with a page to spare. If the stage-6
build exceeds 88, the cuts follow PLAN.md 1.1's order: Appendix B's worked Lagrange example (B.2) first, then the
weight tables of Chapter 9 to `\footnotesize`, then Sidebars 7.1 and 7.2 shortened to half a page each, then one
sidebar of a methods chapter moved to Appendix B. If the author declines E21 or E22, the total is about 86.

# PLAN: Numerical Derivatives with Numerixx

Stage 1 output (editor, 2026-10-23). Awaiting the author's approval; once approved it wins over `BRIEF.md`. Rules:
`CLAUDE.md` in this folder. Everything here was checked against the ground truth in CLAUDE.md's order of precedence:
`include/numerixx/deriv.hpp` and `math.hpp`, then `tests/`, then `docs/design.md` and the Phase documents, then the
tutorial, the techniques survey and the README. Section 11 lists where those sources disagree.

**State of the library.** `BRIEF.md` describes commit `ae17769`. Verified: `include/`, `tests/` and the design
documents are byte-identical between `ae17769` and `HEAD` (`1ffb38f`); the two later commits only added this folder.
The brief's summary of the library is correct, with the two precisions noted in section 11 (D-5).

**Revision 2** (2026-10-23, after the stage-1 reviews by the engineer reader and the numerics reviewer): section 12
records the decision on every finding; the changes are folded into sections 1 to 11.

**Revision 3** (2026-10-23, stage 2, editor): the thread's data are verified (A1 done), and the evidence changed
four things, each recorded in section 12.2 and folded into sections 1.3, 4.3, 5 and 6: the Shomate coefficients are
the NIST Chemistry WebBook's own fit and Chase 1998 is cited for the tabulated data only; the junction of the two
Shomate ranges at 500 K is a range break with a jump in $H$, not a kink, and is used as such (E19, F18, DS20); the
default-step accuracies of DS6, DS7 and DS15 are medians over eleven points, not one point's luck; the inner solver
of E20 is bisection, whose output is quantised at its tolerance. The thread's function lives in one shared header,
`code/shomate.hpp`, quoted with `\snippetheader`.

**Revision 4** (2026-10-23, stage 4, editor): the integration. Section 12.4 records what the authors' stage-3 reports
and the full build changed: the discrepancies D-14 to D-17 in section 11, the corrections to S5, S13, 1.4 (Chapters
3 and 4) and 4.3 folded in, the shared T0 rows, the generated hand tables, the macros added, and the layout test,
gallery and template data set removed.

How to read this plan: section 1 is what each author writes; sections 2 to 6 are the shared inventories (notation,
sidebars, figures with their data sets, examples, references); section 7 says who owns what; section 8 fixes style
beyond CLAUDE.md; section 9 is the schedule; section 10 resolves the brief's open questions and lists what the author
must decide; section 11 records the discrepancies between the sources.

---

## 1. Chapter outline and page targets

### 1.1 Budget

Layout: 114 mm text block, 10 pt, about 350 to 400 words per page of running text; a text-width figure with its
margin caption about 0.4 page, a margin figure 0.2, a diagram 0.3, a sidebar 0.6 to 0.8, a listing of 6 to 12 lines
0.3, a full-width table 0.5. The word budget per chapter below is the running text left after those items, at 375
words per page. Nominal total 82 pages; acceptable band at integration 74 to 88; hard limits 60 and 90 (CLAUDE.md).
If the integrated build exceeds 88, the editor cuts in this order: Appendix B derivations, the weight tables of
Chapter 9 to `\footnotesize`, S13 and S14 shortened to half a page, then one sidebar per methods chapter moved to
Appendix B.

| Block | Item | Pages | Actual at stage 6 |
|---|---|---|---|
| Front matter | Title page | 1 | 1 |
| | About this document (preface, how to read, AI statement) | 2 | 2 |
| | Contents (tocdepth 1) | 2 | 3 |
| | List of sidebars | 1 | 1 |
| Part pages | Part I, Part II | 2 | 2 |
| Part I | 1 Why numerical derivatives | 5 | 6 |
| | 2 Finite differences from Taylor's theorem | 6 | 5 |
| | 3 Rounding and the choice of step | 9 | 9 |
| | 4 Higher-order and one-sided stencils | 7 | 6 |
| | 5 Second to fourth derivatives | 6 | 5 |
| | 6 Ridders' extrapolation | 7 | 8 |
| | 7 Using the library | 9 | 13 |
| Part II | 8 API reference | 5 | 7 |
| | 9 Tables | 5 | 4 |
| | 10 Accuracy, cost, limits and guarantees | 5 | 5 |
| Appendices | A Floating-point essentials | 2 | 2 |
| | B Derivations | 4 | 4 |
| Back matter | Notation and glossary | 1 | 3 |
| | References | 2 | 1 |
| | Quick reference (the last page) | 1 | 1 |
| **Sum** | | **82** | **88** |

The last column is the stage-6 close-out build (2026-10-24): title page, vi, 81 numbered pages, at the top of the
acceptable band. Chapter 7 carries the round-1 additions (E1b's checks, E21 to E23, the CMake file, the decision
table's evidence) and is the chapter to cut first if round 2 adds length; its listings are kept whole by
`\needspace` guards, which cost a few blank lines per page break.

At the stage-7 close (2026-10-24, after the round-2 fixes) the counts are the same, block for block: 88 pages
(title page, vi, 81), with Chapter 7 at 13 and Chapter 8 at 7. The author has said that a slightly longer book is
acceptable; the close-out trims (12.6, S7-6) removed three pages that held two or three lines each, not length.

Front matter 6, part pages 2, Part I 49, Part II 15, appendices 6, back matter 4: 82.

### 1.2 Chapter files

| Chapter | File | Author |
|---|---|---|
| Front matter | `chapters/00-frontmatter.tex` | editor |
| 1 | `chapters/01-why.tex` | fundamentals |
| 2 | `chapters/02-taylor.tex` | fundamentals |
| 3 | `chapters/03-rounding.tex` | fundamentals |
| 4 | `chapters/04-stencils.tex` | methods |
| 5 | `chapters/05-orders.tex` | methods |
| 6 | `chapters/06-ridders.tex` | methods |
| 7 | `chapters/07-using.tex` | practice |
| 8 | `chapters/08-api.tex` | practice |
| 9 | `chapters/09-tables.tex` | practice |
| 10 | `chapters/10-guarantees.tex` | practice |
| A | `chapters/a1-floating-point.tex` | fundamentals |
| B | `chapters/a2-derivations.tex` | methods |
| Back matter | `chapters/99-backmatter.tex` (notation, glossary) | editor |
| Quick reference | `chapters/qr-quick-reference.tex`, input after the references as the book's last page | practice |

Figures are `figures/<chapter-id>-<name>.tex` and examples `code/<chapter-id>-<name>.cpp`, where `<chapter-id>` is
the chapter file's prefix (`03-rounding`). Labels: `ch:`, `sec:`, `fig:`, `tab:`, `sb:`, `eq:`, `lst:` with the ids in
sections 3 to 5.

### 1.3 The thread example

One engineering example runs through Part I: **the heat capacity of a gas from its enthalpy**, $c_p = \dd H / \dd T$,
where $H(T)$ is a molar enthalpy in the Shomate form (a quartic in $T/1000$ plus a $1/T$ term, as in the NIST-JANAF
tables), returned in J/mol for $T$ in K, so that $\dd H/\dd T$ is $c_p$ in J/(mol K) directly. Why this one:
- a derivative engineers compute, with a closed form for $c_p$ and for every higher derivative, so the error of every
  method can be measured honestly (the evidence rule);
- $|x| = T$ is in the hundreds or thousands of kelvin, so the scaling of the default step with $\max(|x|, 1)$ is
  visible, and steps carry a unit (a step of a few millikelvin);
- the fit is valid on a stated range, so its lower end is a realistic edge for the one-sided rules ("evaluable is not
  valid"), and $\dd c_p / \dd T = H''$ is a meaningful second derivative;
- Newton's method for the temperature with a given enthalpy (an energy balance) is a natural use of a derivative
  inside an iteration (Chapter 7).

The choice (verified at stage 2, decision A1): nitrogen gas, the Shomate range 500 to 2000 K of the NIST Chemistry
WebBook, evaluated at $T_0 = \SI{1000}{K}$, which is interior to the range; its lower end, 500 K, is the domain edge
used in Chapters 4 and 7. The coefficients are the WebBook's own fit ("new parameter fit January 2009") to the
NIST-JANAF tables of Chase 1998, which themselves contain no coefficients: the text cites the WebBook for the
equations and the coefficients and Chase 1998 for the data behind them (section 6). $H(T)$ is the enthalpy relative
to 298.15 K, in J/mol; the book says so once. The fit's $c_p$ at 1000 K is 32.69173 J/(mol K) against the JANAF
table's 32.697, so the closed form is the exact answer key for the *function* but describes nitrogen to about four
figures only: the text never presents the closed-form digits as nitrogen's heat capacity. The generator and the
examples share one definition, the header `code/shomate.hpp` (`shomate_H`, `shomate_cp`, `shomate_derivative` up to
order 6, the three ranges, the junction), quoted with `\snippetheader`; the generator checks it against the evidence
reviewer's exact evaluation at every run. The book writes $H(T)$ and $c_p$; never $h$ for the enthalpy, which is the
step. Chapter 1 says plainly that the closed form is the answer key: the book differentiates a function whose
derivative it knows so that every error can be measured, and Chapter 7 adds a function with no formula at hand (E20).
The junction of the two lower ranges at 500 K is the range break of the smoothness pitfall (E19, F18, DS20): the two
published fits do not meet, $H$ jumps by 0.77 J/mol there and $\dd c_p / \dd T$ by 17 %, so the piecewise $H$ is
discontinuous rather than kinked (decision S2-2 in section 12.2); the lesson "evaluable is not smooth" covers both.

The "method" figures (order ladders, Ridders' coverage, float against double) use $\mathrm{e}^x$ at $x = 1$ and the
test functions, because every derivative of $\mathrm{e}^x$ equals the function, so errors are comparable across
orders; the chapter text says so when it switches. The thread function appears in every Part I chapter at least once,
as a listing, a figure or a quoted number.

### 1.4 Chapters

Each entry: learning goal, sections, page target with the word budget, and the figures (F), sidebars (S), examples (E)
and tables (T) it uses, by the ids of sections 3 to 5.

#### Front matter (editor, 6 pages)

Title page (`main.tex`). "About this document": the purpose (a textbook that teaches the methods and a reference for
the module), the audience, how to read (Part I in order, Part II when needed, sidebars optional, the quick-reference
page at the very end for the call, 1.5 for the method choice), which version of the library it describes (the commit on `dev` at release), and "How this document was
prepared" (the AI statement, section 9, stage 4). Contents, list of sidebars.

#### Chapter 1: Why numerical derivatives (5 pages, about 1250 words)

Learning goal: the reader knows where derivatives arise in engineering practice when no formula is at hand, when a
numerical derivative is the right tool and when it is not, has seen the library compute one, and leaves with a
method choice and the one check every call needs.

Sections: 1.1 Derivatives without a formula (property models, simulators, implicit functions, sensitivities, Newton's
method). 1.2 The alternatives in one page (analytic, automatic differentiation, complex step; when each fits; the
numerical derivative's place: any black-box `fn(T)`). 1.3 What the library offers (one function, four methods, orders
1 to 4, any real type, no exceptions, an error estimate from Ridders). 1.4 A first example (the thread: $c_p$ at
$T_0$ with the default step, then with Ridders and its estimate; the closed form is the answer key, stated as such).
1.5 If you read nothing else (half a page: T0, the five-row method choice; the headline accuracy from T14, about ten
digits for $f'$ at the default step in double and about four for $f''''$, as `\val` macros; the headline cost from
T15; the minimal `if (!d)` check; where each of these is explained).

Uses: F1 (margin diagram), S1, E1 (two snippets), T0 (the five-row method choice, text width), the generated values
`\valCpExact`, `\valCpDefaultRelError`, `\valCpRiddersRelError`, `\valCpRiddersEstimate`, `\valCpRiddersEvals`,
`\valHeadlineFirstDigits`, `\valHeadlineFourthDigits`.

#### Chapter 2: Finite differences from Taylor's theorem (6 pages, about 1650 words)

Learning goal: the reader understands where the forward and the central difference and their error terms come from
(Taylor's theorem), knows what the truncation error and the accuracy order are, and can predict how the error changes
when the step is halved while truncation dominates.

Sections: 2.1 The secant and the tangent (the forward difference, its error $h f''/2$; the library's name for it,
`Forward<1>`, given at once with the margin note "P is the accuracy order, not the number of points"). 2.2 Centring
cancels a term (the central difference, `Central<2>`, $h^2 f'''/6$, the even expansion; the sign of the leading term
stated once: the central difference overestimates $f'$ where $f''' > 0$). 2.3 Order of accuracy: what $\order(h^2)$
means (asymptotic: a leading term $C h^P f^{(P+N)}$, the next $h^{P+2}$ for central rules; halving the step quarters
the error while truncation dominates; measured on the thread, F2 and E2). 2.4 Exactness on polynomials (a quadratic
for the central difference) and a first look at the next term. 2.5 What this chapter did not say (rounding: pointer
to Chapter 3; more points: Chapter 4).

Uses: F2, S2, E2, T-margin (the two rules side by side: tag, nodes, weights, error term; the only weight table in the
chapter; complete tables are in Chapter 9, `\cref{tab:weights-central-first}`). Equations numbered: the two
differences and their error terms (referred to from Chapters 3, 4, 6 and Appendix B).

#### Chapter 3: Rounding and the choice of step (9 pages, about 1750 words)

Learning goal: the reader understands why a smaller step is not always better, can read the V-shaped error curve, knows
the default step, what accuracy to expect from it, when to pass a step of their own, and what changes when the
function itself is less accurate than machine epsilon.

The chapter stays with the central difference: its model is written concretely as $E(h) = h^2 |f'''| / 6 +
\epsmach |f| / h$; the general constants $\ampl$ and $\errc$ are introduced in Chapter 4 and the order N in Chapter 5.

Sections: 3.1 Floating-point numbers in brief (a number has about 16 digits; machine epsilon; rounding of $f$; link to
Appendix A). 3.2 The two errors and the V (F3; the concrete model; the model curve over the measured one; the model
gives an absolute error, and "relative accuracy $\epsmach^{2/3}$" presupposes $|f| \sim s |f'|$, which the thread
satisfies to a factor stated as a `\val`). 3.3 The best step and the best accuracy ($h^\ast \sim \epsmach^{1/3}$,
error $\sim \epsmach^{2/3}$: about 10 to 11 of 16 digits). 3.4 The default step ($\epsmach^{1/3} \scale$ for this
rule, with $\scale = \max(|x|, 1)$: relative for large $|x|$, absolute near 0; why the constant is dropped: S5; the
resolution of $x$ and `non_finite`). 3.5 When the default step is wrong: the scale of $f$ (two cases. Large $|x|$
with a narrow feature, F4 and E4: the fix is to replace $s$ by the scale $w$ on which $f$ varies, $h = \epsmach^{1/3}
w$, and the figure shows that step landing at the V's minimum. Small $|x|$ in the user's units, F4b and E4: a mole
fraction of $10^{-6}$ or a parameter of $10^{-4}$ has $s = 1$, so the default step is as large as $x$ or larger and
can cross a physical bound such as zero; the same fix. Or use Ridders with an initial step on $f$'s scale, Chapter 6).
3.6 When $f$ is less accurate than $\epsmach$ (a solver output converged to $10^{-8}$, a simulation that iterates to
convergence; a smooth fit is inaccurate but not noisy, its derivative is the fit's, decision S4-6:
the noise level $\noise$ replaces $\epsmach$ in the rounding branch, the best step becomes about $\noise^{1/3} s$ and
the accuracy about $\noise^{2/3}$, so a $10^{-8}$ solver gives about five digits; Ridders' estimate assumes
$\epsmach$ and becomes optimistic; F17 and E17; pass the step). 3.7 Float, double and more digits (one paragraph with
F5 in the margin). 3.8 A checklist for the step (the scale of $f$ in your units, bounds within a step, the accuracy of
$f$, the order, when to let Ridders choose).

Uses: F3, F4, F4b (margin), F5 (margin), F17 (margin), S3, S4, S5, E3, E4 (two snippets), E17. Values:
`\valCpDefaultStep`, `\valCpDefaultRelError`, `\valCpBestStep`, `\valCpBestRelError`, `\valCpRelativeFactor`,
`\valScaleDefaultRelError`, `\valScaleFitRelError`, `\valSmallDefaultStep`, `\valSmallFitRelError`, `\valNoise...`,
`\valEpsDouble`, `\valEpsFloat`, `\valDigitsLostCentralTwo`.

#### Chapter 4: Higher-order and one-sided stencils (7 pages, about 1350 words)

Learning goal: the reader knows that more points buy a higher order, what that costs in rounding and smoothness, when
to use `Central<P>` with P above 2, and how to differentiate at the edge of a domain with `Forward<P>` and
`Backward<P>`.

Sections: 4.1 Stencils: nodes, weights and the moment conditions (the four-point central rule derived in half a page;
the rest in S6 and Appendix B). 4.2 The order ladder (F6: `Central<2>` to `Central<10>` on one plot; each pair of
points gains two orders; the rounding branches nearly coincide, and a higher order buys accuracy by letting the step
grow, so each minimum slides up the common rounding line; diminishing returns for double, F8). 4.3 The price (the
general error model is introduced here: $E(h) = \errc h^{\aorder} |f^{(\aorder+1)}| + \ampl \epsmach |f| / h$ with
$\errc$ and $\ampl$ defined from the weights; the third term that every rule with more than two nodes carries, the
rounding of the points, $(\ampl/2) \epsmach |x f'| / h$, which is why the tests budget $(\errc + 2\ampl)$ and why the
two-node rules escape it, S4; smoothness over a wider stencil: the error term needs $f$ smooth across the whole
stencil; the step that grows with the order). 4.4 One-sided rules (F7, F9 diagram): every rule gains one order per
node, and a central rule reaches the same order with one node fewer because antisymmetry supplies its even moments
(decision S4-5); their larger constants, the edge of a domain (E6 on the thread's validity range: a central rule samples outside
the fit; `Forward<2>` does not; "evaluable is not valid"). 4.5 Choosing P (a short rule: `Central<2>` or `<4>` for
double inside an iteration, higher only for wide types; one sentence on the caps 10 and 6 with a pointer to Chapter
10.4; pointer to the decision table in Chapter 7 and the full tables in Chapter 9).

Uses: F6, F7 (margin), F8 (margin), F9 (diagram), S6, S7, E5, E6, T-margin (nodes and evaluations per family, as the
header documents them: `Central<P>` has $\aorder + \dorder - 1$ nodes, one-sided rules $\aorder + \dorder$).

#### Chapter 5: Second to fourth derivatives (6 pages, about 1250 words)

Learning goal: the reader can compute $f''$ to $f''''$, knows why the accuracy falls with the order and by how much,
knows the float warning and what to do about it, and understands the `order` argument.

Sections: 5.1 The rules for N = 2 to 4 (the three-point $f''$ from the moment conditions with $M_2 = 2 d$; Chapra's
tables as the reference; `Central<P>` keeps the meaning "error $\order(h^P)$" for every N, so the node count is
$\aorder + \dorder - 1$). 5.2 Why accuracy falls with the order (rounding grows as $1/h^N$; the model gives
$\epsmach^{\aorder/(\aorder+\dorder)}$: F10 the V-curves for N = 1 to 4, F11 the accuracy against N; S8; the
relative figure presupposes $|f| \sim s^N |f^{(N)}|$, and the thread's factor $|H| / (s^2 |H''|)$ is quoted as a
`\val`, so E7's model error is about $2 \cdot 10^{-7}$ relative, not $\epsmach^{1/2}$ itself; at a stationary point
relative accuracy is lost while absolute accuracy stays). 5.3 The default step per order and its sensitivity
(misjudging the scale by 2 costs up to 16 for $f''''$). 5.4 The float warning (what triggers it: any type with
epsilon above $10^{-7}$; how to silence it per call or with `NUMERIXX_NO_PRECISION_WARNING`; E8; the values
`\valFloatFourthCentralRelError`, `\valFloatFourthRiddersRelError`). 5.5 Ridders for higher orders (one paragraph:
same tableau, larger initial step, pointer to Chapter 6). 5.6 The thread: $\dd c_p/\dd T$ (E7). The three-point
$f''$ from the actual offsets is S9 only, referred to from 5.1 in parentheses.

Uses: F10, F11, S8, S9, E7, E8. Tables: none (Chapter 9).

#### Chapter 6: Ridders' extrapolation (7 pages, about 1050 words)

Learning goal: the reader understands Richardson extrapolation and Ridders' adaptive form in one picture, can read the
error estimate and knows when not to trust it, knows what the method does at poles and domain edges, and what it costs.

The text carries what a user needs; the stopping rules (the two-level stall, exact agreement, the rounding floor's
definition), the untrusted-result factors and the level caps are in S11 and Chapter 10, not in the running text.

Sections: 6.1 Richardson's idea (two steps, one cancellation; `Central<4>` as the one-level case, back-reference to
S7). 6.2 The tableau in one picture (F12 the generated tableau of the thread function, F13 convergence by column; the
difference form in S10). 6.3 Reading `abs_error` (the estimate is the change of the best entry from its parent, never
below the rounding of the values it came from; the method stops when more levels only add rounding; the result's
estimate also covers its distance to the latest diagonal entry; S11 gives the rules and the coverage statistics). 6.4
When not to trust it (a non-finite value: the edge, restarts, F15 and E10; a pole between samples or aliasing, S12;
the estimate is not a bound and assumes $f$ accurate to about $\epsmach$: a noisier $f$ makes it optimistic,
back-reference to 3.6; a kink inside the span, forward reference to 7.8). 6.5 Cost (2 or 4 evaluations per level;
typical and maximal counts from F14 and `\valRidders...`; the caps in Chapter 10). 6.6 Stencil or Ridders (why only
Ridders estimates its error; inside an iteration a stencil, for a reported value Ridders).

Uses: F12 (generated table), F13, F14 (scatter plus margin histogram), F15 (margin), S10, S11, S12, E9, E10. Values:
`\valRiddersCoverageFirst`, `\valRiddersCoverageSecond`, `\valRiddersMedianRatioFirst`, `\valRiddersWorstRatioFirst`,
`\valRiddersEvalsMedianFirst`, `\valRiddersEvalsMaxFirst`, the same for `Second`, `\valRiddersEdgeLimit`,
`\valCpRidders...`.

#### Chapter 7: Using the library (9 pages, about 1300 words)

Learning goal: the reader can call every form of `derivative`, handle its result, choose a method for a given job,
check a result that comes without an estimate, use the library with their own number types and at compile time, and
avoid the known pitfalls.

Sections: 7.1 The call (the two overloads, M and N as template arguments, the defaults; E16 is quoted in Chapter 8,
here E1 is recalled). 7.2 The result: `std::expected<Estimate<T>, error>` (E11: `if (d)`, `d->value`, `*d`, the
conversion and its traps, `abs_error` as an optional, `value_or`, `transform`/`and_then`; `error::non_finite` and the
cases that produce it; a `switch` needs a `default`). 7.3 Choosing a method (T1 the decision table, full width: rows
are situations, including the small-$|x|$ and the noisy-$f$ rows; columns the method, the step, the cost, the
accuracy to expect and where the evidence is). 7.4 Worked examples (E12 Newton's method for $T$ from $H$, F16; E13 a
sensitivity of a rate to temperature and to a parameter by capturing the other arguments in a lambda, the scaled
sensitivity $\dd \ln k / \dd \ln E_a$ that engineers report, and a loop that builds a Jacobian column by column with
its cost stated, gradients staying deferred; E20 a function with no formula at hand: a friction factor from an
implicit equation solved by an inner iteration, differentiated with respect to the Reynolds number, the solver
tolerance as the noise level of 3.6 and the implicit-function derivative as the answer key). 7.5 Is my derivative good
enough? (E18: repeat with $2h$ and $h/2$ or with `Central<4>` and count the agreeing digits; Ridders as an independent
check; absolute rather than relative error where $f'$ is near 0; a physical sanity check such as $c_p > 0$). 7.6 User
and multiprecision types (E14; `nxx::real`; S13 expression templates; `T(*r)` and `-> T`). 7.7 At compile time (E15;
S14 overflow guards and constant evaluation; the one GCC exception). 7.8 Pitfalls (the step scale at large and small
$|x|$; a kink or range break inside the stencil's span, E19 and F18: a central rule straddling it returns a wrong value
without failing, higher P buys nothing, Ridders' base differences then converge like $\order(h)$ so its estimate stays
large, "evaluable is not smooth"; aliasing at large $|x|$; noisy functions and the floor; domain edges and the minimum
distance; the fourth derivative in float; `auto` on the result; passing one step for several orders). 7.9 What the
library does not do (bare tables are out of scope: fit or interpolate first, and a linear interpolant has a
piecewise-constant derivative; deferred features named as such: gradients, several orders from shared evaluations, a
tolerance for Ridders, user stencils, noise-aware steps, complex step and automatic differentiation).

Uses: F16 (margin), F18 (margin), T1 (full width), S13, S14, E11, E12, E13, E14, E15, E18, E19, E20. Values:
`\valNewtonIterations`, `\valSensitivity...`, `\valJacobianEvals`, `\valSolver...`, `\valCheck...`, `\valKink...`,
`\valMpCentralSixRelError`, `\valMpRiddersRelError`.

#### Chapter 8: API reference (5 pages, about 1100 words)

Learning goal: every public name, with its parameters, preconditions, return value and failure cases, in one place.

Sections: 8.1 Synopsis (namespace `nxx::deriv`; what is public, what `detail` is; stability: `error` may gain values).
8.2 `derivative` (the two overloads and their deprecated twins; template parameters M, N, T, F and their constraints
`method`, `valid_order`, `callable` in words; parameters; precondition `h > 0` finite as an `assert`, "well above
ulp(x)" for stencils; return; the list of `non_finite` cases from the `error` documentation; exceptions pass through;
constant evaluation and its exception; fn is called only at finite points, once per node). 8.3 `Estimate<T>` (members,
the conversion and where it does not apply, `==`). 8.4 `error`, `order`. 8.5 The method tags (`Central<P>`,
`Forward<P>`, `Backward<P>`, `Ridders`: nodes and evaluations per N, the caps, the default steps, as a table). 8.6
`NUMERIXX_NO_PRECISION_WARNING`. 8.7 `nxx::real` and `nxx::math` (`abs`, `isfinite`, `root_epsilon`; the concept's
requirements in words). 8.8 Rejected calls (what does not compile and why: E16 second snippet).

Uses: E16 (two snippets), T2 (names and signatures, text-width), T3 (methods: nodes, evaluations, default step, cap).
Signatures are given as `\code{}` text and prose, never as a typed listing (CLAUDE.md); the usage snippets are
compiled.

#### Chapter 9: Tables (5 pages, about 350 words)

Learning goal: every number a user may want to look up, generated by the library and checked against the moment
conditions in the generator.

Sections: 9.1 How to read the tables (nodes in units of $h$, integer weights over one denominator; the error constant
$\errc = |M_{\dorder+\aorder}| / ((\dorder+\aorder)!\, d)$; the amplification $\ampl = \sum |w_k| / d$; `Backward<P>`
is `Forward<P>` mirrored with the weights times $(-1)^N$, so it is not tabulated). 9.2 Central rules, N = 1 to 4 (T4
to T7, full width). 9.3 Forward rules, N = 1 to 4 (T8 to T11, full width). 9.4 Default steps and model accuracy (T12:
exponent $1/(\aorder+\dorder)$, the step for double at $|x| \le 1$, the model error $\epsmach^{\aorder/(\aorder+\dorder)}$,
evaluations). 9.5 Ridders' constants (T13: ratio 2, initial step divisors 32 and 4, stall 2, untrusted-result factors
4 and 64, level cap 32, evaluations per level 2 or 4 plus fn(x) for even N, maxima 64, 65, 128, 129).

Uses: T4 to T13, all generated into `data/table-*.tex` (section 4.3).

#### Chapter 10: Accuracy, cost, limits and guarantees (5 pages, about 1400 words)

Learning goal: what accuracy to expect per method, order and type; what it costs; where the methods stop working; what
the library guarantees and under which hypotheses; how that is verified.

Sections: 10.1 Accuracy (T14 the measured matrix: method by order by type, $\mathrm{e}^x$ at 1; the model behind it;
the model gives absolute errors, and the relative figures presuppose $|f| \sim s^N |f^{(N)}|$; what "about" means).
10.2 Ridders' statistics (coverage, median and worst ratio, evaluations, from the random study of Chapter 6, F14, with
its scope stated: nine smooth functions on $[0.2, 3]$, no noise, no singularity within the initial step; the stopping
rules and the untrusted-result factors from S11; the level caps 32 levels, at most 64, 65, 128 and 129 evaluations).
10.3 Cost (T15: evaluations per method and order; typical and maximal for Ridders). 10.4 Limits (the order caps 10
and 6 and why: half a digit for 50-digit types, $\ampl = 27.7$; the domain-edge distance measured in F15; aliasing;
noisy $f$, 3.6; kinks, 7.8; the scale assumption). 10.5 Guarantees (fn called only at finite points; no exceptions
thrown or caught; no allocation; failures are values, in constant expressions too, except an overflow of the final
quotient under GCC; integer weights exact in every type with 24 significant bits or more, else a `static_assert`;
compile-time and run-time results agree under the hypothesis that the compiler neither contracts `a * b + c` into a
fused multiply-add nor uses extended precision, which the build does not pin, verified on the five CI toolchains; the
wording follows decision A5). 10.6 How it is verified (the test files and what each checks; references from closed
forms and published tables; five toolchains in CI; the independent reviews recorded in the Phase documents). 10.7
Deferred features (named, with a pointer to `docs/plan.md`).

Uses: T14 (generated), T15 (generated), back-references to F14 and F15.

#### Appendix A: Floating-point essentials (2 pages, about 650 words)

Learning goal: enough floating-point arithmetic to follow Chapter 3 and the sidebars: representation, machine epsilon
against unit roundoff, ulp, rounding to nearest, the largest finite number, subnormals, why $x + h$ rounds, catastrophic
cancellation. T16: the constants of float and double (generated: epsilon, digits, max, denorm_min). Author:
fundamentals.

#### Appendix B: Derivations (4 pages, about 1400 words)

Learning goal: the proofs the sidebars only cite. B.1 Taylor expansions and the moment conditions (the general
weighted sum, $M_j$, the error term with its constant). B.2 The Lagrange formula for the weights (why it solves the
moment conditions; the worked `Central<4>` example from the header's documentation). B.3 Richardson extrapolation (the
even expansion, the column formula with $4^j - 1$, the difference form, the a posteriori estimate). B.4 The error
model and the best step (the minimisation, the constant $(\dorder \ampl / (\aorder \errc))^{1/(\aorder+\dorder)}$ and
the cost of dropping it). Author: methods.

#### Back matter (editor 3 pages, practice 1 page)

Notation (section 2 as a table) and a glossary of the terms of section 8, item 4; the references (fullbib, two
columns); then the quick reference as the book's last page (`chapters/qr-quick-reference.tex`, practice author,
listed in the contents and pointed to from "How to read"): the include, the call in its four forms (inline `\code`),
the result type, the method tags with cost and accuracy in one line each, the five-row method choice (the same table
as T0 in 1.5, repeated by design: a card repeats), the failure cases, how to check a result without an estimate (7.5
in three lines), the three pitfalls that matter most (step scale in your units, kinks, `auto`).

---

## 2. Notation guide

Mathematics is set in the sans-serif face; variables italic, operators upright (`\fl`, `\ulp`), Euler's number
`\mathrm{e}`. Code names (`Central<4>`, `abs_error`) appear only in `\code{}`, never in mathematics. Primes up to
$f'''$, then $f^{(4)}$ and $f^{(N)}$; the book writes "fourth derivative" in prose and $f''''$ only where the header's
own wording is quoted.

### 2.1 Existing macros (`macros.tex`)

| Symbol | Macro | Meaning |
|---|---|---|
| $\dd$ | `\dd` | the differential, upright |
| $\dd f/\dd x$ | `\deriv{f}{x}` | the first derivative as a quotient |
| $\dd^n f/\dd x^n$ | `\nderiv{f}{x}{n}` | the n-th derivative as a quotient |
| $\mathcal{O}(h^P)$ | `\order` | order of magnitude (big O) |
| $\varepsilon$ | `\epsmach` | machine epsilon: `std::numeric_limits<T>::epsilon()`, the spacing of the numbers just above 1 ($2^{-52}$ for double). The library's rounding model uses $\varepsilon$, not the unit roundoff $\varepsilon/2$; Appendix A says so |
| $\operatorname{ulp}(x)$ | `\ulp` | unit in the last place: the spacing of the numbers at $x$ |
| $\operatorname{fl}(\cdot)$ | `\fl` | the rounded result of an operation |
| $N$ | `\dorder` | the derivative order, 1 to 4 (the template argument N) |
| $P$ | `\aorder` | the accuracy order of a rule: error $\mathcal{O}(h^P)$ (the template argument P) |
| $S$ | `\ampl` | rounding amplification $\sum_k |w_k|$ (equivalently $\sum_k |c_k| / d$ with the integer numerators) |
| $C$ | `\errc` | error constant: leading truncation term $C h^P f^{(P+N)}$, $C = |M_{N+P}| / (N+P)!$ with the rational weights |
| $s$ | `\scale` | $\max(|x|, 1)$, the scale of the default steps |
| | `\code{...}` | a C++ name or expression, monospace |
| | `\nxx` | Numerixx |
| | `\val<Name>` | a generated number (`data/values.tex`); typeset with `\num` or `\SI` (section 8) |

### 2.2 Symbols without a macro (plain mathematics)

| Symbol | Meaning |
|---|---|
| $h$ | the step (for Ridders: the initial step $h_0$; level $i$ has $h_i = h_0 / 2^i$) |
| $x$, $f$, $f(x)$ | the point and the function; `fn` in code |
| $n_k$, $w_k$ | node $k$ (an integer offset, the sample is at $x + n_k h$) and its rational weight: $f^{(N)} \approx \sum_k w_k f_k / h^N$, moments $M_N = N!$. Derivations use $w_k$ |
| $c_k$, $d$ | the integer numerators and their common denominator, $w_k = c_k / d$, as the library stores them and as Chapter 9 tabulates them; with $c_k$ the moments read $M_N = N!\,d$ |
| $x_k$, $f_k$ | $x + n_k h$ and $f(x_k)$ |
| $K$ | the largest $|n_k|$ of a rule |
| $r_i$ | the rounding floor of level $i$ of Ridders' tableau |
| $e_{i,j}$ | the estimate of a tableau entry, $e_{i,j} = \max(|A_{i,j} - A_{i-1,j-1}|, r_i)$; the result's `abs_error` is $\max(e_{\mathrm{best}}, |A_{\mathrm{best}} - A_{i,i}|)$ with $A_{i,i}$ the latest diagonal entry |
| $\xi$ | the intermediate point of Taylor's remainder |
| $T$ (mathematics) | temperature in the thread example; the C++ type parameter is always `\code{T}` |
| $u$ | not used: the book speaks of $\varepsilon$ only |
| $w$ | the scale on which $f$ varies, in the units of $x$; $s$ is the library's guess for it (round 1, T5: Chapter 3 derives the best step with $w$ and introduces $s$ as the guess) |

### 2.3 Proposed macros (editor adds them to `macros.tex` at stage 2; authors mark any others `%% NEW-MACRO`)

| Symbol | Macro | Meaning | First use |
|---|---|---|---|
| $\lvert x \rvert$ | `\abs{x}` | absolute value | ch. 3 |
| $h^\ast$ | `\hbest` | the best step of the error model | ch. 3 |
| $E_{\mathrm{t}}(h)$, $E_{\mathrm{r}}(h)$ | `\etrunc`, `\eround` | truncation and rounding error; $E(h)$ their sum | ch. 3 |
| $D(h)$ | `\Dc{h}` | the base difference (the `Central<2>` rule of order N) at step $h$; $D_i = D(h_i)$ | ch. 2, 6 |
| $A_{i,j}$ | `\Aij{i}{j}` | Ridders' tableau entry: level $i$, column $j$ | ch. 6 |
| $\hat{e}$ | `\errest` | the error estimate `abs_error` | ch. 6 |
| $M_j$ | `\mom{j}` | the $j$-th moment $\sum_k w_k n_k^j$ of a rule | ch. 4 |
| $d$ | `\den` | the common denominator of the integer numerators (italic $d$; the differential is the upright `\dd`) | ch. 4 |
| $\eta$ | `\noise` | the relative accuracy of $f$ when it is worse than $\epsmach$ (a solver tolerance, a fit's residual); it replaces $\epsmach$ in the rounding branch of the model | ch. 3 |
| $\Omega$ | `\fmax` | the largest finite value of the type, `std::numeric_limits<T>::max()` | ch. 7, App. A |
| $n_{\mathrm{f}}$ | `\nevals` | number of evaluations of fn | ch. 6 |
| $H(T)$ | `\enth` | molar enthalpy of the thread, J/mol | ch. 1 |
| $c_p$ | `\cp` | molar heat capacity, J/(mol K) | ch. 1 |
| $R$ | `\Rgas` | the gas constant (sensitivity example) | ch. 7 |

Units: `siunitx` as loaded by the preamble: `\si{\kelvin}`, `\si{\joule\per\mole}`, `\si{\joule\per\mole\per\kelvin}`,
`\SI{\valCpDefaultStep}{\kelvin}`.

---

## 3. Sidebars

At most three per chapter; each self-contained and never required to follow the text. Part II refers to them.

| Id | Title | Ch. | Contents | Why a sidebar |
|---|---|---|---|---|
| S1 `sb:alternatives` | When automatic differentiation or the complex step is the better tool | 1 | Forward-mode AD (exact to rounding, needs a generic fn and a number type), the complex step (one evaluation, needs an analytic fn of a complex argument, a form of AD), contour methods in one sentence; why `nxx::deriv` takes a black-box `fn(T)` and defers these | Most readers only need the conclusion in 1.2 |
| S2 `sb:remainder` | Where the constant 1/6 comes from | 2 | Taylor's theorem with the Lagrange remainder; the subtraction that leaves $h^2 f'''(\xi)/6$; why the next term is $h^4$ | A proof; the text only needs the result |
| S3 `sb:ulp` | Machine epsilon, ulps and the rounding of $x + h$ | 3 | $\varepsilon$ against unit roundoff; ulp; $x + h$ rounds to a neighbour, so the points are not $2h$ apart; the resolution of $x$ and the `non_finite` it causes; the spacing differs on the two sides of a power of two | Floating-point detail a user does not need |
| S4 `sb:actual-distance` | Dividing by the distance actually evaluated | 3 | The two-node rules divide by $\fl(x+h) - \fl(x-h)$ as evaluated, which removes the point-rounding error; that subtraction is exact by Sterbenz's lemma when $\fl(x+h) \le 2 \fl(x-h)$ (roughly $h \le |x|/3$), and otherwise costs one rounding of $\varepsilon/2$ relative, which is harmless; multi-node rules divide by $d h$ instead and why (4x worse otherwise); the point-rounding term $S/2\, \varepsilon |x f'| / h$ they keep | Implementation decision with a lemma |
| S5 `sb:step-constant` | The error model behind the default step, and why its constant is dropped | 3 | Minimising $E(h)$; $h^\ast = (N S \varepsilon |f| / (P C |f^{(P+N)}|))^{1/(P+N)}$; the constant between 1.4 and 2.0 for $f'$ costs at most 1.63 (`Central<10>`, generated as `\valStepCostCentralTenFirst`; the header rounds it to 1.64, D-16); for the higher orders 3.5 and 11 for `Central<2>` ($f''$, $f''''$) and 15 for `Central<4>` and $f''''$ by the model (PHASE4_API.md; `\valStepCost...`), yet half the default was better on test functions, so one convention; proof in App. B.4 | Derivation and a design argument |
| S6 `sb:moments` | The moment conditions and the Lagrange generator | 4 | $M_j = 0$ for $j \ne N$ up to $N+P-1$, $M_N = N!\,d$; the Lagrange formula; exact integer arithmetic; the 24-bit check that keeps weights exact in float; the `static_assert` for narrower types | Algebra the user never sees |
| S7 `sb:richardson-disguise` | `Central<4>` is Richardson extrapolation in disguise | 4 | $(4 D(h/2) - D(h)) / 3$ expands to the four-point rule; fixed-level extrapolation is a stencil; what Ridders adds | Forward link to ch. 6 that would interrupt ch. 4 |
| S8 `sb:fourth-sensitivity` | Why the fourth derivative is so sensitive | 5 | Rounding grows as $1/h^N$; each halving multiplies it by $2^N$; the model error $\varepsilon^{P/(P+N)}$ per N for double and float; why float $f''''$ is noise and $f'''$ still usable; the step sensitivity 16 | Quantitative detail behind one sentence of text |
| S9 `sb:actual-offsets` | The three-point $f''$ from the actual offsets | 5 | $((f_+ - f_0)/h_+ - (f_0 - f_-)/h_-) / ((h_+ + h_-)/2)$ with $h_+ = \fl(x+h) - x$ and $h_- = x - \fl(x-h)$, which are exact for every $0 < h \le |x|$ (Sterbenz); exact for quadratics whatever the rounding of the points; the prototype's 8.7 against 1.4 for float; $h^2$ never formed (underflow) | Implementation decision |
| S10 `sb:difference-form` | Burden's difference form of the tableau | 6 | $A_{i,j} = A_{i,j-1} + (A_{i,j-1} - A_{i-1,j-1}) / (4^j - 1)$; exact when the parents agree; cannot overflow below $\Omega$ where the product form did; only the previous row is kept, filled from `T(0)` | Numerical detail of the implementation |
| S11 `sb:floor-coverage` | How Ridders stops, and how often its estimate holds | 6 | The rules, moved here from the text: $e_{i,j} = \max(|A_{i,j} - A_{i-1,j-1}|, r_i)$ with the floor $r_i$ (why the minimum of noisy differences is biased low); the two-level stall; the stop on exact agreement; the final check against the latest diagonal; the untrusted-result restart (factors 4 and 64); the level cap. Then the random study of F14 with its scope stated (nine smooth functions on $[0.2, 3]$, no noise, no singularity within $h_0$): coverage, median and worst ratio, by N (measured: 99.7 % and 99.9 % for N = 1 and 2, worst ratios 1.8 and 1.4); the tail of the evaluation histogram for N = 2 (restarts of $\log$ and $\sqrt{\cdot}$ where $h_0 = s/4$ exceeds $x$; $1/x$ there has its pole between finite samples and takes more levels instead, without restarting); the design record's figures named as such | Rules and statistics most readers only need the headline of |
| S12 `sb:aliasing` | Aliasing: when the samples lie | 6 | $\sin(64 \pi x)$ at 0 vanishes at every sample of $h = 1/32, 1/64$; the untrusted-result test (at most two correct bits and far above the floor) restarts; exactly periodic samples that lie on a polynomial fool every sampling method | Rare failure mode |
| S13 `sb:expression-templates` | Expression templates and user-defined real types | 7 | What `nxx::real` requires; why `a + b` may not be a `T`; `T const` everywhere, `static_cast<T>` of fn's result; `T(*r)` and `-> T` in `transform`; Boost's arithmetic operator templates do not consider the conversion of an `Estimate<T>` to `T` (its comparisons do compile and are true; D-14) | C++ detail for a minority of readers |
| S14 `sb:constexpr-overflow` | Overflow guards and constant evaluation | 7 | Why GCC rejects overflow and NaN in constant expressions; `sum_fits`, `product_fits`, the input check $h \le \Omega/2Q$, $|x| \le \Omega - Q h$; saturation of the floors; the one exception (the final quotient) | Implementation detail |

Fourteen sidebars, about 9 pages in all, inside the chapter budgets.

---

## 4. Figures and data sets

### 4.1 Conventions

Every plot reads `data/*.csv` with `\addplot table[col sep=comma, comment chars={\#}]` in an `eqtufte` axis. Errors are
plotted as relative errors (dimensionless) on logarithmic axes, so that methods and orders compare; the text quotes
absolute values with units where they help. The default step is marked by a thin `muted` vertical line labelled
"default step", its abscissa taken from the CSV (a one-row file or a column) rather than typed. Curves are labelled
directly in the plot, the one thing to notice in `accent`. Diagrams (no data) are TikZ in the same visual language.
Each caption's first sentence states what to see; the data file is named in `\code{}` at the end.

Widths: text-width `figure` with `\sidecaption` (95 x 45 mm axes as in the layout test), `marginfigure` (53 mm),
`figure*` full width (174 mm) only for F14.

### 4.2 Figures

| Id | Ch. | Shows; what the caption says | Data set | Width |
|---|---|---|---|---|
| F1 `fig:secant` | 1 | A smooth curve with the tangent at $x$ and the secants through $(x, x+h)$ and $(x-h, x+h)$: the central secant is nearly parallel to the tangent | none (diagram) | margin |
| F2 `fig:fd-slopes` | 2 | Relative error of `Forward<1>` and `Central<2>` against $h$ on the thread, large steps only (1e2 to 1e-2 K, where truncation dominates; below 1e-2 K the central curve turns up, which the caption names as Chapter 3's subject, decision S2-14): slopes 1 and 2 on log axes; "halving the step halves one error and quarters the other" | DS1, columns `h, forward1, central2` | text |
| F3 `fig:v-curve` | 3 | The V: `Central<2>` relative error against $h$ over the whole range, the model $E(h)$ dashed, the default step marked; "the two branches, and the default step near the bottom" | DS1, columns `h, central2, model_central2`; DS1b for the marker | text |
| F4 `fig:step-scale` | 3 | The V for a narrow feature at a large $|x|$: the default step $\varepsilon^{1/3} s$ sits far right of the minimum, and the step $\varepsilon^{1/3} w$ built from the feature's width $w$ lands at it; "the default step assumes $f$ varies on the scale $\max(\lvert x \rvert, 1)$; replace $s$ by the scale on which $f$ varies" | DS4; both steps marked from DS4's values | text |
| F4b `fig:step-scale-small` | 3 | The same at a small $|x|$ in the user's units ($y \ln y$ at $y = 10^{-6}$): the default step exceeds $x$ and crosses zero, so the central rule fails, while $\varepsilon^{1/3} y$ gives the expected accuracy; "the 1 in $\max(\lvert x \rvert, 1)$ is in your units" | DS4b | margin |
| F5 `fig:float-double` | 3 | The V of `Central<2>` for float and double on $\mathrm{e}^x$ at 1: the floor sits at about $\varepsilon^{2/3}$ of each type | DS3 | margin |
| F17 `fig:noise` | 3 | The V of `Central<2>` on $\mathrm{e}^x$ at 1 with the function's values perturbed by a relative noise $\eta = 10^{-10}$ and $10^{-8}$ against the clean curve: the rounding branch lifts to $\eta |f| / h$, the minimum moves right to about $(3\eta)^{1/3}$ and up to about $\eta^{2/3}$; "a solver converged to $10^{-8}$ gives about five digits, at a step a thousand times the default" | DS19 | margin |
| F6 `fig:order-ladder` | 4 | Relative error against $h$ for `Central<2>` to `Central<10>`, $\mathrm{e}^x$ at 1: steeper truncation slopes, rounding branches that nearly coincide ($S$ from 1 to 2.3), and minima that slide up the common rounding line as the best step grows; "a higher order buys accuracy by letting the step grow; the minimum falls from about $\varepsilon^{2/3}$ to about $\varepsilon^{10/11}$, three and a half decades" | DS2, columns `h, central2 ... central10` | text |
| F7 `fig:one-sided` | 4 | `Forward<2>` against `Central<2>` and `Forward<4>` against `Central<4>`: same slopes, higher floors (larger $S$) | DS2, columns `h, forward2, central2, forward4, central4` | margin |
| F8 `fig:default-step-accuracy` | 4 | Relative error at the default step against P for Central and Forward, $\mathrm{e}^x$ at 1: diminishing returns beyond 4 or 6 for double | DS7 | margin |
| F9 `fig:stencil-nodes` | 4 | The nodes of `Central<4>`, `Forward<2>` and `Backward<2>` for N = 1 and N = 2 on a line, node 0 filled where the rule uses it | none (diagram) | text |
| F10 `fig:v-curve-orders` | 5 | The V of `Central<2>` for N = 1 to 4, $\mathrm{e}^x$ at 1: the rounding branch steepens as $1/h^N$ and the minimum rises; "why the fourth derivative gives four digits where the first gives eleven" | DS5 | text |
| F11 `fig:order-accuracy` | 5 | Relative error at the default step against N for `Central<2>`, `Central<4>` and Ridders, double and float; the model $\varepsilon^{P/(P+N)}$ as a guide line; the float $f''''$ point marked "warning" | DS6 | text |
| F12 `tab:ridders-tableau` | 6 | The tableau for $c_p$ at $T_0$ from $h_0 = T_0/32$: the levels Ridders used (seven), the entries and, in a second body, their relative errors with the floor; the returned entry in `accent` | DS8 (`data/table-ridders-tableau.tex`, `table-ridders-tableau-errors.tex`) | full width (`table*`): nine to ten columns of twelve-digit numbers measure about 150 mm |
| F13 `fig:ridders-convergence` | 6 | Error of $A_{i,j}$ against level $i$ over 14 to 16 levels, one curve per column $j$, the rounding floor $r_i$ dashed: column 0 descends slowly (its truncation at $h = 1$ K is still far above the floor) while each higher column reaches the floor within a few levels; the floor doubles per level, so the entries then drift up; Ridders' stop (from the counted evaluations) and its returned entry (matched bit for bit against the tableau) marked; "the extrapolation reaches in a few levels what the base difference would need a dozen for; after that, more levels only add rounding" | DS8 (csv) | text |
| F14 `fig:ridders-coverage` | 6 | Scatter of the estimate against the actual error over the random study (N = 1; N = 2 in `muted`), the diagonal drawn: points above it are covered; absolute values (`actual`, `estimate`), since relative ones blow up near a zero of $f^{(N)}$; cases with an actual error of exactly 0 (exact stops on $x^5$) are clipped to $\varepsilon |f^{(N)}|$ and counted in the caption; the caption states the study's scope (nine smooth functions on $[0.2, 3]$, no noise, no singularity within $h_0$); margin histogram of evaluations (F14b `fig:ridders-cost`) with the N = 2 restart tail named (log and sqrt below $x = 0.25$) | DS9 (`ridders-coverage-N1.csv`, `-N2.csv`) | full width; histogram margin |
| F15 `fig:ridders-edge` | 6 | Evaluations and relative error of Ridders for $\log x$ at $x = 2^{-k}$ against $k$: restarts cost about two evaluations per halving; the failure beyond the measured limit | DS10 | margin |
| F16 `fig:newton` | 7 | Newton's residual against iteration for $T$ from $H(T) = H^\ast$ with a `Central<2>` slope: quadratic convergence until rounding | DS11 | margin |
| F18 `fig:kink` | 7 | Relative error of `Central<2>`, `Central<4>` and Ridders (value and estimate) against the distance of $x$ from the range break at 500 K, where the two published Shomate fits do not meet and $H$ jumps (decision S2-2): inside a method's span the error of a straddling stencil grows as $1/h$ whatever P (the data show a jump; for a kink it would be $\order(1)$ or $\order(h)$, said as the general statement), the rule returns it without failing, and only Ridders' estimate shows it; "evaluable is not smooth: the estimate flags it, the stencil does not" | DS20 | margin |

### 4.3 Data sets (`tools/generate_data.cpp`)

All sets below exist (stage 2); the template data set went with the layout test at stage 4, and the two hand tables
(HT) were added then. Common rules: double unless stated; steps $h = 10^{k/4}$ (four per
decade); a result that is not finite is written as `nan` in its cell (pgfplots skips it; the figure may set
`unbounded coords=jump`); reference values from the closed forms coded in the generator or the shared headers (the
Shomate derivatives, $\mathrm{e}^x$, and the closed forms of `tests/test_orders.cpp` ported: `sin_n`, `log_n`,
`runge_n`, plus those of $\sqrt{x}$, $1/x$, $x^5$, $\mathrm{e}^{-x^2}$ and $\tanh x$ to order 4), never from Numerixx;
the random study uses `std::mt19937` with the fixed seed 20261023 written into the CSV header and draws its points
from the generator's raw 32-bit output, so every standard library gives the same points; a float fourth derivative
is computed inside the same pragma the tests use, since the generator builds with warnings as errors. The
"accuracy at the default step" sets (DS6, DS7, DS15) report the median over eleven points $x = 0.5, 0.55, \ldots, 1$
(all with $s = 1$, so the step is the same), because the rounding part of the error at a single point is a matter of
luck by a factor of ten (decision S2-3); the CSVs also carry the smallest and largest value. Output kind
`data/table-<name>.tex`: the body rows of a `tabular`, first line the generator's header as a `%` comment (so
`tools/check.py` accepts it), read with `\input`; numbers in it are wrapped in `\num{}`. The generator reads the
weight tables through `nxx::deriv::detail::make_table<N>(S{})`, as the tests do (decision A3), and checks each
against the moment conditions before writing it, so a table in the book is self-verified; it also checks the thread
function against the evidence reviewer's exact values, Ridders' returned entry against the tableau bit for bit, and
rejects a `\val` macro defined twice, and returns 1 if any check fails. Files are written with LF line ends on every
platform. Macro naming: `<X>RelError` $= \abs{\text{value} - \text{exact}} / \abs{\text{exact}}$; `<X>Actual` the
absolute error; `<X>Estimate` the returned `abs_error` (absolute, in the derivative's unit); `<X>RelEstimate`
$= $ `abs_error` $/ \abs{\text{exact}}$; `<X>Evals` the evaluations of fn; `<X>Value` the returned value to twelve
digits.

| Set | File | Columns | Computation |
|---|---|---|---|
| HT | `table-hand-cube.tex`, `table-hand-calculator.tex`, `values.tex` | the pencil-and-paper tables of 2.3 ($x^3$ at 1, $h = 0.1$ and $0.05$: $f(1 \pm h)$, `Forward<1>` and `Central<2>` values and errors to seven digits) and 3.2 (a six-digit calculator differentiating $\mathrm{e}^x$ at 1 for $h = 0.1$ to $10^{-4}$: every value of $f$ and the quotient rounded to six significant digits, the signed error); `\valHandCalcExact`, `\valHandCalcErrorRatio`, `\valHandCalcLastDiff` | Generated at stage 4 (decision S4-3) so that no computed digit is typed; the cube rows are the library's values printed to seven digits, which equal the exact decimal hand calculation; the calculator is simulated by rounding to six digits |
| E6 | `values.tex` | `\valEdgeForwardTwoRelError`, `\valEdgeForwardFourRelError` | `Forward<2>` and `Forward<4>` at 500 K with the NaN wrapper of `code/04-stencils-edge.cpp`; the generator checks that `Central<2>` fails there |
| DS1 | `step-error-cp.csv` | `h, forward1, central2, central4, model_forward1, model_central2` | Thread $H(T)$ at $T_0$, N = 1; $h$ from 1e2 to 1e-13 K; relative error against the closed-form $c_p$; the model $(C h^P \lvert H^{(P+1)} \rvert + S \varepsilon \lvert H \rvert / h) / c_p$ with $C = 1/2, S = 2$ and $C = 1/6, S = 1$ |
| DS1b | `values.tex` | `\valCpExact`, `\valCpDefaultStep`, `\valCpDefaultRelError`, `\valCpBestStep`, `\valCpBestRelError` (the minimum of DS1's `central2`), `\valCpRiddersRelError`, `\valCpRiddersEstimate`, `\valCpRiddersEvals`, `\valCpRiddersActual`, `\valCpCentralFourRelError`, `\valCpCentralSixRelError`, `\valCpSecondRelError`, `\valCpRiddersSecondRelError` | Thread at $T_0$: the default-step results of E1, E5, E7, E9, counted evaluations |
| DS2 | `step-error-exp.csv` | `h, forward1, forward2, forward4, central2, central4, central6, central8, central10` | $\mathrm{e}^x$ at 1, N = 1, $h$ from 1 to 1e-13 |
| DS3 | `step-error-float.csv` | `h, error_float, error_double` | `Central<2>`, $\mathrm{e}^x$ at 1, $h$ from 1 to 1e-8, float and double; `\valEpsFloat`, `\valEpsDouble`, `\valFloatBestRelError`, `\valDoubleBestRelError` |
| DS4 | `step-scale.csv` | `h, error` and `\valScaleDefaultStep`, `\valScaleDefaultRelError`, `\valScaleFitStep`, `\valScaleFitRelError`, `\valScaleRiddersUserRelError`, `\valScaleRiddersUserEstimate` | $f(x) = \tanh((x - x_0)/w)$, $x_0 = 1000$, $w = 0.5$, at $x = x_0 + 0.2$ (a narrow feature at a large $\lvert x \rvert$), `Central<2>`, $h$ from 10 to 1e-12; the fitted step $h = \varepsilon^{1/3} w$ (the model optimum to within a factor 1.1) and Ridders from $h_0 = w/4$ |
| DS4b | `step-scale-small.csv` | `h, error` and `\valSmallDefaultStep`, `\valSmallDefaultStatus`, `\valSmallFitStep`, `\valSmallFitRelError` | $f(y) = y \ln y$ at $y = 10^{-6}$ (defined for $y > 0$ only): `Central<2>` with the default step $\varepsilon^{1/3}$, which exceeds $y$ and samples below zero (`non_finite`), and with $h = \varepsilon^{1/3} y$; $h$ from $10^{-5}$ to $10^{-20}$ for the curve |
| DS5 | `step-error-orders.csv` | `h, n1, n2, n3, n4` | `Central<2>`, $\mathrm{e}^x$ at 1, N = 1 to 4, $h$ from 1 to 1e-13 |
| DS6 | `order-accuracy.csv` | `N, c2_double, c4_double, rid_double, rid_est_double, rid_evals_double, c2_float, c4_float, rid_float, rid_est_float, rid_evals_float` and `\valFloatFourthCentralRelError`, `\valFloatFourthRiddersRelError`, `\valFloatFourthRiddersRelEstimate`, `\valFloatThirdRiddersRelError`, `\valDoubleFourth...`, and the single point $x = 1$ as `\valFloatFourthRiddersAtOne...`, `\valFloatFourthCentralAtOneRelError` for E8's walk-through | $\mathrm{e}^x$, default steps, N = 1 to 4, double and float: medians over the eleven points (relative error, relative estimate, evaluations) |
| DS7 | `default-step-accuracy.csv` | `P, central, central_min, central_max, forward, forward_min, forward_max` and `\valExp{Central,Forward}<P>RelError` | $\mathrm{e}^x$, N = 1, relative error at the default step for `Central<2..10>` and `Forward<1..6>`: median, smallest and largest over the eleven points (`nan` where a family lacks P) |
| DS8 | `ridders-tableau.csv`, `table-ridders-tableau.tex` | `i, j, A, error, floor` | Thread $c_p$ at $T_0$: base differences $D_i$ = `derivative<Central<2>>(H, T0, h_i)` for $h_i = (T_0/32)/2^i$, $i = 0..15$ (16 levels, so that column 0 reaches its own minimum), the floor $r_i$ from the two values as in the tests' `floor_at`; columns by the difference form computed in the generator (the formula of S10); error against the closed form. Ridders' stop is taken from the evaluation count of a counting wrapper around `derivative<Ridders>(H, T0)` (level = count / 2), and the entry it returned is found by matching its `value` bit for bit against the tableau (both come from `detail::difference` at the same points), not by re-implementing the selection; `\valTableauStopLevel`, `\valTableauBestLevel`, `\valTableauBestColumn`. The `.tex` holds the first six levels as a triangular table with the returned entry flagged |
| DS9 | `ridders-coverage-N1.csv` .. `N4.csv` (one file per N, so a figure needs no row filter) | `fn, name, x, exact, actual, estimate, rel_actual, rel_estimate, evaluations, restarted` | Nine functions with closed-form derivatives to order 4 ($\mathrm{e}^x$, $\sin x$, $\log x$, $\sqrt{x}$, $1/x$, $1/(1+25x^2)$, $x^5$, $\mathrm{e}^{-x^2}$, $\tanh x$), 300 points each uniform in $[0.2, 3]$ (the same points for every N), N = 1 to 4 (3 and 4 for the cost table T15; F14 plots 1 and 2), default steps; absolute and relative actual error and estimate (relative values blow up near a zero of $f^{(N)}$, e.g. $\cos x$ at $\pi/2$, $\mathrm{e}^{-x^2}$'s second derivative at $1/\sqrt{2}$, so F14 plots absolute values or clips); an actual error of exactly 0 is written as 0. `restarted` is detected from a log of the calls: a non-finite value, or a level step that does not decrease (an untrusted result restarts from half the previous initial step and re-evaluates points). Measured: for N = 2, $\log$ and $\sqrt{\cdot}$ restart for $x < 0.25$ (NaN below 0), but $1/x$ does not: the pole lies between finite samples and the tableau simply takes more levels, accurately (the plan's earlier expectation was wrong). Summary macros per N (`First` .. `Fourth`): `\valRiddersCoverage`, `\valRiddersUncovered`, `\valRiddersMedianRatio`, `\valRiddersWorstRatio`, `\valRiddersEvals{Min,Median,Max}`, `\valRiddersExactCases`, `\valRiddersRestarted`, `\valRiddersFailed`; and `\valRiddersCases`, `\valRiddersFunctions`, `\valRiddersPointsPerFunction` |
| DS10 | `ridders-edge.csv` | `k, x, ok, covered, rel_error, rel_estimate, evaluations` | $\log x$ at $x = 2^{-k}$, $k = 2..40$, plus the row $x = 0.01$ of E10 (written with $k = \log_2(1/x)$), Ridders with the default initial step. Three readings of the limit, all measured: `\valRiddersEdgeLimit` ($2^{-34}$, the smallest $x$ that still returns a value; below it the level budget runs out, which matches the derivation of D-6), `\valRiddersEdgeTenDigits` ($2^{-30}$, the smallest $x$ whose result still has ten correct digits; the accuracy degrades as the budget leaves fewer levels for convergence), `\valRiddersEdgeCapExponent` (27, the first $k$ at which all 64 evaluations are used); `\valRiddersEdgeEvalsAtOneHundredth` and `\valRiddersEdgeAtOneHundredth...` from the $x = 0.01$ row |
| DS11 | `newton.csv` | `iteration, T, residual, error` | Newton for $T$ with $H(T) = H(T^\ast)$, $T^\ast$ inside the range, start 200 K below, slope by `derivative(H, T)`; `\valNewtonIterations`, `\valNewtonFinalError` |
| DS12 | `table-weights-central-N1.tex` .. `N4`, `table-weights-forward-N1.tex` .. `N4` | rows: P, nodes, numerators, denominator, $S$, $C$ (reduced fraction), leading error term | From `detail::make_table`, checked against the moment conditions in the generator; Backward stated as the mirror rule |
| DS13 | `table-default-steps.tex` | rows: family, P, N, exponent, step for double at $\lvert x \rvert \le 1$ (`root_epsilon`), model error $\varepsilon^{P/(P+N)}$, evaluations | Public `nxx::math::root_epsilon` and the node counts |
| DS14 | `table-ridders-constants.tex` | the constants of Chapter 9.5 | From `detail::ridders_levels`, `detail::ridders_divisors` and the documented factors |
| DS15 | `accuracy-matrix.csv` (long form: `method, type, N, median, min, max`), `table-accuracy-matrix.tex` | rows: method (`Central<2>`, `Central<4>`, `Central<6>`, `Forward<2>`, `Forward<4>`, `Ridders`); columns: N = 1..4 for float, double, `cpp_bin_float_50` | $\mathrm{e}^x$, default steps, the median relative error over the eleven points; the float $f''''$ cells computed under the pragma and marked `*`; needs Boost.Multiprecision in the generator (decision A2; `\valAccuracyMatrixHasMultiprecision` says whether the cells were computed). The headline digits of 1.5 are $\lfloor -\log_{10} \rfloor$ of the double `Central<2>` medians: `\valHeadlineFirstDigits` (11), `\valHeadlineFourthDigits` (4), and the same for Ridders |
| DS16 | `table-cost.tex` | rows: method; columns: evaluations per N (fixed for stencils), Ridders min/median/max from DS9 and the caps | Node counts and DS9 |
| DS17 | `values.tex` | `\valEpsFloat`, `\valEpsDouble`, `\valMaxDouble`, `\valDenormMinDouble`, `\valDigitsFloat`, `\valDigitsDouble`, `\valUlpAtOneThousand` | `std::numeric_limits` (Appendix A, T16, S3) |
| DS18 | `values.tex` | `\valSensitivityRateRelError`, `\valSensitivityRateEstimate`, `\valSensitivityEaRelError`, `\valSensitivityScaled`, `\valJacobianEvals`, `\valMpCentralSixRelError`, `\valMpRiddersRelError`, `\valMpRiddersEstimate`, `\valHeadlineFirstDigits`, `\valHeadlineFourthDigits`, `\valCpRelativeFactor`, `\valCpSecondRelativeFactor` | The results E13 and E14 print, recomputed by the generator; the headline digits of 1.5 as $-\log_{10}$ of the `Central<2>` entries of DS15 for double; the factors $|H| / (s |H'|)$ and $|H| / (s^2 |H''|)$ at $T_0$ |
| DS19 | `step-error-noise.csv` | `h, clean, noise_ten, noise_eight` and `\valNoise{Ten,Eight}BestStep`, `\valNoise{Ten,Eight}BestRelError`, `\valNoiseEightDefaultRelError`, `\valNoise{Ten,Eight}FitStep`, `\valNoise{Ten,Eight}FitRelError`, `\valNoiseEightModelBest...`, `\valNoiseRidders...` (the single point $x = 1$) and `\valNoiseRiddersEstimateRatio`; the coverage study `\valNoise{Clean,Ten,Eight}Ridders{Coverage,MedianRatio,WorstRatio,MedianRelError}` | `Central<2>` on $\mathrm{e}^x$ at 1, $h$ from 1 to 1e-13, with fn's value multiplied by $1 + \eta u$, $u$ uniform in $[-1, 1]$ from `std::mt19937` seeded per call from the bits of $x$ (`code/noise.hpp`, shared with E17, so the perturbation is identical), $\eta = 10^{-10}$ and $10^{-8}$; the best measured step and error per $\eta$. Ridders on the noisy function: at the single point $x = 1$ the estimate happens to cover the error (ratio 0.17), so the text must not argue from one point; over the 300 random points of the study the coverage falls from 100 % (clean) to about 89 % ($\eta = 10^{-10}$) and 91 % ($10^{-8}$) with worst ratios of 67 and 42: that is the measured form of the header's "optimistic for a noisier fn" (section 3.6, E17) |
| DS20 | `kink.csv` | `delta, central2, central4, ridders, ridders_estimate, ridders_evals` and `\valKinkJump{Enthalpy,Cp,CpSlope,CpSlopeRel}`, `\valKink{Inside,Outside}Delta`, `\valKink{Inside,Outside}{Central,CentralFour}RelError`, `\valKink{Inside,Outside}Ridders...`, and the plan's `\valKinkCentralRelError`, `\valKinkRiddersRelError`, `\valKinkRiddersEstimate` (the inside point) | The published enthalpy of nitrogen across the junction of its two lower Shomate ranges at $x_{\mathrm{k}} = 500$ K (`nitrogen_H_joined` in `code/shomate.hpp`: the lower fit below 500 K, the thread's fit from 500 K on; decision S2-2): a range break where $H$ jumps by 0.77 J/mol, $c_p$ by $-7.5 \cdot 10^{-5}$ J/(mol K) and $\dd c_p/\dd T$ by $-17$ %. $x = x_{\mathrm{k}} + \delta$ for $\delta$ from $10^{-7} s$ to $s$, default steps; relative errors against the closed form of the range $x$ lies in. Measured at $\delta = 1$ mK (inside every default span): `Central<2>` 4.3 (430 %, the jump over $2h$), `Central<4>` 0.041, Ridders $1.9 \cdot 10^{-3}$ with a relative estimate of $6.2 \cdot 10^{-3}$ (0.18 J/(mol K), eleven orders above its usual value: the flag); at $\delta = 50$ K every error is below $10^{-11}$. For a jump the stencil error is $\order(1/h)$, for a kink $\order(1)$ or $\order(h)$; F18's caption says which the data show |
| DS21 | `values.tex`, `step-error-solver.csv` (`h, bisection, fixed_point`) | `\valSolver{Reynolds,Roughness,Tolerance,Friction,FrictionNoise,Exact,ScaleFactor,DefaultStep,FitStep}`, `\valSolverDefaultRelError`, `\valSolverFitRelError`, `\valSolverModelFitRelError`, `\valSolverDeriv{RelError,Actual,Estimate,RelEstimate,Evals,Value,EstimateRatio}` (Ridders), and `\valSolverFixedPoint...` for the comparison | E20's friction factor $f(Re)$ from the Colebrook equation ($\varepsilon/D = 10^{-4}$, $Re = 10^5$) solved by bisection to a relative tolerance $\eta = 10^{-10}$ (decision S2-4: bisection's output is quantised at $\eta$, which is the noise of 3.6; a fixed-point iteration stopped at $\eta$ is smooth between the points where its iteration count changes, so its default-step derivative is accurate to about $\eta$ and misses the lesson, measured as $7.5 \cdot 10^{-11}$ against bisection's $5.5 \cdot 10^{-5}$); `derivative<Ridders>` and `Central<2>` with the default step and with $h = \eta^{1/3} Re$; the answer key by implicit differentiation of the Colebrook equation at the root found to machine precision by Newton's method |

Macro names are letters only (a LaTeX command name cannot contain digits: `CentralFour`, `First`; `N1` only in file
names). `data/README.md` maps every file and macro to its figure or chapter; the editor keeps it current.

### 4.4 Stage-6 extensions (the round-1 decisions, `reviews/R1-decision.md` section 4)

Added by the editor before the authors' revision; every value is in `data/values.tex` and `data/README.md` maps it.

| Set | Macros and files | Computation and use |
|---|---|---|
| HT | `\valHandCalcTruncationHundredth` | $\mathrm{e} \cdot 0.01^2 / 6$, the truncation term alone at $h = 0.01$ (3.2, R1-evidence-03) |
| DS1b | `\valCpTableValue`, `\valCpFitVsTableRelDiff` | The tabulated $c_p$ of nitrogen at 1000 K, 32.697 J/(mol K), a cited constant in the generator (`chase1998`, Part II, p. 1621), and the fit's relative difference from it (1.4, R1-evidence-19) |
| E6 | `\valEdgeRiddersStatus`, `\valEdgeRiddersEvals` (and `\valEdgeRiddersRelError` if a value is returned) | Ridders at 500 K with the NaN wrapper (4.4, R1-evidence-18) |
| DS7 | `\valDefaultStepSpreadAbove`, `\valDefaultStepSpreadBelow` | Over every rule, the largest value over the median and the median over the smallest value (4.2, R1-evidence-14) |
| DS9 | `\valRiddersMargin{Cases,QuarterCovered,SixtyFourthCovered}` (stage 7, R2-numerics-03: the uncovered cases of all four orders together, 14, and how many a quarter and a sixty-fourth of the default initial step cover, 11 and 14; Sidebar 6.2); per N: `\valRiddersZeroError...` (renamed from `RiddersExactCases`: bit-exact results, not observed stops), `\valRiddersZeroOnQuintic...`, `\valRiddersZeroFunctions...`; `\valRiddersTailCases...`, `\valRiddersTailMinEvals...`, `\valRiddersTailRestarts...`, `\valRiddersTailFunctions...` (the functions with their counts, most frequent first, as LaTeX text); where cases are uncovered, `\valRiddersUncoveredMainFunction...`, `...MainCount`, `...MainMinX`, `...MainMaxX`, `...Others` | The histogram's tail is the cases above the median evaluation count (6.5, F14b, T7); the zero-error split (Sidebar 6.2, F14); the uncovered cases' function (Sidebar 6.2, R1-numerics-12). Coverage prints two decimals where one decimal would read 100.0 with a case uncovered (R1-numerics-17e) |
| DS10 | `\valRiddersEdgeEvalsInterior`; `\valRiddersEdgeLimitRelEstimate` | $\log x$ at $x = 1/4$ (the $k = 2$ row), the comparison for $x = 0.01$ (6.4, R1-evidence-10); the relative estimate at $x = 2^{-34}$, the last $x$ that returns a value, 0.10 (6.4, S6-6) |
| DS13 | `\valMisjudgeCost{CentralTwoFirst,CentralTwoFourth,CentralFourFourth}`; `\valCapGainCentralTwelve{First,Fourth}{Double,Multiprecision}`; `\valModelDigitsCentralTwo{First..Fourth}Double` | The model's growth of the total error for a step off by a factor 2 at the balance, $\max((N 2^{-P} + P 2^N), (N 2^P + P 2^{-N})) / (P + N)$ (3.5, Sidebar 5.1, 5.3, T9); the digits a `Central<12>` would gain by the exponent alone (4.4, 9.3, 10.4, T10); $-\log_{10} \varepsilon^{P/(P+N)}$ to one decimal (5.2, T13) |
| DS14 | the untrusted-result row | "estimate at least $\lvert A \rvert / 4$ and more than 64 times the floor" (R1-library-20) |
| DS15 | `\valHeadline{First,Fourth}Digits`, `\valHeadlineRidders{First,Fourth}Digits` now rounded, not floored (values 11, 5, 14, 9); `\valExpRiddersFirstRelError` (Ridders' median relative error for $f'$ in double, 1.3e-14, the figure behind the headline 14; 6.5, S6-6); `\valModelRatioCentralBelowMax`, `\valModelRatioCentralFourthAbove`, `\valModelRatioForwardTwoAbove`, `\valModelRatioForwardFour{Third,Fourth}` | The measured medians against $\varepsilon^{P/(P+N)}$ in double, per family (10.1, T12, T13) |
| DS22 | `\valNoiseEstimate{Friction,Noisy,Clean}` (the estimate, relative to $\lvert f \rvert$), `...K` (the difference order accepted), `...Step`; `\valNoiseEstimateFrictionActual` (the root mean square deviation of the bisection output from the implicit solution over the same points); `\valNoiseEstimateNoisyEta`; `\valNoiseEstimatePoints`, `\valNoiseEstimateEvals` | E21's recipe `book::estimate_noise` (`code/noise.hpp`) at eight points spaced by $\varepsilon^{1/3} \max(\lvert x \rvert, 1)$: on the bisection friction factor at $Re = 10^5$ ($\eta = 10^{-10}$), on the noisy $\mathrm{e}^x$ of E17 at 1 ($\eta = 10^{-8}$) and on the clean $\mathrm{e}^x$ (3.6, 7.4.3, T3) |
| E22 | `\valJacobianBalanceHandMaxRelError`, `\valJacobianBalanceHandEvals`, `\valJacobianBalanceCallMaxRelError`, `\valJacobianBalanceCallEvals`, `\valJacobianBalanceStep{T,N}` | The $2 \times 2$ Jacobian of section 5's E22 by the hand-applied central stencil (two model calls per column) and by one `derivative` call per entry (7.4.2, T2) |

The shared headers gained `code/colebrook.hpp` (the friction factor of E20 and E21: `book::friction`, the
`implicit-function` snippet, `book::friction_fixed_point`, `book::friction_exact`, `book::friction_slope_exact`) and,
in `code/noise.hpp`, the recipe `book::estimate_noise` with the snippet `noise-estimate`. The generator reads both,
so the numbers of DS21 and DS22 come from the same code the examples quote.

---

## 5. Examples (`code/`)

Every listing is a snippet of a compiled file, built as `numerixx_book_<file>` and run as a test. Conventions of
`code/README.md`: markers `// <name>` and `// </name>` in column 0, lines of at most 60 characters, `main` checks what
the text claims and returns 0. Output with `std::println` (the demo builds with it on all five toolchains). Two
shared headers in `code/` (the editor's; the CMake glob builds only `code/*.cpp`): `shomate.hpp` holds the thread's
enthalpy (`book::shomate_H`, `book::shomate_cp`, `book::shomate_derivative`, the ranges `book::nitrogen_low`,
`book::nitrogen`, `book::nitrogen_high`, `book::thread_T0`, `book::junction_T`, `book::nitrogen_H_joined`) and
`noise.hpp` the noisy function (`book::noisy`); every example that uses them includes them with `#include
"shomate.hpp"`, and E1 quotes the definition from the header with `\snippetheader{shomate}{shomate}` and
`\snippetheader{shomate}{shomate-H}` (snippets `shomate`, `shomate-H`, `shomate-cp`, `noise`), so the coefficients
exist in exactly one place (section 8, item 14). The example of a multiprecision type links `boost_multiprecision`
and is skipped when that target is absent (A2). The parameters an example shares with the generator (E13: $A =
10^{10}$ 1/s, $E_a = 80$ kJ/mol, $T = 600$ K, $R = 8.31446261815324$ J/(mol K); E20: $\varepsilon/D = 10^{-4}$,
$Re = 10^5$, $\eta = 10^{-10}$, bisection on $[1, 20]$ in $1/\sqrt{f}$; E19: $\delta = 1$ mK and 50 K) are listed in
`data/README.md`; an example must use the same ones, or its printed numbers differ from the `\val` macros.

| Id | Ch. | File | Snippets; what it demonstrates | What `main` asserts |
|---|---|---|---|---|
| E1 | 1 | `01-why-first.cpp` | The Shomate $H(T)$ in J/mol quoted from the shared header (`\snippetheader{shomate}{shomate}`, `{shomate-H}`); `first` (`derivative(H, T0)` gives $c_p$); `first-ridders` (`derivative<Ridders>(H, T0)` with `abs_error`) | Relative error of the stencil at most $10^{-9}$ (measured $2.4 \cdot 10^{-11}$, the model about $4 \cdot 10^{-11}$; the margin covers libm differences); Ridders' actual error at most 6 times its estimate (the tests' criterion; measured ratio 0.73, 14 evaluations) |
| E2 | 2 | `02-taylor-steps.cpp` | `forward-vs-central`: `Forward<1>` and `Central<2>` at $h = 2$ and 1 K; prints the errors and their ratios | Forward ratio in [1.9, 2.1], central ratio in [3.9, 4.1] (truncation dominates by five orders of magnitude at these steps) |
| E3 | 3 | `03-rounding-step.cpp` | `explicit-step` ($h = 10^{-2}$ and $10^{-8}$ K); `too-small` ($h = 10^{-14}$ K at 1000 K is below the resolution of $x$) | `non_finite` for the tiny step; the default step's error below the $10^{-8}$ step's |
| E4 | 3 | `03-rounding-scale.cpp` | `narrow-feature`: the tanh model at $x_0 + 0.2$: the default step, the step $\varepsilon^{1/3} w$ from the feature's width, Ridders from $h_0 = w/4$; `small-x`: $y \ln y$ at $y = 10^{-6}$ with the default step (fails: it samples below zero) and with $\varepsilon^{1/3} y$ | Fitted-step error below the default's by at least a factor 100 (the model predicts about $10^{4}$); Ridders covered (6x) and error at most $10^{-9}$; `non_finite` for the small-$x$ default and relative error at most $10^{-9}$ for the fitted step |
| E5 | 4 | `04-stencils-orders.cpp` | `higher-order`: `Central<2>`, `<4>`, `<6>` on $H$ at $T_0$ | Each relative error at most $10^{-9}$ |
| E6 | 4 | `04-stencils-edge.cpp` | `validity-edge`: a wrapper returning NaN outside the fit's range; `Central<2>` at the lower end fails, `Forward<2>` and `Forward<4>` succeed | `non_finite` for the central rule; one-sided relative errors at most $10^{-8}$ ($(1/3 + 8)\varepsilon^{2/3}$ is about $3 \cdot 10^{-10}$) |
| E7 | 5 | `05-orders-second.cpp` | `second-derivative`: `derivative<Central<2>, 2>` and `derivative<Ridders, order::second>` for $\dd c_p / \dd T$; `order-enum` | Stencil relative error at most $10^{-6}$ (the model gives about $2 \cdot 10^{-7}$: $4 \varepsilon |H| / (h^2 |H''|)$ with $h = 0.12$ K, the factor $|H| / (s^2 |H''|)$ being about 4 at $T_0$); Ridders covered |
| E8 | 5 | `05-orders-float.cpp` | `float-fourth`: a float fourth derivative inside the silencing pragma (the file builds with warnings as errors, so the pragma is the point); the macro alternative in a comment | `has_value()`; the text quotes DS6 for the numbers |
| E9 | 6 | `06-ridders-basic.cpp` | `ridders-cp`: Ridders on $H$ at $T_0$ with a counting wrapper; prints value, estimate, evaluations | Actual error at most 6 times the estimate; evaluations at most 64 |
| E10 | 6 | `06-ridders-edge.cpp` | `edge-restart`: $\log x$ at 0.01, Ridders restarts (the default initial step 1/32 crosses 0); `Central<2>` with its default step does not | Ridders has a value, covered, evaluations above 16 and at most 64 (measured 24); the stencil's relative error at most $10^{-6}$ (measured $1.2 \cdot 10^{-7}$: the truncation $h^2 / (3 x^2)$ with the step scaled by $s = 1$, not by $x$; the text uses this as the small-$\lvert x \rvert$ cross-reference to 3.5). The plan's earlier $10^{-9}$ was wrong |
| E11 | 7 | `07-using-expected.cpp` | (Round 1: `ridders-check` added, the tolerance test on `abs_error` quoted in 7.2; the `names` snippet is no longer quoted, 8.1 names the using-declarations and the file in prose.) `expected-basics` (`if (d)`, `d.error()`, `*d`, `d->value`, `value_or`); `expected-chain` (`transform`, `and_then`); `estimate-pitfalls` (`auto` deduces `Estimate`; `if (*r)` tests the value; test `has_value()` before `*abs_error`) | `non_finite` for $\log x$ at 0 with the central rule; the chain's result; `std::same_as<decltype(e), Estimate<double>>` |
| E12 | 7 | `07-using-newton.cpp` | `newton`: $T$ from $H(T) = H^\ast$ with the slope from `derivative(H, T)` | Final error at most $10^{-6}$ K; at most 10 iterations |
| E13 | 7 | `07-using-sensitivity.cpp` | `sensitivity`: an Arrhenius rate $k(T) = A \mathrm{e}^{-E_a/(R T)}$ with stated illustrative parameters; $\partial k / \partial T$ and $\partial k / \partial E_a$ by lambdas capturing the rest, Ridders; `scaled`: the scaled sensitivity $\dd \ln k / \dd \ln E_a = (E_a / k)\, \partial k / \partial E_a$; `jacobian`: a loop over the parameters $(A, E_a)$ that fills one Jacobian column per parameter and counts the evaluations | Covered (6x) and relative errors at most $10^{-8}$ against the closed forms $k E_a / (R T^2)$ and $-k / (R T)$; the scaled sensitivity equals $-E_a / (R T)$ to $10^{-8}$; the evaluation count equals the sum over the columns |
| E14 | 7 | `07-using-multiprecision.cpp` | `multiprecision`: `cpp_bin_float_50`, $\mathrm{e}^x$ at 1 with `Central<6>` and Ridders; `T(*r)` and `-> T` | `Central<6>` relative error at most $10^{-40}$ ($\varepsilon^{6/7}$ is about $10^{-43}$); Ridders covered |
| E15 | 7 | `07-using-constexpr.cpp` | `constexpr`: `static_assert(derivative(cube, 2.0, 0.25) == 12.0625)`, a `non_finite` as a constant, Ridders exact on $x^3$ | The `static_assert`s; `main` returns 0 |
| E16 | 8 | `08-api-usage.cpp` | `api-overloads` (default step, explicit step, integer and `order` N, Ridders); `api-rejected` (`requires`-expressions showing that `Central<3>`, order 5, `bool`, an `int` argument do not compile) | The `static_assert`s; `main` calls each accepted form and returns 0 |
| E17 | 3 | `03-rounding-noise.cpp` | (Round 1: the snippet is no longer quoted; 3.6 shows the call as inline code and names the file as the source of the measurement.) `noisy-function`: $\mathrm{e}^x$ with its value perturbed by a deterministic relative noise $\eta = 10^{-8}$ (`book::noisy` from `code/noise.hpp`, the same perturbation as DS19); the default step, the step $\eta^{1/3} s$, and Ridders with its estimate; `noisy-coverage`: Ridders over a handful of points, counting how often the estimate covers the error | The fitted step's error below the default's (measured $1.2 \cdot 10^{-6}$ against $9.4 \cdot 10^{-4}$); the fitted step's relative error at most $10^{-4}$ (the model gives about $\eta^{2/3} = 5 \cdot 10^{-6}$); Ridders' estimate is **not** asserted to fail at $x = 1$ (there it covers, ratio 0.17): the example prints the estimate and the error and counts coverage over its handful of points for the reader at the keyboard, but the chapter text quotes only DS19's coverage macros over 300 points (about 91 % with a worst ratio of 42), never the example's own count |
| E18 | 7 | `07-using-check.cpp` | `check-by-steps`: $c_p$ at $T_0$ with `Central<2>` at $h$, $2h$ and $h/2$, and the agreeing digits counted; `check-by-order`: the same against `Central<4>`; `check-by-ridders`: Ridders as the independent check; `check-absolute`: a point where $f'$ is near 0 (the thread's $H''$ has no zero in range, so $\sin x$ near $\pi/2$), where the absolute error is the one to look at | The three stencil values agree to at least 8 digits with each other and with the closed form; the Ridders estimate covers the stencil's actual error (if not, the example prints why: the stencil's error is not what Ridders estimates); $c_p > 0$ |
| E19 | 7 | `07-using-kink.cpp` | `kink`: the published nitrogen enthalpy across the 500 K range break (`book::nitrogen_H_joined`, DS20) at $\delta = 1$ mK inside the stencil's span: `Central<2>` and `Central<4>` return finite values with a large error, Ridders returns a large estimate; the same at $\delta = 50$ K outside the span, where all agree | Both stencils return a value whose relative error exceeds $10^{-3}$ (measured 4.3 and 0.041); Ridders' estimate exceeds $10^{-3}$ of its value (measured $6.2 \cdot 10^{-3}$, flagged); outside the span every error is at most $10^{-8}$ (measured below $10^{-11}$) |
| E20 | 7 | `07-using-solver.cpp` | `implicit-function`: the Colebrook friction factor $f(Re)$ solved by bisection on $1/\sqrt{f}$ to a relative tolerance $\eta$ (no formula for $f$, so none for $\dd f / \dd Re$; the output is quantised at $\eta$); `implicit-derivative`: `derivative<Ridders>` and `Central<2>` with the default step and with $h = \eta^{1/3} Re$; the answer key by implicit differentiation of the Colebrook equation | Relative error of the fitted-step stencil at most $10^{-4}$ for $\eta = 10^{-10}$ (measured $7 \cdot 10^{-8}$; the model gives about $\eta^{2/3} = 2 \cdot 10^{-7}$); the default-step stencil's error is larger than the fitted step's (measured $5.5 \cdot 10^{-5}$); Ridders' estimate is reported, not asserted (measured: it covers, ratio 0.066). One sentence may contrast a fixed-point solver, whose output is smooth between count changes (`\valSolverFixedPoint...`). **Round 1:** the solver moves to the shared header `code/colebrook.hpp` (`book::friction`), quoted with `\snippetheader{colebrook}{implicit-function}`; the example deletes its own copy and keeps `implicit-derivative` |
| E1b | 1 | `01-why-program.cpp` | `program`: the complete first program of 1.4 (R1-decision T4), quoted whole between markers on the file's first and last lines: `#include "shomate.hpp"`, `#include <numerixx/deriv.hpp>`, `#include <print>`, `main`, the lambda `H` for `book::shomate_H(book::nitrogen, T)`, `T0 = 1000.0`, the call `nxx::deriv::derivative(H, T0)`, `if (d) ... else ...` printing `d->value` or the failure, and `return d ? 0 : 1;`. Every line within 60 characters; the text says once that `book::` is the namespace of the book's example header | Exit status 0 when a value was returned; the printed value is `\valCpDefaultValueLong`'s |
| E21 | 3 | `03-rounding-noise-estimate.cpp` | `noise-estimate-use`: the recipe `book::estimate_noise` (quoted from `code/noise.hpp` with `\snippetheader{noise}{noise-estimate}`) applied to the friction factor of E20 (`book::friction(re, 1e-10)` from `code/colebrook.hpp`) at $Re = 10^5$ and to the noisy $\mathrm{e}^x$ of E17 at 1, with $h = $ `root_epsilon<3, double>() * max(|x|, 1)` and the default eight points; prints the estimate relative to $f$, the difference order used and the cost (eight evaluations). Chapter 3.6 quotes DS22 (`\valNoiseEstimate...`); 7.4.3 quotes `\valNoiseEstimateFriction` | The friction estimate within a factor 10 of $\eta = 10^{-10}$ and the noisy-$\mathrm{e}^x$ estimate within a factor 3 of $\eta / \sqrt{3}$ (a uniform noise of half-width $\eta$ has that root mean square); the order used at least 2 |
| E22 | 7 | `07-using-jacobian.cpp` | (Stage-6 close-out: the book quotes the loop of `jacobian(r, x)` only, de-indented between the markers; the function's header and `return J` stay in the file and the prose names them. `jacobian-by-calls` is named in the prose, not printed.) `jacobian-balances`: a $2 \times 2$ Jacobian of two residuals in the unknowns $T$ (K) and $n$ (mol), $r_1 = n H(T) - Q$ (an energy balance on the thread's enthalpy) and $r_2 = n R T / V - P$ (an ideal-gas pressure), the model as one function returning both residuals, at $T = 1000$ K, $n = 2$ mol, $V = 1$ m$^3$, $R = 8.31446261815324$, with $Q$ and $P$ chosen so that both residuals vanish there; the central stencil applied by hand to the vector residual, column $j$ from the model at $x_j \pm h_j$ with $h_j = $ `root_epsilon<3, double>() * max(|x_j|, 1)` (two model calls per column, four in all); `jacobian-by-calls`: the same by one `derivative` call per entry on a lambda that evaluates the whole model (eight model calls). Closed-form partials as the answer key: $n c_p(T)$, $H(T)$, $n R / V$, $R T / V$. 7.4.2's loop over the parameters of the Arrhenius rate is renamed "parameter sensitivities" and printed (`jacobian` of E13, kept) | Every entry within $10^{-8}$ of its closed form by both ways (the generator measures about $10^{-11}$ and $10^{-10}$, `\valJacobianBalance...`); the call counts 4 and 8 |
| E23 | 7 | `07-using-usertype.cpp` | `usertype-limits`: a strong typedef `Real` over `double` (members `+ - * /`, unary minus, `<=>` and `==`, construction from `int`), with `std::numeric_limits<Real>` specialised to exactly the members the module uses, `is_specialized`, `is_integer`, `is_bounded`, `is_exact`, `radix`, `digits`, `epsilon()`, `max()`, `denorm_min()`, and not inheriting `numeric_limits<double>`; `usertype-call`: `derivative(f, Real(1))` with `Central<2>` and Ridders on $\mathrm{e}^x$. Quoted in Sidebar 7.1 and 8.7 (R1-developer-01) | `static_assert(nxx::real<Real>)`; both calls return a value within $10^{-9}$ of $\mathrm{e}$ |
| CMake | 7 | `getting-started.cmake` (editor-owned; quoted with `\cmakefile{getting-started}`) | The three steps of 7.1 "Getting the library": `FetchContent_Declare` at a pinned commit, `FetchContent_MakeAvailable`, `target_link_libraries(app PRIVATE numerixx::numerixx)`; the one listing the build does not configure (section 8, item 23) | none (not built) |

E14 needs the `boost_multiprecision` target that `tests/CMakeLists.txt` creates (decision A2). E19's function depends
on decision A1 (section 1.3).

---

## 6. References to verify

Every entry enters `refs.bib` with `note = {UNVERIFIED}` and is cited only after the evidence reviewer has checked it
against the work or a library catalogue (authors, title, year, publisher or journal, volume, pages, DOI). Accented
names are written with LaTeX accents in the `.bib`.

| Key | Work | Cited in | For |
|---|---|---|---|
| `burden2016` | R. L. Burden, J. D. Faires, A. M. Burden, *Numerical Analysis*, 10th ed., Cengage, 2016 | 2, 4, 6, B | numerical differentiation from Lagrange interpolation; Richardson extrapolation and its tableau; the a posteriori estimate |
| `chapra2015` | S. C. Chapra, R. P. Canale, *Numerical Methods for Engineers*, 7th ed., McGraw-Hill, 2015, ch. 23 | 4, 5, 9 | the finite-difference formulas up to the fourth derivative, the O(h), O(h^2), O(h^4) labels |
| `ridders1982` | C. J. F. Ridders, Accurate computation of F'(x) and F'(x)F''(x), *Advances in Engineering Software* 4(2), 75-76, 1982 | 6, 10 | the adaptive extrapolation with an error estimate |
| `fornberg1988` | B. Fornberg, Generation of finite difference formulas on arbitrarily spaced grids, *Mathematics of Computation* 51(184), 699-706, 1988 | 4, 9 | published weight tables (the tests' reference) |
| `richardson1911` | L. F. Richardson, The approximate arithmetical solution by finite differences of physical problems involving differential equations, with an application to the stresses in a masonry dam, *Phil. Trans. R. Soc. A* 210, 307-357, 1911 | 6 | the origin of the extrapolation |
| `richardson1927` | L. F. Richardson, J. A. Gaunt, The deferred approach to the limit, *Phil. Trans. R. Soc. A* 226, 299-361, 1927 | 6 | the deferred approach to the limit |
| `higham2002` | N. J. Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed., SIAM, 2002 | 3, A | the rounding model, cancellation, Sterbenz's lemma |
| `goldberg1991` | D. Goldberg, What every computer scientist should know about floating-point arithmetic, *ACM Computing Surveys* 23(1), 5-48, 1991 | 3, A | floating-point basics for the reader |
| `sterbenz1974` | P. H. Sterbenz, *Floating-Point Computation*, Prentice-Hall, 1974 | 3 (S4), 5 (S9) | Sterbenz's lemma (exact subtraction) |
| `sidi2003` | A. Sidi, *Practical Extrapolation Methods: Theory and Applications*, Cambridge University Press, 2003 | 6, B | Richardson extrapolation theory and error estimation |
| `ieee754-2019` | IEEE Std 754-2019, *IEEE Standard for Floating-Point Arithmetic*, 2019 | A | the formats, rounding to nearest, subnormals |
| `lyness1967` | J. N. Lyness, C. B. Moler, Numerical differentiation of analytic functions, *SIAM J. Numer. Anal.* 4(2), 202-210, 1967 | 1 (S1) | the complex step and contour methods |
| `squire1998` | W. Squire, G. Trapp, Using complex variables to estimate derivatives of real functions, *SIAM Review* 40(1), 110-112, 1998 | 1 (S1) | the complex-step derivative |
| `griewank2008` | A. Griewank, A. Walther, *Evaluating Derivatives: Principles and Techniques of Algorithmic Differentiation*, 2nd ed., SIAM, 2008 | 1 (S1) | automatic differentiation |
| `more2012` | J. J. More, S. M. Wild, Estimating derivatives of noisy simulations, *ACM Trans. Math. Software* 38(3), 2012 | 3 (3.6), 7, 10 | derivatives of functions with computational noise; the noise-aware step named as deferred |
| `dennis1996` | J. E. Dennis, R. B. Schnabel, *Numerical Methods for Unconstrained Optimization and Nonlinear Equations*, SIAM Classics in Applied Mathematics 16, 1996 (first published 1983), s. 5.4 and 5.6 | 3 (3.4), B | the finite-difference step scaled by $\max(|x|, \text{typical } x)$, the origin of the default step's scale |
| `nocedal2006` | J. Nocedal, S. J. Wright, *Numerical Optimization*, 2nd ed., Springer, 2006, s. 8.1 | 3, B | the $\sqrt{\varepsilon}$ and $\varepsilon^{1/3}$ analysis of the forward and central difference |
| `gill1983` | P. E. Gill, W. Murray, M. A. Saunders, M. H. Wright, Computing forward-difference intervals for numerical optimization, *SIAM J. Sci. Stat. Comput.* 4(2), 310-321, 1983 | 7 (7.9), 10 | adaptive step selection, named as a deferred feature (optional) |
| `chase1998` | M. W. Chase, Jr., *NIST-JANAF Thermochemical Tables*, 4th ed., *J. Phys. Chem. Ref. Data* Monograph 9, 1998, two parts, pp. 1-1951, ACS and AIP for NIST, Woodbury NY, ISBN 1-56396-831-2 (nitrogen: Part II, p. 1621) | 1, 7 | the tabulated thermochemical data of nitrogen that the Shomate fit represents. **Verified** by the evidence reviewer (2026-10-23); it contains no Shomate coefficients, so it is never cited for them |
| `nist-webbook` | P. J. Linstrom, W. G. Mallard (eds.), *NIST Chemistry WebBook*, NIST Standard Reference Database Number 69, NIST, Gaithersburg MD, DOI 10.18434/T4D303, nitrogen gas-phase thermochemistry page, retrieved 23 October 2026 | 1, 7 | the Shomate equations and the coefficients of the thread (the WebBook's own fit of January 2009). **Verified** (2026-10-23); cited in the WebBook's requested form with the retrieval date |
| `shomate1954` | C. H. Shomate, A method for evaluating and correlating thermodynamic data, *J. Phys. Chem.* 58(4), 368-372, 1954, DOI 10.1021/j150514a018 | 1 (optional) | the origin of the name, only if the text names it; metadata verified in Crossref, content not seen |
| `iso14882-2024` | ISO/IEC 14882:2024, *Programming languages -- C++* | 7, 8 | `std::expected`, `[[deprecated]]`, constant evaluation (optional; only if the practice author cites the standard) |

The design record (`docs/PHASE3_API.md`, `docs/PHASE4_API.md`, `docs/design.md`) is not literature: the book refers to
it by file name in Chapter 10 as "the design record in the repository", not in `refs.bib`. Numerical Recipes is not
cited (nothing in the library comes from it, and the book has no reason to compare).

---

## 7. Ownership

| Owner | Agent | Files |
|---|---|---|
| Editor | `book-editor` | `PLAN.md`, `main.tex`, `preamble.tex` (design changes only with the author's approval), `macros.tex`, `refs.bib` (entries proposed by authors, verified by the evidence reviewer), `chapters/00-frontmatter.tex`, `chapters/99-backmatter.tex`, `chapters/shared-method-summary.tex` (the five rows of T0, input in 1.5 and on the quick-reference page), `tools/generate_data.cpp` and `data/` (every data set of section 4.3), `data/README.md`, `figures/README.md`, `code/README.md`, `CMakeLists.txt` of the book, the temporary `chapters/zz-figure-gallery.tex` (stage 2) and the removal of it, of `chapters/zz-layout-test.tex`, of `code/quickstart.cpp` and of the template data set (stage 4, done), `reviews/R<n>-decision.md` |
| Fundamentals author | `book-author-fundamentals` (Opus, high) | `chapters/01-why.tex`, `02-taylor.tex`, `03-rounding.tex`, `a1-floating-point.tex`; figures `figures/01-why-secant.tex` (F1), `02-taylor-fd-slopes.tex` (F2), `03-rounding-v-curve.tex` (F3), `03-rounding-step-scale.tex` (F4), `03-rounding-step-scale-small.tex` (F4b), `03-rounding-float-double.tex` (F5), `03-rounding-noise.tex` (F17); the table T0 (with the practice author, who repeats it on the quick-reference page); code E1 to E4, E17; sidebars S1 to S5 |
| Methods author | `book-author-methods` (Fable, high) | `chapters/04-stencils.tex`, `05-orders.tex`, `06-ridders.tex`, `a2-derivations.tex`; figures F6 to F15 (`figures/04-stencils-*.tex`, `05-orders-*.tex`, `06-ridders-*.tex`); code E5 to E10; sidebars S6 to S12 |
| Practice author | `book-author-practice` (Opus, high) | `chapters/07-using.tex`, `08-api.tex`, `09-tables.tex`, `10-guarantees.tex`, `qr-quick-reference.tex`; figures F16 (`figures/07-using-newton.tex`), F18 (`figures/07-using-kink.tex`); the tables T1 to T3 (typed, from the header) and the `\input` of T4 to T16; code E11 to E16, E18 to E20; sidebars S13, S14 |
| Reviewers | `book-reader-engineer`, `book-reader-developer`, `book-reviewer-numerics`, `book-reviewer-library`, `book-reviewer-evidence`, Codex | `reviews/R<n>-<reviewer>.md` |
| Authors (revision) | as above | `responses/R<n>-<chapter-id>.md` |

Authors request data sets and macros from the editor by name (sections 2.3 and 4.3 already list them); they do not
edit `generate_data.cpp`, `data/` or `macros.tex`.

---

## 8. Style decisions beyond CLAUDE.md

1. **Part I refers to Part II for tables.** Part I shows at most one small margin table per chapter (the two-point
   rules in Chapter 2, the node counts in Chapter 4) and refers to Chapter 9 for every complete table with `\cref`.
   No table appears twice, with one deliberate exception: the five-row method choice T0 appears in 1.5 ("If you read
   nothing else") and again on the quick-reference page, which repeats by design; the full decision table T1 lives in
   Chapter 7 only, and T0 refers to it with `\cref{tab:method-choice}`.
2. **The decision table (T1).** A `table*` (full width) in Chapter 7: rows are situations ("inside an iteration, fn
   smooth and cheap", "a value to report, with an error estimate", "fn defined on one side only", "fn varies on a
   scale much smaller than |x|", "a second derivative", "fn is the output of an iterative solver (noisy)", "a wide
   type, 50 digits"); columns: method and order to use, the step (default or user), evaluations, accuracy to expect
   (as $\varepsilon^{P/(P+N)}$ and a generated number for double), and where the evidence is (a `\cref` to the figure
   or the test). Every accuracy number in it is a `\val` macro.
3. **Numbers.** `\val` macros hold raw numbers (`1.5e-11`); the text typesets them with `\num{\valX}` or
   `\SI{\valX}{unit}`, never raw. Errors, estimates and steps to two significant digits; coverage to one decimal in
   percent; evaluation counts as integers. The generator formats accordingly. Exact constants (2^-52, 1/6, 4^j - 1)
   are typed in mathematics. Relative errors are dimensionless and are the default in plots and tables; absolute
   values carry units.
4. **Terminology.** "step" (not step size, increment) $h$; "rule" for a formula with its nodes and weights, "stencil"
   for the family (`Central<P>`) and for a rule's node pattern; "node" for an offset $n_k$, "sample" for a value
   $f_k$; "evaluation" for one call of fn; "estimate" for a returned value, "error estimate" for `abs_error`;
   "truncation error" and "rounding error" (never round-off); "machine epsilon" $\varepsilon$; "the default step";
   "the scale $s$"; "accuracy order P" and "derivative order N", always qualified in a sentence that has both;
   "coverage" for the fraction of cases whose estimate is at least the actual error; "level", "column", "entry",
   "tableau", "restart" for Ridders; "non-finite" in prose, `\code{non_finite}` for the enumerator; "the central
   difference" (lower case) for the mathematical rule, `\code{Central<2>}` for the library's method; "noise level"
   $\eta$ for the relative accuracy of $f$ when it is worse than $\varepsilon$; "kink" for a point where $f$ or a
   derivative of $f$ jumps (a range break, a phase boundary, an `if`). Every rule's library tag is named the first time
   the rule appears (Chapter 2), with the margin note "P is the accuracy order, not the number of points: `Central<2>`
   uses two points for $f'$, `Forward<1>` two".
5. **Names.** "Ridders' method" and "Ridders' extrapolation" (apostrophe, no s); "Richardson extrapolation";
   "Taylor's theorem"; "Sterbenz's lemma"; "Newton's method"; "Numerical Recipes" only if the author decides to
   mention it (not planned).
6. **Prototype statistics.** The book does not state the Python prototype's coverage figures as facts about the
   library. Its own coverage study (DS9) is the evidence; the prototype is mentioned once, in S11 and Chapter 10, as
   the design record that chose the constants, by file name.
7. **Measured, proven, rule of thumb.** Each claim is tagged by its phrasing: "is" with an equation reference for a
   proof, "measured ... (\cref{fig:...})" for a generated result, "as a rule" for a heuristic (the default step, the
   edge distance). Derived limits that the tests do not cover (the edge distance) are presented as measured (DS10)
   with the derivation as the explanation.
8. **Code in prose.** `\code{derivative<Central<4>, 2>(fn, x)}`; template arguments without spaces around `<`;
   `fn`, `x`, `h` as the parameter names; results named `d` for a stencil and `r` for Ridders in the examples.
9. **Citations.** Superscript numbers after the punctuation mark; one citation per established method or result at
   its first mention in a chapter; no citation in captions.
10. **Equations.** Numbered only when referred to; `\cref` for every reference; the error-model equation, the two
    differences of Chapter 2, the moment conditions, the tableau formula and the default step are numbered.
11. **Sidebars** are referred to in parentheses at the end of the sentence that motivates them: "(\cref{sb:ulp})".
12. **British spelling list** (beyond CLAUDE.md): behaviour, centre, modelling, analyse, organise, generalise,
    minimise, optimise, initialise, favour, programme only for a schedule (a computer program is "program"), zeros.
13. **Version statement.** The front matter names the commit of `dev` the book describes; the editor fills it at
    release.
14. **The thread's units.** $H$ in J/mol, $T$ in K, $c_p$ in J/(mol K); the Shomate coefficients appear once, in the
    E1 snippet, with the source in a comment and the citation in the text.
15. **Guarantees carry their hypotheses.** A guarantee is stated with the condition it rests on and with how it was
    verified: "compile-time and run-time results agree, provided the compiler neither contracts `a * b + c` into a
    fused multiply-add nor uses extended precision; verified on the five CI toolchains" (decision A5). "Every
    compiler" is never written.
16. **The two lessons of the edge.** "Evaluable is not valid" (a fit sampled outside its range, Chapter 4) and
    "evaluable is not smooth" (a kink inside the stencil's span, Chapter 7) are worded exactly so, once each, and
    the quick-reference page repeats both.
17. **`std::expected` idioms (Chapter 7, from stage 2).** Never suggest `value_or` on the whole result with a
    `double` (the value type is `Estimate<T>`, which no `double` converts to); only `d->abs_error.value_or(...)`.
    `if (*r)` is a `-Wfloat-conversion` error under Clang with `-Werror`, although the header documents it as
    compiling: the examples write `static_cast<bool>(*r)` or compare explicitly, and the text names `if (*r)` only
    as the pitfall, with that caveat (library issue U-2 in section 12.3).
18. **Snippet code.** Every multi-line `if` in a snippet is braced (Clang's `-Wmisleading-indentation` fires on
    column-0 code); a single-line `if (d) std::println(...);` may stay unbraced (round 1, R1-library-16: it builds
    clean on both compilers and reads better). Snippet names may be hyphenated (`\snippet` braces its range ends);
    `tools/check.py` rejects a `\snippet` whose markers are missing or enclose nothing. `return 1;` from `main` stays
    in the examples: they are tests, as 1.4 and the front matter say.
19. **Figures.** `restrict x to domain` (and `y`) on a logarithmic pgfplots axis takes natural-logarithm bounds
    ($\ln 10^{-2} = -4.61$). A margin float (`marginfigure`, `margintable`) must not be placed beside a full-width
    float (`figure*`, `table*`): keep them a page apart (F14b and F15 against F12 and F14 in Chapter 6).
20. **Table bodies.** A generated body is read with `\tablebody{<name>}` on a line of its own between `\midrule`
    and `\bottomrule` (`macros.tex`): LaTeX's `\input` runs file hooks whose tokens open a new row and give
    "Misplaced \noalign".
21. **Listings and the asterisk.** Plex Mono sets `*` raised, so `a * b` in a listing reads as a superscript to
    some eyes. Possible fix `literate={*}{{\textasteriskcentered}}1` in the `bkcode` style: a preamble change, for
    the author to decide (listed in 12.2, S2-19); until then accepted.
22. **The transcribed declarations (round 1, T14).** The `signature` blocks of Chapter 8 are the one typed code in
    the book: the declarations are transcribed from `include/numerixx/deriv.hpp` and `math.hpp` (quoting them would
    need markers in a public header, which the book may not add), and `tools/check.py` fails the build when a
    `\code{...}` line of a block no longer occurs in the headers, comments stripped and whitespace removed. The
    front matter says so.
23. **The CMake file (round 1, T4).** `code/getting-started.cmake`, quoted whole in 7.1 with `\cmakefile`, is the
    one listing the book's build does not configure or run; `check.py` only checks that the file exists. The front
    matter says so. Everything else quoted in the book is a compiled example or a generated table.
24. **"Digits" (round 1, T13).** "About $n$ digits" is the rounded value of $-\log_{10}$ of the error, never its
    floor; the headline macros are rounded by the generator, and the model's digits are quoted to one decimal from
    `\valModelDigits...`. Counts and levels stated as facts are macros, not words.
25. **The two checks (round 1, T1).** Part I reads a result as `d->value` and the estimate as `*r->abs_error`
    (glossed at its first use); the implicit conversion of `Estimate<T>` appears only in 7.2 and Chapter 8. Section
    1.5 and the quick reference teach two checks, "did it work" (`if (d)`) and "is it right" (the silent failures of
    3.5, 3.6 and 7.8, the repeat at $2h$ and $h/2$, and for Ridders a tolerance on `abs_error`).

---

## 9. Schedule of stages 2 to 8

Each stage ends at a checkpoint with the author; commits to `dev` follow an approved stage (CLAUDE.md).

| Stage | Who | Produces | The author approves |
|---|---|---|---|
| 2 Evidence | editor; authors | Editor: the `macros.tex` additions of 2.3; the data sets DS1 to DS21 in `generate_data.cpp` (after the evidence reviewer has verified the Shomate coefficients, A1); `refs.bib` entries marked UNVERIFIED; the CMake change for Boost.Multiprecision (A2); `data/README.md`. Authors: the figure files F1 to F18 against the CSVs; the examples E1 to E20 passing `ctest -R numerixx_book` on the gcc preset. Editor: a temporary `chapters/zz-figure-gallery.tex` that shows every figure with its planned caption, built to a PDF | The figures (gallery PDF), the examples' output, the thread's data source |
| 3 Drafting | authors, in parallel | Each chapter and appendix as a file that passes `tools/check_chapter.sh`, within +-1 page of its target (page count from the wrapped PDF); new macros marked `%% NEW-MACRO`; `refs.bib` requests to the editor | Optional: a reading snapshot of one chapter per author |
| 4 Integration | editor | `main.tex` inputs; merged macros; one voice and terminology (section 8); every `\cref` resolves; the front matter (preface, how to read, the AI statement as the process actually ran: models, roles, rounds) and the back matter (notation, glossary); the layout test and the gallery removed; the full build clean, `tools/check.py` with no errors and every warning cleared or justified; total pages within 74 to 88 | The integrated PDF and the AI statement draft |
| 5 Review round 1 | the five reviewers, in parallel | `reviews/R1-engineer.md`, `R1-developer.md`, `R1-numerics.md`, `R1-library.md`, `R1-evidence.md` (every reference verified or rejected); editor: `reviews/R1-decision.md` | The decision letter's summary; anything contested or scope-changing |
| 6 Revision | authors | Chapters revised per the decision letter; `responses/R1-<chapter-id>.md`; regenerated data if the library changed (none expected) | -- (the editor checks the responses) |
| 7 Review round 2 | targeted reviewers; Codex; Cline with DeepSeek V4.1 Flash 1017 (R1-decision, approved; the model per the author, 2026-10-24, S6-11) | `reviews/R2-*.md` on the revised chapters; `reviews/R2-codex.md` from `codex exec -s read-only --ephemeral --color never -o <scratchpad>/codex_review.md "<prompt>"` naming the PDF, the sources, CLAUDE.md and PLAN.md and asking for ranked findings with evidence, and `reviews/R2-cline.md` from the Cline review with the same brief; `R2-decision.md`; final fixes | The final PDF |
| 8 Release | editor; author | The AI statement made current; the version statement filled; `data/` regenerated and its diff reviewed; commit and push to `dev` | The release |

---

## 10. The open questions of BRIEF.md

### 10.1 Resolved here

| Question | Decision | Reason |
|---|---|---|
| The chapter split and page targets | Section 1: ten chapters as in the brief, two appendices and a quick-reference page, 82 nominal pages | The brief's draft is sound; the targets were set from the figure, sidebar and listing counts so that each chapter keeps at least 1000 words of running text (Chapter 9 excepted, which is tables); after the stage-1 reviews Chapter 6 gave two pages to Chapters 3 and 7 and Chapter 1 gained its summary |
| Part II repeats Part I's tables, or Part I refers to Part II | Part I refers to Part II (style decision 1) | One source of truth for every weight and constant, generated once; the length budget |
| Which examples run through several chapters | The heat capacity from an enthalpy (section 1.3) in every Part I chapter; Newton's method and the sensitivity are single-chapter examples in Chapter 7 | Closed forms for every derivative, a large $\lvert x \rvert$, a natural edge and a natural second derivative |
| The quick-reference page at the front or the back | The book's last page, after the references, with a pointer in "How to read"; its five-row method choice also closes Chapter 1 (1.5) | A reference user opens the back cover; between Appendix B and the notation it was buried (engineer reader, finding 9); the Chapter 1 summary gives the first-time reader the method choice before the theory; moving the page to the front is a one-line change in `main.tex` if the author prefers |
| Sidebar count | Fourteen (section 3) | At most three per chapter, about 9 pages |
| Figures with units | Relative errors on the axes; absolute values with units in the text | Comparability across methods and orders |

### 10.2 Decisions needed from the author

Approved by the author on 2026-10-23, each with the recommended option: A1 NIST Shomate coefficients for nitrogen,
verified before use; A2 the guarded `boost_multiprecision` link in the book's build; A3 the generator reads
`detail::make_table`; A4 82 pages, band 74 to 88; A5 option (i) for the book, option (ii) raised as a separate
library issue.

- **A1. The thread's data.** Recommend: the Shomate coefficients of nitrogen from the NIST Chemistry WebBook (Chase
  1998), verified by the evidence reviewer at stage 2 before any number is generated. Alternative: an illustrative
  polynomial labelled as such, which needs no external data but describes no real substance. The editor fixes the
  range and $T_0$ once the coefficients are verified (provisionally the range containing 1000 K, $T_0 = 1000$ K).
- **A2. A build change in `docs/deriv-book/CMakeLists.txt`.** Recommend: link `boost_multiprecision` (the target
  `tests/CMakeLists.txt` creates; `tests/` is added first, so it exists when the default options are on) to the
  generator and to `07-using-multiprecision.cpp`, guarded by `if (TARGET boost_multiprecision)`; without it, E14 and
  the 50-digit column of T14 are dropped. This touches only the book's build, not the library.
- **A3. The generator reads `detail::make_table`.** Recommend: yes, as the tests do, so that the weight tables of
  Chapter 9 are generated and self-checked rather than typed. The book still describes the tables as internal.
- **A4. Length band.** The plan targets 82 pages with an acceptable 74 to 88 at integration (the brief says 60 to 90).
  Approve the band, or name a tighter target.
- **A5. The compile-time identity claim (numerics finding 2).** The header and `docs/design.md` say the helpers "give
  the same result on every compiler" and the tutorial that compile time and run time agree. The code forms
  `x + T(node) * h`, `c * fk` and `sum + term` (deriv.hpp lines 962 and 970-972), which a compiler may contract into
  fused multiply-adds at run time (GCC's default outside strict ISO mode, Clang's within an expression, on FMA
  hardware), while constant evaluation never contracts; the build pins neither contraction nor extended precision,
  and the tests show the identity only on the five CI toolchains. Two options:
  (i) the book states the guarantee with its hypothesis and its verification ("compile-time and run-time results
  agree provided the compiler neither contracts nor uses extended precision; verified on the five CI toolchains"),
  the library stays as it is, and the design document's "on every compiler" is corrected as a documentation fix;
  (ii) the library pins contraction (for example `#pragma STDC FP_CONTRACT OFF` or the compilers' equivalents in the
  headers, or `-ffp-contract=off` on the interface target), after which the book may state the identity without the
  hypothesis. Option (ii) changes `include/` and the build, so it follows the repository's API-proposal process and
  is outside the book. Recommendation: (i) for the book now, since the book describes the library as it is; raise
  (ii) as a library issue to be decided on its own merits (it costs a little run-time accuracy where FMA would have
  helped, and a contraction that differs between a user's two translation units already makes bit identity across
  them impossible to promise).

---

## 11. Discrepancies between the sources

Recorded per CLAUDE.md ("where they disagree, the higher one wins; report the disagreement"). Code beats tests beats
design documents beats tutorial, survey and README.

| Id | Sources | Disagreement | Resolution for the book |
|---|---|---|---|
| D-1 | `docs/DERIV_TUTORIAL.md` s. 4 against `math.hpp` | The tutorial's `root_epsilon` forms `radix^N` up front (`radix_n`) and uses the textbook Newton update `((N-1) y + a / y^(N-1)) / N`; the code multiplies by the radix one factor at a time (so N = 128 in float no longer overflows, a tested regression) and uses the correction form `y + (a / y^(N-1) - y) / N`, within 0.58 ulp where the textbook form stopped up to 1.8 ulp above the root. The tutorial also says the tests compare with `std::cbrt` and `pow`; they compare with 80-digit references bracketed by neighbours, `cbrt` only for `long double` | Appendix A and S3 describe the code's version; nothing is quoted from the tutorial |
| D-2 | `docs/DERIV_TECHNIQUES.md` s. 2.2 against `deriv.hpp` | Ridders' cost "typically 6 to 12 evaluations" (written before implementation) against "typically 10 to 16 for f' in double" | The book quotes its own counts (DS9) |
| D-3 | `README.md` against the code | "finite differences of orders up to 10" holds for `Central` only (one-sided rules stop at 6); the numbers in its comments ("about 1e-13" for `Central<4>` on sin) are unverified | The book states both caps and uses generated numbers |
| D-4 | `docs/PHASE2_API.md` against Phase 4 and the code | Superseded statements: `Central<P>` has "P evaluations" (now P + N - 1), weights from "closed forms" (now one Lagrange generator), largest numerator 2100 (now 6,222,216 for the fourth derivative), the dropped constant costs "at most 1.64" (now up to 15 for `Central<4>` and f''''), the return type `std::expected<T, error>` (now `Estimate<T>`). `PHASE3_API.md` likewise predates N > 1 | Authors cite the header, `design.md` and `PHASE4_API.md` for current behaviour; PHASE2 and PHASE3 only for the history of a decision |
| D-5 | `BRIEF.md` against the code and tests | "A compile-time warning for a 4th derivative in float": the criterion is any type with epsilon above 1e-7 (`radix^(digits-1) < 10^7`: float, decimal32, bfloat16). "Works for any nxx::real T (float, double, long double, ...)": the derivative tests run on float, double, `cpp_bin_float_50` with and without expression templates and `cpp_dec_float_50`; `long double` is checked only against the concept and in `root_epsilon<3>` | Chapter 5 states the criterion; Chapter 10 says "tested on ..." precisely and "long double satisfies the concept" |
| D-6 | `PHASE3_API.md` and `PHASE4_API.md`; neither in code or tests | Not a disagreement (reclassified after the numerics review): PHASE3 covered N = 1 only (2^-34 from $h_0 = s/32$ and the 32-level cap), PHASE4 keeps it and adds 2^-31 for N >= 2 because $h_0 = s/4$; the same derivation. The distance is derived, not tested | Measured by DS10 and stated as measured (F15); the derivation explains it |
| D-7 | `docs/design.md` "Tests" against `tests/test_stencils.cpp` | `design.md` calls `(C + 2S) eps^(p/(p+1))` "a bound that the measured errors stay 4 to 16 times below"; the test comment says "an estimate, not a bound" | The book calls it an estimate |
| D-8 | `DERIV_TUTORIAL.md` s. 9, `PHASE3_API.md`, `tests/test_ridders.cpp` | The random study's size: "2,700 random cases" (per type) against "5,400" (double and float together); the worst ratio 3.6 is double's, 3.1 float's | Consistent once per type; the book does not quote these as library facts (style decision 6) |
| D-9 | `deriv.hpp` (Ridders documentation) against `PHASE4_API.md` | "13 to 23" typical evaluations for f'' against a prototype mean of 17 and maximum 23; "25 to 45" for f''' and f'''' against maxima of 40 and 41 | No conflict; the book quotes DS9 |
| D-10 | `DERIV_TUTORIAL.md` s. 5 against `deriv.hpp` | The tutorial's `derivative_value` sketch inlines the input check and returns `std::expected<T, error>`; the code's structure is `derivative_with_step` / `stencil_derivative` / `difference`. Labelled a simplification | The book quotes no code from the tutorial |
| D-11 | `demo/deriv_demo.cpp` against `code/quickstart.cpp` | `std::println` against `std::printf` | `std::println` in every example (section 5); the CI builds the demo with it on all five toolchains |
| D-12 | `BRIEF.md` "state at commit ae17769" against `HEAD` | None: the library, tests and design documents are unchanged since `ae17769` | -- |
| D-13 | `math.hpp` file comment and `docs/design.md` ("the same result on every compiler"), `DERIV_TUTORIAL.md` s. 4 and 6, against the code as built | The identity of compile-time and run-time results is not pinned by the code or the build: the sums and products of `deriv.hpp` may be contracted into fused multiply-adds at run time, and the build sets no contraction mode; the tests verify the identity on the five CI toolchains only (numerics finding 2) | Decision A5; until then the book words the guarantee with its hypothesis (style decision 15) |
| D-14 | `deriv.hpp` (the `Estimate` comment, "the operator templates of Boost.Multiprecision") and this plan's S13 ("`operator==` on Boost types ignores the conversion") against the behaviour (practice author, stage 3: GCC 16, Boost 1.92) | Comparisons of an `Estimate<T>` with a Boost number compile and are true; it is the arithmetic operator templates (`*r - e`) whose deduction does not consider the conversion. The header's phrase covers both and is broader than the behaviour | S13 and Chapters 7 and 8 state the arithmetic case only ("the arithmetic operator templates", "its comparisons do accept it"); the header's wording is a library documentation note, U-3 |
| D-15 | `deriv.hpp` (the `detail::noisy` comment: `Central<2>` "about 8e-2 for float, and Ridders a median of 2e-3 with a tenth of the cases more than 100 % off") against DS6 (0.014 and 0.0019, medians over eleven points of $\mathrm{e}^x$) | Different test sets: the header's figures are the design record's over its random cases, the book's are medians on $\mathrm{e}^x$; neither is wrong, and they are not comparable | The book quotes the DS6 values only, says what they were measured on, and does not quote the header's figures (style decision 6) |
| D-16 | `deriv.hpp` and `docs/design.md` ("costs at most a factor 1.64") against the model formula of PHASE4_API.md and Appendix B.4 evaluated from the tables (`\valStepCostCentralTenFirst` = 1.63 for `Central<10>`, $\kappa$ = 1.80) | 1.634 rounds to 1.63, not 1.64; the header's figure is presumably a rounding of $\kappa$ to 1.8 before the cost was formed. The other cost figures agree once rounded (1.12 for `Central<2>` and $f'$; 3.54, 11.2 and 15.1 against the record's 3.5, 11 and 15) | The book quotes the generated values; a library documentation note (U-4) |
| D-17 | `deriv.hpp` (the Ridders comment names "aliasing" without an example; `tests/test_ridders.cpp` tests $\sin(64 \pi x)$ at 0) against the stage-3 draft of 7.8 ("$\sin x$ near $x = 10^4 \pi$ ... is the documented case") | No source documents the $10^4 \pi$ case | 7.8 names the tested case and gives $\sin x$ at a large $\lvert x \rvert$ as the everyday form, without calling it documented |

No disagreement was found between `deriv.hpp` or `math.hpp` and the tests.

---

## 12. Stage-1 review decisions

Reviews: the engineer reader (A1 to A10) and the numerics reviewer (B1 to B12). Each finding: accepted, accepted in a
modified form, or rejected with the reason; the plan sections changed.

| Finding | Decision | Where |
|---|---|---|
| A1 noisy and tabulated functions | Accepted: section 3.6 with the $\eta$-for-$\varepsilon$ rule, F17, DS19, E17; tables out of scope in 7.9; a row in T1 | 1.4 (ch. 3, 7), 2.3, 4, 5 |
| A2 checking a stencil result | Accepted: section 7.5 with E18; three lines on the quick-reference page | 1.4 (ch. 7), 5 |
| A3 small-$\lvert x \rvert$ default step | Accepted: 3.5 covers both cases, F4b, DS4b, a second snippet in E4, the checklist and a T1 row | 1.4 (ch. 3), 4, 5, 8 |
| A4 non-smooth functions | Accepted: pitfall 7.8 with E19, F18, DS20; "evaluable is not smooth" | 1.4 (ch. 7), 4, 5, 8 |
| A5 Chapter 3 density and unearned notation | Accepted: Chapter 3 keeps the concrete central-difference model; $\ampl$ and $\errc$ move to 4.3, N to Chapter 5; float to a paragraph with F5 in the margin; Chapter 3 gains a page from Chapter 6 | 1.4 (ch. 3, 4) |
| A6 Chapter 6 internals; 4.5 and 5.4 | Accepted in a modified form: Chapter 6 cut to 7 pages (not 6: its four figures, the tableau table and three sidebars already take about three pages), the stopping rules, factors and caps moved to S11 and Chapter 10; the caps of 4.5 moved to 10.4; the actual offsets are S9 only | 1.4 (ch. 4, 5, 6, 10), 3 |
| A7 method choice early | Accepted: 1.5 "If you read nothing else" with T0, headline accuracy and cost as `\val` macros, the minimal check; Chapter 1 gains a page | 1.4 (ch. 1), 8 |
| A8 thread motivation; missing examples | Accepted in a modified form: 1.4 states that the closed form is the answer key; E20 is an implicit friction factor with an inner iteration rather than a cubic equation of state (the same lessons, shorter code within 60 characters a line, an answer key by implicit differentiation, no further thermodynamic data to verify); the Jacobian loop and the scaled sensitivity join E13 | 1.3, 1.4 (ch. 7), 5 |
| A9 quick reference buried | Accepted: the book's last page, after the references; pointer kept in "How to read" | 1.1, 1.2, 1.4, 10.1 |
| A10 two meanings of order; tag names | Accepted: tags named at first appearance in Chapter 2 with the margin note; Chapter 2's goal softened | 1.4 (ch. 2), 8 |
| B1 noise and piecewise evidence | Accepted (with A1 and A4): DS19, F17, E17; DS20, F18, E19; the rule in 3.6 and T1 | as A1, A4 |
| B2 compile-time identity and FMA | Not decided by the editor: both options in A5 with a recommendation; D-13 added; Chapter 10.5 and style decision 15 word the guarantee with its hypothesis meanwhile | 10.2, 11, 8 |
| B3 F13 and DS8 | Accepted: 16 levels, the stop from a counting wrapper, the returned entry matched bit for bit, the floor in the CSV; caption rewritten | 4.2, 4.3 |
| B4 the arbitrary user step | Accepted: $h = \varepsilon^{1/3} w$ in E4, DS4 and F4; 3.5 reads "replace $s$ by the scale on which $f$ varies" | 1.4 (ch. 3), 4, 5 |
| B5 the node-rounding term; relative accuracy | Accepted: the third term in 4.3 and S4; the presupposition $\lvert f \rvert \sim s^N \lvert f^{(N)} \rvert$ in 3.2, 5.2 and 10.1 with the thread's factors as `\val` macros; E7's model error corrected to about $2 \cdot 10^{-7}$ | 1.4 (ch. 3, 4, 5, 10), 4.3, 5 |
| B6 F6 caption | Accepted: reworded as suggested | 4.2, 1.4 (ch. 4) |
| B7 F14 scope, zeros, restart tail | Accepted: scope in the caption and S11; zeros clipped and counted; `restarted` column and the N = 2 tail named | 3, 4.2, 4.3 |
| B8 Sterbenz hypotheses | Accepted: S4 and S9 state them and the harmless alternative | 3 |
| B9 $w_k$ ambiguity; the estimate formulas | Accepted: rational $w_k$ in derivations, integer $c_k = d w_k$ in tables; $e_{i,j}$ and `abs_error` given explicitly | 2.1, 2.2 |
| B10 references for the step | Accepted: Dennis and Schnabel, Nocedal and Wright; Gill, Murray, Saunders and Wright as optional | 6 |
| B11 D-6; `\valRiddersEdgeEvalsAtOneHundredth` | Accepted: D-6 reclassified; DS10 gains the $x = 0.01$ row | 11, 4.3 |
| B12 Chapter 2 qualifiers | Accepted: the qualifier, the asymptotic definition and the sign of the leading term | 1.4 (ch. 2) |

Nothing was rejected outright; A6 and A8 were accepted in a modified form for the reasons given.

### 12.2 Stage-2 decisions (editor, 2026-10-23, from the evidence reviewer's verification of the thread's data and
### the generator's first results)

| Id | Finding | Decision | Where |
|---|---|---|---|
| S2-1 | The WebBook's Shomate coefficients are its own fit (January 2009) to Chase 1998; the monograph holds tables, not coefficients. The plan attributed "the Shomate form and the coefficients" to `chase1998` | Cite the WebBook (with the retrieval date, in its requested form) for the equations and the coefficients, Chase 1998 for the data; both entries verified and marked so in `refs.bib`; `shomate1954` added as optional, UNVERIFIED in content. The text says once that $H$ is the enthalpy relative to 298.15 K and never presents the closed-form digits as nitrogen's $c_p$ (the fit is 1.6e-4 below the JANAF table at 1000 K) | 1.3, 6, `refs.bib`, `code/shomate.hpp` |
| S2-2 | The 500 K junction of the two lower ranges is not a kink in $H$: the fits do not meet ($H$ jumps by 0.77 J/mol, $c_p$ by $-7.5 \cdot 10^{-5}$ J/(mol K), $\dd c_p/\dd T$ by $-17$ %). Options: (a) the real piecewise $H$ as a range break with a jump; (b) the piecewise $c_p$, a kink with a small jump; (c) two smooth branches joined by `max` | (a): the pitfall's example is the published function, labelled a range break; the lesson "evaluable is not smooth" covers kinks and jumps alike, and every number is generated (DS20). (b) and (c) are not used: (b)'s jump of $7.5 \cdot 10^{-5}$ already dominates at millikelvin steps, and (c) describes no real substance. F18's caption states that the data show a jump ($\order(1/h)$ for a straddling stencil) and names the kink case ($\order(1)$ or $\order(h)$) as the general statement. **For the author:** this changes the nature of E19/F18 from "kink" to "range break with a jump" | 1.3, 4.3 (DS20), 5 (E19), `code/shomate.hpp` |
| S2-3 | The accuracy at the default step measured at one point ($x = 1$) is luck by a factor of ten in its rounding part (`Central<6>` on $\mathrm{e}^x$ gave $3 \cdot 10^{-16}$ at $x = 1$, $8 \cdot 10^{-15}$ as a median) | DS6, DS7 and DS15 report the median over eleven points in $[0.5, 1]$ (same scale, same step), with the smallest and largest value in the CSV; the headline digits of 1.5 come from the medians (11 for $f'$, 4 for $f''''$, against the plan's "about ten" and "about four"). Single-point values that a walk-through needs are kept under `...AtOne` names | 4.3 (intro, DS6, DS7, DS15) |
| S2-4 | E20's inner solver as planned (fixed-point iteration to a tolerance $\eta$) does not behave as the noise model assumes: its output is a smooth function of $Re$ between the points where the iteration count changes, so the default-step derivative came out accurate to $7.5 \cdot 10^{-11}$, not to $\eta^{2/3}$ | The inner solver is bisection to a relative tolerance $\eta$, whose output is quantised at $\eta$ and shows the noise of 3.6 as modelled (default step $5.5 \cdot 10^{-5}$, fitted step $7 \cdot 10^{-8}$). The fixed-point results are generated too (`\valSolverFixedPoint...`) for one contrasting sentence: a solver's output is a different function, off by up to its tolerance, and the derivative obtained is that function's, which no estimate can see (Ridders on it: error $2.7 \cdot 10^{-11}$, estimate $1.4 \cdot 10^{-12}$). The text states which solver it used | 4.3 (DS21), 5 (E20) |
| S2-5 | At the single point $x = 1$, Ridders' estimate on the noisy function covers the error (ratio 0.17), contrary to E17's planned assertion | E17 does not assert a failure at one point; DS19 gains a coverage study over 300 points (clean 100 %; $\eta = 10^{-10}$: 89.3 %, worst ratio 67; $10^{-8}$: 90.7 %, worst 42), which is what 3.6 and E17 quote | 4.3 (DS19), 5 (E17) |
| S2-6 | DS9's planned "restarted = evaluations above a maximum" is a heuristic; and $1/x$ does not restart for $x < 0.25$ at N = 2 (finite values; the pole between the samples costs levels, not a restart) | Restarts are detected from a log of the calls (a non-finite value, or a level step that does not decrease); the S11 and F14 texts say "log and sqrt restart there" and describe $1/x$ as measured. DS9 extended to N = 3 and 4 for the cost table and split into one file per N | 3 (S11), 4.2 (F14), 4.3 (DS9, DS16) |
| S2-7 | DS10's "smallest x that still succeeds" has two readings, which differ by four octaves | Both are generated: a value is returned down to $2^{-34}$ (the derivation of D-6, now measured), ten correct digits down to $2^{-30}$, the evaluation budget exhausted from $2^{-27}$; F15 shows all three | 4.3 (DS10), 11 (D-6) |
| S2-8 | One source for the thread's coefficients and for the noisy function | The shared headers `code/shomate.hpp` and `code/noise.hpp`, quoted with the new `\snippetheader` macro (`macros.tex`); the "defined in each file" convention of section 5 is withdrawn | 5, `macros.tex`, `CMakeLists.txt` |
| S2-9 | `\lstinputlisting`'s `linerange={#2-/#2}` splits a hyphenated snippet name at its first hyphen and quotes the wrong lines or nothing, silently (all three authors) | `\snippet` (preamble, approved by the author) and `\snippetheader` brace the range ends, `linerange={{#2}-{/#2}}`; `tools/check.py` now rejects a `\snippet` or `\snippetheader` whose file is missing, whose markers are missing or repeated, or whose range is empty | `preamble.tex`, `macros.tex`, `tools/check.py`, 8 (item 18) |
| S2-10 | `\input{data/table-...}` inside a tabular gives "Misplaced \noalign" (the file hooks add tokens after the file) | One form for every body: `\tablebody{<name>}` in `macros.tex`, the TeX primitive `\@@input` behind a readable name; the generator keeps the trailing `\\` on every row; documented in `data/README.md` | `macros.tex`, `data/README.md`, 8 (item 20) |
| S2-11 | The methods author found 12 rows with an actual error of exactly 0 in `ridders-coverage-N1.csv` against `\valRiddersExactCasesFirst` = 11 | The files were from different runs of the same binary, and the runs differed: libm's `pow` (the UCRT's) returned results an ulp apart between runs on the same machine, in the closed forms of $1/x$ and the Runge function, which flipped exact-zero cases (255 to 289 rows per N changed slightly; Numerixx's own values were identical). The generator now forms every power by multiplication, the Runge derivative by complex multiplication, $\sqrt{x}^{\,1-2n}$ as $\sqrt{x}/x^n$, and the step grid $10^{k/4}$ from $\sqrt{\sqrt{10}}$ and exact powers of ten: correctly rounded operations only. The remaining libm calls (exp, sin, log, tanh) are shared by fn and the reference and were identical in every run. `data/` is generated on the gcc preset (GCC 16, UCRT); a regeneration elsewhere may differ in the last bits and in the exact-zero counts, which is why `data/` is committed and its diff reviewed | `tools/generate_data.cpp`, 4.3 |
| S2-12 | E18 (section 7.5) needs its printed numbers as macros | `\valCheck{Steps,Exact,Order}Digits`, `\valCheckRiddersDiff`, `\valCheckStencilActual`, `\valCheckRiddersEstimate`, `\valCheckSine{Derivative,AbsError,RelError,Spread,DigitsDerivative,DigitsScale}`, recomputed as `code/07-using-check.cpp` computes them | 4.3 (DS18), `data/README.md` |
| S2-13 | Corrections from the authors' stage-2 work: F12 needs the full width; E10's planned stencil bound ($10^{-9}$) was wrong (measured $1.2 \cdot 10^{-7}$, the truncation with the step scaled by $s = 1$); E1's measured ratio is 0.73; E17's chapter text quotes DS19's coverage, not the example's count; F18's caption says the data show a jump ($1/h$) | Recorded in 4.2 (F12, F18), 5 (E1, E10, E17) | 4.2, 5 |
| S2-14 | F2's planned range (1e2 to 1e-4 K) shows the central curve turning up below 1e-2 K, against Chapter 2's claim that halving the step quarters the error | The range is cut at 1e-2 K, where truncation still dominates, and the caption points ahead to Chapter 3 ("below 0.01 K the central difference stops improving"); the editor edited `figures/02-taylor-fd-slopes.tex` (the fundamentals author's file) for this, noted here | 4.2 (F2), `figures/02-taylor-fd-slopes.tex` |
| S2-15 | Stage-3 notes from the authors (idioms, Clang warnings, pgfplots, margin floats, listings) | Style decisions 17 to 21 in section 8 | 8 |
| S2-16 | The temporary gallery | `chapters/zz-figure-gallery.tex` shows F1 to F18 (F4b, F14b and the tableau diagram included) in their planned float types with their planned captions and labels, and every snippet of every example; input after the layout test; both go at stage 4 | `main.tex`, `chapters/zz-figure-gallery.tex` |
| S2-17 | The gallery did not build: "TeX capacity exceeded (main memory 5,000,000 words)". Measured per figure: the F14 scatter alone took 4.56 million words (2.6 million above the 1.93 million baseline) reading the two ten-column 2700-row files with an `x expr` per point; every other figure about 20,000 words; the listings 0.6 million | Slimming the columns did not help (4.56 million words again): the cost is pgfplots' record of about 480 words per plotted point. The generator writes `ridders-coverage-plot-N<N>.csv` with every third case (`\valRiddersPlotEvery` = 3, `\valRiddersPlotCases` = 900 per order), the two columns the scatter needs and the exact-zero clipping applied; F14 draws those and its caption says so (all cases enter the statistics and the full files); the editor changed the two `\addplot` lines of `figures/06-ridders-coverage.tex`, noted. The full book will carry the same figures, so this is the structural fix; if a later figure needs more, the options are pgfplots' `external` library (a preamble and build change) or LuaLaTeX, both for the author | `tools/generate_data.cpp`, `figures/06-ridders-coverage.tex`, `data/README.md`, 4.2 (F14) |

### 12.3 Unrelated issues noticed (the library, outside the book; for the author)

| Id | Observation | Where |
|---|---|---|
| U-1 | MSVC `/W4 /WX` reports C4459 (a local hides a global) for the locals `a` and `b` of `root_epsilon` when a user's translation unit has globals of those names | `include/numerixx/math.hpp`, lines 88 and 91 |
| U-2 | The header's `Estimate` documentation says `if (*r)` compiles (and tests the value against 0); under Clang with `-Werror` it is a `-Wfloat-conversion` error. The wording could say "may compile" or name the warning | `include/numerixx/deriv.hpp`, the `Estimate` comment |
| U-3 | The same comment says the conversion is not considered by "the operator templates of Boost.Multiprecision"; the comparisons do accept it (compile and are true), the arithmetic operators do not (D-14). "the arithmetic operator templates" would be exact | `include/numerixx/deriv.hpp`, the `Estimate` comment |
| U-4 | `default_step`'s comment and `docs/design.md` say the dropped constant "costs at most a factor 1.64" for $f'$; the model formula they cite gives 1.63 for `Central<10>` (D-16) | `include/numerixx/deriv.hpp` (the `default_step` comment), `docs/design.md` |
| U-5 | `nxx::real` accepts a type that the code then rejects: the module uses `numeric_limits<T>::radix`, `digits` and `denorm_min()`, which the concept does not require; the tests do not see the gap because their minimal type inherits `numeric_limits<double>` (R1-developer-01: a strong typedef built to the documented requirements passes `static_assert(nxx::real<Real>)` and fails inside the header). Either the concept requires them or the header stops using them | `include/numerixx/math.hpp` (the concept; line 81), `include/numerixx/deriv.hpp` (lines 150, 699, 1009), `tests/test_callables.cpp` (the minimal type) |
| U-6 | The precondition `h > 0` is an `assert`, compiled out under `NDEBUG`; a release build then returns a backward difference for `Forward<1>` with a negative step (fn called below x), while the multi-node rules and Ridders return `non_finite` (R1-developer-04, R1-library-03) | `include/numerixx/deriv.hpp` (the `assert`, line 1337; the two-node `difference`) |
| U-7 | The rationales of the order caps are not supported by the model they cite: "beyond P = 10 the gain is under half a digit even for 50-digit types" (the model's exponent gives 0.7 digit for `Central<12>` at 50 digits, 2.4 for $f''''$) and "beyond 6 the amplification would cost more than the order gains" (by the model `Forward<7>` still gains about 14 %, `Forward<8>` nothing) (R1-numerics-07, R1-evidence-05). The caps are sound design decisions; the stated reasons are not measurements | `include/numerixx/deriv.hpp` (lines 105 and 115), `docs/design.md` (line 60), `docs/PHASE2_API.md` (line 115) |
| U-8 | `std::expected::value()` on the result throws `std::bad_expected_access<error>` on failure (or terminates under `-fno-exceptions`) and converts silently to T through `Estimate`'s conversion, so `double cp = derivative(H, T).value();` compiles without a warning in a library documented as throwing nothing (R1-developer-03). The header's documentation and `docs/design.md` could name it | `include/numerixx/deriv.hpp` (the `Estimate` and `derivative` comments), `docs/design.md` |
| U-9 | The 24 Shomate coefficients reproduced in `code/shomate.hpp` are NIST Standard Reference Data ("(c) U.S. Secretary of Commerce ... All rights reserved" under the Standard Reference Data Act, WebBook footer checked 2026-10-23) in a file of an MIT-licensed repository; factual coefficients with attribution are commonly reproduced, but the terms are NIST's (R1-evidence-26). For the author to confirm; the attribution in the header comment stays. 2026-10-24: the header and the README state that the coefficients are NIST SRD, not covered by the MIT licence; written permission from NIST not yet asked | `docs/deriv-book/code/shomate.hpp` |
| U-10 | The check that a rule's integer weights are exact in T is a `static_assert` in the body of `stencil_derivative`, not a constraint: a requires-expression reports `derivative<Central<10>>` with `std::float16_t` as valid and the call then fails inside the library (R1-library-06) | `include/numerixx/deriv.hpp` (line 1009) |
| U-11 | `nxx::deriv::error` has no `std::formatter` (`std::println("{}", d.error())` does not compile), and `non_finite` covers six causes that a caller cannot tell apart when logging a failure (R1-developer-09) | `include/numerixx/deriv.hpp` (`error`) |
| U-12 | There is no public accessor for a method's default step; a caller who repeats a stencil at $2h$ and $h/2$ must rebuild `root_epsilon<P + N, T>() * max(|x|, 1)` by hand, and may write `* x` (R1-developer-12) | `include/numerixx/deriv.hpp` (`detail::default_step`) |
| U-13 | The public documentation of `error::non_finite` says "h is below the resolution of x" for every rule; the two-node f' rules fail only when both points coincide, every other rule when either side does not move (an internal comment says so). A documentation note (R2-codex-06) | `include/numerixx/deriv.hpp`, lines 42-45 and 1024-1028 |
| U-14 | `tests/test_ridders.cpp` and `tests/test_orders.cpp` bound the multiprecision types' `abs_error` by twice the first run's estimates, regression baselines rather than independent references; the tests' comments and `docs/design.md` say so, the repository's test rule ("never from Numerixx itself") does not (R2-codex-07, R2-cline-05) | `tests/test_ridders.cpp:313`, `tests/test_orders.cpp:632`, the repository's CLAUDE.md |
| U-15 | A kink at x, or so close that no level resolves it, is invisible to every rule and to Ridders' estimate: the central difference returns the mean of the one-sided slopes at every step. A documentation note on the smoothness assumption of `derivative` (R2-codex-01) | `include/numerixx/deriv.hpp`, the documentation of `derivative` |
| U-16 | The comment "the largest intermediate is 12! times 4! (`Central<10>`, n = 4)" names a bound the generator never reaches; the largest intermediate is 12! (the book's 4.1, R1-numerics-17b, R2-cline-11) | `include/numerixx/deriv.hpp:385` |

### 12.4 Stage-4 decisions (editor, 2026-10-23, from the authors' stage-3 reports and the integrated build)

| Id | Finding | Decision | Where |
|---|---|---|---|
| S4-1 | The practice author's two macro blocks (`signature`, `\sigindent`; `\breakcommas`, `\weighttablestart`, `\weighttableend`) and `\qrhead` were marked `%% NEW-MACRO` in the chapters; `\breakcommas` defined the active comma with a global `\gdef` at every use | Merged into `macros.tex`; the active comma is defined once at load time and only prints the ordinary math comma with a break point, `\breakcommas` only sets the comma's mathcode (local to the float's group); removed from the chapters | `macros.tex`, `08-api.tex`, `09-tables.tex`, `qr-quick-reference.tex` |
| S4-2 | T0 in 1.5 and on the quick-reference page differed in wording and in the typesetting of the `\val` macros (raw in Chapter 1) | The five rows live in one editor-owned file, `chapters/shared-method-summary.tex`, input by both tables, so they are identical by construction; both use `\num`; every raw `\val` in the chapters (63, mostly integers in Chapter 6) wrapped in `\num` per style decision 3 | `chapters/shared-method-summary.tex`, `01-why.tex`, `qr-quick-reference.tex`, Chapters 3, 5, 6, App. A |
| S4-3 | The 30 typed-number warnings of `check.py` were the two pencil-and-calculator tables (2.3: $x^3$ at 1; 3.2: a six-digit calculator on $\mathrm{e}^x$) and the prose numbers around them | Generated (HT in 4.3): the cube table is the library's values at $h = 0.1$ and $0.05$ printed to seven digits, which equal the exact decimal hand calculation digit for digit; the calculator is simulated by rounding every value and the quotient to six significant digits; the prose quotes `\valHandCalcExact`, `\valHandCalcErrorRatio` (21, where the draft said 20) and `\valHandCalcLastDiff`, and writes $1/100$ and $1/400$ for $h^2$. Zero warnings remain | `tools/generate_data.cpp`, `02-taylor.tex`, `03-rounding.tex`, 4.3 |
| S4-4 | Macros the authors asked for | Added: `\valAmplForwardSixFourth` (4019; Chapters 9 and 10), `\valEdgeForward{Two,Four}RelError` (E6 at 500 K, 4.4), `\valCpSecondModelRelError` (2.2e-7, 5.2 and 5.6), `\valStepCost{CentralTwoFirst,CentralTwoSecond,CentralTwoFourth,CentralFourFourth,CentralTenFirst}` and `\valStepConstantCentralTenFirst` (S5, App. B.4), `\valMinFloat` (T16 gains a "smallest normal" row), `\valNoiseRiddersCases` (300; 3.6, 6.4, 7.4, 10.2, which had quoted `\valRiddersPointsPerFunction` or, in 6.4, wrongly `\valRiddersCases` = 2700). The $2^{-34}$ and $2^{-31}$ edge limits stay typed with their attribution to the design record (6.4, 10.4); 6,222,216 is a table weight | `tools/generate_data.cpp`, `data/README.md`, the chapters named |
| S4-5 | PLAN 1.4 (Chapter 4) said the one-sided rules "gain one order per point"; every rule does, and the central rules' advantage is one node fewer for the same order, by antisymmetry | The text (4.4) says so; 1.4 corrected | 1.4 |
| S4-6 | The fundamentals author treats a smooth fitted correlation as inaccurate but not noisy (its derivative is the fit's, to $\varepsilon^{2/3}$), against the plan's "a fitted or interpolated table" among the noisy functions; Chapters 6 and 7 had followed the plan ("a solver or a fitted table", T1 "a solver or a fit") | The author is right: noise means irregular errors from one argument to the next, which a smooth fit does not have, whatever its accuracy. 3.6 keeps the paragraph; 6.4, T1 and the pitfall 7.8 now name the solver output (and "any function whose error changes irregularly from one argument to the next") and no longer a fit; 1.4 corrected. Interpolated tables stay in 7.9 as out of scope (a linear interpolant is kinked, a different failure) | 1.4, `06-ridders.tex`, `07-using.tex` |
| S4-7 | `\valCpBestRelError` (2.2e-12) and `\valFloatBestRelError` (5.8e-7) lie below the model's minimum | Expected: the lowest point of a measured V is the luckiest cancellation of the rounding errors (S2-3: luck by a factor of ten at a single point); 3.3 and 3.7 call them so and point the reader to the model's level. No change | 3.3, 3.7 |
| S4-8 | The Colebrook reference requested by the practice author | Entered as `colebrook1939` with the DOI the author believes correct and `note = {UNVERIFIED}`; cited in 7.4 | `refs.bib`, `07-using.tex` |
| S4-9 | Terminology (style decision 4): "stencil" was used in Chapter 1 (1.4, 1.5) before any definition; "node", "sample", "stencil" and "rule" are defined in 4.1's margin note, "node" and "weight" informally in 2.1 | 1.3 gains a one-sentence definition of a stencil (margin note) where the four methods are listed; 2.1 and 4.1 unchanged; the glossary in the back matter defines every term of the decision | `01-why.tex`, `99-backmatter.tex` |
| S4-10 | Repetition between chapters: the noise-coverage figures of DS19 were quoted in 3.6, 6.4, 7.4 and 10.2 | Kept as deliberate cross-references (each quotes the macro and points to 3.6); 6.4's wrong count fixed (S4-4). The list of the nine study functions appears in S11 (self-contained sidebar) and 10.2 (reference): kept | -- |
| S4-11 | Layout of the integrated build (83 pages: title, vi, 76): the notation table of the back matter did not fit below its heading and left a near-empty page; Sidebar 6.3 split after two lines; seven listings split across pages with a few lines stranded (1.4 `shomate-H`, 4.4 `validity-edge`, 7.2 `expected-basics` and `estimate-pitfalls`, 7.4 `newton` and `implicit-derivative`, 7.8 `kink`); Table 9.8 dropped to a page bottom under a gap because `topnumber` = 3 caps the weight tables stacked at a page top | The notation table is two tables; `\needspace` guards before the sidebar and the listings (a guard moves the whole block to the next page, leaving a few blank lines, which is the lesser evil); Chapter 9 sets `topnumber` to 4 and `main.tex` restores 3 after it. Margin floats never share a page with a full-width float (F14b and F15 on pages 35 and 36 against F12 on 32 and F14 on 34; F16 on 42 and F18 on 48 against T1 on 41); Chapter 2's title breaks as planned in the contents and the running heads; no overfull box; sidebars 2.1 and 5.1 still break after a displayed equation or four lines, which is tolerable for a breakable box | `99-backmatter.tex`, `06-ridders.tex`, `01-why.tex`, `04-stencils.tex`, `07-using.tex`, `09-tables.tex`, `main.tex` |
| S4-12 | The appendices were headed "Chapter A" and "Chapter B": `preamble.tex` hard-coded the word "Chapter" in the chapter title format | Changed to `\@chapapp`, which the report class sets to "Appendix" after `\appendix`: a one-token correction of the title format, not a design change, made by the editor and **listed for the author's approval** with the integrated PDF | `preamble.tex` |
| S6-1 | Stage 6, the editor's items before the authors' revision (`reviews/R1-decision.md` section 4; every decision approved by the author on 2026-10-23) | Done: the AI statement rewritten (the coordinating main session on Fable 5.1, the stage-1 plan reviews, the stage-2 verification, the authors' listings and the editor's generator, the transcribed declarations and the CMake file as the two exceptions, round 1, and the two independent reviews to come: OpenAI Codex and Cline, then with MiMo-V2.6-Pro, since 2026-10-24 with DeepSeek V4.1 Flash 1017, S6-11); `refs.bib` per the reference table (`sterbenz1974`, `ieee754-2019` corrected, `more2011` and `martins2003` added, every UNVERIFIED note removed, sections cited at the sentence); `articletitle=true` and the appendix naming (`\crefalias{chapter}{appendix}` in `main.tex`); the glossary and notation of the back matter; `code/colebrook.hpp` and `book::estimate_noise` in `code/noise.hpp`; `code/getting-started.cmake` and `\cmakefile`; the generator's extensions of 4.4; `check.py`'s signature and `\cmakefile` checks; the "restarts" label of F14b removed and `RiddersExactCases` renamed `RiddersZeroError` in its two users (Chapter 6 and the F14 comment); sections 2, 5, 8 (items 18, 22 to 25) and `data/README.md` | this plan, `chapters/00-frontmatter.tex`, `99-backmatter.tex`, `refs.bib`, `preamble.tex`, `main.tex`, `macros.tex`, `code/`, `tools/`, `figures/06-ridders-cost.tex` |
| S4-13 | The log reports `Font shape T1/cmss/b/n undefined` twice, from the tableau table of Chapter 6: `siunitx` in `detect-all` (text) mode takes the minus sign of a negative exponent from the TS1 encoding, which IBM Plex Sans lacks, so every such minus falls back to Computer Modern Sans; the two warnings are the bold accent cells of F12's error table, where the bold fallback does not exist either and medium is used. Cosmetic (the signs are the right size and nearly the right shape) and present since the stage-2 layout | Not changed: a fix is a preamble decision for the author (for example `\sisetup{mode = math}` for `\num`, or loading the TS1 companion of Plex). Noted for the review of the integrated PDF | `preamble.tex` (unchanged) |

### 12.5 Stage-6 close-out decisions (editor, 2026-10-24, from the authors' responses and the re-integrated build)

| Id | Finding | Decision | Where |
|---|---|---|---|
| S6-2 | Two quoted snippet lines exceeded 60 characters (`code/03-rounding-step.cpp` `too-small`, 61; `code/noise.hpp` `noise-estimate`, 62) | Both within 60: `noise-estimate` reflowed (the root-mean-square vector is `rms`, its header comment re-wrapped); `too-small` reads `if (!tiny && tiny.error() == error::non_finite)` under a `using nxx::deriv::error` and a one-line comment, a trivial edit of the fundamentals author's file, noted for the author. `code/getting-started.cmake` reflowed to 60-character comment lines | `code/noise.hpp`, `code/03-rounding-step.cpp`, `code/getting-started.cmake` |
| S6-3 | Chapter 7 over its allowance (14 pages at the authors' hand-in, total 91) | The practice author's three suggestions taken: the CMake reflow; E22's `jacobian-balances` quotes the loop of `jacobian(r, x)` only (the markers moved inside the function, the loop de-indented, the prose names the function and `J`); Sidebar 7.2 shortened (the bounds $Q = \max(2K, d)$, $h \le \Omega/2Q$, $\lvert x \rvert \le \Omega - Qh$ and the helper names `sum_fits`, `product_fits` dropped, every fact kept: pre-arithmetic detection, the exact radix comparison, the saturating floors, the final-quotient exception). Then the `\needspace` guards the author had removed were restored where the integrated layout stranded listing lines: `expected-chain` (two comment lines), `ridders-check` (its last three lines), `check-by-steps`; and nine sentences of 7.4.2, 7.4.3, 7.5, 7.6, 7.7, 7.8 and 7.9 were tightened by a few words each where a paragraph's last line held one to three words (no fact or cross-reference lost; the deferred-features list says "a noise-aware step", as 10.7 does). Chapter 7 is 13 pages (41 to 53) with every listing in one piece; total 88 | `07-using.tex`, `code/07-using-jacobian.cpp` |
| S6-4 | Chapter 4: the margin figure F8 (`fig:one-sided`) and the full-width F7 (`fig:stencil-nodes`, side caption) shared page 26 with their captions overlapping | F8 anchored before the paragraph that cites it; the two now share the page without a collision (the only margin float beside a full-width float in the book; Chapters 6 and 7 are clear as at stage 4) | `04-stencils.tex` |
| S6-5 | Back matter: the first notation table left half of page 77 empty after the authors' revisions moved the break | The split between the two tables moved ($\varepsilon$ to $D_i$ into the first); six glossary entries tightened without a change of meaning. Appendix A.6 repeated the list of `numeric_limits` members that 8.7 and Sidebar 7.1 give; it now points to 8.7 | `99-backmatter.tex`, `a1-floating-point.tex` |
| S6-6 | The two optional macros the authors asked for | Added: `\valExpRiddersFirstRelError` (DS15, 1.3e-14, Ridders' median for $f'$ on $\mathrm{e}^x$ in double), quoted in 6.5 beside the headline digits; `\valRiddersEdgeLimitRelEstimate` (DS10, 0.10, the relative estimate at $x = 2^{-34}$), quoted in 6.4 in place of the typed "a tenth of the value". 7.8 and 10.4 keep "an estimate that has grown to flag the loss" and may quote it in round 2 | `tools/generate_data.cpp`, `06-ridders.tex`, `data/README.md` |
| S6-7 | The reference table (section 6) cites Nocedal and Wright at 3.3 and Dennis and Schnabel at 3.4; the fundamentals author's response omitted the citations | Cited at the sentence, without the section qualifiers: the achemso style prints each qualified citation as an extra numbered entry ("Ref. 14, Sect. 8.1"), and two more would run the reference list to a second page. The sections stay in the table | `03-rounding.tex` |
| S6-8 | The thirteen response files against `R1-decision.md` | Every accepted finding addressed or its non-change justified; the only gap was S6-7. Recorded without change: E22's `jacobian-by-calls` is named in the prose, not printed; the `names` snippet of E11 is no longer quoted (8.1 names the file); T0's Ridders row stays as the decision letter settled it; Chapter 7 names "MSVC and clang-cl of Visual Studio 2026" without a compiler version because `ci.yml` names only the runner | this plan (section 5, E11, E17, E22) |
| S6-9 | The stage-6 layout (88 pages: title, vi, 81) | No overfull box, no undefined reference, `check.py` 0 errors 0 warnings, 24 example tests pass, `data/` diff limited to `values.tex`, `table-ridders-constants.tex` and `README.md`. Guards remain at 3.4 (7), 4.4 (14) and 6.6 (10); Chapter 1 needs none. Tolerated: the sidebars 3.2, 5.1, 6.2, 7.1 and 7.2 break across pages (breakable boxes, as at stage 4); Chapter 3's `noise-estimate-use` listing splits five and four lines over pages 20 and 21 (page 21 is full, so a guard would push the checklist's last item onto a page of its own); half a page blank at the ends of Chapters 5 and 8 and of Appendix B. The `T1/cmss/b/n` font warnings of S4-13 remain | -- |
| S6-10 | The AI statement | Confirmed as it stands apart from S6-11: every sentence true at the close of stage 6 ("the authors revised their chapters"); the last sentence names what is still to come. The version paragraph keeps its TODO for the release commit | `00-frontmatter.tex` |
| S6-11 | The author's change of the round-2 Cline model (2026-10-24): DeepSeek V4.1 Flash 1017 in place of MiMo-V2.6-Pro | "Cline with DeepSeek V4.1 Flash 1017" in the AI statement, in S6-1 above and in the approval note of `reviews/R1-decision.md`; the stage-7 row of the schedule (section 9), which had named Codex only, now lists the Cline review and its output file `reviews/R2-cline.md` | `00-frontmatter.tex`, `reviews/R1-decision.md`, this plan |

### 12.6 Stage-7 decisions (editor, 2026-10-24: the editor's items, the authors' responses and the close-out)

| Id | Finding | Decision | Where |
|---|---|---|---|
| S7-1 | The editor's items of `reviews/R2-decision.md` section 4, before the authors | Done: the AI statement (T5); `refs.bib` (`ieee754-2019` as a book entry with the standard's number as its series; the 17 citation qualifiers moved into the prose beside unqualified `\cite`s, 23 entries on one page); the DS9 margin extension (`\valRiddersMargin{Cases,QuarterCovered,SixtyFourthCovered}` = 14, 11, 14) with PLAN 4.4 and `data/README.md` rows, also for the two stage-6 macros; the figure comment (T4); `check.py` comparing whole declarations (library-03; a mutation test confirmed that one overload losing `[[nodiscard]]` is caught); `getting-started.cmake` with the full hash and "a tag or a full commit hash; CMake 3.25 or later"; Table 1.1 after the two checks and Table 8.2 at the head of 8.5 with `[h]`; the glossary's K and C; `noise.hpp`'s snippet comment cut; `.gitignore` for `tools/__pycache__/` | the files named |
| S7-2 | The thirteen `responses/R2-*.md` against the letter | Every row of section 3 is addressed; no gap. The fundamentals author's probe corrected T9's retry direction (k = 0 means the spacing is too large: retry with a tenth) and added "distrust a level near epsilon for a solver"; recorded in the letter's section 7, the decisions left as written; 7.4.3 and the quick reference make no retry statement and agree with 3.6 | `reviews/R2-decision.md` section 7, `03-rounding.tex` |
| S7-3 | One voice over the round-2 changes | The two-tolerance check reads `res && res->abs_error && *res->abs_error <= atol + rtol * std::abs(res->value)` in 1.5 and in E11's `ridders-check` with the same names (`res`, `atol`, `rtol`), the quick reference with the card's `r`; jump against kink consistent in 6.4, 7.5, 7.8, Table 7.1, Figure 7.2, 1.5 and the quick reference; 10.2 still said the estimate sees "a kink or a jump" as a large estimate, corrected by the editor to "a jump ..., nor a kink at x at all (6.4)" (not in the letter); the sign of C consistent in 2.3, 4.1, 9.1 and the glossary; "about three more digits than Central<2> (14 against 11)" in 6.5, 6.6, 7.3 and 10.3 | `10-guarantees.tex` |
| S7-4 | The `\needspace` guards (the practice author reports the one before `check-by-steps` leaving two thirds of a page empty in the standalone build) | All six kept; in the integrated build none leaves a gap: page 49 is full and the listing starts page 50. Guards: 3.4, 4.4, 6.6, 7.2 (two), 7.5 | -- |
| S7-5 | Layout of the integrated build | Figure 7.1 (margin) shared page 45 with the full-width Table 7.1: anchored after the Newton paragraph, now on page 46 alone. Pages 24 (Figures 4.1, 4.2) and 26 (4.3, 4.4) hold a text-width figure with a side caption and a margin figure without collision. Two listings split across pages are tolerated: `noise-estimate-use` (20/21, page 21 full) and `implicit-function` (47/48, at `double lo = 1.0, hi = 20.0;`); a guard before either would cost a page. Sidebars 2.1, 3.2, 5.1, 6.2, 7.1 and 7.2 break across pages, as before | `07-using.tex` |
| S7-6 | The authors' net additions produced three pages holding two or three lines (a third page of "About", Sidebar 6.3's last two lines on page 41, 7.9's last lines on page 54) | Trimmed by a few words each, the content kept: 6.4's kink sentence, 6.5's digits parenthesis, Sidebar 6.2's margin sentence, Sidebar 6.3's last sentence; 7.3's "other rows", 7.5's physical check, 7.6's two rules and the `-> mp50` sentence, 7.8's aliasing, domain-edge, float-fourth, one-step and value-but-poor pitfalls, 7.9's opening, cache and deferred paragraphs; the front matter's "How to read it", "The evidence", "The version" and the statement. 88 pages (title, vi, 81), every chapter at its stage-6 start page | `06-ridders.tex`, `07-using.tex`, `00-frontmatter.tex` |
| S7-7 | The AI statement at release | Made current: the author's approvals through round 1 and the instruction to close round 2 on the editor's recommendations; round 2 by the five reviewers (30 new findings, the round-1 revisions confirmed), OpenAI Codex (10 findings, compiled counterexamples) and Cline with DeepSeek V4.1 Flash through ClinePass (21 findings); the authors' final fixes; nothing "still to come". The version sentence and the CMake pin wait for the hash | `00-frontmatter.tex` |
| S7-8 | The release pass (stage 8), what remains | The full hash of the final-fixes commit in two places: `chapters/00-frontmatter.tex` line 49 (`\todo{the commit of \code{dev} at release}`) and `code/getting-started.cmake` line 9 (`GIT_TAG`). The "Draft" marker in `main.tex` (`\date{Draft}`) and `preamble.tex` (`\bkRunningDate`, the running head): the author names the release text (a date or version), the latter being a preamble value | `main.tex`, `preamble.tex`, the two files named |

# Milestone 3, Phase 10: the `nxx::roots` chapter through the manual's workflow

Status: approved by the user (2026-11-11), after the critic's review (`reviews/P10-plan-critic.md`) and the user's
decisions on it; Q6 decided: name the Pegasus method only. The second chapter through the workflow of
`docs/sphinx/CLAUDE.md`, to Phase 9's pattern: the roots chapter rewritten to the theory cap and the budgets, with two
figures on the overview, and reviewed in the full round. The shared pages are out of scope; what the chapter repeats
from them is cut and linked. Phase 9 left three items for this phase (`reviews/P9-decision.md` U4 (c) and (d);
`restricted-domain.md`'s overrun), settled below.

## What goes

The pages keep their shape (one job each, problem first, the prototypes' counts, every number asserted); every count
and claim I checked against `roots.hpp` and the tests holds. What goes: about 30% of the words (7,823 now, about
5,440 planned); the preambles; what `basics.md` already says (testing a `std::expected`, "nothing is thrown", the
`switch` default, expression templates, `int` starts, compile time); what the overview says, repeated on the guide
pages (bracket and open methods, Newton's tangent, "every method is called the same way"); the touches-zero case on
four pages; the rounding budgets in prose, which the programs' comments hold; the one-sided-stencil story told twice;
Pegasus's history. Two figures replace prose on the overview; the observer moves to the page whose job it serves.

## The pages

Budgets are prose words by the manual's rule, counted with `python docs/sphinx/count_words.py` and nothing else. They
are targets (the user, 2026-11-11): a page that needs more for content a reader needs may go over, with a one-line
reason here. "Holds" is the theory the page keeps; the rest is the book's.

| Page | Job | Now | Budget | Holds | Leaves to the book |
|---|---|---|---|---|---|
| `index.md` | The problem; a bracket, and why it holds a root; two kinds of method and the call; what `abs_error` covers, and `fn` as computed; where to go | 981 | 650 | A bracket (continuity and a sign change); bisection halves, Brent's method interpolates and falls back; the open methods finish inside the first sign change they cross; `abs_error` bounds the distance to a sign change of `fn` as computed, which rounding in `fn` moves by about its rounding error over the slope | Brent's algorithm and its bound, convergence orders, Dekker's history, the bound's derivation |
| `tutorial.md` | Teach by doing: Kepler's equation from a bracket, checked; the wrong bracket; from a guess with Newton's method | 876 | 620 | Nothing new; the bracket $[M - e, M + e]$ from $\lvert \sin E \rvert \le 1$ | The rounding budget of the check (the program's comment) |
| `choose-a-method.md` | From what you have to a method; measure on your own `fn` | 734 | 520 | Interpolation against halving in a sentence; the two worst cases in the table's rows; multiplicity, defined; a numerical derivative's cost | The orders of convergence, Toms748's algorithm, the teaching methods' stalls |
| `tolerances.md` | When the defaults are wrong: a root near 0, fewer digits for less work | 771 | 520 | $\text{tol}(x)$ and the stop test; the defaults (full precision, no scale assumed); the walk to 0 and `abs_tol`; `rel_tol` | Brent's own tolerance; why interpolation lands on a line's root at once |
| `find-a-bracket.md` | Neither a bracket nor a good guess: `search`, then solve | 511 | 380 | The search starts at the guess's scale; what it cannot find | Nothing |
| `restricted-domain.md` | Keep a solve inside `fn`'s domain | 772 | 560 | NaN as the signal; a bracketing method stays inside, an open one can step out; a failed numerical derivative is a bisection step in a bracket and `not_converged` from a guess; `Forward<2>` only where `fn` is tame near the edge | The stencil's error near an edge (the deriv chapter and the book) |
| `derivatives.md` | Give $f'$ and $f''$ to Newton's and Halley's methods: from `fn`, separately, numerically, from a bracket or a guess | 568 | 420 | Halley's curvature in a clause; 3 and 6 calls per point for the numerical forms; the bound never depends on the derivative | Orders of convergence; why a failed derivative is stored as 0 |
| `other-types.md` | `float`, `long double`, Boost.Multiprecision: the type comes from the start; what changes | 528 | 350 | The open methods' cap grows with the type's bits | Nothing |
| `handle-failures.md` | What each error means and what to try; recover with `or_else`; watch the steps; tell a pole from a root | 697 | 560 | The three errors; the observer; the pole test and when it lies | Why Newton's method cycles |
| `method-of-your-own.md` | Add a method: the protocol, what `solve` checks, the two convergence rules, a bracketing and an open example | 1065 | 650 | `start`, `advance`, `Step`; the checks and `invalid_step`; keep a point $\text{tol}$ inside the ends, and step at least $\text{tol}$ from a guess; `values` | Pegasus's and Chebyshev's theory and history; the overflow-safe arithmetic (the program's comments) |
| `costs.md` | The reference table: evaluations per method on five functions | 293 | 210 | None | None |
| `reference.md` | The API, generated | 27 | as is | | |

Per page, what changes beyond wording:

- **`index.md`.** The opening to about 55 words. Figure 1 replaces the prose on how bisection and Brent's method
  close in, the "52 and 9" sentence included: the text after the figure says what to see, in two sentences. The
  `<brent>` block and output stay. The open-method paragraph to about 90 words, without the tangent (the tutorial
  shows the step). `search` in two sentences and a link; poles and several roots in a clause and a link. "What you
  get back" keeps the two members, the exact-zero rule and the conditioning sentence; figure 2 replaces the cubic's
  paragraph: a caption, two "what to see" sentences and Brent's result; the `<cubic>` block and output stay. "Every
  method works in any real type" becomes a link. The Tutorial, Guide and Reference headings stay, their sentences go
  except "Start with [Choose a method]" under Guide (the user's decision on Q1; the same on `deriv/index.md` and
  `quadrature/index.md`, by the editor in stage 4). A toctree `:caption:` cannot replace the headings: Sphinx 9.2's
  LaTeX builder uses it only as the root document's contents title, and its writer drops it elsewhere. The overview
  and the tutorial are the pages that name the alias `rt`.
- **`tutorial.md`.** The roadmap to one sentence. "Checking the result" to about 170 words: the two-double reference
  and its reason, the actual error within `abs_error` give or take $2 \times 10^{-16}$ from `std::sin`'s rounding,
  the residual in one sentence; the residual's bound and the platform sentences go to the program's comment, where
  they already are. "From a guess": "both elements `double`" as a clause with a link to `derivatives.md`, which
  keeps the `int` pitfall (its `static_assert` stays in `roots_kepler.cpp`).
- **`choose-a-method.md`.** Opens with the three things the choice depends on. The table's two Bisection rows merge;
  the worst cases live in the rows only. The paragraph on how the bracketing methods differ to one sentence; the
  Newton-in-a-bracket paragraph to one sentence and a link to `restricted-domain.md`. "Multiple roots" stays and is
  the one place that defines multiplicity and the touches-zero case; the other pages link. The note on the teaching
  methods to about 35 words, with its one exception. "Measure" keeps the counting lambda and one sentence.
- **`tolerances.md`.** "The defaults" merges into the opening. The walk to 0 stated here once; the Brent and
  multiprecision paragraphs one sentence each; `other-types.md` links here. "`abs_tol` acts at every magnitude" keeps
  its block and output, with one sentence after it: choose `abs_tol` from the scale of the roots you expect, not as
  an accuracy setting; "a method may do better than its guarantee" goes. The stop test is stated here;
  `method-of-your-own.md` links.
- **`find-a-bracket.md`.** "Search, then solve" to about 30 words. "Choosing the guess" keeps the rule and the kelvin
  example, drops "a fiftieth, growing by half" and the second example; the `int` sentence a clause with a link to
  `basics.md`; the touches-zero bullet a link.
- **`restricted-domain.md`.** Three sections kept, cut inside (the user's decision): the two opening paragraphs
  become one; the van der Waals note one sentence; "Prefer a bracket" loses what the overview says; "A numerical
  derivative near the edge" keeps the `<stencil>` and `<one-sided>` blocks and says once, qualified, that `Forward<2>`
  works where `fn` is tame near the edge, as $\ln(x / 10^{-12})$ is, and fails where it is not, as $\sqrt x$ is: the
  asserted counter-case, one sentence, and a bracket as the safeguard. Its anchor stays (deriv's `one-sided.md`,
  `handle-failures.md` and `derivatives.md` link to it).
- **`derivatives.md`.** Phase 9's U4 (d): the `<numerical>` block shows both forms, in the bracket and from the guess
  2, with both counts (21 and 18, asserted already). The `int` pitfall and the "exactly the type" rule live here,
  with a link to `basics.md` for the why; `other-types.md` and `tutorial.md` link. The Halley sentence stays (6 calls
  per point, 30 from 2).
- **`other-types.md`.** The `<expression-templates>` block and its prose go (`derivatives.md` and `basics.md` say
  it); the program keeps its `cpp_dec_float_50` checks unshown. The multiprecision bullets become one sentence each
  with links; the compile-time sentence goes. The `float` sentence keeps its asserted claim, "a `float` `fn` on a
  `double` bracket compiles but is accurate only to `float`", and drops the MSVC clause (settled, Q4).
- **`handle-failures.md`.** The preamble (`basics.md`'s) goes; the error list keeps one clause per cause, each a link.
  New section "Watch the steps", about 60 words: `observer(o)` on the Newton cycle of the "Recover" example, its
  first six steps quoted (0, 1, 0, 1, ...), so the cycle the page describes is shown and asserted. "Check for a pole"
  stays, its heuristic paragraph shortened.
- **`method-of-your-own.md`.** The opening and "Two kinds of method" one paragraph of about 55 words. The page names
  "the Pegasus method" without its authors: "Dowell and Jarratt (1972)" is not in `docs/roots/design.md` or the
  book's `refs.bib` (the 1971 Illinois paper is listed in `PHASE6_API.md` and the header only), so under the
  Citations rule the name leaves the page (Q6). The "four habits" paragraph to one sentence (true values of $f$ in
  the report; a sign test without a product) and a pointer to the program's comments. "What `solve` checks" to three
  sentences and the link to the concept's entry; the stop test a link to `tolerances.md`. The `inside = 0` paragraph
  to about 70 words. "An open method" to about 170: the Chebyshev block, `values = 3`, the two call forms, the
  `min_step = 0` failure as the second convergence rule. "Watching the steps" goes (to `handle-failures.md`).
- **`costs.md`.** The provenance paragraph to two sentences (the prototypes; the FMA and libm caveats, linked). The
  tables and the Steffensen footnote stay.
- **Every page.** `namespace rt = nxx::roots;` moves into each program's includes block, as Phase 9 did with `dv`.
  "Bracket" is defined on each page that uses it, as a clause of at most 12 words. No page says "below".

No page is merged, split or dropped (settled: `method-of-your-own.md` stays one page).

## One voice: the terms

Each defined in a sentence or a clause at first use on each page that uses it, and one name for one thing across the
modules (deriv's list in `PHASE9_PLAN.md` holds for the shared terms):

- *root*: a solution of $f(x) = 0$. *bracket*: an interval at whose ends `fn` has opposite signs. *guess*: a starting
  point. *sign change*: what every result verifies. *evaluation*: one call of `fn`.
- *bracketing method*: keeps a bracket and shrinks it; *open method*: steps from a guess, and finishes inside the
  first sign change it crosses. Never "closed method", "iterative method" or "derivative method" for these.
- *bound*: `abs_error`, always; never "estimate" (deriv's and quadrature's `abs_error` are estimates, and
  `basics.md` draws the line). *tolerance*: `abs_tol`, `rel_tol`, and $\text{tol}(x)$ for their sum at $x$.
- The methods by their prose names: Brent's method, bisection, Newton's method, Halley's method, the secant method,
  Steffensen's method, the Illinois method, regula falsi, the Anderson-Bjorck method, the Pegasus method; `Toms748`
  by its code name, with "the algorithm of Alefeld, Potra and Shi" at its first use on `choose-a-method.md` only.
  *Ridders' method* on a roots page is `nxx::roots::Ridders`, linked at its first use; no roots page mentions the
  deriv one, and no deriv page the roots one (settled, Q5). *The teaching methods*: Ridders', regula falsi, Illinois
  and Steffensen's.
- *interpolate*: fit a line or curve through the points so far and step to its zero (`index.md` and
  `choose-a-method.md`). *multiplicity* and the touches-zero case on `choose-a-method.md` only.
- *pole*: a point where $f$ grows without bound; *domain* and its *edge*, as deriv (never "boundary").
- *step* and *report* on `method-of-your-own.md`: a `Step` is the method's report; *point*: the method's estimate in
  it. *observer*: `observer(o)`.
- $f$ in mathematics, `fn` in code and in prose about the parameter; $f'$ and $f''$, not "the slope" except where the
  tutorial explains Newton's step.

## The figures

Each as `docs/CLAUDE.md` describes: a program (`examples/<name>.cpp`, a test) that holds the table in a raw string and
asserts what the page, caption and alt text say; a script (`plots/<name>.py`) that draws only that table. Both use
`+ - * /` on functions of `+ - * /`, so the tables are the same on every toolchain that does not fuse a multiplication
and an addition. Colours: blue `#3987e5`, orange `#d95926`, text grey. The numbers come from a probe against the
header (GCC 16); the programs confirm or correct them.

1. **How the methods close in** (`index.md`, "Two ways to start"): `roots_closing_in`. The bracket's width after each
   evaluation, for bisection and Brent's method on $x^2 - 2$ over $[1, 2]$, on a logarithmic axis, from the steps
   `observer(o)` reports. What to see: bisection's width halves with every evaluation, a straight line from 1 at the
   second to $8.9 \times 10^{-16}$ at the 52nd; Brent's falls off a cliff, 1, 0.67, 0.086, $5 \times 10^{-3}$,
   $4.8 \times 10^{-3}$, $2.4 \times 10^{-7}$, $4.2 \times 10^{-14}$, $6.7 \times 10^{-16}$ at evaluations 2 to 9,
   its last width the `abs_error` it returns. The one fact: bisection gains one bit per evaluation whatever $f$ is;
   an interpolating method gains digits fast once the bracket is small, and the bracket is why it cannot fail. It
   replaces words: the overview's sentences on halving, interpolating and "52 and 9" become the figure's two "what to
   see" sentences. The flat step at evaluations 5 to 6 is inherent (Wallis's cubic shows the same) and gets one
   clause within the cap: a step can land next to the end it replaces and narrow the bracket little; the next makes
   up for it. `main` asserts the counts (52, 9), that bisection's width halves at every step, that Brent's last width
   is its `abs_error` and below $10^{-15}$, and that both roots are within their bounds of $\sqrt 2$. Reference:
   $\sqrt 2$ to 20 digits. No third series.
2. **A function that cancels near its root** (`index.md`, "What you get back"): `roots_cancellation`. The computed
   values of $(x - 1)^3$ written out, $x^3 - 3x^2 + 3x - 1$ by Horner's rule, and as a product, at 201 points
   $x = 1 + k \cdot 10^{-7}$ across $[1 - 10^{-5}, 1 + 10^{-5}]$, on a linear axis in units of $10^{-16}$. What to
   see: the product's values one smooth curve through 0 at 1; the written-out form's values noise across the window,
   0 or a few multiples of $1.1 \times 10^{-16}$ of either sign, crossing zero dozens of times, agreeing with the
   product only outside about $7 \times 10^{-6}$. The one fact: no method can find the root better than `fn` computes
   its sign, so write $f$ so that it does not cancel. It replaces the cubic's paragraph (about 130 words). `main`
   asserts three things: more than 30 sign changes of the written-out form in the window (the probe gives 41 on
   this grid, 37 on one shifted by a third of a step), exactly one for the product, and the farthest wrong sign or
   zero on the grid beyond $4 \times 10^{-6}$ (5.5 and $6.6 \times 10^{-6}$ on the two grids; the page's "about
   $7 \times 10^{-6}$" is `roots_overview.cpp`'s, asserted on its finer grid, and the text links it). A decimal
   step, not a dyadic one: an exact $x$ would sample the rounding on a lattice. The alt text states the asserted
   counts only.

Considered and not proposed: the walk to 0 on `tolerances.md`, the search's points on `find-a-bracket.md`,
$\lvert f \rvert$ near a pole on `handle-failures.md`: a sentence says each.

## The programs

| Program | Status | `main` asserts | Independent reference |
|---|---|---|---|
| `roots_overview.cpp` | kept; `rt` into the includes block | as now: $\sqrt 2$ within one ulp and within the bound, the factored root within its bound, the expanded form's noise within about $7 \times 10^{-6}$, Brent's zero $4 \times 10^{-6}$ off with `abs_error` the tolerance there | $\sqrt 2$ as two doubles (mpmath); 1 |
| `roots_closing_in.cpp` | new, figure 1 | as figure 1 | $\sqrt 2$ |
| `roots_cancellation.cpp` | new, figure 2 | as figure 2 | 1; the rounding analysis in `roots_overview.cpp` |
| `roots_kepler.cpp` | kept | as now: both roots within `abs_error` + $2.3 \times 10^{-16}$, the bound at most $2\,\text{tol}$, the residual's bound, `no_sign_change` on $[2, 3]$, the `int` pair rejected | $E$ to 100 digits (mpmath and the Bessel series) |
| `roots_choose_a_method.cpp` | kept | as now: the four counts on Wallis's cubic | the prototypes' counts |
| `roots_tolerances.cpp` | kept | as now: 1049, 41, 3, 1076, 51, 21, 52; every result within its bound; the minimum and the negative tolerance | exact roots of $x - c$; $\sqrt 2$ |
| `roots_find_a_bracket.cpp` | kept | as now: the prototype's bracket and 27 evaluations, the root within its bound, `search(fn, 0)` rejected | the prototype; Wallis's root (mpmath) |
| `roots_restricted_domain.cpp` | kept | as now: Newton's lowest point and `non_finite`, $e$ within its bound, the search's 29 evaluations and bracket, $v$ against mpmath, the pole without NaN, `not_converged` from the guess and the root in the bracket, the default step, the $\sqrt x$ counter-case from a guess and in a bracket, one root at 280 K | $e$; mpmath; closed forms |
| `roots_derivatives.cpp` | changed: the `<numerical>` block shows the guess form too | as now: every result within its bound, the counts 7, 6, 5, 21, 18, 30, 8 | Wallis's root (mpmath) |
| `roots_other_types.cpp` | kept; the `<expression-templates>` block no longer shown | as now: four types within their bounds, the `float` rounding case, the cap 213, the expression-template rejections | $\sqrt 2$ to 70 digits |
| `roots_handle_failures.cpp` | changed: an observed Newton run quotes its first six steps; asserts the iterates alternate 0, 1 and the run ends `not_converged` after 213 evaluations | as now, plus the cycle | mpmath; the three pole cases |
| `roots_method_of_your_own.cpp` | changed: the `<observe>` block and its output go; its assertions (the stall after nine evaluations, the five reports at 3, observed equals unobserved) stay as checks; the citation in its comment: Q6 | as now, less the quoted observer output | Wallis's root (OEIS A007493); the `static_assert` overflow cases |
| `roots_costs.cpp` | kept | as now: every cell, every result within its bound | the prototypes; closed forms |

Thirteen programs: eight kept, three changed, two new.

## The reader's tasks (`docs-reader`, from the pages alone)

1. $\cos x = x$: Brent's method in $[0, 1]$; Newton's method from 0.5 with $f'$ in the pair; Newton's method from
   0.5 with `derivative_of` (the guess form). Each within its bound of the Dottie number 0.7390851332151607, and
   say which page showed each call.
2. $e^{20 - 5000 / T} = 1$ from the guess 300 with `search` and Brent's method chained; check against $T = 250$.
3. The root of $x - 10^{-200}$ in $[0, 1]$ by bisection with the defaults and with `abs_tol = 10^{-12}`, counting
   evaluations; say why the counts differ and what the second result is good for; then Brent's method on the same.
4. $\ln x = -x$ from the guess 5 by Newton's method with $f'$ given (it leaves the domain), then recover with `or_else`,
   a search and Brent's method; check against the omega constant 0.5671432904097839.
5. Bisection written as a bracketing method of your own, run through `solve` on $x^2 - 2$ over $[1, 2]$: within its
   bound of $\sqrt 2$, and its count against the library's 52; then make it report a point outside its bracket and
   name the error you get.

## Decisions on the plan review

`reviews/P10-plan-critic.md`, with the user's decisions of 2026-11-11:

1. Budgets: adopted (the user), as targets; about 5,440 in all.
2. The counter into the repository: adopted (the main session): `docs/sphinx/count_words.py`, named in
   `docs/sphinx/CLAUDE.md`.
3. Drop figure 1: rejected (the user); kept on two conditions, met above: it replaces the overview's prose, and the
   flat step gets one clause within the cap.
4. Figure 2 replaces the cubic's paragraph, 201 points, three assertions: adopted.
5. `method-of-your-own.md` at 650, the cuts listed, no split: adopted; the Pegasus citation: Q6.
6. `restricted-domain.md`: keep three sections, cut inside, keep the `<one-sided>` block and the counter-case, 560:
   the user; the critic's larger option (drop the block) not taken. The main session's note stands: the merge would
   have kept the anchor; the sections stay for the reason that they are two ideas.
7. `other-types.md` at 350, the expression-templates block out, the compile-time sentence out, the MSVC clause out:
   adopted.
8. `tolerances.md`: adapted. The block and output of "`abs_tol` acts at every magnitude" stay (a reader should see
   the three digits), its prose becomes the one sentence, and "does better than its guarantee" goes.
9. Cross-page repetition: adopted, each case as listed under the pages.
10. `choose-a-method.md`: adopted (the merged rows, no preamble, the note at 35 words).
11. `handle-failures.md`: adopted, the observer section at about 60 words and six steps quoted.
12. `index.md` at 650 and `tutorial.md` at 620, with the cuts named: adopted.
13. `find-a-bracket.md` at 380: adopted.
14. `derivatives.md` at 420: adopted.
15. Q1: the user: headings kept, sentences dropped except "Start with [Choose a method]", on all three overviews in
    this phase; the `:caption:` option checked and found not to work (above).
16. The plan's own text: adopted.
17. Q5: settled, no rename.
18. Theory cap: adopted. Pegasus's history leaves the page, the residual arithmetic leaves the tutorial, the worst
    cases stay in `choose-a-method.md`'s table rows only.

Settled by the critic's and the editor's agreement: Q2 (no split), Q4 (the MSVC clause dropped), Q5 (no rename).

## Question settled at approval

6. **The Pegasus citation:** the page names the method only (the user, 2026-11-11); the program's comment may keep
   "Dowell and Jarratt (1972)" as a comment.

## Evidence gaps and notes (none contradicts the header)

- `other-types.md`: the `float` sentence, settled as above.
- `handle-failures.md`: "its steps go from 0 to 1 and back to 0, for ever" is asserted only as `not_converged` after
  213 evaluations; the observed run will assert the cycle itself.
- `restricted-domain.md`: "From a guess, a one-sided difference can avoid the problem" is as broad as the sentence
  Codex faulted in deriv (P9-codex, finding 1); the rewrite says it once, qualified, with the counter-case.
- `tolerances.md` states the stop test as `Criteria`'s comment does; since Phase 7 `solve` also stops where the point
  is within $\text{tol}$ of both ends (`bracketing_method`'s comment). The header's brief and the page agree; no
  change.

## Steps

As `docs/sphinx/CLAUDE.md`: the user's approval; evidence (the two figure programs and scripts, the three changed
programs, each passing with GCC, the figures shown to the user); draft; edit (with the two overview edits outside the
chapter); the review round and the decision letter; revision; release on the five toolchains with the Codex review.
Commits to `dev` at the end of each stage the user approves.

## Success criteria

- [x] Plan approved, after the critic's review (2026-11-11)
- [ ] Evidence: the two figures and the changed programs pass with GCC; the user has seen the figures
- [ ] Draft and edit: every page within its budget or with the editor's reason here; `numerixx_docs` builds clean
- [ ] Review round, decision letter and revision (`reviews/P10-*.md`)
- [ ] Release: all five toolchains, site and PDF, CI; the Codex review and its adopted findings

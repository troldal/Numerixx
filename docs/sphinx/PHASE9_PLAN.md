# Milestone 3, Phase 9: the `nxx::deriv` chapter through the manual's workflow

Status: complete (2026-11-10, 3c2bc95). Approved by the user 2026-11-09, after the critic's review of the plan
(`reviews/P9-plan-critic.md`) and the user's decisions on it; the reviews, decisions and responses of every stage are
in `reviews/P9-*.md`.

The pilot of the workflow in `docs/sphinx/CLAUDE.md`: the chapter is rewritten to the theory cap and the budgets, with
figures where they make a concept easier to see, and reviewed in the new round. Roots and quadrature follow in phases
of their own, to this plan's pattern. The shared pages (`basics.md`, `getting-started.md`, `core.md`) are out of
scope; what the chapter repeats from them is cut and linked. `deriv/book.md` stays until Milestone 5.

## What the chapter does well, and what goes

The pages have the shape the rules ask for (one job per page, problem first, every number asserted), and the
overview's one idea, truncation against rounding, is the model `docs/design.md` names. They keep their functions, which
print the same on every toolchain, and their programs. What goes: a quarter of the words (5,832 prose words now, about
4,360 planned); the repetition of `basics.md`, of the roots chapter and of each other; the formulas beyond what a
reader needs to choose a step (the step and error for every $P$ and $N$, the constant of the truncation error), which
the book holds; book links on every page; and two names for one thing. Two figures and one table replace prose, a
code loop and its output.

## The pages

Budgets are prose words (code, output, captions and alt text excluded). "Holds" is the theory the page keeps; the rest
is the book's.

| Page | Job | Now | Budget | Holds | Leaves to the book |
|---|---|---|---|---|---|
| `index.md` | The problem; the one idea and the default step; two kinds of method; the call; the result; where to go | 960 | 700 | The central difference; its truncation error falls fourfold when the step halves, its rounding error grows as $\varepsilon \lvert f \rvert / h$; the default step $\varepsilon^{1/3} \max(\lvert x \rvert, 1)$; a stencil and Ridders' method in a sentence each | Taylor's theorem, the constants, the error expansions, Richardson's tableau and its stopping rules |
| `tutorial.md` | Teach by doing: $dE/dM$ and $d^2E/dM^2$ of Kepler's $E(M)$, a function computed by a root solve | 773 | 700 | Nothing new | As the overview |
| `choose-a-method.md` | From what you need to a method; measure on your own `fn`; more digits in a multiprecision type | 720 | 550 | The accuracy order $P$: the error falls as $h^P$ | Why a higher order can be worse, beyond one sentence |
| `step.md` | When the default step is wrong: a step on `fn`'s scale, or Ridders' method | 521 | 450 | The default step assumes `fn` varies on the scale $\max(\lvert x \rvert, 1)$; for another scale, $\varepsilon^{1/3}$ times it; Ridders' method starts at a 32nd of the scale and halves | The balanced step for every $P$ and $N$, and its constant |
| `one-sided.md` | `fn` defined on one side of $x$: NaN outside, a one-sided method, near the edge, a formula that goes on | 775 | 500 | A one-sided method evaluates `fn` at $x$ and on one side only; Ridders' method restarts below a step where `fn` is NaN | The one-sided tableau and its floor |
| `higher-derivatives.md` | Ask for the $N$-th derivative; how accuracy falls with $N$; not in `float` | 532 | 400 | Rounding is divided by $h^N$; the default step grows with $N$; the digits each order keeps in a `double` | The step and error per order, the stencils' tables |
| `derivative-of.md` | A derivative as a callable; its copy of `fn` | 446 | 300 | None | None |
| `handle-failures.md` | Recognise a failure and a result that is not accurate enough; what to try next | 650 | 500 | `non_finite` and its causes; the estimate assumes `fn` accurate to about $\varepsilon$; the methods' flag, an estimate of at least a quarter of the value | Why the estimate can be optimistic |
| `costs.md` (now `accuracy.md`) | The reference table: evaluations and accuracy per method on three functions | 327 | 200 | None | None |
| `book.md` | The download and the PDF note | 128 | 60 | | |
| `reference.md` | The API, generated | | as is | | |

Per page, what changes beyond wording:

- **`index.md`.** Figure 1 replaces the six-step loop, its output and the prose that walks the table; the code shows
  the default call and the Ridders call. The truncation sentence loses its constant. "What you get back" shrinks to
  the two members of `Estimate`, the rule that only non-finite values are detected (said once in the chapter, here;
  the other pages link) and a link to `basics.md`. The stencil orders in one sentence, without the term "accuracy
  order": `Central<4>` and higher use more points and are more accurate on a smooth function. The Tutorial, Guide and
  Reference sections are the toctrees and one line each. The chapter's one book link.
- **`tutorial.md`.** Opens with the problem. "A simple method" becomes "a stencil"; Ridders' method is redefined in a
  clause; the root-solve rounding in "Checking it" becomes one sentence and the bound. The second derivative stays as
  the call and its check; the sentence on rounding divided by $h^2$ goes to Higher derivatives (linked). It quotes
  bounds, not output, since its function goes through `std::sin`.
- **`choose-a-method.md`.** The table's one-sided, near-the-edge and fixed-cost rows merge into "at the edge" and
  "near the edge". The higher-orders paragraph becomes two sentences and a link to Costs. "Measure on your own
  function" keeps the counting lambda, which no other page shows, with one sentence of commentary. The multiprecision
  section keeps its example and drops the rules, which `basics.md` holds. The only page that defines $P$.
- **`step.md`.** Figure 2 carries why the default step fails on the peak; the "sixteenth of the width" paragraph and
  the "output follows the last example" aside go. Of the formulas, only $\varepsilon^{1/3}$ times the scale for the
  central difference stays, with "`default_step` gives the default for any method and order"; Ridders' start of
  $\max(\lvert x \rvert, 1) / 32$ stays. The code blocks and their output stay.
- **`one-sided.md`.** The near-edge table (below) replaces the "Near the edge" code, output and most of its prose.
  The advice stays in three sentences: try both; keep the smaller estimate only when they agree; an estimate is not a
  proof, with the $\sqrt{1 - x}$ counter-case in one sentence, asserted. The one-sided-stencil paragraph and the note
  on the $10^{-9}$ scale become one sentence each; the note stays here, and `handle-failures.md`'s warning becomes a
  link to "A formula that goes on". Links to the roots chapter's "A numerical derivative near the edge".
- **`higher-derivatives.md`.** The table stays. The $\varepsilon^{2/(2+N)}$ sentence and the Ridders-cost paragraph
  go (Costs holds the counts); the `order` enumeration in a clause; silencing the warning is the macro plus one
  sentence on the pragma.
- **`derivative-of.md`.** The callable example and the `std::ref` state example. The Newton and Halley paragraph
  becomes one link to `roots/derivatives.md#numerically`; the `float` note goes (Higher derivatives says it).
  `roots/derivatives.md`'s note on the copy of `fn` becomes a link to this page's "Its copy of `fn`" (the user's
  decision; the one edit outside the chapter).
- **`handle-failures.md`.** Links to `basics.md` for how to test a result and for "test for success, not for a
  particular error". The warning becomes a link to `one-sided.md`. The 7% and 70% case about 40 words shorter. The
  causes list and the note on a derivative near 0 and a noisy `fn` stay.
- **`costs.md`.** Renamed from `accuracy.md`, with its program, so that the three modules name their reference table
  alike (`roots/costs.md`, `quadrature/costs.md`, toctree entry "Costs"); the title keeps the form "nxx::deriv:
  evaluations and accuracy per method". Three links change (`index.md`, `choose-a-method.md` twice). The "three things
  stand out" bullets go: the table shows them, and the pages they repeat link to it. One sentence keeps the Ridders
  methods' typical and maximum evaluations, linked to the reference, since no other page states them now.
- **`book.md`.** The download sentence and the PDF note.
- **Every page.** `namespace dv = nxx::deriv;` moves into each program's includes block, which the page shows as "The
  program needs"; the sentence "the examples write `dv`" stays only in the overview and the tutorial (which names `rt`
  too).

No page is merged or split: each guide page is one task.

## One voice: the terms

Each defined in a sentence or a clause at first use on each page that uses it (readers land anywhere), and one name
for one thing across the modules:

- *central difference*: `Central<2>` for $f'$. *stencil*: `Central<P>`, `Forward<P>`, `Backward<P>`, a fixed formula
  with a step; never "plain finite difference", "simple method", "rule" or "formula" as its name.
- *accuracy order* $P$ (the error falls as $h^P$), used on `choose-a-method.md` only; *derivative order* $N$, or "the
  second template argument"; never a bare "order" where both could be meant.
- *Ridders' method*: `Ridders`. *The one-sided Ridders methods*: `ForwardRidders` and `BackwardRidders`. *The Ridders
  methods*: all three.
- *truncation error*: the error of the formula itself; *rounding error*: from the rounding of `fn`'s values; *machine
  epsilon* $\varepsilon$, $2.2 \times 10^{-16}$ for a `double`.
- *step* $h$; *the default step*; *a step of your own*; *scale*: $\max(\lvert x \rvert, 1)$, on which the default step
  assumes `fn` varies.
- *estimate*: `abs_error`, never a bound; *error*: the actual error, against a closed form.
- *evaluation*: one call of `fn`; *one-sided*; *domain* and its *edge* (never "boundary").
- $f$ in mathematics, `fn` in code and in prose about the parameter.

## The figures

Each as `docs/CLAUDE.md` describes: a program (`examples/<name>.cpp`, a test) that holds the table in a raw string and
asserts what the page, caption and alt text say; a script (`plots/<name>.py`) that draws only that table. The data are
computed with `+ - * /` and `std::sqrt` at decimal steps $h = 1, 0.5, 0.2, 0.1, \ldots$, each the correctly rounded
quotient of two integers, so the same on every toolchain (stage 2, approved by the user 2026-11-10: dyadic steps
$2^{-k}$ make $2 \pm h$ exact, and the rounded differences repeat or vanish, drawing a flat bottom in figure 1 and
false dips to $3 \times 10^{-16}$ in figure 2; `deriv_trade_off.cpp` records why). Smaller
steps to the right, so that the story reads left to right. Colours: blue `#3987e5`, orange `#d95926`, references in
the text grey. The numbers below come from a probe; the programs confirm or correct them.

1. **Error against step** (`index.md`, "The idea, and its catch"): `deriv_trade_off`. The relative error of the
   central difference of $\sqrt x$ at 2 for $h$ from 1 to $10^{-13}$, on logarithmic axes; the default step as a
   vertical line; Ridders' method's error as a horizontal line. The reader sees a V: the error falls fourfold with
   every halving down to about $5 \times 10^{-6}$ (truncation), then rises raggedly, about as $1/h$ (rounding); the
   default step, $1.2 \times 10^{-5}$, stands at the bottom; Ridders' line lies about a hundred times below it. The
   one fact: the best step is a balance, and the default step finds it. `main` asserts the fourfold fall on the
   truncation branch, the minimum within a factor 8 of the default step with the default step's error within ten times
   it, the rise below $10^{-7}$, and Ridders' error within its estimate and 50 to 500 times below the curve's minimum.
   Reference: $1 / (2\sqrt 2)$.
2. **The valley moves with the function's scale** (`step.md`, "When the default step is too large"): `deriv_scale`.
   The same curve for $\sqrt x$ at 2 (scale 1) and for the peak $x / (x^2 + a^2)$, $a = 10^{-4}$, at
   $2 \times 10^{-4}$ (scale $10^{-4}$), $h$ from 1 to $10^{-14}$, the default step marked. The reader sees two V's of
   one shape, the peak's about four decades to the left; the default step sits at the bottom of one and high on the
   truncation slope of the other. The one fact: the right step depends on the function's scale, not on $x$ alone.
   `main` asserts that the minima's steps differ by a factor between $2^{10}$ and $2^{16}$, that the peak's error at
   the default step is between $10^{-4}$ and $10^{-3}$ (the page's "three or four digits"), and that the scaled step
   $\varepsilon^{1/3} a$ gives an error within ten times the peak's minimum. References: both closed forms.

## The near-edge table (`one-sided.md`)

$1 / (1 + x)$ with NaN below 0, at $d = 10^{-2}$, $10^{-6}$ and $10^{-10}$ from the edge: for Ridders' method and
`ForwardRidders`, the relative error, the estimate and the evaluations, printed by `deriv_one_sided.cpp` and quoted
as output. The reader sees that far from the edge both work and Ridders' method is the more accurate; near it,
Ridders' estimate grows to about $10^{-5}$ and its cost to 64 evaluations, the most it may make, as it restarts below
each step that met NaN, while the one-sided method stays near $2 \times 10^{-13}$ and 10 evaluations. `main` asserts
that both estimates cover their errors at each $d$, the counts (64 at $10^{-10}$; at most 10 for `ForwardRidders`),
and which method is the more accurate at $10^{-2}$ and at $10^{-10}$; the central difference's failure at $10^{-6}$
is asserted already. Reference: $-1 / (1 + d)^2$.

## The programs

| Program | Status | `main` asserts | Independent reference |
|---|---|---|---|
| `deriv_overview.cpp` | changed: the steps loop and its output go; keeps the default call and the Ridders call | the default step's error within ten times the best of figure 1's curve; Ridders' error within its estimate, 50 to 500 times smaller than the central difference's, the estimate 1 to 10 times the error | $1 / (2\sqrt 2)$ |
| `deriv_trade_off.cpp` | new, figure 1 | as figure 1 | $1 / (2\sqrt 2)$ |
| `deriv_kepler.cpp` | kept | as now: the ten decimals, the central difference's bound, Ridders' coverage and gain, 8 to 20 evaluations, $f''$ within its estimate below $10^{-10}$ | the closed forms at the solved $E$ |
| `deriv_choose_a_method.cpp` | kept | as now: each stencil's bound, two digits gained, Ridders within its estimate, the 50-digit result within its estimate near $10^{-46}$, the stencil's two thirds of the digits | $-2x / (1 + x^2)^2$; $1 / (2\sqrt 2)$ in `cpp_bin_float_50` |
| `deriv_step.cpp` | kept | as now: the default step and its error, the scaled step's ten digits, Ridders' 36 and 14 evaluations with coverage, `non_finite` at $10^{-300}$ | $(a^2 - x^2) / (x^2 + a^2)^2$ |
| `deriv_scale.cpp` | new, figure 2 | as figure 2 | both closed forms |
| `deriv_one_sided.cpp` | changed: the near-edge table replaces the "near" block; keeps the edge and kink blocks, the $\sqrt{1 - x}$ and vanishing-polynomial cases | as now, plus the table's claims | $-1 / (1 + x)^2$; $1 - 1 / (2\sqrt{1 - x})$ |
| `deriv_higher.cpp` | kept | as now: the orders' bounds and digits, 3, 4, 5 evaluations, Ridders a hundred times better at every order, the `float` result's digits | $-1/4$, $1/4$, $-3/8$, $3/4$ |
| `deriv_derivative_of.cpp` | kept | as now: bit for bit `derivative<Ridders>`, within the estimates, the copy's and the original's counts | $1 / (2\sqrt x)$ |
| `deriv_handle_failures.cpp` | kept | as now: `non_finite` at $10^{-12}$ and $10^{-300}$, 7% off within its estimate and below the flag, the scaled step's $10^{-12}$, 65% to 75% off with the flag raised | $1 / (2\sqrt x)$ |
| `deriv_accuracy.cpp` | renamed `deriv_costs.cpp` | as now: the table, every Ridders estimate present and covering | the three closed forms |

Every program moves its `dv` alias into its includes block. Eleven programs: nine kept, changed or renamed, two new.

## The reader's tasks (`docs-reader`, from the pages alone)

1. The derivative of $x^x$ at 1.5 with an error estimate, checked against $x^x (\ln x + 1)$; say whether the estimate
   covered the error.
2. The second derivative of $\ln(1 + x^2)$ at 0.3 with an error estimate, and a test, from the estimate, of whether it
   is accurate to a relative $10^{-8}$; checked against the closed form.
3. The derivative of $\sqrt{1 - x^2}$ at 0.999, where the function is defined only up to 1 (NaN beyond): choose the
   method, check against $-x / \sqrt{1 - x^2}$.
4. The derivative of $e^{-(x/10^{-3})^2}$ at $5 \times 10^{-4}$: the default step and a step of your own, each against
   the closed form; say which page told you what to do.
5. Solve $x e^x = 2$ by Newton's method with a derivative from `nxx::deriv`, and check the root against Brent's
   method in a bracket.

## Decisions on the plan review

`reviews/P9-plan-critic.md`, with the user's decisions of 2026-11-09:

1. Drop figure 4 (`deriv_orders`): adopted (the user). The page keeps its table.
2. Figure 3 (`deriv_edge`) as a table: adopted (the user), printed by `deriv_one_sided.cpp` rather than a program of
   its own, since it belongs to that page's one function.
3. Figures 1 and 2 trimmed: adopted (the user). Figure 1 drops the marked Ridders steps and the step-reading
   machinery; figure 2 drops the flat-top and "no correct digit" assertions, which the text does not use.
4. The steps code block goes: adopted.
5. Budgets and their savings: adopted (the user), about 4,360 words. Adapted in one point: the tutorial keeps its
   second-derivative section, as the call and its check, because it is the tutorial's last step and the reader's task
   2 depends on a worked $f''$; only its theory sentence moves.
6. `derivative-of.md` at 300 words, the roots note a link to it: adopted (the user). The alternative, dropping the
   page into the overview and a header comment, rejected by the user.
7. `costs.md` to 200 words: adapted. The bullets and the closing paragraph go, but one sentence keeps the Ridders
   methods' typical and maximum evaluations, which finding 5 removes from `higher-derivatives.md`; one page must state
   them, and the reference table is the place.
8. `one-sided.md`: adopted; the note stays there, the warning in `handle-failures.md` becomes a link.
9. Book links on the overview and `book.md` only, `book.md` at 60 words: adopted (the user).
10. `choose-a-method.md` rows merged, commentary to one sentence: adopted.
11. The `dv` alias into the programs' includes blocks: adopted; the overview and the tutorial still name it.
12. Terms: adapted. "Accuracy order" is used on `choose-a-method.md` only, as proposed. "Scale" is still defined in a
    clause on each page that uses it, because `docs/sphinx/CLAUDE.md` requires a term defined at first use on each
    page; the definition is the formula itself, so it costs a few words.
13. The plan's own text: adopted; the figure descriptions are what to see, the data and the one fact; the tableau
    paragraph is gone.
14. Theory cap: adopted; the constant $h^2 / 6 \lvert f''' \rvert$ leaves the overview, the $P$ and $N$ step formulas
    leave the three pages, Ridders' starting step stays.

## Questions settled at approval

1. The rename `accuracy.md` to `costs.md` (and `deriv_accuracy.cpp` to `deriv_costs.cpp`): adopted (the user,
   2026-11-09).
2. `basics.md`'s "the plain finite differences" for what the headers and the chapter call stencils: left for the
   shared pages' phase (the user, 2026-11-09).

## Stage 4: the edit

The author's three departures from the plan, checked: (1) the costs table shows `Central<10>` less accurate than
`Central<2>` on $1/x$ at 2, so the overview says that `Central<4>` and `Central<6>` are the more accurate ones on a
smooth function, not "`Central<4>` and higher"; (2) the Ridders methods' typical and maximum evaluation counts per
order are the header's documented properties, which no program of the manual computes, so `costs.md` now points to
the reference for them instead of repeating them, and the overview's "10 to 16" links to `nxx::deriv::Ridders`;
(3) the near-edge sentences (64 evaluations, the one-sided method's 10) are asserted by `deriv_one_sided.cpp`. Edits
for wording: `step.md`'s figure moved after the default step's output, which it explains (and which left most of a
PDF page blank before it); `one-sided.md` names the methods instead of "a central method", says what fails at
$10^{-6}$ (the central difference, asserted), and "can be the more accurate one" for $\sqrt{1 - x}$, as asserted at
0.999; the page's roots link is about Newton's method, not a root solve in general; `choose-a-method.md` defines
"scale" and "domain" at first use; `handle-failures.md` introduces the one error value in a full sentence. Prose
words by the editor's count (inline code and formulas as one word each): index 622, tutorial 573, choose-a-method 486,
step 358, one-sided 442, higher-derivatives 352, derivative-of 261, handle-failures 397, costs 135, book 27, all
within budget.

## Steps

As `docs/sphinx/CLAUDE.md`: the user's approval; evidence (the two figure programs and scripts, the changed programs,
each passing with GCC, the figures shown to the user); draft; edit; the review round and the decision letter;
revision; release on the five toolchains with the Codex review. Commits to `dev` at the end of each stage the user
approves.

## Success criteria

- [ ] Plan approved, after the critic's review
- [ ] Evidence: the two figures and the changed programs pass with GCC; the user has seen the figures
- [ ] Draft and edit: every page within its budget or with the editor's reason here; `numerixx_docs` builds clean
- [ ] Review round, decision letter and revision (`reviews/P9-*.md`)
- [ ] Release: all five toolchains, site and PDF, CI; the Codex review and its adopted findings

# Milestone 3, Phase 4: a site organised by module, in a teaching voice

Status: approved by the user (2026-10-29), with the guide kept as seven pages; revised after the reviews of the plan
(`simplicity-critic`, `docs-reviewer`, `numerics-reviewer`;
the record is at the end). It replaces the planned Phase 4 (the roots explanations), whose two topics move into the
roots overview.

## Why

The site follows Diataxis at the top level: Tutorials, How-to, Reference and Explanation, each holding every module's
pages. Numerixx is a collection of mostly independent modules (two now; `ROADMAP.md` plans seven more: integrate,
optimize, interpolate, lsq, multiroots, poly and later others), and a reader comes for one of them, as to a chapter of
the GSL manual, a Boost library or a SciPy subpackage. Split four ways, a module has no place where it can be seen
whole, and each top-level section grows into a long mixed list.

The pages also read as annotated references: dense, every statement qualified, the reasons left out. They should teach,
with more descriptive text, while the book (one book with a part per module, `ROADMAP.md` 6) stays the place for the
full theory.

## The structure

Module first; Diataxis's four kinds of writing kept inside each module.

```
index.md                 Introduction: what Numerixx is; the modules by task, with links; how the site is organised
getting-started.md       Adding Numerixx to a project; a first program (moved from tutorials/)
basics.md                Working with Numerixx: what every module shares, explained once
core.md                  Reference: nxx::real and nxx::math (moved from reference/), under the same caption
roots/
  index.md               Overview: the problem, brackets and guesses, the call, what a result means, where to go
  tutorial.md            Kepler's equation (rewritten in the new voice)
  choose-a-method.md     \
  tolerances.md           |
  find-a-bracket.md       |  The guide: the seven how-to pages, moved and rewritten, one task each
  restricted-domain.md    |
  derivatives.md          |
  other-types.md          |  (its generic parts move to basics.md)
  handle-failures.md     /
  costs.md               Reference: evaluations per method (as now)
  reference.md           Reference: the API (Breathe, as now)
deriv/
  index.md               A short page: what the module does, the book, the reference; the full chapter is Phase 5
  book.md                The book (moved from explanation/derivatives.md)
  reference.md           The API (moved from reference/deriv.md)
examples/                Unchanged: a page shows code only from these programs, built and run as tests
```

- **Navigation.** `index.md` holds one toctree per area (Start here; Roots of functions; Numerical derivatives); each
  module's `index.md` holds that module's pages in its own toctree, so the sidebar shows one entry per module and opens
  the one being read. At nine modules that is nine entries, not fifty. (To check: Lumina collapses nested toctrees.)
- **The guide stays seven pages**, moved and rewritten, not merged into one. The voice rules lengthen every section
  (the sample passage doubles); separate pages keep each program feeding its own page, keep links page-level, and keep
  each page's git history through the move. A page that grows past about 1000 words is split.
- **`basics.md`** holds what both modules share, so it is said once: real types (`float`; `long double` on the tested
  toolchains; Boost.Multiprecision; expression templates and "return T"; `nxx::real` for a type of one's own; reading
  a constraint failure in a compiler message); results as `std::expected` (testing before reading, `and_then`,
  `or_else`, `value_or`, formatting the error; Numerixx itself never throws, while fn's exceptions and `value()`'s
  `std::bad_expected_access` do); what each module's result means (`nxx::roots::Root::abs_error` is a bound and a
  `T`; `nxx::deriv::Estimate::abs_error` is an estimate and an optional, and `*r` on a derivative tests its value
  against 0); compile-time use. Its code comes from its own programs (`basics_*.cpp`) with roots examples; deriv
  examples come with the deriv chapter (Phase 5). What is particular to roots stays in the roots pages: `abs_tol` in
  a multiprecision type, the open methods' iteration cap, `search(fn, 0)` deducing `int`, `Bracket{1.0, 2.0}` being
  a `double` bracket.
- **No "how it works" page.** The roots overview takes the two explanation topics (bracketing and open methods; what a
  result means); conditioning at a multiple root, exact zeros and poles go there or into the failures guide. A page of
  its own is deferred until a reader reports a gap; the book's part on roots (`ROADMAP.md` 6) comes later.
- **The reference.** Breathe directives appear only in `roots/reference.md`, `deriv/reference.md` and `core.md`; every
  other page links with cross-reference roles. The build's checks count rendered declarations across the site, so a
  directive on a teaching page would render a declaration twice and fail the build.

## The voice

Rules for every page, beside the existing ones (code only from tested programs; every number asserted by its
program or linked to what asserts it; the reference generated from the headers):

1. **Problem first.** A page opens with what the reader is trying to do and what goes wrong without the feature, in
   plain words, before any code.
2. **Explain the code.** Before a block, what it does and why; after it, what it prints and what that means. Output
   quoted on a page comes from the program: it keeps the expected output in a raw string that `main` compares with
   what it printed, and the page includes that string. Output that differs between toolchains (a count through libm,
   the last digits of a libm result) is described with the bounds the program asserts, not quoted.
3. **Edge cases out of the main line**, in a `note` box (MyST's `{note}` directive) or a link to the reference,
   unless ignoring them gives a wrong result or a solve that does not end: those stay in the text (the walk to 0 in a
   multiprecision type; `abs_tol` acting at every magnitude).
4. **Define a term at its first use on each page** (bracket, sign change, tolerance, conditioning), in a sentence:
   readers land anywhere.
5. **One name for one thing.** f in mathematics, `fn` in code and in prose about the parameter; mathematics typeset
   ($x^2 - 2$), code in code (`x * x - 2.0`). "The default" only where C++ has one (`solve` has no default method:
   say "the recommended method").
6. **The book's tutorial voice** (`docs/deriv-book/CLAUDE.md`): second person, short sentences, a paragraph where it
   explains better than a list. An overview fits on about one or two screens; overview and tutorial pages end with
   what to read next.

## Build changes

- MyST `dollarmath` with MathJax, for the formulas (tol(x), Kepler's equation, Newton's step, the conditioning terms).
  MathJax loads from a CDN, so a site read offline shows the TeX; accepted. The reference keeps the header comments'
  plain-text formulas.
- `myst_heading_anchors = 3` (planned as 2; the links to `basics.md`'s level-3 sections need 3), so that links to a
  section (`basics.md#results`) resolve under `-W`.
- No `colon_fence`: the backtick form of `{note}` needs no extension.

## Phases

- **Phase 4a, the move:** one commit that only moves files (`git mv`), fixes links, toctrees and the programs'
  first-line comments, and gets a green build; no prose changes, so that git follows every page and a review needs only
  the build. The site is not hosted, so no outside link breaks; this phase settles the paths, which become public at
  first hosting.
- **Phase 4b, the rewrite:** `index.md`, `basics.md` and its programs, the roots chapter in the new voice (overview,
  tutorial, the seven guide pages, costs), deriv's short `index.md`, the build changes; programs keep their checks and
  gain those of new claims.
- **Phase 5:** the `nxx::deriv` chapter in the same shape (the old Phase 5), with deriv examples in `basics.md`.

## Reviews

Phase 4a: the build and its checks. Phase 4b: `docs-reviewer` (with its rule updated, below), `numerics-reviewer`, a
docs-only reader (tasks: a table of Kepler roots over M; van der Waals from a guess; Newton with a numerical
derivative; a root in `cpp_dec_float_50` with expression templates), then Codex.

## Success criteria

- [ ] 4a: the move; `numerixx_docs` builds with no warning and its checks pass; the programs' tests pass
- [ ] 4b: the rewritten pages follow the voice rules; every number asserted or linked; all programs pass on the five
  toolchains, locally and in CI; the docs-only reader's code compiles and gives correct results; the reviews and Codex,
  with the adopted findings recorded here

## Changes to other files

- `docs/plan.md`, Milestone 3: Phases 4a, 4b and 5 as above; "module first, with Diataxis's four kinds of writing
  inside each module".
- `CLAUDE.md`, Layout: `docs/sphinx/` by module; Breathe directives only on reference pages; a program may feed more
  than one page and asserts every claim on each.
- `.claude/agents/docs-reviewer.md`: "flag a page that mixes them" becomes "each page has one main job within its
  module; an overview may explain and show code"; it also checks the voice rules.
- `ROADMAP.md`: "every module milestone ends with its chapter on the site" (no fixed list of pages until a second
  chapter shows the shape).
- The programs keep their file names (their test names come from them); their first-line comments name their new
  pages. `PHASE2_PLAN.md` and `PHASE3_PLAN.md` stay at `docs/sphinx/` (where `exclude_patterns` matches them) and are
  not updated: they are the record of their phases.

## Deferred

- A "how it works" page per module: when a reader reports a gap the overview leaves.
- A fixed chapter shape in `ROADMAP.md` and `CLAUDE.md`: after the second chapter.
- Figures (a bracket closing, Newton's tangents), generated by programs as the book's are.
- Hosting and versioning, `linkcheck`, the install section, redirects for paths changed after hosting.

## Appendix A: the sample, `roots/index.md`

Every number is asserted by an existing program or by the new `roots_overview.cpp` (the expanded cubic: measured
2026-10-29, Brent returns a point 4.2e-6 from 1 with `abs_error` 4.4e-16, the factored form within 1e-16; the computed
sign changes lie within about 5e-6 of 1). Phase 4b measured the noise width on a grid of step 1e-10 as 7.2e-6 (a
running error analysis bounds it by 8.2e-6), and the page says "about 7e-6".

````markdown
# Roots of functions

Many problems end in an equation that has to be solved numerically: the position of a planet in its orbit at a given
time, the volume a real gas fills at a given pressure from an accurate equation of state, the temperature at which a
liquid boils from a vapour-pressure correlation. Each can be written as

$$f(x) = 0$$

for some function $f$, and a solution, a *root* of $f$, is found by evaluating $f$ at a sequence of points that close in
on it. `nxx::roots` does this for a function of one real variable. You give it the function, as a lambda or any other
callable, and a place to start; it gives back the root and a bound on its error. [What you get
back](#what-you-get-back) says exactly what that bound covers.

## Two ways to start

The surest start is a *bracket*: an interval $[a, b]$ at whose ends $f$ has opposite signs. If $f$ is continuous on
$[a, b]$, it must cross zero somewhere inside, so the bracket holds a root. A *bracketing method* keeps such an
interval at every step and shrinks it, so it always ends at a sign change of $f$. The simplest, bisection, halves the
interval with every evaluation; Brent's method, the one to start with, fits a curve through the last few points and
jumps to where the curve crosses zero, bisecting instead whenever the jump would leave the bracket or the bracket is
not shrinking fast enough. To pin the root of $x^2 - 2$ in $[1, 2]$ to a few units in the last place of a `double`,
bisection needs 52 evaluations and Brent's method 9 ([costs](costs.md) compares every method):

```{literalinclude} ../examples/roots_overview.cpp      (Brent on x * x - 2.0 over [1, 2], and its output)
```

If the bracket holds several roots, you get one of them; if $f$ has a pole inside, as $\tan x$ does at $\pi/2$, the
sign change may be the pole, and the method returns that instead ([handle failures](handle-failures.md) shows how to
tell).

Often you know only roughly where the root is. Then start from a *guess* with an *open method*. Newton's method, which
also needs the derivative $f'$ (yours, or one that `nxx::deriv` computes), replaces $f$ by its tangent line at the
current point and steps to where that line crosses zero. Near a root where $f$ crosses zero at an angle this is very
fast; far from one it can wander. Numerixx's open methods therefore watch for the first step across which $f$ changes
sign and finish inside that bracket, so a root they return carries the same guarantee as a bracketing method's. Unlike
a bracketing method they can fail, when their steps never cross a sign change, and the sign change they cross first
need not be the root nearest your guess.

When you have neither a bracket nor a good guess, `search` looks outward from a rough guess until $f$ changes sign and
returns that bracket, which you then pass to `solve` ([find a bracket](find-a-bracket.md)). It cannot find a root where
$f$ touches zero without crossing it, as $x^2$ does.

Every method is called the same way, `solve<Method>(fn, start)`, where the start is a `Bracket` or a guess, and works
in any real type: `float`, `double`, or a multiprecision type ([Working with Numerixx](../basics.md)).

## What you get back

Every solve returns either a root or the reason this solve did not find one, in a `std::expected`. The root has two
parts: `value`, the root itself, and `abs_error`, a bound on its distance to a point where $f$, as your program
computes it, changes sign. With the default tolerances the bound is a few units in the last place: for a `double`, a
relative error of at most about $10^{-15}$.

The bound is about $f$ as computed. Rounding in $f$ moves its sign change away from the exact root, by about $f$'s
rounding error divided by its slope, and no method can see that. For an $f$ computed in a few steps and crossing zero
steeply, the shift is smaller than the bound: about $2 \times 10^{-16}$ for Kepler's equation in the
[tutorial](tutorial.md). Where the terms of $f$ cancel, or its slope at the root is small, the shift can be far larger.
Written out as $x^3 - 3x^2 + 3x - 1$, $(x - 1)^3$ changes sign, as computed, anywhere within about $5 \times 10^{-6}$
of 1, and Brent's method returns a point $4 \times 10^{-6}$ from 1 with a bound of $4 \times 10^{-16}$; written as
`(x - 1) * (x - 1) * (x - 1)`, it is exact.

## Where to go next

The [tutorial](tutorial.md) solves Kepler's equation step by step. [Choose a method](choose-a-method.md) says when a
method other than Brent's does better, and the other guide pages each solve one task.
````

## Appendix B: a passage before and after (the tolerances page)

Before (`how-to/roots-tolerances.md`, today):

> A root at exactly 0, where fn is exactly 0, costs nothing extra: a bracket that contains 0 in its interior is split
> there first. The cost falls on a sign change near 0 that fn does not evaluate to exactly 0. With the default
> `abs_tol`, the method approaches a tiny root down to its own magnitude, and a jump or a root that fn misses at 0 down
> to the type's smallest numbers. Set `abs_tol` to a distance from 0 that is close enough:

After (`roots/tolerances.md`):

> A method stops when the bracket around the root is narrow enough, and "narrow enough" is set by two tolerances added
> together: a relative one, `rel_tol`, which asks for a number of significant digits, and an absolute one, `abs_tol`,
> which asks for a distance. The method stops once the bracket is within $\text{abs\_tol} + \text{rel\_tol}\,|x|$ of
> the root, so the looser of the two decides. The defaults ask for all but the last digit or so of the type, whatever
> the size of the root, so you rarely need to change them.
>
> The exception is a root near zero. With the defaults, a method must first close in on the root's own scale before
> the relative tolerance can stop it. For the root of $x - 10^{-300}$ in $[0, 1]$, bisection halves its way down
> through a thousand powers of two: about a thousand evaluations, almost all of them spent getting near zero rather
> than on digits. If a distance of $10^{-12}$ is close enough for you, say so: [code block], and the same search takes
> 41. The result is then within $2 \times 10^{-12}$ of the root, which for a root at $10^{-300}$ says only that it is
> near zero, and the same $10^{-12}$ applies to every root: one at $10^{-9}$ keeps only about three digits. Choose
> `abs_tol` from the scale of the roots you expect.
>
> In a multiprecision type the walk is far longer, about $2 \times 10^9$ halvings in `cpp_bin_float_50`, long enough
> that a solve seems never to end: there, always set `abs_tol`. A function that jumps at 0, or misses a root at 0, is
> walked down to the type's smallest numbers, and Brent's method walks it too.
>
> ```{note}
> A root at exactly 0, where your function returns exactly 0, is found at once: a bracket that contains 0 is split
> there first.
> ```

## Reviews of this plan

- `simplicity-critic`: the move in its own commit (adopted, Phase 4a); the seven guide pages kept, not merged
  (adopted); no "how it works" page or fixed chapter shape (adopted); no deriv examples in `basics.md` before the deriv
  chapter (adopted); no `colon_fence` (adopted); Breathe directives only on reference pages, the reviewer agent's
  rule, the programs' comments, relative paths, anchors, the plan files' location (adopted). Not adopted: folding
  `costs.md` into choose-a-method (it describes, so it is reference), and flat sidebar captions (nested toctrees scale
  to nine modules).
- `docs-reviewer`: the sample's opening claimed a bound on the distance to the true root (fixed); the tolerances rewrite
  dropped the multiprecision walk and `abs_tol` at every magnitude, so rule 3 now keeps result-changing caveats in the
  text (fixed); a terminology rule, and `basics.md` saying what each module's result means (added); the sample's code
  block, its counts' links, the "default" (fixed); output quoted only from the program (rule 2); nested toctrees;
  `basics.md`'s own programs. Adopted, except one guide page (see above).
- `numerics-reviewer`: a bracket "guarantees" a root and "cannot lose" it, against poles and several roots (fixed);
  "a reason why there is none" (fixed); the conditioning note's "same order as the bound" (replaced with asserted
  figures, the expanded cubic measured); the walk's cause and `abs_tol`'s result (fixed); open methods can fail and
  need not find the nearest root (added); `search` returns the bracket to the caller (fixed); Brent's bisection
  condition and worst case, bisection's 52 as ends plus halvings, Newton needing f', "a few units in the last place",
  the intro's examples (fixed); `long double` as "the tested toolchains", "Numerixx itself never throws" (adopted).
  To check in 4b: a long compile-time walk against MSVC's constexpr step limit before `basics.md` promises
  compile-time use of every call.

## Phase 4b: what was done, and the reviews

Done as planned, with these additions:

- **Quoted output.** `examples/printed.hpp` sends a program's stdout to a file and reads it back; the program compares
  it with raw strings, which the page includes (`:start-after: <name>_output = R"(`). Output through libm
  (`std::sin`, `std::log`) is printed to fewer digits than it holds, or described with the bounds the program asserts
  (the tutorial).
- **`basics.md`** with `basics_types.cpp` and `basics_results.cpp`; the generic parts of other-types (`long double` by
  platform, compile time) moved there.
- **A long constant-expression solve** (the plan's open question): Brent's 1076 evaluations on a jump at 0, Bisection's
  1049 and Brent's 156 on (x - 1)^3 compile with cl, clang-cl, GCC and Clang at their default limits; `basics.md`
  names each compiler's limit flag instead of promising every call.

Reviews of the rewritten pages (`docs-reviewer`, `numerics-reviewer`, a docs-only reader), every finding reproduced
or checked, all adopted unless noted:

- False statements fixed: the derivative's error estimate described as always present (empty for the stencils);
  "the Boost.Multiprecision types are not literal types" (narrowed to the two types on the site); Halley's 6 calls per
  point split as 3 + 3 (it is 1 + 2 + 3; the program's comment too); "each chapter holds a tutorial and a guide" (only
  roots does); sqrt(2) "to the last digit" (one ulp off; now asserted within an ulp against a two-double reference);
  the search's points described as cumulative steps.
- Halley dropped as the method for a root of odd multiplicity (it only ties bisection at m = 3, worse above); the
  teaching-methods note names Ridders' win on x^20 - 1 among the value-only methods.
- The pole test called a heuristic, unreliable where fn cancels (reproduced: on the expanded cubic it happened to call
  none of 90 results a pole, but the three values are rounding noise there).
- Undefined behaviour in `roots_other_types.cpp`'s check (`r->value` read before `r` was tested), fixed.
- New assertions for numbers that were stated unchecked: the walk of about 2e9 halvings in `cpp_bin_float_50`
  (`min_exponent`), the iteration cap 212, the stencil step 6.06e-6, which Boost types use expression templates,
  Kepler's `abs_error` above the 2.3e-16 shift, the x^3 not_converged/abs_tol claim, a deduced pair of
  `cpp_dec_float_50` expressions rejected, a `float` fn with a `double` bracket 3.5e-8 off, van der Waals' single root
  at 280 K and 4 MPa (negative discriminant) and its critical temperature 304 K.
- Added from the reader's report: the tutorial shows its reference constants; `abs_error` at an exact zero; the
  residual's bound derived; the van der Waals three-root note; a using-declaration per function and printing a Boost
  number with `str()`; `derivative_of`'s copy of fn when counting calls; the FMA trigger and remedy; the one-sided
  derivative (`Forward<2>`) near a domain edge, asserted; the multiprecision includes on `basics.md`; how to add the
  library in CMake and the tested compilers.
- Voice: f in mathematics and `fn` in code where they were mixed; terms defined at first use on each page; text after
  code blocks that had none.
- Not adopted: quoting getting-started's output (its second line goes through `std::exp`); an install section (still
  deferred).


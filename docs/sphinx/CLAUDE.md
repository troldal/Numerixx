# The user manual

The rules, roles and process for the user manual (`docs/sphinx/`, built as a site and a PDF), decided by the user on
2026-11-09. The root `CLAUDE.md` describes the layout and the build targets; `docs/CLAUDE.md` says how a figure is
made; `docs/design.md`, Documentation, says what goes in the manual and what in the book.

## Reader and purpose

A developer who knows C++ and calculus, but not numerical analysis, and wants a correct result quickly. The manual
stands alone: it documents every module in full and never requires the book. It is a manual, not a textbook: easy to
understand, but brief.

## Ground truth

What a page says about the library must be true of it. In order of precedence:
1. `include/numerixx/*.hpp` (the code and its comments);
2. the tests in `tests/`;
3. `docs/design.md`, the module's `docs/<module>/design.md` and its phase documents.
Where they disagree, the higher one wins; report the disagreement to the user instead of choosing silently. A page
describes the library as it is: no promised features.

## Evidence rules

- **Code.** A page shows code only through `literalinclude` from a program in `examples/`, built and run as the test
  `numerixx_docs_example_<name>` on all five toolchains. Its `main` checks every claim the page makes about it, against
  a closed form or an independent value, never against Numerixx itself.
- **Output.** Output a page quotes is a raw string in the program, compared with what it printed
  (`examples/printed.hpp`). Output that differs between toolchains is described by the bounds the program asserts.
- **Numbers.** Every computed number on a page is asserted by its program, or linked to what asserts it. Exact
  constants (machine epsilon, a rule's weights) may be typed.
- **Figures.** Drawn from data a tested program checks (below, and `docs/CLAUDE.md`).
- **Citations.** Only sources already verified in a module's `design.md` or the book's `refs.bib`. A new source is
  checked against the work or a library catalogue (authors, title, year, venue, pages, DOI) before a page cites it,
  and added to the module's `design.md`.
- **Reference.** Generated from the headers; Breathe directives only on the reference pages.

## Structure

By module (`deriv/`, `roots/`, `quadrature/`), with `getting-started.md`, `basics.md` and `core.md` for what the
modules share. A module's chapter has an overview, a tutorial, guide pages and a reference page. Each page has one
job, stated in its first paragraph: the overview introduces the problem and the methods, the tutorial teaches by
doing, a guide page solves one task, the reference describes.

## Language

Easy to understand, and brief.

- **The theory cap.** Only the theory a reader needs to choose a method, a step or a tolerance, and to read the result
  and its error estimate. Proofs, derivations, convergence and error analysis and history belong to the book.
- **Budgets.** About 800 words for an overview, 1,000 for a tutorial, 600 for a guide page; the plan sets each page's.
  An overrun needs the editor's reason in the plan. Code blocks, output, captions and alt text do not count; notes
  and headings do, and an inline code span or formula counts as one word (Phase 9: three counting methods gave
  budgets three meanings).
- **Cut what the reader would not miss.** No sentence repeats another page; link to it instead. No preamble, no
  summary of what was just said.
- **Plain words.** One idea per sentence, at most about 25 words; active voice; second person. A term of numerical
  analysis only where the page needs it, defined in a sentence at its first use on each page (readers land anywhere).
- **Problem first.** A page opens with what the reader is trying to do, and what goes wrong without the feature.
- **Explain the code.** Before a block, what it does and why; after it, what it prints and what that means.
- **Edge cases out of the main line**, in a `{note}` or a link to the reference, unless ignoring them gives a wrong
  result or a call that does not end.
- **One name for one thing.** f in mathematics, `fn` in code and in prose about the parameter; mathematics typeset,
  code in code. "The default" only where C++ has one.
- **For screen and print.** The whole is "this manual"; a chapter in the PDF (`getting-started`, `basics`, `core`, a
  module's overview) is a "chapter", a tutorial or guide page a "section". Never "site", "page", "click" or "scroll"
  (`_ext/numerixx_breathe.py` fails the build on them, but "call site"); "the figure", never "below".
- An overview and a tutorial end with what to read next. British spelling (behaviour, organise, centre). ASCII
  sources.

## Figures

A plot or graph wherever it shows something words cannot, or makes a concept or a method easier to see, as in the
quadrature chapter (`PHASE8_PLAN.md`: a weighted sum, where the default method cuts, the tanh-sinh map, the rounding
floor, evaluations against digits, the infinite-interval maps). The plan lists each figure with what the reader should
see in it, its program and its script; the critic and the editor judge whether it earns its place. Each has a caption
that names what it shows, alt text that states what it shows with its numbers, and text that says what to see. It is
drawn by matplotlib at build time from data its program checks: the process, the model to follow and the checks are in
`docs/CLAUDE.md`. A figure replaces words, so a page with figures keeps its budget.

## Files and ownership

| Path | Owner |
|---|---|
| `PHASE<n>_PLAN.md` | editor |
| the phase's pages, `examples/<module>_*.cpp`, `plots/<module>_*.py` | author |
| `reviews/P<n>-<role>.md` | each reviewer (the main session saves the reports of agents that cannot write) |
| `reviews/P<n>-decision.md` | editor |
| `reviews/P<n>-response.md` | author |
| `conf.py`, `_ext/`, `requirements.txt`, `index.md`, `manual.md` | the main session, with the user's approval |

`CLAUDE.md`, `PHASE*.md` and `reviews/` are excluded from the build (`conf.py`).

## Roles

Custom agents in `.claude/agents/` (repository root); the main session coordinates them and talks to the user.

| Agent | Model, effort | Role |
|---|---|---|
| `docs-editor` | Fable, high | The phase plan; one voice and one set of terms across modules; budgets and the theory cap; the edit; the decision letter on each review round |
| `docs-author` | Opus, high | The example programs, the figures' programs and scripts, and the pages, to the approved plan |
| `docs-critic` | Sonnet, high | Argues for cuts: sentences, theory over the cap, repetition across pages, hard words, figures that do not earn their place |
| `docs-reviewer` | Opus, high | Every statement, number and example against the headers and the tests; the structure and the language rules |
| `numerics-reviewer` | Opus, high | Every mathematical and numerical claim |
| `docs-reader` | Opus, high | Learns the library only from the pages, does the plan's tasks, checks the results against independent values |

Plus the repository's independent review by Codex at the end of the phase.

## Process: stages, each ending at a checkpoint with the user

A phase rewrites or writes one chapter (or a set of shared pages). Commit to `dev` at the end of each stage the user
approves.

1. **Plan** (editor): `PHASE<n>_PLAN.md` with each page's job, budget and the theory it holds (and what it leaves to
   the book); the terms each page defines; each example program, with what it asserts and its independent reference;
   each figure, with what to see, its program and its data; the reader's tasks; the open questions for the user. The
   critic reviews the plan for scope; the user approves it. Nothing else is written before that.
2. **Evidence** (author): the example programs and the figures' programs and scripts, before the prose that uses
   them; each program passes with GCC. The user sees the figures.
3. **Draft** (author): the pages, to the plan; `numerixx_docs` builds with no warning and the examples pass.
4. **Edit** (editor): voice, budgets, terms, cross-references and consistency with the other modules; the PDF chapter
   read once (`numerixx_docs_pdf`).
5. **Review** (critic, `docs-reviewer`, `numerics-reviewer`, `docs-reader`, in parallel): `reviews/P<n>-<role>.md`.
   The editor's `reviews/P<n>-decision.md` adopts, adapts or rejects each finding with a reason; the user sees a
   summary and decides anything contested.
6. **Revision** (author): `reviews/P<n>-response.md`, each decision with what was done; then a second, targeted round
   by the reviewers whose findings were major (`reviews/P<n>-<role>-2.md`).
7. **Release**: all five toolchains, the site and the PDF; commit, push, CI; then the Codex review
   (`reviews/P<n>-codex.md`), its findings reproduced before they are acted on. A change to `include/` (a header
   comment) is proposed to the user first.

## Checks

- One example: `cmake --build --preset gcc --target numerixx_docs_example_<name>`, then
  `ctest --preset gcc -R numerixx_docs_example_<name>`.
- The site: `cmake --build --preset gcc --target numerixx_docs`; the PDF: `--target numerixx_docs_pdf`.
- Done means: no Doxygen or Sphinx warning, every check of `_ext/` passing, every example passing on all five
  toolchains, every page within its budget or with the editor's reason, and every decision of the letter done.

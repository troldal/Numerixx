# Milestone 3, Phase 6: one source, two manuals (online and PDF)

Status: done (2026-11-03), at the user's request: steps 1 to 3 and the independent review.

The pages in `docs/sphinx/` become the source of two manuals: the HTML site, as now, and a PDF built by Sphinx's LaTeX
builder and latexmk. Nothing is written twice. The pages stop speaking of a site, and the PDF gets its own root and
its own layout settings.

## What a trial build showed (2026-11-03)

Sphinx 9.2's `latex` builder, then `latexmk -pdf` (MiKTeX), on the current pages: 100 pages, no LaTeX error. MyST,
`dollarmath` (real LaTeX, better than MathJax), the code blocks and Breathe's reference all render; the checks in
`_ext/numerixx_breathe.py` pass under the LaTeX builder too. The problems:

1. **Structure.** `index.md` is the root: its two headings become chapters 1 and 2, and its toctrees, at the end of
   chapter 2, put every other page inside chapter 2. The table of contents has two chapters.
2. **Wording.** 35 lines say "this site", "this page", "guide pages" or "the page quotes", which reads wrongly in
   print.
3. **Blank pages.** The `manual` class opens every chapter on a right-hand page.
4. **Code width.** The examples are written to 120 columns (743 lines are over 95); at the default size the PDF wraps
   every line over about 95 columns, with a continuation mark.
5. **Overfull lines.** About six sources run past the margin: two paragraphs of `index.md` and the tables (the
   module table, for one).
6. **The book link.** `{download}` on `deriv/book.md` prints as plain code text in the PDF and leads nowhere.
7. **No target.** Only HTML is built; the title page has no version.

## The design

**Structure, verified in a spike.** A tag, `-t pdf`, given only to the PDF build, switches `conf.py`:

```python
if tags.has("pdf"):
    root_doc = "manual"
    exclude_patterns += ["index.md"]
else:
    exclude_patterns += ["manual.md", "introduction.md"]
```

- `manual.md` (PDF only): the title and one toctree, `introduction`, `getting-started`, `basics`, `core`,
  `roots/index`, `deriv/index`, `quadrature/index`. Each becomes a chapter.
- `introduction.md` (PDF only): `# Introduction` and a MyST `{include}` of `index.md` between two comment markers
  (`% intro-start`, `% intro-end`), so the home page's text is written once and its toctrees stay out.
- `index.md`: unchanged but for the two markers.

The spike gave seven chapters (Introduction, Getting started, Working with Numerixx, nxx::real and nxx::math, Roots of
functions, Numerical derivatives, Integrals of functions), the general index with the C++ names, and no new warning.
Rejected, also tried: `manual.md` with `index` in its toctree and `index.md`'s toctrees wrapped in `{only} html`.
Sphinx inlines the toctrees before it removes the `only` node, so every page vanished from the PDF.

**LaTeX settings in `conf.py`.** `latex_documents` (root `manual`, `numerixx.tex`, title "Numerixx", the `manual`
class); `latex_elements` with `openany,oneside` (no blank pages), the paper size, margins and a code font size chosen
so that a 120-column line fits; `latex_show_pagerefs = True` (a link prints its page, as "Results (page 4)": a printed
reader cannot click); `release` from CMake's `PROJECT_VERSION` for the title page. `{tabularcolumns}` on the tables
that run over.

**Wording: one vocabulary for both.** "This manual" for the whole. A page that is a chapter in the PDF (`basics`,
`core`, `getting-started`, a module's overview) is a "chapter", as `index.md` already calls a module's folder; a
tutorial or guide page is a "section". So:

| Now | Becomes |
|---|---|
| How this site is organised | How this manual is organised |
| on this site, that this site quotes | in this manual, that this manual quotes |
| This page shows / is about / explains | This section shows ... (in `basics.md`: This chapter explains) |
| a guide of pages that each solve one task | a guide whose sections each solve one task |
| the other guide pages each solve one task | the other sections of the guide each solve one task |
| the output the page quotes | the output the text quotes |
| at the end of the page; on this page | at the end of this section; in this section |

The full list is the 35 lines that `grep -niE "\bsite\b|\bpages?\b"` finds in the pages. The voice rules in
`PHASE4_PLAN.md` gain one rule: write for both media, so no "site", "page", "click" or "scroll". A check in
`_ext/numerixx_breathe.py` (a `source-read` hook, a few lines) warns on "this site", "on this site" and "this page",
so `-W` keeps the rule.

**The book.** `deriv/book.md` keeps `{download}` inside `{only} html`; the PDF's `{only} latex` says the book is the
file `numerixx-deriv-book.pdf` beside this manual, and the PDF target copies `docs/deriv-book/book.pdf` there.

**Build.** `numerixx_docs` runs one Sphinx pass more: after the HTML, `-t pdf -b latex` into
`build/<preset>/docs/public/latex`, with `-E -W`. It needs no TeX, so CI's docs job checks the PDF's source with
every warning and check, unchanged. A new target `numerixx_docs_pdf`, added only where latexmk is found (as
`numerixx_book`), depends on `numerixx_docs` and runs `latexmk -pdf` there, writing
`build/<preset>/docs/public/pdf/numerixx.pdf` and the book beside it. `numerixx_docs_all` adds it where it exists. Not
in CI, like the book's PDF: TeX is a large install.

**Project texts.** `CLAUDE.md` (Layout, Build and test), `docs/plan.md` (Milestone 3's description) and the
`docs-reviewer` and `docs-reader` agents say "documentation" or "manual" for the site, and name the PDF target; the
reviewer checks both outputs.

## Steps

One commit each, the user deciding the next.

1. **Wording.** The table above across the pages, the voice rule, the `source-read` check, the project texts.
   `docs-reviewer` on the changed sentences.
2. **The PDF build.** `manual.md`, `introduction.md`, the markers, the tag, `latex_documents`, the book's `only` pair,
   the CMake targets. The LaTeX builder drops the toctree captions ("Guide", "Reference"), so a module's guide has no
   visible edge in the PDF, while the overviews say "the other sections of the guide": make the grouping visible or
   reword (step 1's review).
3. **Layout.** Paper, margins, code font, `tabularcolumns`, page references, the version; then a read of the whole
   PDF (by eye, and `docs-reader` on the PDF), and Codex on the phase.

## Success criteria

- [x] `numerixx_docs` builds the site and the LaTeX source with no warning, locally and in CI (b57129c: run
  41689684905); the site changed beyond the wording only where step 2 and 3 say (the overviews' last sections, the
  deriv table's two rows, the concepts' requires-expressions, two stale sentences)
- [x] `numerixx_docs_pdf` builds the PDF (95 pages) with no LaTeX error; its contents list the seven chapters and the
  index; every link to a chapter or a section shows its page (145), a C++ name by design not (step 3)
- [x] No "site" or "page" wording left in the pages; the check fails a build that adds one (tried on a scratch copy)
- [x] No code line wraps (the check, and the widest lines read in the PDF); no overfull line over 5 pt in the latexmk
  log (one `\vbox` of 1.2 pt)
- [x] The book link works in the HTML and the PDF names the file beside it
- [x] The examples' tests pass on all five toolchains (run 41689684905; the three re-wrapped programs compare their
  output with the quoted text)
- [x] The reviews and Codex, with the adopted findings recorded here

## Step 1, done

The 35 lines reworded; voice rule 7 in `PHASE4_PLAN.md`; the `source-read` check, which fails the build on "page",
"site" (but "call site"), "click" or "scroll", tried on a scratch copy with added web-only words, one broken across a
line, and a "call site" (each reported with its line, the call site not); `CLAUDE.md`, `docs/plan.md` and the two docs
agents. `docs-reviewer`: nothing major; adopted the broader check, "this tutorial" and "All the output quoted in this
section" where "section" would read as a subsection, the step page's link to where its output is, "throughout this
manual", the aliases sentence and two re-wraps; the guide's edge in the PDF goes to step 2.

## Step 2, done

As designed, with the author's name in `latex_documents`. One finding of the build: each module overview's toctrees
sat in its last section, "Where to go next", so the PDF put the tutorial, the guide and the reference all inside it
(5.3, 5.3.1, ...). The overviews now end with three sections, Tutorial, Guide and Reference, each with its sentence
and an uncaptioned toctree: in the PDF, 5.3 Tutorial, 5.4 Guide with a numbered section per task, 5.5 Reference,
which also gives the guide the visible edge step 1's review asked for; on the site, the same three headings replace
the toctree captions, which the sidebar never showed. `numerixx_docs` runs the LaTeX pass after the HTML one;
`numerixx_docs_pdf` builds a 116-page PDF and puts `numerixx-deriv-book.pdf` beside it.

The user's check of the PDF's tables, before the commit: `tabulary` sizes a column by its average text, and code it
cannot break runs under the next column, unreported. All eleven tables read in the PDF:
- the quadrature overview's "Four methods" and the roots guide's method table clipped their first or second column
  (`GaussKronrod<1`, `nxx::roots::Bisectio`): fixed with `{tabularcolumns}` widths;
- the deriv guide's method table had two rows broken in both the site and the PDF since Phase 5: the `|x|` in
  `$\max(|x|, 1)$` split the cells. Now `\lvert x \rvert`; a cell count of every table row against its header finds
  no other;
- the roots costs table "From a bracket" (12 columns) clips the headers "Bracket", "Brent" and "Toms748", words TeX
  cannot hyphenate: it does not fit the line at the body's size and padding, so it goes to step 3 with the layout;
- the other seven are fine.

## Step 3, done

- **Layout** (`conf.py`): A4, 2 cm margins, `openany,oneside` (no blank pages: 116 pages became 95), code at 8 pt
  (`\footnotesize`), the release (2.0.0, from CMake) and the author on the title page, a smaller ragged index for
  names such as `std::formatter<nxx::quadrature::Failure<T>>::format`, and the default LaTeX chapter heading, since
  fncychap's Bjarne set "NXX::REAL" in capitals.
- **Page references** (`_ext/numerixx_latex.py`, a LaTeX writer): `latex_show_pagerefs` put "(page n)" after every C++
  name too, template parameters included ("T (page 87)"), which pushed the quadrature signatures 96 pt past the
  margin. Only links to a chapter or a section keep it, the same-page ones included (145 in all); a C++ name's anchor
  (`_CPPv...` or a Doxygen ID) or code text marks the others.
- **Code width.** 8 pt fits 115 columns: the user chose to re-wrap the three wider example lines (116, 116, 119) in
  `quadrature_method_of_your_own.cpp`, `quadrature_high_precision.cpp` and `deriv_one_sided.cpp`, outputs unchanged.
  Code inside an API description is at 7 pt (the user's choice, over re-wrapping `include/`) and fits 122 columns at
  its indent; the patch in `_ext/numerixx_breathe.py` puts a concept's requires-expression on its own line, which
  brings `estimating_method` from 128 to 87 columns. A check fails the build on a wider line in either place; it
  found the reference's code, which a check of the examples alone had missed.
- **Tables**: `{tabularcolumns}` on the introduction's table; the 12-column roots costs table at `\footnotesize` with
  3 pt padding.
- **Primes**: Sphinx writes f'' as two closing quotes, which TeX joined into a closing double quote; the writer keeps
  them apart.
- **Text**: "Getting started" named two module chapters and `basics.md` "both modules' headers" (both since Milestone
  4); the includes sentence now precedes its block.
- **The read of the whole PDF** (`docs-reader`): it works as a manual; the guides are clean. Its findings 1 to 5 are
  fixed above; the reference's text from the headers and the page-break polish (items 6 and 7) are deferred in
  `docs/plan.md`, at the user's choice.
- Overfull boxes: none but one `\vbox` of 1.2 pt.

## The independent review (Codex, on 3149643..b57129c)

No blocking defect; four edge cases in the extensions, each reproduced, none reached by the current pages:
1. The code-width check counted a tab as one column; LaTeX sets it to the next multiple of 8 (15 tabs and an "x":
   16 counted, 121 rendered). Adopted: tabs expanded before measuring.
2. A parsed-literal block in the reference would get the 122-column allowance without the smaller type, and a doctest
   block is neither checked nor resized. Not adopted: neither occurs, and MyST makes no doctest blocks.
3. A link whose text is code lost its page even when it points to a section (`` [`Results`](basics.md#results) ``,
   reproduced in a scratch build). Adopted: the target alone decides; the PDF keeps its 145 page references.
4. The "call site" exemption missed "call" and "site" on two lines, and the check read link destinations
   (`.../sites/...`). Adopted: "call" is looked for before "site" across whitespace, and destinations are blanked,
   keeping line numbers.

Found correct: the writer's context and groups against Sphinx 9.2 and Breathe 5.1, the `visit_Text` rewrite, the
root switch and the include between markers (no link into an excluded root, no duplicate label), the CMake order and
cleaning, the case without latexmk, the copied book, the wording, and all 224 code blocks within their limits.

## The user's decisions (2026-11-03)

Plan approved; step 1 first. A4 paper; the title page's author is Kenneth Troldal Balslev; the book's PDF goes beside
the manual's, as above.

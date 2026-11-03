# Milestone 3, Phase 7: plots in the manual

Status: done (2026-11-03, 20fe902, CI run 41716617342), at the user's request: one illustrative plot in the quadrature
chapter, on the site and in the PDF, drawn by matplotlib.

## The plot

In `quadrature/choose-a-method.md`, "Measure the cost on your own function": where the default (`GaussKronrod<15>`,
1350 evaluations) and `TanhSinh` (131) evaluate 1/sqrt x on [0, 1], as the number of evaluations in each decade of x
(the user's choice). The default cuts the panel at 0 again and again, about 100 evaluations in every decade down to
1e-16; `TanhSinh` spends one or two per decade, its smallest point near 1e-176.

## The data (the user's choice)

The manual's rule holds: every number comes from a program built and run as a test. A new program,
`examples/quadrature_sampling.cpp`, records each x the two methods evaluate, counts them per decade, prints the table as
CSV and compares it with the table it holds in a raw string, as the programs compare the output their pages quote; it
also asserts what the page says of it. The plot script, `docs/sphinx/plots/quadrature_sampling.py`, reads the table
from that raw string, so the docs build compiles nothing and CI's docs job is unchanged. A raw string rather than a CSV
file: the examples also run under node (Emscripten), where a file of the source tree cannot be read.

The decade of x is counted by multiplying by 10 until x >= 1: plain IEEE arithmetic, the same on every target, where
log10 from a libm could differ at a decade's edge.

## Rendering

- Sphinx: `matplotlib.sphinxext.plot_directive`, SVG on the site and PDF in the PDF manual (`plot_formats`), no source
  link or format links; MyST's `{plot}` directive on the page, with alt text and a caption.
- `requirements.txt`: matplotlib 3.11.2 and its dependencies, pinned, each released before 2026-10-20 (two weeks before
  this phase), checked to resolve with uv for Python 3.14 to exactly the pin set.
- One figure for both media: a transparent background, grey axes and text, and two series colours from the dataviz
  reference palette, blue `#3987e5` and orange `#d95926`, which pass its validator on both the light and the dark
  surface. Log scales on both axes; a legend and direct labels.

## As built

- The page calls the directive inside `{eval-rst}`: matplotlib's directive inserts reStructuredText into the
  document, which MyST's own directive machinery does not support (`insert_input`). `plot_basedir` is the manual's
  folder, so that the page names `plots/quadrature_sampling.py`: under MyST the directive does not resolve a path from
  the page's folder, and a leading / is the drive's root on Windows.
- The figure has a legend and a note for TanhSinh's two points beyond the axis, not direct labels.
- In print a figure stays where the text puts it (`figure_align` H), beside its paragraph; matplotlib's directive
  takes no `:name:`, so the text cannot cite it by number.

## The review (`docs-reviewer`)

Nothing major; the table and every claim checked independently (all 1481 points as hex floats, the decades in exact
rationals). Adopted: the alt text said "to 1" for "to 0.1"; "three or four cuts" became "about 3.3 cuts", asserted
through the 40 cuts of 30 evaluations, the k-th inside [0, 0.1 / 2^k]; the default's 135 in [0.1, 1) said and
asserted; the range check moved before `decade_of`, which would loop forever on 0; the legend's totals from the
table; a darker grey (`#757470`, 4.7:1 on white, 3.7:1 on the dark site); `decade_of`'s comment; "decade" defined;
"evaluations" throughout; the figure kept beside its paragraph in print; the caption in the page's notation. Not
adopted: reproducible SVG and PDF files (matplotlib writes the build's date into them), which nothing needs.

## Success criteria

- [x] The program passes on all five toolchains (CI, and locally)
- [x] `numerixx_docs` builds with no warning, locally and in CI; the plot shows on the site (SVG) and in the PDF
  (vector), read in both
- [x] The page's text about the plot is asserted by the program
- [x] `docs-reviewer` on the page, the program and the script

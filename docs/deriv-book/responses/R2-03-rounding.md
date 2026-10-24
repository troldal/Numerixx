# R2-03-rounding: responses of the fundamentals author

Chapter 3, `chapters/03-rounding.tex`, and the examples `code/03-rounding-noise-estimate.cpp` (E21) and
`code/03-rounding-step.cpp` (E3), revised against `reviews/R2-decision.md` (section 3, row `03-rounding.tex`, +9 and
-8 lines). Wrapped PDF: 9 pages before and after.

| Finding | Decision | What changed |
|---|---|---|
| T9 (engineer-02, numerics-02) | accept, with one correction | **Snippet.** `noise-estimate-use` now ends with the last step and its guard: `double eta = 1e-10; // the solver's tolerance, as fallback`, then `if (e.k > 0) { eta = e.sigma / std::abs(f(Re)); }`, then `double const h_eta = std::cbrt(eta) * Re; // the step`. The guard is inside the markers. **Example.** It now also differentiates at `h_eta` and asserts a relative error of at most 1e-4 against the implicit derivative; it measures 5.5e-7, against a model eta^(2/3) of 1.2e-7. **Text after the listing.** Three sentences: k = 0 and no value when no three orders agree, so retry with a tenth of the spacing or fall back to the solver's tolerance, as the listing does; a spacing on the scale of f itself can make the orders agree by accident at a wrong level, which the full procedure guards against with a sign test and a spacing adjustment, so use the default step, check k, and distrust a level near eps for a solver; and the header `docs/deriv-book/code/noise.hpp` is under the repository's MIT licence |
| library-05 | accept | `explicit-step` in `code/03-rounding-step.cpp`: the two-line `if (d)` is braced |

**The correction to T9.** The letter says "the spacing is so small that the smooth part has not dropped below the noise ... so retry with ten times the spacing". The direction is reversed. With the book's own `estimate_noise`, on sin(50x)(1 + 1e-8 u) at x = 1:
- h = 0.1 settles at k = 1 with the wrong level 2.2 |f|;
- h = 1e-2 and 1e-3 report k = 0;
- h = 1e-4 and the default step report 4.7e-9 and 7.8e-9, against eta/sqrt(3) = 5.8e-9.

On the clean exp(-x^2) at 2.5, h = 0.1 to 1e-3 report k = 0 and h = 1e-4 reports the rounding level. So k = 0 means that the spacing is too large for the smooth part to fall below the noise within eight points, and the retry is with a tenth of the spacing. Ten times the spacing leads into the wrong-level case. The text says so. The probe is `fund/probe.cpp` in the session's scratchpad.

**"Distrust a level near eps".** The same probe found the second failure on the bisection friction factor. At Re = 1e5 with h = 0.61 exactly (the default step is 0.6055), and at Re = 1e6 with the default step, the recipe settles at k = 3 with about 7e-17 to 9e-17 relative, the rounding level, although the solver's noise is about 1e-10. The quantised output seems to repeat its values at those eight points. At the default steps of Re = 2e4, 5e4, 1e5 and 2e5, and at multiples 0.25 to 10 of the default step at Re = 1e5, it gives 2.6e-11 to 6.1e-11. A level near eps is therefore implausible for an iterative solver, and the text tells the reader to distrust it. DS22's quoted values are unaffected.

**Paid for in the chapter.** Each change was compiled and checked against the 9-page limit:
- the "Ridders' method does not help here" paragraph opens with "Ridders' estimate assumes ...";
- 3.7 is three sentences: float, its measured error, the 50-digit type;
- the units paragraph of 3.5 is one sentence;
- "A smooth but inaccurate f is not noisy" is one sentence;
- the "Two cautions" paragraph is shorter by a line;
- checklist items 5 and 6 are one and two lines.

The editor had already trimmed the noise header's comments. The content of every trimmed passage is kept.

**Examples.** Built with g++ 16 and the book's flags, both exit 0.
- **`03-rounding-noise-estimate`:**
  - `noise of f: 4.3e-11 relative (k = 2), step 35.1`
  - `derivative at that step: relative error 5.5e-07`
  - `noisy exp: 6.2e-09 relative (k = 2), against eta / sqrt(3) = 5.8e-09`
- **`03-rounding-step`:** output unchanged.

`ctest --preset gcc -R "numerixx_book_0[123]"` passes 7/7.

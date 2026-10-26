# Language pass, stage B: practice author (Chapters 7 to 10, quick reference, Table 1.1 rows)

PLAN.md section 12.9, plus the approved factual fix D-22. Files rewritten: `chapters/07-using.tex`,
`chapters/08-api.tex`, `chapters/09-tables.tex`, `chapters/10-guarantees.tex`, `chapters/qr-quick-reference.tex`,
`chapters/shared-method-summary.tex`. The figure files `figures/07-using-{newton,kink}.tex` hold no captions; both
captions are in Chapter 7 and were rewritten there. Code, `data/` and every other author's files are untouched.
Nothing committed; the full book was not built.

## Checks

| File | `check_language_pass.py 3192e8a` | `check_chapter.sh` | Pages, standalone (before -> after) |
|---|---|---|---|
| `07-using.tex` | identical content sequences | COMPILED OK | 13 -> 14 |
| `08-api.tex` | identical content sequences; code list shows the four D-22 changes | COMPILED OK | 7 -> 8 |
| `09-tables.tex` | identical content sequences | COMPILED OK | 4 -> 4 |
| `10-guarantees.tex` | identical content sequences | COMPILED OK | 6 -> 7 |
| `qr-quick-reference.tex` | identical content sequences | COMPILED OK | 1 -> 1 |
| `shared-method-summary.tex` | identical content sequences | via Chapter 1 and the quick reference: COMPILED OK | Chapter 1 stays 7 |

`python tools/check.py`: 0 errors, 0 warnings. The two errors that the stage-A build reported were the D-22 names in
Chapter 8.

Word counts (the script's rough count): 07 4812 -> 5596; 08 2780 -> 3119; 09 1122 -> 1175; 10 3003 -> 3231;
quick reference 683 -> 728; shared rows 102 -> 122.

### Review lines, each read

- Margin notes. Chapter 7: 1 -> 8 (continuous integration, bisection, fixed-point iteration, concept, `constexpr`,
  fused multiply-add with extended precision, spline). Chapter 8: 0 -> 3 (lvalue and rvalue, aggregate, translation
  unit). Each stands on a page without a margin figure or side caption (checked in the standalone PDFs). Chapters 9
  and 10 and the quick reference gain none: their glosses are in the text (dyadic, doctest executable, build tests,
  regression baseline, flagged, non-deduced, completing a type, `x == x`, system header, radix, `\fmax`, saturate,
  span, stationary point, piecewise, pole, optimistic, interpolant, dangling).
- Inline maths added, all spelled-out steps or glosses: 07 `h^2` and `1/h` (why halving the step quarters or doubles
  the error), `H^\ast` and `\enth(T)` (what the residual is), `\lambda` (no closed-form solution for it),
  `\noise = \epsmach^{2/3}` (the noisy-function recipe behind the step of a mixed partial), `ab + c` (the
  fused-multiply-add gloss), `f' = 0` (stationary point), `x` (resolution of `x` as the spacing of the numbers),
  `x = g(x)` and `x \leftarrow g(x)` (fixed-point gloss), `f` (a repeated reference). 08 `\den` (named as the common
  denominator), `f'` and `f''''` ("for f' to f''''" after the Ridders caps), `x` (ulp gloss). 10 `f` and `x`
  (repeated references). No inline maths dropped anywhere.
- Code added: 07 `Central<2>` (what the comparison with `Central<4>` measures), `atol` and `rtol` (defined, as in
  Chapter 1), `constexpr` (gloss), `x` (deduction). 08 `T`, `fn`, `h`, `int` (non-deduction and `valid_order`
  spelled out). 10 `fn` (repeated reference). Quick reference `Forward` and `Backward` (which way each never
  samples). No code dropped except the D-22 names below.
- Prose numbers: none added, none dropped, in every file.

## Two representative sentences

1. Chapter 7, 7.2. Before: "A step on the wrong scale (...), a function less accurate than $\epsmach$ (...) or a
   break within a step of $x$ (...) give a wrong value without a failure; \cref{sec:using-check} shows how to check
   one." After: "Each of these gives a wrong value without a failure: a step on the wrong scale (...), a function less
   accurate than $\epsmach$ (...), and a break within a step of $x$ (...). \Cref{sec:using-check} shows how to check
   such a value."
2. Chapter 8, 8.2 Failures. Before: "For the extrapolation methods, specifically: \code{fn(x)} is not finite where it
   is needed (an even N with \code{Ridders}, every N with the one-sided forms; after one evaluation), or no tableau
   reached two finite levels." After: "The extrapolation methods have two failures of their own. In the first,
   \code{fn(x)} is not finite where it is needed: for an even N with \code{Ridders}, and for every N with the
   one-sided forms. This failure comes after one evaluation. In the second, no tableau reached two finite levels."

## D-22, the factual fix (separate from the language)

`include/numerixx/deriv.hpp` (lines 1740, 1771, 1794, 1815) now requires `nxx::detail::callable<F, T>`; the concept
is in `nxx::detail` in `math.hpp` (lines 52 to 71). Changed in `chapters/08-api.tex`, each a separate edit made before
the language edits:
- line 53, the signature block of the overload with a step: `&& detail::callable<F, T>` -> `&& nxx::detail::callable<F, T>`;
- line 59, the signature block of the overload with the default step: the same;
- line 83, 8.2 Template parameters, item F: `detail::callable<F, T>` -> `nxx::detail::callable<F, T>`;
- line 311, 8.8 Rejected calls: `detail::callable` -> `nxx::detail::callable`.

`check_language_pass.py` lists exactly these four as code dropped and added; its strict sequences are identical. I
checked the rest of what Chapter 8 transcribes against the header: `default_step`, `Estimate`, `error`, `order`, the
method tags, `real`, `abs`, `isfinite` and `root_epsilon` are unchanged since 3192e8a. Only the spelling of the
constraint changed in the signatures. No other chapter of mine names the concept. None of my chapters holds a commit
or version reference, so there was nothing to refresh there.

## Compromises and judgement calls (for the editor)

- Chapter 9, 9.5: the last paragraph reads "With the evaluations per level, it gives the maximum evaluations in the
  last row" rather than the fuller wording. With the fuller wording, Table 9.10 missed page 4 by one line and took a
  page of its own; with this wording the chapter is four pages again. Section 9.4 also stays closer to the original
  sentence shape for the same reason.
- Chapter 10, 10.6: Table 10.4 now stands at the head of the section in the source, with a comment. Below the longer
  paragraph it no longer fitted on page 6, and it left a third of a page empty. The environment and label sequences
  are unchanged.
- Table 7.1: the cells stay terse, as table cells. Only three were reworded ("covers the error in", "accurate only
  to", "every method that samples across the break").
- Chapter 7, 7.7: "How the library avoids that arithmetic" had an unclear referent. It now reads "the arithmetic that
  overflows", which is what the sidebar it points to is about. The same paragraph gains "it sets no option for them"
  after "The build pins neither". I checked that `CMakeLists.txt`, `CMakePresets.json`, `tests/CMakeLists.txt` and
  `docs/deriv-book/CMakeLists.txt` contain no contraction or precision flag.
- Chapter 7, kink paragraph: "At the break itself they return that side's derivative if $f(x)$ belongs to it" is now
  "each of these returns the derivative of its own side, if $f(x)$ belongs to that side". This keeps the condition.
  A shorter reading would have dropped it.
- Chapter 7, Newton: "Each iteration costs three evaluations" now adds "two for the slope, and one for the residual".
  I checked this against `code/07-using-newton.cpp`.
- Chapter 8, The step: "For a stencil, \code{h} well above $\ulp(x)$" (a fragment) is now "\code{h} needs to be well
  above", which is the header's word. I used "needs", not "must", since nothing enforces it.
- Chapter 10, 10.1: "(\code{BackwardRidders} is the mirror of \code{ForwardRidders})" is now "\code{BackwardRidders}
  is not among them, since it is the mirror". I checked this against `data/table-accuracy-matrix.tex`, which has seven
  rows and no `BackwardRidders`.
- `shared-method-summary.tex` is marked editor-owned in its own comment. I rewrote it because the brief asked me to.
  Both of its tables (Table 1.1 and the quick reference) still set cleanly. The quick reference stays one page.

## Noticed outside the pass (not changed)

- Chapter 10, 10.6 and Table 10.4 say "The library's tests are one doctest executable and four build tests" and list
  the files in `tests/`. Since roots Phase 1 (02fffe5), the same executable also holds `test_roots_bracketing.cpp`,
  which the table does not list. Should the table say "the tests of `nxx::deriv`"? That is a content change, so it
  is the editor's or the author's call.
- After D-22, `math.hpp` also declares `nxx::detail::callable`. Table 8.1's caption ("math.hpp those in nxx and
  nxx::math") still holds for the public names, and 8.1 says that everything in a `detail` namespace is internal. So
  I did not add where the concept is declared.

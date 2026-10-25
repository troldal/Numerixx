# P5-08-api: the Phase 5 update (PLAN.md 12.8)

Author: `book-author-practice` (Opus 5.5). File: `chapters/08-api.tex`. `tools/check_chapter.sh`: `COMPILED OK`,
7 pages standalone (as before; about 21 lines longer); `tools/check.py`, with its whole-declaration check of the
signature blocks against `deriv.hpp`: 0 errors, 0 warnings. Example `code/08-api-usage.cpp` unchanged, builds and
passes (`ctest --preset gcc -R numerixx_book_08`).

| Item of 12.8 | What changed |
|---|---|
| 8.1, Table 8.1 | The row `Ridders` is now "`Ridders`, `ForwardRidders`, `BackwardRidders`: classes, empty method tags" |
| 8.2, template parameter M | Lists the two new tags, and names the two groups as the header's concepts do: the stencils and the extrapolation methods |
| 8.2, parameter h | "the initial step of an extrapolation method: `Ridders` or one of its one-sided forms" |
| 8.2, default steps | Adds "for `ForwardRidders` and `BackwardRidders` it is max(\|x\|, 1)/32 for every order" (header, `default_step`) |
| 8.2, result | `abs_error` "filled by the three extrapolation methods" |
| 8.2, failures | fn(x) not finite "where it is needed (an even N with `Ridders`, every N with the one-sided forms; after one evaluation)"; a non-finite level restarts "at half that step", an untrusted result "at half the smallest step it tried" (the header's restart, `deriv.hpp` lines 1511-1517 and 1576-1583) |
| 8.2, cost | The one-sided methods: fn(x) once and N per level, never on the other side of x, at most `\valForwardRiddersMaxEvals{First..Fourth}` (33, 65, 97, 129) |
| 8.3, `Estimate` | `abs_error` of the extrapolation methods; for the one-sided forms the rounding floor "as the tableau amplifies it" (the header's ForwardRidders comment, "covers the rounding that the tableau amplifies") |
| 8.5, signature block | `struct ForwardRidders {};` and `struct BackwardRidders {};` transcribed; `check.py` finds both in the header |
| 8.5, Table 8.2 | Two rows: `ForwardRidders`, `Forward<1>` of order N at h_i = h/2^i, nodes x + k h_i, k = 0..N, fn(x) once plus N per level, at most 33..129 (macros), default step s/32, estimate yes; `BackwardRidders`, mirrored, x - k h_i, the same. The table's `\tabcolsep` is 4 pt and its first column 26.5 mm so that the tag names fit without a break (the total width is unchanged within the 114 mm text block) |
| 8.5, `exact_weights` | "The extrapolation methods are checked against the table they extrapolate: `Ridders` against `Central<2>`, the one-sided forms against the `Forward<1>` table of order N" (PHASE5_API.md, `detail::base`) |
| 8.5, the algorithm | "`ForwardRidders` the `Forward<1>` rule of order N, and `BackwardRidders` runs `ForwardRidders` on f(-t) at -x (`sec:ridders-onesided`)" |
| 8.6, the warning | Scope: "A derivative above the first, by any method, in a type whose machine epsilon is above 1e-7"; the message quoted exactly from `deriv.hpp` lines 1810-1811 and 1831-1832 |
| 8.6, the numbers | The single point x = 1 for f'''' replaced by the DS6 float medians on e^x for f' to f'''' and three methods: `\valFloat{First..Fourth}CentralRelError`, `\valFloat{First..Fourth}RiddersRelError`, `\valFloat{First..Fourth}ForwardRiddersRelError`, over `\valAccuracyMatrixPoints` points |
| unrelated | 8.8: a stray blank line split the last sentence into two paragraphs ("A callable that cannot" / "be called"); joined |

Not mine, for the methods author: the quoted snippet `float-fourth` (`code/05-orders-float.cpp`) still says in its
comment that the warning "says the result is mostly rounding noise"; the message now says "keep few correct
digits". 8.6 quotes that snippet.

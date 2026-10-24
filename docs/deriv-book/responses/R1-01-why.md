# R1-01-why: responses of the fundamentals author

Chapter 1, `chapters/01-why.tex`, revised against `reviews/R1-decision.md` (section 3, row `01-why.tex`). Wrapped
PDF (`check_chapter.sh`): 5 pages before, about 5.3 now (budget +0.5).

| Finding | Decision | What changed |
|---|---|---|
| T11 (developer-07, library-05; developer-03 for 1.3) | accept | 1.3, the bullet "No exceptions from the library": a failure is a value; the library throws nothing itself; an exception thrown by `fn` passes through unchanged, and `derivative` is not `noexcept`. The `.value()` trap itself belongs to 7.2, 8.3 and the quick reference (practice author) |
| T4 (engineer-04, developer-11) | modified | 1.4 now opens with a complete program, the new example E1b `code/01-why-program.cpp` (snippet `program`, the whole file: the includes, `main`, the lambda `H`, the call, `if (d) ... else ...`, `return d ? 0 : 1;`), quoted as the first listing. One paragraph says that the first include is the book's own header and that `book::` is the namespace of the book's examples, before the `shomate` and `shomate-H` snippets. "Getting the library" is 7.1 (practice author) |
| developer-15 | accept | The first listing has an `else` that reports the failure; the paragraph after it says what `if (d)` tests, that only then `d->value` is read, and that the `else` branch handles the failure |
| T1 (engineer-01, developer-02, library-02; engineer-09) | accept | 1.5's "The one check" is now two paragraphs. "Did it work?": `if (d)`, then `d->value`, with an `else`. "Is it right?": the three silent failures in one line each (scale, 3.5; solver output or noise, 3.6; kink or range break, 7.8), the repeat at 2h and h/2 (7.5), and for Ridders the tolerance test `*r->abs_error <= tol * std::abs(r->value)` before `r->value` is used (7.2) |
| engineer-01 (caption) | accept | Table 1.1's caption glosses w ("the distance in x over which f changes appreciably", 3.5) and eta ("how many digits of f you can trust, for example 1e-8 for a solver converged to a relative tolerance of 1e-8 on its output", 3.6). The rows are the editor's shared file `chapters/shared-method-summary.tex`, unchanged |
| engineer-09 (house idiom) | modified | Part I's idiom in Chapter 1: `d->value` and `*r->abs_error`, with a margin gloss at the first `*r->abs_error` ("`abs_error` is a `std::optional`, empty for the stencils; the `*` in `*r->abs_error` reads it"). The implicit conversion does not appear in Chapter 1 |
| library-16 | modified | The two-line `if (r)` of the `first-ridders` snippet (`code/01-why-first.cpp`) is braced; the program's `if`/`else` is braced |
| evidence-07 | accept | Sidebar 1.1: `squire1998` after the complex-step formula; "closely related to forward-mode automatic differentiation" cites `martins2003`; the history reads "The complex step is a special case of Lyness and Moler's earlier contour-integral methods, which also give higher derivatives of analytic functions [lyness1967]; Squire and Trapp brought the simple first-derivative form back into use" |
| evidence-19 | modified | 1.4: "at 1000 K the tables give `\valCpTableValue` J/(mol K) [chase1998, Part II, p. 1621], from which the fit differs by `\valCpFitVsTableRelDiff` relative"; the typed "about four significant figures" is gone |
| T13 (numerics-17h) | accept | 1.5 quotes the headline digits only as macros (`\valHeadlineFirstDigits`, `\valHeadlineFourthDigits`), now rounded (11 and 5); no digit count is written in words in the chapter |

Examples: `code/01-why-program.cpp` (new; prints `cp = 32.6917299992 J/(mol K)`, the value of
`\valCpDefaultValueLong`, exit 0); `code/01-why-first.cpp` (brace only; output and exit status unchanged). Both build
with the book's warning flags under GCC 16 and pass `ctest --preset gcc -R numerixx_book_01`.

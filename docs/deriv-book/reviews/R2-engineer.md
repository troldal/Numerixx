# R2-engineer

Reviewer: book-reader-engineer, round 2 (targeted). I read `reviews/R1-decision.md`, the responses that cite my findings, and the revised build (`build/gcc/docs/deriv-book/pdf/main.pdf`, 88 pages, commit 5432abf). I re-read §1.4–1.5, Chapter 3, §7.1–7.9 and the quick reference, and skimmed the rest for regressions. Page numbers are the printed ones.

(The reviewer agent has no file-writing tool; the main session saved this review verbatim from the reviewer's report.)

## Status of the round-1 findings

| Id | Status | Comment |
|---|---|---|
| R1-engineer-01 | resolved | §1.5 now has "Did it work?" and "Is it right?". The three silent failures are listed with pointers, plus the 2h/h/2 repeat and the Ridders tolerance test. Table 1.1's caption explains w and η in plain words. One layout issue remains: R2-engineer-03. |
| R1-engineer-02 | resolved | §7.4.2 "The Jacobian of a set of balances" (p. 46–47) uses an energy balance and an ideal-gas pressure in T and n. It states the cost of the library's way (m × n calls) and prints the hand-applied central difference with its caveats. This is what I needed. The by-calls variant is only referenced, which is acceptable. |
| R1-engineer-03 | partly | §3.6 (p. 19–21) now warns that a residual tolerance is not η and prints a compiled `estimate_noise` run on the friction factor (4.3 × 10⁻¹¹ against the actual 4.9 × 10⁻¹¹). §7.4.3 points to it. What is still missing is in R2-engineer-02. |
| R1-engineer-04 | resolved | §1.4 is now a complete program with `#include`, `main` and `else`, and `book::` is explained. §7.1 "Getting the library" gives a FetchContent snippet and names the tested toolchains with versions. The CMake pin is noted under 25. |
| R1-engineer-05 | resolved | §3.3 gives the digits lost and the formation-enthalpy order of magnitude. Checklist question 4 in §3.8 covers it. |
| R1-engineer-06 | resolved | w is the true scale and s is "the library's guess for w" throughout Chapter 3 and in the notation table. Chapter 3 now reads cleanly. |
| R1-engineer-07 | resolved | The storage and test details moved into Sidebar 4.1. §4.3 keeps the practical statement. Sidebar 3.1 no longer uses d. |
| R1-engineer-08 | resolved | §7.4.3 (p. 48) leads with "A solver can also be smooth but wrong by its tolerance". It says the estimate did not cover the error and gives the three actions. |
| R1-engineer-09 | resolved | §7.2 states the house idiom (`d->value`, `*r->abs_error`). §1.4 has the margin gloss. §7.5 reads `->value`. |
| R1-engineer-10 | resolved | §7.9 covers the nested call for a mixed partial, its cost, the step and the rule of thumb (about ε^{4/9}), plus caching. A two-line snippet would help, but this is enough to act on. |
| R1-engineer-11 | resolved | §3.5 has the K against °C example: about 15 times too small a step, about one digit lost. |
| R1-engineer-12 | resolved | §6.4, §7.4.3 and §7.5 now lead with their conclusions. "Digits of \|f\|/s" is gone. |
| R1-engineer-13 | resolved | §3.5 and §7.9 now agree on splines. |
| R1-engineer-14 | resolved | "Monadic" is replaced. The C++ terms are in the glossary. |
| R1-engineer-15 | resolved | §3.5 has the margin note on `root_epsilon`. The quick reference has "Your own step" with both lines of code. |
| R1-engineer-16 | resolved | §7.6 gives the two Boost rules in the main text. §3.1 points to Sidebar 3.3. |
| R1-engineer-17 | resolved | §5.4 is prose only; the listing is in §8.6. |
| R1-engineer-18 | resolved | Cross-references now read "Appendix A" and "Appendix B". |
| R1-engineer-19 | partly | The header phrases are gone. I accept the reason for keeping the data paths in the captions. But the main text of §5.3 (p. 30, `05-orders.tex` line 146) now reads "yet in the design record's prototype (`docs/PHASE4_API.md`)": a repository document cited in Part I main text, which was not there before. Repository citations in Sidebars 3.2, 5.2 and 6.2 are fine. |
| R1-engineer-20 | resolved | The "example checks" sentences are gone from §2.3, §4.4 and §5.6. |
| R1-engineer-21 | resolved | §2.2 has the margin gloss on O(h²). |
| R1-engineer-22 | resolved | §6.4 now reads "Ridders returns 100.000 000 000, accurate to 5.0 × 10⁻¹⁴ relative". |
| R1-engineer-23 | resolved | The glossary entries are added. |
| R1-engineer-24 | resolved | Fig. 6.4 (p. 37) now shows x beside k, with labelled regions and one vertical line dropped. Its labels are now very small at margin width, but legible on screen. |
| R1-engineer-25 | partly | "UNVERIFIED" is gone from the references. The version TODO (p. ii) remains, as planned for release. Related: the §7.1 CMake snippet pins `GIT_TAG aeec7e6`, a commit older than this revision. Tie the pin to the release commit named in "The version", so readers do not copy a stale pin. |

## New findings

### R2-engineer-01: the Ridders acceptance check is relative only and gives no guide to `tol`
- **Location:**
  - `chapters/01-why.tex` §1.5 (p. 7)
  - `chapters/07-using.tex` §7.2, the `accept` listing (p. 43–44)
  - the quick reference, "Two checks" (p. 81)
- **Severity:** minor
- **Problem:** the book now tells me to "use `r->value` only if `*r->abs_error <= tol * std::abs(r->value)` for the relative tolerance `tol` you need". This has two gaps.
  - **No guide to `tol`.** The example uses `1e-8` without saying why, and nothing says how to choose it. For instance: "for a result reported to n significant figures, `tol = 10^-n`", or "a value inside a Newton iteration needs no check".
  - **The test rejects good results near zero.** It is purely relative, so it rejects a perfectly good derivative near f' = 0. The §7.2 listing itself differentiates cos x at 0, where the exact answer is 0. That contradicts §7.5's advice: "Where f' is near zero, judge the absolute error".
- **Suggested change:**
  - Add one sentence on choosing `tol`.
  - Either make the check `*r->abs_error <= tol * std::max(std::abs(r->value), scale)`, with the scale named (for example |f|/max(|x|, 1), as §7.5 suggests), or add "near f' = 0, compare with an absolute tolerance instead (Section 7.5)" in §1.5, §7.2 and the quick reference.

### R2-engineer-02: the measured-noise recipe stops one step short and does not say where the code is
- **Location:** `chapters/03-rounding.tex` lines 335–345 (§3.6, p. 20–21); §7.4.3 (p. 48).
- **Severity:** minor
- **Problem:** after `estimate_noise` reports 4.3 × 10⁻¹¹, the text compares it with the actual noise and moves on. It never shows the last step, the one I need: `eta = e.sigma / std::abs(f(Re))`, then `h = std::cbrt(eta) * Re`. Nor does it say whether to round η up for safety, or what to do when it returns NaN ("none settled").
  - `estimate_noise` is "kept in the book's example header", but the book never says where that header is in the repository or whether I may copy it into my own code. Unlike `nxx::deriv`, this is a function I must copy.
- **Suggested change:**
  - Add the two lines that turn `e.sigma` into the step, and one sentence for the NaN case: "spacing too small or too large; retry with 10 times the spacing, or fall back to the solver's tolerance".
  - Give the repository path of the example header and state its licence (the repository's MIT licence, if that applies).

### R2-engineer-03: Table 1.1 splits the "Is it right?" list
- **Location:** `chapters/01-why.tex` §1.5 (p. 6–7).
- **Severity:** editorial
- **Problem:** the float lands between the second and third bullets of "Is it right?", and the third bullet (kinks) is pushed to p. 7 after the table. The three silent failures, the most important half-page of the chapter, read as two plus an afterthought.
- **Suggested change:** place Table 1.1 before "Did it work?" or after the list, or keep the list together with `\needspace`.

### R2-engineer-04: a repository document cited in the Part I main text
- **Location:** `chapters/05-orders.tex` line 146 (§5.3, p. 30).
- **Severity:** editorial
- **Problem:** "yet in the design record's prototype (`docs/PHASE4_API.md`), on test functions that vary faster than s, half the default step was better". This is the kind of internal reference R1-engineer-19 asked to keep out of the main text. Sidebar 3.2 already carries this citation.
- **Suggested change:** in §5.3 write "yet on test functions that vary faster than s, half the default step was better for every order (Sidebar 3.2)", and leave the repository path to the sidebar.

## Verdict

The revision fixes the four things that stopped me in round 1:
- a complete first program and how to get the library into my build;
- a genuine Jacobian of a set of balances, with the library's way and the cheaper hand-written way side by side;
- a measured-η recipe for solver outputs;
- a §1.5 that names the silent failures.

The w/s split makes Chapter 3 read cleanly, and moving the developer detail into sidebars makes Chapters 4 and 5 easier for me. The four new findings are small. With R2-engineer-01 and -02 fixed, I would keep this book at my desk and hand it to new engineers in my group.

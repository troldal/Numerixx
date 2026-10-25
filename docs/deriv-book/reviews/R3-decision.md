# R3-decision: the targeted review of the Phase 5 update (round 3)

Editor, 2026-10-25. Reviews: `R3-numerics.md` (10 findings), `R3-evidence.md` (9), `R3-codex.md` (4, OpenAI Codex,
saved by the main session). The main session reproduced the two majors and three Codex items before this letter
(the poles, the f''' and f'''' edge limit, codex-02 and codex-04). The editor applied every accepted fix in the files
concerned, the authors' files included, as in the release pass (PLAN 12.8, P5-15 to P5-19).

**Counts.** 23 findings: 19 accepted, 2 accepted in a modified form, 2 deferred, 0 rejected.

## Themes, decided once

**T1. Poles and the edge limit** (numerics-02, evidence-02, codex-01). Accepted. The non-finite edge limit stands as
stated. The tests' poles at 1e-3 and 1e-4 lie far inside it, so they show only that the method now recovers where it
did not; no limit for poles is tagged as measured. Closer than the limits, a pole does not give `non_finite`: the
method returns a wrong value.

New generated evidence (`r3_poles`, a sweep of Ridders on 1/(t - a) at 0):
- a = +-c 2^-k for c = 1.25, 1.5, 1.75 and k = 4 to 100, f' to f''''; these are never sample points, so the pole
  always lies between samples;
- 2328 results; none `non_finite`, and all 1586 uncovered results flagged (4 abs_error >= |value|).

The prose is rewritten in 6.4, 7.8 (at equal length), 10.4 and the quick reference: test `abs_error` near a pole.

**T2. The f''' and f'''' limit** (numerics-02). Accepted. Derived, the limit is 2^-34 s for f', 2^-31 s for f'' (its
rule reaches h) and 2^-30 s for f''' and f'''' (their rules reach 2h).

Measured on log x in half halvings (`r3_edge_orders`): a value is returned down to 2^-31.5 for f'' and 2^-30.5 for
f''' and f'''' (`\valRiddersEdgeLimitExponent{Second,Third,Fourth}`).

Corrected in 6.4 and 10.4. The quick reference names the f' figure only. This is a library documentation note, U-18.

**T3. F19 against the flag** (numerics-04, evidence-03, evidence-07). Accepted.
- DS23 gains `fr<N>_estimate_value` (abs_error / |value|). F19 plots it dashed, with a thin line at 1/4, so "flagged"
  is visible.
- The caption and 6.7 and 10.4 say "to the end of the study, k = 45". The generator now fails if a flagged run ends
  earlier.

## Findings

| Finding | Decision | What changed |
|---|---|---|
| numerics-01 (major) one-sided f'''' with a pole ahead | Accepted | `r3_poles` (b): `ForwardRidders`, a pole (1 + u) 2^-e ahead of 0 and the same pole behind, e = 4 to 30, 4000 draws per order. f' and f'' always covered; f''' 3 uncovered (1 neither covered nor flagged, 1.02 times the estimate); f'''' 10 uncovered, all 10 neither covered nor flagged, 57 to 65 times off, up to 5.0 times the estimate; behind, none. 10.2 carries the measured sentence with the advice to repeat at a second initial step; 6.7 has one clause pointing to it (its page is full) |
| numerics-02 (major) | Accepted: T1, T2 | 6.4, 10.4, 7.8, the quick reference; U-18 |
| numerics-03 the silent band | Modified | `r3_band`: `ForwardRidders` f' on log x, sqrt x and cbrt x at four points per halving, k = 2 to 45: 12 of 528 results neither covered nor flagged (log 1, sqrt 7, cbrt 4), each within 2.2 times its estimate. Quoted in 10.4. 6.7 and F19 keep the log reading, which they present as log's; no room in Chapter 6 for more |
| numerics-04 F19's frame | Accepted: T3 | F19 and its caption |
| numerics-05 B.3's "so" | Accepted | B.3 now says it at the same length: the bound to first order, and the alternation is what makes the recursion exact |
| numerics-06 "since the rules reach N h" | Accepted | 6.7: "chosen on a prototype sweep (the design record)" |
| numerics-07 "until the levels run out" | Accepted | 6.7: "beyond its reach the levels run out first" |
| numerics-08 untagged statements | Modified | 10.1: "(measured)" after the factors. The derived amplification in 6.7 is not added (length) |
| numerics-09 Figure 5.2's box | Accepted | An accent outline that steps up between N = 2 and 4, leaving out the double point at N = 4 and clearing the N = 4 markers |
| numerics-10 "relative error 0" | Accepted | Generated wording `\valAliasingResult` ("the double nearest to 64 pi" when the relative error is 0) in Sidebar 6.3 |
| evidence-01 (major) the version pin | Deferred | The release step: `00-frontmatter.tex` and `code/getting-started.cmake` pinned together, at or after f20696f |
| evidence-02 (major) | Accepted: T1 | as T1 |
| evidence-03 | Accepted: T3 | as T3 |
| evidence-04 the AI statement | Accepted | One sentence: updated for Phase 5 in the same way, drafted and reviewed by Claude, with an independent review by OpenAI Codex. Cline is not named for the update |
| evidence-05 the every-power case uncited | Accepted | Burden et al., Section 4.2, in the prose of 6.7 and B.3, as the main session confirmed the section treats an expansion in every power of h |
| evidence-06 the staged `book.pdf` | Deferred | Not the book's sources; the main session raises it with the author |
| evidence-07 "to k = 45" | Accepted: T3 | as T3 |
| evidence-08 a worst ratio of "1.0" | Accepted | `ratio_text`: more digits where two would print 1.0 for a ratio above 1. Changed values: `\valRiddersWorstRatioThird` 1.0 to 1.001 (Table 10.2) and `\valForwardRiddersWorstRatioSecond` 1.0 to 1.01 |
| evidence-09 wording and comments | Accepted | 6.7 "from an initial step of about 16 K"; Table 9.10's caption names the three transcribed rules (in the margin, to keep the page); 0.074 in `code/06-ridders-onesided.cpp` and PLAN E24; `code/05-orders-float.cpp`'s comment no longer claims the single-point macros are quoted |
| codex-01 (major) | Accepted: T1 | as T1 |
| codex-02 the mirror's sign | Accepted | 6.7 "and negates the value for odd N"; Table 8.2 "value negated for odd N" |
| codex-03 smoothness | Accepted | 6.7 "Where f is smooth on both sides, use Ridders"; 10.4 and the quick reference "every method that samples across it" |
| codex-04 the actual offsets | Accepted | 6.7: the forward f' and f'' over the actual distance and offsets (Sidebars 3.1, 5.2); 8.5: "as do the f' and f'' levels of the one-sided Ridders forms" |

## Length

The fixes cost four pages (94). The editor cut them back to 90 without removing a measured fact, mostly where a
chapter repeated another or a caption:
- Chapter 9: the transcribed-rules sentence moved into Table 9.10's margin caption.
- Chapter 7: 7.8's addition was reworded to fit the paragraph's last line.
- Chapter 8: Table 8.2's mirrored row was cut to two lines.
- Chapter 6:
  - Sidebar 6.2: the exact-zero sentence goes (Figure 6.3's caption gives the counts), and so does the closing
    sentence on the tests' factor (10.2 states it);
  - 6.7: the pole-ahead figures moved to 10.2;
  - 6.4: the measured exponents moved to 10.4;
  - Sidebar 6.3: the untrusted test's definition (in 6.2) and the estimate's value go;
  - 6.6: the last sentence shortened.

## Found during the fixes, for the author

- **The repository's reorganisation (1b4dda4)** moved the deriv documents to `docs/deriv/`. The book's 16 path
  references in Chapters 3 to 6, 9 and 10 and Appendix B now read `docs/deriv/...`. `docs/deriv-book/CLAUDE.md`'s
  ground-truth list still names `docs/design.md` and `docs/PHASE2_API.md` to `docs/PHASE5_API.md`. CLAUDE.md is not
  the editor's to change.
- **U-17** (D-18) and **U-18**, worded in PLAN 12.3.

# Phase 9, second round: numerics-reviewer on the revised deriv chapter (uncommitted tree after 9a04d32)

`numerics-reviewer` (Opus, high), 2026-11-10. Saved by the main session as the agent returned it. No shell: "verified"
means derived or read from the code and the quoted output; "suspected" needs a run, named in each case.

**Short answer to the main session's question:** "within a few times its estimate" holds on the evidence that exists,
but the evidence is narrower than the wording. The sentence should give the measured factor, about twice, and say what
it was measured on. "Usually flagged" should go (finding 1).

## Findings, most severe first

**1. Moderate, verified from the recorded data: `deriv.hpp`:307-308 and `docs/deriv/design.md`:191-193.** "Usually
flagged ..., otherwise within a few times its estimate" is supported only as a measured statement on f'.
- The evidence: the book's R3 study (`docs/deriv-book/reviews/R3-numerics.md`:105-121, `reach.cpp`, 4192 calls) ran
  `ForwardRidders` for N = 1 to 4 on log, sqrt, cbrt and 1/t at x = c 2^-k. Every result beyond the reach that was
  not flagged was within 2.2 times its estimate; the worst sqrt at 2^-39: 0.33 off, estimate 0.15. Others: log at
  2^-37 (1.16 times), cbrt (up to about 1.6 times), sqrt at 1e-11 (0.62 times, so covered). For N = 2 to 4 no
  unflagged results at all.
- "Usually" is not supported: for sqrt, every distance from 2^-30 down to about 2^-37.4 gives an unflagged result,
  most covered (1e-11 is one); only below about 2^-39 is the flag raised. The data show a dichotomy: each result is
  either flagged, or within about twice its estimate.
- The example is weak: sqrt at 1e-11 (7% off, estimate 11%) is covered, so it does not show the "few times".
- Fix (header and design.md): "closer, the levels run out and the result is wrong: either flagged by an estimate of at
  least |value| / 4, or, as measured for f' of sqrt, cbrt and log near 0, within about twice its estimate (sqrt at
  1e-11 from 0: 7% off, estimate 11%; at 2^-39: 33% off, estimate 15%)."
- Test: a case at x = 2^-39 (exact in double) asserting error > estimate, error <= 3 * estimate and
  4 * estimate < |value|, with the numbers from the prototype (finding 2).

**2. Minor, verified (rule): the new test at `tests/test_onesided.cpp`:301-312.** Labelled a regression baseline, but
the file's policy does not allow one here, and it pins the wrong point.
- CLAUDE.md allows a baseline only "for which no independent value exists". The file's header (lines 3-4) says
  `.scratch/onesided/ridders5.py` reproduces the C++ results for double "exactly on 28,000 cases"; line 415 uses a
  baseline only "because the prototype does not model these types"; the 1e-12 case at line 290 cites the prototype.
  For double an independent value exists: the 5% bound and the recorded 6.8% and 11% should come from the prototype's
  run of sqrt at 1e-11, cited as at line 256.
- "The bounds leave room for libm" does not apply: the only libm call is `std::sqrt`, correctly rounded. If room is
  needed, it is for FMA contraction.
- At 1e-11 the error is 0.62 times the estimate, so `error <= 2 * abs_error` never exercises "a few times"; the case
  in finding 1 would.
- The other bounds are sound: 4 * abs_error / |value| = 0.47 (margin 2.1); 6.8% against the 5% floor (margin 1.36).

**3. Minor, verified: `docs/deriv/design.md`:195-197 contradicts the header and the manual.** "where both sides exist,
`Ridders` is better, also near a NaN edge, whose NaN guides its restarts", against `deriv.hpp`:301-302 ("either method
can be the more accurate: try both") and the manual's table (`ForwardRidders` 1e-13 against Ridders' 2e-10 at 1e-10).
- Fix: "where both sides exist, `Ridders` is usually better; close to a NaN edge either can be (the manual's table at
  1e-10)".

**4. Minor, from the recorded data (R3-numerics-01, not rerun), for the user: inside the reach, f'''' with a pole
ahead of x.** `ForwardRidders`, N = 4, on 1/(t - a) at 0, a pole ahead with a from 1e-9 to 1e-1 (inside the 2^-30
reach): 9 of 4000 results about 60 times the derivative, neither flagged nor covered, error up to 4.86 times the
estimate. `docs/deriv-book/chapters/10-guarantees.tex`:118-127 carries it; `deriv.hpp` and `PHASE5_API.md`:280 do not.
- Fix: one clause in the header, e.g. "for f'''' a pole just ahead of x, inside the span the rule samples, can give a
  wrong value that is not flagged (measured); repeat at a second step".

**5. Minor, suspected: `one-sided.md`:43-44, "on a function that changes fast near its edge, such as sqrt(1 - x) near
1, Ridders' method is the more accurate one".** Asserted at 0.999 only. Closer, Ridders' restarts exhaust its 64
evaluations (as on 1/(1 + x) at 1e-10), and `BackwardRidders` should win. Run Ridders against `BackwardRidders` on
sqrt(1 - x) at 1 - 1e-6, 1e-8, 1e-9. Fix: "Ridders' method can be the more accurate one".

**6. Minor, suspected: the sum-of-estimates rule on the page's cancelling example.** On sqrt(1 - x^2) at 0.999 both
estimates are 1.6e-10 to two digits; the values differ by 2.7e-10, less than the sum, so the rule keeps the smaller
estimate. If `BackwardRidders`' is the smaller, the rule keeps the result whose estimate misses (by 1.7 times, about
1.2e-11 relative). Not a wrong result in practice, and the next sentence covers it. Run: print both estimates to 4
digits; if so, assert the rule's pick on the factored form only, or add "the rule cannot see this".

**7. Minor, verified: `step.md`:66, "called as `default_step<Method, N>(1.0)`".** `default_step` deduces T from its
argument, so `1.0` gives the double step: in `cpp_bin_float_50`, `default_step<M, N>(1.0) * a` gives about 1e-11
accuracy instead of about 1e-34, silently. Fix: "`default_step<Method, N>(1.0)` (`T(1)` in another type `T`)".
Assertion 11 in `deriv_step.cpp` is exact and correct as written.

**8. Editorial, verified: the note at `one-sided.md`:56** leaves out the header's "from the default step"; with a step
on fn's scale the methods keep about 12 digits at 1e-11 (`handle-failures.md`:59-60). Fix: "...lose accuracy from
their default step where fn changes...".

**9. Editorial, verified: `deriv_one_sided.cpp`:214, "loses up to about 70 eps relative in 1 - x^2".** fl(x*x) in
[0.5, 1) rounds by at most 2^-54; 1 - fl(x^2) is exact (Sterbenz); so 1 - x^2 = 1.999e-3 is off by up to 2.78e-14
relative = 125 eps, about 63 eps after the square root. Fix: "up to about 125 eps in 1 - x^2 (about 60 in its square
root)".

**10. Editorial, verified: `deriv_kepler.cpp`:72 and 74 still say "far below 1e-9"**; the page now says "below".

**11. Editorial, verified: `deriv.hpp`:309 is 134 columns** (`.clang-format` ColumnLimit 120); reflow.

## Checked and correct

- N1/G1: `step.md`:47 "to the right, at smaller steps"; the alt text agrees.
- N2/U1: the old "flagged" promise is gone from the header and design.md; `handle-failures.md`:62-64 agrees, and
  `deriv_handle_failures.cpp`:68-69 asserts 6 to 8% off, covered, with 4 * estimate < |value|. The log case in
  design.md (0.224 against 0.193; 4 * 0.193 = 0.772 < 0.776) matches the book's record.
- N3/G3, `choose-a-method.md`:23-25: `Central<10>` on 1/x at 2 gives 8.4e-11 relative by its truncation term (table
  9e-11, assertion 4). `Central<8>` gives about 7e-12 on 1/x and 1.5e-13 on sqrt x: worse than `<4>` and `<6>`,
  better than `<2>`, so "the highest orders can be worse" on `index.md`:57-58 holds.
- N4/G3, `index.md`:56-58: "smooth on the scale max(|x|, 1)" is right; assertions 3, 5 and 6 match the table.
- N5/G6: the program uses sqrt(1 - x) itself; 1 - x is exact, g_exact good to about 2 ulps; assertion 13 compares
  Ridders against `BackwardRidders` at 0.999.
- Assertion 14: `Backward<2>` at 1 - 1e-7 with h = 6.06e-6 gives -452 against -1581: 71% off, "more than 10%" holds.
- N6/G5: the one-sided stencils advice is limited to fn that changes fast near its edge; the mechanism is right; the
  `Forward<2>` and Newton exception is right.
- N7/G4: correct for every method (Ridders gives 1/32 times the scale), apart from finding 7.
- N8/G8: "below 1e-9", the 3e-10 attributed to the solve; the bound 2.64e-10 + 1.3e-11 holds.
- G9: `default_step<Central<2>>(M)` at M = 1 is 6.06e-6, inside (6e-6, 6.1e-6).
- N9/G10: the 1e-9 is the header's 2^-30 = 9.3e-10, now linked to the `ForwardRidders` reference.
- N10/G7: the alt text matches the quoted table (0.0353, 5.1e-12, 4.09e-12, 2.52e-3); assertion 11 pins it.
- N11: the `cpp_bin_float_50` model 3.0e-34 and the measured 5.7e-35 agree; the bound 1e-32 is about 64% of the
  digits, still "about two thirds".
- U3: "can make an estimate optimistic" is right ("both" would be false: Ridders' 4.9e-12 is covered by 1.6e-10); the
  cancellation (125 eps in 1 - x^2 against about 1.5 eps for (1 - x)(1 + x)) explains `BackwardRidders`' miss by
  1.66; assertions 26 and 27 are deterministic on the five toolchains (only + - * / and correctly rounded sqrt, no FMA
  contraction there; a contracting target such as aarch64 could fail 26, as the file's header states); the reference
  -at / sqrt((1 - at)(1 + at)) is good to about 4e-16.
- G12, the sum rule: if both estimates cover, |v1 - v2| <= e1 + e2, so a larger difference proves one misses. On the
  table it keeps Ridders at 1e-2, `ForwardRidders` at 1e-6 and 1e-10: the more accurate each time. Assertion 25 is
  legitimate since |f'| = 1 within 2e-10.
- D6 and C9: assertions 21, 23 and 24 match the table (64 evaluations, error 2e-10, estimate 8e-6, NaN counts growing
  toward the edge).
- G11: the kink estimate 1.2e-8 is in the quoted output; assertion 10 is consistent.
- The 1e-12 flagged cases (lines 289-299) and the cost 1 + 32N are unchanged and consistent with the header's caps.

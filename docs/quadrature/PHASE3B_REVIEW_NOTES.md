# Phase 3b: review notes on draft 1 (not yet folded in)

Working notes, 2026-11-02, so that the review round can continue in another session. Fold them into
`PHASE3B_API.md` (draft 2, "Decisions from the reviews") and delete this file.

Status of the round:
- **Simplicity review:** done (below).
- **Usage sketch:** done (below). Its finding 1 has been reproduced (below).
- **Numerics review:** still running when these notes were written. Rerun it on draft 1 or on draft 2. The numerics
  findings must also weigh simplicity findings 6 and 7.

## Simplicity review

1. **API "Comments":** say "`PHASE3_API.md` Comments stand, plus:", and list only what is new.
   - GaussKronrod's "until TanhSinh" becomes "or use TanhSinh".
   - TanhSinh's comment gains the cost of constant evaluation for a T without a table.
   - Shorten the member comments to Romberg's form.
   - Drop the `state_t` sentence.
   - The lead-in "three comments" contradicts the five bullets that follow it.
2. **Drop Q4 and Q5:** they are not decisions.
   - D1 is a clarification of the approved spec, which the prototype confirms. Its alternative is a known false
     success.
   - Q4 is narrow, and its alternative is wrong: W needs the weights of every one of level k's nodes beyond the
     cut, not only the nodes evaluated.
3. **Drop the step that adds the printing example's output** to the API document.
4. **Tests:** merge and drop the duplicates.
   - Merge "evaluations per level match B" with the prototype's counts, recorded per level through the observer.
   - The protocol properties go into `test_quadrature_protocol.cpp`'s existing lines, as Romberg's do: observed
     equals unobserved, last report equals the result, named and as a value, reversed ends.
   - One `iterative_method` `static_assert`, for double.
   - Drop the duplicate D1 test bullet.
   - Keep TanhSinh's own tests:
     - level 0 has no current and is not ready;
     - `ready` from level 4;
     - 10 steps, only the last exhausted;
     - no point inside: one report, exhausted, fn never called.
5. **Q3, the prototype:** keep it, but limit it to a double-only loop of A to E that reads the C++ tables.
   - It covers the count cases and D1's case.
   - It needs only Python floats, no mpmath.
   - Reject the alternative of no count tests: the counts are the only check of the cut and pruning logic.
6. **The state:** drop "the slivers' share of the floor" (it can be derived from W, y_e and h).
   - **This affects results:** V_k in the recurrence V_k = V_(k-1)/2 + new terms must be the pure trapezoid
     sum.
   - The sliver term(h W, y_e) is added only to the reported value.
   - V_(k-1) in the estimate is the previous reported value.
   - Optional: keep V and L1 as plain T in the state, with each level's `Sum` local, as `RombergState` does.
7. **Unspecified, and it affects the counts: the pruning threshold at level k.** Which L1: L1_(k-1), the running
   sum, or the sum after the side's inner nodes? Recommended: eps L1_(k-1), which is known before the level starts,
   the same for both sides, and independent of the order of evaluation. D should name the choice.
8. **Trims:**
   - A and B's loop bullets restate G.
   - Steps 1 and 2's test bullets repeat `PHASE3_API.md`'s.
   - Drop the "No Raw Loops" paragraph.
   - Shorten Deferred.
9. **Q2, three steps:** defensible; keep.

**Already as simple as it should be:**
- Q1 and C. "At most one new node beyond the last kept node" was checked for every kind of cut.
- Every per-side field is needed: the geometric pair, the distinct-position pair, which can differ near an end of
  magnitude 1, and W.
- B's table; D2.

## Usage sketch

Its sketch programs were in a session scratchpad and are not kept.

1. **High: the `or_else` fallback turns TanhSinh's honest failures into false successes.** Reproduced with today's
   `GaussKronrod<15>` in double:

   | case | the default's result | error | against its estimate |
   |---|---|---|---|
   | x^-0.99 at rel_tol 1e-4 | ok, 99.4681531307 ± 0.00991 | 0.532 | 54 times |
   | x^-0.9 at the default | ok ± 1.43e-7 | 7.05e-7 | 5 times |
   | (1 - x)^-0.9 at the default | `not_converged`, best ± 0.0362 | 0.234 | 6.5 times |
   | x^-0.999 at the default | `not_converged`, best 509.09 ± 0.907 | 491 | 541 times |

   Proposed:
   - The comment says the fallback is for interior features only.
   - The tested example falls back only on `not_converged`. When both methods fail, it keeps the best with the
     smaller estimate.
   - Add a test that the documented pattern gives no success on x^-0.99 at 1e-4.
   - This revisits `PHASE3_API.md` Q9's fallback: a decision for the user.
2. **"Fails honestly on interior features" overstates.** `PHASE3_API.md` F records 5 and 3 false successes among the
   250 interior cases. Say "usually".
3. **An absent best is not always divergence:** x^-0.999 converges too slowly for T's range. Mention the
   substitution x = u^k, or a wider type, and that the default is no better.
4. **Substitution at an end that is not 0:**
   - name both forms, u = b - x and u = x - a;
   - warn that composing fn(b - u) does not help, because of rounding: fn must be rewritten in u;
   - give a worked example: 1/sqrt(1 - x^2) on [-1, 1], split at 0, with each half as
     ∫_0^1 du / sqrt(u (2 - u)).
   - This is the strongest argument for the deferred (x, xc) form.
5. **A success's `abs_error` can exceed the request** (derived, to confirm with the prototype). The slivers join
   F, and so 2F.
   - fn = 1 on [1, 1 + 1e-9] passes with `abs_error` ≈ 2.2e-7 relative, at the default too.
   - Update the comments of `Criteria` and `Step::attainable`, and add a line to TanhSinh's.
6. **Observing is underspecified.** TanhSinh's comment should say:
   - report i is level i (Romberg's report i is level i + 1);
   - a later report with no current means the estimate is unbounded;
   - `ready` from report 4, `exhausted` at report 9;
   - the evaluations per level vary, roughly doubling.

   Other points:
   - The summary's `print` should be the `Printer` that prints a dash.
   - A failure in `start` produces no report, which is the common case for `non_finite` (sqrt(x - 0.5)): advise
     recording fn's arguments.
   - 1/x runs to the cap: show an observer that throws to stop it.
7. **Cost guidance:**
   - TanhSinh's minimum, 70 to 140 evaluations, is below the default's 150. On the smooth and endpoint battery it
     needs 2.7 times fewer evaluations.
   - Say: for the fewest evaluations on a smooth fn, TanhSinh; below about 100, `GaussLegendre<N>`.
   - Give one line on why the default is not TanhSinh: interior features.
8. **Compile time and constant evaluation:**
   - The compile time, 0.5 to 1.5 s per translation unit, applies to every use. Give it its own bullet, with the
     advice to instantiate in one translation unit.
   - "+ - * /" excludes sqrt, log and pow on Clang and MSVC: say "rational fn, or your own constexpr functions".
   - Give the constexpr steps per evaluation.
   - Run the constexpr 1/x test on all five toolchains.
9. **Multiprecision:**
   - "Full precision" means `rel_tol = 0` (stopping at 2F); otherwise drop the word.
   - Keep `-> mp` in the example, since expression templates are rejected without it.
   - Add the first call's cost: 14 µs per node, 154 ms for all levels in `cpp_bin_float_50`, 343 ms in
     `cpp_dec_float_50`. A type that changes its precision sees its first tables.
10. **Comments that will go stale:**
    - the file's @brief and its paragraph on the iterative methods;
    - `iterative_method`'s first sentence;
    - `integrate`'s `@tparam M` lists;
    - `Step::attainable`, `Step::ready` and `Step::exhausted`;
    - `Criteria`.
11. **Combinators:**
    - Document that the state is independent of fn's type, of fixed size and allocates nothing.
    - Show a `SplitAt<Inner>` estimating combinator in the guide.

**Works well:**
- no new surface;
- the state;
- the failures map cleanly onto the protocol;
- level 0 is reported;
- the tests;
- D's two points.

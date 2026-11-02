# Phase 3b: review notes on draft 1 (not yet folded in)

Working notes, 2026-11-02, so that the review round can continue in another session. Fold them into
`PHASE3B_API.md` (draft 2, "Decisions from the reviews") and delete this file.

Status of the round:
- **Simplicity review:** done (below).
- **Usage sketch:** done (below). Its finding 1 has been reproduced (below).
- **Numerics review:** done (below; the reviewer derived its findings by hand, without running scripts). Its
  findings 1 and 4 cover simplicity findings 6 and 7. Rerun the numerics review on draft 2.

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

## Numerics review

1. **The running sum and the slivers are not separated (C).** As written, the slivers are counted twice. If V_k
   includes S_k = term(h_k W_k, y_e), then V_(k-1)/2 carries S_(k-1)/2, whose weights are all in W_k again.
   - The error is d_k = (d_(k-1) + S_(k-1))/2, which tends to S.
   - fn = 1 on [1, 1 + 1e-9], which `PHASE3_API.md` says is exact, would be off by about 2.2e-7 relative. The
     result is covered by the floor, so it is not a false success, but the value is wrong.
   - **Fix:** the state keeps the raw sum R_k = R_(k-1)/2 + the new terms. The reported value is V_k = R_k + the
     slivers, and |V_k - V_(k-1)| uses the previous reported value. The prototype does the same. This matches
     simplicity finding 6.
2. **D1's evidence is wrong; its conclusion is right.**
   - x^-0.99 in double overflows below max^(-1/0.99) ≈ 4.3e-312, which is subnormal. The mass below that point is
     0.0770, which is `PHASE3_API.md` E's "0.077". So Phase 3's prototype table went into the subnormals (A records
     5.5e-318), and the cut there was a tail cut.
   - Moving the table's end to min made x^-0.99 a table-end case. D1 is what keeps B's "no change in any result"
     true.
   - **Soundness:** for x^-alpha the terms are log-concave on [t_e, ∞) wherever they decrease at t_e >= 0.88. So
     h Te / (1 - r) bounds the mass beyond the table; where they do not decrease, the term is unbounded.
   - **Correct:** fn is finite down to 4.3e-312, not "down to min". 100 min^0.01 = 0.0838.
   - **Fix:** replace the evidence with this argument and these figures. Q5's prototype runs the table to min, and
     also the old subnormal table, to show the case moved from a tail cut to a table-end cut. Q5 then becomes a
     clarification (simplicity finding 2).
3. **The changes of a side's cut omit the table-end kind, which can become a rounding cut.** The counterexample is
   x^-0.99 on [0, 1e-40]: a table end at level 0 (t = 6), then rounding at level 1 (t = 6.5, since hw comp
   underflows).
   - The full set of changes:
     - rounding to rounding, pruning or tail;
     - pruning to tail, and tail to pruning;
     - table end to table end, pruning, tail or rounding.
   - "Cut at t_last + 2^-(k-1)" is false for a table end, but "at most one new node beyond" still holds.
   - **Fix:** list the four kinds and their changes. W starts at 0 and is accumulated at every level whatever the
     kind; it stays 0 for a table-end side. Keep the data on the last kept point and the nearest point at a
     different position for every kind.
4. **The pruning details are open, and one reading is unsound.** The exact counts depend on them.
   - (a) Level 0 drops the run of outermost small terms, taken from the outside in. Cutting at the first small term
     going outward is a false success. Counterexample: (x - x2)/sqrt(x) on [0, 1] with x2 = 7.065e-4, the point at
     t = 2. The term there is exactly 0, the terms beyond are not small, and 2.5e-5 is missed.
   - (b) Level 0's threshold uses the L1 of the whole sweep; L1_0 is then formed from the kept terms.
   - (c) The test at level k needs a stated L1: L1_(k-1), recommended by the simplicity review; the running sum; or
     L1_k.
   - (d) The order of summation: the sides, and the inner nodes before or after the node beyond.
   - **Fix:** pin all four in D2.
5. **The constant-evaluation cost is understated.**
   - Nodes per level: 7, 7, 13, 26, 53, 105, 210, 420, 840, 1680.
   - Totals at 275 steps a node:

     | levels | nodes | steps |
     |---|---|---|
     | to 4 | 106 | 29k |
     | to 5 | 211 | 58k |
     | to 7 | 841 | 231k |
     | to 9 | 3,361 | 924k |

   - Clang's budget is per top-level evaluation, so the table's generation shares it with the run. A user's literal
     class type costs several times 275 steps a node. A run past level 7 likely exceeds the default.
   - The test "1/x in double in constant evaluation" runs 5,164 evaluations in one evaluation: about 203 steps each
     under Clang's default.
   - **To measure:** a constexpr wrapper of double, at level 4 and to the cap; the 1/x test on Clang, clang-cl,
     Emscripten and cl.
6. **A: "replaces A's first bullet" drops too much.** Keep Phase 1's `Interval` with its scale order, and the
   constant-evaluation refusals (`checked_term`, `fits`, `compensated`) that give `non_finite`. Add them to the
   failures list. Replace only the validation and the orientation.
7. **The state must hold the centre's term.** It counts as each side's innermost kept node: for "a side with no
   kept node: h |w(0) y(0)|", and for Tp when the outermost kept node is t = h_k.
8. **The tests:**
   - "Evaluations per level match B" needs the rule:
     - level 0: 1, plus per side each integer t up to and including an evaluated cut node;
     - level k: per side the odd multiples of 2^-k up to t_last, plus 1 if the node beyond is in the table and
       strictly inside.
   - "The last report equals the result" holds only for a < b.

**Noted, not raised against 3b:** a single node with fn exactly 0 can turn a tail, table-end or rounding side into
a pruning cut, and so drop its end term. This is `PHASE3_API.md` C's approved rule.

**Checked and correct:**
- B's fields, and that the loop reproduces E bit for bit; the edge cases;
- level 0's 13 evaluations at most;
- 3,361 nodes;
- at most one new node beyond, by induction;
- the incremental W equals the suffix sum up to rounding;
- the end terms need only their own level;
- D2's halving;
- the cost arithmetic.

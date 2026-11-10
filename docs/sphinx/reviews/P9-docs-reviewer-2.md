# Phase 9, second round: docs-reviewer on the revised deriv chapter (working tree against 9a04d32)

`docs-reviewer` (Opus, high), 2026-11-10. Saved by the main session as the agent returned it. Probes in the session
scratchpad (`p9-review/docs-reviewer-2/`: `probe1.cpp` to `probe4.cpp`, built with `-lstdc++exp`).

**No major findings.** Every first-round finding the letter adopted or adapted is done, and done correctly; U1 to U4
are done correctly. The revision introduced four minor problems and some editorial ones.

**Builds:** the GCC 16 build is up to date and warning-free (the examples with `-Wall -Wextra -Wpedantic -Wshadow
-Wconversion -Wsign-conversion -Werror`). All 12 deriv ctest entries pass. `numerixx_tests` passes
`deriv:*,orders:*,ridders:*,onesided:*` (34 cases, 575 assertions). `numerixx_docs` (site and LaTeX) builds, and the
`_ext` checks pass.

## Minor

**1. one-sided.md:9, narrower than the header and handle-failures.md.** "The methods can only tell that a point lies
outside the domain if `fn` returns NaN there": the revision dropped "(or an infinity)", which no letter item asked
for. The header (deriv.hpp:272): "It restarts below a step at which fn is not finite"; handle-failures.md:13: "`fn`
returned NaN or an infinity". probe1, fn +inf below 0: `Ridders` at 1e-6 succeeds (-0.99999799998).
- Fix: restore "returns NaN (or an infinity) there", or "a value that is not finite".

**2. one-sided.md:78-79, the kink section no longer says what to do.** The revision dropped its last sentence, "So
return NaN outside the domain, or use a one-sided method"; G11 changed only the first sentence. Fix: put it back after
"...and gets -1."

**3. derivative-of.md:27-33 and the `<own-step>` block: nothing after the code, and a lifetime trap.** The block is
good (qualified, `<dv::Ridders>`, asserted at 2 within its estimate; rerun: error 1.1e-14, estimate 4.0e-14).
- No text follows the block; the rules ask for an explanation after a block (for Ridders' method `h` is the initial
  step; the code comment says so, the prose does not).
- The lambda is `[&]`; the next sentence says `derivative_of` "can outlive `fn`". The reader's own lambda cannot: one
  returned from a function dangles. The page's capture warning is about `derivative_of`'s copy, three paragraphs on.
  The header comment at deriv.hpp:1894 has the same `[&]`.
- Fix: "For Ridders' method `h` is the initial step. The lambda refers to `fn` and `h`, so unlike `derivative_of` it
  must not outlive them; capture by value (`[fn, h]`) to keep copies." `[=]` compiles and works (probe3, `make()`).

**4. one-sided.md:47, "trust neither until you have checked at another step": no means given.** Fix: "...checked at
another initial step, the third argument of `derivative` ([Choose a step](step.md))".

## Editorial

**5. tutorial.md:63, "a few ulps".** "ulp" appears nowhere else and is not defined (the wording came from this
reviewer's first round). Fix: "Each E(M) is off from the root solve by a few units in its last digit, which adds at
most about 3e-10."

**6. one-sided.md:47-48, "cancels" and "optimistic".** "A formula often cancels" uses an undefined term: "...cancels
near an edge (subtracts nearly equal numbers and loses digits), as 1 - x^2 does near 1". "Optimistic" is what
handle-failures.md:69 calls "too small": one name for one thing. The sentence itself is true and asserted (probe1 at
0.999: `BackwardRidders` error 2.66e-10, estimate 1.63e-10; Ridders error 4.9e-12, estimate 1.58e-10; on the factored
form both cover). Narrowing "both estimates" was right: across x = 0.9 to 0.99999 the cancelling form makes Ridders'
estimate miss at 0.998 (5.4 times) and `BackwardRidders`' at 0.999. The sum rule picks Ridders at 0.999, the right one.

**7. one-sided.md:3-5 and 27.** The first paragraph no longer states the page's job ("One-sided methods avoid that."
was cut, unasked): restore it, or end with "One-sided methods evaluate `fn` on one side only." Line 27, "with an
estimate", is the term's first use on the page: "with an error estimate, `abs_error`". Lines 32-33 dropped the clause
saying the table's error and estimate are relative (harmless here, |f'| about 1).

**8. one-sided.md:58-59, the note's polynomial.** "on a polynomial that vanishes at the points Ridders' method tries"
holds only because fn is NaN below 0: the restarts land the steps on the zeros (probe1: without the NaN, Ridders
returns 1 with an estimate of 2.2e-16). Fix: "on a polynomial, NaN below 0, whose zeros lie at the points Ridders'
method tries after its restarts".

**9. one-sided.md:52, `Forward<2>` only.** For an upper edge the stencil is `Backward<2>`: "`Forward<2>` or
`Backward<2>` is enough".

**10. deriv.hpp:309, a new line of 133 columns** (the U1 edit; ColumnLimit 120); reflow. Out of scope: the `@code`
block below the U2 fix (deriv.hpp:1900-1902) still writes unqualified `derivative_of`, `Ridders` and `order::second`;
some header `@code` blocks qualify fully (13, 69, 285, 315), others do not (1657, 1725).

**11. tests/test_onesided.cpp:301-304, a tolerance justified wrongly.** "Room for libm" is not the reason: `sqrt` is
correctly rounded and `ForwardRidders` uses only `+ - * /`. The margin is against FMA contraction or a later change
of the method: "The bounds leave room for a target that fuses a * b + c, and for changes of the tableau".

**12. basics.md:163 and 233.** "the Ridders methods" without a gloss, while the stencils have one: "the Ridders
methods (`Ridders`, `ForwardRidders`, `BackwardRidders`)". At 233, "One more exception" now sits two sentences from
the claim it qualifies; "A constant can also differ in one more way:" would be clearer.

## Checked and found correct

- D1 to D18 as the letter decided them: step.md "to the right, at smaller steps"; the lambda in the program and
  asserted; `default_step<Method, N>(1.0)` times the scale right for the stencils and Ridders, asserted, as is
  `default_step<Ridders, 2>(2.0) == 0.5`; index.md lists `<4>` to `<10>`, links Costs, claims asserted in
  `deriv_costs.cpp`; the kink sentence matches the header; "restarts below every step that reaches past the edge"
  matches deriv.hpp:272, NaN counts asserted; the Newton exception agrees with roots/restricted-domain.md; `copied` is
  non-const; the pragma names match deriv.hpp:1788; sqrt(1 - x) itself asserted; the Kepler default step range
  asserted; no prose line over 120 columns in the deriv pages or basics.md.
- one-sided.md in depth: the 1e-10 numbers quoted and asserted; "Ridders' method is the more accurate one" on
  sqrt(1 - x) holds at every distance from 1e-1 to 1e-10 (probe1); the sum rule holds on all three rows; the
  cancelling miss and cure asserted and reproduced; the note's 1e-9 linked; the kink's 1.2e-8 quoted.
- The header comments: ForwardRidders' new wording holds. probe2 swept sqrt, log, cbrt, x^0.1 and 1/x at 2^-30 to
  2^-52, N = 1 to 4: 828 results, 215 not flagged; 11 of the 215 have an error above the estimate, all N = 1, at worst
  2.22 times. sqrt at 1e-11: 6.84% off, estimate 10.9% of the value; log at 2^-37: error 0.2242, estimate 0.1931,
  unflagged; all with the full 33 evaluations (probe4). The qualified `derivative_of` lambda compiles in user code,
  also for N = 2 and captured by value (probe3).
- The new test passes, is labelled as a regression baseline with its commit, and measures against the closed form.
- basics.md: `if (*d)` is false for a zero derivative and true for 3 (probe3); the overflow sentence no longer splits
  the FMA sentences.
- No other regression in the diff against 9a04d32; the `<own-step>` block renders correctly.

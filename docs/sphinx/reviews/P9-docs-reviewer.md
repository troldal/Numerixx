# Phase 9, stage 5: docs-reviewer on the deriv chapter (commit d00357e)

`docs-reviewer` (Opus, high), 2026-11-10. Saved by the main session as the agent returned it. Probes in the session
scratchpad (`p9-review/docs-reviewer/`: `probe.cpp` to `probe4.cpp`, and the rendered figures).

One major finding: step.md says the peak's curve lies to the left, but the figure draws it to the right. The others
are minor or editorial.

**Builds and tests.** All 11 `numerixx_docs_example_deriv_*` programs rebuilt from clean with GCC 16: no warning, all
pass under ctest. `numerixx_docs` builds with no warning, and the `_ext` checks pass. Both figures rendered and looked
at.

## Major

**1. step.md:47: false statement about figure 2.** The page says "the peak's lies about four decades to the left".
`plots/deriv_scale.py:53` puts smaller steps on the right (`set_xlim(max*2, min/2)`), so the peak's lowest error, at
h = 1e-9, is to the right of the square root's, at 5e-6. The alt text says "about four decades further right", and
the rendered figure agrees. The plan's figure 2 text predates the decision to put smaller steps on the right.
- Fix: "...but the peak's lies about four decades to the right, at smaller steps".

## Minor

**2. derivative-of.md:25: a code line that does not compile as written.** `[&](double x) { return derivative(fn, x,
h); }`: in user code `derivative` is unqualified and ADL does not find it (GCC: "'derivative' was not declared in this
scope; did you mean 'nxx::deriv::derivative'?", probe2.cpp, `-DUNQUALIFIED`). It also drops `<Method, N>`, so it
silently becomes `Central<2>` for f'. A complete line a reader will copy, not from an example program.
- Fix: `[&](double x) { return dv::derivative<Method, N>(fn, x, h); }`, ideally held in `deriv_derivative_of.cpp`.

**3. step.md:65-66: "default_step gives the default step, which you can scale the same way" can give a wrong step.**
`default_step<M, N>(x)` is already eps^(1/(P+N)) max(|x|, 1) (deriv.hpp:1648, 1678). Multiplied by fn's scale a it
is |x| times too large when |x| > 1 (x = 50, a = 0.1: 50 times the intended step). For the Ridders methods it is
max(|x|, 1)/32 or /4, not a root of epsilon. The page's example has x = 2e-4, so it does not show the problem.
- Fix: "...gives the default step; divide it by max(|x|, 1) and multiply by fn's own scale". For a stencil this is
  `default_step<M, N>(0.0) * a`.

**4. index.md:54-55: the Central<4> and Central<6> claim.**
- True: on a smooth function `Central<4>` and `Central<6>` beat `Central<2>`; in probe.cpp on all 10 smooth functions
  tried (1/x, sqrt, exp, sin, log, 1/(1+x^2), atan, x^3, exp at 10, cos). The costs and choose-a-method tables agree.
- Problems: the sentence is neither asserted by `deriv_overview.cpp` nor linked to Costs. It implies `Central<8>` and
  `Central<10>` are not more accurate than `Central<2>`; in fact `Central<8>` beats it on 9 of the 10, `Central<10>`
  on 6 of 10 (worse on 1/x and exp at 10, about even on log at 3). What holds is narrower: orders 8 and 10 lose to 4
  and 6 when a singularity is near x. "Central<4> to Central<10>" invites `Central<5>`; only even P is allowed
  (deriv.hpp:123).
- Fix: "`Central<4>`, `<6>`, `<8>` and `<10>` use more points; on a smooth function `Central<4>` and `Central<6>` are
  the more accurate, and higher orders can be worse ([Costs](costs.md))".

**5. one-sided.md:71: "Ridders' method sees a kink at 0" contradicts the header.** deriv.hpp:1719 says a kink at x
itself is invisible even to Ridders' estimate, and the page's own result (-2 with an estimate of 1.2e-8) shows it.
- Fix: "At the kink at 0, Ridders' method returns -2, the mean of the two slopes, with a small estimate: nothing
  tells it the kink is there."

**6. one-sided.md:41: "Far from the edge... Closer, it restarts below each step that met NaN".** Implies no restarts
at 1e-2, but Ridders restarts there too (its first step of 1/32 reaches past 0). Evaluations at NaN points
(probe3.cpp):

| Distance | Evaluations | At NaN points |
|---|---|---|
| 1e-2 | 16 | 2 |
| 1e-6 | 36 | 15 |
| 1e-10 | 64 | 29 |

The table itself is right: 64 is the cap, 32 levels of 2 (`ridders_levels`, deriv.hpp:1246); ForwardRidders' 10
evaluations and every estimate covering its error are asserted.
- Fix: "It restarts below every step that reaches past the edge, and the closer x is, the more often: its estimate
  grows, and at 1e-10 it makes all 64 evaluations it may."

**7. one-sided.md:50-52: contradicts the linked roots section.** The page says "Avoid the one-sided stencils there...
Newton's method with a numerical derivative meets the same problem", but
`roots/restricted-domain.md#a-numerical-derivative-near-the-edge` recommends `Forward<2>` near the edge: sound for
Newton's method, because the derivative's accuracy does not affect the root's bound.
- Fix: say that for Newton's method a one-sided stencil is enough near the edge, and link there.

**8. derivative-of.md, "Its copy of fn": a compile trap not mentioned.** A callable copied from a non-const fn (the
`Counter`, or a `mutable` lambda) cannot be called through a const result: `auto const copied =
dv::derivative_of(counter); copied(2.0);` fails with "no match for call to '(const
nxx::deriv::detail::DerivativeOf<...>) (double)'" (probe2.cpp, `-DCONSTCOPY`); deriv.hpp:1888 documents it. The
section writes `auto copied` without saying why, right after `auto const slope`.
- Fix: one clause: "call it through a non-`const` variable".

**9. higher-derivatives.md:71: the pragma is not named.** "wrap it in your compiler's pragma for deprecation
warnings" does not give the warning's name. The header does (deriv.hpp:1787): `-Wdeprecated-declarations` for GCC
and Clang, C4996 for MSVC. Fix: name them.

**10. higher-derivatives.md:31: term misuse.** "the central difference keeps about 10 digits for f', 8 for f''...":
in the plan's terms the central difference is `Central<2>` for f' only; line 26 calls the f'' rule "the second
difference". Fix: "`Central<2>` keeps about...".

**11. choose-a-method.md:19: an unlinked number.** "Typically 13 to 15 digits" for Ridders is neither asserted nor
linked. It is true (13.2 to 15.0 digits across 10 functions, 13.2 to 14.7 in the costs table). Fix: link to Costs.

## Editorial

- **12. one-sided.md:43.** The page names sqrt(1 - x), but `deriv_one_sided.cpp:141-151` asserts on sqrt(1 - x) + x.
  The claim holds for sqrt(1 - x) itself: Ridders 2.3e-14 against BackwardRidders 1.9e-13 at 0.999 (probe.cpp).
  Name the function the program uses or assert the one the page names.
- **13. tutorial.md:67.** "Each value of E(M) comes from a root solve, but its rounding adds at most about 3e-10":
  the 3e-10 is the root solve's error, up to 1.6e-15 per value, divided by 2h (kepler.cpp:72-75); not "rounding",
  and "but" does not fit. Suggested: "Each E(M) is off by a few ulps from the root solve, which adds at most about
  3e-10."
- **14. tutorial.md:40.** "About 6e-6 here" is not asserted by kepler.cpp. True: eps^(1/3) at M = 1, the value
  `deriv_step.cpp:84` asserts at another x.
- **15. index.md:59.** `fn` is first used without being introduced; make line 5 "the function, `fn`".
- **16. derivative-of.md:26.** "passes such a function to Newton's method" follows the own-step lambda, but
  roots/derivatives.md passes `derivative_of(f)`. Suggested: "passes `derivative_of(fn)` to Newton's method".
- **17. derivative-of.md:30.** "So the result outlives fn" should be "can outlive fn".
- **18. Leftover rewraps.** Prose source lines over 120 columns at choose-a-method.md:5 (150) and index.md:48 (123).

## Unrelated, out of this phase's scope

- **basics.md:230-236.** "One more exception: an overflow..." was inserted between "...fused multiply-add at run
  time." and "Compilers for x86-64 and WebAssembly do not by default", so "do not by default" now reads as if it were
  about `nxx::quadrature` refusing a value.
- **basics.md:162-163.** Still "Ridders' methods" and "the plain finite differences" (deferred, plan question 2).

## Checked and found correct

**API statements against deriv.hpp:** orders 1 to 4 as an integer or an `order`, N = 1 and `Central<2>` as defaults
(lines 1744-1780); `derivative(fn, x, h)` takes the initial step for the Ridders methods; `derivative_of` bit-for-bit
the same (asserted); `Estimate`'s two members, `abs_error` empty for the stencils; `error::non_finite` and its
causes (points that coincide, a step not positive and finite, overflow; Ridders and `Central<4>` return it for
h = 1e-300 at 1); Ridders' initial step max(|x|, 1)/32 for f' (/4 for f'' to f''''); ForwardRidders starts at 1/32
and needs fn to vary on a scale of at least about 2^-30 ("about 1e-9"); the untrusted-result flag, an estimate of at
least |value|/4 (line 1563), for all three Ridders methods; the Ridders methods fail only when fn(x) is not finite or
no tableau reaches two finite levels; the precision warning (N >= 2, a coarse type, every method, none from
`derivative_of`, the macro defined alike in every file; tests/precision_warning.cpp); stencil costs P + N - 1
centrally, P + N one-sided; costs.md's pointer to the reference (counts in the Ridders and ForwardRidders comments,
lines 277-278 and 308-309); `nxx::math::root_epsilon` exists and resolves.

**Numbers:** every quoted output matches its raw string, and every prose number is pinned by output or an assert:
6.0e-12 and the factor of about a hundred; the first ten decimals 1.0373620219; the 14, 36 and 14 evaluations;
3.4e-4; 6.1e-10 and 4.6e-11; 0.2500806; the 7% and 70% cases; the 1e-46 estimate and its 40 correct digits;
Forward<2> at about 10 digits and ForwardRidders at about 13; 64 evaluations at 1e-10. The tolerance justifications
check out where recomputed (the Forward<2> model 2h^2, the Kepler bound 2.64e-10, the h values in deriv_higher, rho
in deriv_trade_off).

**Figures:** figure 1's alt text and caption match the table (3.5e-2, 4.1e-12 at 5e-6, 2.5e-3 at 1e-13, the default
step 1.21e-5, Ridders 5.84e-14 below the whole curve); its asserts cover the fourfold fall, the minimum within a
factor 8, about two thirds of the digits, and the rounding bound. Figure 2's alt text matches the table (the peak's
error at or above 1 for h >= 2e-4; its lowest 1.23e-11 at 1e-9; its default step 6.06e-6, error 3.43e-4). The
decimal steps: both recorded artefacts of powers of two reproduced (sqrt at 2: the same value at six of seven steps
2^-16 to 2^-22; the peak: 3.1e-16 at 2^-31, 2^-34 and 2^-45).

**Links and structure:** every anchor resolves; no references to `accuracy.md` or `deriv_accuracy` remain outside the
plans and reviews; the edited note in roots/derivatives.md is correct; Breathe directives only on reference.md; all
code through `literalinclude` except the inline line in finding 2; the overview and the tutorial end with what to
read next; no forbidden screen-only words.

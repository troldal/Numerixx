# Phase 9, stage 5: numerics-reviewer on the deriv chapter (commit d00357e)

`numerics-reviewer` (Opus, high), 2026-11-10. Saved by the main session as the agent returned it. The reviewer had no
shell: every result comes from re-deriving it or reading the code; one is marked suspected and needs a run.

I read every page in `docs/sphinx/deriv/`, the 11 `deriv_*.cpp` programs, the two plot scripts, and checked them
against `deriv.hpp`, `math.hpp`, `tests/test_onesided.cpp` and `docs/deriv/design.md`.

The mathematics is sound. Every error formula, constant and tolerance I re-derived holds. The problems are one false
sentence about figure 2, one false statement in a header comment that the manual's own program contradicts, and some
wording that could mislead.

## Findings, most severe first

**1. Major, verified: `step.md`:47 gets figure 2's direction wrong.** The page says the peak's curve "lies about four
decades to the left". Both plot scripts reverse the x axis (smaller steps to the right). The peak's lowest error is at
1e-9 and the square root's at 5e-6, so the peak's curve lies to the right. The alt text (line 42) says "about four
decades further right". `PHASE9_PLAN.md`:129 has the same mistake.
- Fix: "lies about four decades to the right, at smaller steps", or "at steps about ten thousand times smaller".

**2. Major, verified, in the ground truth, not the page: the `ForwardRidders` header comment (`deriv.hpp`:306-307)
is contradicted by a manual program.** The comment says that closer than 2^-30 max(|x|, 1), the result "is wrong but
flagged by an estimate of at least |value| / 4".
- `deriv_handle_failures.cpp`:68-69 asserts the opposite for sqrt(x) at 1e-11: the result is 6 to 8% off and its
  estimate is below |value|/4.
- `test_onesided.cpp`:289-299 tests the flag only at 1e-12.
- `handle-failures.md`:61-63 states the gap correctly ("The result at 1e-11 did not raise that flag").
- A user who trusts the header would accept a result that is 7% off.
- Fix: propose to the user a header comment like "closer, the result can be wrong without a flag (sqrt(x) at 1e-11:
  7% off); much closer (1e-12) it is flagged by an estimate of at least |value| / 4", plus a test at 1e-11.

**3. Minor, verified: `choose-a-method.md`:24-25 blames the failure of higher orders only on functions "that change
fast".** The costs table shows a high order losing accuracy on slowly varying functions too:
- `Central<10>` on 1/x at 2: 9e-11, worse than `Central<2>` (4e-11).
- `Central<10>` on sqrt(x) at 2: 1e-12, against 2e-14 for `Central<4>`.
- Cause: the default step eps^(1/(P+1)) assumes derivatives of size f/s^k. A function with a singularity at distance
  s has f^(k) of about k!/s^(k+1). For `Central<10>`, the step 0.0377 * 2 gives a truncation error of about
  (5!)^2/11! * 11! * (0.0377)^10, about 8e-11 relative, which matches the table.
- Fix: "...so it can be less accurate even on a smooth function (`Central<10>` on 1/x at 2, [Costs]), and on a
  function that changes fast it can lose every digit."

**4. Minor, verified: `index.md`:55 says "on a smooth function `Central<4>` and `Central<6>` are more accurate".** The
peak x/(x^2 + a^2) is smooth (analytic on the real line), yet `Central<4>` and `Central<6>` lose every digit on it
(costs table, 1e+00). Fix: "on a function smooth on the scale max(|x|, 1)", as `choose-a-method.md`:18 already says.

**5. Minor (evidence rule): `one-sided.md`:43-44 names sqrt(1 - x) near 1, but `deriv_one_sided.cpp`:141-151 asserts
the claim for g = sqrt(1 - x) + x.**
- Adding x leaves the derivative's shape alone, but raises |g| from 0.03 to 1.03, about 30 times more absolute
  rounding, so it changes which method wins.
- Suspected: the claim probably holds for sqrt(1 - x) itself too. Confirm by running the comparison at 0.999.
- Fix: assert the named function, or name the asserted one.

**6. Minor, verified: `one-sided.md`:50, "Avoid the one-sided stencils there: their step does not shrink near the
edge", reads as general advice.** It is justified only where fn changes fast near its edge (the asserted case is
`Backward<2>` on g at 1 - 1e-7, whose step 6e-6 is 60 times the distance to the edge).
- On the page's own function 1/(1 + x), `Forward<2>`'s points all stay inside the domain, and its error does not
  depend on the distance to the edge: about 1e-10 there as at 0.
- Fix: "Near an edge where fn changes fast, such as sqrt(x) near 0, avoid the one-sided stencils: their step does not
  adapt to fn's scale, and they have no estimate to warn you."

**7. Minor, verified: `step.md`:65-66 says `default_step` "gives the default step, which you can scale the same
way".** `default_step<M, N>(x)` returns eps^(1/(P+N)) max(|x|, 1). A reader at |x| > 1 who writes
`default_step<M, N>(x) * a` gets a step max(|x|, 1) times too large. It works on the page only because
max(|x|, 1) = 1 at x = 2e-4.
- Fix: "`default_step<M, N>(1.0)` gives eps^(1/(P+N)); multiply it by fn's scale."

**8. Minor, verified: `tutorial.md`:67-68.**
- "Far below 1e-9" is asserted only as < 1e-9. The program's own bound is 2.64e-10 + 1.3e-11, about 2.8e-10, only 3.6
  times below. Fix: "below 1e-9", or assert a tighter bound.
- "Its rounding adds at most about 3e-10": the 1.6e-15 behind it is the root solve's error (Brent's tolerance plus
  the shift), not rounding. Fix: "the solve's error adds at most about 3e-10".

**9. Minor (evidence rule): the threshold "about 1e-9 max(|x|, 1)" in the `one-sided.md`:75-76 note is asserted
nowhere.** The tests check 1e-8 (accurate) and 1e-12 (flagged); the manual checks 1e-11 (7% off). The number comes
from the header's 2^-30 = 9.3e-10.
- Fix: link the note to `nxx::deriv::ForwardRidders` in the reference, or add a check at 1e-9.
- Same class, editorial: the tutorial's "about 6e-6" default step is not asserted. It is eps^(1/3), close to an exact
  constant.

**10. Editorial, verified: the alt text of `index.md`:21 says the error falls "in a straight line ... to 4.1e-12 at
5e-6".** The truncation line h^2/32 gives 7.8e-13 at 5e-6, so the point there is 5 times above the line, already
rounding-dominated. The straight part ends near 1e-5 to 2e-5, as the text at line 28 correctly says. Fix: "falls in a
straight line to about 1e-5, reaches its lowest error, 4.1e-12, at 5e-6, then rises...".

**11. Editorial: three assertions are much looser than the claims they pin.**
- `deriv_choose_a_method.cpp`:114 asserts 1e-40 < error < 1e-30 for `cpp_bin_float_50` (eps = 2^-167 = 5.3e-51).
  The model error is eps^(2/3), about 3e-34, so the upper bound is about 3000 times the estimate; 1e-30 is 60% of the
  digits, not two thirds. Suggest < 1e-32.
- `deriv_kepler.cpp`:89 asserts 8 to 20 calls for "about a dozen".
- `deriv_higher.cpp`:101 allows 5 to 7 digits for "5 or 6" (the printed output pins these, so no harm).

## Checked and correct

- Central difference f': truncation h^2 f'''/6. For sqrt(x) at 2, f'''/f' = 3/16, so the relative error is h^2/32;
  the next term is 0.109 h^2, which moves the halving ratio by 1 + 0.082 h^2. The ratio is within 2.5% of 4 at the
  pairs (0.2, 0.1) down to (2e-4, 1e-4), with rounding about 1% at 1e-4. The table's ratios are 4.03, 3.99, 3.99 and
  4.01.
- Rounding "about eps |f| / h": the bound is eps |f| / (2h). The program's rho = 3.2e-16 (2^-52 / (2h * 0.354)) is
  sound. The stencil divides by the actual distance (`stencil_derivative`, two-node branch), as the comment assumes;
  the midpoint shift is about 4e-17.
- Default step eps^(1/3) max(|x|, 1) = 1.21e-5 at 2. The true optimum (3 eps |f| / (2 |f'''|))^(1/3) is 1.9e-5,
  consistent with design.md's omitted constant. The best error, 4.1e-12 (11.4 of 16 digits), and the model
  eps^(2/3) = 3.7e-11 both support "about two thirds". The default step's error, 6.0e-12, is "about as good". At
  1e-13 the error is 2.5e-3, "fewer than four digits".
- Ridders on sqrt(x) at 2: 5.8e-14, 103 times below the central difference; its estimate is 3.8 times the error.
- Peak: f' = -1.2e7 and f = 4000. g'''/g' = -0.56 (from Re[1/(t - i)]), so the truncation is 0.0933 (h/a)^2: 3.43e-4
  at h = 6.06e-6, which matches the output and "three or four digits". The scaled step eps^(1/3) a = 6.1e-10 gives
  truncation 3.4e-12 and a rounding bound of 1.9e-10; the observed 4.6e-11 is within 4 times the curve's lowest
  error, so "about 10 digits". The two curves' lowest-error steps are 3.7 decades apart, "about four".
- Ridders on the peak: initial step /32 with halving, as in `ridders_divisors` and the header. 13 digits, the
  estimate covers the error, and 36 and 14 evaluations are asserted.
- 1/(1 + x^2) at 0.5, via f^(n) = (-1)^n n! Im[(x - i)^-(n+1)]: f''' = 3.69, f^(5) = -21.6, f^(7) = -1110. Bounds
  3.5e-11, 3.4e-13 and 4.7e-13 match the output (3.4e-11, 2.0e-13, 4.9e-13); "about two digits" holds (factors 170
  and 69). The tolerances are 2.5 to 3.6 times the bounds.
- `cpp_bin_float_50`: sqrt(2)/4 = 0.35355339059327376220042218105242451964|2417968..., so all 40 printed digits are
  correct (no rounding tie); the estimate is 6.5e-47.
- Higher orders, 1/x at 2: steps eps^(1/(2+N)) * 2 = 1.21e-5, 2^-12, 1.48e-3 and 4.92e-3. Truncation constants
  (relative): f'' h^2/4 = 1.5e-8; f''' 1.25 h^2 = 2.74e-6; f'''' 1.25 h^2 = 3.0e-5. All match the table, and "10, 8,
  5 or 6, 4 or 5 digits" holds. In float, f'' = 0.2500806 (3.2e-4, truncation 3.5e-4), "three or four digits". The
  derivative values 0.25 and 0.75 are right.
- The warning criterion (eps > 1e-7) matches the header's "7 or fewer significant digits"; `derivative_of` avoids the
  deprecated overloads (`deriv.hpp`:1837).
- Kepler: dE/dM = 1/(1 - e cos E) and d^2E/dM^2 = -e sin E (dE/dM)^3. The propagation bound for the solve's error,
  2 * 1.6e-15 / (2h) = 2.64e-10, is sound.
- One-sided at 0: `Forward<2>` truncation -h^2 f'''/3 = +2h^2 = 7.3e-11 plus rounding (S = 4) gives the observed
  1.14e-10, "about 10 digits"; `ForwardRidders` 1.65e-13, "about 13".
- Near-edge table: the 64 evaluations match the cap `ridders_levels` (2 per level, 32 levels); both estimates cover
  their errors at each distance; the central difference fails at 1e-6; "keep the smaller estimate only if they
  agree" works on the 1e-10 row. Kink: (f(h) - f(-h))/2h = -2 + O(h), the mean of the slopes, as the header
  documents. The vanishing-estimate case is asserted.
- Handle failures: sqrt(x) at 1e-11 has derivative 1.5811e5. The result 1.4730e5 is 6.8% off, with an estimate (11%)
  that covers it. The x/8 step gives an estimate of 3.2e-12 relative. At 1e-12, 70% off with the flag raised.
- `derivative_of`: identical results bit for bit; the 2-call count; 1/(2 sqrt(x)) to 15 decimals.
- Costs table: every claim on the page matches its rows; "ten thousand times faster" is 1/a. The "10 to 16" and "13
  to 15 digits" figures match the header and the table.
- Tolerance comments in `deriv_trade_off.cpp`, `deriv_scale.cpp`, `deriv_step.cpp`, `deriv_one_sided.cpp`,
  `deriv_higher.cpp` and `deriv_handle_failures.cpp`: each a sound estimate, each tolerance within 1 to 6 times it,
  apart from item 11. Every reference is a closed form, not a Numerixx result.

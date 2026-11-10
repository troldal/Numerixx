# Phase 9, stage 7: the independent review by Codex (e70ed13)

OpenAI Codex CLI 0.166.0, model gpt-6.1-sol, reasoning effort high, read-only sandbox, 2026-11-10, on the phase from
8f565e3 to e70ed13; told not to read the earlier reviews first (it did not). Saved by the main session as Codex wrote
it, with the reproductions below.

I found one major documentation issue and two minor issues. No files were modified. I did not read `docs/sphinx/reviews/`.

1. **Major -- the Newton recommendation is too broad.**  
   [docs/sphinx/deriv/one-sided.md:52](docs/sphinx/deriv/one-sided.md:52) says `Forward<2>` or `Backward<2>` "is enough" when the derivative only steers Newton's method.

   Counterexample: \(f(x)=\sqrt{x}-10^{-6}\), starting at \(x=2\times10^{-12}\). The root is \(10^{-12}\). With its default step \(h\approx6.05545\times10^{-6}\), `Forward<2>` gives approximately **525.049**, against the true derivative **353553.391**. Newton therefore proposes **-7.86905x10^-10**, where `fn` returns NaN; Numerixx returns `non_finite`. With the exact derivative, the next point is **8.28427x10^-13**, inside the domain and across the root.

   **Correction:** qualify the recommendation: these stencils can suffice when their step resolves the function's scale. Near a rapidly changing edge, retain the advice to scale the step or use an adaptive method; a bracket also safeguards Newton.

2. **Minor -- "two thirds" is stated for all stencils.**  
   [docs/sphinx/deriv/choose-a-method.md:75](docs/sphinx/deriv/choose-a-method.md:75) says a stencil keeps about two thirds of the type's digits. That describes `Central<2>`, which is the only stencil checked by this part of its program.

   The implementation's model is \(\varepsilon^{P/(P+1)}\) for a first derivative. In the example's 50-digit type, that predicts roughly **33, 40 and 43 digits** for `Central<2>`, `<4>` and `<6>`. A reproduction using 168-bit arithmetic on the same square-root example gave relative errors approximately **5.68x10^-35, 1.10x10^-41 and 5.51x10^-44**.

   **Correction:** say "The central difference keeps about two thirds..."; higher accuracy orders retain more.

3. **Minor -- the Kepler reference is not fully independent of Numerixx.**  
   [docs/sphinx/examples/deriv_kepler.cpp:43](docs/sphinx/examples/deriv_kepler.cpp:43) obtains `E` through Numerixx's root solver, then uses it to construct both derivative references. The uncertainty discussion at line 65 makes the current numerical comparisons reasonable, but the oracle still depends on Numerixx, contrary to the evidence rule.

   Independently solving Kepler's equation with 85-digit arithmetic gives:
   
   - \(E'(1)=1.0373620218936458704990245805\ldots\)
   - \(E''(1)=-0.5567130326685877794636349315\ldots\)

   **Correction:** retain the displayed closed-form calculation for teaching, but check the results against independently recorded derivative constants, with their rounding allowance.

What I checked and found correct:

- **Both figures:** independently reproduced all **40** trade-off rows and **43** scale rows to their displayed precision. The scripts use the checked tables; captions, alt text, minima and default-step markers agree with them.
- **Near-edge table:** reproduced all six rows' displayed errors, estimates and evaluation counts, including Ridders' 64 evaluations at \(10^{-10}\).
- **The revised "about twice" warning:** fair as an empirical description of the named first-derivative samples, not a general guarantee. At \(2^{-39}\), the reproduced error/estimate ratio is **2.22365**. The new test correctly checks unflagged behaviour against the closed form, with a looser factor-three coverage check.
- **The fourth-derivative pole clause:** its failure mechanism and magnitude are correct. For example, \(a\approx4.29445668859717\times10^{-6}\) produces about **61.65 times** the correct derivative, unflagged, with error **4.627 times** its estimate. I verified representative cases, not the original nine-of-4000 sample.
- **`derivative_of`:** the revised reference-capture lifetime explanation, value-capture example, copying, `std::ref`, constness and warning behaviour agree with the implementation.
- The derivative-order syntax, default-step formulas, stencil evaluation counts, optional estimates, cancellation examples and three `basics.md` edits agree with the code.

The existing clang-cl derivative tests passed **34 cases and 573 assertions**; the Kepler example also passed. I did not rebuild or run the examples that write `.out` files; those were checked through source inspection and in-memory numerical reproduction.

## Reproductions (main session)

- Finding 1 reproduced (GCC 16): `Forward<2>` on sqrt(x) - 1e-6 (NaN below 0) at 2e-12 gives 525.05 against the
  exact 353553.39; Newton's first step goes to -7.87e-10, and `solve<Newton>` with `derivative_of<Forward<2>>` from
  2e-12 returns `non_finite`. The roots chapter's own case, ln(x / 1e-12) from 1e-13, converges to 1e-12, so its
  example is right but its sentence "From a guess, a one-sided difference avoids the problem" is as broad as ours.
- Finding 2 follows from the default step eps^(1/(P+1)) in `deriv.hpp`: a stencil of order P keeps about P/(P+1) of
  the digits, two thirds only for `Central<2>`; the costs table shows `Central<4>` beating `Central<2>` in double.
- Finding 3 reproduced: `deriv_kepler.cpp` builds its references on E from `nxx::roots`. mpmath at 60 digits,
  independently: E = 1.49870113351784831405798549726, E' = 1.0373620218936458704990245805,
  E'' = -0.556713032668587779463634931485, the same as Codex's.

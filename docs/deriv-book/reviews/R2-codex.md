# R2-codex

Independent review by OpenAI Codex (codex-cli 0.161.0, `codex exec -s read-only --ephemeral`), round 2, on commit 5432abf; the main session saved the output verbatim.

Reviewed commit **5432abf**. I made no file changes.

## Ranked findings

1. **R2-codex-01 — major — Kinks can pass both recommended accuracy checks**

   **Locations:** [06-ridders.tex:208](C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/06-ridders.tex:208), [07-using.tex:243](C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/07-using.tex:243), and the quick-reference page.

   The book says a kink produces a large Ridders estimate, and that a kink or jump “does not survive” comparison at \(h,2h,h/2\). Neither is assured.

   A counterexample compiled directly against the header with GCC 16 is
   \[
   f(t)=t+\lvert t+2^{-40}\rvert,\qquad x=0.
   \]
   The derivative exists and equals **2**. Default `derivative<Ridders>` returns:

   - value: **1.000000000067909**
   - `abs_error`: **2.1500912161798169e-10**

   Thus a result wrong by approximately 50% passes the suggested relative tolerance check at \(10^{-8}\).

   For the step-comparison claim, \(g(t)=t+|t|\) at zero gives exactly **1** with `Central<2>` at \(h=0.125,0.25,0.5\), although the derivative does not exist. Ridders also returns 1, with an estimate of approximately machine epsilon.

   **Suggested change:** Present these checks as heuristics requiring a smooth, sufficiently resolved function. Restrict the large-estimate observation to the measured nitrogen example. Add a counterexample and recommend checking branch boundaries and one-sided derivatives when smoothness is uncertain.

2. **R2-codex-02 — minor — The minimum-error multiplier is inverted**

   **Location:** [a2-derivations.tex:153](C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/a2-derivations.tex:153).

   Write the model as \(E=T+R=ah^P+b/h^N\). At its minimum,
   \[
   PT=NR,\qquad E(h_*)=(1+N/P)R.
   \]
   The book instead gives \((1+P/N)R\).

   For \(N=1,P=2\), its expression gives \(3R\), while the correct value is \(3R/2\). Chapter 3’s special-case calculation uses the correct factor.

   **Suggested change:** Replace \(1+P/N\) with \(1+N/P\). The preceding optimum-step formula and the subsequent epsilon exponent are correct.

3. **R2-codex-03 — minor — A positive error constant is used as a signed coefficient**

   **Locations:** [04-stencils.tex:39](C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/04-stencils.tex:39), [02-taylor.tex:134](C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/02-taylor.tex:134), [99-backmatter.tex:31](C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/99-backmatter.tex:31).

   The book defines
   \[
   C=\frac{|M_{N+P}|}{(N+P)!}
   \]
   and then calls \(+Ch^P f^{(N+P)}\) the leading truncation term. The signed coefficient is \(M_{N+P}/(N+P)!\).

   For `Central<4>`, \(M_5=-4\), so
   \[
   D_hf-f'=-h^4f^{(5)}/30+\cdots.
   \]
   For \(f(t)=t^5\) at zero, the rule returns \(-4h^4\), whereas the general statement predicts \(+4h^4\). The worked derivation immediately below and the generated table correctly retain the negative sign.

   **Suggested change:** Distinguish the signed coefficient from the positive magnitude \(C\), consistently across the introductory definition, derivation and glossary.

4. **R2-codex-04 — minor — Figure 7.2 applies the two-point jump formula to other methods**

   **Location:** [07-using.tex:362](C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/07-using.tex:362), Figure 7.2.

   The caption describes a plateau \(J/(2hc_p)\) whenever the break lies inside a method’s span. That coefficient applies to `Central<2>`, not generally to other rules.

   For `Central<4>` and a jump at \(a\), with \(x=a+\delta\), its weights give the jump contribution
   \[
   \begin{cases}
   7J/(12h),&0<\delta<h,\\
   -J/(12h),&h<\delta<2h.
   \end{cases}
   \]
   Both regions are inside its span. The CSV exhibits the corresponding two levels: approximately **0.04092** relative error at \(\delta=0.28117\), and **0.005843** at \(\delta=0.5\), a factor of seven apart. An adaptive Ridders calculation also has no single fixed-\(h\) plateau formula.

   **Suggested change:** Limit \(J/(2h)\) to the two-point rule. Explain the general contribution as \(J/h\) times the sum of weights on one side of the jump, and identify the second `Central<4>` plateau.

5. **R2-codex-05 — minor — The stated order of divisions omits an important branch**

   **Locations:** [05-orders.tex:212](C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/05-orders.tex:212), [08-api.tex:210](C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/08-api.tex:210).

   Both passages say the other rules divide by \(dh\), then by \(h\), \(N-1\) times. In [deriv.hpp:733](C:/Dev/XLThermo/Numerixx/include/numerixx/deriv.hpp:733), the implementation reverses that order when \(dh>1>h\): divisions by \(h\) come first, and division by \(dh\) last.

   This changes observable results near underflow. The regression in [test_orders.cpp:289](C:/Dev/XLThermo/Numerixx/tests/test_orders.cpp:289) uses `Central<10>`, \(N=4\), \(h=19\,2^{-22}\). Dividing by \(453600h\) first produces zero; the implemented order preserves a normal result near \(2.6\times10^{-308}\).

   **Suggested change:** Describe the conditional ordering. Retain the correct explanation that the implementation avoids forming \(h^N\).

6. **R2-codex-06 — minor — Appendix A generalises the two-node resolution behaviour**

   **Location:** [a1-floating-point.tex:74](C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/a1-floating-point.tex:74).

   The appendix says that when only one of \(x+h,x-h\) moves, “the result is finite but mostly rounding”. The multi-node implementation instead rejects the call if **either** side fails to move: [deriv.hpp:1028](C:/Dev/XLThermo/Numerixx/include/numerixx/deriv.hpp:1028).

   A compiled counterexample uses \(f(t)=t\), \(x=-1\), \(h=\texttt{0x1.8p-54}\). Here \(x+h\ne x\), but \(x-h=x\):

   - `Central<2>` for \(f'\) succeeds and returns exactly **1**.
   - `Central<4>` for \(f'\) returns `non_finite`.
   - `Central<2>` for \(f''\) returns `non_finite`.

   Thus neither the universal success claim nor “mostly rounding” describes this example.

   **Suggested change:** State the different resolution checks for two-node first derivatives and multi-node rules.

7. **R2-codex-07 — minor — The test references are not all independent of Numerixx**

   **Location:** [10-guarantees.tex:213](C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/10-guarantees.tex:213).

   The assertion that reference values come from closed forms, published tables or an independent prototype, “never from Numerixx itself”, is too broad.

   [test_ridders.cpp:313](C:/Dev/XLThermo/Numerixx/tests/test_ridders.cpp:313) explicitly identifies its multiprecision error-estimate references as coming from the **first C++ run**. [test_orders.cpp:632](C:/Dev/XLThermo/Numerixx/tests/test_orders.cpp:632) likewise bounds higher-derivative estimates using that first run.

   These are useful regression baselines, but they do not independently validate the magnitude of the reported error estimate.

   **Suggested change:** Distinguish independent derivative-value oracles from implementation-derived regression ceilings for `abs_error`.

8. **R2-codex-08 — minor — Continuity alone does not justify the stated next-power remainder**

   **Location:** [a2-derivations.tex:36](C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/a2-derivations.tex:36).

   The argument replaces \(f^{(N+P)}(\xi_k)\) by its limiting value at \(x\), then claims an \(O(h^{P+1})\) remainder. Continuity supplies \(o(h^P)\); the stronger remainder needs an additional regularity assumption.

   For example, \(f(t)=\operatorname{sgn}(t)|t|^{7/2}\) is \(C^3\), with \(f'''(0)=0\). Its central first difference at zero is \(h^{5/2}\). With \(N=1,P=2\), the displayed leading term vanishes, but the remaining error is **not** \(O(h^3)\).

   **Suggested change:** Explicitly assume a locally Lipschitz \(f^{(N+P)}\), or a continuous additional derivative, for this statement. Otherwise write \(o(h^P)\). The appendix’s general opening assumption should identify which regularity this step uses.

9. **R2-codex-09 — minor — The Sterbenz proof asserts a false factor-of-two condition**

   **Location:** [05-orders.tex:205](C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/05-orders.tex:205).

   The passage claims that, for every \(0<h\le|x|\), both pairs \((x_+,x)\) and \((x,x_-)\) have the same sign and lie within a factor of two.

   Even with exactly represented binary values, \(x=1,h=0.75\) gives \(x_-=0.25\): the second pair differs by a factor of four. At \(h=x\), one point is zero. The cited hypothesis therefore does not establish the claimed interval.

   **Suggested change:** Supply a proof covering the remaining interval separately, including signs and endpoints, or restrict the displayed Sterbenz argument to a range where its hypotheses hold. This finding concerns the argument; it does not invalidate the binary actual-offset technique itself.

10. **R2-codex-10 — editorial — “Every rule is exact in float” overstates exact coefficient storage**

    **Location:** [09-tables.tex:77](C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/09-tables.tex:77), Table 9.4 caption.

    Integer numerators fitting the significand means those integers are represented exactly. It does not make the finite-difference rule exact: rational weights such as \(1/6\), function evaluations, products and divisions still round. The book’s own float accuracy measurements demonstrate this distinction.

    **Suggested change:** Say that all stored integer numerators and denominators are represented exactly in `float`, avoiding coefficient-conversion error.

## Checked and found correct

- **All 44 stencil-table rows:** independently checked integer moment conditions, first nonzero error moments, signs and error constants. Compile-time comparisons also confirmed their nodes, integer numerators and denominators against the header.
- **Public API:** names, overloads, supported orders, constraints, return types, `std::expected` usage, explicit-step preconditions, precision warning and user-defined-type requirements agree with the implementation, apart from the specific descriptions reported above.
- **Default steps and costs:** stencil exponents and scaling, Ridders’ initial-step divisors, the 32-level budget, and maximum evaluation counts **64, 65, 128, 129** agree with the code.
- **Ridders’ principal description:** the Richardson recurrence, rounding floors, two-level stall, final check, restart logic and shared budget describe the library implementation rather than a generic algorithm.
- **Generated evidence:** the coverage CSVs reproduce the reported statistics. Across 2,700 cases per derivative order, the uncovered counts are **9, 4, 1, 0**; median evaluation counts are **14, 17, 32, 29**. I found no additional mismatch between the quoted generated values and the data they reference.
- **Engineering examples:** the Shomate units and coefficients match the [NIST source](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Mask=1). The enthalpy derivatives, implicit Colebrook differentiation, sensitivities and Jacobian calculations are consistent.
- **Noise analysis:** the central-difference noise scaling and the distinction between a simplified noise-estimation recipe and the published method are sound. The cited rate is supported by [Moré and Wild’s paper](https://www.mcs.anl.gov/papers/P1785.pdf).
- **Executable checks:** all **24 current example executables** passed. The library executable passed **24 test cases and 430 assertions**. `tools/check.py` reported **0 errors and 0 warnings**.
- **PDF:** read the extracted text of the 88-page build and visually inspected selected plots, dense reference tables and the quick-reference page. Those inspected pages were readable and cleanly laid out.

These execution checks used the existing binaries; I did not rebuild or regenerate files. Full-text citation verification remains incomplete for inaccessible sources, particularly the original Ridders article; I could verify its metadata and abstract, but not certify the paper-specific step-ratio and stopping-rule attribution.

## Verdict

The book has a strong technical foundation: the coefficient tables, API coverage, engineering examples and generated evidence largely withstand independent checking. It needs revision before serving as a dependable engineering reference, chiefly because its practical guidance sometimes turns useful diagnostics into assurances they cannot provide. Correcting that major issue, the optimum-error multiplier, and the smaller proof and implementation-description discrepancies would materially strengthen an otherwise careful and useful book.
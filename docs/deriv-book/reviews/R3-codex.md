# R3-codex

Independent review by OpenAI Codex (`codex exec -s read-only --ephemeral`) of the Phase 5 update (uncommitted changes on 56d27a1); saved verbatim by the main session.

**Ranked findings**

1. **R3-codex-01 — major — The pole claim conflates numerical reach with the domain-edge failure boundary.**  
   Location: [10-guarantees.tex:177](/C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/10-guarantees.tex:177), also [06-ridders.tex:198](/C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/06-ridders.tex:198) and [qr-quick-reference.tex:48](/C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/qr-quick-reference.tex:48).

   “The limit holds for a pole … too” extends the preceding boundary for returning a value to finite poles. These have different failure behaviour: [deriv.hpp:1583](/C:/Dev/XLThermo/Numerixx/include/numerixx/deriv.hpp:1583) returns the last untrusted result when the budget expires. Consequently, a pole closer than \(2^{-34}s\) can still produce a successful result.

   For example, an in-memory evaluation of the header’s first-derivative arithmetic for \(f(t)=1/(t-2^{-100})\) at zero returns approximately \(3.68935\times10^{20}\), with estimate \(5.90296\times10^{21}\); the exact derivative is \(-2^{200}\), giving relative error effectively 1. The same calculation reproduces the stored central-log reach readings. The generator’s eight pole cases use only distances \(10^{-3}\) and \(10^{-4}\), so they do not establish the claimed boundary.

   **Suggested change:** retain the measured pole improvement, but distinguish the log/NaN success boundary from finite-pole accuracy and fallback behaviour. Attribute the pole-limit statement to the design record and qualify its meaning.

2. **R3-codex-02 — minor — The backward construction omits the odd-order sign correction.**  
   Location: [06-ridders.tex:281](/C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/06-ridders.tex:281).

   Applying ForwardRidders to \(g(t)=f(-t)\) at \(-x\) estimates \((-1)^N f^{(N)}(x)\). BackwardRidders additionally negates the returned value for odd \(N\), explicitly implemented at [deriv.hpp:1640](/C:/Dev/XLThermo/Numerixx/include/numerixx/deriv.hpp:1640) and checked by [test_onesided.cpp:328](/C:/Dev/XLThermo/Numerixx/tests/test_onesided.cpp:328). Without that step, the description gives the wrong sign for \(f(t)=t\).

   **Suggested change:** add “then negates the value for odd derivative orders; the error estimate is unchanged.”

3. **R3-codex-03 — minor — The method-selection guidance needs a smoothness qualification.**  
   Location: [06-ridders.tex:299](/C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/06-ridders.tex:299).

   “Where \(f\) exists on both sides, use Ridders” conflicts with the section’s stated use of one-sided methods at a kink. Both branches of the test function at [test_onesided.cpp:153](/C:/Dev/XLThermo/Numerixx/tests/test_onesided.cpp:153) exist: its left and right slopes at 1 are 2 and 3, and the tests verify that BackwardRidders and ForwardRidders return those respective slopes. Central Ridders instead targets their mean.

   Related unchanged statements in [10-guarantees.tex:198](/C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/10-guarantees.tex:198) and [qr-quick-reference.tex:85](/C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/qr-quick-reference.tex:85) say a nearby break gives a wrong value from “every method”, although a one-sided method sampling the desired smooth branch can be correct. The updated Table 7.1 already qualifies this as sampling “across the break”.

   **Suggested change:** recommend central Ridders when the desired branch is smooth and valid on both sides; restrict the blanket break warning to methods sampling across it.

4. **R3-codex-04 — minor — The arithmetic reference omits the one-sided second-derivative exception.**  
   Location: [08-api.tex:236](/C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/08-api.tex:236).

   The paragraph identifies `Central<2>` for \(f''\) as using actual offsets, then assigns nominal-step division to “the other rules”. ForwardRidders for \(N=2\) also uses actual offsets: [deriv.hpp:1344](/C:/Dev/XLThermo/Numerixx/include/numerixx/deriv.hpp:1344) passes the actual adjacent spacings to `second_difference_of`; BackwardRidders inherits this through reflection.

   **Suggested change:** explicitly include the second-derivative levels of both one-sided extrapolation methods among the actual-offset cases.

**Checked and found correct**

- The \(2^j-1\) divisor, one-order gain per column, alternating extrapolation coefficients, and carried-floor recurrence in §6.7 and Appendix B.3 are correct.
- Method names and declarations, default initial step \(s/32\), cached \(f(x)\), one-sided sampling through restarts, and maxima \(33,65,97,129\) agree with the implementation and tests. The revised restart descriptions agree with the code.
- Warning scope, epsilon threshold, quoted message, and suppression mechanisms agree with the header. The float table markings and Figure 5.2 annotation match the expanded scope.
- No incorrect macro selection was found in the newly quoted numerical values. The generator computes the stated error ratios, coverage, costs, and pole references as described.
- All four reach readings recomputed from the CSV match their macros:

  | Derivative order | First full budget | Covered through \(k\) | Flagged from \(k\), through 45 |
  |---|---:|---:|---:|
  | 1 | 25 | 36 | 38 |
  | 2 | 25 | 34 | 35 |
  | 3 | 26 | 33 | 34 |
  | 4 | 26 | 32 | 33 |

  The sole uncovered, unflagged point is correctly reported at \(k=37\): error \(0.2242128743\), estimate \(0.1931093837\), relative to the exact derivative. Figure 6.5 acknowledges this exception; its central-log comparison and cutoff match the CSV.
- The new nitrogen example reproduces the quoted errors, estimates, and evaluation counts. Its exact initial step is 15.625 K, correctly described as “about 16 K”.
- Existing binaries passed all 35 library test cases, all 25 book examples, and the quickstart. `check.py` reported zero errors and warnings against the existing build log. I did not rebuild or regenerate data, and modified no files.

**Verdict**

The update’s core mathematics, implementation descriptions, generated evidence, and new reach figure are sound. I recommend revision before acceptance: clarify the pole-limit claim, include the backward sign correction and actual-offset exception, and reconcile the smoothness guidance with the newly available one-sided methods. These are local explanatory fixes; I found no numerical-data defect requiring changes to the new example or reach dataset.
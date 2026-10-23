# R1-numerics: review of the numerical and mathematical claims (round 1)

Reviewer: `book-reviewer-numerics` (Fable), 2026-10-23, on the integrated build (commit `aeec7e6`, 83-page PDF).

Read: every chapter source in `chapters/` (01 to 10, a1, a2, qr, 00, 99, `shared-method-summary.tex`), `macros.tex`,
`PLAN.md` (sections 1 to 12), the PDF (`build/gcc/docs/deriv-book/pdf/main.pdf`, all figures and tables), the ground
truth `include/numerixx/deriv.hpp` and `math.hpp`, `tests/test_ridders.cpp`, `test_stencils.cpp`, `test_orders.cpp`,
`test_math.cpp` (the root_epsilon claims), `docs/design.md`, `PHASE2_API.md` (the cap rationale), `PHASE3_API.md`,
`PHASE4_API.md`, `tools/generate_data.cpp`, every file in `data/`, every figure in `figures/`, the shared headers and
the examples in `code/`.

Recomputed (scratchpad `r1num/checks.py` and `tail.py`, Python with `fractions`): every weight, moment, remainder
constant and remainder sign the text derives, in exact rational arithmetic; the Shomate thread (H, cp, H'', H''' at
1000 K; the jumps at the 500 K junction) from the published coefficients; the error model's best step, the constant
kappa and the cost of dropping it for every rule the book quotes; the step-sensitivity factors; the model gain of the
orders beyond the caps; the generator's largest intermediate; the statistics of the random study (coverage, uncovered
cases, median and worst ratio, zero-error cases, restarts, evaluation counts, the histogram) from the four
`ridders-coverage-N*.csv`; the minima of the measured V curves; the model-to-measured ratios of Table 10.1; Newton's
convergence ratios; the edge figure's evaluation counts and the 2^-34 derivation.

Severity as in `reviews/README.md`: major = wrong or misleading; minor = imprecise or unsupported; editorial.
Findings are ranked, most severe first.

## Findings

### R1-numerics-01 (major). The sign of the remainder in Eq. (4.3)

Location: `chapters/04-stencils.tex` line 47 (Eq. 4.3, PDF page 21); referred to from 4.3 (line 180) and Appendix
B.1 (`a2-derivations.tex` line 38).

Problem: Eq. (4.3) reads f'(x) = [f(x-2h) - 8f(x-h) + 8f(x+h) - f(x+2h)]/(12h) - (h^4/30) f^(5)(xi). The sign of
the remainder is wrong; it is + (h^4/30) f^(5)(xi).

Evidence. (i) The chapter's own derivation one line later: M_5 = -4, so the weighted sum is h f' + (M_5/5!) h^5
f^(5) = h f' - (h^5/30) f^(5); divided by h, the rule is f' - (h^4/30) f^(5), hence f' = rule + (h^4/30) f^(5).
(ii) Exact check on f = x^5 at x = 0 with h = 1 (`checks.py`, section 3): the rule gives -4 while f'(0) = 0 and
h^4 f^(5)/30 = 4, so f' = rule + h^4 f^(5)/30. (iii) Table 9.1's own convention (9.1: "the rule's value minus
f^(N)(x) is about +-C h^P f^(N+P), with the sign given in the last column") lists -(1/30) h^4 f^(5) for P = 4,
that is rule - f' = -h^4 f^(5)/30: the opposite of Eq. (4.3). (iv) `tests/test_stencils.cpp` line 123 asserts
`derivative<Central<4>>(x^5, 1, 0.25) == 5 - 4/256`, that is rule = f' - h^4 f^(5)/30. (v) Burden et al., the
five-point midpoint formula: f'(x0) = [...]/(12h) + (h^4/30) f^(5)(xi). Eq. (5.1), on the facing page, is correct
(rule - f'' = + h^2 f^(4)/12, so f'' = rule - h^2 f^(4)/12, as Table 9.2 says), which makes the inconsistency
visible to a careful reader.

Correction: change the minus to a plus before h^4/30 in Eq. (4.3). Chapter 2 makes a point of the sign of the
error term (2.2: "the central difference overestimates f' where f''' > 0"), so this is not cosmetic. Appendix B.1's
sentence "That is the form in Eq. (4.3)" is also affected by R1-numerics-04.

### R1-numerics-02 (major). The tail of the N = 2 evaluation histogram is not the restarts

Location: `chapters/06-ridders.tex` lines 207-209 (caption of Fig. 6.5, page 36) and lines 216-217 (6.5: "the
tail of the second-derivative histogram is the restarts of log x and sqrt(x) at points below 0.25"); the label
"restarts" in `figures/06-ridders-cost.tex` line 405.

Problem: the tail is attributed to the 11 restarts; the data say otherwise.

Evidence (`tail.py`, from `data/ridders-coverage-N2.csv`): the histogram of evaluations for N = 2 is 7 (4 cases),
9 (3), 11 (293), 15 (531), 17 (1537), 19 (288), 21 (44). The bar at 19 holds 288 cases, none of them restarted,
from e^(-x^2) (117), the Runge function (57), 1/x (50), log x (35) and sqrt(x) (29), at x from 0.23 to 2.99. The
bar at 21 holds 44 cases of which 11 restarted (every restart takes exactly 21 evaluations). So the tail of 332
cases is the ordinary, slower convergence of five functions; the restarts are a thirtieth of it and invisible at
the figure's scale. Sidebar 6.2's statement "11 for f'', which are the cases of log x and sqrt(x) below x = 0.25"
is correct; only the attribution of the histogram's tail is wrong.

Correction: caption and text: "the bars at 19 and 21 are the slower convergence of e^(-x^2), the Runge function,
1/x, log x and sqrt(x), whose higher derivatives are large relative to f''; the 11 restarts of log x and sqrt(x)
below x = 0.25 all take 21 evaluations and are a small part of the last bar". Move or remove the "restarts" label
in the figure.

### R1-numerics-03 (major). "13 and 8 exact stops, all on x^5"

Location: `chapters/06-ridders.tex` lines 292-293 (Sidebar 6.2) and lines 135-136 (caption of Fig. 6.3: "An
actual error of exactly 0 (an exact stop; 13 and 8 cases)").

Problem: the counts are of cases whose returned value equals the closed form bit for bit (`actual == 0`,
`generate_data.cpp` line 1027), not of exact stops (the stopping rule "an entry equals its parent", which the
generator does not log). For f' the 13 cases lie on seven functions, not on x^5, and their evaluation counts show
ordinary stalls, not early exact stops.

Evidence (`tail.py`): N = 1 zero-error cases (function, x, evaluations): exp (3.0, 14); sin (0.339, 12); log
(0.996, 10), (0.848, 14), (0.659, 14), (0.992, 14); 1/x (0.29, 16); Runge (0.601, 14); x^5 (0.525, 12), (1.095,
12), (0.43, 12); e^(-x^2) (2.562, 16), (2.498, 14). An exact stop on a polynomial costs 4 to 8 evaluations
(`tests/test_ridders.cpp` lines 200-206); 10 to 16 evaluations are the usual two-level stall. For N = 2 all 8
cases are on x^5 (11 evaluations, one with 9), consistent with exact stops, but still not logged as such.

Correction: Sidebar 6.2: "13 and 8 cases with an actual error of exactly 0: for f'' all on x^5, for f' a bit-exact
coincidence spread over seven functions". Fig. 6.3's caption: drop "an exact stop" or write "a bit-exact result".
If the book wants to count exact stops, the generator must log the stop reason.

### R1-numerics-04 (minor). Appendix B.1's single-point remainder does not cover Eq. (4.3)

Location: `chapters/a2-derivations.tex` lines 34-38 (page 68).

Problem: the text says the intermediate-value theorem collapses the remainder to one point "if the products
w_k n_k^(N+P) all have one sign" and that "that is the form in Eq. (4.3) and Eq. (5.1)". For Central<4> the
products are -8/3, 2/3, 2/3, -8/3 (`checks.py`, section 3): mixed signs, so the argument does not deliver (4.3).
The same holds for every central rule for f' above order 2 (Central<6>: + - + + - +; Central<8> and Central<10>
likewise). It does deliver (5.1) (products 1, 0, 1). The last sentence of the paragraph ("When the products do not
share a sign, as for some one-sided rules") therefore misstates where the hypothesis fails.

Correction: the single-xi form of (4.3) is true, but by the interpolation error formula, not by this argument:
f'(x0) = p'(x0) + f^(5)(xi) prod_(j != 0) (x0 - x_j) / 5!, where the derivative of f^(5)(xi(x)) drops out because
the node product vanishes at x0, and prod_(j != 0) (x0 - x_j) = 4h^4 gives + h^4 f^(5)(xi)/30 (Burden et al., the
derivation of the five-point formulas). Either cite that and say that the Taylor argument gives the asymptotic form
only, or state (4.3) with the asymptotic remainder C h^4 f^(5)(x) + O(h^6); and say "for the central rules above
order 2 and for some one-sided rules" in the last sentence.

### R1-numerics-05 (minor). "The double column agrees with [the model] to within about a factor of ten"

Location: `chapters/10-guarantees.tex` lines 36-38 (page 61).

Problem: not true of the one-sided rules at the higher derivative orders.

Evidence (`checks.py`, section 9, from `data/accuracy-matrix.csv`): model eps^(P/(P+N)) over the measured median:
Forward<4> at N = 3: 1.1e-9 / 3.6e-8 = 0.03 (the measured error is 33 times the model); at N = 4: 1.5e-8 / 5.5e-7
= 0.03 (37 times); Forward<2> at N = 4: 0.17 (6 times); Central<6> at N = 4: 0.18 (5.6 times). The central rules
for f' to f''' lie within a factor 6.3 of the model, on the good side. The excess of the one-sided rules is their
amplification (S = 208 and 709 for Forward<4> at N = 3 and 4, Table 9.8), which the model error drops by
construction (9.4 says so).

Correction: "for the central rules the double column lies within a factor of about six of the model for f' to
f''' and up to six times above it for f''''; the one-sided rules lie above it by up to their amplification S, 33
to 37 times for Forward<4> at N = 3 and 4".

### R1-numerics-06 (minor). The step-sensitivity factors 2^N

Location: `chapters/05-orders.tex` lines 104-108 (Sidebar 5.1) and 141-143 (5.3); `03-rounding.tex` line 277.

Problem: "misjudging that scale by a factor 2 costs up to 2^N in error: a factor 2 for f' with Central<2> and up
to 16 for f''''". 2^N is the growth of the rounding term alone when the step is halved; when the step is doubled
the truncation term grows by 2^P, which for f' with Central<2> is 4, not 2, and which 3.5 states ("a factor of two
in w costs at most a factor of four in the error of Central<2>"). The two chapters contradict each other. By the
model with both terms at their balance, the total error grows by at most 1.67 (f') and 5.5 (f'''') for Central<2>
(`checks.py`, section 6: E(2h*)/E(h*) = (N 2^P + P 2^-N)/(P+N) and E(h*/2)/E(h*) = (N 2^-P + P 2^N)/(P+N);
for Central<4> and f'''' both are 8.0).

Correction: "halving the step multiplies the rounding term by 2^N and doubling it multiplies the truncation term
by 2^P; at the model's balance the total error grows by at most 1.7 for f' and 5.5 for f'''' with Central<2> (the
header's 16 is the growth of the rounding term alone)". Make 3.5 and Sidebar 5.1 agree.

### R1-numerics-07 (minor). The rationales of the order caps are not supported by the model

Location: `chapters/10-guarantees.tex` lines 126-129 (10.4); `chapters/04-stencils.tex` lines 235-236 (4.4);
`chapters/09-tables.tex` lines 82-86 (9.3).

Problem. (a) "Central<P> stops at P = 10: beyond it the gain is under half a digit even for 50-digit types." By the
book's own model the gain of Central<12> over Central<10> is 50.3 (12/13 - 10/11) = 0.70 digit for f' (0.68 with
the constants S and C in the model minimum) and 2.4 digits for f'''' (`checks.py`, section 4). (b) "The one-sided
rules go up to order 6; beyond that the amplification ... would cost more than the order gains": by the model
Forward<7> still gains 14 % (E_min 6.4e-13 against 7.5e-13 for Forward<6>, with S = 46 against 27.7) and
Forward<8> gains nothing; neither costs more. Both sentences repeat the header's comments (`deriv.hpp` lines 105
and 115, `PHASE2_API.md` line 115) with the phrasing of a proof (style decision 7), and nothing in the book or the
design record measures them.

Correction: present the caps as design decisions with the model's figures: "beyond P = 10 the model gains 0.7
digit per two orders for f', even at 50 digits, and more for the higher derivatives, and the library stops there";
"beyond P = 6 the amplification grows about as fast as the order pays: by the model order 7 gains a seventh of a
decade and order 8 nothing".

### R1-numerics-08 (minor). The "factor 4" of dividing by the outer-node distance

Location: `chapters/04-stencils.tex` lines 175-176 ("dividing by the outer distance would correct the outer nodes
only and enlarge the error of the others, by a factor 4 for Central<4>"); `chapters/03-rounding.tex` lines 102-104
(Sidebar 3.1: "would correct the rounding of those two points only and, for Central<4>, make the error four times
larger").

Problem: both describe the effect wrongly. With the points rounded by delta_k, the rule's point-rounding error is
(f'/h) (delta_-2 - 8 delta_-1 + 8 delta_1 - delta_2)/12. Dividing instead by 3 (x_2 - x_-2) = 12h (1 + (delta_2 -
delta_-2)/(4h)) gives f' + (f'/h) [4 (delta_-2 - delta_2)/12 + 8 (delta_1 - delta_-1)/12]: the outer nodes' term is
quadrupled (its coefficient rises from 1/12 to 4/12), not corrected, and the inner nodes' term is unchanged, not
enlarged. The design record's "4x" is right; the book's two paraphrases are not.

Correction: "dividing by the distance between the outer nodes would not remove their rounding but quadruple it
(their coefficient rises from 1/12 to 4/12) and leave the inner nodes' untouched; only for the two-node rules does
the actual distance remove the point rounding".

### R1-numerics-09 (minor). "No singularity within h_0" is false for N >= 2

Location: `chapters/06-ridders.tex` lines 131-132 (caption of Fig. 6.3) and 148-149 (6.3);
`chapters/10-guarantees.tex` lines 56-57 (10.2).

Problem: the scope statement is made for the whole study. For N = 2 the initial step is s/4 = 0.25, so log x,
sqrt(x) and 1/x at x < 0.25 have their singularity within h_0; for N = 3 and 4 the levels reach x +- 2h_0 =
x +- 0.5. These are exactly the 11 (N = 2) and 61 (N = 3, 4) restarts the sidebar reports (`tail.py`: the restarted
cases are log x and sqrt(x) with x up to 0.247 for N = 2 and up to 0.495 for N = 3 and 4). Sidebar 6.2 qualifies
"for f'"; the caption and 10.2 do not.

Correction: "no singularity within the initial step for f'; for f'' to f'''' the initial step s/4 reaches past zero
for x below 0.25 (0.5 for the span of f''' and f''''), which are the restarts of Table 10.2".

### R1-numerics-10 (minor). The floor formula in Sidebar 6.2 divides by d twice

Location: `chapters/06-ridders.tex` line 261.

Problem: the rounding of the points is written (1/2) eps (|x| + 2h) |f'| sum_k |w_k| / (d h^N). With w_k = c_k/d
(the notation of 4.1 and the back matter) this is sum_k |c_k| / (d^2 h^N). The code uses the integer numerators:
`spread = (eps/2) (|x| + 2h) * sum |c_k|`, divided by `d h` and then by h, N - 1 times (`deriv.hpp` lines
1151-1163).

Correction: sum_k |c_k| / (d h^N), or sum_k |w_k| / h^N.

### R1-numerics-11 (minor). The implemented method is not Ridders' published method; say so

Location: `chapters/06-ridders.tex` lines 7-11 (introduction) and 6.1; Sidebar 6.2, "Where the constants come
from" (lines 279-284).

Problem: the chapter attributes the method to Ridders (1982) and the sidebar says the constants "were chosen with
a prototype", but the book never says how the implementation differs from the paper, which a reader who knows
Ridders' routine (or its Numerical Recipes descendant) will ask: Ridders halves by 1.4, uses the product form of
the tableau, stops when the error grows, and has no rounding floor, final check, restart or default step
(`PHASE3_API.md` lines 167-169). In the book "Ridders' method" means Ridders' idea with the library's rules, and
every statistic quoted is of the latter.

Correction: one sentence in 6.1 or at the head of Sidebar 6.2: "Ridders' paper halves the step by 1.4 and stops
when the error grows; the library uses the ratio 2, a rounding floor, a two-level stall, a final check and
restarts, chosen on a prototype (Sidebar 6.2): the idea is Ridders', the rules are the library's."

### R1-numerics-12 (minor). Where the estimate fails is not said

Location: `chapters/06-ridders.tex` 6.3 lines 140-149 and Sidebar 6.2.

Problem: the study reports 9 and 4 uncovered cases but not what they are. They are all e^(-x^2) at x between 1.6
and 2.9 (ratios 1.01 to 1.81 for f', 1.06 to 1.36 for f''; `checks.py`, section 1), and the single one for f''' is
the Runge function at x = 0.54 (ratio 1.001). In the Gaussian tail the derivatives grow geometrically (f^(k) is
about (2x)^k f), so the even expansion is only marginally asymptotic at h_0 = s/32 (a_4 h_0^2 / a_2 is about
(2x)^2 h_0^2 / 20, about 0.2 at x = 2.5), and the change from the parent undershoots the entry's error. The design
record's exp(10x) case for f''' (96.3 % covered before the point-rounding floor) is the same mechanism. A
numerical analyst expects the book to name the class of function that makes the estimate tight, because it is
what the reader should watch for (a steep exponential, a rate with a large activation energy, a Gaussian tail) and
it is the complement of the noise and kink cases the chapter does name.

Correction: a sentence in 6.3 or Sidebar 6.2: "every uncovered case of the study is e^(-x^2) at x above 1.6, where
the derivatives grow like (2x)^k and the first levels are not yet in the asymptotic regime; a function whose
derivatives grow geometrically on the scale of h_0 is where the estimate is tightest, and a smaller initial step
restores the margin".

### R1-numerics-13 (minor). The rounding of an f that is computed by cancellation

Location: `chapters/03-rounding.tex` 3.1 lines 25-27 and 3.6.

Problem: the model takes each value of f as off by about eps |f|. When f is itself computed by cancellation (a
difference of large terms, an enthalpy of formation plus a sensible-heat integral, an alternating sum), its
rounding error is eps times the largest intermediate, k eps |f| with k >> 1, and the floor of the V rises by k.
3.1 says "sometimes less, sometimes more" and 3.3 discusses a larger |H|, but the intermediate case (smooth, not
noisy, yet off by far more than eps |f|) is never named, although it is common in property models. The remedy is
3.6's with eta = k eps.

Correction: one sentence in 3.1 or 3.6: "if f is computed by cancellation its error is eps times the largest
intermediate, k eps |f| with k >> 1; treat k eps as the noise level eta of Section 3.6".

### R1-numerics-14 (minor). 4^j - 1 rounds in binary too

Location: `chapters/06-ridders.tex` lines 250-251 (Sidebar 6.1: "The factors 4^j are built by repeated
multiplication in the real type, exact in binary; in a decimal type 4^j - 1 rounds at large j").

Problem: 4^j is exact in binary, but 4^j - 1 needs 2j bits: it rounds in float from j = 13 (26 bits > 24) and in
double from j = 27 (54 bits > 53), both inside the 32-level cap. `PHASE3_API.md` line 153 notes the float case. The
conclusion ("by at most half an epsilon, which does not matter") holds in every type.

Correction: "4^j is exact; 4^j - 1 rounds once 2j exceeds the significand (j >= 13 in float, j >= 27 in double, far
beyond the cap in a 50-digit type), by at most half an epsilon, which does not matter".

### R1-numerics-15 (editorial). The factor 4^j in 6.3

Location: `chapters/06-ridders.tex` lines 108-109.

Problem: "the parent's error is larger than the entry's by about a factor 4^j, so the difference between them is
nearly the parent's whole error". The factor 4^j relates the parent A_(i-1,j-1) to the entry's left neighbour
A_(i,j-1) (Appendix B.3); the entry A_(i,j) itself is of a higher order in h. The conclusion is right.

Correction: "the parent's error is 4^j times that of the entry's left neighbour, and the entry's own error is
smaller still, so the difference between entry and parent is nearly the parent's whole error".

### R1-numerics-16 (editorial). The 24 against 14 evaluations at the edge

Location: `chapters/06-ridders.tex` lines 159-161.

Problem: "after 24 evaluations against 14 in the interior: a restart costs about two evaluations per halving". Of
the 10 extra evaluations only 4 are the two failed levels; the converging tableau takes 10 levels (20 evaluations)
against 7 at T_0, because log x at 0.01 varies on the scale x and the restart step 1/128 is 0.8 x. The per-halving
cost of two is right and is what Fig. 6.4 measures (`checks.py`, section 11: two evaluations per unit of k, and
nine levels after the k - 4 failing ones at every k, which also gives the cap at k = 27 and the limit k = 34).

Correction: "after 24 evaluations: four for the two failed levels and twenty for a tableau that starts at 1/128,
nearly the size of x, and needs ten levels; Fig. 6.4 shows that each further halving of x adds one failed level,
two evaluations".

### R1-numerics-17 (editorial). Small slips

(a) `03-rounding.tex` line 65: 2.71855 and 2.71801 "share their first three digits": four (2.718), which with six
digits leaves the two the sentence counts.

(b) `04-stencils.tex` lines 102-103 (Sidebar 4.1): "The largest intermediate is 12! * 4!": the largest intermediate
of the generator for Central<10> and N = 4 is 12! = 4.8e8 (the end node's product of node differences); 12! * 4! =
1.1e10 is a bound that is not reached (`checks.py`, section 5). "Far inside long long" stands.

(c) `06-ridders.tex` line 90: "needs 13 levels to reach its best": level 13 is the fourteenth level; "13 halvings".

(d) `04-stencils.tex` lines 116-117 (caption of Fig. 4.1): "the minimum falls from about eps^(2/3) to about
eps^(10/11), three and a half decades": the two model values are 3.8 decades apart (3.7e-11 to 5.9e-15); 3.5
decades is the gap of the measured minima (2.6e-12 to 8.2e-16, `checks.py`, section 8). Say which is meant.

(e) `10-guarantees.tex` Table 10.2 and Sidebar 6.2: f''' shows coverage "100.0" with 1 uncovered and a worst ratio
of "1.0" (it is 1.001, `checks.py`, section 1); print 99.96, or say "one case 0.1 % over its estimate".

(f) `a1-floating-point.tex` line 83: Sterbenz's lemma needs a and b of the same sign (a/2 <= b <= 2a with a, b > 0);
"two nearby numbers" implies it, say it.

(g) `03-rounding.tex` lines 322-324: "their spread estimates eta |f|": the k-th difference of independent noise has
a spread sqrt(C(2k, k)) times the noise (about 1.4, 2.4 and 4.5 times for k = 1, 2, 3); More and Wild divide by
that factor. "Estimates eta |f| up to a known factor".

(h) `05-orders.tex` line 68 against line 72, `10-guarantees.tex` line 39 and `01-why.tex` line 168: the model's
eps^(1/3) = 6.1e-6 is called "five digits" (5.2) and the measured 1.2e-5 "4 digits" by the generator's floor (4.9).
Either "about five" for both, or say that the headline figure floors.

(i) `06-ridders.tex` line 199 and `07-using.tex` line 322: "eleven orders above its usual value": 0.0062 against
2.7e-13 (the estimate 50 K from the break) is 10.4 orders; against the thread's 8.2e-14 at T_0 it is 10.9. "Ten
orders" for the first comparison.

(j) `06-ridders.tex` line 221: "Ridders costs three to seven times as much for f'": 14 against 2 to 6 evaluations
is 2.3 to 7 times.

(k) `06-ridders.tex` line 164 and `07-using.tex` line 340: the limit 2^-34 is where a value is still returned,
with 0.4 % relative error; the text says so two clauses earlier, but the quick-reference and 10.4 readers should
not take "needs the nearest bad point 2^-34 s away" for ten-digit accuracy. Consider "about 2^-30 s for full
accuracy, 2^-34 s for any value at all".

## Checked and correct

The following were verified and need not be re-checked in round 2 unless the text changes.

Chapter 2. Eqs. (2.1) to (2.4) and the derivation in Sidebar 2.1 (the mean of f'''(xi+) and f'''(xi-) and the
intermediate-value step; the even series). The backward difference's error -h f''/2 and the central difference as
the mean of the two. The signs of H'' (+0.005744) and H''' (-5.8e-6) at 1000 K, recomputed from the coefficients,
and the stated over- and under-estimation. The x^3 hand table (exact: 3.31/0.31, 3.01/0.01, 3.1525/0.1525,
3.0025/0.0025; forward error 3h + h^2, central h^2). Exactness on quadratics; the degree P + N - 1. The Richardson
combination (4 D(h/2) - D(h))/3 = f' + O(h^4). Fig. 2.1's slopes 1 and 2 (four and eight decades over four
decades of h in `step-error-cp.csv`). The definition of O(h^P) and the factor-of-10 statement (sqrt(10) = 3.2).

Chapter 3. The six-digit calculator table recomputed (errors +0.0046, +0.00022, +0.0017, -0.018; ratio 21).
Eq. (3.1) and its derivation (2 eps |f| over 2h). Eq. (3.2): h* = (3 eps |f| / |f'''|)^(1/3) and E(h*) = (3/2)
eps |f| / h*; with the scale assumption E(h*)/|f'| = 1.04 eps^(2/3). eps^(1/3) = 6.1e-6, eps^(2/3) = 3.7e-11,
eps^(1/2) = 1.5e-8. The thread's factors |H|/(s |cp|) = 0.66 and |H|/(s^2 |H''|) = 3.7 (exact: 0.6565, 3.736).
Sidebar 3.1: the distance error eps |x| and the relative error eps |x| / (2h) = eps^(2/3)/2 at the default step;
Sterbenz for h <= |x|/3; the single rounding above that; the code divides by the evaluated distance (`difference`).
Sidebar 3.2: h*, kappa = (N S / (P C))^(1/(P+N)), the cost formula (N kappa^-P + P kappa^N)/(P+N) (both terms at
the ratio N : P), 7 * 3^(1/3)/9 = 1.12; kappa lies in [1.44, 2.00] over the first-derivative rules with Central<2>
and Forward<1> as the extremes, and the cost peaks at 1.634 for Central<10>; the higher-order kappas 2.63 (f''),
2.40 (f'''') and the costs 3.54, 11.2, 15.1 (`checks.py`, section 7, all agree with `values.tex`). Sidebar 3.3:
ulp(x) = eps 2^e, ulp(1000) = 1.1e-13, the spacing below 1, the unit roundoff, the both-sides check in the code
(`deriv.hpp` line 1028). The default-step numbers (0.0061 K; 1e-14 K below the spacing). The narrow-feature
analysis: h^2 |f'''| / (6 |f'|) at the default step reproduces the measured 2.8e-5 (2.77e-5). The small-|x| case
(6.1e-6 > 1e-6; non-finite). The noise model (eta^(1/3) s, eta^(2/3), (eta/eps)^(1/3) = 356; 4.6e-6 for eta =
1e-8) and the coverage figures 100.0, 89.3, 90.7 %. The float figures (eps^(2/3) = 2.4e-5, step 0.0049). The
50-digit eps = 2^-167 = 5.3e-51.

Chapter 4. Eq. (4.2); the Central<4> derivation (M_1 = 2a + 4b = 1, M_3 = 2a + 16b = 0, a = 2/3, b = -1/12,
M_5 = -4, C = 1/30; the magnitude of (4.3) is right, only its sign is not). The Lagrange example w_(-1) = -2/3.
Exactness for degree P + N - 1 and the test's Central<4> on x^4. The order ladder: 2^P per halving, S from 1.00 to
2.28, the model minima 3.7e-11 and 5.9e-15, the measured medians 5.8e-12 and 1.6e-15, the diminishing returns (1.9
then 1.6 decades). The third term (S/2) eps |x f'| / h and its absence for the two-node rules. The test budget
(C + 2S) eps^(P/(P+1)) and "4 to 16 times below", quoted correctly from `test_stencils.cpp`. The one-sided moment
count (P + 1 conditions, order P), the even moments by antisymmetry, S and C of Forward<1> to <6> and Central<2> to
<10> (Tables 9.1 and 9.5 agree with exact recomputation). The thread's steps 0.74 K and 5.8 K and the span of three
steps. E6's numbers and the claim that Central<2> fails at 500 K with the NaN wrapper.

Chapter 5. Eq. (5.1) with its sign (rule - f'' = + h^2 f^(4)/12), the f''' and f'''' rules against Tables 9.3 and
9.4, the constants 1/6, 1/12, 1/4, 1/6, Eq. (5.2), the exponents 2/3, 1/2, 2/5, 1/3 and their values, the
relative-factor caveat and 4 eps^(1/2) x 3.7 = 2.2e-7 (recomputed 2.23e-7), the coarse criterion (radix^(digits-1)
< 10^7 is eps > 1e-7: float, decimal32 and bfloat16 yes, double no, as `detail::coarse` computes it), the margin
note on sin x / x, Ridders' per-order facts (s/4; 2 or 4 evaluations per level; f(x) once, kept across restarts;
the median 29 for f''''), Sidebar 5.2 (exact offsets for h <= |x| by Sterbenz, exactness on quadratics, h^2 never
formed, the 1.5e-154 bound, the order of the divisions in `divide_by_powers`).

Chapter 6. Eqs. (6.1) to (6.3) against `detail::ridders` (lines 1244-1252 and 1270-1286). Table 6.1 recomputed
from `ridders-tableau.csv`: the floor at level 0 (4.66e-15 relative) and its doubling, column 0 gaining a factor 4
per level, column 1 starting at 3.9e-9, column 2 at 2.4e-13 at level 2, the stop after level 6 by the two-level
stall (best e_(4,4) = 8.15e-14 relative; levels 5 and 6 pinned at their floors 1.49e-13 and 2.99e-13), abs_error =
e_(4,4) = 2.66e-12 (above |A_(4,4) - A_(6,6)| = 2.0e-12), the actual error 1.95e-12, the base difference's minimum
5.8e-13 at level 13, column 2 at the floor by level 3. Every statistic of the random study for all four orders,
recomputed exactly from the CSVs (coverage 99.67, 99.85, 99.96, 100.00 %; uncovered 9, 4, 1, 0; median ratios
0.141, 0.059, 0.015, 0.0041; worst 1.81, 1.36, 1.001, 0.33; 13, 8, 3, 12 zero-error cases; 0, 11, 61, 61 restarts;
evaluations 8/14/18, 7/17/21, 12/32/44, 9/29/41); the parity of the counts (even for f', odd for f''); the 11
restarts of f'' being log x and sqrt(x) below 0.25, and 1/x not restarting. The edge figure: two evaluations per
halving, the cap first reached at k = 27 (23 failing levels plus the nine the tableau needs), ten digits to k = 30,
a value to k = 34, and the derivation of 2^-34 (k - 4 failing levels, at least two left of 32) and of 2^-31 for the
initial step s/4. The stencil's 1.2e-7 at x = 0.01 (h^2 / (3 x^2) = 1.22e-7). The untrusted test as "two correct
bits". Sidebar 6.1 (4^7 2^115 > 2^128; the difference form exact for equal parents; two rows filled from T(0)).
Sidebar 6.2's rules against the code line by line (the pre-checks, the estimate, the best entry, the stall
counter's i > 0 guard, the exact stop, the final check, the untrusted test with 4 and 64, the fallback, the
resolution exception, f(x) once). Sidebar 6.3 (sin(+-2 pi) and sin(+-pi) at the first two levels, +-1 at the third
with D = 128, the restart from 1/64; cos(64 pi x) for f'' stopping at level 1 on exact agreement). The kink
numbers: J/(2h)/cp = 4.3 and (7/12) J/(h cp) = 0.041 recomputed; the jumps 0.768 J/mol, -7.5e-5 J/(mol K) and
-17 % recomputed exactly from the published coefficients.

Chapter 7. Newton: quadratic convergence (e_(i+1)/e_i^2 = 8.8e-5 = |H''/(2 H')| at every step), the final residual
at eps |H| = 4.8e-12, three evaluations per iteration, the delta term of an inexact slope. The sensitivity closed
forms and -E_a/(RT) = -16.04. The Colebrook steps 0.61 and 46.4, the rounding-branch estimate eta |lambda| /
(h |lambda'|) = 9e-5 against the measured 5.5e-5, eta^(2/3) = 2.2e-7. Section 7.5's checks (J/(2h) doubling as
the step halves; 11.3 and 6.1 digits recomputed from the values). The multiprecision figures (eps^(1/7) = 6.6e-8,
eps^(6/7) against 7.2e-45, eps^(2/3) against 1.5e-35). The constexpr example (12.0625 = 12 + h^2 f'''/6; Ridders'
column 1 removes it exactly). The pitfalls' numbers against `values.tex`.

Chapters 8 to 10. Default steps, nodes and evaluation counts per family and order, the caps 64, 65, 128, 129, the
failure list, the precondition wording, the coarse criterion, root_epsilon within an ulp with N to 2000 against
80-digit references (`test_math.cpp` lines 44 and 59). The weight tables: Central<4> for f', Central<2> for f'' to
f'''', Forward<2> and Backward<2> for f'''' recomputed exactly; every row satisfies the moment conditions by the
generator's check. Table 9.9's exponents, steps and model errors (spot-checked); Table 9.10; Tables 10.2 and 10.3
against the CSVs; the headline digits (11, 4, 13, 9 by floor); the guarantees with their hypotheses (the FMA and
extended-precision wording of decision A5); "the estimate is not a bound" wherever abs_error is explained.

Appendices. A: eps = 2^(1-p), the ulp statements, the unit roundoff, the subnormal floor (`value_rounding`), the
overflow policy. B.1: Eq. (B.1), (B.2), the symmetry argument (M_j vanishes for j - N odd), the rounding term; B.2:
the Lagrange formula, uniqueness, m = P + N and P + N - 1, the worked example and Chapra's row; B.3: the induction,
the two forms of (B.4), the a posteriori estimate and why the floor is needed; B.4: (B.5) to (B.7) and every figure
quoted.

Notation and sidebars. `macros.tex`, the notation tables and the glossary agree with the chapters' usage (N, P,
S, C, s, w, eta, D(h), A_(i,j), r_i, e_(i,j), e-hat, Omega); the one slip is Sidebar 6.2's sum |w_k| / (d h^N)
(R1-numerics-10). Every sidebar (2.1, 3.1 to 3.3, 4.1, 4.2, 5.1, 5.2, 6.1 to 6.3, 7.1, 7.2) is self-contained, and
the running text that refers to one does not depend on it.

## Verdict

The mathematics of the book is sound and its evidence is unusually well tied to the library: every derivation I
recomputed (the Taylor remainders, the moment conditions and the Lagrange weights, the Richardson tableau and its
estimate, the error model with its best step, constant and cost, the Ridders rules against the code) came out as
stated, with one sign error in a numbered equation (Eq. 4.3) and one proof that does not cover the case it is
applied to (B.1 for the central rules above order 2). The quantitative claims are generated and, where I
recomputed them, exact; the two places where the text misreads its own data are the attribution of the N = 2
histogram's tail to the eleven restarts and the "exact stops, all on x^5". The remaining findings are precision of
language (the sensitivity factors, the "factor of ten" of Table 10.1, the caps' rationales, the "factor 4") and
three caveats a numerical analyst would expect (how the method differs from Ridders' paper, which functions make the
estimate tight, and cancellation inside f). Nothing found changes the book's advice to the reader. With findings 01
to 03 corrected and the minor ones addressed, the book withstands scrutiny on its numerics.

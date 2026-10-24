# R2-cline

Independent review by Cline 3.0.73 in plan mode (read-only) with the model DeepSeek V4.1 Flash through ClinePass (`cline -p -P cline-pass -m cline-pass/deepseek-v4.1-flash --thinking high`), round 2, on commit 5432abf. The run hit its 1650-second limit while printing; the main session extracted the review text from the run log (the streamed thinking lines removed) and saved it here verbatim. It contains the findings in several passes (Appendices A and B, Chapters 9 and 10, the back matter and refs.bib; then Chapters 7 and 8 against the headers) and a few of the model's working notes between them; the editor reproduces every finding before acting on it. Finding ids: R2-cline-NN in the order of appearance, assigned by the editor in the decision letter.

I have completed the verification. Here are my findings (all arithmetic done with exact-fraction Python against the files).

## Findings

**A. `docs/deriv-book/chapters/a2-derivations.tex`**

1. **a2:152-153 — MAJOR.** The factor in `E(h_best)` has N and P swapped, contradicting the next paragraph.
   - Verbatim: `152: At the minimum the two terms stand in the ratio $\dorder : \aorder$ (truncation to rounding), so` / `153: $E(\hbest) = (1 + \aorder / \dorder)\, \ampl \epsmach \abs{f} / {\hbest}^{\dorder}$, which is of order`
   - Counter-evidence: with T:R = N:P, T = (N/P)R, so E = R(1 + N/P) = R(P+N)/P. a2:166-167 says the opposite: `since those values are the fractions $\dorder / (\aorder + \dorder)$ and $\aorder / (\aorder + \dorder)$ of` / `$E(\hbest)$` (T/E = N/(P+N), R/E = P/(P+N)). Concrete: Central<2>, f' (N=1, P=2, S=1, C=1/6): h*=3^(1/3)ε^(1/3), R = 0.6934 ε^(2/3), T = 0.3467 ε^(2/3), E = 1.0401 ε^(2/3) = 1.5·R, not 3·R. The generator agrees with the N:P reading: `tools/generate_data.cpp:1510-1511` "with the two terms in the ratio N : P at the balance" and `:1499` `(N * std::pow(k, -P) + P * std::pow(k, N)) / (P + N)`.
   - Fix: write `(1 + \dorder / \aorder)`.

2. **a2:178-179 — MAJOR.** False superlative: 15.1 is not the largest cost.
   - Verbatim: `178: \num{\valStepCostCentralTwoFourth}; for \code{Central<4>} and $f''''$ it is \num{\valStepCostCentralFourFourth}, the` / `179: largest of the library's rules. These agree with the design record's figures of 3.5, 11 and 15`
   - Counter-evidence (same formula, `generate_data.cpp:1495-1502`): Central<10>, N=4: S = 623872/14175 = 44.0125, C = 59/277200 → κ = 2.2452, cost = (4κ^-10 + 10κ^4)/14 = **18.15**; Central<8>, N=4 = 17.72; Central<6>, N=4 = 16.88 — all > 15.1. The generator only evaluates the cost for the five rules the book quotes (`:1503-1507`), so the claim was never computed.
   - Fix: delete "the largest of the library's rules" (or say "the largest of the four rules the design record quotes").

3. **a2:180-181 — MINOR.** "1.64 … is a rounding of κ = 1.8" is not a true statement.
   - Verbatim: `180: (\code{docs/PHASE4_API.md}); the 1.64 that the library's documentation gives for \code{Central<10>} is a rounding` / `181: of $\kappa = \num{\valStepConstantCentralTenFirst}$, where the formula gives \num{\valStepCostCentralTenFirst
- `abs(mp(*r) - e)` ✓
- `r->value - e` ✓
- `*r == e` ✓ (Boost's comparisons DO accept Estimate, as chapter 8 line 149 claims)
- `*r < e` ✓
- `mp(*r) - e` ✓
- `derivative(gl, mpet(1))` with a lambda returning an expression (no `-> T`) compiles ✓
- `derivative(g2, mpet(1))` with `-> mpet` ✓

So the claim "it
}.
`
   - Counter-evidence: 1.64/1.63 is the *cost* at κ = 1.7975 (my recomputation: 1.634), not a rounding of κ; the documentation it refers to (`include/numerixx/deriv.hpp:521`, `docs/design.md:112-113`)
 says "For f' it lies between 1.4 and 2.0 and costs at most a factor 1.64" — i.e. 1.64 is the cost, and the documentation never names Central<10>.
   - Fix
: "…is the cost at κ = 1.8 rounded to three digits, where the formula gives 1.63".

**B. `docs/deriv-book/chapters/a1-floating-point.tex`**

4. **a1:73 — MINOR.** "up to ulp(x)/2" is false when x
 +
 h crosses a power of two.
   - Verbatim: `73: $\
fl(x + h)$, which differs from $x + h
$ by up to $\ulp(x) / 2$. Two consequences follow, both met in`
   - Counter-evidence: x = 1 − 2⁻⁵³ (ulp(x) =
 
2⁻⁵³), h = 2⁻⁵²: x + h = 1 + 2⁻⁵³ is exactly the midpoint of 1
 and
 1 + 2⁻⁵², ties-to-even gives fl
(x+h) = 1,
 error = 2⁻⁵³ = ulp(x) = 2 × ulp(x)/2. The library

H
's own comment notes the two-sided spacing at a power of
 the radix (`deriv.hpp:1025`).
   -
 Fix: "differs from
 x + h by up to half the spacing at the result (up to ulp(x) when x + h crosses a power of two)".

**C
.
 `docs/deriv-book/chapters/10-guarantees.tex`**

5. **10

Hmm, but there's one place: `detail::weighted_sum` caches values in a `WeightedSum<T>` struct with T members — it stores `static_cast<T>(...)` presumably ✓ (T, not auto).

So the book's claim "give a
:213-216 — MINOR.** "never from nxx itself" is falsified by the tests.
   - Verbatim: `216: prototype, never from \n

Let me actually RUN a test: compile and run a small
xx{} itself, and each tolerance is justified by an error estimate in a comment. The stencils,`
   - Counter-evidence: `tests/test_ridders.cpp:313` `// References from the first C++ run (relative estimates): bin50 exp 4.4e-
47, log 3.7e-46; dec50 exp 4.1e-46,`; `tests/test_orders.cpp:633` `// the estimate bounded per type by twice the first C++ run, as for f' (docs/

Alternative: I can check at compile time with constexpr! Boost's cpp_bin_float is constexpr-usable? Hmm, `exp` for multiprecision is constexpr in recent Boost? Probably yes for c
design
.md: the decimal prototype does` and `:634` `// not model these types' guard digits). The first run gave relative estimates of 2.7e-43, 1.6e-39 and 8.1e-36 for`.
   - Fix: "…or from an independent prototype — for the multiprecision estimates the prototype does not model, from the library's own first run".

**D. `docs/deriv-book/chapters/09-tables.tex
`**

6. **09

```
constexpr mpet one(1);
constexpr auto d = nxx::deriv::der
:78-79 — EDITORIAL.** Caption overstates: "every rule is exact in float".
   -
 Verbatim: `78: 6,222,216, is the centre weight of \code{Central<10>}; every numerator fits in 24 bits, so every rule is exact in` / `79: float.\label{tab:weights-central-fourth}}`
   - Counter-evidence: what is exact is the weights; the float column of `tab:accuracy-matrix` is 5.0e-6 (Central<2>, f') to 0.014 (f''''). Compare the careful wording at `10:191-192
static_assert(d->value > mpet("2.7"));
```
If the lambda returns a dangling expression, the constexpr evaluation would fail. Hmm, in a constant expression, dangling references are...
`:
 "The integer weights are exact in every type with 24 significant bits or more".
   - Fix: "…so every weight is exact in float".

**E. `docs/deriv-book/chapters/99-back
matter.tex`**

7. **99:29 — EDITORIAL.** "the rule spans K h on each side of x" is false for the one-sided rules.
   - Verbatim: `29: $K$ & the largest $\abs{n_k}$ of a rule; the rule spans $K h$ on each side of $x$\\`
   - Counter-evidence: `Forward<P>` has nodes 0…P+N−1 (`deriv.hpp

Hmm, this is getting deep
:474-478`), so it spans K h above
 x
 and 0 below; the header defines K as
 a radius:
 `

Actually
der
iv.hpp:538
` "With K = max(|lo|, |hi
|
) the largest
 node
".
   - Fix: "the rule reaches at
 most
 K h from x".


**F. `docs
/
deriv-book/data/README.md`**

8.
 **
README
:13-14 — EDITOR
IAL.** The stated formatting
 convention om
its the documented exception.
  
 - Verbatim: `
13: Errors
,

And
 estimates
 and steps have two significant digits, counts are integers
, coverage is a percentage with one`
 /
 `14: decimal; typeset them with \num
{\valX} or \SI{\valX
}{unit}.`
  
 - Counter
-evidence: `data
/values.tex:152
 \

If
newcommand
{\val
RiddersCoverageThird
}{
99.96
}`
 (two decimals),
 which
 the
 generator documents: `tools/generate_data.cpp:1040-104

Hmm, let me check the actual
1
` "Coverage in
 percent
 to one decimal, as
 the style decision
 says, except where one
 decimal would read 
100.0
 / while a case is uncovered: then
 two".
   - Fix: add the exception

Let me grep the header for `
 to the README sentence


.

**G
. `docs/deriv-book/refs.bib`**

9. **refs.b
ib:200 (with `:212`) — EDITORIAL.** The entry comment's instruction is not followed at any use site.
   - Verbatim: `200: % inexact slope on Newton's method in chapter 5; cite the section at the sentence.` and `212: % Metadata verified; section 8.1 (finite-difference derivative approximations, from p. 194) from the
 contents.`
   - Counter-evidence: `chapters/03-rounding.tex:162` `standard in numerical optimisation.\cite{dennis1996}` and `chapters/07-using.tex:156` `$h
$
 adds a term of about $\delta$ times the error to each step,\cite{dennis199
6
}` carry no section; `chapters/03-rounding.tex:133` `digits.\cite{n
oc
ed
al2006}` likewise. R1-decision row 25 required "the section numbers move from the notes to the citing sentences".
   - Fix: `\cite[Sect.~5.4]{dennis1996
}` at

Hmm, actually, let me test
 the two sites and `\cite
[Sect.~8.1]{noced
al
2006}` at 3.3.

10. **refs.bib:255-264 — EDIT
OR
IAL.** `nist-webbook` has no `year` field, so the entry prints with no date (the retrieval
 date is only inside `note
`).
   -
 Fix: add
 `year = {2026
}` (the retrieval
 year
) or cite
 the WebBook
 as a
 database
 release
.

*
Cross
-scope note (verified, outside
 my chapters
):* `
deriv.hpp:385
`
 "
The largest

Alternative: use a constexpr evaluation test. `
 intermediate is 12! times 4! (`Central<10>`, n =
 
4)" disagrees with `ch
apters/04-stenc
ils
.tex:103-104` "The
 largest intermediate it reaches is $12!$ …
 $
12
! \
cdot 4!$ is a bound that is
 not reached
". My recom
putation
 gives
 max
 |
n!·
poly
[n]|
 = 7,
111,104 and max
 |q| = 12! = 479,001
,600,
 so 12!·4! is indeed never reached and the header is the loose one.

## Areas checked and found correct (no findings)

- **Weight tables T4–T11 (all 8 bodies):** nodes, numerators, common denominator, S

Let me try:
```
constexpr auto one = mpet(1);
constexpr auto dd = nxx::deriv::derivative
 and

 C = |
M
```
Hmm, if it
_(N+P)|/((
N+P)!·d) and the sign of the leading
 term
 reproduce
 my exact
-f
raction Lagrange
 computation
 in every
 one

Actually, let me
 of the 44 rows (
e
.g
. Central<2>/f': {
−1,1}/2, S=1.
00, C=1/
6
, +; Central<4>/f': {1,−8,8,−1}/12, 1.50, 1/30, −; Central<6>/f'' {2,−27,270,−490,270,−27,2}/180, 1/560, +; Central<10>/f'''' 13 nodes, d=453600, S=44.0, 59/277200
, +; Forward<1..

Let me grep.
6>/f' d = 1,2,6,12,60,60, S = 2.00…27.7, C = 1/2…1/7, alternating signs). **Largest numerator is 6,222,216**, the centre weight of Central<10>/f'''' (all values < 2²⁴ = 16,777,216; max denominator 453,600), so `09:78`, `10:191-193` and the 24-bit claim hold.
- **`table-default-steps.tex`:** all 44 rows match eps^(1/(P+N)), eps^(P/(P+N)) (double) and the counts P+N−1 (central) / P+N (forward).
- **`table-cost.tex`:** stencil counts = P+N−1 / P+N exactly; Ridders caps 64/65/128/129 = 32 levels × 2/2/4/4 (+1 for even N, `ridders_divisors{32,4,4,4}`, `ridders_levels=32`); medians/ranges recomputed from the CSVs (8/14/18, 7/17/21, 12/32/44, 9/29/41).
- **`table-accuracy-matrix.tex`:** all 72 cells equal the medians in `data/accuracy-matrix.csv`; headline digits
 11/5/14/9 and the model ratios
 6.3, 5.6, 5.8, 32, 37 all reproduce.
- **`table-ridders-constants.tex`:** ratio 2, divisors 32/4, stall 2, factors 4 and 64, cap 32, evaluations 2,2,4,4 and f(x) once 0,1,0,1 — all match `deriv.hpp` (`factor*4`, `stall == 2`, `/
4`, `/64`, `ridders_levels`, `ridders_levels`).
- **`tab:ridders-statistics` (ch10) and the edge limits:** recomputed from `ridders-coverage-N
*.csv`: coverage 99.67
/99.85/99
.96/100.00, uncovered 9/4/1/0, median ratio 0.1406/0.0592/0.01513/0.004081, worst 1.807/1.359/1.001/0.3291, evals and restarts 0/11/61/61. `ridders-edge.csv` confirms the return limit at k=34 (rel_estimate 0.10), failure at k=35, ten digits to k=30 (9.2e-14 vs 1.16e-10 at k=31) and 64 evaluations from k=27; PHASE3_API.md:207-208 / PHASE4_API.md:300 give 2⁻³⁴ and 2⁻³¹ as the book states.
- **Appendix A:** eps = 2^(1−p) (2⁻⁵²), digits 24/53, ulp(1) = 2.2e-16, below 1 = 1.1e-16, ulp(1000) = 2⁻⁴³ = 1.1e-13, min/denorm_min/max, cpp_bin_float_50 168 bits / 5.3e-51 = 2⁻¹⁶⁷, subnormal spacing = denorm_min, `value_rounding` = eps|f| but never below denorm_min (`deriv.hpp:695-702`), the pre-arithmetic overflow checks with the final-quotient exception, and Sterbenz as stated (one sign, Thm 4.3.1) — all verified.
- **Appendix B (other than finding 1–3):** the moment conditions and C = |M_(N+P)|/((N+P)!d); the single-point form (Central<4> product 4, m = 5 → +h⁴f⁽⁵⁾/30, consistent with the table's −1/30); products 1,0,1 for the three-point f'' and −8/3,2/3,2/3,−8/3 for Central<4>; w_k = N![t^N]L_k and its moment argument; the symmetry argument (M_j = 0 for j−N odd, one node fewer); the Richardson induction (column j starts at h^(2j+2)) and the a posteriori estimate; h_best, κ = (N S/(P C))^(1/(P+N)), κ ∈ [1.4422, 2] for f' (max at Forward<1>, min at Central<2>), the costs 1.12/1.63/3.54/11.2/15.1 and eq. (derivations-cost) — all correct.
- **Chapter 10 guarantees:** no `throw`/`try` and no static/global state in the headers; no allocation; fn once per node
 and
 never
 across
 x
 for one
-sided rules
; Ridders ≤ 64/65/128/129 (test counts 64/64); failures as
 values with the single final-quotient exception; exact weights (`static
_
assert` + `holds
_integers`); exact
ness on polynomials with the 2
⁵³ bit budget;
 compile
-time = run-time (
no `-ffp-contract`/`-
ff
ast-math` anywhere in
 the project's CMake; five
 CI toolchains GCC
/Clang/cl
/clang-cl/E
mscripten with `-Werror`/`/
WX`); `root_epsilon` within an
 ulp
 (
80-digit refs, N to
 200
0;
 
60-digit for multiprecision); types tested (float, double, bin50 with and without
 ET
,
 dec
50; long double in
 `
root
_
epsilon` and
 the `no
isy` check
 only). The deferred list matches `docs
/plan.md:49
-
50`, PHASE
2:226
, PHASE3:265
-272, PHASE4
:418-422 and DER
IV
_TECHNIQUES
.md:162
-
165. The st
encil tolerance
 claim
 `(C + 2
S) eps^(P/(
P

H
+N
))` and "
4 to 16
 times below it
" are verb
atim in
 `tests/test_stencils.cpp:202-204`; the "six times its estimate / twice the
 prototype" claim matches `tests/test_ridders.cpp:3-8, 87
-88`; "
one doctest executable and three build tests
" matches `tests/CMakeLists.txt
:29

H
-54
`.

- **Quick reference:** every statement checked
 against
 the header and the macros — call forms and
 defaults,
 `
Est
imate` result
 type,
 `
d
.value()` throwing, `if
 (d)`, the failure list and the `
h
 ≤
 0`
 assert, Central P
 =
 2…10 with
 P
+N−1 evaluations, Forward
/
Backward P = 1…
6 with P+N, R
idd
ers' median 
14 / caps 64-

OK
65
-128-129 / 14 digits, the edge distances 2⁻³
⁴
 and 2⁻³⁰, the two-step check
, `root_epsilon
<3,double>()` and `std::cbrt(eta)`, the pitfalls, and every cross-reference
 label (`ch
:api`, `ch:tables`, `tab:accuracy-matrix`, `tab:cost`, `sec:guar
antees-guarantees`, `tab:method-summary`, `tab:
method-choice`) resolves.

- **refs.bib:** no entry
 carries
 `note
 = {UNVERIFIED}`
; all 25 keys are
 plausible
 and the risky
 ones agree
 with the R1 evidence
 table
's
 Crossref/catalogue records (ridders1982 Adv. Eng. Softw. 4(2):75-76 1982; richardson191
1 210:307-357; richardson


1927 226:299-361; ly
ness
196
7 
4(2
):
202-210
; squ
ire
1998 40(
1):110-112
; fornberg198
8 51(184):699-706
;
 chapra
2015 7th ed
.
 ISBN 978-0
-
07-339792-4; more201
2 TOMS 38
(3) art. 19
; gill
1983 4
(
2):310-321; dennis1996

Hmm, hold
 SI
AM Classics 16
; noced
al
2006; sterbenz
197
4 with
 Thm 4.3.
1
/p
. 
138 and ISBN 0-13-322495-
3
; shomate
1954 58(4):368-372; chase199
8
 Part II
 p
. 162
1; nist-we
bbook DOI
 
10.184
34/T4D303). Every
 cited key exists
 in the file; the two
 uncited entries say
 so
. Use sites
 in my chapters
 support their claims (Fornberg/Chapra for the published rows at `
09:14-
15`, `10:
214
-215`, `a
2:102-104
`
; Burden §
1
.1/§4
.
1/§4.2;
 Ster
benz Thm 4.3.1 at `
a1:83
`; Higham/Goldberg/IEEE 754 at `a1:6-7
`; more2012 for the
 deferred noise-aware step at `10:258`, matching DERIV_TECHNIQUES.md §2.4 and the R1 note "Appropriate at 
3.
6 l. 294, 7.9 and 10.7").
- **Back matter:** all notation and glossary entries check out against the header (moments, C, S, eps/
ulp

OK. I'll report it as a minor finding with careful wording.

H
/fmax/scale/noise/E(h)/h_best/D_i/A_ij/r_i/estimate, aggregate, aliasing, cancellation, coverage, kink, level/column/tableau, lvalue, non-finite, restart, Richardson, rounding floor, structured binding, unit roundoff, untrusted result) except item 7 above. There is no index in the book and `PLAN.md:66-67,95` does not call for one, so nothing to check there.

Fine. Let me now also check the claim in chapter 7 line 284-285: "In arithmetic, write \code{mp50(*r)} or \code{r->value}, because Boost's operators do not see the conversion of an \code{Estimate}; this holds for every Boost type." ✓ verified (the error message shows no match for operator- with Estimate operand) ✓.

And chapter 8 line 149: "its comparisons do accept it" ✓ verified (*r == e compiles).

Good. Now let me wrap up with a few remaining checks:

1. Chapter 7 §7.3 (Choosing a method) table row: "$f$ varies on a scale $w$ much smaller than $|x|$, or $|x| \ll 1$ in your units & \code{Central<2>}, or \code{Ridders} from $h = w/4$ & $\epsmach^{1/3} w$ & 2 & ...". The "Evaluations" column says 2 for a row offering Ridders too. Minor.

2. Let me check the numbers in tab:method-choice against the referenced macros once more — done.

3. Let me verify the claim in chapter 7 line 373: "is far above its usual size (\num{\valKinkOutsideRiddersRelEstimate} of the value away from the break)" = 2.7e-13 ✓.

4. Chapter 7 line 396-402 (edge pitfall) ✓.

5. Chapter 7 line 386-389 (aliasing) ✓.

Now let me re-examine something I noticed: chapter 7 line 47: "The one reason so far is \code{error::non_finite}." ✓

And line 51: "The central rule at $x = 0$ samples $\log(-h)$, which is NaN, so the first call fails." ✓ (h = 6.055e-6 at x=0 since scale(0) = 1).

Hmm — one thing: at x = 0, the default step for Central<2> f' = root_epsilon<3,double>() × max(0,1) = 6.055e-6 ✓. The snippet says "samples log(-h): NaN" ✓.

Now let me check chapter 8's §8.2 statement: "One call of each form, on the thread's enthalpy at T0 (code/08-api-usage.cpp)" and the snippet `<api-overloads>`: `auto const b = derivative(H, T0, 0.01);` — hmm, `

Hmm, in the snippet: `auto const d = derivative<Forward<2>>(H, T0, 0.1);` — Forward<2> at h=0.1 for f' of H at 1000 K: the error is h²|H'''|/(3 cp) ≈ 6e-10 per the code comment ✓ and the check tolerance is 1e-7 ✓.

Hmm, `Forward<2>` samples x, x+h, x+2h — all ≥ T0 = 1000 K, inside the fit's range [500, 2000] ✓ fine.

Now, let me double check the claim in ch8 line 227-230 once more: "On $\mathrm{e}^x$ at 1, \code{Central<2>} is \num{\valFloatFourthCentralAtOneRelError} off and Ridders' method \num{\valFloatFourthRiddersAtOneRelError}, with a relative estimate of \num{\valFloatFourthRiddersAtOneRelEstimate} that says the result is rough (measured; \cref{sec:orders-float})." = 0.0067, 9.7e-5, 0.28. ✓ values exist.

Hmm, "The result is the same with or without the warning." ✓ (the deprecated overload just forwards).

Now, one more important check: chapter 8's §8.2 statement (line 86-88): "Preconditions. h > 0 and h finite, checked by an assert: a bad step is a bug in the caller. Under NDEBUG the check is compiled out and the behaviour unspecified (Forward<1> with a negative step evaluates below x); validate a computed step, or use the overload without a step, which has no precondition."

Hmm — "the behaviour unspecified". Actually with NDEBUG and h<0, the behaviour is... well, for a stencil with h<0, the nodes are still x + k h; the "Forward" direction is inverted. For Ridders with h<0? h_i = h/2^i → negative; `h_i > 0` fails → the level is treated as the resolution → the tableau ends → probably returns non_finite. So "behaviour unspecified" is a fair description. OK.

Now let me also double check one of the trickiest claims: chapter 8 line 107-108: "or, for a rule with more nodes, $x + h$ or $x - h$ rounds to $x$ or neighbouring points round to the same value" ✓ matches the code (line 1028: `if (T(x + h) == x || T(x - h) == x) return non_finite;` and line 1226-1229 for Ridders).

Hmm, but for a rule with more nodes, the code only checks `x + h` and `x - h` (not x + 2h vs x + h etc.) in the stencil path (line 1028) — the comment at line 1024-1027 says both sides are checked. For Ridders (line 1226-1229) it also checks 2h for N ≥ 3. The book says "x + h or x - h rounds to x or neighbouring points round to the same value" — hmm, "neighbouring points round to the same value" is vague but consistent-ish. For a Central<4> stencil (4+ nodes) at h near 0.6 ulp, x+h and x+2h could round to the same value while x+h != x and x-h != x. Then the code would NOT detect it (it only checks x±h against x) — but the book says the failure is detected ("error::non_finite when ... h is below the resolution of x: ... or, for a rule with more nodes, x + h or x - h rounds to x or neighbouring points round to the same value").

Hmm — is that a claim that the library returns non_finite in that case? The sentence structure: "error::non_finite when no finite estimate exists: ... or h is below the resolution of x: the points of a two-node rule coincide, or, for a rule with more nodes, x + h or x - h rounds to x or neighbouring points round to the same value."

Hmm — the header's own doc for the enum says only "or h is below the resolution of x". The book expands. Is "neighbouring points round to the same value" a case that yields non_finite? For Central<4> with h = 0.6 ulp(x): x+h → x + ulp; x+2h → x + ulp (rounds to the same). Then the stencil computes a value from a set of points where two coincide — is that detected? The check `T(x+h) == x || T(x-h) == x` fails to catch it. Then the weighted sum uses f(x+h) and f(x+2h) which are equal → the result is a finite but wrong value. Hmm! So the book's claim would be wrong in that case.

But hold on — test_stencils.cpp:153-154 has:
```
static_assert(derivative<Backward<2>>(power<1>, -1.0, 0x1.8p-54) == non_finite);
static_assert(derivative<Central<4>>(power<1>, -1.0, 0x1.8p-54) == non_finite);
```
0x1.8p-54 = 1.5 * 2^-54 = 0.75 ulp at 1. For Central<4> at x = -1 with h = 0.75 ulp: x + h = -1 + 1.5*2^-54; the spacing above -1 is 2^-53 = 2*2^-54, so 1.5*2^-54 = 0.75 * 2^-53 → rounds to -1 + 2^-53 (hmm, 0.75 ulp rounds to 1 ulp). And x - h = -1 - 1.5*2^-54 → below -1 the spacing is 2^-54, so 1.5*2^-54 = 1.5 ulp → rounds to -1 - 2*2^-54 = -1 - 2^-53. Hmm. So x+h = -1 + 2^-53 ≠ x, x-h = -1 - 2^-53 ≠ x. And x + 2h = -1 + 3*2^-54 → 1.5*2^-53 → rounds to -1 + 2*2^-53 = -1 + 2^-52 (since 1.5 ulp rounds to 2 ulp). So x+2h ≠ x+h ✓. Hmm, so why is it non_finite?

Hmm, so the test says Central<4> at x=-1, h=0.75ulp is non_finite. Per the code, the only checks for size>2 are `T(x+h)==x || T(x-h)==x`. Since neither holds, we go on... Hmm, then how does it return non_finite? 

Hmm, maybe via `detail::fits(x, h, q)`: q = max(2*k_max, denominator) = max(2*2, 12) = 12 for Central<4> f' (nodes -2..2, denominator 12). fits requires h <= max/(2q) and |x| <= max - q h. With h = 1.5*2^-54 ≈ 8.3e-17, that's tiny, so fits is true. Hmm.

So how does it fail? Let me re-read: maybe the check is in the else-branch of the size==2 case:
```
if (T(x + h) == x || T(x - h) == x) return std::unexpected(error::non_finite);
```
Hmm, for x = -1 (a power of two), spacing above -1 is 2^-52? No wait. -1 = -2^0. The spacing (ulp) of numbers in [1,2) is 2^-52, and in [0.5,1) is 2^-53. For x = -1: ulp(-1) = 2^-52 for values just above -1 in magnitude... hmm, careful with signs.

x = -1. Numbers just below -1 (i.e., -1.0000001) have magnitude in (1,2] → spacing 2^-52. Numbers just above -1 (i.e., -0.99999) have magnitude in [0.5,1) → spacing 2^-53.

So h = 0x1.8p-54 = 1.5 * 2^-54 = 0.75 * 2^-53.
x + h = -1 + 0.75*2^-53 → the representable values near -1 from above are -1 + k*2^-53 → 0.75 rounds to 1 → x + h = -1 + 2^-53.
x - h = -1 - 0.75*2^-53 → representable values below -1 are -1 - k*2^-52 → 0.75*2^-53 = 0.375*2^-52 → rounds to 0 → x - h = -1 = x! 

So x - h == x ✓ → non_finite ✓. That explains the test.

OK so my hypothetical (both x+h and x-h differ from x but coincide with each other) requires h > ulp/2 on both sides... but the spacings differ by a factor 2 at a power of two, so if h > ulp_above/2 then h could still be < ulp_below/2? Hmm: at x=-1, ulp above (toward 0) = 2^-53, ulp below = 2^-52. If h = 0.3*2^-52 = 0.6*2^-53: x + h = -1 + 0.6*2^-53 → rounds to -1 + 2^-53 (0.6 → 1) ≠ x ✓. x - h = -1 - 0.3*2^-52 → rounds to -1 (0.3 → 0) == x → detected ✓.

For both to differ from x, we need h ≥ 0.5*ulp_above and h ≥ 0.5*ulp_below; since ulp_below = 2*ulp_above, h ≥ ulp_above suffices for both. Then x+h and x+2h: h ≥ ulp_above = 2^-53 means 2h ≥ 2^-52 = ulp_below; x+2h rounds to at least -1+2^-52, while x+h rounds to at least -1+2^-53... they could coincide if x+h rounds to -1+2^-52 (when h is just under 1.5*2^-53?) hmm: h = 1.4*2^-53: x+h = -1 + 1.4*2^-53 → rounds to -1 + 2^-53 (1.4 → 1)? Hmm, 1.4 ulp_above rounds to 1 ulp_above = 2^-53. x+2h = -1 + 2.8*2^-53 = -1 + 1.4*2^-52 → rounds to -1 + 2^-52 (1.4 → 1). So x+h = -1+2^-53 ≠ x+2h = -1+2^-52 ✓. Hmm.

Try h = 1.6*2^-53: x+h = -1 + 1.6*2^-53 → rounds to -1 + 2*2^-53 = -1+2^-52. x+2h = -1+3.2*2^-53 = -1+1.6*2^-52 → rounds to -1 + 2*2^-52 = -1+2^-51. Hmm, so x+h = -1+2^-52 and x+2h = -1+2^-51, differ ✓.

Hmm, x-h: -1 - 1.6*2^-53 = -1 - 0.8*2^-52 → rounds to -1 - 2^-52 ✓ ≠ x.

So with h ≥ 1.5 ulp_above, all four points are distinct. Hmm, so maybe the "neighbouring points round to the same value" case can't actually arise for Central<4> at a power of two, because the rounding is monotone: x - h < x - ... hmm. Actually rounding is monotone non-decreasing, and for h ≥ 0.5 ulp on both sides, x-h < x < x+h. Then x-2h ≤ x-h and x+2h ≥ x+h, with equality only if the rounding collapses them. Rounding collapses x-h and x-2h only if they're within half an ulp of each other, i.e., h is below the local spacing. Since h ≥ 0.5 ulp(x) ≥ 0.5 × spacing/2... hmm, at a power of two, spacing below = 2×spacing above = 2^-52 and h ≥ 2^-53 = 0.5 × 2^-52 → x-2h and x-h differ by 2h ≥ 2^-52 = one ulp below ✓ so they're distinct? Rounding: x-2h and x-h differ by 2h ≥ ulp, so they round to different values (monotonicity: if a < b and b - a ≥ ulp, they round to different values... not exactly, but if b-a ≥ ulp then fl(a) < fl(b)? Hmm, fl is monotone non-decreasing, and fl(b) - fl(a) ≥ b - a - ulp ≥ 0; could be 0 if b-a is exactly ulp and both round down... hmm, e.g., a = 0, b = ulp: fl(0)=0, fl(ulp)=ulp, differ. a = 0.6ulp, b = 1.6ulp: fl = 1ulp, 2ulp differ. Hmm, a = 0.4ulp, b = 1.4ulp: fl(a) = 0, fl(b) = 1ulp differ. So yes, if b - a ≥ ulp, then fl(b) > fl(a) (since fl(a) ≤ a + ulp/2 and fl(b) ≥ b - ulp/2 ≥ a + ulp/2 > fl(a)... careful with the ≥/>, but generally true).

Hmm OK so it seems the "neighbouring points round to the same value" scenario is essentially the same as "x ± h rounds to x" for these rules. The book's phrasing might be paraphrasing the code's comment for Ridders (which checks x + 2h vs x + h for N ≥ 3). Hmm, that's the Ridders case, not the stencil case. So the book's failure list is a bit loose but not demonstrably wrong. I'll skip it.

OK. I think I have enough. Let me also do a final check of a couple of items in the scope list I haven't touched:

- "table-accuracy-matrix.tex" ✓ read.
- "accuracy-matrix.csv" ✓ read (partially; the tex matches the CSV medians for Central<2>: float 5.04e-6 → 5.0e-6 ✓, 1.378e-4 → 1.4e-4 ✓, 3.976e-4 → 4.0e-4 ✓, 0.013606 → 0.014 ✓, double 5.8296e-12 → 5.8e-12 ✓, 8.136e-9 → 8.1e-9 ✓, 1.385e-7 → 1.4e-7 ✓, 1.164e-5 → 1.2e-5 ✓, bin50 5.985e-35 → 6.0e-35 ✓). ✓ Consistent.

Hmm, chapter 7 line 106 says the model eps^(2/3) = 3.7e-11 "measured median 5.8e-12" ✓ = the accuracy-matrix Central<2> double N=1 median ✓.

- "getting-started.cmake" ✓ read (GIT_TAG aeec7e6 exists ✓).
- "code/README.md" ✓ read.
- "values.tex" ✓.
- "newton.csv" ✓, "kink.csv" ✓.
- figures 07-using-newton.tex, 07-using-kink.tex ✓.

Hmm, one thing about the newton figure: `xtick distance=1, xmin=-0.2, xmax=4.4` and the data goes to iteration 4 ✓. And `ytick={1e-12,1e-8,1e-4,1,1e4}` ✓.

Hmm, the figure's node "quadratic: digits double" — hmm, "digits double" = "the correct digits double". OK.

The kink figure: `unbounded coords=jump` — used because central2 at delta=0.005 is 6.99e-12 and the axis ymin=1e-15 ✓ fine (no unbounded coords actually). And ridders_evals is not plotted ✓.

Hmm, one issue: the kink figure's y-axis is log and includes values down to 1e-15; the ridders_estimate for delta=5 is 1.3e-13 ✓ within range.

Now let me think about whether there are any *major* issues I've missed. Let me re-read the chapter 8 text I haven't scrutinized: lines 44-60, 130-160, 158-213, 215-231, 232-260, 261-271.

Let me re-read line 46-60 and check the "signature" environment usage: the two declarations are printed with `\code{...}` and `\sigindent`. The second overload's line: `\sigindent\code{derivative(F &&fn, T const &x);}` — hmm, missing the return type! Let me look again:

```
\begin{signature}
\code{template<typename M = Central<2>, auto N = 1, real T, typename F>}\\
\sigindent\code{requires detail::method<M> && detail::valid_order<N>}\\
\sigindent\sigindent\code{&& detail::callable<F, T>}\\
\code{[[nodiscard]] constexpr std::expected<Estimate<T>, error>}\\
\sigindent\code{derivative(F &&fn, T const &x, std::type_identity_t<T> const &h);}\\[1ex]
\code{template<typename M = Central<2>, auto N = 1, real T, typename F>}\\
\sigindent\code{requires detail::method<M> && detail::valid_order<N>}\\
\sigindent\sigindent\code{&& detail::callable<F, T>}\\
\code{[[nodiscard]] constexpr std::expected<Estimate<T>, error>}\\
\sigindent\code{derivative(F &&fn, T const &x);}
\end{signature}
```
Both have the return type ✓. Good — I misread earlier.

Hmm, in the actual header the return type is on the same line as `[[nodiscard]] constexpr std::expected<...>` then the function name on the next line ✓ matches the book's layout.

OK.

Now — the book says (line 22) "function template: two overloads, each with a deprecated twin" ✓ (4 total).

Hmm, chapter 8 line 46: "Two overloads, with a step and with the default step." ✓

Now let me look for any remaining issue with the "Estimates" column in tab:api-methods: "Central<P> ... no; Forward<P> ... no; Backward<P> ... no; Ridders ... yes" ✓.

Hmm, the "Default step" column for Ridders: "$s/32$ ($N = 1$), $s/4$ ($N \ge 2$)" ✓.

And the Nodes column for Ridders: "Central<2> of order N at $h/2^i$: $\pm h_i$ for $N \le 2$, $\pm h_i, \pm 2h_i$ for $N \ge 3$, and $x$ for even N". Hmm — for N=2, the base rule is the three-point second difference using x, x±h — so the nodes include x ✓ ("and x for even N"). ✓ For N=4 the table is Central<2> of order 4 = nodes -2,-1,0,1,2 ✓.

Hmm, "± h_i for N ≤ 2": for N=2 the rule is the second difference at x±h plus x ✓.

OK, that's consistent.

Now let me double check the book's claim (ch8 line 210-212): "The two-node rules divide by the distance between the points actually evaluated, \code{Central<2>} for $f''$ by its actual offsets (\cref{sb:actual-distance,sb:actual-offsets}), the other rules by $\den h$ and then $N - 1$ times by $h$."

✓ matches the code: `detail::difference` divides by dx; the N==2 Central<2> branch uses actual offsets; others use `divide_by_powers<N, ...>(value, denominator*h, h)` → divides by (den*h) then N-1 times by h ✓ (the doc says "the evaluator divides once by denominator * h, then N - 1 times by h" ✓).

Hmm, wait: "divide once by denominator*h, then N-1 times by h" — that's N divisions total, giving h^N ✓ correct.

OK. Now let me spot-check chapter 7's §7.4.1 (Newton) numbers vs newton.csv once more:
- "Started at \SI{\valNewtonStart}{K}" = 800 ✓ (newton.csv iteration 0 T=800 ✓).
- "for the enthalpy of \SI{1000}{K}, it converges in \num{\valNewtonIterations} iterations" = 4 ✓.
- "to within \SI{\valNewtonFinalError}{K}" = 1.1e-13 ✓ (row 4 error = 1.1369e-13 ✓).
- "the final residual, \SI{\valNewtonFinalResidual}{\joule\per\mole}, is the rounding of $\enth$ itself" = 3.6e-12 ✓ (row 4 residual -3.638e-12 ✓).
- "Each iteration costs three evaluations of $\enth$" ✓.
- "Ridders' method would cost \num{\valCpRiddersEvals} per slope" = 14 ✓.

Hmm — the code's comment says "The text claims at most 10 steps and 1e-6 K" and checks `s->steps <= 10 && error <= 1e-6`. The text says "converges in 4 iterations to within 1.1e-13 K" ✓ consistent (4 ≤ 10, 1.1e-13 ≤ 1e-6).

Hmm, one thing: the chapter's Newton snippet uses `derivative(H, T)` — the *default* overload. But note that the text says "the method is the default Central<2>" ✓.

Hmm — is there an issue with the Newton iteration: the code computes `dT = (H(T) - H_target)/cp->value` and then `T -= dT`, and the convergence test is `|dT| < 1e-9`. The text says "stop when the step is below 1e-9 K" (code comment) — the text (ch7 line 143-144) says "A failed slope and a failure to converge are both returned as errors, so a caller cannot take an unconverged temperature for a solution" ✓.

Now let me check chapter 7's claim at line 154-155: "the final residual, 3.6e-12 J/mol, is the rounding of H itself". ✓

OK. Let me also verify the sensitivity example's claim (line 177-179) about evaluations: "The cost is the sum over the parameters: 20 evaluations, of which 4 for A (the rate is linear in A, so the first extrapolated entry agrees exactly with its parent and the method stops) and 16 for Ea."

Hmm — "the first extrapolated entry agrees exactly with its parent". Hmm, for a linear function, the Central<2> difference is exact at every level, so A[i][0] = exact for all i, and then A[i][1] = A[i][0] + (A[i][0] - A[i-1][0])/3 = A[i][0] exactly (since the change is 0). Hmm! But the `change == 0` check uses `change = distance(row[j], prev[j-1])` = |A[i][j] - A[i-1][j-1]|. For j=1: |A[1][1] - A[0][0]| = |exact - exact| = 0 ✓ → exact → break, best = A[1][1] = exact. So it stops after 2 levels = 4 evaluations ✓. So "the first extrapolated entry agrees exactly with its parent" — A[1][1] vs A[0][0]: they agree ✓. OK, the phrasing is a bit off (the extrapolated entry agrees with its *diagonal parent*, which is what the code compares) but fine.

Hmm, hold on: for Ridders, the check for the level j = 1 is `row[1]` vs `prev[0]` — and the code's comment says "An entry equal to its parent ends the tableau". And the loop is `for (std::size_t j = 1; j <= i; ++j)` where i is the level index. At level i=1, j=1: change = |row[1] - prev[0]| ✓. So yes, after 2 levels (i=0 and i=1) → 2 levels × 2 evaluations = 4 ✓.

Good, matches valJacobianEvalsA = 4.

Now — one more: the text says "of which 4 for A ... and 16 for Ea" and valJacobianEvals = 20 = 4+16 ✓.

OK. I'm now fairly confident about the findings list. Let me also double-check the "m × n calls" issue by re-reading the sentence in context (ch7 line 185-187):

"The library differentiates a function of one argument with one value, so it supports one way today: one \code{derivative} call per entry, on a lambda that evaluates the whole model with one unknown varied, at $m \times n$ calls, \num{\valJacobianBalanceCallEvals} here (snippet \code{jacobian-by-calls} in \code{code/07-using-jacobian.cpp})."

Hmm, m and n: "the Jacobian of several residuals" — m residuals, n unknowns. m×n = 4 derivative calls. But valJacobianBalanceCallEvals = 8 is the *model* call count. So the sentence conflates. Hmm — unless "m × n calls" is meant as "m×n derivative calls" and the macro is wrong... but the code asserts library_calls == 8. So the text's "at m × n calls, 8 here" is inconsistent. Hmm, hmm. Actually maybe the author intended "2 m n calls" (two model calls per entry). I'll report it.

Hmm, wait. Let me reconsider: maybe m and n in the book's notation are not "residuals × unknowns". Let me re-read: line 181-184: "Newton's method for several unknowns needs the Jacobian of several residuals, which a model usually returns from one call. Take an energy balance r1 = ... and an ideal-gas pressure r2 = ... in the unknowns T and n, at ... The library differentiates a function of one argument with one value, so it supports one way today: one derivative call per entry, on a lambda that evaluates the whole model with one unknown varied, at m × n calls, 8 here".

So m = 2 residuals, n = 2 unknowns → m×n = 4 entries. The number of derivative calls = 4, and the number of model calls = 8. The text says "at m × n calls, 8 here" → 8 ≠ 4. So there's a factor-2 mismatch. Given "one derivative call per entry ... at m×n calls", the phrase implies m×n derivative calls, and 8 is wrong for that (should be 4). OR the phrase means model calls and then it should be "2 m n". Either way, inconsistent. Good finding (minor).

Hmm, but careful: maybe "m × n" refers to something else, e.g. m = number of derivative calls... no.

Hmm, also compare line 186-187 with line 187's "(\num{\valJacobianBalanceCallEvals} here ...)". And later line 192: "Here that is \num{\valJacobianBalanceHandEvals} model calls" = 4 — explicitly "model calls". So in the hand case, the text says 4 model calls; in the library case it says 8 calls (and 8 is the model-call count per the code). So the text's "m × n calls" for the library case should be "2 m n model calls". Confirmed inconsistency. Good.

Now, let me write up. Let me also decide on the kink caption finding and the sensitivity "below its estimate" finding, the check-by-ridders "to within it" finding, the member-function-pointer finding, the "T is always deduced" finding, the "-> mp50 matters" finding, the "m×n" finding, and the ch8/ch7 numeric_limits list inconsistency.

Hmm, let me reconsider the sensitivity one. Text (ch7 line 168-170): "Against the closed forms ∂k/∂T = k Ea/(R T²) and ∂k/∂Ea = −k/(R T) the measured relative errors are 1.1e-13 and 1.2e-13, each below its estimate."

The estimates: valSensitivityRateEstimate = 3.9e-12 (absolute), valSensitivityRateRelEstimate = 1.4e-13 (relative!). Hmm! There IS a relative-estimate macro: valSensitivityRateRelEstimate = 1.4e-13 and valSensitivityEaRelEstimate = 3.8e-13. So the chapter could have used those. The measured relative errors are 1.1e-13 and 1.2e-13; the relative estimates are 1.4e-13 and 3.8e-13. So "each below its estimate" is TRUE if "its estimate" = the relative estimate (1.1e-13 < 1.4e-13 ✓; 1.2e-13 < 3.8e-13 ✓). And also true for absolute (3.2e-12 < 3.9e-12 ✓; 2.6e-14 < 8.3e-14 ✓). So the claim holds either way! Good — no finding. I'll drop it.

Hmm, let me double check valSensitivityEaRelEstimate = 3.8e-13 vs measured 1.2e-13 ✓ below ✓. And the code's check is `covered(dk_dT, exact_T) && covered(dk_dEa, exact_Ea)` with the 6× criterion ✓. So no finding. Good.

Hmm, but there's a subtlety: the text says "the measured relative errors are ... each below its estimate" — the natural reading is that each relative error is below its error estimate, which mixes relative and absolute. Since |value| = 29.0 for the first, the absolute estimate 3.9e-12 corresponds to a relative 1.34e-13 ≈ 1.4e-13 = valSensitivityRateRelEstimate ✓. So the claim is fine. Drop it.

Now the check-by-ridders one: I'll keep it but phrase it as "the printed numbers imply a difference of 1e-11, 3.7× the estimate quoted, so the claim cannot be checked from the text; the library's own test criterion is 6× the estimate". Hmm. Actually, let me reconsider once more whether to include it. The evidence: chapter line 260-261 says "the difference |d − r| = 7.7e-10 J/(mol K) is the stencil's error to within it; the actual error is 7.8e-10". The values I computed: |d − exact| = 7.7506e-10 (displayed 7.8e-10) and |d − r| = 7.7e-10 (2 s.f.). So |d−r| and |d−exact| differ by at most 1.0e-11 (and at least ~0). The estimate is 2.7e-12. So the text's claim "to within it" is only *possibly* true (it requires the true difference to be ≤ 2.7e-12, i.e. r within 2.7e-12 of exact). Given valCpRiddersActual = 2.0e-12 for the same call, it's probably true. Hmm. So the finding is weak: the claim is probably true but unverifiable from the printed numbers, and the library's own stated criterion (code comment line 99, "to within 6 times that estimate") is looser. Hmm — the text says "to within it", which is stricter than the code's comment. That IS a mismatch between the book and the code comment. I'll report it as minor with both pieces of evidence, noting the reader can't verify it from the two-digit numbers. That's honest and useful.

Hmm, but hold on — is the code comment's "6 times" a *guarantee* for all cases, or just the criterion the tests use? The comment says: "But when Ridders is covered (|r - exact| <= 6 * its estimate, the tests' criterion), the triangle inequality makes |d - r| the stencil's error to within 6 times that estimate." So the code's claim: if Ridders is covered within 6×, then |d−r| approximates the stencil error within 6× the estimate. The book says "within it" (1×). So the book is stricter than the code. But since Ridders' actual error is typically ≪ its estimate (the book says the median ratio is 0.14 for f' in the random study), "within it" is often true. Hmm — in this very case, the ratio is 2.0e-12/2.7e-12 = 0.74, so the true difference is 0.74 × the estimate ≤ 1× ✓ TRUE. Hmm! So the claim holds for this case.

Hmm. So it's not a factual error for this case; it's just not verifiable from the printed numbers. Hmm. I'll report it as "editorial/minor: the claim cannot be checked from the two-significant-digit numbers quoted (they imply a difference of 1e-11, 3.7× the estimate); the library's own criterion is 6× the estimate". Hmm — but the implied difference from rounded numbers isn't the true one. I need to be careful not to overstate. Let me phrase: "the reader cannot reproduce the claim from the quoted numbers: 7.8e-10 − 7.7e-10 = 1e-11 ≫ 2.7e-12". And note the true margin (2.0e-12 from valCpRiddersActual) is within the estimate only by a factor 1.35, while the library's own comment claims only 6×. Severity: minor/editorial.

OK. Now, one more area to check: chapter 7's §7.8 pitfall "The step's scale" claims: "at a small $|x|$ such as a mole fraction of $10^{-6}$, the step is larger than $x$ and the result is \code{non_finite} or crosses a physical bound" — hmm, at x = 1e-6, the default step is root_epsilon<3,double>() × max(1e-6, 1) = 6.055e-6 > x ✓ "the step is larger than x" ✓. And valSmallDefaultStep = 6.1e-6 ✓, valSmallDefaultStatus = non-finite ✓. So "the result is non_finite" — for a function defined on x>0 (like log), sampling x-h < 0 gives NaN → non_finite ✓. OK.

Hmm, the pitfall text mentions "crosses a physical bound" ✓.

OK, now let me consider one more possible finding: chapter 7 line 245-246: "The library has no accessor for its default step, so the check rebuilds it, \code{root_epsilon<P + N, T>() * max(|x|, 1)}". ✓ verified: there's no public accessor (detail::default_step is internal) ✓.

And chapter 8 line 95: "there is no accessor for it" ✓.

Now let me also double check chapter 7's statement (line 316): "Comparing the result with a value, as in \code{derivative(cube, 2.0, 0.25) == 12.0625}, compares \code{value} only, through the conversion, and is false when the result is a failure" ✓ verified (static_assert(!(bad == 2.0)) passed).

Alright. Let me also quickly check the front matter/quick-reference isn't in scope. Right.

One more thing to verify: the scope mentions "the claim that N is 1..4 and that a size_t/int/order is accepted". Chapter 7 line 28-29: "The derivative order N is 1 by default; it is an integer from 1 to 4 or a value of the enumeration order, such as order::second." Chapter 8 line 68-70 ✓. The header's valid_order accepts any integral type except bool, plus `order` ✓. Hmm — the book says "an integer from 1 to 4" and "as a value of any integer type except bool". ✓ And `size_t`: `std::integral<size_t>` ✓ and `N >= 1 && N <= 4` in its own type ✓. So `derivative<Central<2>, std::size_t(2)>` works ✓. Let me not test; the header comment confirms.

Hmm, hold on — one potential issue: `derivative<Central<2>, std::size_t(2)>(...)`: `valid_order<size_t(2)>` → `std::integral<size_t> && !same_as<size_t, bool> && 2 >= 1 && 2 <= 4` ✓ true. Then `static_cast<int>(N)` ✓. Fine.

Now let me write the final report. Let me organize by file:line with severity.

Let me collect the findings:

**Chapter 7**

1. `07-using.tex:187` (and 185-186): "at $m \times n$ calls, \num{\valJacobianBalanceCallEvals} here" — the macro is 8 (model calls), but m×n = 4 (derivative calls). Evidence: text line 185-187; code/07-using-jacobian.cpp:129 `hand_calls == 4 && library_calls == 8`; values.tex:454 `\valJacobianBalanceCallEvals = 8`. Severity: minor. Fix: "at $2mn$ model calls" or "at $m\times n$ calls, \num{...} model evaluations".

2. `07-using.tex:260-261`: "the difference $|d - r| = 7.7e-10$ ... is the stencil's error to within it; the actual error is 7.8e-10" — from the printed numbers the difference is 1e-11, 3.7× the estimate 2.7e-12 quoted in the same sentence; the library's own comment states the guarantee as 6× (code/07-using-check.cpp:99), and the test only checks `|diff - actual| <= 6 * estimate`. Severity: minor. Fix: quote more digits or state "to within a few times the estimate".

3. `07-using.tex:363-366` (figure caption in the chapter): "the error is a plateau of about $J / (2 h \cp)$ ... and \num{\valKinkInsideCentralFourRelError} for \code{Central<4>} ($h = 0.37$ K), so it falls as $1/h$ across the rules". Computed: J/(2h cp) for Central<4> = 0.7684/(2×0.37005×29.5814) = 0.0351, not the quoted 0.041; the exact plateau is 7J/(12h cp) = 0.04095 (Central<4> has two nodes below the break). And 4.289/0.04095 = 105 ≠ h ratio 122, so it does not fall exactly as 1/h across the rules. Severity: minor.

4. `07-using.tex:32-34` (footnote) and `08-api.tex:78`: "or a member function of \code{T} itself" / "so does a pointer to a member function of \code{T} itself" — only a *const* member function works; verified `detail::callable<R (R::*)(), R>` is false and `detail::callable<R (R::*)() const, R>` is true; the library's own test uses `constexpr Real square() const` (tests/test_callables.cpp:38). Severity: minor. Fix: "a pointer to a const member function of T".

5. `07-using.tex:286`: "give a lambda an explicit \code{-> mp50}: it matters for a type with expression templates" — no failure demonstrated: `derivative([](mpet const& x){ return exp(x); }, mpet(1))` is well-formed (verified with cpp_bin_float<50> et_on), and every call in the library casts immediately (deriv.hpp:804-805, 865-866, 963, 1031, 1200). Severity: minor/editorial. Fix: say it is a safe habit rather than that it matters.

**Chapter 8**

6. `08-api.tex:270`: "\code{T} is always deduced from \code{x}" — false: T can be given explicitly as the third template argument; verified `derivative<Central<2>, 1, float>(sq, 1.0)` compiles and returns a value. Severity: minor. Fix: "T is deduced from x, not from h; pass it third if you must".

7. `08-api.tex:248`: "The algorithms also read \code{radix}, \code{digits} and \code{denorm_min()}" — the list omits epsilon() and max(), which the concept also does not check and which the chapter's own snippet (line 252 → code/07-using-usertype.cpp:54-56) lists; chapter 7's sidebar (07-using.tex:292) gives all five. Severity: editorial.

Hmm, also is that right? The concept doesn't check epsilon()/max() existence. Right — `is_bounded` is just a bool constant. So a type could define is_bounded=true and no max() → the concept passes and the library fails. The book's line 243 says is_bounded true means "max() is finite and epsilon() positive" — which is about semantics, not existence. OK, editorial.

8. `08-api.tex:117-119`: "at most \num{\valRiddersMaxEvalsFirst}, ... times" ✓ no issue.

Hmm, anything else in chapter 8? Let me think about whether any of chapter 8's transcribed declarations mismatched. I found none. I should say so explicitly.

Hmm, let me double check one more: the book's §8.2 says "Each has a twin with the additional constraint \code{detail::noisy<N, T>}". ✓ And §8.5's static_assert messages: the book quotes ``\code{Central<P> needs an even order P from 2 to 10}'' and ``\code{Forward<P> needs an order P from 1 to 6}'' ✓ both verbatim from deriv.hpp:219, 234 ✓. And "this rule's integer weights are not exact in T" ✓ deriv.hpp:1010 (the full message is "this rule's integer weights are not exact in T: use a type with more significant digits" — the book quotes the prefix, fine).

Hmm, one more: chapter 8 line 181 says the tag's static_assert messages are ``\code{Central<P> needs an even order P from 2 to 10}'', ``\code{Forward<P> needs an order P from 1 to 6}''. It omits Backward's ("Backward<P> needs an order P from 1 to 6") — fine.

Now let me also double check §8.4's claim: "P is the accuracy order, the $\order(h^P)$ of the error, not the number of points." ✓

And §8.4: "A tag with an order outside its range is rejected by the constraint of \code{derivative}; completing the type itself fails a \code{static_assert}". ✓ verified by the 08-api-usage static_asserts compiling (Central<3>, Central<12>, Forward<7> rejected) ✓ — hmm, actually the example's static_asserts compile ✓ which confirms both claims.

Hmm, does `!accepted<Central<3>, 1>` hold? The example compiles in CI (presumably) and my -fsyntax-only of 08-api-usage.cpp passed ✓. So yes ✓.

Now, other possible findings in chapter 8:

- Line 40: "The enumeration \code{error} may gain values in a later version." ✓
- Line 55-56 of ch7: "\code{error} may gain values, so a \code{switch} on it needs a \code{default}." ✓

Hmm, what about chapter 8 line 156: "\code{r.value()} is not the member but the accessor of \code{std::expected}, which throws \code{std::bad_expected_access} when there is no result (or terminates without exceptions)". ✓ verified by the code's static_assert `std::same_as<decltype(r.value()), Estimate<double> const &>` ✓.

OK. Another: chapter 8 line 101-102: "An \code{Estimate<T>} (\cref{sec:api-estimate}) whose \code{abs_error} is filled by Ridders' method and empty for the stencils, or \code{error::non_finite}." ✓

Hmm — "error::non_finite" is the only error ✓.

Now let me also consider: does chapter 8 mention that `derivative` is `[[nodiscard]]`? Yes, in the signature ✓.

Does it mention that Estimate's operator== is a friend (so `Estimate == Estimate` works but not `Estimate == T`)? Line 152-154 ✓.

Hmm, what about `Estimate<T>`'s `operator T()` being `constexpr` — the book shows `constexpr operator T() const;` ✓.

OK. Let me also verify chapter 8's claim (line 111-113) that the untrusted result "is nevertheless returned, as a value with its large abs_error and not as a failure, when the step reaches the resolution of x or the levels run out during the restarts" ✓ (code lines 1314-1319: `if (!untrusted || resolution) return result; fallback = result; ... if (fallback) return *fallback;`) ✓.

Good.

Now let me think about whether there are any *major* findings I've missed. The scope asks to be adversarial about anything the target reader (an engineer) would misunderstand or be unable to act on. Let me re-scan chapter 7 for actionable issues:

- Line 243-251: the "Repeat with 2h and h/2" check. "While truncation dominates, halving the step quarters the error of Central<2>; while rounding dominates, it doubles it. Values that agree to many digits at h, 2h and h/2 are accurate to about that many." ✓ Reasonable.

- Line 250-251: "A kink or a jump inside the span does not survive this check: a jump $J$ adds $J / (2h)$ to the central difference, which doubles as the step halves." ✓ correct (J/(2h) doubles when h halves) ✓.

- Line 253-255: "Use a rule of higher order. derivative<Central<4>>(H, T0) has its own default step and a much smaller error, so its agreement with Central<2>, 10.6 digits here, measures the latter's error." ✓

- Line 263-268 ✓ (discussed).

- Line 270-271: "Finally, a physical check costs nothing: a heat capacity is positive..." ✓

- §7.6: ✓

- §7.7: ✓

- §7.8: ✓

Hmm, what about line 38-40: "A stencil calls fn once per node, so Central<2> costs two evaluations for f' and three for f''; Ridders' method costs two or four per level, with a median of 14 for f' and 32 for f''' in the random study of ch:ridders (tab:cost)." ✓ verified.

Hmm — "A stencil calls fn once per node" — hmm, for Central<2> f'' the nodes include 0 (x itself) → 3 evaluations ✓.

OK. And line 30-31: "The number type T is deduced from x, and a wider step converts to it, so 0.01 is a valid step for a 50-digit x" ✓.

Hmm, "and a wider step converts to it" — a *wider* step (double) converts to a narrower type... hmm, the step's type is `std::type_identity_t<T>`, so for a 50-digit T, a double step converts to T ✓ "wider step" meaning "the step you write is a double, which is wider"? Hmm, a double is *narrower* than a 50-digit type! So "a wider step converts to it" is backwards: the step is *narrower* (double) and converts to the wider T. Hmm! Let me re-read: "The number type T is deduced from x, and a wider step converts to it, so 0.01 is a valid step for a 50-digit x".

Hmm. I think the author means "a step of a wider/different type converts to T" — but for a 50-digit x, 0.01 is a double, which is *narrower*. Hmm, "wider step" could mean "a step written in a wider type" — no.

Hmm, hold on: maybe "wider" here means "a step with a wider *range* of types available" — no. I think the intended meaning is "a step of a different type converts to it" (i.e., the conversion is allowed). Since 0.01 is a double and T is a 50-digit type, calling the step "wider" is wrong — a double is narrower than cpp_bin_float_50 (53 bits vs 168). So "a wider step" is misleading/wrong.

Hmm, is it? For the *float* case: "for a float x write 1e-3f, since the narrowing of a double is diagnosed under -Wconversion" — there, a double step → float is a *narrowing* ✓. So the book's own next sentence says the double→float conversion narrows. So "a wider step converts to it" contradicts the next sentence for float (double is wider than float, so it narrows). Hmm — for float, a double IS wider, and it converts (narrowing). For a 50-digit type, a double is narrower and it converts (widening). So "a wider step converts to it" is true for float (the double step is wider) but false/misleading for the 50-digit case in the same sentence. Hmm, the sentence says "a wider step converts to it, so 0.01 is a valid step for a 50-digit x" — the "so" makes the 50-digit case the example of the "wider step" claim. So it's inconsistent: for a 50-digit T, 0.01 is a *narrower* value.

Hmm, this is the sort of thing an adversarial reviewer should catch, but is it a real error? The sentence is trying to say "a step of a different (double) type converts to T". "Wider" is wrong for the 50-digit example. I'd flag it as editorial/minor: "wider" is the wrong word; a double step is narrower than a 50-digit T (and the next sentence calls the same double→float conversion a narrowing).

Hmm, hmm. Let me reconsider: maybe the author means "a wider step type" as in "a step given in a type other than T"... Ambiguous. I'll report as editorial (low confidence on intent, but the wording is at best confusing given the next sentence).

Hmm, that might be too nitpicky. The task says "Be adversarial about anything the target reader (an engineer) would misunderstand". I'll include it briefly as editorial.

Hmm, I'd rather include the stronger ones. Let me count my findings: 7 solid. Let me also reconsider a few possible issues I may have skimmed:

- ch7 line 55-56: "It has no std::formatter: print static_cast<int>(d.error()), or a word chosen in that switch." ✓ — hmm, `d.error()` returns `error`; `static_cast<int>` ✓. Fine.

- ch7 line 68-70: "Third, r.value() is the accessor of std::expected, not the estimate's member: it compiles, converts to double and works until the first failure, where it throws std::bad_expected_access (or terminates a program built without exceptions)." ✓ verified by the static_assert in the code.

- ch7 line 71-72: "Fourth, d.value_or(0.0) does not compile, since no double converts to an Estimate<double>; d->abs_error.value_or(0.0) does." Hmm! `d` here — in the snippet, `r` is the expected; the prose says `d.value_or(0.0)`. Hmm, is `d` a typo for `r`? Let me re-read the snippet:

```
// <estimate-pitfalls>
auto const cosine = [](double x) { return std::cos(x); };
auto const r = derivative<Ridders>(cosine, 0.0);
if (r) { ... }
// r.value() is std::expected's accessor: ...
```
The prose (line 71-72) refers to "d.value_or(0.0)" and "d->abs_error.value_or(0.0)". But the snippet's variable is `r`, not `d`! Hmm. So the prose uses `d` where the snippet (and the rest of the paragraph, "r.value()") uses `r`. Hmm, and earlier in the chapter (line 73) the prose says "The rest of the book reads a result as d->value and Ridders' estimate as *r->abs_error". Hmm, so `d` is used as a generic name for a result elsewhere in the book. But within this paragraph about the snippet, using `d` for the variable named `r` is inconsistent — the reader looking at the listing sees `r`. Hmm, is that a finding? It's a small editorial inconsistency: the paragraph says "d.value_or(0.0) does not compile ... d->abs_error.value_or(0.0) does" while the listing binds the result to `r`. Hmm, also `r->abs_error.value_or(0.0)` — hmm, does that compile? `r->abs_error` is `std::optional<double>`, `.value_or(0.0)` ✓ compiles.

And does `r.value_or(0.0)` compile? `std::expected<Estimate<double>, error>::value_or(U&&)` — hmm, `value_or` requires the argument to be convertible to `Estimate<double>`; `0.0` (double) is not → so it doesn't compile ✓. The book says "since no double converts to an Estimate<double>" ✓ correct.

So the only issue is the variable name `d` vs `r`. Let me check: is it `d` in the prose? Line 71: "Fourth, \code{d.value_or(0.0)} does not compile, since no \code{double} converts to an \code{Estimate<double>}; \code{d->abs_error.value_or(0.0)} does." Hmm, and the preceding sentences use "r" ("Second, if (*r)..."; "Third, r.value() is the accessor..."). So yes, within the same paragraph the result object is called `r` in sentences 2 and 3 and `d` in sentence 4. Since the listing names it `r`, the `d` in sentence 4 is a slip. Minor editorial finding. Worth reporting? It's the kind of thing the target reader would notice ("there is no d in the listing"). I'll include it as editorial.

Hmm, hold on. Let me double check the actual text: line 71-72: "Fourth, \code{d.value_or(0.0)} does not compile, since no \code{double} converts to an \code{Estimate<double>}; \code{d->abs_error.value_or(0.0)} does. And \code{int n = *r} truncates, with a \code{-Wfloat-conversion} warning on GCC and Clang."

Yes — sentence 4 uses `d`, then the next sentence uses `*r`. Inconsistent. Minor editorial.

OK, that's another one. Let me verify by re-reading the exact lines from my earlier read: 

"65: First, \code{auto} deduces \code{Estimate<double>}, and \code{std::println} does not apply the conversion: print
66: \code{e.value}. Second, \code{if (*r)} tests the value against zero, not whether \code{r} holds a result: the
67: derivative of $\cos x$ at 0 is exactly 0, a good result that tests false. Clang with \code{-Wconversion -Werror}
68: rejects that line; GCC accepts it without a warning. Third, \code{r.value()} is the accessor of \code{std::expected},
69: not the estimate's member: it compiles, converts to \code{double} and works until the first failure, where it throws
70: \code{std::bad_expected_access} (or terminates a program built without exceptions). Read \code{r->value} after the
71: test. Fourth, \code{d.value_or(0.0)} does not compile, since no \code{double} converts to an \code{Estimate<double>};
72: \code{d->abs_error.value_or(0.0)} does. And \code{int n = *r} truncates, with a \code{-Wfloat-conversion} warning on
73: GCC and Clang. The rest of the book reads a result as \code{d->value} and Ridders' estimate as \code{*r->abs_error}."

Hmm, line 73 says "The rest of the book reads a result as \code{d->value} and Ridders' estimate as \code{*r->abs_error}." — here `d` is a generic name and `r` a Ridders result. So `d` in line 71 is... hmm, in the listing the Ridders result is `r`. So line 71's `d` should be `r`. Hmm, but maybe the author intends `d` generically = "the derivative result". In line 73 the generic name for a stencil result is `d`. In the listing, `r` is the Ridders result. So line 71's `d.value_or(0.0)` is a slip (should be `r`). Yes, minor editorial. I'll report it with the verbatim quote and the listing evidence.

Hmm, is it possible that the snippet had a `d` earlier? The snippet `<estimate-pitfalls>` has only `cosine`, `r`, `e`, `v`. So `d` is undefined in the listing. Confirmed editorial slip.

Now let me do a final check of a couple of chapter-7 §7.5 numbers I haven't reproduced: valCheckOrderDigits = 10.6 (digits(d, c4)) — requires Central<4>. Let me compute with Python: Central<4> at T0 with h = root_epsilon<5,double>()*1000 = 0.370048.

digits(d, c4) = -log10(|d - c4|/|c4|). Let me compute. My prediction: |d - c4| ≈ 7.75e-10 (d's error) → 10.63 ✓ = 10.6 ✓.

Let me verify quickly. Also valCheckRiddersDiff: I can't compute Ridders without reimplementing. Skip.

Let me compute the Central<4> digit count to be thorough.

Also — I realize I should double check the claim in ch7 line 249: "Here they agree to \num{\valCheckStepsDigits} digits, and the value agrees with the closed form to \num{\valCheckExactDigits} (measured)." Both 10.6 ✓ (my computed min 10.5534 → 10.6 ✓, and 10.6251 → 10.6 ✓).

Hmm, interesting: valCheckOrderDigits = 10.6 too, but the text at line 253-255 says "so its agreement with Central<2>, \num{\valCheckOrderDigits} digits here, measures the latter's error" ✓.

Let me compute Central<4> digits to confirm 10.6.

Let me run one more Python check and then write up.

Tags: v1.0.0 and v1.1.0-legacy exist, but no v2.0.0 tag. So the getting-started snippet pins a commit (aeec7e6) which is a book-development commit rather than a library release. Hmm — is that a finding? The snippet comment says "# a pinned commit or tag". The book says "With FetchContent, at a pinned commit or tag". Pinning the book's own stage-4 commit (aeec7e6) is odd: a reader fetching that commit gets the library as of that commit, which includes the book sources — fine. But it's not a release. Also note the repo's own VERSION is 2.0.0 and tags exist only for v1.0.0/v1.1.0-legacy. So there is no v2 tag to pin. Pinning a commit is the reasonable choice. Hmm, but pinning a commit that is *inside the book's own development history* (stage 4) rather than the latest library commit is slightly odd — a reader would get the library before review round 1's fixes. Hmm! That could matter: if the book's claims were fixed in stage 6 (5432abf), the pinned commit aeec7e6 predates those. But the library code (include/) presumably didn't change in stage 5-6 (only the book). Let me check whether include/ changed between aeec7e6 and HEAD.

If include/ changed, then pinning aeec7e6 gives a library that may differ from the one the book documents → a real (minor) finding. Let me check with git diff --stat aeec7e6..HEAD -- include/.

Let me run that.

Also let me double check: does the repo have a `main.cpp` for the getting-started snippet? Not needed.

Let me check the diff.

Hmm, but the snippet's comment says "# a pinned commit or tag" and the GIT_TAG is `aeec7e6` (a short SHA). That works with FetchContent ✓.

Hmm — one more consideration: the book says the library "ships no install rule and no find_package configuration, so a CMake project adds it as a subdirectory". The snippet uses FetchContent ✓. Fine.

Hmm, PLAN.md line 804 mentions "the CMake file (round 1, T4)" — item 23. Let me read that item to see if the CMake file has known caveats. Let me read PLAN.md around line 800-810.

Actually, let me check whether PLAN.md says the CMake listing is *not* checked against the repo (i.e., it might contain a fake URL). The snippet has `GIT_REPOSITORY https://github.com/troldal/Numerixx.git`. Is that the real remote? Let me check `git remote -v`.

Let me check.

Let me verify check.py actually does this (for the chapter-8 signature blocks), because chapter 8 claims they're "checked at every build" (per my task description: "it says they are transcribed from the header and checked at every build"). Let me check check.py.

Hmm — chapter 8's own text (line 8-9) says: "It describes the headers as they are, in their own wording where it summarises a comment of the header." And the chapter's header comment says "Every statement is checked against include/numerixx/deriv.hpp and include/numerixx/math.hpp." So the "checked at every build" is check.py's signature check. Let me read check.py to confirm what it checks — and specifically whether it checks *every* line of chapter 8's signature blocks (e.g., the `template<real T> struct Estimate` line vs the header's `template<real T>\nstruct Estimate`). Let me read check.py.

This matters: if check.py normalizes whitespace and strips comments, then the book's signature lines must occur verbatim (modulo whitespace) in the header. Let me verify a couple of the chapter-8 signature lines against the header text with the same normalization. E.g.:

- "template<typename M = Central<2>, auto N = 1, real T, typename F>" → header line 1400: `template<typename M = Central<2>, auto N = 1, real T, typename F>` ✓ (identical).
- "requires detail::method<M> && detail::valid_order<N>" + "&& detail::callable<F, T>" → header line 1401: `requires detail::method<M> && detail::valid_order<N> && detail::callable<F, T>` → with whitespace removed: "requiresdetail::method<M>&&detail::valid_order<N>" then "&&detail::callable<F,T>" — the concatenation of the two book lines (whitespace-stripped) = "requiresdetail::method<M>&&detail::valid_order<N>&&detail::callable<F,T>" which occurs in the header's line ✓ (if check.py joins the block's lines before searching).
- "[[nodiscard]] constexpr std::expected<Estimate<T>, error>" → header 1402 ✓.
- "derivative(F &&fn, T const &x, std::type_identity_t<T> const &h);" → header 1403: `derivative(F &&fn, T const &x, std::type_identity_t<T> const &h)` ✓ (the book adds a semicolon; with whitespace stripped the book's string has a trailing ';' which wouldn't occur in the header!). Hmm! That would fail check.py unless check.py strips semicolons too. Hmm.

Let me read check.py to see how it normalizes. If it strips ';' and whitespace, then the book's lines must still occur as substrings. Let me read.

This is worth checking because the chapter claims the declarations are verified at every build; if check.py can't match the book's lines (e.g. because of the added ';' or because the book reorders), then the claim is false. But if the build passes, it's fine. I can't run the build (it writes files). Let me read check.py.

So for chapter 8's signature blocks, each `\code{...}` line must appear in the header text with whitespace removed and trailing ';' removed.

Now let me verify chapter 8's signature lines against this check. Let me extract the \code lines from chapter 8's signature blocks and test each against the normalized headers using Python — I can do this without modifying anything (a read-only Python one-liner). That's a great verification of the "checked at every build" claim: if any line fails, the build would fail (a major finding); if all pass, I confirm the claim.

Let me write a Python script that replicates check.py's logic and reports the lines that fail. This is read-only.

Let me do it carefully:

```python
import pathlib, re
HERE = pathlib.Path('docs/deriv-book')
def normalised_code(text):
    text = re.sub(r'/\*.*?\*/', ' ', text, flags=re.S)
    text = re.sub(r'//[^\n]*', ' ', text)
    return re.sub(r'\s+', '', text)
headers = ''
for name in ('deriv.hpp','math.hpp'):
    headers += (pathlib.Path('include/numerixx')/name).read_text(encoding='utf-8', errors='replace') + '\n'
nh = normalised_code(headers)
for path in sorted(HERE.glob('chapters/*.tex')):
    text = path.read_text(encoding='ascii', errors='replace')
    for block in re.finditer(r'\\begin\{signature\}(.*?)\\end\{signature\}', text, flags=re.S):
        for piece in re.split(r'\\\\', block.group(1)):
            if piece.strip().startswith('//') or piece.strip().startswith('namespace'):
                continue
            code = normalised_code(piece).rstrip(';')
            if code and code not in nh:
                print('FAIL', path.name, repr(piece.strip()[:90]))
```
Note: the regex for the block — check.py's line 65 is truncated but likely `r'\\begin\{signature\}(.*?)\\end\{signature\}'`. And the split on `\\\\` (i.e. `\\` in the .tex). Note `\sigindent` is a macro, so the pieces are like `\code{...}\sigindent\code{...}` — hmm, the split is on `\\` (a literal double backslash in the tex = line break), so within a piece there may be `\sigindent\code{...}` which is not a `\code{...}` line... Hmm, check.py's truncated lines 65-67 probably extract `\code{...}` with a findall. Let me look at the actual lines 65-67. The read was truncated ("...[truncated 236 chars]..."). Let me re-read those lines specifically.

Let me get lines 63-73 verbatim.

Hmm, but note: the check only verifies that each `\code{...}` line *occurs* somewhere in the concatenated headers. That means a line could occur in a different place/order than the book shows (e.g., the book's ordering or an attribution could be wrong while each line matches). Since I read the header directly, I verified the ordering and context too. Fine.

One thing the check can't catch: `\code{enum class order : int { first = 1, second = 2, third = 3, fourth = 4 };}` — normalized → "enumclassorder:int{first=1,second=2,third=3,fourth=4}" — the header has `enum class order : int {\n first = 1, /**< f' */\n second = 2, ...}` → after comment stripping and whitespace removal: "enumclassorder:int{first=1,second=2,third=3,fourth=4}" ✓ occurs. Good.

And `\code{enum class error { non_finite };}` → "enumclasserror{non_finite}" ✓ occurs in the header ✓.

And `\code{template<real T> struct Estimate}` → "template<realT>structEstimate" ✓ occurs (header line 66-67) ✓.

And `\code{T value;}` → "Tvalue" ✓; `\code{std::optional<T> abs_error;}` → "std::optional<T>abs_error" ✓; `\code{constexpr operator T() const;}` → "constexproperatorT()const" ✓; `\code{friend constexpr bool}` → "friendconstexprbool" ✓; `\code{operator==(Estimate const &, Estimate const &) = default;}` → ✓.

And the math ones ✓, and the methods ✓ (`template<int P> struct Central;` → "template<intP>structCentral" ✓ occurs in the header line 216-217 ✓).

Great. So chapter 8's declarations are all verified.

Now let me finalize. But first, let me reconsider the "half-ulp" and a few other statements once more, and let me re-examine one thing in chapter 8 that I want to be sure about: line 94-95: "For a stencil of accuracy order P, $h = \epsmach^{1/(P+N)} \max(|x|, 1)$, computed as \code{root_epsilon<P + N, T>()} times the scale".

Hmm — the header computes `detail::make_table<N>(M{}).order + N` as the exponent, where `make_table<N>(M{}).order` is the table's *order* field. Is the table's `order` == P? For Central<P>, the table's order is P ✓ (the accuracy order). Let me double check the Table struct's `order` field. Let me grep for `order` in the Table struct (line 338-460).

Let me check that `Table::order` is the accuracy order P.

Also, chapter 8 line 94: "For a stencil of accuracy order P" — and tab:default-steps' "exponent 1/(P+N)" ✓ consistent.

Let me verify Table::order.

Also one more: chapter 8 line 198: "the $P + N - 1$ integers around 0" ✓.

Let me grep the Table struct.

Hmm, but is 6222216 < 2^24 = 16777216 ✓ and even < 2^23 = 8388608 ✓. The header says "must fit 24 bits" and "keeps it exact in float" (float has 24 bits) ✓. And `holds_integers<float>(6222216)` → 2^24 > 6222216 ✓ true ✓.

OK. Now, I'm done with the deep checks. Let me compile the final report.

Let me reconsider each finding's severity and make sure the evidence is verbatim.

FINDINGS:

### Chapter 7 (07-using.tex)

**F1. 07-using.tex:186-187 — minor.** "at $m \times n$ calls, \num{\valJacobianBalanceCallEvals} here".
- Evidence: line 185-187 verbatim: "so it supports one way today: one \code{derivative} call per entry, on a lambda that evaluates the whole model with one unknown varied, at $m \times n$ calls, \num{\valJacobianBalanceCallEvals} here (snippet \code{jacobian-by-calls} in \code{code/07-using-jacobian.cpp})."
- Counter-evidence: with m = 2 residuals and n = 2 unknowns, m×n = 4 derivative calls; the macro is 8 = the number of *model* calls (values.tex:454 `\newcommand{\valJacobianBalanceCallEvals}{8}`; code/07-using-jacobian.cpp:129 `bool const ok = ... && hand_calls == 4 && library_calls == 8;`, where `calls` counts model invocations).
- Fix: "at $2mn$ model calls, \num{\valJacobianBalanceCallEvals} here" (or "$m\times n$ calls, each two model evaluations").

**F2. 07-using.tex:71 — editorial.** "Fourth, \code{d.value_or(0.0)} does not compile... \code{d->abs_error.value_or(0.0)} does."
- Evidence: the listing the paragraph discusses (code/07-using-expected.cpp:88-94, snippet `<estimate-pitfalls>`) binds the result to `r`; sentences 2 and 3 of the same paragraph use `r` (`if (*r)`, `r.value()`), and the next sentence uses `*r`.
- Fix: use `r` (or rename the snippet's variable).

**F3. 07-using.tex:260-261 — minor.** "the difference $|d - r| = \SI{\valCheckRiddersDiff}{\joule\per\mole\per\kelvin}$ is the stencil's error to within it; the actual error is \SI{\valCheckStencilActual}{...}"
- Counter-evidence: the two macros are 7.7e-10 and 7.8e-10 (values.tex:216-217), so the quoted numbers imply a difference of 1e-11 = 3.7 × the estimate 2.7e-12 quoted in the same sentence (values.tex:218). The library's own criterion is 6× the estimate: code/07-using-check.cpp:99 "the triangle inequality makes |d - r| the stencil's error to within 6 times that estimate" and line 101 `std::abs(diff - actual) <= 6 * *r->abs_error`.
- My reproduction: |d − exact| = 7.7506e-10 (so 7.8e-10) and Ridders' actual error from DS1b is 2.0e-12, so the true margin is 2.0e-12 vs the 2.7e-12 estimate — the claim holds, but only by a factor 1.35, and cannot be checked from the printed digits.
- Fix: quote one more digit, or say "to within a few times its estimate (the tests allow six)".

**F4. 07-using.tex:363-366 — minor.** Figure caption formula.
- Evidence: "While the break lies inside a method's span, the error is a plateau of about $J / (2 h \cp)$, $J$ the jump in $\enth$: \num{\valKinkInsideCentralRelError} for \code{Central<2>} ($h = \SI{\valKinkCentralStep}{K}$) and \num{\valKinkInsideCentralFourRelError} for \code{Central<4>} ($h = \SI{\valKinkCentralFourStep}{K}$), so it falls as $1/h$ across the rules, whatever their order."
- Counter-evidence (computed from code/shomate.hpp's coefficients and the default steps): with J = 0.7684 J/mol, cp = 29.5814 J/(mol K), h(Central<2>) = 0.0030277, h(Central<4>) = 0.370048: J/(2h cp) = 4.2896 ✓ (matches 4.3) but 0.0351 for Central<4>, not the quoted 0.041; the exact plateau for Central<4> is 7J/(12h cp) = 0.04095 (two of its four nodes lie below the break). The ratio of the two plateaus is 105, while the ratio of the steps is 122, so it does not fall exactly as 1/h across the rules.
- Fix: state the plateau as $J/(2h c_p)$ for the two-node rules and $7J/(12h c_p)$ for \code{Central<4>}, or drop the general formula.

**F5. 07-using.tex:32-34 (and 08-api.tex:78) — minor.** "or a member function of T itself" / "so does a pointer to a member function of \code{T} itself".
- Counter-evidence: verified by compiling: `static_assert(nxx::deriv::detail::callable<R (R::*)() const, R>)` passes and `static_assert(!nxx::deriv::detail::callable<R (R::*)(), R>)` passes; `std::invocable<R (R::*)(), R const&>` is false. The library's own test uses a const member: tests/test_callables.cpp:38 `[[nodiscard]] constexpr Real square() const { return Real(v * v); }` and line 29 "a member-function pointer that std::invocable<F &, Real const &> accepts".
- Fix: "a pointer to a const member function of T".

**F6. 07-using.tex:286 — minor/editorial.** "And give a lambda an explicit \code{-> mp50}: it matters for a type with expression templates, such as \code{number<cpp_bin_float<50>, et_on>}".
- Counter-evidence: `requires { nxx::deriv::derivative([](mpet const& x){ return exp(x); }, mpet(1)); }` is satisfied for mpet = `number<cpp_bin_float<50>, et_on>` (compiled, Boost 1.92, GCC 16), and every call in the library casts the result immediately: deriv.hpp:804-805, 865-866, 963, 1031, 1200 (`T const f_lo = static_cast<T>(std::invoke(fn, x_lo));` etc.), so no unevaluated expression outlives the call.
- Fix: say it is a safe habit (as the same sentence already does for cpp_bin_float_50) rather than that it matters.

**F7. 07-using.tex:29-30 — editorial.** "The number type T is deduced from \code{x}, and a wider step converts to it, so \code{0.01} is a valid step for a 50-digit \code{x}".
- Counter-evidence: `0.01` is a double, which is *narrower* than a 50-digit T (53 bits vs 168: values.tex:316-317 `\valEpsMultiprecision{5.3e-51}`, `\valDigitsMultiprecision{168}`), and the next sentence calls the same double→float conversion a narrowing ("since the narrowing of a \code{double} is diagnosed under \code{-Wconversion}").
- Fix: "a step of another type converts to it".

### Chapter 8 (08-api.tex)

**F8. 08-api.tex:270 — minor.** "\code{T} is always deduced from \code{x}".
- Counter-evidence: T is the third template parameter and can be given explicitly; verified `static_assert(nxx::deriv::derivative<nxx::deriv::Central<2>, 1, float>(sq, 1.0).has_value());` compiles (x is a double, converted to float).
- Fix: "T comes from x, not from h; give it third only if you mean to".

**F9. 08-api.tex:248 — editorial.** "The algorithms also read \code{radix}, \code{digits} and \code{denorm_min()}, which the concept does not check".
- Counter-evidence: the library also reads `epsilon()` and `max()` (deriv.hpp:698, 554, 573, 589, 605, 624, 639, 661-662, 680, 1156), and neither is checked by the concept; the chapter's own snippet (08-api.tex:252 → code/07-using-usertype.cpp:54-56) and Chapter 7's sidebar (07-using.tex:292) list all five.
- Fix: add epsilon() and max().

Hmm, F9 — is it a real problem? "The algorithms also read radix, digits and denorm_min(), which the concept does not check: a type without them passes the constraints and then fails to compile inside the header." A reader could think max() and epsilon() ARE checked by the concept (they're not). Since the paragraph's point is "what the concept does not check", the omission is misleading. Keep as editorial.

Now, should I report anything as MAJOR? Let me think about which of these an engineer would be materially harmed by:
- F5 (member function pointer): an engineer writing `derivative(&Model::H, T0)` with a non-const H gets a wall of template errors. The book says it works. That's actionable and wrong → I'd call it minor-to-major; I'll say minor (it's one word: const).
- F1: an engineer counting calls would be off by 2×. Minor.
- F4: a wrong formula in a figure caption. Minor.
- F8: minor.

No major findings. Hmm. Let me double check a couple of "major" candidates I might have dismissed too quickly:

(a) The claim in ch7 line 55: "It has no \code{std::formatter}: print \code{static_cast<int>(d.error())}" ✓ fine.

(b) ch8 line 217-219: "A fourth derivative in a type whose machine epsilon is above $10^{-7}$, tested with integers as $\beta^{\text{digits} - 1} < 10^7$". Hmm — the code: `for (int i = 1; i < digits; ++i) { power *= radix; if (power >= 1e7) return false; } return true;` So power ends at radix^(digits-1) (i goes 1..digits-1, that's digits-1 multiplications) ✓. And returns true iff radix^(digits-1) < 1e7 ✓ (it returns false as soon as power >= 1e7; if it never reaches 1e7, true). ✓ The book's formula matches ✓.

Hmm, one subtlety: for a type with digits = 1 (not possible for a real type)... skip.

(c) ch8 line 220-221 quotes the warning verbatim: "the 4th derivative in a type with 7 or fewer significant digits is mostly rounding noise: use double, or silence this warning at the call or with NUMERIXX_NO_PRECISION_WARNING" ✓ matches deriv.hpp:1452-1453 ✓ (verbatim, modulo the line break).

(d) ch8 line 222: "The warning is reported at the caller's line, so a pragma there silences it (\code{-Wdeprecated-declarations} on GCC and Clang, C4996 on MSVC), and a system-header include does not hide it." ✓ matches deriv.hpp:193-194 and 1436-1437 ✓.

(e) ch7 line 44: "The library itself throws nothing: a failure is a value. An exception thrown by fn or by the operations of T passes through unchanged, and derivative is not noexcept." ✓

(f) ch7 line 75-77: "A step on the wrong scale, a function less accurate than epsmach or a break within a step of x give a wrong value without a failure" ✓.

(g) ch8 §8.2 line 89-91 — the half-ulp claim ✓ verified above.

(h) One more: ch7 line 30: "a wider step converts to it" (F7).

(i) ch7 §7.7 line 326-329: "The compiler and the processor compute the same value, provided the compiler neither contracts a * b + c into a fused multiply-add nor uses extended precision; the build pins neither, and the identity is verified on the five toolchains of the continuous integration (\cref{sec:guarantees-guarantees})."

Hmm — "the build pins neither" ✓ (no -ffp-contract=off in the CMake, presumably). Let me verify quickly that the build doesn't set -ffp-contract. Let me grep CMakeLists for fp-contract / fast-math. Hmm, PLAN.md D-13 (line 903) confirms: "The identity of compile-time and run-time results is not pinned by the code or the build... the build sets no contraction mode" ✓. So the claim matches the repo's own record ✓.

Hmm, and the code 07-using-constexpr.cpp:62 checks `run_time && *run_time == *cp` at run time ✓ (bit-for-bit). Fine.

(j) ch7 line 17-19 toolchains ✓ verified against ci.yml.

(k) ch8 line 38: "Both are header-only, need C++23 and no other library" ✓.

(l) ch8 line 39-40: "Everything in a namespace detail is internal: it may change without notice, and the book reads it only to generate the tables of ch:tables." Hmm — the book *also* reads `detail::noisy` and `detail::method` etc. in §8.2 and the examples use `nxx::deriv::detail::callable`? Hmm: chapter 8's §8.2 text mentions `detail::method<M>`, `detail::valid_order<N>`, `detail::callable<F, T>`, `detail::noisy<N, T>` — that's reading detail names in the text (not "to generate the tables"). Hmm! The sentence says "the book reads it only to generate the tables of \cref{ch:tables}" — but chapter 8 names several `detail` entities in the reference (and 08-api-usage.cpp uses `nxx::math::root_epsilon` — that's not detail; and 07-using-jacobian.cpp uses `nxx::math::root_epsilon` too). Hmm, does any *example* use `nxx::deriv::detail::...`? Let me check: grep the code dir for "detail::".

If an example uses a detail name, that contradicts "the book reads it only to generate the tables". Let me check.

Also the chapter's own signature block includes `requires detail::method<M> && detail::valid_order<N> && detail::callable<F, T>` — that's a transcription of the declaration (necessary), and §8.2's prose explains them. So the sentence "the book reads it only to generate the tables" is at least incomplete: chapter 8 quotes `detail::` names in its transcribed declarations and explains `detail::callable`, `detail::noisy`, `detail::valid_order`, `detail::method`. Hmm. Is that a finding? It's a self-description inconsistency: the book does name internal names outside the tables. An engineer reading "Everything in a namespace detail is internal: it may change without notice, and the book reads it only to generate the tables" might conclude that `detail::callable` etc. are not part of the documented API and skip them — but chapter 8 then relies on them for the rejection behaviour (§8.6). Hmm, mild. I could report as editorial. Let me check whether any code example uses a detail name (that would strengthen it).

Let me grep.

(m) One more: ch8 line 265: "A type too narrow for a rule's integer weights is rejected by a \code{static_assert} inside the library instead (\cref{sec:api-methods}), which a \code{requires}-expression does not see." ✓ verified.

Let me do the grep for detail:: in the code dir.

Hmm, but the sentence might be intended as "the book reads [the internals] only to generate the tables" meaning the *other* chapters. Ugh. I'll include it as a low-severity note since the reader could be misled into thinking `detail::` names are not documented/stable... Actually the sentence already says they may change without notice, which is the actionable part. Hmm. I'll include it as editorial but flagged as "the chapter's own §8.2 does the opposite".

Hmm, let me reconsider. The task says "Do not restate the book's own claims as findings" and "Report ONLY findings you have personally verified". This one is verified (verbatim lines exist). I'll include it, ranked last.

Now — no fast-math/ffp-contract flags in the build ✓ (grep found nothing) so "the build pins neither" ✓ correct.

Let me also double check one claim I want to be sure about: ch7 line 326-327 "The compiler and the processor compute the same value, provided the compiler neither contracts a*b+c into a fused multiply-add nor uses extended precision; the build pins neither, and the identity is verified on the five toolchains of the continuous integration". ✓ OK.

Now let me write the final answer. I should be concise but dense, with file:line, verbatim line, counter-evidence, severity, fix.

Before that, let me do a final sanity check on two numbers I haven't verified: valCpDefaultRelError = 2.4e-11 (used in ch7 line 157: "a δ near 2.4e-11 is far too small to slow the convergence"). My computed relative error of Central<2> at the default step on the thread's cp: 2.362e-11 → 2.4e-11 ✓ (values.tex:12 = 2.4e-11 ✓). And DS1b's CpDefaultRelError ✓.

Hmm, and ch7 line 156: "A slope that is off by a relative δ adds a term of about δ times the error to each step, \cite{dennis1996}". Dennis & Schnabel chapter 5 discusses the effect of an inexact derivative on Newton's method ✓ (per refs.bib's comment "the effect of an inexact slope on Newton's method in chapter 5; cite the section at the sentence"). Hmm! refs.bib says "cite the section at the sentence" — i.e., the citation should name the section. The chapter cites `\cite{dennis1996}` without a section. Hmm — refs.bib's own instruction says "Sections and theorems are cited at the citing sentence, not in the entry" (line 7) and the dennis1996 comment says "The scaled step is in section 5.4, the effect of an inexact slope on Newton's method in chapter 5; cite the section at the sentence." So the entry asks for a section reference at the citing sentence; the chapter doesn't give one (e.g. "chapter 5" or "§5.4.2"). Hmm! That's a traceability nitpick. Is it worth reporting? The refs.bib comment is a directive to the authors. The chapter says "A slope that is off by a relative δ adds a term of about δ times the error to each step,\cite{dennis1996}". The claim is supported by Dennis & Schnabel's chapter 5 (they show that the Newton step with an approximate derivative gives a linear term of order δ). ✓ supported. But the section should be cited per the book's own rule. I'll mention as editorial.

Hmm, also `\cite{more2012,gill1983}` for the deferred "noise-aware step" — more2012 ("Estimating derivatives of noisy simulations") supports it ✓. gill1983 ("Computing forward-difference intervals for numerical optimization") — it's about choosing the difference interval given the function's noise/accuracy; supports the noise-aware step claim ✓ (refs.bib doesn't say it's unverified). Hmm — is gill1983 relevant to a *noise-aware step*? Gill et al. (1983) derive the optimal interval h = (something)^(1/2) using an estimate of the relative error in f. Yes, that's a noise-aware step ✓ supported.

And the citation is attached to "a noise-aware step;" while the next list item (complex step / AD) has no citation — but those are cited in Chapter 1's sidebar (squire1998, martins2003, griewank2008) ✓ so not a defect here.

And colebrook1939 ✓ supports the Colebrook equation (the book says "in its common form" because the printed constants weren't seen — honest).

And nocedal2006, griewank2008, higham2002 are not cited anywhere in chapters 7-8 (verified by grep) — I'll note that as "not used in these chapters, so nothing to check; higham2002/griewank2008 are cited elsewhere".

Hmm, hold on: is that right? Let me double check by grepping the whole book for these keys, so my statement "not used in ch7/8" is precise and I can say where they are used. Let me do that quickly.

Let me also double check the "five toolchains" claim and the kink/nist-webbook citation ✓.

One more: ch7 line 203 cites colebrook1939 for the equation, but the *solver* (bisection) isn't cited — fine.

Let me run the final grep for citations across the book.

The four citations in scope:
- dennis1996 (ch7:156) — supports "a slope off by a relative δ adds about δ times the error to each step" ✓ (D&S ch. 5 on inexact derivatives in Newton's method). The refs.bib entry asks for the section at the sentence ("cite the section at the sentence"); the chapter cites the book alone. Editorial.
- colebrook1939 (ch7:203) — supports the Colebrook equation ✓ (with the "common form" hedge documented in refs.bib).
- nist-webbook (ch7:356) — supports the two Shomate ranges of nitrogen ✓.
- more2012, gill1983 (ch7:432) — support a "noise-aware step" ✓ (Moré & Wild 2012 §4 treats the central difference with noise; Gill et al. 1983 computes forward-difference intervals).

So the citations are OK, with the one editorial note about the section.

Now let me write the final report. I'll structure it:

1. Scope summary + what I verified (declarations all match; check.py's signature check reproduced → 0 failures).
2. Findings, ordered by severity.
3. Areas with nothing wrong.

Let me get exact line numbers and verbatim text for each finding. Let me re-verify the verbatim lines:

F1: 07-using.tex:185-187. Let me quote line 187 exactly: "be called with a \code{T const} lvalue..." no wait, that's ch8. Let me get ch7 lines 185-187 exactly:

185: "one argument with one value, so it supports one way today: one \code{derivative} call per entry, on a lambda that"
186: "evaluates the whole model with one unknown varied, at $m \times n$ calls,"
187: "\num{\valJacobianBalanceCallEvals} here (snippet \code{jacobian-by-calls} in"

Hmm, from my earlier read:
"184: $V = \SI{1}{\metre\cubed}$, with $Q$ and $P$ such that both vanish there. The library differentiates a function of
185: one argument with one value, so it supports one way today: one \code{derivative} call per entry, on a lambda that
186: evaluates the whole model with one unknown varied, at $m \times n$ calls,
187: \num{\valJacobianBalanceCallEvals} here (snippet \code{jacobian-by-calls} in"

✓ So the verbatim line 186-187 is "evaluates the whole model with one unknown varied, at $m \times n$ calls, \num{\valJacobianBalanceCallEvals} here".

F2: 07-using.tex:71: "test. Fourth, \code{d.value_or(0.0)} does not compile, since no \code{double} converts to an \code{Estimate<double>};" and 72: "\code{d->abs_error.value_or(0.0)} does. And \code{int n = *r} truncates, with a \code{-Wfloat-conversion} warning on".

F3: 07-using.tex:260-261: "belongs to Ridders' value, not to the stencil's. But when it is small, the difference $\abs{d - r} = \SI{\valCheckRiddersDiff}{\joule\per\mole\per\kelvin}$" / "is the stencil's error to within it; the actual error is \SI{\valCheckStencilActual}{\joule\per\mole\per\kelvin}."

F4: 07-using.tex:363-366.

F5: 07-using.tex:33: "function, a lambda, a function object, or a member function of T itself. For a member function of your model," and 08-api.tex:78: "does a pointer to a member function of \code{T} itself. For a member function of your own model, capture the".

Hmm, 08-api.tex line 77-79:
"77:     lvalue through \code{std::invoke}, with a \code{T const} lvalue, so a callable that rejects rvalues works, and so
78:     does a pointer to a member function of \code{T} itself. For a member function of your own model, capture the
79:     object in a lambda."
✓.

F6: 07-using.tex:286-287: "Boost type. And give a lambda an explicit \code{-> mp50}: it matters for a type with expression templates, such as" / "\code{number<cpp_bin_float<50>, et_on>}, and a safe habit for \code{cpp_bin_float_50} (\cref{sb:expression-templates})."

F7: 07-using.tex:29-30: "from \code{x}, and a wider step converts to it, so \code{0.01} is a valid step for a 50-digit \code{x}; for a".

F8: 08-api.tex:270: "be called with a \code{T const} lvalue..." no. Line 269-271:
"269: A type given explicitly as the second template argument, as in
270: \code{derivative<Ridders, double>(fn, x)}, does not compile either, because that argument is the order; \code{T} is always deduced from \code{x}. A callable that cannot
271: be called with a \code{T const} lvalue, or whose result does not convert to \code{T}, fails \code{detail::callable}."
✓ So the verbatim claim is on line 270.

F9: 08-api.tex:248: "algorithms also read \code{radix}, \code{digits} and \code{denorm_min()}, which the concept does not check: a type"

F10: 08-api.tex:39-40: "lists the public names. Everything in a namespace \code{detail} is internal: it may change without notice, and the" / "book reads it only to generate the tables of \cref{ch:tables}."

F11 (editorial): 07-using.tex:156 cite without section (per refs.bib:200 "cite the section at the sentence").

Let me double-check F11's exact wording: "A slope that is off by a relative" (155) "$\\delta$ adds a term of about $\\delta$ times the error to each step,\\cite{dennis1996} and a $\\delta$ near" (156).

OK. Also, I should report the positive verification results: (a) all chapter-8 declarations match the headers, including every detail the task asked about; (b) check.py's signature check reproduces with 0 failures; (c) all ten chapter-7 examples and 08-api-usage.cpp compile cleanly with -Wall -Wextra (GCC 16, C++23) — and I should note 07-using-multiprecision.cpp compiles with Boost 1.92 standalone; (d) the CSVs/values.tex reproduce from first principles for the check section, the kink, the Jacobian steps, the sensitivity closed forms, the multiprecision epsilon, the default steps; (e) the diagnostic claims (-Wconversion on clang, -Wfloat-conversion on gcc) hold; (f) the deprecated-twin/noisy rule and the float16 requires-expression claim hold; (g) citations support their claims.

Also I should mention the areas checked with nothing wrong: chapter 8's §8.3 Estimate semantics, §8.4 method ranges/static_asserts, §8.5 warning rule, §8.7 real/math, §8.6 rejected calls.

Let me now also double check one loose end: the task asks about "the `noisy`/deprecation rule" — I verified. And "`root_epsilon<N,T>`'s N>=2 requirement" — book line 256: "\code{root_epsilon<N, T>()} returns $\epsmach^{1/N}$ for $N \ge 2$ (a \code{static_assert})" ✓ matches deriv.hpp:80 `static_assert(N >= 2);` ✓. Note the book's signature block shows `template<int N, real T> constexpr T root_epsilon();` ✓ matches.

Hmm, one more thing: the book says root_epsilon is "accurate to within an ulp" ✓ (header line 69 "Accurate to within an ulp") ✓.

Also: "math::abs/isfinite" ✓.

And "the constexpr claims (which operations are rejected by constant evaluation)" — ch7 line 320-322: "The one exception is an overflow of the final quotient, for example a fourth derivative with a tiny step: at run time it is non_finite, but GCC rejects it at compile time." ✓ matches the header's comment and test_stencils. And ch7 line 332-340 (sidebar) ✓ — let me double check the sidebar's claim "the input checks keep every point $x + k h$ at least about $K h$ below $\fmax$, with $K$ the largest node (a NaN $x$ fails them too)" ✓ matches `fits`'s doc (deriv.hpp:536-543) ✓. And "every difference and partial sum is tested before it is formed, by comparing $a/\beta + b/\beta$ with $\fmax/\beta$" ✓ matches sum_fits (line 661-680) ✓. And "the rounding floors of Ridders' method saturate at $\fmax$" ✓ (saturating_sum/saturating_quotient) ✓. And "The exception is the final quotient, with the repeated divisions by $h$ for $\dorder \ge 2$" ✓ (divide_by_powers<N, false>) ✓.

Good.

One more: ch7 line 313-314: "With a step that is a power of two, $x^3$ gives exact results in binary arithmetic, which makes the truncation error of \cref{ch:taylor} visible to the compiler" ✓ (12.0625 = 12 + h² with h=0.25 ✓ verified by arithmetic).

Alright, I'm confident. Let me write the report now.

Actually, hold on. Let me reconsider F3 once more, because I want to be accurate about whether the claim is false or just unverifiable. Let me try to compute Ridders' value for cp at T0 exactly by reimplementing the algorithm in Python. It's ~40 lines. Hmm. The algorithm:
- N=1, f0 not used (N odd).
- level: difference(fn, x-h, x+h) → value and rounding_floor.
- tableau: A[i][j] = A[i][j-1] + (A[i][j-1] - A[i-1][j-1])/(4^j - 1).
- best entry: min estimate, estimate = max(change, level_floor).
- abs_error = max(best_error, |best - diagonal|).
- untrusted → restart, etc.

That's a fair amount of work (the floor computations involve value_rounding, etc.). The level_floor affects only the estimate, not the value. The value: the tableau stops when two levels stall or an entry equals its parent.

Honestly, the value r can be approximated: for cp (a smooth function), Ridders converges to cp + tiny. The book's DS1b says the Ridders value is 32.6917300000 with an actual error of 2.0e-12. So r = exact ± 2.0e-12. And |d - exact| = 7.7506e-10, so |d - r| = 7.7506e-10 ∓ 2.0e-12 → 7.7306e-10 (if r > exact) or 7.7706e-10 (if r < exact).

Hmm, and the macro valCheckRiddersDiff = 7.7e-10. If |d−r| = 7.7706e-10, the generator's `sig()` would print 7.8e-10, not 7.7e-10. So r > exact (r − exact = +2.0e-12)? Hmm, wait: I need to be careful with signs.

d = 32.69172999922494 (below exact by 7.7506e-10).
If r = exact + 2.0e-12 = 32.691730000002, then d − r = −7.7706e-10 → |d−r| = 7.7706e-10 → sig → 7.8e-10 ✗ (macro says 7.7e-10).
If r = exact − 2.0e-12 = 32.691729999998, then d − r = −7.7306e-10 → |d−r| = 7.7306e-10 → sig → 7.7e-10 ✓.

So r = exact − 2.0e-12 → |r − exact| = 2.0e-12 = valCpRiddersActual ✓✓ consistent!

Then | |d−r| − |d−exact| | = |7.7306e-10 − 7.7506e-10| = 2.0e-12, and the estimate is 2.7e-12 → 2.0e-12 ≤ 2.7e-12 ✓ the claim "to within it" is TRUE.

So F3 is a false positive in substance; the only issue is that the reader can't verify it from the printed digits. Hmm. I'll report it as editorial at most, or maybe just note it. Hmm.

Given the instruction "Report ONLY findings you have personally verified", and this claim turns out to be true (just not checkable from the quoted 2-s.f. numbers), I should either drop it or report it as an editorial remark about the two-digit macros making the claim unverifiable. I'll include it as editorial with a clear statement that the claim itself is consistent with the data (I derived r = exact − 2.0e-12), but the printed numbers suggest otherwise (7.8 − 7.7 = 0.1e-10 = 1e-11 ≫ 2.7e-12). That's an honest, useful note for the editor. Hmm, I could also mention the library's own 6× criterion. OK, keep as editorial.

Hmm, but hold on — I derived r's sign from the assumption that valCpRiddersActual corresponds to the same call. Both valCpRiddersActual = 2.0e-12 (DS1b) and valCheckRiddersEstimate = 2.7e-12 (= valCpRiddersEstimate) are for `derivative<Ridders>(H, T0)` ✓ so they're the same call ✓. Good, that's solid evidence.

OK, writing the report now. I'll aim for dense but complete.

Let me also mention explicitly the items from the task's checklist where I found nothing wrong, so the other reviewer knows what was covered.

One more quick check I should do: the task mentions "the claim about generic code and the default-step overload" (ch8 line 97-98: "use the overload without a step in generic code that varies N (tab:default-steps)" ✓ matches the header's "@ref ... generic code that varies N should use this overload rather than pass one h for every order" ✓ — the header says it at line 1417 ✓). ✓ No finding.

And "the claim that N is 1..4 and that a size_t/int/order is accepted" ✓.

Let me write it.
I have completed my review. Here are my findings.

## Verification performed (so the next reviewer knows what was covered)

Ground truth read in full: `include/numerixx/deriv.hpp` (all 1477 lines, incl. `detail`), `include/numerixx/math.hpp`, `tests/test_callables.cpp`, `test_math.cpp`, `test_deriv.cpp`, `precision_warning.cpp`, `docs/deriv-book/CLAUDE.md`, `tools/check.py`, `CMakeLists.txt`, `CMakePresets.json`, `.github/workflows/ci.yml`, `refs.bib`, `values.tex`, `newton.csv`, `kink.csv`, `accuracy-matrix.csv`, `table-*.tex`, `code/07-using-*.cpp`, `08-api-usage.cpp`, `colebrook.hpp`, `noise.hpp`, `shomate.hpp`, `getting-started.cmake`, both chapter-7 figures.

Checks I ran (read-only): all ten chapter-7 examples + `08-api-usage.cpp` compile clean under `g++ 16 -std=c++23 -Wall -Wextra` (`-fsyntax-only`); `07-using-multiprecision.cpp` compiles against the vendored Boost 1.92 (standalone); I reproduced `check.py`'s `signature_problems()` logic in Python against the headers → **0 failures**, i.e. every transcribed line of ch. 8's `signature` blocks does occur in the headers; I re-derived the data with Python from the Shomate coefficients and the header's own step formula (matches `values.tex` for `Check*`, `Kink*`, `JacobianBalanceStep*`, `Sensitivity*`, `EpsMultiprecision`, `RootEps*`, `CpDefaultRelError`, `Newton*`).

**Chapter 8's declarations: no mismatches found.** Signature, `M`/`N`/`T`/`F` order and defaults, `requires` clause, `std::type_identity_t<T>` for `h`, `[[nodiscard]]`/`[[deprecated]]`, the `noisy` rule (`precision_warning && static_cast<int>(N)==4 && coarse<T>()`, `β^(digits-1) < 10^7`), `enum class error { non_finite }`, `enum class order : int {first=1…fourth=4}`, `Estimate`'s members/aggregate/`operator T()`/defaulted `==` and its documented traps, `Central<2..10 even>`, `Forward/Backward<1..6>` and both static_assert texts, `Ridders`' cost/caps (`64, 65, 128, 129`), `real`, `math::abs/isfinite`, `root_epsilon` (`static_assert(N>=2)`, ulp-accurate, β⁻ᴺ scaling), "24 significant bits" (largest numerator is 6 222 216, per `deriv.hpp:386`), the "requires-expression cannot see the weights static_assert" claim (verified with `std::float16_t`: the `requires`-expression is satisfied and the real call fails with *"this rule's integer weights are not exact in T"*), `assert(h>0 && isfinite(h))`, `NDEBUG`/`Forward<1>` below x, the half-ulp/`dx == 0` story, the Ridders failure/restart/`untrusted` rule (`abs_error >= |best|/4 && abs_error/64 > floor`), the `-Wconversion`/`-Wfloat-conversion` claims (Clang warns on `if (*r)`, GCC does not; both warn on `int n = *r`), `std::println` not applying the conversion, `r.value()` returning the `Estimate`, `*r == 2.0` comparing the value while `r == 2.0` is false on a failure, designated initialisers + structured bindings, Boost arithmetic operators rejecting `Estimate` while its comparisons accept it, and the C++23 `<cmath>` non-`constexpr` claim (cppreference: `constexpr since C++26`; GCC accepts as an extension, Clang rejects). **Nothing to report there.**

Findings, ordered by severity (none is major):

- **`docs/deriv-book/chapters/07-using.tex:186-187` — minor.** *"evaluates the whole model with one unknown varied, at $m \times n$ calls, \num{\valJacobianBalanceCallEvals} here"* — with m = 2 residuals and n = 2 unknowns, m×n = 4, but the macro is 8, which is the **model**-call count. Counter-evidence: `values.tex:454` `\newcommand{\valJacobianBalanceCallEvals}{8}`; `code/07-using-jacobian.cpp:129` `bool const ok = worst(J_hand) <= 1e-8 && worst(J) <= 1e-8 && hand_calls == 4 && library_calls == 8;` (and the same paragraph calls 4 "model calls" for the hand loop). Fix: "at $2mn$ model calls, \num{\valJacobianBalanceCallEvals} here".

- **`07-using.tex:33` and `08-api.tex:78` — minor.** *"or a member function of T itself"* / *"so does a pointer to a member function of \code{T} itself"* — only a **const** member function satisfies `detail::callable`. Counter-evidence (compiled): `static_assert(nxx::deriv::detail::callable<R (R::*)() const, R>);` passes and `static_assert(!nxx::deriv::detail::callable<R (R::*)(), R>);` passes; `std::invocable<R (R::*)(), R const&>` is false, because `fn` is invoked with a `T const` lvalue. The library's own test uses a const member: `tests/test_callables.cpp:38` `[[nodiscard]] constexpr Real square() const { return Real(v * v); }` (and line 29 says so). Fix: "a pointer to a **const** member function of T".

- **`07-using.tex:363-366` — minor (figure caption).** *"the error is a plateau of about $J / (2 h \cp)$, $J$ the jump in $\enth$: \num{\valKinkInsideCentralRelError} for \code{Central<2>} ($h = \SI{\valKinkCentralStep}{K}$) and \num{\valKinkInsideCentralFourRelError} for \code{Central<4>} ($h = \SI{\valKinkCentralFourStep}{K}$), so it falls as $1/h$ across the rules, whatever their order."* Counter-evidence (computed from `code/shomate.hpp`'s coefficients and the default steps): J = 0.7684 J/mol, c_p = 29.5814 J/(mol K), h₂ = 0.0030277 K, h₄ = 0.370048 K → J/(2h₂c_p) = 4.2896 (matches 4.3) but J/(2h₄c_p) = 0.0351, **not** the quoted 0.041; the Central<4> plateau is 7J/(12h₄c_p) = 0.04095, because two of its four nodes (`x−h`, `x−2h`) lie below the break. The two plateaus then stand in ratio 105 while the steps stand in ratio 122, so it does not fall exactly as 1/h "whatever their order". Fix: give the constant per rule, or say "about J/(2h c_p) for the two-node rules".

- **`08-api.tex:270` — minor.** *"\code{derivative<Ridders, double>(fn, x)}, does not compile either, because that argument is the order; \code{T} is always deduced from \code{x}."* Counter-evidence (compiled): `T` is the third template parameter and can be given explicitly — `static_assert(nxx::deriv::derivative<nxx::deriv::Central<2>, 1, float>(sq, 1.0).has_value());` passes with a `double` x converted to `float`. Fix: "T comes from x (never from h); give it third only if you mean to convert".

- **`07-using.tex:286-287` — minor.** *"And give a lambda an explicit \code{-> mp50}: it matters for a type with expression templates, such as \code{number<cpp_bin_float<50>, et_on>}"* — no failure is demonstrated: `requires { nxx::deriv::derivative([](mpet const& x){ return exp(x); }, mpet(1)); }` is satisfied for `et_on` (compiled), and the library casts every call immediately, so no expression object outlives the call: `deriv.hpp:804-805` `T const f_lo = static_cast<T>(std::invoke(fn, x_lo));`, likewise 865-866, 963, 1031, 1200. Fix: keep it as the "safe habit" the same sentence already recommends for `cpp_bin_float_50`.

- **`07-using.tex:71-72` — editorial.** *"Fourth, \code{d.value_or(0.0)} does not compile, since no \code{double} converts to an \code{Estimate<double>}; \code{d->abs_error.value_or(0.0)} does."* — the listing this paragraph annotates binds the result to `r`, not `d` (`code/07-using-expected.cpp:88` `auto const r = derivative<Ridders>(cosine, 0.0);`, snippet `<estimate-pitfalls>`), and the two preceding sentences of the same paragraph use `r` (`if (*r)`, `r.value()`). Fix: use `r`. (The claim itself is right: `value_or` needs a `double`→`Estimate<double>` conversion.)

- **`07-using.tex:260-261` — editorial.** *"the difference $\abs{d - r} = \SI{\valCheckRiddersDiff}{\joule\per\mole\per\kelvin}$ is the stencil's error to within it; the actual error is \SI{\valCheckStencilActual}{\joule\per\mole\per\kelvin}$."* The claim is consistent with the data (my computation gives |d−exact| = 7.7506e-10 and, from `\valCpRiddersActual` = 2.0e-12, |d−r| = 7.7306e-10, margin 2.0e-12 ≤ the 2.7e-12 estimate), but the two macros are printed to two significant digits (7.7e-10, 7.8e-10), so a reader subtracting them gets 1e-11 ≈ 3.7 × the estimate — the claim looks false as printed. The library's own criterion is looser: `code/07-using-check.cpp:99` *"the triangle inequality makes |d - r| the stencil
error: run timed out after 1650s


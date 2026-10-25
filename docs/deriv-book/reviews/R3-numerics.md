# R3-numerics: targeted review of the Phase 5 update (round 3)

Reviewer: numerics (Opus), 2026-10-25, on the uncommitted update over `56d27a1` (`git diff 56d27a1 -- docs/deriv-book`)
and the 90-page build `build/gcc/docs/deriv-book/pdf/main.pdf`.

Scope: the changed passages only. 6.7 and F19 (`figures/06-ridders-onesided-reach.tex`,
`data/forward-ridders-reach.csv`); B.3's new paragraph; the restart rule and the pole sentence in 6.4, Sidebars 6.2
and 6.3; 5.4 (and the sidebar paragraph on p. 30) and Figure 5.2; 10.1, 10.2 and 10.4; the method rows of 7.1, 8.2,
8.5, Table 8.2, Table 9.10, Table 10.1, Table 10.3, T0 and the quick reference. Ground truth: `deriv.hpp` (lines
262-330, 1370-1585), `docs/PHASE5_API.md`, `tests/test_onesided.cpp`, `tests/test_ridders.cpp`, PLAN.md 12.8.

Recomputed (scratchpad `r3num/`, GCC 16, `-O2`, against `include/` as it stands):
- `b3.py` (exact `fractions`): the coefficients of the tableau with column factor 2^j (and 4^j) to j = 12, their
  signs and sums, the floor recursion against sum |gamma_k| r_k on 200 random floor vectors, and the cancellation of
  h^1 ... h^j by column j;
- `reach.cpp`: `ForwardRidders` N = 1 to 4 on log, sqrt, cbrt and 1/t at x = c 2^-k (c = 1, 1.25, 1.5, 1.75;
  k = 2 to 46), and `Ridders` N = 1 to 4 on 1/(t - a) at 0 for a = +-c 2^-k (k = 4 to 44): 4192 calls;
- `edge.cpp`: `Ridders` N = 1 to 4 on log at x = 2^-k, k = 28 to 36 in half steps;
- `fpole.cpp`, `fpole2.cpp`, `fpole3.cpp`, `one.cpp`: `ForwardRidders` and `Ridders` at poles ahead of and behind
  x = 0, a = 10^U(-9, -1), 400 and 4000 draws per order, and a dense grid a = 2^(l/16).

The book's own readings reproduce exactly: k = 37 for log f' (0.224 off, estimate 0.193, not flagged by 0.5 %),
9.38e-9 and 0.00142 at k = 30, and the reach exponents of DS23.

## Findings

### R3-numerics-01 (major). One-sided f'''' with a pole ahead of x: wrong by a factor 60, neither covered nor flagged

Location: `06-ridders.tex` 6.7, "The side has a price" paragraph; `10-guarantees.tex` 10.2, the `ForwardRidders`
paragraph, and 10.4 "Domain edges".

Problem. The book's one-sided coverage (99.9 to 99.96 %, worst 1.7) comes from the nine smooth functions of the
study, which have no feature within the forward span (x >= 0.2, nearest feature at 0 or at +-i/5); 6.7 says so
("finite above every x"). It says nothing about a pole *ahead* of x, inside the span the forward rule samples, which
is the natural case at a correlation's end (a singularity just beyond the range). `PHASE5_API.md` lists it under
"Caveats to document": "2 of 400 poles ahead of x and inside the initial span are uncovered for N = 4, by up to 4.9
times". The book does not carry it, and the failure is worse than "uncovered": the result is not flagged either.

Evidence (`fpole3.cpp`): `ForwardRidders` on 1/(t - a) at x = 0, a = 10^U, U uniform on (-9, -1), 4000 draws per
order:

| N | uncovered | worst error / estimate | neither covered nor flagged |
|---|---|---|---|
| 1, 2, 3 | 0 | - | 0 |
| 4 | 10 | 4.86 | 9 |

The nine silent f'''' cases have a = 0.01756, 0.001097, 0.001098, 6.845e-5, 4.282e-6, 2.679e-7, 1.678e-8,
1.681e-8, 1.045e-9 (successive ratios 16: a resonance of the nodes x + k h_i, k <= 4, with the pole), relative error
58 to 66 (the value is about 60 times the derivative), relative estimate about 12 to 15, so 4 abs_error < |value|.
They are narrow: at a = 1.04e-9 and 1.06e-9 the error is 0.013 and 0.012, covered. Every one lies inside the
documented reach (a >= 2^-30), so the rule of 6.7 ("keep f's nearest feature at least about 2^-30 s from x") does
not protect against it. The same draws with the pole behind x (all N), and central `Ridders` (all N, 4000 draws
each), gave no uncovered case.

Suggested change. One measured sentence in 6.7 after the coverage figures, and the same in 10.2: "With a pole ahead
of x, inside the span the rule samples, f'''' can be wrong by far more than its estimate: on 1/(t - a) at 0, a pole
10^-9 to 10^-1 ahead, 9 of 4000 results were about 60 times the derivative with an estimate 4.3 to 4.9 times below
the error, too small to flag them (measured); f' to f''' and poles behind x had none. Near a known singularity,
check a one-sided f'''' at a second initial step." The worst ratio stays within the factor 6 of the library's tests,
which is worth saying in the same breath.

### R3-numerics-02 (major). "The edge limit holds for a pole too": not shown by the cases quoted, and one halving optimistic for f''' and f''''

Location: `06-ridders.tex` 6.4, "A pole between the samples" paragraph ("so the edge limit above holds for a pole
too: ..."); `10-guarantees.tex` 10.4 "Domain edges" ("The limit holds for a pole ... (measured)"); the quick
reference ("Near a bad point (a NaN or a pole)"); 7.8 ("a NaN or a pole").

Problem (a), the evidence. The limit is 2^-34 s for f' and 2^-31 s for f'' to f''''. The quoted cases put the pole
at 10^-3 = 2^-10 and 10^-4 = 2^-13.3, twenty octaves inside it: they show that the fast restart handles a pole, not
where it stops handling one. In 10.4 the sentence is tagged "(measured)", while the claim it supports is derived from
the restart rule (PHASE5_API.md item C: "the limit then holds for finite poles too", with no measurement near the
limit).

Problem (b), the f''' and f'''' limit, for NaN and pole alike. `Central<2>` for N = 3, 4 reaches x +- 2h. From
h0 = s/4 the 32 levels end at h = 2^-33 s; two finite levels need 2 * 2^-32 s < d, so d > 2^-31 s, and the last
dyadic distance that works is 2^-30 s, not 2^-31 s. The f'' rule reaches +-h, and its limit 2^-31 s is right.
Measured (`edge.cpp`, `reach.cpp`):

| `Ridders`, log x at x = 2^-k | k = 30 | k = 30.5 | k = 31 | k = 32 |
|---|---|---|---|---|
| f'' | 8.2e-5 off | 9.5e-4 | 6.8e-3, estimate 0.157 | non_finite |
| f''' | 0.012 off, estimate 0.24 | 0.066, estimate 0.68 | non_finite (35 evaluations) | non_finite |
| f'''' | 0.014 off, estimate 0.27 | 0.078, estimate 0.77 | non_finite (36 evaluations) | non_finite |

| `Ridders`, 1/(t - a) at 0 | largest uncovered abs(a) | every uncovered one flagged? |
|---|---|---|
| f' | 3.64e-11 = 1.25 * 2^-35 | yes |
| f'' | 2.91e-10 = 1.25 * 2^-33 | yes |
| f''', f'''' | 5.82e-10 = 1.25 * 2^-31 (error 3.4, estimate 2.57); at 2^-31, error 1.0, estimate 0.0015 | yes |

So for f' and f'' the pole limit agrees with the NaN limit; for f''' and f'''' both are 2^-30 s, and at the stated
2^-31 s a pole gives a result 100 % to 340 % off.

Problem (c), what "holds" means. Beyond the NaN limit the result is `non_finite`, loud. Beyond the pole limit it is
a finite, wrong value; in all 1312 pole calls of `reach.cpp` every uncovered result was flagged (4 abs_error >=
|value|), so the reader is protected only by testing `abs_error`.

Suggested change. 6.4: "That restart descends as the one after a non-finite value does, so a pole is handled down to
about the same distance (derived); on the tests' poles, 10^-3 ahead and 10^-4 behind, 8 of 8 results are covered ...
(measured). Closer than that, Ridders returns a wrong value rather than `non_finite`, flagged by an estimate of at
least a quarter of it." In 6.4 and 10.4 give the higher-order limit as 2^-31 s for f'' and 2^-30 s for f''' and
f'''' (their reach is 2h), or measure it in the generator as for f'. The f''' and f'''' figure predates the update,
but the update now extends it to poles.

### R3-numerics-03 (minor). The silent band beyond the one-sided reach is wider than "one distance"

Location: F19 caption ("At one distance for f' (the circle) it is neither covered nor flagged"); 6.7 ("At k = 37
between them ..."); 10.4.

Problem. The single silent case is a property of log at dyadic x. On other functions with an edge at 0 the gap
between "covered" and "flagged" is a band about two halvings wide. `reach.cpp`, `ForwardRidders` f' at x = c 2^-k:

| function | x where neither covered nor flagged | relative error / relative estimate |
|---|---|---|
| log | 2^-37 only | 0.224 / 0.193 |
| sqrt | 1.8e-12 to 5.5e-12 (7 points, 2^-39 to 2^-37.4) | 0.14 to 0.33 / 0.13 to 0.15 |
| cbrt | 3.6e-12 to 6.4e-12 (4 points) | 0.16 to 0.27 / 0.16 to 0.17 |
| 1/t | none | - |

N = 2 to 4: none on any of the four. The worst case, sqrt at 2^-39, is 2.2 times its estimate, within the tests'
factor 6, and the rule "at least about 2^-30 s" keeps clear of every one by seven halvings.

Suggested change. Keep the log numbers, and add in 6.7 (and 10.4): "between the two lies a narrow band, one point for
log and about two halvings for sqrt and the cube root, where f' is 15 % to 33 % off with an estimate of 13 % to 19 %:
neither covered nor flagged, but within about twice the estimate (measured)". The caption's "At one distance"
becomes "At one distance for log".

### R3-numerics-04 (minor). F19 plots the estimate in the wrong frame for the "flagged" claim

Location: F19 (`figures/06-ridders-onesided-reach.tex`) and its caption; 6.7.

Problem. The dashed curves are abs_error / |exact|; the flag is abs_error / |value| >= 1/4. Beyond the reach the
value tends to 0, so in the plotted frame the f'''' estimate falls to 3.5e-11 at k = 45 while the error is 1: the
figure shows an estimate wrong by ten orders of magnitude, and the caption tells the reader the opposite
("flagged"). The reader can only compute abs_error / |value|; the figure's source comment acknowledges the mismatch.

Suggested change. Plot abs_error / |value| for the dashed curves (indistinguishable from the present ones wherever
the result is good, since value ~ exact there) and draw the flag threshold 1/4 as a thin horizontal line, so that
"flagged" is visible as dashed above the line. If the frame stays, say in the caption: "estimates relative to the
exact value; relative to the returned value, which tends to 0, they exceed 1/4 from k = 38 (f') and 33 (f'''')".

### R3-numerics-05 (minor). B.3: the alternation makes the recursion exact; it is not what makes the bound hold

Location: `a2-derivations.tex`, the new paragraph: "whose coefficients alternate in sign, so the rounding it inherits
is at most sum |gamma_k| r_k, which is what the recursion ... computes".

Problem. Both statements are right (verified exactly: for the factors 2^j and 4^j to j = 12 the coefficients
gamma_k have the sign (-1)^k and sum to 1; the floor recursion equals sum |gamma_k| r_k exactly on 200 random r;
column j cancels h^1 to h^j and leaves h^(j+1)), but the "so" is misplaced. |sum gamma_k e_k| <= sum |gamma_k| r_k
for |e_k| <= r_k is the triangle inequality, whatever the signs. The alternation is needed for the next step: the
recursion adds the parents' floors with the positive weights 2^j / (2^j - 1) and 1 / (2^j - 1), while the
coefficients combine as 2^j gamma_k - gamma'_{k-1}; the two agree in magnitude only because those terms have the
same sign. The bound is also first order: it takes r_k as a bound on the base rounding (r_k is the model floor) and
ignores the rounding of the extrapolation itself.

Suggested change. "... a combination sum_k gamma_k D_k of the base values, so the rounding it inherits is at most
sum_k |gamma_k| r_k, to first order. The coefficients alternate in sign, so 2^j gamma_k and -gamma'_{k-1} never
cancel, and the recursion r_{i,j} = (2^j r_{i,j-1} + r_{i-1,j-1}) / (2^j - 1) computes that sum exactly."
Optional, derived, and useful for 6.7's "the smaller divisors amplify rounding more": the newest level's coefficient
tends to prod 1/(1 - 2^-m) = 3.46 against prod 1/(1 - 4^-m) = 1.45 centrally (the header's "about 3.5 times rather
than 1.4"), and sum |gamma_k| to 8.25 against 1.97.

### R3-numerics-06 (minor). "s/32 for every order, since the rules reach N h"

Location: 6.7, second paragraph.

Problem. The reach N h does not imply one divisor for every N (if anything it argues for a smaller step as N grows).
The divisor was chosen by a sweep: PHASE5_API.md, "s/32 and s/64 are the starts with no gross errors" for N = 1 to
3, and s/32 has "the fewest gross errors" for N = 4 (5, against 9 at s/16 and 34 at s/64).

Suggested change. "The default initial step is s/32 for every order, chosen on a prototype sweep (the design record):
the start with the fewest gross errors for every N."

### R3-numerics-07 (editorial). "restarts ... until the levels run out"

Location: 6.7, last paragraph: "a one-sided base sees no non-finite value, only untrusted results, and restarts at
half the smallest step it tried until the levels run out".

Problem. Within the reach the restarts end in a converged tableau (log f' uses its whole budget from k = 25 and is
still accurate to 6.8e-13 at k = 28); the levels run out only beyond it. As written, every near-edge call exhausts
the budget.

Suggested change. "... and restarts at half the smallest step it tried; closer than its reach, the levels run out
first."

### R3-numerics-08 (editorial). Two untagged statements

Location: 10.1, "`ForwardRidders`, whose columns gain one order rather than two, loses more with each order: ... 14,
212, 1261 and 3007 times Ridders'"; 6.7, "The smaller divisors amplify rounding more".

Problem. Chapter 10 promises a tag on every statement; this one has none. The factors are measured (the median of
the one-sided errors over the median of Ridders' on the same 2700 cases, `generate_data.cpp` line 1351); the
mechanism is derived (R3-numerics-05).

Suggested change. "(measured)" after the factors in 10.1; optionally the derived amplification of -05 in 6.7.

### R3-numerics-09 (editorial). Figure 5.2's box touches a double point and cuts the boxed markers

Location: `figures/05-orders-order-accuracy.tex`, the accent rectangle (1.85, 1.5e-5) to (4.06, 3.5e-2); PDF p. 30.

Problem. The lower edge at 1.5e-5 meets the marker of `Central<2>` double f'''' (1.2e-5, a factor 1.25 below), so
that point reads as boxed; the right edge at N = 4.06 runs through the N = 4 float markers. The lowest float point to
enclose is `Central<4>` f'' at 2.0e-5, so a rectangle has a factor 1.3 of room at N = 2 and none at N = 4.

Suggested change. Colour the float markers for N >= 2 in the accent instead of boxing them (the caption then reads
"in the accent"), or end the box at N = 4.2 with a sloped lower edge (from 1.6e-5 at N = 2 to 2.5e-4 at N = 4,
above the double point).

### R3-numerics-10 (editorial). Sidebar 6.3: "relative error 0"

Location: Sidebar 6.3, "(measured: relative error 0, estimate 2.2e-16, 20 evaluations in all)".

Problem. The reference is 64 pi rounded to double; 0 means the result is that double, not 64 pi, whose relative
distance from it is up to 2^-53 = 1.1e-16.

Suggested change. "(measured: the double nearest 64 pi, estimate 2.2e-16 relative, 20 evaluations in all)".

## Checked and correct (round 3)

- **B.3, the divisor 2^j - 1.** With A_{i,j-1} = f^(N) + c h_i^j + O(h_i^(j+1)) and the parent at 2 h_i,
  (2^j A_{i,j-1} - A_{i-1,j-1}) / (2^j - 1) cancels h^j: verified exactly to j = 6 on D(h) = sum b_m h^m with
  independent b_m (first surviving power j + 1, constant term 1). The `Forward<1>` rule of order N (nodes 0 to N,
  binomial weights, largest 1, 2, 3, 6) has an expansion in every power of h for generic f.
- **B.3 and 6.7, the carried floor.** The coefficients alternate and sum to 1 for both factors; the recursion equals
  sum |gamma_k| r_k exactly; the code (`propagated_floor`, `deriv.hpp` 1387; factor 2^j at column j, 1501; the
  entry's carried floor in the estimate and the untrusted test, 1528 and 1572) is what the book describes.
- **The restart rule.** `h_start = h_i / 2` after a non-finite level (1516) and after an untrusted result (1581),
  h_i being the smallest step tried: "at half the smallest step it tried" (6.4, 6.7, Sidebars 6.2 and 6.3, 8.2,
  10.2, Table 9.10) and "at half that step" (8.2, non-finite) are right. Table 10.2's caption is right that a restart
  after an untrusted result continues the halving and cannot be told apart in the call log.
- **6.7 and F19, against the CSV.** Ten digits for f' to k = 29 (4.0e-11; 9.4e-9 at 30); covered unbroken from k = 2
  to 36, 34, 33, 32; flagged from 38, 35, 34, 33, contiguous to 45 for every N; one silent point for log, f' at 37;
  central f' more accurate than one-sided f' at every k where it returns a value, non-finite from 35; the axis labels
  (2^-10 ~ 10^-3 ... 2^-40 ~ 10^-12); the cap line at k = 25, with f' accurate for the four halvings after it.
- **Costs.** 1 + 32 N = 33, 65, 97, 129, reached (DS23, k >= 26 for N = 3, 4); medians 10, 17, 22, 29; Table 8.2's
  nodes x + k h_i, k = 0 ... N; f(x) once for every N; `BackwardRidders` as the mirror with the sign of odd N.
- **Study figures.** 99.96, 99.96, 99.96, 99.9 % (1, 1, 1, 3 of 2700), worst ratios 1.3, 1.0, 1.6, 1.7, medians
  1.5e-13 to 3.5e-6, factors 14 to 3007 as ratios of medians on the same cases; the scope "finite above every x" is
  true of the study (x in [0.2, 3.0]).
- **E24.** Initial step 500/32 = 15.6 K, "about 16 K"; Ridders non_finite after 64 evaluations.
- **The float warning.** eps > 10^-7 if and only if beta^(digits - 1) < 10^7: float, decimal32, bfloat16 in, double
  out; the quoted message matches `deriv.hpp`; 5.4's medians (5.0e-6, 1.4e-4, 4.0e-4, 0.014 for `Central<2>`;
  2.7e-5 to 1.9e-3 for Ridders, estimates 5.3e-4 and 0.28; 0.0022 to 5.7 one-sided, "more than the value itself")
  match `values.tex` and Table 10.1, where every float cell above f' carries the asterisk.
- **The pole cases quoted.** 8 of 8 covered, 1.6e-15 to 3.8e-9, 32 to 88 evaluations (generator, the tests' two
  poles); my runs agree in kind (central Ridders covered for every abs(a) >= 7e-10, every N).
- **Tables 9.10 and 8.2, T0, the quick reference, 7.1:** every one-sided constant and count agrees with the header.

## Verdict

The update is sound in its mathematics: the 2^j - 1 divisor and the carried floor are derived correctly (B.3 needs
one "so" moved), the restart rule is described as the code implements it, and every reach and coverage number
reproduces from the library. Two things should change before release. The design record's own caveat for the
one-sided methods, a pole ahead of x for f'''', is missing, and it is a silent factor-60 error well inside the
documented reach (01). The new claim that the edge limit holds for poles rests on cases twenty octaves inside the
limit, and the f''' and f'''' limit it inherits is one halving optimistic, for NaN and pole alike (02). The rest is a
wider silent band than F19 suggests (03), a figure frame that hides the flag it describes (04), two reasons to
re-tag (05, 06) and wording.

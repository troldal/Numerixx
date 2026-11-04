# Milestone 3, Phase 8: the quadrature chapter's figures

Status: approved (2026-11-03, figures 1 to 5 as proposed and a figure 0 at the start, at the user's request), in
progress.

Six figures in the quadrature chapter, each because it makes a concept or the working of a method easier to see than
the text alone does; the test figure of Phase 7 goes. Each is made as `docs/CLAUDE.md` describes: a program built and
run as a test (`examples/<name>.cpp`) that holds the figure's table in a raw string and asserts every claim the text,
the caption and the alt text make; a script (`plots/<name>.py`) that draws only that table; the directive on the
page. Every table is computed with `+ - * /` and `std::sqrt` only, or in Boost's software arithmetic, so that it is the
same on every toolchain.

## The figures

0. **A weighted sum** (`index.md`, "A weighted sum, and an error estimate", after the rule's formula):
   `quadrature_weighted_sum`. $4/(1 + x^2)$ on $[0, 1]$, the overview's own example, with the area under it shaded,
   and `GaussLegendre<5>` as five bars side by side from 0, the $i$-th of width $w_i$ and height $f(x_i)$, its node
   marked on the curve: the rule replaces the area by bars whose total area is almost the same. The nodes are recorded
   through `fn`; each weight is the rule applied to the function that is 1 at its node and 0 elsewhere (public API
   only). The program asserts that each bar contains its node (the separation theorem of Gauss rules) and that the
   bars' total is the rule's result.
1. **Where the default method cuts** (`index.md`, after the paragraph on adaptive Gauss-Kronrod):
   `quadrature_panels`. The peak of "Set tolerances", $1/((x - 0.3)^2 + 10^{-4})$ on $[0, 1]$, above two rows of the
   final panels' ends, at the default request (6 bisections, 16 panels, the narrowest 1/80 wide) and at $10^{-12}$
   (16, 26, 1/160); away from the peak the panels keep the initial width 0.1. The panels are rebuilt from the points
   `fn` receives: the first 150 are the 10 initial panels, each further 30 the two halves of a bisected panel, the one
   that contains them all, which the program checks.
2. **Why `TanhSinh` copes with a singularity at an end** (`rough-integrands.md`, "A singularity at an end"):
   `quadrature_tanh_sinh_map`. Left, the map $x = (1 + \tanh(\sinh t))/2$ onto $[0, 1]$, the library's own ($c = 1$),
   with points at equal steps in $t$ crowding towards both ends; right, on a logarithmic scale, the integrand in $t$
   for $f = 1/\sqrt x$, finite and falling off double-exponentially on both sides, so that equal steps in $t$ converge
   fast. `exp` is computed in the program from `+ - * /` (a series after halving), so that no libm's last bits reach
   the table; the program checks the map against the points `TanhSinh` evaluates.
3. **The rounding floor** (`tolerances.md`, "More digits"): `quadrature_floor`. The same peak at requests from
   $10^{-2}$ to $10^{-16}$: the relative estimate below the request until the floor, then at the floor whatever the
   request (690 evaluations from $10^{-13}$ on), and the actual error, from an independent value of
   $100(\arctan 70 + \arctan 30)$, two to three orders of magnitude below the estimate. The page's table stays.
4. **Evaluations against digits** (`high-precision.md`, "Ask for the digits you want"): `quadrature_digits`.
   $1/(1 + x^2)$ on $[0, 1]$ in `cpp_bin_float_50`, requests of 5 to 48 digits: the default at 150 evaluations up to
   20 digits, then 390 to 18,990; `TanhSinh` at 172 up to 25 and 325 up to 48. The text gains that below about 20
   digits the default costs no more.
5. **Where the infinite-interval maps look** (`infinite-intervals.md`, "Where the maps look"): `quadrature_line`.
   `TanhSinh`'s 133 points on the whole line (from $1/(1 + x^2)$), with the peaks at $x_0 = 10$ and $x_0 = 40$
   marked: about 2 apart near 10, about 11 near 40, wider than the peak.

The Phase 7 figure (`quadrature_sampling`, in "Measure the cost on your own function") is deleted with figure 2, which
explains the same thing through its cause; its counts stay in the text and in Costs. `docs/CLAUDE.md` then names
figure 1 as the example to follow.

## Steps

One commit per figure, in the order 1, 0, 2, 3, 4, 5, each after:

- the program on the five toolchains (`ctest -R <name>`; WSL from a copy);
- `numerixx_docs` and `numerixx_docs_pdf` with no warning, and the figure read on the site and in the PDF, on the
  light and the dark surface;
- `docs-reviewer` on the page, the program and the script, its adopted findings recorded here.

After the last, CI and the independent review (Codex), as for every phase.

## As built, and the reviews

**Figure 1.** The narrowest panels are four of 1/80 at the default request and eight of 1/160 at $10^{-12}$, together
covering $[0.275, 0.325]$ (the probe had suggested two). `docs-reviewer`: nothing major; every number confirmed by a
reconstruction of its own (each block's centre and outer nodes), and the table unchanged at -O0 to -O3, with
`-mfma -ffp-contract=fast` and `-march=native`. All ten findings adopted: "where $f$ changes fastest" (its slope is 0
at the peak) became "around the peak"; "whatever the request" asserted through the longest run, at a request of 0 (690
evaluations, 18 cuts), of which every run is a prefix; the caption says what a tick marks; "which asks for a relative
error"; the figure's paragraph split from the program's; "a relative error of $10^{-12}$" throughout and
`rel_tol = 1e-12` in the label, from the table; shorter sentences; one return code per claim, and "each cut a panel that
ends at 0.3" asserted; the default request as `Criteria<double>{}.rel_tol`.

**Figure 0.** The weights come from the rule applied to an indicator of each node, which the reviewer showed to be
exact: on $[0, 1]$ every point is $c$ or $c \pm d$, bit-identical between calls, and a single nonzero term is summed
without rounding. The rule's error, $-1.37 \times 10^{-8}$, from numpy and from the closed-form nodes in 50 digits.
`docs-reviewer`: nothing major, all eleven findings adopted: the error asserted to two digits, in
$[-1.4, -1.35] \times 10^{-8}$; the figure 2.6 inches high, so that it fits on the chapter's first page of the PDF rather than leaving 40% of it
blank; the shaded area at alpha 0.4, visible on the dark site; the paragraph no longer repeats the caption, says why a
bar containing its node matters (a Riemann sum with uneven bars), says what the sum gains and misses, and bridges to the
error estimate ("Here the integral is known; in general it is not"); the alt text's dot "on the bar's top edge"; the
theorem named (Chebyshev-Markov-Stieltjes, Szego 3.41.1); the legend in the text's order, $f(x)$ throughout; the
program's names `integrand` and `Rule`.

## Success criteria

- [x] Figure 1, where the default method cuts (five toolchains and WSL GCC locally; `docs-reviewer`; df9851c: run
  41722681387)
- [x] Figure 0, a weighted sum (five toolchains and WSL GCC locally; `docs-reviewer`)
- [ ] Figure 2, the tanh-sinh map; the Phase 7 figure deleted
- [ ] Figure 3, the rounding floor
- [ ] Figure 4, evaluations against digits
- [ ] Figure 5, where the infinite-interval maps look
- [ ] CI passes; the independent review and its adopted findings

# Language pass, stage B: methods author (Chapters 4 to 6, Appendix B)

PLAN.md section 12.9. Files rewritten: `chapters/04-stencils.tex`, `chapters/05-orders.tex`,
`chapters/06-ridders.tex`, `chapters/a2-derivations.tex`. The figure files `figures/04-*`, `05-*` and `06-*` hold no
captions (every caption is in the chapters and was rewritten there); their plot labels are unchanged. The code in
`code/`, Chapters 1 and 2 and every other author's files are untouched. Nothing committed.

## Checks

| File | `check_language_pass.py 3192e8a` | `check_chapter.sh` | Pages (standalone, before -> after) |
|---|---|---|---|
| `04-stencils.tex` | RESULT: identical content sequences (headings 6, environments 24, displays 7, `\val` 37, snippets 2, cite keys 3, labels 17, cref targets 42, inputs 4) | COMPILED OK, 0 overfull boxes | 6 -> 7 |
| `05-orders.tex` | RESULT: identical content sequences (headings 7, environments 16, displays 5, `\val` 63, snippets 2, cite keys 2, labels 14, cref targets 39, inputs 2) | COMPILED OK, 0 overfull boxes | 5 -> 6 |
| `06-ridders.tex` | RESULT: identical content sequences (headings 12, environments 26, displays 5, `\val` 143, snippets 3, cite keys 7, labels 21, cref targets 47, inputs 7) | COMPILED OK, 0 overfull boxes | 10 -> 11 |
| `a2-derivations.tex` | RESULT: identical content sequences (headings 5, environments 14, displays 11, `\val` 10, cite keys 7, labels 12, cref targets 32) | COMPILED OK, 0 overfull boxes | 4 -> 5 |

`python tools/check.py` (no log): 0 errors, 0 warnings. ASCII scan of the four files: clean. The full book was not
built, as instructed. The standalone PDFs of Chapters 4 to 6 were inspected page by page for margin-float collisions;
none.

### Review lines, each read

Every added inline formula is a step the text already implies; no number, claim or `\code` name was dropped.

Chapter 4
- Margin notes 1 -> 1 (the existing node/sample/stencil/rule note, split into sentences). Glosses in the text:
  moments and moment conditions (emphasised at definition), amplification, decade ("a factor of 100"), V-shaped
  curve, `consteval` (inside Sidebar 4.1, where margin notes are not used), the 24 bits as float's significand.
- Inline maths added: $n_k h$, $h^{\dorder} f^{(\dorder)}$, $j = \dorder$ (why $\mom{\dorder} = \dorder!$ gives the
  derivative term); $n$, $-n$ (why antisymmetry kills the even moments); $12 b = -1$ (solving the two conditions);
  $\sum_k f_k L_k^{(\dorder)}(0)$ (why the derivative of $p$ gives the weights); $\dorder! = 1$, $-4$ (the worked
  weight $-4/6$); $2^{\aorder}$ (why the slope is the order); $E(h)$ (minimising the model); the four-sample
  expansion of Sidebar 4.2 ($4 \Dc{h'}$, $8 (f(x + h') - f(x - h'))$, $4h'$, $\Dc{h} = \Dc{2h'}$,
  $f(x + 2h') - f(x - 2h')$, $12h'$); $(-1)^{\dorder}$, $h$, $-h$, $h^{\dorder}$ (why the mirror image multiplies
  by $(-1)^{\dorder}$); single symbols.
- Code added: `Central<2>` (naming the rule that `\valErrcCentralTwo` belongs to).
- Prose numbers added: 10 ("The library stops at order 10", for "stops there"), 100 (two decades), 24 (float's
  significand), 3 and 6 (the division by 3 in Sidebar 4.2; "$-4$ over 6").

Chapter 5
- Margin notes 1 -> 3. New: "deprecated" (a C++ term the float warning depends on) and "levels" and "tableau" of
  Ridders' method, which Section 5.5 uses one chapter before Chapter 6 defines them. Neither margin note sits beside
  a margin float.
- Inline maths added: $2/24 = 1/12$ (the error constant of the three-point rule); $-1$ (the symmetric pair);
  $\aorder = 2$, $2/(2 + \dorder)$, $\dorder = 1$ (the exponents of `Central<2>`); $2^4$ (the factor sixteen);
  $1, 2, 1$ (why $\ampl = 4$); $\abs{f^{(\aorder + \dorder)}}$, $\abs{f} / \scale^{\aorder + \dorder}$,
  $h^{\aorder + \dorder} = \epsmach \scale^{\aorder + \dorder}$ (how balancing the two terms gives
  \cref{eq:orders-step}); $h_+ = h_- = h$ (the offset form reduces to the three-point rule); single symbols.
- Prose numbers added: 1, 1, 1, 2, 4 ("the factor is about 1", "far below 1", the absolute weights 1, 2, 1).
- Sidebar 5.2 now states Sterbenz's lemma in one sentence before applying it.

Chapter 6
- Margin notes 0 -> 0 (Section 6.4 and 6.7 have margin floats beside them; glosses go in the text). Glosses in the
  text: tableau, two-level stall ("a stop after two levels without progress"), covered (the estimate at least as
  large as the actual error), untrusted, resolution of $x$ ("where the step no longer moves the points"),
  default-constructible.
- Inline maths added: $(h/2)^2 = h^2 / 4$; $D_i + (D_i - D_{i-1}) / 3 = (4 D_i - D_{i-1}) / 3$ and
  $h = h_{i-1}$ (column 1 is Richardson's combination); $\Aij{i}{j-1}$, $\Aij{i-1}{j-1}$, $j - 1$, $h^{2j}$ (why
  the parent's error is $4^j$ times the left neighbour's); $f''' = 2/x^3$, $f' = 1/x$, $(h^2/6) f''' / f'$ (where
  the stencil's $h^2 / 3x^2$ comes from); $c \abs{h} = c \abs{-h}$ (why the central difference of $c \abs{t}$ is
  zero); $(-1)^{\dorder}$, $f(-t)$ (why `BackwardRidders` negates odd orders); $2^j$ (why the one-sided divisor is
  $2^j - 1$); $4^7 \cdot 2^{115} = 2^{129}$, $2^{128}$ (the float overflow of the product form);
  $1/4 = 2^{-2}$ (two correct bits); $64 \pi / 32 = 2\pi$ (the aliasing samples); $\scale = 1$, $h_0 = 1/4$ (why
  $\log x$ and $\sqrt{x}$ restart below 0.25); single symbols.
- Code added: `ForwardRidders`, `BackwardRidders` (repeated as sentence subjects after a split).
- Prose numbers added: 4 ("with the factor $4^2 = 16$ in place of 4").

Appendix B
- Margin notes 0 -> 1: little-o, which B.1 uses without definition (Chapter 2 glosses only big-O). Glosses in the
  text: a posteriori ("the one computed from the results themselves"), the intermediate-value theorem, the triangle
  inequality (named).
- Inline maths added, all steps of the proofs: the moment as the factor of $h^j$; $\mom{\dorder} / \dorder! = 1$;
  the limit of $f^{(\dorder + \aorder)}(\xi_k)$ and $\sum_k w_k n_k^{\dorder + \aorder}$ (the asymptotic remainder);
  the weighted-mean argument for one intermediate point; the product rule at $t = 0$ and $h f'(x) - p'(0)$;
  $4 / 5! = 1/30$; the pair contribution $w_n n^j + w_{-n} (-n)^j = w_n n^j (1 + (-1)^{\dorder + j})$; why
  $L_k$ is 1 and 0 at the nodes, and why two interpolants agree; the coefficient $t^1$; $(2 h_i)^{2j} = 4^j h_i^{2j}$;
  $4^j = (4^j - 1) + 1$ (the line of algebra between the two forms); $\Aij{i}{j} - f^{(\dorder)}$ is
  $\order(h_i^{2j+2})$; convexity of $h^{\aorder}$ and $1/h^{\dorder}$; the ratio $\dorder : \aorder$ from
  multiplying the minimum condition by $h$; ${\hbest}^{\dorder}$ proportional to
  $\epsmach^{\dorder/(\aorder + \dorder)}$; the values of $\aorder$, $\ampl$, $\errc$ that give $3^{1/3}$ and
  $(2 \cdot 2)^{1/2}$.
- Prose numbers added: 0 (six times: "tends to 0", "weight 0", "is 0", "difference is zero") and 1.

## Word count (`wc -w` on the source; the script's prose estimate in brackets)

| File | Before (3192e8a) | After |
|---|---|---|
| `04-stencils.tex` | 3262 (2960) | 3818 (3484) |
| `05-orders.tex` | 2507 (2297) | 2883 (2620) |
| `06-ridders.tex` | 4971 (4762) | 5711 (5455) |
| `a2-derivations.tex` | 2251 (2027) | 2803 (2511) |
| Total | 12991 (12046) | 15215 (14070) |

## Two representative sentences

Before (6.3): "The method keeps the entry with the smallest estimate, and stops when two further levels produce no
smaller one, because each level doubles the floor, so more levels only add rounding (on the nitrogen example, the
last two levels of \cref{tab:ridders-tableau})."

After: "The method keeps the entry with the smallest estimate. It stops when two further levels produce no smaller
one. The reason is that each level doubles the floor, so more levels only add rounding. On the nitrogen example,
these are the last two levels of \cref{tab:ridders-tableau}."

Before (4.1): "A rule that is to be exact for the odd function $f(x + t) = t$ and give zero for even functions has
antisymmetric weights, [...], so every even moment vanishes by itself. The two odd moments that remain give two
conditions, $\mom{1} = 2a + 4b = 1$ and $\mom{3} = 2a + 16b = 0$, whose solution is $a = \tfrac23$ and
$b = -\tfrac{1}{12}$:"

After: "We want a rule that is exact for the odd function $f(x + t) = t$ and gives zero for even functions. Such a
rule has antisymmetric weights: [...]. With these weights every even moment vanishes by itself, because the terms for
the nodes $n$ and $-n$ cancel. Two odd moments remain, and they give two conditions: $\mom{1} = 2a + 4b = 1$ and
$\mom{3} = 2a + 16b = 0$. Subtracting the first condition from the second gives $12 b = -1$. So the solution is
$a = \tfrac23$ and $b = -\tfrac{1}{12}$, and the rule with its error is"

## Compromises and points for the editor

1. **Section 6.7 breathes.** Its four dense paragraphs are now six: the one-sided forms; how the one-sided tableau
   differs (with the floor recursion explained as adding the floors with the entry's weights); the default step and
   cost; the price; the nitrogen edge example; the reach. No `\needspace` or float placement changed in any of the
   four files.
2. **One rewording forced by layout (4.5).** "recorded in the design record of the repository
   (`docs/deriv/PHASE2_API.md`)" left the unbreakable path overfull by 15 pt in the new line breaks. It now reads
   "The design record of the repository (`docs/deriv/PHASE2_API.md`) records them." Same content.
3. **Style decision 19 (Chapter 6).** The margin figures F14b (cost) and F15 (one-sided reach) are still a page or
   more away from the full-width F12 and F14 in the standalone build, but the chapter is a page longer: please
   confirm at the full build, where the parity differs.
4. **A forward gloss in Chapter 5.** Section 5.5 uses "tableau" and "level" before Chapter 6 defines them. The new
   margin note glosses both, naming the next chapter without a `\cref` (none may be added).
5. **Terminology kept as written.** "Parent" means the left neighbour's parent, $\Aij{i-1}{j-1}$, in 6.2, and "the
   entry's parent" in 6.3 and B.3 is the same entry. The pass makes this explicit in 6.3 by naming both entries; it
   does not rename anything.
6. **Signposting** added at the openings of 4.1, 4.2, 4.3, 5.2, 5.3, 6.2, 6.4 (a list of the four situations its
   paragraphs treat), 6.5 and Appendix B (its four questions, one per section). No closing summaries: 4.5 and 6.6
   already close their chapters' arguments, and adding one to Chapter 6 would sit after the sidebars.
7. **Typed numbers kept.** The prototype figures in sidebars (8.7 and 1.4 in 5.2; 3.5, 11, 15, 1.6 and 1.2 to 1.8
   in B.4; 2.3 to 7 in 6.5; 6,222,216 in 4.1) are unchanged, as before.
8. **No disagreement between the sources** was found in the passages touched. The glosses checked against
   `include/numerixx/deriv.hpp`: the 24-bit limit ("every type with at least 24 significant bits", line 150), the
   untrusted test (`abs_error >= |best| / 4` and `abs_error / 64 > best_floor`, line 1559), the two-level stall
   (line 1541), `abs_error` as `std::optional<T>` (line 82), Ridders' initial steps $\scale / 32$ and $\scale / 4$
   (line 267), and the deprecated overloads (line 1796).

No new data sets or references are needed.

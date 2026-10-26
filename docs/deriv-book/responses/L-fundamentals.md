# Language pass, stage B: fundamentals author (Chapter 3, Appendix A)

PLAN.md section 12.9. Files rewritten: `chapters/03-rounding.tex`, `chapters/a1-floating-point.tex`. The figure files
`figures/03-rounding-*.tex` hold no captions (all five captions are in the chapter and were rewritten there); their
plot labels are unchanged. Chapters 1 and 2, the code and every other author's files are untouched. Nothing committed.

## Checks

| File | `check_language_pass.py 3192e8a` | `check_chapter.sh` | Pages (standalone) |
|---|---|---|---|
| `03-rounding.tex` | RESULT: identical content sequences (headings 12, environments 40, displays 9, `\val` 54, snippets 6, cite keys 8, labels 20, cref targets 49, inputs 6) | COMPILED OK, 0 overfull boxes | 11 |
| `a1-floating-point.tex` | RESULT: identical content sequences (headings 7, environments 10, displays 3, `\val` 21, cite keys 4, labels 8, cref targets 9) | COMPILED OK, 0 overfull boxes | 3 |

`python tools/check.py` (no log): 0 errors, 0 warnings. The full book was not built, as instructed.

### Review lines, each read

Chapter 3
- Margin notes 3 -> 3. No new `\mnote`: margin figures stand beside most of Sections 3.5 and 3.6, so new glosses
  go in the text (rule of 12.9). Glosses in the text: decade, ulp, `fl` (rounded value), enthalpy of formation (a
  constant that shifts every value), Arrhenius factor, flash calculation, residual, tear stream, "optimistic"
  (estimate smaller than the actual error), "the 1 is in your units". NaN, relative error, number type and correct
  digits were glossed in Chapter 1.
- Inline maths dropped: `\epsmach \abs{f} / h`, which now ends the spelled-out step
  `2 \epsmach \abs{f} / (2h) = \epsmach \abs{f} / h`.
- Inline maths added. All are spelled-out steps:
  - the division by $2h$ (above);
  - $h^3 = 3\epsmach|f|/|f'''|$ (solving the derivative condition);
  - $|f'''| = 3\epsmach|f|/(\hbest)^3$ (where $E(\hbest) = \frac32 \epsmach|f|/\hbest$ comes from: the truncation
    term is half the rounding term);
  - $\frac32 \cdot 3^{-1/3}$ (the constant that $\approx \epsmach^{2/3}$ drops; it is 1.04);
  - $x + h \le 2(x - h)$ and $3h \le x$ (why Sterbenz's condition means $h \le |x|/3$);
  - $\aorder = 2$, $\dorder = 1$, twice (giving back $\epsmach^{1/3}\scale$, and $7 \cdot 3^{1/3}/9$ from the cost
    formula);
  - $\epsmach^{1/3}\scale$ and $\mathrm{e}^x$, which are repeated references;
  - $1/h$ (why a step 15 times too small costs about one digit);
  - $\noise u$ (why the noise estimate is near $\noise/\sqrt3$: the root-mean-square of $\noise u$ for $u$ uniform on
    $[-1, 1]$);
  - $T$ (the Arrhenius scale $\Rgas T^2/E_a$ is the change in $T$ over which the exponent changes by about 1);
  - $h^2$ (truncation falls as $h^2$, hence the factor 100 in the hand table);
  - single symbols $f$, $f'$, $h$, $k$, $x$.
- Code added: `double` (in "So a `double` cannot hold every real number").
- Prose numbers added: 1, 1, 1, 10 ("a factor of 10" in the step, "close to 1", "changes by about 1", "1 in
  whatever unit").

Appendix A
- Margin notes 0 -> 0.
- Inline maths added:
  - $2^{1-p}$ and $p$ (why $\epsmach = 2^{1-p}$: the worth of the last significand digit just above 1);
  - $1 \le m < 2$ (the definition of a normal number, from the display in A.1);
  - $\circ$ and $\fl$ (naming the symbols of the display);
  - $x + h$ (signposting).
- Code added: `double`.
- Prose numbers added: 0, 1, 2 ("exponent is 0", "the next larger number than 1", "base 2").

## Word count (non-comment source words, `wc -w`; the script's prose estimate in brackets)

| File | Before (3192e8a) | After |
|---|---|---|
| `03-rounding.tex` | 4233 (4002) | 5046 (4750) |
| `a1-floating-point.tex` | 962 (907) | 1254 (1190) |

## Two representative sentences

Before (3.3): "Setting the derivative of \cref{eq:error-model} with respect to $h$ to zero,
$h \abs{f'''} / 3 - \epsmach \abs{f} / h^2 = 0$, gives the \emph{best step}"

After: "The bottom of the V is where the model is smallest, and there its derivative with respect to $h$ is zero.
Differentiating \cref{eq:error-model} with respect to $h$ and setting the result to zero gives
$h \abs{f'''} / 3 - \epsmach \abs{f} / h^2 = 0$, so $h^3 = 3 \epsmach \abs{f} / \abs{f'''}$. This is the
\emph{best step}, with the model error that goes with it:" (and after the display: "The second formula follows from
putting $\abs{f'''} = 3 \epsmach \abs{f} / (\hbest)^3$ into the model: at the best step, the truncation term is half
the rounding term.")

Before (3.6): "A solver's tolerance is the noise level only when its convergence test is on the output itself,
relative to it; a flash that converges a residual (an energy balance, equal fugacities) in its own units, or a
simulator that converges tear streams, has a tolerance that is not $\noise$: measure $\noise$, as shown below."

After: "First, a solver's tolerance is the noise level only when its convergence test is on the output itself,
relative to the output. Some solvers test something else. A flash may converge a residual in its own units, such as
an energy balance or the equality of fugacities; a residual is the amount by which an equation fails to hold. A
simulator may converge its tear streams, the recycle streams whose values it guesses and then iterates. For such
solvers the tolerance is not $\noise$. Measure $\noise$ instead, as shown below."

## Compromises and points for the editor

1. **Layout guards (Chapter 3).** In the standalone build the new wording split two listings badly. The
   `narrow-feature` listing left an empty strip of its frame at the foot of a page, and `noise-estimate-use` left
   its first comment line alone on a page. Each now has a `\needspace` guard
   (`\needspace{15\baselineskip}` and `\needspace{13\baselineskip}`), with a comment saying why. The lead-in
   sentence stays at the page foot and the listing moves whole to the next page. Please confirm this at the full
   build: page breaks there may differ.
2. **`\hbest^3`** is a double superscript (the macro carries one already). The spelled-out step writes
   `(\hbest)^3`.
3. **Two sentences now name what they refer to.** I checked both against the sources.
   - Section 3.6 ended "differentiates the friction factor with this rule". After the new Ridders sentences,
     "this rule" could be read as Ridders' method. It now reads "with the noise-adjusted step of this section",
     which is what `07-using.tex` (Section 7, lines 229 to 236) does: $\noise^{1/3} s$.
   - "If no three orders agree" now reads "If the estimates of no three successive orders agree". This follows
     `code/noise.hpp`, which accepts the first $k$ whose $\sigma_k$, $\sigma_{k+1}$ and $\sigma_{k+2}$ agree.
   - "(order 2, 8 evaluations)" became "It settles at order 2, and uses 8 evaluations of $f$". "Settles" is the
     header's word.
4. **A loose claim, kept as written (content is frozen).** Section 3.5, "A small $|x|$": "a mole fraction of
   $10^{-6}$, a strain of $10^{-4}$ or a time constant in seconds of a fast process all have $s = 1$, so the
   default step is $\epsmach^{1/3}$, as large as the variable itself or larger." For `double`,
   $\epsmach^{1/3} \approx 6 \times 10^{-6}$. That exceeds the mole fraction but is about 16 times smaller than the
   strain. "As large as or larger" is true of the first example only. It may be worth a factual fix outside the
   pass (for example "comparable to the variable, or larger", or "a large fraction of the variable"). Not changed
   here.
5. **Glosses in the text, not in the margin.** These are flash calculation, tear stream, residual, Arrhenius factor
   and enthalpy of formation. They are written as definitions with no numbers. "The temperature dependence of a
   rate constant" for the Arrhenius factor and "finds the phases of a mixture in equilibrium" for a flash are
   standard definitions, not claims about the library.
6. **Signposting** was added at the openings of 3.1, 3.2, 3.3, 3.5 and 3.8, and of Appendix A. I added no closing
   summaries: the checklist of 3.8 already closes the chapter.
7. **No disagreement between the sources** was found in the passages touched.

No new data sets or references are needed.

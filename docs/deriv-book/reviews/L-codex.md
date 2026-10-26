# L-codex: independent review of the language pass (OpenAI Codex, gpt-6.1-sol, read-only)

Run 2026-10-26 by the main session on the uncommitted pass against 3192e8a; prompt: changes of meaning only.
The review follows verbatim; the decisions are at the end.

1. **L-codex-01 — major — [03-rounding.tex:196](/C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/03-rounding.tex:196)**

   **Original:** “makes the default step smaller than the model's best step, by that factor when the guess is right, $w = \scale$.”

   **New:** “makes the default step smaller than the model's best step. When the guess is right, $w = \scale$, it is smaller by that factor.”

   **Difference:** The sentence split makes “smaller” unconditional; the condition now qualifies only the factor. When the true scale is much smaller than the guessed scale, the default step can instead be much **larger** than the model’s best step. The implementation uses the guessed scale directly ([deriv.hpp:591](/C:/Dev/XLThermo/Numerixx/include/numerixx/deriv.hpp:591)). This changes the reader’s diagnosis of an unsuitable step.

   **Suggested:** “When the scale guess is right, $w = \scale$, dropping the constant makes the default step smaller than the model’s best step by that factor.”

2. **L-codex-02 — major — [01-why.tex:208](/C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/01-why.tex:208)**

   **Original:** “$f$ has a kink or a range break within a few steps of $x$”

   **New:** The same statement, followed by: “A kink is a point where $f$ or one of its derivatives jumps.”

   **Difference:** The added definition includes discontinuities of $f$ and jumps in arbitrary higher derivatives. The original book distinguishes a jump in $f$ from a kink: Chapter 6 identifies a kink as a jump in $f'$ while $f$ remains continuous. That distinction matters because the book says Ridders’ estimate can flag a jump but cannot reveal a kink at $x$. The new definition makes those diagnostic instructions misleading.

   **Suggested:** “A kink is a point where $f$ stays continuous but its slope jumps.”

3. **L-codex-03 — minor — [shared-method-summary.tex:16](/C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/shared-method-summary.tex:16)**

   **Original:** “$f$ accurate to $\noise > \epsmach$ only (a solver output)”

   **New:** “$f$ is accurate only to $\noise > \epsmach$, as the output of a solver is”

   **Difference:** An example has become an unqualified assertion about solver outputs. A solver can converge to machine precision; its output also need not have the irregular error assumed by the noise model. This wording appears in both Chapter 1’s summary and the quick reference, where readers choose their step.

   **Suggested:** “$f$ is accurate only to $\noise > \epsmach$, as can happen with solver output”

4. **L-codex-04 — minor — [02-taylor.tex:232](/C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/02-taylor.tex:232)**

   **Original:** “cancels the $h^2$ term.”

   **New:** “cancels the $h^2$ term: the subtraction leaves $3 f'(x)$ plus terms of order $h^4$, and the division by 3 leaves $f'(x)$.”

   **Difference:** The added algebraic explanation drops the remaining error after division. Dividing by 3 leaves $f'(x)+O(h^4)$, not the exact derivative. For example, for $f(x)=x^5$ at zero, the combination is $-h^4/4$, whereas the derivative is zero. The unchanged displayed equation correctly retains the remainder.

   **Suggested:** “The subtraction leaves $3f'(x)$ plus terms of order $h^4$. Dividing by 3 therefore leaves $f'(x)$ plus an error of order $h^4$.”

5. **L-codex-05 — minor — [08-api.tex:85](/C:/Dev/XLThermo/Numerixx/docs/deriv-book/chapters/08-api.tex:85)**

   **Original:** “\code{fn} is called as a named lvalue through \code{std::invoke}”

   **New:** The same statement, with the added definition: “An rvalue is a temporary value.”

   **Difference:** Rvalues also include expressions referring to existing objects, such as `std::move(existing_object)`. Equating rvalues with temporaries gives a false explanation of the callable requirement. The actual constraint is invocation through `F&`, regardless of whether the caller supplied a temporary ([math.hpp:67](/C:/Dev/XLThermo/Numerixx/include/numerixx/math.hpp:67)).

   **Suggested:** Replace the glossary note with: “Here, the named reference \code{fn} is an lvalue, even when the caller passed a temporary.”

I read these chapters in full and found their rewritten prose faithful:

- Chapter 4, *Higher-order and one-sided stencils*
- Chapter 5, *Second to fourth derivatives*
- Chapter 6, *Ridders’ extrapolation*
- Chapter 7, *Using the library*
- Chapter 9, *Tables*
- Chapter 10, *Accuracy, cost, limits and guarantees*, allowing the approved terminology change
- Appendix A, *Floating-point essentials*
- Appendix B, *Derivations*

The rewrite is largely faithful, but it is not yet entirely meaning-preserving. Two changes affect how readers interpret step selection and failure diagnostics; three smaller additions introduce misleading generalisations or explanations. I completed the requested word-diff comparison across all chapter files and checked the relevant library claims against both headers. No files were modified.

## Decisions (main session, 2026-10-26)

All five reproduced against 3192e8a and the headers, all accepted:

- L-codex-01: 3.3 now says dropping the constant changes the step by that factor, and that the default step is smaller than the best step by it when the guess is right.
- L-codex-02: the margin gloss in 1.5 now defines a kink as a point where f is continuous but its slope jumps, as 6.4 does.
- L-codex-03: Table 1.1 (and the quick reference) says "as the output of a solver can be".
- L-codex-04: 2.4 says the division leaves f'(x) with an error of order h^4.
- L-codex-05: the rvalue gloss in 8.2 now names a temporary or an expression that may be moved from, such as std::move(v).

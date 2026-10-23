# R1-engineer: review by the target reader (practising process engineer)

Reviewer: book-reader-engineer, round 1. I read the integrated build (`build/gcc/docs/deriv-book/pdf/main.pdf`, 83 pages; page numbers below are the printed ones). I read it as a process engineer who knows calculus, writes Python and Excel and a little C++, and needs derivatives of property models and flowsheet calculations. I did not check the mathematics or the code; where something looked wrong to me as a reader, I say so.

(The reviewer agent has no file-writing tool; the main session saved this review verbatim from the reviewer's report.)

## Findings (ranked by severity)

### R1-engineer-01: "If you read nothing else" gives the wrong one check
- **Location:** `chapters/01-why.tex` lines 144–180 (§1.5, p. 6), Table 1.1; same content on the quick-reference page (p. 76).
- **Severity:** major
- **Problem:** §1.5 says: "**The one check.** Test the result before you use it, as in the listings above: `if (d)`, then `d->value`."
  - To me "the one check" means "the check that the answer is right". But `if (d)` only tells me a number came back.
  - The book itself shows that the dangerous cases return a number without complaint:
    - a central rule straddling the 500 K range break is "4.3 off, relative" (§7.8, p. 48);
    - the default step on a narrow feature gives 2.8 × 10⁻⁵ instead of 5.7 × 10⁻¹² (§3.5, p. 16);
    - a solver output with η = 10⁻¹⁰ gives 5.5 × 10⁻⁵ (Table 7.1).
  - None of these is named in §1.5. Rows 4 and 5 of Table 1.1 mention a "scale w" and a noise "η", but at page 6 I do not know what either means, how to find it, or that ignoring it gives a silently wrong number.
  - A reader who really reads nothing else leaves with false confidence. That is worse than leaving with nothing.
- **Suggested change:** split "The one check" into two short paragraphs:
  - *Did it work:* `if (d)`.
  - *Is it right:* three silent failures in one line each, with pointers:
    - the step assumes f varies on max(|x|, 1) in your units (§3.5);
    - a kink or range break within a few steps gives a wrong value (§7.8);
    - a solver output needs a larger step (§3.6).
  - Then add "if in doubt, repeat at 2h and h/2 (§7.5)".
  - In Table 1.1's caption, add one plain sentence each for w ("the distance in x over which f changes appreciably") and η ("how many digits of f you can trust, e.g. 10⁻⁸ for a solver tolerance of 10⁻⁸").

### R1-engineer-02: the Jacobian of a set of balances is promised but not delivered
- **Location:** `chapters/01-why.tex` line 6 (p. 2); `chapters/07-using.tex` lines 152–157 (§7.4.2, p. 43); `code/07-using-sensitivity.cpp` lines 85–108; §7.9 (p. 49).
- **Severity:** major
- **Problem:**
  - Chapter 1 opens with "a Jacobian is the slope of a set of balances", and Newton's method appears in the list of uses. That is exactly my daily need: a Newton or Broyden solve of a flowsheet or flash, m residuals in n unknowns.
  - §7.4.2 then says: "A Jacobian is the same thing done once per parameter: a loop over the parameters ... (snippet `jacobian` in `code/07-using-sensitivity.cpp`)". Three problems with that:
    - The snippet is not printed in the book. I am sent to a source file, for the one listing I most wanted.
    - What the snippet computes is `J[0], J[1]` for the scalar `rate(A, Ea, T)`. That is a 1×2 gradient, not a Jacobian of a set of balances.
    - §7.9 then says "gradients, Jacobians and Hessians" are deferred, with "a loop of calls, as in Section 7.4.2" as the stopgap. So I am told a Jacobian is both shown and deferred.
  - The real question is never answered. My model evaluates all m residuals (or a flash returns T, the vapour fraction and every composition) in one expensive call, but `derivative` takes a scalar function. How do I get column j without calling the whole model m times per evaluation, and what does that cost?
- **Suggested change:** print the listing in §7.4.2 and make it a genuine m×n Jacobian of a small vector function, for example two residuals of an energy and a mass balance. Show one of two things:
  - (a) one `derivative` call per (row, column), with its cost of m × n × evaluations stated; or
  - (b) a fixed-step stencil applied by hand to the vector, using the weights of Table 9.1, with the cost of 2n model calls for `Central<2>`.

  State plainly which of these the library supports today. Rename the present example "a gradient" or "parameter sensitivities".

### R1-engineer-03: I cannot obtain the noise level η for my own solver
- **Location:** `chapters/03-rounding.tex` lines 318–325 (§3.6, p. 18–19); Table 1.1 row 5 (p. 6); Table 7.1 row 5 (p. 41); §7.4.3 (p. 43–44).
- **Severity:** major
- **Problem:**
  - The cure for a noisy f is clear: "pass h = η^{1/3} s ... with η the solver's tolerance or your best estimate of the accuracy of f".
  - In the friction-factor example η is the bisection tolerance on the output itself, so it is easy. Real solvers are not like that:
    - a flash converges on a residual (fugacity equality, an energy balance) in its own units;
    - a simulator converges tear streams;
    - none of these tolerances is "the relative accuracy of the output".
  - The fallback is one paragraph: "Evaluate f at a handful of points ... form their successive differences ... their spread estimates η|f|. Moré and Wild turn this idea into a reliable procedure. The library does not do it for you". There is no code, no number of points, no spacing, and no worked case.
  - A solver output is one of the five situations the book says "cover most uses", and this is the only one I could not carry out after reading.
- **Suggested change:**
  - Add a compiled 10–15-line snippet that estimates η from a difference table at, say, 8 points spaced η-scale apart. Run it on the friction factor and show that it recovers about 10⁻¹⁰. Quote how many evaluations it costs.
  - Add one sentence warning that a residual tolerance is not η, and that η should be measured when the solver's convergence test is not on the output.

### R1-engineer-04: no complete first program, and no "how do I get it into my build"
- **Location:** `chapters/01-why.tex` lines 112–131 (§1.4, p. 4–5); `chapters/07-using.tex` lines 14–16 (§7.1, p. 39); `chapters/08-api.tex` §8.1 (p. 51); §10.5 (p. 64).
- **Severity:** major
- **Problem:**
  - Every listing is a fragment. The first call, `nxx::deriv::derivative(H, T0)`, appears without `#include <numerixx/deriv.hpp>`, without `<print>`, and without `main`.
  - `book::shomate_H(book::nitrogen, T)` uses a namespace `book::` that is never explained.
  - The include line first appears in §7.1 (p. 39) and the `#include` listing in §8.1 (p. 51).
  - Nowhere does the book tell me:
    - how to obtain the library (a CPM or FetchContent line, or a copy of the headers);
    - which compiler versions work. C++23, `std::expected` and `std::println` rule out many installed compilers. §10.5 names the CI toolchains (GCC, Clang, MSVC, clang-cl, Emscripten) but no versions.
  - As someone who writes "a little C++", I would get stuck before my first derivative.
- **Suggested change:**
  - Make the first listing in §1.4 a complete, compiled program: includes, a lambda for H(T), `main`, the call, `if (d)`, print.
  - Put a short "Getting the library" box in §7.1: a CMake snippet (FetchContent or CPM plus `target_link_libraries(... numerixx::numerixx)`) and the minimum compiler versions tested.
  - Say once that `book::` is the namespace of the book's own example header.

### R1-engineer-05: a large |f| (enthalpy of formation) is mentioned once and missing from the checklist
- **Location:** `chapters/03-rounding.tex` lines 132–138 (§3.3, p. 14) and 352–366 (§3.8, p. 19); §10.1 (p. 61).
- **Severity:** minor
- **Problem:**
  - "Had the enthalpy been referred to 0 K, or included an enthalpy of formation, |H| would be larger and the derivative proportionally less accurate". This is the single most relevant remark in the book for my work.
  - Simulator enthalpies normally include formation enthalpies: hundreds of kJ/mol against a c_p of tens of J/(mol K). Pressures in Pa with a small slope are similar.
  - The remark sits in the middle of a paragraph about ε^{2/3}, and the five-question checklist in §3.8 does not ask it. I nearly missed it.
- **Suggested change:**
  - Add a sixth checklist question: "Is |f| much larger than |x f'| (an enthalpy with a formation term, an absolute pressure)? Expect to lose log₁₀(|f| / (s|f'|)) digits; subtract a constant reference from f before differentiating if you can."
  - Give one number for a realistic case, for example water vapour with its formation enthalpy, as a generated value.

### R1-engineer-06: the letter s changes meaning between §3.3 and §3.4
- **Location:** `chapters/03-rounding.tex` lines 118–125 against 146–150 (p. 14–15); §3.5 (p. 16).
- **Severity:** minor
- **Problem:**
  - In §3.3, s is the scale on which f actually varies: "Suppose f varies on a length scale s, so that each derivative is about the previous one divided by s".
  - In §3.4, s = max(|x|, 1), the library's assumption.
  - In §3.5 the true scale becomes w.
  - I read §3.4 as "the library computes the scale of f", and only in §3.5 understood that s is a guess. The whole lesson of the chapter (s is a guess, w is the truth) hangs on this distinction.
- **Suggested change:** use w (or another symbol) for the true scale from §3.3 onwards, and introduce s only in §3.4 as "the library's guess at w".

### R1-engineer-07: developer and test detail in the Part I main text; `d` used before it is defined
- **Location:** `chapters/03-rounding.tex` line 102 (Sidebar 3.1, p. 14); `chapters/04-stencils.tex` lines 52–60 (§4.1, p. 21) and 170–178 (§4.3, p. 23).
- **Severity:** minor
- **Problem:**
  - Sidebar 3.1 says "Rules with more nodes divide by dh". At p. 14, d has not been defined; it first appears on p. 21.
  - §4.1, in the main text, then explains how the library stores weights as integer numerators, and that "The tests use this with a point and a step that are short binary fractions ... so that `Central<4>` on x⁴ returns exactly 4 at compile time and at run time."
  - §4.3 says "That is why the tests budget (C + 2S) ε^{P/(P+1)} for every multi-node rule, an estimate rather than a bound, below which the measured errors stay by a factor of 4 to 16."
  - These are the first sentences in Part I where I lost the thread. I do not know what "the tests budget" means, why I should care, or what d is.
  - The equations of §4.1 itself, (4.1)–(4.3), I could follow with effort.
- **Suggested change:**
  - Move the integer-numerator storage, the binary-fraction test and the "tests budget" sentence into Sidebar 4.1 or Chapter 10.
  - In Sidebar 3.1, write "divide by a fixed multiple of h" instead of "dh".
  - In §4.3, keep only the practical statement: rules with more than two nodes carry an extra rounding term from the points, about (S/2) ε|x f'|/h.

### R1-engineer-08: "The solver matters" ends without telling me what to do
- **Location:** `chapters/07-using.tex` lines 186–190 (§7.4.3, p. 44).
- **Severity:** minor
- **Problem:**
  - "A fixed-point iteration stopped at the same tolerance is smooth between the values of Re at which its iteration count changes, but off by up to its tolerance: Ridders' method on it was 2.7 × 10⁻¹¹ off with a relative estimate of 1.4 × 10⁻¹² (measured). A method differentiates the function it is given; no estimate sees how far that is from the one you meant."
  - I had to work out for myself that the estimate failed to cover the error by a factor of about 20; the text does not say so.
  - I also could not tell which kind of solver I have (bisection-like noisy, or fixed-point-like smooth but biased), or what to change: tighten the inner tolerance? fix the iteration count? warm-start?
  - This matters to me: I would use the flash solver as it is, and need to know whether to trust the result.
- **Suggested change:**
  - State the conclusion first: "a solver can be smooth but wrong by its tolerance; then no derivative method, and no estimate, can see it."
  - Then give the action: converge the inner solver tighter than the derivative accuracy you need, ideally to η ≤ ε × 10³; keep the iteration count fixed across the stencil's points if your solver allows; or differentiate the converged equations implicitly.
  - Say explicitly that the estimate did not cover the error here.

### R1-engineer-09: too many ways to read the result; the book does not recommend one
- **Location:** `chapters/01-why.tex` §1.4 listing (`*r->abs_error`, p. 5); `chapters/07-using.tex` lines 51–59 (§7.2, p. 40–41); §7.5 listing (`digits(*d2, *d)`, p. 45); §7.6 (`mp50(*r)`, p. 46).
- **Severity:** minor
- **Problem:**
  - Ways to read the result that appear in the book:
    - `d->value`;
    - `*d` (an `Estimate`, which converts);
    - `double v = *d`;
    - `auto e = *r` (a trap);
    - `*r->abs_error`;
    - `r->abs_error.has_value()`;
    - `mp50(*r)`.
  - §7.2 lists "three traps" of the conversion, and then §7.5 uses the conversion anyway (`digits(*d2, *d)`).
  - On p. 5, `*r->abs_error` meets me before I know that `abs_error` is an optional. I read it as a multiplication.
  - As an occasional C++ user, I want one safe idiom and the reason for it.
- **Suggested change:**
  - State one house rule in §1.4 and keep to it in every Part I listing: read `d->value`, and for Ridders `d->abs_error.value()`.
  - Mention the implicit conversion only in §7.2 and Chapter 8, as an alternative with its traps.
  - Add a margin note at the first `*r->abs_error` ("an optional; `*` reads it").

### R1-engineer-10: mixed partial derivatives and several outputs of one expensive call
- **Location:** §7.4.2 (p. 42–43), §7.9 (p. 49), §10.7 (p. 64–65).
- **Severity:** minor
- **Problem:**
  - Property work needs mixed partials, such as ∂²A/∂T∂V or ∂²G/∂T∂P for an equation of state. The book lists "mixed partials and Hessians" as deferred (§10.7) but gives no interim advice: may I nest two `derivative` calls, and with what accuracy and cost? Or should I use a four-point cross formula?
  - Equally, there is nothing on differentiating several outputs of one expensive call (a flash returning T and the phase fractions) without repeating the call per output.
- **Suggested change:**
  - Add a short paragraph to §7.9 or §7.4.2 with the recommended interim approach for a mixed partial (a nested call, or the cross formula), its cost in evaluations, and its expected accuracy as a generated value.
  - Add one sentence on caching the expensive call so that several outputs share the evaluations.

### R1-engineer-11: "the 1 is in your units" needs an engineering example of units
- **Location:** `chapters/03-rounding.tex` line 258 (§3.5, p. 17) and the Fig. 3.3 caption.
- **Severity:** minor
- **Problem:**
  - The y ln y example at y = 10⁻⁶ is excellent. But the general lesson ("the same physical function gets a different default step depending on whether I work in K or °C, Pa or bar, mol or kmol") is left for me to infer.
  - Example: a property near 20 °C gets s = 20 and a step about 15 times smaller than the same property at 293 K. I had to work that out myself.
- **Suggested change:** add two sentences with one units example (K against °C, or Pa against bar) and its consequence in digits, as a generated value or with "about".

### R1-engineer-12: numbers bury the conclusions in several paragraphs
- **Location:** `chapters/06-ridders.tex` lines 156–166 (§6.4, p. 35); `chapters/07-using.tex` lines 176–184 (§7.4.3, p. 44) and 216–222 (§7.5, p. 45).
- **Severity:** minor
- **Problem:**
  - Several paragraphs carry four to six measured numbers each, with the conclusion at the end or only implied.
  - Example: "Ridders returns 100.000 000 000 to 5.0 × 10⁻¹⁴ relative, covered by its estimate, after 24 evaluations against about 14 in the interior: a restart costs about two evaluations per halving."
  - In §7.5: "which is 11.3 digits of |f|/s and only 6.1 of f'". I did not understand "digits of |f|/s".
  - I trust the numbers, but on a working read I want the sentence I should remember first.
- **Suggested change:**
  - Lead each such paragraph with its practical conclusion, then give the evidence.
  - In §7.5, replace "digits of |f|/s" with plain words: "the absolute error is tiny, about 10⁻¹² in units of f per unit x, even though only six digits of the small slope agree".

### R1-engineer-13: contradictory guidance on tabulated or interpolated data
- **Location:** `chapters/03-rounding.tex` lines 276–277 (§3.5, p. 18); `chapters/07-using.tex` lines 358–362 (§7.9, p. 49).
- **Severity:** minor
- **Problem:**
  - §3.5: "For an interpolated table it [the scale w] is the spacing of the table's nodes". This reads as "differentiate an interpolated table with h on that scale".
  - §7.9: "a linear interpolant has a piecewise-constant derivative and a kink at every node, so a smooth fit or a spline is the place to start".
  - Steam tables, vendor curves and plant data are everyday inputs for me. I came away unsure whether to differentiate my interpolant numerically or not.
- **Suggested change:**
  - In §3.5, say "for a smooth interpolant (a cubic spline) it is the node spacing", and point to §7.9.
  - In §7.9, add: "with a spline, use the spline's own derivative when it offers one; numerically, keep the stencil within one interval or accept the kink error of §7.8".

### R1-engineer-14: C++ vocabulary assumed in the Part I chapter 7
- **Location:** `chapters/07-using.tex` line 48 (§7.2, p. 40); §8.3 (p. 54); front matter "comfortable with calculus and with C++".
- **Severity:** minor
- **Problem:**
  - "the standard monadic operations work on it" is jargon to me; I know the word only from functional programming talk.
  - Other terms such as "aggregate", "designated initialisers", "structured bindings", "lvalue" and "requires-expression" are fine in Chapter 8 but are not glossed anywhere.
  - The front matter promises a reader "comfortable with C++". The target engineer may write C++ only occasionally.
- **Suggested change:**
  - Replace "monadic operations" with "member functions `transform` and `and_then`, which apply a function only when there is a value".
  - Either add the C++ terms to the glossary or keep them out of Part I.

### R1-engineer-15: no code for "pass ε^{1/3} w"; `root_epsilon` appears unexplained
- **Location:** `chapters/03-rounding.tex` §3.5 listing (p. 17); `chapters/qr-quick-reference.tex` (p. 76); Table 1.1 (p. 6).
- **Severity:** minor
- **Problem:**
  - The quick reference and Table 1.1 tell me to pass h = ε^{1/3} w or η^{1/3} s, but not how to write that in code.
  - The §3.5 listing introduces `root_epsilon<3, double>()` without a word of explanation. I guessed it is the cube root of ε, and it is only confirmed in §8.7 (p. 55).
- **Suggested change:**
  - Add one line to the quick reference: `h = nxx::math::root_epsilon<3, double>() * w; // eps^(1/3) * w`, plus `std::cbrt(eta) * std::max(std::abs(x), 1.0)` for the noisy case.
  - Add a margin note at the first use of `root_epsilon` on p. 17.

### R1-engineer-16: two sidebars are in the wrong place or not optional
- **Location:** `chapters/03-rounding.tex` line 203 (Sidebar 3.3, p. 16); `chapters/07-using.tex` §7.6 (Sidebar 7.1, p. 46).
- **Severity:** minor
- **Problem:**
  - I could skip Sidebars 2.1, 3.2, 4.1, 4.2, 5.1, 5.2, 6.1, 6.2, 6.3 and 7.2 without losing the thread, and I did. They are truly optional. Two are not:
    - Sidebar 7.1 is required for anyone who uses Boost.Multiprecision. The main text says "Two details of the snippet matter ... `-> mp50` on the lambda and `mp50(*r)` ... Sidebar 7.1 explains both", and the main text never says what to do.
    - Sidebar 3.3 ("Machine epsilon, ulps and the rounding of x + h") sits after §3.4. I expected it at §3.1, where ε is introduced.
- **Suggested change:**
  - Give the two rules in the §7.6 main text in one sentence ("give the lambda an explicit `-> mp50` return type, and write `mp50(*r)` or `r->value` in Boost arithmetic") and leave the why to the sidebar.
  - Move Sidebar 3.3 next to §3.1, or point to it from there.

### R1-engineer-17: the float-warning pragma listing is heavy for Part I
- **Location:** `chapters/05-orders.tex` §5.4, line 148 onwards (p. 29).
- **Severity:** minor
- **Problem:** a 20-line listing of compiler-specific `#pragma` blocks shows how to silence a warning for a fourth derivative in `float`. The book's own advice is never to compute that ("Cure: double"). For an engineer this is the least useful listing in Part I, and it breaks the flow between §5.3 and §5.5.
- **Suggested change:** keep the two sentences on what the warning means and the cure ("use double"). Move the pragma listing to §8.6, which already points back to it.

### R1-engineer-18: appendices are cross-referenced as "Chapter A" and "Chapter B"
- **Location:** `chapters/03-rounding.tex` line 28 (p. 12, "Chapter A has the details"); Sidebar 3.2 (p. 15, "The derivation in full is in Chapter B").
- **Severity:** editorial
- **Problem:** the references read "Chapter A" and "Chapter B". The contents page calls them "Appendix A" and "Appendix B", and for a moment I looked for chapters I had missed.
- **Suggested change:** fix the `\cref` name for appendices so they render as "Appendix A".

### R1-engineer-19: references to artefacts the reader has not seen
- **Location:** `chapters/05-orders.tex` line 146 ("as the header warns"); `chapters/06-ridders.tex` line 156 ("the header's own example"); Sidebar 6.2, §9.5 and §10.2–10.7 ("the design record", `docs/PHASE3_API.md`); every figure caption ("Data `data/step-error-cp.csv`").
- **Severity:** editorial
- **Problem:**
  - I have not read the header or the design record, so "as the header warns" asks me to trust a document I have not seen.
  - The data-file path at the end of every caption is noise when reading for understanding. It belongs to the evidence trail, not to the reader's caption.
- **Suggested change:**
  - In Part I, replace "as the header warns" and "the header's own example" with the plain statement.
  - Move the data-file paths from the captions to a list of figures in the back matter, or set them in a smaller, separate line.

### R1-engineer-20: "The example checks ..." sentences describe the test harness, not the method
- **Location:** `chapters/04-stencils.tex` line 254 (p. 25); `chapters/05-orders.tex` line 190 (p. 30); §2.3 (p. 10–11, "the example checks that each lies within 0.1").
- **Severity:** editorial
- **Problem:** "The example checks that the two results compare equal, value and estimate alike, which the defaulted == of Estimate does." These tell me what the book's CI asserts, not what I learn. They read as filler.
- **Suggested change:** delete them, or keep them only where the check carries a lesson, and say what that lesson is.

### R1-engineer-21: O(h²) used one page before it is defined
- **Location:** `chapters/02-taylor.tex` line 95 (§2.2, p. 8) against line 140 (§2.3, p. 9).
- **Severity:** editorial
- **Problem:** "its error is −h f''(x)/2 + O(h²)" comes before §2.3 explains the O notation. I knew it, but the target reader may not.
- **Suggested change:** add a margin gloss at line 95 ("terms of order h² and smaller; §2.3").

### R1-engineer-22: ambiguous sentence in §6.4
- **Location:** `chapters/06-ridders.tex` line 159 (p. 35).
- **Severity:** editorial
- **Problem:** "Ridders returns 100.000 000 000 to 5.0 × 10⁻¹⁴ relative". I first read it as "returns values from 100 to 5 × 10⁻¹⁴".
- **Suggested change:** "Ridders returns 100.000 000 000, accurate to 5.0 × 10⁻¹⁴ relative and covered by its estimate, after 24 evaluations ...".

### R1-engineer-23: glossary gaps
- **Location:** Notation and glossary (p. 72–74).
- **Severity:** editorial
- **Problem:** the glossary lacks several terms that the text uses as technical words:
  - "cancellation";
  - "amplification" (S);
  - "rounding floor";
  - "moment conditions";
  - "Richardson extrapolation";
  - "aliasing";
  - "unit roundoff";
  - "untrusted result".

  I looked up "floor" while reading §6.3 and did not find it.
- **Suggested change:** add these entries, one line each, with the section where each is defined.

### R1-engineer-24: Figure 6.4 is hard to read
- **Location:** `chapters/06-ridders.tex` line 179 (Fig. 6.4, p. 35).
- **Severity:** editorial
- **Problem:** two stacked panels share four vertical lines labelled "cap", "10 digits" and "last value", plus a curve, and the x-axis is k with x = 2⁻ᵏ. I needed the caption and the text side by side to decode it. The message ("Ridders survives an edge down to about 10⁻¹⁰ in x") does not jump out.
- **Suggested change:** put a second axis or tick labels in x (10⁻⁸, 10⁻¹⁰), label the three regions directly ("works", "ten digits lost", "fails"), and drop one of the vertical lines.

### R1-engineer-25: visible draft markers
- **Location:** front matter, "The version" (p. ii, "[TODO: the commit of dev at release]"); References (p. 75, "UNVERIFIED" on 19 of 21 entries).
- **Severity:** editorial
- **Problem:** as a reader, "UNVERIFIED" against almost every reference undermines the book's claim to evidence on its first read. I understand that it is a process marker for this round.
- **Suggested change:** none beyond the planned release steps. Make sure both disappear before release.

## What works

- **The V.** Fig. 3.1 together with the six-digit pocket-calculator table (p. 13) is the clearest explanation of why a smaller step is not better that I have read. I could explain it to a colleague now.
- **Real engineering cases.** These are my problems, not toy one-liners:
  - c_p from a Shomate enthalpy;
  - y ln y at a mole fraction of 10⁻⁶;
  - the 500 K range break between two published fits;
  - Newton on an energy balance;
  - the Colebrook friction factor;
  - Arrhenius sensitivities with the scaled sensitivity d ln k / d ln E_a.
- **Actionable lessons.** "Evaluable is not valid" with the NaN wrapper (§4.4) is something I will put into my property wrappers tomorrow. "Evaluable is not smooth" with Fig. 7.2 explains a class of bugs I have actually had.
- **Practical chapter 7.** §7.5 "Is my derivative good enough?" (2h and h/2, a higher order, Ridders, the absolute error near zero, the physical sign check) and §7.8 Pitfalls in symptom/cure form are exactly what I need at the desk.
- **Honesty.** The book says that Ridders' estimate is not a bound, quantifies its coverage (99.7 %), and shows it failing under noise (90.7 %, worst ratio 42). It separates "measured", "derived" and "rule of thumb". I trust it more for that.
- **Tables 1.1 and 7.1** with an evidence column, and the one-page quick reference, are what I would photocopy.
- **The sidebars** are almost all truly optional and clearly boxed. I skipped ten of them without losing the thread.
- **The plots** help, notably Figs. 2.1, 3.1–3.5, 5.2, 6.3, 7.1 and 7.2: each caption says what to see.
- **Tone.** Nothing read as patronising. The prose is plain and the voice consistent; the fault, where there is one, is density rather than condescension.

## Overall verdict

Yes, I would keep this book at my desk: Chapters 1, 3 and 7 and the quick reference, with Chapters 2, 5 and 6 read once and Part II consulted rarely (I would never use the Chapter 9 weight tables). After reading it I can choose between a stencil and Ridders, I know when the default step is wrong and why, and I know how to sanity-check a result. Four things stop it from being the book I reach for first:
- a complete, compilable first program with instructions for getting the library into a CMake build and the supported compiler versions;
- a real Jacobian of a residual vector, with its cost, rather than a gradient left in a source file;
- a measured-η recipe in code for solver outputs whose tolerance is not on the output;
- a §1.5 that names the silent failures instead of offering `if (d)` as "the one check".

With those, and with the developer detail of §4.1 and §4.3 moved into sidebars, it would be the reference I hand to new engineers in my group.

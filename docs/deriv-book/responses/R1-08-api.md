# R1-08-api: response to review round 1

Author: `book-author-practice` (Opus 5.5), stage 6. File: `chapters/08-api.tex`; it now quotes E8's `float-fourth`
(`code/05-orders-float.cpp`, the methods author's, unchanged by me) and E23's `usertype-limits`
(`code/07-using-usertype.cpp`, new). `tools/check_chapter.sh`: `COMPILED OK`, 7 pages standalone (6.35 full pages,
against 5.35 before). `tools/check.py`, including its new check of the `signature` blocks against the headers: no
error, no warning.

| Finding | Decision | What changed |
|---|---|---|
| developer-04, library-03 | accept | 8.2 Preconditions: "checked by an assert ... Under NDEBUG the check is compiled out and the behaviour unspecified (Forward<1> with a negative step evaluates below x); validate a computed step, or use the overload without a step, which has no precondition" (U-6) |
| library-01 | accept | 8.2 Preconditions: "a two-node rule fails only when its two points coincide, and near a power of two a step can move x on one side only, so that Central<2> divides by a half-ulp distance and returns a finite value that is mostly rounding"; Failures: "the points of a two-node rule coincide, or, for a rule with more nodes, x + h or x - h rounds to x or neighbouring points round to the same value" |
| library-07 | accept | 8.2 Failures: an untrusted result (abs_error at least a quarter of |value| and far above rounding) makes Ridders restart, "and is nevertheless returned, as a value with its large abs_error and not as a failure, when the step reaches the resolution of x or the levels run out during the restarts. Test abs_error before relying on such a value" |
| developer-06 | accept | 8.2, template parameter T: "the conversion is implicit, so for a float x write 1e-3f, since a double step is diagnosed under -Wconversion" |
| developer-10 | accept | 8.2, template parameter F: "a pointer to a member function of T itself. For a member function of your own model, capture the object in a lambda" |
| developer-12 | accept | 8.2 Default steps: "there is no accessor for it" (U-12) |
| evidence-23 | accept | The caps are `\valRiddersMaxEvals{First..Fourth}` in 8.2 Cost and in Table 8.2's Ridders row |
| T11 (developer-03) | accept | 8.2: "neither overload is noexcept"; 8.3: "`r.value()` is not the member but the accessor of std::expected, which throws std::bad_expected_access when there is no result (or terminates without exceptions); read `r->value` after testing r" (U-8) |
| library-10, developer-18 | accept | 8.3: abs_error "saturates at max() where the estimate would overflow (a large fn over a tiny step), and then says only that the method could not estimate its error"; "comparing with a T, as `*r == 2.0` or `r == 2.0` on the std::expected, compares the value only, and the latter is false when r holds an error" |
| library-22, developer-17 | accept | 8.3: if (*r) "diagnosed by Clang under -Wconversion, not by GCC"; int n = *r "truncates, with a -Wfloat-conversion warning" |
| library-06 | accept | 8.5: the weight check "is in the body of the evaluator, not in a constraint: a requires-expression reports such a call as valid (for Central<10> with std::float16_t, say), and the call then fails to compile inside the library" (U-10). 8.8: "A type too narrow for a rule's integer weights is rejected by a static_assert inside the library instead (Section 8.5), which a requires-expression does not see" |
| engineer-17 | modified | 8.6 quotes the pragma listing (`\snippet{05-orders-float}{float-fourth}`) with this call's own numbers: Central<2> `\valFloatFourthCentralAtOneRelError`, Ridders `\valFloatFourthRiddersAtOneRelError` with a relative estimate of `\valFloatFourthRiddersAtOneRelEstimate` (the AtOne macros of DS6; the E8 call is e^x at 1, developer-16). The pointer to 5.4 stays |
| developer-01 | accept | 8.7: the concept's conditions, then "The algorithms also read radix, digits and denorm_min(), which the concept does not check: a type without them passes the constraints and then fails to compile inside the header", and the compiled minimal type quoted (`usertype-limits` of E23: a wrapper of double without a default constructor, numeric_limits specialised with exactly is_specialized, is_integer, is_bounded, is_exact, radix, digits, epsilon(), max(), denorm_min(), not inheriting numeric_limits<double>; E23 static_asserts `nxx::real<Real>` and checks both calls on e^x against e to 1e-9) (U-5) |
| library-11 | accept | "One call of each form, on the thread's enthalpy" |
| library-18 | accept | Table 8.1 is placed with `[h]` at the head of 8.1, before the synopsis text; the lead-in "Two overloads, with a step and with the default step." ends with a full stop |
| T14 (library-15, evidence-28) | modified | The `signature` blocks are unchanged and pass `check.py`'s new check against the headers |
| engineer-19 | modified | 8.2 Default steps: "The header advises passing a step ..." is now the plain instruction "Pass a step when ..." |

Not changed: the names snippet of E11 is no longer quoted in 8.1 (a sentence names the using-declarations and the
file), to make room for the two listings the round added; Chapter 7 no longer quotes it either, and Chapter 1's
complete program shows the include.

Length: 6.35 pages against 5.35, +1.0 where the plan allows +0.5; the moved pragma listing (20 lines) and E23's
member list (18 lines) are most of it.

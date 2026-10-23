# R1-library: review of the integrated book against the code

Reviewer: book-reviewer-library, round 1. Subject: the integrated book at `aeec7e6` (`chapters/`, `code/`, the PDF
`build/gcc/docs/deriv-book/pdf/main.pdf`, built at 21:33:11, after the last source change at 21:33:02). Ground truth:
`include/numerixx/deriv.hpp` and `math.hpp`, then `tests/`, then `docs/design.md` and the PHASE documents.

Method: every chapter and every example read against the headers and the tests; all 20 examples compiled and run with
GCC 16 and with Clang 22 (llvm-mingw, libc++), both with the book's flags `-Wall -Wextra -Wpedantic -Wshadow
-Wconversion -Wsign-conversion -Werror`; the Part II pages and the Chapter 7 listings inspected in the PDF; each claim
below that is not a plain reading of the header checked with a scratch program in
`C:\Users\kenne\AppData\Local\Temp\claude\C--Dev-XLThermo-Numerixx\5424a43d-ee2e-47d6-89e3-20a17882b9f6\scratchpad\r1lib`
(`claims1.cpp`, `cx1.cpp`, `boost1.cpp`, `half.cpp`, `neg.cpp`, `second.cpp`, `warn.cpp`).

Summary: no major finding; 10 minor, 12 editorial. The statements about the API are accurate almost everywhere; the
findings are edge cases stated more strongly than the code supports, a listing that reads a result before testing it,
a few test claims that overstate the tests, and gaps in Part II's completeness.

## Findings

### R1-library-01 (minor): "a step below the resolution of x gives non_finite" is not true of the central difference

Location: `chapters/03-rounding.tex` 212-214 (Sidebar 3.3); `chapters/a1-floating-point.tex` 75-76; `chapters/08-api.tex`
85-86 ("a step below the resolution of x is detected and returns non_finite"); echoed in `chapters/07-using.tex` 44-45.

Problem. Sidebar 3.3 says: "If h is below half a spacing of x, fl(x + h) = x: the step has no effect, and the library
reports error::non_finite, because no difference can be formed"; Appendix A.4 says the same. For the two-node rules
the library detects only coinciding points (`deriv.hpp` 801, `if (dx == T(0))`): `Central<2>` for f' fails only when
*both* x + h and x - h round back to x. At a power of the radix the spacing below x is half the spacing above, so a
step below half the upper spacing moves x - h but not x + h, a difference is formed over a distance of half an ulp,
and the result is a finite value that is mostly rounding. The sidebar's next sentence explains this asymmetry, but
draws the conclusion only for the multi-node rules.

Evidence (`claims1.cpp`, fn = 3x, h = 0.4 eps or 0.4 ulp):

```
x=1 h=8.88e-17: x+h==x true x-h==x false | C2 4 B1 4 F1 non_finite
x=-1 h=8.88e-17: x+h==x false x-h==x true | C2 4 B1 non_finite F1 4
x=2 h=1.78e-16: x+h==x true x-h==x false | C2 4 B1 4 F1 non_finite
x=1.5 h=8.88e-17: x+h==x true x-h==x true | C2 non_finite B1 non_finite F1 non_finite
```

`derivative<Central<2>>(fn, 1.0, 0.4 * eps)` returns 4 for a derivative of 3, without failing. The header claims no
more than this: non_finite when "h is below the resolution of x", with the precondition "h well above ulp(x)"
(`deriv.hpp` 1395).

Suggested change. Sidebar 3.3: "If h is below half a spacing of x on both sides, both points round back to x and the
library reports non_finite. Near a power of two a step can move x on one side only; the central difference then
divides by the half-ulp distance that remains and returns a finite value that is mostly rounding. Hence the
precondition: a step well above ulp(x)." Appendix A.4 and 8.2 accordingly ("detected when x + h and x - h both round
back to x").

### R1-library-02 (minor): the pitfalls listing reads the result before testing it

Location: `code/07-using-expected.cpp` 87-90, snippet `estimate-pitfalls` (Section 7.2, PDF p. 48, book p. 41).

Problem. The listing writes `auto const e = *r;` and `double const v = *r;` and only then `if (r) ...`. Dereferencing
an empty `std::expected` is undefined behaviour; it works here only because Ridders on cos at 0 succeeds. Two
paragraphs earlier the text says "The usual pattern is to test the result, then read it", and Section 1.5 calls the
test "the one check". A reader who copies the pitfall listing copies the unchecked access.

Suggested change. Put the reads inside the test:

```
auto const r = derivative<Ridders>(cosine, 0.0);
if (r) {
    auto const e = *r;     // an Estimate<double>, not a double
    double const v = *r;   // the value, by the conversion
    ...
}
```

### R1-library-03 (minor): the precondition is an assert; in a release build a wrong step is not caught, and "Forward<P> never below x" then fails

Location: `chapters/08-api.tex` 84-87; `chapters/10-guarantees.tex` 163-166 (guarantee "Where fn is called");
`chapters/qr-quick-reference.tex` 30.

Problem. The book says the precondition `h > 0` is "checked with assert". `assert` (`deriv.hpp` 1337) is compiled out
under `NDEBUG`, so a release build does not check it, and the behaviour then depends on the rule. Chapter 10 states
"Forward<P> never below x, Backward<P> never above" as a guarantee without its hypothesis (style decision 15).

Evidence (`neg.cpp`, `-DNDEBUG`, fn = x^2 at x = 1, h = -0.25, and the lowest point fn was called at):

```
Forward<1> h=-0.25: 1.75 lowest point 0.75
Forward<2> h=-0.25: non_finite lowest point 1e+09
```

`Forward<1>` with a negative step silently evaluates below x and returns a backward difference; the multi-node rules
return non_finite (their point-order check, `deriv.hpp` 1051); Ridders returns non_finite (`h_i > 0` fails).

Suggested change. 8.2: "checked with assert in builds without NDEBUG (and as a compile error in a constant
expression); a release build does not check it". 10.5: "Forward<P> never below x and Backward<P> never above, for a
step h > 0 (the precondition)". The quick reference may keep "a bug (assert)".

### R1-library-04 (minor): cpp_bin_float_50 has no expression templates, but the listing and the text say the `-> mp50` matters for it

Location: `code/07-using-multiprecision.cpp` 21-22 (comment "-> mp50: a number, not an expression template");
`chapters/07-using.tex` 238-239 ("Two details of the snippet matter with such types, -> mp50 on the lambda and mp50(*r)
in the arithmetic"), 252-253.

Problem. `cpp_bin_float_50` is `number<cpp_bin_float<50>>`, whose default is `et_off` for a backend without an
allocator; `exp(x)` of it returns a `cpp_bin_float_50`, not an expression. The `-> mp50` is harmless but does nothing
for this type; it matters for `number<cpp_bin_float<50>, et_on>`, the type the tests use for the expression-template
case (`tests/test_deriv.cpp` 13 and 95, where that lambda wraps its result in `bin50_et(...)`). The second detail is
right for both: `*r - e` does not compile with or without expression templates.

Evidence (`boost1.cpp`, Boost 1.92, GCC 16): `expression_template_default<cpp_bin_float<50>>::value` prints 0
(`et_off`); `static_assert(std::same_as<decltype(exp(x)), cpp_bin_float_50>)` passes;
`std::same_as<decltype(exp(y)), number<cpp_bin_float<50>, et_on>>` is false; `abs(*r - e)` fails to compile for both
types ("no match for operator-"); `*r < e` and `*r == e` compile for both.

Suggested change. The comment: "-> mp50: a number, even for a type with expression templates (cpp_bin_float_50 has
none; number<..., et_on> has)". The text: "mp50(*r) in the arithmetic is needed for every Boost type; -> mp50 on the
lambda matters only for a type with expression templates, and is a safe habit". The sidebar's "may return an
expression" can stay.

### R1-library-05 (minor): "derivative never throws"

Location: `chapters/07-using.tex` 36.

Problem. The library throws nothing, but `derivative` is not `noexcept` and passes on whatever `fn` or `T` throws
(`deriv.hpp` 1380; the book says so correctly in 8.2 and 10.5). "derivative never throws" invites wrapping it in a
`noexcept` function, which terminates the program when fn throws.

Suggested change. "derivative itself never throws: a failure is a value. An exception thrown by fn passes through
unchanged (Section 8.2)."

### R1-library-06 (minor): the weight check is a static_assert in the body, not a constraint

Location: `chapters/08-api.tex` 163-167 (8.5) and 237-242 (8.8, "What the constraints reject does not compile: the
error is reported at the call"); `chapters/10-guarantees.tex` 131-132.

Problem. The check that a rule's integer weights are exact in T (`deriv.hpp` 1009) is a `static_assert` inside
`stencil_derivative`, reached only when the body is instantiated. A requires-expression, the tool 8.8 itself uses,
reports such a call as valid, and generic code that selects a method with `if constexpr (requires { ... })` then fails
to compile inside the library. Part II should say which rejections are constraints and which are hard errors.

Evidence (`half.cpp`, GCC 16): `nxx::real<std::float16_t>` and `nxx::real<std::bfloat16_t>` hold; a concept
`requires(X x) { derivative<M>(fn, x); }` is true for `Central<10>` with `std::float16_t`, while the call itself fails
with "this rule's integer weights are not exact in T" at `deriv.hpp` 1009 (`Central<8>` with `std::bfloat16_t`
likewise).

Suggested change. 8.5: "This check is a static_assert in the body, not a constraint: a requires-expression still
reports the call as valid." 8.8: "The constraints reject the cases below; a type too narrow for a rule's weights is
rejected by a static_assert instead (Section 8.5)."

### R1-library-07 (minor): Part II's account of Ridders' outcomes is incomplete

Location: `chapters/08-api.tex` 98-104 (8.2, Failures).

Problem. 8.2 says an untrusted result is returned "if the levels run out during restarts after an untrusted result".
The code returns it in two more cases: when a tableau ends at the resolution of x, its result is returned even if
untrusted (`deriv.hpp` 1314, `if (!untrusted || resolution) return result;`), and when a later tableau ends at the
resolution of x before two successful levels, the earlier untrusted result is returned (1299 `break`, then 1318
`if (fallback) return *fallback;`). The first case is tested: the Heaviside step at 1 with h = 2^-51 returns a value
with `abs_error >= |value| / 4` after 2 levels (`tests/test_ridders.cpp` 296-300). Sidebar 6.2 states the resolution
case; Part II, the reference, does not. A user who reads 8.2 would not expect a value whose estimate exceeds a quarter
of it.

Suggested change. Add to 8.2: "A value is also returned, with its large abs_error, when the step reaches the resolution
of x (no restart can go lower). An abs_error of a quarter of |value| or more marks a result the method itself does not
trust."

### R1-library-08 (minor): what the tests check is overstated in three places

Location: `chapters/04-stencils.tex` 177-178; `chapters/10-guarantees.tex` 214-215 (Table 10.5, `test_stencils.cpp`)
and 183-184.

Problem and evidence.
- 4.3: "the tests budget (C + 2S) eps^(P/(P+1)) for every multi-node rule". `tests/test_stencils.cpp` 212-221 applies
  it to ten rules for f', among them the two-node `Central<2>` and `Forward<1>`, and not to `Central<8>`,
  `Forward<3>`, `Forward<5>`, `Backward<1>` or `Backward<3>` to `Backward<5>`; `test_orders.cpp` 614-627 extends it
  to fourteen (rule, N) pairs. Table 10.5 says "every stencil".
- 10.5: `root_epsilon` "verified against references computed to 80 digits, for float and double and N from 2 to
  2000". `tests/test_math.cpp` 47-66 checks selected N: 2, 3, 5, 7, 11, 42, 297 and 2000 for double, 3, 5, 7, 11, 87,
  128 and 146 for float, plus the multiprecision types against 60-digit references. The range 2 to 5000 at 0.58 ulp
  is a measurement in the header's comment (`math.hpp` 107-108), not a test.

Suggested change. 4.3: "the tests budget (C + 2S) eps^(P/(P+N)) for a selection of rules and orders, two-node rules
included". Table 10.5: "stencils of each family with their default step against ...". 10.5: "for selected N from 2
to 2000 (float and double, against 80-digit references) and for the multiprecision types".

### R1-library-09 (minor): the cost of misjudging the scale is stated differently in Chapters 3 and 5

Location: `chapters/05-orders.tex` 104-108 and 141-142; `chapters/03-rounding.tex` 277-278.

Problem. Chapter 3: "a factor of two in w costs at most a factor of four in the error of Central<2>". Sidebar 5.1:
"misjudging that scale by a factor 2 costs up to 2^N in error: a factor 2 for f' with Central<2>". The sidebar's own
previous sentence gives the other side: a step twice too large multiplies the truncation term by 2^P, which for
`Central<2>` and f' is 4, and for `Central<4>` and f' is 16. The term-wise bound is max(2^P, 2^N); the header's "up
to 16" (`deriv.hpp` 1415) is the f'''' case of `Central<2>`, where 2^N = 16 exceeds 2^P = 4.

Suggested change. Sidebar 5.1 and 5.3: "costs up to max(2^P, 2^N) in one of the two terms: 4 for f' with Central<2>,
16 for f'''' (the figure the header states)". Chapter 3 stays.

### R1-library-10 (minor): Part II omits two properties of Estimate that the book's own listings rely on

Location: `chapters/08-api.tex` 135-140 (8.3).

Problem.
- Comparing an `Estimate<T>` (or the `std::expected` holding one) with a T compares the value only, through the
  conversion: the defaulted `==` takes two Estimates, and a T does not convert to an Estimate. Section 7.7's listing
  depends on it (`static_assert(derivative(cube, 2.0, 0.25) == 12.0625)`, `code/07-using-constexpr.cpp` 22), as do the
  tests (`one == 1.0`, `tests/test_ridders.cpp` 168). 8.3 says only that "without [the defaulted ==], == would compare
  only the values", which suggests that `==` never compares the value alone.
- `abs_error` can be `std::numeric_limits<T>::max()`: the rounding floors and the change from the parent saturate
  there (`deriv.hpp` 589, 603-608). The tests return `{0, max}` for a constant over a tiny step
  (`tests/test_orders.cpp` 331-338; `test_ridders.cpp` 139). A reader should know that such an estimate carries no
  information.

Suggested change. 8.3: "Comparing with a T (`*r == 2.0`, or the expected with a T) compares the value only, through
the conversion. abs_error saturates at max() where the error estimate would overflow (a large fn over a tiny step);
such a result says only that the method could not estimate its error."

### R1-library-11 (editorial): "Every accepted form" is a selection

Location: `chapters/08-api.tex` 116; `code/08-api-usage.cpp` 1-3.

Problem. The listing shows eight calls; it does not show Ridders with an initial step, N = 3 or 4, an integer type
other than int for N, or a member-pointer callable, all of which are accepted.

Suggested change. "One call of each form" or "the common forms".

### R1-library-12 (editorial): the f'' rule from the actual offsets is not bit-identical to the textbook form at exact points

Location: `chapters/05-orders.tex` 40-42 ("the value is the same whenever the points are exact").

Problem. With exact points the two forms are equal in exact arithmetic, but they round differently unless h is a
power of two. Evidence (`second.cpp`, e^x at x = 1 + k/256, k = 0 to 199, `-ffp-contract=off`): with h = 2^-10 none of
200 results differ from `((f(x - h) - 2 f(x)) + f(x + h)) / h / h`; with h = 3 * 2^-12 (still exact points) 131 of 200
differ, at x = 1 by 1.5e-13 relative.

Suggested change. "the two forms agree in exact arithmetic whenever the points are exact".

### R1-library-13 (editorial): the node-count margin table describes Ridders incompletely

Location: `chapters/04-stencils.tex` 74 (Table 4.1).

Problem. The Ridders row gives its nodes as "+-h_0 / 2^i". For f''' and f'''' a level also samples +-2 h_i, and for
even N x itself once (`deriv.hpp` 1196-1200, 1229). The column also changes meaning, from a count for the stencils to
a set for Ridders.

Suggested change. "+-h_i (N <= 2), +-h_i and +-2h_i (N >= 3); x once for even N", or a pointer to Table 8.2, which is
exact.

### R1-library-14 (editorial): "the factor 1 instead of 4" is attached to the wrong constant

Location: `chapters/06-ridders.tex` 279-284 (Sidebar 6.2).

Problem. The list of constants names "the factor 1 on the floor" and then, as evidence, "the factor 1 instead of 4
left poles between samples uncovered". The latter is about the untrusted test's factor 4 on |value| (`deriv.hpp`
1310; `docs/PHASE3_API.md` 164), not about the floor; the reader will attach it to the floor.

Suggested change. "...; an untrusted test at |A| (factor 1) instead of |A| / 4 left poles between samples
uncovered; ...".

### R1-library-15 (editorial): Part II's declarations are typed, not quoted

Location: `chapters/08-api.tex` 47-58, 121-129, 144-147, 156-161, 213-220 (the `signature` blocks).

Problem. CLAUDE.md: "No listing is typed into a chapter." The declarations of 8.2 to 8.7 are typed. I compared each
with the header and all are exact today (`deriv.hpp` 41-47, 66-83, 92-97, 216-274, 1400-1475; `math.hpp` 30-39,
49-78), but nothing keeps them so; quoting from `include/` would need markers in a public header.

Suggested change. Keep the blocks, and add to `tools/check.py` a check that each code line of a `signature` block
occurs, whitespace-normalised, in `deriv.hpp` or `math.hpp`; or record the exception to the listing rule in PLAN.md
section 8.

### R1-library-16 (editorial): style decision 18 is not followed in ten quoted snippets

Location: unbraced `if`s in quoted snippets: `01-why-first.cpp` 25 and 31-32 (a two-line statement);
`02-taylor-steps.cpp` 28; `03-rounding-step.cpp` 24-25 (two lines); `04-stencils-edge.cpp` 32-33 (two lines);
`07-using-newton.cpp` 33 and 37; `07-using-solver.cpp` 33-36 (an if/else over four lines); `07-using-expected.cpp` 90;
`07-using-check.cpp` 49 and 54; `07-using-multiprecision.cpp` 27. (The unquoted snippets `constexpr-failure`,
`check-absolute` and `jacobian` have more.)

Problem. PLAN.md section 8, item 18: "Every if in a snippet is braced". None of them warns: all 20 examples build
without a diagnostic under GCC 16 and Clang 22 with `-Werror`, so this is style, not a build risk. Several snippets
also `return 1;` from a `main` the reader never sees (`forward-vs-central`, `check-by-steps`, `multiprecision`).

Suggested change. Brace them, or narrow the rule to multi-line `if`s; in the listings, print and stop rather than
`return 1;`.

### R1-library-17 (editorial): Appendix B says its cost figures are the header's, rounded

Location: `chapters/a2-derivations.tex` 165.

Problem. "These are the figures the header and the design record state, rounded." For `Central<10>` and f' the book
gives `\valStepCostCentralTenFirst` = 1.63, the header and `docs/design.md` 1.64 (`deriv.hpp` 521; PLAN.md D-16 and
U-4); 1.634 does not round to 1.64.

Suggested change. "These agree with the figures of the header and the design record (3.5, 11 and 15), except that
both give 1.64 for Central<10>, where the model gives 1.63."

### R1-library-18 (editorial): Table 8.1 separates the declarations from their lead-in

Location: PDF pp. 58-59 (book pp. 51-52); `chapters/08-api.tex` 21-46.

Problem. The float `tab:api-names` is set after "Declarations. Two overloads, with a step and with the default
step:", so the colon is followed by the table of public names, and the signatures start on the next page.

Suggested change. Place the table at the top of the section, or before 8.2's heading, or end the lead-in without a
colon.

### R1-library-19 (editorial): why exact types are rejected

Location: `chapters/07-using.tex` 246-247 (Sidebar 7.1).

Problem. "Exact types such as rationals are rejected, since their epsilon is 0." The concept rejects them by
`is_exact` (and `cpp_rational` also by `is_bounded`), `math.hpp` 31-32; epsilon is not tested. The test comment gives
both properties (`tests/test_math.cpp` 20).

Suggested change. "...are rejected (numeric_limits marks them exact): their epsilon is 0 and they have no largest
value, so neither the step nor the overflow checks would mean anything."

### R1-library-20 (editorial): the Ridders constants table

Location: `chapters/09-tables.tex` 145-149 and Table 9.10 (`data/table-ridders-constants.tex`;
`tools/generate_data.cpp` 1386-1391).

Problem. The row "Untrusted result: estimate at least |A| / 4 and above the floor by: 4, 64" repeats the 4 and reads
as "floor + 64"; the condition is `abs_error / 64 > floor`, more than 64 times the floor. The rows for the ratio 2,
the stall of 2 and the factors 4 and 64 are typed in the generator (they are literals in `detail::ridders`, with no
named constant), while 9.1 says every body was written "from the library itself".

Suggested change. Row: "Untrusted result: estimate at least |A| / 4 and more than 64 times the floor". 9.5: "the cap
and the initial steps are read from the library; the other constants are literals of its code, transcribed".

### R1-library-21 (editorial): "Central<2>, 2" in Table 7.1

Location: `chapters/07-using.tex` 97 (Table 7.1, row "A second derivative").

Problem. `\code{Central<2>, 2}` and `\code{Ridders, 2}` are not code anyone writes; the reader has to infer the
template argument list.

Suggested change. `derivative<Central<2>, 2>` and `derivative<Ridders, 2>`, or "Central<2> or Ridders, N = 2".

### R1-library-22 (editorial): `int n = *r` is a compile error under the book's own flags

Location: `chapters/07-using.tex` 59; `chapters/08-api.tex` 139.

Problem. The text gives the Clang caveat for `if (*r)` (style decision 17) but not for `int n = *r`, which with
`-Wconversion -Werror` fails on both compilers (`cx1.cpp -DINTN`: GCC "conversion from 'double' to 'int' may change
value [-Werror=float-conversion]", Clang "implicit conversion turns floating-point number into integer"). GCC accepts
`if (*r)` under the same flags; only Clang rejects it, as the book says.

Suggested change. "And int n = *r truncates (a -Wfloat-conversion warning on GCC and Clang)."

## Checked and correct

Names and declarations
- Part II lists every public name: `derivative` (two overloads, each with a deprecated twin), `Estimate<T>` (members,
  `operator T`, defaulted friend `==`), `error::non_finite`, `order` (`: int`, first to fourth), `Central<P>`,
  `Forward<P>`, `Backward<P>`, `Ridders`, `NUMERIXX_NO_PRECISION_WARNING`, `nxx::real`, `nxx::math::abs`, `isfinite`,
  `root_epsilon`. The signature blocks match the header token for token (finding 15 is about keeping them so).
- Template parameters: `M` defaulted to `Central<2>`; `auto N = 1`; `real T` deduced from x alone, the step
  `std::type_identity_t<T> const &`; `detail::method` closed; `detail::valid_order` accepts any integer type but bool,
  and `order`, compared in its own type (`0x1'0000'0001ULL` rejected, `tests/test_orders.cpp` 199);
  `detail::callable` is `std::invocable<F &, T const &>` with a result convertible to T, called through `std::invoke`
  with a T const lvalue (`test_callables.cpp`).
- Valid orders: `Central<P>` even 2 to 10; `Forward<P>`, `Backward<P>` 1 to 6; the static_assert texts quoted in 8.5
  are the header's; an out-of-range tag fails the constraint, not a hard error (`test_stencils.cpp` 108-113; E16).
- `derivative<Ridders, double>(fn, x)` does not compile (GCC and Clang: "invalid explicitly-specified argument for
  template parameter 'N'").
- `real<T>`: the six conditions and the operator requirements; the types the tests accept and reject (8.7) are those
  of `tests/test_math.cpp` 17-21. `abs`, `isfinite` (`x == x && abs(x) <= max()`), `root_epsilon` (static_assert
  N >= 2, scaling into (radix^-N, 1], Newton from 1) as described. CMake target `numerixx::numerixx`, `cxx_std_23`.

Results, failures, costs
- `Estimate`: aggregate, structured bindings, implicit conversion; not considered by `auto`, `std::format` and Boost's
  arithmetic operators, while Boost's comparisons accept it (`boost1.cpp`, et_off and et_on); `if (*r)` tests the
  value (false for cos at 0 by Ridders, whose value is exactly 0, E11) and is a Clang `-Wfloat-conversion` error, not
  a GCC one; `d.value_or(0.0)` does not compile on GCC or Clang ("argument has to be convertible to value_type");
  `abs_error` empty for every stencil, always filled by Ridders.
- The failure list of 8.2 is the header's (`deriv.hpp` 43-44); fn(x) non-finite for even N costs Ridders one
  evaluation (`test_orders.cpp` 499-505); a NaN or infinite x gives non_finite without a call, with or without a step
  (`claims1.cpp`: "NaN x with h: false false false calls 0").
- Evaluations: `Central<P>` P + N - 1, one-sided P + N; Ridders 2 or 4 per level plus fn(x) once for even N, at most
  64, 65, 128 and 129 (cap 32, `deriv.hpp` 1085); `sin x / x` at 0 fails for f'' and works for f' and f''' (central and
  Ridders), as 5.1's margin note says.
- Default steps: eps^(1/(P+N)) max(|x|, 1) via `root_epsilon<P + N, T>()`; Ridders s/32 for f', s/4 for f'' to f''''
  (`ridders_divisors`); the counts and exponents of the default-step table.
- Ridders: the tableau in Burden's difference form, two arrays of 32 filled from T(0); the estimate max(change, floor),
  raised to the distance from the latest diagonal entry; stops on exact agreement, after a two-level stall, at the
  resolution of x, at the cap; restarts from half the failing step, and from half the initial step for an untrusted
  result (|A| / 4, 64 times the floor); the floors of f', f'' and f''' / f'''' including the point-rounding term.
  `log x` at 0.01: 24 evaluations, two restarts; `sin(64 pi x)` at 0: exactly 64 pi after 24 evaluations;
  `cos(64 pi x)` for f'' at 0: 0 with the floor as its estimate after 5 evaluations, the uncaught case (`claims1.cpp`).
- The two-node rules divide by the evaluated distance, `Central<2>` for f'' uses the actual offsets, the other rules
  divide by d h and then by h, never by h^N; both sides are checked for the multi-node rules; neighbouring points must
  differ.

Warning, constexpr, types, exceptions
- The precision warning: the criterion radix^(digits-1) < 10^7 (float, `std::float16_t` and `std::bfloat16_t` noisy,
  double not); any method, N = 4 only; reported at the caller's line on GCC and Clang, also when the header is
  included as a system header; the message quoted in 8.6 is the header's; `-Wdeprecated-declarations` / C4996; the
  macro turns it off (`warn.cpp`).
- Constant evaluation: E15's static_asserts hold; a failure is a value at compile time; a fourth derivative with a
  tiny step is non_finite at run time, a GCC compile error ("'(-2.0e+240 / 9.9999999999999998e-121)' is not a
  constant expression") and non_finite under Clang (`cx1.cpp -DOVF`), as Sidebar 7.2 says; compile-time and run-time
  results agree bit for bit in E15 on GCC; the guarantee carries its hypothesis (style decision 15).
- Exceptions and memory: no `throw`, `try` or `catch` in the headers; no allocation in the library.
- Test claims: the moment conditions for every rule and order and the published rows of Fornberg and Chapra
  (`test_stencils.cpp`, `test_orders.cpp`); `Central<4>` on x^4 exactly 4 at compile and at run time; the
  static_asserts on overflowing, non-finite and saturating cases; exactness on polynomials with the bit budget
  (`fits_53_bits`); coverage within 6 times the estimate and the estimate within twice the prototype's; the types
  `derivative` is tested on; three build tests for the warning and one doctest executable.
- Deferred features: each item of 7.9 and 10.7 is in `docs/plan.md`, a PHASE document or `DERIV_TECHNIQUES.md`.

Listings
- All 36 `\snippet` and `\snippetheader` uses name existing markers; the PDF is newer than every source, so each
  listing shown is the compiled text.
- All 20 examples build with no diagnostic under GCC 16 and Clang 22 with the book's warning flags and `-Werror`, and
  all exit 0. Their printed numbers agree with `data/values.tex` (for example cp = 32.6917299992, and Ridders' 2.7e-12
  after 14 evaluations, in E1; 4 Newton steps to 1.1e-13 K; the Jacobian's 4 + 16 evaluations; sensitivity errors of
  3.2e-12 and 2.6e-14, below their estimates 3.9e-12 and 8.3e-14; 7.2e-45, and 22 evaluations for Ridders, with
  cpp_bin_float_50).

## Verdict

The book describes `nxx::deriv` faithfully. I found nothing that would lead a user to wrong code in normal use: every
name, parameter, constraint, default, cost and constant in Part II matches the header, the examples compile without a
diagnostic on two compilers and check what the text claims, and the claims about the tests are, apart from finding 8,
what the tests do. The ten minor findings should be fixed before release, most with a sentence or two: edge cases
stated more strongly than the code supports (01, 03, 05), one listing that contradicts the idiom the chapter teaches
(02), one Boost fact (04), Part II's coverage of what a requires-expression, a comparison with a T and Ridders'
outcomes actually do (06, 07, 10), and two places where Part I overstates the tests or disagrees with itself (08,
09). The editorial findings are for the editor to decide. With those changes Part II is complete and exact, and I
expect to sign off on the library content in round 2.

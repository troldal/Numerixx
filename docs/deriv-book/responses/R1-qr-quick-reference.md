# R1-qr-quick-reference: response to review round 1

Author: `book-author-practice` (Opus 5.5), stage 6. File: `chapters/qr-quick-reference.tex`. `tools/check_chapter.sh`:
`COMPILED OK`, one page, with room left; "Where to look" did not need to go. The method-choice rows are the editor's
shared `chapters/shared-method-summary.tex`, unchanged.

| Finding | Decision | What changed |
|---|---|---|
| T1 (engineer-01, developer-02) | accept | New block "Two checks": "Did it work? if (d). Is it right? A wrong step scale, a noisy f or a break near x return a wrong value without failing: if in doubt, repeat at 2h and h/2. For Ridders, also require `*r->abs_error <= tol * std::abs(r->value)`" |
| T11 (developer-03) | accept | "The result" names `d.value()` as std::expected's accessor that throws on a failure ("test, then read d->value"); the last pitfall is "auto and .value()" |
| developer-04 | accept | Failures: "A step <= 0 is a bug (an assert, gone under NDEBUG)" |
| engineer-15 | accept | New block "Your own step": `h = nxx::math::root_epsilon<3, double>() * w;` for eps^(1/3) w, and `h = std::cbrt(eta) * std::max(std::abs(x), 1.0);` for a noisy f |
| evidence-23 | accept | Ridders' caps by `\valRiddersMaxEvals{First..Fourth}` |
| T12 (numerics-17k, evidence-13) | accept | The Ridders line: "Near a bad point it returns a value down to 2^-`\valRiddersEdgeLimitExponent` max(|x|, 1) away, ten digits only down to 2^-`\valRiddersEdgeTenDigitsExponent` max(|x|, 1) (for f')" |
| T13 | accept | Rounded headline digits by macro: Central<2> about `\valHeadlineFirstDigits` digits for f' and `\valHeadlineFourthDigits` for f'''' (11, 5); Ridders about `\valHeadlineRiddersFirstDigits` (14) |

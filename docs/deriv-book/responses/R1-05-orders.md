# R1-05-orders: response to review round 1

Author: `book-author-methods` (Fable 5.1), stage 6. File: `chapters/05-orders.tex`; also `code/05-orders-float.cpp`.
Compiled with `tools/check_chapter.sh` (`COMPILED OK`, 5 pages standalone, with the last page now half empty: the
pragma listing has gone to Chapter 8); `tools/check.py` reports no warning in the file; E8 rebuilt and run (exit 0).

| Finding | Decision | What changed |
|---|---|---|
| T9 (numerics-06, library-09) | accept | Sidebar 5.1, last paragraph: "a step a factor 2 below the right one multiplies the rounding term by 2^N, and a step a factor 2 above it multiplies the truncation term by 2^P ... at the model's balance misjudging that scale by a factor 2 grows the total error by at most `\valMisjudgeCostCentralTwoFirst` for f' and `\valMisjudgeCostCentralTwoFourth` for f'''' with Central<2>, and `\valMisjudgeCostCentralFourFourth` for Central<4> and f''''; the factor 16 that the library's documentation quotes for f'''' is the growth of the rounding term alone, 2^N". 5.3 quotes the same two Central<2> macros and names 2^N as the rounding term's growth; the former "costs up to 2^N in error: a factor 2 for f' ... 16 for f''''" is gone from both places, so 3.5, Sidebar 5.1 and 5.3 agree by construction |
| T13 (numerics-17h, evidence-22) | modified | 5.2: "roughly eleven, eight, six and five digits" is now the generated `\valModelDigitsCentralTwo{First,Second,Third,Fourth}Double` (10.4, 7.8, 6.3, 5.2); the headline sentence quotes the rounded `\valHeadlineFourthDigits` and `\valHeadlineFirstDigits` (5 and 11). Figure 5.1's caption: "gives about `\valHeadlineFourthDigits` digits where the first gives `\valHeadlineFirstDigits`". Sidebar 5.1: "the difference between about `\valModelDigitsCentralTwoFirstDouble` digits and about `\valModelDigitsCentralTwoFourthDouble`" |
| engineer-17 | modified (approved) | 5.4 keeps what the warning means, its criterion and the cure ("The cure is to use double") in a few sentences, says that the warning can be silenced per call or per build and that the result is the same, and points to 8.6 for the quoted call and its numbers. The `\snippet{05-orders-float}{float-fourth}` is removed from this chapter; the practice author places it in 8.6 |
| developer-16 | modified | `code/05-orders-float.cpp`: the function is now e^x at `1.0f`, the call the generator makes for the `\valFloatFourth...AtOne` macros of DS6; it prints the relative errors against e and Ridders' estimate, absolute and relative. Output: Central<2> 6.7e-03, Ridders 9.7e-05 with estimate 7.7e-01 (2.8e-01 relative), which are `\valFloatFourthCentralAtOneRelError`, `\valFloatFourthRiddersAtOneRelError`, `\valFloatFourthRiddersAtOneEstimate` and `\valFloatFourthRiddersAtOneRelEstimate`. The file's header comment says so |
| engineer-19 | modified | 5.3: "as the header warns" removed (the sentence ends at "for every order"); Sidebar 5.1: "the figure the header states" is now "the factor 16 that the library's documentation quotes" |
| library-12 | accept | 5.1: "the two forms agree in exact arithmetic whenever the points are exact" |
| evidence-21 | accept | Sidebar 5.2 names `docs/PHASE4_API.md` at the 8.7 against 1.4 figures ("in the prototype of the design record"); 5.3 attributes the half-step result to "the design record's prototype (`docs/PHASE4_API.md`)" |
| engineer-20 | accept | 5.6: "The example checks that the two results compare equal ..." deleted; the lead-in to the `order-enum` snippet says that the two calls return equal results, value and estimate alike |
| T5 (engineer-06) | modified | Checked: 5.2, 5.3 and Sidebar 5.1 say "the default step assumes that f varies on the scale s", which is the assumption as worded in the decision; no change needed |
| numerics-17f (via Sidebar 5.2) | accept | The Sterbenz sentence now says "lie within a factor 2 of each other and have the same sign", with `\cite[Thm~4.3.1]{sterbenz1974}` |
| Citations | accept | `\cite[ch.~23]{chapra2015}` at Chapra's tables |

Not changed, and why: 5.5's "a median of `\valRiddersEvalsMedianFourth` evaluations over the random study" stays as the
study's figure (not the e^x points'), as the sentence says.

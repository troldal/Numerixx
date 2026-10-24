# R1-10-guarantees: response to review round 1

Author: `book-author-practice` (Opus 5.5), stage 6. File: `chapters/10-guarantees.tex`. `tools/check_chapter.sh`:
`COMPILED OK`, 5 full pages standalone (4.6 before, +0.4 against the plan's +0.3). `tools/check.py`: no warning.

| Finding | Decision | What changed |
|---|---|---|
| T12 (numerics-05, evidence-11) | accept | 10.1: the "factor of ten" is gone. "Against the double column, the central rules lie below the model by up to a factor of `\valModelRatioCentralBelowMax` for f' to f''' and up to `\valModelRatioCentralFourthAbove` times above it for f''''; Forward<2> lies up to `\valModelRatioForwardTwoAbove` times above it, and Forward<4> `\valModelRatioForwardFourThird` and `\valModelRatioForwardFourFourth` times above it for f''' and f'''', by its larger amplification S, which the model drops (measured)" |
| T13 | accept | 10.1's headline digits are the rounded macros (now 11 and 5; Ridders 14 and 9) |
| T12 (numerics-09, evidence-12) | accept | 10.2: "For f' no singularity lies within the initial step; for f'' to f'''' the initial step s/4 (and the span 2 h0 of f''' and f'''') reaches past zero for x below 0.25 (0.5), which are the restarts of Table 10.2" |
| numerics-17e | accept | Table 10.2 prints the generator's two-decimal coverage for f''' (99.96) |
| T12 (numerics-17j, evidence-09) | accept | 10.3: Ridders' `\valRiddersEvalsMedianFirst` evaluations against a stencil's 2 to 6 buy about `\valHeadlineRiddersFirstDigits` digits against `\valHeadlineFirstDigits` for Central<2>, as many as Central<4> and none over Central<6> (Table 10.1), and an estimate |
| T10 (numerics-07, evidence-05) | modified | 10.4 Orders: "The caps on P are design decisions (docs/PHASE2_API.md)", Central<12>'s gains by the `\valCapGainCentralTwelve...` macros, and the one-sided amplification "grows about as fast as a higher order pays"; "under half a digit" is gone (U-7) |
| T12 (numerics-17k, evidence-13) | accept | 10.4 Domain edges: a value down to k = `\valRiddersEdgeLimitExponent` "with an estimate grown to flag the loss, but ten correct digits only down to k = `\valRiddersEdgeTenDigitsExponent`"; "For ten digits, keep the nearest bad point at least 2^-`\valRiddersEdgeTenDigitsExponent` max(|x|, 1) away"; the design record's derivation of 2^-34 and 2^-31 stays, the second named as not measured |
| library-03 | accept | 10.5 "Where fn is called": "for a step h > 0 (the precondition, an assert that NDEBUG removes), Forward<P> never below x and Backward<P> never above" |
| evidence-23 | accept | 10.5's caps are the `\valRiddersMaxEvals...` macros |
| developer-13 | accept | 10.5, new guarantee "No shared state": "The headers hold no global, static or thread-local variable other than compile-time constants, and fn is called on the calling thread, one call at a time. Concurrent calls of derivative are therefore as safe as the concurrent calls of fn they make. Verified by inspection" (the only namespace-scope variables in the headers are `inline constexpr` constants) |
| library-08 | accept | 10.5: root_epsilon "verified against references computed to 80 digits for selected N from 2 to 2000 in float and double, and for the multiprecision types against 60-digit references". Table 10.5: "stencils of each family with their default step against (C + 2S) eps^(P/(P+1))" |
| T6 (evidence-15, library-08) | modified | 10.6 takes the test budget from 4.3: "The stencils, for example, are held to (C + 2S) eps^(P/(P+N)) relative on e^x, truncation plus the rounding of the samples and of the points, for a selection of rules and orders, two-node rules included; the tests' comment calls it an estimate, not a bound, and reports measured errors 4 to 16 times below it on its three points" (attributed, not generated) |
| evidence-20 | accept | 10.6: "Each phase was reviewed by specialised review agents, Phase 4 also by another vendor's model (OpenAI Codex); the findings and decisions are recorded in the design record ..." |
| evidence-25 (citation form) | accept | 10.6 cites "the published tables of Fornberg\cite{fornberg1988} and Chapra\cite[ch.~23]{chapra2015}" |

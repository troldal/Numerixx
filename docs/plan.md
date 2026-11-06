# Plan

## Milestone 1: numerical derivatives of a scalar function (`nxx::deriv`)

One phase per step; the user picks the next phase and approves its API before it is written. A box is ticked when the
phase's `static_assert` and run-time tests pass on GCC, Clang, MSVC, clang-cl and Emscripten, locally and in CI.

### Phase 1: central difference

- [x] `derivative(fn, x, h)` and `derivative(fn, x)` with the default step, constexpr (29ba464)
- [x] Tests on `float`, `double` and Boost.Multiprecision (`cpp_bin_float_50` with and without expression templates,
  `cpp_dec_float_50`) (29ba464)
- [x] Non-finite results return `error::non_finite` (29ba464)
- [x] CI passes on all five toolchains (run 39437081443)

### Phase 2: stencils

API approved: `docs/deriv/PHASE2_API.md`.

- [x] Higher-order central stencils, `Central<P>` for even P up to 10 (56cffb1)
- [x] Forward and backward stencils, `Forward<P>` and `Backward<P>` for P up to 6, for points where fn can be evaluated
  on one side only (56cffb1)
- [x] `root_epsilon<N, T>` replaces `cbrt_epsilon` (56cffb1)
- [x] CI passes on all five toolchains (run 39486718887)

### Phase 3: Ridders' extrapolation

API approved: `docs/deriv/PHASE3_API.md`.

- [x] `derivative<Ridders>` with an error estimate, restarting below failing or untrusted steps (1023248)
- [x] One return type for every method, `std::expected<Estimate<T>, error>`, implicitly convertible to T (1023248)
- [x] CI passes on all five toolchains (run 39605999468)

### Phase 4: derivatives of order 1 to 4

API approved: `docs/deriv/PHASE4_API.md`.

- [x] `derivative<M, N>` for N = 1 to 4, as an integer or an `enum class order`, with every method (d39e31b)
- [x] One Lagrange-formula generator for every table; Phase 2 tables reproduced bit for bit (d39e31b)
- [x] Ridders for N = 1 to 4 on the `Central<2>` rule of each order (d39e31b)
- [x] No overflow before the final quotient (fixes a constant-evaluation failure of `Central<4>` on a constant 1e308)
  (d39e31b)
- [x] A compile-time warning for a float 4th derivative, silenceable at the call or by `NUMERIXX_NO_PRECISION_WARNING`
  (d39e31b)
- [x] CI passes on all five toolchains (run 39910979076)

### Phase 5: one-sided Ridders, and the float warning from f''

API approved: `docs/deriv/PHASE5_API.md`.

- [x] `derivative<ForwardRidders, N>` and `derivative<BackwardRidders, N>` for N = 1 to 4, on a `Forward<1>` base with
  the rounding floor carried through the tableau, and `BackwardRidders` as the forward method on the mirrored function
  (f20696f)
- [x] One extrapolation driver on the base stencil; `Ridders` unchanged bit for bit except after untrusted results
  (f20696f; 28,000 cases compared with the previous header)
- [x] After an untrusted result every Ridders method restarts below the smallest step tried (fixes poles 1e-3 and 1e-4
  from x) (f20696f)
- [x] The low-precision warning for every derivative above the first in float, every method (f20696f)
- [x] An independent Python prototype that reproduces the C++ results for double bit for bit, as the source of the
  test references (`.scratch/onesided/ridders5.py`, 28,000 cases)
- [x] CI passes on all five toolchains (run 40138064893)

## Milestone 2: roots of a real function of one variable (`nxx::roots`)

Architecture approved: `docs/roots/BRIEF.md` (2026-10-25). Each phase gets its API approved before it is written
(`docs/roots/PHASE<n>_API.md`), as in Milestone 1.

- Phase 1: bracketing core, `Bracket`, `Root`, `error`, `Criteria`, `solve<M>(fn, Bracket[, criteria])`; `Bisection`
  and `Brent`. API approved: `docs/roots/PHASE1_API.md`.
  - [x] `roots.hpp` with `solve<Bisection>` and `solve<Brent>`, constexpr, the validation, zero first and the defaults
    of the specification; `nxx::detail::callable` shared from `math.hpp`
  - [x] Tests on float, double, long double and Boost.Multiprecision (`cpp_bin_float_50` with and without expression
    templates, `cpp_dec_float_50`), with the prototype's evaluation counts reproduced exactly in double
  - [x] CI passes on all five toolchains (02fffe5: run 40192833145)
  - [x] The independent review (Codex) and its adopted findings: five, all reproduced and fixed
    (`docs/roots/PHASE1_API.md`)
  - [x] A demo, `numerixx_roots_demo`: Bisection against Brent, Kepler's equation, tolerances, a bracket around 0,
    types, compile time, failures and composition (added after the review, at the user's request)
- Phase 2: open methods and derivatives (fn returning a tuple-like, or separate derivative callables); `Newton` and
  `Halley` on a guess or a bracket, `Secant`; the entry point renamed `solve`. API approved:
  `docs/roots/PHASE2_API.md`.
  - [x] Newton, Halley and Secant from a guess with the hand-over, Newton and Halley on a bracket as Brent's loop with
    a derivative step, both derivative forms, `not_converged`; Phase 1's results unchanged
  - [x] Tests in double (the prototype's counts exactly), float, long double and the Boost.Multiprecision types; the
    constraints; a demo section
  - [x] CI passes on all five toolchains (9c00c98 with 8b3ff81: run 40477325328)
  - [x] The independent review (Codex) and its adopted findings: eight, all reproduced and fixed
    (`docs/roots/PHASE2_API.md`)
- Phase 3: numerical derivatives, `nxx::deriv::derivative_of<M, N>(fn)`, and derivative callables returning
  `std::expected`. API approved: `docs/roots/PHASE3_API.md`.
  - [x] `derivative_of` (a decayed copy of fn, called with the callable's constness, no precision warning); `dfn` and
    `d2fn` returning `std::expected`, a failed derivative counting as one that is not finite
  - [x] Tests: `derivative_of` against `derivative` bit for bit; the equivalences and the partial failure; the
    prototype's `Central<2>` counts exactly; the types; constant evaluation; the constraints; a demo section (five
    toolchains and WSL GCC and Clang, locally)
  - [x] CI passes on all five toolchains (2292381: run 40499775174)
  - [x] The independent review (Codex, and Kimi K3.1 Pro through Cline) and its adopted findings: six, all
    reproduced and fixed (`docs/roots/PHASE3_API.md`)
- Phase 4: `Toms748`, `AndersonBjorck`. API approved: `docs/roots/PHASE4_API.md`.
  - [x] `Toms748` (Algorithm 4.2 as printed) and `AndersonBjorck` (with the bisection guard), sharing
    `detail::enclose` (the paper's modified `bracket`); every operation on fn's values checked
  - [x] Tests: the prototype's counts exactly; accuracy; the placement fallbacks; constant evaluation; the types; a
    demo section (five toolchains and WSL GCC and Clang, locally)
  - [x] CI passes on all five toolchains (8951e7c: run 40511081791)
  - [x] The independent reviews (Codex, and Kimi K3.1 Pro through Cline): no defect in the header; two test
    findings fixed (`docs/roots/PHASE4_API.md`)
- Phase 5: searching, `search(fn, x0)` (the brief's `expand`, renamed). API approved: `docs/roots/PHASE5_API.md`.
  - [x] `search(fn, x0)`: the outward search (first step |x0| / 50, growth 3/2), domain edges probed (at most
    `digits` probes), `no_sign_change` when none is found
  - [x] Tests: the prototype's counts and brackets exactly; the domain edge; the composition; the edges; constant
    evaluation; the types; a demo section (five toolchains and WSL GCC and Clang, locally)
  - [x] CI passes on all five toolchains (ff58d32: run 40519529656)
  - [x] The independent reviews (Codex, and Kimi K3.1 Pro through Cline): six findings, all reproduced and fixed
    (`docs/roots/PHASE5_API.md`)
- Phase 6: `Ridders`, `RegulaFalsi`, `Illinois`, `Steffensen` (the user chose all four). API approved:
  `docs/roots/PHASE6_API.md`.
  - [x] `Ridders` (two `enclose` steps, the cancellation-free point, `detail::sqrt`), `RegulaFalsi` and `Illinois`
    (`false_position<M>`, guarded), `Steffensen` (Secant's step shared, the hand-over at z)
  - [x] Tests: the prototype's counts exactly; the roots Steffensen finds; the sqrt; constant evaluation; the types; the
    demo (five toolchains and WSL GCC and Clang, locally)
  - [x] CI passes on all five toolchains (840e04b)
  - [x] The independent review (Codex; Kimi K3.1 Pro incomplete): three findings, all reproduced and fixed
    (`docs/roots/PHASE6_API.md`)
- Phase 7: the method protocol, so that a method, the library's or a user's, is a type with `start` and `advance`
  (`docs/open-closed/ANALYSIS.md` section 10, step 2; the model is quadrature's Phase 3a). `solve` owns the loop,
  checks each step's bracket and sign change, applies `Criteria`'s stop test, forms `abs_error`, and ends every method
  (each bracket strictly inside the last); the observer moves to `math.hpp`; a latent bound above 2 tol in types
  without subnormals is fixed. The spike: `docs/open-closed/roots-spike/`.
  - [x] API approved: `docs/roots/PHASE7_API.md` (2026-11-06, Q1 to Q19 as recommended, draft 2; the one-call concepts
    deferred)
  - [x] Step 0: `nxx::observer` in `math.hpp`, quadrature's using-declarations, `iterative_method`'s state by
    reference; quadrature's tests unchanged, with `static_assert`s added for the clause and the shared observer
  - [ ] Step 1: on a bracket: the public types, the concepts, the loop, the fix, the nine bracketing methods stepped,
    the new forms; `tests/test_roots_protocol.cpp`
  - [ ] Step 2: from a guess: the four open methods stepped, the cap in the loop; the documentation and the guide page
  - [ ] CI passes on all five toolchains
  - [ ] The independent review (Codex) and its adopted findings

## Milestone 3: user documentation (Doxygen, Breathe, Sphinx)

A manual for library users, beside the API reference: Doxygen reads the headers and stays the gate that every public
entity is documented; Sphinx with Breathe and MyST (Markdown) renders the site, and from Phase 6 the same pages as a
PDF, so they are worded for both. The manual is organised by module (the user's choice, 2026-10-29, replacing Diataxis
at the top level chosen 2026-10-27), with Diataxis's four kinds of writing inside each module: an overview, a tutorial,
guide pages, reference. The deriv book stays the long explanation of `nxx::deriv` and is offered as its PDF, not
repeated. `numerixx_docs_all` builds the whole documentation: the site and the PDF manual, Doxygen's internal HTML and
the book. Each phase gets its plan approved before it is written, and the pages get the reviews of a phase. A box is
ticked when the docs build with no warning locally and in CI, and the examples' tests pass on all five toolchains.

- Phase 1: spike. Render `roots.hpp`, `deriv.hpp` and `math.hpp` through Doxygen XML, Breathe and Sphinx in a scratch
  folder, and list what the C++ domain misrenders or rejects (concepts, `requires` clauses, constrained `auto`,
  `std::expected` return types, deduction guides, overloads of `solve<M>`), with the workaround for each or the
  fallback (Doxygen HTML for the reference, Sphinx for the guides).
  - [x] Findings reported (Doxygen 1.18.0, Sphinx 9.2.0, Breathe 5.1.0, myst-parser 5.1.0, Furo 2025.12.19): Breathe
    parses `<requiresclause>` but drops it (breathe#482, open since 2020), although Sphinx accepts it since 5.2;
    `doxygenconcept` shows no definition; `doxygennamespace` omits concepts; one overload is selected only by its full
    parameter list. A patch in `conf.py` fixes the first two, and with it `sphinx-build -W` reports nothing. MkDocs
    was rejected: MkDoxy 1.2.8 reads neither clauses nor concepts, and MkDocs is unmaintained
  - [x] Toolchain chosen by the user (2026-10-27): Sphinx, Breathe and MyST for the guides and the full reference
- Phase 2: skeleton. Plan approved: `docs/sphinx/PHASE2_PLAN.md` (2026-10-27). `docs/sphinx/` (`conf.py`, pinned
  `requirements.txt`, MyST pages); `numerixx_docs` becomes the Sphinx site (Doxygen XML, then Sphinx `-E -W` through
  `uv run`), `numerixx_docs_internal` stays Doxygen HTML.
  - [x] The Breathe patch in `conf.py`, with Breathe pinned: the `requires` clause appended to the template prefix
    (`SphinxRenderer.create_template_prefix`), and a concept's `<initializer>` shown as a code block (the
    `Node_compounddefType` entry of `SphinxRenderer.node_handlers`)
  - [x] Doxygen: `EXCLUDE_SYMBOLS = *::detail *::detail::*` (deriv's `detail` leaks into the public docs now) and
    `GENERATE_DEPRECATEDLIST = NO`; Sphinx: `breathe_show_include = False`, names without a target ignored by
    `nitpick_ignore_regex`
  - [x] One reference page per namespace, with the concept and the `std::formatter` specialisations; no suppressed
    warning (`clang-cl`: no warning from Doxygen or Sphinx)
  - [x] Checks that fail the build: every class, concept and declaration in the XML rendered exactly once, and every
    clause and concept definition given exactly once, by Doxygen ID of what Breathe rendered (ad96830, after the
    reviews; the first version read the pages' directives)
  - [x] Examples as tests: pages show code only through `literalinclude` from standalone programs, each built and run
    as a test on all five toolchains
  - [x] The header `@code` blocks covered by tests in the module test files that assert their claims; the
    `deriv.hpp:14` comment corrected (it gave relative figures as absolute)
  - [x] A CI job (Doxygen 1.18.0 pinned, uv, Sphinx `-W`) (06267f3: run 40572439269)
  - [x] The reviews of the phase (`numerics-reviewer`) and the independent review (Codex): ten findings, all
    reproduced or checked and fixed (ad96830, the `deriv.hpp:310` comment after the user's approval;
    `docs/sphinx/PHASE2_PLAN.md`)
  - [x] After the phase, at the user's request: the site in Diataxis parts (`tutorials/`, `reference/`,
    `explanation/`, with `how-to/` from Phase 3), the book's PDF offered by the explanation part, the target
    `numerixx_docs_all`, and the patch and checks moved from `conf.py` to the local extension
    `docs/sphinx/_ext/numerixx_breathe.py`
- Phase 3: `nxx::roots`, the tutorial (Kepler's equation), seven how-to guides (choose a method, tolerances, find a
  bracket, a restricted domain, derivatives, other types, failures) and a reference page of costs, typed from the
  independent prototypes' counts and asserted by its program (the program that writes tables is dropped: a generated
  count would be Numerixx's own). Plan approved: `docs/sphinx/PHASE3_PLAN.md` (2026-10-27), with two comment fixes in
  `roots.hpp`.
  - [x] Pages, programs and tests; docs build clean (8422f4a: run 40584146349)
  - [x] Reviews: `docs-reviewer`, numerics, a docs-only reader whose four programs compiled and ran correctly; then
    the independent review (Codex). Findings reproduced and fixed (dbab9db: run 40586125853), and the header comment
    fixes the user approved (`docs/sphinx/PHASE3_PLAN.md`)
- Phase 4: the site organised by module, in a teaching voice (the user's request, 2026-10-29). Plan approved:
  `docs/sphinx/PHASE4_PLAN.md`, with the guide kept as seven pages. It absorbs the old Phase 4 (the roots explanations
  move into the roots overview).
  - [x] 4a: the pages moved into the module layout, links and toctrees fixed, no prose changed; docs build clean
    (ff7ab73: run 40965441740)
  - [x] 4b: `index.md`, `basics.md`, the roots chapter and deriv's short page in the new voice; docs build clean;
    programs pass on all five toolchains (83afe10: run 41019005701)
  - [x] 4b reviews: `docs-reviewer`, numerics, a docs-only reader (fixed in 83afe10); then the independent review
    (Codex): seven findings, all adopted (`docs/sphinx/PHASE4_PLAN.md`)
- Phase 5: the `nxx::deriv` chapter in the same shape, from `docs/deriv/DERIV_TUTORIAL.md`, with deriv examples in
  `basics.md`.
  - [x] Pages, examples and tests; docs build clean; programs pass on all five toolchains (4393031: run 41023892369)
  - [x] Reviews as in Phase 4b: `docs-reviewer`, numerics and a docs-only reader (fixed in 4393031); Codex: eight
    findings, all adopted (`docs/sphinx/PHASE5_PLAN.md`)
- Phase 6: one source, two manuals: the site and a PDF through Sphinx's LaTeX builder, with wording for both media.
  Plan approved: `docs/sphinx/PHASE6_PLAN.md` (2026-11-03; A4, the author's name, the book beside the PDF).
  - [x] Step 1: wording for both media, the voice rule and its check, the project texts (36e59c4: run 41677577233)
  - [x] Step 2: the PDF build (`manual.md`, `introduction.md`, the `pdf` tag, `numerixx_docs_pdf`); each overview's
    "Where to go next" split into Tutorial, Guide and Reference, so that the PDF shows them (c648bbf: run 41682779179)
  - [x] Step 3: layout, a read of the whole PDF, reviews and Codex (b57129c: run 41689684905; Codex: four edge cases,
    three adopted, `docs/sphinx/PHASE6_PLAN.md`)
- Phase 7: plots in the manual, drawn by matplotlib at build time from data a tested program checks; the first, in
  the quadrature guide, where the default and `TanhSinh` evaluate 1/sqrt x. Plan approved: `docs/sphinx/PHASE7_PLAN.md`
  (2026-11-03).
  - [x] The figure, its program and its script; built and read on the site and in the PDF; the program passing on the
    five toolchains, locally and in CI; `docs-reviewer` (20fe902: run 41716617342)
- Phase 8: the quadrature chapter's figures, six where a figure makes a concept or a method easier to see: a weighted
  sum, where the default cuts, the tanh-sinh map, the rounding floor, evaluations against digits, and the
  infinite-interval maps; the Phase 7 figure goes. Plan approved: `docs/sphinx/PHASE8_PLAN.md` (2026-11-03).
  - [x] The six figures, one commit each, with their programs on the five toolchains and `docs-reviewer` (a6c2085: run
    41819573354)
  - [x] CI passes; the independent review (Fable in place of Codex): no defect, two wording findings adopted
    (`docs/sphinx/PHASE8_PLAN.md`)

Open decisions: hosting (GitHub Pages from CI, or Read the Docs) and versioning (one site per release, or `dev` only),
with `linkcheck` and the install section of getting-started;
whether the headers' `@code` blocks move to `@snippet` from tests (a change to `include/`); whether to offer the
`requires` patch to Breathe upstream.

## Milestone 4: integrals of a function of one variable (`nxx::quadrature`)

Architecture approved: `docs/quadrature/BRIEF.md` (2026-10-30, D1 to D14 as recommended). Each phase gets its API
approved before it is written (`docs/quadrature/PHASE<n>_API.md`), as in Milestones 1 and 2.

- Phase 1: the core, `Integral`, `Failure`, `error`, `Criteria`, `integrate<M>(fn, a, b[, criteria])`; the fixed rule
  `GaussLegendre<N>` with the arithmetic-only node generator and compile-time tables; `Romberg`.
  - [x] API approved: `docs/quadrature/PHASE1_API.md` (2026-10-30, Q1 to Q7 as recommended)
  - [x] `quadrature.hpp`, constexpr; tests on float, double, long double and Boost.Multiprecision; docs build clean
    (all five presets and `numerixx_docs` locally, 2026-10-30)
  - [x] CI passes on all five toolchains (5f38606: run 41051017542)
  - [x] The independent review (Codex) and its adopted findings (six findings, all adopted:
    `docs/quadrature/PHASE1_API.md`)
- Phase 2: adaptive `GaussKronrod<15>`, the default (the heap, the floor, final panels, stagnation, the cap, the
  initial mesh, compensated sums; the node cache for types without a table; the README example and the demo).
  - [x] API approved: `docs/quadrature/PHASE2_API.md` (2026-10-31, Q1 to Q12 as recommended)
  - [x] `GaussKronrod<15>` and the default method, constexpr; tests in every type; docs build clean
    (all five presets, `numerixx_docs` and the demo locally, 2026-10-31)
  - [x] CI passes on all five toolchains (2098416: run 41240105813)
  - [x] The independent review (Codex) and its adopted findings (three findings, all adopted:
    `docs/quadrature/PHASE2_API.md`)
- Phase 3a: the method protocol, so that a method, the library's or a user's, is a type and adding one edits nothing
  (`docs/open-closed/ANALYSIS.md`: the user's decisions D1, open to users, and D4, quadrature first). The concepts
  `fixed_rule`, `iterative_method` and `estimating_method`; `integrate` owns the loop and the test against the
  criteria (`Step`); the method value first; the ends ordered by `integrate`; `qd::observer` for each step; Romberg
  and Gauss-Kronrod stepped, results unchanged except Romberg's floor on reversed ends.
  - [x] API approved: `docs/quadrature/PHASE3A_API.md` (2026-11-02, Q1 to Q15 as recommended, draft 5)
  - [x] The protocol in `quadrature.hpp`, constexpr; `tests/test_quadrature_protocol.cpp`; the existing tests
    unchanged; docs build clean (GCC 13 and Clang 18 locally, the 2,438-result battery identical, 2026-11-02)
  - [x] CI passes on all five toolchains and the docs build (c2398c9: run 41467982100, PR #19)
  - [x] The independent review (Claude Fable in place of Codex, 2026-11-02): no defect; six low findings, all
    adopted: three in the docs and tests (PR #19), and the three small fixes to `include/` after the user's approval
    (`PHASE3A_API.md`, "The independent review")
- Phase 3b: exp (internal) and `TanhSinh` on finite intervals, with the endpoint handling and the end-term test, as
  an iterative method of Phase 3a's protocol.
  - [x] API approved: `docs/quadrature/PHASE3_API.md` (2026-11-01, Q1 to Q10 as recommended; no K31, closing
    Phase 2's Q7), amended by `PHASE3A_API.md` G (2026-11-02)
  - [x] API approved: `docs/quadrature/PHASE3B_API.md` (2026-11-02, Q1 to Q5 as recommended, draft 2): TanhSinh
    through the protocol, in three steps
  - [x] Step 1, `detail::exp`: against mpmath within 1.5 ulp in every type, constexpr; all five presets locally and
    the docs build (2026-11-02; Emscripten with `-Wno-error=unused-template`, since the local emsdk 6.0.8 flags
    `test_quadrature_protocol.cpp`, which CI's 6.0.10 does not)
  - [x] Step 2, the tables in pieces: against mpmath within the a-priori bound, the counts per level, compile-time
    and run-time tables bit-identical, a level read at once; the budgets per toolchain in `design.md`; the
    prototype in `docs/quadrature/tanh-sinh-prototype/`; all five presets locally (Emscripten as in step 1) and the
    docs build (2026-11-02)
  - [x] Step 3, `TanhSinh`; tests in every type; docs build clean (all five presets locally, Emscripten 6.0.13 with
    no workaround, and the docs build, 2026-11-02; the numerics review of the implementation, nine findings, all
    adopted)
  - [x] CI passes on all five toolchains (f551f10: run 41528481205)
  - [x] The independent review (Codex) and its adopted findings (three findings, all adopted with the user's
    approval: the power law's products below min(), thresholds below min(), the search for x_min overflowing;
    `PHASE3B_API.md`, "The independent review")
- Phase 4: `TanhSinh` on half lines and the line (exp-sinh, sinh-sinh), infinite ends opt-in per method.
  - [x] API approved: `docs/quadrature/PHASE4_API.md` (2026-11-03, Q1 to Q8 as recommended, draft 3)
  - [x] Step 1, the table (internal), with Phase 3b's latent fix to `exp_constants`: against mpmath within the
    derived bound (largest ratio 0.40 double, 0.46 `cpp_bin_float_50`), the counts per level and side the
    prototype's, compile-time and run-time tables bit-identical; all five presets locally and the docs build
    (2026-11-03)
  - [x] Step 2, the infinite ends: the opt-in, `TanhSinh`'s maps and sides; tests in every type; docs build clean
    (all five presets locally and the docs build, 2026-11-03; the kept `es.py` with the C++ table agrees bit for bit on
    26 counted cases; the numerics review of the implementation: no defect, its corrections to comments and tests
    adopted)
  - [x] CI passes on all five toolchains (cc7e889: run 41637250176)
  - [x] The independent review (Codex) and its adopted findings (three findings, all adopted with the user's
    approval, all in a type whose max() is small against 1 / eps: the centre overflowing, 1 / eps not fitting, the
    radix wording; `PHASE4_API.md`, "The independent review")
- Phase 5: the module's chapter on the documentation site, in the shape and voice of the roots and deriv chapters
  (`docs/sphinx/PHASE4_PLAN.md`, "The voice"), and a pass over the Doxygen comments of `quadrature.hpp` that puts each
  entity's purpose first (the user's request). Plan approved 2026-11-03, after the plan's reviews
  (`simplicity-critic`, `docs-reviewer`, `numerics-reviewer`), with the user's decisions: the header keeps one sentence
  per limit with its remedy, and its measured figures (the catalogue of TanhSinh's limits, the timings) move to
  `docs/quadrature/design.md` and to the guide pages whose programs assert them; the pass covers the `detail`
  namespaces too, so that the implementation can be read from its comments; "Numerixx never throws" gains its
  exception, `std::bad_alloc` from an allocation (`GaussKronrod` always, `TanhSinh` in a type without a table), in
  `basics.md` and `docs/design.md` rather than in the code.

  The chapter, `docs/sphinx/quadrature/`, each page's program `docs/sphinx/examples/quadrature_<page>.cpp`:
  ```
  index.md               Overview: integrals with no closed form; a weighted sum of samples, and an error estimate
                         from two answers; the four methods and when each fits; the call; what you get back (an
                         Integral whose abs_error is an estimate, not a bound; a Failure and its best); what the
                         module does not do (sampled data, complex or multidimensional integrands, principal values);
                         where to go
  tutorial.md            The period of a pendulum at large amplitude, 4 sqrt(L/g) K(k), parametrised by the
                         complementary modulus k' = cos(theta0/2): K in its smooth form, 1 / sqrt(k'^2 + k^2 cos^2 phi),
                         by the default method; checked against pi / (2 AGM(1, k')), within abs_error plus the AGM's
                         bound; K in its algebraic form, whose singularity at t = 1 is moved to 0 by u = 1 - t with fn
                         written in u, by TanhSinh; k' toward 0 (amplitude toward 180 degrees), where K grows as
                         ln(4/k'); at k' = 0 the integral diverges, and TanhSinh returns not_converged with no best
  choose-a-method.md     \  From what you have to a method; an integral inside a root solve, an optimiser or another
                          |   integral (GaussLegendre<N> smooth in fn's parameters; Newton on F' = fn; inner tolerance
                          |   tighter than the outer; an inner Failure as NaN); TanhSinh's compile-time cost
  tolerances.md           | rel_tol and why sqrt(eps) (said here, linked from elsewhere); the floor; abs_tol, which
                          |   saves evaluations near 0; read abs_error against your request
  rough-integrands.md     | An endpoint singularity: TanhSinh; at an end other than 0, substitute with fn written in u;
                          |   the default's limit; kinks, steps and peaks inside: split at the known points; what no
                          |   estimate sees
  infinite-intervals.md   | Infinite ends with TanhSinh: the maps' scale and centre, slow decay, factors that overflow,
                          |   a large number written for infinity, oscillation
  high-precision.md       | A multiprecision type: rel_tol (the default keeps about half the digits), fn returning T
                          |   with no double literals, pi in T, the first call's table build; TanhSinh in
                          |   cpp_bin_float_50 against the AGM in that type
  handle-failures.md      | What each error means, a failure's best (meaningless for a divergent integral), recovering;
                          |   observer() to watch the steps
  method-of-your-own.md  /  The three concepts; one fixed rule with a parameter, checked by a static_assert
  costs.md               Reference: evaluations and accuracy per method on 4 / (1 + x^2), sqrt x and 1 / sqrt x on
                         [0, 1], and 1 / (1 + x^2) on the line (TanhSinh); only rows the prototypes reproduce
  reference.md           The API (as now; its interim introduction replaced)
  ```

  Outside the chapter: `index.md` gains the module's row and the toctree entry `quadrature/index`, and its opening
  no longer calls integrals planned; `basics.md` names three modules, states the `std::bad_alloc` exception
  (as `docs/design.md` does), adds an interval's ends where it says the type
  comes from the start, and in "What a result holds" gains an integral (an `Integral`'s `abs_error` is an estimate and
  a `T`; the error side is a `Failure` with a best, not a bare `error`) and the collision of `Criteria` and `error`
  with `nxx::roots`' under using-directives (said there once, linked from elsewhere).

  Examples: quoted output only in `double` and from integrands built from `+ - * /`, `std::sqrt`, comparisons,
  `std::abs` and min/max, which IEEE 754 makes the same on every toolchain, under the FMA caveat of `basics.md`
  ("Compile time"), repeated on each page that quotes; integrands through libm, and the amplitudes whose k is not exact
  or from a square root, are described by the bounds their programs assert. References are closed forms or, for K, the
  AGM (stopped at a - b <= 4 eps a then one more step, since `while (a != b)` can cycle between adjacent numbers;
  never at k' = 0), never Numerixx. The counts in `costs.md` come from the kept prototypes
  (`docs/quadrature/tanh-sinh-prototype/`, Phase 2's `final2.py`), and its program asserts them.

  Doxygen: every entity of `quadrature.hpp`, public and `detail`, opens with its purpose; for the public ones (the
  namespace, `error` and its enumerators, `Integral`, `Failure`, `Criteria`, `Step`, the three concepts, `state_t`,
  `Observer` and `observer`, the four methods and their members, the six forms of `integrate`, the two formatters) what
  the caller uses it for and when, in the caller's terms; for a `detail` entity why the implementation needs it and
  which decision of `design.md` it carries out; then what it does; then its limits. Function bodies gain comments where
  a step's reason is not plain from the code. The `@brief` is one sentence stating the purpose. Known cases: the briefs
  of `GaussKronrod`, `GaussLegendre`, `TanhSinh` and `Step` say what each is, not when to use it; `Romberg`'s purpose
  comes last, and "equally spaced data" misleads, since `integrate` takes only fn; the method choice sits in one
  overload's comment. Comments only, no code token changes; every `@tparam` and `@param` kept (WARN_NO_PARAMDOC fails
  the build); "a < b" spaced; no line starting with a number or "- " unless a list is meant. Each `@code` block has a
  test that asserts its claims; the missing ones are added.

  Steps, each built on the five toolchains, locally and in CI, with the docs build clean:
  - [x] 5a: the Doxygen pass, shown to the user as a diff before it is committed (`include/`): every entity, `detail`
    included, purpose first; code tokens unchanged (compared with comments stripped); TanhSinh's measured figures in
    `design.md`, "TanhSinh's measured limits"; three comment defects fixed (the heap's tie-break described backwards,
    `exp_sinh_start`'s comment above `inverse_epsilon`, `KronrodState`'s undocumented constants); the header's
    examples run as written in `tests/test_quadrature_examples.cpp`, and `Midpoint{.panels = 1000}` in the protocol
    tests. All five presets locally, both Doxygen builds clean (2026-11-03); CI (cdca6e1: run 41655023180)
  - [x] 5b: the chapter's pages and programs, and the `index.md` and `basics.md` changes. Done as planned, with
    `std::bad_alloc` stated in `basics.md`, the site's `index.md`, the README and `docs/design.md`; the prototypes
    behind `costs.md` kept in `docs/quadrature/costs-prototype/` (`costs.py` prints every count it reproduces); the
    concepts' definitions rendered with their default argument (Doxygen 1.18 drops the "=" in a concept's initializer;
    the extension rebuilds the template head from the parameter list). All five presets locally and the docs build
    (2026-11-03); CI (d0e4a10: run 41663101675)
  - [x] 5b reviews, as in Milestone 3 Phase 5: `docs-reviewer`, `numerics-reviewer`, `docs-reader` (tasks: the
    probability that a standard normal variable exceeds 3; an integrand with a kink at a known interior point; a time
    from an integral, solved for its upper limit with `nxx::roots`; an integral in `cpp_bin_float_50` to 40 digits; a
    fixed rule of one's own), then Codex on the commit; the adopted findings recorded in this entry. The reviews (all
    findings reproduced or derived, all adopted): the tutorial's accuracy loss near 180 degrees is the coarse doubles
    near pi/2 (the sliver fl(pi/2) misses, about 6e-17 / k', and the rounding of the points there), not the rounding of
    cos, which the first draft blamed; the default steps around an isolated non-finite value anywhere, not only near an
    end; a divergent integral can succeed at a loose request (1/x at rel_tol 0.2: ok, 9.33 after 150 evaluations), and
    the default's best on 1/x is the same at every request, so the pages teach the real signs (no best from TanhSinh, a
    result growing as the singular end moves); x^51 e^-x fails on inf * 0 past 1.1e6, not on x^51 overflowing first;
    sin x / x on [0, inf) gives no best in the C++ (the prototype's best below its error, recorded in `design.md`, was
    never asserted); TanhSinh's prototype reads the header's nodes, which `costs.md` now says; the floor, the slivers'
    2e-8, the AGM's error budget, cpp_bin_float_50's epsilon (2^-167) and the Simpson bound corrected; the reader's
    gaps closed (the mesh's panel boundaries and a kink at a tight request, the Boost header, a `Criteria` variable,
    `qd` at file scope); "the default method" throughout, since C++ has one. The tutorial stays one page at about
    1,450 words with its directives, above the 1,000 of the voice rules: one continuous problem.
    The independent review (Codex, on cdca6e1 and d0e4a10): nine findings, all reproduced or derived and adopted. In
    the docs: an inner integral that fails at one outer point is stepped around, the outer one succeeding (asserted,
    with the advice to count inner failures); Boost's run-time-precision mpfr_float is not an nxx::real
    (is_specialized false), so the high-precision page no longer names it; GaussLegendre<N> takes no criteria;
    constant evaluation also refuses sums beyond max / 2 that run time accepts (basics.md); a half line from 1e17
    stops because the map's centre a + 1 rounds onto a, before any evaluation (also in design.md); slow decay needs a
    wider exponent range, not more digits (asserted in cpp_bin_float_50); a split kink's error "within a few ulps",
    and a tolerance derived. In quadrature.hpp, after the user's approval, comments only: Neumaier's sum errs by eps
    |sum| plus a second-order term, not one rounding whatever the terms (2^100, 1, 2^-100, -2^100, -1 sums to 0);
    TanhSinh's slow-decay remedy and far-half-line cause as above; a bad node is avoided because the halves are
    sampled afresh, not because it becomes an end.
- Milestone 4 is complete: the module, its chapter and the comments of its header.

## Deferred

### Documentation

- The same purpose-first Doxygen pass for `roots.hpp`, `deriv.hpp` and `math.hpp`, if the user wants it after
  `quadrature.hpp`'s.
- In the quadrature chapter: combinators over methods (a substitution as a method), until a reader asks.
- The reference in print (Milestone 3 Phase 6's read of the PDF, items 6 and 7): `@code` examples in the headers that
  lose their indentation (on the site too); header comments that cite internal files (`PHASE4_API.md`,
  `docs/<module>/design.md`, `costs.py`); ASCII math and curly-quoted code in reference prose; `@tparam` on a concept
  rendered as "tparam T"; declarations split from their template line at page breaks, near-empty pages, and the
  hyphenated headers of the roots costs table. Changes to `include/` need the user's approval.

### Derivatives (`nxx::deriv`)

- Gradients and Jacobians, until a module needs them.
- For `derivative_of` (`docs/roots/PHASE3_API.md`): a step or scale (a lambda around `derivative(fn, x, h)` does it
  now), and a second call operator that takes fn(x) from the caller, so that Halley's numerical f'' need not evaluate
  it again.

Customization (reviewed 2026-10-25 with a usage sketch and a simplicity review; every item can be added without
breaking the API, so none needs a hook now). For `nxx::quadrature` the open method protocol is Phase 3a
(`docs/quadrature/PHASE3A_API.md`, approved 2026-11-02); the same for deriv and roots follows it
(`docs/open-closed/ANALYSIS.md`, section 10):
- An open method protocol, so that a user's method works in `derivative<M>`: method objects with a const member
  template, which also carry parameters. Until a library algorithm is generic in the method (gradients, Jacobians,
  `derivatives<M, N>`) and a user needs a method of their own in it. A spike compiled it in constant expressions on
  GCC, Clang and MSVC.
- fn's noise level (fn accurate to eta rather than eps), attached to fn, e.g. `with_noise(fn, 1e-10)`, for every
  method: Ridders' floor and the stencils' default step. The only item with a numerical case: at relative noise 1e-8
  Ridders' estimate covered the error in 90.7% of the book's cases, against 100% on clean fn, and a caller cannot fix
  that from outside. Until a caller needs a reliable estimate on noisy fn; a prototype study first (coverage, and the
  untrusted test that a larger floor weakens, on the pole set).
- Ridders' cost and tolerance (an evaluation cap, an early stop at a tolerance), returning the best estimate so far
  rather than an error. Until a measured need that a passed step cannot meet.
- Stencils from arbitrary nodes (`Stencil<-1, 0, 1, 2>`), from the existing generator. Until a node set no built-in
  family gives; the accuracy order of symmetric sets is not nodes - N, which the default step must account for.
- Not planned: a scale parameter (pass `default_step<M, N>(x) / max(|x|, 1) * w` for fn's scale w), exposing
  Ridders' tuned constants, observers (wrap fn), non-uniform points (run-time weights).

Complex numbers (2026-10-25; `docs/deriv/DERIV_TECHNIQUES.md` 2.5 and 2.7):
- A complex value for a real argument (f: R -> C): the same algorithms, with a value type separate from the argument
  type, the modulus in the floors and estimates, and overflow guards per component. Until a caller differentiates a
  complex-valued fn.
- A complex argument for an analytic fn: the same stencils and Ridders along a real step, and the Cauchy contour
  method, which gives every order to near machine precision from one set of samples. Needs a generic fn and a complex
  type; until a module needs derivatives of analytic functions (a root finder for analytic, non-polynomial fn, or
  Taylor coefficients).
- Non-analytic fn of a complex argument (Wirtinger derivatives, the 2x2 real Jacobian): with gradients and Jacobians.
- The complex step (one evaluation, about eps for f'): needs fn written for complex arguments; Boost.Math has
  `complex_step_derivative`. With the contour method, if ever.
- Polynomials (planned) get exact derivatives from Horner's scheme and synthetic division, not from `nxx::deriv`. Their
  complex type (a concept for real and multiprecision types, with `abs` and `isfinite`) is the foundation the items
  above would reuse.

### Build and tests

- Clang 23 (23.1.4, WSL, 2026-11-03) fails `test_quadrature_tanh_sinh.cpp`: `ts_reciprocal_result`, 1/x on [0, 1] run
  to the cap in constant evaluation, exceeds the default constexpr step limit ("hit maximum step limit" in
  `nodes.each`, `quadrature.hpp`), where CI's Clang 22 passes. Before CI moves to Clang 23: the test's budget or the
  library's step count.
- With Clang, a class-type real (the tests' `Wrapped` and `Capped`) makes TanhSinh build its levels in constant
  evaluation during code generation: 31 s and 68 s of a test file's compile, where `double` costs under a second
  (`-ftime-trace`, 2026-11-03). A user with a type of their own pays it in every file that integrates with it. A
  question for `include/`, with the user's approval.
- CI runs clang-cl only on pull requests, pushes to master and manual runs (2026-11-03): ccache cannot keep its
  dependency files (`CMakeLists.txt`), so it is the one slow job. Run it on `dev` by hand
  (`gh workflow run ci.yml --ref dev`) when a change could differ there.

## History

- 2026-10-20: the first rewrite was abandoned; master went back to the 1.x code (`v1.1.0-legacy`) and Numerixx 2
  restarted on `dev`.

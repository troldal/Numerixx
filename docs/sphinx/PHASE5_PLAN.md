# Milestone 3, Phase 5: the `nxx::deriv` chapter

Status: done (2026-10-30), at the user's request after Phase 4b. The shape and the voice are those approved in
`PHASE4_PLAN.md` ("Phase 5: the `nxx::deriv` chapter in the same shape"); this document records the pages and their
programs.

## The pages

```
deriv/
  index.md               Overview: the problem; the finite difference and the choice of step; stencils and Ridders'
                         method; the call; what you get back; where to go
  tutorial.md            The slope of Kepler's E(M), a function computed by a root solve: dE/dM and d2E/dM2, checked
                         against their closed forms
  choose-a-method.md     \
  step.md                 |
  one-sided.md            |  The guide, one task per page
  higher-derivatives.md   |
  derivative-of.md        |
  handle-failures.md     /
  accuracy.md            Reference: evaluations and accuracy per method on three functions
  book.md                The book (as now)
  reference.md           The API (as now)
```

`basics.md` gains a derivative example in "What a result holds" (an `Estimate`, its optional `abs_error`).

## The examples

Each page's program is `docs/sphinx/examples/deriv_<page>.cpp`. Quoted output must be the same on every toolchain, so
the functions whose output is quoted use only `+ - * /` and `std::sqrt`, which IEEE 754 rounds correctly: 1/x,
sqrt(x), x / (x^2 + a^2), 1/(1 + x) with NaN below 0. The tutorial's function goes through `std::sin`, so its page
describes the results by the bounds its program asserts, as the roots tutorial does. References are closed forms.

## Source material

`docs/deriv/DERIV_TUTORIAL.md` (the guided tour of the implementation), `docs/deriv/design.md` and the phase documents
for the numbers measured on the prototypes (Ridders' coverage, medians, costs). The book stays the place for the
theory; the overview gives the one idea a user needs, the trade-off between truncation and rounding.

## Reviews

As in Phase 4b: `docs-reviewer`, `numerics-reviewer`, `docs-reader` (tasks: the derivative of a function given by a
root solve; a derivative at the edge of a domain; f'' with an error estimate; a derivative passed to Newton's method),
then Codex on the commit.

## Success criteria

- [x] The pages and programs; every quoted output compared by its program, every number asserted or linked; the docs
  build clean; all programs pass on the five toolchains, locally and in CI
- [x] The reviews and Codex, with the adopted findings recorded here

## What was done, and the reviews

Done as planned. The examples use 1/x, sqrt x, 1/(1 + x^2), the peak x / (x^2 + a^2) and 1/(1 + x) with NaN below 0,
whose output is the same on every toolchain (checked on all five); the step page computes its scaled step with
`nxx::math::root_epsilon`, since `std::cbrt` need not be correctly rounded. `basics.md` gained a derivative example and
the `int` derivative case; `choose-a-method.md` a section on a multiprecision type.

Reviews (`docs-reviewer`, `numerics-reviewer`, a docs-only reader run from the new `docs-reader` agent's instructions),
each finding reproduced or derived, all adopted unless noted:

- The trust test: the page taught "estimate below a quarter of the value" as the test before using a result. That is
  only the methods' own give-up flag (`ForwardRidders` on sqrt at 1e-11 is 7% off with an estimate of 11%, and passes
  it; on a derivative that is 0 it misfires). The page now teaches a test against the accuracy the reader needs, with
  the flag described as such, and a note on derivatives near 0.
- Only non-finite values are detected: a formula that goes on past the domain, or a kink among the points, gives a
  wrong result with a small estimate (Ridders at a kink at 0: -2, the mean of the slopes, estimate 1e-8). Said in
  `one-sided.md` with an example, in `handle-failures.md` and the overview.
- A NaN at a point fails a stencil, while Ridders' methods restart with smaller steps (overview, tutorial, failures).
- Near an edge: on 1/(1 + x) `ForwardRidders` is more accurate at 1e-6, on sqrt(1 - x) near 1 `Ridders` is (the
  reader's case, now asserted); "try both, keep the smaller estimate", and a one-sided stencil's default step does not
  shrink with the distance to the edge (`Backward<2>` 71% off at 1e-7 from it, asserted). The method table has a
  near-edge row.
- "Forward<2> gives fewer digits than the central difference would" was false on that function at 1 (5.2e-12 against
  2.8e-11): dropped.
- Numbers and tolerances: the digits per order (4 or 5 for f''''), the error models in the programs' comments (eps^(1/2)
  = 1.5e-8, eps^(2/5) = 5.5e-7), tolerances from derived bounds (choose-a-method: 2e-10, 3e-12, 2e-12; f''': 4e-6),
  the tutorial's quoted digits and f'' bound (below 1e-10) asserted, "about a dozen evaluations", "13 to 15 digits"
  for Ridders, "about a sixteenth" of the peak's width; Ridders given a step starts from it (14 evaluations, asserted).
- Gaps: a derivative in `cpp_bin_float_50`, Halley's f'' through `derivative_of`, the macro on the command line, the
  `dv` and `rt` aliases named, the deriv header on `basics.md`, terms defined at first use, one wording for the FMA
  caveat.
- MSVC's C4459 on a program's global `a`, which a local of `deriv.hpp` hides: the program's name changed; not a header
  defect.
- Header comments, proposed to the user rather than changed (`include/` needs approval): `deriv.hpp:278` and
  `roots.hpp:1730` start a line with "129." and "0.", which Doxygen reads as an ordered list (the Ridders reference
  loses "129." and a paragraph renders as a list item); `ForwardRidders`' "where fn is defined on both sides, use
  `Ridders`, also near a domain edge" is not borne out on 1/(1 + x) at 1e-6.

The independent review (Codex, on 4393031; eight findings, each reproduced or checked, all adopted):

- "Keep the result with the smaller estimate" near an edge can pick a wrong result: on u (1 - (64u)^2) (1 - (128u)^2),
  u = t - 1/64, NaN below 0, Ridders at 1/64 returns 0 with an estimate of 0 (its steps land on the polynomial's
  zeros), ForwardRidders 1, the derivative. The page now says to compare, keep the smaller estimate only when the two
  agree, and that an estimate is not a proof; the case is asserted in `deriv_one_sided.cpp`.
- `derivative_of(std::ref(fn))`: the original must outlive the callable, and the copy does not extend what `fn`
  captures by reference.
- The `dv` alias named on `higher-derivatives.md` and `derivative-of.md`.
- The accuracy table's coverage check now requires the Ridders rows to have an estimate.
- The tutorial's ten decimals are checked on the computed slope, not only the reference; "70% off" checked as
  65% to 75%; the Ridders cost ranges linked to `nxx::deriv::Ridders`.
- A step "so that x + h rounds to x" does not always fail (the two-point rules need only their two points to differ:
  h = 0x1.8p-54 at 1 succeeds): the pages now say "the method's points round to the same number".
- "Fewer digits than `Central<2>`" for `Forward<2>` is "typically"; the accuracy table compares "the main methods".


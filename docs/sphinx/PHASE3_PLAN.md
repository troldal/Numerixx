# Milestone 3, Phase 3: the `nxx::roots` pages

Status: approved by the user (2026-10-27) after the reviews of the plan (`simplicity-critic`, `numerics-reviewer`,
`docs-reviewer`; the record is at the end), with the comment fixes, and split: the two explanations are Phase 4.
Milestone 3 is in `docs/plan.md`; Phase 2 built the site, its checks and its Diataxis parts
(`docs/sphinx/PHASE2_PLAN.md`). The ground truth for every statement is `include/numerixx/roots.hpp`, its tests and
`docs/roots/design.md`.

## Scope

The pages a user of `nxx::roots` needs beyond the generated reference: a tutorial, how-to guides and a reference page of
costs. The explanations below are Phase 4's; until then the pages link to the reference. `include/` changes only by the
comment fixes below, which the user approved.

## Rules for every page

- Each page does one Diataxis job. A page links to the reference for what a declaration does (MyST cross-reference
  roles such as `{cpp:struct}`, checked by `nitpicky` and `-W`) and does not restate it.
- Code appears only through `literalinclude` from a program in `docs/sphinx/examples/`, named after its page
  (`roots-kepler.md`, `roots_kepler.cpp`), built and run as a test.
- Every number on a page is asserted by that page's program against an independent source cited in a comment (a
  closed form, mpmath, a prototype count recorded in `docs/roots/design.md`, or a derivation in the comment), or is
  a link to the reference or test that asserts it. No output is typed unchecked; no evaluation count is stated for an
  fn that calls libm, whose counts differ between platforms.
- Formulas in plain text, as the header comments write them; the math extension stays deferred.

## Pages

**Tutorial** (`tutorials/roots-kepler.md`, `roots_kepler.cpp`): Kepler's equation E - e sin E = M, e = 0.5, M = 1,
from where getting-started ends. Write fn as a lambda capturing e and M; derive the bracket [M - e, M + e] from
|E - M| = e |sin E| <= e; solve with `Brent` and read `value` and `abs_error` by member name (not a structured
binding); print the actual error against the reference beside `abs_error`, and check by substitution; give a bracket
without a sign change and handle `no_sign_change`; solve from the guess M with `Newton`, f and f' as a
`std::pair<double, double>` (the pitfall `std::pair{f, 1}` deduces `pair<double, int>`). f' = 1 - e cos E >= 0.5, so
the root is unique. Reference: E = 1.498701133517848314057985497256239901590695552547389812626730414735991
(mpmath's findroot at 100 digits, and the Bessel series E = M + sum 2/n J_n(n e) sin(n M) as a second, independent
computation; they agree to 1e-101); the
check is |value - ref| <= abs_error + c, c = 2.3e-16 (sin's 2 ulps times e, the subtraction's half ulp of 1, over
f' = 0.964) plus half an ulp of the root for the double-held reference. No counts.

**How-to guides** (`how-to/`), each a task and its program:

| Page | Task |
|---|---|
| `roots-choose-a-method.md` | Pick a method: a bracket or a guess, f' available or not, a multiple root, the cost of fn; measure the evaluations on your own fn with a counting wrapper. Links to the costs table (reference) and to the reference for each method |
| `roots-tolerances.md` | Set `abs_tol` and `rel_tol`: tol(x) = abs_tol + rel_tol \|x\| acts at every magnitude (abs_tol = 1e-12 caps a root at 1e-9 at about 3 digits); the walk to a sign change near 0 that fn misses, and which methods walk (Bisection on x - 1e-300 over [0, 1]: 1049 evaluations, 41 with abs_tol = 1e-12); fewer digits for less cost (Bisection on x^2 - 2: 52, 21 with rel_tol = 1e-6); values below the floors are raised, negative ones are `non_finite`; the bound abs_error <= 2 tol(value) |
| `roots-find-a-bracket.md` | Find a bracket from a guess with `search` and chain it into `solve` with `and_then`; `search(fn, 0)` deduces `int` |
| `roots-restricted-domain.md` | Solve where fn is defined only on part of the line: return NaN outside it; `search`, then a bracketing method, which never leaves the bracket; an open method gives `non_finite` when an iterate leaves the domain; a numerical derivative's stencil near the edge; van der Waals, as in the demo |
| `roots-derivatives.md` | Give f' (and f''): fn returning a `std::pair` (a `std::tuple` for Halley; `std::array` in a sentence), separate callables, `nxx::deriv::derivative_of`; the calls each costs, from the Phase 3 prototype's `Central<2>` counts only |
| `roots-other-types.md` | Solve in `float`, `long double` and Boost.Multiprecision, and in a constant expression (its own section). `long double` is double on MSVC, x87 extended on GCC and Clang on x86-64, binary128 on Emscripten: take tolerances from `numeric_limits`. With expression templates (`cpp_dec_float_50`, `number<cpp_bin_float<50>, et_on>`) fn returns T explicitly, and `std::pair<T, T>`; maths functions by ADL. The walk to 0 is about 2e9 halvings in `cpp_bin_float_50`: set `abs_tol`. Constant expressions: fn of `+ - * /` only, built-in types only (Boost's are not literal types, and `std::sin` is not constexpr) |
| `roots-handle-failures.md` | React to a failed solve: test success first, and recover by what to try next (a wider bracket or `search` after `no_sign_change`; a bracket after `not_converged`; the domain after `non_finite`); a `switch` on the error has a `default`, since causes may move between values; `value_or`, `or_else`, `std::format` of the error; `.value()` throws `std::bad_expected_access`, fn's exceptions pass through; check a result for a pole with \|fn(value)\| |

**Reference** (`reference/roots-costs.md`, `roots_costs.cpp`): evaluations per method, typed from the independent
prototypes' counts and asserted exactly by the program, which also checks each value within `abs_error + c` of its
closed-form root and `abs_error <= 2 tol(value)`. One evaluation is one call of fn (giving f, with f' for Newton and
f'' for Halley); valid for x86-64 and WebAssembly builds without FMA contraction, and for fn written as the
prototypes write it.

- Bracket table, every bracketing method plus Newton and Halley (analytic, and Newton with `derivative_of`): x^2 - 2
  on [1, 2], x^3 - 2x - 5 on [2, 3], x^20 - 1 on [0, 5], (x - 1)^3 on [0, 3], 1/x - 1 on [0.5, 3]. Sources:
  `docs/roots/design.md` Evidence, Phases 1 to 4 and 6; Bisection on 1/x - 1 is 54, derived (half the width
  2.5 * 2^-(k+1) against tol near 1, about 2^-51) and given by `.scratch/roots1/proto3.py` (rerun 2026-10-27).
- Guess table: x^2 - 2 from 1 with Newton 7, Halley 5, Secant 9, Steffensen 10 (Steffensen's count and root depend
  on fn's scale: a footnote), and Newton with `derivative_of` 21.

**Explanation** (`explanation/`; Phase 4), no programs unless a point needs one:

- `roots-sign-changes.md` (bracketing and open methods): a result is a verified sign change of fn as evaluated, or a
  point where fn is exactly 0 (an end, zero first, a midpoint, an iterate), which may be an even root or another root
  than the bracket's sign change; an even root cannot be bracketed; zero first evaluates fn(0) in a bracket that
  contains 0; the hand-over, and so which root an open method finds; why Brent is the default and which methods are
  teaching examples, each linked to its reference entry; the sources as `docs/roots/design.md` gives them, with how
  each was checked.
- `roots-accuracy.md` (what a result means): `abs_error` bounds the distance to a sign change of fn as evaluated, at
  most 2 tol(value), and tol(value) where fn is exactly 0; conditioning, fn's rounding over |f'| at a simple root
  and (m! delta / |f^(m)|)^(1/m) at a root of multiplicity m; the default tolerances and why they assume no scale; a
  pole is returned as a root, or gives `non_finite` where a point lands on it, depending on the method.

## Navigation

`index.md` gets a "How-to guides" toctree with the caption "How-to: roots" (one per module, as Phase 4 adds deriv's),
the intro sentence names the how-to guides again, the costs page joins the reference toctree, and getting-started's
last line points to the roots tutorial.

## Examples

`tests/CMakeLists.txt` replaces its single example target by a `file(GLOB ... CONFIGURE_DEPENDS)` over
`docs/sphinx/examples/*.cpp`: each an executable and a test `numerixx_docs_example_<name>`, with the tests' warning
flags, `boost_multiprecision` (header-only) and, on MinGW GCC, `stdc++exp`. The CI docs job compiles no C++.

## Comment fixes to propose (`include/`)

- `roots.hpp:1900` and `:1930`: the bracket overload's `@brief` and `@tparam M` omit `Ridders`, `RegulaFalsi` and
  `Illinois`, which `detail::bracketing` accepts.
- `roots.hpp:1970` and `:1987`: the guess overload's `@brief` and `@tparam M` omit `Steffensen`.

## Reviews of the pages

`docs-reviewer`, `numerics-reviewer`, and a docs-only reader (`usage-sketcher`, given only the rendered pages) who
writes client code for four tasks that must compile and give correct results: a table of Kepler roots over M; van
der Waals with `search` and NaN outside the domain; Newton with `derivative_of`; a root in `cpp_dec_float_50`. Then
the independent review (Codex) on the phase's commit.

## Success criteria

- [ ] The pages, their programs and the costs table; `numerixx_docs` builds with no warning and its checks pass
- [ ] Every number on a page asserted by its program against an independent source, or linked to what asserts it
- [ ] Every program passes as a test on all five toolchains, locally and in CI
- [ ] The docs-only reader's code compiles and gives correct results, or each failure is traced to a page and fixed
- [ ] The reviews and the independent review, with the adopted findings recorded here

## Changes to `docs/plan.md` and `CLAUDE.md`

- Plan, Phase 3: the pages above; the costs table typed and asserted, not generated (the program that writes tables
  is dropped: a generated count would be Numerixx's own, and Phase 2's rule is that a count is derived again or the
  prototype rerun, never copied from Numerixx); the explanations as above.
- `CLAUDE.md`, Layout: the examples' glob.

## The split (user decision, 2026-10-27)

Phase 3: the tutorial, the seven how-to guides and the costs page (nine pages, nine programs). Phase 4: the two
explanations. Phase 5: `nxx::deriv`'s pages.

## Deferred

- A generated table, and a table with libm fns (cos x - x, e^x - 2), whose counts differ between platforms: if a
  reader asks.
- A migration how-to (Boost.Math's `toms748_solve`, SciPy's `brentq`, GSL's driver loop): when asked.
- A per-method history page: when asked. The math extension: when a formula needs more than plain text.

## Reviews of this plan

- `numerics-reviewer`: "every result is a verified sign change" and "no method can find an even root" are false (an
  exact zero is returned wherever it lands: x^2 (x - 1) on [-1, 2] returns 0 by zero first); counts hold only without
  FMA contraction and for the prototypes' forms of fn; the table's counts must be asserted, with Bisection on 1/x - 1
  derived (54, reproduced by the prototype); conditioning at a multiple root goes as the 1/m power; poles can give
  `non_finite`; Kepler's c = 2.3e-16 plus half an ulp; the tolerance claims and the `long double` table. All adopted.
- `docs-reviewer`: a failures guide instead of telling the errors apart (the header warns that causes may move); the
  expression-template and constexpr rules in the types guide; a restricted-domain guide; the costs table as
  reference; no `search` in the tutorial, which derives its bracket; navigation, names, cross-reference roles; the
  comment fixes above; the docs-only reader's tasks. All adopted.
- `simplicity-critic`: no generator (typed, asserted counts); one examples glob in `tests/CMakeLists.txt` linking
  Boost for all; the methods explanation merged into the sign-changes page; one return form shown in the derivatives
  guide; no Ridders-derivative counts; formulas as plain text. Adopted. Not adopted: cutting the errors guide (the
  docs review shows the reader needs the recovery it teaches) and not splitting (the pages grew from the other
  reviews; the user decides).

## Reviews of the pages (8422f4a)

CI passed on 8422f4a (run 40584146349). Then `docs-reviewer`, `numerics-reviewer`, a docs-only reader and the
independent review (Codex); every finding was reproduced or checked before it was acted on.

- The docs-only reader wrote the four programs (Kepler roots over M; van der Waals with `search` and NaN; Newton with
  `derivative_of`; a root in `cpp_dec_float_50`) from the rendered pages alone. All four compiled under the tests'
  warning flags and gave correct results. Its gaps, adopted where a page was wrong: the `mp` alias outside the shown
  snippet; the rest (a table of call shapes, printing Boost values, an install section) deferred below.
- The pole check (all four reviews): the page taught |fn(value)| against its size at the ends, the rule
  `docs/roots/PHASE1_API.md` E rejected; reproduced, the resonance (x^2 - 2) / ((x^2 - 2)^2 + 1e-24) on [1, 2] is
  called a pole at sqrt(2). Fixed: the local test E describes (|fn| at value against value +- 1024 abs_error), asserted
  on the pole 1/(x - 1/3) and E's three cases.
- The method table (all but the reader): Brent "never much slower than bisection" and Ridders "cheapest on (x - 1)^3"
  contradicted the costs table (156 and 86 against Bisection's 54). Fixed; the header says the same (below).
- `search` "cannot find an even root" (Codex): it returns {1, 1} for (x - 1)^2 from 1, where fn is exactly 0.
  Reproduced; fixed, with the guess's scale and the root nearest the guess.
- The van der Waals text credited the root to the edge approach (docs review): the upper side finds it, and without the
  NaN the same code returns the pole at b. Fixed, and the pole asserted.
- A bracket keeps fn's own points inside the domain, not a numerical derivative's nodes (three reviews); the ends must
  be finite. Fixed.
- The tutorial printed an "actual error" of 0, against the reference rounded to a double (Codex, docs review): the
  reference is now held as two doubles (the error prints 8.0e-17). Its `int` pitfall was wrong for the shown
  expression (numerics, Codex): reworded, and `std::pair{f, 1}`'s rejection asserted. "Two iterates have opposite
  signs": fn's values at them.
- `abs_tol = 1e-12` "comes out" at 3 digits: it is the guarantee; Brent returns this linear root almost exactly
  (error 2e-25, reproduced). Reworded. The walk: a tiny root down to its own magnitude; 2e9 halvings is
  `cpp_bin_float_50`'s.
- Unchecked numbers (Codex, docs review): numerical Halley's 6 calls per point (now 30 from 2, the Phase 3
  prototype's) and Newton's first iterate -19.9 (now recorded) asserted. `roots_costs.cpp`'s numerical rows now check
  the `abs_error` contract too.
- The other-types page: expression templates refer to fn's locals or by-value parameters, `cpp_bin_float_50` has them
  off; the constexpr maths functions; the FMA condition on "the constant gives what a run gives" (not reproduced with
  GCC `-mfma -ffp-contract=fast`, but stated); 80-bit x87. Fixed.
- Tolerance derivations in the programs' comments (numerics): ln(x / 1e-12) c 2e-28 (was 7e-28), van der Waals 5e-19
  (2e-18), the pole's c 0, the tolerances header. Fixed; all within their bounds on every toolchain.
- Inline code in prose (`search(fn, 0.0)`, a type name): allowed; the rule is that code blocks come only from programs.

Proposed to the user (changes to `include/`), from the same findings:
- `roots.hpp`, `Ridders`: "the cheapest method here on roots of odd multiplicity" is false (Bisection and Halley 54
  on (x - 1)^3, Ridders 86); and `solve`'s "Choosing a method" recommends Ridders there and lists Steffensen among the
  bracket overload's teaching examples.
- `roots.hpp`, `solve`: "where fn may have poles, check |fn(value)|" teaches the unsound test; say how |fn| changes
  near value.
- `roots.hpp`, `Criteria`: "about 2e9 in a Boost.Multiprecision type" is `cpp_bin_float_50`'s (`cpp_dec_float_50`:
  about 2e8).

Deferred from the reader: a table of call shapes (fn, fn and f', fn, f' and f'' against a bracket or a guess);
printing a Boost value (`str(50)`); a complete CMake and `main.cpp` (with the install section, at the first release).

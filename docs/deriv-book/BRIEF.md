# Brief: the approved approach and the starting material

Written before the work started, from the author's request and decisions (2026-10-23). Stage 1 turns this into
`PLAN.md`; where they differ, the approved `PLAN.md` wins.

## The request

A textbook-style document for `nxx::deriv` that describes the methods and how to use them and serves as a reference;
for readers comfortable with calculus but not experts (practising engineers); in LaTeX, in the layout of the EquiCore
design paper; with sidebars for complex topics a user does not need; with plots where they make a point; written by
several authors, an editor, and reviewers who mix the target audience and peer reviewers; solid enough to withstand
scrutiny. Decisions: see `CLAUDE.md`.

## What the library offers (state at commit ae17769; check against the code)

- `nxx::deriv::derivative<M, N>(fn, x[, h])`: the N-th derivative, N = 1 to 4 (an integer or `order::first` ...
  `order::fourth`), by method M: `Central<P>` (even P up to 10), `Forward<P>` and `Backward<P>` (P up to 6), or
  `Ridders`. P is the accuracy order, O(h^P), not the number of points.
- Returns `std::expected<Estimate<T>, error>`; `Estimate{value, abs_error}` converts to T; only Ridders fills
  `abs_error`. `error::non_finite` is the one failure. Nothing throws.
- Works for any `nxx::real` T (float, double, long double, Boost.Multiprecision with or without expression templates,
  user types), in constant expressions too.
- Default steps: `eps^(1/(P+N)) * max(|x|, 1)` for a stencil; Ridders starts at `max(|x|, 1) / 32` for f' and `/ 4`
  for f'' to f''''.
- Weights from one Lagrange-formula generator (exact integers); Chapra's and Fornberg's tables reproduced.
- A compile-time warning for a 4th derivative in float; `NUMERIXX_NO_PRECISION_WARNING`.
- Overflow handled before it happens; failures are values even at compile time (one documented exception).

## Draft outline (for stage 1 to refine; page counts are rough)

**Front matter:** title, About this document (preface, how to read, the AI statement), contents, list of sidebars.

**Part I, Fundamentals**
1. Why numerical derivatives (4): where derivatives come from in engineering practice when no formula is at hand;
   the alternatives (analytic, automatic differentiation, complex step) and when each fits; what the library offers;
   a first example.
2. Finite differences from Taylor's theorem (6): forward and central differences; truncation error and the order of
   accuracy; what O(h^2) means in practice.
3. Rounding and the choice of step (7): floating-point numbers in brief; why a smaller step is not always better (the
   V-shaped error curve, as a plot); the best step and the default step; scaling with |x|; what accuracy to expect.
4. Higher-order and one-sided stencils (7): more points for more accuracy; the price in rounding; one-sided rules at
   the edge of a domain; tables of the common rules.
5. Second to fourth derivatives (6): the rules (Chapra's tables), why accuracy falls with the order (a plot), the
   float warning.
6. Ridders' extrapolation (8): Richardson's idea; the tableau; the error estimate and when to trust it; restarts at
   poles and domain edges; cost (plots of convergence and of the estimate against the actual error).
7. Using the library (8): calling it, error handling with `std::expected`, choosing a method (a decision table),
   user and multiprecision types, compile-time use, pitfalls (step scale, aliasing, noisy functions, domain edges),
   worked engineering examples.

**Part II, Reference**
8. API reference (6): every public name, parameter, precondition, return value and failure.
9. Tables (6): weights of every rule for N = 1 to 4, error constants, rounding amplification, default steps.
10. Accuracy, cost and guarantees (5): what accuracy to expect per method, order and type; evaluation counts;
   Ridders' constants and limits; the overflow and constant-evaluation guarantees; how it is verified.

**Appendices:** floating-point essentials; derivations (Taylor expansions, the Lagrange formula and the moment
conditions, Richardson extrapolation); notation; references.

## Sidebar candidates

Machine epsilon, ulps and the rounding of x + h; dividing by the actual distance (Sterbenz's lemma); the error model
behind the default step and why its constant is dropped; the moment conditions and the Lagrange generator; why the
fourth derivative is so sensitive; Burden's difference form of the tableau; the rounding floor and the coverage
statistics of the estimate; aliasing; overflow guards and constant evaluation; expression templates and user types;
why only Ridders estimates its error.

## Figure candidates (each from generated data)

Error against step for the central difference (truncation and rounding branches); the same for `Central<P>`, P = 2 to
10; accuracy against derivative order for `Central<2>` and Ridders; one-sided against central near a domain edge;
Ridders' tableau convergence for one case; estimate against actual error over a random study (coverage); evaluation
counts; float against double.

## Example candidates (each a compiled file in code/)

A first derivative with the default step; error handling; a heat capacity cp = dh/dT from an enthalpy function; a
sensitivity of a model output to a parameter; Newton's method with a numerical derivative; a curvature or an
inflection point with the second derivative; a derivative at the edge of a domain with `Forward<P>`; a multiprecision
derivative; a compile-time derivative with `static_assert`.

## Literature candidates (every one to be verified by the evidence reviewer before it is cited)

- R. L. Burden, J. D. Faires, A. M. Burden, *Numerical Analysis* (numerical differentiation; Richardson
  extrapolation).
- S. C. Chapra, R. P. Canale, *Numerical Methods for Engineers*, ch. 23 (finite-difference formulas up to the 4th
  derivative).
- C. J. F. Ridders, "Accurate computation of F'(x) and F'(x)F''(x)", Advances in Engineering Software, 1982.
- B. Fornberg, "Generation of finite difference formulas on arbitrarily spaced grids", Mathematics of Computation,
  1988.
- L. F. Richardson, on the deferred approach to the limit (1911; and with J. A. Gaunt, 1927).
- N. J. Higham, *Accuracy and Stability of Numerical Algorithms* (floating point, rounding analysis).
- D. Goldberg, "What every computer scientist should know about floating-point arithmetic", ACM Computing Surveys,
  1991.
- A. Sidi, *Practical Extrapolation Methods*, 2003.
- IEEE 754-2019 (floating-point arithmetic).
- On the alternatives, if discussed: J. N. Lyness and C. B. Moler on the complex step (1967); W. Squire and G. Trapp
  (1998); A. Griewank and A. Walther, *Evaluating Derivatives* (automatic differentiation).

## Open questions for stage 1

- The exact chapter split and page targets within 60 to 90 pages.
- Whether Part II repeats the tables of Part I or Part I refers to Part II.
- Which examples run through several chapters (one engineering example as a thread is recommended).
- Whether a short "quick reference" page (a one-page summary of the API and the method choice) belongs at the front
  or the back.

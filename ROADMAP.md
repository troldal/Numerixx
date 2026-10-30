# Roadmap

Status: proposal (2026-10-29), not yet in `docs/plan.md`. The modules and books after `nxx::deriv` and `nxx::roots`,
in the suggested order. Each milestone starts with its open decisions, which are the user's; its phases and checkboxes
go into `docs/plan.md` when it starts. Every module milestone ends with its pages on the documentation site (tutorial,
how-to guides, reference, explanation).

Boundaries that hold for all of it:
- No linear algebra library in Numerixx and no dependency for one (Eigen and Blaze considered, 2026-10-28). The dense
  kernels that least squares and systems need stay internal, implemented from the papers.
- Code only from Numerixx 1.x and Boost.Math; every other method from papers and textbooks.
- Dense problems of modest size; large sparse problems are out of scope.

## Done

- Milestone 1: `nxx::deriv`, derivatives of a scalar function.
- Milestone 2: `nxx::roots`, roots of a real function of one variable.
- Milestone 3: user documentation, the site organised by module (Phases 1 to 5). Open: hosting and versioning of the
  site.

## 4. `nxx::quadrature`: integrals of a function of one variable (current)

Architecture approved: `docs/quadrature/BRIEF.md` (2026-10-30), which refines the list below (infinite intervals for
tanh-sinh only, `GaussKronrod<15>` alone at first, the fixed rule returning a number); the working name
`nxx::integrate` became `nxx::quadrature` there.

- Fixed Gauss-Legendre rules, the nodes as constexpr data.
- Adaptive Gauss-Kronrod, the general-purpose default.
- Romberg (Richardson extrapolation, the counterpart of Ridders).
- Tanh-sinh, for endpoint singularities.
- Exp-sinh and sinh-sinh, for half-infinite and infinite intervals.
- Every result carries an error estimate, in `std::expected`.
- Deferred: oscillatory integrands, Cauchy principal values, integrals in several dimensions.

## 5. `nxx::optimize`: minimum of a function of one variable

- Bracketing a minimum.
- Golden-section search.
- Brent's minimiser.
- Reuses the bracket and tolerance decisions of `nxx::roots`.

## 6. The combined book, part on root finding

- One book on Numerixx with a part per module, each readable on its own.
- The groundwork shared by all parts moves out of the deriv book into an opening part: floating point, rounding,
  error and tolerances, user types, how to read a result.
- The deriv book becomes the part on derivatives; the roots part is written.
- The book agents are generalised to work per part.

## 7. `nxx::interpolate`: interpolation of data in one dimension

- Linear interpolation.
- Cubic splines (natural, clamped, not-a-knot).
- PCHIP and Akima/makima, which keep the shape of the data.
- Barycentric rational interpolation (Floater-Hormann).
- Exact derivatives and integrals of the interpolants.
- Open: the API pattern for data held in `std::span`s, new to the library.
- Deferred: interpolation in two dimensions, general B-splines, smoothing splines (after least squares).

## 8. `nxx::lsq`: least squares

Open before the milestone: where the shared dense kernels live (a module may not use another module's `detail`, so
either a shared internal header or the `detail` of `math.hpp`; a change to `docs/design.md`), and whether linear least
squares is public API.

- Phase 1: the internal kernels (QR with column pivoting, Givens updates, triangular solves, an overflow-safe scaled
  norm) and linear least squares.
- Phase 2: Gauss-Newton with a line search, and Levenberg-Marquardt in More's trust-region form (1978); finite-
  difference Jacobians and gradients added to `nxx::deriv`.
- Phase 3: Powell's dogleg.
- Phase 4: weights, robust losses (Huber, soft-L1, Cauchy; Triggs et al. 2000), and fit statistics (covariance,
  standard errors, rank).
- Phase 5: bound constraints (trust-region reflective, Coleman and Li 1996, or a projected Levenberg-Marquardt).
- Phase 6: the large-residual correction (NL2SOL style; Dennis, Gay and Welsch 1981).
- Test references: the NIST Statistical Reference Datasets (linear and nonlinear regression, certified values) and the
  More-Garbow-Hillstrom test set (1981).
- Sources: More (1978); More, Garbow and Hillstrom, User Guide for MINPACK-1 (ANL-80-74), not the MINPACK code;
  Dennis and Schnabel; Nocedal and Wright; Bjorck, Numerical Methods for Least Squares Problems (2nd ed.).
- Deferred: geodesic acceleration, variable projection, orthogonal distance regression, general constraints.

## 9. `nxx::multiroots`: systems of nonlinear equations

- Newton's method with a line search.
- Broyden's method.
- Powell's hybrid method.
- All on the kernels and Jacobians of Milestone 8.

## 10. `nxx::poly`: polynomials

- `Polynomial<T>`, a value type: Horner evaluation with derivatives and an error bound, arithmetic, division, the
  derivative and the integral.
- Stable closed forms for quadratics and cubics.
- Real roots isolated (Sturm sequences or Descartes' rule) and refined by `nxx::roots`.
- All complex roots by Aberth-Ehrlich, then polished.
- Fitting through linear least squares.
- Open: the complex-type concept for user types, which complex-valued `nxx::deriv` would later reuse.

## Later, as demand arises

- Parts of the combined book on integration, least squares and the other modules, once each module has settled.
- The items deferred in `docs/plan.md`: the open method protocol, `with_noise`, complex derivatives and the rest.

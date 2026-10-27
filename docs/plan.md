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
- Phase 4: `Toms748`, `AndersonBjorck`.
- Phase 5: searching, `expand(fn, x0)`.
- Phase 6, at the user's choice: `Ridders`, `RegulaFalsi`, `Illinois`, `Steffensen`.

## Deferred

### Derivatives (`nxx::deriv`)

- Gradients and Jacobians, until a module needs them.
- For `derivative_of` (`docs/roots/PHASE3_API.md`): a step or scale (a lambda around `derivative(fn, x, h)` does it
  now), and a second call operator that takes fn(x) from the caller, so that Halley's numerical f'' need not evaluate
  it again.

Customization (reviewed 2026-10-25 with a usage sketch and a simplicity review; every item can be added without
breaking the API, so none needs a hook now):
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

## History

- 2026-10-20: the first rewrite was abandoned; master went back to the 1.x code (`v1.1.0-legacy`) and Numerixx 2
  restarted on `dev`.

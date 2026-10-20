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

- [ ] Higher-order central stencils
- [ ] Forward and backward stencils, for points where fn can be evaluated on one side only

### Phase 3: Ridders' extrapolation

- [ ] Ridders' method with an error estimate

### Phase 4: second derivatives

- [ ] Second derivatives of a scalar function

## Deferred

- Gradients and Jacobians, until a module needs them.
- A derivative as a callable (`derivative_of(fn)`), until there is a use for it.

## History

- 2026-10-20: the first rewrite was abandoned; master went back to the 1.x code (`v1.1.0-legacy`) and Numerixx 2
  restarted on `dev`.

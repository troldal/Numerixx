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

API approved: `docs/PHASE2_API.md`.

- [x] Higher-order central stencils, `Central<P>` for even P up to 10 (56cffb1)
- [x] Forward and backward stencils, `Forward<P>` and `Backward<P>` for P up to 6, for points where fn can be evaluated
  on one side only (56cffb1)
- [x] `root_epsilon<N, T>` replaces `cbrt_epsilon` (56cffb1)
- [x] CI passes on all five toolchains (run 39486718887)

### Phase 3: Ridders' extrapolation

API approved: `docs/PHASE3_API.md`.

- [x] `derivative<Ridders>` with an error estimate, restarting below failing or untrusted steps (1023248)
- [x] One return type for every method, `std::expected<Estimate<T>, error>`, implicitly convertible to T (1023248)
- [x] CI passes on all five toolchains (run 39605999468)

### Phase 4: derivatives of order 1 to 4

API approved: `docs/PHASE4_API.md`.

- [x] `derivative<M, N>` for N = 1 to 4, as an integer or an `enum class order`, with every method (d39e31b)
- [x] One Lagrange-formula generator for every table; Phase 2 tables reproduced bit for bit (d39e31b)
- [x] Ridders for N = 1 to 4 on the `Central<2>` rule of each order (d39e31b)
- [x] No overflow before the final quotient (fixes a constant-evaluation failure of `Central<4>` on a constant 1e308)
  (d39e31b)
- [x] A compile-time warning for a float 4th derivative, silenceable at the call or by `NUMERIXX_NO_PRECISION_WARNING`
  (d39e31b)
- [x] CI passes on all five toolchains (run 39910979076)

### Phase 5: one-sided Ridders, and the float warning from f''

API approved: `docs/PHASE5_API.md`.

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

## Deferred

- Gradients and Jacobians, until a module needs them.
- A derivative as a callable (`derivative_of(fn)`), until there is a use for it.

## History

- 2026-10-20: the first rewrite was abandoned; master went back to the 1.x code (`v1.1.0-legacy`) and Numerixx 2
  restarted on `dev`.

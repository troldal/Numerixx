# Migrating from Numerixx 1.x

Numerixx 2 is a rewrite (see [docs/redesign/PLAN.md](docs/redesign/PLAN.md)). The old API is preserved at two tags,
so existing code can stay on it and migrate on its own schedule:

| Tag | Commit | Contents |
|---|---|---|
| `v1.0.0` | `5de1e07` | the last master of Numerixx 1 (C++20) |
| `v1.1.0-legacy` | `8528e94` | the tip of the `dev-reorg` development branch (C++23; reworked optimize module with `fminimize`/`fmaximize`, reworked interpolate module, the `.result()` API, `mdiff`) |

```cmake
CPMAddPackage(NAME Numerixx GITHUB_REPOSITORY troldal/Numerixx GIT_TAG v1.1.0-legacy)
```

Both old trees need their original dependencies:

- `v1.0.0` finds gcem, tl-expected, Blaze, LAPACK and Boost unconditionally (`find_package(... REQUIRED)`), and
  OpenMP with clang-cl.
- `v1.1.0-legacy` needs Boost always, and builds its vendored Blaze, which needs LAPACK and BLAS, whenever
  interpolate, optimize or multiroots is enabled.

## Build and headers

| Numerixx 1.x | Numerixx 2 |
|---|---|
| vcpkg dependencies (Boost, Blaze, LAPACK, gcem, tl-expected, ...) | CPM, fetched automatically; the scalar modules need only the standard library |
| flat include directories: `#include <Roots.hpp>`, `<Deriv.hpp>`, ... | one include root: `#include <numerixx/roots.hpp>`, `<numerixx/deriv.hpp>`, ... |
| targets `numerixx::roots`, `numerixx::deriv`, ...; umbrella `numerixx::all` | targets `numerixx::roots`, `numerixx::deriv`, ...; umbrella `numerixx::numerixx` |
| `v1.1.0-legacy`: modules selected with `NUMERIXX_ROOTS`, `NUMERIXX_DERIV`, ... | the scalar modules are always defined; `NUMERIXX_WITH_LINALG` adds `linalg` and `multiroots` (Eigen), `NUMERIXX_WITH_FXT` adds `pipes` (FXT) |
| `numerixx::func` (`v1.0.0`) | removed; a few function adaptors live in `nxx::fn` (`counted`, `catching`, `extend_linearly`) |
| C++20 (`v1.0.0`) or C++23 (`v1.1.0-legacy`) | C++23 |

## API

The modules arrive phase by phase (see the roadmap in [PLAN.md](docs/redesign/PLAN.md) §8); this table grows with
them.

| Numerixx 1.x | Numerixx 2 |
|---|---|
| `fsolve<Bisection>(f, {lo, hi})` | `roots::bisection{}(f, {lo, hi})`, or the facade `roots::solve(f, {lo, hi})` |
| `fdfsolve<Newton>(f, df, x0)` | `roots::newton{}.with_derivative(df)(f, x0)`; safeguarded: `roots::solve(f, df, x0)` |
| `fdfsolve<Secant>(f, df, x0)` | the derivative-free `roots::secant{}(f, x0)`: no derivative argument |
| the tolerance argument, `fsolve<S>(f, {lo, hi}, eps, maxiter)` and `fdfsolve<S>(f, df, x0, eps, maxiter)` | a criterion that names its test, and `.with_budget(maxiter)`. `v1.1.0-legacy`'s `eps` was both parts of a mixed test (width, or step, ≤ eps·x + eps/2), so `eps = 1e-10` becomes `roots::bisection{nxx::width_tol{5e-11, nxx::rel_tolerance{1e-10}}}` or `roots::newton{nxx::x_tol{5e-11, nxx::rel_tolerance{1e-10}}}` (the width test scales with min(\|lo\|, \|hi\|), not the estimate). `v1.1.0-legacy` multiplied `eps` by the signed estimate, so a bracketing solve of a root below −1/2 never met its tolerance and always ran `maxiter` iterations; Numerixx 2 stops when the width meets the tolerance; `v1.0.0`'s was a residual bound, \|f(x)\| < eps, now the opt-in `nxx::f_tol{eps}`, which does not bound the error in x. One number is absolute (`width_tol{1e-10}`); the relative part is always named, so two bare numbers (`width_tol{a, b}`, `make(a, b)`) do not compile. Run-time values go through `make()`: `width_tol<double>::make(a)`, or `make(a, *rel)` with `rel = rel_tolerance<double>::make(r)`, or `width_tol{*tol, *rel}` from a validated `tolerance<double>` (DESIGN §6.2, §12.24). A number, a validated `tolerance<T>` or a part is not a criterion: `bisection{1e-10}`, `brent{*tol}` and `s.with_stop(1e-10)` do not compile, and the reason names the criterion to wrap it in, `brent{nxx::width_tol{*tol}}` (DESIGN §6.6, §7.2) |
| `.result()`, `.result<T>()` (`v1.1.0-legacy`) | the returned `std::expected` itself: check it (`if (r)`, then `r->x`) or use `value_or`; a failure never arrives as a plain value |
| `.result(fn)` (`v1.1.0-legacy`) | `transform(fn)`, a member of `std::expected` or an FXT pipe |
| `search<...>(f, bounds)` | `roots::expand`, `roots::scan` or `roots::subdivide`; the result is the input of every bracketing solver |
| `fminimize` / `fmaximize` (`v1.1.0-legacy`) | `optimize::minimize` / `optimize::maximize`, returning `extremum{x, fx}` |
| `diff<ALGO>(f, x)` | `deriv::diff(f, x, stencil)`, for example `deriv::diff(f, x, deriv::central_1_4)` |
| `mdiff` (`v1.1.0-legacy`) | `deriv::mixed` |
| `multisolve<MultiNewton>(...)` | `multiroots::newton` or the facade `multiroots::solve(F, x0)` |
| `polysolve(p)` | `poly::roots(p)` |
| `derivativeOf(f)` | `deriv::derivative_of(f)` (`v1.0.0`'s `derivativeOf` discarded `f`; `v1.1.0-legacy`'s keeps it) |
| `integralOf(f)` | `integrate::integral_of(f)` (both 1.x versions discarded `f`) |

There is no compatibility layer: it would have to reproduce behaviour that Numerixx 2 removes on purpose.

## Results that change

Some 1.x results were wrong, so Numerixx 2 returns different values or an error:

- second and mixed derivatives with the default step: `v1.1.0-legacy` used √ε for every stencil, which is O(1)
  wrong for second derivatives and for its `mdiff`; `v1.0.0` used ε^(1/3), which leaves errors around 1e-5 relative.
  Numerixx 2 chooses the step per stencil;
- solvers that ran out of iterations and returned their last iterate as a success now return an error that carries
  the best estimate; `nxx::best(r)` and `nxx::best_x(r)` read the solution or that estimate in one call;
- Newton on a function without a real root (for example x² + 1) could return a non-finite or arbitrary "root"; it now
  fails with the reason and the best estimate;
- quadratics with complex roots no longer return NaN as a success;
- bracketing without a sign change, and poles, are errors instead of "roots". A pole fails with
  `sign_change_not_root`; its best estimate's x locates the pole, and a detected pole's estimate carries no enclosure, so a
  `first_of` chain does not rank it above another solver's estimate (a bracketing failure next to an undetected pole
  still carries its enclosure: DESIGN §6.7 and §6.10, known limits until phase 3);
- a NaN or infinite bracket end, which `v1.1.0-legacy`'s `validateBounds` does not check (it rejects only equal ends;
  `numerixx/roots/impl/Common.impl.hpp` lines 41-47), fails in-band with `non_finite_input`;
  equal ends fail in-band with `invalid_input` instead of throwing `NumerixxError`.

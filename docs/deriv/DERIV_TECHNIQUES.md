# Numerical differentiation: techniques for `nxx::deriv`

Survey of the techniques for the derivative of a scalar function, judged against the project's constraints: any
`nxx::real` (including Boost.Multiprecision with expression templates), only `+ - * /`, comparisons and `nxx::math`,
constexpr, no exceptions, Emscripten, and code only from 1.x and Boost.Math. Recommendations at the end.

## 1. Error model that drives every choice

A first-derivative formula of accuracy order p with step h has

- truncation error `C_p * h^p * |f^(p+1)|`, and
- rounding error `S_p * eps * |f| / h`, where `S_p` is the sum of the absolute stencil weights.

Minimising the sum gives `h* = (S_p * eps / (p * C_p))^(1/(p+1))` (times the scale on which f varies) and a best error
of order `eps^(p/(p+1))`. So a method of order p loses about `D / (p+1)` of the type's D significant digits:

| order p | evals (central) | `C_p` | `S_p` | digits lost, double (D = 15.7) | digits lost, 50-digit type |
|---------|-----------------|--------|-------|-------------------------------|----------------------------|
| 2       | 2               | 1/6    | 1     | 5.2                           | 16.7                       |
| 4       | 4               | 1/30   | 3/2   | 3.1                           | 10.0                       |
| 6       | 6               | 1/140  | 11/6  | 2.2                           | 7.1                        |
| 8       | 8               | 1/630  | 2.08  | 1.7                           | 5.6                        |

The constants are the Taylor remainders of the standard central stencils. For p = 2 and p = 4 they reproduce the
steps Boost.Math uses, `(3 eps)^(1/3)` and `(11.25 eps)^(1/5)`.

Two consequences:
- For `double`, going past order 4 or 6 gains less than one digit. For multiprecision types the gain grows with D,
  which is a general argument for offering higher orders in a library aimed at user types.
- Higher order is not always better: it needs f smooth on the scale `h*`, which grows with p (the Boost.Math docs
  give an example where order 2 gives 5 digits and orders 4 to 8 give none). That is the case for an error estimate.

## 2. Techniques

### 2.1 Fixed finite-difference stencils (central, forward, backward)

Weights `w_k` on points `x + k h` such that `sum_k w_k k^j = m! * delta(j, m)` for j below the stencil size, where m is
the derivative order. The weights are rational; with a common denominator they are small integers, as the current
Phase 1 code and Boost.Math already use.

- **Fit:** only `+ - * /`; constexpr; any `nxx::real`; no new callback requirements. Integer numerators over an
  integer denominator are exact in every type, so `T(1)` style construction from `int` is all that is needed.
- **One-sided stencils** are needed where fn exists on one side only (a domain boundary, a phase boundary). A
  first-derivative forward stencil of order p needs p + 1 points, including `f(x)`, and has larger `C_p` and `S_p`
  than the central one (forward order 2 is `(-3 f0 + 4 f1 - f2) / 2h`, `S = 4` against `1`), so its default step and
  error are worse. Open variants skip `f(x)` (for example the 1.x rule on `x + h, ..., x + 4h`), for fn that cannot be
  evaluated at x itself.
- **Second derivatives** come from the same moment conditions with m = 2: `(f(x+h) - 2 f(x) + f(x-h)) / h^2` has
  truncation `h^2 f''''/12`, rounding `4 eps |f| / h^2`, best step about `eps^(1/4)` and best error about
  `eps^(1/2)`. The current divisor decision generalises: with `hp = (x+h) - x` and `hm = x - (x-h)`, use the
  non-uniform three-point formula `2 ((fp - f0)/hp - (f0 - fm)/hm) / (hp + hm)`.
- **Default step:** `eps^(1/(p+1)) * max(|x|, 1)`, optionally with the constant `(S_p / (p C_p))^(1/(p+1))`. This
  needs an n-th root of epsilon in `nxx::math`; `cbrt_epsilon` generalises directly (same exact radix scaling, Newton
  iteration for `y^n = a`).
- **Weights, tables or algorithm:** Fornberg (1988) computes weights for any derivative order on any node set with a
  short recursion using only `+ - * /`, so it could even run in constexpr. For a handful of fixed stencils, integer
  tables are simpler and exact; the moment conditions above give a `static_assert` check on every table that is
  independent of Numerixx (a closed form, as the test rules require). Fornberg becomes worthwhile only if arbitrary
  orders or non-uniform nodes are wanted.
- **Licensing:** the weights are textbook facts. Boost.Math's `finite_difference.hpp` (BSL-1.0) has central orders
  1, 2, 4, 6, 8 with error estimates, but it calls `std::pow`, `sqrt` and `nextafter`, is not constexpr, and has no
  one-sided or second-derivative stencils, so it is a reference for coefficients and error formulas rather than code
  to copy.

### 2.2 Richardson extrapolation and Ridders' method

The central difference `D(h)` has an error expansion in even powers of h. Richardson extrapolation combines `D(h)`,
`D(h/c)`, ... to cancel those terms. With a fixed number of levels it is just a higher-order stencil (one level with
c = 2 is the 4-point stencil on `+-h/2, +-h`), so on its own it adds nothing to 2.1.

Ridders (1982) makes it adaptive: start from a fairly large step, shrink it geometrically, extrapolate each new
`D(h)` through a Neville tableau, and take as the error estimate the change between neighbouring tableau entries.
Stop when the estimate starts to grow (rounding has taken over) and return the entry with the smallest estimate.

- **Fit:** only `+ - * /` and comparisons; a fixed-size tableau (`std::array`) keeps it constexpr and heap-free; works
  for first and second derivatives, since both central formulas have even-power expansions.
- **Value:** the only technique here that gives a usable error estimate and tolerates a poor initial step, which
  addresses the failure mode in section 1. Cost is typically 6 to 12 evaluations.
- **API consequence:** the result carries an error estimate, so `std::expected<result, error>` with a value and an
  error, and possibly a new error such as "no estimate below the initial one".
- **Licensing:** the widely copied implementation is Numerical Recipes' `dfridr`, which is excluded, along with its
  constants (shrink factor and stopping safety factor). Implement from Ridders' paper and the general Richardson
  theory (Burden et al., *Numerical Analysis*; Sidi, *Practical Extrapolation Methods*, 2003), and choose and justify
  our own constants.

### 2.3 Adaptive step selection without extrapolation

Algorithms that probe fn to estimate the scale on which it varies, then pick the step from the error model:
Stepleman and Winarsky (1979) for central differences; Gill, Murray, Saunders and Wright (1983) for forward
differences in optimisation. They fix the weak point of the default step (the assumption that fn varies on the scale
`max(|x|, 1)`) at a cost of extra evaluations, but give no error estimate. Ridders covers the same need and gives an
estimate, so these are deferred unless Ridders proves too expensive.

### 2.4 Noise-aware differences

When fn is computed by an iterative solver or a simulation, its error is far above `eps` and the step must be set from
the noise level. Moré and Wild estimate the noise (*Estimating computational noise*, SIAM J. Sci. Comput. 2011) and
derive near-optimal steps from it (*Estimating derivatives of noisy simulations*, ACM TOMS 2012). A cheap first step
in that direction is to let the caller state fn's relative accuracy (instead of `eps`) in the default-step formula.
Worth remembering for any client whose fn wraps a solver; deferred.

### 2.5 Complex step

`f'(x) ~= Im f(x + i h) / h` (Lyness and Moler 1967; Squire and Trapp 1998). No subtraction, so h can be tiny and
the result is accurate to about `eps` with one evaluation. Boost.Math has `complex_step_derivative`.

- **Fit, poor:** fn must accept a complex argument and be the restriction of an analytic function. That excludes
  `abs`, `max`, comparisons on the argument and most iterative solvers inside fn, and it needs a complex type for T:
  `std::complex<T>` is unspecified for types other than the standard floating types, so multiprecision needs Boost's
  complex adaptors. It also changes the callback contract from `fn(T)` to a generic fn.
- **Dominated** by forward-mode automatic differentiation (Martins, Sturdza and Alonso 2003 show the method is a form
  of it): if fn has to be generic anyway, dual numbers give the exact derivative and handle comparisons.
- Deferred; not recommended for `nxx::deriv`.

### 2.6 Automatic differentiation

Forward mode with dual numbers (or reverse mode) gives derivatives exact to rounding, but it is a different module:
it needs fn templated on its argument type and a number type with overloads of every elementary function. Boost.Math
already provides `autodiff` (BSL-1.0) for users who want it. Out of scope for `nxx::deriv`; a candidate for a
separate module only if a concrete use appears.

### 2.7 Contour integrals and Taylor coefficients for analytic functions

The Cauchy integral on a circle, discretised by the trapezoidal rule (Lyness and Moler 1967; Fornberg, ACM TOMS 1981;
Bornemann, Found. Comput. Math. 2011), gives many derivatives at once and high derivatives stably. Same analyticity
and complex-argument requirements as the complex step, plus `sin`/`cos` of the nodes and radius selection, which
`nxx::math` lacks. Deferred.

### 2.8 Global approximation and sampled data

Chebyshev interpolation differentiated term by term (Boost.Math `chebyshev_transform::prime`) is excellent for smooth
fn on an interval but belongs to an approximation module. Differentiating sampled or noisy data (Savitzky-Golay,
Lanczos smoothing differentiators, Boost.Math `discrete_lanczos_derivative`, splines) is a different problem: the
input is data, not a callable. Both deferred to their own modules.

## 3. Licensing notes from the legacy code

- 1.x `Order1CentralRichardson` and the 3- and 5-point stencils are standard textbook formulas (cited from Chapra) and
  are safe.
- 1.x `Order1ForwardRichardson` (and its backward twin) evaluates at `x + h, ..., x + 4h` without `f(x)`, with the
  coefficients 22, 62 and 52 over 12. Four nodes and accuracy order 3 give four moment conditions with a unique
  solution, so these weights follow from the node choice, as Richardson extrapolation in Burden et al., *Numerical
  Analysis*, derives them; they can be carried over.
- Ridders: implement from the paper, not from Numerical Recipes (2.2).

## 4. Recommendation

The existing plan already has the right shape; this survey mainly confirms it and fills in details.

1. **Phase 2, stencils (2.1).** Central first-derivative stencils of orders 2, 4, 6 (8 is cheap to add but gains
   little for `double`; include it only if multiprecision use justifies it) and forward/backward stencils, as
   integer weight tables checked by `static_assert` against the moment conditions, with one evaluator. Add an n-th
   root of epsilon to `nxx::math` for the default steps. The 1.x one-sided rules that skip `f(x)` are worth keeping as
   an option for when fn cannot be evaluated at x itself.
2. **Phase 3, Ridders (2.2)** with an error estimate, from the paper, with our own justified constants. This is the
   method to recommend to users who do not know their function's scale.
3. **Phase 4, second derivatives (2.1)** from the same tables (m = 2), with the non-uniform three-point form for the
   basic case. Ridders extends to it with a different starting formula.

Deferred, with the reason:
- Adaptive step selection (2.3): Ridders covers it with an error estimate.
- Noise-aware steps (2.4): until a fn with known computational noise needs it; a caller-stated accuracy is the cheap
  version.
- Fornberg's general weights (2.1): until arbitrary orders or non-uniform nodes are needed.
- Complex step (2.5), automatic differentiation (2.6), contour methods (2.7): they need a generic, analytic fn and
  complex or dual number types; not a fit for a `fn(T)` callback.
- Chebyshev and sampled-data differentiation (2.8): separate modules.

## Sources

- Fornberg, *Generation of finite difference formulas on arbitrarily spaced grids*, Math. Comp. 51 (1988) 699-706,
  <https://www.ams.org/mcom/1988-51-184/S0025-5718-1988-0935077-0/>
- Ridders, *Accurate computation of F'(x) and F'(x)F''(x)*, Adv. Eng. Software 4(2) (1982) 75-76.
- Stepleman and Winarsky, *Adaptive numerical differentiation*, Math. Comp. 33 (1979) 1257-1264,
  <https://www.ams.org/journals/mcom/1979-33-148/S0025-5718-1979-0537969-8>
- Gill, Murray, Saunders and Wright, *Computing forward-difference intervals for numerical optimization*, SIAM J. Sci.
  Stat. Comput. 4 (1983) 310-321 (cited in
  <https://www.nag.com/numeric/nl/nagdoc_27.2/flhtml/e04/e04xaf.html>).
- Moré and Wild, *Estimating derivatives of noisy simulations*, ACM TOMS 38(3) (2012),
  <https://optimization-online.org/2010/11/2787/>
- Lyness and Moler, *Numerical differentiation of analytic functions*, SIAM J. Numer. Anal. 4 (1967) 202-210.
- Squire and Trapp, *Using complex variables to estimate derivatives of real functions*, SIAM Review 40 (1998)
  110-112.
- Martins, Sturdza and Alonso, *The complex-step derivative approximation*, ACM TOMS 29(3) (2003) 245-262,
  <https://mdolab.engin.umich.edu/bibliography/Martins2003a.html>
- Fornberg, *Numerical differentiation of analytic functions*, ACM TOMS 7(4) (1981) 512-526.
- Bornemann, *Accuracy and stability of computing high-order derivatives of analytic functions by Cauchy integrals*,
  Found. Comput. Math. 11 (2011), <https://arxiv.org/abs/0910.1841>
- Boost.Math numerical differentiation docs,
  <https://www.boost.org/doc/libs/latest/libs/math/doc/html/math_toolkit/diff.html>, and source,
  <https://github.com/boostorg/math/blob/develop/include/boost/math/differentiation/finite_difference.hpp>

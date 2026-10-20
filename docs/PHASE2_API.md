# Phase 2 API proposal: stencils

Status: revised after review (simplicity-critic, usage-sketcher, numerics-reviewer, spike). Nothing in `include/` has
changed.

## Summary

Each stencil is a named constexpr value, passed to `derivative` as a non-type template parameter that defaults to the
Phase 1 central difference, so existing calls do not change:

```cpp
using namespace nxx::deriv;
auto const f = [](double x) { return std::exp(x); };

derivative(f, 1.0);                 // Phase 1, unchanged: central2, default step
derivative<central4>(f, 1.0);       // 4 evaluations, error ~ 1.5 eps^(4/5)
derivative<forward2>(f, 0.0, 1e-5); // never evaluates fn below x
```

## Public API (`deriv.hpp`)

```cpp
namespace nxx::deriv {

    // Central rules of accuracy order 2, 4 and 6: p evaluations.
    inline constexpr detail::stencil<2> central2 = /* nodes {-1, 1},             weights {-1, 1}, den 2   */;
    inline constexpr detail::stencil<4> central4 = /* nodes {-2, -1, 1, 2},      weights {1, -8, 8, -1}, den 12 */;
    inline constexpr detail::stencil<6> central6 = /* nodes {-3, ..., 3} \ {0},  weights {-1, 9, -45, 45, -9, 1}, den 60 */;

    // One-sided rules of accuracy order 2 and 4, for fn defined on one side of x only: p + 1 evaluations, including
    // fn(x). forward evaluates at x + k h, backward at x - k h, k >= 0.
    inline constexpr detail::stencil<3> forward2  = /* nodes {0, 1, 2},          weights {-3, 4, -1}, den 2 */;
    inline constexpr detail::stencil<5> forward4  = /* nodes {0, 1, 2, 3, 4},    weights {-25, 48, -36, 16, -3}, den 12 */;
    inline constexpr detail::stencil<3> backward2 = /* nodes {-2, -1, 0},        weights {1, -4, 3}, den 2 */;
    inline constexpr detail::stencil<5> backward4 = /* nodes {-4, -3, -2, -1, 0}, weights {3, -16, 36, -48, 25}, den 12 */;

    /// The derivative of fn at x by stencil S with step h. Preconditions: h > 0, and h well above ulp(x).
    template<auto S = central2, real T, typename F>
    requires detail::is_stencil<decltype(S)> && /* F as in Phase 1 */
    [[nodiscard]] constexpr std::expected<T, error> derivative(F &&fn, T const &x, std::type_identity_t<T> const &h);

    /// The same with the default step eps^(1/(p+1)) * max(|x|, 1), p the accuracy order of S.
    template<auto S = central2, real T, typename F>
    requires /* as above */
    [[nodiscard]] constexpr std::expected<T, error> derivative(F &&fn, T const &x);

} // namespace nxx::deriv
```

The table type is internal:

```cpp
namespace nxx::deriv::detail {
    // f'(x) ~= sum_k weights[k] * fn(x + nodes[k] * h) / (denominator * h), with error O(h^order).
    template<std::size_t N>
    struct stencil {
        std::array<int, N> nodes;
        std::array<int, N> weights;
        int denominator;
        int order;
    };
}
```

`error` keeps its one value, `non_finite`; its documentation now also covers a node `x + k h` that is not finite.

## Public API (`math.hpp`)

```cpp
/// The n-th root of machine epsilon (N >= 2): exact scaling by powers of radix^N, then Newton's iteration from 1,
/// which decreases monotonically to the root.
template<int N, real T>
constexpr T root_epsilon();
```

It replaces `cbrt_epsilon<T>()`, which is removed (nothing is released); its callers become `root_epsilon<3, T>()`.
The default steps use N = 3, 5 and 7.

## Behaviour

- **Divisor.** Two-node rules keep the Phase 1 formula, `(f_hi - f_lo) / (x_hi - x_lo)` with the points actually
  evaluated, so `central2` is unchanged bit for bit. Rules with more nodes divide `sum_k w_k f_k` by `T(den) * h` in a
  single division: dividing by the outer-node distance would cancel the outer-node rounding only for two-node rules
  and enlarges it for the others (by 4x for `central4`, 11x for `central6`). Node rounding then adds about
  `S_p / 2 * eps * |x f'| / h`, the same order as the rounding of fn itself.
- **Failures before arithmetic**, as Phase 1 now does (`design.md`): the input check generalises to
  `|x| + K h <= max` with K the largest |node| plus one (h <= max / 2K), so every node and the divisor are finite and fn
  is only called at finite points; fn's values are checked before they are combined; the final result is checked last.
- **Arithmetic.** Weights enter as `T(w)`; fn values go through `static_cast<T>` as in Phase 1; the sum accumulates in
  a `T`. Only `+ - * /`.
- **Default step** `root_epsilon<p + 1, T>() * max(|x|, 1)`, without the order-dependent constant
  `(S_p / (p C_p))^(1/(p+1))`, as in Phase 1. The constants lie between 1.4 and 1.8 for these rules and cost at most a
  factor 1.5 in error.
- **One-sided rules** never call fn on the other side of x: h > 0, every node has one sign, and rounding is monotone.

## Tests

- **Tables** (`static_assert`, with a moment helper in the test file, independent of the header): for each named
  stencil, `M_0 = 0`, `M_1 = den`, `M_j = 0` for `2 <= j <= order` and `M_(order+1) != 0`, where
  `M_j = sum_k w_k * n_k^j`.
- **Exactness** (`static_assert`): a rule of order p is exact for polynomials of degree <= p; with dyadic x and
  power-of-two h the result is exact in binary floating point. Degree p + 1 gives the exact truncation term, for
  example `-4 h^4` for `central4` on `x^5`.
- **Run time**, per stencil: exp and sin at a few points against `std::exp`/`std::cos`. Tolerances come from the
  error model including the rounding factor `S_p`: about `1.5 eps^(4/5)` for `central4`, `1.8 eps^(6/7)` for
  `central6`, `4 eps^(2/3)` for `forward2` and `(32/3) eps^(4/5)` for `forward4`. A one-sided case where fn is NaN on
  the excluded side.
- **Errors:** `non_finite` from fn, and from a node that overflows (x near `max`) with a saturating fn.
- **Multiprecision:** `cpp_bin_float_50` with and without expression templates, `cpp_dec_float_50`.
- **`root_epsilon`:** the existing `cbrt_epsilon` tests, ported, plus N = 5 and 7 against `pow`.

## Breaking changes

- `cbrt_epsilon` is removed in favour of `root_epsilon<3, T>`.
- An explicit `derivative<double>(fn, x)` no longer compiles, since the first template parameter is now the stencil.
  No code in the repository does this.

## Return type, and how Phase 3 fits

Fixed stencils keep returning `std::expected<T, error>`; the error-estimating method of Phase 3 gets its own function
and result type, so adding it breaks nothing:

```cpp
derivative<central4>(fn, x)    -> std::expected<T, error>              // Phase 2, unchanged by Phase 3
ridders(fn, x[, h0])           -> std::expected<estimate<T>, error>    // Phase 3 (sketch, to be proposed then)

template<real T>
struct estimate {
    T value;
    T abs_error; // not `error`, which is already the enum in nxx::deriv
};
```

Why not return the struct from `derivative` now:
- A fixed stencil has no honest value for the error: bounding its truncation error costs extra evaluations (Boost.Math
  spends one or two), and the rounding part alone is smallest exactly when truncation is largest, so it misleads.
- The struct would have no current use in Phase 2 (two-uses rule), and every caller would write `d->value` instead of
  `*d`, also inside `transform` and `and_then` chains.
- Ridders is not a stencil: it adapts the step, takes an initial step and possibly a tolerance, and its number of
  evaluations varies, so it does not fit `derivative<S>` as one more template argument.

If a second module later returns a value with an error estimate (numerical integration, say), `estimate<T>` moves to
the shared core, with an alias left in `nxx::deriv` so the move is non-breaking.

## Decisions from the review

- **The table type is internal** (simplicity-critic): no current use needs user-written stencils. The moment checks
  live in the tests, and each table stores its order instead of the header computing it.
- **No order-1 rules** (simplicity-critic): `forward1` only pays when fn(x) is already known, which this API cannot
  accept; `forward2` costs one more call and gives `eps^(2/3)` instead of `eps^(1/2)`.
- **Backward tables are literals**, not built by a mirror helper; the moment tests catch typos.
- **`template<auto S = central2, ...>`** with `requires detail::is_stencil<decltype(S)>`: the spike built it on
  clang-cl, `cl`, GCC 16 and Emscripten, at run time and in constant expressions, with `double` and `cpp_bin_float_50`
  with expression templates. A non-stencil argument is rejected by the constraint, not by a hard error.
  `template<detail::stencil S>` also works everywhere, but needs the internal type in the public signature.

## Alternatives considered

- **Stencil as a function argument**, `derivative(fn, x, central4, h)`: puts h third in one overload and fourth in
  the other, or forces the stencil last and loses the defaulted step.
- **Separate functions** (`central`, `forward`, `backward` with an order parameter): three entry points instead of
  one.
- **Weights computed by Fornberg's algorithm**: unnecessary for seven fixed tables (`DERIV_TECHNIQUES.md` 2.1).

## Deferred

All non-breaking:
- A public stencil type for user-written rules, with in-library checks.
- Choosing a stencil at run time (an enum and a dispatcher); today it takes a `switch` over `derivative<S>` calls.
- Distinguishing errors (fn outside its domain, overflow); today all are `non_finite`.
- An error estimate for fixed stencils at the cost of extra evaluations (as in Boost.Math), as a separate function
  returning `estimate<T>`, if Ridders proves too expensive for some use.
- `forward1`/`backward1`, together with an overload that accepts a known fn(x).
- `central8`, open one-sided rules that skip fn(x), one-sided orders 3 and 5 and up.
- The default-step constant per stencil (changes results, so tolerances need rechecking).
- A default step aware of a domain boundary; near one, the caller passes h.
- Second derivatives (Phase 4): a derivative-order field in the internal table.

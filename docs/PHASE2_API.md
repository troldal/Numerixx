# Phase 2 API proposal: stencils

Status: approved on 2026-10-20. Second revision: the first was reviewed by simplicity-critic, usage-sketcher,
numerics-reviewer and a compile spike; this one replaces the named constants (`central4`, `forward2`) by stencil
families `Central<P>`, `Forward<P>` and `Backward<P>`, whose weights are computed from closed forms.

## Naming

Functions, concepts and enumerations follow the standard library (`derivative`, `isfinite`, `real`, `error`,
`error::non_finite`); classes, that is types that can have member functions, are PascalCase (`Central<P>`,
`Estimate<T>`).

## Summary

A stencil is a type that names a family and an accuracy order. It is the first template argument of `derivative`, and
defaults to the Phase 1 central difference, so existing calls do not change:

```cpp
using namespace nxx::deriv;
auto const f = [](double x) { return std::exp(x); };

derivative(f, 1.0);                     // Phase 1, unchanged: Central<2>, default step
derivative<Central<4>>(f, 1.0);         // 4 evaluations, error ~ 1.5 eps^(4/5)
derivative<Forward<2>>(f, 0.0, 1e-5);   // never evaluates fn below x
```

## Public API (`deriv.hpp`)

```cpp
namespace nxx::deriv {

    /// Central rule of accuracy order P (even, 2 <= P <= 10) on the nodes -P/2 ... P/2 without 0: P evaluations.
    template<int P> struct Central;

    /// One-sided rules of accuracy order P (1 <= P <= 6), for fn defined on one side of x only: P + 1 evaluations,
    /// including fn(x). Forward evaluates at x + k h, Backward at x - k h, k = 0 ... P.
    template<int P> struct Forward;
    template<int P> struct Backward;

    /// The derivative of fn at x by stencil S with step h. Preconditions: h > 0, and h well above ulp(x).
    template<typename S = Central<2>, real T, typename F>
    requires detail::stencil<S> && /* F as in Phase 1 */
    [[nodiscard]] constexpr std::expected<T, error> derivative(F &&fn, T const &x, std::type_identity_t<T> const &h);

    /// The same with the default step eps^(1/(P+1)) * max(|x|, 1).
    template<typename S = Central<2>, real T, typename F>
    requires /* as above */
    [[nodiscard]] constexpr std::expected<T, error> derivative(F &&fn, T const &x);

} // namespace nxx::deriv
```

The three class templates are empty tags: a `static_assert` on P is all they contain. Their tables live in `detail`,
where nothing outside the library can depend on them:

```cpp
namespace nxx::deriv::detail {
    // f'(x) ~= sum_k numerators[k] * fn(x + nodes[k] * h) / (denominator * h), with error O(h^order).
    template<std::size_t N>
    struct Table {
        std::array<int, N> nodes;
        std::array<int, N> numerators;
        int denominator;
        int order;
    };

    template<typename S> inline constexpr auto table = /* computed for Central<P>, Forward<P>, Backward<P> */;
    template<typename S> concept stencil = /* S is a specialisation of Central, Forward or Backward */;
}
```

## Public API (`math.hpp`)

```cpp
/// The n-th root of machine epsilon (N >= 2): exact scaling by powers of radix^N, then Newton's iteration from 1,
/// which decreases monotonically to the root.
template<int N, real T>
constexpr T root_epsilon();
```

It replaces `cbrt_epsilon<T>()`, which is removed (nothing is released); its callers become `root_epsilon<3, T>()`.
The default steps use N = P + 1: 3 to 11 for the central rules, 2 to 7 for the one-sided ones.

## The weights

The tables are computed at compile time from closed forms, in exact rational arithmetic (`long long` and `std::gcd`),
then put over one common denominator:

- **Central, P = 2m**, nodes +-1 ... +-m: the coefficient of `fn(x + k h) - fn(x - k h)` is
  `(-1)^(k+1) (m!)^2 / (k (m-k)! (m+k)!)`, k = 1 ... m.
- **Forward, P**, nodes 0 ... P: `w_0 = -(1 + 1/2 + ... + 1/P)`, and `w_k = (-1)^(k+1) C(P, k) / k`, k = 1 ... P.
- **Backward**: Forward mirrored, nodes and weights negated.

Checked in exact rational arithmetic for every allowed P: each satisfies the moment conditions `M_0 = 0`, `M_1 = 1`,
`M_j = 0` for 2 <= j <= P, `M_(P+1) != 0`, and they reproduce the tables of the first revision. The results:

| Stencil      | Numerators / denominator                       | S = sum abs(w) | Error constant C |
|--------------|------------------------------------------------|----------------|------------------|
| `Central<2>` | -1, 1 / 2                                      | 1              | 1/6              |
| `Central<4>` | 1, -8, 8, -1 / 12                              | 1.5            | 1/30             |
| `Central<6>` | -1, 9, -45, 45, -9, 1 / 60                     | 1.83           | 1/140            |
| `Central<8>` | 3, -32, 168, -672, 672, -168, 32, -3 / 840     | 2.08           | 1/630            |
| `Central<10>`| -2, 25, -150, 600, -2100, ... (antisym.) / 2520| 2.28           | 1/2772           |
| `Forward<1>` | -1, 1 / 1                                      | 2              | 1/2              |
| `Forward<2>` | -3, 4, -1 / 2                                  | 4              | 1/3              |
| `Forward<3>` | -11, 18, -9, 2 / 6                             | 6.67           | 1/4              |
| `Forward<4>` | -25, 48, -36, 16, -3 / 12                      | 10.7           | 1/5              |
| `Forward<5>` | -137, 300, -300, 200, -75, 12 / 60             | 17.1           | 1/6              |
| `Forward<6>` | -147, 360, -450, 400, -225, 72, -10 / 60       | 27.7           | 1/7              |

The truncation error is `C h^P f^(P+1)`; the rounding error is about `S eps |f| / h`.

**Why the caps (Central P <= 10, one-sided P <= 6).** The integers are not the limit: the largest numerator is 2100,
exact even in `float`, and the intermediate factorials (up to 20!) fit in `long long`. The caps are numerical:
- `Central<10>` loses about D/11 of a type's D digits; beyond it the gain is under half a digit even for a 50-digit
  type, while fn must be smooth over a wider and wider stencil.
- One-sided rounding amplification grows fast (S = 27.7 for P = 6, against 2.28 for `Central<10>`).

A `static_assert` in the class template rejects other orders with a clear message. Raising a cap is non-breaking.

**Order 1.** `Forward<1>` and `Backward<1>` come for free with the family, so the first revision's reason to drop
them (one more table each) is gone. They cost one evaluation less than order 2 but give `eps^(1/2)` instead of
`eps^(2/3)`; their documentation says so.

## Behaviour

- **Divisor.** Two-node rules (`Central<2>`, `Forward<1>`, `Backward<1>`) keep the Phase 1 formula,
  `(f_hi - f_lo) / (x_hi - x_lo)` with the points actually evaluated, so `Central<2>` is unchanged bit for bit. Rules
  with more nodes divide `sum_k n_k f_k` by `T(den) * h` in a single division: dividing by the outer-node distance
  would cancel the outer-node rounding only for two-node rules and enlarges it for the others (4x for `Central<4>`).
  Node rounding then adds about `S / 2 * eps * |x f'| / h`, the same order as the rounding of fn itself.
- **Failures before arithmetic**, as Phase 1 does (`design.md`): with `Q = max(hi - lo, den)`, where `lo` and `hi`
  are the outermost nodes, the input check is `h <= max / (2 Q)` and `|x| <= max - Q h`. Every node lies within
  `(hi - lo) h <= Q h` of x and the divisor is at most `Q h <= max / 2`, so all are finite and fn is only called at
  finite points; the check itself cannot overflow. For `Central<2>` (Q = 2) this is exactly the Phase 1 check. `x + h == x` (h below the resolution of
  x) is `non_finite`, as `dx == 0` is in Phase 1. fn's values are checked before they are combined; the final result is
  checked last.
- **Arithmetic.** Numerators enter as `T(n)`; fn values go through `static_cast<T>` as in Phase 1; the sum accumulates
  in a `T`. Only `+ - * /`.
- **Default step** `root_epsilon<P + 1, T>() * max(|x|, 1)`, without the order-dependent constant
  `(S / (P C))^(1/(P+1))`, as in Phase 1. The constants K lie between 1.4 and 2.0 for these stencils; leaving one out
  multiplies the best error by `(K^-P + P K) / (P + 1)`, at most 1.64 (for `Central<10>`).
- **One-sided rules** never call fn on the other side of x: h > 0, every node has one sign, and rounding is monotone.

## Tests

- **Tables** (`static_assert`, with a moment helper in the test file, independent of the header): for every allowed P
  of each family, `M_0 = 0`, `M_1 = den`, `M_j = 0` for `2 <= j <= P` and `M_(P+1) != 0`, where
  `M_j = sum_k n_k * node_k^j`; and the order stored in the table is P.
- **Exactness** (`static_assert`): a rule of order P is exact for polynomials of degree <= P; with dyadic x and
  power-of-two h the result is exact in binary floating point. Degree P + 1 gives the exact truncation term, for
  example `-4 h^4` for `Central<4>` on `x^5`.
- **Run time**, per family at several orders (Central 2, 4, 6, 10; Forward and Backward 1, 2, 4, 6): exp and sin at a
  few points against `std::exp`/`std::cos`. Tolerances come from the error model including S, for example about
  `1.5 eps^(4/5)` for `Central<4>` and `10.7 eps^(4/5)` for `Forward<4>`. A one-sided case where fn is NaN on the
  excluded side.
- **Errors:** `non_finite` from fn, from a node that would overflow (x near `max`), and from h below the resolution
  of x, also in constant expressions.
- **Rejected orders:** `Central<3>`, `Central<12>`, `Forward<0>`, `Forward<7>` do not satisfy `detail::stencil` (a
  `requires` check in the tests).
- **Multiprecision:** `cpp_bin_float_50` with and without expression templates, `cpp_dec_float_50`.
- **`root_epsilon`:** the existing `cbrt_epsilon` tests, ported, plus N = 2, 5, 7 and 11 against `pow`.

## Breaking changes

All pre-release:
- `cbrt_epsilon` is removed in favour of `root_epsilon<3, T>`.
- An explicit `derivative<double>(fn, x)` no longer compiles, since the first template parameter is now the stencil.
  No code in the repository does this.

## Return type, and how Phase 3 fits

Fixed stencils keep returning `std::expected<T, error>`; the error-estimating method of Phase 3 gets its own function
and result type, so adding it breaks nothing:

```cpp
derivative<Central<4>>(fn, x)  -> std::expected<T, error>               // Phase 2, unchanged by Phase 3
ridders(fn, x[, h0])           -> std::expected<Estimate<T>, error>     // Phase 3 (sketch, to be proposed then)

template<real T>
struct Estimate {
    T value;
    T abs_error;
};
```

Why not return the struct from `derivative` now:
- A fixed stencil has no honest value for the error: bounding its truncation error costs extra evaluations (Boost.Math
  spends one or two), and the rounding part alone is smallest exactly when truncation is largest, so it misleads.
- The struct would have no current use in Phase 2 (two-uses rule), and every caller would write `d->value` instead of
  `*d`, also inside `transform` and `and_then` chains.
- Ridders is not a stencil: it adapts the step, takes an initial step and possibly a tolerance, and its number of
  evaluations varies, so it does not fit `derivative<S>` as one more template argument.

If a second module later returns a value with an error estimate (numerical integration, say), `Estimate<T>` moves to
the shared core, with an alias left in `nxx::deriv` so the move is non-breaking.

## Decisions from the reviews

- **Families instead of named constants** (second revision): `Central<P>` reads as the parameterised family it is,
  which `central4` only imitated. A family has several current uses, so it passes the two-uses rule, and closed forms
  make the generality cost a few lines. The stencil becomes a type template parameter, which is simpler than the
  class-type non-type parameter of the first revision (the first revision's spike confirmed that one on clang-cl,
  `cl`, GCC 16 and Emscripten; a type parameter needs no such check).
- **The tables are internal** (simplicity-critic): no current use needs user-written stencils. The moment checks live
  in the tests.
- **Divisor and failure checks** (numerics-reviewer): as in Behaviour above.
- **Tolerances include S** (numerics-reviewer): without it the tolerance for `Forward<4>` would be 10x too tight.

## Alternatives considered

- **Named constants** `central4`, `forward2` (first revision): names that mix letters and digits suggest a family
  without being one.
- **Stencil as a function argument**, `derivative(fn, x, Central<4>{}, h)`: puts h third in one overload and fourth in
  the other, or forces the stencil last and loses the defaulted step.
- **Separate functions** (`central<P>`, `forward<P>`, `backward<P>`): three entry points instead of one.
- **Fornberg's general algorithm** for the weights: it handles arbitrary nodes, which nothing needs; the closed forms
  cover the three families with less code.

## Deferred

All non-breaking:
- A public stencil type for user-written rules, with in-library checks.
- Choosing a stencil at run time (an enum and a dispatcher); today it takes a `switch` over `derivative<S>` calls.
- Distinguishing errors (fn outside its domain, overflow); today all are `non_finite`.
- An overload of the one-sided rules that accepts a known fn(x), saving one evaluation.
- Higher orders than the caps; open one-sided rules that skip fn(x).
- The default-step constant per stencil (changes results, so tolerances need rechecking).
- A default step aware of a domain boundary; near one, the caller passes h.
- An error estimate for fixed stencils at the cost of extra evaluations (as in Boost.Math), as a separate function
  returning `Estimate<T>`, if Ridders proves too expensive for some use.
- Second derivatives (Phase 4): a derivative-order parameter on the families, or families of their own, with the
  same evaluator.

# Phase 1 API proposal: the core, the fixed Gauss-Legendre rule, and Romberg

Status: approved by the user on 2026-10-30, Q1 to Q7 as recommended (draft 2). One review round (numerics-reviewer,
simplicity-critic, usage-sketcher) is folded in ("Decisions from the reviews"). Architecture:
`docs/quadrature/BRIEF.md` (approved 2026-10-30, D1 to D14 as recommended).

## Summary

```cpp
#include <numerixx/quadrature.hpp>
namespace qd = nxx::quadrature;

auto const f = [](double x) { return std::exp(x); };

// The fixed rule: a number, smooth in parameters inside fn; no estimate, no criteria.
std::expected<double, qd::error> const g = qd::integrate<qd::GaussLegendre<10>>(f, 0.0, 1.0); // e - 1 to a few ulps

// Romberg: an Integral with an error estimate, or a Failure that may carry the best estimate.
auto const r = qd::integrate<qd::Romberg>(f, 0.0, 1.0);                     // 17 evaluations, error 3.3e-14
if (r) { double const v = r->value; double const e = r->abs_error; }      // e = 3.4e-10, an estimate
auto const s = qd::integrate<qd::Romberg>(f, 0.0, 1.0, {.rel_tol = 1e-14}); // 65 evaluations
auto const t = qd::integrate<qd::Romberg>(f, 1.0, 0.0);                     // oriented: 1 - e

// Compile time: .value(), whose throw on an empty result is a compile error; compared with a tolerance (Q6).
constexpr double nine = qd::integrate<qd::GaussLegendre<7>>([](double x) { return x * x; }, 0.0, 3.0).value();
static_assert(nxx::math::abs(nine - 9.0) <= 4 * 9.0 * std::numeric_limits<double>::epsilon());
```

- A new header, `include/numerixx/quadrature.hpp`, namespace `nxx::quadrature`. It includes `math.hpp` only.
- Two methods, `GaussLegendre<N>` and `Romberg`, one function name with two overloads (BRIEF.md 3.1).
- No default method yet: `GaussKronrod<15>` becomes the default in Phase 2, when it exists (adding a default template
  argument breaks nobody).
- Everything is constexpr, and nothing throws.

## Public API (`quadrature.hpp`)

```cpp
namespace nxx::quadrature {

    /// Why integrate() returned no result. More values may be added, and a listed cause may move to a new value: a
    /// switch should have a default. Formattable: std::format("{}", e) writes the enumerator's name.
    enum class error {
        non_finite,    ///< An end is NaN or infinite (for an infinite interval use TanhSinh, Phase 4), a tolerance
                       ///< is negative or not finite, fn returned a value that is not finite, or a sum overflowed.
        not_converged, ///< The tolerance was not met within the method's cap (Romberg: 16 levels).
    };

    /// An integral and an estimate of its error. The estimate is not a bound: no method that samples fn at finitely
    /// many points can bound the error (a narrow peak between the points is missed). It is never below the rounding
    /// floor (C). Read it by member name: members may be added.
    template<real T>
    struct Integral
    {
        T value;
        T abs_error;
        friend constexpr bool operator==(Integral const &, Integral const &) = default;
    };

    /// Why a method with an estimate returned no Integral, and for not_converged the best estimate it reached.
    /// Formattable: std::format("{}", f) writes the reason's name.
    template<real T>
    struct Failure
    {
        error reason;
        std::optional<Integral<T>> best; ///< present for not_converged, empty for non_finite
        friend constexpr bool operator==(Failure const &, Failure const &) = default;
    };

    /// When to stop: when the error estimate is at most max(abs_tol, rel_tol |value|), or at the rounding floor (C).
    template<real T>
    struct Criteria
    {
        T abs_tol = T(0);                              ///< for integrals near 0
        T rel_tol = nxx::math::root_epsilon<2, T>();   ///< sqrt(eps): 1.5e-8 in double
        friend constexpr bool operator==(Criteria const &, Criteria const &) = default;
    };

    /// The N-point Gauss-Legendre rule, exact for polynomials of degree 2N - 1 up to the rounding of its nodes and
    /// weights. 1 <= N <= 64 (Q1).
    template<int N>
    struct GaussLegendre
    {
    };

    /// Romberg's method: the trapezoidal rule on halved steps, extrapolated by Richardson (Romberg 1955).
    struct Romberg
    {
    };

    /// The fixed rule: N evaluations of fn, no estimate, no criteria.
    template<typename M, real T, typename F>
        requires detail::fixed_rule<M> && nxx::detail::callable<F, T>
    [[nodiscard]] constexpr std::expected<T, error> integrate(F &&fn, T const &a, T const &b);

    /// A method with an estimate (Romberg; Gauss-Kronrod and TanhSinh in later phases): an Integral, or a Failure.
    template<typename M, real T, typename F>
        requires std::same_as<M, Romberg> && nxx::detail::callable<F, T>
    [[nodiscard]] constexpr std::expected<Integral<T>, Failure<T>>
    integrate(F &&fn, T const &a, T const &b, Criteria<std::type_identity_t<T>> const &criteria = {});

} // namespace nxx::quadrature

template<> struct std::formatter<nxx::quadrature::error>;                        // the enumerator's name
template<nxx::real T> struct std::formatter<nxx::quadrature::Failure<T>>;        // the reason's name (Q7)
```

- **Two overloads, one per kind of method,** each constrained to its methods, so `integrate<GaussLegendre<10>>(fn,
  a, b, {.rel_tol = 1e-12})` does not compile (no criteria for the fixed rule), and `integrate<Romberg>` returns a
  `Failure` where the fixed rule returns a bare `error`. The second overload's constraint is `std::same_as<M,
  Romberg>` until Phase 2 adds a second method with an estimate; a named concept waits for it.
- **`Failure`'s formatter writes the reason only** (BRIEF.md 3.3: "formattable through its reason"): it works for
  every T, while writing `best` would need a formatter for T, which Boost's types do not have.
- **fn** as in the other modules: `nxx::detail::callable<F, T>`, invoked with a `T const` lvalue; with the
  expression-template rule of `math.hpp` (`nxx::detail::value_of`: for a class type T, fn returns T or an arithmetic
  type, not an expression), which this module inherits. For a built-in T anything convertible is accepted, so a fn returning
  double for T = float narrows silently, as in the other modules; documented.
- **Both ends deduce T,** so `integrate<M>(fn, 0.0, 1)` does not compile; as roots' `Bracket{1, 2.0}`, which catches
  an integer end and an expression-template end at compile time (the brief's usage sketch counted it a strength).
- **`GaussLegendre<N>` is unconstrained as a type**, as `nxx::deriv::Central<P>`; `detail::fixed_rule` checks
  1 <= N <= 64, so `GaussLegendre<0>` fails at the call with the constraint's name.

## Behaviour

### A. Validation at the boundary, in this order

1. **Tolerances** (Romberg). Negative or not finite: `non_finite`. Then `rel_tol` is clamped to [0, 1], so that
   `rel_tol |value|` stays finite (constant evaluation rejects an overflow); `abs_tol` needs no clamp, since the test is
   a max of finite values. No lower clamp: the floor (C) is the effective minimum.
2. **Ends.** NaN or infinite: `non_finite` (BRIEF.md D13: infinite intervals are `TanhSinh`'s, Phase 4).
3. **a == b:** the result is 0 (fixed rule) or `{0, 0}` (Romberg) without calling fn; 0 is exact, so is its error.
4. **Orientation:** b < a gives the negated integral, through the signed scale (B).

### B. Nodes and scale

- **Nodes:** `c + hw u` for u in (-1, 1), with `c = a / 2 + b / 2` and `hw = b / 2 - a / 2`, which cannot overflow
  ([-max, max] has hw = max) and round once each in radix 2 above the subnormal range. Romberg evaluates its ends as
  fn(a) and fn(b), not as c -+ hw.
- **The scale of each term** is formed before the sum, not applied to it: the fixed rule's term is `s_i fn(x_i)` with
  `s_i = (b - a) (w_i / 2)` where b - a is finite, else `hw w_i`; Romberg's level k adds `h_k fn(x_j)` with
  `h_k = (b - a) / 2^k` (or `2 hw / 2^k`). Two cases the numerics review found: draft 1's `hw (fn(a) + fn(b))`
  overflows the sum of the values on fn = 1e308 over [0, 1e-10] (true integral 1e298); and `hw = (b - a) / 2` rounds
  a subnormal width ([0, 3 denorm_min]: hw = 1.5 denorm_min rounds to 2, a 33 % error in the result), while
  `(b - a)` is exact there and `w_i / 2` exact in radix 2.
- **What no map avoids, documented:** an interval narrow against its position (c far larger than hw) rounds each node
  by about eps |c|, so fn's argument errs by that and the result by about eps |c| times the integral of |f'|; this is
  fn's conditioning, as in roots, not part of the estimate. On an interval a few ulps wide a Gauss node can round onto
  an end, where the fixed rule then evaluates fn.

### C. Sums, the rounding floor and the result's error

- **Value sums are Neumaier-compensated** (`detail::Sum<T>`, shared with Phase 2): Romberg's level k adds 2^(k-1)
  values, where a plain running sum drifts by about n eps L1, and the fixed rule's N terms would err by up to about
  N eps L1, above its weights' contribution (D). The sums of |fn| for L1 are plain: their terms are positive, so they
  err by about n eps relative, against a floor of 50 eps.
- **The rounding floor** is `50 max(eps L1, denorm_min)`, L1 the integral of |fn| by the same rule on the same
  evaluations (Romberg: the trapezoid sums of |fn|; its extrapolation's weights are positive, Bauer, Rutishauser and
  Stiefel 1963, with sum |c_j| below about 1.97, so Romberg's own rounding is at most about 14 eps L1 and 50 is
  conservative for it). QUADPACK's constant (BRIEF.md 3.6), the module's one, shared with Phase 2. The denorm_min term
  ends an integral below about min / eps (1e-320 e^x), where eps L1 and rel_tol |value| both underflow to 0. With the
  floor at 50 eps, a request below about 1e-14 relative is raised to it in double.
- **`abs_error` = max(estimate, floor)**, so an integral that is exactly 0 (x^3 on [-1, 1]) reports the floor, 5.6e-15
  in the prototype, not 0. The fixed rule reports no error.
- **Finiteness is checked once per sum, not per value:** with every scale finite and nonzero, a non-finite value of fn
  makes its sum non-finite (Neumaier's correction too), so each of Romberg's level sums and L1, and the fixed rule's
  result, is checked once: `non_finite`. That covers fn's own NaN and infinities and an overflow of the sum (1e308 sin x
  on [0, 2 pi]) alike, and keeps N branches out of the fixed rule's loop, which runs inside root finders.
- **Samples that are all 0** give {0, 0} even for a positive integral ((x (1 - x) prod_(j=1..15) (16 x - j))^2 on
  [0, 1] vanishes at every point through level 4); an exact-zero result is a sampling artefact, documented with the
  narrow-peak warning.

### D. GaussLegendre<N>

- **The rule:** nodes x_i, the zeros of P_N, and weights `w_i = 2 / ((1 - x_i)(1 + x_i) P_N'(x_i)^2)`; the integral
  `sum s_i fn(c + hw x_i)` (B). N evaluations of fn.
- **The generator** (BRIEF.md 3.4), arithmetic only, with no pi and no cos: the number of sign changes in the
  sequence P_0(x), ..., P_N(x) (a Sturm sequence, from the three-term recurrence; a zero term is skipped and the last
  nonzero sign kept) counts the zeros of P_N strictly above x, so bisection on that count isolates each zero of the
  positive half in a bracket of its own, from [0, 1]; then Newton on P_N, safeguarded by bisection inside the bracket,
  converges to it; x_i = -x_(N+1-i), and 0 for odd N. Three requirements from the numerics review, each a defect of the
  spike's code:
  - the isolation bisection is capped at `digits log2(radix) + 8` steps, after which the bracket is used as it is (a
    type with few digits, such as `std::float16_t`, could otherwise miscount forever);
  - the bracket's orientation comes from the zero's index, not from the computed sign of P_N at its end: P_N > 0 above
    x_1 and its sign on (x_(k+1), x_k) is (-1)^k (a computed P_N(lo) of exactly 0 turned the spike's orientation
    around and duplicated a node);
  - Newton's stagnation stop (the step no longer shrinking) applies only once the step is below sqrt(eps) |x|, and the
    result is asserted strictly inside its bracket (the spike's stop could fire on two growing steps near an extremum
    of P_N).
- **Cost, for double** (a spike, 2026-10-30; smallest passing constexpr limit, clang-cl / cl / GCC 16, against the
  defaults 1,048,576 steps, 1,048,576 steps and 2^25 operations):

  | N | Sturm count (chosen) | Tricomi start (needs pi and cos in T) | Interlacing over n = 2 ... N |
  |---|---|---|---|
  | 10 | 4.9k / 4.2k / 22k | 7.8k / 3.7k / 54k | 14k / 11k / 45k |
  | 30 | 41k / 36k / 237k | 19k / 14k / 76k | 251k / 211k / 929k |
  | 64 | 200k / 189k / 1.2M | 60k / 48k / 237k | 2.2M / 1.8M / 8.1M, over the default |
  | 128 | 878k / 783k / 5.4M | 200k / 168k / 783k | 16M / 13.5M / 63M, over |

  The Sturm count guarantees its brackets, at about N^2 log N; Tricomi's start is cheaper but needs pi and cos in T,
  two helpers for this one use (BRIEF.md 3.10); interlacing over every lower degree is N^3 and stops at N = 48 on
  clang-cl. N_max = 64 keeps every compiler below a fifth of its limit; 128 fits at 84 % of clang-cl's (Q1). The
  phase measures long double (GCC's 80 bits, Emscripten's binary128) against the same limits.
- **Accuracy of the computed rule, Sturm route** (the spike's measurement at 11 values of N, `out_followup2_gcc.txt`;
  the phase measures every N against the tests' independent references): nodes within 0.97 eps relative (double,
  N = 50), 0.89 (float), 1.07 (`cpp_bin_float_50`, N = 64), and within 0.44 eps absolute; weights up to about 360 eps
  relative at the outer nodes (the amplification -2x / (1 - x^2) of a node's rounding, from the Legendre equation),
  but their absolute errors sum to at most 9.4 eps for N <= 64 (1.42 eps at N = 7, 5.7 at N = 20, 9.39 at N = 64, in
  double).
- **Tables:** for `float`, `double` and `long double`, named (not `std::is_floating_point_v`, which also admits the
  unmeasured `std::float16_t` and `std::bfloat16_t`), each rule is a namespace-scope variable template
  `detail::gauss_legendre_table<N, T>`, `inline constexpr`, computed at compile time once per program; for every other
  T the generator runs at the start of each call (about a millisecond for N = 20 in `cpp_bin_float_50`, measured in
  the phase); a cache for those is Phase 2's (BRIEF.md D10).
- **Accuracy:** the truncation error of the N-point rule on [a, b] is `(b - a)^(2N+1) (N!)^4 / ((2N + 1) ((2N)!)^3)
  f^(2N)(xi)` (textbook; N = 1 gives the midpoint rule's (b - a)^3 / 24). The rounding error of the computed rule,
  per node, is about `|hw| sum_i (|dw_i| |f_i| + w_i |f'_i| (|hw| |dx_i| + eps (|c| + |hw x_i|)))` plus the summation's
  and fn's own rounding; draft 1's `16 eps (max |f| + max |f'|) |hw|` was too tight (it omitted the rounding of c and
  of the nodes: fn = x - 1e8 on [1e8, 1e8 + 7 2^-26] with N = 1 errs by 7.8e-16 against an allowed 1.9e-22) and too
  loose (x^127 on [0, 1] at N = 64: 1024 eps allowed, about 1 eps actual). The tests compute it per case from the
  reference nodes and weights. Polynomial results are not exact to the ulp: in the spike GL7 gave 9 - 1 ulp for x^2
  on [0, 3] (on its interlacing table; the Sturm table differs in the last bits), and even Boost's correctly rounded
  table misses x^3 on [0, 1] by half an ulp at N = 7 (Q6).
- **Why a fixed rule has no estimate and no criteria:** its use is inside a root finder, an optimiser or a nested
  integral, where a result that changes smoothly with a parameter in fn matters more than an estimate (Lyness 1983);
  an estimate would cost N - 1 more evaluations on every call. `Romberg` or the default method give one (BRIEF.md D4).

### E. Romberg

- **Levels:** T_0 = (b - a) (fn(a) + fn(b)) / 2, formed as `h_0 fn(a) / 2 + h_0 fn(b) / 2` (B); level k >= 1 adds
  the 2^(k-1) new midpoints u_j = -1 + (2j - 1) / 2^(k-1) (exact: at most 16 bits at the cap), T_k = T_(k-1) / 2 +
  h_k sum fn(c + hw u_j). Two rows of the Richardson table, R(k, m) = R(k, m-1) + (R(k, m-1) - R(k-1, m-1)) /
  (4^m - 1); 4^m - 1 is exact in double and in radix 10 at every m up to the cap, and in float up to m = 12 (beyond,
  the correction is at rounding level).
- **The estimate:** |R(k, k) - R(k-1, k-1)|, the difference of the last two diagonal entries. It estimates the error
  of the previous entry, so it is pessimistic for the returned R(k, k) on smooth fn (30 to 10^4 times the true error
  in the prototype), and it stays honest where Richardson's assumptions fail: for sqrt x, whose expansion has a term
  in h^1.5, it is about 1.83 times the error, while the subdiagonal |R(k, k) - R(k, k-1)| would be 1.83 / (4^k - 1)
  times it, a factor 140 too small at k = 4 (numerics review). The result is R(k, k).
- **Stop:** at level k >= 4 (17 evaluations; 1.x returned 0 for sin^2 on [0, 2 pi] after 2) when
  estimate <= max(abs_tol, rel_tol |R(k, k)|, floor).
- **The cap:** 16 levels, 2^16 + 1 = 65,537 evaluations; then `not_converged` with best = {R(16, 16),
  max(estimate, floor)}. Fixed, not from T: in the prototype smooth fn need at most 9 levels in double and 14 at 50
  digits for 1e-45 (16,385 evaluations); what reaches the cap is a fn Romberg suits badly (sqrt x and a kink at 1/5 at
  rel_tol = 1e-14). Draft 1's 20 levels let such a misuse cost a million evaluations. Raising the cap later turns
  failures into successes and breaks nobody; lowering it would. Q3.
- **No stagnation rule:** noise in fn below the floor ends at the floor, and noise above it is invisible to the
  levels before it shows (the prototype: e^x (1 + 1e-10 sin(1e7 x)) at rel_tol = 1e-12 stops at level 5 reporting
  3.1e-13, true error 2.1e-11). No estimate from samples sees noise it has not sampled, and a rule for estimates that
  stop shrinking would rarely fire on noise that falls slowly from level to level; at a tight tolerance such a fn can
  run to the cap. The documentation says that fn's own accuracy bounds the result's, and that the estimate assumes fn
  accurate to about eps. Q4.
- **Cost:** 2^k + 1 evaluations at level k.

## Prototype results

`scratchpad/romberg/proto.py` (session scratchpad, not kept): Romberg as in E, in double (Python floats), Neumaier
sums, floor 50 eps L1; references from mpmath at 60 digits. Its map (always `b/2 - a/2`, `a/2 + b/2`, the scale
applied after the sum) matches B on every interval below; the phase aligns it before the asserted counts depend on
it.

| Case | rel_tol | Status | Level | Evaluations | Estimate | Actual error |
|---|---|---|---|---|---|---|
| e^x on [0, 1] | 1.5e-8 | ok | 4 | 17 | 3.4e-10 | 3.3e-14 |
| e^x on [0, 1] | 1e-14 | ok | 6 | 65 | 1.9e-14 | 7.7e-17 |
| e^x on [1, 0] | 1.5e-8 | ok, 1 - e | 4 | 17 | 3.4e-10 | 3.3e-14 |
| 1 / (1 + x^2) on [0, 1] | 1.5e-8 | ok | 5 | 33 | 2.9e-9 | 1.2e-11 |
| 1 / (1 + x^2) on [0, 1] | 1e-14 | ok | 8 | 257 | 8.7e-15 | 3.1e-17 |
| sin on [0, pi] | 1.5e-8 | ok | 5 | 33 | 5.4e-9 | 1.3e-12 |
| sin^2 on [0, 2 pi] | 1.5e-8 | ok | 7 | 129 | 5.7e-9 | 1.4e-12 |
| sin^2 on [0, 2 pi] | 1e-14 | ok | 9 | 513 | 3.5e-14 | 1.2e-16 |
| x^2 on [0, 3] | any | ok (the floor) | 4 | 17 | 1.0e-13 | 0 |
| x^3 on [-1, 1] (= 0) | any | ok (the floor) | 4 | 17 | 5.6e-15 | 0 |
| \|x - 1/5\| on [0, 1] | 1.5e-8 | ok | 14 | 16,385 | 2.9e-9 | 5.9e-10 |
| \|x - 1/5\| on [0, 1] | 1e-14 | not_converged | 16 | 65,537 | 1.8e-10 | 3.7e-11 |
| e^(-1e4 (x - 1/2)^2) on [0, 1] | 1.5e-8 | ok | 12 | 4,097 | 1.2e-12 | 2.9e-16 |
| cos(1000 x) on [0, 1] | 1.5e-8 | ok | 14 | 16,385 | 3.4e-14 | 3.5e-18 |
| 1e-300 e^x on [0, 1] | 1.5e-8 | ok | 4 | 17 | 3.4e-310 | 3.3e-314 |
| sqrt x on [0, 1] | 1.5e-8 | ok | 16 | 65,537 | 7.5e-9 | 4.1e-9 |
| sqrt x on [0, 1] | 1e-14 | not_converged | 16 | 65,537 | 7.5e-9 | 4.1e-9 |
| 1 / sqrt x on [0, 1] | any | non_finite (fn(0)) | 0 | 2 | | |
| e^x (1 + 1e-12 sin(1e7 x)) | 1e-14 | ok, false | 6 | 65 | 1.9e-14 | 2.1e-13 |
| e^x (1 + 1e-10 sin(1e7 x)) | 1e-12 | ok, false | 5 | 33 | 3.1e-13 | 2.1e-11 |
| e^x on [0, 1], 50 digits (mpmath) | 1e-45 | ok | 11 | 2,049 | | 2.7e-51 |
| 1 / (1 + x^2), 50 digits | 1e-45 | ok | 14 | 16,385 | | 1.3e-51 |

The two "false" rows are fn noise far above eps (documented, not detectable, and not tested: a test would pin a known
wrong success to a platform's `sin`). Draft 1 had a kink at 1/3, which is degenerate: 2^k / 3 mod 1 alternates
between 1/3 and 2/3, both giving the trapezoid error (2/9) h^2 exactly, which Richardson's first column removes.

## Tests (`tests/test_quadrature_gauss.cpp`, `tests/test_quadrature_romberg.cpp`; test cases `quadrature:...`)

- **References** from closed forms (e - 1, pi / 4, 2, pi, 9, 0, 17/50 for the kink at 1/5, 2/3, sin(1000) / 1000,
  sqrt(pi) erf(50) / 100 for the peak), at 50 digits, cross-checked with mpmath; never from Numerixx.
- **The generator** against independent nodes and weights: Boost.Math's 115-digit tables for N = 7, 10, 15, 20, 25,
  30, and mpmath at 60 digits (Golub-Welsch, `mp.eigsy` on the Jacobi matrix) for a sample of other N (1, 2, 3, 5, 40,
  50, 64), copied into the test with their sources; nodes within a stated number of eps, weights within their
  amplification bound. Without a reference, at every N from 1 to 64: the weights sum to 2, and x^(2N-2) and x^(2N-1)
  on [0, 1] are integrated within the per-node rounding bound of D. Draft 1's 100-digit run of the generator as the
  reference for other N was Numerixx checking itself (the test rule; a systematic fault such as a duplicated node
  would pass).
- **GaussLegendre accuracy:** `|value - ref| <=` the truncation bound of D (from a stated bound on |f^(2N)|) plus the
  per-node rounding bound; float, long double and both 50-digit types on a sample of N.
- **Romberg accuracy:** `|value - ref| <=` the closed-form bound of Bauer, Rutishauser and Stiefel (as Stoer and
  Bulirsch state it), `|I - R(k, k)| <= (b - a)^(2k+3) 2^(-k(k+1)) |B_(2k+2)| / (2k + 2)! max |f^(2k+2)|`, at the
  level k the case stops at (from the prototype), plus about 14 eps L1 of rounding (C). For e^x on [0, 1] at k = 4 it
  is 5.4e-14 against the prototype's 3.3e-14. Where the bound is useless (cos(1000 x), the peak), only the honesty
  check. The honesty check is separate and labelled: `|value - ref| <= abs_error`, the estimate covering the error, on
  every smooth case. Draft 1 used it as the accuracy test, which a Richardson defect losing three digits would pass
  (the estimate is 30 to 10^4 times the error). The evaluation counts exactly in double for fn of `+ - * /` (the
  prototype as the independent implementation, as roots' was), and the stop level.
- **Edges:** reversed ends (negated, both methods); a == b (0, fn never called); NaN and infinite ends; a negative,
  NaN and huge tolerance; abs_tol = 0 with an integral of 0 (the floor ends it); fn not finite at an end and inside;
  1e308 sin x on [0, 2 pi] (`non_finite`); fn = 1e308 on [0, 1e-10] (1e298, no overflow); [-max, max] with fn = 0 and
  fn = 1e-300; [0, 3 denorm_min] with fn = 1e300 (the subnormal width); 1e-320 e^x on [0, 1] (the denorm_min floor
  ends it); sqrt x and the kink at 1e-14 (`not_converged` at 65,537 evaluations, best present and within its own
  estimate); the all-zero samples (documented behaviour, {0, 0}).
- **Types:** double for the edges; float, long double, `cpp_bin_float_50` with and without expression templates and
  `cpp_dec_float_50` (radix 10) for a few accuracy cases with fn of `+ - * /`; `static_assert`s for the built-in
  types only (Boost's are not literal).
- **Constant evaluation:** `static_assert` for N = 1, 2, 7 and 64 (the smallest, the first even, the summary's case,
  the budget limit) within the rounding bound, and Romberg on x^2 over [0, 3] within its floor. The run-time tests of
  every N already put each double table through constant evaluation (`inline constexpr`), on every toolchain.
- **fn's exceptions pass through** unchanged; the count of evaluations is checked with a counting fn.

## Documentation

- `quadrature.hpp`: Doxygen comments for every public name: what an estimate is and is not, the floor and its
  effective minimum, the cost, the noise limit, the all-zero samples, an interval narrow against its position, the
  infinite-interval pointer to Phase 4, the cap, the silent narrowing of a built-in fn result.
- `docs/quadrature/design.md`: the decisions in force, the sources, the prototype table.
- The documentation site: a minimal `docs/sphinx/quadrature/reference.md` with the Breathe directives, since the docs
  build fails unless every public declaration is rendered; it says that `Criteria` here and in `nxx::roots` differ in
  their defaults (sqrt(eps) here, full precision there) and collide under two using-directives. The module's chapter
  (overview, tutorial, guides) at the end of the milestone, as Milestone 3 did for roots and deriv.
- `docs/plan.md`: Phase 1's checkboxes. CLAUDE.md already lists the module (`quadrature:*`).
- No README example and no demo yet: both wait for the default method (Phase 2). Q5.

## Decisions for the user

- **Q1. GaussLegendre's generator and range: the Sturm count, 1 <= N <= 64.** Recommended: yes. Alternatives:
  N <= 128 (84 % of clang-cl's constexpr budget per rule); Tricomi's start with pi and cos internal (cheaper at
  compile time, two more helpers).
- **Q2. Compile-time tables for float, double and long double, named; other T per call (the cache in Phase 2).**
  Recommended: yes (BRIEF.md D10).
- **Q3. Romberg's cap: 16 levels (65,537 evaluations), fixed.** Recommended: yes (changed from draft 1's 20; raising
  it later breaks nobody). Alternative: 20 levels, for 50-digit full precision on less smooth fn.
- **Q4. No stagnation rule for Romberg; fn's noise documented as a limit no estimate sees.** Recommended: yes.
- **Q5. No README example or demo until Phase 2.** Recommended: yes.
- **Q6. No exact-equality checks of quadrature results** in tests, examples or documentation: a tolerance from the
  rounding bound of D; constant-evaluation examples use `.value()`, whose throw on an empty result is a compile
  error, not `*`, which is undefined behaviour there. Recommended: yes.
- **Q7. `Failure` formattable through its reason** (BRIEF.md 3.3). Recommended: yes, a formatter that writes the
  reason's name. Alternative: no formatter; callers format `f.reason` (the simplicity review: no current use).

## Decisions from the reviews

From the numerics review (all adopted):
- The generator's three defects (an uncapped isolation loop, the bracket's orientation from a computed sign, a
  stagnation stop on growing steps): requirements in D.
- The rounding bound: per node, with the rounding of c and of the nodes, replacing a formula that was too tight and too
  loose (D).
- The accuracy figures were the interlacing route's; replaced by the Sturm route's (D).
- The generator's references: independent (Boost's tables, mpmath), not a 100-digit run of the generator (Tests).
- The kink at 1/3 is degenerate: 1/5 (prototype table).
- Subnormal widths and an intermediate overflow: the scale per term (B); the floor's denorm_min term (C).
- 4^m - 1 is not exact in float beyond m = 12 (E).
- All-zero samples give {0, 0} (C, documented).
- Romberg's accuracy test: the closed-form bound, not the estimate (Tests).
- Confirmed: the Sturm count's direction and its handling of zero terms; the weight formula's amplification; the
  truncation formula; the map's overflow cases; the level formula and Richardson; the diagonal estimate (the
  subdiagonal would be a false convergence on sqrt x); the floor after Richardson (positive weights).

From the simplicity review:
- Adopted: `std::same_as<M, Romberg>` instead of a concept named for one method (Romberg is not adaptive); four
  `static_assert` instantiations instead of 64; no test of the noisy row; finiteness once per sum; no `abs_tol`
  clamp; the stop test without the redundant max; the cap at 16; plain L1 sums; the tables named for three types.
- Not adopted: the map as a single formula for the scale (the numerics review's subnormal case needs `b - a` where it
  is finite; the nodes use the single formula); dropping `Failure`'s formatter, which the approved brief includes
  (put to the user as Q7).

From the usage sketch:
- Adopted: `.value()` in constant-evaluation examples (Q6); the infinite-end pointer in `error::non_finite`'s comment;
  the two `Criteria`'s defaults on the reference page; `Failure` formattable (Q7); the cap lowered (Q3).
- Not adopted: `b` as `std::type_identity_t<T>` (strict deduction catches integer and expression-template ends; the
  brief's sketch counted it a strength).
- Raised separately, not for this phase: fn returning `std::expected` (an integral as fn of a root finder or of an
  outer integral now needs `.value_or(nan)`, and the cause is lost). It concerns every module's `callable`
  (`math.hpp`), as roots' separate derivatives already accept it.

## Deferred

- The default method (Phase 2), the node cache (Phase 2), Gauss-Kronrod (Phase 2), `TanhSinh` (Phases 3, 4).
- N above 64; a Romberg cap above 16; compile-time tables for C++23's extended floating-point types.
- Romberg on equally spaced data (an interface over a span of samples), the use the literature leaves it.
- fn returning `std::expected`, for every module (raised with the user).
- Everything in BRIEF.md section 6.

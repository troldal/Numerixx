# Phase 4 API proposal: derivatives of order 1 to 4

Status: approved by the user on 2026-10-23 and implemented.

The user has decided:
- the derivative order is a template argument after the method, written as an integer or as a scoped enum `order`;
- the highest order is 4, as in Chapra, *Numerical Methods for Engineers*, ch. 23;
- the low-precision case gets a compiler warning, with an opt-out.

Two review rounds (numerics-reviewer, simplicity-critic, usage-sketcher) and the prototypes in `.scratch/phase4` are
folded in.

## Summary

```cpp
using namespace nxx::deriv;
auto const f = [](double x) { return std::exp(x); };

auto const d1 = derivative(f, 1.0);                           // as now: Central<2>, f'
auto const d2 = derivative<Ridders, 2>(f, 1.0);               // f'' to about 1e-12, with an error estimate
auto const d3 = derivative<Central<2>, order::third>(f, 1.0); // f''' with Chapra's 4-point O(h^2) formula
auto const e3 = derivative<Central<4>, order::third>(f, 1.0); // f''' with the more accurate 6-point O(h^4) formula
auto const d4 = derivative<Backward<2>, 4>(f, 1.0);           // Chapra's (3f(x) - 14f(x-h) + 26f(x-2h) - ...) / h^4
```

- Every method works for every order.
- `derivative<M>` is `derivative<M, 1>`, so no existing call changes. Its results are bit-identical, except where
  the Phase 2 code fails today (the overflow fix below).

## Public API (`deriv.hpp`)

```cpp
namespace nxx::deriv {

    /// The derivative orders, for readability: derivative<Ridders, order::second> is derivative<Ridders, 2>.
    enum class order : int { first = 1, second = 2, third = 3, fourth = 4 };

    /// The N-th derivative of fn at x by method M with step h: for a stencil the step, for Ridders the initial step.
    /// N is an integer or an order, 1 <= N <= 4. Precondition: h > 0 and h finite.
    template<typename M = Central<2>, auto N = 1, real T, typename F>
    requires detail::method<M> && detail::valid_order<N> && detail::callable<F, T>
    [[nodiscard]] constexpr std::expected<Estimate<T>, error>
        derivative(F &&fn, T const &x, std::type_identity_t<T> const &h);

    /// The same with the default step: eps^(1/(P+N)) * max(|x|, 1) for a stencil of order P; for Ridders
    /// max(|x|, 1) / 32 for N = 1 and max(|x|, 1) / 4 for N >= 2.
    template<typename M = Central<2>, auto N = 1, real T, typename F>
    requires detail::method<M> && detail::valid_order<N> && detail::callable<F, T>
    [[nodiscard]] constexpr std::expected<Estimate<T>, error> derivative(F &&fn, T const &x);

    // The same two overloads, additionally constrained by detail::noisy<N, T> and marked [[deprecated]]
    // (low-precision warning, below).

} // namespace nxx::deriv
```

- **`order`, the scoped enum.**
  - A scoped enum does not convert to `int` implicitly, so the order parameter is `auto N`.
  - `detail::valid_order<N>` accepts any integer type except `bool`, plus `order`, with a value from 1 to 4. Generic
    code can therefore pass an `int` or a `std::size_t` from an index sequence. 0, 5, `true` and `2.0` are rejected.
  - The order is converted to `int` once, and everything below the public overloads uses `int`.
  - `derivative<M, 2>` and `derivative<M, order::second>` are two instantiations of the public wrapper with one
    implementation.
  - The spike confirmed all of this on GCC, Clang, `cl`, clang-cl and Emscripten.
  - Being scoped, `first` ... `fourth` do not enter `nxx::deriv`, and adding `fifth` later breaks nothing.
- **The method tags, `Estimate` and `error` are unchanged.**
  - Their documentation gives the nodes for each order.
  - `Estimate::value` becomes "the estimate of the derivative".
- **`error::non_finite`:** the same cases as now, plus fn(x) not finite.
  - fn(x) is evaluated by every rule with a node at 0: the one-sided rules, the central rules for even N, and
    Ridders for N = 2 and 4.
  - So sin(x)/x at 0 works for f' and f''' with the central rules, and fails for f''.
- **The default steps depend on the method and on N.** Generic code that varies N should use the default-step
  overload rather than pass one h for every order.

### The methods for each order

**P is the accuracy order, not the number of points:** the error of `Central<P>`, `Forward<P>` and `Backward<P>` is
O(h^P), Chapra's O(h), O(h^2) and O(h^4) labels. For the first derivative a central rule happens to use P points; for
order N it uses P + N - 1, a one-sided rule P + N. The documentation of the tags says so, with this mapping for the
third derivative (each sum over the denominator times h^3):

| Call | Nodes | Formula | Chapra |
|---|---|---|---|
| `derivative<Central<2>, 3>` | +-1, +-2 | (f(x+2h) - 2f(x+h) + 2f(x-h) - f(x-2h)) / 2h^3 | centred, O(h^2) |
| `derivative<Central<4>, 3>` | +-1, +-2, +-3 | (-f(x+3h) + 8f(x+2h) - 13f(x+h) + 13f(x-h) - 8f(x-2h) + f(x-3h)) / 8h^3 | centred, O(h^4) |
| `derivative<Forward<1>, 3>` | 0 ... 3 | (f(x+3h) - 3f(x+2h) + 3f(x+h) - f(x)) / h^3 | forward, O(h) |
| `derivative<Forward<2>, 3>` | 0 ... 4 | (-3f(x+4h) + 14f(x+3h) - 24f(x+2h) + 18f(x+h) - 5f(x)) / 2h^3 | forward, O(h^2) |

The accuracy order is kept as the parameter, rather than a number of points, because it decides the error and keeps
one meaning across orders: `Central<2>` is the O(h^2) rule for every N.

| Method | Nodes | Evaluations |
|---|---|---|
| `Central<P>`, even P <= 10 | P + N - 1 symmetric nodes: -m ... m without 0 for odd N, with 0 for even N | P + N - 1 |
| `Forward<P>`, P <= 6 | 0 ... P + N - 1 | P + N |
| `Backward<P>`, P <= 6 | -(P + N - 1) ... 0 | P + N |
| `Ridders` | the `Central<2>` nodes for N, at h_start / 2^i | 2 per level for N = 1, 2 and 4 per level for N = 3, 4, plus fn(x) once for N = 2, 4; at most 32 levels |

Rounding amplification S (the sum of the absolute weights) and the model error at the default step, for double
with f and its derivatives of order 1:

| N | `Central<2>` ... `Central<10>`: S, error | `Forward<1>` ... `Forward<6>`: S, error |
|---|---|---|
| 1 | 1 ... 2.3, 4e-11 ... 1e-14 | 2 ... 28, 4e-8 ... 1e-12 |
| 2 | 4 ... 6.8, 6e-8 ... 6e-13 | 4 ... 192, 3e-5 ... 4e-10 |
| 3 | 3 ... 10.7, 2e-6 ... 1e-11 | 8 ... 956, 1e-3 ... 4e-8 |
| 4 | 16 ... 44, 1e-4 ... 3e-10 | 16 ... 4019, 1e-2 ... 2e-6 |

- The Phase 2 order caps are kept for every N: one rule is simpler than a cap per N.
- The high one-sided orders are poor for f''' and f'''', but they are what works exactly at a domain edge.

### Low-precision warning

- **What triggers it:** N = 4 in a type whose epsilon is above 1e-7, that is float, decimal32 or bfloat16, not
  double. There f'''' is mostly rounding noise:
  - `Central<2>` about 8e-2 for float;
  - Ridders' median error 2e-3, with 10% of the cases more than 100% off.
  - f''' is still usable (Ridders about 2e-4) and gets no warning.
- **How it is tested:** with integers, `radix^(digits - 1) < 10^7`. That holds in any radix, and needs no constexpr
  `epsilon()`, which some user types lack.
- **The mechanism: a deprecated overload, not a deprecated helper.**
  - The two public overloads have twins constrained additionally by `detail::noisy<N, T>` and marked
    `[[deprecated("...")]]`. Subsumption picks the twin for noisy calls.
  - The spike checked this on all five toolchains, both with ordinary include paths and as a system header
    (`-isystem`, `/external:I` with `/external:W0`, `/imsvc`). In every case:
    - exactly one warning, reported at the user's call;
    - none for double, and none for float f''';
    - a pragma around the call silences it.
  - A deprecated helper called inside the body was the first idea. It is reported inside the header instead, so a
    pragma at the call cannot silence it, and GCC and `cl` suppress it entirely when the library is a system header,
    which is how CMake gives it to `find_package` users.
- **Opt-out (user's choice: a macro).** Defining `NUMERIXX_NO_PRECISION_WARNING` turns `detail::noisy` off.
  - Since the warning can now also be silenced per call, the macro is an extra rather than the only way out.
  - It must be defined the same way in every translation unit, best on the command line: overload resolution would
    otherwise differ between them.
  - Our float f'''' tests silence the warning per call, so the tests do not need the macro.
- **Tests:**
  - a CTest that compiles a float f'''' with warnings as errors and is expected to fail (`WILL_FAIL`);
  - the same with the per-call pragma, and with the macro, expected to pass;
  - a `static_assert` that the warning's condition is false for double.

## Algorithms

### Tables: one generator

The weights of every rule come from the Lagrange formula (Burden et al., *Numerical Analysis*, numerical
differentiation):

```
w_k = N! [t^N] prod_{j != k} (t - n_j) / (n_k - n_j)
```

This is the N-th derivative at 0 of the polynomial that interpolates fn at the nodes, the same rules Chapra derives
from Taylor series.

- **One constexpr function over (nodes, N).** Each family only supplies its nodes, which removes the `Backward`
  mirroring: for even N the backward weights are the forward ones mirrored, not negated.
- **Exact rational arithmetic** (long long, reduced with gcd), over one common denominator, as now.
  - The largest intermediate is 12! = 4.8e8 (`Central<10>`, N = 4), times N! = 24, far inside long long.
  - The largest numerator is 6,222,216 and the largest denominator 453,600, both below 2^24, so the weights are
    exact in float.
  - A `static_assert` in the generator checks every value against 2^24 before it is narrowed to `int`.
- **It replaces the three Phase 2 closed forms** (`make_table`). The `Table` size parameter is renamed, since N now
  means the order. `docs/design.md` ("Fornberg's algorithm is not needed for three fixed families") is updated.
- **Phase 2 tables:** a constexpr spike reproduces every one of them exactly. The implementation step checks all of
  them once and records the result in the commit, so no second implementation is kept in the tests.
- **Test references** for every table:
  - the moments, M_j = 0 for j != N up to N + P - 1, M_N = N! denominator, and M_{N+P} != 0, so the order is exactly
    P. They are a Taylor check, independent of the interpolation formula.
  - a reduced denominator (the gcd of the numerators and the denominator is 1), which with the moments fixes a table
    uniquely;
  - Fornberg (1988), as now, for derivatives 0 to 4;
  - the user's Chapra row, `Backward<2>`, N = 4.
  - A test that the backward weights are (-1)^N times the mirrored forward weights.

### Evaluation

- **The general rule:** `sum_k numerators[k] * fn(x + nodes[k] h)`, then `sum / (denominator * h)`, then N - 1 more
  divisions by h.
  - h^N is never formed: h^4 is below the smallest normal double for h < 1.2e-77.
  - The repeated divisions cannot overflow or underflow spuriously. If h >= 1, each intermediate is at most the one
    before; otherwise the final value is the largest.
- **No overflow before the final quotient.** Before forming each term, `|f_k| <= max / |numerators[k]|` is checked,
  then `sum_fits` on the partial sum. A failing check is `non_finite`.
  - This fixes a Phase 2 bug. `derivative<Central<4>>` of the constant 1e308 at 0 with h = 1 forms -8e308, then
    inf - inf: GCC rejects the overflow in a constant expression, and Clang the NaN. The run-time result is
    `non_finite`; the true answer is 0. With the larger numerators of N >= 2 (up to 6.2e6) it is easy to reach.
  - Results that are finite today are unchanged, bit for bit. The cost is two comparisons per node.
  - The only remaining compile-time exception is an overflow of the final quotient, which GCC rejects. It is
    documented as now, and becomes easier to reach for N = 4 with a tiny user step.
- **Actual offsets for the two-node rules (N = 1, as now) and for `Central<2>` with N = 2.** Every other rule divides
  by the nominal h, as the multi-node rules do now.
  - N = 2 uses `((f+ - f0) / h+ - (f0 - f-) / h-) / ((h+ + h-) / 2)` with h+ = (x + h) - x and h- = x - (x - h).
    - Before each subtraction, `sum_fits(f+, -f0)`, `sum_fits(f0, -f-)` and the finiteness of d+ and d- are checked.
    - The form is exact for quadratics whenever the offsets are exact, which by Sterbenz's lemma holds for
      0 < h <= |x| in any radix.
  - The evidence, from the Ridders prototype, nominal h against actual offsets:
    - N = 1: coverage 98.0% against 99.3%, and the worst case 22 times the estimate against 3.6.
    - N = 2: coverage equal, but the worst case 8.7 against 1.4 times for float, and 2.5 against 1.5 for double.
    - N = 3 and 4 do well with the nominal h, also near zeros of f, where node rounding would show first: sin near
      pi and cos(5x) near its zeros are 99.4 to 99.7% covered, worst 1.7 times.
- **Overflow and resolution checks:** `fits(x, h, q)` with q = max(2K, denominator), and x +- h != x, as now. Ridders
  uses q = 2 for N = 1, 2 and q = 4 for N = 3, 4, whose bases reach x +- 2h.
- **Default step** `eps^(1/(P+N)) * max(|x|, 1)`, the f' rule with P + N in place of P + 1.
  - It comes from the error model C h^P |f^(P+N)| + S eps |f| / h^N, with the constant (N S / (P C))^(1/(P+N))
    dropped.
  - That costs more as N grows. At the model optimum, the default step's error is 1.1 times the optimum for
    `Central<2>` with N = 1, 3.5 times with N = 2, 2.5 with N = 3, 11 with N = 4, and 15 for `Central<4>` with N = 4.
    `docs/design.md` ("at most 1.64") is updated with this.
  - The constant is still left out. On nine test functions, half the default step was better for every order, N = 1
    included: by 1.2 to 1.8 times for `Central<2>`, up to 200 times for `Central<8>`. These functions vary faster
    than the assumed scale max(|x|, 1), and the model's larger step would make them worse.
  - So it is one convention for all N, and the documentation states it with the sensitivity: for N = 4, misjudging
    fn's scale by a factor 2 costs up to 16.

### Ridders

The Phase 3 algorithm and constants (`PHASE3_API.md`, steps 2 to 8: ratio 2, a stall of two levels, the
untrusted-result factors 4 and 64, 32 levels). Three things are new.
- **The base difference is the `Central<2>` rule for N.** Its error expansion is even for every N, so the factors
  4^j - 1 apply.
  - N = 1 and 2: the actual-offset forms above.
  - N = 3: `(f(x+2h) - 2f(x+h) + 2f(x-h) - f(x-2h)) / 2h^3`, which is f''' + h^2 f^(5)/4 + ...
  - N = 4: `(f(x+2h) - 4f(x+h) + 6f(x) - 4f(x-h) + f(x-2h)) / h^4`, which is f'''' + h^2 f^(6)/6 + ...
- **The rounding floor** is `eps sum_k |w_k| |f_k| / denominator`, with eps applied to |f_k| before the weight, then
  divided by h N times.
  - For N = 2 it is written from the actual offsets.
  - Every division saturates at max, and every sum through `sum_fits`. One saturation is not enough: at x = 0,
    fn = 1e300 with h = 1e-300 overflows the first division, and fn = 1 overflows only the last with h = 1e-170 for
    N = 2, 1e-130 for N = 3 and 1e-100 for N = 4.
- **fn(x)** (N = 2, 4) is evaluated once per call and kept across restarts. If it is not finite, the result is
  `non_finite` after one evaluation.
- **The initial step** is max(|x|, 1) / 32 for N = 1 (as now) and max(|x|, 1) / 4 for N = 2 to 4.
  - The base's rounding error grows as 1/h^N, so a larger start pays.
  - N = 2: s/4 against s/8 is 1.6 times as accurate in double (median 3.9e-13 against 6.2e-13), for 1 more
    evaluation; float is also better. Against s/32 it is 7 times.
  - N = 3 and 4: s/4 against s/8 is 1.6 and 2.6 times as accurate, for 3 and 2 more evaluations.
  - The steps are a constexpr array indexed by N.
- **One driver,** `detail::ridders<N>`, with the base chosen by `if constexpr` inside, each base returning its value
  and its floor. The implementation step compares N = 1 against the current code on random cases once, as a
  regression check, and records the result.
- **Cost:** at most 64, 65, 128 and 129 evaluations for N = 1 to 4.

## Prototype results (Ridders)

Python, numpy `float64` and `float32`, and `decimal` for 50 and 100 digits, from the default initial steps. Random
cases: nine functions, x in [0.2, 3], 300 each, with closed-form derivatives of every order.

| N, initial step | Covered (double) | Worst actual / estimate | Median relative error | Evaluations, mean / max |
|---|---|---|---|---|
| 1, s/32 (Phase 3) | 99.3% | 3.6 | 9.5e-15 | 14 / 22 |
| 2, s/4 | 99.7% | 1.5 | 3.9e-13 | 17 / 23 |
| 3, s/4 | 99.6% | 2.0 | 3.1e-11 | 32 / 40 |
| 4, s/4 | 100% | 0.94 | 7.4e-10 | 30 / 41 |

- **float:**
  - N = 2: 99.4% covered, worst 1.4 times the estimate.
  - N = 3: 99.5% covered, worst 5.3 times, median error 1.6e-4.
  - N = 4: 99.6% covered, worst 1.2 times, median 2.1e-3. The estimate reports the noise: for exp at 1 it is 29 times
    the value.
- **Poles:** all 400 poles 1/(x - a), a in (1/128, 1/32), are covered for every N.
- **Near zeros of f** (N = 3, 4): 99.4 to 99.7% covered, at most 44 evaluations. sin at pi for N = 4 gives -1.9e-16
  against 1.2e-16, with an estimate of 4.5e-13, after 17 evaluations.

| Case (double) | Relative error | Relative estimate | Evaluations |
|---|---|---|---|
| N = 2: exp at 1 (from s/8) | 4.0e-13 | 3.6e-12 | 13 |
| N = 3: exp at 1 | 3.4e-11 | 7.2e-10 | 28 |
| N = 3: x^5 at 2 | 0 | 4.9e-14 | 12 |
| N = 3: log at 1e-4 (edge, restarts) | 3.9e-10 | 2.0e-8 | 45 |
| N = 3: 1/(x - 0.01) at 0 (untrusted restart) | 2.2e-11 | 7.5e-9 | 108 |
| N = 4: exp at 1 | 5.4e-10 | 2.1e-7 | 25 |
| N = 4: log at 1e-4 (edge, restarts) | 5.1e-8 | 4.6e-6 | 46 |
| N = 4: 1/(x - 0.01) at 0 (untrusted restart) | 2.3e-8 | 3.3e-6 | 109 |
| N = 3: exp at 1, 50 / 100 digits | 3.0e-41 / 2.8e-87 | 6.4e-39 / 5.6e-85 | 48 / 64 |
| N = 4: exp at 1, 50 / 100 digits | 1.3e-37 / 1.5e-81 | 5.6e-35 / 1.6e-78 | 45 / 65 |

For N = 2 the Phase 3 constants were rechecked:
- A one-level stall lets float's worst case reach 3.3 times the estimate.
- The untrusted factor 1 instead of 4 leaves 32 of the 400 poles uncovered.
- Factors of 16 to 256 over the floor give identical results.

For N = 3 and 4 only the coverage was measured.

## Behaviour to document

- **Accuracy** is about eps^(P/(P+N)) for a rule of order P.
  - Ridders' median is about 1e-14, 4e-13, 3e-11 and 7e-10 relative for N = 1 to 4 in double.
  - Use the stencils inside iterations and Ridders for a reported value.
  - For float, f''' is good to about 1e-4 and f'''' is mostly noise (the warning).
- **The default steps assume fn varies on the scale max(|x|, 1).** The sensitivity grows with N (above).
  - At large |x| a periodic fn aliases: sin at fl(10000 pi) from the default start. Ridders for f' has the same
    limit.
  - Pass a step on fn's own scale.
- **Exactly periodic samples fool every sampling method.** If the samples of the first two levels lie on a
  polynomial, Ridders stops on exact agreement: cos(64 pi x) at 0 for N = 2, sin(48 pi x) for N = 3, cos(48 pi x) for
  N = 4.
- **Domain edges.**
  - Ridders needs the nearest bad point at least about 2^-31 max(|x|, 1) away for N >= 2, and 2^-34 for N = 1.
  - Exactly at an edge, use `Forward<P>` or `Backward<P>`.
- **Cost:** the method table. Ridders for N = 3 and 4 typically takes 25 to 45 evaluations, and up to 110 after
  restarts.
- **The final quotient** can overflow for N = 4 with a tiny user step (any nonzero sum over h^4 for h below about
  1e-81). That is `non_finite` at run time, and GCC rejects it at compile time.

## Tests

For every N, as in Phase 3, plus:
- **Tables:** moments, a reduced denominator, Fornberg's rows, the Chapra row, and the (-1)^N mirror. All values are
  below 2^24.
- **Polynomials** (`static_assert` and run time, double):
  - Every rule is exact for degree up to N + P - 1, with dyadic x and a power-of-two h. For each case claimed exact,
    a `static_assert` checks the bit budget `sum_k |numerators[k]| |a + n_k b|^d < 2^53` (x = a 2^-s, h = b 2^-s).
    Float only for the low orders.
  - Ridders gives x^(N+1) exactly at level 1. x^(N+2) needs column 1: at x = 2 with h = 1/8 it is exact for
    N = 2, 3, 4 (7, 12 and 13 evaluations).
  - x^(N+4) needs column 2, which divides by 15. It is tested against a tolerance derived from the rounding of the /3
    and /15 steps on the h^4 term (x^8 for N = 4 is 2.6e-11 off).
- **The actual offsets:** fn(t) = t at x = 1 with h = 0x1.4p-53 gives f'' = 0 exactly, where a nominal-h formula
  gives 5.8e15.
- **Never forming h^N:** exact values with steps outside the range of h^N, per N.
  - N = 2: 2^-600 and 2^600.
  - N = 3: fn = (2^300 t)^3 with h = 2^-600, and fn = (t / 2^260)^3 with h = 2^600.
  - N = 4: fn = (2^100 t)^4 with h = 2^-300, giving f'''' = 24 * 2^400, and its mirror.
- **Overflow:**
  - the Phase 2 bug above, the constant 1e308 for `Central<4>`, now `non_finite` in a constant expression too;
  - a constant 1e308 for N = 4;
  - f0 = -max with f+ = max for N = 2;
  - the floor saturation cases above, per N, at x = 0;
  - fn(x) not finite, at a cost of one evaluation for Ridders N = 2 and 4.
- **Accuracy and coverage** against closed forms, double and float, with tolerances from the model or from the
  prototype as in Phase 3.
  - float N = 3 gets its own bound of 10 times the estimate, since the prototype's worst case is 5.3.
  - Multiprecision for N = 3 and 4: `cpp_bin_float_50` with and without expression templates, and
    `cpp_dec_float_50`, with saturation inputs chosen per type.
- **The constraint:**
  - `derivative<M, 0>`, `derivative<M, 5>`, `derivative<M, true>` and `derivative<M, 2.0>` are rejected;
  - `derivative<M, order::second>` equals `derivative<M, 2>`;
  - `order` does not convert to `int` implicitly.
- **The warning:** the compile tests above.
- **N = 1:** the existing tests, unchanged.

## Decisions from the reviews

- **First round, kept:**
  - the overflow guards of the N = 2 form;
  - a floor saturated at every step;
  - the statement of the unequal-offsets term (0 in binary unless x +- h straddle a binade);
  - the aliasing at large |x|, which is documented;
  - x^4 and x^6, the actual-offset case, and the documentation of the default steps, fn(x) and which method suits
    which job.
- **Second round, adopted:**
  - the per-term overflow check (the Phase 2 bug);
  - q = 4 for Ridders N = 3, 4;
  - the cost of the dropped step constant, now stated;
  - the h^N test cases per N;
  - the floor-saturation cases per N;
  - the integer warning criterion;
  - the deprecated overload instead of the deprecated helper;
  - the bit-budget check;
  - the float N = 3 bound;
  - the 2^24 check;
  - one generator over nodes;
  - one-off checks instead of a kept copy of the old code;
  - one Ridders step for N >= 2;
  - a constexpr step table;
  - the `Forward<1>` and `Backward<1>` wording, whose N = 1 path stays as now.
- **Not adopted:**
  - **Nominal h for N = 2** (simplicity-critic): float's worst case is 8.7 against 1.4 times the estimate.
  - **Divided differences for N = 3, 4** (numerics-reviewer, conditional on node-rounding failures): the runs near
    zeros of f showed none.
  - **A default-step constant per rule:** the test functions favour smaller steps, not larger ones.
  - **Dropping the warning** (simplicity-critic): the user wants it, and the per-call silencing removes the main
    objection.
  - `second_derivative`, the reduced method set, a `default_step<M, N>` accessor, an "estimate below |value|"
    helper.

## Changes during implementation

From the implementation review (numerics-reviewer), each reproduced first:
- **The term check is `|f_k| >= max / |c|`, not `>`.** fl(max / 6) rounds up, and 6 times it ties to infinity, so a
  constant max / 6 under `Central<2>` for f'''' still overflowed in a constant expression. The same applies to the new
  `saturating_product`.
- **The one-sided quotients of the second difference** are pre-checked (`|a| / max < b / 2`, conservative), since an
  overflow there is not the final quotient. A step function of +-1e300 with h = 1e-10 is now `non_finite` in a constant
  expression too.
- **Ridders' floor for N = 3 and 4 includes the rounding of the points**, `eps/2 (|x| + 2h) |f'| sum |w_k| /
  (den h^N)`, with f' from the values at x +- h, which the nominal step leaves in. Without it, f''' of exp(a x),
  a = 10 to 50, at x in [1, 3] was covered in 96.3 to 96.7% of the cases; with it in 99.0 to 99.3%, and f'''' in
  100%. The random study stays at 100%, and float's worst case for f''' drops from 5.3 to 2.0 times the estimate, so
  float f''' needs no bound of its own.
- **Neighbouring points of a multi-node stencil must differ**: at h near 0.6 ulp(x), x + h and x + 2h round to the same
  value, and `Central<2>` for f'''' returned a finite wrong value. Ridders had the check already.
- **A stencil whose integer weights T cannot hold exactly** (radix^digits not above the largest table value, half
  precision under `Central<10>`) is a `static_assert`, not a silently inexact rule (or, beyond the range, undefined
  behaviour). (Changed after Phase 4: a constraint of `derivative`, `detail::exact_weights`, so that a
  requires-expression sees it; Ridders is checked against the `Central<2>` table it uses.)
- **An integer order is compared in its own type**, so `0x1'0000'0001ULL` is rejected rather than wrapped into 1.
From an independent review by another model (OpenAI Codex), each reproduced first:
- **The term check is exact** (`product_fits`): `|f_k| >= max / |c|` also rejected |f| = max under a weight of 1, where
  the product cannot overflow; Phase 3 returned max / 12 for a `Central<4>` case that the check turned into
  `non_finite`. Now a value below fl(max / |c|) fits, above it does not, and the tie is decided by the scaled product,
  exact like `sum_fits`. The comparison of 180,000 first-derivative results with the Phase 3 header, with values up to
  max, is identical except for Ridders on subnormal values (next point).
- **The rounding of a value is at least denorm_min** for a nonzero subnormal value, where rounding is absolute: with
  `eps |f|` alone, Ridders on a fn with subnormal values reported a result 192 denorm_min off with an estimate of 0.
  For normal values nothing changes.
- **The divisions are ordered to avoid a spurious underflow**: with `denominator * h > 1 > h`, dividing by
  `denominator * h` first shrank a subnormal sum to 0 (`Central<10>` for f'''' returned 0 for a normal 2.6e-308); the
  divisions by h now come first there, and only a result within that factor of max can overflow spuriously, which is
  `non_finite` rather than a silent 0.
- **Breaking, for an unusual call:** T given explicitly as the second template argument, `derivative<Ridders,
  double>(...)`, no longer compiles, since that argument is now the order. T is deduced from x as before; no call in
  the repository or the documentation spelled it.

## Deferred

- **Several orders from shared evaluations,** for Newton's method for a minimum, Halley's method or curvature:
  `derivatives<M, N>` returning all orders up to N. Nothing here blocks it.
- **Orders above 4:** the constraint and the enum widen without a break.
- **A default-step constant**, or a scale parameter for fn, for all N at once.
- **Mixed partials and Hessians**, with gradients and Jacobians (`plan.md`).

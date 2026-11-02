# Phase 3a API proposal: the method protocol

Status: draft 1, for review (numerics-reviewer, simplicity-critic, usage-sketcher), then the user. The user's
decisions on `docs/open-closed/ANALYSIS.md` (2026-11-02): the methods open to users (D1), quadrature first (D4).
Architecture: `docs/quadrature/BRIEF.md` (approved 2026-10-30), whose D1 ("methods as a closed set of empty tags")
this phase amends. Phase 3 (`PHASE3_API.md`, TanhSinh, approved 2026-11-01, not yet implemented) becomes Phase 3b and
adds `TanhSinh` through this protocol. The evidence is a spike on a patched copy of `quadrature.hpp` (Spike, below).

## Summary

```cpp
#include <numerixx/quadrature.hpp>
namespace qd = nxx::quadrature;

// Unchanged: every call that compiles today compiles and returns the same bits.
auto const r = qd::integrate(fn, 0.0, 1.0);
auto const s = qd::integrate<qd::Romberg>(fn, 0.0, 1.0, {.rel_tol = 1e-12});
auto const g = qd::integrate<qd::GaussLegendre<20>>(fn, 0.0, 1.0);

// New: a method of one's own is a type with one const member, used like the library's.
struct Midpoint
{
    int panels = 100;

    template<nxx::real T, typename F>
    constexpr std::expected<T, qd::error> apply(F &fn, T const &a, T const &b) const; // a fixed rule
};
auto const m = qd::integrate<Midpoint>(fn, 0.0, 1.0);
auto const n = qd::integrate(fn, 0.0, 1.0, Midpoint{.panels = 1000}); // a method with parameters, by value

// A method may compose others through the public integrate(): a substitution over any estimating method.
auto const t = qd::integrate(fn, 0.0, 1.0, {.rel_tol = 1e-12}, SquareAtLowerEnd<qd::GaussKronrod<15>>{});
```

- Two public concepts, `fixed_rule` and `estimating_method`, each a single const member: `apply(fn, a, b)` and
  `estimate(fn, a, b, criteria)`.
- `GaussLegendre<N>`, `Romberg` and `GaussKronrod<15>` implement them, by forwarding to the `detail` functions they
  are dispatched to today; `integrate` calls the member instead of naming the methods. Results are bit for bit
  unchanged.
- `integrate` takes the method as an optional last argument, so that a method with parameters can be passed.
- The closed lists `detail::gauss_points`, `detail::fixed_rule` and `detail::estimating` are removed, and so is the
  `if constexpr` dispatch on `Romberg`. Adding a method, the library's or a user's, adds a type and edits nothing.
- No new types. Everything constexpr; nothing throws; no ADL (the protocol is a member call).

## Public API (changes to `quadrature.hpp`)

```cpp
namespace nxx::quadrature {

    /// A fixed rule: a method M whose const member apply(fn, a, b) returns the integral from a to b, or
    /// error::non_finite. integrate() has checked the ends first (contract: B).
    template<typename M, typename F, typename T>
    concept fixed_rule = std::copy_constructible<M> && requires(M const &m, F &fn, T const &a, T const &b) {
        { m.apply(fn, a, b) } -> std::same_as<std::expected<T, error>>;
    };

    /// A method with an error estimate: a method M whose const member estimate(fn, a, b, criteria) returns the
    /// integral from a to b with an estimate of its error, or a Failure. integrate() has checked the ends and the
    /// criteria first (contract: B).
    template<typename M, typename F, typename T>
    concept estimating_method = std::copy_constructible<M>
        && requires(M const &m, F &fn, T const &a, T const &b, Criteria<T> const &criteria) {
               { m.estimate(fn, a, b, criteria) } -> std::same_as<std::expected<Integral<T>, Failure<T>>>;
           };

    template<int N>
    struct GaussLegendre
    {
        template<real T, typename F>
        requires(N >= 1 && N <= 64)
        constexpr std::expected<T, error> apply(F &fn, T const &a, T const &b) const;
    };

    struct Romberg
    {
        template<real T, typename F>
        constexpr std::expected<Integral<T>, Failure<T>>
            estimate(F &fn, T const &a, T const &b, Criteria<T> const &criteria) const;
    };

    template<int N>
    requires(N == 15)
    struct GaussKronrod
    {
        template<real T, typename F>
        constexpr std::expected<Integral<T>, Failure<T>>
            estimate(F &fn, T const &a, T const &b, Criteria<T> const &criteria) const;
    };

    template<typename M, real T, typename F>
    requires nxx::detail::callable<F, T> && fixed_rule<M, std::remove_reference_t<F>, T>
    [[nodiscard]] constexpr std::expected<T, error> integrate(F &&fn, T const &a, T const &b, M const &m = M{});

    template<typename M = GaussKronrod<15>, real T, typename F>
    requires nxx::detail::callable<F, T> && estimating_method<M, std::remove_reference_t<F>, T>
    [[nodiscard]] constexpr std::expected<Integral<T>, Failure<T>> integrate(F &&fn, T const &a, T const &b,
        Criteria<std::type_identity_t<T>> const &criteria = {}, M const &m = M{});
}
```

- **The members** are declared in the structs, so that each method's Doxygen comment shows its protocol, and defined
  after `detail`, just before `integrate`. A qualified name such as `detail::interval_of` in a member defined inside
  the class would be looked up where the class is defined, before `detail` exists.
- **`GaussLegendre<N>` stays an unconstrained struct.** The range 1 to 64 moves from `detail::gauss_points` to a
  requires-clause on its member, so `GaussLegendre<0>` and `GaussLegendre<65>` are still types that `integrate`
  rejects. A constrained struct would make naming them a hard error, which breaks the existing `!integrable<...>`
  tests. `GaussKronrod<N>` keeps its `requires(N == 15)` on the type (Phase 2, Q11).
- **One function per kind serves both call forms** (Q3). `M` is named (`integrate<Romberg>(fn, a, b)`, with `m`
  value-initialised) or deduced from a method value (`integrate(fn, a, b, {}, Romberg{})`). A fixed rule has no
  criteria, so its method value comes fourth: `integrate(fn, a, b, Midpoint{.panels = 1000})`. With no method named
  or passed, `M` falls back to `GaussKronrod<15>` as today. The fixed overload has no default for `M`, so `{}` in
  fourth place cannot select it unless `M` is named.
- **The constraints check fn first** (`nxx::detail::callable`, unchanged), so a bad fn is reported as today, then the
  method. Each concept takes F and T: a method may accept fewer of them (a method that needs exp in T, a method for
  one T). F is the reference-stripped type that `integrate` passes as an lvalue.
- **The documentation** (Doxygen) of `integrate`, `Integral` and `Criteria` changes as in B: what was said of every
  method becomes what is said of the library's methods, and the `@tparam M` lists read "a `fixed_rule`
  (`GaussLegendre<N>`, 1 <= N <= 64, or a user's)" and "an `estimating_method` (`GaussKronrod<15>`, the default,
  `Romberg`, or a user's)".

## Behaviour

### A. What `integrate` does

Unchanged, in the same order (design.md, "Validation, in order"); then the method's member:
- **Fixed rule:** the ends (NaN or infinite: `non_finite`); a == b: 0 without calling fn; then `m.apply(fn, a, b)`.
- **Estimating:** the tolerances (negative or not finite: `non_finite`), `rel_tol` clamped to at most 1; the ends;
  a == b: `{0, 0}` without calling fn; then `m.estimate(fn, a, b, Criteria<T>{abs_tol, rel_tol})`, the criteria
  validated and clamped. The validated pair already went to `detail::romberg` and `detail::gauss_kronrod` as two
  arguments. It now goes as the type the caller wrote, so the protocol adds no type (Q6).
- **The result** is returned as the member gave it, after a debug check of the contract's checkable part (B, Q5).

### B. The contract

What every method may assume, and what it must give. The library's three methods meet it, and so does TanhSinh's
approved specification. It goes into the concepts' Doxygen comments, into `design.md`, and into the reference page.

- **Given to the method:** a and b finite, a != b, in the caller's order (b < a asks for the negated integral, Q4).
  For an estimating method, criteria with both members finite and non-negative and `rel_tol <= 1`; `abs_tol` is not
  clamped (it only enters a maximum). fn as `integrate` received it, as a non-const lvalue: call it with
  `std::invoke(fn, x)` on a `T const` lvalue and convert the result with `T(...)` (the module's rule on fn's result,
  `docs/design.md`).
- **Required of a fixed rule:** a finite T, or `error::non_finite`.
- **Required of an estimating method:** an `Integral` whose value and `abs_error` are finite, `abs_error` not negative,
  when its own stop test was met; otherwise a `Failure`: `non_finite` without a best, or `not_converged` with a best
  that is finite in the same way, or without one.
- **Required in meaning, not checkable:** `abs_error` estimates |value - exact| (an estimate, not a bound, as
  `Integral` says); success means the method's own stop test against `max(abs_tol, rel_tol |value|)` (or its own
  floor) was met. fn's exceptions pass through, the method's own too: the library catches nothing.
- **Checked:** `integrate` asserts the checkable part (finite value; finite, non-negative `abs_error`, in the result
  and in a best; no best with `non_finite`), as `roots.hpp` asserts its unreachable states. Under `NDEBUG` nothing is
  checked; in a constant expression a failed assert is a compile error.
- **No longer universal:** `Integral` says "abs_error is never below the rounding floor", and `Criteria` says "the
  rounding floor ... is the effective minimum of every request". Both are properties of the library's methods. The
  rounding floor needs the method's own L1, so `integrate` cannot enforce it for another method. The wording becomes
  "for the library's methods".

### C. The library's methods

- `GaussLegendre<N>::apply` is the body of today's fixed `integrate` after its checks: `detail::interval_of`, then
  `detail::with_gauss_legendre<N, T>` and `detail::apply`.
- `Romberg::estimate` is `detail::romberg(fn, a, b, criteria.abs_tol, criteria.rel_tol)`.
- `GaussKronrod<N>::estimate` is `detail::with_gauss_kronrod<T>` around `detail::gauss_kronrod`.
- Nothing inside `detail` changes. The members add one inlined call per integration.

### D. Phase 3b (TanhSinh) through the protocol

`PHASE3_API.md` stays approved as it is, except for its "Public API" bullet, which becomes:
- **`TanhSinh` implements `estimate`**, so the existing estimating overload of `integrate` takes it. No list, no
  dispatch and no overload changes. The approved text said that `detail::estimating<M>` admits `TanhSinh`.
- "`TanhSinh` stays a plain struct: a later option goes into `Criteria` or a new tag" gains a third way: a later
  option can be a defaulted data member (`TanhSinh{.levels = 9}`), which leaves `qd::TanhSinh` and `TanhSinh{}`
  valid.

Its behaviour (A to F), the decisions Q1 to Q10 and the tests are unchanged. `docs/plan.md` splits Phase 3 into 3a
(this protocol) and 3b (TanhSinh, `PHASE3_API.md`).

### E. What a method of one's own cannot use yet

The building blocks that make the library's methods robust stay in `detail`: `Interval` (the overflow-free map and
the order of each term's products), `Sum` (Neumaier), `rounding_floor`, the overflow checks for constant evaluation,
and the node tables. A user's method works without them, as the spike's do, but it does not inherit their guarantees.
Publishing them is deferred until a method needs one (ANALYSIS.md, 6.4).

## Spike

Kept in `docs/open-closed/quadrature-spike/`, beside the roots spike:
- `quadrature.patch.py` writes this proposal's changes into a copy of `quadrature.hpp` (concepts, members, entry
  points; the closed lists removed).
- `dump.cpp` prints every result of a battery bit for bit: 8 functions (e^x, sqrt x, 1/sqrt x, log x, a kink, a peak
  of width 1e-3, sin 50x, NaN on half the interval) on 6 intervals, both orientations, by `GaussLegendre<1>`, `<7>`,
  `<20>` and `<64>`, `GaussKronrod<15>` at the default, 1e-12 and 0, and Romberg at the default, 1e-10 and an
  oversized tolerance, in float, double, long double, `cpp_bin_float_50` and `cpp_dec_float_50`. It also covers edges
  (the constants 1 and 1/4 on [-max, max], a == b, a NaN end, a negative tolerance) and three results computed in
  constant expressions. Built once against today's header and once against the patched one: the 2,438 lines are
  identical (`cmp`). They hold 919 successes, 553 failures with or without a best and 180 fixed-rule errors, the
  rest fixed-rule values; hexadecimal for the built-in types, max_digits10 digits for the Boost types.
- `user_methods.cpp`: methods outside `nxx`, on the public protocol alone. A fixed rule with a run-time parameter
  (`Midpoint{.panels}`) and a combinator over any estimating method (`SquareAtLowerEnd<Inner>`, x = a + (b - a) u^2),
  in constant expressions and at run time. Also the concepts' rejections (a wrong return type, a non-const member,
  the two kinds kept apart, `GaussLegendre<0>`, `<65>`, `int`). Built with `-Wall -Wextra -Wpedantic -Wshadow
  -Wconversion`, no warnings. On 1/sqrt x over [0, 1] at rel_tol 1e-12: `GaussKronrod<15>` 2,370 evaluations,
  error 3.0e-13; through the substitution 150 evaluations, exact to the last bit. `SquareAtLowerEnd<Romberg>` is
  `non_finite`, as Romberg is documented to be: it samples the ends, where the substituted integrand is 0 times inf.
- The existing quadrature tests (`test_quadrature_gauss.cpp`, `_romberg.cpp`, `_kronrod.cpp`) built unchanged against
  the patched header pass exactly as against today's. That means 25 of 26 test cases in both builds: the same one
  fails in both, the `GaussLegendre<2>` long double node baseline, which was recorded on the project's toolchains
  and does not hold for GCC 13's 80-bit long double.
- **Orientation** (Q4): with today's header, `integrate(fn, b, a) == -integrate(fn, a, b)` held bit for bit in 3,000
  cases for each of `GaussLegendre<64>`, Romberg and `GaussKronrod<15>`. The compensated sums hide the reversed
  summation order, but nothing proves they always will.
- Built with GCC 13.3 and Boost 1.83 (the container's; the project pins 1.92, whose download the container's network
  refuses), not on the project's five toolchains. The phase itself builds and tests on all five, as every phase does.

## Tests (`tests/test_quadrature_protocol.cpp`, test cases `quadrature:protocol:*`)

- **The existing suite unchanged**, and passing bit for bit (its exact counts and values).
- **Both call forms agree:** `integrate<M>(fn, a, b[, c])` equals `integrate(fn, a, b, [c,] M{})` for the three
  library methods, in double and `cpp_bin_float_50`; `static_assert` for double.
- **A user fixed rule with a parameter** (the composite midpoint rule) against closed forms. On x^2 over [0, 3] with
  3 panels the rule's error, (b - a) h^2 f'' / 24 = 1/4, is exact for a quadratic, so the value is exactly 8.75 in
  double, and -8.75 reversed. `static_assert` with run-time twins.
- **A user combinator over a library method** (`SquareAtLowerEnd<GaussKronrod<15>>`) on 1/sqrt x over [0, 1] = 2 and
  log x over [0, 1] = -1, within `abs_error` and the tolerance; in double and `cpp_bin_float_50`.
- **The library's validation reaches a user method:** a NaN end, a negative tolerance and a == b return what they
  return today, with a fn that fails the test if it is called.
- **The concepts:** a wrong return type, a non-const member and `int` are rejected; the two kinds stay apart
  (`Romberg` is not a `fixed_rule`, `GaussLegendre<7>` not an `estimating_method`); criteria for a fixed rule still
  do not compile (the existing `with_criteria`).
- Not tested: the assert of B (doctest has no death tests, and Emscripten aborts the whole run).

## Documentation

- `quadrature.hpp`: the file comment (a method is a type, and a user's works the same way); the concepts' comments
  with the contract of B; each member's comment; `integrate`'s `@tparam M`, `@param m` and an example with a method
  value; `Integral` and `Criteria` as in B.
- `docs/quadrature/design.md`: "The interface" rewritten for the protocol and the contract; BRIEF.md D1 marked as
  amended (the set is open).
- `docs/sphinx/quadrature/reference.md`: the two concepts. The guide page on writing a method, with
  `user_methods.cpp` turned into an example built as a test, goes into the milestone's chapter on the site.
- `docs/plan.md`: Phase 3 split into 3a and 3b; the "Customization" item on an open method protocol points here for
  quadrature.

## Decisions for the user

- **Q1. Member names `apply` (fixed rule) and `estimate` (estimating method).** Recommended: yes. The name says which
  kind a method is, and a method cannot satisfy both kinds by accident. Alternative: `integrate` for both, told apart
  by the arguments and the return type. That reads naturally but invites calling the member directly, which skips the
  validation.
- **Q2. Concept names `fixed_rule` and `estimating_method`, with the parameters `<M, F, T>`.** Recommended: yes. They
  follow the standard library's naming and reuse the names already in `detail`. Alternative: `<M, T>` alone, checked
  with a stand-in fn; that is simpler, but it cannot express a method that accepts only some fn.
- **Q3. The method value as the last argument, defaulted, so one function per kind.** Recommended: yes. It adds no
  overload, keeps today's argument order and keeps `integrate<M>` the usual form. One new acceptance:
  `integrate<GaussLegendre<10>>(fn, a, b, {})` now compiles (the `{}` makes the method), which is harmless. The
  designated `{.rel_tol = ...}` is still rejected. Alternatives: separate overloads that take the method value (two
  more overloads), or the method first (`integrate(m, fn, a, b)`, unlike `solve`'s planned form).
- **Q4. The orientation stays with the method:** it receives a and b in the caller's order and returns the integral
  from a to b. Recommended: yes. All three methods already do so, TanhSinh's specification does too, and a rule
  written with a signed step (the spike's midpoint rule) gets it for free. Alternative: `integrate` orders the ends
  and negates. That is simpler for a user's method, and agreed bit for bit in the 9,000 cases above, but it changes
  summation orders inside `GaussLegendre` and Romberg, which would need a proof or an accepted change in the last
  bit.
- **Q5. The checkable part of the contract asserted in `integrate`.** Recommended: yes (debug builds and constant
  evaluation; free under `NDEBUG`). Alternatives: no check; or a run-time check that turns a broken result into
  `non_finite`, which would hide a broken method instead of showing it.
- **Q6. The validated criteria passed as `Criteria<T>`, no new type.** Recommended: yes. The contract states what
  is validated. Alternative: a `Tolerances<T>` type whose construction validates, at the cost of one more public
  type.
- **Q7. The protocol as Phase 3a, before TanhSinh (3b).** Recommended: yes. It is the cheapest moment: three methods,
  one dispatch, and TanhSinh becomes the first method added without touching existing code.
- **Q8. `Integral`'s and `Criteria`'s statements about the rounding floor restricted to the library's methods.**
  Recommended: yes (B). Alternative: require the floor of every method. That is unenforceable, and it would make a
  user's method violate the documentation silently.

## Deferred

- The building blocks (`Interval`, `Sum`, `rounding_floor`, the overflow checks, the tables) as public API, each when
  a method needs it (E).
- `Adaptive<Pair>`: the globally adaptive loop over any embedded pair (`ANALYSIS.md` 6.3). Its constants were
  measured for G7-K15, so it waits for a second pair and a numerics review.
- The same protocol for `nxx::roots` and `nxx::deriv` (`ANALYSIS.md` 10), after this phase.

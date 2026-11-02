# Phase 3a API proposal: the method protocol

Status: draft 2, after one review round (numerics-reviewer, simplicity-critic, usage-sketcher; "Decisions from the
reviews", below), for the user. The user's decisions on `docs/open-closed/ANALYSIS.md` (2026-11-02): the methods
open to users (D1), quadrature first (D4). Architecture: `docs/quadrature/BRIEF.md` (approved 2026-10-30), whose D1
("methods as a closed set of empty tags") this phase amends. Phase 3 (`PHASE3_API.md`, TanhSinh, approved
2026-11-01, not yet implemented) becomes Phase 3b and adds `TanhSinh` through this protocol. The evidence is a spike
on a patched copy of `quadrature.hpp` (Spike, below).

## Summary

```cpp
#include <numerixx/quadrature.hpp>
namespace qd = nxx::quadrature;

// Unchanged: every call that compiles today compiles, with the same results (Romberg's floor on reversed ends aside, Q4).
auto const r = qd::integrate(fn, 0.0, 1.0);
auto const s = qd::integrate<qd::Romberg>(fn, 0.0, 1.0, {.rel_tol = 1e-12});
auto const g = qd::integrate<qd::GaussLegendre<20>>(fn, 0.0, 1.0);

// New: a method is a type with one const member, used like the library's.
struct Midpoint
{
    int panels = 100;

    template<nxx::real T, typename F>
    constexpr std::expected<T, qd::error> apply(F &fn, T const &a, T const &b) const; // a fixed rule; a < b
};
static_assert(qd::fixed_rule<Midpoint, double>);

auto const m = qd::integrate<Midpoint>(fn, 0.0, 1.0);
// A method with parameters is passed by value, first, for either kind:
auto const n = qd::integrate(Midpoint{.panels = 1000}, fn, 0.0, 1.0);
auto const t = qd::integrate(SquareAtLowerEnd<qd::GaussKronrod<15>>{}, fn, 0.0, 1.0, {.rel_tol = 1e-12});
```

- Two public concepts, `fixed_rule` and `estimating_method`, each a single const member: `apply(fn, a, b)` and
  `estimate(fn, a, b, criteria)`.
- `GaussLegendre<N>`, `Romberg` and `GaussKronrod<15>` implement them by forwarding to the `detail` functions they
  are dispatched to today.
- `integrate` validates, puts the ends in order (b < a negates), calls the member, and checks the result's finiteness.
  It takes a method either named (`integrate<M>(fn, ...)`, as today) or as a value, first
  (`integrate(m, fn, ...)`).
- The closed lists `detail::gauss_points`, `detail::fixed_rule` and `detail::estimating` are removed, and so is the
  `if constexpr` dispatch on `Romberg`. Adding a method, the library's or a user's, adds a type and edits nothing.
- No new types. Everything constexpr; nothing throws; no ADL (the protocol is a member call).

## Public API (changes to `quadrature.hpp`)

```cpp
namespace nxx::quadrature {

    // ... error, Integral, Failure, Criteria as today; then detail; then:

    /// A fixed rule: M's const member apply(fn, a, b) returns the integral from a to b, or an error (contract: B).
    template<typename M, typename T, typename F = T(T const &)>
    concept fixed_rule = requires(M const &m, F &fn, T const &a, T const &b) {
        { m.apply(fn, a, b) } -> std::same_as<std::expected<T, error>>;
    };

    /// A method with an error estimate: M's const member estimate(fn, a, b, criteria) returns the integral from a
    /// to b with an estimate of its error, or a Failure (contract: B).
    template<typename M, typename T, typename F = T(T const &)>
    concept estimating_method = requires(M const &m, F &fn, T const &a, T const &b, Criteria<T> const &criteria) {
        { m.estimate(fn, a, b, criteria) } -> std::same_as<std::expected<Integral<T>, Failure<T>>>;
    };

    template<int N>
    struct GaussLegendre
    {
        template<real T, typename F>
        requires(N >= 1 && N <= 64)
        constexpr std::expected<T, error> apply(F &fn, T const &a, T const &b) const { /* today's fixed integrate */ }
    };

    struct Romberg { /* estimate(fn, a, b, criteria) const: detail::romberg */ };

    template<int N>
    requires(N == 15)
    struct GaussKronrod { /* estimate(fn, a, b, criteria) const: detail::gauss_kronrod */ };

    // A method as a value, first: any method, so generic code writes one call.
    template<typename M, real T, typename F>
    requires nxx::detail::callable<F, T> && fixed_rule<M, T, std::remove_reference_t<F>>
    [[nodiscard]] constexpr std::expected<T, error> integrate(M const &m, F &&fn, T const &a, T const &b);

    template<typename M, real T, typename F>
    requires nxx::detail::callable<F, T> && estimating_method<M, T, std::remove_reference_t<F>>
    [[nodiscard]] constexpr std::expected<Integral<T>, Failure<T>>
        integrate(M const &m, F &&fn, T const &a, T const &b, Criteria<std::type_identity_t<T>> const &criteria = {});

    // A method named, as today: integrate(M{}, fn, ...), for an M that value-initialises.
    template<typename M, real T, typename F>
    requires std::default_initializable<M> && nxx::detail::callable<F, T> && fixed_rule<M, T, std::remove_reference_t<F>>
    [[nodiscard]] constexpr std::expected<T, error> integrate(F &&fn, T const &a, T const &b);

    template<typename M = GaussKronrod<15>, real T, typename F>
    requires std::default_initializable<M> && nxx::detail::callable<F, T>
             && estimating_method<M, T, std::remove_reference_t<F>>
    [[nodiscard]] constexpr std::expected<Integral<T>, Failure<T>>
        integrate(F &&fn, T const &a, T const &b, Criteria<std::type_identity_t<T>> const &criteria = {});
}
```

- **The three methods move** from the top of the header to just after `detail`, before `integrate`, with their
  Doxygen comments, and define their members inline. Nothing in `detail` names a method once the lists are gone. The
  site's reference page orders the declarations itself. Draft 1 kept them in place and defined the members out of
  line, a pattern the library does not use, with repeated requires-clauses that MSVC and Doxygen are known to
  mismatch (simplicity review).
- **The concepts take `<M, T, F>`**, the order of `integrate<M, T, F>`. F defaults to a plain function of T, so
  `static_assert(qd::fixed_rule<MyRule, double>)` beside a method's definition is the usual check (usage sketch).
  F is the type `integrate` passes, reference removed: `integrate` already holds F and T, so this is the simplest
  check, not a feature. A combinator's inner method sees another fn; it checks against the outer F as a stand-in.
  No `std::copy_constructible`: `integrate` never copies a method, and adding the requirement later would break
  users, while dropping it later would not.
- **The named forms keep today's signatures exactly** and forward to the value forms with `M{}`. They require
  `std::default_initializable<M>`, so a method without a default is rejected softly, the same on every compiler.
  A braced list can never become a method: `integrate<M>(fn, a, b, {1e-6})` is criteria or nothing.
- **The value forms put the method first**, so one call, `integrate(m, fn, a, b)`, works for both kinds. Generic code
  takes any method without a branch on its kind. A method given with criteria: `integrate(m, fn, a, b, criteria)`.
- **`GaussLegendre<N>` stays an unconstrained struct**, its range a requires-clause on the member, so that
  `GaussLegendre<0>` and `<65>` are types that `integrate` rejects (the existing `!integrable<...>` tests).
- **A method is of one kind.** A type with both members makes the named call ambiguous; said in the concepts'
  comment, not prevented.

## Behaviour

### A. What `integrate` does

In order:
- The tolerances, for an estimating method (negative or not finite: `non_finite`), then `rel_tol` clamped to at most
  1. The ends (NaN or infinite: `non_finite`). a == b: 0, or `{0, 0}`, without calling fn. All as today.
- **The ends in order** (Q4): the method receives `(min, max)`; where b < a, `integrate` negates the value, and a
  best's value. Today each method orients itself; Gauss-Legendre and Gauss-Kronrod give the exact negation, so they
  are unchanged, and Romberg's `abs_error` on reversed ends changes where its plain L1 sum, summed in the other
  order, moves the floor (Spike).
- **The method's member**, with the criteria validated and clamped, as `Criteria<T>` (Q6).
- **The result checked** (Q5), in every build: a fixed rule's value must be finite; an `Integral` (the result, or a
  best) must have a finite value and a finite, non-negative `abs_error`, and a best may come only with
  `not_converged`. Anything else is returned as `non_finite` without a best. The library's methods never trip it
  (numerics review: no counterexample in any of the four, TanhSinh's specification included).

### B. The contract

Written once, in the concepts' Doxygen comments (the reference page shows it through Breathe; `design.md` records the
decision and points there).

- **Given to the method:** a < b, both finite. They may be extreme: b - a may overflow ([-max, max]), the ends may
  be subnormal, so that halves round ([D, 3D], D = denorm_min), adjacent (no T strictly between them), or of
  opposite signs, so that b - a rounds. For an estimating method, criteria with both members finite and non-negative,
  `rel_tol <= 1`, `abs_tol` unclamped (it enters only a maximum), and possibly both 0. fn as `integrate` received
  it, as an lvalue (const if the caller's was): call it with `std::invoke` and an argument of type T, and convert the
  result with `T(...)` (the module's rule on fn's result, `docs/design.md`).
- **Required of every method:**
  - fn is evaluated only at points of [a, b]. A point formed as a + h t can round outside, so clamp it.
  - It terminates for every criteria it can be given, both tolerances 0 included: a cap, or a floor of its own.
  - A finite result (checked, A), or an error: `non_finite` without a best, or `not_converged` with or without one.
  - Its parameters are its own preconditions (Q9): checked by the method, or made unrepresentable when invalid.
  - fn's exceptions, and the method's own, pass through: the library catches nothing.
- **Required in meaning, not checkable:** an estimating method's `abs_error` estimates |value - exact| (an
  estimate, not a bound, as `Integral` says), covering the rounding of its own sum; success means its own stop test
  against `max(abs_tol, rel_tol |value|)`, or its floor, was met.
- **The library's methods only** (Q8): the rounding floor (abs_error never below 50 max(eps L1, denorm_min) for
  a != b, which also corrects today's wording, false at a == b), `abs_tol = 0` working, a best with every
  `not_converged` (Romberg and GaussKronrod; TanhSinh omits it where its estimate is unbounded), fn never evaluated
  at an end (GaussLegendre, TanhSinh), and `non_finite` in constant evaluation where an overflow or a NaN would
  otherwise be a compile error. `Integral`, `Failure`, `Criteria`, `error` and `integrate`'s comments say so.
- **Infinite ends** stay outside this contract (Q10): a method written to it never receives one. Phase 4
  (exp-sinh, sinh-sinh) makes them opt-in per method.

### C. The library's methods

- `GaussLegendre<N>::apply` is the body of today's fixed `integrate` after its checks: `detail::interval_of`, then
  `detail::with_gauss_legendre<N, T>` and `detail::apply`.
- `Romberg::estimate` is `detail::romberg(fn, a, b, criteria.abs_tol, criteria.rel_tol)`.
- `GaussKronrod<N>::estimate` is `detail::with_gauss_kronrod<T>` around `detail::gauss_kronrod`.
- Nothing inside `detail` changes. Each method still handles either order of its ends, though it now receives a < b;
  the code that does it stays, untouched.

### D. Phase 3b (TanhSinh) through the protocol

`PHASE3_API.md` stays approved as it is, except:
- **Its "Public API" bullet:** `TanhSinh` implements `estimate`, so the existing estimating forms of `integrate` take
  it; no list, no dispatch, no overload changes (it said that `detail::estimating<M>` admits `TanhSinh`).
- **Its A:** "b < a negates" is now done by `integrate`, so TanhSinh receives a < b.

Its behaviour otherwise, the decisions Q1 to Q10 and the tests are unchanged. `docs/plan.md` splits Phase 3 into 3a
(this protocol) and 3b (TanhSinh, `PHASE3_API.md`).

### E. What a method of one's own cannot use yet

The building blocks that make the library's methods robust stay in `detail`: `Interval` (the overflow-free map and
the order of each term's products), `Sum` (Neumaier), `rounding_floor`, the overflow checks for constant evaluation,
and the node tables. A user's method works without them, as the spike's do, but does not inherit their guarantees.
The guide page gives the recipe for a floor (50 eps L1 from the method's own samples, `Criteria`'s), so that a user's
estimating method can end at `rel_tol = 0`. Publishing the blocks waits until a method needs one (ANALYSIS.md, 6.4).

## Spike

`docs/open-closed/quadrature-spike/` (README there): `quadrature.patch.py` writes this proposal into a copy of
`quadrature.hpp`, with `normalize` for Q4's recommendation or without it for the alternative. Built with GCC 13.3 and
Boost 1.83 (the container's; the project pins 1.92, whose download the container's network refuses), not on the
project's five toolchains; the phase itself builds and tests on all five.

- **The existing quadrature tests** (`test_quadrature_gauss.cpp`, `_romberg.cpp`, `_kronrod.cpp`), built unchanged
  against both variants, pass exactly as against today's header: 25 of 26 test cases in each build. The one failure
  is the same in all, the `GaussLegendre<2>` long double node baseline, recorded on the project's toolchains, which
  does not hold for GCC 13's 80-bit long double.
- **Bit for bit** (`dump.cpp`): 8 integrands (e^x, sqrt x, 1/sqrt x, log x, a kink, a peak of width 1e-3, sin 50x,
  NaN on half the interval) on 6 intervals, three of them reversed, by `GaussLegendre<1>`, `<7>`, `<20>` and `<64>`,
  `GaussKronrod<15>` at the default, 1e-12 and 0, and Romberg at the default, 1e-10 and an oversized tolerance, in
  float, double, long double, `cpp_bin_float_50` and `cpp_dec_float_50`; edges (1 and 1/4 on [-max, max], a == b, a
  NaN end, a negative tolerance) and three results computed in constant expressions: 2,438 results. Draft 1:
  identical. Draft 2 without ordering the ends: identical. Draft 2 with it (Q4): DUMP_NORMALIZE (run pending).
- **Orientation today** (`orient.cpp`): 3,000 intervals (positive, negative, across 0, up to 2e6 wide) and three
  integrand families, each integrated forwards and reversed, the whole result compared (status, value, abs_error,
  best). Gauss-Legendre<64> and Gauss-Kronrod<15> (at the default and at rel_tol 0): the exact negation every time.
  Romberg: the exact negation in value and status every time, a different `abs_error` in 6 runs at the default and
  368 at rel_tol 0 (374 of the 6,000 runs), from the floor's plain L1 sum.
  Draft 1's claim, "held bit for bit", compared values only.
- **Methods of one's own** (`user_methods.cpp`, built with `-Wall -Wextra -Wpedantic -Wshadow -Wconversion`, no
  warnings, with GCC 13.3 and with Clang 18, which needs `-U__cpp_concepts -D__cpp_concepts=202002L` for libstdc++
  13's `<expected>`): a fixed rule with a run-time parameter (`Midpoint{.panels}`), clamped and refusing a width that
  overflows, so fn stays in [-max, max] on [-max, max]; a combinator over any estimating method
  (`SquareAtLowerEnd<Inner>`, x = a + (b - a) u^2); in constant expressions and at run time. A spy method records
  that it receives the ends in order, `rel_tol = 5` as 1 and `abs_tol = 1e300` unclamped, and that a NaN end, a
  negative tolerance and a == b never call it. A method returning a NaN value or a negative `abs_error` gives
  `non_finite`, with and without `NDEBUG`. The concepts reject a wrong return type, a non-const member, the other
  kind, a type of both kinds when named, `GaussLegendre<0>`, `<65>`, `int`, and a braced list as a method.
- **The substitution** at rel_tol 1e-12, against `GaussKronrod<15>` alone: on 1/sqrt x over [0, 1] (which it turns
  into the constant 2, its best case) 150 evaluations, exact, against 2,370 and an error of 3.0e-13; reversed, the
  same, negated; on cos x / sqrt x (g = 2 cos u^2; reference 1.80904847580054416295 from mpmath, by quadrature and by
  Fresnel's C) 150 evaluations, an error of 2.2e-16, against 2,400 and 2.1e-13. Over Romberg it is `non_finite`, as
  Romberg is documented to be: it samples u = 0, where the substituted integrand is 0 times inf.

## Tests (`tests/test_quadrature_protocol.cpp`, test cases `quadrature:protocol:*`)

- **The existing suite unchanged**, passing with its exact counts and values, except any check of Romberg's
  `abs_error` on reversed ends (none in the spike's run; Q4).
- **Both forms agree:** `integrate<M>(fn, a, b[, c]) == integrate(M{}, fn, a, b[, c])` for the three library
  methods, one `static_assert` each in double, with run-time twins.
- **Reversal is exact for every method:** `integrate(fn, b, a) == -integrate(fn, a, b)`, the whole result, for the
  three library methods and the user methods, including Romberg at rel_tol 0 on a case of `orient.cpp` that differs
  today.
- **What a method receives** (the spy): the ends in order; `rel_tol = 5` as 1; `abs_tol = 1e300` unclamped; never
  called on a NaN end, a negative tolerance or a == b.
- **A user fixed rule with a parameter** (the composite midpoint rule) against closed forms: on x^2 over [0, 3] with
  3 panels the rule's error, (b - a) h^2 f'' / 24 = 1/4, is exact for a quadratic and every node is dyadic, so the
  value is exactly 8.75 in double, -8.75 reversed; `static_assert` with run-time twins. On [-max, max] it is
  `non_finite` and fn, checked by the existing `Confined` callback, is never called outside the interval.
- **A user combinator over a library method** (`SquareAtLowerEnd<GaussKronrod<15>>`), in double and
  `cpp_bin_float_50` (one of them with expression templates): 1/sqrt x over [0, 1] = 2, where g(u) = 2u fl(1 /
  fl(sqrt(fl(u^2)))) is 2 within about 3 ulps per evaluation, so |value - 2| <= 8 eps (K15's 15 weights, the
  compensated sum); cos x / sqrt x over [0, 1] against the mpmath reference within `abs_error` and the tolerance.
- **The checked result:** a method returning a NaN value, a negative `abs_error`, or a best with `non_finite` gives
  `non_finite`.
- **The concepts:** a wrong return type, a non-const member, the other kind, `int`, a type of both kinds when named;
  criteria for a fixed rule still do not compile (the existing `with_criteria`); a braced list is never a method.
  No test for a method that cannot be value-initialised beyond `!std::default_initializable` in the named forms.

## Documentation

- `quadrature.hpp`: the file comment (a method is a type, and a user's works the same way; both forms); the concepts'
  comments: the contract of B, "a method is of one kind", the usual causes of a failed concept (a missing `const`,
  the exact return type, the wrong kind for the call) and `static_assert(qd::fixed_rule<MyRule, double>)` beside a
  method; each member's comment; `integrate`'s `@tparam M`, `@param m` and an example of each form; `Integral`,
  `Failure`, `Criteria` and `error::non_finite` as in B (the last: "or a method reported its own failure").
- `docs/quadrature/design.md`: "The interface" rewritten for the protocol, pointing to the concepts for the contract;
  BRIEF.md D1 marked as amended (the set is open).
- `docs/sphinx/quadrature/reference.md`: the two concepts and the value forms.
- The guide page on writing a method, for the milestone's chapter on the site, with `user_methods.cpp` as an
  example built as a test: written as an engineer writes it (T(...) only where an expression-template type needs
  it, `std::invoke` for member pointers), the clamp, a cap and the floor's recipe (E), state through a pointer member
  (the member is const), the reversal check.
- `docs/plan.md`: Phase 3 split into 3a and 3b; the "Customization" item on an open method protocol points here for
  quadrature.

## Decisions for the user

- **Q1. Member names `apply` (fixed rule) and `estimate` (estimating method).** Recommended: yes. The name says which
  kind a method is. Alternative: `integrate` for both, told apart by arguments and return type, which invites calling
  the member directly and skipping the validation.
- **Q2. Concepts `fixed_rule<M, T, F = T(T const &)>` and `estimating_method<M, T, F = ...>`, with no copy
  requirement.** Recommended: yes (changed from draft 1's `<M, F, T>` and `std::copy_constructible`). Alternative:
  `<M, T>` alone, checked with a stand-in fn in `detail` that users cannot name.
- **Q3. A method value goes first, `integrate(m, fn, a, b[, criteria])`, beside today's named forms, unchanged.**
  Recommended: yes (changed from draft 1, which put it last). One call for both kinds, so generic code needs no
  branch; a braced list never becomes a method; a method that cannot be value-initialised fails softly. Cost: four
  overloads instead of two (the named two forward in one line each). Alternative: draft 1's defaulted last
  parameter: two overloads, but an estimating method by value needs `{}` before it, and `integrate<M>(fn, a, b,
  {1e-6})` silently builds a user's M whose first member is a floating type.
- **Q4. `integrate` puts the ends in order and negates.** Recommended: yes (changed from draft 1). Every method
  receives a < b: the contract is simpler, a method tied to an end knows which, and two mistakes the usage sketch
  made that nothing catches (a loop `x < b` from a, an `if (b < a) swap`) cannot happen. Reversal becomes exact for
  every method. Cost: Romberg's `abs_error` on reversed ends changes where the floor's L1 sum rounds differently
  (values and statuses never, in 6,000 runs; DUMP_NORMALIZE). Alternative: the method orients itself, bit for bit with
  today, with the orientation as the first rule of the contract and a reversal check in the guide.
- **Q5. The result's finiteness checked in every build, a broken result returned as `non_finite`.** Recommended: yes
  (changed from draft 1's assert). `integrate`'s comment promises `non_finite` where a value is not finite; under
  `NDEBUG` an assert would return a user's NaN as a success. It costs a few comparisons per call, and it is tested.
  Alternatives: the assert of draft 1; or no check, the contract as documentation only (simplicity review).
- **Q6. The validated criteria passed as `Criteria<T>`, no new type.** Recommended: yes.
- **Q7. The protocol as Phase 3a, before TanhSinh (3b).** Recommended: yes.
- **Q8. The statements that hold only for the library's methods say so** (B: the floor, `abs_tol = 0`, a best with
  every `not_converged`, never at an end, constant evaluation), and the floor's wording is corrected for a == b.
  Recommended: yes. Alternative: require them of every method; the floor needs the method's own L1, so `integrate`
  could not check it.
- **Q9. A method's parameters are its own preconditions** (checked by the method, asserted as the library asserts
  its preconditions, or made unrepresentable: a template parameter, an unsigned type); no new enumerator now, and
  `non_finite`'s comment says a method may report its own failure. Recommended: yes; an enumerator when a library
  method has a run-time parameter. Alternative: `error::invalid_method` now.
- **Q10. Finite ends only, and infinite ends opt-in per method in Phase 4.** Recommended: yes. A method written to
  this contract then never receives an infinite end. Phase 4 chooses the opt-in (a member, or a third concept).

## Decisions from the reviews

From the numerics review:
- Adopted: fn evaluated only in [a, b] required of every method, and the spike's methods clamp and refuse an
  overflowing width (B; draft 1's midpoint rule called fn at +inf on [-max, max], and its substitution at b + 2 on
  [-1, 2^53 + 2]); the orientation of methods tied to an end, answered by Q4; termination for every criteria, both
  tolerances 0 included, and parameters as the method's preconditions (B, Q9); the edge inputs listed under "Given",
  and constant evaluation under "the library's methods only" (B); Q8 extended (`integrate`'s @return, `abs_tol = 0`,
  `Failure`'s best) and the a == b inconsistency of the floor's wording fixed; the tests' tolerances derived (8 eps
  for the substitution of 1/sqrt x), a spy for what a method receives and a method, not fn, that must not be called,
  an expression-template instantiation, a reversed combinator, the domain check; the spike's presentation (1/sqrt x
  named the best case, cos x / sqrt x added; smoothness only in x - a; "an lvalue"); `orient.cpp` compares whole
  results, over intervals of both signs, which found Romberg's `abs_error`.
- Confirmed: the checked part holds for all four methods, the edges and constant evaluation; the members pass the same
  arguments; `rel_tol <= 1` and an unclamped `abs_tol`; the midpoint and substitution values; Q4's derivation.
From the simplicity review:
- Adopted: the methods moved after `detail` with inline members (no out-of-line constrained definitions); no
  `std::copy_constructible`; the contract written once, in the header; "a method is of one kind"; the tests trimmed
  (one `static_assert` per library method for the two forms; a spy for the clamped criteria); "an lvalue"; Q2's
  reason restated.
- Not adopted: Q5's "no check", in favour of the release check above (the usage sketch showed a NaN returned as a
  success under `NDEBUG`); dropping the second integrand of the substitution (the numerics review asked for a case
  that is not the best case).
From the usage sketch:
- Adopted: the method value first, for both kinds (Q3; draft 1 rejected an estimating method by value without
  criteria, and built a user's method from a braced tolerance); the ends ordered by `integrate` (Q4); the release
  check (Q5); parameters (Q9); infinite ends recorded (Q10); `<M, T, F>` with F defaulted (Q2); no copy requirement,
  and state through a pointer member in the guide; the concepts' comments on the usual diagnostics and a
  `static_assert` beside a method; the floor's recipe and the guide's style (E, Documentation).
- Not adopted: a precedence between the kinds for a type of both (a sentence instead); publishing `rounding_floor`
  now (a recipe instead, until a method needs the block).

## Deferred

- The building blocks (`Interval`, `Sum`, `rounding_floor`, the overflow checks, the tables) as public API, each when
  a method needs it (E).
- Infinite ends as an opt-in (Q10): Phase 4.
- `Adaptive<Pair>`: the globally adaptive loop over any embedded pair (`ANALYSIS.md` 6.3); its constants were measured
  for G7-K15, so it waits for a second pair and a numerics review.
- The same protocol for `nxx::roots` and `nxx::deriv` (`ANALYSIS.md` 10), after this phase, with the method value
  first there too.

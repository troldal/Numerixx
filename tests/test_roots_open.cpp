// Newton, Halley and Secant in nxx::roots (Phase 2, docs/roots/PHASE2_API.md): from a guess, with the hand-over to
// the bracket loop at the first sign change, and Newton and Halley safeguarded by a bracket; both derivative forms.
//
// Evaluation counts are asserted only where no libm is involved (+ - * / in double), against an independent Python
// implementation of the specification (.scratch/roots2/proto2.py), assuming no FMA contraction in fn or in the steps
// (true on x64 and wasm, the CI's targets). Every other case checks accuracy (roots_check.hpp).

#include "roots_check.hpp"

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <format>
#include <limits>
#include <tuple>
#include <utility>

using namespace roots_check;

namespace {

    // Each case as fn returning (f, f', f''): Newton uses the first two, Secant the first, Brent and Bisection the
    // first. The expressions are those of the prototype, in the same order, so that the counts match.
    constexpr auto sq = [](double x) { return std::tuple{ x * x - 2.0, 2.0 * x, 2.0 }; };
    constexpr auto wal = [](double x) { return std::tuple{ (x * x - 2.0) * x - 5.0, 3.0 * x * x - 2.0, 6.0 * x }; };
    constexpr auto triple = [](double x) {
        return std::tuple{ (x - 1.0) * (x - 1.0) * (x - 1.0), 3.0 * (x - 1.0) * (x - 1.0), 6.0 * (x - 1.0) };
    };
    constexpr auto x20 = [](double x) {
        double const x2 = x * x;
        double const x4 = x2 * x2;
        double const x8 = x4 * x4;
        double const x16 = x8 * x8;
        return std::tuple{ x8 * x8 * x4 - 1.0, 20.0 * x16 * x2 * x, 380.0 * x16 * x2 };
    };
    constexpr auto cycle = [](double x) { return std::tuple{ (x * x - 2.0) * x + 2.0, 3.0 * x * x - 2.0, 6.0 * x }; };
    constexpr auto reciprocal = [](double x) { return std::tuple{ 1.0 / x - 1.0, -1.0 / (x * x), 2.0 / (x * x * x) }; };
    constexpr auto cubed = [](double x) { return std::tuple{ x * x * x, 3.0 * x * x, 6.0 * x }; };
    constexpr auto line = [](double x) { return std::tuple{ x - 1.0, 1.0, 0.0 }; };
    constexpr auto steep = [](double x) { return std::tuple{ 1e300 * (x - 1.0), 1e300, 0.0 }; };

    // Wilkinson's W10 = (x - 1) ... (x - 10), expanded and evaluated by Horner's scheme with its derivatives: a simple
    // root at 7 whose computed values are rounding noise within about 1e-10 of it.
    constexpr auto w10 = [](double x) {
        constexpr std::array<double, 11> c{
            1.0, -55.0, 1320.0, -18150.0, 157773.0, -902055.0, 3416930.0, -8409500.0, 12753576.0, -10628640.0, 3628800.0
        };
        double p = 0.0;
        double dp = 0.0;
        double d2 = 0.0;
        for (double const a : c) {
            d2 = d2 * x + 2.0 * dp;
            dp = dp * x + p;
            p = p * x + a;
        }
        return std::tuple{ p, dp, d2 };
    };

    constexpr double cycle_root = -1.769292354238631415240409; // x^3 - 2x + 2 (mpmath)
    constexpr double kepler_root = 1.862086686874532271834067; // E - 0.9 sin E = 1, with 0.9 as a double (mpmath)

    template<typename F>
    bool sign_change_within(F const &fn, rt::Root<double> const &r)
    {
        // The bracket's other end is at value - abs_error or value + abs_error: fn there has the sign opposite to fn at
        // value (or one is 0). The ends value -+ abs_error themselves may share a sign, since a noisy fn can change
        // sign several times in between.
        double const at = std::get<0>(fn(r.value));
        double const lo = std::get<0>(fn(r.value - r.abs_error));
        double const hi = std::get<0>(fn(r.value + r.abs_error));
        return at == 0.0 || lo == 0.0 || hi == 0.0 || (lo < 0.0) != (at < 0.0) || (hi < 0.0) != (at < 0.0);
    }

} // namespace

// Constant evaluation: each method from a guess and on a bracket, and the overflow checks (no infinity is formed).
static_assert(nxx::math::abs(rt::solve<rt::Newton>(sq, 1.0)->value - sqrt2) <= 1e-15);
static_assert(nxx::math::abs(rt::solve<rt::Halley>(sq, 1.0)->value - sqrt2) <= 1e-15);
static_assert(nxx::math::abs(rt::solve<rt::Secant>(sq, 1.0)->value - sqrt2) <= 1e-15);
static_assert(nxx::math::abs(rt::solve<rt::Newton>(wal, rt::Bracket{ 2.0, 3.0 })->value - wallis) <= 1e-15);
static_assert(nxx::math::abs(rt::solve<rt::Halley>(wal, rt::Bracket{ 2.0, 3.0 })->value - wallis) <= 1e-15);
static_assert(rt::solve<rt::Newton>(line, 1e308) == std::unexpected(rt::error::not_converged)); // a step above max / 2
static_assert(rt::solve<rt::Secant>(steep, 1e7)->value == 1.0); // the ratio first: f (x - xp) would overflow
static_assert(rt::solve<rt::Halley>([](double x) { return std::tuple{ 1e308 * (x - 1.0), 1e308, 1e308 }; }, 2.0)->value
              == 1.0); // f'' / f' and t f'' / (2 f') near max
static_assert(
    nxx::math::abs(
        rt::solve<rt::Newton>([](float x) { return std::pair{ x * x - 2.0f, 2.0f * x }; }, 1.0f)->value - 1.41421356f)
    <= 1e-6f);
static_assert(
    nxx::math::abs(
        rt::solve<rt::Halley>([](long double x) { return std::tuple{ x * x - 2.0L, 2.0L * x, 2.0L }; }, 1.0L)->value
        - 1.41421356237309504880L)
    <= 1e-15L);
static_assert(rt::solve<rt::Newton>(cycle, 0.0) == std::unexpected(rt::error::not_converged));

// The constraints: each method accepts only its starts and the values it needs.
namespace {

    constexpr auto value_only = [](double x) { return x * x - 2.0; };
    constexpr auto pair_fn = [](double x) { return std::pair{ x * x - 2.0, 2.0 * x }; };

    template<typename M, typename F, typename S>
    concept solvable = requires(F f, S s) { rt::solve<M>(f, s); };

    template<typename M, typename F, typename D, typename S>
    concept solvable_with = requires(F f, D d, S s) { rt::solve<M>(f, d, s); };

    template<typename M, typename F, typename D, typename D2, typename S>
    concept solvable_with_two = requires(F f, D d, D2 d2, S s) { rt::solve<M>(f, d, d2, s); };

    using Value = decltype(value_only);
    using Pair = decltype(pair_fn);
    using Triple = decltype(sq);
    using B = rt::Bracket<double>;

} // namespace

static_assert(
    solvable<rt::Newton, Pair, double> && solvable<rt::Newton, Pair, B> && solvable<rt::Secant, Value, double>);
static_assert(solvable<rt::Halley, Triple, double> && solvable<rt::Halley, Triple, B>);
static_assert(solvable<rt::Brent, Pair, B> && solvable<rt::Bisection, Triple, B>); // a tuple fn's first element
static_assert(!solvable<rt::Newton, Value, double>); // no derivative
static_assert(!solvable<rt::Halley, Pair, double>); // no second derivative
static_assert(!solvable<rt::Bisection, Value, double>); // a bracketing method from a guess
static_assert(!solvable<rt::Secant, Value, B>); // Secant on a bracket is Brent
static_assert(solvable_with<rt::Newton, Value, Value, double> && solvable_with<rt::Newton, Value, Value, B>);
static_assert(!solvable_with<rt::Newton, Pair, Value, double>); // a mixed form
static_assert(!solvable_with<rt::Halley, Value, Value, double>); // Halley needs d2fn
static_assert(solvable_with_two<rt::Halley, Value, Value, Value, double>);
static_assert(!solvable_with_two<rt::Newton, Value, Value, Value, double>);
// Elements exactly T: a reference (std::tie of fn's locals) or another type is rejected, not converted.
static_assert(!solvable<rt::Newton, double (*)(double), double>);
static_assert(!solvable<rt::Newton, std::tuple<double &, double &> (*)(double const &), double>);
static_assert(!solvable<rt::Newton, std::pair<float, float> (*)(double const &), double>);
static_assert(!solvable<rt::Newton, std::pair<double, double> const &(*)(double const &), double>);
static_assert(solvable<rt::Newton, std::array<double, 2> (*)(double const &), double>);

TEST_CASE("roots: open methods, evaluation counts in double")
{
    // The prototype's counts (proto2.py), for the methods from a guess and on a bracket.
    struct Case
    {
        char const *name;
        int newton;
        int halley;
        int secant;
    };
    auto const check = [](auto const &fn, double x0, int newton, int halley, int secant) {
        CHECK(counted<rt::Newton>(fn, x0).second == newton);
        CHECK(counted<rt::Halley>(fn, x0).second == halley);
        if (secant > 0) CHECK(counted<rt::Secant>(fn, x0).second == secant);
    };
    check(sq, 1.0, 7, 5, 9); // the early hand-over: 1.5 is past the root
    check(sq, 10.0, 9, 7, 13); // the late hand-over: monotone from above, crossing only at the minimum step
    check(sq, 1e-3, 17, 12, 25);
    check(sq, 1.4142135623730951, 2, 2, 4); // a guess at the root: one minimum step past it
    check(wal, 2.0, 6, 5, 8);
    check(wal, 0.0, 21, 12, 15);
    check(wal, 10.0, 11, 8, 15);
    check(triple, 2.0, 89, 53, 127); // a triple root: linear convergence
    check(x20, 0.5, 22, 11, 44);
    check(x20, 2.0, 20, 11, 28);
    check(cycle, -2.0, 7, 5, 9);
    check(reciprocal, 0.5, 7, 2, 10);
    check(steep, 1e7, 2, 2, 4);
    check(w10, 7.2, 21, 19, 23);
    check(w10, 6.8, 22, 16, 23);

    auto const bracket = [](auto const &fn, rt::Bracket<double> b, int newton, int halley) {
        CHECK(counted<rt::Newton>(fn, b).second == newton);
        CHECK(counted<rt::Halley>(fn, b).second == halley);
    };
    bracket(sq, { 1.0, 2.0 }, 8, 6);
    bracket(wal, { 2.0, 3.0 }, 7, 6);
    bracket(triple, { 0.0, 3.0 }, 91, 54);
    bracket(x20, { 0.0, 5.0 }, 11, 11);
    bracket(cycle, { -2.0, -1.0 }, 8, 6);
    bracket(reciprocal, { 0.5, 3.0 }, 9, 6);
    bracket(cubed, { -1.0, 2.0 }, 3, 3); // zero first
    bracket(line, { 0.0, 1e308 }, 7, 7);
    bracket(w10, { 6.7, 7.6 }, 28, 24);
    CHECK(counted<rt::Brent>(w10, rt::Bracket{ 6.7, 7.6 }).second == 23); // Brent on a tuple fn: its first element
}

TEST_CASE("roots: open methods, accuracy in double")
{
    // The conditioning terms c, derived as in Phase 1. x^2 - 2: 7.9e-17. Wallis: 2e-16. (x - 1)^3, x - 1, 1e300 (x - 1)
    // and x^20 - 1 near 1: x - 1 is exact by Sterbenz, so the sign change is at 1 (x^20 - 1: about 19 half-ulp
    // roundings of x^20 near 1, 2.1e-15, over f' = 20: 1.1e-16). x^3 - 2x + 2 near -1.77: about 11.5 half-ulp roundings
    // (1.3e-15) over |f'| = 7.39: 1.8e-16. 1/x - 1: fl(1/x) is 1 up to half an ulp below 1: 1.2e-16.
    for (double const x0 : { 1.0, 10.0, 1e-3 }) {
        CHECK(brackets(rt::solve<rt::Newton>(sq, x0), sqrt2, 8e-17));
        CHECK(brackets(rt::solve<rt::Halley>(sq, x0), sqrt2, 8e-17));
        CHECK(brackets(rt::solve<rt::Secant>(sq, x0), sqrt2, 8e-17));
    }
    for (double const x0 : { 2.0, 0.0, 10.0 }) {
        CHECK(brackets(rt::solve<rt::Newton>(wal, x0), wallis, 2e-16));
        CHECK(brackets(rt::solve<rt::Halley>(wal, x0), wallis, 2e-16));
        CHECK(brackets(rt::solve<rt::Secant>(wal, x0), wallis, 2e-16));
    }
    CHECK(brackets(rt::solve<rt::Newton>(triple, 2.0), 1.0, 0.0));
    CHECK(brackets(rt::solve<rt::Halley>(triple, rt::Bracket{ 0.0, 3.0 }), 1.0, 0.0));
    CHECK(brackets(rt::solve<rt::Secant>(x20, 0.5), 1.0, 1.1e-16));
    CHECK(brackets(rt::solve<rt::Newton>(cycle, -2.0), cycle_root, 1.8e-16));
    CHECK(brackets(rt::solve<rt::Halley>(reciprocal, 3.0), 1.0, 1.2e-16));
    CHECK(brackets(rt::solve<rt::Newton>(line, rt::Bracket{ 0.0, 1e308 }), 1.0, 0.0));

    // libm cases, with 2-ulp budgets: cos x - x (1.4e-16, Phase 1), e^x - 2 (4.5e-16), x e^x - 1 (1.3e-16), and Kepler:
    // 2 ulps of sin near 0.96 times 0.9, the product's and the difference's half ulps (3.6e-16) over
    // f' = 1 - 0.9 cos E = 1.26: 2.9e-16.
    auto const dottie_fn = [](double x) { return std::tuple{ std::cos(x) - x, -std::sin(x) - 1.0, -std::cos(x) }; };
    auto const ln2_fn = [](double x) { return std::tuple{ std::exp(x) - 2.0, std::exp(x), std::exp(x) }; };
    auto const omega_fn = [](double x) {
        return std::tuple{ x * std::exp(x) - 1.0, (1.0 + x) * std::exp(x), (2.0 + x) * std::exp(x) };
    };
    auto const kepler_fn = [](double x) {
        return std::tuple{ x - 0.9 * std::sin(x) - 1.0, 1.0 - 0.9 * std::cos(x), 0.9 * std::sin(x) };
    };
    for (double const x0 : { 0.0, 1.0, 3.0 }) {
        CHECK(brackets(rt::solve<rt::Newton>(dottie_fn, x0), dottie, 1.4e-16));
        CHECK(brackets(rt::solve<rt::Halley>(dottie_fn, x0), dottie, 1.4e-16));
        CHECK(brackets(rt::solve<rt::Secant>(dottie_fn, x0), dottie, 1.4e-16));
        CHECK(brackets(rt::solve<rt::Newton>(kepler_fn, x0), kepler_root, 2.9e-16));
        CHECK(brackets(rt::solve<rt::Halley>(kepler_fn, x0), kepler_root, 2.9e-16));
        CHECK(brackets(rt::solve<rt::Secant>(kepler_fn, x0), kepler_root, 2.9e-16));
    }
    for (double const x0 : { 0.0, 5.0 }) {
        CHECK(brackets(rt::solve<rt::Newton>(ln2_fn, x0), ln2, 4.5e-16));
        CHECK(brackets(rt::solve<rt::Halley>(ln2_fn, x0), ln2, 4.5e-16));
        CHECK(brackets(rt::solve<rt::Secant>(ln2_fn, x0), ln2, 4.5e-16));
        CHECK(brackets(rt::solve<rt::Newton>(omega_fn, x0 / 2.5), omega, 1.3e-16));
        CHECK(brackets(rt::solve<rt::Halley>(omega_fn, x0 / 2.5), omega, 1.3e-16));
    }
    CHECK(brackets(rt::solve<rt::Newton>(dottie_fn, rt::Bracket{ 0.0, 1.0 }), dottie, 1.4e-16));
    CHECK(brackets(rt::solve<rt::Halley>(kepler_fn, rt::Bracket{ 0.0, 3.141592653589793 }), kepler_root, 2.9e-16));

    // W10: the bound describes the computed fn, whose sign changes lie up to about 1e-10 from 7 (a Horner error bound
    // gives c = 7e-8, too loose to catch a regression), so the sign change itself is checked.
    for (auto const &r : { rt::solve<rt::Newton>(w10, 7.2),
             rt::solve<rt::Halley>(w10, 6.8),
             rt::solve<rt::Secant>(w10, 7.2),
             rt::solve<rt::Newton>(w10, rt::Bracket{ 6.7, 7.6 }) }) {
        REQUIRE(r.has_value());
        CHECK(sign_change_within(w10, *r));
        CHECK(r->abs_error <= 2.0 * default_tol(r->value));
    }
}

TEST_CASE("roots: the two derivative forms give the same result")
{
    // Bit for bit: the separate callables are joined into the tuple form, evaluated in the same order.
    auto const f = [](double x) { return (x * x - 2.0) * x - 5.0; };
    auto const df = [](double x) { return 3.0 * x * x - 2.0; };
    auto const d2f = [](double x) { return 6.0 * x; };
    for (double const x0 : { 2.0, 0.0, 10.0 }) {
        CHECK(rt::solve<rt::Newton>(f, df, x0) == rt::solve<rt::Newton>(wal, x0));
        CHECK(rt::solve<rt::Halley>(f, df, d2f, x0) == rt::solve<rt::Halley>(wal, x0));
    }
    CHECK(rt::solve<rt::Newton>(f, df, rt::Bracket{ 2.0, 3.0 }) == rt::solve<rt::Newton>(wal, rt::Bracket{ 2.0, 3.0 }));
    CHECK(rt::solve<rt::Halley>(f, df, d2f, rt::Bracket{ 2.0, 3.0 })
          == rt::solve<rt::Halley>(wal, rt::Bracket{ 2.0, 3.0 }));

    // A pair, a tuple of three and an array give Newton the same values.
    auto const as_pair = [](double x) { return std::pair{ (x * x - 2.0) * x - 5.0, 3.0 * x * x - 2.0 }; };
    auto const as_array = [](double x) { return std::array{ (x * x - 2.0) * x - 5.0, 3.0 * x * x - 2.0 }; };
    CHECK(rt::solve<rt::Newton>(as_pair, 0.0) == rt::solve<rt::Newton>(wal, 0.0));
    CHECK(rt::solve<rt::Newton>(as_array, 0.0) == rt::solve<rt::Newton>(wal, 0.0));

    // A tuple fn on Brent and Bisection: its first element, so the value-only fn's result.
    CHECK(rt::solve<rt::Brent>(wal, rt::Bracket{ 2.0, 3.0 }) == rt::solve<rt::Brent>(f, rt::Bracket{ 2.0, 3.0 }));
    CHECK(
        rt::solve<rt::Bisection>(wal, rt::Bracket{ 2.0, 3.0 }) == rt::solve<rt::Bisection>(f, rt::Bracket{ 2.0, 3.0 }));

    // Each separate callable is evaluated once per point.
    int nf = 0;
    int nd = 0;
    auto const cf = [&nf](double x) {
        ++nf;
        return x * x - 2.0;
    };
    auto const cd = [&nd](double x) {
        ++nd;
        return 2.0 * x;
    };
    REQUIRE(rt::solve<rt::Newton>(cf, cd, 1.0).has_value());
    CHECK(nf == 7);
    CHECK(nd == 7);
}

TEST_CASE("roots: open methods, failures")
{
    // Newton's cycle between 0 and 1 never changes sign: the cap, 4 x 53 iterations after the guess's evaluation.
    auto const [c, nc] = counted<rt::Newton>(cycle, 0.0);
    CHECK(c == std::unexpected(rt::error::not_converged));
    CHECK(nc == 1 + 212);
    // Halley follows Newton there (its correction points against Newton's step), Secant escapes.
    CHECK(rt::solve<rt::Halley>(cycle, 0.0) == std::unexpected(rt::error::not_converged));
    CHECK(brackets(rt::solve<rt::Secant>(cycle, 0.0), cycle_root, 1.8e-16));

    // Divergence: 1/x - 1 from 3 moves away from the root.
    CHECK(counted<rt::Newton>(reciprocal, 3.0)
          == std::pair{ std::expected<rt::Root<double>, rt::error>{ std::unexpected(rt::error::not_converged) }, 10 });
    CHECK(counted<rt::Secant>(reciprocal, 3.0).first == std::unexpected(rt::error::not_converged));

    // A multiple root at 0 from a guess with the default abs_tol: no sign change before the cap; a bracket solves it.
    CHECK(rt::solve<rt::Newton>(cubed, 1.0) == std::unexpected(rt::error::not_converged));
    CHECK(rt::solve<rt::Secant>(cubed, 1.0) == std::unexpected(rt::error::not_converged));
    CHECK(rt::solve<rt::Newton>(cubed, rt::Bracket{ -1.0, 2.0 }) == rt::Root<double>{ 0.0, dmin });

    // No usable step: f' = 0 at the guess; Secant's two values equal.
    CHECK(rt::solve<rt::Newton>(sq, 0.0) == std::unexpected(rt::error::not_converged));
    auto const flat = [](double x) { return x < 5.0 ? -1.0 : 1.0; };
    CHECK(rt::solve<rt::Secant>(flat, 0.0) == std::unexpected(rt::error::not_converged));

    // A guess or a value of fn that is not finite.
    CHECK(rt::solve<rt::Newton>(sq, dnan) == std::unexpected(rt::error::non_finite));
    CHECK(rt::solve<rt::Secant>(sq, dinf) == std::unexpected(rt::error::non_finite));
    auto const ln = [](double x) { return std::tuple{ std::log(x / 1e-12), 1.0 / x, -1.0 / (x * x) }; };
    CHECK(rt::solve<rt::Newton>(ln, 1e-11) == std::unexpected(rt::error::non_finite)); // the first step is below 0
    CHECK(brackets(rt::solve<rt::Newton>(ln, 1e-13), 1e-12, 2.5e-28)); // from below it converges

    CHECK(std::format("{}", rt::error::not_converged) == "not_converged");
}

TEST_CASE("roots: safeguarded methods, edges")
{
    // A derivative that is not finite makes the step unusable: bisection steps, not an error, at an end and at 0.
    auto const sqrt_fn = [](double x) { return std::pair{ std::sqrt(x) - 0.5, 0.5 / std::sqrt(x) }; }; // f'(0) = inf
    auto const cbrt_fn = [](double x) { return std::pair{ std::cbrt(x) - 0.5, 1.0 / (3.0 * std::cbrt(x * x)) }; };
    // sqrt x - 0.5 changes sign at 0.25, where sqrt is exact and correctly rounded (IEEE): c = 0. cbrt is a libm
    // function, not correctly rounded: glibc's returns 0.5 + 2^-53 one ulp below 0.125, so the computed fn changes
    // sign about 1.1e-16 from 0.125. Its 2-ulp budget, 2.2e-16 over f' = 1 / (3 cbrt(x)^2) = 1.33: c = 1.7e-16.
    CHECK(brackets(rt::solve<rt::Newton>(sqrt_fn, rt::Bracket{ 0.0, 1.0 }), 0.25, 0.0));
    CHECK(brackets(
        rt::solve<rt::Newton>(cbrt_fn, rt::Bracket{ -1.0, 1.0 }), 0.125, 1.7e-16)); // f'(0) = inf at zero first

    // A pole inside a bracket is returned as a root, as for Brent (documented): 1/(x - 1/3) on [0, 1].
    double const third = 1.0 / 3.0;
    auto const pole = [third](double x) {
        double const u = x - third;
        return std::pair{ 1.0 / u, -1.0 / (u * u) };
    };
    auto const p = rt::solve<rt::Newton>(pole, rt::Bracket{ 0.0, 1.0 });
    CHECK((p && nxx::math::abs(p->value - third) <= p->abs_error));

    // A hand-over across 0: atan x from 1, Newton crosses 0 and the bracket loop finds atan's exact zero there.
    auto const atan_fn = [](double x) { return std::pair{ std::atan(x), 1.0 / (1.0 + x * x) }; };
    CHECK(rt::solve<rt::Newton>(atan_fn, 1.0) == rt::Root<double>{ 0.0, dmin });
}

TEST_CASE("roots: open methods in float, long double and multiprecision types")
{
    // x^2 - 2 and the Wallis cubic with c = eps |ref| (fn's rounding, 0.35 and 0.41 eps relative, as for double).
    auto const check = [](auto zero) {
        using T = decltype(zero);
        T const eps = std::numeric_limits<T>::epsilon();
        T const tiny = T(std::numeric_limits<T>::denorm_min());
        T const s2 = from_string<T>(sqrt2_70);
        T const w = from_string<T>(wallis_70);
        auto const sq_t = [](T const &x) { return std::tuple<T, T, T>{ T(x * x - 2), T(2 * x), T(2) }; };
        auto const wal_t = [](T const &x) {
            return std::tuple<T, T, T>{ T(T(x * x - 2) * x - 5), T(T(3 * x) * x - 2), T(6 * x) };
        };
        for (auto const &r : { rt::solve<rt::Newton>(sq_t, T(1)),
                 rt::solve<rt::Halley>(sq_t, T(1)),
                 rt::solve<rt::Secant>(sq_t, T(1)),
                 rt::solve<rt::Newton>(sq_t, rt::Bracket<T>{ T(1), T(2) }) })
            CHECK(brackets(r, s2, T(eps * s2), T(2 * T(tiny + T(2 * eps * s2)))));
        for (auto const &r : { rt::solve<rt::Newton>(wal_t, T(0)),
                 rt::solve<rt::Halley>(wal_t, T(10)),
                 rt::solve<rt::Secant>(wal_t, T(2)),
                 rt::solve<rt::Halley>(wal_t, rt::Bracket<T>{ T(2), T(3) }) })
            CHECK(brackets(r, w, T(eps * w), T(2 * T(tiny + T(2 * eps * w)))));
    };
    check(float{});
    check(static_cast<long double>(0));
    check(bin50{});
    check(bin50_et{});
    check(dec50{});
}

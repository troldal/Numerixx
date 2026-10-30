// Romberg of nxx::quadrature (Phase 1, docs/quadrature/PHASE1_API.md E), and the module's result types.
//
// Accuracy: |value - ref| within the bound of Bauer, Rutishauser and Stiefel (1963), as Stoer and Bulirsch state it,
//   |I - R(k, k)| <= (b - a)^(2k+3) 2^(-k(k+1)) |B_(2k+2)| / (2k + 2)! max |f^(2k+2)|,
// at the level k at which the case stops (from its evaluations, 2^k + 1), plus the rounding: about 14 eps L1 from
// Romberg's own arithmetic (PHASE1_API.md C: the extrapolation's weights are positive with sum below 1.97; the
// allowance covers radix 10, where w / 2 and 2^-k round too) and fn's own rounding, rho L1, through the same weights,
// 1.97 rho L1. Where that bound is useless (cos 1000 x, a narrow peak) and where Richardson's assumptions fail (a kink,
// sqrt x), only the honesty check. The honesty check is labelled as such: |value - ref| <= abs_error, the estimate
// covering the error; it would pass a defect that lost three digits (the estimate is 30 to 10^4 times the error), so it
// never stands in for the bound. Evaluation counts are asserted exactly in double for fn of + - * / and sqrt
// (correctly rounded), against an independent Python implementation of the specification that mirrors the header's
// operation order; its counts are recorded in docs/quadrature/design.md, the prototype itself is not kept. No count is
// asserted where a libm function's last bits could move a stop.

#include "quadrature_check.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <format>
#include <functional>
#include <limits>
#include <stdexcept>
#include <utility>

using namespace quadrature_check;

namespace {

    // |B_2m| / (2m)! for 2m = 2, 4, ..., 34 (mpmath, 20 digits): the Bernoulli numbers of the bound, up to the cap's
    // level 16.
    constexpr double bernoulli_over_factorial[] = {
        8.3333333333333333333e-2,  1.3888888888888888889e-3,  3.3068783068783068783e-5,  8.2671957671957671958e-7,
        2.0876756987868098979e-8,  5.2841901386874931848e-10, 1.3382536530684678833e-11, 3.3896802963225828668e-13,
        8.5860620562778445641e-15, 2.174868698558061873e-16,  5.5090028283602295152e-18, 1.3954464685812523341e-19,
        3.5347070396294674717e-21, 8.9535174270375468504e-23, 2.2679524523376830603e-24, 5.7447906688722024453e-26,
        1.4551724756148649019e-27,
    };

    // The level k at which a run stopped, from its evaluations 2^k + 1.
    int level_of(int evaluations)
    {
        int k = 0;
        while ((1 << k) + 1 < evaluations) ++k;
        REQUIRE((1 << k) + 1 == evaluations);
        return k;
    }

    // The bound of Bauer, Rutishauser and Stiefel at level k on an interval of the given width, for
    // max |f^(2k+2)| = exp(log_max_derivative); 0 for a polynomial of degree below 2k + 2 (log = -inf).
    double truncation(int k, double width, double log_max_derivative)
    {
        double const m = static_cast<double>(k);
        double const bound = std::exp((2.0 * m + 3.0) * std::log(width) - m * (m + 1.0) * std::log(2.0)
                                      + std::log(bernoulli_over_factorial[k]) + log_max_derivative);
        REQUIRE(std::isfinite(bound)); // a bound that overflowed would pass every check
        return bound;
    }

    // log max |f^(2k+2)| for 1 / (1 + x^2): at most log (2k + 2)! (f = Im 1 / (x - i), |x - i| >= 1).
    double runge_derivative(int k)
    { return std::lgamma(2.0 * k + 3.0); }

    // Romberg's rounding: 14 eps L1 of its own, and fn's rho L1 through weights whose absolute values sum below 1.97.
    double rounding(double eps, double rho, double l1)
    {
        double const bound = (14.0 * eps + 1.97 * rho) * l1;
        REQUIRE(std::isfinite(bound));
        return bound;
    }

    // A run with its evaluations counted.
    template<typename T, typename F>
    std::pair<std::expected<qd::Integral<T>, qd::Failure<T>>, int> counted(F fn, T a, T b, qd::Criteria<T> c = {})
    {
        int n = 0;
        auto const r = qd::integrate<qd::Romberg>(Counted<F>{ fn, &n }, a, b, c);
        return { r, n };
    }

    constexpr double pi = 3.141592653589793;
    constexpr double no_derivative = -std::numeric_limits<double>::infinity(); // the truncation of a low polynomial
    constexpr auto runge = [](double x) { return 1.0 / (1.0 + x * x); };
    constexpr auto kink = [](double x) { return x < 0.2 ? 0.2 - x : x - 0.2; };
    // (x (1 - x) prod_(j=1..15) (16 x - j))^2: 0 at every multiple of 1/16, so at every sample through level 4.
    constexpr auto zero_samples = [](double x) {
        double p = x * (1.0 - x);
        for (int j = 1; j <= 15; ++j) p *= 16.0 * x - j;
        return p * p;
    };
    // Values near max of alternating signs, 0 elsewhere: every value and L1 finite, but the extrapolation overflows,
    // in the estimate at level 1 and in R(2, 2) at level 2 (the prototype's cases).
    constexpr auto estimate_overflow = [](double x) {
        if (x == 0.0 || x == 1.0) return -0.9 * dmax;
        return x == 0.5 ? 0.9 * dmax : 0.0;
    };
    constexpr auto diagonal_overflow = [](double x) {
        if (x == 0.0 || x == 1.0) return -0.3 * dmax;
        if (x == 0.5) return -0.7 * dmax;
        return (x == 0.25 || x == 0.75) ? dmax : 0.0;
    };

    // The floor of x^2 on [0, 3] at level 4: 50 eps L1, L1 the trapezoid sum of |fn|, 9 + 9 / 512 (the trapezoid
    // rule's error h^2 / 2 at h = 3 / 16). Every sample and sum is exact in double (dyadic points, values of a few
    // bits), so the floor is that double exactly; so is x^3's on [-1, 1], L1 = 2 (1/4 + h^2 / 4) = 0.5078125 at h = 1/8
    // (the trapezoid rule's error h^2 (f'(1) - f'(0)) / 12 on each half).
    constexpr double square_floor = 50.0 * deps * (9.0 + 9.0 / 512);
    constexpr double cube_floor = 50.0 * deps * 0.5078125;

} // namespace

// x^2 on [0, 3] in a constant expression: R(k, k) is exact for it from k = 1 on, so the value is within its rounding,
// and the estimate is the floor (run-time twins in "quadrature: romberg: accuracy in double"). .value(), whose throw on
// an empty result is a compile error (Q6).
static_assert(nxx::math::abs(qd::integrate<qd::Romberg>(square, 0.0, 3.0).value().value - 9.0) <= 14.0 * deps * 9.0);
static_assert(qd::integrate<qd::Romberg>(square, 0.0, 3.0).value().abs_error == square_floor);
static_assert(nxx::math::abs(qd::integrate<qd::Romberg>([](float x) { return x * x; }, 0.0f, 3.0f).value().value - 9.0f)
              <= 14.0f * std::numeric_limits<float>::epsilon() * 9.0f);
static_assert(
    nxx::math::abs(qd::integrate<qd::Romberg>([](long double x) { return x * x; }, 0.0L, 3.0L).value().value - 9.0L)
    <= 14.0L * std::numeric_limits<long double>::epsilon() * 9.0L);
// The edges and the failures are values in constant expressions too.
static_assert(qd::integrate<qd::Romberg>(square, 2.0, 2.0) == qd::Integral<double>{ 0.0, 0.0 });
static_assert(qd::integrate<qd::Romberg>(square, dnan, 1.0)
              == std::unexpected(qd::Failure<double>{ qd::error::non_finite, std::nullopt }));
static_assert(qd::integrate<qd::Romberg>(square, 0.0, 1.0, { .abs_tol = -1.0 })
              == std::unexpected(qd::Failure<double>{ qd::error::non_finite, std::nullopt }));

// The independent review's constant-evaluation cases (PHASE1_API.md, "Independent review"), run-time twins in
// "quadrature: romberg: constant evaluation of the edges".
namespace {
    constexpr auto zero_fn = [](double) { return 0.0; };
    constexpr auto nan_fn = [](double) { return dnan; };
    constexpr auto inf_fn = [](double x) { return x == 0.5 ? dinf : x; };
    constexpr auto max_fn = [](double) { return dmax; };
    constexpr auto large_fn = [](double) { return 0.8 * dmax; };
    constexpr auto wrapped_one = [](Wrapped const &) { return Wrapped(1); };
} // namespace
// A zero fn on [-max, max], where b - a overflows and is not formed: 0 with the floor 50 denorm_min, at level 4.
static_assert(qd::integrate<qd::Romberg>(zero_fn, -dmax, dmax) == qd::Integral<double>{ 0.0, 50.0 * dmin });
// A value of fn that is not finite, at an end and inside, and a sum that overflows: non_finite, as at run time.
static_assert(qd::integrate<qd::Romberg>(nan_fn, 0.0, 1.0)
              == std::unexpected(qd::Failure<double>{ qd::error::non_finite, std::nullopt }));
static_assert(qd::integrate<qd::Romberg>(inf_fn, 0.0, 1.0)
              == std::unexpected(qd::Failure<double>{ qd::error::non_finite, std::nullopt }));
static_assert(qd::integrate<qd::Romberg>(max_fn, 0.0, 3.0)
              == std::unexpected(qd::Failure<double>{ qd::error::non_finite, std::nullopt }));
// Conservative in constant evaluation: T_0 = 0.8 max lies above max / 2, so the extrapolation is refused there, where
// run time gives 0.8 max.
static_assert(qd::integrate<qd::Romberg>(large_fn, 0.0, 1.0)
              == std::unexpected(qd::Failure<double>{ qd::error::non_finite, std::nullopt }));
// A real type without a default constructor: 1 on [0, 1], every sample and sum exact (dyadic weights), so 1 exactly
// with the floor 50 eps.
static_assert(qd::integrate<qd::Romberg>(wrapped_one, Wrapped(0), Wrapped(1)).value().value == Wrapped(1));
static_assert(
    qd::integrate<qd::Romberg>(wrapped_one, Wrapped(0), Wrapped(1)).value().abs_error == Wrapped(50.0 * deps));

// The defaults, and the overload's constraint: only Romberg takes criteria, and both ends deduce T.
static_assert(qd::Criteria<double>{} == qd::Criteria<double>{ 0.0, nxx::math::root_epsilon<2, double>() });
static_assert(qd::Criteria<double>{}.rel_tol > 1.49e-8 && qd::Criteria<double>{}.rel_tol < 1.50e-8);
template<typename A, typename B>
concept romberg_integrable = requires(A a, B b) { qd::integrate<qd::Romberg>(square, a, b); };
static_assert(romberg_integrable<double, double> && !romberg_integrable<double, int>);

TEST_CASE("quadrature: romberg: evaluation counts in double")
{
    // The prototype's counts (docs/quadrature/design.md), for fn of + - * / and sqrt.
    struct Case
    {
        char const *name;
        std::pair<std::expected<qd::Integral<double>, qd::Failure<double>>, int> run;
        bool ok;
        int evaluations;
    };
    auto const sqrt_fn = [](double x) { return std::sqrt(x); };
    Case const cases[] = {
        { "1/(1+x^2) [0,1]", counted(runge, 0.0, 1.0), true, 33 },
        { "1/(1+x^2) [0,1] 1e-14", counted(runge, 0.0, 1.0, { .rel_tol = 1e-14 }), true, 257 },
        { "1/(1+x^2) [0,1] rel 0", counted(runge, 0.0, 1.0, { .rel_tol = 0.0 }), true, 257 },
        { "1/(1+x^2) [1,0]", counted(runge, 1.0, 0.0), true, 33 },
        { "1/(1+x^2) rel 1e300", counted(runge, 0.0, 1.0, { .rel_tol = 1e300 }), true, 17 },
        { "1/(1+x^2) abs 1e300", counted(runge, 0.0, 1.0, { .abs_tol = 1e300 }), true, 17 },
        { "x^2 [0,3]", counted(square, 0.0, 3.0), true, 17 },
        { "x^2 [0,3] 1e-14", counted(square, 0.0, 3.0, { .rel_tol = 1e-14 }), true, 17 },
        { "x^2 [3,0]", counted(square, 3.0, 0.0), true, 17 },
        { "x^3 [-1,1]", counted(cube, -1.0, 1.0), true, 17 },
        { "|x-1/5| [0,1]", counted(kink, 0.0, 1.0), true, 16385 },
        { "|x-1/5| [0,1] 1e-14", counted(kink, 0.0, 1.0, { .rel_tol = 1e-14 }), false, 65537 },
        { "sqrt x [0,1]", counted(sqrt_fn, 0.0, 1.0), true, 65537 },
        { "sqrt x [0,1] 1e-14", counted(sqrt_fn, 0.0, 1.0, { .rel_tol = 1e-14 }), false, 65537 },
        { "1e308 [0,1e-10]", counted([](double) { return 1e308; }, 0.0, 1e-10), true, 17 },
        { "1e-310 [0,1e10]", counted([](double) { return 1e-310; }, 0.0, 1e10), true, 17 },
        { "0 [-max,max]", counted([](double) { return 0.0; }, -dmax, dmax), true, 17 },
        { "1e-300 [-max,max]", counted([](double) { return 1e-300; }, -dmax, dmax), true, 17 },
        { "1e300 [0,3dmin]", counted([](double) { return 1e300; }, 0.0, 3.0 * dmin), true, 17 },
        { "zero samples [0,1]", counted(zero_samples, 0.0, 1.0), true, 17 },
        { "nan inside", counted([](double x) { return x == 0.5 ? dnan : x; }, 0.0, 1.0), false, 3 },
        { "inf at an end", counted([](double x) { return x == 1.0 ? dinf : x; }, 0.0, 1.0), false, 2 },
        { "estimate overflows", counted(estimate_overflow, 0.0, 1.0), false, 3 },
        { "denorm_min [-max,max]", counted([](double) { return dmin; }, -dmax, dmax), true, 17 },
        { "0 inside [D,3D], 1e300 outside",
          counted([](double x) { return x < dmin || x > 3.0 * dmin ? 1e300 : 0.0; }, dmin, 3.0 * dmin), true, 17 },
        { "R(2,2) overflows", counted(diagonal_overflow, 0.0, 1.0), false, 5 },
    };
    for (Case const &c : cases) {
        INFO(c.name);
        CHECK(c.run.first.has_value() == c.ok);
        CHECK(c.run.second == c.evaluations);
    }
}

TEST_CASE("quadrature: romberg: accuracy in double")
{
    // Each case: fn, the interval, the criteria, the reference (mpmath), log max |f^(2k+2)| on the interval as a
    // function of k, fn's rounding rho and L1. The truncation bound at the level reached plus the rounding; then,
    // separately, the honesty check. Where the truncation bound at the level reached is negligible (e^x at 1e-14,
    // k = 6: 8e-24; 1/(1+x^2) at 1e-14, k = 8: 1e-20; the polynomials), the check is the rounding allowance alone,
    // about 100 times the errors measured (7.7e-17 and 3.1e-17): a worst-case bound on rounding, not a sharp one.
    struct Case
    {
        char const *name;
        std::function<double(double)> fn;
        double a;
        double b;
        qd::Criteria<double> criteria;
        char const *ref;
        std::function<double(int)> log_derivative;
        double rho;
        double l1;
    };
    auto const exp_fn = [](double x) { return std::exp(x); };
    auto const e_derivative = [](int) { return 1.0; };                  // log e
    auto const flat = [](int) { return no_derivative; };                // a polynomial of degree below 2k + 2
    auto const sin2_derivative = [](int k) { return (2.0 * k + 1.0) * std::log(2.0); }; // (1 - cos 2x) / 2
    Case const cases[] = {
        // 1 / (1 + x^2) = pi / 4: rho = 2 eps (three roundings).
        { "1/(1+x^2)", runge, 0.0, 1.0, {}, pi_over_4, runge_derivative, 2.0 * deps, 0.7854 },
        { "1/(1+x^2) 1e-14", runge, 0.0, 1.0, { .rel_tol = 1e-14 }, pi_over_4, runge_derivative, 2.0 * deps, 0.7854 },
        // e^x = e - 1: std::exp within 2 ulps. At k = 4 the bound is 5.4e-14 (the prototype's error: 3.3e-14).
        { "e^x", exp_fn, 0.0, 1.0, {}, e_minus_1, e_derivative, 2.0 * deps, 1.7183 },
        { "e^x 1e-14", exp_fn, 0.0, 1.0, { .rel_tol = 1e-14 }, e_minus_1, e_derivative, 2.0 * deps, 1.7183 },
        // sin on [0, fl(pi)] = 1 - cos(fl(pi)) = 2 - 7.5e-33; std::sin within 2 ulps.
        { "sin", [](double x) { return std::sin(x); }, 0.0, pi, {}, "1.999999999999999999999999999999992501201",
          [](int) { return 0.0; }, 2.0 * deps, 2.0 },
        // sin^2 on [0, 2 fl(pi)] = pi to 1e-33; rho = 3 eps (a libm sin and a product). 1.x returned 0 here after 2
        // levels (the samples at 0, pi and 2 pi vanish); the minimum of 4 levels prevents it.
        { "sin^2",
          [](double x) {
              double const s = std::sin(x);
              return s * s;
          },
          0.0, 2.0 * pi, {}, "3.141592653589793238462643383279502884197", sin2_derivative, 3.0 * deps, pi },
        { "sin^2 1e-14",
          [](double x) {
              double const s = std::sin(x);
              return s * s;
          },
          0.0, 2.0 * pi, { .rel_tol = 1e-14 }, "3.141592653589793238462643383279502884197", sin2_derivative,
          3.0 * deps, pi },
        // Polynomials of degree below 2k + 2, integrated exactly by R(k, k): the rounding alone.
        { "x^2", square, 0.0, 3.0, {}, "9", flat, deps / 2, 9.0 },
        { "x^3", cube, -1.0, 1.0, { .rel_tol = 1e-14 }, "0", flat, deps, 0.5 },
        { "x^3 abs 0", cube, -1.0, 1.0, { .abs_tol = 0.0, .rel_tol = 0.0 }, "0", flat, deps, 0.5 },
        // Reversed ends: the negated integral, at the same level.
        { "e^x reversed", exp_fn, 1.0, 0.0, {},
          "-1.718281828459045235360287471352662497757247093699959574966967627724077", e_derivative, 2.0 * deps,
          1.7183 },
        // 1e-300 e^x: the bound scaled by 1e-300 (log e - 300 log 10); rho = 3 eps.
        { "1e-300 e^x", [](double x) { return 1e-300 * std::exp(x); }, 0.0, 1.0, {},
          "1.71828182845904527841886960948e-300", [](int) { return 1.0 - 300.0 * std::log(10.0); }, 3.0 * deps,
          1.72e-300 },
    };
    for (Case const &c : cases) {
        INFO(c.name);
        auto const [r, n] = counted(c.fn, c.a, c.b, c.criteria);
        REQUIRE(r.has_value());
        int const k = level_of(n);
        double const error = error_of(r->value, c.ref);
        CHECK(error <= truncation(k, std::abs(c.b - c.a), c.log_derivative(k)) + rounding(deps, c.rho, c.l1));
        CHECK(error <= r->abs_error); // the honesty check
    }

    // The floors, exact in double (see square_floor): an integral that is exactly 0 reports its floor, not 0; and the
    // run-time twin of the constant-evaluation case.
    auto const x2 = qd::integrate<qd::Romberg>(square, 0.0, 3.0);
    REQUIRE(x2.has_value());
    CHECK(x2->abs_error == square_floor);
    auto const x3 = qd::integrate<qd::Romberg>(cube, -1.0, 1.0, { .abs_tol = 0.0, .rel_tol = 1e-14 });
    REQUIRE(x3.has_value());
    CHECK(x3->abs_error == cube_floor);
}

TEST_CASE("quadrature: romberg: the estimate covers the error (honesty check, not the accuracy test)")
{
    // Where the bound is useless: cos 1000 x (whose 2k+2-th derivative is 1000^(2k+2)) and a narrow peak
    // e^(-1e4 (x - 1/2)^2); and where Richardson's assumptions fail: the kink at 1/5 and sqrt x. References at 40
    // digits (mpmath).
    struct Case
    {
        char const *name;
        std::expected<qd::Integral<double>, qd::Failure<double>> r;
        char const *ref;
    };
    Case const cases[] = {
        { "cos 1000x", qd::integrate<qd::Romberg>([](double x) { return std::cos(1000.0 * x); }, 0.0, 1.0),
          "0.0008268795405320025602558874291092181412127" },
        { "peak",
          qd::integrate<qd::Romberg>([](double x) { return std::exp(-1e4 * (x - 0.5) * (x - 0.5)); }, 0.0, 1.0),
          "0.01772453850905516027298167483341145182798" },
        { "kink", qd::integrate<qd::Romberg>(kink, 0.0, 1.0), "0.34" },
        { "sqrt", qd::integrate<qd::Romberg>([](double x) { return std::sqrt(x); }, 0.0, 1.0),
          "0.6666666666666666666666666666666666666667" },
    };
    for (Case const &c : cases) {
        INFO(c.name);
        REQUIRE(c.r.has_value());
        CHECK(error_of(c.r->value, c.ref) <= c.r->abs_error);
    }
}

TEST_CASE("quadrature: romberg: other types")
{
    // 1 / (1 + x^2) on [0, 1] with the default criteria, sqrt(eps) of each type, within the bound at the level the
    // run stopped at (fn's rho = 2 eps, L1 = pi / 4), and x^2 on [0, 3] within its rounding.
    auto const check = []<typename T>(T) {
        double const eps = eps_of<T>();
        auto const fn = [](T const &x) { return T(T(1) / T(T(1) + T(x * x))); };
        auto const [r, n] = counted(fn, T(0), T(1));
        REQUIRE(r.has_value());
        int const k = level_of(n);
        CHECK(error_of(r->value, pi_over_4)
              <= truncation(k, 1.0, runge_derivative(k)) + rounding(eps, 2.0 * eps, 0.7854));
        CHECK(error_of(r->value, pi_over_4) <= static_cast<double>(wide(r->abs_error))); // the honesty check
        auto const s = qd::integrate<qd::Romberg>([](T const &x) { return T(x * x); }, T(0), T(3));
        REQUIRE(s.has_value());
        CHECK(error_of(s->value, "9") <= rounding(eps, eps / 2, 9.0));
    };
    check(0.0f);
    check(0.0L);
    check(bin50{});
    check(bin50_et{});
    check(dec50{});

    // A 50-digit request in cpp_bin_float_50: rel_tol 1e-45 (the prototype at 50 digits needed 16,385 evaluations).
    auto const fn = [](bin50 const &x) { return bin50(1 / (1 + x * x)); };
    auto const [r, n] = counted(fn, bin50(0), bin50(1), { .rel_tol = bin50("1e-45") });
    REQUIRE(r.has_value());
    int const k = level_of(n);
    double const eps = eps_of<bin50>();
    CHECK(error_of(r->value, pi_over_4) <= truncation(k, 1.0, runge_derivative(k)) + rounding(eps, 2.0 * eps, 0.7854));
}

TEST_CASE("quadrature: romberg: constant evaluation of the edges")
{
    // The run-time twins of the independent review's static_asserts.
    qd::Failure<double> const invalid{ qd::error::non_finite, std::nullopt };
    CHECK(qd::integrate<qd::Romberg>(zero_fn, -dmax, dmax) == qd::Integral<double>{ 0.0, 50.0 * dmin });
    CHECK(qd::integrate<qd::Romberg>(nan_fn, 0.0, 1.0) == std::unexpected(invalid));
    CHECK(qd::integrate<qd::Romberg>(inf_fn, 0.0, 1.0) == std::unexpected(invalid));
    CHECK(qd::integrate<qd::Romberg>(max_fn, 0.0, 3.0) == std::unexpected(invalid));
    // At run time the conservative case gives its value: 0.8 max, every sample and sum exact (power-of-2 weights).
    auto const large = qd::integrate<qd::Romberg>(large_fn, 0.0, 1.0);
    REQUIRE(large.has_value());
    CHECK(error_of(large->value, bin100(0.8 * dmax)) <= 0.0);
    auto const w = qd::integrate<qd::Romberg>(wrapped_one, Wrapped(0), Wrapped(1));
    REQUIRE(w.has_value());
    CHECK(w->value == Wrapped(1));
    CHECK(w->abs_error == Wrapped(50.0 * deps));
}

TEST_CASE("quadrature: romberg: points stay inside the interval")
{
    // Subnormal widths, where c + hw u would leave the interval without the clamp (on [D, 3D], D = denorm_min,
    // c = hw = 2D), and [-max, max]; both orientations. Romberg's ends are fn(a) and fn(b) themselves.
    double const D = dmin;
    for (auto const &[a, b] : { std::pair{ D, 3.0 * D }, std::pair{ 3.0 * D, D }, std::pair{ 0.0, 3.0 * D },
                                std::pair{ D, 2.0 * D }, std::pair{ -3.0 * D, D }, std::pair{ -dmax, dmax },
                                std::pair{ dmax, -dmax }, std::pair{ 0.0, dmax } }) {
        INFO("[" << a << ", " << b << "]");
        int outside = 0;
        auto const fn = Confined{ [](double x) { return x / 2; }, std::min(a, b), std::max(a, b), &outside };
        (void)qd::integrate<qd::Romberg>(fn, a, b);
        (void)qd::integrate<qd::Romberg>(fn, a, b, { .rel_tol = 0.0 });
        CHECK(outside == 0);
    }
    // The review's counterexample: 0 on [D, 3D], 1e300 outside it; the integral 0, with the floor 50 denorm_min.
    auto const r = qd::integrate<qd::Romberg>([D](double x) { return x < D || x > 3.0 * D ? 1e300 : 0.0; }, D, 3.0 * D);
    CHECK(r == qd::Integral<double>{ 0.0, 50.0 * dmin });
}

TEST_CASE("quadrature: romberg: edges")
{
    // a == b: {0, 0}, fn not called; NaN and infinite ends and invalid tolerances: non_finite, fn not called.
    int n = 0;
    auto const c = Counted{ square, &n };
    auto const same = qd::integrate<qd::Romberg>(c, 2.0, 2.0);
    REQUIRE(same.has_value());
    CHECK(same->value == 0.0); // defined as 0 (A3), not computed
    CHECK(same->abs_error == 0.0);
    qd::Failure<double> const invalid{ qd::error::non_finite, std::nullopt };
    for (auto const &[a, b] : { std::pair{ dnan, 1.0 }, std::pair{ 0.0, dnan }, std::pair{ -dinf, 1.0 },
                                std::pair{ 0.0, dinf } })
        CHECK(qd::integrate<qd::Romberg>(c, a, b) == std::unexpected(invalid));
    for (qd::Criteria<double> const criteria :
         { qd::Criteria<double>{ .abs_tol = -1e-12 }, qd::Criteria<double>{ .rel_tol = -1e-12 },
           qd::Criteria<double>{ .abs_tol = dnan }, qd::Criteria<double>{ .rel_tol = dnan },
           qd::Criteria<double>{ .abs_tol = dinf }, qd::Criteria<double>{ .rel_tol = dinf } })
        CHECK(qd::integrate<qd::Romberg>(c, 0.0, 1.0, criteria) == std::unexpected(invalid));
    CHECK(n == 0);

    // fn not finite at an end, and inside: non_finite with no best estimate.
    for (auto const &r : { qd::integrate<qd::Romberg>([](double x) { return 1.0 / std::sqrt(x); }, 0.0, 1.0),
                           qd::integrate<qd::Romberg>([](double x) { return x == 0.5 ? dnan : x; }, 0.0, 1.0),
                           qd::integrate<qd::Romberg>([](double x) { return x == 0.25 ? -dinf : x; }, 0.0, 1.0) })
        CHECK(r == std::unexpected(invalid));

    // Every value finite, but a sum or the extrapolation overflows: non_finite. 1e308 sin x on [0, 2 pi], whose sum of
    // |fn| overflows at level 2; values near max whose extrapolation overflows in the estimate or in R(2, 2), which the
    // check of the estimate catches (R(k-1, k-1) is finite, so a non-finite R(k, k) makes it non-finite).
    auto const sine =
        qd::integrate<qd::Romberg>([](double x) { return 1e308 * std::sin(x); }, 0.0, 6.283185307179586);
    CHECK(sine == std::unexpected(invalid));
    CHECK(qd::integrate<qd::Romberg>(estimate_overflow, 0.0, 1.0) == std::unexpected(invalid));
    CHECK(qd::integrate<qd::Romberg>(diagonal_overflow, 0.0, 1.0) == std::unexpected(invalid));

    // The order of each term's products (detail::Interval::term), on constant fn: the rounding alone (rho = 0).
    // fn = 1e308 on [0, 1e-10]: 1e298 (1e308 fl(1e-10)), no overflow, since each term is scaled before the sum.
    auto const huge = qd::integrate<qd::Romberg>([](double) { return 1e308; }, 0.0, 1e-10);
    REQUIRE(huge.has_value());
    CHECK(error_of(huge->value, "1.00000000000000004741126094494e298") <= rounding(deps, 0.0, 1e298));
    // fn = 1e-310 on [0, 1e10]: fl(1e-310) 1e10, the scale (b - a) 2^-k formed first, so nothing is subnormal.
    auto const subnormal_fn = qd::integrate<qd::Romberg>([](double) { return 1e-310; }, 0.0, 1e10);
    REQUIRE(subnormal_fn.has_value());
    CHECK(error_of(subnormal_fn->value, "9.99999999999996944932750289769e-301") <= rounding(deps, 0.0, 1e-300));
    // [-max, max], b - a overflows and the terms are hw (w fn): fn = 0 gives 0 with the floor's denorm_min term,
    // 50 denorm_min exactly; fn = 1e-300 gives 2 max fl(1e-300) (L1 = 3.6e8).
    auto const zero = qd::integrate<qd::Romberg>([](double) { return 0.0; }, -dmax, dmax);
    REQUIRE(zero.has_value());
    CHECK(error_of(zero->value, "0") <= 0.0);
    CHECK(zero->abs_error == 50.0 * dmin);
    auto const wide_tiny = qd::integrate<qd::Romberg>([](double) { return 1e-300; }, -dmax, dmax);
    REQUIRE(wide_tiny.has_value());
    CHECK(error_of(wide_tiny->value, "359538626.972463141629054847463") <= rounding(deps, 0.0, 3.6e8));
    // The constant denorm_min on [-max, max]: 2 max denorm_min = 1.78e-15, with (hw y) w formed for |y| <= 1, where
    // hw (w y) underflowed (then not_converged with 3.3e-20). Every sample and sum is exact here: the value exactly.
    auto const smallest = qd::integrate<qd::Romberg>([](double) { return dmin; }, -dmax, dmax);
    REQUIRE(smallest.has_value());
    CHECK(error_of(smallest->value, bin100(2) * bin100(dmax) * bin100(dmin)) <= 0.0);

    // [0, 3 denorm_min] with fn = 1e300: 1.48e-23, fn's value multiplied first, where (b - a) 2^-k would round
    // 1.5 denorm_min to 2 and 3/8 denorm_min to 0 by level 3.
    auto const subnormal = qd::integrate<qd::Romberg>([](double) { return 1e300; }, 0.0, 3.0 * dmin);
    REQUIRE(subnormal.has_value());
    CHECK(error_of(subnormal->value, "1.48219693752373971035210123429e-23") <= rounding(deps, 0.0, 1.5e-23));

    // 1e-320 e^x on [0, 1]: eps L1 underflows to 0, and so does rel_tol |value|; the floor's denorm_min term,
    // 50 denorm_min, ends it at the first level allowed, 17 evaluations, whatever libm's last bits.
    {
        auto const [r, evaluations] = counted([](double x) { return 1e-320 * std::exp(x); }, 0.0, 1.0);
        REQUIRE(r.has_value());
        CHECK(evaluations == 17);
        CHECK(r->abs_error == 50.0 * dmin);
        CHECK(error_of(r->value, "1.71826269914134988937977026597e-320") <= r->abs_error);
    }

    // The cap: sqrt x and the kink at 1/5 at rel_tol 1e-14 end not_converged after 65,537 evaluations, with the best
    // estimate present and within its own error estimate.
    for (auto const &[r, ref] :
         { std::pair{ qd::integrate<qd::Romberg>(kink, 0.0, 1.0, { .rel_tol = 1e-14 }), "0.34" },
           std::pair{ qd::integrate<qd::Romberg>([](double x) { return std::sqrt(x); }, 0.0, 1.0, { .rel_tol = 1e-14 }),
                      "0.6666666666666666666666666666666666666667" } }) {
        REQUIRE(!r.has_value());
        CHECK(r.error().reason == qd::error::not_converged);
        REQUIRE(r.error().best.has_value());
        CHECK(error_of(r.error().best->value, ref) <= r.error().best->abs_error);
    }

    // Samples that all vanish: 0, with the floor's denorm_min term as the estimate, although the integral is positive.
    auto const zeros = qd::integrate<qd::Romberg>(zero_samples, 0.0, 1.0);
    REQUIRE(zeros.has_value());
    CHECK(error_of(zeros->value, "0") <= 0.0);
    CHECK(zeros->abs_error == 50.0 * dmin);
}

// Emscripten's default build compiles throw but aborts on it rather than unwinding to a catch, so the test runs
// everywhere else.
#ifndef __EMSCRIPTEN__
TEST_CASE("quadrature: romberg: fn's exceptions pass through")
{
    auto const throwing = [](double x) -> double {
        if (x > 0.6) throw std::runtime_error("fn failed");
        return x;
    };
    CHECK_THROWS_AS((void)qd::integrate<qd::Romberg>(throwing, 0.0, 1.0), std::runtime_error);
}
#endif

// The output in the examples of the formatters' comments (quadrature.hpp).
TEST_CASE("quadrature: error and Failure format as their names")
{
    CHECK(std::format("{}", qd::error::non_finite) == "non_finite");
    CHECK(std::format("{}", qd::error::not_converged) == "not_converged");
    CHECK(std::format("[{:>15}]", qd::error::not_converged) == "[  not_converged]");
    auto const r = qd::integrate<qd::Romberg>(kink, 0.0, 1.0, { .rel_tol = 1e-14 });
    REQUIRE(!r.has_value());
    CHECK(std::format("no integral: {}", r.error()) == "no integral: not_converged");
    CHECK(std::format("{}", qd::Failure<bin50>{ qd::error::non_finite, std::nullopt }) == "non_finite");
    CHECK(std::format("{:<12}|", qd::Failure<double>{ qd::error::non_finite, std::nullopt }) == "non_finite  |");
}

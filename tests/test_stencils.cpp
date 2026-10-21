#include <numerixx/deriv.hpp>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>

namespace mp = boost::multiprecision;

using bin50 = mp::cpp_bin_float_50;
using bin50_et = mp::number<mp::cpp_bin_float<50>, mp::et_on>;
using dec50 = mp::cpp_dec_float_50;

using nxx::deriv::Backward;
using nxx::deriv::Central;
using nxx::deriv::derivative;
using nxx::deriv::Forward;

namespace {

    template<typename S>
    constexpr auto table = nxx::deriv::detail::make_table(S{});

    // The moment M_j = sum_k w_k n_k^j of a table, computed here, independently of the closed forms in the library.
    template<std::size_t N>
    constexpr long long moment(nxx::deriv::detail::Table<N> const &t, int j)
    {
        long long m = 0;
        for (std::size_t k = 0; k < N; ++k) {
            long long power = 1;
            for (int i = 0; i < j; ++i) power *= t.nodes[k];
            m += t.numerators[k] * power;
        }
        return m;
    }

    // By Taylor's theorem, a rule of order p has M_0 = 0, M_1 = denominator, M_j = 0 for 2 <= j <= p and
    // M_(p+1) != 0. Its nodes increase strictly.
    template<typename S>
    constexpr bool has_order(int p)
    {
        constexpr auto t = table<S>;
        if (t.order != p) return false;
        for (std::size_t k = 1; k < t.nodes.size(); ++k)
            if (!(t.nodes[k - 1] < t.nodes[k])) return false;
        if (moment(t, 0) != 0 || moment(t, 1) != t.denominator) return false;
        for (int j = 2; j <= p; ++j)
            if (moment(t, j) != 0) return false;
        return moment(t, p + 1) != 0;
    }

    template<typename S, std::size_t N>
    constexpr bool has_weights(std::array<int, N> const &nodes, std::array<int, N> const &numerators, int denominator)
    {
        constexpr auto t = table<S>;
        if constexpr (t.nodes.size() != N)
            return false;
        else
            return t.nodes == nodes && t.numerators == numerators && t.denominator == denominator;
    }

    template<int N>
    constexpr auto power = [](double x) {
        double y = 1;
        for (int i = 0; i < N; ++i) y *= x;
        return y;
    };

    constexpr double dmax = std::numeric_limits<double>::max();
    constexpr auto half = [](double x) { return x / 2; };
    constexpr auto infinite = [](double) { return std::numeric_limits<double>::infinity(); };
    constexpr auto not_a_number = [](double) { return std::numeric_limits<double>::quiet_NaN(); };
    constexpr auto non_finite = std::unexpected(nxx::deriv::error::non_finite);

    // |a - b| <= tol * |b|
    template<typename T>
    bool close(std::expected<T, nxx::deriv::error> const &a, T const &b, T const &tol)
    {
        using std::abs;
        return a.has_value() && abs(*a - b) <= tol * abs(b);
    }

} // namespace

// Every table of every allowed order satisfies the moment conditions of its order.
static_assert(has_order<Central<2>>(2) && has_order<Central<4>>(4) && has_order<Central<6>>(6));
static_assert(has_order<Central<8>>(8) && has_order<Central<10>>(10));
static_assert(has_order<Forward<1>>(1) && has_order<Forward<2>>(2) && has_order<Forward<3>>(3));
static_assert(has_order<Forward<4>>(4) && has_order<Forward<5>>(5) && has_order<Forward<6>>(6));
static_assert(has_order<Backward<1>>(1) && has_order<Backward<2>>(2) && has_order<Backward<3>>(3));
static_assert(has_order<Backward<4>>(4) && has_order<Backward<5>>(5) && has_order<Backward<6>>(6));

// Spot checks against the published tables (Fornberg 1988, Table 1).
static_assert(has_weights<Central<2>, 2>({ -1, 1 }, { -1, 1 }, 2));
static_assert(has_weights<Central<4>, 4>({ -2, -1, 1, 2 }, { 1, -8, 8, -1 }, 12));
static_assert(
    has_weights<Central<8>, 8>({ -4, -3, -2, -1, 1, 2, 3, 4 }, { 3, -32, 168, -672, 672, -168, 32, -3 }, 840));
static_assert(has_weights<Forward<1>, 2>({ 0, 1 }, { -1, 1 }, 1));
static_assert(has_weights<Forward<4>, 5>({ 0, 1, 2, 3, 4 }, { -25, 48, -36, 16, -3 }, 12));
static_assert(has_weights<Backward<2>, 3>({ -2, -1, 0 }, { 1, -4, 3 }, 2));

// Only the stencil families in their allowed orders are accepted, and the orders outside are rejected by the
// constraint, not by a hard error.
template<typename S>
concept usable = requires { derivative<S>(power<2>, 1.0); };
static_assert(usable<Central<2>> && usable<Central<10>> && usable<Forward<1>> && usable<Backward<6>>);
static_assert(!usable<Central<0>> && !usable<Central<3>> && !usable<Central<12>>);
static_assert(!usable<Forward<0>> && !usable<Forward<7>> && !usable<Backward<0>> && !usable<Backward<7>>);
static_assert(!usable<double>);

// A rule of order p is exact for polynomials of degree <= p. With x = 1 and h = 1/4 every value is a short binary
// fraction and the evaluator divides once, so the results are exact. Degree p + 1 leaves exactly the truncation term
// M_(p+1) / ((p+1)! den) h^p f^(p+1).
static_assert(derivative<Central<4>>(power<4>, 1.0, 0.25) == 4.0);
static_assert(derivative<Central<6>>(power<6>, 1.0, 0.25) == 6.0);
static_assert(derivative<Central<10>>(power<10>, 1.0, 0.25) == 10.0);
static_assert(derivative<Forward<4>>(power<4>, 1.0, 0.25) == 4.0);
static_assert(derivative<Backward<6>>(power<6>, 1.0, 0.25) == 6.0);
static_assert(derivative<Central<4>>(power<5>, 1.0, 0.25) == 5.0 - 4.0 / 256); // -h^4/30 * 120
static_assert(derivative<Central<6>>(power<7>, 1.0, 0.25) == 7.0 + 36.0 / 4096); // +h^6/140 * 5040
static_assert(derivative<Forward<4>>(power<5>, 1.0, 0.25) == 5.0 - 24.0 / 256); // -h^4/5 * 120
static_assert(derivative<Backward<4>>(power<5>, 1.0, 0.25) == 5.0 - 24.0 / 256); // the same, mirrored
static_assert(derivative<Forward<2>>(power<3>, 1.0, 0.25) == 3.0 - 2.0 / 16); // -h^2/3 * 6
static_assert(derivative<Forward<1>>(power<2>, 1.0, 0.25) == 2.0 + 0.25); // +h/2 * 2
static_assert(derivative<Backward<1>>(power<2>, 1.0, 0.25) == 2.0 - 0.25); // -h/2 * 2

// Failures are error values in constant expressions too.
static_assert(derivative<Central<4>>(infinite, 1.0, 0.5) == non_finite);
static_assert(derivative<Forward<2>>(not_a_number, 1.0, 0.5) == non_finite);
static_assert(derivative<Central<6>>(power<3>, 1.0, 1e-300) == non_finite); // h is below the resolution of x
static_assert(derivative<Forward<1>>(power<3>, 1.0, 1e-300) == non_finite);
// Central<4> has Q = max(2 max|node|, den) = 12: a step above max / 24 is rejected before any point is formed.
static_assert(derivative<Central<4>>(half, 0.0, dmax / 20) == non_finite);
// A large step that passes: the points are exact multiples of max / 32, and only the sum and the divisor 12 h round,
// so the result is within a few ulp of 0.5.
static_assert(nxx::math::abs(*derivative<Central<4>>(half, 0.0, dmax / 32) - 0.5) < 1e-15);

// Regression: with Q = max(hi - lo, den), Forward<1> had Q = 1 and no margin. At x one ulp below max with
// h = 1.5 ulp(max), max - Q h ties down to x, so the check passed, and x + h tied up to infinity: fn was called at
// infinity and a bounded fn gave a successful 0. With Q = max(2 max|node|, den) = 2 the step is rejected.
namespace {
    constexpr auto saturating = [](double x) { return x < dmax ? 1.0 : 2.0; };
} // namespace
static_assert(derivative<Forward<1>>(saturating, dmax - 0x1p+971, 0x1.8p+971) == non_finite);
static_assert(derivative<Backward<1>>(saturating, -(dmax - 0x1p+971), 0x1.8p+971) == non_finite);

// Regression: at x = -1 the spacing is 2^-53 above x and 2^-52 below it, so h = 0x1.8p-54 moves x up but not down.
// Checking only x + h == x let Backward<2> evaluate coinciding points and return a successful 0 (true value 1).
static_assert(derivative<Backward<2>>(power<1>, -1.0, 0x1.8p-54) == non_finite);
static_assert(derivative<Central<4>>(power<1>, -1.0, 0x1.8p-54) == non_finite);

// Regression: a non-finite x with the default step gave a NaN or infinite step, and a NaN step tripped the h > 0
// precondition (an abort in debug builds, a compile error in constant expressions) instead of returning non_finite.
namespace {
    constexpr double qnan = std::numeric_limits<double>::quiet_NaN();
    constexpr double inf = std::numeric_limits<double>::infinity();
} // namespace
static_assert(derivative(power<1>, qnan) == non_finite);
static_assert(derivative<Central<6>>(power<1>, qnan) == non_finite);
static_assert(derivative<Forward<2>>(power<1>, qnan) == non_finite);
static_assert(derivative<Backward<1>>(power<1>, qnan) == non_finite);
static_assert(derivative(power<1>, inf) == non_finite && derivative<Forward<2>>(power<1>, -inf) == non_finite);

namespace {

    // The default step on f = x^(p+1) at x = 1, where the truncation error is exactly M_(p+1) / den * h^p
    // (f^(p+1) = (p+1)!). Rounding adds about S eps (1 + (p+1)/2) / h: fn's values are near 1, and rounding the
    // points contributes S/2 eps |x f'| / h with |x f'| = p + 1. Twice that is the tolerance.
    template<typename S>
    constexpr bool default_step_matches()
    {
        constexpr auto t = table<S>;
        constexpr int p = t.order;
        double const eps = std::numeric_limits<double>::epsilon();
        double const h = nxx::math::root_epsilon<p + 1, double>();
        double hp = 1;
        for (int i = 0; i < p; ++i) hp *= h;
        double const truncation = double(moment(t, p + 1)) / t.denominator * hp;
        double weights = 0;
        for (int const w : t.numerators) weights += w < 0 ? -w : w;
        double const rounding = weights / t.denominator * eps * (1 + (p + 1) / 2.0) / h;
        auto const d = derivative<S>(power<p + 1>, 1.0);
        return d.has_value() && nxx::math::abs(*d - (p + 1 + truncation)) <= 2 * rounding;
    }

} // namespace

static_assert(default_step_matches<Central<2>>() && default_step_matches<Central<4>>());
static_assert(default_step_matches<Central<6>>() && default_step_matches<Central<10>>());
static_assert(default_step_matches<Forward<1>>() && default_step_matches<Forward<2>>());
static_assert(default_step_matches<Forward<4>>() && default_step_matches<Forward<6>>());
static_assert(default_step_matches<Backward<2>>() && default_step_matches<Backward<6>>());

TEST_CASE("deriv: stencils with the default step")
{
    // h = eps^(1/(p+1)) for |x| <= 1. For exp, every derivative equals f, so the relative error is about
    // C eps^(p/(p+1)) from truncation plus S eps^(p/(p+1)) from rounding fn's values, and up to S/2 more from rounding
    // the points: (C + 2S) eps^(p/(p+1)), with C the error constant and S the sum of the absolute weights. This is an
    // estimate, not a bound (the products and partial sums round too, about 0.75 S more at worst); the measured errors
    // lie 4 to 16 times below it.
    double const eps = std::numeric_limits<double>::epsilon();
    auto const check = [eps]<typename S>(S, double c, double s) {
        int const p = table<S>.order;
        double const tol = (c + 2 * s) * std::pow(eps, double(p) / (p + 1));
        auto const exp = [](double x) { return std::exp(x); };
        for (double const x : { -1.0, 0.5, 1.0 }) CHECK(close(derivative<S>(exp, x), std::exp(x), tol));
    };
    check(Central<2>{}, 1.0 / 6, 1.0);
    check(Central<4>{}, 1.0 / 30, 1.5);
    check(Central<6>{}, 1.0 / 140, 11.0 / 6);
    check(Central<10>{}, 1.0 / 2772, 2.28);
    check(Forward<1>{}, 1.0 / 2, 2.0);
    check(Forward<2>{}, 1.0 / 3, 4.0);
    check(Forward<4>{}, 1.0 / 5, 32.0 / 3);
    check(Forward<6>{}, 1.0 / 7, 27.7);
    check(Backward<2>{}, 1.0 / 3, 4.0);
    check(Backward<6>{}, 1.0 / 7, 27.7);

    // The compile-time cases above, at run time.
    CHECK(default_step_matches<Central<2>>());
    CHECK(default_step_matches<Central<4>>());
    CHECK(default_step_matches<Central<6>>());
    CHECK(default_step_matches<Central<10>>());
    CHECK(default_step_matches<Forward<1>>());
    CHECK(default_step_matches<Forward<2>>());
    CHECK(default_step_matches<Forward<4>>());
    CHECK(default_step_matches<Forward<6>>());
    CHECK(default_step_matches<Backward<2>>());
    CHECK(default_step_matches<Backward<6>>());
}

TEST_CASE("deriv: one-sided stencils never evaluate across x")
{
    // fn is NaN on one side of 1: a central rule fails there, the one-sided rule on the defined side does not.
    auto const right = [](double x) { return x < 1.0 ? std::numeric_limits<double>::quiet_NaN() : x * x * x; };
    auto const left = [](double x) { return x > 1.0 ? std::numeric_limits<double>::quiet_NaN() : x * x * x; };
    CHECK(derivative<Central<4>>(right, 1.0) == non_finite);
    // Forward<2> on x^3 with h = eps^(1/3) = 6e-6: truncation -2 h^2 = 7e-11, rounding of fn about S eps / h = 1.5e-10
    // and of the points about S/2 eps |x f'| / h = 2.2e-10 (S = 4), so 1.5e-10 relative to 3.
    CHECK(close(derivative<Forward<2>>(right, 1.0), 3.0, 5e-10));
    CHECK(close(derivative<Backward<2>>(left, 1.0), 3.0, 5e-10));
    // Order 4 is exact for x^3, so only rounding remains: with S = 32/3 and h = eps^(1/5) = 7.4e-4, about
    // S eps / h + S/2 eps 3 / h = 8e-12, so 2.7e-12 relative.
    CHECK(close(derivative<Forward<4>>(right, 1.0), 3.0, 1e-11));
    CHECK(close(derivative<Backward<4>>(left, 1.0), 3.0, 1e-11));
}

TEST_CASE("deriv: stencils call fn once per node")
{
    // For fn(x) = x the values are the rounded points, so only point rounding remains: about S/2 eps / h, 4e-14 for
    // Central<6> (h = eps^(1/7)) and 1.6e-12 for Forward<4> (h = eps^(1/5)).
    int calls = 0;
    auto const fn = [&calls](double x) {
        ++calls;
        return x;
    };
    CHECK(close(derivative<Central<6>>(fn, 1.0), 1.0, 2e-13));
    CHECK(calls == 6);
    calls = 0;
    CHECK(close(derivative<Forward<4>>(fn, 1.0), 1.0, 1e-11));
    CHECK(calls == 5);
    calls = 0;
    CHECK(derivative<Backward<1>>(fn, 1.0, 0.5) == 1.0);
    CHECK(calls == 2);
}

TEST_CASE("deriv: stencil failures")
{
    int calls = 0;
    auto const counted = [&calls](double x) {
        ++calls;
        return x;
    };
    // Forward<4> has Q = 12: x = 0.99 max with h = max / 100 fails |x| <= max - 12 h, so x + 4h (which would
    // overflow) is never formed and fn is not called.
    CHECK(derivative<Forward<4>>(counted, 0.99 * dmax, dmax / 100) == non_finite);
    CHECK(derivative<Central<4>>(counted, 1.0, 1e-300) == non_finite);
    CHECK(calls == 0);

    // The weighted sum overflows although every value is finite.
    auto const huge = [](double x) { return x < 0.0 ? -dmax : dmax; };
    CHECK(derivative<Central<4>>(huge, 0.0, 0.5) == non_finite);

    // The regressions above, at run time: no point reaches infinity and fn is not called.
    auto const saturating_counted = [&calls](double x) {
        ++calls;
        return x < dmax ? 1.0 : 2.0;
    };
    calls = 0;
    CHECK(derivative<Forward<1>>(saturating_counted, dmax - 0x1p+971, 0x1.8p+971) == non_finite);
    CHECK(derivative<Backward<2>>(counted, -1.0, 0x1.8p-54) == non_finite);
    CHECK(derivative<Central<4>>(counted, -1.0, 0x1.8p-54) == non_finite);
    CHECK(calls == 0);

    // A non-finite x with the default step: no step is computed and fn is not called.
    CHECK(derivative(counted, qnan) == non_finite);
    CHECK(derivative<Central<6>>(counted, qnan) == non_finite);
    CHECK(derivative<Forward<2>>(counted, qnan) == non_finite);
    CHECK(derivative<Backward<1>>(counted, qnan) == non_finite);
    CHECK(derivative(counted, inf) == non_finite);
    CHECK(derivative<Forward<2>>(counted, -inf) == non_finite);
    CHECK(calls == 0);
}

TEST_CASE("deriv: stencils on multiprecision types")
{
    // Central<6>: (C + 2S) eps^(6/7) with C = 1/140 and S = 11/6, as for double. Forward<4>: (1/5 + 64/3) eps^(4/5).
    // cpp_dec_float carries guard digits beyond its epsilon, so its rounding is negligible and only truncation is
    // left: 2C eps^(p/(p+1)), with the factor 2 for the spread of exp over the stencil and the estimate's slack.
    auto const check = []<typename T>(T, bool guard_digits) {
        T const eps = std::numeric_limits<T>::epsilon();
        auto const exp = [](T const &x) { return T(mp::exp(x)); };
        T const tol6 = (guard_digits ? T(2) / T(140) : T(3.68)) * T(mp::pow(eps, T(6) / T(7)));
        T const tol4 = (guard_digits ? T(2) / T(5) : T(21.6)) * T(mp::pow(eps, T(4) / T(5)));
        CHECK(close(derivative<Central<6>>(exp, T(1)), T(mp::exp(T(1))), tol6));
        CHECK(close(derivative<Forward<4>>(exp, T(1)), T(mp::exp(T(1))), tol4));
        CHECK(close(derivative<Backward<4>>(exp, T(1)), T(mp::exp(T(1))), tol4));
    };
    check(bin50{}, false);
    check(bin50_et{}, false);
    check(dec50{}, true);
}

#include <numerixx/deriv.hpp>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <doctest/doctest.h>

#include <cmath>
#include <limits>

namespace mp = boost::multiprecision;

using bin50 = mp::cpp_bin_float_50;
using bin50_et = mp::number<mp::cpp_bin_float<50>, mp::et_on>;
using dec50 = mp::cpp_dec_float_50;

using nxx::deriv::derivative;

namespace {

    constexpr auto square = [](double x) { return x * x; };
    constexpr auto cube = [](double x) { return x * x * x; };

    // |a - b| <= tol * |b|
    template<typename T>
    bool close(std::expected<nxx::deriv::Estimate<T>, nxx::deriv::error> const &a, T const &b, T const &tol)
    {
        using std::abs;
        return a.has_value() && abs(T(a->value - b)) <= tol * abs(b);
    }

} // namespace

// The central difference is exact for quadratics, and has error h^2 for x^3: with steps that are powers of two the
// results are exact in binary floating point.
static_assert(derivative(square, 3.0, 0.5) == 6.0);
static_assert(derivative(cube, 2.0, 0.25) == 12.0 + 0.0625);
static_assert(nxx::math::abs(*derivative(cube, 2.0) - 12.0) < 1e-8);

// Failures are error values in constant expressions too: they are detected before the arithmetic that would produce
// an infinity or a NaN.
namespace {

    constexpr double dmax = std::numeric_limits<double>::max();
    constexpr auto half = [](double x) { return x / 2; };
    constexpr auto infinite = [](double) { return std::numeric_limits<double>::infinity(); };
    constexpr auto not_a_number = [](double) { return std::numeric_limits<double>::quiet_NaN(); };
    constexpr auto non_finite = std::unexpected(nxx::deriv::error::non_finite);

} // namespace

static_assert(derivative(infinite, 1.0, 0.5) == non_finite); // inf - inf is never formed
static_assert(derivative(not_a_number, 1.0, 0.5) == non_finite);
static_assert(derivative(cube, 1.0, 1e-300) == non_finite); // x + h == x - h == x: no division by zero
static_assert(derivative(half, 0.0, dmax) == non_finite); // the distance between x +- h would overflow
static_assert(derivative(half, dmax, 1.0) == non_finite); // x + h == x - h == max
static_assert(derivative(half, -dmax * 0.9, dmax / 5) == non_finite); // x - h would overflow
static_assert(derivative(half, 0.0, dmax / 8) == 0.5); // large but representable steps still work

// Integer arguments are not real numbers.
template<typename X>
concept differentiable_at = requires(X x) { derivative(square, x); };
static_assert(differentiable_at<double> && !differentiable_at<int>);

TEST_CASE("deriv: central difference with the default step")
{
    // h = eps^(1/3) * max(|x|, 1): truncation (h^2 / 6 * f''') and rounding (eps * f / h) errors are both about
    // eps^(2/3), 4e-11 for double, when f varies on the scale of max(|x|, 1).
    auto const exp = [](double x) { return std::exp(x); };
    auto const sin = [](double x) { return std::sin(x); };
    CHECK(close(derivative(exp, 1.0), std::exp(1.0), 1e-9));
    CHECK(close(derivative(sin, 0.5), std::cos(0.5), 1e-9));
    CHECK(close(derivative(cube, 1e6), 3e12, 1e-9)); // the step is relative for large |x|

    // exp varies faster than the scale |x| = 20: h = 1.2e-4 and the relative error is h^2 / 6 = 2.4e-9.
    CHECK(close(derivative(exp, -20.0), std::exp(-20.0), 5e-9));

    // For |x| < 1 the step is absolute, and so is the accuracy: the error of x^3 is h^2 = 3.7e-11 at any small x.
    double const h = nxx::math::root_epsilon<3, double>();
    CHECK(close(derivative(cube, 1e-6), 3e-12 + h * h, 1e-9));
    CHECK(close(derivative(cube, 0.0), h * h, 1e-9));
}

TEST_CASE("deriv: float")
{
    auto const sin = [](float x) { return std::sin(x); };
    CHECK(close(derivative(sin, 0.5f), std::cos(0.5f), 1e-3f));
}

TEST_CASE("deriv: multiprecision types")
{
    // eps^(2/3) is about 1e-33 for 50 decimal digits.
    auto const exp_bin = [](bin50 const &x) { return mp::exp(x); };
    CHECK(close(derivative(exp_bin, bin50(1)), bin50(mp::exp(bin50(1))), bin50("1e-30")));

    auto const exp_et = [](bin50_et const &x) { return bin50_et(mp::exp(x)); };
    CHECK(close(derivative(exp_et, bin50_et(1)), bin50_et(mp::exp(bin50_et(1))), bin50_et("1e-30")));

    auto const exp_dec = [](dec50 const &x) { return dec50(mp::exp(x)); };
    CHECK(close(derivative(exp_dec, dec50(1)), dec50(mp::exp(dec50(1))), dec50("1e-30")));

    // The step converts to the argument's type.
    auto const square_bin = [](bin50 const &x) { return bin50(x * x); };
    CHECK(derivative(square_bin, bin50(3), 0.5) == bin50(6));
}

TEST_CASE("deriv: fn is called twice")
{
    int calls = 0;
    auto fn = [&calls](double x) {
        ++calls;
        return x;
    };
    CHECK(derivative(fn, 1.0) == 1.0);
    CHECK(calls == 2);
}

TEST_CASE("deriv: a non-finite estimate is an error")
{
    auto const log = [](double x) { return std::log(x); }; // NaN to the left of 0
    CHECK(derivative(log, 0.0) == std::unexpected(nxx::deriv::error::non_finite));

    auto const inf = [](double) { return std::numeric_limits<double>::infinity(); };
    CHECK(derivative(inf, 1.0) == std::unexpected(nxx::deriv::error::non_finite));

    auto const huge = [](double x) {
        return x < 0.0 ? -std::numeric_limits<double>::max() : std::numeric_limits<double>::max();
    };
    CHECK(derivative(huge, 0.0, 0.5) == std::unexpected(nxx::deriv::error::non_finite)); // the difference overflows

    // Before the fix, the distance between the points overflowed to infinity and the result was a successful 0.
    CHECK(derivative(half, 0.0, dmax) == std::unexpected(nxx::deriv::error::non_finite));

    // fn is not called at a point that would overflow, nor when h is below the resolution of x.
    int calls = 0;
    auto const counted = [&calls](double x) {
        ++calls;
        return x;
    };
    CHECK(derivative(counted, dmax, 1.0) == std::unexpected(nxx::deriv::error::non_finite));
    CHECK(derivative(counted, 1.0, 1e-300) == std::unexpected(nxx::deriv::error::non_finite));
    CHECK(calls == 0);
}

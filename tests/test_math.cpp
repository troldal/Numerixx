#include <numerixx/math.hpp>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <doctest/doctest.h>

#include <cmath>
#include <limits>

namespace mp = boost::multiprecision;

using bin50 = mp::cpp_bin_float_50;
using bin50_et = mp::number<mp::cpp_bin_float<50>, mp::et_on>;
using dec50 = mp::cpp_dec_float_50;

static_assert(nxx::real<float> && nxx::real<double> && nxx::real<long double>);
static_assert(nxx::real<bin50> && nxx::real<bin50_et> && nxx::real<dec50>);
static_assert(!nxx::real<int> && !nxx::real<unsigned>);

static_assert(nxx::math::abs(-2.5) == 2.5 && nxx::math::abs(2.5) == 2.5 && nxx::math::abs(0.0) == 0.0);
static_assert(nxx::math::isfinite(0.0) && nxx::math::isfinite(std::numeric_limits<double>::max()));
static_assert(nxx::math::isfinite(std::numeric_limits<double>::denorm_min()));
static_assert(!nxx::math::isfinite(std::numeric_limits<double>::infinity()));
static_assert(!nxx::math::isfinite(-std::numeric_limits<double>::infinity()));
static_assert(!nxx::math::isfinite(std::numeric_limits<double>::quiet_NaN()));

// Usable in constant expressions; the cube of the result is epsilon to within a few rounding errors.
static_assert(nxx::math::cbrt_epsilon<double>() > 6.05e-6 && nxx::math::cbrt_epsilon<double>() < 6.06e-6);

namespace {

    // |a - b| <= k * eps * |b|
    template<typename T>
    bool close(T const &a, T const &b, int k)
    {
        using std::abs;
        return abs(a - b) <= T(k) * std::numeric_limits<T>::epsilon() * abs(b);
    }

} // namespace

TEST_CASE("math: cbrt_epsilon matches cbrt(epsilon)")
{
    CHECK(close(nxx::math::cbrt_epsilon<float>(), std::cbrt(std::numeric_limits<float>::epsilon()), 2));
    CHECK(close(nxx::math::cbrt_epsilon<double>(), std::cbrt(std::numeric_limits<double>::epsilon()), 2));
    CHECK(close(nxx::math::cbrt_epsilon<long double>(), std::cbrt(std::numeric_limits<long double>::epsilon()), 2));

    // Standalone Boost.Multiprecision has no cbrt; pow(eps, 1/3) is correct to a few eps.
    CHECK(close(
        nxx::math::cbrt_epsilon<bin50>(), bin50(mp::pow(std::numeric_limits<bin50>::epsilon(), 1 / bin50(3))), 4));
    CHECK(close(nxx::math::cbrt_epsilon<bin50_et>(),
        bin50_et(mp::pow(std::numeric_limits<bin50_et>::epsilon(), 1 / bin50_et(3))),
        4));
    CHECK(close(
        nxx::math::cbrt_epsilon<dec50>(), dec50(mp::pow(std::numeric_limits<dec50>::epsilon(), 1 / dec50(3))), 4));
}

TEST_CASE("math: cbrt_epsilon is the same at compile time and at run time")
{
    constexpr double at_compile_time = nxx::math::cbrt_epsilon<double>();
    double const at_run_time = nxx::math::cbrt_epsilon<double>();
    CHECK(at_compile_time == at_run_time);
}

TEST_CASE("math: abs and isfinite on multiprecision types")
{
    CHECK(nxx::math::abs(bin50(-3)) == bin50(3));
    CHECK(nxx::math::isfinite(bin50(1)));
    CHECK(!nxx::math::isfinite(std::numeric_limits<bin50>::infinity()));
    CHECK(!nxx::math::isfinite(std::numeric_limits<bin50>::quiet_NaN()));
    CHECK(nxx::math::abs(dec50(-3)) == dec50(3));
    CHECK(!nxx::math::isfinite(std::numeric_limits<dec50>::infinity()));
}

#include <numerixx/math.hpp>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
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
// Exact and unbounded: epsilon is 0 and max is not a bound, so root_epsilon and isfinite would be meaningless.
static_assert(!nxx::real<mp::cpp_rational>);

static_assert(nxx::math::abs(-2.5) == 2.5 && nxx::math::abs(2.5) == 2.5 && nxx::math::abs(0.0) == 0.0);
static_assert(nxx::math::isfinite(0.0) && nxx::math::isfinite(std::numeric_limits<double>::max()));
static_assert(nxx::math::isfinite(std::numeric_limits<double>::denorm_min()));
static_assert(!nxx::math::isfinite(std::numeric_limits<double>::infinity()));
static_assert(!nxx::math::isfinite(-std::numeric_limits<double>::infinity()));
static_assert(!nxx::math::isfinite(std::numeric_limits<double>::quiet_NaN()));

namespace {

    // |a - b| <= k * eps * |b|
    template<typename T>
    constexpr bool close(T const &a, T const &b, int k)
    { return nxx::math::abs(T(a - b)) <= T(k) * std::numeric_limits<T>::epsilon() * nxx::math::abs(b); }

} // namespace

// References: eps^(1/N) = 2^(-52/N) for double (2^(-23/N) for float), computed to 90 digits independently and correctly
// rounded. root_epsilon is within an ulp of the root (the last Newton update rounds; measured within 0.4 eps) and the
// reference within half an ulp, so 2 eps covers both.
static_assert(nxx::math::root_epsilon<2, double>() == 0x1p-26); // exact
static_assert(close(nxx::math::root_epsilon<3, double>(), 0x1.965fea53d6e3dp-18, 2));
static_assert(close(nxx::math::root_epsilon<5, double>(), 0x1.8406003b2ae5cp-11, 2));
static_assert(close(nxx::math::root_epsilon<7, double>(), 0x1.7c6a1f29e2ce6p-8, 2));
static_assert(close(nxx::math::root_epsilon<11, double>(), 0x1.3545598e5c232p-5, 2));
static_assert(close(nxx::math::root_epsilon<3, float>(), 0x1.428a3p-8f, 2));

TEST_CASE("math: root_epsilon matches the exact roots")
{
    // The same cases as the static_asserts above, at run time.
    CHECK(nxx::math::root_epsilon<2, double>() == 0x1p-26);
    CHECK(close(nxx::math::root_epsilon<3, double>(), 0x1.965fea53d6e3dp-18, 2));
    CHECK(close(nxx::math::root_epsilon<5, double>(), 0x1.8406003b2ae5cp-11, 2));
    CHECK(close(nxx::math::root_epsilon<7, double>(), 0x1.7c6a1f29e2ce6p-8, 2));
    CHECK(close(nxx::math::root_epsilon<11, double>(), 0x1.3545598e5c232p-5, 2));
    CHECK(close(nxx::math::root_epsilon<3, float>(), 0x1.428a3p-8f, 2));
    // long double differs by platform; std::cbrt is within an ulp for the exact argument eps.
    CHECK(close(nxx::math::root_epsilon<3, long double>(), std::cbrt(std::numeric_limits<long double>::epsilon()), 2));

    // cpp_bin_float_50 has 168 bits, eps = 2^-167; cpp_dec_float_50 has eps = 10^-49. References to 60 digits.
    CHECK(close(nxx::math::root_epsilon<3, bin50>(),
        bin50("1.74849169850398606480436442988966643043193186599291719119095e-17"),
        2));
    CHECK(close(nxx::math::root_epsilon<3, bin50_et>(),
        bin50_et("1.74849169850398606480436442988966643043193186599291719119095e-17"),
        2));
    CHECK(close(nxx::math::root_epsilon<7, bin50>(),
        bin50("6.58088632627613745718899600281324816540078282821252444693168e-8"),
        2));
    CHECK(close(nxx::math::root_epsilon<3, dec50>(),
        dec50("4.64158883361277889241007635091944657655134912501124363765069e-17"),
        2));
    CHECK(close(nxx::math::root_epsilon<11, dec50>(),
        dec50("3.51119173421513132134864589708621504202017817462297471566355e-5"),
        2));
}

TEST_CASE("math: root_epsilon is the same at compile time and at run time")
{
    constexpr double at_compile_time = nxx::math::root_epsilon<7, double>();
    double const at_run_time = nxx::math::root_epsilon<7, double>();
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

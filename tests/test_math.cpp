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

    // lo <= a <= hi, for the two floating-point neighbours lo < hi of an exact root: a is within an ulp of the root.
    template<typename T>
    constexpr bool brackets(T const &a, T const &lo, T const &hi)
    { return lo <= a && a <= hi; }

} // namespace

// References: eps^(1/N) = 2^(-52/N) for double (2^(-23/N) for float), computed to 80 digits with Python's decimal,
// independently of Numerixx; lo and hi are the two values of the type on either side of it, so brackets() enforces
// the documented accuracy of an ulp exactly.
static_assert(nxx::math::root_epsilon<2, double>() == 0x1p-26); // exact
static_assert(brackets(nxx::math::root_epsilon<3, double>(), 0x1.965fea53d6e3cp-18, 0x1.965fea53d6e3dp-18));
static_assert(brackets(nxx::math::root_epsilon<5, double>(), 0x1.8406003b2ae5cp-11, 0x1.8406003b2ae5dp-11));
static_assert(brackets(nxx::math::root_epsilon<7, double>(), 0x1.7c6a1f29e2ce6p-8, 0x1.7c6a1f29e2ce7p-8));
static_assert(brackets(nxx::math::root_epsilon<11, double>(), 0x1.3545598e5c232p-5, 0x1.3545598e5c233p-5));
static_assert(brackets(nxx::math::root_epsilon<3, float>(), 0x1.428a2ep-8f, 0x1.428a30p-8f));
static_assert(brackets(nxx::math::root_epsilon<5, float>(), 0x1.51cb44p-5f, 0x1.51cb46p-5f));
static_assert(brackets(nxx::math::root_epsilon<7, float>(), 0x1.a402fep-4f, 0x1.a40300p-4f));
static_assert(brackets(nxx::math::root_epsilon<11, float>(), 0x1.e0bb8cp-3f, 0x1.e0bb8ep-3f));
// Regression: radix^N used to be formed up front and overflowed for large N (2^128 for float), which constant
// evaluation rejects.
static_assert(brackets(nxx::math::root_epsilon<128, float>(), 0x1.c40ab4p-1f, 0x1.c40ab6p-1f));
static_assert(brackets(nxx::math::root_epsilon<2000, double>(), 0x1.f6db004386a06p-1, 0x1.f6db004386a07p-1));
// Regression: the update ((N - 1) y + a / y^(N-1)) / N rounded a sum of size N y, and its stopping test fired up to
// 1.8 ulps above the root: 1.64 ulps for float and N = 87 (0x1.aa45c6p-1), 1.77 for N = 146, 1.47 for double and
// N = 42, 1.59 for N = 297, the worst cases of N = 2 ... 200 (float) and 2 ... 300 (double).
static_assert(brackets(nxx::math::root_epsilon<87, float>(), 0x1.aa45c2p-1f, 0x1.aa45c4p-1f));
static_assert(brackets(nxx::math::root_epsilon<146, float>(), 0x1.cb0966p-1f, 0x1.cb0968p-1f));
static_assert(brackets(nxx::math::root_epsilon<42, double>(), 0x1.b21b3aa923da0p-2, 0x1.b21b3aa923da1p-2));
static_assert(brackets(nxx::math::root_epsilon<297, double>(), 0x1.c57c8beb88756p-1, 0x1.c57c8beb88757p-1));

TEST_CASE("math: root_epsilon matches the exact roots")
{
    // The same cases as the static_asserts above, at run time.
    CHECK(nxx::math::root_epsilon<2, double>() == 0x1p-26);
    CHECK(brackets(nxx::math::root_epsilon<3, double>(), 0x1.965fea53d6e3cp-18, 0x1.965fea53d6e3dp-18));
    CHECK(brackets(nxx::math::root_epsilon<5, double>(), 0x1.8406003b2ae5cp-11, 0x1.8406003b2ae5dp-11));
    CHECK(brackets(nxx::math::root_epsilon<7, double>(), 0x1.7c6a1f29e2ce6p-8, 0x1.7c6a1f29e2ce7p-8));
    CHECK(brackets(nxx::math::root_epsilon<11, double>(), 0x1.3545598e5c232p-5, 0x1.3545598e5c233p-5));
    CHECK(brackets(nxx::math::root_epsilon<3, float>(), 0x1.428a2ep-8f, 0x1.428a30p-8f));
    CHECK(brackets(nxx::math::root_epsilon<5, float>(), 0x1.51cb44p-5f, 0x1.51cb46p-5f));
    CHECK(brackets(nxx::math::root_epsilon<7, float>(), 0x1.a402fep-4f, 0x1.a40300p-4f));
    CHECK(brackets(nxx::math::root_epsilon<11, float>(), 0x1.e0bb8cp-3f, 0x1.e0bb8ep-3f));
    CHECK(brackets(nxx::math::root_epsilon<128, float>(), 0x1.c40ab4p-1f, 0x1.c40ab6p-1f));
    CHECK(brackets(nxx::math::root_epsilon<2000, double>(), 0x1.f6db004386a06p-1, 0x1.f6db004386a07p-1));
    CHECK(brackets(nxx::math::root_epsilon<87, float>(), 0x1.aa45c2p-1f, 0x1.aa45c4p-1f));
    CHECK(brackets(nxx::math::root_epsilon<146, float>(), 0x1.cb0966p-1f, 0x1.cb0968p-1f));
    CHECK(brackets(nxx::math::root_epsilon<42, double>(), 0x1.b21b3aa923da0p-2, 0x1.b21b3aa923da1p-2));
    CHECK(brackets(nxx::math::root_epsilon<297, double>(), 0x1.c57c8beb88756p-1, 0x1.c57c8beb88757p-1));
    // long double differs by platform; std::cbrt is within an ulp for the exact argument eps, and an ulp is at most
    // eps relative, so 2 eps covers both errors.
    CHECK(close(nxx::math::root_epsilon<3, long double>(), std::cbrt(std::numeric_limits<long double>::epsilon()), 2));

    // cpp_bin_float_50 has 168 bits, eps = 2^-167; cpp_dec_float_50 has eps = 10^-49. References to 60 digits, far
    // beyond either type, so the error is root_epsilon's alone: within an ulp, which is at most eps relative.
    CHECK(close(nxx::math::root_epsilon<3, bin50>(),
        bin50("1.74849169850398606480436442988966643043193186599291719119095e-17"),
        1));
    CHECK(close(nxx::math::root_epsilon<3, bin50_et>(),
        bin50_et("1.74849169850398606480436442988966643043193186599291719119095e-17"),
        1));
    CHECK(close(nxx::math::root_epsilon<7, bin50>(),
        bin50("6.58088632627613745718899600281324816540078282821252444693168e-8"),
        1));
    CHECK(close(nxx::math::root_epsilon<3, dec50>(),
        dec50("4.64158883361277889241007635091944657655134912501124363765069e-17"),
        1));
    CHECK(close(nxx::math::root_epsilon<11, dec50>(),
        dec50("3.51119173421513132134864589708621504202017817462297471566355e-5"),
        1));
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

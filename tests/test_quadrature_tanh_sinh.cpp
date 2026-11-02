// TanhSinh of nxx::quadrature (Phase 3b, docs/quadrature/PHASE3B_API.md, with the numerics of PHASE3_API.md). Step 1:
// detail::exp, the exponential in T from + - * / on the domain the node tables need.
//
// exp() is checked against mpmath at 60 digits within the derived 1.48 ulp (PHASE3_API.md B: the final rounding of
// 1 + p, 0.5 ulp; Horner's scheme, 0.73; the rounding of r, 0.25); the test allows 1.5. The references are in
// quadrature_tanh_sinh_reference.hpp, from the script in docs/quadrature/design.md.

#include "quadrature_check.hpp"
#include "quadrature_tanh_sinh_reference.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <span>

using namespace quadrature_check;

namespace {

    namespace ref = tanh_sinh_reference;

    // x = k / 2^s, exact in every type tested (the script keeps k within the type's digits).
    template<typename T>
    T dyadic(long long k, int s)
    { return T(T(k) / T(1LL << s)); }

    // |v / radix^e - m| for a reference value m radix^e, 1 <= m < radix, measured in bin100: in units of eps it is the
    // error in ulps of the reference's binade. A radix-2 value is split exactly by frexp first, since
    // cpp_bin_float_50's smallest values lie below cpp_bin_float_100's range.
    template<typename T>
    double scaled_error(T const &v, char const *mantissa, long long exponent)
    {
        bin100 const m(mantissa);
        using std::frexp, std::ldexp;
        if constexpr (std::numeric_limits<T>::radix == 2) {
            int e = 0;
            T const f = frexp(v, &e); // v = f 2^e, 1/2 <= f < 1; std::frexp or Boost's, found by lookup
            return static_cast<double>(bin100(abs(bin100(ldexp(wide(f), static_cast<int>(e - exponent))) - m)));
        }
        else {
            bin100 const scale = pow(bin100(10), static_cast<int>(exponent));
            return static_cast<double>(bin100(abs(bin100(wide(v) / scale) - m)));
        }
    }

    // exp() at the script's points of its domain, in ulps, within the derived 1.48 (the test allows 1.5); 0 at and
    // below x_min, and normal just above it.
    template<typename T>
    double check_exp(std::span<ref::ExpPoint const> points)
    {
        qd::detail::ExpConstants<T> const constants = qd::detail::exp_constants<T>();
        double const eps = eps_of<T>();
        double worst = 0.0;
        for (ref::ExpPoint const &point : points) {
            T const x = dyadic<T>(point.k, point.s);
            double const ulps = scaled_error(qd::detail::exp(x, constants), point.mantissa, point.exponent) / eps;
            INFO("x = " << point.k << " / 2^" << point.s);
            CHECK(ulps <= 1.5);
            worst = std::max(worst, ulps);
        }
        CHECK(qd::detail::exp(constants.x_min, constants) == T(0));
        CHECK(qd::detail::exp(T(constants.x_min - T(1)), constants) == T(0));
        T const inside = T(constants.x_min * T(T(1) - T(T(64) * T(std::numeric_limits<T>::epsilon()))));
        CHECK(qd::detail::exp(inside, constants) >= T((std::numeric_limits<T>::min)()));
        return worst;
    }

    constexpr qd::detail::ExpConstants<double> exp_double = qd::detail::exp_constants<double>();

} // namespace

// exp() at both ends of its domain in double, within its 1.5 ulp of mpmath's value plus the half ulp of rounding that
// to a double: e^-708 = 3.307553003638408e-308 in [2^-1022, 2^-1021), an ulp 2^-1074, just above min(); e^6.5 =
// 665.1416330443618 in [2^9, 2^10), an ulp 2^-43 (the double table needs t up to 6.56); and 0 at x_min.
static_assert(nxx::math::abs(qd::detail::exp(-708.0, exp_double) - 3.307553003638408e-308) <= 2.0 * 0x1p-1074);
static_assert(nxx::math::abs(qd::detail::exp(6.5, exp_double) - 665.1416330443618) <= 2.0 * 0x1p-43);
static_assert(qd::detail::exp(exp_double.x_min, exp_double) == 0.0);

TEST_CASE("quadrature: tanh_sinh: exp against mpmath")
{
    double const worst_float = check_exp<float>(ref::exp_float);
    double const worst_double = check_exp<double>(ref::exp_double);
    double worst_long = 0.0;
    if constexpr (std::numeric_limits<long double>::min_exponent == -16381)
        worst_long = check_exp<long double>(ref::exp_long_double); // x87 and binary128 share the domain
    else worst_long = check_exp<long double>(ref::exp_double);
    double const worst_bin50 = check_exp<bin50>(ref::exp_bin50);
    double const worst_bin50_et = check_exp<bin50_et>(ref::exp_bin50);
    double const worst_dec50 = check_exp<dec50>(ref::exp_dec50);
    MESSAGE("exp, largest error in ulps: float " << worst_float << ", double " << worst_double << ", long double "
                                                 << worst_long << ", cpp_bin_float_50 " << worst_bin50 << " ("
                                                 << worst_bin50_et << " with expression templates), cpp_dec_float_50 "
                                                 << worst_dec50);
    // The compile-time constants equal the run-time ones (the static_asserts above use them).
    qd::detail::ExpConstants<double> const runtime = qd::detail::exp_constants<double>();
    CHECK(runtime.c1 == exp_double.c1);
    CHECK(runtime.c2 == exp_double.c2);
    CHECK(runtime.x_min == exp_double.x_min);
}

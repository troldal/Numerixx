// Maths helpers (DESIGN §6.1): thin wrappers over the ADL idiom, so built-in floats use std:: and multiprecision types
// find their own functions. abs and isfinite are also usable in constant expressions (comparisons only there); sqrt
// is run time only, because <cmath> is not constexpr. The power-of-two constants are exact on every platform.
#pragma once

#include <numerixx/config.hpp>
#include <numerixx/core/scalar.hpp>

#include <cmath>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <type_traits>

NXX_BEGIN_HEADER

namespace nxx::math
{
    template<real T>
    constexpr T abs(const T& x) noexcept
    {
        if consteval { return x < T(0) ? -x : x; }
        else {
            using std::abs;    // block scope: unqualified lookup stops here, so this never recurses into nxx::math::abs
            return abs(x);
        }
    }

    template<real T>
    constexpr bool isfinite(const T& x) noexcept
    {
        if consteval {    // comparisons only: NaN compares unequal to itself, infinities lie outside [lowest, max]
            return x == x && x <= (std::numeric_limits<T>::max)() && x >= std::numeric_limits<T>::lowest();
        }
        else {
            using std::isfinite;
            return isfinite(x);
        }
    }

    template<real T>
    constexpr bool isnan(const T& x) noexcept
    { return !(x == x); }

    template<real T>
    T sqrt(const T& x) noexcept
    {
        using std::sqrt;
        return sqrt(x);
    }

    // 2^e, exact for every e the type can represent.
    template<class P>
    constexpr P pow2(int e) noexcept
    {
        P r(1);
        if (e >= 0)
            for (int i = 0; i < e; ++i) r *= P(2);
        else
            for (int i = 0; i < -e; ++i) r /= P(2);
        return r;
    }

    // eps^(num/den) rounded to a power of two: root_eps<double>(1, 2) == 0x1p-26 == sqrt(eps) rounded.
    template<class P>
    constexpr P root_eps(int num, int den) noexcept
    { return pow2<P>(-((std::numeric_limits<P>::digits - 1) * num) / den); }

    // The midpoint of [a, b] without overflow, also for a = -1.7e308, b = 1.7e308.
    template<real T>
    constexpr T midpoint(const T& a, const T& b) noexcept
    {
        if constexpr (std::is_floating_point_v<T>)
            return std::midpoint(a, b);
        else if ((a < T(0)) != (b < T(0)))
            return a / T(2) + b / T(2);
        else
            return a + (b - a) / T(2);
    }
}    // namespace nxx::math

NXX_END_HEADER

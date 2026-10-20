// Numerixx: the real-number concept and the constexpr helpers the algorithms need.
//
// The helpers use only comparisons and + - * /, so they work in constant expressions, give the same result on every
// compiler, and accept user-defined types such as boost::multiprecision::cpp_bin_float_50.

#pragma once

#include <concepts>
#include <limits>

namespace nxx {

    /// A real-number type: a built-in floating-point type, or a user type that specialises std::numeric_limits and
    /// provides the arithmetic operators.
    template<typename T>
    concept real = std::numeric_limits<T>::is_specialized && !std::numeric_limits<T>::is_integer
                   && std::constructible_from<T, int> && std::totally_ordered<T> && requires(T const a, T const b) {
                          { a + b } -> std::convertible_to<T>;
                          { a - b } -> std::convertible_to<T>;
                          { a * b } -> std::convertible_to<T>;
                          { a / b } -> std::convertible_to<T>;
                          { -a } -> std::convertible_to<T>;
                      };

    namespace math {

        template<real T>
        constexpr T abs(T const &x)
        {
            if (x < T(0)) return T(-x);
            return x;
        }

        /// False for infinities and NaN (a NaN is the only value that compares unequal to itself).
        template<real T>
        constexpr bool isfinite(T const &x)
        { return x == x && math::abs(x) <= (std::numeric_limits<T>::max)(); }

        /// The cube root of machine epsilon.
        ///
        /// Epsilon is first scaled into (radix^-3, 1] by powers of radix^3, which is exact. Newton's iteration for the
        /// cube root of the scaled value then starts at 1, above the root, and decreases until rounding stops it.
        template<real T>
        constexpr T cbrt_epsilon()
        {
            T const radix = static_cast<T>(std::numeric_limits<T>::radix);
            T const radix3 = radix * radix * radix;

            T a = std::numeric_limits<T>::epsilon();
            T scale = T(1);
            while (a * radix3 <= T(1)) {
                a = a * radix3;
                scale = scale * radix;
            }

            T y = T(1);
            while (true) {
                T const next = (T(2) * y + a / (y * y)) / T(3);
                if (!(next < y)) break;
                y = next;
            }
            return y / scale;
        }

    } // namespace math
} // namespace nxx

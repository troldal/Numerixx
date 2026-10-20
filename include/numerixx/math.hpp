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
    ///
    /// The type must be bounded and inexact, like a floating-point type: the algorithms rely on a finite max() and a
    /// positive epsilon(). Exact types such as rationals are rejected.
    //
    // Constructible from int, because constants are written T(1), T(2): a double literal would silently lose
    // precision in a wider type. The operators' results need only convert to T, because Boost.Multiprecision with
    // expression templates returns an expression object from a + b, not a T.
    template<typename T>
    concept real = std::numeric_limits<T>::is_specialized && !std::numeric_limits<T>::is_integer
                   && std::numeric_limits<T>::is_bounded && !std::numeric_limits<T>::is_exact
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
            if (x < T(0)) return T(-x); // T(...) evaluates a possible expression template
            return x;
        }

        /// False for infinities and NaN (a NaN is the only value that compares unequal to itself).
        // max is parenthesised so that a max macro, as defined by <windows.h>, cannot expand there.
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

            // Multiplying by a power of the radix changes only the exponent, so a and scale stay exact, and
            // eps = a / scale^3 gives cbrt(eps) = cbrt(a) / scale. For double, 2^-52 becomes a = 0.5 with scale = 2^17.
            T a = std::numeric_limits<T>::epsilon();
            T scale = T(1);
            while (a * radix3 <= T(1)) {
                a = a * radix3;
                scale = scale * radix;
            }

            // Newton's iteration for y^3 = a. From y = 1, at or above the root since a <= 1, every iterate stays above
            // the root (by the inequality of arithmetic and geometric means) and decreases; in floating point it stops
            // decreasing once it reaches the root to within rounding. The scaling put the root in (1/radix, 1], so this
            // takes a few steps; started from 1 against eps itself, it would take dozens.
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

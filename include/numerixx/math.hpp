/**
 * @file
 * @brief The real-number concept and the constexpr helpers shared by the Numerixx algorithms.
 *
 * The helpers use only comparisons and + - * /, so they work in constant expressions, give the same result on every
 * compiler, and accept user-defined types such as boost::multiprecision::cpp_bin_float_50.
 */

#pragma once

#include <concepts>
#include <limits>

/** @brief Numerixx: numerical algorithms for built-in and user-defined real types. */
namespace nxx {

    /**
     * @brief A real-number type: a built-in floating-point type, or a user type such as a Boost.Multiprecision number.
     *
     * The type specialises std::numeric_limits as bounded, inexact and non-integer, so it has a finite max() and a
     * positive epsilon(). It is constructible from int and totally ordered, and its + - * / and unary minus return
     * something convertible to the type, so expression templates are allowed. Exact types such as rationals are
     * rejected.
     *
     * @tparam T The candidate type.
     */
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

    /** @brief Constexpr helpers for any nxx::real, built from comparisons and + - * / only. */
    namespace math {

        /**
         * @brief The absolute value of @p x.
         * @param x The value.
         * @return |x|.
         */
        template<real T>
        constexpr T abs(T const &x)
        {
            if (x < T(0)) return T(-x); // T(...) evaluates a possible expression template
            return x;
        }

        /**
         * @brief Whether @p x is finite.
         * @param x The value.
         * @return false for infinities and NaN, true otherwise.
         */
        // A NaN is the only value that compares unequal to itself. max is parenthesised so that a max macro, as defined
        // by <windows.h>, cannot expand there.
        template<real T>
        constexpr bool isfinite(T const &x)
        { return x == x && math::abs(x) <= (std::numeric_limits<T>::max)(); }

        /**
         * @brief The N-th root of machine epsilon, the building block of finite-difference steps.
         *
         * Accurate to within an ulp, and the same at compile time and run time.
         *
         * @tparam N The root, N >= 2.
         * @tparam T The real type whose epsilon is taken.
         * @return eps^(1/N) for T.
         */
        // Epsilon is first scaled into (radix^-N, 1] by powers of radix^N, which is exact. Newton's iteration for the
        // N-th root of the scaled value then starts at 1, above the root, and decreases until rounding stops it.
        template<int N, real T>
        constexpr T root_epsilon()
        {
            static_assert(N >= 2);
            T const radix = static_cast<T>(std::numeric_limits<T>::radix);
            T radix_n = T(1);
            for (int i = 0; i < N; ++i) radix_n = radix_n * radix;

            // Multiplying by a power of the radix changes only the exponent, so a and scale stay exact, and
            // eps = a / scale^N gives eps^(1/N) = a^(1/N) / scale. For double and N = 3, 2^-52 becomes a = 0.5 with
            // scale = 2^17.
            T a = std::numeric_limits<T>::epsilon();
            T scale = T(1);
            while (a * radix_n <= T(1)) {
                a = a * radix_n;
                scale = scale * radix;
            }

            // Newton's iteration for y^N = a: y <- ((N - 1) y + a / y^(N-1)) / N. From y = 1, at or above the root
            // since a <= 1, every iterate stays above the root (by the inequality of arithmetic and geometric means)
            // and decreases; in floating point it stops decreasing once it reaches the root to within rounding. The
            // scaling put the root in (1/radix, 1], so this takes a few steps; started from 1 against eps itself, it
            // would take dozens.
            T y = T(1);
            while (true) {
                T y_n1 = T(1);
                for (int i = 1; i < N; ++i) y_n1 = y_n1 * y;
                T const next = (T(N - 1) * y + a / y_n1) / T(N);
                if (!(next < y)) break;
                y = next;
            }
            return y / scale;
        }

    } // namespace math
} // namespace nxx

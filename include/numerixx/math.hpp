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
#include <type_traits>

/** @brief Numerixx: numerical algorithms for built-in and user-defined real types. */
namespace nxx {

    /**
     * @brief A real-number type: a built-in floating-point type, or a user type such as a Boost.Multiprecision number.
     *
     * The type specialises std::numeric_limits as bounded, inexact and non-integer, with the members the library
     * reads: radix (at least 2) and digits (at least 1) as constant expressions, epsilon() and max() convertible to
     * the type, and denorm_min() from which the type can be constructed. It is constructible from int, copyable and
     * totally ordered, and its + - * / and unary minus return something convertible to the type, so expression
     * templates are allowed. Exact types such as rationals are rejected.
     *
     * @tparam T The candidate type.
     */
    // Constructible from int, because constants are written T(1), T(2): a double literal would silently lose precision
    // in a wider type. The operators' results need only convert to T, because Boost.Multiprecision with expression
    // templates returns an expression object from a + b, not a T. The numeric_limits members are those the headers
    // read: radix and digits in constant expressions (the exact-weights and float-warning constraints of derivative),
    // radix at least 2 because root_epsilon and the overflow guards scale by it; epsilon() and max() are
    // copy-initialised into T, while denorm_min() is read as T(denorm_min()), so a type that inherits
    // numeric_limits<double> and constructs from double only explicitly is accepted.
    template<typename T>
    concept real = std::numeric_limits<T>::is_specialized && !std::numeric_limits<T>::is_integer
                   && std::numeric_limits<T>::is_bounded && !std::numeric_limits<T>::is_exact
                   && std::constructible_from<T, int> && std::copyable<T> && std::totally_ordered<T>
                   && requires(T const a, T const b) {
                          requires std::numeric_limits<T>::radix >= 2 && std::numeric_limits<T>::digits >= 1;
                          { std::numeric_limits<T>::epsilon() } -> std::convertible_to<T>;
                          { (std::numeric_limits<T>::max)() } -> std::convertible_to<T>;
                          requires std::constructible_from<T, decltype(std::numeric_limits<T>::denorm_min())>;
                          { a + b } -> std::convertible_to<T>;
                          { a - b } -> std::convertible_to<T>;
                          { a * b } -> std::convertible_to<T>;
                          { a / b } -> std::convertible_to<T>;
                          { -a } -> std::convertible_to<T>;
                      };

    /** @cond INTERNAL */
    namespace detail {

        /**
         * @brief A function the algorithms can evaluate at a T: invocable with a T const lvalue, with a result
         * convertible to T.
         *
         * F & in the constraint, because fn is called as a named reference, an lvalue, through std::invoke with a T
         * const lvalue: member pointers work, and so does a callable that rejects rvalues. Spelled
         * nxx::detail::callable in the modules, whose own detail namespaces would otherwise hide it.
         *
         * @tparam F The callable type.
         * @tparam T The real type.
         */
        template<typename F, typename T>
        concept callable =
            std::invocable<F &, T const &> && std::convertible_to<std::invoke_result_t<F &, T const &>, T>;

    } // namespace detail
    /** @endcond */

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
        constexpr bool isfinite(T const &x) { return x == x && math::abs(x) <= (std::numeric_limits<T>::max)(); } // NOLINT

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

            // Multiplying by a power of the radix changes only the exponent, so scaled and scale stay exact, and eps =
            // scaled / scale^N gives eps^(1/N) = scaled^(1/N) / scale. For double and N = 3, 2^-52 becomes scaled = 0.5
            // with scale = 2^17. Each step tries scaled * radix^N one factor of radix at a time and stops as soon as
            // the value passes 1, so no intermediate exceeds radix: forming radix^N itself would overflow for large N
            // (2^128 for float), which constant evaluation rejects. For such N no step is taken, and scaled stays eps.
            T scaled = std::numeric_limits<T>::epsilon();
            T scale = T(1);
            while (true) {
                T trial = scaled;
                for (int i = 0; i < N && trial <= T(1); ++i) trial = trial * radix;
                if (!(trial <= T(1))) break;
                scaled = trial;
                scale = scale * radix;
            }

            // Newton's iteration for y^N = scaled: y <- y + (scaled / y^(N-1) - y) / N. From y = 1, at or above the
            // root since scaled <= 1, every iterate stays above the root (by the inequality of arithmetic and geometric
            // means) and decreases; in floating point it stops decreasing once it reaches the root to within rounding.
            // The scaling put the root in (1/radix, 1], so this takes a few steps; started from 1 against eps itself,
            // it would take dozens. y^(N-1) cannot underflow: if eps was scaled, then radix^N <= 1/eps and y > 1/radix,
            // so y^(N-1) > radix^(1-N) >= radix eps; otherwise scaled = eps and y^(N-1) >= eps^((N-1)/N) >= eps.
            //
            // The update is written as a correction to y. Near the root scaled / y^(N-1) - y is exact (the operands are
            // within a factor 2), and the rounding errors of scaled / y^(N-1), up to about N/2 ulps from the N - 1
            // products, are divided by N, so the result is within an ulp: at most 0.58 ulps for float and double, N = 2
            // ... 5000, against 60-digit roots. The textbook form ((N - 1) y + scaled / y^(N-1)) / N rounds a sum of
            // size N y, worth up to an ulp of y, so its stopping test could fire 1.8 ulps above the root.
            T y = T(1);
            while (true) {
                T y_n1 = T(1);
                for (int i = 1; i < N; ++i) y_n1 = y_n1 * y;
                T const next = T(y + T(scaled / y_n1 - y) / T(N));
                if (!(next < y)) break;
                y = next;
            }
            return y / scale;
        }

    } // namespace math
} // namespace nxx

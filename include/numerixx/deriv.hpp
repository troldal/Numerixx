// Numerixx: numerical derivatives by finite differences.

#pragma once

#include <numerixx/math.hpp>

#include <cassert>
#include <concepts>
#include <expected>
#include <limits>
#include <type_traits>

namespace nxx::deriv {

    enum class error {
        /// No finite estimate: fn returned an infinite or NaN value, the difference or the quotient overflowed, x + h
        /// or x - h would overflow, or h is below the resolution of x.
        non_finite
    };

    namespace detail {

        // The central difference has truncation error O(h^2) and rounding error O(eps / h). A step of
        // eps^(1/3) * max(|x|, 1) balances the two, scaled so that the step is relative for large |x|.
        template<real T>
        constexpr T default_step(T const &x)
        {
            T const scale = math::abs(x) < T(1) ? T(1) : math::abs(x);
            return math::cbrt_epsilon<T>() * scale;
        }

    } // namespace detail

    /// The derivative of fn at x by the central difference (fn(x + h) - fn(x - h)) / 2h, with error O(h^2).
    ///
    /// fn is called twice, and only at finite points. The divisor is the distance between the two points actually
    /// evaluated, since x + h and x - h are rounded. Precondition: h > 0.
    ///
    /// The function never throws by itself; an exception thrown by fn or by T's operations passes through unchanged.
    /// Failures are checked before the arithmetic that would produce them, so they are also error values in constant
    /// expressions, except an overflow of the final difference or quotient, which GCC rejects at compile time.
    template<real T, typename F>
    requires std::invocable<F &, T const &> && std::convertible_to<std::invoke_result_t<F &, T const &>, T>
    [[nodiscard]] constexpr std::expected<T, error> derivative(F &&fn, T const &x, std::type_identity_t<T> const &h)
    {
        assert(h > T(0));
        // With |x| + 2h <= max, the points and their distance are finite, so nothing overflows before fn is called.
        T const largest = (std::numeric_limits<T>::max)();
        if (!(h <= largest / T(4) && math::abs(x) <= largest - T(2) * h)) return std::unexpected(error::non_finite);
        T const xp = x + h;
        T const xm = x - h;
        T const dx = xp - xm;
        if (dx == T(0)) return std::unexpected(error::non_finite); // h is below the resolution of x
        T const fp = static_cast<T>(fn(xp));
        T const fm = static_cast<T>(fn(xm));
        if (!math::isfinite(fp) || !math::isfinite(fm)) return std::unexpected(error::non_finite);
        T const d = (fp - fm) / dx;
        if (!math::isfinite(d)) return std::unexpected(error::non_finite); // the difference or the quotient overflowed
        return d;
    }

    /// The central difference with the default step, eps^(1/3) * max(|x|, 1).
    ///
    /// The step assumes that fn varies on a scale of max(|x|, 1). The error is then about eps^(2/3) times the size of
    /// fn''' on that scale; pass a step when fn varies faster, or when the derivative is much smaller than that error.
    template<real T, typename F>
    requires std::invocable<F &, T const &> && std::convertible_to<std::invoke_result_t<F &, T const &>, T>
    [[nodiscard]] constexpr std::expected<T, error> derivative(F &&fn, T const &x)
    { return nxx::deriv::derivative(fn, x, detail::default_step(x)); }

} // namespace nxx::deriv

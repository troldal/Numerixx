// Numerixx: numerical derivatives by finite differences.

#pragma once

#include <numerixx/math.hpp>

#include <cassert>
#include <concepts>
#include <expected>
#include <type_traits>

namespace nxx::deriv {

    enum class error {
        non_finite ///< The estimate is infinite or NaN: fn returned a non-finite value, or the difference overflowed.
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
    /// fn is called twice. The divisor is the distance between the two points actually evaluated, since x + h and
    /// x - h are rounded. Precondition: h > 0.
    template<real T, typename F>
    requires std::invocable<F &, T const &> && std::convertible_to<std::invoke_result_t<F &, T const &>, T>
    [[nodiscard]] constexpr std::expected<T, error> derivative(F &&fn, T const &x, std::type_identity_t<T> const &h)
    {
        assert(h > T(0));
        T const xp = x + h;
        T const xm = x - h;
        T const fp = static_cast<T>(fn(xp));
        T const fm = static_cast<T>(fn(xm));
        T const d = (fp - fm) / (xp - xm);
        if (!math::isfinite(d)) return std::unexpected(error::non_finite);
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

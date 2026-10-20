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
        //
        // In detail, the error is about h^2/6 |f'''| + eps |f| / h, smallest at h = (3 eps |f| / |f'''|)^(1/3) with
        // an error of order eps^(2/3). f and f''' are unknown, so the step assumes that f varies on the scale
        // max(|x|, 1) and drops the constant 3^(1/3).
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
    //
    // By Taylor's theorem, (f(x + h) - f(x - h)) / 2h = f'(x) + h^2/6 f'''(x) + O(h^4): the even terms cancel.
    //
    // F & in the constraint, because fn is called as a named reference, an lvalue. h is a non-deduced context
    // (type_identity_t), so T comes from x alone and a double step converts to a multiprecision T instead of making
    // deduction fail.
    template<real T, typename F>
    requires std::invocable<F &, T const &> && std::convertible_to<std::invoke_result_t<F &, T const &>, T>
    [[nodiscard]] constexpr std::expected<T, error> derivative(F &&fn, T const &x, std::type_identity_t<T> const &h)
    {
        assert(h > T(0)); // a non-positive step is a bug in the caller, not a run-time failure
        // With |x| + 2h <= max, the points and their distance are finite, so nothing overflows before fn is called.
        // The check cannot overflow itself: 2h <= max/2, and max - 2h only subtracts.
        T const largest = (std::numeric_limits<T>::max)();
        if (!(h <= largest / T(4) && math::abs(x) <= largest - T(2) * h)) return std::unexpected(error::non_finite);
        T const xp = x + h;
        T const xm = x - h;
        // Divide by the distance between the points actually evaluated, not by 2h: x + h and x - h are rounded, and at
        // x = 1 with h = 1e-5 their spacing can differ from 2h by about 1e-11 relative, as much as the error aimed for.
        T const dx = xp - xm;
        if (dx == T(0)) return std::unexpected(error::non_finite); // h is below the resolution of x
        // fn may return an expression template, and every intermediate is a T const, never auto: an unevaluated
        // expression would hold references to its operands, which can dangle.
        T const fp = static_cast<T>(fn(xp));
        T const fm = static_cast<T>(fn(xm));
        // Checked before subtracting, so inf - inf is never formed: constant evaluation rejects arithmetic that
        // produces a NaN (and on GCC an overflow), so the failure would be a compile error instead of an error value.
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
    {
        // Qualified, so that argument-dependent lookup cannot pick a derivative() from the namespace of T or F.
        return nxx::deriv::derivative(fn, x, detail::default_step(x));
    }

} // namespace nxx::deriv

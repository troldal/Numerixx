// Intervals (DESIGN §6.2): bracket<T>, a finite interval with lo < hi. The integration domain interval<T> (any
// orientation) arrives with integrate in phase 7.
#pragma once

#include <numerixx/core/error.hpp>
#include <numerixx/core/math.hpp>
#include <numerixx/core/refined.hpp>
#include <numerixx/core/scalar.hpp>

#include <expected>

NXX_BEGIN_HEADER

namespace nxx
{
    template<real T>
    class bracket
    {
        T lo_;
        T hi_;

    public:
        using value_type      = T;
        using nxx_refined_tag = void;

        constexpr bracket(detail::trust_me, T lo, T hi) noexcept : lo_(lo), hi_(hi) {}

        // A literal: bracket{2.0, 1.0} or bracket{0.0, inf} does not compile. A run-time value does not compile either
        // ("not a constant expression"): use make(), or pass {lo, hi} straight to a solver.
        consteval bracket(T lo, T hi) : lo_(lo), hi_(hi)
        {
            if (!(lo < hi) || !math::isfinite(lo) || !math::isfinite(hi))
                detail::literal_violates_invariant("a bracket needs finite lo < hi");
        }

        // Re-orders the endpoints; equal or non-finite endpoints are an error.
        static constexpr auto make(T a, T b) noexcept -> std::expected<bracket, errc>
        {
            if (!math::isfinite(a) || !math::isfinite(b) || a == b) return std::unexpected(errc::invalid_input);
            return a < b ? bracket { detail::trust_me {}, a, b } : bracket { detail::trust_me {}, b, a };
        }

        constexpr T lo() const noexcept { return lo_; }
        constexpr T hi() const noexcept { return hi_; }
        constexpr T width() const noexcept { return hi_ - lo_; }                       // may overflow to inf on extreme brackets
        constexpr T half_width() const noexcept { return hi_ / T(2) - lo_ / T(2); }    // finite for every finite bracket
        constexpr T midpoint() const noexcept { return math::midpoint(lo_, hi_); }

        friend constexpr bool operator==(const bracket&, const bracket&) = default;
    };

    template<real T>
    bracket(T, T) -> bracket<T>;
}    // namespace nxx

NXX_END_HEADER

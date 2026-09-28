// Errors and results (DESIGN §6.3): every failure is a value in std::expected, with a stable code, the algorithm that
// failed, what it cost, the best estimate so far and the user's own callback error, unchanged.
#pragma once

#include <numerixx/config.hpp>

#include <concepts>
#include <cstdint>
#include <expected>
#include <optional>
#include <type_traits>
#include <utility>

NXX_BEGIN_HEADER

namespace nxx
{
    enum class errc : std::uint8_t {
        // 1..31: input errors, found before iterating
        invalid_input = 1,
        no_sign_change,
        non_finite_input,
        dimension_mismatch,
        not_increasing,
        leading_coefficient_zero,
        out_of_domain,
        // 32..63: numerical errors, found while iterating
        budget_exhausted = 32,
        evaluations_exhausted,
        stalled,
        non_finite_value,
        zero_derivative,
        singular,
        line_search_failed,
        diverged,
        sign_change_not_root,
        local_minimum,
        // 64..: the user's callback said no
        callback_failed = 64
    };

    constexpr bool is_input_error(errc e) noexcept { return std::to_underlying(e) < 32; }
    constexpr bool is_budget_error(errc e) noexcept { return e == errc::budget_exhausted || e == errc::evaluations_exhausted; }

    // Open set of algorithm ids: roots 1-39, optimize 40-59, multiroots 60-79, integrate 80-99, deriv 100-109,
    // poly 110-119, users 200 and up. Each module names its ids (nxx::roots::algos::brent, ...).
    enum class algo : std::uint8_t { none = 0, user_first = 200 };

    enum class stop_reason : std::uint8_t { exact_zero, criterion, resolution_limit, algorithm };

    struct counters
    {
        std::uint32_t iterations  = 0;
        std::uint32_t evaluations = 0;    // in units of calls of the user's function (cost_of, D33)

        friend constexpr counters operator+(counters a, counters b) noexcept
        { return { a.iterations + b.iterations, a.evaluations + b.evaluations }; }
        friend constexpr bool operator==(counters, counters) = default;
    };

    // The cause of a plain callback, which cannot fail: takes no space.
    struct none
    {
        friend constexpr bool operator==(none, none) = default;
    };

    template<class UE>
    using cause_slot = std::conditional_t<std::is_same_v<UE, none>, none, std::optional<UE>>;

    // Per-evaluation and per-step error.
    template<class UE = none>
    struct fault
    {
        errc                  code {};
        std::uint32_t         evals = 0;    // evaluations consumed by the failing step
        NXX_NO_UNIQUE_ADDRESS cause_slot<UE> cause {};

        friend constexpr bool operator==(const fault&, const fault&) = default;
    };

    template<class Est>
    struct solution : Est
    {
        counters    used {};
        algo        by = algo::none;
        stop_reason how {};

        friend constexpr bool operator==(const solution&, const solution&) = default;
    };

    template<class Est, class UE = none>
    struct failure
    {
        using estimate_type = Est;
        using cause_type    = UE;

        errc                  code {};
        algo                  where = algo::none;
        counters              used {};
        std::optional<Est>    best {};                    // best estimate so far; nullopt only if nothing was evaluated
        NXX_NO_UNIQUE_ADDRESS cause_slot<UE> cause {};    // the user's callback error, unchanged

        friend constexpr bool operator==(const failure&, const failure&) = default;
    };

    template<class Est, class UE = none>
    using result = std::expected<solution<Est>, failure<Est, UE>>;

    // The solution's x, or the failure's best->x, or nothing.
    template<class R>
    constexpr auto best_x(const R& r)
    {
        using X = std::remove_cvref_t<decltype(r->x)>;
        if (r) return std::optional<X>(r->x);
        if (r.error().best) return std::optional<X>(r.error().best->x);
        return std::optional<X> {};
    }

    // nxx::is_fatal(e): whether a user's callback error stops first_of's fall-through. Customise it with a function
    // is_fatal(const YourError&) in your error type's namespace (found by ADL); the default is false.
    namespace detail::fatal
    {
        void is_fatal() = delete;    // poison pill: unqualified calls below find only ADL candidates

        struct is_fatal_fn
        {
            template<class E>
            constexpr bool operator()(const E& e) const noexcept
            {
                if constexpr (requires {
                                  { is_fatal(e) } -> std::convertible_to<bool>;
                              })
                    return is_fatal(e);
                else
                    return false;
            }
        };
    }    // namespace detail::fatal

    inline constexpr detail::fatal::is_fatal_fn is_fatal {};
}    // namespace nxx

NXX_END_HEADER

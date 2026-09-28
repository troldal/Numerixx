// Run-time solver chains (DESIGN §6.10), opt-in: the only header that uses std::function, included by neither
// <numerixx/core.hpp> nor <numerixx/numerixx.hpp>.
//
// any_solver<F, Est, UE> erases the type of any curried solver (a bound solver, a then_t, a first_of_t, another run-time
// chain) whose call const F& -> result<Est, UE> matches exactly, for one fixed callable type F such as
// std::function<double(double)>. first_of over a range of them has the same laziness, merge and cost accounting as
// the static first_of, and the result is itself an any_solver.
//
// Heap: std::function may allocate when a solver is wrapped or copied; a call with an existing F never allocates.
// Pass f as an F: a lambda that is not already an F is converted on every call.
#pragma once

#include <numerixx/config.hpp>
#include <numerixx/core/compose.hpp>
#include <numerixx/core/error.hpp>

#include <concepts>
#include <functional>
#include <optional>
#include <ranges>
#include <type_traits>
#include <utility>
#include <vector>

NXX_BEGIN_HEADER

namespace nxx
{
    // A curried solver for callables of type F: a copyable value s with s(f) -> R for const F& f. The call is checked
    // before copyability: copying a chain that holds an any_solver asks whether its copyable box converts to an
    // any_solver, and checking copy_constructible first would make that question depend on itself (Clang rejects it).
    template<class S, class F, class R>
    concept curried_solver_for =
        std::invocable<const S&, const F&> && std::same_as<std::invoke_result_t<const S&, const F&>, R> && std::copy_constructible<S>;

    // Never empty: there is no default constructor and there are no move operations (a move is a copy), because a
    // moved-from std::function is unspecified and may be empty.
    template<class F, class Est, class UE = none>
    class any_solver
    {
    public:
        using function_type = F;
        using estimate_type = Est;
        using cause_type    = UE;
        using result_type   = result<Est, UE>;

        template<class S>
            requires(!std::same_as<std::remove_cvref_t<S>, any_solver> && curried_solver_for<std::remove_cvref_t<S>, F, result<Est, UE>>)
        any_solver(S&& s) : impl_(std::forward<S>(s))    // implicit: every matching solver value is an any_solver
        {}

        // Deleted with a reason only for a curried solver of the wrong kind (callable with F, other result type): a
        // catch-all would make any_solver a viable conversion target for every type and overload sets ambiguous.
        template<class S>
            requires(!std::same_as<std::remove_cvref_t<S>, any_solver> && std::invocable<const std::remove_cvref_t<S>&, const F&> &&
                     !curried_solver_for<std::remove_cvref_t<S>, F, result<Est, UE>>)
        any_solver(S&&) NXX_DELETE("nxx::any_solver<F, Est, UE>: needs a copyable curried solver `const F& -> result<Est, UE>` "
                                   "with exactly this F, Est and UE");

        any_solver(const any_solver&)            = default;
        any_solver& operator=(const any_solver&) = default;    // replaces the whole value (strong guarantee)
        ~any_solver()                            = default;

        [[nodiscard]] result_type operator()(const F& fn) const { return impl_(fn); }

    private:
        std::function<result_type(const F&)> impl_;
    };

    template<class T>
    inline constexpr bool is_any_solver_v = false;
    template<class F, class Est, class UE>
    inline constexpr bool is_any_solver_v<any_solver<F, Est, UE>> = true;

    namespace detail
    {
        template<class Policy, class F, class Est, class UE>
        class runtime_first_of    // owns its alternatives; never changed after construction
        {
            NXX_NO_UNIQUE_ADDRESS copyable_box<Policy> policy_;
            std::vector<any_solver<F, Est, UE>>        alts_;

        public:
            runtime_first_of(Policy policy, std::vector<any_solver<F, Est, UE>> alts) : policy_(std::move(policy)), alts_(std::move(alts))
            {}

            [[nodiscard]] result<Est, UE> operator()(const F& fn) const
            {
                using R    = result<Est, UE>;
                using Fail = failure<Est, UE>;
                if (alts_.empty()) return R { std::unexpect, Fail { errc::invalid_input, algo::none, {}, std::nullopt, {} } };
                std::optional<Fail> failed;    // merged so far: the last code and cause, the best estimate, the total cost
                for (const auto& s : alts_) {
                    R r = s(fn);
                    if (r) {    // lazy: later alternatives never run; the success pays for the failed attempts
                        if (failed) r->used = r->used + failed->used;
                        return r;
                    }
                    const bool go_on = std::invoke(*policy_, r.error());
                    failed           = failed ? nxx::detail::merge(*failed, std::move(r).error()) : std::move(r).error();
                    if (!go_on) break;
                }
                return R { std::unexpect, *std::move(failed) };
            }
        };

        template<class Rng>
        inline constexpr bool any_solver_range_v =
            std::ranges::input_range<Rng> && is_any_solver_v<std::remove_cvref_t<std::ranges::range_reference_t<Rng>>>;

        template<class P, class Rng>
        auto runtime_chain(P policy, Rng alts)
        {
            using S     = std::remove_cvref_t<std::ranges::range_reference_t<Rng>>;
            using Chain = runtime_first_of<P, typename S::function_type, typename S::estimate_type, typename S::cause_type>;
            if constexpr (std::is_same_v<Rng, std::vector<S>>)
                return S { Chain { std::move(policy), std::move(alts) } };
            else
                return S { Chain { std::move(policy), std::ranges::to<std::vector<S>>(alts) } };
        }
    }    // namespace detail

    // first_of over a run-time range (std::vector, std::span, std::array, ...) of any_solver. The chain owns a copy (it
    // is curried, so borrowing the range would dangle). An empty range fails with errc::invalid_input when called: run-time
    // configuration can be empty, so this illegal state is representable and reported in-band (DESIGN §6.10, FLAG).
    // Same shape as the one-argument first_of(S) (one by-value parameter, a plain class head), so this constrained
    // overload is the more specialised one.
    template<class Rng>
        requires detail::any_solver_range_v<Rng>
    auto first_of(Rng alts)
    { return detail::runtime_chain(continue_unless_fatal {}, std::move(alts)); }

    template<class P, class Rng>
        requires detail::any_solver_range_v<Rng>
    auto first_of_with(P policy, Rng alts)
    { return detail::runtime_chain(std::move(policy), std::move(alts)); }
}    // namespace nxx

NXX_END_HEADER

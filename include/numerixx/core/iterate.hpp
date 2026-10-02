// The solver protocol and the single bounded-iteration driver (DESIGN §6.6, §6.7). Every iterative algorithm runs
// through nxx::iterate: one loop that owns the budget, the stop criterion, the observer and the best estimate.
#pragma once

#include <numerixx/core/criteria.hpp>
#include <numerixx/core/error.hpp>

#include <concepts>
#include <cstdint>
#include <expected>
#include <functional>
#include <optional>
#include <type_traits>
#include <utility>

NXX_BEGIN_HEADER

namespace nxx
{
    // What prepare() produces: the callable (a std::reference_wrapper for one-shot solves), the validated input and the
    // evaluations prepare() already spent. Exposed for manual stepping only.
    template<class F, class In>
    struct problem
    {
        F             f;
        In            in;
        std::uint32_t nfev0 = 0;
    };

    namespace detail
    {
        template<class A, class P>
        using init_result_t = decltype(std::declval<const A&>().init(std::declval<const P&>()));
        template<class A, class P>
        using state_t = typename init_result_t<A, P>::value_type;
        template<class A, class P>
        using init_failure_t = typename init_result_t<A, P>::error_type;
        template<class A, class P>
        using estimate_t = std::remove_cvref_t<decltype(std::declval<const A&>().estimate(std::declval<const state_t<A, P>&>()))>;

        // Whether e is a better failure payload than best: the family's better_than (found by ADL), else merit_of.
        template<class Est>
        constexpr bool better(const Est& e, const Est& best)
        {
            if constexpr (requires {
                              { better_than(e, best) } -> std::convertible_to<bool>;
                          })
                return better_than(e, best);
            else
                return merit_of(e) < merit_of(best);
        }
    }    // namespace detail

    // The protocol, all const: init(p) -> expected<S, failure>, step(p, s) -> expected<S, fault>, view(s) for the stop
    // criteria, estimate(s) on success, best(s) on failure, intrinsic(s) for the algorithm's own stops, options() for
    // the stop criterion, budget and observer, and s.nfev, the evaluations so far.
    template<class A, class P>
    concept iterative_solver_for = requires(const A& a, const P& p, const detail::state_t<A, P>& s) {
        { A::id } -> std::convertible_to<algo>;
        a.options();
        a.init(p);
        a.step(p, s);
        a.view(s);
        a.estimate(s);
        a.best(s);
        { a.intrinsic(s) } -> std::same_as<std::optional<stop_reason>>;
        { s.nfev } -> std::convertible_to<std::uint32_t>;
    };

    namespace detail
    {
        // One step; a fault becomes a failure that carries the best estimate and the failing step's evaluations.
        template<class A, class P, class S, class FEst>
        constexpr auto advance(const A& alg, const P& p, const S& s, std::uint32_t k, const FEst& best)
            -> std::expected<S, init_failure_t<A, P>>
        {
            using Fail = init_failure_t<A, P>;
            auto next  = alg.step(p, s);
            if (!next) {
                const auto& e = next.error();
                return std::unexpected(Fail { e.code, A::id, counters { k, s.nfev + e.evals }, best, e.cause });
            }
            return *std::move(next);
        }

        // A success, after the solver's optional post-condition (bracketing methods: the pole check, DESIGN §7.2):
        // finish(p, solution) -> std::optional<failure>, a failure when the post-condition does not hold. The result is
        // built once, in place; passing a std::expected through a by-value hook made GCC 16 report a false
        // -Wmaybe-uninitialized.
        template<class R, class A, class P, class SEst>
        constexpr R succeed(const A& alg, const P& p, solution<SEst> sol)
        {
            if constexpr (requires { alg.finish(p, sol); }) {
                if (auto fail = alg.finish(p, sol)) return R { std::unexpect, *std::move(fail) };
            }
            return R { std::in_place, std::move(sol) };
        }
    }    // namespace detail

    template<class A, class P>
        requires iterative_solver_for<A, P>
    constexpr auto iterate(const A& alg, const P& p)
    {
        using Fail = detail::init_failure_t<A, P>;    // failure<FEst, UE>
        using FEst = typename Fail::estimate_type;    // root_estimate, also for searchers
        using SEst = detail::estimate_t<A, P>;        // sign_bracket for searchers, else FEst
        using R    = std::expected<solution<SEst>, Fail>;

        auto first = alg.init(p);
        if (!first) return R { std::unexpect, std::move(first).error() };
        auto s    = *std::move(first);
        FEst best = alg.best(s);
        if (auto how = alg.intrinsic(s))
            return detail::succeed<R>(alg, p, solution<SEst> { alg.estimate(s), counters { 0, s.nfev }, A::id, *how });

        const auto&         opt = alg.options();
        const std::uint32_t n   = opt.budget.value();
        for (std::uint32_t k = 1; k != 0 && k <= n; ++k) {
            auto next = detail::advance(alg, p, s, k, best);
            if (!next) return R { std::unexpect, std::move(next).error() };
            if (const FEst e = alg.best(*next); detail::better(e, best)) best = e;    // best iterate on every exit
            opt.observe(alg.view(*next));                                             // logging lives here
            if (auto how = alg.intrinsic(*next))
                return detail::succeed<R>(alg, p, solution<SEst> { alg.estimate(*next), counters { k, next->nfev }, A::id, *how });
            switch (opt.stop(alg.view(s), alg.view(*next), counters { k, next->nfev })) {
                case verdict::converged:
                    return detail::succeed<R>(
                        alg,
                        p,
                        solution<SEst> { alg.estimate(*next), counters { k, next->nfev }, A::id, stop_reason::criterion });
                case verdict::stalled:
                    return R { std::unexpect, Fail { errc::stalled, A::id, counters { k, next->nfev }, best, {} } };
                case verdict::exhausted:
                    return R { std::unexpect, Fail { errc::evaluations_exhausted, A::id, counters { k, next->nfev }, best, {} } };
                case verdict::proceed:
                    break;
            }
            s = *std::move(next);    // the only mutation: a local
        }
        return R { std::unexpect, Fail { errc::budget_exhausted, A::id, counters { n, s.nfev }, best, {} } };
    }
}    // namespace nxx

NXX_END_HEADER

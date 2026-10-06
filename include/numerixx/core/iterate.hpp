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
    }    // namespace detail

    // nxx::better_than(a, b): whether failure estimate a is a better best estimate than b (DESIGN §6.6). It calls the
    // better_than(const Est&, const Est&) that ADL finds next to the estimate type: a hidden friend or a function in the
    // estimate type's namespace. It must be a strict weak order, so the best estimate of a chain does not depend on how
    // its alternatives are grouped or folded (§6.7). Without one it is not invocable (std::is_invocable_v is false); a
    // better_than in another namespace, or a member function, counts as missing. There is no fallback.
    namespace detail::better_cpo
    {
        void better_than() = delete;    // poison pill: unqualified calls below find only ADL candidates

        template<class E>
        inline constexpr bool found_v = requires(const E& a, const E& b) {
            { better_than(a, b) } -> std::convertible_to<bool>;
        };

        struct better_than_fn
        {
            // No deleted sibling: without found_v<E> the call is simply not invocable, and the compiler names found_v.
            template<class E>
                requires found_v<E>
            constexpr bool operator()(const E& a, const E& b) const noexcept(noexcept(static_cast<bool>(better_than(a, b))))
            { return static_cast<bool>(better_than(a, b)); }
        };
    }    // namespace detail::better_cpo

    inline constexpr detail::better_cpo::better_than_fn better_than {};

    namespace detail
    {
        // Whether the failure estimate type Est has an order that nxx::better_than can call.
        template<class Est>
        inline constexpr bool has_better_than_v = better_cpo::found_v<std::remove_cvref_t<Est>>;

        // Whether e is a better failure payload than best (DESIGN §6.6): the family's better_than. A <= order instead
        // of a < order trips the precondition in assert builds.
        template<class Est>
        constexpr bool better(const Est& e, const Est& best)
        {
            NXX_EXPECTS(!nxx::better_than(e, e));
            return nxx::better_than(e, best);
        }
    }    // namespace detail

    // The protocol, all const: init(p) -> expected<S, failure>, step(p, s) -> expected<S, fault>, view(s) for the stop
    // criteria, estimate(s) on success, best(s) on failure, intrinsic(s) for the algorithm's own stops, options() for
    // the stop criterion, budget and observer, and s.nfev, the evaluations so far. The failure estimate type (the
    // estimate_type of init's failure) must have an order, found by ADL through nxx::better_than (DESIGN §6.6).
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
        requires detail::has_better_than_v<typename detail::init_failure_t<A, P>::estimate_type>;
    };

    namespace detail
    {
        // A step never reports invalid_input or non_finite_input (DESIGN §6.3, §6.7): input codes mean "rejected before
        // iterating". nxx::evaluate already maps them when a callback returns a Numerixx fault (§6.4), and the
        // library's steps produce neither code; this is the backstop for a user-written step that returns one
        // directly. It becomes non_finite_value with the fault's evaluations and cause. Stage 2 of then and
        // warm_fallback maps its failure the same way, with its cost and cause: its input is stage 1's output, not the
        // caller's (§6.10). Not noexcept: moving a user cause may throw, and the library is exception-neutral (D10).
        template<class E>    // a fault<UE>, or a failure<Est, UE> of a stage 2
        constexpr E step_fault(E e)
        {
            if (e.code == errc::invalid_input || e.code == errc::non_finite_input) e.code = errc::non_finite_value;
            return e;
        }

        // alg.step(p, s) with step_fault applied: the driver (advance) and steps_view both step through it, so
        // steps_view yields what the driver sees.
        template<class A, class P, class S>
        constexpr auto checked_step(const A& alg, const P& p, const S& s)
        {
            auto next = alg.step(p, s);
            if (!next) return decltype(next) { std::unexpect, nxx::detail::step_fault(std::move(next).error()) };
            return next;
        }

        // One step; a fault becomes a failure that carries the best estimate and the failing step's evaluations.
        template<class A, class P, class S, class FEst>
        constexpr auto advance(const A& alg, const P& p, const S& s, std::uint32_t k, const FEst& best)
            -> std::expected<S, init_failure_t<A, P>>
        {
            using Fail = init_failure_t<A, P>;
            auto next  = nxx::detail::checked_step(alg, p, s);
            if (!next) {
                const auto& e = next.error();
                return std::unexpected(Fail { e.code, A::id, counters { k, s.nfev + e.evaluations }, best, e.cause });
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

// Combinators (DESIGN §6.10): solver values built from solver values. They are named class templates, not closures, so
// they are copy-assignable whenever their parts are copy-constructible.
//
//   first_of(s1, s2, ...)        the first success; on total failure the last code and cause, the best estimate over
//                                all attempts and the total cost. Falls through unless the user's error is_fatal.
//   first_of_with(policy, s...)  the same, with policy(const failure&) -> bool deciding whether to continue.
//   then(s1, s2, ...)            staging: stage 2 gets the function and stage 1's value (search -> bracketing solver,
//                                bracketing solver -> open method).
//   warm_fallback(s1, s2)        restart s2 from s1's best estimate when s1 fails.
#pragma once

#include <numerixx/config.hpp>
#include <numerixx/core/callable.hpp>
#include <numerixx/core/error.hpp>
#include <numerixx/core/iterate.hpp>

#include <expected>
#include <functional>
#include <type_traits>
#include <utility>

NXX_BEGIN_HEADER

namespace nxx
{
    // The default first_of policy: input and numerical errors try the next alternative; a callback error stops the
    // chain only if nxx::is_fatal says so.
    struct continue_unless_fatal
    {
        template<class Est, class UE>
        constexpr bool operator()(const failure<Est, UE>& e) const noexcept
        {
            if constexpr (std::is_same_v<UE, none>)
                return true;
            else
                return !(e.cause && nxx::is_fatal(*e.cause));
        }
    };

    namespace detail
    {
        // Two failures of a chain: the later code and cause, the better best estimate, the total cost.
        template<class Fail>
        constexpr Fail merge(const Fail& earlier, Fail later)
        {
            later.used = earlier.used + later.used;
            if (earlier.best && (!later.best || nxx::detail::better(*earlier.best, *later.best))) later.best = earlier.best;
            return later;
        }

        template<class R>
        constexpr R add_cost(R r, counters extra)
        {
            if (r)
                r->used = r->used + extra;
            else
                r.error().used = r.error().used + extra;
            return r;
        }

        // A stage-1 failure as the failure type of stage 2: the same type, or a failure without a cause (a plain
        // callback) widened to stage 2's cause type (a fallible derivative in stage 2).
        template<class To, class From>
        constexpr To rebind_failure(From f)
        {
            if constexpr (std::is_same_v<To, From>)
                return f;
            else {
                static_assert(std::is_same_v<typename From::estimate_type, typename To::estimate_type> &&
                                  std::is_same_v<typename From::cause_type, none>,
                              "nxx::then: the stages report different callback error types; map one with .transform_error so they agree");
                return To { f.code, f.where, f.used, f.best, {} };
            }
        }

        template<class S, class... A>
        inline constexpr bool callable_v = std::is_invocable_v<const S&, const A&...>;
        template<class S, class... A>
        using result_t = std::invoke_result_t<const S&, const A&...>;
        template<class...>
        inline constexpr bool always_false = false;
    }    // namespace detail

    template<class Policy, class S1, class S2>
    class first_of_t
    {
        NXX_NO_UNIQUE_ADDRESS detail::copyable_box<Policy> policy_;
        detail::copyable_box<S1>                           s1_;
        detail::copyable_box<S2>                           s2_;

    public:
        constexpr first_of_t(Policy policy, S1 s1, S2 s2) : policy_(std::move(policy)), s1_(std::move(s1)), s2_(std::move(s2)) {}

        template<class... A>
        constexpr auto operator()(const A&... a) const
        {
            if constexpr (!detail::callable_v<S1, A...> || !detail::callable_v<S2, A...>) {
                static_assert(detail::always_false<A...>,
                              "nxx::first_of: an alternative is not callable with these arguments; did you forget .on(input)?");
            }
            else if constexpr (!std::is_same_v<detail::result_t<S1, A...>, detail::result_t<S2, A...>>) {
                static_assert(detail::always_false<A...>,
                              "nxx::first_of: every alternative must return the same std::expected<solution<Est>, failure<Est, UE>>; "
                              "adapt the odd one with .transform/.transform_error");
                return std::invoke(*s1_, a...);
            }
            else {
                auto r1 = std::invoke(*s1_, a...);
                if (r1 || !std::invoke(*policy_, r1.error())) return r1;    // lazy: later alternatives never run
                auto r2 = std::invoke(*s2_, a...);
                if (r2) {
                    r2->used = r2->used + r1.error().used;    // success pays for the failed attempts
                    return r2;
                }
                return decltype(r1) { std::unexpect, detail::merge(r1.error(), std::move(r2).error()) };
            }
        }
    };

    template<class S>
    constexpr S first_of(S s)
    { return s; }

    template<class S1, class S2, class... Ss>
    constexpr auto first_of(S1 s1, S2 s2, Ss... ss)
    {
        if constexpr (sizeof...(Ss) == 0)
            return first_of_t<continue_unless_fatal, S1, S2> { {}, std::move(s1), std::move(s2) };
        else {
            using Rest = decltype(nxx::first_of(std::move(s2), std::move(ss)...));
            return first_of_t<continue_unless_fatal, S1, Rest> { {}, std::move(s1), nxx::first_of(std::move(s2), std::move(ss)...) };
        }
    }

    template<class P, class S>
    constexpr S first_of_with(P, S s)
    { return s; }

    template<class P, class S1, class S2, class... Ss>
    constexpr auto first_of_with(P policy, S1 s1, S2 s2, Ss... ss)
    {
        if constexpr (sizeof...(Ss) == 0)
            return first_of_t<P, S1, S2> { std::move(policy), std::move(s1), std::move(s2) };
        else {
            using Rest = decltype(nxx::first_of_with(policy, std::move(s2), std::move(ss)...));
            return first_of_t<P, S1, Rest> { policy, std::move(s1), nxx::first_of_with(policy, std::move(s2), std::move(ss)...) };
        }
    }

    template<class S1, class S2>
    class then_t
    {
        detail::copyable_box<S1> s1_;
        detail::copyable_box<S2> s2_;

    public:
        constexpr then_t(S1 s1, S2 s2) : s1_(std::move(s1)), s2_(std::move(s2)) {}

        template<class F>
        constexpr auto operator()(const F& fn) const
        {
            auto r1  = std::invoke(*s1_, fn);
            using V1 = typename decltype(r1)::value_type;
            if constexpr (!std::is_invocable_v<const S2&, const F&, const V1&>) {
                static_assert(detail::always_false<F>,
                              "nxx::then: stage 2 cannot start from stage 1's result (a bracketing solver needs a bracket or a "
                              "search result; use .from_enclosure() or put a searcher first)");
                return r1;
            }
            else {
                using R2 = std::invoke_result_t<const S2&, const F&, const V1&>;
                if (!r1) return R2 { std::unexpect, detail::rebind_failure<typename R2::error_type>(std::move(r1).error()) };
                return detail::add_cost(std::invoke(*s2_, fn, *r1), r1->used);
            }
        }
    };

    template<class S1, class S2>
    constexpr auto then(S1 s1, S2 s2)
    { return then_t<S1, S2> { std::move(s1), std::move(s2) }; }

    template<class S1, class S2, class S3, class... Ss>
    constexpr auto then(S1 s1, S2 s2, S3 s3, Ss... ss)
    { return nxx::then(nxx::then(std::move(s1), std::move(s2)), std::move(s3), std::move(ss)...); }

    template<class S1, class S2>
    class warm_fallback_t
    {
        detail::copyable_box<S1> s1_;
        detail::copyable_box<S2> s2_;

    public:
        constexpr warm_fallback_t(S1 s1, S2 s2) : s1_(std::move(s1)), s2_(std::move(s2)) {}

        template<class F>
        constexpr auto operator()(const F& fn) const
        {
            auto r1 = std::invoke(*s1_, fn);
            if (r1 || !r1.error().best) return r1;
            auto r2 = std::invoke(*s2_, fn, *r1.error().best);
            static_assert(std::is_same_v<decltype(r1), decltype(r2)>, "nxx::warm_fallback: both stages must return the same result type");
            if (r2) {
                r2->used = r2->used + r1.error().used;
                return r2;
            }
            return decltype(r1) { std::unexpect, detail::merge(r1.error(), std::move(r2).error()) };
        }
    };

    template<class S1, class S2>
    constexpr auto warm_fallback(S1 s1, S2 s2)
    { return warm_fallback_t<S1, S2> { std::move(s1), std::move(s2) }; }
}    // namespace nxx

NXX_END_HEADER

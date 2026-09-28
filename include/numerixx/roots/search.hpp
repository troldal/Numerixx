// Bracket search (DESIGN §7.2): expand grows a start window until the function changes sign. Its success value, a
// sign_bracket with samples, is the input of every bracketing solver, so then(expand{}.on(w), brent{}) composes without
// re-evaluating the ends. It has no configurable stop criterion: it stops at the first sign change (intrinsic) or when
// the budget (60) runs out, failing with the best sample.
//
// Growth: geometric when the window is positive (lo / g, hi * g, g = 1.6), mirrored when it is negative, additive
// otherwise; only the end with the smaller |f| moves (one evaluation per step); on success the tightest sign-changing
// pair among the last samples is returned.
//
// Deferred to phase 3: a guess as start (with limits), the domain-edge backtrack on NaN, scan and subdivide.
#pragma once

#include <numerixx/roots/bracket.hpp>

#include <concepts>
#include <cstdint>
#include <expected>
#include <limits>
#include <optional>

NXX_BEGIN_HEADER

namespace nxx::roots
{
    template<real T>
    struct expand_state
    {
        T             lo, flo, hi, fhi;
        T             prev, fprev;    // where the end that moved last was before
        int           side;           // -1: lo moved last, +1: hi moved last, 0: nothing moved yet
        std::uint32_t nfev;
    };

    template<class Opt = nxx::options<never>>
    class expand : public search_facade
    {
        Opt opt_;

    public:
        static constexpr algo      id    = algos::expand;
        static constexpr view_kind views = view_kind::point;
        template<class In>
        static constexpr bool accepts_v = detail::window_input_v<In>;
        template<class F, class In>
        static constexpr bool callable_v = std::is_invocable_v<const F&, const detail::bracket_scalar_t<In>&>;

        constexpr expand()
            requires std::same_as<Opt, nxx::options<never>>
            : opt_ { never {}, max_iterations { 60 } }
        {}

        constexpr expand(nxx::detail::from_options_t, Opt o) : opt_(std::move(o)) {}

        constexpr const Opt& options() const noexcept { return opt_; }

        template<class O2>
        constexpr auto rebuild(O2 o) const
        { return expand<O2> { nxx::detail::from_options, std::move(o) }; }

        template<class F, class In>
            requires detail::window_input_v<In>
        constexpr auto prepare(const F& fn, const In& in) const
        {
            using T = detail::bracket_scalar_t<In>;
            using R = std::expected<problem<F, bracket<T>>, root_failure<F, T>>;
            auto b  = detail::to_window<F>(id, in);
            if (!b) return R { std::unexpect, std::move(b).error() };
            return R { problem<F, bracket<T>> { fn, *b, 0 } };
        }

        template<class F, real T>
        constexpr auto init(const problem<F, bracket<T>>& p) const -> std::expected<expand_state<T>, root_failure<F, T>>
        {
            using Fail = root_failure<F, T>;
            auto flo   = nxx::evaluate_sample(p.f, p.in.lo());
            if (!flo)
                return std::unexpected(Fail { flo.error().code, id, counters { 0, flo.error().evals }, std::nullopt, flo.error().cause });
            const std::uint32_t one = cost_of(p.f);
            auto                fhi = nxx::evaluate_sample(p.f, p.in.hi());
            if (!fhi)
                return std::unexpected(Fail { fhi.error().code,
                                              id,
                                              counters { 0, one + fhi.error().evals },
                                              root_estimate<T> { p.in.lo(), *flo, detail::unknown<T>(), std::nullopt },
                                              fhi.error().cause });
            return expand_state<T> { p.in.lo(), *flo, p.in.hi(), *fhi, p.in.lo(), *flo, 0, one + cost_of(p.f) };
        }

        template<class F, real T>
        constexpr auto step(const problem<F, bracket<T>>& p, const expand_state<T>& s) const
            -> std::expected<expand_state<T>, fault<callback_error_t<F, T>>>
        {
            using UE       = callback_error_t<F, T>;
            const T top    = (std::numeric_limits<T>::max)();
            const T bottom = std::numeric_limits<T>::lowest();
            // Endpoints stay finite: growth saturates at +-max (an infinite endpoint would make a sign_bracket that a
            // bracketing successor cannot split). w is the width; 2 * (hi/2 - lo/2) equals hi - lo whenever that is
            // finite and is inf otherwise, which the clamp turns into the largest finite step.
            const auto clamp       = [&](const T& v) { return v > top ? top : (v < bottom ? bottom : v); };
            const T    g           = T(8) / T(5);
            const T    w           = T(2) * (s.hi / T(2) - s.lo / T(2));
            const bool lo_at_limit = !(s.lo > bottom);
            const bool hi_at_limit = !(s.hi < top);
            if (lo_at_limit && hi_at_limit) return std::unexpected(fault<UE> { errc::stalled, 0, {} });    // nothing left to search
            bool move_lo = math::abs(s.flo) < math::abs(s.fhi);
            if (move_lo && lo_at_limit) move_lo = false;
            if (!move_lo && hi_at_limit) move_lo = true;
            expand_state<T> n = s;
            if (move_lo) {
                n.prev  = s.lo;
                n.fprev = s.flo;
                n.side  = -1;
                n.lo    = clamp(s.lo > T(0) ? s.lo / g : (s.hi < T(0) ? s.lo * g : s.lo - g * w));
                auto y  = nxx::evaluate_sample(p.f, n.lo);
                if (!y) return std::unexpected(y.error());
                n.flo = *y;
            }
            else {
                n.prev  = s.hi;
                n.fprev = s.fhi;
                n.side  = +1;
                n.hi    = clamp(s.hi < T(0) ? s.hi / g : (s.lo > T(0) ? s.hi * g : s.hi + g * w));
                auto y  = nxx::evaluate_sample(p.f, n.hi);
                if (!y) return std::unexpected(y.error());
                n.fhi = *y;
            }
            n.nfev += cost_of(p.f);
            return n;
        }

        template<real T>
        constexpr root_estimate<T> best(const expand_state<T>& s) const noexcept
        {
            return math::abs(s.flo) <= math::abs(s.fhi) ? root_estimate<T> { s.lo, s.flo, s.hi - s.lo, std::nullopt }
                                                        : root_estimate<T> { s.hi, s.fhi, s.hi - s.lo, std::nullopt };
        }

        template<real T>
        constexpr point_view<T> view(const expand_state<T>& s) const noexcept
        {
            const root_estimate<T> e = best(s);
            return point_view<T> { e.x, e.fx };
        }

        // Valid only on a sign change (the intrinsic stop), which is when the driver calls it; on any other state it
        // would forge a sign_bracket without a sign change (a contract violation, checked in debug builds).
        template<real T>
        constexpr sign_bracket<T> estimate(const expand_state<T>& s) const noexcept
        {
            NXX_EXPECTS(detail::opposite(s.flo, s.fhi));
            if (s.side < 0 && detail::opposite(s.flo, s.fprev))
                return sign_bracket<T> { nxx::detail::trust_me {}, s.lo, s.flo, s.prev, s.fprev };
            if (s.side > 0 && detail::opposite(s.fprev, s.fhi))
                return sign_bracket<T> { nxx::detail::trust_me {}, s.prev, s.fprev, s.hi, s.fhi };
            return sign_bracket<T> { nxx::detail::trust_me {}, s.lo, s.flo, s.hi, s.fhi };
        }

        template<real T>
        constexpr std::optional<stop_reason> intrinsic(const expand_state<T>& s) const noexcept
        {
            if (detail::opposite(s.flo, s.fhi)) return stop_reason::algorithm;
            return std::nullopt;
        }
    };

    expand() -> expand<>;
}    // namespace nxx::roots

NXX_END_HEADER

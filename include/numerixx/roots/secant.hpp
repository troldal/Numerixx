// The secant method (DESIGN §7.2): derivative-free. The second point is x0 + 2^-10 max(|x0|, 1); a flat secant is
// errc::stalled; a non-finite step is errc::diverged. With .with_projection(p), every proposed iterate is projected
// before it is evaluated, and an iterate pinned at the edge fails with errc::stalled. Default criterion
// step_tol<7, 10>, budget 50.
//
// Deferred to phase 3: the progress window (stall and divergence over a window of steps) and the step-length cap.
#pragma once

#include <numerixx/roots/bracket.hpp>

#include <concepts>
#include <cstdint>
#include <expected>
#include <optional>

NXX_BEGIN_HEADER

namespace nxx::roots
{
    template<real T>
    struct secant_state
    {
        T             x0, f0, x1, f1;
        T             step;    // length of the last proposed step, before projection: what the stop criteria see
        std::uint32_t nfev;
    };

    template<class Opt = nxx::options<step_tol<7, 10>>>
    class secant : public open_facade
    {
        Opt opt_;

        using stop_type = typename Opt::stop_type;

    public:
        static constexpr algo      id       = algos::secant;
        static constexpr view_kind views    = view_kind::point;
        static constexpr bool      projects = true;
        template<class In>
        static constexpr bool accepts_v = detail::open_input_v<In>;
        template<class F, class In>
        static constexpr bool callable_v = std::is_invocable_v<const F&, const detail::open_scalar_t<In>&>;

        constexpr secant()
            requires std::same_as<Opt, nxx::options<step_tol<7, 10>>>
            : opt_ { step_tol<7, 10> {}, max_iterations { 50 } }
        {}

        constexpr explicit secant(stop_type stop)
            requires(criterion_for_v<stop_type, view_kind::point> && std::same_as<Opt, nxx::options<stop_type>>)
            : opt_ { stop, max_iterations { 50 } }
        {}

        template<class C>
            requires(is_criterion_v<C> && !criterion_for_v<C, view_kind::point>)
        explicit secant(C) NXX_DELETE("width_tol needs a bracketing method (the view has no enclosure()); use x_tol or step_tol");

        constexpr secant(nxx::detail::from_options_t, Opt o) : opt_(std::move(o)) {}

        constexpr const Opt& options() const noexcept { return opt_; }

        template<class O2>
        constexpr auto rebuild(O2 o) const
        { return secant<O2> { nxx::detail::from_options, std::move(o) }; }

        template<class F, class In>
            requires detail::open_input_v<In>
        constexpr auto prepare(const F& fn, const In& in) const
        { return detail::prepare_open(id, fn, in); }

        template<class F, real T>
        constexpr auto init(const problem<F, detail::open_start<T>>& p) const -> std::expected<secant_state<T>, root_failure<F, T>>
        {
            using Fail       = root_failure<F, T>;
            const T       x0 = opt_.project(p.in.x0);
            T             f0 {};
            std::uint32_t used = 0;
            if (p.in.fx0 && x0 == p.in.x0)
                f0 = *p.in.fx0;    // seeded by a previous stage: no re-evaluation
            else {
                auto y = nxx::evaluate(p.f, x0);
                if (!y) return std::unexpected(Fail { y.error().code, id, counters { 0, y.error().evals }, std::nullopt, y.error().cause });
                f0   = *y;
                used = cost_of(p.f);
            }
            if (f0 == T(0)) return secant_state<T> { x0, f0, x0, f0, T(0), used };    // exact zero: intrinsic stop at once

            // The second point: x0 + h, or x0 - h where that overflows or the projection pins x0 + h to x0.
            const T h    = math::pow2<T>(-10) * (std::max)(math::abs(x0), T(1));
            const T up   = x0 + h;
            const T down = x0 - h;
            T       x1   = math::isfinite(up) ? opt_.project(up) : x0;
            if (x1 == x0 && math::isfinite(down)) x1 = opt_.project(down);
            const root_estimate<T> first { x0, f0, detail::unknown<T>(), std::nullopt };
            if (x1 == x0) return std::unexpected(Fail { errc::stalled, id, counters { 0, used }, first, {} });
            auto y1 = nxx::evaluate(p.f, x1);
            if (!y1) return std::unexpected(Fail { y1.error().code, id, counters { 0, used + y1.error().evals }, first, y1.error().cause });
            return secant_state<T> { x0, f0, x1, *y1, math::abs(x1 - x0), used + cost_of(p.f) };
        }

        template<class F, real T>
        constexpr auto step(const problem<F, detail::open_start<T>>& p, const secant_state<T>& s) const
            -> std::expected<secant_state<T>, fault<callback_error_t<F, T>>>
        {
            using UE = callback_error_t<F, T>;
            if (s.f1 == s.f0) return std::unexpected(fault<UE> { errc::stalled, 0, {} });    // flat secant
            const T dx = s.x1 - s.x0;
            const T df = s.f1 - s.f0;
            T       delta;
            if (math::isfinite(df)) {
                delta = s.f1 * dx / df;
                if (!math::isfinite(delta)) delta = dx / df * s.f1;    // the product overflowed (roots near 1e300)
            }
            else    // f1 - f0 overflowed: the ratio of the halves is the same and representable
                delta = dx * ((s.f1 / T(2)) / (s.f1 / T(2) - s.f0 / T(2)));
            const T proposed = s.x1 - delta;
            if (!math::isfinite(proposed)) return std::unexpected(fault<UE> { errc::diverged, 0, {} });
            const T x2 = opt_.project(proposed);
            if (x2 != proposed && x2 == s.x1) return std::unexpected(fault<UE> { errc::stalled, 0, {} });    // pinned at the edge
            auto f2 = nxx::evaluate(p.f, x2);
            if (!f2) return std::unexpected(f2.error());
            return secant_state<T> { s.x1, s.f1, x2, *f2, math::abs(proposed - s.x1), s.nfev + cost_of(p.f) };
        }

        template<real T>
        constexpr point_view<T> view(const secant_state<T>& s) const noexcept
        { return point_view<T> { s.x1, s.f1, s.step }; }

        template<real T>
        constexpr root_estimate<T> estimate(const secant_state<T>& s) const noexcept
        { return root_estimate<T> { s.x1, s.f1, math::abs(s.x1 - s.x0), std::nullopt }; }

        template<real T>
        constexpr root_estimate<T> best(const secant_state<T>& s) const noexcept
        {
            if (math::abs(s.f0) < math::abs(s.f1)) return root_estimate<T> { s.x0, s.f0, math::abs(s.x1 - s.x0), std::nullopt };
            return estimate(s);
        }

        template<real T>
        constexpr std::optional<stop_reason> intrinsic(const secant_state<T>& s) const noexcept
        {
            if (s.f1 == T(0)) return stop_reason::exact_zero;
            return std::nullopt;
        }
    };

    template<class C>
        requires is_criterion_v<C>
    secant(C) -> secant<nxx::options<C>>;
}    // namespace nxx::roots

NXX_END_HEADER

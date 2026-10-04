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
        T             step;    // the larger of the proposed and the projected step: what the stop criteria see
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
            requires(stop_criterion_for_v<stop_type, view_kind::point> && std::same_as<Opt, nxx::options<stop_type>>)
            : opt_ { stop, max_iterations { 50 } }
        {}

        template<class C>
            requires(criterion_for_v<C, view_kind::point> && nxx::detail::guard_only_v<C>)
        explicit secant(C) NXX_DELETE("min_iterations only guards another criterion: combine it with a convergence test "
                                      "using && (your_test && min_iterations{n})");

        template<class C>
            requires(is_criterion_v<C> && !criterion_for_v<C, view_kind::point>)
        explicit secant(C) NXX_DELETE("width_tol needs a bracketing method (the view has no enclosure()); use x_tol or step_tol");

        constexpr secant(nxx::detail::from_options_t, Opt o) : opt_(std::move(o)) {}

        constexpr const Opt& options() const noexcept { return opt_; }

        // Only options whose stop criterion can stop this solver: rebuild is public, so it must not be a way around the
        // constructors and with_stop (a bare min_iterations guard would report success without testing accuracy).
        template<class O2>
            requires nxx::detail::stop_allowed_v<secant, typename O2::stop_type>
        constexpr auto rebuild(O2 o) const
        { return secant<O2> { nxx::detail::from_options, std::move(o) }; }

        template<class F, class In>
            requires detail::open_input_v<In>
        constexpr auto prepare(const F& fn, const In& in) const
        { return detail::prepare_open(id, fn, in); }

        template<class F, real T>
        constexpr auto init(const problem<F, detail::open_start<T>>& p) const -> std::expected<secant_state<T>, root_failure<F, T>>
        {
            using Fail = root_failure<F, T>;
            const T x0 = opt_.project(p.in.x0);
            // A projection can leave the reals (clamp_to{inf, inf}, or a custom one): f is never evaluated at a non-finite
            // point, so no exact zero or criterion can be reported there (DESIGN §7.2, before each evaluation).
            if (!math::isfinite(x0)) return std::unexpected(Fail { errc::non_finite_input, id, counters {}, std::nullopt, {} });
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
            // An exact zero: the intrinsic test stops at once. No step was taken, so its length is unknown (inf). A real
            // step has a finite length unless |proposed - x1| overflows, and then inf is still an honest uncertainty.
            if (f0 == T(0)) return secant_state<T> { x0, f0, x0, f0, detail::unknown<T>(), used };

            // The second point: x0 + h, or x0 - h where that overflows or the projection pins x0 + h to x0 or sends it off
            // the reals (a domain whose outside is marked with NaN or inf).
            const T h    = math::pow2<T>(-10) * (std::max)(math::abs(x0), T(1));
            const T up   = x0 + h;
            const T down = x0 - h;
            const T a    = math::isfinite(up) ? opt_.project(up) : x0;    // x0 where x0 + h overflows
            T       x1   = a;
            if ((x1 == x0 || !math::isfinite(x1)) && math::isfinite(down)) x1 = opt_.project(down);
            const root_estimate<T> first { x0, f0, detail::unknown<T>(), std::nullopt };
            // Neither neighbour is usable: diverged if the projection sent one of them off the reals, stalled if it only
            // pinned them to x0, whichever side it was (DESIGN §7.2).
            if (x1 == x0 || !math::isfinite(x1)) {
                const errc why = math::isfinite(a) && math::isfinite(x1) ? errc::stalled : errc::diverged;
                return std::unexpected(Fail { why, id, counters { 0, used }, first, {} });
            }
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
            if (!math::isfinite(x2)) return std::unexpected(fault<UE> { errc::diverged, 0, {} });            // projected off the reals
            if (x2 != proposed && x2 == s.x1) return std::unexpected(fault<UE> { errc::stalled, 0, {} });    // pinned at the edge
            auto f2 = nxx::evaluate(p.f, x2);
            if (!f2) return std::unexpected(f2.error());
            // The criteria see the larger of the proposed and the actual step: a projection that moves the point further
            // (to a far finite value) must not look like convergence, and one that pins it must not either.
            const T moved = math::abs(x2 - s.x1);
            return secant_state<T> { s.x1, s.f1, x2, *f2, (std::max)(math::abs(proposed - s.x1), moved), s.nfev + cost_of(p.f) };
        }

        template<real T>
        constexpr point_view<T> view(const secant_state<T>& s) const noexcept
        { return point_view<T> { s.x1, s.f1, s.step }; }

        template<real T>
        constexpr root_estimate<T> estimate(const secant_state<T>& s) const noexcept
        {
            // |x1 - x0|, or unknown when the solve stopped at its start (an exact zero at x0): no step was taken. x1 == x0
            // alone does not say that, because a real step can round to 0.
            const T uncertainty = math::isfinite(s.step) ? math::abs(s.x1 - s.x0) : detail::unknown<T>();
            return root_estimate<T> { s.x1, s.f1, uncertainty, std::nullopt };
        }

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

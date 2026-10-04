// Newton's method (DESIGN §7.2, D13) with an explicit derivative source:
//   (a) a callable df:                          newton{}.with_derivative(df)
//   (b) a derivative policy with bind(f):       newton{}.with_derivative(deriv::numeric{}) (works in curried chains)
//   (c) a structural .derivative() on f itself: a polynomial or a spline
// Recognised structurally, so roots does not include deriv. Without a source (ready_v false), calling it does not
// compile, with the reason in open_facade. A zero derivative is errc::zero_derivative; a non-finite step is
// errc::diverged; the projection works as for secant. The error type of the solve is the common cause of f and f'
// (DESIGN §6.4). Default criterion step_tol<3, 5>, budget 30.
//
// Deferred to phase 3: the progress window and the step-length cap.
#pragma once

#include <numerixx/roots/bracket.hpp>

#include <concepts>
#include <cstdint>
#include <expected>
#include <optional>
#include <type_traits>

NXX_BEGIN_HEADER

namespace nxx::roots
{
    namespace detail
    {
        template<class D, class F>
        inline constexpr bool derivative_available_v = !std::is_same_v<D, no_derivative> || has_derivative_source_v<F>;

        template<class D, class F>
        constexpr auto bind_derivative(const D& d, const F& fn)
        {
            if constexpr (requires { d.bind(fn); })
                return d.bind(fn);    // a policy: a function of f
            else if constexpr (std::is_same_v<D, no_derivative>)
                return nxx::derivative_source(fn);    // f's own .derivative()
            else
                return d;    // a callable
        }

        // The derivative is boxed, so a problem (and a steps_view over it) stays copy-assignable when df is a lambda with
        // captures.
        template<real T, class DF>
        struct newton_start
        {
            T                             x0;
            std::optional<T>              fx0;
            nxx::detail::copyable_box<DF> df;
        };
    }    // namespace detail

    template<real T>
    struct newton_state
    {
        T             x, fx, dx;
        T             step;    // the larger of the proposed and the projected step: what the stop criteria see
        std::uint32_t nfev;
    };

    template<class Opt = nxx::options<step_tol<3, 5>>>
    class newton : public open_facade
    {
        Opt opt_;

        using stop_type = typename Opt::stop_type;
        using D         = typename Opt::derivative_type;

        template<class F, real T>
        using df_t = decltype(detail::bind_derivative(std::declval<const D&>(), std::declval<const F&>()));
        template<class F, real T>
        using cause_t = common_cause_t<callback_error_t<F, T>, callback_error_t<df_t<F, T>, T>>;

    public:
        static constexpr algo      id              = algos::newton;
        static constexpr view_kind views           = view_kind::point;
        static constexpr bool      uses_derivative = true;
        static constexpr bool      projects        = true;
        template<class In>
        static constexpr bool accepts_v = detail::open_input_v<In>;
        template<class F, class In>
        static constexpr bool callable_v = std::is_invocable_v<const F&, const detail::open_scalar_t<In>&>;
        template<class F>
        static constexpr bool ready_v = detail::derivative_available_v<D, F>;

        constexpr newton()
            requires std::same_as<Opt, nxx::options<step_tol<3, 5>>>
            : opt_ { step_tol<3, 5> {}, max_iterations { 30 } }
        {}

        constexpr explicit newton(stop_type stop)
            requires(stop_criterion_for_v<stop_type, view_kind::point> && std::same_as<Opt, nxx::options<stop_type>>)
            : opt_ { stop, max_iterations { 30 } }
        {}

        template<class C>
            requires(criterion_for_v<C, view_kind::point> && nxx::detail::guard_only_v<C>)
        explicit newton(C) NXX_DELETE("min_iterations only guards another criterion: combine it with a convergence test "
                                      "using && (your_test && min_iterations{n})");

        template<class C>
            requires(is_criterion_v<C> && !criterion_for_v<C, view_kind::point>)
        explicit newton(C) NXX_DELETE("width_tol needs a bracketing method (the view has no enclosure()); use x_tol or step_tol");

        constexpr newton(nxx::detail::from_options_t, Opt o) : opt_(std::move(o)) {}

        constexpr const Opt& options() const noexcept { return opt_; }

        // Only options whose stop criterion can stop this solver: rebuild is public, so it must not be a way around the
        // constructors and with_stop (a bare min_iterations guard would report success without testing accuracy).
        template<class O2>
            requires nxx::detail::stop_allowed_v<newton, typename O2::stop_type>
        constexpr auto rebuild(O2 o) const
        { return newton<O2> { nxx::detail::from_options, std::move(o) }; }

        template<class F, class In>
            requires(detail::open_input_v<In> && detail::derivative_available_v<D, F>)
        constexpr auto prepare(const F& fn, const In& in) const
        {
            using T    = detail::open_scalar_t<In>;
            using DF   = df_t<F, T>;
            using UE   = cause_t<F, T>;
            using Fail = failure<root_estimate<T>, UE>;
            using R    = std::expected<problem<F, detail::newton_start<T, DF>>, Fail>;
            auto p0    = detail::prepare_open(id, fn, in);
            if (!p0) return R { std::unexpect, detail::widen<UE>(p0.error()) };
            return R { problem<F, detail::newton_start<T, DF>> {
                fn,
                detail::newton_start<T, DF> { p0->in.x0,
                                              p0->in.fx0,
                                              nxx::detail::copyable_box<DF> { detail::bind_derivative(*opt_.derivative, fn) } },
                0 } };
        }

        template<class F, real T, class DF>
        constexpr auto init(const problem<F, detail::newton_start<T, DF>>& p) const
            -> std::expected<newton_state<T>, failure<root_estimate<T>, cause_t<F, T>>>
        {
            using UE   = cause_t<F, T>;
            using Fail = failure<root_estimate<T>, UE>;
            const T x0 = opt_.project(p.in.x0);
            // As in secant: a projection that leaves the reals is rejected before f is evaluated there (DESIGN §7.2).
            if (!math::isfinite(x0)) return std::unexpected(Fail { errc::non_finite_input, id, counters {}, std::nullopt, {} });
            if (p.in.fx0 && x0 == p.in.x0)
                return newton_state<T> { x0, *p.in.fx0, detail::unknown<T>(), detail::unknown<T>(), 0 };    // seeded
            auto y = nxx::evaluate(p.f, x0);
            if (!y) {
                const fault<UE> e = detail::widen<UE>(y.error());
                return std::unexpected(Fail { e.code, id, counters { 0, e.evals }, std::nullopt, e.cause });
            }
            return newton_state<T> { x0, *y, detail::unknown<T>(), detail::unknown<T>(), cost_of(p.f) };
        }

        template<class F, real T, class DF>
        constexpr auto step(const problem<F, detail::newton_start<T, DF>>& p, const newton_state<T>& s) const
            -> std::expected<newton_state<T>, fault<cause_t<F, T>>>
        {
            using UE              = cause_t<F, T>;
            const std::uint32_t c = cost_of(*p.in.df);
            auto                d = nxx::evaluate(*p.in.df, s.x);
            if (!d) return std::unexpected(detail::widen<UE>(d.error()));
            if (*d == T(0)) return std::unexpected(fault<UE> { errc::zero_derivative, c, {} });
            const T proposed = s.x - s.fx / *d;
            if (!math::isfinite(proposed)) return std::unexpected(fault<UE> { errc::diverged, c, {} });
            const T x1 = opt_.project(proposed);
            if (!math::isfinite(x1)) return std::unexpected(fault<UE> { errc::diverged, c, {} });           // projected off the reals
            if (x1 != proposed && x1 == s.x) return std::unexpected(fault<UE> { errc::stalled, c, {} });    // pinned at the edge
            auto y = nxx::evaluate(p.f, x1);
            if (!y) {
                fault<UE> e = detail::widen<UE>(y.error());
                e.evals += c;
                return std::unexpected(e);
            }
            // The criteria see the larger of the proposed and the actual step: a projection that moves the point further
            // (to a far finite value) must not look like convergence, and one that pins it must not either.
            const T moved = math::abs(x1 - s.x);
            return newton_state<T> { x1, *y, moved, (std::max)(math::abs(proposed - s.x), moved), s.nfev + c + cost_of(p.f) };
        }

        template<real T>
        constexpr point_view<T> view(const newton_state<T>& s) const noexcept
        { return point_view<T> { s.x, s.fx, s.step }; }

        template<real T>
        constexpr root_estimate<T> estimate(const newton_state<T>& s) const noexcept
        { return root_estimate<T> { s.x, s.fx, s.dx, std::nullopt }; }

        template<real T>
        constexpr root_estimate<T> best(const newton_state<T>& s) const noexcept
        { return estimate(s); }

        template<real T>
        constexpr std::optional<stop_reason> intrinsic(const newton_state<T>& s) const noexcept
        {
            if (s.fx == T(0)) return stop_reason::exact_zero;
            return std::nullopt;
        }
    };

    template<class C>
        requires is_criterion_v<C>
    newton(C) -> newton<nxx::options<C>>;
}    // namespace nxx::roots

NXX_END_HEADER

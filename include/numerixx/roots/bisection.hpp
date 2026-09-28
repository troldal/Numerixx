// Bisection (DESIGN §7.2): one evaluation per step, overflow-safe midpoint, exact zero and unsplittable brackets as
// intrinsic stops, and the pole check. Default criterion floored_width{}, budget 200: floored_width{} is achievable in
// every T (DESIGN §3.5), and cpp_bin_float_50 (168 bits) needs about 165 halvings from a unit bracket.
//
// Deferred to phase 3: the representation-space midpoint for built-in floats (a root at 1e-200 needs about 665
// value-space halvings and exhausts the default budget).
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
    struct bisection_state
    {
        sign_bracket<T> b;
        std::uint32_t   nfev;
    };

    template<class Opt = nxx::options<floored_width>>
    class bisection : public bracketing_facade
    {
        Opt opt_;

        using stop_type = typename Opt::stop_type;

    public:
        static constexpr algo      id    = algos::bisection;
        static constexpr view_kind views = view_kind::enclosure;
        template<class In>
        static constexpr bool accepts_v = detail::bracket_input_v<In>;
        template<class F, class In>
        static constexpr bool callable_v = std::is_invocable_v<const F&, const detail::bracket_scalar_t<In>&>;

        constexpr bisection()
            requires std::same_as<Opt, nxx::options<floored_width>>
            : opt_ { floored_width {}, max_iterations { 200 } }
        {}

        constexpr explicit bisection(stop_type stop)
            requires(criterion_for_v<stop_type, view_kind::enclosure> && std::same_as<Opt, nxx::options<stop_type>>)
            : opt_ { stop, max_iterations { 200 } }
        {}

        template<class C>
            requires(is_criterion_v<C> && !criterion_for_v<C, view_kind::enclosure>)
        explicit bisection(C) NXX_DELETE("x_tol and step_tol compare successive iterates; bracketing methods converge on the "
                                         "enclosure: use width_tol{abs[, rel]} or floored_width{}");

        constexpr bisection(nxx::detail::from_options_t, Opt o) : opt_(std::move(o)) {}

        constexpr const Opt& options() const noexcept { return opt_; }

        template<class O2>
        constexpr auto rebuild(O2 o) const
        { return bisection<O2> { nxx::detail::from_options, std::move(o) }; }

        template<class F, class In>
            requires detail::bracket_input_v<In>
        constexpr auto prepare(const F& fn, const In& in) const
        { return detail::prepare_bracketing(id, fn, in); }

        template<class F, real T>
        constexpr auto init(const problem<F, sign_bracket<T>>& p) const -> std::expected<bisection_state<T>, root_failure<F, T>>
        { return bisection_state<T> { p.in, p.nfev0 }; }

        template<class F, real T>
        constexpr auto step(const problem<F, sign_bracket<T>>& p, const bisection_state<T>& s) const
            -> std::expected<bisection_state<T>, fault<callback_error_t<F, T>>>
        {
            const T m  = math::midpoint(s.b.lo(), s.b.hi());
            auto    fm = nxx::evaluate_sample(p.f, m);
            if (!fm) return std::unexpected(fm.error());
            return bisection_state<T> { s.b.narrowed(m, *fm), s.nfev + cost_of(p.f) };
        }

        template<real T>
        constexpr enclosure_view<T> view(const bisection_state<T>& s) const noexcept
        {
            const root_estimate<T> e = s.b.best();
            return enclosure_view<T> { e.x, e.fx, s.b.lo(), s.b.hi() };
        }

        template<real T>
        constexpr root_estimate<T> estimate(const bisection_state<T>& s) const noexcept
        { return s.b.best(); }

        template<real T>
        constexpr root_estimate<T> best(const bisection_state<T>& s) const noexcept
        { return s.b.best(); }

        template<real T>
        constexpr std::optional<stop_reason> intrinsic(const bisection_state<T>& s) const noexcept
        {
            if (s.b.has_exact_zero()) return stop_reason::exact_zero;
            const T m = math::midpoint(s.b.lo(), s.b.hi());
            if (!(s.b.lo() < m && m < s.b.hi())) return stop_reason::resolution_limit;    // unsplittable
            return std::nullopt;
        }

        template<class F, real T>
        constexpr auto finish(const problem<F, sign_bracket<T>>& p, const solution<root_estimate<T>>& sol) const
        { return detail::pole_check<root_failure<F, T>>(p.in, sol); }
    };

    template<class C>
        requires is_criterion_v<C>
    bisection(C) -> bisection<nxx::options<C>>;
}    // namespace nxx::roots

NXX_END_HEADER

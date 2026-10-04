// Brent's method, zeroin (DESIGN §7.2, §6.8): inverse quadratic interpolation and secant steps safeguarded by
// bisection. Its tolerance is a width criterion, width_tol{abs[, rel]} or floored_width{} (the default):
// tol1 = max(threshold / 2, 2 eps |b|), where threshold is the tolerance's enclosure form (width_tol: abs + rel
// min(|lo|, |hi|)) and 2 eps |b| is Brent's resolution floor (a smaller step would not move b). Its intrinsic stop
// (|c - b| / 2 <= tol1) therefore meets the threshold whenever the threshold is attainable, and reports
// stop_reason::criterion exactly then, as every width criterion guarantees (DESIGN §9.3). With a tolerance below the
// floor it stops at the floor, with width <= 4 eps |b|: stop_reason::criterion if that width still meets the threshold,
// stop_reason::resolution_limit otherwise. The external stop defaults to never{}, so there are not two sources of
// truth. The intrinsic test runs before the external stop, so with_stop adds an early exit (f_tol) or a failure
// (max_evaluations) and cannot tighten the tolerance: with_stop and rebuild reject a stop criterion that contains a
// width criterion, which goes to the constructor instead (brent{width_tol{1e-12}}). It bisects while a sample is
// infinite (log(x) on [0, 2]) and runs the pole check. Budget 100.
#pragma once

#include <numerixx/roots/bracket.hpp>

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <expected>
#include <limits>
#include <optional>
#include <type_traits>

NXX_BEGIN_HEADER

namespace nxx::roots
{
    namespace detail
    {
        // A tolerance for Brent: a width criterion with a threshold for an enclosure. False, not a hard error, for any
        // other W: each test is asked only once the one before it holds, because an && in a variable template's
        // initializer does not stop the instantiation of its later operands (W::applies_to was formed for a double).
        template<class W>
        inline constexpr bool width_tolerance_v = [] {
            if constexpr (criterion_for_v<W, view_kind::enclosure>)
                return requires(const W& w) { w.threshold(1.0, 2.0); } && std::remove_cvref_t<W>::applies_to == view_kind::enclosure;
            else
                return false;
        }();
    }    // namespace detail

    template<real T>
    struct brent_state
    {
        T             a, fa, b, fb, c, fc, d, e;
        std::uint32_t nfev;
    };

    template<class Tol = floored_width, class Opt = nxx::options<never>>
    class brent : public bracketing_facade
    {
        Tol tol_;
        Opt opt_;

        // The tolerance for the current enclosure [min(b, c), max(b, c)].
        template<real T>
        constexpr T threshold(const brent_state<T>& s) const noexcept
        { return tol_.threshold((std::min)(s.b, s.c), (std::max)(s.b, s.c)); }

        template<real T>
        constexpr T tol1(const brent_state<T>& s) const noexcept
        { return (std::max)(threshold(s) / T(2), T(2) * std::numeric_limits<T>::epsilon() * math::abs(s.b)); }

        // (c - b) / 2, overflow-safe: on [-1.7e308, 1.7e308] the plain difference is inf. The plain form stays the default,
        // so ordinary brackets round exactly as before.
        template<real T>
        static constexpr T half_step(const T& c, const T& b) noexcept
        {
            const T h = (c - b) / T(2);
            return math::isfinite(h) ? h : c / T(2) - b / T(2);
        }

        template<real T>
        static constexpr void normalise(brent_state<T>& s) noexcept
        {
            if (!detail::opposite(s.fb, s.fc)) {
                s.c  = s.a;
                s.fc = s.fa;
                s.d = s.e = s.b - s.a;
            }
            if (math::abs(s.fc) < math::abs(s.fb)) {
                s.a  = s.b;
                s.b  = s.c;
                s.c  = s.a;
                s.fa = s.fb;
                s.fb = s.fc;
                s.fc = s.fa;
            }
        }

    public:
        static constexpr algo      id    = algos::brent;
        static constexpr view_kind views = view_kind::enclosure;
        // The tolerance decides convergence, so with_stop and rebuild reject a width criterion (solver_facade).
        static constexpr bool internal_tolerance = true;
        template<class In>
        static constexpr bool accepts_v = detail::bracket_input_v<In>;
        template<class F, class In>
        static constexpr bool callable_v = std::is_invocable_v<const F&, const detail::bracket_scalar_t<In>&>;

        constexpr brent()
            requires(std::same_as<Tol, floored_width> && std::same_as<Opt, nxx::options<never>>)
            : tol_ {},
              opt_ { never {}, max_iterations { 100 } }
        {}

        constexpr explicit brent(Tol tol)
            requires(detail::width_tolerance_v<Tol> && std::same_as<Opt, nxx::options<never>>)
            : tol_(tol),
              opt_ { never {}, max_iterations { 100 } }
        {}

        template<class C>
            requires(is_criterion_v<C> && !detail::width_tolerance_v<C>)
        explicit brent(C) NXX_DELETE("brent's tolerance is a width criterion: width_tol{abs[, rel]} or floored_width{} "
                                     "(x_tol and step_tol compare successive iterates; bracketing methods converge on the enclosure)");

        template<class R>
            requires((std::is_arithmetic_v<R> || real<R>) && !std::is_convertible_v<R, Tol>)
        explicit brent(R) NXX_DELETE("a tolerance is a criterion, not a number: write brent{nxx::width_tol{1e-10}}");

        // Constrained as rebuild is, so no options that carry a width criterion build a brent, not even through the
        // detail key (DESIGN §6.8). The constraint spells brent<Tol, Opt>, not the injected-class-name: in the deduction
        // guide that cl builds from this constructor, a constraint on the injected name is never satisfied.
        constexpr brent(nxx::detail::from_options_t, Opt o, Tol tol)
            requires nxx::detail::stop_allowed_v<brent<Tol, Opt>, typename Opt::stop_type>
            : tol_(tol),
              opt_(std::move(o))
        {}

        constexpr const Opt& options() const noexcept { return opt_; }
        constexpr const Tol& tolerance() const noexcept { return tol_; }

        // Only options whose stop criterion can stop this solver: rebuild is public, so it must not be a way around the
        // constructors and with_stop (a bare min_iterations guard would report success without testing accuracy), nor
        // take a width criterion, which the tolerance's intrinsic test would preempt (with_stop's rule, DESIGN §6.8).
        template<class O2>
            requires nxx::detail::stop_allowed_v<brent, typename O2::stop_type>
        constexpr auto rebuild(O2 o) const
        { return brent<Tol, O2> { nxx::detail::from_options, std::move(o), tol_ }; }

        template<class F, class In>
            requires detail::bracket_input_v<In>
        constexpr auto prepare(const F& fn, const In& in) const
        { return detail::prepare_bracketing(id, fn, in); }

        template<class F, real T>
        constexpr auto init(const problem<F, sign_bracket<T>>& p) const -> std::expected<brent_state<T>, root_failure<F, T>>
        {
            brent_state<T> s { p.in.lo(), p.in.flo(), p.in.hi(), p.in.fhi(), p.in.lo(), p.in.flo(), T(0), T(0), p.nfev0 };
            s.d = s.e = s.b - s.a;
            normalise(s);
            return s;
        }

        template<class F, real T>
        constexpr auto step(const problem<F, sign_bracket<T>>& p, const brent_state<T>& s0) const
            -> std::expected<brent_state<T>, fault<callback_error_t<F, T>>>
        {
            brent_state<T> s              = s0;    // a local copy: the state value itself is never mutated
            const T        t1             = tol1(s);
            const T        xm             = half_step(s.c, s.b);
            const bool     finite_samples = math::isfinite(s.fa) && math::isfinite(s.fb) && math::isfinite(s.fc);
            if (finite_samples && math::isfinite(s.e) && math::abs(s.e) >= t1 && math::abs(s.fa) > math::abs(s.fb)) {
                T       pp, q, r;
                const T sr = s.fb / s.fa;
                if (s.a == s.c) {    // secant
                    pp = T(2) * xm * sr;
                    q  = T(1) - sr;
                }
                else {    // inverse quadratic interpolation
                    q  = s.fa / s.fc;
                    r  = s.fb / s.fc;
                    pp = sr * (T(2) * xm * q * (q - r) - (s.b - s.a) * (r - T(1)));
                    q  = (q - T(1)) * (r - T(1)) * (sr - T(1));
                }
                if (pp > T(0)) q = -q;
                pp         = math::abs(pp);
                const T m1 = T(3) * xm * q - math::abs(t1 * q);
                const T m2 = math::abs(s.e * q);
                if (T(2) * pp < (std::min)(m1, m2)) {
                    s.e = s.d;
                    s.d = pp / q;
                }
                else {
                    s.d = xm;
                    s.e = s.d;
                }
            }
            else {    // bisection, also while a sample is infinite
                s.d = xm;
                s.e = s.d;
            }
            s.a     = s.b;
            s.fa    = s.fb;
            s.b     = math::abs(s.d) > t1 ? s.b + s.d : s.b + (xm > T(0) ? t1 : -t1);
            auto fb = nxx::evaluate_sample(p.f, s.b);
            if (!fb) return std::unexpected(fb.error());
            s.fb = *fb;
            s.nfev += cost_of(p.f);
            normalise(s);
            return s;
        }

        template<real T>
        constexpr enclosure_view<T> view(const brent_state<T>& s) const noexcept
        { return enclosure_view<T> { s.b, s.fb, (std::min)(s.b, s.c), (std::max)(s.b, s.c) }; }

        template<real T>
        constexpr root_estimate<T> estimate(const brent_state<T>& s) const noexcept
        {
            if (s.b == s.c) return root_estimate<T> { s.b, s.fb, T(0), std::nullopt };
            const sign_bracket<T> enc = s.b < s.c ? sign_bracket<T> { nxx::detail::trust_me {}, s.b, s.fb, s.c, s.fc }
                                                  : sign_bracket<T> { nxx::detail::trust_me {}, s.c, s.fc, s.b, s.fb };
            return root_estimate<T> { s.b, s.fb, enc.width(), enc };
        }

        template<real T>
        constexpr root_estimate<T> best(const brent_state<T>& s) const noexcept
        { return estimate(s); }

        template<real T>
        constexpr std::optional<stop_reason> intrinsic(const brent_state<T>& s) const noexcept
        {
            if (s.fb == T(0)) return stop_reason::exact_zero;
            if (math::abs(half_step(s.c, s.b)) <= tol1(s))    // Brent's own test: the next step would be below resolution
                return math::abs(s.c - s.b) <= threshold(s) ? stop_reason::criterion : stop_reason::resolution_limit;
            return std::nullopt;
        }

        template<class F, real T>
        constexpr auto finish(const problem<F, sign_bracket<T>>& p, const solution<root_estimate<T>>& sol) const
        { return detail::pole_check<root_failure<F, T>>(p.in, sol); }
    };

    template<class C>
        requires is_criterion_v<C>
    brent(C) -> brent<C>;
}    // namespace nxx::roots

NXX_END_HEADER

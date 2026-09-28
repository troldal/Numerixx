// What every 1-D root solver shares (DESIGN §6.3, §6.5, §6.6, §7.2): the estimate type, the sign bracket, the views
// the stop criteria see, the accepted inputs and their validation, the per-iterate projection clamp_to, and the pole
// check of the bracketing methods.
#pragma once

#include <numerixx/core.hpp>

#include <algorithm>
#include <cstdint>
#include <expected>
#include <limits>
#include <optional>
#include <type_traits>
#include <utility>

NXX_BEGIN_HEADER

namespace nxx::roots
{
    // Algorithm ids (open enum, roots 1-39).
    namespace algos
    {
        inline constexpr algo bisection { 1 };
        inline constexpr algo brent { 4 };
        inline constexpr algo secant { 9 };
        inline constexpr algo newton { 11 };
        inline constexpr algo expand { 13 };
    }    // namespace algos

    namespace detail
    {
        // Opposite signs, or a zero at either end: comparisons, not products (no overflow, no underflow).
        template<real T>
        constexpr bool opposite(const T& a, const T& b) noexcept
        { return (a <= T(0) && b >= T(0)) || (a >= T(0) && b <= T(0)); }

        // The uncertainty of an estimate that has neither an enclosure nor a step.
        template<real T>
        constexpr T unknown() noexcept
        { return std::numeric_limits<T>::infinity(); }
    }    // namespace detail

    template<real T>
    struct root_estimate;

    // A bracket with endpoint samples of opposite sign (or an exact zero). Samples may be +-inf. The only transition is
    // narrowed(m, fm), which keeps the half that still changes sign.
    template<real T>
    class sign_bracket
    {
        T lo_;
        T flo_;
        T hi_;
        T fhi_;

    public:
        using value_type = T;

        constexpr sign_bracket(nxx::detail::trust_me, T lo, T flo, T hi, T fhi) noexcept : lo_(lo), flo_(flo), hi_(hi), fhi_(fhi) {}

        constexpr T    lo() const noexcept { return lo_; }
        constexpr T    hi() const noexcept { return hi_; }
        constexpr T    flo() const noexcept { return flo_; }
        constexpr T    fhi() const noexcept { return fhi_; }
        constexpr T    width() const noexcept { return hi_ - lo_; }
        constexpr bool has_exact_zero() const noexcept { return flo_ == T(0) || fhi_ == T(0); }

        constexpr bracket<T> as_bracket() const noexcept { return bracket<T> { nxx::detail::trust_me {}, lo_, hi_ }; }

        // The endpoint with the smaller |f|, with this bracket as its enclosure.
        constexpr root_estimate<T> best() const noexcept;

        constexpr sign_bracket narrowed(const T& m, const T& fm) const noexcept
        {
            NXX_EXPECTS(lo_ < m && m < hi_);
            if (fm == T(0) || detail::opposite(flo_, fm)) return sign_bracket { nxx::detail::trust_me {}, lo_, flo_, m, fm };
            return sign_bracket { nxx::detail::trust_me {}, m, fm, hi_, fhi_ };
        }

        friend constexpr bool operator==(const sign_bracket&, const sign_bracket&) = default;
    };

    // The result of every 1-D root solver, so chains of different solvers type-check.
    template<real T>
    struct root_estimate
    {
        T                              x;
        T                              fx;
        T                              uncertainty;    // enclosure width (a bound) or |x_k - x_{k-1}| (an indicator); inf if unknown
        std::optional<sign_bracket<T>> enclosure;      // with samples: a bracketing successor starts without re-evaluating

        // x and f(x) are both required: an open method starts from the estimate's fx without evaluating f again, so a
        // value-initialised fx (root_estimate{1.0}) would pass for an exact zero that was never evaluated.
        constexpr root_estimate(T                              x_,
                                T                              fx_,
                                T                              uncertainty_ = detail::unknown<T>(),
                                std::optional<sign_bracket<T>> enclosure_   = std::nullopt)
            : x(x_),
              fx(fx_),
              uncertainty(uncertainty_),
              enclosure(std::move(enclosure_))
        {}

        friend constexpr bool operator==(const root_estimate&, const root_estimate&) = default;
    };

    template<real T>
    constexpr root_estimate<T> sign_bracket<T>::best() const noexcept
    {
        return math::abs(flo_) <= math::abs(fhi_) ? root_estimate<T> { lo_, flo_, width(), *this }
                                                  : root_estimate<T> { hi_, fhi_, width(), *this };
    }

    // The better of two failure payloads (DESIGN §6.7), a strict weak order, so the best estimate of a chain does not
    // depend on how its alternatives are grouped or folded: an estimate with a sign-changing enclosure beats one without
    // (the enclosure is a guarantee, a small |f(x)| is not); between two enclosures the narrower one (nested, hence newer)
    // wins; ties, and two estimates without enclosures, go by the smaller |f(x)|.
    template<real T>
    constexpr bool better_than(const root_estimate<T>& a, const root_estimate<T>& b) noexcept
    {
        if (a.enclosure.has_value() != b.enclosure.has_value()) return a.enclosure.has_value();
        if (a.enclosure) {
            const T wa = a.enclosure->width();
            const T wb = b.enclosure->width();
            if (wa < wb) return true;
            if (wb < wa) return false;
        }
        return math::abs(a.fx) < math::abs(b.fx);
    }

    // What the stop criteria of open methods see.
    //
    // distance() is the length of the step as the method PROPOSED it, before a projection shortened it: a clamped step
    // that stops just short of the edge says nothing about convergence, so it must not satisfy x_tol or step_tol (the
    // iterate would be reported as a root at the box edge). Without a projection it is |x_k - x_{k-1}|.
    template<real T>
    class point_view
    {
        T x_;
        T fx_;
        T step_;

    public:
        constexpr point_view(T x, T fx, T step = std::numeric_limits<T>::infinity()) noexcept : x_(x), fx_(fx), step_(step) {}
        constexpr T x() const noexcept { return x_; }
        constexpr T fx() const noexcept { return fx_; }
        constexpr T distance(const point_view& prev) const noexcept { return math::isfinite(step_) ? step_ : math::abs(x_ - prev.x_); }
        constexpr T residual() const noexcept { return math::abs(fx_); }
        constexpr T scale() const noexcept { return math::abs(x_); }
    };

    template<real T>
    class enclosure_bounds
    {
        T lo_;
        T hi_;

    public:
        constexpr enclosure_bounds(T lo, T hi) noexcept : lo_(lo), hi_(hi) {}
        constexpr T lo() const noexcept { return lo_; }
        constexpr T hi() const noexcept { return hi_; }
        constexpr T width() const noexcept { return hi_ - lo_; }
    };

    // What the stop criteria of bracketing methods see: no distance(), because successive iterates say nothing about
    // an enclosure (x_tol does not apply).
    template<real T>
    class enclosure_view
    {
        T x_;
        T fx_;
        T lo_;
        T hi_;

    public:
        constexpr enclosure_view(T x, T fx, T lo, T hi) noexcept : x_(x), fx_(fx), lo_(lo), hi_(hi) {}
        constexpr T                   x() const noexcept { return x_; }
        constexpr T                   fx() const noexcept { return fx_; }
        constexpr T                   residual() const noexcept { return math::abs(fx_); }
        constexpr T                   scale() const noexcept { return math::abs(x_); }
        constexpr enclosure_bounds<T> enclosure() const noexcept { return { lo_, hi_ }; }
    };

    template<class F, real T>
    using root_failure = failure<root_estimate<T>, callback_error_t<F, T>>;

    // Per-iterate projection (D20): applied to every proposed iterate before the function is evaluated. An iterate
    // pinned at the edge makes the solve fail with errc::stalled, never "the root is the boundary".
    template<real T>
    struct clamp_to
    {
        T lo;
        T hi;

        constexpr T operator()(const T& x) const noexcept { return x < lo ? lo : (hi < x ? hi : x); }
    };

    template<real T>
    clamp_to(T, T) -> clamp_to<T>;

    namespace detail
    {
        // Accepted inputs of the bracketing methods (DESIGN §6.5). kind: 0 bracket, 1 sampled (sign_bracket or a
        // search result), 2 braced {lo, hi} or std::pair, 3 the result of bracket<T>::make.
        template<class In>
        struct bracket_input
        {
            static constexpr bool value = false;
        };
        template<real T>
        struct bracket_input<bracket<T>>
        {
            static constexpr bool value = true;
            static constexpr int  kind  = 0;
            using scalar                = T;
        };
        template<real T>
        struct bracket_input<sign_bracket<T>>
        {
            static constexpr bool value = true;
            static constexpr int  kind  = 1;
            using scalar                = T;
        };
        template<real T>
        struct bracket_input<solution<sign_bracket<T>>>
        {
            static constexpr bool value = true;
            static constexpr int  kind  = 1;
            using scalar                = T;
        };
        template<real T>
        struct bracket_input<std::pair<T, T>>
        {
            static constexpr bool value = true;
            static constexpr int  kind  = 2;
            using scalar                = T;
        };
        template<real T>
        struct bracket_input<std::expected<bracket<T>, errc>>
        {
            static constexpr bool value = true;
            static constexpr int  kind  = 3;
            using scalar                = T;
        };

        template<class In>
        inline constexpr bool bracket_input_v = bracket_input<std::remove_cvref_t<In>>::value;
        template<class In>
        using bracket_scalar_t = typename bracket_input<std::remove_cvref_t<In>>::scalar;

        // Start windows of the searchers: a bracket-like input without samples.
        template<class In>
        inline constexpr bool window_input_v = [] {
            if constexpr (bracket_input_v<In>)
                return bracket_input<std::remove_cvref_t<In>>::kind != 1;
            else
                return false;
        }();

        // Accepted inputs of the open methods: a guess of a real type, or a root estimate (also a solution).
        template<class In>
        struct open_input
        {
            static constexpr bool value = real<In>;
            using scalar                = In;
        };
        template<real T>
        struct open_input<root_estimate<T>>
        {
            static constexpr bool value = true;
            using scalar                = T;
        };
        template<real T>
        struct open_input<solution<root_estimate<T>>>
        {
            static constexpr bool value = true;
            using scalar                = T;
        };

        template<class In>
        inline constexpr bool open_input_v = open_input<std::remove_cvref_t<In>>::value;
        template<class In>
        using open_scalar_t = typename open_input<std::remove_cvref_t<In>>::scalar;

        // Validates a window: invalid endpoints fail in-band with invalid_input (or the errc make() reported), at zero cost.
        template<class F, class In>
        constexpr auto to_window(algo id, const In& in)
            -> std::expected<bracket<bracket_scalar_t<In>>, root_failure<F, bracket_scalar_t<In>>>
        {
            using T            = bracket_scalar_t<In>;
            using Fail         = root_failure<F, T>;
            constexpr int kind = bracket_input<std::remove_cvref_t<In>>::kind;
            if constexpr (kind == 0)
                return in;
            else if constexpr (kind == 2) {
                auto b = bracket<T>::make(in.first, in.second);
                if (!b) return std::unexpected(Fail { b.error(), id, {}, std::nullopt, {} });
                return *b;
            }
            else {
                static_assert(kind == 3);
                if (!in) return std::unexpected(Fail { in.error(), id, {}, std::nullopt, {} });
                return *in;
            }
        }

        // Samples both ends of a window: the prepare() of every bracketing method. Samples may be +-inf.
        template<class F, real T>
        constexpr auto sample(algo id, const F& fn, const bracket<T>& b) -> std::expected<problem<F, sign_bracket<T>>, root_failure<F, T>>
        {
            using Fail = root_failure<F, T>;
            auto flo   = nxx::evaluate_sample(fn, b.lo());
            if (!flo)
                return std::unexpected(Fail { flo.error().code, id, counters { 0, flo.error().evals }, std::nullopt, flo.error().cause });
            const std::uint32_t one = cost_of(fn);
            auto                fhi = nxx::evaluate_sample(fn, b.hi());
            if (!fhi)
                return std::unexpected(Fail { fhi.error().code,
                                              id,
                                              counters { 0, one + fhi.error().evals },
                                              root_estimate<T> { b.lo(), *flo, unknown<T>(), std::nullopt },
                                              fhi.error().cause });
            const std::uint32_t used = one + cost_of(fn);
            if (!detail::opposite(*flo, *fhi)) {
                const root_estimate<T> best = math::abs(*flo) <= math::abs(*fhi)
                                                  ? root_estimate<T> { b.lo(), *flo, unknown<T>(), std::nullopt }
                                                  : root_estimate<T> { b.hi(), *fhi, unknown<T>(), std::nullopt };
                return std::unexpected(Fail { errc::no_sign_change, id, counters { 0, used }, best, {} });
            }
            return problem<F, sign_bracket<T>> { fn, sign_bracket<T> { nxx::detail::trust_me {}, b.lo(), *flo, b.hi(), *fhi }, used };
        }

        template<class F, class In>
        constexpr auto prepare_bracketing(algo id, const F& fn, const In& in)
            -> std::expected<problem<F, sign_bracket<bracket_scalar_t<In>>>, root_failure<F, bracket_scalar_t<In>>>
        {
            using T = bracket_scalar_t<In>;
            if constexpr (bracket_input<std::remove_cvref_t<In>>::kind == 1)
                return problem<F, sign_bracket<T>> { fn, static_cast<const sign_bracket<T>&>(in), 0 };    // already sampled
            else {
                auto b = detail::to_window<F>(id, in);
                if (!b) return std::unexpected(std::move(b).error());
                return detail::sample(id, fn, *b);
            }
        }

        // The start of an open method: the guess, and f(guess) when a previous stage already evaluated it.
        template<real T>
        struct open_start
        {
            T                x0;
            std::optional<T> fx0;
        };

        template<class F, class In>
        constexpr auto prepare_open(algo id, const F& fn, const In& in)
            -> std::expected<problem<F, open_start<open_scalar_t<In>>>, root_failure<F, open_scalar_t<In>>>
        {
            using T    = open_scalar_t<In>;
            using Fail = root_failure<F, T>;
            if constexpr (real<std::remove_cvref_t<In>>) {
                if (!math::isfinite(in)) return std::unexpected(Fail { errc::non_finite_input, id, {}, std::nullopt, {} });
                return problem<F, open_start<T>> { fn, open_start<T> { in, std::nullopt }, 0 };
            }
            else {
                const root_estimate<T>& e = in;
                if (!math::isfinite(e.x)) return std::unexpected(Fail { errc::non_finite_input, id, {}, std::nullopt, {} });
                return problem<F, open_start<T>> { fn,
                                                   open_start<T> { e.x, math::isfinite(e.fx) ? std::optional<T>(e.fx) : std::nullopt },
                                                   0 };
            }
        }

        template<class UE, class UE2>
        constexpr fault<UE> widen(const fault<UE2>& f) noexcept
        {
            if constexpr (std::is_same_v<UE, UE2>)
                return f;
            else
                return fault<UE> { f.code, f.evals, {} };
        }

        template<class UE, class Est, class UE2>
        constexpr failure<Est, UE> widen(const failure<Est, UE2>& f)
        {
            if constexpr (std::is_same_v<UE, UE2>)
                return f;
            else
                return failure<Est, UE> { f.code, f.where, f.used, f.best, {} };
        }

        // The pole check of the bracketing methods (DESIGN §7.2): after a criterion or resolution-limit stop, a residual
        // that grew while the bracket shrank means the sign change is a pole (tan on [1, 2]), not a root. The failure,
        // or nothing when the solution stands (the driver's finish protocol, core/iterate.hpp).
        //
        // The initial level is the larger of the FINITE initial samples: an endpoint may sample +-inf (log(x) on
        // [0, 2]), and an infinite level would switch the check off exactly where a pole sits at an end (1/x on
        // [-1, 0]). A success whose f(x) is not finite is never accepted (a pole reached exactly). Known limit: a large
        // finite endpoint sample (1/x on [-1e-20, 1]) still raises the level and hides the pole.
        template<class Fail, real T>
        constexpr auto pole_check(const sign_bracket<T>& initial, const solution<root_estimate<T>>& sol) -> std::optional<Fail>
        {
            const auto pole = [&] {
                return std::optional<Fail>(
                    Fail { errc::sign_change_not_root, sol.by, sol.used, static_cast<const root_estimate<T>&>(sol), {} });
            };
            if (!math::isfinite(sol.fx)) return pole();
            if (sol.how == stop_reason::exact_zero || !sol.enclosure) return std::nullopt;
            const T    flo       = math::abs(initial.flo());
            const T    fhi       = math::abs(initial.fhi());
            const bool lo_finite = math::isfinite(flo);
            const bool hi_finite = math::isfinite(fhi);
            if (!lo_finite && !hi_finite) return std::nullopt;    // both ends infinite: only the finiteness test applies
            const T before = !lo_finite ? fhi : (!hi_finite ? flo : (std::max)(flo, fhi));
            const T now    = (std::min)(math::abs(sol.enclosure->flo()), math::abs(sol.enclosure->fhi()));
            if (!(now > before)) return std::nullopt;
            return pole();
        }
    }    // namespace detail
}    // namespace nxx::roots

NXX_END_HEADER

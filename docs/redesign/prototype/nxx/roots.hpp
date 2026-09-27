// Prototype of PLAN_v1 6.5/6.6/6.9/7.2: 1-D root solvers as immutable values with init/step.
// Does NOT include deriv (plan 5.2 DAG): derivative sources are recognised structurally.
#pragma once
#include "core.hpp"

namespace nxx::roots {

template<real T> struct root_estimate {
    T x;
    T fx;
    std::optional<bracket<T>> enclosure;
};
template<real T> constexpr T merit_of(const root_estimate<T>& e) noexcept { return math::abs(e.fx); }

namespace detail {
    template<real T> constexpr bool opposite(const T& a, const T& b) noexcept {   // comparisons, no product
        return (a <= T(0) && b >= T(0)) || (a >= T(0) && b <= T(0));
    }
}

// sign_bracket: bracket + endpoint samples of opposite sign (or an exact zero). Only transition: narrowed().
template<real T> class sign_bracket {
    T lo_, flo_, hi_, fhi_;
public:
    constexpr sign_bracket(nxx::detail::trust_me, T lo, T flo, T hi, T fhi) noexcept
        : lo_(lo), flo_(flo), hi_(hi), fhi_(fhi) {}
    constexpr T lo() const noexcept { return lo_; }
    constexpr T hi() const noexcept { return hi_; }
    constexpr T flo() const noexcept { return flo_; }
    constexpr T fhi() const noexcept { return fhi_; }
    constexpr T width() const noexcept { return hi_ - lo_; }
    constexpr bracket<T> as_bracket() const noexcept { return bracket<T>{nxx::detail::trust_me{}, lo_, hi_}; }
    constexpr bool has_exact_zero() const noexcept { return flo_ == T(0) || fhi_ == T(0); }
    constexpr root_estimate<T> best() const noexcept {
        return math::abs(flo_) <= math::abs(fhi_) ? root_estimate<T>{lo_, flo_, as_bracket()}
                                                   : root_estimate<T>{hi_, fhi_, as_bracket()};
    }
    // keeps the half that still changes sign; m must be strictly inside
    constexpr sign_bracket narrowed(T m, T fm) const noexcept {
        if (fm == T(0)) return sign_bracket{nxx::detail::trust_me{}, lo_, flo_, m, fm};
        if (detail::opposite(flo_, fm)) return sign_bracket{nxx::detail::trust_me{}, lo_, flo_, m, fm};
        return sign_bracket{nxx::detail::trust_me{}, m, fm, hi_, fhi_};
    }
};
template<real T> constexpr T merit_of(const sign_bracket<T>& b) noexcept { return merit_of(b.best()); }

// Views: what stop criteria see (core.hpp protocol).
template<real T> struct point_view {
    T x_, fx_;
    constexpr T x() const noexcept { return x_; }
    constexpr T fx() const noexcept { return fx_; }
    constexpr T distance(const point_view& p) const noexcept { return math::abs(x_ - p.x_); }
    constexpr T residual() const noexcept { return math::abs(fx_); }
    constexpr T scale() const noexcept { return math::abs(x_); }
};
template<real T> struct bracket_view : point_view<T> {
    bracket<T> enc;
    constexpr bracket<T> enclosure() const noexcept { return enc; }
};

template<class F, real T> using root_failure = failure<root_estimate<T>, callback_error_t<F, T>>;

namespace detail {
    template<class F, real T>
    constexpr auto make_sign_bracket(const F& f, const bracket<T>& b) -> std::expected<sign_bracket<T>, root_failure<F, T>> {
        using Fail = root_failure<F, T>;
        auto fa = evaluate(f, b.lo());
        if (!fa) return std::unexpected(Fail{fa.error().code, algo::none, {0, 1}, std::nullopt, fa.error().cause});
        auto fb = evaluate(f, b.hi());
        if (!fb) return std::unexpected(Fail{fb.error().code, algo::none, {0, 2},
                                             root_estimate<T>{b.lo(), *fa, std::nullopt}, fb.error().cause});
        const sign_bracket<T> sb{nxx::detail::trust_me{}, b.lo(), *fa, b.hi(), *fb};
        if (!opposite(*fa, *fb)) return std::unexpected(Fail{errc::no_sign_change, algo::none, {0, 2}, sb.best(), {}});
        return sb;
    }
    template<class UE, class UE2> constexpr fault<UE> widen(const fault<UE2>& f) {   // fault<none> -> fault<UE>
        if constexpr (std::is_same_v<UE, UE2>) return f; else return fault<UE>{f.code, {}};
    }
}

// ================================================================ bisection
template<real T> struct bisection_state { sign_bracket<T> b; std::uint32_t nfev; };

template<class Stop = floored_width>
class bisection : public solver_facade {
    Stop stop_{};
    max_iterations budget_ = 100;
public:
    static constexpr algo id = algo::bisection;
    constexpr bisection() = default;
    constexpr explicit bisection(Stop s, max_iterations b = 100) noexcept : stop_(s), budget_(b) {}
    constexpr const Stop& stop() const noexcept { return stop_; }
    constexpr max_iterations budget() const noexcept { return budget_; }

    template<class F, real T>
    constexpr auto prepare(const F& f, const bracket<T>& b) const -> std::expected<problem<F, sign_bracket<T>>, root_failure<F, T>> {
        auto sb = detail::make_sign_bracket(f, b);
        if (!sb) { auto e = sb.error(); e.where = id; return std::unexpected(e); }
        return problem<F, sign_bracket<T>>{f, *sb, 2};
    }
    template<class F, real T>
    constexpr auto prepare(const F& f, const sign_bracket<T>& sb) const -> std::expected<problem<F, sign_bracket<T>>, root_failure<F, T>> {
        return problem<F, sign_bracket<T>>{f, sb, 0};
    }
    template<class Self, class F, real T>
    void operator()(this const Self&, const F&, const T&)
        NXX_DELETE("bisection needs a bracket: nxx::bracket{lo, hi}, a sign_bracket or a search result");
    using solver_facade::operator();

    template<class F, real T>
    constexpr auto init(const problem<F, sign_bracket<T>>& p) const -> std::expected<bisection_state<T>, root_failure<F, T>> {
        return bisection_state<T>{p.in, p.nfev0};
    }
    template<class F, real T>
    constexpr auto step(const problem<F, sign_bracket<T>>& p, const bisection_state<T>& s) const
        -> std::expected<bisection_state<T>, fault<callback_error_t<F, T>>> {
        const T m = s.b.lo() + (s.b.hi() - s.b.lo()) / T(2);
        auto fm = evaluate(p.f, m);
        if (!fm) return std::unexpected(fm.error());
        return bisection_state<T>{s.b.narrowed(m, *fm), s.nfev + 1};
    }
    template<real T> constexpr bracket_view<T> view(const bisection_state<T>& s) const noexcept {
        const auto e = s.b.best();
        return bracket_view<T>{{e.x, e.fx}, s.b.as_bracket()};
    }
    template<real T> constexpr root_estimate<T> estimate(const bisection_state<T>& s) const noexcept { return s.b.best(); }
    template<real T> constexpr root_estimate<T> best(const bisection_state<T>& s) const noexcept { return s.b.best(); }
    template<real T> constexpr std::optional<stop_reason> intrinsic(const bisection_state<T>& s) const noexcept {
        if (s.b.has_exact_zero()) return stop_reason::exact_zero;
        const T m = s.b.lo() + (s.b.hi() - s.b.lo()) / T(2);
        if (!(s.b.lo() < m && m < s.b.hi())) return stop_reason::resolution_limit;   // unsplittable
        return std::nullopt;
    }
};
bisection() -> bisection<>;

// ================================================================ brent (zeroin)
template<real T> struct brent_state { T a, fa, b, fb, c, fc, d, e; std::uint32_t nfev; };

template<class Stop = never>
class brent : public solver_facade {
    Stop stop_{};
    max_iterations budget_ = 100;
    template<real T> constexpr T xtol() const noexcept {   // Brent's internal tolerance is configuration (D6)
        if constexpr (requires { stop_.abs_tolerance(); }) return T(stop_.abs_tolerance()); else return T(0);
    }
    template<real T> constexpr T tol1(const brent_state<T>& s) const noexcept {
        return T(2) * std::numeric_limits<T>::epsilon() * math::abs(s.b) + xtol<T>() / T(2);
    }
    template<real T> static constexpr void normalise(brent_state<T>& s) noexcept {
        if (!detail::opposite(s.fb, s.fc)) { s.c = s.a; s.fc = s.fa; s.d = s.e = s.b - s.a; }
        if (math::abs(s.fc) < math::abs(s.fb)) { s.a = s.b; s.b = s.c; s.c = s.a; s.fa = s.fb; s.fb = s.fc; s.fc = s.fa; }
    }
public:
    static constexpr algo id = algo::brent;
    constexpr brent() = default;
    constexpr explicit brent(Stop s, max_iterations b = 100) noexcept : stop_(s), budget_(b) {}
    constexpr const Stop& stop() const noexcept { return stop_; }
    constexpr max_iterations budget() const noexcept { return budget_; }

    template<class F, real T>
    constexpr auto prepare(const F& f, const bracket<T>& b) const -> std::expected<problem<F, sign_bracket<T>>, root_failure<F, T>> {
        auto sb = detail::make_sign_bracket(f, b);
        if (!sb) { auto e = sb.error(); e.where = id; return std::unexpected(e); }
        return problem<F, sign_bracket<T>>{f, *sb, 2};
    }
    template<class F, real T>
    constexpr auto prepare(const F& f, const sign_bracket<T>& sb) const -> std::expected<problem<F, sign_bracket<T>>, root_failure<F, T>> {
        return problem<F, sign_bracket<T>>{f, sb, 0};
    }
    template<class F, real T>
    constexpr auto init(const problem<F, sign_bracket<T>>& p) const -> std::expected<brent_state<T>, root_failure<F, T>> {
        brent_state<T> s{p.in.lo(), p.in.flo(), p.in.hi(), p.in.fhi(), p.in.lo(), p.in.flo(), T(0), T(0), p.nfev0};
        s.d = s.e = s.b - s.a;
        normalise(s);
        return s;
    }
    template<class F, real T>
    constexpr auto step(const problem<F, sign_bracket<T>>& p, const brent_state<T>& s0) const
        -> std::expected<brent_state<T>, fault<callback_error_t<F, T>>> {
        brent_state<T> s = s0;   // a local copy; the state value itself is never mutated
        const T t1 = tol1(s), xm = (s.c - s.b) / T(2);
        if (math::abs(s.e) >= t1 && math::abs(s.fa) > math::abs(s.fb)) {
            T pp, q, r;
            const T sr = s.fb / s.fa;
            if (s.a == s.c) { pp = T(2) * xm * sr; q = T(1) - sr; }
            else {
                q = s.fa / s.fc; r = s.fb / s.fc;
                pp = sr * (T(2) * xm * q * (q - r) - (s.b - s.a) * (r - T(1)));
                q = (q - T(1)) * (r - T(1)) * (sr - T(1));
            }
            if (pp > T(0)) q = -q;
            pp = math::abs(pp);
            const T min1 = T(3) * xm * q - math::abs(t1 * q), min2 = math::abs(s.e * q);
            if (T(2) * pp < math::min(min1, min2)) { s.e = s.d; s.d = pp / q; }
            else { s.d = xm; s.e = s.d; }
        } else { s.d = xm; s.e = s.d; }
        s.a = s.b; s.fa = s.fb;
        s.b = math::abs(s.d) > t1 ? s.b + s.d : s.b + (xm > T(0) ? t1 : -t1);
        auto fb = evaluate(p.f, s.b);
        if (!fb) return std::unexpected(fb.error());
        s.fb = *fb;
        s.nfev += 1;
        normalise(s);
        return s;
    }
    template<real T> constexpr bracket_view<T> view(const brent_state<T>& s) const noexcept {
        return bracket_view<T>{{s.b, s.fb}, bracket<T>{nxx::detail::trust_me{}, math::min(s.b, s.c), math::max(s.b, s.c)}};
    }
    template<real T> constexpr root_estimate<T> estimate(const brent_state<T>& s) const noexcept {
        if (s.b == s.c) return root_estimate<T>{s.b, s.fb, std::nullopt};
        return root_estimate<T>{s.b, s.fb, bracket<T>{nxx::detail::trust_me{}, math::min(s.b, s.c), math::max(s.b, s.c)}};
    }
    template<real T> constexpr root_estimate<T> best(const brent_state<T>& s) const noexcept { return estimate(s); }
    template<real T> constexpr std::optional<stop_reason> intrinsic(const brent_state<T>& s) const noexcept {
        if (s.fb == T(0)) return stop_reason::exact_zero;
        if (math::abs((s.c - s.b) / T(2)) <= tol1(s)) return stop_reason::algorithm;
        return std::nullopt;
    }
};
brent() -> brent<>;

// ================================================================ projection support
struct no_projection { template<class X> constexpr X operator()(const X& x) const noexcept { return x; } };
template<real T> struct clamp_to {   // D20: per-iterate projection
    T lo, hi;
    constexpr T operator()(const T& x) const noexcept { return x < lo ? lo : (hi < x ? hi : x); }
};
template<real T> clamp_to(T, T) -> clamp_to<T>;

// ================================================================ secant (derivative-free)
template<real T> struct secant_state { T x0, f0, x1, f1; std::uint32_t nfev; };

template<class Stop = default_step>
class secant : public solver_facade {
    Stop stop_{};
    max_iterations budget_ = 100;
public:
    static constexpr algo id = algo::secant;
    constexpr secant() = default;
    constexpr explicit secant(Stop s, max_iterations b = 100) noexcept : stop_(s), budget_(b) {}
    constexpr const Stop& stop() const noexcept { return stop_; }
    constexpr max_iterations budget() const noexcept { return budget_; }

    template<class F, real T>
    constexpr auto prepare(const F& f, const T& x0) const -> std::expected<problem<F, T>, root_failure<F, T>> {
        if (!math::isfinite(x0)) return std::unexpected(root_failure<F, T>{errc::non_finite_input, id, {}, std::nullopt, {}});
        return problem<F, T>{f, x0, 0};
    }
    template<class F, real T>
    constexpr auto prepare(const F& f, const root_estimate<T>& e) const { return prepare(f, e.x); }
    template<class F, real T>
    constexpr auto init(const problem<F, T>& p) const -> std::expected<secant_state<T>, root_failure<F, T>> {
        using Fail = root_failure<F, T>;
        auto f0 = evaluate(p.f, p.in);
        if (!f0) return std::unexpected(Fail{f0.error().code, id, {0, 1}, std::nullopt, f0.error().cause});
        const T x1 = p.in + math::pow2<T>(-10) * math::max(math::abs(p.in), T(1));
        auto f1 = evaluate(p.f, x1);
        if (!f1) return std::unexpected(Fail{f1.error().code, id, {0, 2}, root_estimate<T>{p.in, *f0, std::nullopt}, f1.error().cause});
        return secant_state<T>{p.in, *f0, x1, *f1, 2};
    }
    template<class F, real T>
    constexpr auto step(const problem<F, T>& p, const secant_state<T>& s) const
        -> std::expected<secant_state<T>, fault<callback_error_t<F, T>>> {
        if (s.f1 == s.f0) return std::unexpected(fault<callback_error_t<F, T>>{errc::stalled, {}});   // flat secant
        const T x2 = s.x1 - s.f1 * (s.x1 - s.x0) / (s.f1 - s.f0);
        auto f2 = evaluate(p.f, x2);
        if (!f2) return std::unexpected(f2.error());
        return secant_state<T>{s.x1, s.f1, x2, *f2, s.nfev + 1};
    }
    template<real T> constexpr point_view<T> view(const secant_state<T>& s) const noexcept { return {s.x1, s.f1}; }
    template<real T> constexpr root_estimate<T> estimate(const secant_state<T>& s) const noexcept { return {s.x1, s.f1, std::nullopt}; }
    template<real T> constexpr root_estimate<T> best(const secant_state<T>& s) const noexcept {
        return math::abs(s.f0) < math::abs(s.f1) ? root_estimate<T>{s.x0, s.f0, std::nullopt} : estimate(s);
    }
    template<real T> constexpr std::optional<stop_reason> intrinsic(const secant_state<T>& s) const noexcept {
        if (s.f1 == T(0)) return stop_reason::exact_zero;
        return std::nullopt;
    }
};
secant() -> secant<>;

// ================================================================ newton (explicit derivative source, D13)
struct no_derivative {};
template<class F> concept has_member_derivative = requires(const F& f) { f.derivative(); };

namespace detail {
    // what can serve as f' for f: (a) a callable df, (b) a policy with bind(f) (e.g. deriv::numeric),
    // (c) nothing, if f itself has .derivative() (polynomial, spline)
    template<class D, class F> constexpr auto bind_derivative(const D& d, const F& f) {
        if constexpr (requires { d.bind(f); }) return d.bind(f);
        else if constexpr (std::is_same_v<D, no_derivative>) return nxx::detail::unref(f).derivative();
        else return d;
    }
    template<class D, class F, class T>
    inline constexpr bool has_derivative_source_v =
        !std::is_same_v<D, no_derivative> || has_member_derivative<std::remove_cvref_t<decltype(nxx::detail::unref(std::declval<const F&>()))>>;
}
template<real T, class DF> struct newton_input { T x0; DF df; };
template<real T> struct newton_state { T x, fx; std::uint32_t nfev; };

template<class Stop = default_step, class D = no_derivative, class Proj = no_projection>
class newton : public solver_facade {
    Stop stop_{};
    max_iterations budget_ = 100;
    NXX_NO_UNIQUE_ADDRESS D d_{};
    NXX_NO_UNIQUE_ADDRESS Proj proj_{};
    template<class, class, class> friend class newton;
    constexpr newton(Stop s, max_iterations b, D d, Proj p) noexcept : stop_(s), budget_(b), d_(d), proj_(p) {}
public:
    static constexpr algo id = algo::newton;
    constexpr newton() = default;
    constexpr explicit newton(Stop s, max_iterations b = 100) noexcept : stop_(s), budget_(b) {}
    constexpr const Stop& stop() const noexcept { return stop_; }
    constexpr max_iterations budget() const noexcept { return budget_; }
    template<class D2> constexpr auto with_derivative(D2 d) const { return newton<Stop, D2, Proj>{stop_, budget_, d, proj_}; }
    template<class P2> constexpr auto with_projection(P2 p) const { return newton<Stop, D, P2>{stop_, budget_, d_, p}; }

    template<class F, real T>
        requires detail::has_derivative_source_v<D, F, T>
    constexpr auto prepare(const F& f, const T& x0) const {
        using DF = decltype(detail::bind_derivative(d_, f));
        using UE = common_cause_t<callback_error_t<F, T>, callback_error_t<DF, T>>;
        using R = std::expected<problem<F, newton_input<T, DF>>, failure<root_estimate<T>, UE>>;
        if (!math::isfinite(x0)) return R{std::unexpect, failure<root_estimate<T>, UE>{errc::non_finite_input, id, {}, std::nullopt, {}}};
        return R{problem<F, newton_input<T, DF>>{f, {x0, detail::bind_derivative(d_, f)}, 0}};
    }
    template<class F, real T>
        requires detail::has_derivative_source_v<D, F, T>
    constexpr auto prepare(const F& f, const root_estimate<T>& e) const { return prepare(f, e.x); }

    template<class Self, class F, class In>
        requires(!detail::has_derivative_source_v<D, F, In>)
    void operator()(this const Self&, const F&, const In&)
        NXX_DELETE("newton needs a derivative: .with_derivative(df), .with_derivative(deriv::numeric{}), "
                   "a callable with .derivative(), or use secant");
    using solver_facade::operator();

    template<class F, real T, class DF>
    constexpr auto init(const problem<F, newton_input<T, DF>>& p) const {
        using UE = common_cause_t<callback_error_t<F, T>, callback_error_t<DF, T>>;
        using R = std::expected<newton_state<T>, failure<root_estimate<T>, UE>>;
        auto fx = evaluate(p.f, p.in.x0);
        if (!fx) return R{std::unexpect, failure<root_estimate<T>, UE>{fx.error().code, id, {0, 1}, std::nullopt,
                                                                       detail::widen<UE>(fx.error()).cause}};
        return R{newton_state<T>{p.in.x0, *fx, 1}};
    }
    template<class F, real T, class DF>
    constexpr auto step(const problem<F, newton_input<T, DF>>& p, const newton_state<T>& s) const
        -> std::expected<newton_state<T>, fault<common_cause_t<callback_error_t<F, T>, callback_error_t<DF, T>>>> {
        using UE = common_cause_t<callback_error_t<F, T>, callback_error_t<DF, T>>;
        auto d = evaluate(p.in.df, s.x);
        if (!d) return std::unexpected(detail::widen<UE>(d.error()));
        if (*d == T(0)) return std::unexpected(fault<UE>{errc::zero_derivative, {}});
        const T proposed = s.x - s.fx / *d;
        const T x1 = proj_(proposed);
        if (x1 != proposed && x1 == s.x) return std::unexpected(fault<UE>{errc::stalled, {}});   // pinned at the edge
        auto fx = evaluate(p.f, x1);
        if (!fx) return std::unexpected(detail::widen<UE>(fx.error()));
        return newton_state<T>{x1, *fx, s.nfev + 2};
    }
    template<real T> constexpr point_view<T> view(const newton_state<T>& s) const noexcept { return {s.x, s.fx}; }
    template<real T> constexpr root_estimate<T> estimate(const newton_state<T>& s) const noexcept { return {s.x, s.fx, std::nullopt}; }
    template<real T> constexpr root_estimate<T> best(const newton_state<T>& s) const noexcept { return estimate(s); }
    template<real T> constexpr std::optional<stop_reason> intrinsic(const newton_state<T>& s) const noexcept {
        if (s.fx == T(0)) return stop_reason::exact_zero;
        return std::nullopt;
    }
};
newton() -> newton<>;

// ================================================================ expand_out (searcher: succeeds with a sign_bracket)
template<real T> struct expand_state { T lo, flo, hi, fhi; std::uint32_t nfev; };

class expand_out : public solver_facade {
    max_iterations budget_ = 60;
public:
    static constexpr algo id = algo::expand_out;
    constexpr expand_out() = default;
    constexpr explicit expand_out(max_iterations b) noexcept : budget_(b) {}
    constexpr never stop() const noexcept { return {}; }   // not configurable: estimate() is only valid on a sign change
    constexpr max_iterations budget() const noexcept { return budget_; }

    template<class F, real T>
    constexpr auto prepare(const F& f, const bracket<T>& b) const -> std::expected<problem<F, bracket<T>>, root_failure<F, T>> {
        return problem<F, bracket<T>>{f, b, 0};
    }
    template<class F, real T>
    constexpr auto init(const problem<F, bracket<T>>& p) const -> std::expected<expand_state<T>, root_failure<F, T>> {
        using Fail = root_failure<F, T>;
        auto fa = evaluate(p.f, p.in.lo());
        if (!fa) return std::unexpected(Fail{fa.error().code, id, {0, 1}, std::nullopt, fa.error().cause});
        auto fb = evaluate(p.f, p.in.hi());
        if (!fb) return std::unexpected(Fail{fb.error().code, id, {0, 2}, root_estimate<T>{p.in.lo(), *fa, std::nullopt}, fb.error().cause});
        return expand_state<T>{p.in.lo(), *fa, p.in.hi(), *fb, 2};
    }
    template<class F, real T>
    constexpr auto step(const problem<F, bracket<T>>& p, const expand_state<T>& s) const
        -> std::expected<expand_state<T>, fault<callback_error_t<F, T>>> {
        const T w = (s.hi - s.lo) * T(1.6);
        expand_state<T> n = s;
        if (math::abs(s.flo) < math::abs(s.fhi)) {
            n.lo = s.lo - w;
            auto f = evaluate(p.f, n.lo);
            if (!f) return std::unexpected(f.error());
            n.flo = *f;
        } else {
            n.hi = s.hi + w;
            auto f = evaluate(p.f, n.hi);
            if (!f) return std::unexpected(f.error());
            n.fhi = *f;
        }
        n.nfev += 1;
        return n;
    }
    template<real T> constexpr point_view<T> view(const expand_state<T>& s) const noexcept { return {best(s).x, best(s).fx}; }
    template<real T> constexpr sign_bracket<T> estimate(const expand_state<T>& s) const noexcept {
        return sign_bracket<T>{nxx::detail::trust_me{}, s.lo, s.flo, s.hi, s.fhi};
    }
    template<real T> constexpr root_estimate<T> best(const expand_state<T>& s) const noexcept {
        return estimate(s).best();
    }
    template<real T> constexpr std::optional<stop_reason> intrinsic(const expand_state<T>& s) const noexcept {
        if (detail::opposite(s.flo, s.fhi)) return stop_reason::algorithm;
        return std::nullopt;
    }
};

}    // namespace nxx::roots

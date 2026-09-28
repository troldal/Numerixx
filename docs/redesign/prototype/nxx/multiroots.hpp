// Prototype of PLAN_v1 7.5 multiroots: damped Newton (Armijo backtracking) as an immutable solver value,
// generic over any dense_vector (in-house vec or Eigen) through linalg::vector_traits.
#pragma once
#include "linalg.hpp"

namespace nxx::multiroots {

template<class V> using traits = linalg::vector_traits<V>;
template<class V> using scalar_of = typename traits<V>::scalar;

template<class V> struct system_estimate {
    V x;
    V fx;
    scalar_of<V> merit;   // 0.5 * ||F||^2
};
template<class V> constexpr auto merit_of(const system_estimate<V>& e) noexcept { return e.merit; }

template<class V> struct nd_view {
    V x_, fx_;
    constexpr scalar_of<V> distance(const nd_view& p) const { const V d = x_ - p.x_; return traits<V>::norm_inf(d); }
    constexpr scalar_of<V> residual() const { return traits<V>::norm_inf(fx_); }
    constexpr scalar_of<V> scale() const { return traits<V>::norm_inf(x_); }
};

template<class V> struct nd_state { V x, fx; scalar_of<V> merit; std::uint32_t nfev; };

struct forward_difference {};   // default Jacobian source: n extra evaluations per step
struct no_projection { template<class V> constexpr V operator()(const V& x) const { return x; } };

template<class Stop = default_step, class Jac = forward_difference, class Proj = no_projection>
class newton : public solver_facade {
    Stop stop_{};
    max_iterations budget_ = 50;
    NXX_NO_UNIQUE_ADDRESS Jac jac_{};
    NXX_NO_UNIQUE_ADDRESS Proj proj_{};
    template<class, class, class> friend class newton;
    constexpr newton(Stop s, max_iterations b, Jac j, Proj p) : stop_(s), budget_(b), jac_(j), proj_(p) {}

    template<class F, class V>
    constexpr auto eval(const F& f, const V& x) const -> std::expected<V, fault<callback_error_t<F, V>>> {
        auto y = evaluate(f, x);
        if (y && !traits<V>::all_finite(*y)) return std::unexpected(fault<callback_error_t<F, V>>{errc::non_finite_value, {}});
        return y;
    }
    template<class F, class V>
    constexpr auto jacobian(const F& f, const V& x, const V& fx, std::uint32_t& nfev) const
        -> std::expected<typename traits<V>::matrix, fault<callback_error_t<F, V>>> {
        using T = scalar_of<V>;
        if constexpr (std::is_same_v<Jac, forward_difference>) {
            auto J = traits<V>::zero_matrix(x);
            const std::size_t n = traits<V>::size(x);
            for (std::size_t j = 0; j < n; ++j) {
                V xh = x;
                const T h0 = math::root_eps<T>(1, 2) * math::max(math::abs(x[j]), T(1));
                xh[j] = x[j] + h0;
                const T h = xh[j] - x[j];
                auto fj = eval(f, xh);
                ++nfev;
                if (!fj) return std::unexpected(fj.error());
                for (std::size_t i = 0; i < n; ++i) J(i, j) = ((*fj)[i] - fx[i]) / h;
            }
            return J;
        } else {
            ++nfev;
            return typename traits<V>::matrix(std::invoke(jac_, x));
        }
    }
public:
    static constexpr algo id = algo::multi_newton;
    constexpr newton() = default;
    constexpr explicit newton(Stop s, max_iterations b = 50) : stop_(s), budget_(b) {}
    constexpr const Stop& stop() const noexcept { return stop_; }
    constexpr max_iterations budget() const noexcept { return budget_; }
    template<class J2> constexpr auto with_jacobian(J2 j) const { return newton<Stop, J2, Proj>{stop_, budget_, j, proj_}; }
    template<class P2> constexpr auto with_projection(P2 p) const { return newton<Stop, Jac, P2>{stop_, budget_, jac_, p}; }

    template<class F, linalg::dense_vector V>
    constexpr auto prepare(const F& f, const V& x0) const
        -> std::expected<problem<F, V>, failure<system_estimate<V>, callback_error_t<F, V>>> {
        if (!traits<V>::all_finite(x0))
            return std::unexpected(failure<system_estimate<V>, callback_error_t<F, V>>{errc::non_finite_input, id, {}, std::nullopt, {}});
        return problem<F, V>{f, x0, 0};
    }
    template<class F, class V>
    constexpr auto init(const problem<F, V>& p) const -> std::expected<nd_state<V>, failure<system_estimate<V>, callback_error_t<F, V>>> {
        using Fail = failure<system_estimate<V>, callback_error_t<F, V>>;
        auto fx = eval(p.f, p.in);
        if (!fx) return std::unexpected(Fail{fx.error().code, id, {0, 1}, std::nullopt, fx.error().cause});
        return nd_state<V>{p.in, *fx, traits<V>::sq_norm(*fx) / 2, 1};
    }
    template<class F, class V>
    constexpr auto step(const problem<F, V>& p, const nd_state<V>& s) const -> std::expected<nd_state<V>, fault<callback_error_t<F, V>>> {
        using T = scalar_of<V>;
        using Fault = fault<callback_error_t<F, V>>;
        std::uint32_t nfev = s.nfev;
        auto J = jacobian(p.f, s.x, s.fx, nfev);
        if (!J) return std::unexpected(J.error());
        const V rhs = -s.fx;
        auto dx = traits<V>::solve(*J, rhs);
        if (!dx) return std::unexpected(Fault{dx.error(), {}});   // singular Jacobian
        const T eps = std::numeric_limits<T>::epsilon();
        const T xs = math::max(T(1), traits<V>::norm_inf(s.x));
        T lam(1);
        for (int ls = 0; ls < 30; ++ls, lam /= T(2)) {
            const V proposed = s.x + lam * *dx;
            const V xn = proj_(proposed);
            const V moved = xn - s.x;
            if (!(xn == proposed) && traits<V>::norm_inf(moved) == T(0)) return std::unexpected(Fault{errc::stalled, {}});
            auto fn = eval(p.f, xn);
            ++nfev;
            if (!fn) { if (fn.error().code == errc::non_finite_value) continue; return std::unexpected(fn.error()); }
            const T mn = traits<V>::sq_norm(*fn) / 2;
            const bool armijo = mn <= (T(1) - T(2) * T(1e-4) * lam) * s.merit;
            const bool rounding_level = traits<V>::norm_inf(moved) <= T(16) * eps * xs;   // A's lesson: no false stall
            if (armijo || rounding_level) return nd_state<V>{xn, *fn, mn, nfev};
        }
        return std::unexpected(Fault{errc::line_search_failed, {}});
    }
    template<class V> constexpr nd_view<V> view(const nd_state<V>& s) const { return {s.x, s.fx}; }
    template<class V> constexpr system_estimate<V> estimate(const nd_state<V>& s) const { return {s.x, s.fx, s.merit}; }
    template<class V> constexpr system_estimate<V> best(const nd_state<V>& s) const { return estimate(s); }
    template<class V> constexpr std::optional<stop_reason> intrinsic(const nd_state<V>& s) const noexcept {
        if (s.merit == scalar_of<V>(0)) return stop_reason::exact_zero;
        return std::nullopt;
    }
};
newton() -> newton<>;

template<class V> struct box {   // per-iterate projection onto [lo, hi] componentwise (FlashHS box limits)
    V lo, hi;
    constexpr V operator()(V x) const {
        for (std::size_t i = 0; i < traits<V>::size(x); ++i) x[i] = x[i] < lo[i] ? lo[i] : (hi[i] < x[i] ? hi[i] : x[i]);
        return x;
    }
};

}    // namespace nxx::multiroots

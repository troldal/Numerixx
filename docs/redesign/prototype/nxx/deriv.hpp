// Prototype of PLAN_v1 7.1 (deriv) and 6.12 (derivative_of): stencils as integer data, steps as powers of two.
#pragma once
#include "core.hpp"

namespace nxx::deriv {

template<int Order, int Accuracy, std::size_t Points> struct stencil {
    static constexpr int order = Order, accuracy = Accuracy;
    std::array<int, Points> offset, weight;
    int denominator;
};
inline constexpr stencil<1, 2, 2> central_1_2{{-1, 1}, {-1, 1}, 2};
inline constexpr stencil<1, 4, 4> central_1_4{{-2, -1, 1, 2}, {1, -8, 8, -1}, 12};
inline constexpr stencil<2, 2, 3> central_2_2{{-1, 0, 1}, {1, -2, 1}, 1};
inline constexpr stencil<1, 1, 2> forward_1_1{{0, 1}, {-1, 1}, 1};

// step<T>: h = factor * max(|x|, floor), then h = (x + h) - x. optimal() is resolved per stencil at the call,
// because step<T> alone does not know the stencil (the plan's `step::optimal()` cannot be computed eagerly).
template<real T> class step {
    enum class kind : std::uint8_t { optimal, relative, absolute };
    kind k_ = kind::optimal;
    T factor_ = T(0), floor_ = T(1);
    constexpr step(kind k, T f, T fl) noexcept : k_(k), factor_(f), floor_(fl) {}
public:
    constexpr step() = default;
    static constexpr step optimal() noexcept { return {}; }
    static constexpr step relative(tolerance<primal_t<T>> factor, tolerance<primal_t<T>> floor = tolerance<primal_t<T>>{1.0}) noexcept {
        return step{kind::relative, T(factor.value()), T(floor.value())};
    }
    static constexpr step absolute(tolerance<primal_t<T>> h) noexcept { return step{kind::absolute, T(h.value()), T(0)}; }
    template<int O, int A> constexpr T resolve(const T& x) const noexcept {
        T h = k_ == kind::absolute ? factor_
            : (k_ == kind::optimal ? math::root_eps<T>(1, O + A) : factor_) * math::max(math::abs(x), floor_);
        const T xh = x + h;
        return xh - x;
    }
};

// Plan 7.1 writes `diff(f, x, const stencil<O,A,N>& s = central_1_2, ...)`: the default cannot deduce O, A, N,
// so `diff(f, x)` fails to compile (see neg/neg_diff_default.cpp). Default template arguments fix it.
template<class F, real T, int O = 1, int A = 2, std::size_t N = 2>
constexpr auto diff(const F& f, T x, const stencil<O, A, N>& s = central_1_2, step<T> h = step<T>::optimal())
    -> std::expected<T, failure<T, callback_error_t<F, T>>> {
    using Fail = failure<T, callback_error_t<F, T>>;
    const T hh = h.template resolve<O, A>(x);
    T acc(0);
    for (std::size_t i = 0; i < N; ++i) {
        if (s.weight[i] == 0) continue;
        auto y = evaluate(f, x + T(s.offset[i]) * hh);
        if (!y) return std::unexpected(Fail{y.error().code, algo::none, {0, static_cast<std::uint32_t>(i + 1)}, std::nullopt, y.error().cause});
        acc += T(s.weight[i]) * *y;
    }
    T den = T(s.denominator);
    for (int k = 0; k < O; ++k) den *= hh;
    return acc / den;
}

struct optimal_step {};   // untyped "optimal": derivative_of does not know T until it is called

template<class F, class S = std::remove_cvref_t<decltype(central_1_2)>, class H = optimal_step>
class derivative_of_t {
    F f_;
    S s_;
    H h_;
public:
    constexpr derivative_of_t(F f, S s, H h) : f_(std::move(f)), s_(s), h_(h) {}
    template<real T>
    constexpr auto operator()(const T& x) const {   // const call (fixes dev-reorg); returns expected, never throws
        if constexpr (std::is_same_v<H, optimal_step>) return diff(f_, x, s_, step<T>::optimal());
        else return diff(f_, x, s_, h_);
    }
};
template<class F, class S = std::remove_cvref_t<decltype(central_1_2)>, class H = optimal_step>
constexpr auto derivative_of(F f, S s = central_1_2, H h = {}) { return derivative_of_t<F, S, H>{std::move(f), s, h}; }

// A derivative *policy*: a function of f, so a curried chain (built before f is known) can still say
// "Newton with a numeric derivative". roots::newton recognises anything with .bind(f).
template<class S = std::remove_cvref_t<decltype(central_1_2)>, class H = optimal_step>
struct numeric {
    S s = central_1_2;
    H h = {};
    template<class F> constexpr auto bind(const F& f) const { return derivative_of(f, s, h); }
};

}    // namespace nxx::deriv

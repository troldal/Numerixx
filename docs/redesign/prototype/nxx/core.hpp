// Prototype of PLAN_v1 section 6 (core abstractions). Written from the plan text, not from the spikes.
// Includes only the standard library (plan D18).
#pragma once
#include <array>
#include <cmath>
#include <compare>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <functional>
#include <iterator>
#include <limits>
#include <optional>
#include <ranges>
#include <type_traits>
#include <utility>

// ---------------------------------------------------------------- config.hpp (plan 5.3)
#if defined(_MSC_VER)
#  define NXX_NO_UNIQUE_ADDRESS [[msvc::no_unique_address]]
#else
#  define NXX_NO_UNIQUE_ADDRESS [[no_unique_address]]
#endif

// Plan 5.3: "= delete(reason) where supported, wrapped in -Wc++26-extensions suppression".
#if defined(__cpp_deleted_function) && __cpp_deleted_function >= 202403L
#  if defined(__clang__) && !defined(NXX_NO_DELETE_PRAGMA)
#    define NXX_DELETE(reason)                                                                     \
        _Pragma("clang diagnostic push") _Pragma("clang diagnostic ignored \"-Wc++26-extensions\"") \
            = delete (reason) _Pragma("clang diagnostic pop")
#  else
#    define NXX_DELETE(reason) = delete (reason)
#  endif
#else
#  define NXX_DELETE(reason) = delete
#endif

namespace nxx {

// ---------------------------------------------------------------- 6.1 scalars
template<class T> struct scalar_traits {};
template<class T>
    requires(std::numeric_limits<T>::is_specialized && !std::numeric_limits<T>::is_integer &&
             !std::numeric_limits<T>::is_exact)
struct scalar_traits<T> {
    using primal_type = T;
    static constexpr bool is_real = true;
    static constexpr T epsilon() noexcept { return std::numeric_limits<T>::epsilon(); }
    static constexpr T primal(const T& x) noexcept { return x; }
};
template<class T>
inline constexpr bool is_real_v = requires { requires scalar_traits<std::remove_cvref_t<T>>::is_real; };
template<class T>
concept real = is_real_v<T> && std::regular<T> && std::totally_ordered<T> && requires(const T a, const T b) {
    { a + b } -> std::convertible_to<T>;
    { a - b } -> std::convertible_to<T>;
    { a * b } -> std::convertible_to<T>;
    { a / b } -> std::convertible_to<T>;
    { -a } -> std::convertible_to<T>;
};
template<real T> using primal_t = typename scalar_traits<T>::primal_type;

namespace math {
namespace adl {
    void abs() = delete;
    void isfinite() = delete;
    void sqrt() = delete;
    struct abs_fn {
        template<class T>
        constexpr T operator()(const T& x) const noexcept {
            if constexpr (std::is_floating_point_v<T>) {
                if consteval { return x < T(0) ? -x : x; } else { return std::fabs(x); }
            } else if constexpr (requires { { abs(x) } -> std::convertible_to<T>; }) {
                return abs(x);
            } else {
                return x < T(0) ? -x : x;
            }
        }
    };
    struct isfinite_fn {
        template<class T>
        constexpr bool operator()(const T& x) const noexcept {
            if constexpr (std::is_floating_point_v<T>) {
                if consteval {   // comparisons only (plan 6.1)
                    return x == x && x <= std::numeric_limits<T>::max() && x >= std::numeric_limits<T>::lowest();
                } else {
                    return std::isfinite(x);
                }
            } else if constexpr (requires { { isfinite(x) } -> std::convertible_to<bool>; }) {
                return isfinite(x);
            } else {
                return x == x;
            }
        }
    };
    struct sqrt_fn {   // run time only; the plan forbids <cmath> on constexpr paths
        template<class T>
        T operator()(const T& x) const noexcept {
            if constexpr (std::is_floating_point_v<T>) return std::sqrt(x);
            else return sqrt(x);
        }
    };
}    // namespace adl
inline constexpr adl::abs_fn abs{};
inline constexpr adl::isfinite_fn isfinite{};
inline constexpr adl::sqrt_fn sqrt{};
template<class T> constexpr T max(const T& a, const T& b) noexcept { return a < b ? b : a; }
template<class T> constexpr T min(const T& a, const T& b) noexcept { return b < a ? b : a; }
template<class P> constexpr P pow2(int e) noexcept {   // exact
    P r(1);
    if (e >= 0) for (int i = 0; i < e; ++i) r *= P(2);
    else        for (int i = 0; i < -e; ++i) r /= P(2);
    return r;
}
template<class P> constexpr P root_eps(int num, int den) noexcept {   // eps^(num/den) rounded to a power of two
    return pow2<P>(-((std::numeric_limits<P>::digits - 1) * num) / den);
}
}    // namespace math

// ---------------------------------------------------------------- 6.3 errors and results
enum class errc : std::uint8_t {
    invalid_input = 1, no_sign_change, non_finite_input, dimension_mismatch, out_of_domain,
    budget_exhausted = 32, evaluations_exhausted, stalled, non_finite_value, zero_derivative, singular,
    line_search_failed, bracket_collapsed, diverged,
    callback_failed = 64
};
constexpr bool is_input_error(errc e) noexcept { return std::to_underlying(e) < 32; }
constexpr bool is_budget_error(errc e) noexcept {
    return e == errc::budget_exhausted || e == errc::evaluations_exhausted;
}
enum class algo : std::uint8_t { none, bisection, brent, secant, newton, expand_out, multi_newton, user = 200 };
enum class stop_reason : std::uint8_t { exact_zero, criterion, resolution_limit, algorithm };

struct counters {
    std::uint32_t iterations = 0, evaluations = 0;
    friend constexpr counters operator+(counters a, counters b) noexcept {
        return {a.iterations + b.iterations, a.evaluations + b.evaluations};
    }
    friend constexpr bool operator==(counters, counters) = default;
};

struct none { friend constexpr bool operator==(none, none) = default; };
template<class UE> using cause_slot = std::conditional_t<std::is_same_v<UE, none>, none, std::optional<UE>>;

template<class UE = none> struct fault {   // per-evaluation / per-step error
    errc code{};
    NXX_NO_UNIQUE_ADDRESS cause_slot<UE> cause{};
};

template<class Est> struct solution : Est {
    counters used{};
    algo by = algo::none;
    stop_reason how{};
};
template<class Est, class UE = none> struct failure {
    using estimate_type = Est;
    using cause_type = UE;
    errc code{};
    algo where = algo::none;
    counters used{};
    std::optional<Est> best{};
    NXX_NO_UNIQUE_ADDRESS cause_slot<UE> cause{};
};
template<class Est, class UE = none> using result = std::expected<solution<Est>, failure<Est, UE>>;

template<class R> constexpr auto best_x(const R& r) {   // optional: the solution's x, or the failure's best->x
    using X = decltype(r->x);
    if (r) return std::optional<X>(r->x);
    if (r.error().best) return std::optional<X>(r.error().best->x);
    return std::optional<X>{};
}

// ---------------------------------------------------------------- 6.2 refined types
namespace detail {
    inline void literal_violates_invariant(const char*) noexcept {}   // NOT constexpr: reaching it = compile error
    struct trust_me { explicit constexpr trust_me() = default; };

    template<class Tag, class T> class refined {
        T v_;
    public:
        using value_type = T;
        constexpr refined(trust_me, T v) noexcept : v_(v) {}
        consteval refined(T v) : v_(v) {   // invalid literal: compile error
            if (!Tag::check(v)) literal_violates_invariant(Tag::message);
        }
        template<class B> requires std::same_as<B, bool>
        refined(B) NXX_DELETE("a bool is not a numeric refinement");
        static constexpr auto make(T v) noexcept -> std::expected<refined, errc> {
            if (!Tag::check(v)) return std::unexpected(errc::invalid_input);
            return refined{trust_me{}, v};
        }
        constexpr T value() const noexcept { return v_; }
        friend constexpr auto operator<=>(const refined&, const refined&) = default;
    };
}    // namespace detail

namespace tag {
    struct abs_tolerance {
        static constexpr const char* message = "tolerance must be finite and > 0";
        template<class T> static constexpr bool check(const T& v) { return v > T(0) && math::isfinite(v); }
    };
    struct rel_tolerance {
        static constexpr const char* message = "relative tolerance must be finite and in [0, 1)";
        template<class T> static constexpr bool check(const T& v) { return v >= T(0) && v < T(1); }
    };
}    // namespace tag
template<real T> using tolerance     = detail::refined<tag::abs_tolerance, T>;
template<real T> using rel_tolerance = detail::refined<tag::rel_tolerance, T>;

class max_iterations {   // 1 .. 2^32-1
    std::uint32_t n_;
    constexpr max_iterations(detail::trust_me, std::uint32_t n) noexcept : n_(n) {}
public:
    template<std::integral I> requires(!std::same_as<I, bool>)
    consteval max_iterations(I n) : n_(static_cast<std::uint32_t>(n)) {
        if (n < I(1) || static_cast<unsigned long long>(n) > 0xFFFF'FFFFull)
            detail::literal_violates_invariant("max_iterations must be in [1, 2^32)");
    }
    template<class B> requires std::same_as<B, bool>
    max_iterations(B) NXX_DELETE("a bool is not an iteration count");
    static constexpr auto make(long long n) noexcept -> std::expected<max_iterations, errc> {
        if (n < 1 || n > 0xFFFF'FFFFll) return std::unexpected(errc::invalid_input);
        return max_iterations{detail::trust_me{}, static_cast<std::uint32_t>(n)};
    }
    constexpr std::uint32_t value() const noexcept { return n_; }
    friend constexpr bool operator==(max_iterations, max_iterations) = default;
};

template<real T> class bracket {   // finite, lo < hi
    T lo_, hi_;
public:
    constexpr bracket(detail::trust_me, T lo, T hi) noexcept : lo_(lo), hi_(hi) {}
    consteval bracket(T lo, T hi) : lo_(lo), hi_(hi) {
        if (!(lo < hi) || !math::isfinite(lo) || !math::isfinite(hi))
            detail::literal_violates_invariant("bracket requires finite lo < hi");
    }
    static constexpr auto make(T a, T b) noexcept -> std::expected<bracket, errc> {
        if (!math::isfinite(a) || !math::isfinite(b) || a == b) return std::unexpected(errc::invalid_input);
        return a < b ? bracket{detail::trust_me{}, a, b} : bracket{detail::trust_me{}, b, a};
    }
    constexpr T lo() const noexcept { return lo_; }
    constexpr T hi() const noexcept { return hi_; }
    constexpr T width() const noexcept { return hi_ - lo_; }
    friend constexpr bool operator==(const bracket&, const bracket&) = default;
};
template<real T> bracket(T, T) -> bracket<T>;

// ---------------------------------------------------------------- 6.4 callables and evaluation
namespace detail {
    template<class R> struct unwrap { using value = R; using error = none; static constexpr int kind = 0; };
    template<class V, class E> struct unwrap<std::expected<V, E>> {
        using value = V; using error = E; static constexpr int kind = 1; };
    template<class V, class UE> struct unwrap<std::expected<V, fault<UE>>> {
        using value = V; using error = UE; static constexpr int kind = 2; };
#if defined(NXX_FIX_UNWRAP_FAILURE)   // prototype fix: also unwrap failure<Est, UE> (what 6.12 says derivative_of returns)
    template<class V, class Est, class UE> struct unwrap<std::expected<V, failure<Est, UE>>> {
        using value = V; using error = UE; static constexpr int kind = 3; };
#endif
    template<class F> constexpr decltype(auto) unref(const F& f) noexcept {
        if constexpr (requires { typename F::type; f.get(); }) return f.get(); else return (f);
    }
}    // namespace detail
template<class F, class X>
using callback_error_t = typename detail::unwrap<std::remove_cvref_t<std::invoke_result_t<const F&, const X&>>>::error;

template<class X, class F>
constexpr auto evaluate(const F& f, const X& x) noexcept(std::is_nothrow_invocable_v<const F&, const X&>)
    -> std::expected<X, fault<callback_error_t<F, X>>> {
    using U = detail::unwrap<std::remove_cvref_t<std::invoke_result_t<const F&, const X&>>>;
    using UE = typename U::error;
    auto finite = [](const X& y) {
        if constexpr (real<X>) return math::isfinite(y);
        else return true;   // vector types: the caller checks
    };
    if constexpr (U::kind == 0) {
        X y = std::invoke(f, x);
        if (!finite(y)) return std::unexpected(fault<UE>{errc::non_finite_value, {}});
        return y;
    } else {
        auto r = std::invoke(f, x);
        if (!r) {
            if constexpr (U::kind == 1) return std::unexpected(fault<UE>{errc::callback_failed, cause_slot<UE>(r.error())});
            else if constexpr (U::kind == 2) return std::unexpected(r.error());
            else return std::unexpected(fault<UE>{r.error().code, r.error().cause});
        }
        X y = *std::move(r);
        if (!finite(y)) return std::unexpected(fault<UE>{errc::non_finite_value, {}});
        return y;
    }
}

// Combining the error types of two callbacks (f and f') -- the plan is silent on this.
namespace detail {
    template<class A, class B> struct common_cause {
        static_assert(std::is_same_v<A, B>, "nxx: f and its derivative report different error types; "
                                            "map one with .transform_error so they agree");
        using type = A;
    };
    template<class A> struct common_cause<A, none> { using type = A; };
    template<class B> struct common_cause<none, B> { using type = B; };
    template<> struct common_cause<none, none> { using type = none; };
}
template<class A, class B> using common_cause_t = typename detail::common_cause<A, B>::type;

// ---------------------------------------------------------------- 6.8 stop criteria
enum class verdict : std::uint8_t { proceed, converged, stalled, exhausted };
struct criterion_base;
template<class T> inline constexpr bool is_criterion_v = std::is_base_of_v<criterion_base, T>;
template<class A, class B> struct any_of_t;
template<class A, class B> struct all_of_t;

struct criterion_base {   // hidden friends, bool-variable-template constraints (clang-cl safe)
    template<class A, class B> requires(is_criterion_v<A> && is_criterion_v<B>)
    friend constexpr auto operator||(const A& a, const B& b) noexcept { return any_of_t<A, B>{{}, a, b}; }
    template<class A, class B> requires(is_criterion_v<A> && is_criterion_v<B>)
    friend constexpr auto operator&&(const A& a, const B& b) noexcept { return all_of_t<A, B>{{}, a, b}; }
};

// Views expose: x(), fx(), distance(prev) (primal), residual() (primal), scale() (primal); bracketing
// views add enclosure(). Criteria only talk to that protocol, so they are shared by 1-D and N-D.
template<class A, class B> struct any_of_t : criterion_base {   // first non-proceed wins
    A a; B b;
    template<class V> constexpr verdict operator()(const V& p, const V& n, counters c) const {
        if (auto v = a(p, n, c); v != verdict::proceed) return v;
        return b(p, n, c);
    }
};
template<class A, class B> struct all_of_t : criterion_base {   // both converge; failure wins
    A a; B b;
    template<class V> constexpr verdict operator()(const V& p, const V& n, counters c) const {
        const auto va = a(p, n, c), vb = b(p, n, c);
        if (va == verdict::stalled || va == verdict::exhausted) return va;
        if (vb == verdict::stalled || vb == verdict::exhausted) return vb;
        return (va == verdict::converged && vb == verdict::converged) ? verdict::converged : verdict::proceed;
    }
};

template<real T> struct x_tol : criterion_base {   // |dx| <= abs + rel*|x|
    tolerance<T> abs;
    rel_tolerance<T> rel = rel_tolerance<T>{detail::trust_me{}, T(0)};
    constexpr x_tol(tolerance<T> a) noexcept : abs(a) {}
    constexpr x_tol(tolerance<T> a, rel_tolerance<T> r) noexcept : abs(a), rel(r) {}
    constexpr T abs_tolerance() const noexcept { return abs.value(); }
    template<class V> constexpr verdict operator()(const V& p, const V& n, counters) const {
        return n.distance(p) <= abs.value() + rel.value() * n.scale() ? verdict::converged : verdict::proceed;
    }
};
template<real T> x_tol(T) -> x_tol<T>;
template<real T> x_tol(T, T) -> x_tol<T>;

template<real T> struct f_tol : criterion_base {   // |f| <= abs
    tolerance<T> abs;
    constexpr f_tol(tolerance<T> a) noexcept : abs(a) {}
    template<class V> constexpr verdict operator()(const V&, const V& n, counters) const {
        return n.residual() <= abs.value() ? verdict::converged : verdict::proceed;
    }
};
template<real T> f_tol(T) -> f_tol<T>;

template<real T> struct width_tol : criterion_base {   // compiles only for views with enclosure()
    tolerance<T> abs;
    constexpr width_tol(tolerance<T> a) noexcept : abs(a) {}
    template<class V> constexpr verdict operator()(const V&, const V& n, counters) const {
        static_assert(requires { n.enclosure(); },
                      "nxx::width_tol needs a bracketing method (the view has no enclosure()); use x_tol");
        if constexpr (requires { n.enclosure(); })
            return n.enclosure().width() <= abs.value() ? verdict::converged : verdict::proceed;
        else return verdict::proceed;
    }
};
template<real T> width_tol(T) -> width_tol<T>;

struct floored_width : criterion_base {   // absolute floor (roots at 0 terminate); default for bracketing methods
    template<class V> constexpr verdict operator()(const V&, const V& n, counters) const {
        const auto e = n.enclosure();
        using T = std::remove_cvref_t<decltype(e.lo())>;
        const T eps = std::numeric_limits<T>::epsilon();
        const T floor = T(4) * eps * math::max(T(1), math::min(math::abs(e.lo()), math::abs(e.hi())));
        return e.width() <= floor ? verdict::converged : verdict::proceed;
    }
};
struct default_step : criterion_base {   // |dx| <= 8 eps max(1, |x|), default for open methods
    template<class V> constexpr verdict operator()(const V& p, const V& n, counters) const {
        using T = std::remove_cvref_t<decltype(n.scale())>;
        const T eps = std::numeric_limits<T>::epsilon();
        return n.distance(p) <= T(8) * eps * math::max(T(1), n.scale()) ? verdict::converged : verdict::proceed;
    }
};
struct min_iterations : criterion_base {   // use with &&
    std::uint32_t n;
    constexpr min_iterations(max_iterations m) noexcept : n(m.value()) {}
    template<class V> constexpr verdict operator()(const V&, const V&, counters c) const {
        return c.iterations >= n ? verdict::converged : verdict::proceed;
    }
};
struct max_evaluations : criterion_base {
    std::uint32_t n;
    constexpr max_evaluations(max_iterations m) noexcept : n(m.value()) {}
    template<class V> constexpr verdict operator()(const V&, const V&, counters c) const {
        return c.evaluations >= n ? verdict::exhausted : verdict::proceed;
    }
};
struct never : criterion_base {
    template<class V> constexpr verdict operator()(const V&, const V&, counters) const noexcept { return verdict::proceed; }
};

// ---------------------------------------------------------------- problem, better(), observer
template<class F, class In> struct problem {   // f is a std::reference_wrapper for one-shot calls
    F f;
    In in;
    std::uint32_t nfev0 = 0;
};
struct no_observer { template<class V> constexpr void operator()(const V&) const noexcept {} };

namespace detail {
    template<class Est> constexpr bool better(const Est& e, const Est& best) {   // per family, by ADL hook
        return merit_of(e) < merit_of(best);
    }
    template<class A, class P> using state_t = typename decltype(std::declval<const A&>().init(std::declval<const P&>()))::value_type;
    template<class A, class P> using init_failure_t = typename decltype(std::declval<const A&>().init(std::declval<const P&>()))::error_type;
}

// ---------------------------------------------------------------- 6.6 solver protocol as a concept
template<class A, class P>
concept iterative_solver_for = requires(const A& a, const P& p, const detail::state_t<A, P>& s) {
    { A::id } -> std::convertible_to<algo>;
    a.init(p);
    a.step(p, s);
    a.view(s);
    a.estimate(s);
    a.best(s);
    { a.intrinsic(s) } -> std::same_as<std::optional<stop_reason>>;
    { s.nfev } -> std::convertible_to<std::uint32_t>;
};

// ---------------------------------------------------------------- 6.7 the single bounded-iteration driver
// Deviation from the plan's signature: success estimate (estimate) and failure estimate (best) may differ,
// because searchers succeed with a sign_bracket but fail with a root_estimate (plan 6.3).
template<class A, class P, class Stop, class Obs = no_observer>
    requires iterative_solver_for<A, P>
constexpr auto iterate(const A& alg, const P& p, const Stop& stop, max_iterations budget, const Obs& observe = {}) {
    using Fail = detail::init_failure_t<A, P>;
    using FEst = typename Fail::estimate_type;
    using SEst = std::remove_cvref_t<decltype(alg.estimate(std::declval<const detail::state_t<A, P>&>()))>;
    using R = std::expected<solution<SEst>, Fail>;
    auto first = alg.init(p);
    if (!first) return R{std::unexpect, std::move(first).error()};
    auto s = *std::move(first);
    FEst best = alg.best(s);
    if (auto how = alg.intrinsic(s)) return R{solution<SEst>{alg.estimate(s), {0, s.nfev}, A::id, *how}};
    for (std::uint32_t k = 1; k <= budget.value(); ++k) {
        auto next = alg.step(p, s);
        if (!next) return R{std::unexpect, Fail{next.error().code, A::id, {k, s.nfev}, best, next.error().cause}};
        if (const FEst e = alg.best(*next); detail::better(e, best)) best = e;   // best iterate on EVERY exit
        std::invoke(observe, alg.view(*next));
        if (auto how = alg.intrinsic(*next))
            return R{solution<SEst>{alg.estimate(*next), {k, next->nfev}, A::id, *how}};
        switch (stop(alg.view(s), alg.view(*next), counters{k, next->nfev})) {
            case verdict::converged:
                return R{solution<SEst>{alg.estimate(*next), {k, next->nfev}, A::id, stop_reason::criterion}};
            case verdict::stalled:   return R{std::unexpect, Fail{errc::stalled, A::id, {k, next->nfev}, best, {}}};
            case verdict::exhausted: return R{std::unexpect, Fail{errc::evaluations_exhausted, A::id, {k, next->nfev}, best, {}}};
            case verdict::proceed:   break;
        }
        s = *std::move(next);   // the only mutation: a local
    }
    return R{std::unexpect, Fail{errc::budget_exhausted, A::id, {budget.value(), s.nfev}, best, {}}};
}

// ---------------------------------------------------------------- solver facade (deducing this, no CRTP)
template<class S, class In> struct bound {   // solver.on(input): a value f -> result
    S solver;
    In in;
    template<class F> constexpr auto operator()(const F& f) const { return solver(f, in); }
};

struct solver_facade {
    template<class Self, class F, class In>
        requires requires(const Self& s, const F& f, const In& in) { s.prepare(std::cref(f), in); }
    constexpr auto operator()(this const Self& self, const F& f, const In& in) {
        auto p = self.prepare(std::cref(f), in);
        using R = decltype(iterate(self, *p, self.stop(), self.budget()));
        if (!p) return R{std::unexpect, std::move(p).error()};
        return iterate(self, *p, self.stop(), self.budget());
    }
    template<class Self, class In>
    constexpr auto on(this const Self& self, In in) { return bound<Self, In>{self, std::move(in)}; }
};

// ---------------------------------------------------------------- 6.9 lazy iteration: steps_view
// An input range of expected<state, fault>. The first element is init (its failure mapped to a fault);
// it ends after the first error or intrinsic stop, otherwise it is infinite: bound it with views::take.
template<class A, class P>
class steps_view : public std::ranges::view_interface<steps_view<A, P>> {
    using S = detail::state_t<A, P>;
    using UE = typename detail::init_failure_t<A, P>::cause_type;
public:
    using element = std::expected<S, fault<UE>>;
    constexpr steps_view(A a, P p) : alg_(std::move(a)), p_(std::move(p)) {}

    class iterator {
        const steps_view* parent_ = nullptr;
        std::optional<element> cur_{};
        bool done_ = false;
        friend class steps_view;
    public:
        using value_type = element;
        using difference_type = std::ptrdiff_t;
        using iterator_concept = std::input_iterator_tag;
        iterator() = default;
        constexpr const element& operator*() const { return *cur_; }
        constexpr iterator& operator++() {
            if (!*cur_ || parent_->alg_.intrinsic(**cur_)) done_ = true;
            else cur_ = parent_->alg_.step(parent_->p_, **cur_);
            return *this;
        }
        constexpr void operator++(int) { ++*this; }
        friend constexpr bool operator==(const iterator& it, std::default_sentinel_t) { return it.done_; }
    };
    constexpr iterator begin() const {
        iterator it;
        it.parent_ = this;
        auto s = alg_.init(p_);
        if (s) it.cur_.emplace(*std::move(s));
        else it.cur_.emplace(std::unexpect, fault<UE>{s.error().code, s.error().cause});
        return it;
    }
    constexpr std::default_sentinel_t end() const noexcept { return {}; }
private:
    A alg_;
    P p_;
};

// ---------------------------------------------------------------- 6.10 combinators
namespace detail {
    template<class Fail> constexpr Fail merge(const Fail& a, Fail b) {   // last code/cause, BEST estimate, total cost
        b.used = a.used + b.used;
        if (a.best && (!b.best || better(*a.best, *b.best))) b.best = a.best;
        return b;
    }
    template<class R> constexpr R add_cost(R r, counters extra) {
        if (r) r->used = r->used + extra; else r.error().used = r.error().used + extra;
        return r;
    }
}

template<class S> constexpr auto first_of(S s) { return s; }
template<class S1, class S2, class... Ss>
constexpr auto first_of(S1 s1, S2 s2, Ss... ss) {   // unconstrained variadic: clang-cl safe
    return [s1, rest = first_of(std::move(s2), std::move(ss)...)]<class... Args>(const Args&... a) {
        using R1 = std::invoke_result_t<const S1&, const Args&...>;
        using R2 = std::invoke_result_t<const decltype(rest)&, const Args&...>;
        if constexpr (!std::is_same_v<R1, R2>) {
            static_assert(std::is_same_v<R1, R2>, "nxx::first_of: every alternative must return the same "
                          "std::expected<solution<Est>, failure<Est, UE>>; adapt the odd one with .transform/.transform_error");
            return R1{std::unexpect};
        } else {
            auto r1 = std::invoke(s1, a...);
            if (r1) return r1;   // lazy: later alternatives never run
            auto r2 = std::invoke(rest, a...);
            if (r2) { r2->used = r2->used + r1.error().used; return r2; }   // success pays for failed attempts
            return R1{std::unexpect, detail::merge(r1.error(), std::move(r2).error())};
        }
    };
}

template<class S1, class S2>
constexpr auto then(S1 s1, S2 s2) {   // Kleisli with environment: stage 2 gets f AND stage 1's value
    return [s1, s2]<class F>(const F& f) {
        auto r1 = s1(f);
#if defined(NXX_THEN_CONTRACT)   // the plan's 6.10 code has no contract check; this is the proposed fix
        using V1 = typename decltype(r1)::value_type;
        if constexpr (!std::is_invocable_v<const S2&, const F&, const V1&>) {
            static_assert(std::is_invocable_v<const S2&, const F&, const V1&>,
                          "nxx::then: stage 2 cannot start from stage 1's result (a bracketing solver needs a bracket or a "
                          "search result; use brent{}.from_enclosure() or put a searcher first)");
            return r1;
        } else
#endif
        {
            using R2 = decltype(s2(f, *r1));
            if (!r1) return R2{std::unexpect, std::move(r1).error()};
            return detail::add_cost(s2(f, *r1), r1->used);
        }
    };
}
template<class S1, class S2, class S3, class... Ss>
constexpr auto then(S1 s1, S2 s2, S3 s3, Ss... ss) { return then(then(std::move(s1), std::move(s2)), std::move(s3), std::move(ss)...); }

template<class S1, class S2>
constexpr auto warm_fallback(S1 s1, S2 s2) {   // restart s2 from s1's failure best
    return [s1, s2]<class F>(const F& f) {
        auto r1 = s1(f);
        if (r1 || !r1.error().best) return r1;
        auto r2 = s2(f, *r1.error().best);
        static_assert(std::is_same_v<decltype(r1), decltype(r2)>, "nxx::warm_fallback: both stages must return the same result type");
        if (r2) { r2->used = r2->used + r1.error().used; return r2; }
        return decltype(r1){std::unexpect, detail::merge(r1.error(), std::move(r2).error())};
    };
}

}    // namespace nxx

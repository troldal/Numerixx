// Stop criteria (DESIGN §6.8): pure functions of (previous view, next view, counters). Each criterion declares the view
// kinds it applies to, so a solver can reject one that makes no sense for it (x_tol on a bracketing method) when it is
// configured, with a reason.
//
// Views expose x(), fx(), residual() and scale(). Point and system views add distance(prev); enclosure views add
// enclosure() (with lo() and hi()) but not distance(), because successive iterates say nothing about an enclosure.
#pragma once

#include <numerixx/config.hpp>
#include <numerixx/core/error.hpp>
#include <numerixx/core/math.hpp>
#include <numerixx/core/refined.hpp>
#include <numerixx/core/scalar.hpp>

#include <algorithm>
#include <cstdint>
#include <expected>
#include <limits>
#include <type_traits>
#include <utility>

NXX_BEGIN_HEADER

namespace nxx
{
    enum class verdict : std::uint8_t { proceed, converged, stalled, exhausted };

    // A bit set: point = open methods, enclosure = bracketing methods, system = N-D solvers.
    enum class view_kind : std::uint8_t { point = 1, enclosure = 2, system = 4 };

    constexpr view_kind operator|(view_kind a, view_kind b) noexcept
    { return static_cast<view_kind>(std::to_underlying(a) | std::to_underlying(b)); }
    constexpr view_kind operator&(view_kind a, view_kind b) noexcept
    { return static_cast<view_kind>(std::to_underlying(a) & std::to_underlying(b)); }
    inline constexpr view_kind all_views = view_kind::point | view_kind::enclosure | view_kind::system;

    struct criterion_base;

    template<class T>
    inline constexpr bool is_criterion_v = std::is_base_of_v<criterion_base, std::remove_cvref_t<T>>;

    // Whether criterion C applies to a solver whose views are of kind K.
    template<class C, view_kind K>
    inline constexpr bool criterion_for_v = [] {
        if constexpr (is_criterion_v<C>)
            return std::to_underlying(std::remove_cvref_t<C>::applies_to & K) != 0;
        else
            return false;
    }();

    template<class A, class B>
    struct any_of_t;
    template<class A, class B>
    struct all_of_t;

    // Hidden friends with bool-variable-template constraints: clang-cl safe (DESIGN §5.3).
    struct criterion_base
    {
        using nxx_refined_tag = void;

        template<class A, class B>
            requires(is_criterion_v<A> && is_criterion_v<B>)
        friend constexpr auto operator||(const A& a, const B& b) noexcept
        { return any_of_t<A, B> { {}, a, b }; }

        template<class A, class B>
            requires(is_criterion_v<A> && is_criterion_v<B>)
        friend constexpr auto operator&&(const A& a, const B& b) noexcept
        { return all_of_t<A, B> { {}, a, b }; }
    };

    // a || b: the first verdict that is not "proceed" wins.
    template<class A, class B>
    struct any_of_t : criterion_base
    {
        static constexpr view_kind applies_to = A::applies_to & B::applies_to;
        A                          a;
        B                          b;

        template<class V>
        constexpr verdict operator()(const V& prev, const V& next, counters c) const
        {
            if (const verdict v = a(prev, next, c); v != verdict::proceed) return v;
            return b(prev, next, c);
        }
    };

    // a && b: converged only when both converge; a failure verdict of either wins.
    template<class A, class B>
    struct all_of_t : criterion_base
    {
        static constexpr view_kind applies_to = A::applies_to & B::applies_to;
        A                          a;
        B                          b;

        template<class V>
        constexpr verdict operator()(const V& prev, const V& next, counters c) const
        {
            const verdict va = a(prev, next, c);
            const verdict vb = b(prev, next, c);
            if (va == verdict::stalled || va == verdict::exhausted) return va;
            if (vb == verdict::stalled || vb == verdict::exhausted) return vb;
            return (va == verdict::converged && vb == verdict::converged) ? verdict::converged : verdict::proceed;
        }
    };

    namespace detail
    {
        template<class T>
        constexpr bool mixed_tolerance_ok(const T& abs, const T& rel) noexcept
        { return tag::abs_tolerance::check(abs) && tag::rel_tolerance::check(rel) && (abs > T(0) || rel > T(0)); }
    }    // namespace detail

    // |x_k - x_{k-1}| <= abs + rel * |x_k|. Open methods and systems only.
    template<real T>
    class x_tol : public criterion_base
    {
        T abs_;
        T rel_;

        constexpr x_tol(detail::trust_me, T a, T r) noexcept : abs_(a), rel_(r) {}

    public:
        static constexpr view_kind applies_to = view_kind::point | view_kind::system;

        constexpr x_tol(tolerance<T> a) noexcept : abs_(a.value()), rel_(T(0)) {}

        // x_tol{0.0, 1e-8} (purely relative) is legal; x_tol{0.0, 0.0} is not.
        consteval x_tol(T a, T r) : abs_(a), rel_(r)
        {
            if (!detail::mixed_tolerance_ok(a, r))
                detail::literal_violates_invariant("x_tol needs abs >= 0, 0 <= rel < 1, and abs > 0 or rel > 0");
        }

        static constexpr auto make(T a, T r) noexcept -> std::expected<x_tol, errc>
        {
            if (!detail::mixed_tolerance_ok(a, r)) return std::unexpected(errc::invalid_input);
            return x_tol { detail::trust_me {}, a, r };
        }

        // The same, from validated role types (a configuration holds abs_tolerance and rel_tolerance, so the roles
        // cannot be swapped); only the joint invariant abs > 0 || rel > 0 is left to check.
        static constexpr auto make(abs_tolerance<T> a, rel_tolerance<T> r) noexcept -> std::expected<x_tol, errc>
        { return make(a.value(), r.value()); }

        constexpr T abs() const noexcept { return abs_; }
        constexpr T rel() const noexcept { return rel_; }
        // Computed in the scalar type of the problem, which may differ from T (a double tolerance on a float or
        // multiprecision problem).
        template<real U>
        constexpr U threshold(const U& x) const noexcept
        { return U(abs_) + U(rel_) * math::abs(x); }

        template<class V>
        constexpr verdict operator()(const V& prev, const V& next, counters) const
        { return next.distance(prev) <= threshold(next.x()) ? verdict::converged : verdict::proceed; }
    };

    template<real T>
    x_tol(T) -> x_tol<T>;
    template<real T>
    x_tol(T, T) -> x_tol<T>;

    // |dx| <= 2^-ceil(p * Num / Den) * max(|x|, 1), p = digits of T: the open-method defaults (Newton step_tol<3, 5>,
    // secant step_tol<7, 10>). The absolute floor at scale 1 lets a root at 0 terminate.
    template<int Num, int Den>
    struct step_tol : criterion_base
    {
        static_assert(Num > 0 && Den > 0 && Num < Den,
                      "nxx::step_tol<Num, Den> needs 0 < Num / Den < 1 (Num == Den asks for a step below one ulp)");
        static constexpr view_kind applies_to = view_kind::point | view_kind::system;

        template<real T>
        static constexpr T threshold(const T& x) noexcept
        {
            constexpr int e = (std::numeric_limits<T>::digits * Num + Den - 1) / Den;
            return math::pow2<T>(-e) * (std::max)(math::abs(x), T(1));
        }

        template<class V>
        constexpr verdict operator()(const V& prev, const V& next, counters) const
        { return next.distance(prev) <= threshold(next.x()) ? verdict::converged : verdict::proceed; }
    };

    // hi - lo <= abs + rel * min(|lo|, |hi|): every point of the enclosure, the returned x included, is within
    // tolerance. Bracketing methods only.
    template<real T>
    class width_tol : public criterion_base
    {
        T abs_;
        T rel_;

        constexpr width_tol(detail::trust_me, T a, T r) noexcept : abs_(a), rel_(r) {}

    public:
        static constexpr view_kind applies_to = view_kind::enclosure;

        constexpr width_tol(tolerance<T> a) noexcept : abs_(a.value()), rel_(T(0)) {}

        consteval width_tol(T a, T r) : abs_(a), rel_(r)
        {
            if (!detail::mixed_tolerance_ok(a, r))
                detail::literal_violates_invariant("width_tol needs abs >= 0, 0 <= rel < 1, and abs > 0 or rel > 0");
        }

        static constexpr auto make(T a, T r) noexcept -> std::expected<width_tol, errc>
        {
            if (!detail::mixed_tolerance_ok(a, r)) return std::unexpected(errc::invalid_input);
            return width_tol { detail::trust_me {}, a, r };
        }

        // The same, from validated role types (a configuration holds abs_tolerance and rel_tolerance, so the roles
        // cannot be swapped); only the joint invariant abs > 0 || rel > 0 is left to check.
        static constexpr auto make(abs_tolerance<T> a, rel_tolerance<T> r) noexcept -> std::expected<width_tol, errc>
        { return make(a.value(), r.value()); }

        constexpr T abs() const noexcept { return abs_; }
        constexpr T rel() const noexcept { return rel_; }

        // The threshold at a point x, and for an enclosure (min(|lo|, |hi|)), in the scalar type of the problem.
        template<real U>
        constexpr U threshold(const U& x) const noexcept
        { return U(abs_) + U(rel_) * math::abs(x); }
        template<real U>
        constexpr U threshold(const U& lo, const U& hi) const noexcept
        { return U(abs_) + U(rel_) * (std::min)(math::abs(lo), math::abs(hi)); }

        template<class V>
        constexpr verdict operator()(const V&, const V& next, counters) const
        {
            const auto e = next.enclosure();
            return e.hi() - e.lo() <= threshold(e.lo(), e.hi()) ? verdict::converged : verdict::proceed;
        }
    };

    template<real T>
    width_tol(T) -> width_tol<T>;
    template<real T>
    width_tol(T, T) -> width_tol<T>;

    // w <= max(2^(1 - bits), 4 eps) * max(1, min(|lo|, |hi|)): the default tolerance of bracketing methods. The
    // absolute floor at scale 1 makes roots at 0 terminate. bits defaults to the digits of T.
    class floored_width : public criterion_base
    {
        int bits_ = 0;    // 0: the digits of T

    public:
        static constexpr view_kind applies_to = view_kind::enclosure;

        constexpr floored_width() noexcept = default;

        consteval explicit floored_width(int bits) : bits_(bits)
        {
            if (bits < 1) detail::literal_violates_invariant("floored_width needs bits >= 1");
        }

        constexpr int bits() const noexcept { return bits_; }

        template<real T>
        constexpr T factor() const noexcept
        {
            const int b = bits_ == 0 ? std::numeric_limits<T>::digits : bits_;
            return (std::max)(math::pow2<T>(1 - b), T(4) * std::numeric_limits<T>::epsilon());
        }

        template<real T>
        constexpr T threshold(const T& x) const noexcept
        { return factor<T>() * (std::max)(T(1), math::abs(x)); }

        template<real T>
        constexpr T threshold(const T& lo, const T& hi) const noexcept
        { return factor<T>() * (std::max)(T(1), (std::min)(math::abs(lo), math::abs(hi))); }

        template<class V>
        constexpr verdict operator()(const V&, const V& next, counters) const
        {
            const auto e = next.enclosure();
            return e.hi() - e.lo() <= threshold(e.lo(), e.hi()) ? verdict::converged : verdict::proceed;
        }
    };

    // |f(x)| <= abs. Opt-in; applies to every view (minimisers, phase 4, will reject it).
    template<real T>
    class f_tol : public criterion_base
    {
        T abs_;

    public:
        static constexpr view_kind applies_to = all_views;

        constexpr f_tol(tolerance<T> a) noexcept : abs_(a.value()) {}
        constexpr T abs() const noexcept { return abs_; }

        template<class V>
        constexpr verdict operator()(const V&, const V& next, counters) const
        {
            using U = std::remove_cvref_t<decltype(next.residual())>;
            return next.residual() <= U(abs_) ? verdict::converged : verdict::proceed;
        }
    };

    template<real T>
    f_tol(T) -> f_tol<T>;

    // An evaluation budget, for expensive functions: exhausted once n evaluations are spent.
    class max_evaluations : public criterion_base
    {
        std::uint32_t n_;

    public:
        static constexpr view_kind applies_to = all_views;

        constexpr max_evaluations(evaluation_budget n) noexcept : n_(n.value()) {}
        constexpr std::uint32_t value() const noexcept { return n_; }

        template<class V>
        constexpr verdict operator()(const V&, const V&, counters c) const noexcept
        { return c.evaluations >= n_ ? verdict::exhausted : verdict::proceed; }
    };

    // Converged once n iterations have run; combine with &&.
    class min_iterations : public criterion_base
    {
        std::uint32_t n_;

    public:
        static constexpr view_kind applies_to = all_views;

        constexpr min_iterations(max_iterations n) noexcept : n_(n.value()) {}
        constexpr std::uint32_t value() const noexcept { return n_; }

        template<class V>
        constexpr verdict operator()(const V&, const V&, counters c) const noexcept
        { return c.iterations >= n_ ? verdict::converged : verdict::proceed; }
    };

    // Never stops: leaves the decision to the solver's intrinsic test and the budget.
    struct never : criterion_base
    {
        static constexpr view_kind applies_to = all_views;

        template<class V>
        constexpr verdict operator()(const V&, const V&, counters) const noexcept
        { return verdict::proceed; }
    };
}    // namespace nxx

NXX_END_HEADER

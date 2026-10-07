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
#include <concepts>
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

    namespace detail
    {
        // A guard (min_iterations) delays convergence but never establishes it: alone, or under ||, it would stop a
        // solver after n iterations and report success without any accuracy test.
        template<class C>
        inline constexpr bool guard_only_v = [] {
            if constexpr (requires { std::remove_cvref_t<C>::guard_only; })
                return bool(std::remove_cvref_t<C>::guard_only);
            else
                return false;
        }();
    }    // namespace detail

    // Whether C can stop a solver whose views are of kind K: it applies to them and is not a bare guard. Solvers
    // constrain their stop criterion on this, with reasoned deletions for the rest.
    template<class C, view_kind K>
    inline constexpr bool stop_criterion_for_v = criterion_for_v<C, K> && !detail::guard_only_v<C>;

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
        static constexpr bool      guard_only = detail::guard_only_v<A> || detail::guard_only_v<B>;    // either may stop alone
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
        static constexpr bool      guard_only = detail::guard_only_v<A> && detail::guard_only_v<B>;    // both must converge
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
        // Whether C is a width criterion (one that applies to enclosures only: width_tol, floored_width) or contains one
        // under || or &&, at any depth. A solver with its own tolerance (brent) rejects such a stop criterion: its
        // intrinsic test runs first and reports success once its own tolerance holds, so an external width criterion
        // could not tighten it (DESIGN §6.8). Each step is an if constexpr, so no operand is instantiated that the
        // answer does not need.
        template<class C>
        struct contains_width
        {
            static constexpr bool value = [] {
                if constexpr (!is_criterion_v<C>)
                    return false;
                else if constexpr (requires { C::applies_to; })
                    return C::applies_to == view_kind::enclosure;
                else
                    return false;
            }();
        };

        template<class A, class B>
        struct contains_width<any_of_t<A, B>>
        {
            static constexpr bool value = [] {
                if constexpr (contains_width<std::remove_cvref_t<A>>::value)
                    return true;
                else
                    return contains_width<std::remove_cvref_t<B>>::value;
            }();
        };

        template<class A, class B>
        struct contains_width<all_of_t<A, B>>
        {
            static constexpr bool value = [] {
                if constexpr (contains_width<std::remove_cvref_t<A>>::value)
                    return true;
                else
                    return contains_width<std::remove_cvref_t<B>>::value;
            }();
        };

        template<class C>
        inline constexpr bool contains_width_v = contains_width<std::remove_cvref_t<C>>::value;
    }    // namespace detail

    namespace detail
    {
        template<class T>
        constexpr bool mixed_tolerance_ok(const T& abs, const T& rel) noexcept
        { return tag::abs_tolerance::check(abs) && tag::rel_tolerance::check(rel) && (abs > T(0) || rel > T(0)); }

        // abs + rel_part (rel_part = rel * s >= 0), saturated at the largest finite value instead of overflowing to inf.
        // An overflowing sum accepted an infinite width: width_tol{1.7e308, nxx::rel_tolerance{0.5}} on [-1.7e308, 1.7e308] computed
        // 2.55e308 and 3.4e308 as inf, and inf <= inf (brent then stopped before its first step, since tol1 = inf / 2).
        // The test uses halves, which cannot overflow and round like the full sum (scaling by 2 is exact), so it
        // saturates exactly when the rounded sum would overflow, and every other sum is computed as before. An abs that
        // is inf in U (a double tolerance above FLT_MAX on a float problem) saturates too; NaN stays NaN.
        template<real U>
        constexpr U saturating_sum(const U& abs, const U& rel_part) noexcept
        {
            const U top = (std::numeric_limits<U>::max)();
            return abs / U(2) + rel_part / U(2) > top / U(2) ? top : abs + rel_part;
        }
    }    // namespace detail

    // |x_k - x_{k-1}| <= abs + rel * |x_k|. Open methods and systems only.
    //
    // Role-typed (DESIGN §6.2): one number is absolute, x_tol{1e-10}; the relative part is always named, mixed
    // x_tol{1e-10, nxx::rel_tolerance{1e-8}} and purely relative x_tol{0.0, nxx::rel_tolerance{1e-8}}, so the roles
    // cannot be swapped. A validated tolerance<T> takes a relative part too, x_tol{*tol, *rel}, also at run time. Two
    // bare numbers, and a part alone, are deleted with reasons. At run time: make(abs), make(abs, rel_tolerance),
    // make(abs_tolerance, rel_tolerance) and make(tolerance, rel_tolerance). The mixed literal constructor is consteval,
    // and cl 19.51 lacks P2564, so generic code that forwards the parts calls make() (§6.2 FLAG).
    template<real T>
    class x_tol : public criterion_base
    {
        T abs_;
        T rel_;

        constexpr x_tol(detail::trust_me, T a, T r) noexcept : abs_(a), rel_(r) {}

    public:
        static constexpr view_kind applies_to = view_kind::point | view_kind::system;

        constexpr x_tol(tolerance<T> a) noexcept : abs_(a.value()), rel_(T(0)) {}

        // The parts are checked by their own literals (abs finite and >= 0, 0 <= rel < 1); only the joint invariant is
        // left: x_tol{0.0, nxx::rel_tolerance{1e-8}} (purely relative) is legal, x_tol{0.0, nxx::rel_tolerance{0.0}} is
        // not. R is a template so that a bare number never converts into the relative part.
        template<class R>
            requires detail::is_rel_v<R, T>
        consteval x_tol(abs_tolerance<T> a, R r) : abs_(a.value()),
                                                   rel_(r.value())
        {
            if (!(abs_ > T(0) || rel_ > T(0)))
                detail::literal_violates_invariant("x_tol needs abs >= 0, 0 <= rel < 1, and abs > 0 or rel > 0");
        }

        // A validated absolute tolerance with a relative part (DESIGN §6.2, §12.24): tolerance<T> is finite and > 0, so
        // the joint invariant holds and nothing is left to check. constexpr, not consteval, so it takes run-time values
        // and forwards on cl. A is a template so that a bare number never converts into tolerance<T> here, which would
        // make {1e-10, nxx::rel_tolerance{1e-8}} ambiguous with the mixed literal.
        template<class A, class R>
            requires(std::same_as<std::remove_cvref_t<A>, tolerance<T>> && detail::is_rel_v<R, T>)
        constexpr x_tol(A a, R r) noexcept : abs_(a.value()),
                                             rel_(r.value())
        {}

        template<class A, class B>
            requires(detail::is_bare_number_v<A> && detail::is_bare_number_v<B>)
        x_tol(A, B) NXX_DELETE("say which number is relative: x_tol{1e-10, nxx::rel_tolerance{1e-8}}; purely relative: "
                               "x_tol{0.0, nxx::rel_tolerance{1e-8}}; one number is absolute: x_tol{1e-10}");

        template<class R>
            requires detail::is_tolerance_part_v<std::remove_cvref_t<R>>
        x_tol(R) NXX_DELETE("a part alone is not a criterion: x_tol{a} is absolute (x_tol<T>::make(a) for a number a at run time); "
                            "x_tol{0.0, nxx::rel_tolerance{r}} is purely relative (make(0.0, *rel) at run time)");

        // An absolute tolerance: finite and > 0, as tolerance<T>::make checks.
        static constexpr auto make(T a) noexcept -> std::expected<x_tol, errc>
        {
            if (!tag::positive_tolerance::check(a)) return std::unexpected(errc::invalid_input);
            return x_tol { detail::trust_me {}, a, T(0) };
        }

        // The run-time mirror of the literal x_tol{a, nxx::rel_tolerance{r}}: the absolute part is checked in-band
        // (finite and >= 0), with the joint invariant abs > 0 || rel > 0.
        template<class R>
            requires detail::is_rel_v<R, T>
        static constexpr auto make(T a, R r) noexcept -> std::expected<x_tol, errc>
        {
            if (!detail::mixed_tolerance_ok(a, r.value())) return std::unexpected(errc::invalid_input);
            return x_tol { detail::trust_me {}, a, r.value() };
        }

        // The same, from validated role types (a configuration holds abs_tolerance and rel_tolerance).
        template<class R>
            requires detail::is_rel_v<R, T>
        static constexpr auto make(abs_tolerance<T> a, R r) noexcept -> std::expected<x_tol, errc>
        { return make(a.value(), r); }

        // From a validated absolute tolerance: it cannot fail (tolerance<T> is finite and > 0). It returns std::expected
        // like the two-argument make forms; there is no one-argument make for a validated value (DESIGN §6.2).
        template<class A, class R>
            requires(std::same_as<std::remove_cvref_t<A>, tolerance<T>> && detail::is_rel_v<R, T>)
        static constexpr auto make(A a, R r) noexcept -> std::expected<x_tol, errc>
        { return x_tol { a, r }; }

        static void make(T, T) NXX_DELETE("say which number is relative: make(a, *rel) with rel = rel_tolerance<T>::make(r); "
                                          "make(a) for an absolute tolerance");

        constexpr T abs() const noexcept { return abs_; }
        constexpr T rel() const noexcept { return rel_; }
        // Computed in the scalar type of the problem, which may differ from T (a double tolerance on a float or
        // multiprecision problem).
        template<real U>
        constexpr U threshold(const U& x) const noexcept
        { return detail::saturating_sum<U>(U(abs_), U(rel_) * math::abs(x)); }

        template<class V>
        constexpr verdict operator()(const V& prev, const V& next, counters) const
        { return next.distance(prev) <= threshold(next.x()) ? verdict::converged : verdict::proceed; }
    };

    template<real T>
    x_tol(T) -> x_tol<T>;
    template<class A, real T>
        requires detail::is_bare_number_v<A>
    x_tol(A, rel_tolerance<T>) -> x_tol<T>;
    template<real T>
    x_tol(abs_tolerance<T>, rel_tolerance<T>) -> x_tol<T>;
    template<real T>
    x_tol(tolerance<T>, rel_tolerance<T>) -> x_tol<T>;
    // So that the deleted constructors report their reasons, not CTAD.
    template<class A, class B>
        requires(detail::is_bare_number_v<A> && detail::is_bare_number_v<B>)
    x_tol(A, B) -> x_tol<detail::bare_scalar_t<A>>;
    template<class R>
        requires detail::is_tolerance_part_v<R>
    x_tol(R) -> x_tol<typename R::value_type>;

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
    // tolerance. Bracketing methods only. Role-typed like x_tol (DESIGN §6.2): width_tol{1e-10} is absolute,
    // width_tol{1e-10, nxx::rel_tolerance{1e-8}} mixed, width_tol{0.0, nxx::rel_tolerance{1e-8}} purely relative.
    template<real T>
    class width_tol : public criterion_base
    {
        T abs_;
        T rel_;

        constexpr width_tol(detail::trust_me, T a, T r) noexcept : abs_(a), rel_(r) {}

    public:
        static constexpr view_kind applies_to = view_kind::enclosure;

        constexpr width_tol(tolerance<T> a) noexcept : abs_(a.value()), rel_(T(0)) {}

        // The parts are already checked by their literals; only the joint invariant is left.
        template<class R>
            requires detail::is_rel_v<R, T>
        consteval width_tol(abs_tolerance<T> a, R r) : abs_(a.value()),
                                                       rel_(r.value())
        {
            if (!(abs_ > T(0) || rel_ > T(0)))
                detail::literal_violates_invariant("width_tol needs abs >= 0, 0 <= rel < 1, and abs > 0 or rel > 0");
        }

        // A validated absolute tolerance with a relative part (DESIGN §6.2, §12.24): tolerance<T> is finite and > 0, so
        // the joint invariant holds and nothing is left to check. constexpr, not consteval, so it takes run-time values
        // and forwards on cl. A is a template so that a bare number never converts into tolerance<T> here, which would
        // make {1e-10, nxx::rel_tolerance{1e-8}} ambiguous with the mixed literal.
        template<class A, class R>
            requires(std::same_as<std::remove_cvref_t<A>, tolerance<T>> && detail::is_rel_v<R, T>)
        constexpr width_tol(A a, R r) noexcept : abs_(a.value()),
                                                 rel_(r.value())
        {}

        template<class A, class B>
            requires(detail::is_bare_number_v<A> && detail::is_bare_number_v<B>)
        width_tol(A, B) NXX_DELETE("say which number is relative: width_tol{1e-10, nxx::rel_tolerance{1e-8}}; purely relative: "
                                   "width_tol{0.0, nxx::rel_tolerance{1e-8}}; one number is absolute: width_tol{1e-10}");

        template<class R>
            requires detail::is_tolerance_part_v<std::remove_cvref_t<R>>
        width_tol(R) NXX_DELETE("a part alone is not a criterion: width_tol{a} is absolute "
                                "(width_tol<T>::make(a) for a number a at run time); "
                                "width_tol{0.0, nxx::rel_tolerance{r}} is purely relative (make(0.0, *rel) at run time)");

        // An absolute tolerance: finite and > 0, as tolerance<T>::make checks.
        static constexpr auto make(T a) noexcept -> std::expected<width_tol, errc>
        {
            if (!tag::positive_tolerance::check(a)) return std::unexpected(errc::invalid_input);
            return width_tol { detail::trust_me {}, a, T(0) };
        }

        // The run-time mirror of the literal width_tol{a, nxx::rel_tolerance{r}}: the absolute part is checked in-band
        // (finite and >= 0), with the joint invariant abs > 0 || rel > 0.
        template<class R>
            requires detail::is_rel_v<R, T>
        static constexpr auto make(T a, R r) noexcept -> std::expected<width_tol, errc>
        {
            if (!detail::mixed_tolerance_ok(a, r.value())) return std::unexpected(errc::invalid_input);
            return width_tol { detail::trust_me {}, a, r.value() };
        }

        // The same, from validated role types (a configuration holds abs_tolerance and rel_tolerance).
        template<class R>
            requires detail::is_rel_v<R, T>
        static constexpr auto make(abs_tolerance<T> a, R r) noexcept -> std::expected<width_tol, errc>
        { return make(a.value(), r); }

        // From a validated absolute tolerance: it cannot fail (tolerance<T> is finite and > 0). It returns std::expected
        // like the two-argument make forms; there is no one-argument make for a validated value (DESIGN §6.2).
        template<class A, class R>
            requires(std::same_as<std::remove_cvref_t<A>, tolerance<T>> && detail::is_rel_v<R, T>)
        static constexpr auto make(A a, R r) noexcept -> std::expected<width_tol, errc>
        { return width_tol { a, r }; }

        static void make(T, T) NXX_DELETE("say which number is relative: make(a, *rel) with rel = rel_tolerance<T>::make(r); "
                                          "make(a) for an absolute tolerance");

        constexpr T abs() const noexcept { return abs_; }
        constexpr T rel() const noexcept { return rel_; }

        // The threshold at a point x, and for an enclosure (min(|lo|, |hi|)), in the scalar type of the problem.
        template<real U>
        constexpr U threshold(const U& x) const noexcept
        { return detail::saturating_sum<U>(U(abs_), U(rel_) * math::abs(x)); }
        template<real U>
        constexpr U threshold(const U& lo, const U& hi) const noexcept
        { return detail::saturating_sum<U>(U(abs_), U(rel_) * (std::min)(math::abs(lo), math::abs(hi))); }

        template<class V>
        constexpr verdict operator()(const V&, const V& next, counters) const
        {
            const auto e = next.enclosure();
            return e.hi() - e.lo() <= threshold(e.lo(), e.hi()) ? verdict::converged : verdict::proceed;
        }
    };

    template<real T>
    width_tol(T) -> width_tol<T>;
    template<class A, real T>
        requires detail::is_bare_number_v<A>
    width_tol(A, rel_tolerance<T>) -> width_tol<T>;
    template<real T>
    width_tol(abs_tolerance<T>, rel_tolerance<T>) -> width_tol<T>;
    template<real T>
    width_tol(tolerance<T>, rel_tolerance<T>) -> width_tol<T>;
    // So that the deleted constructors report their reasons, not CTAD.
    template<class A, class B>
        requires(detail::is_bare_number_v<A> && detail::is_bare_number_v<B>)
    width_tol(A, B) -> width_tol<detail::bare_scalar_t<A>>;
    template<class R>
        requires detail::is_tolerance_part_v<R>
    width_tol(R) -> width_tol<typename R::value_type>;

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

    // Converged once n iterations have run. A guard: solvers accept it only under && with a convergence test, because
    // alone, or under ||, it would report success after n iterations without testing accuracy.
    class min_iterations : public criterion_base
    {
        std::uint32_t n_;

    public:
        static constexpr view_kind applies_to = all_views;
        static constexpr bool      guard_only = true;    // solvers take it only as `test && min_iterations{n}`

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

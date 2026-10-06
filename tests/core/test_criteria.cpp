// Stop criteria (DESIGN §6.8): the || and && verdict algebra, applies_to and criterion_for_v for every criterion and
// view kind, the thresholds of each criterion on small mock views, the defaults' achievability for float, double and
// long double (>= 4 eps, DESIGN §3.5), and the criteria driving nxx::iterate through a small mock solver.
#include <numerixx/core.hpp>

#include <doctest/doctest.h>

#include <array>
#include <cstdint>
#include <expected>
#include <limits>
#include <optional>
#include <type_traits>
#include <utility>

namespace
{
    using nxx::verdict;
    using nxx::view_kind;

    // What an open method's criteria see: x, f(x), and the distance to the previous iterate.
    template<class T>
    struct mock_point_view
    {
        T x_;
        T fx_;

        constexpr T x() const noexcept { return x_; }
        constexpr T fx() const noexcept { return fx_; }
        constexpr T distance(const mock_point_view& prev) const noexcept { return nxx::math::abs(x_ - prev.x_); }
        constexpr T residual() const noexcept { return nxx::math::abs(fx_); }
        constexpr T scale() const noexcept { return nxx::math::abs(x_); }
    };

    template<class T>
    struct mock_bounds
    {
        T lo_;
        T hi_;

        constexpr T lo() const noexcept { return lo_; }
        constexpr T hi() const noexcept { return hi_; }
    };

    // What a bracketing method's criteria see: an enclosure, and no distance().
    template<class T>
    struct mock_enclosure_view
    {
        T x_;
        T fx_;
        T lo_;
        T hi_;

        constexpr T              x() const noexcept { return x_; }
        constexpr T              fx() const noexcept { return fx_; }
        constexpr T              residual() const noexcept { return nxx::math::abs(fx_); }
        constexpr T              scale() const noexcept { return nxx::math::abs(x_); }
        constexpr mock_bounds<T> enclosure() const noexcept { return { lo_, hi_ }; }
    };

    using pv = mock_point_view<double>;
    using ev = mock_enclosure_view<double>;

    // A criterion that returns a fixed verdict, for the truth tables of || and &&.
    struct scripted : nxx::criterion_base
    {
        static constexpr view_kind applies_to = nxx::all_views;
        verdict                    v          = verdict::proceed;

        template<class V>
        constexpr verdict operator()(const V&, const V&, nxx::counters) const noexcept
        { return v; }
    };

    static_assert(nxx::criterion_for_v<scripted, view_kind::enclosure>);

    constexpr std::array<verdict, 4> all_verdicts { verdict::proceed, verdict::converged, verdict::stalled, verdict::exhausted };

    constexpr bool is_failure(verdict v) { return v == verdict::stalled || v == verdict::exhausted; }

    // a || b: the first verdict that is not proceed wins.
    constexpr verdict expected_or(verdict a, verdict b) { return a != verdict::proceed ? a : b; }

    // a && b: a failure of either wins (a's first); converged only when both converge.
    constexpr verdict expected_and(verdict a, verdict b)
    {
        if (is_failure(a)) return a;
        if (is_failure(b)) return b;
        return (a == verdict::converged && b == verdict::converged) ? verdict::converged : verdict::proceed;
    }

    // ---- A minimal solver over the protocol (DESIGN §6.6): x moves by a fixed step, each step costs `cost` -------------
    struct mock_estimate
    {
        double x;
        double fx;
    };

    constexpr bool better_than(const mock_estimate& a, const mock_estimate& b) noexcept
    { return nxx::math::abs(a.fx) < nxx::math::abs(b.fx); }

    struct mock_state
    {
        double        x;
        std::uint32_t nfev;
    };

    struct mock_problem
    {
    };

    // The failure estimate's order is a customisation point (DESIGN §6.6): mock_estimate has one (above); these have
    // none that nxx::better_than finds. A member function, or a function in a namespace that is not the type's, does
    // not count.
    struct unordered_estimate
    {
        double x;
        double fx;
    };

    struct member_ordered_estimate
    {
        double         x;
        double         fx;
        constexpr bool better_than(const member_ordered_estimate& b) const noexcept { return fx < b.fx; }
    };

    namespace estimate_home
    {
        struct foreign_estimate
        {
            double x;
            double fx;
        };
    }    // namespace estimate_home

    namespace elsewhere
    {
        [[maybe_unused]] constexpr bool better_than(const estimate_home::foreign_estimate& a,
                                                    const estimate_home::foreign_estimate& b) noexcept
        { return a.fx < b.fx; }
    }    // namespace elsewhere

    template<class Stop, class Est = mock_estimate>
    struct mock_solver
    {
        static constexpr nxx::algo id { nxx::algo::user_first };

        nxx::options<Stop> opt;
        std::uint32_t      init_cost = 1;
        std::uint32_t      cost      = 1;

        constexpr const nxx::options<Stop>& options() const noexcept { return opt; }

        constexpr auto init(const mock_problem&) const -> std::expected<mock_state, nxx::failure<Est>>
        { return mock_state { 10.0, init_cost }; }

        constexpr auto step(const mock_problem&, const mock_state& s) const -> std::expected<mock_state, nxx::fault<>>
        {
            return mock_state { s.x / 2.0, s.nfev + cost };    // never reaches 0: only the criterion or the budget stops it
        }

        constexpr pv                              view(const mock_state& s) const noexcept { return { s.x, s.x }; }
        constexpr Est                             estimate(const mock_state& s) const noexcept { return { s.x, s.x }; }
        constexpr Est                             best(const mock_state& s) const noexcept { return { s.x, s.x }; }
        constexpr std::optional<nxx::stop_reason> intrinsic(const mock_state&) const noexcept { return std::nullopt; }
    };

    // nxx::better_than is invocable only on an estimate type with an ADL better_than, and iterative_solver_for requires
    // one for the failure estimate type, so a solver whose estimate has no order is not a solver (DESIGN §6.6).
    template<class Est>
    constexpr bool ordered_v = std::is_invocable_v<decltype(nxx::better_than), const Est&, const Est&>;
    static_assert(ordered_v<mock_estimate> && nxx::detail::has_better_than_v<mock_estimate>);
    static_assert(!ordered_v<unordered_estimate> && !nxx::detail::has_better_than_v<unordered_estimate>);
    static_assert(!ordered_v<member_ordered_estimate> && !nxx::detail::has_better_than_v<member_ordered_estimate>);
    static_assert(!ordered_v<estimate_home::foreign_estimate> && !nxx::detail::has_better_than_v<estimate_home::foreign_estimate>);
    static_assert(!ordered_v<double> && !ordered_v<int>);
    static_assert(nxx::iterative_solver_for<mock_solver<nxx::never>, mock_problem>);
    static_assert(!nxx::iterative_solver_for<mock_solver<nxx::never, unordered_estimate>, mock_problem>);
    static_assert(!nxx::iterative_solver_for<mock_solver<nxx::never, member_ordered_estimate>, mock_problem>);
    static_assert(!nxx::iterative_solver_for<mock_solver<nxx::never, estimate_home::foreign_estimate>, mock_problem>);
    // The CPO is noexcept exactly when the order and its conversion to bool are: an order whose result converts to bool
    // through a conversion that may throw gives a CPO that is not noexcept.
    struct boolish
    {
        bool      value;
        constexpr operator bool() const { return value; }    // NOLINT(google-explicit-constructor): not noexcept
    };
    struct boolish_estimate
    {
        double x;
        double fx;
    };
    [[maybe_unused]] constexpr boolish better_than(const boolish_estimate& a, const boolish_estimate& b) noexcept
    { return { a.fx < b.fx }; }
    static_assert(ordered_v<boolish_estimate>);
    static_assert(noexcept(nxx::better_than(std::declval<const mock_estimate&>(), std::declval<const mock_estimate&>())));
    static_assert(!noexcept(nxx::better_than(std::declval<const boolish_estimate&>(), std::declval<const boolish_estimate&>())));
    // The order is used through the CPO, in constant expressions too.
    static_assert(nxx::better_than(mock_estimate { 0.0, 1.0 }, mock_estimate { 0.0, -2.0 }));
    static_assert(!nxx::better_than(mock_estimate { 0.0, 1.0 }, mock_estimate { 0.0, 1.0 }));

    template<class Stop>
    constexpr mock_solver<Stop> make_mock(Stop stop, nxx::max_iterations budget, std::uint32_t init_cost, std::uint32_t cost)
    { return mock_solver<Stop> { nxx::options<Stop> { stop, budget }, init_cost, cost }; }

    template<class T>
    constexpr T four_eps = T(4) * std::numeric_limits<T>::epsilon();
}    // namespace

TEST_SUITE("core")
{
    TEST_CASE("criteria algebra: || takes the first verdict that is not proceed")
    {
        const pv p0 { 1.0, 1.0 };
        for (const verdict a : all_verdicts)
            for (const verdict b : all_verdicts) {
                const auto either = scripted { {}, a } || scripted { {}, b };
                CHECK(either(p0, p0, nxx::counters {}) == expected_or(a, b));
            }
        static_assert((scripted { {}, verdict::converged } || scripted { {}, verdict::stalled })(pv {}, pv {}, {}) == verdict::converged);
    }

    TEST_CASE("criteria algebra: && converges only when both converge, and a failure verdict of either wins")
    {
        const ev e0 { 1.0, 1.0, 0.0, 2.0 };
        for (const verdict a : all_verdicts)
            for (const verdict b : all_verdicts) {
                const auto both = scripted { {}, a } && scripted { {}, b };
                CHECK(both(e0, e0, nxx::counters {}) == expected_and(a, b));
            }
        static_assert((scripted { {}, verdict::converged } && scripted { {}, verdict::proceed })(pv {}, pv {}, {}) == verdict::proceed);
        static_assert((scripted { {}, verdict::converged } && scripted { {}, verdict::exhausted })(pv {}, pv {}, {}) == verdict::exhausted);
    }

    TEST_CASE("criteria algebra: nesting and the real criteria")
    {
        // min_iterations{3} && x_tol{1e-6}: a small step alone is not enough before the third iteration.
        const auto c = nxx::min_iterations { 3 } && nxx::x_tol { 1e-6 };
        const pv   a { 1.0, 0.1 };
        const pv   b { 1.0 + 1e-8, 0.1 };
        CHECK(c(a, b, nxx::counters { 2, 5 }) == verdict::proceed);
        CHECK(c(a, b, nxx::counters { 3, 5 }) == verdict::converged);
        CHECK(c(a, pv { 2.0, 0.1 }, nxx::counters { 3, 5 }) == verdict::proceed);

        // (never || max_evaluations{10}) || f_tol{1e-8}: the budget verdict comes first.
        const auto d = (nxx::never {} || nxx::max_evaluations { 10 }) || nxx::f_tol { 1e-8 };
        CHECK(d(a, pv { 1.0, 0.0 }, nxx::counters { 1, 10 }) == verdict::exhausted);
        CHECK(d(a, pv { 1.0, 0.0 }, nxx::counters { 1, 9 }) == verdict::converged);
        CHECK(d(a, pv { 1.0, 1.0 }, nxx::counters { 1, 9 }) == verdict::proceed);

        // floored_width{} && f_tol: both the enclosure and the residual.
        const auto w = nxx::floored_width {} && nxx::f_tol { 1e-6 };
        CHECK(w(ev {}, ev { 1.0, 1e-7, 1.0, 1.0 }, nxx::counters {}) == verdict::converged);
        CHECK(w(ev {}, ev { 1.0, 1e-5, 1.0, 1.0 }, nxx::counters {}) == verdict::proceed);
        CHECK(w(ev {}, ev { 1.0, 1e-7, 1.0, 2.0 }, nxx::counters {}) == verdict::proceed);
    }

    TEST_CASE("applies_to: each criterion's view kinds, and the intersection under || and &&")
    {
        static_assert(nxx::x_tol<double>::applies_to == (view_kind::point | view_kind::system));
        static_assert(nxx::step_tol<3, 5>::applies_to == (view_kind::point | view_kind::system));
        static_assert(nxx::width_tol<double>::applies_to == view_kind::enclosure);
        static_assert(nxx::floored_width::applies_to == view_kind::enclosure);
        static_assert(nxx::f_tol<double>::applies_to == nxx::all_views);
        static_assert(nxx::max_evaluations::applies_to == nxx::all_views);
        static_assert(nxx::min_iterations::applies_to == nxx::all_views);
        static_assert(nxx::never::applies_to == nxx::all_views);
        static_assert(std::to_underlying(nxx::all_views) == 7);
        static_assert(((view_kind::point | view_kind::enclosure) & view_kind::enclosure) == view_kind::enclosure);

        using x_or_budget = decltype(nxx::x_tol { 1e-6 } || nxx::max_evaluations { 10 });
        static_assert(x_or_budget::applies_to == (view_kind::point | view_kind::system));
        using width_and_f = decltype(nxx::floored_width {} && nxx::f_tol { 1e-6 });
        static_assert(width_and_f::applies_to == view_kind::enclosure);
        using nothing = decltype(nxx::width_tol { 1e-6 } || nxx::x_tol { 1e-6 });
        static_assert(std::to_underlying(nothing::applies_to) == 0);
        CHECK(true);
    }

    TEST_CASE("criterion_for_v for every criterion and view kind")
    {
        constexpr auto P = view_kind::point;
        constexpr auto E = view_kind::enclosure;
        constexpr auto S = view_kind::system;

        static_assert(nxx::criterion_for_v<nxx::x_tol<double>, P> && !nxx::criterion_for_v<nxx::x_tol<double>, E> &&
                      nxx::criterion_for_v<nxx::x_tol<double>, S>);
        static_assert(nxx::criterion_for_v<nxx::step_tol<3, 5>, P> && !nxx::criterion_for_v<nxx::step_tol<3, 5>, E> &&
                      nxx::criterion_for_v<nxx::step_tol<3, 5>, S>);
        static_assert(nxx::criterion_for_v<nxx::step_tol<7, 10>, P> && !nxx::criterion_for_v<nxx::step_tol<7, 10>, E>);
        static_assert(!nxx::criterion_for_v<nxx::width_tol<double>, P> && nxx::criterion_for_v<nxx::width_tol<double>, E> &&
                      !nxx::criterion_for_v<nxx::width_tol<double>, S>);
        static_assert(!nxx::criterion_for_v<nxx::floored_width, P> && nxx::criterion_for_v<nxx::floored_width, E> &&
                      !nxx::criterion_for_v<nxx::floored_width, S>);
        static_assert(nxx::criterion_for_v<nxx::f_tol<double>, P> && nxx::criterion_for_v<nxx::f_tol<double>, E> &&
                      nxx::criterion_for_v<nxx::f_tol<double>, S>);
        static_assert(nxx::criterion_for_v<nxx::max_evaluations, P> && nxx::criterion_for_v<nxx::max_evaluations, E> &&
                      nxx::criterion_for_v<nxx::max_evaluations, S>);
        static_assert(nxx::criterion_for_v<nxx::min_iterations, P> && nxx::criterion_for_v<nxx::min_iterations, E> &&
                      nxx::criterion_for_v<nxx::min_iterations, S>);
        static_assert(nxx::criterion_for_v<nxx::never, P> && nxx::criterion_for_v<nxx::never, E> && nxx::criterion_for_v<nxx::never, S>);

        // min_iterations is a guard: it applies everywhere, but a solver takes it only under && with a real test.
        using guarded     = decltype(nxx::x_tol { 1e-6 } && nxx::min_iterations { 3 });
        using guard_or    = decltype(nxx::x_tol { 1e-6 } || nxx::min_iterations { 3 });
        using two_guards  = decltype(nxx::min_iterations { 2 } && nxx::min_iterations { 3 });
        using guard_limit = decltype((nxx::x_tol { 1e-6 } && nxx::min_iterations { 3 }) || nxx::max_evaluations { 50 });
        static_assert(!nxx::stop_criterion_for_v<nxx::min_iterations, P> && !nxx::stop_criterion_for_v<nxx::min_iterations, E>);
        static_assert(nxx::stop_criterion_for_v<guarded, P>);
        static_assert(!nxx::stop_criterion_for_v<guard_or, P>);    // || can stop on the guard alone
        static_assert(!nxx::stop_criterion_for_v<two_guards, P>);
        static_assert(nxx::stop_criterion_for_v<guard_limit, P>);
        static_assert(nxx::stop_criterion_for_v<nxx::never, P> && nxx::stop_criterion_for_v<nxx::floored_width, E>);

        // cv-ref qualified criteria, combinations, and non-criteria.
        static_assert(nxx::criterion_for_v<const nxx::floored_width&, E>);
        using x_or_budget = decltype(nxx::x_tol { 1e-6 } || nxx::max_evaluations { 10 });
        static_assert(nxx::criterion_for_v<x_or_budget, P> && !nxx::criterion_for_v<x_or_budget, E>);
        using never_or_budget = decltype(nxx::never {} || nxx::max_evaluations { 10 });
        static_assert(nxx::criterion_for_v<never_or_budget, E> && nxx::criterion_for_v<never_or_budget, P>);
        using nothing = decltype(nxx::width_tol { 1e-6 } && nxx::x_tol { 1e-6 });
        static_assert(!nxx::criterion_for_v<nothing, P> && !nxx::criterion_for_v<nothing, E> && !nxx::criterion_for_v<nothing, S>);
        static_assert(!nxx::criterion_for_v<double, P>);
        static_assert(!nxx::criterion_for_v<nxx::tolerance<double>, E>);
        static_assert(!nxx::criterion_for_v<nxx::max_iterations, P>);

        static_assert(nxx::is_criterion_v<nxx::never> && nxx::is_criterion_v<const nxx::x_tol<double>&>);
        static_assert(!nxx::is_criterion_v<double> && !nxx::is_criterion_v<nxx::tolerance<double>>);
        CHECK(true);
    }

    TEST_CASE("x_tol: |x_k - x_{k-1}| <= abs + rel |x_k|")
    {
        const nxx::x_tol absolute { 1e-6 };
        CHECK(absolute(pv { 1.0, 0.0 }, pv { 1.0 + 5e-7, 0.0 }, {}) == verdict::converged);
        CHECK(absolute(pv { 1.0, 0.0 }, pv { 1.0 + 2e-6, 0.0 }, {}) == verdict::proceed);

        const nxx::x_tol relative { 0.0, 1e-8 };
        CHECK(relative.threshold(1e6) == doctest::Approx(1e-2));
        CHECK(relative(pv { 1e6, 0.0 }, pv { 1e6 + 5e-3, 0.0 }, {}) == verdict::converged);
        CHECK(relative(pv { 1e6, 0.0 }, pv { 1e6 + 2e-2, 0.0 }, {}) == verdict::proceed);
        CHECK(relative(pv { 0.0, 0.0 }, pv { 1e-300, 0.0 }, {}) == verdict::proceed);    // purely relative: no floor at 0

        const nxx::x_tol mixed { 1e-6, 1e-3 };
        CHECK(mixed.threshold(-10.0) == doctest::Approx(1e-6 + 1e-2));
        static_assert(nxx::x_tol { 1e-6 }(pv { 1.0, 0.0 }, pv { 1.0, 0.0 }, {}) == verdict::converged);
    }

    TEST_CASE("step_tol: |dx| <= 2^-ceil(p Num / Den) max(|x|, 1)")
    {
        static_assert(nxx::step_tol<3, 5>::threshold(0.0) == 0x1p-32);    // ceil(53 * 3 / 5) = 32; floor at scale 1
        static_assert(nxx::step_tol<3, 5>::threshold(8.0) == 0x1p-29);
        static_assert(nxx::step_tol<3, 5>::threshold(-8.0) == 0x1p-29);
        static_assert(nxx::step_tol<7, 10>::threshold(1.0) == 0x1p-38);      // ceil(53 * 7 / 10) = 38
        static_assert(nxx::step_tol<3, 5>::threshold(1.0f) == 0x1p-15f);     // ceil(24 * 3 / 5) = 15
        static_assert(nxx::step_tol<7, 10>::threshold(1.0f) == 0x1p-17f);    // ceil(24 * 7 / 10) = 17

        const nxx::step_tol<3, 5> s;
        CHECK(s(pv { 0.0, 1.0 }, pv { 0x1p-33, 1.0 }, {}) == verdict::converged);    // a root at 0 terminates
        CHECK(s(pv { 0.0, 1.0 }, pv { 0x1p-31, 1.0 }, {}) == verdict::proceed);
        CHECK(s(pv { 1e6, 1.0 }, pv { 1e6 + 1e-4, 1.0 }, {}) == verdict::converged);    // relative above 1
        CHECK(s(pv { 1e6, 1.0 }, pv { 1e6 + 1e-3, 1.0 }, {}) == verdict::proceed);
    }

    TEST_CASE("width_tol: hi - lo <= abs + rel min(|lo|, |hi|)")
    {
        const nxx::width_tol absolute { 1e-3 };
        CHECK(absolute(ev {}, ev { 1.0, 0.0, 1.0, 1.0005 }, {}) == verdict::converged);
        CHECK(absolute(ev {}, ev { 1.0, 0.0, 1.0, 1.002 }, {}) == verdict::proceed);

        const nxx::width_tol relative { 0.0, 1e-3 };
        CHECK(relative.threshold(1000.0, 1000.5) == doctest::Approx(1.0));
        CHECK(relative.threshold(-1000.5, -1000.0) == doctest::Approx(1.0));    // min(|lo|, |hi|) = 1000
        CHECK(relative(ev {}, ev { 1000.0, 0.0, 1000.0, 1000.5 }, {}) == verdict::converged);
        CHECK(relative(ev {}, ev { -1000.0, 0.0, -1000.5, -1000.0 }, {}) == verdict::converged);
        CHECK(relative(ev {}, ev { 1000.0, 0.0, 1000.0, 1002.0 }, {}) == verdict::proceed);
        // A purely relative width test cannot be met by an enclosure of 0 (why floored_width has a floor).
        CHECK(relative(ev {}, ev { 0.0, 0.0, -1e-300, 1e-300 }, {}) == verdict::proceed);
    }

    TEST_CASE("width_tol and x_tol: a threshold that overflows saturates, so an infinite width or distance never passes")
    {
        // abs + rel * s overflowed to inf when abs is near max, and so did the width of [-max, max]: inf <= inf said
        // converged although 2 max > 1.5 max. Saturated at max, the threshold is still a lower bound of the exact one,
        // and it stays a constant expression (GCC rejects an overflow there).
        constexpr double         big = (std::numeric_limits<double>::max)();
        constexpr nxx::width_tol wide { big, 0.5 };
        static_assert(wide.threshold(-big, big) == big);
        static_assert(wide.threshold(big) == big);
        CHECK(wide(ev {}, ev { big, 0.0, -big, big }, {}) == verdict::proceed);            // 2 max > 1.5 max
        CHECK(wide(ev {}, ev { big, 0.0, 0.25 * big, big }, {}) == verdict::converged);    // a finite width still passes

        constexpr nxx::x_tol far { 1e308, 0.5 };
        static_assert(far.threshold(big) == big);
        CHECK(far(pv { -big, 0.0 }, pv { big, 0.0 }, {}) == verdict::proceed);    // distance 2 max > 1e308 + max / 2
    }

    TEST_CASE("floored_width: w <= max(2^(1 - bits), 4 eps) max(1, min(|lo|, |hi|))")
    {
        constexpr double         eps = std::numeric_limits<double>::epsilon();
        const nxx::floored_width fw;
        CHECK(fw.factor<double>() == 4.0 * eps);    // bits = digits: 2^(1 - 53) = eps, floored at 4 eps
        CHECK(fw(ev {}, ev { 1.0, 0.0, 1.0, 1.0 + 4.0 * eps }, {}) == verdict::converged);
        CHECK(fw(ev {}, ev { 1.0, 0.0, 1.0, 1.0 + 8.0 * eps }, {}) == verdict::proceed);
        CHECK(fw(ev {}, ev { 0.0, 0.0, -1e-17, 1e-17 }, {}) == verdict::converged);       // absolute floor: a root at 0 terminates
        CHECK(fw(ev {}, ev { 1e6, 0.0, 1e6, 1e6 + 5e-10 }, {}) == verdict::converged);    // 4 ulps of 1e6 <= 4 eps * 1e6
        CHECK(fw(ev {}, ev { 1e6, 0.0, 1e6, 1e6 + 1e-8 }, {}) == verdict::proceed);

        constexpr nxx::floored_width coarse { 10 };
        static_assert(coarse.bits() == 10);
        static_assert(coarse.factor<double>() == 0x1p-9);
        static_assert(coarse.threshold(4.0) == 0x1p-7);
        static_assert(coarse.threshold(-4.0, -8.0) == 0x1p-7);
        static_assert(nxx::floored_width { 100 }.factor<double>() == 4.0 * eps);    // more bits than T has: floored
        CHECK(coarse(ev {}, ev { 0.5, 0.0, 0.5, 0.5 + 0x1p-10 }, {}) == verdict::converged);
    }

    TEST_CASE("f_tol, max_evaluations, min_iterations, never")
    {
        const nxx::f_tol ft { 1e-8 };
        CHECK(ft(pv {}, pv { 1.0, -5e-9 }, {}) == verdict::converged);
        CHECK(ft(pv {}, pv { 1.0, 2e-8 }, {}) == verdict::proceed);
        CHECK(ft(ev {}, ev { 1.0, 5e-9, 0.0, 2.0 }, {}) == verdict::converged);    // every view kind

        const nxx::max_evaluations me { 10 };
        CHECK(me.value() == 10u);
        CHECK(me(pv {}, pv {}, nxx::counters { 3, 9 }) == verdict::proceed);
        CHECK(me(pv {}, pv {}, nxx::counters { 3, 10 }) == verdict::exhausted);
        const auto budget = nxx::evaluation_budget::make(4u);    // a run-time budget goes through make()
        CHECK(budget.has_value());
        if (budget) {
            const nxx::max_evaluations rt_me { *budget };
            CHECK(rt_me(ev {}, ev {}, nxx::counters { 1, 4 }) == verdict::exhausted);
        }
        else {
            FAIL_CHECK("evaluation_budget::make(4) failed");
        }

        const nxx::min_iterations mi { 5 };
        CHECK(mi.value() == 5u);
        CHECK(mi(pv {}, pv {}, nxx::counters { 4, 100 }) == verdict::proceed);
        CHECK(mi(pv {}, pv {}, nxx::counters { 5, 100 }) == verdict::converged);

        const nxx::never nv;
        CHECK(nv(pv {}, pv {}, nxx::counters { 1000, 1000 }) == verdict::proceed);
        CHECK(nv(ev {}, ev {}, nxx::counters {}) == verdict::proceed);
    }

    TEST_CASE_TEMPLATE("the defaults are achievable: step_tol and floored_width thresholds >= 4 eps", T, float, double, long double)
    {
        static_assert(nxx::floored_width {}.factor<T>() >= four_eps<T>);
        static_assert(nxx::floored_width {}.threshold(T(1)) >= four_eps<T>);
        static_assert(nxx::floored_width {}.threshold(T(0)) >= four_eps<T>);    // the absolute floor at 0
        static_assert(nxx::step_tol<3, 5>::threshold(T(1)) >= four_eps<T>);     // Newton's default
        static_assert(nxx::step_tol<7, 10>::threshold(T(1)) >= four_eps<T>);    // secant's default
        static_assert(nxx::step_tol<3, 5>::threshold(T(0)) >= four_eps<T>);
        static_assert(nxx::step_tol<7, 10>::threshold(T(0)) >= four_eps<T>);

        // Relative above 1: still >= 4 eps |x| (four ulps or more of x).
        const T big = T(1e6);
        CHECK(nxx::floored_width {}.threshold(big) >= four_eps<T> * big);
        CHECK(nxx::step_tol<3, 5>::threshold(big) >= four_eps<T> * big);
        CHECK(nxx::step_tol<7, 10>::threshold(big) >= four_eps<T> * big);

        // An enclosure of exactly the default width converges; twice that does not.
        using V        = mock_enclosure_view<T>;
        const T w      = nxx::floored_width {}.factor<T>();
        const T one    = T(1);
        const T one_w  = one + w;
        const T one_2w = one + T(2) * w;
        const V prev   = { one, T(0), one, T(2) };
        const V tight  = { one, T(0), one, one_w };
        const V loose  = { one, T(0), one, one_2w };
        CHECK(nxx::floored_width {}(prev, tight, nxx::counters {}) == verdict::converged);
        CHECK(nxx::floored_width {}(prev, loose, nxx::counters {}) == verdict::proceed);

        // A step of exactly the default threshold converges.
        using P      = mock_point_view<T>;
        const T st   = nxx::step_tol<7, 10>::threshold(one);
        const P from = { one, T(1) };
        const P to   = { one + st, T(1) };
        CHECK(nxx::step_tol<7, 10> {}(from, to, nxx::counters {}) == verdict::converged);
    }

    TEST_CASE("criteria drive nxx::iterate: never || max_evaluations, the budget, min_iterations")
    {
        // init costs 1 evaluation and every step 1: after k steps nfev = 1 + k, so 10 evaluations are spent at k = 9.
        const auto spent =
            nxx::iterate(make_mock(nxx::never {} || nxx::max_evaluations { 10 }, nxx::max_iterations { 100 }, 1, 1), mock_problem {});
        CHECK_FALSE(spent.has_value());
        if (!spent) {
            CHECK(spent.error().code == nxx::errc::evaluations_exhausted);
            CHECK(spent.error().used == nxx::counters { 9, 10 });
            CHECK(spent.error().by == nxx::algo::user_first);
            CHECK(spent.error().best.has_value());
        }

        // Two evaluations per step: 10 are spent after 5 steps (nfev = 0 + 2k).
        const auto spent2 =
            nxx::iterate(make_mock(nxx::never {} || nxx::max_evaluations { 10 }, nxx::max_iterations { 100 }, 0, 2), mock_problem {});
        CHECK((!spent2 && spent2.error().code == nxx::errc::evaluations_exhausted && spent2.error().used == nxx::counters { 5, 10 }));

        // The iteration budget: budget_exhausted after exactly n steps, never success.
        const auto starved = nxx::iterate(make_mock(nxx::never {}, nxx::max_iterations { 5 }, 1, 1), mock_problem {});
        CHECK((!starved && starved.error().code == nxx::errc::budget_exhausted && starved.error().used == nxx::counters { 5, 6 }));
        if (!starved && starved.error().best) {
            CHECK(starved.error().best->x == 10.0 / 32.0);    // the best iterate (smallest |f|) is the last one here
        }

        // min_iterations{3}: converged with stop_reason::criterion at the third iteration.
        const auto three = nxx::iterate(make_mock(nxx::min_iterations { 3 }, nxx::max_iterations { 100 }, 1, 1), mock_problem {});
        CHECK(three.has_value());
        if (three) {
            CHECK(three->how == nxx::stop_reason::criterion);
            CHECK(three->used == nxx::counters { 3, 4 });
            CHECK(three->x == 10.0 / 8.0);
            CHECK(three->by == nxx::algo::user_first);
        }

        // In a constant expression too.
        static_assert(
            nxx::iterate(make_mock(nxx::never {} || nxx::max_evaluations { 10 }, nxx::max_iterations { 100 }, 1, 1), mock_problem {})
                .error()
                .code == nxx::errc::evaluations_exhausted);
    }
}

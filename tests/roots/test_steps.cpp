// Manual stepping (DESIGN §6.9, §9.3): nxx::steps_view is an input range and a view over the solver protocol. Element 0
// is init(p); the range ends after the first error or the first intrinsic stop; it composes with std::views; and it
// yields the same iterates as the driver, up to the driver's stopping point.
#include <numerixx/roots.hpp>

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <functional>
#include <limits>
#include <optional>
#include <random>
#include <ranges>
#include <type_traits>
#include <utility>
#include <vector>

// The iterates are compared bit for bit with the driver's: no fused multiply-adds in this file's functions (Clang).
#if defined(__clang__)
#    pragma clang fp contract(off)
#endif

namespace
{
    namespace r = nxx::roots;
    using est_t = r::root_estimate<double>;

    constexpr auto sq2  = [](double x) { return x * x - 2.0; };    // root sqrt(2)
    constexpr auto dsq2 = [](double x) { return 2.0 * x; };

    enum class step_error { domain };

    // x^2 - 2 for x >= 0, a user error for x < 0.
    constexpr auto g_neg = [](double x) -> std::expected<double, step_error> {
        if (x < 0.0) return std::unexpected(step_error::domain);
        return x * x - 2.0;
    };

    // A newton whose step returns an input code itself, as a user-written step might (DESIGN §6.7).
    using newton_t = decltype(r::newton {}.with_derivative(dsq2));
    struct raw_input_step : newton_t
    {
        template<class P, class S>
        constexpr auto step(const P&, const S&) const -> std::expected<S, nxx::fault<step_error>>
        { return std::unexpected(nxx::fault<step_error> { nxx::errc::invalid_input, 1, step_error::domain }); }
    };

    // ---- Types ---------------------------------------------------------------------------------------------------------
    template<class S, class F, class In>
    using problem_of = typename decltype(std::declval<const S&>().prepare(std::declval<const std::reference_wrapper<const F>&>(),
                                                                          std::declval<const In&>()))::value_type;

    using sq2_t   = std::remove_const_t<decltype(sq2)>;
    using g_neg_t = std::remove_const_t<decltype(g_neg)>;

    using brent_view     = nxx::steps_view<r::brent<>, problem_of<r::brent<>, sq2_t, std::pair<double, double>>>;
    using bisection_view = nxx::steps_view<r::bisection<>, problem_of<r::bisection<>, sq2_t, nxx::bracket<double>>>;
    using secant_view    = nxx::steps_view<r::secant<>, problem_of<r::secant<>, sq2_t, double>>;
    using newton_view    = nxx::steps_view<newton_t, problem_of<newton_t, sq2_t, double>>;
    using newton_g_view  = nxx::steps_view<newton_t, problem_of<newton_t, g_neg_t, double>>;

    static_assert(std::ranges::input_range<brent_view>);
    static_assert(std::ranges::view<brent_view>);
    static_assert(!std::ranges::forward_range<brent_view>);    // the iterator caches its element: single pass
    static_assert(std::ranges::input_range<bisection_view> && std::ranges::view<bisection_view>);
    static_assert(std::ranges::input_range<secant_view> && std::ranges::view<secant_view>);
    static_assert(std::ranges::input_range<newton_view> && std::ranges::view<newton_view>);
    static_assert(std::ranges::input_range<newton_g_view> && std::ranges::view<newton_g_view>);

    // Elements: expected<state, fault<UE>>, UE the user's callback error (none for a plain f).
    static_assert(std::is_same_v<brent_view::element, std::expected<r::brent_state<double>, nxx::fault<nxx::none>>>);
    static_assert(std::is_same_v<std::ranges::range_value_t<brent_view>, brent_view::element>);
    static_assert(std::is_same_v<newton_g_view::element, std::expected<r::newton_state<double>, nxx::fault<step_error>>>);

    // Bounded and transformed views stay views.
    using taken_t = decltype(std::declval<brent_view>() | std::views::take(8));
    static_assert(std::ranges::input_range<taken_t> && std::ranges::view<taken_t>);

    // ---- Helpers -------------------------------------------------------------------------------------------------------
    bool same_bits(double a, double b) { return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b); }

    bool same_estimate(const est_t& a, const est_t& b)
    { return same_bits(a.x, b.x) && same_bits(a.fx, b.fx) && same_bits(a.uncertainty, b.uncertainty) && a.enclosure == b.enclosure; }

    bool same_state(const r::brent_state<double>& a, const r::brent_state<double>& b)
    {
        return same_bits(a.a, b.a) && same_bits(a.fa, b.fa) && same_bits(a.b, b.b) && same_bits(a.fb, b.fb) && same_bits(a.c, b.c) &&
               same_bits(a.fc, b.fc) && same_bits(a.d, b.d) && same_bits(a.e, b.e) && a.nfev == b.nfev;
    }

    // Runs a bracketing solver through the driver (recording the view of every iterate with an observer) and through
    // steps_view, and checks that the view yields the driver's iterates, that the driver stopped at the first element
    // where the solver's intrinsic test or its stop criterion fires, and that the element there is the driver's result.
    template<class S, class F>
    void check_against_driver(const S& solver, const F& fn, double lo, double hi)
    {
        std::vector<std::array<double, 4>> seen;
        const auto                         observed =
            solver.with_observer([&seen](const auto& v) { seen.push_back({ v.x(), v.fx(), v.enclosure().lo(), v.enclosure().hi() }); });
        const auto driver = observed(fn, std::pair { lo, hi });
        const auto p      = solver.prepare(std::cref(fn), std::pair { lo, hi });
        if (!driver || !p) {
            FAIL_CHECK("precondition: the driver succeeds and the problem is valid");
            return;
        }
        const std::uint32_t k_stop = driver->used.iterations;
        CHECK(seen.size() == k_stop);    // the observer sees every iterate the driver takes

        using state_t = std::remove_cvref_t<decltype(*solver.init(*p))>;
        std::optional<state_t> prev;
        std::uint32_t          k = 0;
        for (const auto& st : nxx::steps_view { solver, *p } | std::views::take(static_cast<std::ptrdiff_t>(k_stop) + 1)) {
            if (!st) {
                FAIL_CHECK("an error element on a well-posed problem");
                return;
            }
            if (k >= 1 && k <= seen.size()) {    // the same iterate the driver observed
                const auto  v = solver.view(*st);
                const auto& o = seen[k - 1];
                CHECK(same_bits(v.x(), o[0]));
                CHECK(same_bits(v.fx(), o[1]));
                CHECK(same_bits(v.enclosure().lo(), o[2]));
                CHECK(same_bits(v.enclosure().hi(), o[3]));
            }
            const std::optional<nxx::stop_reason> how = solver.intrinsic(*st);
            const nxx::verdict                    vd =
                prev ? solver.options().stop(solver.view(*prev), solver.view(*st), nxx::counters { k, st->nfev }) : nxx::verdict::proceed;
            if (k < k_stop) {    // the driver did not stop earlier than it had to
                CHECK_FALSE(how.has_value());
                CHECK(vd == nxx::verdict::proceed);
            }
            else {    // ... and stopped here, with this element's estimate
                CHECK((how.has_value() || vd == nxx::verdict::converged));
                CHECK(driver->how == (how ? *how : nxx::stop_reason::criterion));
                CHECK(st->nfev == driver->used.evaluations);
                CHECK(same_estimate(solver.estimate(*st), *driver));
            }
            prev = *st;
            ++k;
        }
        CHECK(k == k_stop + 1);    // init plus one element per iteration
    }

    template<class S>
    void check_against_driver_on_cubics(const S& solver)
    {
        std::mt19937                           gen(20260928u);
        std::uniform_real_distribution<double> coef(0.5, 20.0);    // roots cbrt(c) in (0.79, 2.72), inside [0, 3]
        for (int i = 0; i < 100; ++i) {
            const double c = coef(gen);
            check_against_driver(solver, [c](double x) { return x * x * x - c; }, 0.0, 3.0);
        }
    }
}    // namespace

TEST_SUITE("roots")
{
    TEST_CASE("steps_view: is an input range and a view")
    {
        // The contracts are the static_asserts above; this also builds one of each at run time.
        const auto solver = r::brent {};
        const auto p      = solver.prepare(std::cref(sq2), std::pair { 1.0, 2.0 });
        if (!p) {
            FAIL_CHECK("precondition: [1, 2] brackets sqrt(2)");
            return;
        }
        const nxx::steps_view view { solver, *p };
        static_assert(std::is_same_v<std::remove_const_t<decltype(view)>, brent_view>);
        CHECK((view.begin() != view.end()));
        CHECK((*view.begin()).has_value());
    }

    TEST_CASE("steps_view: take(8) over brent starts at init and ends at the intrinsic stop")
    {
        const auto solver = r::brent {};
        const auto p      = solver.prepare(std::cref(sq2), nxx::bracket { 1.0, 2.0 });
        const auto driver = solver(sq2, nxx::bracket { 1.0, 2.0 });
        if (!p || !driver) {
            FAIL_CHECK("precondition: brent succeeds on x2 - 2 over [1, 2]");
            return;
        }
        const auto init = solver.init(*p);
        if (!init) {
            FAIL_CHECK("brent's init cannot fail on a prepared problem");
            return;
        }

        std::vector<r::brent_state<double>> states;
        std::size_t                         errors = 0;
        for (const auto& st : nxx::steps_view { solver, *p } | std::views::take(8)) {
            if (st)
                states.push_back(*st);
            else
                ++errors;
        }
        CHECK(errors == 0u);
        const std::size_t n_driver = static_cast<std::size_t>(driver->used.iterations) + 1;    // init + one per iteration
        CHECK(states.size() == std::min<std::size_t>(8, n_driver));
        if (states.empty()) return;

        CHECK(same_state(states.front(), *init));    // element 0 is init(p)
        CHECK(states.front().nfev == 2u);            // the two endpoint samples of prepare()
        for (std::size_t i = 0; i + 1 < states.size(); ++i) CHECK_FALSE(solver.intrinsic(states[i]).has_value());
        if (n_driver <= 8) {    // it converged within the window: the last element is the intrinsic stop
            CHECK(solver.intrinsic(states.back()) == std::optional { nxx::stop_reason::criterion });
            CHECK(same_estimate(solver.estimate(states.back()), *driver));
        }

        // Unbounded, it still ends by itself at the intrinsic stop.
        const auto n_all = std::ranges::distance(nxx::steps_view { solver, *p } | std::views::take(1000));
        CHECK(n_all == static_cast<std::ptrdiff_t>(n_driver));
        CHECK(n_all < 1000);

        // A view is a value: begin() starts again from init(p), so a second pass yields the same states.
        const nxx::steps_view view { solver, *p };
        std::size_t           i = 0;
        for (const auto& st : view | std::views::take(8)) {
            if (st && i < states.size()) { CHECK(same_state(*st, states[i])); }
            ++i;
        }
        CHECK(i == states.size());
    }

    TEST_CASE("steps_view: element 0 of bisection is init and each step halves the sign-changing bracket")
    {
        const auto solver = r::bisection {};
        const auto p      = solver.prepare(std::cref(sq2), nxx::bracket { 0.0, 2.0 });
        if (!p) {
            FAIL_CHECK("precondition: [0, 2] brackets sqrt(2)");
            return;
        }
        std::optional<r::bisection_state<double>> prev;
        std::size_t                               n = 0;
        for (const auto& st : nxx::steps_view { solver, *p } | std::views::take(20)) {
            if (!st) {
                FAIL_CHECK("no error on a well-posed problem");
                break;
            }
            if (!prev) {
                CHECK(st->b == p->in);    // init(p): the prepared bracket
                CHECK(st->nfev == 2u);
            }
            else {
                CHECK(st->b.width() == prev->b.width() / 2.0);    // exact: halving [0, 2] in binary
                CHECK(st->b.flo() < 0.0);                         // the sign change is kept
                CHECK(st->b.fhi() > 0.0);
                CHECK(st->nfev == prev->nfev + 1u);
            }
            prev = *st;
            ++n;
        }
        CHECK(n == 20u);    // bisection has no intrinsic stop this early: the view runs past nothing but take(20)
    }

    TEST_CASE("steps_view: yields the driver's iterates up to its stopping point - brent (property)")
    { check_against_driver_on_cubics(r::brent {}); }

    TEST_CASE("steps_view: yields the driver's iterates up to its stopping point - brent with width_tol (property)")
    { check_against_driver_on_cubics(r::brent { nxx::width_tol { 1e-6 } }); }

    TEST_CASE("steps_view: yields the driver's iterates up to its stopping point - bisection (property)")
    { check_against_driver_on_cubics(r::bisection {}); }

    TEST_CASE("steps_view: yields the driver's iterates up to its stopping point - bisection with width_tol (property)")
    { check_against_driver_on_cubics(r::bisection { nxx::width_tol { 1e-6 } }); }

    TEST_CASE("steps_view: a failing init yields one error element and ends")
    {
        // A NaN at the guess: the secant's init fails with non_finite_value after one evaluation.
        const auto nan_at_0 = [](double x) { return x == 0.0 ? std::numeric_limits<double>::quiet_NaN() : x - 1.0; };
        const auto sec      = r::secant {};
        const auto p        = sec.prepare(std::cref(nan_at_0), 0.0);
        if (!p) {
            FAIL_CHECK("prepare only checks that the guess is finite");
            return;
        }
        std::vector<nxx::fault<>> faults;
        std::size_t               oks = 0;
        for (const auto& st : nxx::steps_view { sec, *p } | std::views::take(10)) {
            if (st)
                ++oks;
            else
                faults.push_back(st.error());
        }
        CHECK(oks == 0u);
        CHECK(faults.size() == 1u);
        if (faults.size() == 1) {
            CHECK(faults[0].code == nxx::errc::non_finite_value);
            CHECK(faults[0].evaluations == 1u);
        }
        const auto driver = sec(nan_at_0, 0.0);    // the driver reports the same failure
        if (driver) { FAIL_CHECK("the secant cannot start at a NaN"); }
        else {
            CHECK(driver.error().code == nxx::errc::non_finite_value);
            CHECK(driver.error().used == nxx::counters { 0, 1 });
        }

        // A fallible callback: the fault keeps the user's error.
        const auto nt = r::newton {}.with_derivative(dsq2);
        const auto pg = nt.prepare(std::cref(g_neg), -1.0);
        if (!pg) {
            FAIL_CHECK("prepare only checks that the guess is finite");
            return;
        }
        std::vector<nxx::fault<step_error>> gfaults;
        for (const auto& st : nxx::steps_view { nt, *pg } | std::views::take(10)) {
            if (st) { FAIL_CHECK("g(-1) fails"); }
            else {
                gfaults.push_back(st.error());
            }
        }
        CHECK(gfaults.size() == 1u);
        if (gfaults.size() == 1) {
            CHECK(gfaults[0].code == nxx::errc::callback_failed);
            CHECK(gfaults[0].evaluations == 1u);
            CHECK(gfaults[0].cause == std::optional { step_error::domain });
        }
    }

    TEST_CASE("steps_view: a failing step ends the range after its error element")
    {
        const auto nt = r::newton {}.with_derivative(dsq2);
        const auto p  = nt.prepare(std::cref(sq2), 0.0);    // f'(0) = 0: the first step fails
        if (!p) {
            FAIL_CHECK("precondition: 0 is a finite guess");
            return;
        }
        std::vector<newton_view::element> elems;
        for (const auto& st : nxx::steps_view { nt, *p } | std::views::take(10)) elems.push_back(st);
        CHECK(elems.size() == 2u);
        if (elems.size() == 2) {
            CHECK(elems[0].has_value());
            if (elems[0]) {
                CHECK(elems[0]->x == 0.0);
                CHECK(elems[0]->fx == -2.0);
            }
            CHECK_FALSE(elems[1].has_value());
            if (!elems[1]) {
                CHECK(elems[1].error().code == nxx::errc::zero_derivative);
                CHECK(elems[1].error().evaluations == 1u);    // the derivative evaluation
            }
        }
    }

    TEST_CASE("steps_view: a step's input code becomes non_finite_value, with its evaluations and cause")
    {
        // detail::step_fault (DESIGN §6.3, §6.7), which the driver and steps_view apply to every step's fault: only the
        // two input codes Numerixx's own callables produce are mapped; every other code passes unchanged.
        // tests/usage/test_composition.cpp has the end-to-end rows for nested callables (mapped inside nxx::evaluate),
        // because they need the deriv module.
        using nxx::errc;
        constexpr auto mapped = nxx::detail::step_fault(nxx::fault<int> { errc::invalid_input, 3, 7 });
        static_assert(mapped.code == errc::non_finite_value && mapped.evaluations == 3 && mapped.cause == std::optional { 7 });
        static_assert(nxx::detail::step_fault(nxx::fault<> { errc::non_finite_input, 2, {} }).code == errc::non_finite_value);
        static_assert(nxx::detail::step_fault(nxx::fault<> { errc::non_finite_input, 2, {} }).evaluations == 2);
        static_assert(nxx::detail::step_fault(nxx::fault<> { errc::callback_failed, 1, {} }).code == errc::callback_failed);
        static_assert(nxx::detail::step_fault(nxx::fault<> { errc::zero_derivative, 1, {} }).code == errc::zero_derivative);
        static_assert(nxx::detail::step_fault(nxx::fault<> { errc::out_of_domain, 1, {} }).code == errc::out_of_domain);
        static_assert(!noexcept(nxx::detail::step_fault(std::declval<nxx::fault<int>>())));    // a user cause may throw (D10)

        // A user-written step that returns an input code directly (nxx::evaluate already maps a callback's, DESIGN §6.4):
        // detail::checked_step is the backstop. The code is mapped, the cause and the evaluations are kept, in the view
        // and in the driver alike.
        const auto raw = raw_input_step { r::newton {}.with_derivative(dsq2) };
        const auto p   = raw.prepare(std::cref(g_neg), 1.0);
        if (!p) {
            FAIL_CHECK("prepare only checks that the guess is finite");
            return;
        }
        std::vector<std::expected<r::newton_state<double>, nxx::fault<step_error>>> elems;
        for (const auto& st : nxx::steps_view { raw, *p } | std::views::take(10)) elems.push_back(st);
        CHECK(elems.size() == 2u);
        if (elems.size() == 2 && !elems[1]) {
            CHECK(elems[1].error().code == nxx::errc::non_finite_value);
            CHECK(elems[1].error().evaluations == 1u);
            CHECK(elems[1].error().cause == std::optional { step_error::domain });
        }
        const auto driver = nxx::iterate(raw, *p);
        CHECK_FALSE(driver.has_value());
        if (!driver) {
            CHECK(driver.error().code == nxx::errc::non_finite_value);
            CHECK(driver.error().used == nxx::counters { 1, 2 });    // f(1) in init, then the step's one evaluation
            CHECK(driver.error().cause == std::optional { step_error::domain });
            CHECK(driver.error().best.has_value());
        }
    }

    TEST_CASE("steps_view: no sign change is found by prepare before a view exists")
    {
        // Bracketing methods sample and check the bracket in prepare() (DESIGN §6.5, §6.6), so their init cannot fail and
        // a steps_view never sees this error.
        const auto pb = r::brent {}.prepare(std::cref(sq2), std::pair { 3.0, 4.0 });
        CHECK_FALSE(pb.has_value());
        if (!pb) {
            CHECK(pb.error().code == nxx::errc::no_sign_change);
            CHECK(pb.error().by == r::algos::brent);
            CHECK(pb.error().used == nxx::counters { 0, 2 });
            CHECK(pb.error().best.has_value());
        }
        const auto pi = r::bisection {}.prepare(std::cref(sq2), std::pair { 1.0, 1.0 });
        CHECK_FALSE(pi.has_value());
        if (!pi) {
            CHECK(pi.error().code == nxx::errc::invalid_input);
            CHECK(pi.error().used == nxx::counters {});
        }
    }

    TEST_CASE("steps_view: composes with std::views::transform from states to estimates")
    {
        const auto nt = r::newton {}.with_derivative(dsq2);
        const auto pn = nt.prepare(std::cref(sq2), 3.0);
        if (!pn) {
            FAIL_CHECK("precondition: 3 is a finite guess");
            return;
        }
        const auto to_estimate = [&nt](const auto& elem) {
            return elem.transform([&nt](const auto& state) { return nt.estimate(state); });
        };
        auto estimates    = nxx::steps_view { nt, *pn } | std::views::transform(to_estimate) | std::views::take(4);
        using estimates_t = decltype(estimates);
        static_assert(std::ranges::input_range<estimates_t>);
        static_assert(std::ranges::view<estimates_t>);
        static_assert(std::is_same_v<std::ranges::range_value_t<estimates_t>, std::expected<est_t, nxx::fault<>>>);

        std::vector<double> xs;
        for (const auto& e : estimates) {
            if (e)
                xs.push_back(e->x);
            else
                FAIL_CHECK("no error on the way to sqrt(2) from 3");
        }
        std::vector<double> direct;
        for (const auto& st : nxx::steps_view { nt, *pn } | std::views::take(4))
            if (st) direct.push_back(nt.estimate(*st).x);

        CHECK(xs.size() == 4u);
        CHECK(xs == direct);
        if (xs.size() != 4) return;
        CHECK(xs[0] == 3.0);    // element 0: the guess
        const double root = std::sqrt(2.0);
        CHECK(std::abs(xs[3] - root) < std::abs(xs[2] - root));
        CHECK(std::abs(xs[2] - root) < std::abs(xs[1] - root));
        CHECK(std::abs(xs[1] - root) < std::abs(xs[0] - root));

        // The driver takes the same steps.
        std::vector<double> seen;
        const auto          observed = nt.with_observer([&seen](const auto& v) { seen.push_back(v.x()); });
        const auto          driver   = observed(sq2, 3.0);
        CHECK(driver.has_value());
        if (seen.size() >= 3) {
            CHECK(same_bits(seen[0], xs[1]));
            CHECK(same_bits(seen[1], xs[2]));
            CHECK(same_bits(seen[2], xs[3]));
        }
        else {
            FAIL_CHECK("newton from 3 takes more than three steps");
        }
    }
}

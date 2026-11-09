// Numerixx demo: roots of a function of one variable with nxx::roots, and methods of your own.

#include <numerixx/deriv.hpp>
#include <numerixx/roots.hpp>

#include <array>
#include <cmath>
#include <expected>
#include <format>
#include <limits>
#include <numbers>
#include <optional>
#include <print>
#include <string_view>
#include <tuple>
#include <utility>

namespace {

    namespace rt = nxx::roots;

    template<typename T>
    void show(std::string_view label, std::expected<rt::Root<T>, rt::error> const &r)
    {
        if (r)
            std::println("   {:<30} {:.17g} +- {:.1e}", label, r->value, r->abs_error);
        else
            std::println("   {:<30} error: {}", label, r.error());
    }

    // fn with a count of its evaluations, which the methods do not report.
    template<typename F>
    struct Counted
    {
        F fn;
        int *count;

        auto operator()(double x) const // a value, or a tuple of the value and its derivatives
        {
            ++*count;
            return fn(x);
        }
    };

    void counted_result(std::string_view label, std::expected<rt::Root<double>, rt::error> const &r, int n)
    {
        if (r)
            std::println("   {:<30} {:.17g} +- {:.1e}  ({} evaluations)", label, r->value, r->abs_error, n);
        else
            std::println("   {:<30} error: {}  ({} evaluations)", label, r.error(), n);
    }

    template<typename M, typename F, typename S>
    void count_from(std::string_view label, F const &fn, S const &start, rt::Criteria<double> c)
    {
        int n = 0;
        auto const r = rt::solve<M>(Counted<F>{ fn, &n }, start, c);
        counted_result(label, r, n);
    }

    // The same with the method given as a value, for a method with parameters.
    template<typename M, typename F, typename S>
    void count_with(std::string_view label, M const &m, F const &fn, S const &start, rt::Criteria<double> c = {})
    {
        int n = 0;
        auto const r = rt::solve(m, Counted<F>{ fn, &n }, start, c);
        counted_result(label, r, n);
    }

    template<typename M, typename F>
    void count(std::string_view label, F const &fn, rt::Bracket<double> const &bracket, rt::Criteria<double> c = {})
    { count_from<M>(label, fn, bracket, c); }

    template<typename M, typename F>
    void count(std::string_view label, F const &fn, double guess, rt::Criteria<double> c = {})
    { count_from<M>(label, fn, guess, c); }

    void first_root()
    {
        std::println("1. sqrt(2) as the root of x^2 - 2 in [1, 2], by Brent's method");
        auto const f = [](double x) { return x * x - 2.0; };
        auto const r = rt::solve<rt::Brent>(f, rt::Bracket{ 1.0, 2.0 });
        show("root", r);
        std::println("   {:<30} {:.17g}", "exact", std::numbers::sqrt2);
        // abs_error bounds the distance to the point where fn, as evaluated, changes sign.
        if (r) std::println("   {:<30} {:.1e}", "actual error", std::abs(r->value - std::numbers::sqrt2));
    }

    void methods()
    {
        // Bisection gains one bit per evaluation, whatever fn is. Brent interpolates and is much faster on a smooth fn
        // with a simple root, but slower on a multiple root, where the interpolation converges only linearly. Toms748
        // costs about what Brent costs, with a better worst case; AndersonBjorck wins where fn is strongly curved near
        // one end of the bracket, as ln(x / 1e-12) is near 1e-20.
        std::println("\n2. Bisection, Brent, Toms748 and AndersonBjorck");
        auto const dottie = [](double x) { return std::cos(x) - x; };
        auto const wallis = [](double x) { return (x * x - 2.0) * x - 5.0; };
        auto const triple = [](double x) { return (x - 1.0) * (x - 1.0) * (x - 1.0); };
        count<rt::Bisection>("cos x = x: Bisection", dottie, { 0.0, 1.0 });
        count<rt::Brent>("cos x = x: Brent", dottie, { 0.0, 1.0 });
        count<rt::Bisection>("x^3 - 2x - 5: Bisection", wallis, { 2.0, 3.0 });
        count<rt::Brent>("x^3 - 2x - 5: Brent", wallis, { 2.0, 3.0 });
        count<rt::Bisection>("(x - 1)^3: Bisection", triple, { 0.0, 3.0 });
        count<rt::Brent>("(x - 1)^3: Brent", triple, { 0.0, 3.0 });
        count<rt::Toms748>("cos x = x: Toms748", dottie, { 0.0, 1.0 });
        count<rt::AndersonBjorck>("cos x = x: AndersonBjorck", dottie, { 0.0, 1.0 });
        auto const trace = [](double x) { return std::log(x / 1e-12); };
        count<rt::Brent>("ln(x / 1e-12): Brent", trace, { 1e-20, 1.0 });
        count<rt::Toms748>("ln(x / 1e-12): Toms748", trace, { 1e-20, 1.0 });
        count<rt::AndersonBjorck>("ln(x / 1e-12): AndersonBjorck", trace, { 1e-20, 1.0 });

        // The teaching methods on the triple root: Ridders' exponential fit is the cheapest of all there; false position
        // (RegulaFalsi) stalls at one end, which the bisection guard cuts off, and Illinois halves the stalled value.
        count<rt::Ridders>("(x - 1)^3: Ridders", triple, { 0.0, 3.0 });
        count<rt::RegulaFalsi>("(x - 1)^3: RegulaFalsi", triple, { 0.0, 3.0 });
        count<rt::Illinois>("(x - 1)^3: Illinois", triple, { 0.0, 3.0 });
    }

    void kepler()
    {
        // Kepler's equation E - e sin E = M gives a planet's eccentric anomaly E from its mean anomaly M. For e < 1,
        // f(0) = -M and f(pi) = pi - M bracket the root for every M in (0, pi).
        std::println("\n3. Kepler's equation E - e sin E = M for e = 0.9, a table over M");
        double const e = 0.9;
        for (double const m : { 0.1, 0.5, 1.0, 2.0, 3.0 }) {
            auto const fn = [e, m](double big_e) { return big_e - e * std::sin(big_e) - m; };
            show(std::format("M = {:.1f}", m), rt::solve<rt::Brent>(fn, rt::Bracket{ 0.0, std::numbers::pi }));
        }
    }

    void tolerances()
    {
        // The default criteria give the full precision of the type at every magnitude: a root at 1e-12 keeps its 16
        // digits, with no tolerance to set. A looser tolerance saves evaluations.
        std::println("\n4. Tolerances");
        auto const trace = [](double x) { return std::log(x / 1e-12); };
        count<rt::Brent>("ln(x / 1e-12): default", trace, { 1e-20, 1.0 });
        auto const f = [](double x) { return x * x - 2.0; };
        count<rt::Brent>("x^2 - 2: default", f, { 1.0, 2.0 });
        count<rt::Brent>("x^2 - 2: rel_tol = 1e-6", f, { 1.0, 2.0 }, { .rel_tol = 1e-6 });

        // A sign change near 0 that fn misses is approached down to the type's smallest numbers: set abs_tol there.
        auto const tiny = [](double x) { return x - 1e-300; };
        count<rt::Bisection>("x - 1e-300: default", tiny, { 0.0, 1.0 });
        count<rt::Bisection>("x - 1e-300: abs_tol = 1e-12", tiny, { 0.0, 1.0 }, { .abs_tol = 1e-12 });
    }

    void around_zero()
    {
        // A bracket around 0 is first split at 0, so a root at 0 costs one evaluation more.
        std::println("\n5. A bracket around 0");
        auto const sin = [](double x) { return std::sin(x); };
        count<rt::Brent>("sin x in [-1, 2]", sin, { -1.0, 2.0 });
        // Where fn is not finite at 0 (a removable singularity), that split is skipped.
        auto const sinc = [](double x) { return std::sin(x) / x - 0.5; };
        count<rt::Brent>("sin(x)/x = 0.5 in [-1, 3]", sinc, { -1.0, 3.0 });
    }

    template<typename T>
    void precision(std::string_view type)
    {
        auto const f = [](T x) { return x * x - T(2); };
        if (auto const r = rt::solve<rt::Brent>(f, rt::Bracket{ T(1), T(2) }))
            std::println("   {:<12} abs_error {:.1e}, error {:.1e}",
                type,
                static_cast<double>(r->abs_error),
                static_cast<double>(std::abs(r->value - std::sqrt(T(2)))));
    }

    void types()
    {
        // The error is a few units of the type's precision. MSVC makes long double the same as double.
        std::println("\n6. sqrt(2) in float, double and long double");
        precision<float>("float");
        precision<double>("double");
        precision<long double>("long double");
    }

    void compile_time()
    {
        std::println("\n7. At compile time");
        constexpr auto f = [](double x) { return x * x - 2.0; };
        constexpr auto r = rt::solve<rt::Brent>(f, rt::Bracket{ 1.0, 2.0 });
        static_assert(r.has_value() && nxx::math::abs(r->value - 1.4142135623730951) <= r->abs_error);
        show("constexpr sqrt(2)", r);
        // Bisection reaches a dyadic root exactly: 0.5, 0.25, then 0.375.
        static_assert(
            rt::solve<rt::Bisection>([](double x) { return x - 0.375; }, rt::Bracket{ 0.0, 1.0 })->value == 0.375);
    }

    void errors()
    {
        std::println("\n8. Failures are values, not exceptions");
        auto const f = [](double x) { return x * x - 2.0; };
        show("x^2 - 2 in [2, 3]", rt::solve<rt::Brent>(f, rt::Bracket{ 2.0, 3.0 })); // the same sign at both ends
        auto const log = [](double x) { return std::log(x) + 1.0; };
        show("log x + 1 in [0, 1]", rt::solve<rt::Brent>(log, rt::Bracket{ 0.0, 1.0 })); // log 0 is -inf
        show("x^2 - 2, rel_tol = -1", rt::solve<rt::Brent>(f, rt::Bracket{ 1.0, 2.0 }, { .rel_tol = -1.0 }));

        // A pole changes sign like a root, and is returned as one: check |fn| at the result where fn may have poles.
        auto const pole = [](double x) { return 1.0 / (x - 1.0 / 3.0); };
        auto const p = rt::solve<rt::Brent>(pole, rt::Bracket{ 0.0, 1.0 });
        show("1 / (x - 1/3) in [0, 1]", p);
        if (p) std::println("   {:<30} {:.1e}", "|fn| there", std::abs(pole(p->value)));
    }

    void composition()
    {
        // Results compose through std::expected: the root of one solve sets up the next, and an error at any stage
        // stops the chain.
        std::println("\n9. Composition: the x where x e^x = 1, then the y where y^2 = that x");
        auto const omega = [](double x) { return x * std::exp(x) - 1.0; };
        auto const both = rt::solve<rt::Brent>(omega, rt::Bracket{ 0.0, 1.0 }).and_then([](rt::Root<double> const &w) {
            auto const g = [x = w.value](double y) { return y * y - x; };
            return rt::solve<rt::Brent>(g, rt::Bracket{ 0.0, 1.0 });
        });
        show("sqrt(W(1))", both);
        std::println("   {:<30} {:.17g}", "via std::sqrt", std::sqrt(0.56714329040978387300));
    }

    void open_methods()
    {
        // Newton and Halley take fn's derivatives from fn, which returns them with the value, or from callables of
        // their own; Secant needs none. From a guess, each iterates until two iterates bracket a sign change, then
        // finishes inside that bracket, so every result has a bound like Brent's.
        std::println("\n10. Newton, Halley and Secant");
        auto const fdf = [](double x) { return std::tuple{ x * x - 2.0, 2.0 * x, 2.0 }; }; // f, f', f''
        count<rt::Newton>("x^2 - 2: Newton from 1", fdf, 1.0);
        count<rt::Halley>("x^2 - 2: Halley from 1", fdf, 1.0);
        count<rt::Secant>("x^2 - 2: Secant from 1", fdf, 1.0);
        count<rt::Newton>("x^2 - 2: Newton on [1, 2]", fdf, { 1.0, 2.0 });

        // The derivative as a callable of its own: the same result.
        auto const f = [](double x) { return std::exp(x) - 2.0; };
        auto const df = [](double x) { return std::exp(x); };
        show("e^x = 2: Newton, f and f'", rt::solve<rt::Newton>(f, df, 0.0));

        // Kepler's equation again, now from a guess: E = M is a common start.
        auto const kepler = [](double big_e) {
            return std::pair{ big_e - 0.9 * std::sin(big_e) - 1.0, 1.0 - 0.9 * std::cos(big_e) };
        };
        count<rt::Newton>("Kepler, M = 1: Newton from M", kepler, 1.0);

        // An open method can fail where a bracket cannot: Newton cycles between 0 and 1 on x^3 - 2x + 2.
        auto const cycle = [](double x) { return std::pair{ (x * x - 2.0) * x + 2.0, 3.0 * x * x - 2.0 }; };
        count<rt::Newton>("x^3 - 2x + 2: Newton from 0", cycle, 0.0);
        count<rt::Newton>("x^3 - 2x + 2: Newton, [-2, -1]", cycle, { -2.0, -1.0 });

        // Steffensen: Newton without a derivative, at the price of scale dependence. f(x) is added to x, so scaling f
        // changes the cost and even the root found (1000 (x^2 - 2) from 1: x + f(x) = -999 brackets -sqrt 2).
        auto const sq = [](double x) { return x * x - 2.0; };
        auto const sq1000 = [](double x) { return 1000.0 * (x * x - 2.0); };
        count<rt::Steffensen>("x^2 - 2: Steffensen from 1", sq, 1.0);
        count<rt::Steffensen>("1000 (x^2 - 2): Steffensen", sq1000, 1.0);
        count<rt::Secant>("1000 (x^2 - 2): Secant", sq1000, 1.0);
    }

    void numerical_derivatives()
    {
        // Where f' is inconvenient to derive, nxx::deriv computes it: derivative_of(f) is a callable whose call at x
        // is derivative(f, x). Each point then costs the derivative's calls of f too (Central<2>: 2 more), so for a
        // value-only f, Secant and Brent are usually cheaper. Every call of f is counted here.
        namespace dv = nxx::deriv;
        std::println("\n11. Numerical derivatives (calls of f)");
        auto const report = [](std::string_view label, std::expected<rt::Root<double>, rt::error> const &r, int n) {
            if (r)
                std::println("   {:<30} {:.17g} +- {:.1e}  ({} calls)", label, r->value, r->abs_error, n);
            else
                std::println("   {:<30} error: {}  ({} calls)", label, r.error(), n);
        };
        auto const kepler = [](double big_e) { return big_e - 0.9 * std::sin(big_e) - 1.0; };
        int n = 0;
        Counted<decltype(kepler)> const f{ kepler, &n };
        auto r = rt::solve<rt::Newton>(f, dv::derivative_of(f), 1.0);
        report("Kepler: Newton, Central<2>", r, n);
        n = 0;
        r = rt::solve<rt::Halley>(f, dv::derivative_of(f), dv::derivative_of<dv::Central<2>, 2>(f), 1.0);
        report("Kepler: Halley, Central<2>", r, n);
        n = 0;
        r = rt::solve<rt::Newton>(f, dv::derivative_of<dv::Ridders>(f), 1.0);
        report("Kepler: Newton, deriv Ridders", r, n);
        n = 0;
        r = rt::solve<rt::Secant>(f, 1.0);
        report("Kepler: Secant", r, n);

        // A numerical derivative fails where its stencil leaves f's domain. A failure counts as a derivative that is
        // not finite: from a guess there is no step, on a bracket a bisection step. ln(x / 1e-12) varies on a scale
        // far below the default step's max(|x|, 1) eps^(1/3) = 6e-6.
        auto const ln = [](double x) { return std::log(x / 1e-12); };
        Counted<decltype(ln)> const g{ ln, &n };
        n = 0;
        r = rt::solve<rt::Newton>(g, dv::derivative_of(g), 1e-13);
        report("ln: Newton from 1e-13", r, n);
        n = 0;
        r = rt::solve<rt::Newton>(g, dv::derivative_of(g), rt::Bracket{ 1e-20, 1.0 });
        report("ln: Newton on [1e-20, 1]", r, n);
        n = 0;
        r = rt::solve<rt::Secant>(g, 1e-13);
        report("ln: Secant from 1e-13", r, n);
    }

    void searching()
    {
        // Without a bracket, search(fn, x0) finds one: fn at x0, then at x0 -+ h for h = |x0| / 50 growing by 3/2,
        // until the sign differs from fn(x0). A bracketing method then finishes, through and_then.
        std::println("\n12. Searching for a bracket from a guess");
        auto const show_bracket = [](std::string_view label, std::expected<rt::Bracket<double>, rt::error> const &b) {
            if (b)
                std::println("   {:<30} [{:.6g}, {:.6g}]", label, b->lower, b->upper);
            else
                std::println("   {:<30} error: {}", label, b.error());
        };
        auto const psat = [](double t) { return std::exp(20.0 - 5000.0 / t) - 1.0; };
        show_bracket("e^(20 - 5000/T) = 1 from 300", rt::search(psat, 300.0));
        show("... and Brent on it", rt::search(psat, 300.0).and_then([&](rt::Bracket<double> b) {
            return rt::solve<rt::Brent>(psat, b);
        }));

        // A domain edge: where fn is not finite, that side bisects towards its last finite point.
        show_bracket("ln(x / 1e-12) = 0 from 1", rt::search([](double x) { return std::log(x / 1e-12); }, 1.0));

        // van der Waals is finite between 0 and its co-volume b, with a pole at b, which a search would bracket.
        // Returning NaN outside the domain keeps the search inside it.
        double const b = 4.27e-5;
        auto const vdw = [b](double v) { return 8.314 * 280.0 / (v - b) - 0.364 / (v * v) - 4e6; };
        auto const vdw_domain = [&vdw, b](double v) {
            return v > b ? vdw(v) : std::numeric_limits<double>::quiet_NaN();
        };
        show_bracket("van der Waals from 1e-4", rt::search(vdw, 1e-4));
        show_bracket("... NaN outside v > b", rt::search(vdw_domain, 1e-4));

        // What an outward search cannot find: two roots within one step.
        show_bracket("x^2 = 2 from 100", rt::search([](double x) { return x * x - 2.0; }, 100.0));
    }

    // A method of your own on a bracket: a type with two const members, start and advance, each leaving a Step (the
    // bracket, fn's values at its ends and the method's point) in its state. solve() checks every Step, applies the
    // stop test and forms the bound, so the method has neither; the manual's "A method of your own" has the rules.
    // This one evaluates fn at the fraction `at` of the bracket from its lower end (0.5 is bisection) and keeps the
    // part with the sign change; its point is the end with the smaller |f|.
    struct Section
    {
        double at = 0.5;

        struct State
        {
            double a, fa, b, fb; // the bracket and fn's values at its ends
            rt::Step<double> step;
        };

        static constexpr rt::Step<double> report(State const &s)
        { return rt::Step<double>{ s.a, s.fa, s.b, s.fb, nxx::math::abs(s.fa) < nxx::math::abs(s.fb) ? s.a : s.b }; }

        template<typename F>
        constexpr std::expected<State, rt::error>
            start(F &, rt::Start<double> const &s, rt::Tolerance<double> const &) const
        {
            State state{ s.lower, s.at_lower[0], s.upper, s.at_upper[0], {} };
            state.step = report(state);
            return state;
        }

        template<typename F>
        constexpr std::expected<void, rt::error> advance(F &fn, State &s, rt::Tolerance<double> const &tol) const
        {
            // A weighted average of the ends, which cannot overflow, kept tol inside them: a point that approached the
            // root from one side only would never move the other end, and the bracket would not close.
            double const margin = tol.at(s.step.point);
            double x = (1.0 - at) * s.a + at * s.b;
            if (x < s.a + margin) x = s.a + margin;
            else if (s.b - margin < x) x = s.b - margin;
            if (!(s.a < x && x < s.b)) x = s.a / 2 + s.b / 2;
            if (!(s.a < x && x < s.b)) return std::unexpected(rt::error::not_converged);
            double const fx = fn(x);
            if (!nxx::math::isfinite(fx)) return std::unexpected(rt::error::non_finite);
            if (fx == 0.0) {
                s.step = rt::Step<double>{ s.a, s.fa, x, fx, x };
                return {};
            }
            if ((fx < 0.0) == (s.fa < 0.0)) { // the sign test without a product, which could underflow
                s.a = x;
                s.fa = fx;
            } else {
                s.b = x;
                s.fb = fx;
            }
            s.step = report(s);
            return {};
        }
    };
    static_assert(rt::bracketing_method<Section, double>); // names the member at fault, where solve() would not

    // A method that breaks the contract: its advance leaves the bracket as it was.
    struct Stuck
    {
        struct State
        {
            rt::Step<double> step;
        };

        template<typename F>
        constexpr std::expected<State, rt::error>
            start(F &, rt::Start<double> const &s, rt::Tolerance<double> const &) const
        { return State{ rt::Step<double>{ s.lower, s.at_lower[0], s.upper, s.at_upper[0], s.lower } }; }

        template<typename F>
        constexpr std::expected<void, rt::error> advance(F &, State &, rt::Tolerance<double> const &) const
        { return {}; }
    };

    void own_bracketing_method()
    {
        std::println("\n13. A method of your own on a bracket");
        auto const f = [](double x) { return x * x - 2.0; };
        count<rt::Bisection>("x^2 - 2: Bisection", f, { 1.0, 2.0 });
        count<Section>("x^2 - 2: Section", f, { 1.0, 2.0 }); // named: Section{}, at = 0.5
        count_with("x^2 - 2: Section, at = 0.1", Section{ .at = 0.1 }, f, rt::Bracket{ 1.0, 2.0 });
        auto const dottie = [](double x) { return std::cos(x) - x; };
        count_with("cos x = x: Section, at = 0.1", Section{ .at = 0.1 }, dottie, rt::Bracket{ 0.0, 1.0 });

        // observer(o) before fn: o sees each Step, here at a loose tolerance, so that there are few.
        std::println("   the steps of Section at rel_tol = 1e-2:");
        auto const print = [](rt::Step<double> const &s) {
            std::println("     [{:.6f}, {:.6f}] point {:.6f}", s.lower, s.upper, s.point);
        };
        (void)rt::solve<Section>(rt::observer(print), f, rt::Bracket{ 1.0, 2.0 }, { .rel_tol = 1e-2 });

        // solve() checks each Step: a bracket that does not shrink ends the solve with invalid_step, so even a broken
        // method cannot loop for ever.
        show("x^2 - 2: Stuck", rt::solve<Stuck>(f, rt::Bracket{ 1.0, 2.0 }));

        // A constexpr method runs at compile time, as the library's do.
        constexpr auto r = rt::solve<Section>([](double x) { return x * x - 2.0; }, rt::Bracket{ 1.0, 2.0 });
        static_assert(r.has_value() && nxx::math::abs(r->value - 1.4142135623730951) <= r->abs_error);
    }

    // A method of your own from a guess (values = 2, so fn gives f and f'): Newton's step, at least tol(x) long, until
    // two iterates show a sign change; then Brent's method on them, whose state it holds. solve() counts the steps
    // without a sign change against its cap, so this method ends even where Newton cycles or diverges.
    struct NewtonThenBrent
    {
        static constexpr int values = 2;

        struct State
        {
            double x;              // the latest iterate
            std::array<double, 2> v; // f and f' there
            std::optional<rt::bracketing_state_t<rt::Brent, double>> brent; // after the first sign change
            rt::Step<double> step;
        };

        template<typename F>
        static constexpr std::array<double, 2> values_at(F &fn, double x)
        {
            auto const r = fn(x);
            return { std::get<0>(r), std::get<1>(r) };
        }

        template<typename F>
        constexpr std::expected<State, rt::error> start(F &fn, double const &guess, rt::Tolerance<double> const &) const
        {
            std::array<double, 2> const v = values_at(fn, guess);
            if (!nxx::math::isfinite(v[0])) return std::unexpected(rt::error::non_finite);
            return State{ guess, v, std::nullopt, rt::Step<double>{ guess, v[0], guess, v[0], guess } };
        }

        template<typename F>
        constexpr std::expected<void, rt::error> advance(F &fn, State &s, rt::Tolerance<double> const &tol) const
        {
            if (s.brent) { // Brent's step, and its report
                std::expected<void, rt::error> const r = rt::Brent{}.advance(fn, *s.brent, tol);
                s.step = s.brent->step;
                return r;
            }
            if (s.v[1] == 0.0 || !nxx::math::isfinite(s.v[1])) return std::unexpected(rt::error::not_converged);
            double step = -s.v[0] / s.v[1];
            double const least = tol.at(s.x);
            if (nxx::math::abs(step) < least) step = step < 0.0 ? -least : least;
            double const next = s.x + step;
            std::array<double, 2> const vn = values_at(fn, next);
            if (!nxx::math::isfinite(vn[0])) return std::unexpected(rt::error::non_finite);
            bool const up = s.x < next;
            if (vn[0] != 0.0 && (vn[0] < 0.0) != (s.v[0] < 0.0)) { // a sign change: Brent from here
                rt::Start<double> const ends = up ? rt::Start<double>{ s.x, { s.v[0] }, next, { vn[0] } }
                                                  : rt::Start<double>{ next, { vn[0] }, s.x, { s.v[0] } };
                auto started = rt::Brent{}.start(fn, ends, tol);
                if (!started) return std::unexpected(started.error());
                s.brent = *started;
                s.step = s.brent->step;
                return {};
            }
            // Before it: the last two iterates in order, the latest the point.
            s.step = up ? rt::Step<double>{ s.x, s.v[0], next, vn[0], next }
                        : rt::Step<double>{ next, vn[0], s.x, s.v[0], next };
            s.x = next;
            s.v = vn;
            return {};
        }
    };
    static_assert(rt::open_method<NewtonThenBrent, double>);

    void own_open_method()
    {
        std::println("\n14. A method of your own from a guess");
        // values = 2 gives the method the forms with f' as a separate callable, as for Newton.
        auto const f = [](double x) { return std::exp(x) - 2.0; };
        auto const df = [](double x) { return std::exp(x); };
        show("e^x = 2: own, f and f'", rt::solve<NewtonThenBrent>(f, df, 0.0));
        show("e^x = 2: Newton, f and f'", rt::solve<rt::Newton>(f, df, 0.0));
        // As a value, or with an observer, fn gives both values. Brent's interpolation finishes in two evaluations
        // more than the library's Newton, which keeps its own step inside the bracket.
        auto const kepler = [](double big_e) {
            return std::pair{ big_e - 0.9 * std::sin(big_e) - 1.0, 1.0 - 0.9 * std::cos(big_e) };
        };
        count_with("Kepler, M = 1: own from M", NewtonThenBrent{}, kepler, 1.0);
        count<rt::Newton>("Kepler, M = 1: Newton from M", kepler, 1.0);
        // No root: x^2 + 1 steps from 1 to 0, where f' = 0; from 0.5 it wanders until the cap.
        auto const no_root = [](double x) { return std::pair{ x * x + 1.0, 2.0 * x }; };
        count_with("x^2 + 1: own from 1", NewtonThenBrent{}, no_root, 1.0);
        count_with("x^2 + 1: own from 0.5", NewtonThenBrent{}, no_root, 0.5);
    }

} // namespace

int main()
{
    first_root();
    methods();
    kepler();
    tolerances();
    around_zero();
    types();
    compile_time();
    errors();
    composition();
    open_methods();
    numerical_derivatives();
    searching();
    own_bracketing_method();
    own_open_method();
}

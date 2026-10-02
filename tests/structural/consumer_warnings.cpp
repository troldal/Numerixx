// A consumer translation unit compiled with strict warnings as errors (DESIGN §9.1, §5.3). The global `f` makes
// MSVC's C4459 ("declaration hides global declaration") fire in any Numerixx header that names a parameter `f`,
// so the library names callback parameters `fn` or `func`.
//
// Including the headers is not enough: most warnings (C4459, C4100, C4127, C4702, -Wshadow, -Wconversion, ...) are
// reported where a template is instantiated. So the functions below instantiate the library with `f` in scope, the
// way a consumer does: static chains (first_of, then, warm_fallback, first_of_with), Newton with derivative_of and
// with the numeric policy, steps_view, projection and observer, criteria algebra, fallible callbacks, run-time chains,
// finite differences, the FXT pipes, and float and long double. This is an object library: nothing here runs.
#include <numerixx/core/any_solver.hpp>
#include <numerixx/numerixx.hpp>

#include <cstdint>
#include <expected>
#include <functional>
#include <limits>
#include <ranges>
#include <utility>
#include <vector>

double f(double x);
double f(double x) { return x * x - 2.0; }

namespace consumer
{
    namespace r = nxx::roots;

    double df(double x);
    double df(double x) { return 2.0 * x; }

    // The user's own callback error, with the is_fatal customisation found by ADL.
    enum class eval_error : std::uint8_t { domain, fatal };
    constexpr bool is_fatal(eval_error e) noexcept { return e == eval_error::fatal; }

    int    uses_numerixx();
    double solve_with_the_global_function();
    double three_solver_chain();
    double staged_pipeline();
    double newton_with_derivatives();
    double manual_stepping();
    int    projection_and_observer();
    bool   criteria_algebra();
    int    fallible_callbacks();
    double combinators();
    double run_time_chain();
    double finite_differences();
    double other_scalars();
    double pipes();

    int uses_numerixx() { return nxx::version.major + static_cast<int>(f(2.0)); }

    // The global function itself as the callback, by reference and by pointer (canonical call 1).
    double solve_with_the_global_function()
    {
        const auto by_reference = r::brent {}(f, { 1.0, 2.0 });
        const auto by_pointer   = r::brent {}(&f, { 1.0, 2.0 });
        const auto facade       = nxx::roots::solve(f, { 1.0, 2.0 });
        const auto curried      = r::bisection {}.on({ 1.0, 2.0 })(f);
        return nxx::best_x(by_reference).value_or(0.0) + nxx::best_x(by_pointer).value_or(0.0) + nxx::best_x(facade).value_or(0.0) +
               nxx::best_x(curried).value_or(0.0);
    }

    // Canonical call 9: newton, secant and brent (three state types) with run-time inputs.
    double three_solver_chain()
    {
        const double x0    = 1.0;
        const double lo    = 1.0;
        const double hi    = 2.0;
        const auto   chain = nxx::first_of(r::newton {}.with_derivative(df).on(x0), r::secant {}.on(x0), r::brent {}.on({ lo, hi }));
        const auto   res   = chain(f);
        return res ? res->x : res.error().best ? res.error().best->x : std::numeric_limits<double>::quiet_NaN();
    }

    // expand -> bisection -> newton, and canonical call 10.
    double staged_pipeline()
    {
        const auto pipeline = nxx::then(r::expand {}.on(nxx::bracket { 1.0, 1.2 }),
                                        r::bisection { nxx::width_tol { 1e-4 } },
                                        r::newton {}.with_derivative(df));
        const auto coarse   = nxx::then(r::bisection { nxx::width_tol { 1e-3 } }.with_budget(100).on({ 1.0, 2.0 }), r::secant {});
        return nxx::best_x(pipeline(f)).value_or(0.0) + nxx::best_x(coarse(f)).value_or(0.0);
    }

    // Newton with every derivative source: a callable, derivative_of, the numeric policy (also curried).
    double newton_with_derivatives()
    {
        const auto analytic = r::newton {}.with_derivative(df)(f, 1.0);
        const auto of       = r::newton {}.with_derivative(nxx::deriv::derivative_of(f))(f, 1.0);
        const auto of4 =
            r::newton {}.with_derivative(nxx::deriv::derivative_of(f, nxx::deriv::central_1_4, nxx::deriv::relative { 1e-4, 1.0 }))(f, 1.0);
        const auto numeric = r::newton {}.with_derivative(nxx::deriv::numeric {})(f, 1.0);
        const auto curried = nxx::first_of(r::newton {}.with_derivative(nxx::deriv::numeric {}).on(1.0), r::brent {}.on({ 1.0, 2.0 }))(f);
        return nxx::best_x(analytic).value_or(0.0) + nxx::best_x(of).value_or(0.0) + nxx::best_x(of4).value_or(0.0) +
               nxx::best_x(numeric).value_or(0.0) + nxx::best_x(curried).value_or(0.0);
    }

    // steps_view over brent and newton (canonical call 12), also lazily transformed into estimates.
    double manual_stepping()
    {
        double     sum    = 0.0;
        const auto solver = r::brent {};
        const auto p      = solver.prepare(std::cref(f), nxx::bracket<double>::make(1.0, 2.0));
        if (p)
            for (const auto& st : nxx::steps_view { solver, *p } | std::views::take(20))
                if (st) sum += solver.estimate(*st).x;

        const auto nt = r::newton {}.with_derivative(df);
        const auto pn = nt.prepare(std::cref(f), 3.0);
        if (pn) {
            auto xs =
                nxx::steps_view { nt, *pn } |
                std::views::transform([&nt](const auto& st) { return st.transform([&nt](const auto& s) { return nt.estimate(s); }); }) |
                std::views::take(4);
            for (const auto& e : xs)
                if (e) sum += e->x;
        }
        return sum;
    }

    // Canonical call 3 and an observed, clamped Newton; builders in either order.
    int projection_and_observer()
    {
        int        calls   = 0;
        const auto counter = [&calls](const auto& view) { calls += view.x() > 0.0 ? 1 : 0; };
        const auto clamped = r::secant {}.with_projection(r::clamp_to { 0.0, 3.0 }).with_observer(counter)(f, 1.0);
        const auto newton  = r::newton {}.with_observer(counter).with_derivative(df).with_projection(r::clamp_to { 1.0, 3.0 })(f, 10.0);
        const auto brent   = r::brent {}.with_observer(counter)(f, { 1.0, 2.0 });
        const auto search  = r::expand {}.with_observer(counter)(f, { 1.0, 1.2 });
        const auto length  = clamped.transform([](const auto& s) { return s.x; });
        return calls + (length ? 1 : 0) + (newton ? 1 : 0) + (brent ? 1 : 0) + (search ? 1 : 0) + (nxx::best_x(clamped) ? 1 : 0);
    }

    // Stop criteria and their algebra, on the solvers they apply to.
    bool criteria_algebra()
    {
        const auto a = r::bisection { nxx::width_tol { 1e-10 } || nxx::max_evaluations { 20 } }(f, { 1.0, 2.0 });
        const auto b = r::secant { nxx::x_tol { 1e-12 } && nxx::min_iterations { 3 } }(f, 1.0);
        const auto c = r::newton { nxx::never {} || nxx::max_evaluations { 10 } }.with_derivative(df)(f, 1.0);
        const auto d = r::newton { nxx::x_tol { 1e-12, 1e-10 } || nxx::f_tol { 1e-14 } }.with_derivative(df)(f, 1.0);
        const auto e = r::brent { nxx::width_tol { 1e-12, 1e-10 } }(f, { 1.0, 2.0 });
        const auto g = r::bisection {}.with_stop(nxx::floored_width { 40 } || nxx::f_tol { 1e-9 }).with_budget(200)(f, { 1.0, 2.0 });
        const auto h = r::secant {}.with_stop(nxx::step_tol<1, 2> {} || nxx::f_tol { 1e-12 })(f, 1.0);

        // Run-time tolerances go through make().
        const auto tol = nxx::width_tol<double>::make(1e-9, 0.0);
        const bool i   = tol && r::bisection { *tol }(f, { 1.0, 2.0 }).has_value();
        return a.has_value() && b.has_value() && c.has_value() && d.has_value() && e.has_value() && g.has_value() && h.has_value() && i;
    }

    // Fallible callbacks keep the user's error through solvers, combinators and derivative_of (common cause).
    int fallible_callbacks()
    {
        const auto domain = [](double x) -> std::expected<double, eval_error> {
            if (x < 0.0) return std::unexpected(eval_error::domain);
            return f(x);
        };
        const auto chain = nxx::first_of(r::newton {}.with_derivative(nxx::deriv::derivative_of(domain)).on(0.0),
                                         r::secant {}.on(-1.0),
                                         nxx::then(r::expand {}.on(nxx::bracket { 1.0, 1.2 }), r::brent {}),
                                         r::bisection {}.on({ 0.0, 2.0 }));
        const auto res   = chain(domain);
        const auto bad   = r::bisection {}(domain, { -1.0, 2.0 });
        int        score = res ? 1 : 0;
        if (!bad && bad.error().cause && *bad.error().cause == eval_error::domain) score += 2;
        return score;
    }

    // warm_fallback, first_of_with and nested chains.
    double combinators()
    {
        const auto warm   = nxx::warm_fallback(r::bisection {}.with_budget(5).on({ 1.0, 2.0 }), r::newton {}.with_derivative(df))(f);
        const auto picky  = nxx::first_of_with([](const auto& e) { return !nxx::is_input_error(e.code); },
                                               r::secant {}.on(1.0),
                                               r::brent {}.on({ 1.0, 2.0 }))(f);
        const auto chain  = nxx::first_of(r::secant {}.with_budget(3).on(1.0), r::brent {}.on({ 1.0, 2.0 }));
        auto       copy   = chain;
        copy              = chain;    // chains are copy-assignable values
        const auto nested = nxx::first_of(copy, nxx::then(r::expand {}.on({ 1.0, 1.2 }), r::bisection {}))(f);
        return nxx::best_x(warm).value_or(0.0) + nxx::best_x(picky).value_or(0.0) + nxx::best_x(nested).value_or(0.0);
    }

    // Canonical call 11: a run-time chain over std::function, and a static chain converted to it.
    double run_time_chain()
    {
        using fn_t               = std::function<double(double)>;
        using solver_t           = nxx::any_solver<fn_t, r::root_estimate<double>>;
        const fn_t            fn = f;
        std::vector<solver_t> alternatives { r::secant {}.on(1.0), r::brent {}.on({ 1.0, 2.0 }) };
        const solver_t        chain = nxx::first_of(std::move(alternatives));
        const solver_t        mixed = nxx::first_of(r::newton {}.with_derivative(df).on(1.0), chain);
        const auto            picky = nxx::first_of_with(nxx::continue_unless_fatal {}, std::vector<solver_t> { mixed, chain });
        return nxx::best_x(chain(fn)).value_or(0.0) + nxx::best_x(mixed(fn)).value_or(0.0) + nxx::best_x(picky(fn)).value_or(0.0);
    }

    // One-shot differentiation with every stencil and step kind.
    double finite_differences()
    {
        const auto a = nxx::deriv::diff(f, 1.0);
        const auto b = nxx::deriv::central(f, 1.0, nxx::deriv::relative { 1e-5, 1.0 });
        const auto c = nxx::deriv::diff(f, 1.0, nxx::deriv::central_2_4, nxx::deriv::absolute { 1e-3 });
        const auto d = nxx::deriv::diff(f, 0.0, nxx::deriv::forward_1_1);
        const auto e = nxx::deriv::diff(f, 1.0, nxx::deriv::backward_1_1, nxx::deriv::optimal {});
        const auto g = nxx::deriv::diff(f, 1.0, nxx::deriv::central_2_2);
        return a.value_or(0.0) + b.value_or(0.0) + c.value_or(0.0) + d.value_or(0.0) + e.value_or(0.0) + g.value_or(0.0);
    }

    // float and long double through the same solvers.
    double other_scalars()
    {
        const auto ff  = [](float x) { return x * x - 2.0f; };
        const auto dff = [](float x) { return 2.0f * x; };
        const auto fl  = [](long double x) { return x * x - 2.0L; };

        const auto a = r::brent {}(ff, { 1.0f, 2.0f });
        const auto b = r::newton {}.with_derivative(dff)(ff, 1.0f);
        const auto c = nxx::first_of(r::secant {}.on(1.0f), r::bisection {}.on({ 1.0f, 2.0f }))(ff);
        const auto d =
            nxx::then(r::expand {}.on(nxx::bracket { 1.0f, 1.2f }), r::brent {}, r::newton {}.with_derivative(nxx::deriv::numeric {}))(ff);
        const auto e = r::bisection {}(fl, { 1.0L, 2.0L });
        const auto g = r::newton {}.with_derivative(nxx::deriv::derivative_of(fl))(fl, 1.0L);
        const auto h = nxx::deriv::diff(fl, 1.0L, nxx::deriv::central_1_4);

        const float sf =
            nxx::best_x(a).value_or(0.0f) + nxx::best_x(b).value_or(0.0f) + nxx::best_x(c).value_or(0.0f) + nxx::best_x(d).value_or(0.0f);
        const long double sl = nxx::best_x(e).value_or(0.0L) + nxx::best_x(g).value_or(0.0L) + h.value_or(0.0L);
        return static_cast<double>(sf) + static_cast<double>(sl);
    }

#if __has_include(<fxt/monads/Expected.hpp>)
    // The FXT pipes over results (DESIGN §8.1).
    double pipes()
    {
        using nxx::operator|;
        const auto   x_of  = [](const auto& s) { return s.x; };
        const auto   chain = nxx::first_of(r::newton {}.with_derivative(df).on(0.0), r::brent {}.on({ 1.0, 2.0 }));
        int          seen  = 0;
        const double a     = chain(f) | fxt::tap([&seen](const auto&) { ++seen; }) | fxt::transform(x_of) | fxt::value_or(0.0);
        const double b     = r::expand {}(f, { 1.0, 1.2 }) | fxt::and_then([](const auto& sb) { return r::brent {}(f, sb); }) |
                             fxt::transform(x_of) | fxt::value_or(0.0);
        const double c     = r::bisection {}.with_budget(3)(f, { 1.0, 2.0 }) |
                             fxt::or_else([](const auto&) { return r::brent {}(f, { 1.0, 2.0 }); }) |
                             fxt::match(x_of, [](const auto&) { return 0.0; });
        const auto   d     = r::secant {}(f, 1.0) | fxt::transform_error([](const auto& e) { return e.code; });
        return a + b + c + (d ? d->x : 0.0) + seen;
    }
#else
    double pipes() { return 0.0; }
#endif
}    // namespace consumer

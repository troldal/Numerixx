// A quick tour of Numerixx 2 as it stands after the de-risking spike (DESIGN §10.2): 1-D root finding, numerical
// derivatives, and composing solvers. Numerical failures are values, not exceptions: every result is a std::expected,
// either a solution (x, fx, how it stopped, what it cost) or a failure (an error code, what it cost, the best estimate
// so far and, for a fallible callback, the callback's own error). Numerixx itself does not throw, and it is
// exception-neutral: an exception thrown by your callback propagates out of the solver (DESIGN D10). The opt-in run-time
// chain of section 6 also allocates through std::function and std::vector, so wrapping or copying a solver there can
// throw std::bad_alloc, or whatever the callable's copy constructor throws.
//
// std::printf rather than std::println: with MinGW's libstdc++, std::print needs -lstdc++exp at link time.
#include <numerixx/core/any_solver.hpp>    // opt-in: run-time solver chains (the only header that uses std::function)
#include <numerixx/deriv.hpp>
#include <numerixx/roots.hpp>

#include <cmath>
#include <cstdio>
#include <expected>
#include <functional>
#include <limits>
#include <ranges>
#include <vector>

namespace r = nxx::roots;
namespace d = nxx::deriv;

namespace
{
    const char* name(nxx::errc e)
    {
        switch (e) {
            case nxx::errc::invalid_input:
                return "invalid_input";
            case nxx::errc::no_sign_change:
                return "no_sign_change";
            case nxx::errc::budget_exhausted:
                return "budget_exhausted";
            case nxx::errc::zero_derivative:
                return "zero_derivative";
            case nxx::errc::sign_change_not_root:
                return "sign_change_not_root";
            case nxx::errc::callback_failed:
                return "callback_failed";
            default:
                return "another error";
        }
    }

    const char* name(nxx::stop_reason s)
    {
        switch (s) {
            case nxx::stop_reason::exact_zero:
                return "exact zero";
            case nxx::stop_reason::criterion:
                return "criterion met";
            case nxx::stop_reason::resolution_limit:
                return "resolution limit";
            default:
                return "algorithm";
        }
    }

    // Prints a root-finding result: the solution, or the failure and its best estimate.
    template<class R>
    void report(const char* what, const R& res)
    {
        if (res) {
            std::printf("  %-38s x = %.17g, f(x) = %.2g (%s; %u iterations, %u evaluations)\n",
                        what,
                        res->x,
                        res->fx,
                        name(res->how),
                        res->used.iterations,
                        res->used.evaluations);
            return;
        }
        const auto& err = res.error();
        std::printf("  %-38s failed: %s after %u evaluation%s",
                    what,
                    name(err.code),
                    err.used.evaluations,
                    err.used.evaluations == 1 ? "" : "s");
        if (err.best) std::printf("; best x = %.17g, f(x) = %.2g", err.best->x, err.best->fx);
        std::printf("\n");
    }

    constexpr auto f  = [](double x) { return x * x - 2.0; };    // the root is sqrt(2)
    constexpr auto df = [](double x) { return 2.0 * x; };
}    // namespace

int main()
{
    // Run-time values, as they would come from a file or a user.
    double       lo = 1.0, hi = 2.0;
    double       user_tolerance = 1e-10;
    long long    user_budget    = 40;
    const double nan            = std::numeric_limits<double>::quiet_NaN();

    std::printf("1. One call\n");
    // The bracketing default (Brent, provisionally). {lo, hi} is validated in-band: the ends may come in either order,
    // and equal or non-finite ends give an errc::invalid_input failure, not an exception.
    report("solve(f, {lo, hi})", r::solve(f, { lo, hi }));
    report("solve(f, {hi, lo})", r::solve(f, { hi, lo }));
    report("solve(f, {hi, hi})", r::solve(f, { hi, hi }));

    std::printf("\n2. Choosing a solver and its stop criterion\n");
    // Bracketing solvers converge on their enclosure, so they take a width criterion. x_tol, which compares successive
    // iterates, does not compile with them, and the error says why.
    report("bisection{width_tol{1e-6}}", r::bisection { nxx::width_tol { 1e-6 } }(f, { lo, hi }));
    const auto tight = r::brent { nxx::width_tol { 1e-12 } }(f, { lo, hi });
    report("brent{width_tol{1e-12}}", tight);
    if (tight && tight->enclosure)
        std::printf("  %-38s the root is in [%.17g, %.17g]\n", "", tight->enclosure->lo(), tight->enclosure->hi());

    // Literals are checked at compile time (width_tol{-1.0} does not compile). Run-time values go through make(),
    // which returns a std::expected.
    const auto tol    = nxx::width_tol<double>::make(user_tolerance, 0.0);
    const auto budget = nxx::max_iterations::make(user_budget);
    if (tol && budget) report("brent{run-time tol}.with_budget(40)", r::brent { *tol }.with_budget(*budget)(f, { lo, hi }));

    std::printf("\n3. Open methods, from a guess\n");
    report("newton with f'", r::newton {}.with_derivative(df)(f, 1.0));
    report("newton with a numeric f'", r::newton {}.with_derivative(d::numeric {})(f, 1.0));    // counts f's calls too
    report("secant (derivative-free)", r::secant {}(f, 1.0));
    report("secant kept inside [0, 10]", r::secant {}.with_projection(r::clamp_to { 0.0, 10.0 })(f, 5.0));
    // Open methods stop on successive iterates: x_tol{abs[, rel]}, or the default step_tol.
    report("secant{x_tol{1e-12}}", r::secant { nxx::x_tol { 1e-12 } }(f, 1.0));

    std::printf("\n4. Failures are values\n");
    report("brent on [2, 3] (no sign change)", r::brent {}(f, { 2.0, 3.0 }));
    report("bisection on tan over [1, 2] (a pole)", r::bisection {}([](double x) { return std::tan(x); }, { 1.0, 2.0 }));
    const auto starved = r::newton {}.with_derivative(df).with_budget(3)(f, 100.0);
    report("newton with a budget of 3", starved);
    // best_x takes the solution's x or the failure's best x. A failure's best estimate is only a candidate (here it is
    // far from the root), so check it, for example its f(x) or uncertainty, before accepting it.
    std::printf("  %-38s best_x = %.17g\n", "", nxx::best_x(starved).value_or(nan));

    std::printf("\n5. A fallible callback keeps its own error\n");
    enum class model_error : unsigned char { out_of_range };
    const auto model = [](double x) -> std::expected<double, model_error> {
        if (x < 0.0) return std::unexpected(model_error::out_of_range);
        return std::sqrt(x) - 1.2;
    };
    const auto bad = r::bisection {}(model, { -1.0, 4.0 });
    report("bisection over [-1, 4]", bad);
    if (!bad && bad.error().cause && *bad.error().cause == model_error::out_of_range)
        std::printf("  %-38s the cause is model_error::out_of_range\n", "");
    report("bisection over [0, 4]", r::bisection {}(model, { 0.0, 4.0 }));

    std::printf("\n6. Composing solvers\n");
    // first_of: if one solver fails, the next one tries. The function is supplied later (.on binds the input), so the
    // chain is a value that can be stored, copied and run at compile time. The cost includes the failed attempts.
    constexpr auto chain = nxx::first_of(r::newton {}.with_derivative(df).on(0.0),          // f'(0) = 0: zero_derivative
                                         r::secant {}.with_budget(5).on(0.0),               // needs more than 5 steps
                                         r::bisection {}.on(nxx::bracket { 0.0, 2.0 }));    // succeeds
    static_assert(chain(f).has_value());                                                    // the whole chain runs at compile time too
    report("first_of(newton, secant, bisection)", chain(f));

    // then: each stage starts from the previous stage's result. Grow a window until f changes sign, narrow the
    // enclosure coarsely, then polish with Newton.
    constexpr auto staged =
        nxx::then(r::expand {}.on(nxx::bracket { 2.0, 2.5 }), r::bisection { nxx::width_tol { 1e-4 } }, r::newton {}.with_derivative(df));
    report("then(expand, bisection, newton)", staged(f));

    // A chain assembled at run time, from solvers of different types: any_solver erases a solver's type, for one type F
    // of the function passed to the chain (here std::function<double(double)>). Only that function must be an F; the
    // solvers inside keep their own callables, such as Newton's derivative.
    using fn_t      = std::function<double(double)>;
    using solver_t  = nxx::any_solver<fn_t, r::root_estimate<double>>;
    auto candidates = std::vector<solver_t> { r::newton {}.with_derivative(df).on(0.0), r::brent {}.on({ lo, hi }) };
    report("first_of(run-time vector)", nxx::first_of(candidates)(fn_t { f }));

    std::printf("\n7. Derivatives\n");
    const auto sine = [](double x) { return std::sin(x); };
    if (const auto slope = d::central(sine, 1.0))
        std::printf("  %-38s %.17g (error %.1e)\n", "central(sin, 1)", *slope, *slope - std::cos(1.0));
    if (const auto slope = d::diff(sine, 1.0, d::central_1_4))
        std::printf("  %-38s %.17g (error %.1e)\n", "diff(sin, 1, central_1_4)", *slope, *slope - std::cos(1.0));
    const auto cosine = d::derivative_of(sine);    // a callable x -> expected<double, fault>, with its own cost
    if (const auto slope = cosine(0.5))
        std::printf("  %-38s %.17g (error %.1e)\n", "derivative_of(sin)(0.5)", *slope, *slope - std::cos(0.5));

    std::printf("\n8. Stepping by hand\n");
    // steps_view yields the solver's states one by one; it applies neither the stop criterion nor the budget.
    const auto solver = r::brent {};
    if (const auto p = solver.prepare(std::cref(f), nxx::bracket { 1.0, 2.0 })) {
        int k = 0;
        for (const auto& st : nxx::steps_view { solver, *p } | std::views::take(8)) {
            if (!st) break;
            const auto e = solver.estimate(*st);
            std::printf("  step %d: x = %.17g, uncertainty %.2g\n", k++, e.x, e.uncertainty);
        }
    }
    // An observer sees every iterate of a normal solve.
    int seen = 0;
    report("bisection with an observer", r::bisection {}.with_observer([&seen](const auto&) { ++seen; })(f, { lo, hi }));
    std::printf("  %-38s the observer saw %d iterates\n", "", seen);
    return 0;
}

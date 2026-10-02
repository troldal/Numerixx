// FXT pipes over Numerixx results (DESIGN §8.1) and spike exit criterion 3 (DESIGN §10.2): first_of and then over
// solvers with different state types (newton, secant, brent, bisection, expand) and fallible callbacks, including
// .with_projection/.with_observer values, called through the family facades and curried with .on(), and consumed
// with the FXT pipes of <numerixx/pipes.hpp> (DESIGN §6.10, §6.11).
#include <numerixx/deriv.hpp>
#include <numerixx/pipes.hpp>
#include <numerixx/roots.hpp>

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <expected>
#include <functional>
#include <limits>
#include <optional>
#include <type_traits>
#include <utility>

namespace
{
    enum class demo_error { failed };

    constexpr auto halve(int value) -> std::expected<int, demo_error>
    {
        if (value % 2 != 0) return std::unexpected(demo_error::failed);
        return value / 2;
    }

    namespace r = nxx::roots;

    // The user's own callback error. is_fatal, found by ADL, makes eval_error::fatal stop first_of's fall-through.
    enum class eval_error { domain, fatal };

    [[maybe_unused]] constexpr bool is_fatal(eval_error e) noexcept { return e == eval_error::fatal; }

    // x^2 - 2 on its domain x >= 0 (root sqrt 2); below it the callback fails with the user's error.
    constexpr auto parabola = [](double x) -> std::expected<double, eval_error> {
        if (x < 0.0) return std::unexpected(eval_error::domain);
        return x * x - 2.0;
    };
    // A plain derivative: Newton's cause is still eval_error, the common cause of f and f' (DESIGN §6.4).
    constexpr auto parabola_slope = [](double x) { return 2.0 * x; };

    // The same parabola, but a negative argument is fatal.
    constexpr auto parabola_fatal = [](double x) -> std::expected<double, eval_error> {
        if (x < 0.0) return std::unexpected(eval_error::fatal);
        return x * x - 2.0;
    };

    // log(x) - 1/2 on x > 0 (root e^(1/2)). Newton from 6 proposes 6 - 6 (log 6 - 1/2) = -1.75, outside the domain.
    constexpr auto logarithm = [](double x) -> std::expected<double, eval_error> {
        if (!(x > 0.0)) return std::unexpected(eval_error::domain);
        return std::log(x) - 0.5;
    };
    constexpr auto logarithm_slope = [](double x) { return 1.0 / x; };

    // Plain callbacks: the cause is nxx::none.
    constexpr auto plain        = [](double x) { return x * x - 2.0; };
    constexpr auto plain_slope  = [](double x) { return 2.0 * x; };
    constexpr auto no_real_root = [](double x) { return x * x + 1.0; };

    constexpr double sqrt2        = 1.4142135623730951;
    constexpr double sqrt_e       = 1.6487212707001282;
    constexpr double not_a_number = std::numeric_limits<double>::quiet_NaN();

    bool close_to(double a, double b, double tol) { return std::abs(a - b) <= tol; }

    constexpr auto x_of     = [](const auto& sol) { return sol.x; };
    constexpr auto by_of    = [](const auto& sol) { return sol.by; };
    constexpr auto used_of  = [](const auto& sol) { return sol.used; };
    constexpr auto code_of  = [](const auto& err) { return err.code; };
    constexpr auto where_of = [](const auto& err) { return err.where; };
    constexpr auto cause_of = [](const auto& err) { return err.cause; };
    constexpr auto spent_of = [](const auto& err) { return err.used; };

    using estimate   = r::root_estimate<double>;
    using fallible   = nxx::result<estimate, eval_error>;
    using infallible = nxx::result<estimate>;

    // A result's x, or NaN on failure; its cost, or zero on failure.
    template<class R>
    double x_or_nan(const R& res)
    {
        using nxx::operator|;
        return res | fxt::transform(x_of) | fxt::value_or(not_a_number);
    }

    template<class R>
    nxx::counters used_or_zero(const R& res)
    {
        using nxx::operator|;
        return res | fxt::transform(used_of) | fxt::value_or(nxx::counters {});
    }

    // An observer's record: solvers are immutable values, so the state lives outside, captured by reference.
    struct trace
    {
        int    calls = 0;
        double lo    = std::numeric_limits<double>::infinity();
        double hi    = -std::numeric_limits<double>::infinity();
    };

    auto recorder(trace& t)
    {
        return [&t](const auto& view) {
            ++t.calls;
            t.lo = (std::min)(t.lo, view.x());
            t.hi = (std::max)(t.hi, view.x());
        };
    }

    // Two results of the same computation: same success or failure, same x, f(x), cost and algorithm.
    template<class R1, class R2>
    bool same_outcome(const R1& a, const R2& b)
    {
        if (a.has_value() != b.has_value()) return false;
        if (a) return a->x == b->x && a->fx == b->fx && a->used == b->used && a->by == b->by && a->how == b->how;
        return a.error().code == b.error().code && a.error().where == b.error().where && a.error().used == b.error().used &&
               nxx::best_x(a) == nxx::best_x(b);
    }
}    // namespace

TEST_SUITE("pipes")
{
    TEST_CASE("FXT pipes compose std::expected values")
    {
        using nxx::operator|;

        const std::expected<int, demo_error> start { 8 };
        const auto quarter = start | fxt::and_then(halve) | fxt::and_then(halve) | fxt::transform([](int v) { return v + 1; });
        CHECK(quarter == 3);

        const auto odd = std::expected<int, demo_error> { 3 } | fxt::and_then(halve);
        CHECK_FALSE(odd.has_value());
        CHECK((odd | fxt::value_or(-1)) == -1);
    }

    TEST_CASE("first_of over newton secant and brent with a fallible callback consumed through pipes")
    {
        using nxx::operator|;

        const auto chain = nxx::first_of(r::newton {}.with_derivative(parabola_slope).on(0.0),    // f'(0) = 0: zero_derivative
                                         r::secant {}.on(-1.0),                                   // f(-1) fails: domain
                                         r::brent {}.on({ 0.0, 2.0 }));                           // succeeds
        const auto res   = chain(parabola);
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(res)>, fallible>);

        CHECK((res | fxt::transform(by_of) | fxt::value_or(nxx::algo::none)) == r::algos::brent);
        CHECK(close_to(res | fxt::transform(x_of) | fxt::value_or(not_a_number), sqrt2, 1e-14));

        // The success pays for the failed attempts: newton 1 iteration and 2 evaluations (f(0), f'(0)), secant 1
        // evaluation (the failed f(-1)).
        const auto alone = r::brent {}(parabola, { 0.0, 2.0 });
        const auto spent = alone | fxt::transform(used_of) | fxt::value_or(nxx::counters {});
        CHECK((res | fxt::transform(used_of)) == nxx::counters { spent.iterations + 1, spent.evaluations + 3 });
        CHECK(x_or_nan(res) == x_or_nan(alone));

        // Each alternative's own failure, through the pipes.
        const auto newton_alone = r::newton {}.with_derivative(parabola_slope)(parabola, 0.0);
        CHECK((newton_alone | fxt::transform_error(code_of)) == std::unexpected(nxx::errc::zero_derivative));
        CHECK((newton_alone | fxt::transform_error(spent_of)) == std::unexpected(nxx::counters { 1, 2 }));
        const auto secant_alone = r::secant {}(parabola, -1.0);
        CHECK((secant_alone | fxt::transform_error(code_of)) == std::unexpected(nxx::errc::callback_failed));
        CHECK((secant_alone | fxt::transform_error(cause_of)) == std::unexpected(std::optional { eval_error::domain }));

        // tap sees the solution and passes the result on unchanged.
        int        tapped = 0;
        const auto again  = chain(parabola) | fxt::tap([&tapped](const auto& sol) { tapped += sol.by == r::algos::brent ? 1 : 0; });
        CHECK(tapped == 1);
        CHECK(same_outcome(again, res));
    }

    TEST_CASE("then stages expand bisection and newton over a fallible callback and the pipes stage the same solves")
    {
        using nxx::operator|;

        const auto pipeline = nxx::then(r::expand {}.on(nxx::bracket { 2.0, 2.5 }),       // grows to a sign change
                                        r::bisection { nxx::width_tol { 1e-4 } },         // coarse enclosure
                                        r::newton {}.with_derivative(parabola_slope));    // polish, seeded
        const auto res      = pipeline(parabola);
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(res)>, fallible>);

        CHECK((res | fxt::transform(by_of) | fxt::value_or(nxx::algo::none)) == r::algos::newton);
        CHECK(close_to(res | fxt::transform(x_of) | fxt::value_or(not_a_number), sqrt2, 1e-14));

        // The same stages by hand: the family facades' operator(), chained with fxt::and_then. then charges every
        // stage's cost to the result; and_then does not, so tap adds them up.
        nxx::counters spent {};
        const auto    charge = [&spent](const auto& sol) { spent = spent + sol.used; };
        const auto    staged =
            r::expand {}(parabola, { 2.0, 2.5 }) | fxt::tap(charge) |
            fxt::and_then([](const auto& sb) { return r::bisection { nxx::width_tol { 1e-4 } }(parabola, sb); }) | fxt::tap(charge) |
            fxt::and_then([](const auto& coarse) { return r::newton {}.with_derivative(parabola_slope)(parabola, coarse); }) |
            fxt::tap(charge);
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(staged)>, fallible>);
        CHECK(x_or_nan(staged) == x_or_nan(res));    // same values of f, same path
        CHECK((res | fxt::transform(used_of)) == spent);

        // A failing first stage passes through then unchanged: the user's cause, and where it failed.
        const auto blocked = nxx::then(r::expand {}.on(nxx::bracket { -3.0, -2.0 }), r::brent {})(parabola);
        CHECK((blocked | fxt::transform_error(code_of)) == std::unexpected(nxx::errc::callback_failed));
        CHECK((blocked | fxt::transform_error(where_of)) == std::unexpected(r::algos::expand));
        CHECK((blocked | fxt::transform_error(cause_of)) == std::unexpected(std::optional { eval_error::domain }));
        CHECK((blocked | fxt::transform_error(spent_of)) == std::unexpected(nxx::counters { 0, 1 }));
    }

    TEST_CASE("first_of over then chains and bracketing solvers keeps the cause and a fatal cause stops the chain")
    {
        using nxx::operator|;

        // An alternative that is a then chain, then one that is a curried bisection.
        const auto chain = nxx::first_of(nxx::then(r::expand {}.on(nxx::bracket { -3.0, -2.0 }), r::brent {}),    // domain
                                         nxx::then(r::expand {}.on(nxx::bracket { 2.0, 2.5 }), r::brent {}),      // succeeds
                                         r::bisection {}.on({ 0.0, 2.0 }));
        const auto res   = chain(parabola);
        CHECK((res | fxt::transform(by_of) | fxt::value_or(nxx::algo::none)) == r::algos::brent);
        CHECK(close_to(res | fxt::transform(x_of) | fxt::value_or(not_a_number), sqrt2, 1e-14));
        const auto second = nxx::then(r::expand {}.on(nxx::bracket { 2.0, 2.5 }), r::brent {})(parabola);
        CHECK(second.has_value());
        CHECK((res | fxt::transform(used_of)) == used_or_zero(second) + nxx::counters { 0, 1 });

        // Every alternative fails with the user's error: the last code and cause, where the last one failed, the
        // total cost, and no best estimate (nothing was evaluated successfully).
        const auto all_fail = nxx::first_of(r::bisection {}.on({ -2.0, -1.0 }), r::newton {}.with_derivative(parabola_slope).on(-1.0));
        const auto failed   = all_fail(parabola);
        CHECK((failed | fxt::transform_error(code_of)) == std::unexpected(nxx::errc::callback_failed));
        CHECK((failed | fxt::transform_error(where_of)) == std::unexpected(r::algos::newton));
        CHECK((failed | fxt::transform_error(cause_of)) == std::unexpected(std::optional { eval_error::domain }));
        CHECK((failed | fxt::transform_error(spent_of)) == std::unexpected(nxx::counters { 0, 2 }));
        CHECK_FALSE(nxx::best_x(failed).has_value());

        // match folds both channels into one value.
        const auto verdict = failed | fxt::match([](const auto&) { return std::optional<eval_error> {}; }, cause_of);
        CHECK(verdict == std::optional { eval_error::domain });

        // or_else recovers with another solver of the same value type.
        const auto recovered = failed | fxt::or_else([](const auto&) { return r::brent {}(parabola, { 0.0, 2.0 }); });
        CHECK(close_to(recovered | fxt::transform(x_of) | fxt::value_or(not_a_number), sqrt2, 1e-14));

        // A fatal cause stops the fall-through: brent never runs (its observer is never called) and f is called once.
        int           later   = 0;
        const auto    guarded = nxx::first_of(r::bisection {}.on({ -1.0, 2.0 }),
                                              r::brent {}.with_observer([&later](const auto&) { ++later; }).on({ 0.0, 2.0 }));
        std::uint32_t calls   = 0;
        const auto    stopped = guarded(nxx::fn::counted(parabola_fatal, calls));
        CHECK(later == 0);
        CHECK(calls == 1);
        CHECK((stopped | fxt::transform_error(where_of)) == std::unexpected(r::algos::bisection));
        CHECK((stopped | fxt::transform_error(cause_of)) == std::unexpected(std::optional { eval_error::fatal }));

        // The same chain with a non-fatal cause falls through to brent.
        const auto fell_through = guarded(parabola);
        CHECK(later > 0);
        CHECK((fell_through | fxt::transform(by_of) | fxt::value_or(nxx::algo::none)) == r::algos::brent);
        CHECK(close_to(fell_through | fxt::transform(x_of) | fxt::value_or(not_a_number), sqrt2, 1e-14));
    }

    TEST_CASE("with_projection and with_observer values inside first_of and then")
    {
        using nxx::operator|;

        // Without a projection, Newton from 6 leaves the domain at its first step.
        const auto free_newton = r::newton {}.with_derivative(logarithm_slope);
        const auto escaped     = free_newton(logarithm, 6.0);
        CHECK((escaped | fxt::transform_error(code_of)) == std::unexpected(nxx::errc::callback_failed));
        CHECK((escaped | fxt::transform_error(cause_of)) == std::unexpected(std::optional { eval_error::domain }));
        CHECK(nxx::best_x(escaped) == std::optional { 6.0 });

        // With a clamp, the same start converges; the observer sees every iterate, all inside the box.
        trace      seen;
        const auto clamped =
            r::newton {}.with_derivative(logarithm_slope).with_projection(r::clamp_to { 0.5, 10.0 }).with_observer(recorder(seen));
        const auto chain = nxx::first_of(free_newton.on(6.0), clamped.on(6.0));
        const auto res   = chain(logarithm);
        CHECK((res | fxt::transform(by_of) | fxt::value_or(nxx::algo::none)) == r::algos::newton);
        CHECK(close_to(res | fxt::transform(x_of) | fxt::value_or(not_a_number), sqrt_e, 1e-14));
        CHECK(seen.calls > 0);
        CHECK(seen.lo >= 0.5);
        CHECK(seen.hi <= 10.0);
        CHECK(seen.lo == 0.5);    // the first step was clamped
        // One observed call per iteration of the clamped solve; the chain adds the free Newton's one failed step.
        CHECK((res | fxt::transform([](const auto& sol) { return sol.used.iterations; })) == static_cast<std::uint32_t>(seen.calls) + 1u);
        const int chain_calls = seen.calls;

        // Builders are order-independent: the same configuration in another order gives the same result.
        trace      seen2;
        const auto reordered =
            r::newton {}.with_observer(recorder(seen2)).with_projection(r::clamp_to { 0.5, 10.0 }).with_derivative(logarithm_slope);
        CHECK(same_outcome(reordered(logarithm, 6.0), clamped(logarithm, 6.0)));
        CHECK(seen2.calls == chain_calls);
        CHECK(seen.calls == 2 * chain_calls);    // the clamped solver, called once more, reports to the same record

        // A clamp that excludes the root pins the iterate at the edge: stalled, with the best estimate at the edge.
        const auto pinned = r::newton {}.with_derivative(logarithm_slope).with_projection(r::clamp_to { 2.0, 10.0 }).on(6.0)(logarithm);
        CHECK((pinned | fxt::transform_error(code_of)) == std::unexpected(nxx::errc::stalled));
        CHECK(nxx::best_x(pinned) == std::optional { 2.0 });
        // or_else accepts the best estimate of a stalled solve as a value.
        const auto accepted = pinned | fxt::or_else([](const auto& err) -> fallible {
                                  if (err.code == nxx::errc::stalled && err.best)
                                      return nxx::solution<estimate> { *err.best, err.used, err.where, nxx::stop_reason::algorithm };
                                  return std::unexpected(err);
                              });
        CHECK((accepted | fxt::transform(x_of) | fxt::value_or(not_a_number)) == 2.0);

        // A secant with a projection and an observer as the second stage of then, with run-time bounds.
        const double lo = 1.0;
        const double hi = 3.0;
        trace        polish;
        const auto   staged = nxx::then(r::bisection { nxx::width_tol { 1e-2 } }.on({ lo, hi }),
                                        r::secant {}.with_projection(r::clamp_to { lo, hi }).with_observer(recorder(polish)))(logarithm);
        CHECK((staged | fxt::transform(by_of) | fxt::value_or(nxx::algo::none)) == r::algos::secant);
        CHECK(close_to(staged | fxt::transform(x_of) | fxt::value_or(not_a_number), sqrt_e, 1e-14));
        CHECK(polish.calls > 0);
        CHECK(polish.lo >= lo);
        CHECK(polish.hi <= hi);

        // An observer on a bracketing solver and on a searcher: one call per iteration.
        trace      enclosing;
        trace      searching;
        const auto searched       = nxx::then(r::expand {}.with_observer(recorder(searching)).on(nxx::bracket { 2.0, 2.5 }),
                                              r::brent {}.with_observer(recorder(enclosing)))(parabola);
        const auto searched_alone = r::expand {}(parabola, { 2.0, 2.5 });
        CHECK((searched_alone | fxt::transform([](const auto& sol) { return sol.used.iterations; })) ==
              static_cast<std::uint32_t>(searching.calls));
        CHECK((searched | fxt::transform([](const auto& sol) { return sol.used.iterations; })) ==
              static_cast<std::uint32_t>(searching.calls + enclosing.calls));
        CHECK(close_to(searched | fxt::transform(x_of) | fxt::value_or(not_a_number), sqrt2, 1e-14));
    }

    TEST_CASE("direct calls through the family facades equal the curried calls")
    {
        using nxx::operator|;

        const double lo = 0.0;
        const double hi = 2.0;

        // Bracketing facade: braced run-time bracket, std::pair, bracket<T>::make, curried.
        const auto braced = r::brent {}(parabola, { lo, hi });
        CHECK(same_outcome(braced, r::brent {}.on({ lo, hi })(parabola)));
        CHECK(same_outcome(braced, r::brent {}(parabola, std::pair { lo, hi })));
        CHECK(same_outcome(braced, r::brent {}(parabola, nxx::bracket<double>::make(hi, lo))));    // re-ordered
        CHECK(same_outcome(r::bisection {}(parabola, { lo, hi }), r::bisection {}.on({ lo, hi })(parabola)));

        // Open facade: a guess, curried or not.
        CHECK(same_outcome(r::secant {}(parabola, 1.0), r::secant {}.on(1.0)(parabola)));
        CHECK(same_outcome(r::newton {}.with_derivative(parabola_slope)(parabola, 1.0),
                           r::newton {}.with_derivative(parabola_slope).on(1.0)(parabola)));

        // Search facade: a window, curried or not.
        const auto window  = r::expand {}(parabola, { 2.0, 2.5 });
        const auto curried = r::expand {}.on({ 2.0, 2.5 })(parabola);
        CHECK(window.has_value());
        const auto lo_of = [](const auto& sb) { return sb.lo(); };
        CHECK((window | fxt::transform(lo_of) | fxt::value_or(not_a_number)) ==
              (curried | fxt::transform(lo_of) | fxt::value_or(not_a_number)));
        CHECK(used_or_zero(window) == used_or_zero(curried));
        CHECK((window | fxt::transform([](const auto& sb) { return sb.lo() < sqrt2 && sqrt2 < sb.hi(); }) | fxt::value_or(false)));

        // first_of over uncurried solvers, called with (f, input): bracketing with a bracket, open methods with a guess.
        const auto bracketing = nxx::first_of(r::bisection {}.with_budget(5), r::brent {})(parabola, std::pair { lo, hi });
        CHECK((bracketing | fxt::transform(by_of) | fxt::value_or(nxx::algo::none)) == r::algos::brent);
        CHECK(close_to(bracketing | fxt::transform(x_of) | fxt::value_or(not_a_number), sqrt2, 1e-14));
        const auto starved = r::bisection {}.with_budget(5)(parabola, std::pair { lo, hi });
        CHECK((starved | fxt::transform_error(code_of)) == std::unexpected(nxx::errc::budget_exhausted));
        CHECK((starved | fxt::transform_error(spent_of)) == std::unexpected(nxx::counters { 5, 7 }));
        CHECK((bracketing | fxt::transform(used_of)) == used_or_zero(braced) + nxx::counters { 5, 7 });

        const auto open = nxx::first_of(r::newton {}.with_derivative(parabola_slope).with_budget(1), r::secant {})(parabola, 1.0);
        CHECK((open | fxt::transform(by_of) | fxt::value_or(nxx::algo::none)) == r::algos::secant);
        CHECK(close_to(open | fxt::transform(x_of) | fxt::value_or(not_a_number), sqrt2, 1e-14));

        // Fix f, vary the input (DESIGN §6.14).
        const auto solve_at = std::bind_front(r::brent {}, parabola);
        CHECK(same_outcome(solve_at(std::pair { lo, hi }), braced));
    }

    TEST_CASE("derivative_of with a fallible callback and the numeric policy inside chains")
    {
        using nxx::operator|;

        // derivative_of(parabola) at 0 samples parabola(-h): the derivative fails with the user's cause, unnested, and
        // the failing step costs the one evaluation it spent.
        const auto analytic_free = r::newton {}.with_derivative(nxx::deriv::derivative_of(parabola));
        const auto at_zero       = analytic_free(parabola, 0.0);
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(at_zero)>, fallible>);
        CHECK((at_zero | fxt::transform_error(code_of)) == std::unexpected(nxx::errc::callback_failed));
        CHECK((at_zero | fxt::transform_error(cause_of)) == std::unexpected(std::optional { eval_error::domain }));
        CHECK((at_zero | fxt::transform_error(spent_of)) == std::unexpected(nxx::counters { 1, 2 }));

        // A curried chain: the numeric policy binds f when the chain is called.
        const auto chain = nxx::first_of(analytic_free.on(0.0), r::newton {}.with_derivative(nxx::deriv::numeric {}).on(1.0));
        const auto res   = chain(parabola);
        CHECK((res | fxt::transform(by_of) | fxt::value_or(nxx::algo::none)) == r::algos::newton);
        CHECK(close_to(res | fxt::transform(x_of) | fxt::value_or(not_a_number), sqrt2, 1e-13));

        // Polish a coarse bisection with a derivative_of Newton through and_then.
        const auto polished = r::bisection { nxx::width_tol { 1e-3 } }(parabola, { 1.0, 2.0 }) | fxt::and_then([](const auto& coarse) {
                                  return r::newton {}.with_derivative(nxx::deriv::derivative_of(parabola))(parabola, coarse);
                              });
        CHECK((polished | fxt::transform(by_of) | fxt::value_or(nxx::algo::none)) == r::algos::newton);
        CHECK(close_to(polished | fxt::transform(x_of) | fxt::value_or(not_a_number), sqrt2, 1e-13));
    }

    TEST_CASE("a plain callback has cause nxx::none and the headline chain runs through pipes at compile time")
    {
        using nxx::operator|;

        constexpr auto headline = nxx::first_of(r::newton {}.with_derivative(plain_slope).on(0.0),    // zero_derivative
                                                r::secant {}.with_budget(5).on(0.0),                  // budget_exhausted
                                                r::bisection {}.on(nxx::bracket { 0.0, 2.0 }));       // succeeds
        static_assert((headline(plain) | fxt::transform(by_of) | fxt::value_or(nxx::algo::none)) == r::algos::bisection);

        const auto res = headline(plain);
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(res)>, infallible>);
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(res)>::error_type::cause_type, nxx::none>);
        CHECK(close_to(res | fxt::transform(x_of) | fxt::value_or(not_a_number), sqrt2, 1e-14));

        // The cost is that of the three attempts, each run alone.
        nxx::counters spent {};
        const auto    charge_failure = [&spent](const auto& err) { spent = spent + err.used; };
        const auto    charge_success = [&spent](const auto& sol) { spent = spent + sol.used; };
        const auto    first          = r::newton {}.with_derivative(plain_slope)(plain, 0.0) | fxt::tap_error(charge_failure);
        const auto    second         = r::secant {}.with_budget(5)(plain, 0.0) | fxt::tap_error(charge_failure);
        const auto    third          = r::bisection {}(plain, nxx::bracket { 0.0, 2.0 }) | fxt::tap(charge_success);
        CHECK((first | fxt::transform_error(code_of)) == std::unexpected(nxx::errc::zero_derivative));
        CHECK((second | fxt::transform_error(code_of)) == std::unexpected(nxx::errc::budget_exhausted));
        CHECK((res | fxt::transform(used_of)) == spent);
        CHECK(x_or_nan(res) == x_or_nan(third));

        // No real root: every alternative fails; the chain reports the last code and keeps a best estimate.
        const auto none_found = headline(no_real_root);
        CHECK((none_found | fxt::transform_error(code_of)) == std::unexpected(nxx::errc::no_sign_change));
        CHECK((none_found | fxt::transform_error(where_of)) == std::unexpected(r::algos::bisection));
        CHECK(nxx::best_x(none_found).has_value());
        CHECK((none_found | fxt::match([](const auto&) { return 0; },
                                       [](const auto& err) { return err.code == nxx::errc::no_sign_change ? 1 : 2; })) == 1);

        // Coarse bisection, then secant (canonical call 10), with a run-time bracket.
        const double lo     = 0.0;
        const double hi     = 2.0;
        const auto   staged = nxx::then(r::bisection { nxx::width_tol { 1e-3 } }.with_budget(100).on({ lo, hi }), r::secant {})(plain);
        CHECK((staged | fxt::transform(by_of) | fxt::value_or(nxx::algo::none)) == r::algos::secant);
        CHECK(close_to(staged | fxt::transform(x_of) | fxt::value_or(not_a_number), sqrt2, 1e-14));
    }
}

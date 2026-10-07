// Canonical calls (DESIGN §6.14, normative): calls 1, 2, 3, 9 and 10 (spike exit criterion 8, §10.2) and, from phase 1,
// calls 13, 14 and 15, spelled as in the table, with RUN-TIME brackets, guesses, tolerances and budgets. Run-time scalars come from
// volatile reads; tolerances and budgets go through make() and are dereferenced only after checking. Each module phase adds its calls here.
#include <numerixx/roots.hpp>

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <type_traits>
#include <utility>

namespace r = nxx::roots;

namespace
{
    // Not constant expressions: the calls below cannot validate these at compile time.
    double rt(double v)
    {
        volatile double x = v;
        return x;
    }

    long long rt_ll(long long v)
    {
        volatile long long n = v;
        return n;
    }

    const double     sqrt2 = std::sqrt(2.0);
    constexpr double qnan  = std::numeric_limits<double>::quiet_NaN();

    // A user's strong type (call 3).
    struct Length
    {
        double value;
    };

    // Two results are the same: same outcome, value, cost, algorithm and stop reason (or code and best estimate).
    template<class R>
    bool same_result(const R& a, const R& b)
    {
        if (a.has_value() != b.has_value()) return false;
        if (a)
            return a->x == b->x && a->fx == b->fx && a->uncertainty == b->uncertainty && a->enclosure == b->enclosure &&
                   a->used == b->used && a->by == b->by && a->how == b->how;
        const auto& ea = a.error();
        const auto& eb = b.error();
        return ea.code == eb.code && ea.by == eb.by && ea.used == eb.used && ea.best.has_value() == eb.best.has_value() &&
               (!ea.best || ea.best->x == eb.best->x) && ea.cause == eb.cause;
    }
}    // namespace

TEST_SUITE("usage")
{
    TEST_CASE("call 1: a root in a run-time bracket")
    {
        const auto   f  = [](double x) { return x * x - 2.0; };
        const double lo = rt(1.0);
        const double hi = rt(2.0);

        const auto a = nxx::roots::solve(f, { lo, hi });
        const auto b = r::brent {}(f, { lo, hi });
        static_assert(std::is_same_v<decltype(a), decltype(b)>);
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(a)>, nxx::result<r::root_estimate<double>>>);
        CHECK(a.has_value());
        if (a) {
            CHECK(std::abs(a->x - sqrt2) < 1e-14);
            CHECK(a->by == r::algos::brent);
            CHECK(a->used.evaluations > 2);
        }
        CHECK(same_result(a, b));    // solve(f, {lo, hi}) is brent: a composition, not a second implementation

        // Run-time tolerance and budget, through make().
        const auto tol    = nxx::tolerance<double>::make(rt(1e-10));
        const auto budget = nxx::max_iterations::make(rt_ll(60));
        CHECK(tol.has_value());
        CHECK(budget.has_value());
        if (tol && budget) {
            const auto c = r::brent { nxx::width_tol { *tol } }.with_budget(*budget)(f, { lo, hi });
            CHECK(c.has_value());
            if (c) {
                CHECK(std::abs(c->x - sqrt2) <= 1e-10);
                CHECK(c->uncertainty <= 1e-10 + 8.0 * std::numeric_limits<double>::epsilon() * std::abs(c->x));
            }
            // The same through a run-time bracket from make() and .on().
            const auto d = r::brent { nxx::width_tol { *tol } }.with_budget(*budget).on(nxx::bracket<double>::make(lo, hi))(f);
            CHECK(same_result(c, d));
        }

        // A run-time bracket that is not one fails in-band, at zero cost and with no best estimate.
        const auto same_ends = nxx::roots::solve(f, { lo, lo });
        CHECK_FALSE(same_ends.has_value());
        if (!same_ends) {
            CHECK(same_ends.error().code == nxx::errc::invalid_input);
            CHECK(same_ends.error().used == nxx::counters {});
            CHECK_FALSE(same_ends.error().best.has_value());
        }
        const auto reversed = nxx::roots::solve(f, { hi, lo });    // re-ordered, as bracket<T>::make does
        CHECK(same_result(reversed, a));
    }

    TEST_CASE("call 2: Newton with an analytic derivative")
    {
        const auto   f  = [](double x) { return x * x - 2.0; };
        const auto   df = [](double x) { return 2.0 * x; };
        const double x0 = rt(1.0);

        const auto res = r::newton {}.with_derivative(df)(f, x0);
        CHECK(res.has_value());
        if (res) {
            CHECK(std::abs(res->x - sqrt2) < 1e-15);
            CHECK(res->by == r::algos::newton);
        }

        const auto neg = r::newton {}.with_derivative(df)(f, rt(-1.0));
        CHECK(neg.has_value());
        if (neg) { CHECK(std::abs(neg->x + sqrt2) < 1e-15); }

        // Run-time budget and tolerance, through make().
        const auto budget = nxx::max_iterations::make(rt_ll(30));
        const auto tol    = nxx::tolerance<double>::make(rt(1e-12));
        CHECK(budget.has_value());
        CHECK(tol.has_value());
        if (budget && tol) {
            const auto same_budget = r::newton {}.with_derivative(df).with_budget(*budget)(f, x0);
            CHECK(same_result(same_budget, res));    // 30 is the default

            const auto by_x_tol = r::newton { nxx::x_tol { *tol } }.with_derivative(df).with_budget(*budget)(f, x0);
            CHECK(by_x_tol.has_value());
            if (by_x_tol) {
                CHECK(std::abs(by_x_tol->x - sqrt2) < 1e-12);
                CHECK(by_x_tol->how == nxx::stop_reason::criterion);
            }
        }

        // A run-time budget that is too small: budget_exhausted with the best estimate, never success.
        const auto two = nxx::max_iterations::make(rt_ll(2));
        CHECK(two.has_value());
        if (two) {
            const auto starved = r::newton {}.with_derivative(df).with_budget(*two)(f, x0);
            CHECK_FALSE(starved.has_value());
            if (!starved) {
                CHECK(starved.error().code == nxx::errc::budget_exhausted);
                CHECK(starved.error().used.iterations == 2);
                CHECK(starved.error().best.has_value());
                CHECK(std::abs(nxx::best_x(starved).value_or(qnan) - sqrt2) < 1e-2);
            }
        }
        // solve(f, df, x0) (rtsafe, safeguarded) arrives in phase 3.
    }

    TEST_CASE("call 3: secant with a box clamp, into a user strong type")
    {
        double     lo_seen   = std::numeric_limits<double>::infinity();
        double     hi_seen   = -std::numeric_limits<double>::infinity();
        const auto objective = [&lo_seen, &hi_seen](double x) {
            lo_seen = std::min(lo_seen, x);
            hi_seen = std::max(hi_seen, x);
            return x * x - 2.0;
        };
        const double xmin = rt(0.0);
        const double xmax = rt(3.0);
        const double g    = rt(0.1);

        const auto len =
            r::secant {}.with_projection(r::clamp_to { xmin, xmax })(objective, g).transform([](const auto& s) { return Length { s.x }; });
        static_assert(std::is_same_v<typename std::remove_cvref_t<decltype(len)>::value_type, Length>);
        CHECK(len.has_value());
        if (len) { CHECK(std::abs(len->value - sqrt2) < 1e-12); }
        // The first secant step from 0.1 proposes about 10; the projection clamps it before f is evaluated.
        CHECK(lo_seen >= xmin);
        CHECK(hi_seen <= xmax);
        CHECK(hi_seen == xmax);

        // Without the projection, the same secant evaluates far outside the box.
        hi_seen          = -std::numeric_limits<double>::infinity();
        const auto loose = r::secant {}(objective, g);
        CHECK(loose.has_value());
        CHECK(hi_seen > xmax);
    }

    TEST_CASE("call 3: a pinned iterate is stalled, and best_x accepts the best estimate")
    {
        const auto   objective = [](double x) { return x * x - 2.0; };    // the root sqrt(2) lies outside [0, 1]
        const double xmin      = rt(0.0);
        const double xmax      = rt(1.0);
        const double g         = rt(0.5);

        const auto res = r::secant {}.with_projection(r::clamp_to { xmin, xmax })(objective, g);
        CHECK_FALSE(res.has_value());    // never "the root is the boundary"
        if (!res) {
            CHECK(res.error().code == nxx::errc::stalled);
            CHECK(res.error().by == r::algos::secant);
            CHECK(res.error().best.has_value());
        }
        CHECK(nxx::best_x(res) == std::optional<double>(1.0));    // the best estimate, at the edge
        const double accepted = nxx::best_x(res).value_or(qnan);
        CHECK(accepted == 1.0);

        const auto len = res.transform([](const auto& s) { return Length { s.x }; });
        CHECK_FALSE(len.has_value());
        if (!len) { CHECK(len.error().code == nxx::errc::stalled); }
    }

    TEST_CASE("call 9: a three-solver chain")
    {
        const auto f  = [](double x) { return x * x - 2.0; };
        const auto df = [](double x) { return 2.0 * x; };

        SUBCASE("newton succeeds; nothing else runs")
        {
            const double x0  = rt(1.0);
            const double lo  = rt(1.0);
            const double hi  = rt(2.0);
            const auto   res = nxx::first_of(r::newton {}.with_derivative(df).on(x0), r::secant {}.on(x0), r::brent {}.on({ lo, hi }))(f);
            CHECK(res.has_value());
            if (res) {
                CHECK(res->by == r::algos::newton);
                CHECK(std::abs(res->x - sqrt2) < 1e-15);
            }
            CHECK(same_result(res, r::newton {}.with_derivative(df)(f, x0)));
        }
        SUBCASE("newton fails at f'(0) = 0; a later alternative succeeds and pays for the failed attempt")
        {
            const double x0  = rt(0.0);
            const double lo  = rt(1.0);
            const double hi  = rt(2.0);
            const auto   res = nxx::first_of(r::newton {}.with_derivative(df).on(x0), r::secant {}.on(x0), r::brent {}.on({ lo, hi }))(f);

            const auto nt = r::newton {}.with_derivative(df)(f, x0);
            CHECK_FALSE(nt.has_value());
            if (!nt) {
                CHECK(nt.error().code == nxx::errc::zero_derivative);
                CHECK(nt.error().used == nxx::counters { 1, 2 });    // f(0), then f'(0)
            }
            const auto sc = r::secant {}(f, x0);
            const auto br = r::brent {}(f, { lo, hi });
            CHECK(res.has_value());
            if (res && !nt) {
                CHECK(res->by != r::algos::newton);
                CHECK(std::abs(std::abs(res->x) - sqrt2) < 1e-12);
                if (sc) {
                    CHECK(res->by == r::algos::secant);
                    CHECK(res->x == sc->x);
                    CHECK(res->used == nt.error().used + sc->used);
                }
                else if (br) {
                    CHECK(res->by == r::algos::brent);
                    CHECK(res->used == nt.error().used + sc.error().used + br->used);
                }
            }
        }
        SUBCASE("a non-finite run-time guess: both open methods fail in-band at zero cost; brent succeeds")
        {
            const double x0  = rt(qnan);
            const double lo  = rt(1.0);
            const double hi  = rt(2.0);
            const auto   res = nxx::first_of(r::newton {}.with_derivative(df).on(x0), r::secant {}.on(x0), r::brent {}.on({ lo, hi }))(f);
            CHECK(res.has_value());
            if (res) {
                CHECK(res->by == r::algos::brent);
                CHECK(std::abs(res->x - sqrt2) < 1e-14);
            }
            CHECK(same_result(res, r::brent {}(f, { lo, hi })));
        }
        SUBCASE("everything fails: the last code, the total cost, the best estimate")
        {
            const double x0  = rt(qnan);
            const double lo  = rt(2.0);    // no sign change on [2, 3]
            const double hi  = rt(3.0);
            const auto   res = nxx::first_of(r::newton {}.with_derivative(df).on(x0), r::secant {}.on(x0), r::brent {}.on({ lo, hi }))(f);
            CHECK_FALSE(res.has_value());
            if (!res) {
                CHECK(res.error().code == nxx::errc::no_sign_change);
                CHECK(res.error().by == r::algos::brent);
                CHECK(res.error().used == nxx::counters { 0, 2 });
                CHECK(nxx::best_x(res) == std::optional<double>(2.0));
            }
        }
    }

    TEST_CASE("call 10: coarse bisection, then secant")
    {
        const auto   f  = [](double x) { return x * x - 2.0; };
        const double lo = rt(1.0);
        const double hi = rt(2.0);

        const auto tol = nxx::tolerance<double>::make(rt(1e-3));
        const auto n   = nxx::max_iterations::make(rt_ll(100));
        CHECK(tol.has_value());
        CHECK(n.has_value());
        if (tol && n) {
            const auto res = nxx::then(r::bisection { nxx::width_tol { *tol } }.with_budget(*n).on({ lo, hi }), r::secant {})(f);
            CHECK(res.has_value());
            if (res) {
                CHECK(res->by == r::algos::secant);
                CHECK(std::abs(res->x - sqrt2) < 1e-14);
            }

            // Stage by stage: the coarse enclosure (width 2^-10 <= 1e-3 after 10 halvings), then a seeded secant.
            const auto coarse = r::bisection { nxx::width_tol { *tol } }.with_budget(*n)(f, { lo, hi });
            CHECK(coarse.has_value());
            if (coarse && res) {
                CHECK(coarse->used == nxx::counters { 10, 12 });
                CHECK(coarse->uncertainty <= 1e-3);
                const auto polish = r::secant {}(f, *coarse);    // seeded with x and f(x): no re-evaluation
                CHECK(polish.has_value());
                if (polish) {
                    CHECK(res->x == polish->x);
                    CHECK(res->used == coarse->used + polish->used);    // the total cost of both stages
                }
            }

            // The same pipeline with literals (validated at compile time) gives the same result.
            const auto literal = nxx::then(r::bisection { nxx::width_tol { 1e-3 } }.with_budget(100).on({ lo, hi }), r::secant {})(f);
            CHECK(same_result(literal, res));
        }

        // A run-time budget too small for stage 1: its failure is the pipeline's failure.
        const auto three = nxx::max_iterations::make(rt_ll(3));
        CHECK(three.has_value());
        if (tol && three) {
            const auto res = nxx::then(r::bisection { nxx::width_tol { *tol } }.with_budget(*three).on({ lo, hi }), r::secant {})(f);
            CHECK_FALSE(res.has_value());
            if (!res) {
                CHECK(res.error().code == nxx::errc::budget_exhausted);
                CHECK(res.error().by == r::algos::bisection);
                CHECK(res.error().used == nxx::counters { 3, 5 });
                CHECK(res.error().best.has_value());
            }
        }
    }

    TEST_CASE("call 13: a mixed tolerance literal")
    {
        const auto   f  = [](double x) { return x * x - 2.0; };
        const double lo = rt(1.0);
        const double hi = rt(2.0);

        // The relative part is named, so the roles cannot be swapped (DESIGN §6.2).
        const auto res = r::bisection { nxx::width_tol { 1e-10, nxx::rel_tolerance { 1e-8 } } }(f, { lo, hi });
        CHECK(res.has_value());
        if (res) {
            CHECK(res->how == nxx::stop_reason::criterion);
            CHECK(res->enclosure.has_value());
            if (res->enclosure) {
                const double width = res->enclosure->hi() - res->enclosure->lo();
                CHECK(width <= 1e-10 + 1e-8 * (std::min)(std::abs(res->enclosure->lo()), std::abs(res->enclosure->hi())));
                CHECK(width > 1e-10);    // the relative part counts: 1e-10 alone needs more halvings
            }
            CHECK(std::abs(res->x - sqrt2) <= 1e-10 + 1e-8 * sqrt2);
        }

        // At run time: the relative part through rel_tolerance<T>::make, then make(a, *rel); the same result.
        const auto rel = nxx::rel_tolerance<double>::make(rt(1e-8));
        CHECK(rel.has_value());
        if (rel) {
            const auto tol = nxx::width_tol<double>::make(rt(1e-10), *rel);
            CHECK(tol.has_value());
            if (tol) { CHECK(same_result(r::bisection { *tol }(f, { lo, hi }), res)); }
            CHECK(nxx::width_tol<double>::make(rt(-1e-10), *rel) == std::unexpected(nxx::errc::invalid_input));

            // From a validated tolerance, as a configuration holds it: width_tol{*abs, *rel} needs no check (§6.2).
            const auto abs_tol = nxx::tolerance<double>::make(rt(1e-10));
            CHECK(abs_tol.has_value());
            if (abs_tol) { CHECK(same_result(r::bisection { nxx::width_tol { *abs_tol, *rel } }(f, { lo, hi }), res)); }
        }
    }

    TEST_CASE("call 14: a run-time tolerance")
    {
        const auto   f  = [](double x) { return x * x - 2.0; };
        const double lo = rt(1.0);
        const double hi = rt(2.0);
        const double t  = rt(1e-10);

        bool solved = false;
        if (auto tol = nxx::width_tol<double>::make(t)) {
            const auto res = r::brent { *tol }(f, { lo, hi });
            solved         = true;
            CHECK(res.has_value());
            if (res) {
                CHECK(res->how == nxx::stop_reason::criterion);
                CHECK(std::abs(res->x - sqrt2) <= 1e-10);
                CHECK(same_result(res, r::brent { nxx::width_tol { 1e-10 } }(f, { lo, hi })));    // the literal's result
            }
        }
        CHECK(solved);

        // A tolerance that is not one fails in-band: no solver is built.
        CHECK(nxx::width_tol<double>::make(rt(0.0)) == std::unexpected(nxx::errc::invalid_input));
        CHECK(nxx::width_tol<double>::make(rt(-1e-10)) == std::unexpected(nxx::errc::invalid_input));
        CHECK(nxx::width_tol<double>::make(rt(qnan)) == std::unexpected(nxx::errc::invalid_input));
    }

    TEST_CASE("call 15: the best estimate, on success or failure")
    {
        const auto   f  = [](double x) { return x * x - 2.0; };
        const double lo = rt(1.0);
        const double hi = rt(2.0);

        const auto res = r::brent {}(f, { lo, hi });
        static_assert(std::is_same_v<decltype(nxx::best(res)), std::optional<r::root_estimate<double>>>);
        const auto b = nxx::best(res);
        CHECK(b.has_value());
        if (res && b) { CHECK(*b == static_cast<const r::root_estimate<double>&>(*res)); }    // the solution's estimate

        // A run-time budget too small: the failure's best estimate, through the same call.
        const auto three = nxx::max_iterations::make(rt_ll(3));
        CHECK(three.has_value());
        if (three) {
            const auto starved = r::bisection {}.with_budget(*three)(f, { lo, hi });
            const auto sb      = nxx::best(starved);
            CHECK_FALSE(starved.has_value());
            CHECK(sb.has_value());
            if (!starved) { CHECK(sb == starved.error().best); }
            if (sb) { CHECK(nxx::best_x(starved) == std::optional<double>(sb->x)); }
        }

        // Nothing evaluated: no best estimate.
        CHECK_FALSE(nxx::best(r::brent {}(f, { lo, lo })).has_value());
    }
}

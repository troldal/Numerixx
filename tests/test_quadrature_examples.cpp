// The examples in quadrature.hpp's comments (docs/plan.md, Milestone 4 Phase 5a): each block's calls as written, with
// its free names defined, and the claims its comments make that the other files do not already assert (those are marked
// where they are). The blocks with a method of one's own (Midpoint, SquareAtLowerEnd, Trapezoid) are asserted in
// test_quadrature_protocol.cpp, beside those types. References are closed forms (quadrature_check.hpp); tolerances are
// justified where they are used.

#include "quadrature_check.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <format>
#include <string>
#include <vector>

using namespace quadrature_check;

TEST_CASE("quadrature: the examples of the header comments")
{
    double const inf = dinf;
    auto const fn = [](double x) { return std::exp(x); };
    // The blocks' `use` and `report`, the caller's own: here they record what they are given.
    std::vector<double> used;
    auto const use = [&](double value, double abs_error) { used = { value, abs_error }; };

    // The file comment. The default on e^x over [0, 1]: its estimate on the initial mesh of 10 panels (150
    // evaluations) already meets the default request, so no bisection follows. Its accuracy is asserted in
    // test_quadrature_kronrod.cpp ("K01 e^x" and the smooth cases); here the count, and the honesty check.
    int count = 0;
    auto const r = qd::integrate(Counted<decltype(fn)>{ fn, &count }, 0.0, 1.0); // e - 1, 150 evaluations
    REQUIRE(r.has_value());
    CHECK(count == 150);
    CHECK(error_of(r->value, e_minus_1) <= r->abs_error);
    if (r) use(r->value, r->abs_error);
    CHECK(used == std::vector{ r->value, r->abs_error });
    // The block's `log`, here log_step (a local named log would shadow GCC's built-in): Romberg reports after level 1
    // and after each level until it stops, at level 4 or later, so at least four steps.
    int steps = 0;
    auto const log_step = [&](qd::Step<double> const &) { ++steps; };
    auto const s = qd::integrate<qd::Romberg>(qd::observer([&](qd::Step<double> const &step) { log_step(step); }),
                                              [](double x) { return std::exp(x); }, 0.0, 1.0);
    REQUIRE(s.has_value());
    CHECK(steps >= 4);
    // sqrt(pi): asserted in test_quadrature_infinite.cpp; here the call as written.
    auto const g = qd::integrate<qd::TanhSinh>([](double x) { return std::exp(-x * x); }, -inf, inf); // sqrt(pi)
    REQUIRE(g.has_value());
    char const *const sqrt_pi = "1.772453850905516027298167483341145182797549456122387128213807789852911";
    CHECK(error_of(g->value, sqrt_pi) <= g->abs_error);

    // Criteria's comment, and integrate()'s for a method as a value: Romberg on e^x at rel_tol 1e-12. Its floor,
    // 50 eps L1 = 1.9e-14, lies below the request, so a success has abs_error <= 1e-12 (e - 1); Romberg's estimate
    // overstates the error on smooth fn (design.md, "Romberg"), so the error lies below it.
    auto const romberg = nxx::quadrature::integrate<nxx::quadrature::Romberg>(fn, 0.0, 1.0, { .rel_tol = 1e-12 });
    REQUIRE(romberg.has_value());
    CHECK(romberg->abs_error <= 1e-12 * romberg->value);
    CHECK(error_of(romberg->value, e_minus_1) <= romberg->abs_error);
    CHECK(qd::integrate(qd::Romberg{}, [](double x) { return std::exp(x); }, 0.0, 1.0, { .rel_tol = 1e-12 })
          == romberg);

    // observer()'s comment: print, writing to lines instead of standard output. The last step reported is the one
    // integrate() returns (Romberg's settle changes nothing), so the last line formats the result.
    std::vector<std::string> lines;
    auto const print = [&](qd::Step<double> const &step) {
        if (step.current)
            lines.push_back(std::format("{:.15f} +- {:.3g}", step.current->value, step.current->abs_error));
    };
    auto const printed = qd::integrate<qd::Romberg>(qd::observer(print), fn, 0.0, 1.0, { .rel_tol = 1e-12 });
    REQUIRE(printed.has_value());
    REQUIRE(!lines.empty());
    CHECK(lines.back() == std::format("{:.15f} +- {:.3g}", printed->value, printed->abs_error));
    CHECK(printed == romberg);

    // integrate()'s comments for a fixed rule: GaussLegendre<20> on 1 / (1 + x^2) over [0, 1], pi / 4. The poles at
    // x = +-i map to u = -1 +- 2i on [-1, 1], on the Bernstein ellipse rho = |u + sqrt(u^2 - 1)| = 4.61, so the
    // truncation error, of order rho^-40 = 3e-27, is negligible. The rounding: fn's value within 2 roundings, its
    // product with the weight and the half-width within 2 more, node and weight within an ulp each with |f'| <= 1,
    // so each term within about 6 half-ulps of its size; the terms sum to the integral, 0.785, in a compensated sum,
    // so the result is within 6 / 2 eps 0.785 = 2.4 eps; 4 eps (8.9e-16) bounds it.
    auto const runge = [](double x) { return 1.0 / (1.0 + x * x); };
    auto const legendre =
        qd::integrate(qd::GaussLegendre<20>{}, [](double x) { return 1.0 / (1.0 + x * x); }, 0.0, 1.0);
    REQUIRE(legendre.has_value());
    CHECK(error_of(*legendre, pi_over_4) <= 4.0 * deps);
    CHECK(qd::integrate<qd::GaussLegendre<20>>(runge, 0.0, 1.0) == legendre);
    // `nine`: asserted in test_quadrature_gauss.cpp.

    // integrate()'s comment for a method with an estimate, named: e^x by the default and by Romberg (as above), log x
    // on [0, 1], -1, and x^2 e^-x on [0, inf), 2 (Gamma(3)), by TanhSinh, as written; the honesty check on each.
    auto const t = qd::integrate<qd::TanhSinh>([](double x) { return std::log(x); }, 0.0, 1.0); // -1
    REQUIRE(t.has_value());
    CHECK(std::fabs(t->value + 1.0) <= t->abs_error);
    auto const e = qd::integrate<qd::TanhSinh>([](double x) { return x * x * std::exp(-x); }, 0.0, inf); // 2
    REQUIRE(e.has_value());
    CHECK(std::fabs(e->value - 2.0) <= e->abs_error);
    std::vector<double> reported;
    auto const report = [&](double value, double abs_error) { reported = { value, abs_error }; };
    auto const failed = qd::integrate<qd::TanhSinh>([](double x) { return 1.0 / x; }, 0.0, 1.0);
    if (failed) use(failed->value, failed->abs_error);
    else if (failed.error().best) report(failed.error().best->value, failed.error().best->abs_error);
    // 1/x: not_converged without a best (error's comment), so neither is called.
    REQUIRE(!failed.has_value());
    CHECK(failed.error() == qd::Failure<double>{ qd::error::not_converged, std::nullopt });
    CHECK(used == std::vector{ r->value, r->abs_error });
    CHECK(reported.empty());

    // The fallback block in integrate()'s comment: asserted in test_quadrature_tanh_sinh.cpp. The observer forms:
    // test_quadrature_protocol.cpp ("The observer"). The formatters: test_quadrature_romberg.cpp.
}

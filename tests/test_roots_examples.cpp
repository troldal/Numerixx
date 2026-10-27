// The examples in roots.hpp's comments (docs/sphinx/PHASE2_PLAN.md, "The headers' examples"): each block's calls as
// written, with its free names defined, and the claims its comments make that the other files do not already assert
// (those are marked where they are). Counts are the independent prototypes' or derived here, assuming no FMA
// contraction in fn, as in the other files; accuracy bounds as in roots_check.hpp.

#include "roots_check.hpp"

#include <numerixx/deriv.hpp>

#include <doctest/doctest.h>

#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>

using namespace roots_check;

// Bracket's comment: the deduction from two doubles.
static_assert(std::is_same_v<decltype(rt::Bracket{ 300.0, 400.0 }), rt::Bracket<double>>);

TEST_CASE("roots: the examples of the header comments")
{
    // The file comment and solve's comment for a guess: Newton on f and f' as a std::pair, and Secant on f, from 1,
    // after 7 and 9 evaluations, the prototype's counts that test_roots_open.cpp asserts for a tuple with the same
    // value and derivative. Brent's case of the file comment is in test_roots_bracketing.cpp.
    auto const fdf = [](double x) { return std::pair{ x * x - 2.0, 2.0 * x }; };
    auto const [newton, newton_count] = counted<rt::Newton>(fdf, 1.0);
    CHECK(newton_count == 7);
    CHECK(brackets(newton, sqrt2, 8e-17));
    auto const [secant, secant_count] = counted<rt::Secant>([](double x) { return x * x - 2.0; }, 1.0);
    CHECK(secant_count == 9);
    CHECK(brackets(secant, sqrt2, 8e-17));

    // search's comment: from 300 for exp(20 - 5000/t) - 1, whose root is 5000/20 = 250. The search steps 6 * 1.5^k
    // (k = 0, 1, ...) below and then above 300, every point exact in double; fn > 0 above 250, so the first sign change
    // is below, between 300 - 45.5625 = 254.4375 (fn = 0.42) and 300 - 68.34375 = 231.65625 (fn = -0.80), after
    // 1 + 2 * 6 + 1 = 14 evaluations. |fn| >= 0.42 at every point, so no rounding of exp can change a sign.
    auto const f = [](double t) { return std::exp(20.0 - 5000.0 / t) - 1.0; };
    int searches = 0;
    auto const b = rt::search(Counted<decltype(f)>{ f, &searches }, 300.0);
    REQUIRE(b.has_value());
    CHECK(searches == 14);
    CHECK(*b == rt::Bracket{ 231.65625, 254.4375 });
    // And Brent in it: 5000/t rounds by at most half an ulp of 20 (1.8e-15), 20 - 5000/t and exp(.) - 1 are exact near
    // the root, exp adds at most 2 ulps of 1 (4.4e-16), and |f'(250)| = 5000/250^2 = 0.08, so the sign change is within
    // 2.2e-15 / 0.08 < 2.8e-14 of 250.
    auto const r = rt::search(f, 300.0).and_then([&](rt::Bracket<double> bracket) {
        return rt::solve<rt::Brent>(f, bracket);
    });
    CHECK(brackets(r, 250.0, 2.8e-14));

    // The blocks that make no claim, as written: Criteria's comment, solve's comment for a bracket (psat is an
    // Antoine-like stand-in, 101325 at 373.15 K), and Newton's comment on numerical derivatives (fn is cos x - x).
    auto const fn = [](double x) { return std::cos(x) - x; };
    CHECK(rt::solve<rt::Brent>(fn, rt::Bracket{ 0.0, 1.0 }, { .abs_tol = 1e-12 }).has_value());

    auto const psat = [](double t) { return 101325.0 * std::exp(10.0 * (1.0 - 373.15 / t)); };
    CHECK(rt::solve<rt::Brent>([&](double t) { return psat(t) - 101325.0; }, rt::Bracket{ 300.0, 400.0 }).has_value());
    CHECK(rt::solve<rt::Bisection>(fn, rt::Bracket{ 0.0, 1.0 }, { .rel_tol = 1e-8 }).has_value());
    CHECK(rt::solve<rt::Newton>([](double x) { return std::pair{ x * x - 2.0, 2.0 * x }; }, rt::Bracket{ 1.0, 2.0 })
              .has_value());

    CHECK(rt::solve<rt::Newton>(fn, nxx::deriv::derivative_of(fn), rt::Bracket{ 0.0, 1.0 }).has_value());
    CHECK(rt::solve<rt::Newton>(fn, [&](double x) { return nxx::deriv::derivative(fn, x, 1e-4); }, 1.0).has_value());
    auto const d = nxx::deriv::derivative_of(fn);
    auto const fd = [&](double x) {
        auto const df = d(x);
        return std::pair{ fn(x), df ? df->value : std::numeric_limits<double>::quiet_NaN() };
    };
    CHECK(rt::solve<rt::Newton>(fd, 1.0).has_value());
}

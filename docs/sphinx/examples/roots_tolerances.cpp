// The program of docs/sphinx/how-to/roots-tolerances.md: abs_tol and rel_tol. The page shows the parts between the
// markers; main returns non-zero unless what the page says holds. Counts are the independent prototypes'
// (docs/roots/design.md, Evidence) or derived in a comment, for fn written as here and builds without FMA
// contraction; every fn here is exact near its root, so a root's error is bounded by its abs_error alone.

// clang-format off
// <includes>
#include <numerixx/roots.hpp>

#include <print>
// </includes>
// clang-format on

#include <cmath>
#include <expected>
#include <limits>
#include <utility>

namespace {

    namespace rt = nxx::roots;

    // The number of calls of fn that solve<M> makes, and its result.
    template<typename M, typename F>
    std::pair<std::expected<rt::Root<double>, rt::error>, int> counted(F const &fn, rt::Bracket<double> bracket,
                                                                       rt::Criteria<double> criteria = {})
    {
        int calls = 0;
        auto const r = rt::solve<M>(
            [&](double x) {
                ++calls;
                return fn(x);
            },
            bracket, criteria);
        return {r, calls};
    }

    // The contract every result keeps: abs_error <= 2 tol(value).
    bool within_contract(rt::Root<double> const &r, rt::Criteria<double> const &criteria)
    {
        return r.abs_error <= 2 * (criteria.abs_tol + criteria.rel_tol * std::abs(r.value));
    }

} // namespace

int main()
{
    // clang-format off
    // <near-zero>
    auto const tiny = [](double x) { return x - 1e-300; };
    auto const full = rt::solve<rt::Bisection>(tiny, rt::Bracket{0.0, 1.0});
    auto const enough = rt::solve<rt::Bisection>(tiny, rt::Bracket{0.0, 1.0}, {.abs_tol = 1e-12});
    if (full && enough) {
        std::println("default:        {:.3e} +- {:.1e}", full->value, full->abs_error);
        std::println("abs_tol 1e-12:  {:.3e} +- {:.1e}", enough->value, enough->abs_error);
    }
    // </near-zero>

    // <every-magnitude>
    auto const small = [](double x) { return x - 1e-9; };
    auto const capped = rt::solve<rt::Brent>(small, rt::Bracket{0.0, 1.0}, {.abs_tol = 1e-12});
    if (capped) std::println("{:.6e} +- {:.1e}", capped->value, capped->abs_error);
    // </every-magnitude>

    // <fewer-digits>
    auto const square = [](double x) { return x * x - 2.0; };
    auto const coarse = rt::solve<rt::Bisection>(square, rt::Bracket{1.0, 2.0}, {.rel_tol = 1e-6});
    if (coarse) std::println("{:.10f} +- {:.1e}", coarse->value, coarse->abs_error);
    // </fewer-digits>
    // clang-format on

    if (!full || !enough || !capped || !coarse) return 1;
    rt::Criteria<double> const defaults{};

    // The walk to a sign change near 0: Bisection on x - 1e-300 over [0, 1] takes 1049 evaluations by default, and 41
    // with abs_tol = 1e-12: 2 for the ends, then halvings of [0, 1] until half the width is at most tol, about
    // 1e-12, so 2^-(k+1) <= 1e-12 at k = 39. Brent reaches it in 3, and on a jump at 0 walks as Bisection does, 1076
    // evaluations, where AndersonBjorck takes 51.
    if (!(counted<rt::Bisection>(tiny, {0.0, 1.0}).second == 1049)) return 2;
    if (!(counted<rt::Bisection>(tiny, {0.0, 1.0}, {.abs_tol = 1e-12}).second == 41)) return 3;
    if (!(counted<rt::Brent>(tiny, {0.0, 1.0}).second == 3)) return 4;
    auto const jump = [](double x) { return x < 0.0 ? -1.0 : 1.0; };
    if (!(counted<rt::Brent>(jump, {-1.0, 2.0}).second == 1076)) return 5;
    if (!(counted<rt::AndersonBjorck>(jump, {-1.0, 2.0}).second == 51)) return 6;
    // x - 1e-300 is exact near its root, so the result is within abs_error of 1e-300; with abs_tol = 1e-12 it is at
    // most 2e-12 from it.
    if (!(std::abs(full->value - 1e-300) <= full->abs_error && within_contract(*full, defaults))) return 7;
    rt::Criteria<double> const abs12{.abs_tol = 1e-12};
    if (!(std::abs(enough->value - 1e-300) <= enough->abs_error && within_contract(*enough, abs12))) return 8;

    // abs_tol acts at every magnitude: the root 1e-9 is within abs_error <= 2 (1e-12 + 2 eps 1e-9), about 2e-12, so
    // about 3 significant digits.
    if (!(std::abs(capped->value - 1e-9) <= capped->abs_error && within_contract(*capped, abs12))) return 9;
    if (!(capped->abs_error <= 2.0001e-12)) return 10;

    // Fewer digits: with rel_tol = 1e-6, Bisection on x^2 - 2 stops once 2^-(k+1) <= 1e-6 sqrt(2) = 1.414e-6, at
    // k = 19 (2^-20 = 9.5e-7), so 21 evaluations against 52; the root within abs_error + 8e-17 (x * x rounds by half
    // an ulp of 2, over |f'| = 2.83) of sqrt(2), plus half an ulp of sqrt(2) for the reference in a double.
    if (!(counted<rt::Bisection>(square, {1.0, 2.0}, {.rel_tol = 1e-6}).second == 21)) return 11;
    if (!(counted<rt::Bisection>(square, {1.0, 2.0}).second == 52)) return 12;
    double const sqrt2 = 1.4142135623730950488;
    if (!(std::abs(coarse->value - sqrt2) <= coarse->abs_error + 8e-17 + 1.2e-16)) return 13;
    if (!within_contract(*coarse, {.rel_tol = 1e-6})) return 14;

    // Tolerances below the minimum are raised to it, giving the default result; a negative one is non_finite.
    if (!(rt::solve<rt::Brent>(square, rt::Bracket{1.0, 2.0}, {.abs_tol = 0.0, .rel_tol = 1e-20})
          == rt::solve<rt::Brent>(square, rt::Bracket{1.0, 2.0})))
        return 15;
    if (!(rt::solve<rt::Brent>(square, rt::Bracket{1.0, 2.0}, {.rel_tol = -1.0})
          == std::unexpected(rt::error::non_finite)))
        return 16;
    return 0;
}

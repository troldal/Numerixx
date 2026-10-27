// The program of docs/sphinx/tutorials/roots-kepler.md: Kepler's equation E - e sin E = M solved step by step. The page
// shows the parts between the markers; main returns non-zero unless what the page says about the results holds,
// checked against an independent reference.

// clang-format off
// <includes>
#include <numerixx/roots.hpp>

#include <cmath>
#include <print>
#include <utility>
// </includes>
// clang-format on

#include <limits>

int main()
{
    // E for e = 0.5, M = 1: mpmath's findroot at 100 digits, and the Bessel series E = M + sum 2/n J_n(n e) sin(n M)
    // as a second, independent computation (they agree to 1e-101).
    constexpr double reference = 1.498701133517848314057985497256239901590695552547389812626730414735991;

    // clang-format off
    // <function>
    double const e = 0.5; // the orbit's eccentricity
    double const M = 1.0; // the mean anomaly, in radians
    auto const kepler = [e, M](double E) { return E - e * std::sin(E) - M; };
    // </function>

    // <bracket>
    namespace rt = nxx::roots;
    auto const r = rt::solve<rt::Brent>(kepler, rt::Bracket{M - e, M + e});
    if (r) {
        std::println("E         = {:.17g}", r->value);
        std::println("abs_error = {:.1e}", r->abs_error);
    }
    // </bracket>

    // <check>
    if (r) {
        std::println("actual error = {:.1e}", std::abs(r->value - reference));
        std::println("kepler(E)    = {:.1e}", kepler(r->value));
    }
    // </check>

    // <failure>
    auto const wrong = rt::solve<rt::Brent>(kepler, rt::Bracket{2.0, 3.0});
    if (!wrong) std::println("no root: {}", wrong.error());
    // </failure>

    // <newton>
    auto const kepler_and_slope = [e, M](double E) {
        return std::pair{E - e * std::sin(E) - M, 1.0 - e * std::cos(E)};
    };
    auto const n = rt::solve<rt::Newton>(kepler_and_slope, M);
    if (n) std::println("Newton: E = {:.17g} +- {:.1e}", n->value, n->abs_error);
    // </newton>
    // clang-format on

    // Each check states what must hold, so that a NaN, which fails every comparison, fails it.
    if (!r || !n || wrong) return 1;
    double const eps = std::numeric_limits<double>::epsilon();
    double const bound = 2 * (std::numeric_limits<double>::denorm_min() + 2 * eps * reference);

    // "The actual error is within abs_error": within abs_error + c, c the shift of fn's sign change from the exact
    // root. Near the root sin(E) = 0.997 errs by at most 2 ulps (2.2e-16), times e = 0.5 (exact) 1.1e-16; E - e sin E
    // rounds by half an ulp of 1 (1.1e-16), and - M is exact; the budget 2.2e-16 over f' = 1 - e cos E = 0.964 is
    // 2.3e-16. The reference held in a double adds half an ulp of the root, 1.1e-16.
    double const c = 2.3e-16 + 1.2e-16;
    for (auto const &root : {*r, *n}) {
        if (!(std::abs(root.value - reference) <= root.abs_error + c)) return 2;
        if (!(root.abs_error <= bound)) return 3;
    }
    // "kepler(E) is close to 0": |fn(value)| <= max |f'| (1.5) times the distance to the exact root, plus fn's own
    // rounding (2.2e-16).
    if (!(std::abs(kepler(r->value)) <= 1.5 * (r->abs_error + c) + 2.2e-16)) return 4;
    // [2, 3] lies above the root: fn > 0 at both ends.
    if (!(wrong.error() == rt::error::no_sign_change)) return 5;
    return 0;
}

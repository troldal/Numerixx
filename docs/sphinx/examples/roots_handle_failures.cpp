// The program of docs/sphinx/how-to/roots-handle-failures.md: reacting to a failed solve. The page shows the parts
// between the markers; main returns non-zero unless what the page says holds, against closed forms and mpmath, and
// the independent prototypes' outcomes (docs/roots/design.md, Evidence, Phases 2 and 6).

// clang-format off
// <includes>
#include <numerixx/roots.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <print>
#include <utility>
// </includes>
// clang-format on

int main()
{
    namespace rt = nxx::roots;

    // clang-format off
    // <recover>
    // x^3 - 2x + 2 and its derivative: Newton from 0 cycles between 0 and 1.
    auto const fdf = [](double x) { return std::pair{(x * x - 2.0) * x + 2.0, 3.0 * x * x - 2.0}; };
    auto const r = rt::solve<rt::Newton>(fdf, 0.0).or_else([&](rt::error) {
        return rt::search(fdf, 0.0).and_then([&](rt::Bracket<double> const &b) {
            return rt::solve<rt::Brent>(fdf, b);
        });
    });
    // </recover>

    // <react>
    if (!r) {
        switch (r.error()) {
        case rt::error::no_sign_change:
            std::println("no sign change: widen the bracket, or search from a guess");
            break;
        default:
            std::println("no root: {}", r.error());
            break;
        }
    }
    double const x = r.transform([](rt::Root<double> const &root) { return root.value; })
                         .value_or(std::numeric_limits<double>::quiet_NaN());
    std::println("x = {:.17g}", x);
    // </react>

    // <pole>
    // 1/(x - 1/3) changes sign at its pole, and Brent returns the pole as a root.
    auto const g = [](double t) { return 1.0 / (t - 1.0 / 3.0); };
    auto const p = rt::solve<rt::Brent>(g, rt::Bracket{0.0, 1.0});
    if (p && std::abs(g(p->value)) > std::max(std::abs(g(0.0)), std::abs(g(1.0))))
        std::println("{:.17g} is a pole, not a root", p->value);
    // </pole>
    // clang-format on

    // Newton alone from 0 gives not_converged (Phase 2's prototype: 213 evaluations, its cap).
    if (!(rt::solve<rt::Newton>(fdf, 0.0) == std::unexpected(rt::error::not_converged))) return 1;
    // The recovered root -1.769292354238631415 (mpmath) within abs_error + c: fn's rounding near the root (x * x,
    // * x; - 2 and + 2 exact) about 8e-16 over |f'| = 7.4, 1.1e-16, plus half an ulp of the root, 1.1e-16.
    if (!r || !(std::abs(r->value + 1.769292354238631415) <= r->abs_error + 2.5e-16)) return 2;
    if (!(x == r->value)) return 3;
    // Brent returns the pole (Phase 6's prototype: 55 evaluations), where |g| is far above its size at the ends.
    if (!p || !(std::abs(p->value - 1.0 / 3.0) <= p->abs_error + 1e-16) || !(std::abs(g(p->value)) > 1e10)) return 4;
    return 0;
}

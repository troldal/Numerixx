// The program of docs/sphinx/roots/handle-failures.md: reacting to a failed solve. The page shows the parts
// between the markers and the output in the raw strings; main returns non-zero unless the program prints that output
// and what the page says holds, against closed forms and mpmath, and the independent prototypes' outcomes
// (docs/roots/design.md, Evidence, Phases 2 and 6).

// clang-format off
// <includes>
#include <numerixx/roots.hpp>

#include <cmath>
#include <limits>
#include <print>
#include <utility>
// </includes>
// clang-format on

#include "printed.hpp"

#include <string>

// clang-format off
// <pole-test>
// Whether a result is more likely a pole than a root: |fn| falls away from a pole and grows away from a root.
// Compares |fn| at the result with |fn| a little to each side, 1024 times abs_error away.
template<typename F>
bool looks_like_a_pole(F const &fn, nxx::roots::Root<double> const &r)
{
    double const away = 1024 * r.abs_error;
    double const here = std::abs(fn(r.value));
    return here > std::abs(fn(r.value - away)) && here > std::abs(fn(r.value + away));
}
// </pole-test>
// clang-format on

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;
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
    if (p && looks_like_a_pole(g, *p)) std::println("{:.17g} is a pole, not a root", p->value);
    // </pole>
    // clang-format on

    constexpr std::string_view react_output = R"(
x = -1.7692923542386314
)";
    constexpr std::string_view pole_output = R"(
0.33333333333333309 is a pole, not a root
)";
    if (!(docs::printed() == std::string(docs::quoted(react_output)) + std::string(docs::quoted(pole_output))))
        return 8;

    // Newton alone from 0 gives not_converged (Phase 2's prototype: 213 evaluations, its cap).
    if (!(rt::solve<rt::Newton>(fdf, 0.0) == std::unexpected(rt::error::not_converged))) return 1;
    // The recovered root -1.769292354238631415 (mpmath) within abs_error + c: fn's rounding near the root (x * x,
    // * x; - 2 and + 2 exact) about 8e-16 over |f'| = 7.4, 1.1e-16, plus half an ulp of the root, 1.1e-16.
    if (!r || !(std::abs(r->value + 1.769292354238631415) <= r->abs_error + 2.5e-16)) return 2;
    if (!(x == r->value)) return 3;
    // Brent returns the pole (Phase 6's prototype: 55 evaluations). g's computed sign change lies exactly at the double
    // 1/3, since t - fl(1/3) is exact near it (Sterbenz), so the result is within abs_error of that double with no
    // further term; the test calls it a pole.
    if (!p || !(std::abs(p->value - 1.0 / 3.0) <= p->abs_error) || !looks_like_a_pole(g, *p)) return 4;

    // The test on docs/roots/PHASE1_API.md E's cases, where a test on |fn(value)| alone fails both ways:
    // a resonance whose |fn| at the well-conditioned root sqrt(2) is about 4e8, far above its size at the ends;
    auto const resonance = [](double t) {
        double const d = t * t - 2.0;
        return d / (d * d + 1e-24);
    };
    auto const r1 = rt::solve<rt::Brent>(resonance, rt::Bracket{1.0, 2.0});
    if (!r1 || looks_like_a_pole(resonance, *r1)) return 5;
    // a pole whose |fn| is far below its size at the ends;
    auto const hidden = [](double t) {
        double const d = t - 1.0 / 3.0;
        return 1.0 / d + 1e20 * d * d * d;
    };
    auto const r2 = rt::solve<rt::Brent>(hidden, rt::Bracket{0.0, 1.0});
    if (!r2 || !looks_like_a_pole(hidden, *r2)) return 6;
    // a root at 1/3 whose |fn| at the ends is 2e-49, far below its size near the root.
    auto const narrow = [](double t) {
        double const d = t - 1.0 / 3.0;
        return d * std::exp(-1000.0 * d * d);
    };
    auto const r3 = rt::solve<rt::Brent>(narrow, rt::Bracket{0.0, 1.0});
    if (!r3 || looks_like_a_pole(narrow, *r3)) return 7;
    return 0;
}

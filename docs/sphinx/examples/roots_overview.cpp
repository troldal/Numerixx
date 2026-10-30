// The program of docs/sphinx/roots/index.md: a root from a bracket, and what its error bound covers. The page shows
// the parts between the markers and the output in the raw string; main returns non-zero unless the program prints
// that output and what the page says holds, against closed forms and the rounding analysis in the comments, for fn
// written as here and builds without FMA contraction.

// clang-format off
// <includes>
#include <numerixx/roots.hpp>

#include <print>
// </includes>
// clang-format on

#include "printed.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    // clang-format off
    // <brent>
    namespace rt = nxx::roots;
    auto const r = rt::solve<rt::Brent>([](double x) { return x * x - 2.0; }, rt::Bracket{1.0, 2.0});
    if (r) std::println("root {:.17g}, abs_error {:.1e}", r->value, r->abs_error);
    else std::println("no root: {}", r.error());
    // </brent>

    // <cubic>
    auto const expanded = [](double x) { return ((x - 3.0) * x + 3.0) * x - 1.0; }; // x^3 - 3x^2 + 3x - 1
    auto const factored = [](double x) { return (x - 1.0) * (x - 1.0) * (x - 1.0); };
    auto const e = rt::solve<rt::Brent>(expanded, rt::Bracket{0.0, 3.0});
    auto const f = rt::solve<rt::Brent>(factored, rt::Bracket{0.0, 3.0});
    if (e && f) {
        std::println("expanded: {:.17g}, abs_error {:.1e}", e->value, e->abs_error);
        std::println("factored: {:.17g}, abs_error {:.1e}", f->value, f->abs_error);
    }
    // </cubic>
    // clang-format on

    constexpr std::string_view brent_output = R"(
root 1.4142135623730949, abs_error 6.7e-16
)";
    constexpr std::string_view cubic_output = R"(
expanded: 0.99999581735634269, abs_error 4.4e-16
factored: 0.99999999999999989, abs_error 3.3e-16
)";
    if (!(docs::printed() == std::string(docs::quoted(brent_output)) + std::string(docs::quoted(cubic_output))))
        return 1;

    // Each check states what must hold, so that a NaN, which fails every comparison, fails it.
    if (!r || !e || !f) return 2;
    double const eps = std::numeric_limits<double>::epsilon();
    auto const within_contract = [&](rt::Root<double> const &root) {
        return root.abs_error <= 2 * (std::numeric_limits<double>::denorm_min() + 2 * eps * std::abs(root.value));
    };

    // sqrt(2) within abs_error + c: x * x rounds by at most half an ulp of 2 and the subtraction is exact, so the
    // computed sign change lies within 2.2e-16 / |f'| = 2.2e-16 / 2.83 < 8e-17 of sqrt(2); plus half an ulp of sqrt(2)
    // for the reference held in a double, 1.2e-16. "A relative error of at most about 1e-15": abs_error <= 2 tol, that
    // is 4 eps |value| = 8.9e-16 relative.
    if (!(std::abs(r->value - 1.4142135623730950488) <= r->abs_error + 8e-17 + 1.2e-16 && within_contract(*r)))
        return 3;
    // "To within one unit in the last place": an ulp of sqrt(2) in a double is eps = 2.2e-16. sqrt(2) is held as the
    // nearest double and the rest (mpmath), so that the difference is not rounded: value - hi is exact (Sterbenz).
    double const sqrt2_hi = 1.4142135623730951;
    double const sqrt2_lo = -9.66729331345291e-17;
    if (!(std::abs((r->value - sqrt2_hi) - sqrt2_lo) < eps)) return 7;

    // The factored form computes its sign exactly: x - 1 is exact for x in [0.5, 2] (Sterbenz) and correctly signed
    // elsewhere, and the product of three numbers of one sign keeps it (it underflows only within 1.7e-108 of 1,
    // where no double other than 1 lies). Its only sign change is at 1, so the result is within abs_error of 1.
    if (!(std::abs(f->value - 1.0) <= f->abs_error && within_contract(*f))) return 4;

    // The expanded form: where the computed value has the wrong sign, measured on a grid of step 1e-10 within 1e-5 of
    // 1, and the number of sign changes there. A running error analysis bounds the rounding of Horner's rule near 1:
    // below 1, x - 3 and t * x round by up to 2.2e-16 each (an ulp of 2 and of 1, halved), + 3 is exact (Sterbenz),
    // * x rounds by 1.1e-16 and - 1 is exact, so at most 5.5e-16; above 1, 3.3e-16. |(x - 1)^3| exceeds that beyond
    // (5.5e-16)^(1/3) = 8.2e-6 below 1 and 6.9e-6 above: no computed sign change lies farther from 1. The page says
    // "within about 7e-6", the grid's farthest wrong sign, and "over and over".
    double farthest = 0.0;
    int changes = 0;
    bool positive = expanded(1.0 - 1e-5) > 0.0;
    for (int k = -100000; k <= 100000; ++k) {
        double const x = 1.0 + k * 1e-10;
        double const y = expanded(x);
        if (!((y > 0.0) == (x > 1.0)) || y == 0.0) farthest = std::max(farthest, std::abs(x - 1.0));
        if (!((y > 0.0) == positive)) ++changes;
        positive = y > 0.0;
    }
    if (!(farthest >= 6.5e-6 && farthest <= 7.5e-6 && farthest <= 8.2e-6 && changes > 1000)) return 5;
    // Brent's point: about 4e-6 from 1, inside the noise and within the rounding bound, with abs_error about 4e-16.
    double const distance = std::abs(e->value - 1.0);
    if (!(distance >= 4e-6 && distance <= 4.5e-6 && within_contract(*e)))
        return 6;
    return 0;
}

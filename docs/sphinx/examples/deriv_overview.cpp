// The program of docs/sphinx/deriv/index.md: the central difference of sqrt at 2 with the default step, and Ridders'
// method. The page shows the parts between the markers and the output in the raw strings; main returns non-zero unless
// the program prints that output and what the page says holds, against the closed form f'(2) = 1 / (2 sqrt 2).
// std::sqrt is correctly rounded (IEEE 754), so the output is the same on every target that does not contract
// a * b + c into a fused multiply-add (x86-64 without -mfma, WebAssembly). The figure on the same page, the error
// against the step, is deriv_trade_off.cpp's.

// clang-format off
// <includes>
#include <numerixx/deriv.hpp>

#include <cmath>
#include <print>

namespace dv = nxx::deriv;
// </includes>
// clang-format on

#include "printed.hpp"

#include <algorithm>
#include <limits>
#include <string>
#include <vector>

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    // clang-format off
    // <default>
    auto const fn = [](double x) { return std::sqrt(x); };
    double const exact = std::sqrt(2.0) / 4.0; // f'(2) = 1 / (2 sqrt 2)
    auto const d = dv::derivative(fn, 2.0);    // the central difference, with the default step
    if (d) std::println("default step: relative error {:.1e}", std::abs(d->value - exact) / exact);
    // </default>

    // <ridders>
    auto const r = dv::derivative<dv::Ridders>(fn, 2.0);
    if (r) std::println("value {:.17g}, abs_error {:.1e}, actual error {:.1e}", r->value, *r->abs_error,
                        std::abs(r->value - exact));
    // </ridders>
    // clang-format on

    constexpr std::string_view default_output = R"(
default step: relative error 6.0e-12
)";
    constexpr std::string_view ridders_output = R"(
value 0.35355339059329444, abs_error 8.0e-14, actual error 2.1e-14
)";
    std::string const expected = std::string(docs::quoted(default_output)) + std::string(docs::quoted(ridders_output));
    if (!(docs::printed() == expected)) return 1;

    // Each check states what must hold, so that a NaN, which fails every comparison, fails it. exact is within half an
    // ulp of 1 / (2 sqrt 2), a relative 1.1e-16, far below every error compared.
    if (!d || !r || !r->abs_error) return 2;

    // The default step, eps^(1/3) max(|x|, 1) = 1.2e-5 at x = 2, gives about the best error of the figure's curve
    // (deriv_trade_off.cpp): within ten times the smallest error of the central difference at the steps 1, 0.5, 0.2,
    // 0.1, ... down to 1e-13, each the double nearest its decimal value.
    double const step = dv::default_step<dv::Central<2>>(2.0);
    if (!(step > 1.2e-5 && step < 1.22e-5)) return 3;
    std::vector<double> hs;
    double p = 1.0; // 10^j, exact
    for (int j = 0; j < 13; ++j, p *= 10.0) hs.insert(hs.end(), {1.0 / p, 5.0 / (10.0 * p), 2.0 / (10.0 * p)});
    hs.push_back(1.0 / p);
    double best = std::numeric_limits<double>::infinity();
    for (double const h : hs) {
        auto const s = dv::derivative(fn, 2.0, h);
        if (!s) return 4;
        best = std::min(best, std::abs(s->value - exact) / exact);
    }
    if (!(std::abs(d->value - exact) / exact <= 10.0 * best)) return 5;
    // Ridders: its estimate covers the actual error, about 1e-13 relative or better.
    if (!(std::abs(r->value - exact) <= *r->abs_error && std::abs(r->value - exact) <= 1e-13 * exact)) return 6;
    // "About a hundred times smaller than the central difference's" error, and an estimate "a few times larger" than
    // the actual error: between 50 and 500 times, and between 1 and 10 times.
    double const ratio = std::abs(d->value - exact) / std::abs(r->value - exact);
    if (!(ratio > 50 && ratio < 500 && *r->abs_error < 10 * std::abs(r->value - exact))) return 7;
    return 0;
}

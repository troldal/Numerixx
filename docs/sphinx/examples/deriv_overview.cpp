// The program of docs/sphinx/deriv/index.md: the central difference of sqrt at 2 for a range of steps, the default
// step, and Ridders' method. The page shows the parts between the markers and the output in the raw strings; main
// returns non-zero unless the program prints that output and what the page says holds, against the closed form
// f'(2) = 1 / (2 sqrt 2). std::sqrt is correctly rounded (IEEE 754), so the output is the same on every target that
// does not contract a * b + c into a fused multiply-add (x86-64 without -mfma, WebAssembly).

// clang-format off
// <includes>
#include <numerixx/deriv.hpp>

#include <cmath>
#include <initializer_list>
#include <print>
// </includes>
// clang-format on

#include "printed.hpp"

#include <limits>
#include <string>

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    // clang-format off
    // <steps>
    namespace dv = nxx::deriv;
    auto const fn = [](double x) { return std::sqrt(x); };
    double const exact = std::sqrt(2.0) / 4.0; // f'(2) = 1 / (2 sqrt 2)
    for (double const h : {1e-1, 1e-3, 1e-5, 1e-7, 1e-9, 1e-11}) {
        auto const d = dv::derivative(fn, 2.0, h);
        if (d) std::println("h = {:.0e}   relative error {:.1e}", h, std::abs(d->value - exact) / exact);
    }
    // </steps>

    // <default>
    auto const d = dv::derivative(fn, 2.0); // the default step
    if (d) std::println("default step: relative error {:.1e}", std::abs(d->value - exact) / exact);
    // </default>

    // <ridders>
    auto const r = dv::derivative<dv::Ridders>(fn, 2.0);
    if (r) std::println("value {:.17g}, abs_error {:.1e}, actual error {:.1e}", r->value, *r->abs_error,
                        std::abs(r->value - exact));
    // </ridders>
    // clang-format on

    constexpr std::string_view steps_output = R"(
h = 1e-01   relative error 3.1e-04
h = 1e-03   relative error 3.1e-08
h = 1e-05   relative error 5.1e-12
h = 1e-07   relative error 1.4e-09
h = 1e-09   relative error 3.1e-08
h = 1e-11   relative error 8.2e-06
)";
    constexpr std::string_view default_output = R"(
default step: relative error 6.0e-12
)";
    constexpr std::string_view ridders_output = R"(
value 0.35355339059329444, abs_error 8.0e-14, actual error 2.1e-14
)";
    std::string const expected = std::string(docs::quoted(steps_output)) + std::string(docs::quoted(default_output))
                                 + std::string(docs::quoted(ridders_output));
    if (!(docs::printed() == expected)) return 1;

    // Each check states what must hold, so that a NaN, which fails every comparison, fails it.
    if (!d || !r || !r->abs_error) return 2;
    double const eps = std::numeric_limits<double>::epsilon();

    // The error model for each step: truncation h^2 / 6 max |f'''| over [2 - h, 2 + h], with f''' = 3/8 x^(-5/2)
    // largest at 2 - h, plus rounding: each sqrt is within half an ulp, eps/2 sqrt(2 + h), and the difference of the
    // two over the distance dx = 2h (the points' own rounding is divided out, as the code divides by the actual
    // distance) adds eps sqrt(2 + h) / 2h; all relative to f'(2) = 0.3536, plus 1.2e-16 for exact's own rounding.
    for (double const h : {1e-1, 1e-3, 1e-5, 1e-7, 1e-9, 1e-11}) {
        auto const s = dv::derivative(fn, 2.0, h);
        if (!s) return 3;
        double const truncation = h * h / 6.0 * 0.375 * std::pow(2.0 - h, -2.5);
        double const rounding = eps * std::sqrt(2.0 + h) / (2.0 * h);
        if (!(std::abs(s->value - exact) <= 1.0001 * (truncation + rounding) + 1.2e-16 * exact)) return 4;
    }
    // "The truncation error falls as h^2" at the large steps: at 1e-1 and 1e-3 the error is the truncation term to
    // within 15% (the next term of the series is smaller by a factor 0.109 h^2 at x = 2).
    for (double const h : {1e-1, 1e-3}) {
        double const err = std::abs(dv::derivative(fn, 2.0, h)->value - exact);
        double const model = h * h / 6.0 * 0.375 * std::pow(2.0, -2.5);
        if (!(err > 0.85 * model && err < 1.15 * model)) return 5;
    }
    // The default step, eps^(1/3) max(|x|, 1) = 1.2e-5 at x = 2, gives about the best error of the table: within
    // ten times the best measured (5.1e-12 relative at 1e-5).
    if (!(std::abs(d->value - exact) / exact <= 5.1e-11)) return 6;
    // Ridders: its estimate covers the actual error, about 1e-13 relative or better.
    if (!(std::abs(r->value - exact) <= *r->abs_error && std::abs(r->value - exact) <= 1e-13 * exact)) return 7;
    // "About a hundred times smaller than the central difference's" error, and an estimate "a few times larger" than
    // the actual error: between 50 and 500 times, and between 1 and 10 times.
    double const ratio = std::abs(d->value - exact) / std::abs(r->value - exact);
    if (!(ratio > 50 && ratio < 500 && *r->abs_error < 10 * std::abs(r->value - exact))) return 8;
    return 0;
}

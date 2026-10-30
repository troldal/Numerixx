// The program of docs/sphinx/deriv/step.md: a function that varies on a scale far below max(|x|, 1), where the default
// step is too large, x / (x^2 + a^2) with a = 1e-4 at x = 2e-4, whose derivative is (a^2 - x^2) / (x^2 + a^2)^2. The
// page shows the parts between the markers and the output in the raw string; main returns non-zero unless the program
// prints that output and what the page says holds. fn and the steps use only + - * / (nxx::math::root_epsilon, not
// std::cbrt, which need not be correctly rounded), so the output is the same on every target that does not contract
// a * b + c into a fused multiply-add (x86-64 without -mfma, WebAssembly).

// clang-format off
// <includes>
#include <numerixx/deriv.hpp>

#include <cmath>
#include <print>
// </includes>
// clang-format on

#include "printed.hpp"

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;
    namespace dv = nxx::deriv;

    // clang-format off
    // <scale>
    // A peak of width a = 1e-4: f changes on the scale of a, far below max(|x|, 1) = 1.
    double const a = 1e-4;
    auto const fn = [a](double x) { return x / (x * x + a * a); };
    double const x = 2e-4;
    double const exact = (a * a - x * x) / ((x * x + a * a) * (x * x + a * a)); // -1.2e7
    // </scale>

    // <default>
    auto const coarse = dv::derivative(fn, x); // h = 6e-6, about a sixteenth of the peak's width
    if (coarse) std::println("default step:   relative error {:.1e}", std::abs(coarse->value - exact) / -exact);
    // </default>

    // <scaled>
    // The step for f's own scale: the cube root of epsilon, times the scale a.
    double const h = nxx::math::root_epsilon<3, double>() * a;
    auto const scaled = dv::derivative(fn, x, h);
    if (scaled) std::println("h = {:.1e}:     relative error {:.1e}", h, std::abs(scaled->value - exact) / -exact);
    // </scaled>

    // <ridders>
    auto const r = dv::derivative<dv::Ridders>(fn, x);
    if (r) std::println("Ridders:        relative error {:.1e}, abs_error {:.1e}", std::abs(r->value - exact) / -exact,
                        *r->abs_error / -exact);
    // </ridders>

    // <ridders-step>
    // Given a step, Ridders' method starts from it: a / 8 is on the peak's own scale.
    auto const started = dv::derivative<dv::Ridders>(fn, x, a / 8);
    if (started) std::println("Ridders from a / 8: relative error {:.1e}", std::abs(started->value - exact) / -exact);
    // </ridders-step>
    // clang-format on

    constexpr std::string_view output = R"(
default step:   relative error 3.4e-04
h = 6.1e-10:     relative error 4.6e-11
Ridders:        relative error 3.3e-14, abs_error 1.9e-13
Ridders from a / 8: relative error 5.0e-14
)";
    if (!(docs::printed() == docs::quoted(output))) return 1;

    // Each check states what must hold, so that a NaN fails it. exact is computed in a few roundings, about 4 eps
    // relative, far below every error checked.
    if (!coarse || !scaled || !r || !r->abs_error) return 2;
    // The default step, eps^(1/3) max(|x|, 1) = 6.06e-6, is about a sixteenth of the peak's width: the truncation
    // error, h^2 / 6 |f'''| / |f'| = 0.0933 (h / a)^2 = 3.4e-4. The page says "three or four digits".
    double const step = dv::default_step<dv::Central<2>, 1>(x);
    if (!(step > 6e-6 && step < 6.1e-6)) return 3;
    double const coarse_error = std::abs(coarse->value - exact) / -exact;
    if (!(coarse_error > 1e-4 && coarse_error < 1e-3)) return 4;
    // The scaled step, 6.06e-10: truncation 0.0933 (h / a)^2 = 3.4e-12, rounding about 1.5 eps |f| / (h |f'|) with
    // |f| = 4000 and |f'| = 1.2e7, 1.9e-10 relative; the page says "about 10 digits".
    double const scaled_error = std::abs(scaled->value - exact) / -exact;
    if (!(scaled_error < 3e-10)) return 5;
    // Ridders finds the scale itself: "13 to 15 digits", its estimate covering the actual error and below 1e-12
    // relative.
    if (!(std::abs(r->value - exact) <= *r->abs_error && *r->abs_error < 1e-12 * -exact)) return 6;
    // "At the cost of 36 evaluations", against 2 for the central difference.
    int calls = 0;
    (void)dv::derivative<dv::Ridders>(
        [&](double t) {
            ++calls;
            return fn(t);
        },
        x);
    if (!(calls == 36)) return 7;
    // From a / 8: 14 evaluations, as accurate.
    calls = 0;
    auto const from_step = dv::derivative<dv::Ridders>(
        [&](double t) {
            ++calls;
            return fn(t);
        },
        x, a / 8);
    if (!(calls == 14 && started && from_step == started)) return 8;
    if (!(std::abs(started->value - exact) <= *started->abs_error && *started->abs_error < 1e-12 * -exact)) return 9;
    return 0;
}

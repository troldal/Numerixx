// The program of docs/sphinx/deriv/higher-derivatives.md: derivatives of order 1 to 4 of 1/x at 2, whose derivatives
// are -1/4, 1/4, -3/8 and 3/4, by Central<2> and Ridders, and f'' in float. The page shows the parts between the
// markers and the output in the raw strings; main returns non-zero unless the program prints that output and what the
// page says holds. fn uses only + - * /, so the output is the same on every target that does not contract a * b + c
// into a fused multiply-add (x86-64 without -mfma, WebAssembly).

// clang-format off
// <silence>
// Turns off the low-precision warning for every derivative above the first in float, in this program. It must come
// before the header, and be defined alike in every file of the program.
#define NUMERIXX_NO_PRECISION_WARNING
// </silence>

// <includes>
#include <numerixx/deriv.hpp>

#include <cmath>
#include <print>
// </includes>
// clang-format on

#include "printed.hpp"

#include <array>
#include <string>

namespace {

    namespace dv = nxx::deriv;

    // The relative error of derivative<M, N> of 1/x at 2.
    template<typename M, int N>
    double relative_error(double exact)
    {
        auto const d = dv::derivative<M, N>([](double x) { return 1.0 / x; }, 2.0);
        return d ? std::abs(d->value - exact) / std::abs(exact) : -1.0;
    }

} // namespace

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    // clang-format off
    // <orders>
    auto const fn = [](double x) { return 1.0 / x; };
    auto const second = dv::derivative<dv::Central<2>, 2>(fn, 2.0);           // f'', as an integer
    auto const fourth = dv::derivative<dv::Ridders, dv::order::fourth>(fn, 2.0); // f'''', as an order
    if (second) std::println("f''(2)   = {:.10f}", second->value);
    if (fourth) std::println("f''''(2) = {:.10f}, abs_error {:.1e}", fourth->value, *fourth->abs_error);
    // </orders>

    // <float>
    auto const in_float = dv::derivative<dv::Central<2>, 2>([](float x) { return 1.0f / x; }, 2.0f);
    if (in_float) std::println("f''(2) in float = {}", in_float->value);
    // </float>
    // clang-format on

    // The table of the page, by the same calls: relative errors of Central<2> and Ridders for f' to f''''.
    std::array<double, 4> const exact{-0.25, 0.25, -0.375, 0.75};
    std::array<double, 4> const central{relative_error<dv::Central<2>, 1>(exact[0]),
                                        relative_error<dv::Central<2>, 2>(exact[1]),
                                        relative_error<dv::Central<2>, 3>(exact[2]),
                                        relative_error<dv::Central<2>, 4>(exact[3])};
    std::array<double, 4> const ridders{relative_error<dv::Ridders, 1>(exact[0]),
                                        relative_error<dv::Ridders, 2>(exact[1]),
                                        relative_error<dv::Ridders, 3>(exact[2]),
                                        relative_error<dv::Ridders, 4>(exact[3])};
    for (int n = 0; n < 4; ++n)
        std::println("N = {}   Central<2> {:.1e}   Ridders {:.1e}", n + 1, central[static_cast<std::size_t>(n)],
                     ridders[static_cast<std::size_t>(n)]);

    constexpr std::string_view orders_output = R"(
f''(2)   = 0.2500000037
f''''(2) = 0.7500000013, abs_error 1.4e-06
)";
    constexpr std::string_view float_output = R"(
f''(2) in float = 0.2500806
)";
    constexpr std::string_view table_output = R"(
N = 1   Central<2> 3.7e-11   Ridders 1.6e-15
N = 2   Central<2> 1.5e-08   Ridders 3.3e-13
N = 3   Central<2> 2.7e-06   Ridders 3.0e-11
N = 4   Central<2> 2.9e-05   Ridders 1.8e-09
)";
    std::string const expected = std::string(docs::quoted(orders_output)) + std::string(docs::quoted(float_output))
                                 + std::string(docs::quoted(table_output));
    if (!(docs::printed() == expected)) return 1;

    // Each check states what must hold, so that a NaN fails it.
    if (!second || !fourth || !fourth->abs_error || !in_float) return 2;
    // The accuracy falls with the order. Central<2> at its default step eps^(1/(2+N)) max(|x|, 1) has a model error of
    // order eps^(2/(2+N)): 3.7e-11, 1.5e-8, 5.5e-7 and 6.1e-6 for N = 1 to 4. The bounds for 1/x at 2, truncation plus
    // rounding: N = 1, h = 1.2e-5, 3.7e-11 + 3.7e-11; N = 2, h = 2^-12, h^2 / 4 = 1.5e-8 + 1.5e-8; N = 3,
    // h = 1.48e-3, 1.25 h^2 = 2.74e-6 + 1.4e-7 (values) + 1.4e-7 (points); N = 4, h = 4.92e-3, 1.25 h^2 = 3.0e-5 +
    // 4e-6. Ridders is far more accurate at every order, and its estimate covers its error.
    if (!(central[0] < 4e-10 && central[1] < 6e-8 && central[2] < 4e-6 && central[3] < 6e-5)) return 3;
    // "About 10 digits for f', 8 for f'', 5 or 6 for f''' and 4 or 5 for f''''" (the constants the default step
    // leaves out cost up to a digit at the higher orders).
    if (!(central[0] < 1e-9 && central[1] < 1e-7 && central[2] < 1e-5 && central[2] > 1e-7 && central[3] < 1e-4
          && central[3] > 1e-6))
        return 8;
    // "3, 4 and 5 evaluations for Central<2>" for f'' to f'''': P + N - 1 points.
    int calls = 0;
    auto const counted = [&](double t) {
        ++calls;
        return 1.0 / t;
    };
    std::array<int, 3> const nodes{(calls = 0, (void)dv::derivative<dv::Central<2>, 2>(counted, 2.0), calls),
                                   (calls = 0, (void)dv::derivative<dv::Central<2>, 3>(counted, 2.0), calls),
                                   (calls = 0, (void)dv::derivative<dv::Central<2>, 4>(counted, 2.0), calls)};
    if (!(nodes == std::array<int, 3>{3, 4, 5})) return 9;
    for (std::size_t n = 0; n < 4; ++n) {
        if (!(central[n] > 0.0 && ridders[n] >= 0.0 && ridders[n] < central[n] / 100)) return 4;
        if (n > 0 && !(central[n] > central[n - 1] && ridders[n] > ridders[n - 1])) return 5;
    }
    if (!(std::abs(fourth->value - 0.75) <= *fourth->abs_error)) return 6;
    // float: f'' is good to only about 3 or 4 digits, as the warning says. Model: eps_f^(1/2) = 3.5e-4 for Central<2>.
    if (!(std::abs(in_float->value - 0.25f) / 0.25f < 1e-3f)) return 7;
    return 0;
}

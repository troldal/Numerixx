// The program of docs/sphinx/deriv/one-sided.md: a function defined for x >= 0 only, 1 / (1 + x) with NaN below 0,
// whose derivative is -1 / (1 + x)^2: at the edge 0 and at 1e-6 from it. The page shows the parts between the markers
// and the output in the raw strings; main returns non-zero unless the program prints that output and what the page
// says holds. fn uses only + - * /, so the output is the same on every target that does not contract a * b + c into a
// fused multiply-add (x86-64 without -mfma, WebAssembly).

// clang-format off
// <includes>
#include <numerixx/deriv.hpp>

#include <cmath>
#include <limits>
#include <print>
// </includes>
// clang-format on

#include "printed.hpp"

#include <string>

namespace {

    namespace dv = nxx::deriv;

    // The relative error of a result, or -1 for a failure.
    template<typename R>
    double relative_error(R const &d, double exact)
    { return d ? std::abs(d->value - exact) / std::abs(exact) : -1.0; }

} // namespace

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    // clang-format off
    // <edge>
    // 1 / (1 + x), defined for x >= 0 only: NaN below 0.
    auto const fn = [](double x) {
        if (x < 0.0) return std::numeric_limits<double>::quiet_NaN();
        return 1.0 / (1.0 + x);
    };
    auto const central = dv::derivative(fn, 0.0);
    auto const forward = dv::derivative<dv::Forward<2>>(fn, 0.0);
    auto const forward_ridders = dv::derivative<dv::ForwardRidders>(fn, 0.0);
    if (!central) std::println("Central<2>:     {}", central.error());
    if (forward) std::println("Forward<2>:     {:.15f}", forward->value);
    if (forward_ridders) std::println("ForwardRidders: {:.15f}, abs_error {:.1e}", forward_ridders->value,
                                      *forward_ridders->abs_error);
    // </edge>

    // <near>
    // Near the edge, not at it: Ridders' steps cross 0, and it restarts with smaller ones.
    auto const near = dv::derivative<dv::Ridders>(fn, 1e-6);
    if (near) std::println("Ridders at 1e-6: {:.15f}, abs_error {:.1e}", near->value, *near->abs_error);
    // </near>

    // <kink>
    // A formula that goes on below 0 with another slope: no NaN tells the methods where the edge is.
    auto const kink = [](double x) { return x >= 0.0 ? 1.0 / (1.0 + x) : 1.0 - 3.0 * x; };
    auto const both_sides = dv::derivative<dv::Ridders>(kink, 0.0);
    auto const one_side = dv::derivative<dv::ForwardRidders>(kink, 0.0);
    if (both_sides) std::println("Ridders:        {:.9f}, abs_error {:.1e}", both_sides->value, *both_sides->abs_error);
    if (one_side) std::println("ForwardRidders: {:.9f}, abs_error {:.1e}", one_side->value, *one_side->abs_error);
    // </kink>
    // clang-format on

    constexpr std::string_view edge_output = R"(
Central<2>:     non_finite
Forward<2>:     -0.999999999886047
ForwardRidders: -0.999999999999835, abs_error 4.8e-12
)";
    constexpr std::string_view near_output = R"(
Ridders at 1e-6: -0.999997999984771, abs_error 9.3e-10
)";
    constexpr std::string_view kink_output = R"(
Ridders:        -1.999999991, abs_error 1.2e-08
ForwardRidders: -1.000000000, abs_error 4.8e-12
)";
    std::string const expected = std::string(docs::quoted(edge_output)) + std::string(docs::quoted(near_output))
                                 + std::string(docs::quoted(kink_output));
    if (!(docs::printed() == expected)) return 1;

    // Each check states what must hold, so that a NaN fails it.
    if (!forward || !forward_ridders || !forward_ridders->abs_error || !near || !near->abs_error) return 2;
    // At the edge, the central methods fail: their points reach below 0, where fn is NaN.
    if (!(central == std::unexpected(dv::error::non_finite))) return 3;
    if (!(dv::derivative<dv::Ridders>(fn, 0.0) == std::unexpected(dv::error::non_finite))) return 4;
    // Forward<2>: f'(0) = -1 to its model's accuracy, h^2 / 3 |f'''| = h^2 / 3 * 6 with h = eps^(1/3) = 6.06e-6, about
    // 7e-11, plus rounding about S eps / h = 4 * 2.2e-16 / 6e-6 = 1.5e-10; the page says "about 10 digits".
    if (!(relative_error(forward, -1.0) < 3e-10)) return 5;
    // ForwardRidders: its estimate covers the actual error, below 1e-12 (the page: "about 13 digits").
    if (!(std::abs(forward_ridders->value + 1.0) <= *forward_ridders->abs_error && *forward_ridders->abs_error < 1e-11
          && std::abs(forward_ridders->value + 1.0) < 1e-12))
        return 6;
    // At 1e-6, the central difference's default step, 6.06e-6, still reaches below 0.
    if (!(dv::derivative(fn, 1e-6) == std::unexpected(dv::error::non_finite))) return 8;
    // Near the edge: Ridders' estimate covers its error. -1 / (1 + 1e-6)^2 is held to within 2 ulps.
    double const exact_near = -1.0 / ((1.0 + 1e-6) * (1.0 + 1e-6));
    if (!(std::abs(near->value - exact_near) <= *near->abs_error + 4.5e-16 && *near->abs_error < 1e-8)) return 7;
    // "On this function ForwardRidders is more accurate at 1e-6" than Ridders.
    auto const forward_near = dv::derivative<dv::ForwardRidders>(fn, 1e-6);
    if (!forward_near || !(std::abs(forward_near->value - exact_near) < std::abs(near->value - exact_near))) return 12;
    // "On sqrt(1 - x) near 1, Ridders' method can be the more accurate one", and a one-sided stencil with its default
    // step fails there: g(x) = sqrt(1 - x) + x, defined up to 1, g'(x) = 1 - 1 / (2 sqrt(1 - x)). At 0.999 Ridders
    // beats BackwardRidders; at 1 - 1e-7, Backward<2>'s step, 6e-6, spans sixty times the distance to the edge, and
    // its result is more than 10% off, with no estimate.
    auto const g = [](double x) {
        if (x > 1.0) return std::numeric_limits<double>::quiet_NaN();
        return std::sqrt(1.0 - x) + x;
    };
    double const at = 0.999;
    double const g_exact = 1.0 - 0.5 / std::sqrt(1.0 - at);
    auto const g_central = dv::derivative<dv::Ridders>(g, at);
    auto const g_backward = dv::derivative<dv::BackwardRidders>(g, at);
    if (!g_central || !g_backward
        || !(std::abs(g_central->value - g_exact) < std::abs(g_backward->value - g_exact)))
        return 13;
    double const close = 1.0 - 1e-7;
    auto const g_stencil = dv::derivative<dv::Backward<2>>(g, close);
    double const g_close = 1.0 - 0.5 / std::sqrt(1.0 - close);
    if (!g_stencil || !(std::abs(g_stencil->value - g_close) > 0.1 * std::abs(g_close))) return 14;
    // "On a function whose values happen to vanish at the points a method tries, [the estimate] can even be 0 for a
    // wrong result" (Codex's case): u (1 - (64 u)^2) (1 - (128 u)^2), u = t - 1/64, NaN below 0, at 1/64, where the
    // derivative is 1. Ridders' first step 1/32 reaches below 0; its next steps, 1/64 and 1/128, land on zeros of the
    // polynomial, so it returns 0 with an estimate of 0. ForwardRidders returns 1.
    auto const vanishing = [](double t) {
        if (t < 0.0) return std::numeric_limits<double>::quiet_NaN();
        double const u = t - 1.0 / 64.0;
        return u * (1.0 - (64.0 * u) * (64.0 * u)) * (1.0 - (128.0 * u) * (128.0 * u));
    };
    auto const fooled = dv::derivative<dv::Ridders>(vanishing, 1.0 / 64.0);
    auto const right = dv::derivative<dv::ForwardRidders>(vanishing, 1.0 / 64.0);
    if (!fooled || !fooled->abs_error || !(fooled->value == 0.0 && *fooled->abs_error == 0.0)) return 15;
    if (!right || !right->abs_error || !(std::abs(right->value - 1.0) <= *right->abs_error)) return 16;
    // The kink: Ridders returns about the mean of the slopes -3 and -1, wrong, with a small estimate; ForwardRidders
    // sees only the right-hand formula and returns -1 within its estimate.
    if (!both_sides || !both_sides->abs_error || !one_side || !one_side->abs_error) return 9;
    if (!(std::abs(both_sides->value + 2.0) < 1e-6 && *both_sides->abs_error < 1e-6)) return 10;
    if (!(std::abs(one_side->value + 1.0) <= *one_side->abs_error)) return 11;
    return 0;
}

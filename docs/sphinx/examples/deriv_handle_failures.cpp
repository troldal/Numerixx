// The program of docs/sphinx/deriv/handle-failures.md: a failed derivative, and results whose estimates say they are
// not accurate enough: sqrt x at 1e-11 and 1e-12 by ForwardRidders, whose derivative is 1 / (2 sqrt x). The page shows
// the parts between the markers and the output in the raw strings; main returns non-zero unless the program prints that
// output and what the page says holds. std::sqrt is correctly rounded (IEEE 754), so the output is the same on every
// target that does not contract a * b + c into a fused multiply-add (x86-64 without -mfma, WebAssembly).

// clang-format off
// <includes>
#include <numerixx/deriv.hpp>

#include <cmath>
#include <print>
// </includes>
// clang-format on

#include "printed.hpp"

#include <expected>
#include <string>

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;
    namespace dv = nxx::deriv;

    // clang-format off
    // <failure>
    auto const fn = [](double x) { return std::sqrt(x); };
    auto const central = dv::derivative(fn, 1e-12); // its points reach below 0, where sqrt is NaN
    if (!central) std::println("no derivative: {}", central.error());
    // </failure>

    // <trust>
    // Whether a result is accurate enough for your purpose: it has an estimate, and the estimate is at most the
    // relative error you can accept.
    auto const good_enough = [](dv::Estimate<double> const &d, double relative) {
        return d.abs_error && *d.abs_error <= relative * std::abs(d.value);
    };
    auto const near_zero = dv::derivative<dv::ForwardRidders>(fn, 1e-11);
    if (near_zero && !good_enough(*near_zero, 1e-8))
        std::println("not good enough: {:.6e}, abs_error {:.1e}", near_zero->value, *near_zero->abs_error);
    // A step on the function's own scale, a fraction of x:
    auto const scaled = dv::derivative<dv::ForwardRidders>(fn, 1e-11, 1e-11 / 8);
    if (scaled && good_enough(*scaled, 1e-8))
        std::println("with h = x / 8:  {:.6e}, abs_error {:.1e}", scaled->value, *scaled->abs_error);
    // </trust>
    // clang-format on

    constexpr std::string_view failure_output = R"(
no derivative: non_finite
)";
    constexpr std::string_view trust_output = R"(
not good enough: 1.472966e+05, abs_error 1.6e+04
with h = x / 8:  1.581139e+05, abs_error 5.0e-07
)";
    if (!(docs::printed() == std::string(docs::quoted(failure_output)) + std::string(docs::quoted(trust_output))))
        return 1;

    // Each check states what must hold, so that a NaN fails it.
    if (!(central == std::unexpected(dv::error::non_finite))) return 2;
    if (!near_zero || !near_zero->abs_error || !scaled || !scaled->abs_error) return 3;
    // At 1e-11 the result is 7% off ("off by 7%"), its estimate covers that and fails the 1e-8 test, though it is below
    // a quarter of the value; the scaled step's result is within its estimate, near 1e-12 relative. The reference
    // 0.5 / sqrt(1e-11) is within an ulp (1.2e-16 relative), far below these.
    double const exact = 0.5 / std::sqrt(1e-11);
    double const off = std::abs(near_zero->value - exact) / exact;
    if (!(off > 0.06 && off < 0.08 && std::abs(near_zero->value - exact) <= *near_zero->abs_error)) return 4;
    if (!(*near_zero->abs_error < std::abs(near_zero->value) / 4)) return 7;
    if (!(std::abs(scaled->value - exact) <= *scaled->abs_error + 1e-10 && *scaled->abs_error / exact < 1e-11))
        return 5;
    // At 1e-12 the method gives up: 70% off, and its estimate at least a quarter of the value, its own flag.
    auto const gave_up = dv::derivative<dv::ForwardRidders>(fn, 1e-12);
    double const exact12 = 0.5 / std::sqrt(1e-12);
    if (!gave_up || !gave_up->abs_error || !(std::abs(gave_up->value - exact12) / exact12 > 0.6)
        || !(*gave_up->abs_error >= std::abs(gave_up->value) / 4))
        return 8;
    // A step below the resolution of x is an error too: x + h rounds to x.
    if (!(dv::derivative(fn, 1.0, 1e-300) == std::unexpected(dv::error::non_finite))) return 6;
    return 0;
}

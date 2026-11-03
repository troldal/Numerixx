// The program of docs/sphinx/quadrature/index.md: the integral of 4 / (1 + x^2) over [0, 1], which is pi, by the
// default method and by the other three, and the failure on 1/x over [0, 1], whose integral diverges. The page shows
// the parts between the markers and the output in the raw strings; main returns non-zero unless the program prints that
// output and what the page says holds. The integrands use only + - * /, which IEEE 754 rounds correctly, and the
// methods only arithmetic, so the output is the same on every target that does not contract a * b + c into a fused
// multiply-add (x86-64 without -mfma, WebAssembly).

// clang-format off
// <includes>
#include <numerixx/quadrature.hpp>

#include <cmath>
#include <numbers>
#include <print>
// </includes>
// clang-format on

#include "printed.hpp"

#include <string>

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    // clang-format off
    // <first>
    namespace qd = nxx::quadrature;
    int evaluations = 0;
    auto const fn = [&evaluations](double x) {
        ++evaluations;
        return 4.0 / (1.0 + x * x);
    };
    double const pi = std::numbers::pi; // the exact integral
    auto const r = qd::integrate(fn, 0.0, 1.0);
    if (r) {
        std::println("value     = {:.17g}", r->value);
        std::println("abs_error = {:.1e}", r->abs_error);
        std::println("error     = {:.1e}", std::abs(r->value - pi));
        std::println("{} evaluations", evaluations);
    }
    // </first>

    // <methods>
    evaluations = 0;
    auto const g = qd::integrate<qd::GaussLegendre<10>>(fn, 0.0, 1.0);
    if (g) std::println("GaussLegendre<10>  error {:.1e}  {:3} evaluations", std::abs(*g - pi), evaluations);
    evaluations = 0;
    auto const t = qd::integrate<qd::TanhSinh>(fn, 0.0, 1.0);
    if (t) std::println("TanhSinh           error {:.1e}  {:3} evaluations", std::abs(t->value - pi), evaluations);
    evaluations = 0;
    auto const m = qd::integrate<qd::Romberg>(fn, 0.0, 1.0);
    if (m) std::println("Romberg            error {:.1e}  {:3} evaluations", std::abs(m->value - pi), evaluations);
    // </methods>

    // <failure>
    auto const d = qd::integrate([](double x) { return 1.0 / x; }, 0.0, 1.0);
    if (!d) std::println("no integral: {}", d.error());
    // </failure>
    // clang-format on

    constexpr std::string_view first_output = R"(
value     = 3.1415926535897931
abs_error = 3.5e-14
error     = 0.0e+00
150 evaluations
)";
    constexpr std::string_view methods_output = R"(
GaussLegendre<10>  error 2.5e-13   10 evaluations
TanhSinh           error 4.4e-16  121 evaluations
Romberg            error 4.8e-11   33 evaluations
)";
    constexpr std::string_view failure_output = R"(
no integral: not_converged
)";
    std::string const expected = std::string(docs::quoted(first_output)) + std::string(docs::quoted(methods_output))
                                 + std::string(docs::quoted(failure_output));
    if (!(docs::printed() == expected)) return 1;

    // Each check states what must hold, so that a NaN, which fails every comparison, fails it.
    if (!r || !g || !t || !m || d) return 2;
    // "The estimate is larger than the actual error", for each method with one: the honesty the page describes.
    if (!(std::abs(r->value - pi) <= r->abs_error && std::abs(t->value - pi) <= t->abs_error
          && std::abs(m->value - pi) <= m->abs_error))
        return 3;
    // "Romberg's estimate overstates its error": by more than a hundred times here (1.2e-8 against 4.8e-11).
    if (!(m->abs_error > 100.0 * std::abs(m->value - pi))) return 4;
    // "Every request was met": the default request, sqrt(eps) relative (1.5e-8), for the three with an estimate.
    double const request = 1.5e-8 * pi;
    if (!(r->abs_error <= request && t->abs_error <= request && m->abs_error <= request)) return 5;
    // "Romberg's estimate, 1.2e-8"; "abs_error ... about 1e-14 relative", the default's floor 50 eps L1 with L1 = pi.
    if (!(m->abs_error >= 1.15e-8 && m->abs_error < 1.25e-8)) return 7;
    if (!(r->abs_error / pi > 5e-15 && r->abs_error / pi < 2e-14)) return 8;
    // "Kept cutting panels near 0 until its estimate stopped improving": the mesh's 150 evaluations, then 32
    // bisections of 30 evaluations each without a new minimum of the excess, Gauss-Kronrod's stagnation rule.
    int divergent = 0;
    auto const counted = qd::integrate([&divergent](double x) { ++divergent; return 1.0 / x; }, 0.0, 1.0);
    if (!(counted == d && divergent == 150 + 32 * 30)) return 9;
    // The divergent integral: not_converged, and the default's best, which the page calls meaningless, is a finite
    // number (31.5 here) where the integral is infinite.
    if (!(d.error().reason == qd::error::not_converged && d.error().best && std::isfinite(d.error().best->value)))
        return 6;
    return 0;
}

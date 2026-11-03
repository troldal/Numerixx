// The program of docs/sphinx/quadrature/choose-a-method.md: the cost of two methods on 1 / sqrt x, measured by
// counting evaluations in fn; an upper limit found by Newton's method on an integral; and a double integral by nesting.
// The page shows the parts between the markers and the output in the raw strings; main returns non-zero unless the
// program prints that output and what the page says holds. The integrands use only + - * / and std::sqrt, which IEEE
// 754 rounds correctly, and the methods only arithmetic, so the output is the same on every target that does not
// contract a * b + c into a fused multiply-add (x86-64 without -mfma, WebAssembly).

// clang-format off
// <includes>
#include <numerixx/quadrature.hpp>
#include <numerixx/roots.hpp>

#include <cmath>
#include <limits>
#include <print>
#include <utility>
// </includes>
// clang-format on

#include "printed.hpp"

#include <string>

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    // clang-format off
    // <cost>
    namespace qd = nxx::quadrature;
    int evaluations = 0;
    auto const fn = [&evaluations](double x) {
        ++evaluations;
        return 1.0 / std::sqrt(x);
    };
    auto const k = qd::integrate(fn, 0.0, 1.0);
    if (k) std::println("GaussKronrod<15>  error {:.1e}  {} evaluations", std::abs(k->value - 2.0), evaluations);
    evaluations = 0;
    auto const t = qd::integrate<qd::TanhSinh>(fn, 0.0, 1.0);
    if (t) std::println("TanhSinh          error {:.1e}  {} evaluations", std::abs(t->value - 2.0), evaluations);
    // </cost>

    // <limit>
    namespace rt = nxx::roots;
    double const nan = std::numeric_limits<double>::quiet_NaN();
    // F(x), the integral of sqrt t from 0 to x, or NaN where there is none
    auto const F = [nan](double x, double rel_tol) {
        auto const i = qd::integrate([](double s) { return std::sqrt(s); }, 0.0, x, {.rel_tol = rel_tol});
        return i ? i->value : nan;
    };
    for (double const rel_tol : {1.5e-8, 1e-14}) {
        // F(x) = 2/3 and F'(x) = sqrt x, as a pair for Newton's method
        auto const fdf = [&F, rel_tol](double x) { return std::pair{F(x, rel_tol) - 2.0 / 3.0, std::sqrt(x)}; };
        auto const r = rt::solve<rt::Newton>(fdf, 2.0);
        if (r) std::println("rel_tol {:.1e}: x = {:.17g}, abs_error {:.1e}, error {:.1e}", rel_tol, r->value,
                            r->abs_error, std::abs(r->value - 1.0));
    }
    // </limit>

    // <nested>
    long count = 0;
    // the inner integral over x, for a given y, with a tighter tolerance than the outer one
    auto const inner = [&count, nan](double y) {
        auto const i = qd::integrate([&count, y](double x) {
            ++count;
            return 1.0 / (1.0 + x + y);
        }, 0.0, 1.0, {.rel_tol = 1e-12});
        return i ? i->value : nan;
    };
    auto const outer = qd::integrate(inner, 0.0, 1.0);
    double const exact = 0.5232481437645478365168072249348708416404711360; // 3 ln 3 - 4 ln 2
    if (outer)
        std::println("{:.16f}, error {:.1e}, {} evaluations", outer->value, std::abs(outer->value - exact), count);
    // </nested>
    // clang-format on

    constexpr std::string_view cost_output = R"(
GaussKronrod<15>  error 1.4e-08  1350 evaluations
TanhSinh          error 0.0e+00  131 evaluations
)";
    constexpr std::string_view limit_output = R"(
rel_tol 1.5e-08: x = 0.99999999970610376, abs_error 4.4e-16, error 2.9e-10
rel_tol 1.0e-14: x = 0.99999999999999978, abs_error 4.4e-16, error 2.2e-16
)";
    constexpr std::string_view nested_output = R"(
0.5232481437645479, error 0.0e+00, 22500 evaluations
)";
    std::string const expected = std::string(docs::quoted(cost_output)) + std::string(docs::quoted(limit_output))
                                 + std::string(docs::quoted(nested_output));
    if (!(docs::printed() == expected)) return 1;

    // Each check states what must hold, so that a NaN, which fails every comparison, fails it.
    if (!k || !t || !outer) return 2;
    // "Ten times fewer evaluations", and "both meet the request", sqrt(eps) relative = 1.5e-8 of 2 (the output above:
    // 1,350 against 131).
    if (!(std::abs(k->value - 2.0) <= 3e-8 && std::abs(t->value - 2.0) <= 3e-8)) return 3;
    // "The default's estimate covers its error" (1.4e-8 against 2.1e-8).
    if (!(std::abs(k->value - 2.0) <= k->abs_error)) return 4;
    // The root's abs_error is about F as computed: with the default request, F(1) is off by 2.9e-10 (its estimate
    // 5.1e-9), and the root, where F' = 1, by the same; the root's abs_error, 4.4e-16, does not see it.
    auto const f1 = qd::integrate([](double s) { return std::sqrt(s); }, 0.0, 1.0);
    if (!(f1 && std::abs(f1->value - 2.0 / 3.0) > 1e-10 && std::abs(f1->value - 2.0 / 3.0) <= f1->abs_error)) return 5;
    // "Each solve took six integrals" (a Newton step costs one value of fdf; six from the guess 2).
    int integrals = 0;
    auto const counted = [&](double x) {
        ++integrals;
        return std::pair{F(x, 1e-14) - 2.0 / 3.0, std::sqrt(x)};
    };
    auto const again = rt::solve<rt::Newton>(counted, 2.0);
    if (!(again && integrals == 6)) return 6;
    // A failed inner integral reaches the outer method as NaN, which it reports as non_finite: an inner integrand
    // that is NaN for y > 0.5.
    auto const broken = qd::integrate([nan](double y) {
        auto const i = qd::integrate([y, nan](double x) { return y > 0.5 ? nan : x + y; }, 0.0, 1.0);
        return i ? i->value : nan;
    }, 0.0, 1.0);
    if (!(!broken && broken.error().reason == qd::error::non_finite)) return 7;
    // "An isolated failure can be stepped around" (the independent review): the inner integral fails at the one outer
    // point y = 0.35, a node of the default's mesh, and the outer integral succeeds, 1, ignoring it; counting the
    // failures in F shows it.
    int failures = 0;
    auto const isolated = qd::integrate([nan, &failures](double y) {
        auto const i = qd::integrate([y, nan](double x) { return y == 0.35 ? nan : 1.0 + 0.0 * x; }, 0.0, 1.0);
        if (!i) ++failures;
        return i ? i->value : nan;
    }, 0.0, 1.0);
    if (!(isolated && isolated->value == 1.0 && failures == 1)) return 8;
    // The nested integral: 150 inner integrals of 150 evaluations each, both on their initial mesh; within the
    // outer request, and the outer estimate covers the error.
    if (!(count == 150L * 150L && std::abs(outer->value - exact) <= outer->abs_error)) return 8;
    return 0;
}

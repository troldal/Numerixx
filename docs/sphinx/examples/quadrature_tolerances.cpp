// The program of docs/sphinx/quadrature/tolerances.md: more digits and the rounding floor on a peak, abs_tol for an
// integral that is 0, a result whose abs_error exceeds the request, and a noisy function. The page shows the parts
// between the markers and the output in the raw strings; main returns non-zero unless the program prints that output
// and what the page says holds. The integrands use only + - * / and integer operations, which are exact or correctly
// rounded, and the methods only arithmetic, so the output is the same on every target that does not contract a * b + c
// into a fused multiply-add (x86-64 without -mfma, WebAssembly).

// clang-format off
// <includes>
#include <numerixx/quadrature.hpp>

#include <cmath>
#include <initializer_list>
#include <print>
// </includes>
// clang-format on

#include "printed.hpp"

#include <bit>
#include <cstdint>
#include <string>

namespace {

    namespace qd = nxx::quadrature;

    // Noise in [-1, 1] from a hash of x's bits, as in the tests: it never resolves, however close the points.
    double hash(double x)
    {
        std::uint64_t h = std::bit_cast<std::uint64_t>(x) * 0x9E3779B97F4A7C15ull;
        h ^= h >> 29;
        h *= 0xBF58476D1CE4E5B9ull;
        h ^= h >> 32;
        return static_cast<double>(h >> 11) / 4503599627370496.0 - 1.0;
    }

    // A braced value without a name is the first member, abs_tol.
    static_assert(qd::Criteria<double>{ 1e-12 }.abs_tol == 1e-12);
    static_assert(qd::Criteria<double>{ .rel_tol = 1e-12 }.abs_tol == 0.0);

} // namespace

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    // 100 (atan 70 + atan 30), the closed form, from mpmath at 60 digits.
    double const exact = 309.39869151241494108699839806789614552002783882631;

    // clang-format off
    // <more>
    int evaluations = 0;
    auto const peak = [&evaluations](double x) {
        ++evaluations;
        return 1.0 / ((x - 0.3) * (x - 0.3) + 1e-4);
    };
    for (double const rel_tol : {1.5e-8, 1e-10, 1e-12, 1e-14}) {
        evaluations = 0;
        auto const r = qd::integrate(peak, 0.0, 1.0, {.rel_tol = rel_tol});
        if (r) std::println("rel_tol {:<7g}: {:3} evaluations, abs_error {:.1e}, error {:.1e}", rel_tol, evaluations,
                            r->abs_error, std::abs(r->value - exact));
    }
    // </more>

    // <variable>
    qd::Criteria<double> const criteria{.rel_tol = 1e-12};
    auto const result = qd::integrate(peak, 0.0, 1.0, criteria);
    // </variable>

    // <zero>
    auto const odd = [&evaluations](double x) {
        ++evaluations;
        return 1.0 / ((x - 0.3) * (x - 0.3) + 1e-4) - 1.0 / ((x + 0.3) * (x + 0.3) + 1e-4);
    };
    for (double const abs_tol : {0.0, 1e-8, 1e-6}) {
        evaluations = 0;
        auto const r = qd::integrate(odd, -1.0, 1.0, {.abs_tol = abs_tol});
        if (r) std::println("abs_tol {:<5g}: {:4} evaluations, value {:.1e}, abs_error {:.1e}", abs_tol, evaluations,
                            r->value, r->abs_error);
    }
    // </zero>

    // <sliver>
    double const b = 1.0 + 1e-9;
    auto const s = qd::integrate<qd::TanhSinh>([](double) { return 1.0; }, 1.0, b);
    if (s) std::println("abs_error {:.1e}, relative {:.1e}", s->abs_error, s->abs_error / s->value);
    // </sliver>

    // <noisy>
    auto const noisy = [](double x) { return 1.0 / (1.0 + x) * (1.0 + 1e-10 * hash(x)); };
    double const ln2 = 0.69314718055994530941723212145817656807550013436026;
    auto const n = qd::integrate<qd::TanhSinh>(noisy, 0.0, 1.0, {.rel_tol = 1e-12});
    if (n) std::println("abs_error {:.1e}, error {:.1e}", n->abs_error, std::abs(n->value - ln2));
    // </noisy>
    // clang-format on

    constexpr std::string_view more_output = R"(
rel_tol 1.5e-08: 330 evaluations, abs_error 5.4e-07, error 2.3e-13
rel_tol 1e-10  : 390 evaluations, abs_error 8.7e-09, error 1.7e-13
rel_tol 1e-12  : 630 evaluations, abs_error 1.3e-10, error 1.7e-13
rel_tol 1e-14  : 690 evaluations, abs_error 6.6e-12, error 5.7e-14
)";
    constexpr std::string_view zero_output = R"(
abs_tol 0    : 1470 evaluations, value -5.6e-14, abs_error 1.3e-11
abs_tol 1e-08:  960 evaluations, value -5.2e-13, abs_error 9.2e-09
abs_tol 1e-06:  720 evaluations, value 3.0e-13, abs_error 8.1e-07
)";
    constexpr std::string_view sliver_output = R"(
abs_error 2.0e-16, relative 2.0e-07
)";
    constexpr std::string_view noisy_output = R"(
abs_error 2.1e-13, error 2.3e-12
)";
    std::string const expected = std::string(docs::quoted(more_output)) + std::string(docs::quoted(zero_output))
                                 + std::string(docs::quoted(sliver_output)) + std::string(docs::quoted(noisy_output));
    if (!(docs::printed() == expected)) return 1;

    // Each check states what must hold, so that a NaN, which fails every comparison, fails it.
    double const eps = std::numeric_limits<double>::epsilon();
    // "The default request is sqrt(eps), 1.5e-8 in double", and the defaults of Criteria.
    if (!(qd::Criteria<double>{}.rel_tol == std::sqrt(eps) && qd::Criteria<double>{}.abs_tol == 0.0)) return 2;
    // The peak. "Every result is far more accurate than requested" (the estimate overstates); "a request below the
    // floor ends at the floor": 50 eps L1 with L1 the integral, 309.4, is 3.4e-12, and the default and TanhSinh accept
    // at twice it, 6.9e-12; so rel_tol 1e-14 and 0 give the same result.
    auto const at = [&](double rel_tol) { return qd::integrate(peak, 0.0, 1.0, {.rel_tol = rel_tol}); };
    for (double const rel_tol : {1.5e-8, 1e-10, 1e-12, 1e-14}) {
        auto const r = at(rel_tol);
        if (!(r && std::abs(r->value - exact) <= r->abs_error)) return 3;
    }
    // "The default request already gave 15 digits here": relative error below 1e-15 (7.3e-16).
    auto const fifteen = at(1.5e-8);
    if (!(std::abs(fifteen->value - exact) / exact < 1e-15)) return 11;
    // The variable holds the same request as the braced list.
    if (!(result && result == at(1e-12))) return 13;
    // "A request of 0 gives the same result as 1e-14": both below the floor. (1e-13 gives it too, since the estimate
    // fell from 1.3e-10 to 6.6e-12 in two bisections; quadrature_floor.cpp asserts that for the section's figure.)
    auto const floor14 = at(1e-14);
    if (!(floor14 == at(0.0))) return 4;
    double const floor = 50.0 * eps * exact;
    // "Just under twice the floor, here 3.4e-12".
    if (!(floor14->abs_error > 1.8 * floor && floor14->abs_error <= 2.0 * floor && floor > 3.3e-12 && floor < 3.5e-12))
        return 5;
    // "About 1e-14 relative": 6.6e-12 / 309.4 = 2.1e-14.
    if (!(floor14->abs_error / exact > 5e-15 && floor14->abs_error / exact < 5e-14)) return 6;
    // The integral that is 0 (odd about 0): "ends at the floor without abs_tol"; each value is within its abs_error of
    // the exact 0.
    for (double const abs_tol : {0.0, 1e-8, 1e-6}) {
        auto const r = qd::integrate(odd, -1.0, 1.0, {.abs_tol = abs_tol});
        if (!(r && std::abs(r->value) <= r->abs_error)) return 7;
    }
    // The sliver: "the request was 1.5e-8 relative"; "the value is exact": the integral of 1 is b - 1, exact in double.
    if (!(s && s->abs_error / s->value > 1.5e-8 && s->value == b - 1.0)) return 8;
    // The noisy function: "ok, with an abs_error about a tenth of its actual error" (2.1e-13 against 2.3e-12), and at
    // the default request "the same function is fine": error far below 1.5e-8 ln 2.
    if (!(n && n->abs_error < 0.2 * std::abs(n->value - ln2))) return 9;
    auto const d = qd::integrate<qd::TanhSinh>(noisy, 0.0, 1.0);
    if (!(d && std::abs(d->value - ln2) <= d->abs_error && std::abs(d->value - ln2) <= 1.5e-8 * ln2)) return 10;
    // "On this function the recommended method notices the noise and gives up instead" at 1e-12.
    auto const k = qd::integrate(noisy, 0.0, 1.0, {.rel_tol = 1e-12});
    if (!(!k && k.error().reason == qd::error::not_converged)) return 12;
    return 0;
}

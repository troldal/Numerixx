// The data of the figure in docs/sphinx/quadrature/tolerances.md ("More digits"): the default method on the peak
// 1 / ((x - 0.3)^2 + 1e-4) over [0, 1] at every relative tolerance from 1e-2 to 1e-16, with its estimate and its actual
// error, relative to the integral, against the request and the rounding floor. The program prints them as CSV and
// compares them with the table in the raw string, which docs/sphinx/plots/quadrature_floor.py plots; main returns
// non-zero unless the two agree and what the page says of them holds (PHASE8_PLAN.md). The integrand and the method
// use only + - * /, which IEEE 754 rounds correctly, so the table is the same on every target that does not contract
// a * b + c into a fused multiply-add (x86-64 without -mfma, WebAssembly).

#include <numerixx/quadrature.hpp>

#include "printed.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <print>
#include <string_view>
#include <vector>

namespace {

    namespace qd = nxx::quadrature;

    double peak(double x) { return 1.0 / ((x - 0.3) * (x - 0.3) + 1e-4); }

    struct Run {
        double rel_tol;
        std::size_t evaluations;
        double value;
        double abs_error;
    };

} // namespace

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    // 100 (atan 70 + atan 30), the closed form, from mpmath at 60 digits (as in quadrature_tolerances.cpp).
    double const exact = 309.39869151241494108699839806789614552002783882631;
    double const eps = std::numeric_limits<double>::epsilon();
    // The rounding floor, 50 eps times the integral of |f|, which is the integral itself (f > 0); relative to it,
    // 50 eps.
    double const floor = 50.0 * eps;

    std::vector<Run> runs;
    for (double const rel_tol : std::array{1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13,
                                           1e-14, 1e-15, 1e-16}) {
        std::size_t evaluations = 0;
        auto const r = qd::integrate([&evaluations](double x) {
            ++evaluations;
            return peak(x);
        }, 0.0, 1.0, {.rel_tol = rel_tol});
        if (!r) return 2;
        runs.push_back({.rel_tol = rel_tol, .evaluations = evaluations, .value = r->value, .abs_error = r->abs_error});
    }

    // Relative to the integral: the request, the estimate, the actual error, the floor.
    std::println("rel_tol,evaluations,estimate,error,floor");
    for (Run const &r : runs)
        std::println("{:.0e},{},{:.3g},{:.3g},{:.3g}", r.rel_tol, r.evaluations, r.abs_error / exact,
                     std::abs(r.value - exact) / exact, floor);

    constexpr std::string_view floor_table = R"(
rel_tol,evaluations,estimate,error,floor
1e-02,150,0.00116,4.52e-06,1.11e-14
1e-03,180,0.000615,2.27e-06,1.11e-14
1e-04,210,6.58e-05,1.42e-08,1.11e-14
1e-05,270,4.9e-07,3.06e-12,1.11e-14
1e-06,270,4.9e-07,3.06e-12,1.11e-14
1e-07,330,1.73e-09,7.35e-16,1.11e-14
1e-08,330,1.73e-09,7.35e-16,1.11e-14
1e-09,360,8.79e-10,0,1.11e-14
1e-10,390,2.8e-11,5.51e-16,1.11e-14
1e-11,510,7.73e-12,5.51e-16,1.11e-14
1e-12,630,4.29e-13,5.51e-16,1.11e-14
1e-13,690,2.13e-14,1.84e-16,1.11e-14
1e-14,690,2.13e-14,1.84e-16,1.11e-14
1e-15,690,2.13e-14,1.84e-16,1.11e-14
1e-16,690,2.13e-14,1.84e-16,1.11e-14
)";
    if (!(docs::printed() == docs::quoted(floor_table))) return 1;

    // What the page says.
    // From 1e-2 to 1e-12 every estimate meets its request, and the evaluations grow from 150 to 630, never falling.
    for (std::size_t i = 0; i <= 10; ++i)
        if (!(runs[i].abs_error <= runs[i].rel_tol * std::abs(runs[i].value))) return 3;
    if (!(runs[0].evaluations == 150 && runs[10].evaluations == 630)) return 4;
    for (std::size_t i = 1; i < runs.size(); ++i)
        if (!(runs[i].evaluations >= runs[i - 1].evaluations)) return 5;
    // Two more bisections, 60 evaluations, take the estimate from 4.3e-13 (at 1e-12) to 2.1e-14; from 1e-13 on, the
    // same run whatever the request: 690 evaluations and the same estimate, just under twice the floor.
    if (!(runs[11].evaluations - runs[10].evaluations == 60)) return 12;
    if (!(std::abs(runs[10].abs_error / exact - 4.3e-13) < 0.05e-13 &&
          std::abs(runs[11].abs_error / exact - 2.1e-14) < 0.05e-14))
        return 13;
    for (std::size_t i = 11; i < runs.size(); ++i)
        if (!(runs[i].evaluations == 690 && runs[i].abs_error == runs[11].abs_error)) return 6;
    // At 1e-13 that estimate still meets the request; from 1e-14 on it is above it, and the method accepts it as
    // within twice the floor, the default's attainable error.
    if (!(runs[11].abs_error <= runs[11].rel_tol * std::abs(runs[11].value))) return 7;
    for (std::size_t i = 12; i < runs.size(); ++i)
        if (!(runs[i].abs_error > runs[i].rel_tol * std::abs(runs[i].value))) return 14;
    double const at_floor = runs[11].abs_error / exact;
    if (!(at_floor > 1.8 * floor && at_floor <= 2.0 * floor)) return 8;
    // At every request the actual error is at least a hundred times smaller than the estimate.
    for (Run const &r : runs)
        if (!(100.0 * std::abs(r.value - exact) < r.abs_error)) return 9;
    // From 1e-7 on it is within four units in the last place of the integral, 2^-44 for a number in [256, 512); at
    // 1e-9 it is 0, the double nearest the integral.
    double const ulp = std::ldexp(1.0, -44);
    for (std::size_t i = 5; i < runs.size(); ++i)
        if (!(std::abs(runs[i].value - exact) <= 4.0 * ulp)) return 10;
    if (!(runs[7].value == exact)) return 11;
    return 0;
}

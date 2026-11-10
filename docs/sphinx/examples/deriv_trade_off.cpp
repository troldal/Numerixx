// The data of the figure in docs/sphinx/deriv/index.md ("The idea, and its catch"): the relative error of the central
// difference of sqrt x at 2 for the steps 1, 0.5, 0.2, 0.1, ... down to 1e-13, with the default step and its error,
// and the error and estimate of Ridders' method, against the closed form f'(2) = 1 / (2 sqrt 2). The program prints
// them as CSV and compares them with the table in the raw string, which docs/sphinx/plots/deriv_trade_off.py plots;
// main returns non-zero unless the two agree and what the page says of them holds (PHASE9_PLAN.md). Each step is the
// correctly rounded quotient of two integers; fn is std::sqrt, which IEEE 754 rounds correctly, and the methods use
// only + - * /, so the table is the same on every target that does not contract a * b + c into a fused multiply-add
// (x86-64 without -mfma, WebAssembly).
//
// Decimal steps rather than powers of two: with h = 2^-k, 2 + h and 2 - h are exact and the computed difference is a
// multiple of 2^(k-53), so one lucky multiple can repeat over several halvings; the curve then shows a flat stretch
// (the same result at six of the seven steps 2^-16 to 2^-22) that is an artefact of those steps, not the rounding error
// a reader meets.

#include <numerixx/deriv.hpp>

#include "printed.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <print>
#include <string_view>
#include <vector>

namespace dv = nxx::deriv;

namespace {

    struct Point {
        double h;
        double error; // relative to f'(2)
    };

    // The steps 1, 0.5, 0.2, 0.1, 0.05, 0.02, ... down to 10^-last, each the double nearest its decimal value: 10^j is
    // exact for j <= 22, and a quotient of two exact integers is correctly rounded.
    std::vector<double> steps(int last)
    {
        std::vector<double> hs;
        double p = 1.0; // 10^j
        for (int j = 0; j < last; ++j, p *= 10.0) {
            hs.push_back(1.0 / p);
            hs.push_back(5.0 / (10.0 * p));
            hs.push_back(2.0 / (10.0 * p));
        }
        hs.push_back(1.0 / p);
        return hs;
    }

    // The point of the table at the step nearest h.
    Point const &at(std::vector<Point> const &curve, double h)
    {
        return *std::ranges::min_element(curve, {}, [h](Point const &p) { return std::abs(std::log(p.h / h)); });
    }

} // namespace

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    auto const fn = [](double x) { return std::sqrt(x); };
    // 1 / (2 sqrt 2) = sqrt(2) / 4: the correctly rounded sqrt(2), divided exactly by 4, is within half an ulp of the
    // closed form, a relative 1.1e-16, below every error the program compares.
    double const exact = std::sqrt(2.0) / 4.0;
    auto const relative = [exact](double value) { return std::abs(value - exact) / exact; };

    std::vector<Point> curve;
    for (double const h : steps(13)) {
        auto const d = dv::derivative(fn, 2.0, h);
        if (!d) return 2;
        curve.push_back({.h = h, .error = relative(d->value)});
    }
    double const step = dv::default_step<dv::Central<2>>(2.0);
    auto const central = dv::derivative(fn, 2.0);
    auto const ridders = dv::derivative<dv::Ridders>(fn, 2.0);
    if (!central || !ridders || !ridders->abs_error) return 3;
    double const central_error = relative(central->value);
    double const ridders_error = relative(ridders->value);
    double const ridders_estimate = *ridders->abs_error / exact;

    std::println("h,error");
    for (Point const &p : curve) std::println("{:.6g},{:.3g}", p.h, p.error);
    std::println("");
    std::println("default_step,default_error,ridders_error,ridders_estimate");
    std::println("{:.3g},{:.3g},{:.3g},{:.3g}", step, central_error, ridders_error, ridders_estimate);

    constexpr std::string_view trade_off_table = R"(
h,error
1,0.0353
0.5,0.00803
0.2,0.00126
0.1,0.000313
0.05,7.81e-05
0.02,1.25e-05
0.01,3.13e-06
0.005,7.81e-07
0.002,1.25e-07
0.001,3.13e-08
0.0005,7.81e-09
0.0002,1.25e-09
0.0001,3.12e-10
5e-05,8.08e-11
2e-05,1.07e-11
1e-05,5.1e-12
5e-06,4.09e-12
2e-06,5.58e-11
1e-06,1.11e-10
5e-07,2.22e-10
2e-07,4.25e-10
1e-07,1.44e-09
5e-08,5.9e-10
2e-08,4.84e-09
1e-08,2.61e-08
5e-09,5.31e-09
2e-09,2.42e-08
1e-09,3.13e-08
5e-10,3.13e-08
2e-10,3.45e-07
1e-10,1.22e-06
5e-11,4.36e-06
2e-11,7.51e-06
1e-11,8.2e-06
5e-12,2.32e-05
2e-12,4.15e-05
1e-12,6e-05
5e-13,0.000374
2e-13,0.000393
1e-13,0.00252

default_step,default_error,ridders_error,ridders_estimate
1.21e-05,6.04e-12,5.84e-14,2.27e-13
)";
    if (!(docs::printed() == docs::quoted(trade_off_table))) return 1;

    // What the page says. The error model: the truncation error, relative to f'(2), is h^2 / 32 (1 + 0.109 h^2 + ...)
    // (h^2 f''' / 6 + h^4 f^(5) / 120 over f', with f''' / f' = 3/16 and f^(5) / f' = 105/256 at 2); the rounding
    // error is at most rho / h with rho = 3.2e-16: each sqrt is within half an ulp, 2^-53, of sqrt(2 +- h), their
    // difference within 2^-52, divided by the actual distance of the points, which is within 3.3e-16 of 2h and exact
    // (Sterbenz), and by f'(2) = 0.354. That the points round moves their midpoint by at most 1.7e-16, which shifts
    // the result by |f''| 1.7e-16, a relative 4e-17: negligible.
    double const rho = 3.2e-16;

    // "The error falls fourfold when the step halves": at each halving from 0.2 down to 1e-4 the ratio is 4 within
    // 2.5%. The truncation's next term changes the ratio by a factor 1 + 0.082 h^2, at most 1.004 here; the rounding,
    // at most rho / h, is at most 1% of the truncation at 1e-4 (3.2e-12 against 3.1e-10), so 2% on the ratio.
    for (double const h : {0.2, 0.02, 0.002, 0.0002}) {
        double const ratio = at(curve, h).error / at(curve, h / 2.0).error;
        if (!(ratio > 0.975 * 4.0 && ratio < 1.025 * 4.0)) return 4;
    }
    // It falls at every step from 1 down to the smallest error of the curve, which is within a factor 8 of the default
    // step, 1.2e-5 (eps^(1/3) times the scale 2); the default step's own error is within ten times that smallest one.
    auto const best = std::ranges::min_element(curve, {}, &Point::error);
    for (auto p = curve.begin(); p != best; ++p)
        if (!(std::next(p)->error < p->error)) return 5;
    if (!(best->h < 8.0 * step && best->h > step / 8.0)) return 6;
    if (!(central_error <= 10.0 * best->error)) return 7;
    if (!(step > 1.2e-5 && step < 1.22e-5)) return 8;
    // "About two thirds of the type's digits": the smallest error keeps 10 to 13 of a double's 16 (the model,
    // eps^(2/3), is 3.7e-11).
    if (!(best->error < 1e-10 && best->error > 1e-13)) return 14;
    // Then rounding takes over: below 1e-7 every error is within the rounding bound rho / h (with the truncation, at
    // most 3.2e-16 there, and 1.1e-16 for exact), and every error below 1e-9 is more than a thousand times the
    // smallest.
    for (Point const &p : curve) {
        if (p.h < 1e-7 && !(p.error <= rho / p.h + p.h * p.h / 32.0 + 1.1e-16 + 1.1e-16)) return 9;
        if (p.h < 1e-9 && !(p.error > 1000.0 * best->error)) return 10;
    }
    // "Two or three digits are left" at 1e-13: an error between 1e-3 and 1e-2.
    if (!(curve.back().error > 1e-3 && curve.back().error < 1e-2)) return 11;
    // Ridders' method: its estimate covers its error, which is 50 to 500 times below the curve's smallest error; the
    // estimate is 1 to 10 times the error.
    if (!(ridders_error <= ridders_estimate && ridders_estimate <= 10.0 * ridders_error)) return 12;
    if (!(best->error > 50.0 * ridders_error && best->error < 500.0 * ridders_error)) return 13;
    return 0;
}

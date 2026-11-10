// The data of the figure in docs/sphinx/deriv/step.md ("When the default step is too large"): the relative error of the
// central difference for the steps 1, 0.5, 0.2, 0.1, ... down to 1e-14 on two functions, sqrt x at 2, which varies on
// the scale max(|x|, 1), and the peak x / (x^2 + a^2), a = 1e-4, at 2e-4, which varies on the scale a, with each
// function's default step and its error; against the closed forms 1 / (2 sqrt 2) and (a^2 - x^2) / (x^2 + a^2)^2. The
// program prints them as CSV and compares them with the table in the raw string, which docs/sphinx/plots/deriv_scale.py
// plots; main returns non-zero unless the two agree and what the page says of them holds (PHASE9_PLAN.md). The
// functions use only + - * / and std::sqrt, which IEEE 754 rounds correctly, and so do the methods and the default
// step (nxx::math::root_epsilon), so the table is the same on every target that does not contract a * b + c into a
// fused multiply-add (x86-64 without -mfma, WebAssembly).
//
// The steps are decimal, as in deriv_trade_off.cpp and for the same reason: with powers of two, one lucky multiple of
// the computed difference's quantum repeats over several halvings, here an error of 3e-16 on the peak at 2^-31, 2^-34
// and 2^-45, far below the rounding error at any other step, which says nothing about the choice of a step.

#include <numerixx/deriv.hpp>

#include "printed.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <print>
#include <span>
#include <string_view>
#include <vector>

namespace dv = nxx::deriv;

namespace {

    constexpr double a = 1e-4; // the peak's width
    constexpr double peak_x = 2e-4;

    double root(double x) { return std::sqrt(x); }
    double peak(double x) { return x / (x * x + a * a); }

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

    struct Curve {
        std::vector<double> errors; // relative, one per step
        double default_step;
        double default_error;
    };

    // The central difference of fn at x for every step, and with the default step, relative to exact; nothing if a
    // call fails.
    bool trace(double (*fn)(double), double x, double exact, std::span<double const> hs, Curve &c)
    {
        auto const relative = [exact](double value) { return std::abs(value - exact) / std::abs(exact); };
        for (double const h : hs) {
            auto const d = dv::derivative(fn, x, h);
            if (!d) return false;
            c.errors.push_back(relative(d->value));
        }
        auto const d = dv::derivative(fn, x);
        if (!d) return false;
        c.default_step = dv::default_step<dv::Central<2>>(x);
        c.default_error = relative(d->value);
        return true;
    }

    // The index of a curve's smallest error, the first if several.
    std::size_t best(Curve const &c)
    { return static_cast<std::size_t>(std::ranges::min_element(c.errors) - c.errors.begin()); }

} // namespace

int main(int, char **argv)
{
    if (!docs::capture(argv[0])) return 100;

    // The references. sqrt(2) / 4 is within half an ulp of 1 / (2 sqrt 2), a relative 1.1e-16. For the peak, -1.2e7:
    // a^2 and x^2 are each within half an ulp, their difference -3 a^2 and the squared sum 25 a^4 then within about 3
    // ulps, a relative 1e-15 in all; both far below the smallest error compared, 4e-12.
    double const root_exact = std::sqrt(2.0) / 4.0;
    double const peak_exact = (a * a - peak_x * peak_x) / ((peak_x * peak_x + a * a) * (peak_x * peak_x + a * a));
    std::vector<double> const hs = steps(14);
    Curve slow;
    Curve fast;
    if (!trace(root, 2.0, root_exact, hs, slow) || !trace(peak, peak_x, peak_exact, hs, fast)) return 2;

    std::println("h,sqrt,peak");
    for (std::size_t i = 0; i < hs.size(); ++i)
        std::println("{:.6g},{:.3g},{:.3g}", hs[i], slow.errors[i], fast.errors[i]);
    std::println("");
    std::println("function,default_step,default_error");
    std::println("sqrt,{:.3g},{:.3g}", slow.default_step, slow.default_error);
    std::println("peak,{:.3g},{:.3g}", fast.default_step, fast.default_error);

    constexpr std::string_view scale_table = R"(
h,sqrt,peak
1,0.0353,1
0.5,0.00803,1
0.2,0.00126,1
0.1,0.000313,1
0.05,7.81e-05,1
0.02,1.25e-05,1
0.01,3.13e-06,1
0.005,7.81e-07,1
0.002,1.25e-07,1.02
0.001,3.13e-08,1.09
0.0005,7.81e-09,1.37
0.0002,1.25e-09,1.49
0.0001,3.12e-10,0.167
5e-05,8.08e-11,0.0274
2e-05,1.07e-11,0.00383
1e-05,5.1e-12,0.00094
5e-06,4.09e-12,0.000234
2e-06,5.58e-11,3.73e-05
1e-06,1.11e-10,9.33e-06
5e-07,2.22e-10,2.33e-06
2e-07,4.25e-10,3.73e-07
1e-07,1.44e-09,9.33e-08
5e-08,5.9e-10,2.33e-08
2e-08,4.84e-09,3.73e-09
1e-08,2.61e-08,9.33e-10
5e-09,5.31e-09,2.28e-10
2e-09,2.42e-08,2.58e-11
1e-09,3.13e-08,1.23e-11
5e-10,3.13e-08,7.73e-11
2e-10,3.45e-07,9.07e-11
1e-10,1.22e-06,1.31e-10
5e-11,4.36e-06,3.29e-10
2e-11,7.51e-06,1.19e-09
1e-11,8.2e-06,1.24e-09
5e-12,2.32e-05,3.95e-09
2e-12,4.15e-05,1.09e-08
1e-12,6e-05,1.47e-09
5e-13,0.000374,3.64e-08
2e-13,0.000393,7.99e-08
1e-13,0.00252,2.04e-07
5e-14,0.00284,8.81e-07
2e-14,0.0101,1.42e-06
1e-14,0.00539,1.42e-06

function,default_step,default_error
sqrt,1.21e-05,6.04e-12
peak,6.06e-06,0.000343
)";
    if (!(docs::printed() == docs::quoted(scale_table))) return 1;

    // What the page says. The two valleys: the smallest errors lie at steps a factor 2^10 to 2^16 apart, the peak's to
    // the left (smaller), as a function that varies on a scale 10^4 times smaller predicts: the balanced step
    // eps^(1/3) times the scale.
    double const apart = hs[best(slow)] / hs[best(fast)];
    if (!(apart > 1024.0 && apart < 65536.0)) return 3;
    // The default step: 1.2e-5 for sqrt at 2, at the bottom of its valley (its error within ten times the smallest);
    // 6.1e-6 for the peak, whose scale max(|x|, 1) is 1, on its truncation slope: an error between 1e-4 and 1e-3,
    // "three or four digits" (the truncation, 0.0933 (h / a)^2 = 3.4e-4, as deriv_step.cpp derives).
    if (!(slow.default_error <= 10.0 * slow.errors[best(slow)])) return 4;
    if (!(fast.default_step > 6e-6 && fast.default_step < 6.1e-6)) return 5;
    if (!(fast.default_error > 1e-4 && fast.default_error < 1e-3)) return 6;
    // The step on the peak's own scale, eps^(1/3) a = 6.1e-10, gives an error within ten times the smallest of the
    // peak's curve: its truncation, 3.4e-12, and its rounding, at most about 1.5 eps |f| / (h |f'|) = 1.9e-10 relative
    // (|f| = 4000, |f'| = 1.2e7), put it in the valley's bottom.
    double const scaled_step = nxx::math::root_epsilon<3, double>() * a;
    auto const scaled = dv::derivative(peak, peak_x, scaled_step);
    if (!scaled) return 7;
    double const scaled_error = std::abs(scaled->value - peak_exact) / std::abs(peak_exact);
    if (!(scaled_error <= 10.0 * fast.errors[best(fast)])) return 8;
    return 0;
}

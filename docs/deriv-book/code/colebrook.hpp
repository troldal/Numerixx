// The friction-factor function of Chapter 7 (E20, DS21) and of the noise-estimate recipe of Chapter 3 (E21, DS22):
// the Darcy friction factor f of a rough pipe from the Colebrook equation, 1/sqrt(f) = -2 log10(eps_r / 3.7 +
// 2.51 / (Re sqrt f)), in its common form (C. F. Colebrook, J. Inst. Civ. Eng. 11(4), 133-156, 1939). Written in
// x = 1/sqrt(f) the equation is G(x) = x + 2 log10(a + b x / Re) = 0 with a = eps_r / 3.7 and b = 2.51; G increases
// with x, so bisection on [1, 20] (f from 1 to 1/400) finds the root. The generator and the examples include this
// one file, so that the solver exists in exactly one place. Parameters (data/README.md): eps_r = 1e-4, Re = 1e5,
// eta = 1e-10.

#pragma once

#include <cmath>
#include <numbers>

namespace book {

// clang-format off
// Snippet lines stay within 60 characters (code/README.md).
// <implicit-function>
// The Darcy friction factor f(Re): the Colebrook equation
// in x = 1/sqrt(f), relative roughness 1e-4, is G(x) = 0.
// G grows with x: bisect on [1, 20] to the tolerance eta.
inline double friction(double Re, double eta)
{
    auto const G = [Re](double x) {
        return x
            + 2 * std::log10(1e-4 / 3.7 + 2.51 * x / Re);
    };
    double lo = 1.0, hi = 20.0;
    while (hi - lo > eta * hi) {
        double const mid = (lo + hi) / 2;
        if (G(mid) > 0) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    double const x = (lo + hi) / 2;
    return 1 / (x * x);
}
// </implicit-function>
// clang-format on

inline constexpr double colebrook_a = 1e-4 / 3.7;
inline constexpr double colebrook_b = 2.51;

// Fixed-point iteration x <- -2 log10(a + b x / Re) from x = 7 until the step is below eta |x|: the output is
// smooth in Re except where the iteration count changes (the comparison of decision S2-4, not the book's example).
inline double friction_fixed_point(double Re, double eta)
{
    double x = 7.0;
    for (int i = 0; i < 1000; ++i) {
        double const next = -2 * std::log10(colebrook_a + colebrook_b * x / Re);
        bool const done = std::abs(next - x) <= eta * std::abs(next);
        x = next;
        if (done) break;
    }
    return 1 / (x * x);
}

// The answer key: the root to machine precision by Newton's method on G, from x = 7.
inline double friction_exact(double Re)
{
    double x = 7.0;
    for (int i = 0; i < 100; ++i) {
        double const u = colebrook_a + colebrook_b * x / Re;
        double const g = x + 2 * std::log10(u);
        double const dg = 1 + 2 / std::numbers::ln10 * (colebrook_b / Re) / u;
        double const next = x - g / dg;
        if (next == x) break;
        x = next;
    }
    return 1 / (x * x);
}

// df/dRe by implicit differentiation, df/dRe = -G_Re / G_f, with G written in f: G = f^(-1/2) + 2 log10(u),
// u = a + b / (Re sqrt f).
inline double friction_slope_exact(double Re)
{
    double const f = friction_exact(Re);
    double const u = colebrook_a + colebrook_b / (Re * std::sqrt(f));
    double const c = 2 / std::numbers::ln10;
    double const G_f = -0.5 * std::pow(f, -1.5) * (1 + c * colebrook_b / (Re * u));
    double const G_Re = -c * colebrook_b / (u * Re * Re * std::sqrt(f));
    return -G_Re / G_f;
}

} // namespace book

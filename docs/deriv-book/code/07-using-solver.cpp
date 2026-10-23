// E20 (Chapter 7, section 7.4; values DS21): a function with no formula at hand. The Darcy friction factor f of a
// rough pipe follows from the Colebrook equation, 1/sqrt(f) = -2 log10(eps_r / 3.7 + 2.51 / (Re sqrt(f))), which has
// no closed-form solution, so it is solved by an inner iteration; here bisection to a relative tolerance eta (PLAN.md,
// decision S2-4), whose output is quantised at eta: the noise of section 3.6. Differentiated with respect to the
// Reynolds number with the default step, with the step eta^(1/3) Re, and by Ridders. Parameters as in
// data/README.md (the generator uses the same): eps_r = 1e-4, Re = 1e5, eta = 1e-10, bisection on [1, 20] in
// x = 1/sqrt(f). The answer key is the implicit-function derivative df/dRe = -G_Re / G_f of the Colebrook equation
// at its root found to machine precision by Newton's method.

#include <numerixx/deriv.hpp>

#include <cmath>
#include <numbers>
#include <print>

using nxx::deriv::derivative;
using nxx::deriv::Ridders;

// clang-format off
// <implicit-function>
// The Darcy friction factor f(Re) of a pipe with relative
// roughness 1e-4, from the Colebrook equation written in
// x = 1/sqrt(f): G(x) = x + 2 log10(1e-4/3.7 + 2.51 x/Re)
// = 0. G grows with x, so bisection on [1, 20] finds the
// root, to the relative tolerance eta.
double friction(double Re, double eta)
{
    auto const G = [Re](double x) {
        return x
            + 2 * std::log10(1e-4 / 3.7 + 2.51 * x / Re);
    };
    double lo = 1.0;
    double hi = 20.0;
    while (hi - lo > eta * hi) {
        double const mid = (lo + hi) / 2;
        if (G(mid) > 0)
            hi = mid;
        else
            lo = mid;
    }
    double const x = (lo + hi) / 2;
    return 1 / (x * x);
}
// </implicit-function>
// clang-format on

namespace {

    constexpr double colebrook_a = 1e-4 / 3.7;
    constexpr double colebrook_b = 2.51;

    // The root to machine precision by Newton's method on G, from x = 7.
    double friction_exact(double Re)
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

    // df/dRe = -G_Re / G_f, with G written in f: G = f^(-1/2) + 2 log10(u), u = a + b / (Re sqrt(f)), where a and b
    // are colebrook_a and colebrook_b.
    double friction_slope_exact(double Re)
    {
        double const f = friction_exact(Re);
        double const u = colebrook_a + colebrook_b / (Re * std::sqrt(f));
        double const c = 2 / std::numbers::ln10;
        double const G_f = -0.5 * std::pow(f, -1.5) * (1 + c * colebrook_b / (Re * u));
        double const G_Re = -c * colebrook_b / (u * Re * Re * std::sqrt(f));
        return -G_Re / G_f;
    }

    double relative(double value, double exact) { return std::abs(value - exact) / std::abs(exact); }

} // namespace

int main()
{
    bool ok = true;
    // clang-format off
// <implicit-derivative>
double const Re = 1e5, eta = 1e-10;
auto const f = [eta](double re) {
    return friction(re, eta);
};
// The default step assumes f accurate to machine epsilon.
auto const d = derivative(f, Re);
// f is accurate to about eta: the step eta^(1/3) Re.
double const h = std::cbrt(eta) * Re;
auto const dh = derivative(f, Re, h);
// Ridders, whose estimate assumes machine epsilon too.
auto const r = derivative<Ridders>(f, Re);
// </implicit-derivative>
    // clang-format on
    if (!d || !dh || !r || !r->abs_error) {
        std::println("FAILED: a method did not return a value");
        return 1;
    }
    double const exact = friction_slope_exact(Re);
    std::println("f = {:.6f}, df/dRe = {:.6e} (implicit differentiation)", friction_exact(Re), exact);
    std::println("default step:  {:.6e}, relative error {:.1e}", d->value, relative(*d, exact));
    std::println("step {:.0f}:     {:.6e}, relative error {:.1e}", h, dh->value, relative(*dh, exact));
    std::println("Ridders:       {:.6e}, relative error {:.1e}, relative estimate {:.1e}, ratio {:.2f}", r->value,
        relative(*r, exact), *r->abs_error / std::abs(exact), std::abs(r->value - exact) / *r->abs_error);
    // The fitted step: the model gives about eta^(2/3) = 2e-7, the generator measures 7e-8 (DS21); the text claims
    // at most 1e-4. The default step's error is larger (measured 5.5e-5): the rounding branch, eta |f| / h, at a
    // step 75 times too small. Ridders' estimate is reported, not checked (it covers here, ratio 0.066, but it
    // assumes f accurate to eps, section 3.6).
    ok = ok && relative(*dh, exact) <= 1e-4 && relative(*d, exact) > relative(*dh, exact);
    if (!ok) std::println("FAILED: a claim of the text does not hold");
    return ok ? 0 : 1;
}

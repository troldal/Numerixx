// E20 (Chapter 7, section 7.4; values DS21): a function with no formula at hand. The Darcy friction factor f of a
// rough pipe follows from the Colebrook equation, 1/sqrt(f) = -2 log10(eps_r / 3.7 + 2.51 / (Re sqrt(f))), which has
// no closed-form solution, so it is solved by an inner iteration; here bisection to a relative tolerance eta (PLAN.md,
// decision S2-4), whose output is quantised at eta: the noise of section 3.6. Differentiated with respect to the
// Reynolds number with the default step, with the step eta^(1/3) Re, and by Ridders. Parameters as in
// data/README.md (the generator uses the same): eps_r = 1e-4, Re = 1e5, eta = 1e-10, bisection on [1, 20] in
// x = 1/sqrt(f). The answer key is the implicit-function derivative df/dRe = -G_Re / G_f of the Colebrook equation
// at its root found to machine precision by Newton's method. The solver and the answer key are those of the shared
// header colebrook.hpp (book::friction, book::friction_exact, book::friction_slope_exact), which the generator
// and E21 use too.

#include <numerixx/deriv.hpp>

#include "colebrook.hpp"

#include <cmath>
#include <print>

using nxx::deriv::derivative;
using nxx::deriv::Ridders;

namespace {

    double relative(double value, double exact) { return std::abs(value - exact) / std::abs(exact); }

} // namespace

int main()
{
    bool ok = true;
    // clang-format off
// <implicit-derivative>
double const Re = 1e5, eta = 1e-10;
auto const f = [eta](double re) {
    return book::friction(re, eta);
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
    double const exact = book::friction_slope_exact(Re);
    std::println("f = {:.6f}, df/dRe = {:.6e} (implicit differentiation)", book::friction_exact(Re), exact);
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

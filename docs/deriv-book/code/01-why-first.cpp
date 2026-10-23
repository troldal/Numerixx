// E1 (Chapter 1, section 1.4): the first example of the book. The heat capacity of nitrogen at 1000 K as the
// derivative of its enthalpy, first with the default stencil, then with Ridders' method and its error estimate. The
// enthalpy comes from code/shomate.hpp (quoted in the chapter with \snippetheader); its closed-form cp is the answer
// key, so main can check what the text claims. The same numbers are DS1b of the generator (\valCpDefaultRelError,
// \valCpRidders...).

#include "shomate.hpp"

#include <numerixx/deriv.hpp>

#include <cmath>
#include <print>

int main()
{
    // clang-format off
// <first>
// H(T) in J/mol for T in K, so dH/dT is cp in J/(mol K)
auto const H = [](double T) {
    return book::shomate_H(book::nitrogen, T);
};
double const T0 = 1000.0; // K, inside the fit's range

auto const d = nxx::deriv::derivative(H, T0);
if (d) std::println("cp = {:.10f} J/(mol K)", d->value);
// </first>

// <first-ridders>
using nxx::deriv::Ridders;
auto const r = nxx::deriv::derivative<Ridders>(H, T0);
if (r) std::println("cp = {:.12f} +- {:.1e} J/(mol K)",
                    r->value, *r->abs_error);
// </first-ridders>
    // clang-format on

    // The answer key: the closed-form cp of the same fit.
    double const exact = book::shomate_cp(book::nitrogen, T0);
    if (T0 != book::thread_T0 || !d || !r) return 1;
    double const d_rel = std::abs(d->value - exact) / exact;
    double const r_actual = std::abs(r->value - exact);
    std::println("exact cp = {:.12f} J/(mol K)", exact);
    std::println("default step: relative error {:.1e}", d_rel);
    std::println("Ridders: actual error {:.1e}, estimate {:.1e}, ratio {:.2f}", r_actual, *r->abs_error,
        r_actual / *r->abs_error);

    // The stencil: the error model gives about 4e-11 relative at the default step (measured 2.4e-11); 1e-9 leaves
    // room for a different libm. Ridders: covered within the tests' factor 6 (measured ratio 0.73, 14 evaluations).
    bool const ok = d_rel <= 1e-9 && r_actual <= 6.0 * *r->abs_error;
    return ok ? 0 : 1;
}

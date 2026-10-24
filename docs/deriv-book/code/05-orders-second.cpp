// E7 (Chapter 5): the second derivative of the nitrogen example, d cp / dT = H''(T0), with the three-point central rule and
// with Ridders. The stencil's relative error is at most 1e-6 (the model gives about 2e-7: 4 eps |H| / (h^2 |H''|)
// with h = 0.12 K, the factor |H| / (s^2 |H''|) being about 4 at T0; measured 6.5e-8); Ridders' estimate covers its
// error. Quoted with \snippet{05-orders-second}{second-derivative} and \snippet{05-orders-second}{order-enum}.

#include "shomate.hpp"

#include <numerixx/deriv.hpp>

#include <cmath>
#include <print>

int main()
{
    using namespace nxx::deriv;
    auto const H = [](double T) { return book::shomate_H(book::nitrogen, T); };
    double const T0 = book::nitrogen_T0;
    // clang-format off
// <second-derivative>
// d cp / dT = H''(T0): the three-point rule with its
// default step (eps^(1/4) T0, about 0.12 K), and Ridders
// with an error estimate.
auto const d = derivative<Central<2>, 2>(H, T0);
auto const r = derivative<Ridders, order::second>(H, T0);
// </second-derivative>
    // clang-format on
    // clang-format off
// <order-enum>
// order::second and the integer 2 are the same template
// argument; the enumeration only reads better in a call.
auto const a = derivative<Central<2>, order::second>(H, T0);
auto const b = derivative<Central<2>, 2>(H, T0);
// </order-enum>
    // clang-format on
    auto const &by_enum = a;
    auto const &by_int = b;
    double const exact = book::shomate_derivative(book::nitrogen, T0, 2);
    std::println("H'' exact     {:.10e} J/(mol K^2)", exact);
    if (!d || !r) return 1;
    double const stencil_error = std::abs(d->value - exact) / std::abs(exact);
    double const ridders_actual = std::abs(r->value - exact);
    std::println("Central<2>    {:.10e}  rel. error {:.1e}", d->value, stencil_error);
    std::println("Ridders       {:.10e}  error {:.1e}  estimate {:.1e}", r->value, ridders_actual, *r->abs_error);
    bool const same = by_enum == by_int; // Estimate's == compares value and abs_error
    bool const ok = stencil_error <= 1e-6 && ridders_actual <= *r->abs_error && same;
    return ok ? 0 : 1;
}

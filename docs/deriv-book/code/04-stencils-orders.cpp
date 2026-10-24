// E5 (Chapter 4): higher-order central stencils on the nitrogen enthalpy. Central<2>, Central<4> and Central<6>
// all give cp = dH/dT at T0 with their default steps; the relative error of each against the closed form is at most
// 1e-9 (measured 2.4e-11, 2.1e-13 and 3.3e-14; the model about 4e-11 for Central<2>). Quoted with
// \snippet{04-stencils-orders}{higher-order}.

#include "shomate.hpp"

#include <numerixx/deriv.hpp>

#include <cmath>
#include <print>

namespace {

    // The relative error of an estimate against the closed form, or 1 when there is no estimate.
    double rel_error(std::expected<nxx::deriv::Estimate<double>, nxx::deriv::error> const &d, double exact)
    { return d ? std::abs(d->value - exact) / std::abs(exact) : 1.0; }

} // namespace

int main()
{
    using namespace nxx::deriv;
    auto const H = [](double T) { return book::shomate_H(book::nitrogen, T); };
    double const T0 = book::nitrogen_T0;
    // clang-format off
// <higher-order>
// cp = dH/dT at T0 with the default step of each rule: more
// points, a higher order, and a larger default step.
auto const d2 = derivative<Central<2>>(H, T0);
auto const d4 = derivative<Central<4>>(H, T0);
auto const d6 = derivative<Central<6>>(H, T0);
// </higher-order>
    // clang-format on
    double const cp = book::shomate_cp(book::nitrogen, T0);
    double const e2 = rel_error(d2, cp);
    double const e4 = rel_error(d4, cp);
    double const e6 = rel_error(d6, cp);
    std::println("cp exact      {:.12f} J/(mol K)", cp);
    if (d2) std::println("Central<2>    {:.12f}  rel. error {:.1e}", d2->value, e2);
    if (d4) std::println("Central<4>    {:.12f}  rel. error {:.1e}", d4->value, e4);
    if (d6) std::println("Central<6>    {:.12f}  rel. error {:.1e}", d6->value, e6);
    // Each at most 1e-9: the model gives about 4e-11 for Central<2> and less for the higher orders; the margin
    // covers rounding differences between platforms (the nitrogen function needs no libm).
    bool const ok = e2 <= 1e-9 && e4 <= 1e-9 && e6 <= 1e-9;
    return ok ? 0 : 1;
}

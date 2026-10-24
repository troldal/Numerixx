// E3 (Chapter 3, sections 3.2 to 3.4): a smaller step is not always better. On the nitrogen example, a step of 1e-2 K is
// near the bottom of the V (F3), 1e-8 K is far down its rounding branch, and the default step beats the latter. A
// step below the spacing of the doubles at 1000 K cannot move T0 at all, and the library reports non_finite.

#include "shomate.hpp"

#include <numerixx/deriv.hpp>

#include <cmath>
#include <print>

int main()
{
    auto const H = [](double T) { return book::shomate_H(book::nitrogen, T); };
    double const T0 = book::nitrogen_T0; // 1000 K
    double const cp = book::shomate_cp(book::nitrogen, T0);
    using nxx::deriv::derivative;
    using nxx::deriv::error;

    // clang-format off
// <explicit-step>
// the step, in K, as the third argument
for (double const h : {1e-2, 1e-8}) {
    auto const d = derivative(H, T0, h);
    if (d) {
        std::println("h = {:.0e} K: rel. error {:.1e}",
                     h, std::abs(d->value - cp) / cp);
    }
}
// </explicit-step>

// <too-small>
// T0 + 1e-14 K rounds back to T0: there is no difference
auto const tiny = derivative(H, T0, 1e-14);
if (!tiny && tiny.error() == error::non_finite) {
    std::println("h = 1e-14 K: non_finite");
}
// </too-small>
    // clang-format on

    // What the text claims: the step 1e-14 K is below the spacing of the doubles at T0, so non_finite; the default
    // step's error is below that of the 1e-8 K step (DS1: 2.4e-11 against 9.5e-7).
    double const spacing = std::nextafter(T0, 2 * T0) - T0;
    auto const d_default = derivative(H, T0);
    auto const d_small = derivative(H, T0, 1e-8);
    if (!d_default || !d_small) return 1;
    double const e_default = std::abs(d_default->value - cp) / cp;
    double const e_small = std::abs(d_small->value - cp) / cp;
    std::println("spacing of the doubles at {} K: {:.1e} K", T0, spacing);
    std::println("default step: relative error {:.1e}", e_default);
    bool const ok = 1e-14 < spacing / 2 && !tiny && tiny.error() == nxx::deriv::error::non_finite
                    && e_default < e_small;
    return ok ? 0 : 1;
}

// <program>
// The heat capacity of nitrogen at 1000 K: cp = dH/dT.
#include "shomate.hpp" // the book's H(T), namespace book

#include <numerixx/deriv.hpp>

#include <print>

int main()
{
    // H(T) in J/mol for T in K, so dH/dT is in J/(mol K)
    auto const H = [](double T) {
        return book::shomate_H(book::nitrogen, T);
    };
    double const T0 = 1000.0; // K, inside the fit's range

    auto const d = nxx::deriv::derivative(H, T0);
    if (d) {
        std::println("cp = {:.10f} J/(mol K)", d->value);
    } else {
        std::println("no estimate: error::non_finite");
    }
    return d ? 0 : 1;
}
// </program>

// E2 (Chapter 2, section 2.3): the accuracy order measured on the nitrogen example. Halving the step halves the error of the
// forward difference (order 1) and quarters that of the central difference (order 2), as long as truncation
// dominates: at steps of 1 and 2 K the rounding error is about five orders of magnitude below the truncation error.

#include "shomate.hpp"

#include <numerixx/deriv.hpp>

#include <cmath>
#include <print>

int main()
{
    auto const H = [](double T) { return book::shomate_H(book::nitrogen, T); };
    double const T0 = book::nitrogen_T0; // 1000 K

    // clang-format off
// <forward-vs-central>
using nxx::deriv::derivative;
using nxx::deriv::Central, nxx::deriv::Forward;
double const cp = book::shomate_cp(book::nitrogen, T0);

// steps of 2 K and 1 K: the second is half the first
auto const f2 = derivative<Forward<1>>(H, T0, 2.0);
auto const f1 = derivative<Forward<1>>(H, T0, 1.0);
auto const c2 = derivative<Central<2>>(H, T0, 2.0);
auto const c1 = derivative<Central<2>>(H, T0, 1.0);
if (!f2 || !f1 || !c2 || !c1) return 1;

// absolute errors in J/(mol K), and their ratio
double const ef2 = std::abs(f2->value - cp);
double const ef1 = std::abs(f1->value - cp);
double const ec2 = std::abs(c2->value - cp);
double const ec1 = std::abs(c1->value - cp);
std::println("forward: {:.2e} then {:.2e}, ratio {:.3f}",
             ef2, ef1, ef2 / ef1);
std::println("central: {:.2e} then {:.2e}, ratio {:.3f}",
             ec2, ec1, ec2 / ec1);
// </forward-vs-central>
    // clang-format on

    // Order 1 halves, order 2 quarters; the next terms of the expansions shift the ratios by less than 1e-3 here.
    double const forward_ratio = ef2 / ef1;
    double const central_ratio = ec2 / ec1;
    bool const ok = forward_ratio >= 1.9 && forward_ratio <= 2.1 && central_ratio >= 3.9 && central_ratio <= 4.1;
    return ok ? 0 : 1;
}

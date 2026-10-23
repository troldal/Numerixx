// E19 (Chapter 7, section 7.8; figure F18, data set DS20): "evaluable is not smooth". The published enthalpy of
// nitrogen joins two Shomate ranges at 500 K, and the two fits do not meet: H jumps by about 0.77 J/mol there (a range
// break; book::nitrogen_H_joined in code/shomate.hpp, PLAN.md decision S2-2). At 1 mK above the break, inside the span
// of every default step, both stencils return a finite but wrong value without failing, and Ridders returns a value
// with a large error estimate; at 50 K above it, outside every span, all three agree with the closed form. The answer
// key is the closed form of the range that x lies in. The points are those of the generator (data/README.md).

#include <numerixx/deriv.hpp>

#include "shomate.hpp"

#include <cmath>
#include <print>

using nxx::deriv::Central;
using nxx::deriv::derivative;
using nxx::deriv::Ridders;

namespace {

    double relative(double value, double exact) { return std::abs(value - exact) / std::abs(exact); }

} // namespace

int main()
{
    bool ok = true;
    // clang-format off
// <kink>
// H of nitrogen as published: the fit for 100 to 500 K
// below 500 K, the fit for 500 to 2000 K from 500 K on.
auto const H = [](double T) {
    return book::nitrogen_H_joined(T);
};
double const T = book::junction_T + 1e-3; // 1 mK above
auto const c2 = derivative(H, T);              // succeeds
auto const c4 = derivative<Central<4>>(H, T);  // succeeds
auto const r = derivative<Ridders>(H, T);      // succeeds
if (c2 && c4 && r) {
    std::println("{:.4f} {:.4f} {:.4f} +- {:.2f}",
        c2->value, c4->value, r->value, *r->abs_error);
}
// </kink>
    // clang-format on
    if (!c2 || !c4 || !r || !r->abs_error) {
        std::println("FAILED: a method did not return a value");
        return 1;
    }
    double const exact = book::nitrogen_derivative_joined(T, 1);
    std::println("1 mK above: exact {:.4f}; relative errors {:.1e} (Central<2>), {:.1e} (Central<4>), "
                 "{:.1e} (Ridders), Ridders' relative estimate {:.1e}",
        exact, relative(*c2, exact), relative(*c4, exact), relative(*r, exact), *r->abs_error / std::abs(exact));
    // Inside the spans: the generator measures 4.3 (Central<2>, the jump over 2h), 0.041 (Central<4>) and a Ridders
    // estimate of 6.2e-3 of the value, eleven orders above its usual size (DS20). The text claims more than 1e-3 for
    // each, and Ridders' estimate covering its own error.
    ok = ok && relative(*c2, exact) > 1e-3 && relative(*c4, exact) > 1e-3;
    ok = ok && *r->abs_error > 1e-3 * std::abs(r->value) && std::abs(r->value - exact) <= *r->abs_error;

    // clang-format off
// <kink-outside>
double const T_far = book::junction_T + 50.0; // 50 K above
auto const f2 = derivative(H, T_far);
auto const f4 = derivative<Central<4>>(H, T_far);
auto const rf = derivative<Ridders>(H, T_far);
// </kink-outside>
    // clang-format on
    if (!f2 || !f4 || !rf || !rf->abs_error) {
        std::println("FAILED: a method did not return a value");
        return 1;
    }
    double const exact_far = book::nitrogen_derivative_joined(T_far, 1);
    std::println("50 K above: {:.10f} {:.10f} {:.10f} +- {:.1e}; exact {:.10f}", f2->value, f4->value, rf->value,
        *rf->abs_error, exact_far);
    // Outside every span, the usual accuracy: the generator measures below 1e-11 (DS20); the text claims 1e-8.
    ok = ok && relative(*f2, exact_far) <= 1e-8 && relative(*f4, exact_far) <= 1e-8 && relative(*rf, exact_far) <= 1e-8;

    if (!ok) std::println("FAILED: a claim of the text does not hold");
    return ok ? 0 : 1;
}

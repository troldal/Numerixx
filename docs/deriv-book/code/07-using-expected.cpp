// E11 (Chapter 7, section 7.2): the result of derivative(), a std::expected<Estimate<T>, error>. How to test it,
// read it, chain it, and the traps of the conversion of Estimate<T> to T. The thread's enthalpy H(T) comes from the
// shared header; T0 = 1000 K.

// clang-format off
// <names>
#include <numerixx/deriv.hpp>

using nxx::deriv::derivative; // the one function
using nxx::deriv::Central;    // the methods
using nxx::deriv::Ridders;
using nxx::deriv::Estimate;   // a value, an error estimate
using nxx::deriv::error;      // why there is no value
// </names>
// clang-format on

#include "shomate.hpp"

#include <cmath>
#include <concepts>
#include <print>

namespace {

    double H(double T) { return book::shomate_H(book::nitrogen, T); }

    double relative(double value, double exact) { return std::abs(value - exact) / std::abs(exact); }

} // namespace

int main()
{
    bool ok = true;
    double const T0 = book::thread_T0;

    // clang-format off
// <expected-basics>
auto const ln = [](double x) { return std::log(x); };
auto const d = derivative(ln, 0.0); // samples log(-h): NaN
if (d) {
    std::println("slope {}", d->value);
} else if (d.error() == error::non_finite) {
    std::println("no finite estimate at x = 0");
}

auto const g = derivative(ln, 2.0); // succeeds
if (g) {
    double const a = g->value;  // the member
    double const b = *g;        // the conversion to double
    // A stencil has no error estimate: abs_error is empty.
    double const e = g->abs_error.value_or(0.0);
    std::println("{} {} {}", a, b, e);
}
// </expected-basics>
    // clang-format on
    ok = ok && !d && d.error() == error::non_finite;
    ok = ok && g && g->value == double(*g) && !g->abs_error.has_value();
    ok = ok && relative(g->value, 0.5) <= 1e-9; // the model gives about eps^(2/3) = 4e-11

    // clang-format off
// <expected-chain>
// cp in J/(kg K) from cp in J/(mol K): M = 28.0134 g/mol
auto const cp = derivative(H, T0).transform(
    [](double c) { return c / 28.0134e-3; });

// d ln cp / dT = H'' / H' in 1/K: two calls that can each
// fail, combined; the first failure is the result.
auto const slope = derivative(H, T0).and_then(
    [&](double c) {
        return derivative<Central<2>, 2>(H, T0).transform(
            [c](double c2) { return c2 / c; });
    });
// </expected-chain>
    // clang-format on
    double const cp_exact = book::shomate_cp(book::nitrogen, T0);
    double const slope_exact = book::shomate_derivative(book::nitrogen, T0, 2) / cp_exact;
    std::println("cp = {:.6f} J/(kg K), d ln cp/dT = {:.6e} 1/K (exact {:.6e})", cp.value_or(0.0),
        slope.value_or(0.0), slope_exact);
    // cp: Central<2> at the default step, relative error about 2.4e-11 (DS1b); 1e-9 covers libm differences.
    ok = ok && cp && relative(*cp, cp_exact / 28.0134e-3) <= 1e-9;
    // H'' by Central<2> at the default step: about 2e-7 relative by the model (PLAN.md, E7); 1e-5 leaves room.
    ok = ok && slope && relative(*slope, slope_exact) <= 1e-5;

    // clang-format off
// <estimate-pitfalls>
auto const cosine = [](double x) { return std::cos(x); };
auto const r = derivative<Ridders>(cosine, 0.0);
auto const e = *r;     // an Estimate<double>, not a double
double const v = *r;   // the value, by the conversion
if (r) std::println("a result: {}", v); // true: r holds one
// if (*r) would test v != 0, not whether r holds a value.
// Test has_value() before reading the error estimate.
if (r->abs_error.has_value()) {
    std::println("{} +- {:.1e}", e.value, *e.abs_error);
}
// </estimate-pitfalls>
    // clang-format on
    static_assert(std::same_as<decltype(e), Estimate<double> const>);
    // cos is even about 0, so cos(h) - cos(-h) is exactly 0 at every level: the value is exactly 0, and *r tests
    // false although r holds a value (written out as a cast: Clang's -Wconversion rejects if (*r) itself).
    ok = ok && r && v == 0.0 && !static_cast<bool>(*r) && e.abs_error.has_value();

    if (!ok) std::println("FAILED: a claim of the text does not hold");
    return ok ? 0 : 1;
}

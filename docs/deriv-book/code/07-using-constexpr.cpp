// E15 (Chapter 7, section 7.7): derivative() in constant expressions. Results and failures are values at compile time
// as at run time; the static_asserts are the checks of this example. With steps that are powers of two, x^3 gives
// exact results in binary floating point: the central difference is f' + h^2 f'''/6 = 12 + 0.0625 at x = 2,
// h = 1/4, and Ridders' first column removes the h^2 term, so it returns 12 exactly.

#include <numerixx/deriv.hpp>

#include "shomate.hpp"

#include <limits>
#include <print>

using nxx::deriv::derivative;
using nxx::deriv::error;
using nxx::deriv::Forward;
using nxx::deriv::Ridders;

// clang-format off
// <constexpr>
constexpr auto cube = [](double x) { return x * x * x; };
// The truncation error h^2 f'''/6 = 0.0625, exactly:
static_assert(derivative(cube, 2.0, 0.25) == 12.0625);
// Ridders' extrapolation removes it:
static_assert(derivative<Ridders>(cube, 2.0, 0.0625)
              == 12.0);
// </constexpr>

// <constexpr-failure>
// The thread's fit, NaN outside its range: a failure is a
// value at compile time too.
constexpr auto H_fit = [](double T) {
    if (T < 500.0)
        return std::numeric_limits<double>::quiet_NaN();
    return book::shomate_H(book::nitrogen, T);
};
constexpr auto at_edge = derivative(H_fit, 500.0);
static_assert(!at_edge);
static_assert(at_edge.error() == error::non_finite);
// The one-sided rule samples at and above 500 K only:
static_assert(derivative<Forward<2>>(H_fit, 500.0));
// </constexpr-failure>

// <constexpr-value>
// cp at 1000 K, computed by the compiler (default step),
// against the closed form:
constexpr auto cp = derivative(H_fit, 1000.0);
constexpr double exact =
    book::shomate_cp(book::nitrogen, 1000.0);
static_assert(cp);
static_assert(nxx::math::abs(cp->value - exact)
              < 1e-9 * exact);
// </constexpr-value>
// clang-format on

int main()
{
    // The same call at run time. The two agree bit for bit when the compiler neither contracts a * b + c into a
    // fused multiply-add nor uses extended precision (PLAN.md, style decision 15): true of the CI toolchains.
    double const T = 1000.0;
    auto const run_time = derivative(H_fit, T);
    std::println("compile time {:.15f}, run time {:.15f}", cp->value, run_time ? run_time->value : 0.0);
    bool const ok = run_time && *run_time == *cp;
    if (!ok) std::println("FAILED: compile-time and run-time results differ");
    return ok ? 0 : 1;
}

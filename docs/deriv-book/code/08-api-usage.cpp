// E16 (Chapter 8, sections 8.2 and 8.8): every accepted form of the call, on the thread's enthalpy at T0 = 1000 K,
// and the calls that do not compile, as requires-expressions that are false. main calls each accepted form and checks
// that it succeeds; the static_asserts check the rejected ones.

#include <numerixx/deriv.hpp>

#include "shomate.hpp"

#include <cmath>
#include <print>

using nxx::deriv::Backward;
using nxx::deriv::Central;
using nxx::deriv::derivative;
using nxx::deriv::Forward;
using nxx::deriv::order;
using nxx::deriv::Ridders;

// clang-format off
// <api-rejected>
constexpr auto sq = [](double x) { return x * x; };

// Whether derivative<M, N>(sq, x) compiles for x of type X.
template<typename M, auto N, typename X = double>
concept accepted =
    requires(X x) { derivative<M, N>(sq, x); };

static_assert(accepted<Central<4>, 1>);
static_assert(!accepted<Central<3>, 1>);  // P must be even
static_assert(!accepted<Central<12>, 1>); // P at most 10
static_assert(!accepted<Forward<7>, 1>);  // P at most 6
static_assert(!accepted<Central<2>, 5>);  // N from 1 to 4
static_assert(!accepted<Central<2>, static_cast<order>(5)>);
static_assert(!accepted<Central<2>, true>); // a bool
static_assert(!accepted<double, 1>);      // not a method
static_assert(!accepted<Central<2>, 1, int>); // x not real
// </api-rejected>
// clang-format on

namespace {

    double H(double T) { return book::shomate_H(book::nitrogen, T); }

    // |value - exact| <= tol |exact|
    bool close(double value, double exact, double tol) { return std::abs(value - exact) <= tol * std::abs(exact); }

} // namespace

int main()
{
    double const T0 = book::thread_T0;
    // clang-format off
// <api-overloads>
// f' by Central<2> with its default step
auto const a = derivative(H, T0);
// the same rule with the step h = 0.01 K
auto const b = derivative(H, T0, 0.01);
// another method: a template argument
auto const c = derivative<Central<4>>(H, T0);
auto const d = derivative<Forward<2>>(H, T0, 0.1);
// f'': the order as an integer or as an order
auto const e = derivative<Central<2>, 2>(H, T0);
auto const g =
    derivative<Backward<2>, order::second>(H, T0);
// Ridders: f' and f'', each with an error estimate
auto const r1 = derivative<Ridders>(H, T0);
auto const r2 = derivative<Ridders, order::second>(H, T0);
// </api-overloads>
    // clang-format on
    double const cp = book::shomate_cp(book::nitrogen, T0);
    double const cp_slope = book::shomate_derivative(book::nitrogen, T0, 2);
    auto const show = [](auto const &x) { return x ? x->value : std::nan(""); };
    std::println("cp:     {:.12f} {:.12f} {:.12f} {:.12f} {:.12f} (exact {:.12f})", show(a), show(b), show(c), show(d),
        show(r1), cp);
    std::println("dcp/dT: {:.9e} {:.9e} {:.9e} (exact {:.9e})", show(e), show(g), show(r2), cp_slope);
    bool ok = a && b && c && d && e && g && r1 && r2;
    // Relative tolerances from the error models, with room for libm differences (H''' = -5.8e-6 J/(mol K^3) at
    // T0): Central<2> at the default step about 2.4e-11 (DS1b); at h = 0.01 K the truncation h^2 |H'''| / (6 cp) is
    // 3e-12 and the rounding eps |H| / (h cp) 1.5e-11; Central<4> about 2e-13 (DS1b); Forward<2> at h = 0.1 K the
    // truncation h^2 |H'''| / (3 cp), 6e-10.
    ok = ok && close(*a, cp, 1e-9) && close(*b, cp, 1e-9) && close(*c, cp, 1e-9) && close(*d, cp, 1e-7);
    // f'' at the default steps: Central<2> about 2e-7 (E7); Backward<2> eps^(2/4) times a larger constant.
    ok = ok && close(*e, cp_slope, 1e-5) && close(*g, cp_slope, 1e-4);
    // Ridders: covered by 6 times the estimate, the tests' criterion.
    ok = ok && r1->abs_error && std::abs(r1->value - cp) <= 6 * *r1->abs_error;
    ok = ok && r2->abs_error && std::abs(r2->value - cp_slope) <= 6 * *r2->abs_error;
    ok = ok && !a->abs_error && !e->abs_error; // the stencils have no error estimate
    if (!ok) std::println("FAILED: an accepted form did not succeed as the text says");
    return ok ? 0 : 1;
}

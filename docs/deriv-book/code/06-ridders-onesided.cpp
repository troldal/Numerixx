// E24 (Chapter 6): Ridders' one-sided form at the edge of a domain. The nitrogen fit is valid on [500, 2000] K, and
// the wrapper of E6 (code/04-stencils-edge.cpp) returns NaN outside that range. At the lower end, 500 K, central
// Ridders samples below the range at every level and returns non_finite after its level cap (64 evaluations);
// ForwardRidders samples at 500 K and above only and returns cp and dcp/dT with error estimates, from its default
// initial step 500 K / 32. main checks: Ridders non_finite; each one-sided actual error within 6 times its estimate
// (the tests' criterion; measured 0.074 and 0.013 times); relative errors at most 1e-10 for f' and 1e-7 for f''
// (measured 2.1e-13 and 4.5e-9; f'' is the larger because |H| / (s^2 |H''|) is about 6 at 500 K). The printed
// numbers are those of \valEdgeForwardRidders... and \valEdgeRidders... (the same wrapper and default steps).
// Quoted with \snippet{06-ridders-onesided}{onesided-edge}.

#include "shomate.hpp"

#include <numerixx/deriv.hpp>

#include <cmath>
#include <limits>
#include <print>

int main()
{
    using namespace nxx::deriv;
    // E6's wrapper: NaN outside the fit's range, so that a method sampling there fails instead of returning a
    // plausible number.
    auto const H = [](double T) {
        using book::nitrogen;
        if (T < nitrogen.T_min || T > nitrogen.T_max) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return book::shomate_H(nitrogen, T);
    };
    // clang-format off
// <onesided-edge>
// H is NaN below 500 K: Ridders samples there and fails.
double const T = book::nitrogen.T_min; // 500 K
auto const c = derivative<Ridders>(H, T); // non_finite
auto const r1 = derivative<ForwardRidders>(H, T);
auto const r2 = derivative<ForwardRidders, 2>(H, T);
// </onesided-edge>
    // clang-format on
    if (c || c.error() != error::non_finite || !r1 || !r2) return 1;
    std::println("Ridders: {}", c.error());
    std::println("cp      = {:.10f} +- {:.1e}", r1->value, *r1->abs_error);
    std::println("dcp/dT  = {:.14f} +- {:.1e}", r2->value, *r2->abs_error);

    // The evaluations of each call, by the same calls through a counting wrapper (the function is deterministic).
    int calls = 0;
    auto const counted = [&calls, &H](double t) {
        ++calls;
        return H(t);
    };
    auto const count = [&calls](auto const &result) {
        (void)result;
        int const n = calls;
        calls = 0;
        return n;
    };
    int const n_c = count(derivative<Ridders>(counted, T));
    int const n_1 = count(derivative<ForwardRidders>(counted, T));
    int const n_2 = count(derivative<ForwardRidders, 2>(counted, T));

    double const cp = book::shomate_cp(book::nitrogen, T);
    double const cp2 = book::shomate_derivative(book::nitrogen, T, 2);
    double const a1 = std::abs(r1->value - cp);
    double const a2 = std::abs(r2->value - cp2);
    double const e1 = a1 / std::abs(cp);
    double const e2 = a2 / std::abs(cp2);
    std::println("initial step {} K", default_step<ForwardRidders>(T));
    std::println("Ridders        non_finite after {} evaluations", n_c);
    std::println("f'   exact {:.12f}, rel. error {:.1e}, rel. estimate {:.1e}, ratio {:.3f}, {} evaluations", cp, e1,
        *r1->abs_error / std::abs(cp), a1 / *r1->abs_error, n_1);
    std::println("f''  exact {:.12e}, rel. error {:.1e}, rel. estimate {:.1e}, ratio {:.3f}, {} evaluations", cp2, e2,
        *r2->abs_error / std::abs(cp2), a2 / *r2->abs_error, n_2);
    bool const ok = a1 <= 6 * *r1->abs_error && a2 <= 6 * *r2->abs_error && e1 <= 1e-10 && e2 <= 1e-7;
    return ok ? 0 : 1;
}

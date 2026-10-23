// E12 (Chapter 7, section 7.4; figure F16, data set DS11): Newton's method for the temperature at which the thread's
// enthalpy takes a given value, an energy balance, with the slope dH/dT = cp from derivative(). The same iteration as
// the generator's DS11: H_target = H(1000 K), start at 800 K, stop when the step is below 1e-9 K.

#include <numerixx/deriv.hpp>

#include "shomate.hpp"

#include <cmath>
#include <expected>
#include <print>

// clang-format off
// <newton>
// H(T) - H(298.15 K) of nitrogen in J/mol, T in K.
double H(double T)
{ return book::shomate_H(book::nitrogen, T); }

struct Newton
{
    double T;  // the last iterate, K
    int steps; // 20 means: not converged
};

// The temperature at which H(T) = H_target, by Newton's
// method from T: the slope dH/dT = cp from derivative().
std::expected<Newton, nxx::deriv::error>
solve_T(double H_target, double T)
{
    Newton n{ T, 0 };
    while (n.steps < 20) {
        auto const cp = nxx::deriv::derivative(H, n.T);
        if (!cp) return std::unexpected(cp.error());
        double const dT = (H(n.T) - H_target) / cp->value;
        n.T -= dT;
        ++n.steps;
        if (std::abs(dT) < 1e-9) break; // K
    }
    return n;
}
// </newton>
// clang-format on

int main()
{
    double const T_star = book::thread_T0;
    // clang-format off
// <newton-call>
auto const s = solve_T(H(1000.0), 800.0);
if (s && s->steps < 20) {
    std::println("T = {:.9f} K, {} steps", s->T, s->steps);
}
// </newton-call>
    // clang-format on
    if (!s) {
        std::println("FAILED: no slope");
        return 1;
    }
    double const error = std::abs(s->T - T_star);
    std::println("error {:.1e} K, residual {:.1e} J/mol", error, H(s->T) - H(T_star));
    // Quadratic convergence from 200 K away: the generator measures 4 steps and an error of about 1e-13 K (DS11),
    // the rounding of H (eps |H| / cp = 1.5e-13 K). The text claims at most 10 steps and 1e-6 K.
    bool const ok = s->steps <= 10 && error <= 1e-6;
    if (!ok) std::println("FAILED: Newton did not converge as the text says");
    return ok ? 0 : 1;
}

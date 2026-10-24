// E12 (Chapter 7, section 7.4; figure F16, data set DS11): Newton's method for the temperature at which the thread's
// enthalpy takes a given value, an energy balance, with the slope dH/dT = cp from derivative(). The same iteration as
// the generator's DS11: H_target = H(1000 K), start at 800 K, stop when the step is below 1e-9 K.

#include <numerixx/deriv.hpp>

#include "shomate.hpp"

#include <cmath>
#include <expected>
#include <print>

// H(T) - H(298.15 K) of nitrogen in J/mol, T in K.
double H(double T)
{ return book::shomate_H(book::nitrogen, T); }

// clang-format off
// <newton>
// Why solve_T returns no temperature.
enum class solve_error { no_slope, no_convergence };

struct Newton
{
    double T;  // the solution, K
    int steps; // the iterations it took
};

// The temperature at which H(T) = H_target, by Newton's
// method from T: the slope dH/dT = cp from derivative().
std::expected<Newton, solve_error>
solve_T(double H_target, double T)
{
    for (int steps = 1; steps <= 20; ++steps) {
        auto const cp = nxx::deriv::derivative(H, T);
        if (!cp) {
            return std::unexpected(solve_error::no_slope);
        }
        double const dT = (H(T) - H_target) / cp->value;
        T -= dT;
        if (std::abs(dT) < 1e-9) { // K
            return Newton{ T, steps };
        }
    }
    return std::unexpected(solve_error::no_convergence);
}
// </newton>
// clang-format on

int main()
{
    double const T_star = book::thread_T0;
    // clang-format off
// <newton-call>
auto const s = solve_T(H(1000.0), 800.0);
if (s) {
    std::println("T = {:.9f} K, {} steps", s->T, s->steps);
} else if (s.error() == solve_error::no_convergence) {
    std::println("no convergence in 20 steps");
}
// </newton-call>
    // clang-format on
    if (!s) {
        std::println("FAILED: no temperature");
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

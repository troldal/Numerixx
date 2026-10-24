// E21 (Chapter 3, section 3.6): measuring the noise level of a function whose accuracy is not known. The recipe
// book::estimate_noise of code/noise.hpp (a difference table at eight points, after More and Wild) is applied to the
// friction factor of Chapter 7 (book::friction of code/colebrook.hpp, bisection to the relative tolerance 1e-10, at
// Re = 1e5) and to the noisy exp of E17 (eta = 1e-8, at x = 1), with the points spaced by the default step of
// Central<2>, eps^(1/3) max(|x|, 1). The parameters are those of the generator (DS22, data/README.md), so the printed
// numbers are \valNoiseEstimateFriction and \valNoiseEstimateNoisy.

#include "colebrook.hpp"
#include "noise.hpp"

#include <numerixx/deriv.hpp>
#include <numerixx/math.hpp>

#include <cmath>
#include <print>

int main()
{
    // clang-format off
// <noise-estimate-use>
// f(Re) from a solver converged to the tolerance 1e-10
auto const f = [](double re) {
    return book::friction(re, 1e-10);
};
double const Re = 1e5;
double const h = nxx::deriv::default_step(Re);
auto const e = book::estimate_noise(f, Re, h); // 8 values
double eta = 1e-10; // the solver's tolerance, as fallback
if (e.k > 0) {      // the estimate settled
    eta = e.sigma / std::abs(f(Re));
}
double const h_eta = std::cbrt(eta) * Re; // the step
// </noise-estimate-use>
    // clang-format on
    std::println("noise of f: {:.1e} relative (k = {}), step {:.1f}", eta, e.k, h_eta);

    // The step from the measured noise gives the accuracy of section 3.6, about eta^(2/3) (the answer key is the
    // implicit derivative of code/colebrook.hpp).
    auto const d = nxx::deriv::derivative(f, Re, h_eta);
    if (!d) return 1;
    double const slope = book::friction_slope_exact(Re);
    double const rel_d = std::abs(d->value - slope) / std::abs(slope);
    std::println("derivative at that step: relative error {:.1e}", rel_d);

    // The noisy exp of E17: a uniform relative noise of half-width eta has the root mean square eta / sqrt(3).
    double const eta_g = 1e-8;
    auto const g = [eta_g](double x) { return book::noisy(std::exp(x), x, eta_g); };
    double const x = 1.0;
    double const hx = nxx::deriv::default_step(x);
    auto const n = book::estimate_noise(g, x, hx);
    if (e.k == 0 || n.k == 0) return 1;
    double const rel_f = e.sigma / std::abs(f(Re));
    double const rel_g = n.sigma / std::abs(g(x));
    std::println("noisy exp: {:.1e} relative (k = {}), against eta / sqrt(3) = {:.1e}", rel_g, n.k,
        eta_g / std::sqrt(3.0));

    // The friction factor's estimate within a factor 10 of the solver's tolerance (DS22 measures 4.3e-11 against a
    // root mean square deviation from the exact root of 4.9e-11); the noisy exp's within a factor 3 of
    // eta / sqrt(3); a difference order of at least 2, so that the smooth part of f has dropped out; and the
    // derivative at the step from the measured noise within 1e-4 (the model gives about eta^(2/3), 1.2e-7).
    bool const ok = rel_f >= 1e-11 && rel_f <= 1e-9 && rel_g >= eta_g / std::sqrt(3.0) / 3
                    && rel_g <= 3 * eta_g / std::sqrt(3.0) && e.k >= 2 && n.k >= 2 && rel_d <= 1e-4;
    return ok ? 0 : 1;
}

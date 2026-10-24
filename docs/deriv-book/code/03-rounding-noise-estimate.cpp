// E21 (Chapter 3, section 3.6): measuring the noise level of a function whose accuracy is not known. The recipe
// book::estimate_noise of code/noise.hpp (a difference table at eight points, after More and Wild) is applied to the
// friction factor of Chapter 7 (book::friction of code/colebrook.hpp, bisection to the relative tolerance 1e-10, at
// Re = 1e5) and to the noisy exp of E17 (eta = 1e-8, at x = 1), with the points spaced by the default step of
// Central<2>, eps^(1/3) max(|x|, 1). The parameters are those of the generator (DS22, data/README.md), so the printed
// numbers are \valNoiseEstimateFriction and \valNoiseEstimateNoisy.

#include "colebrook.hpp"
#include "noise.hpp"

#include <numerixx/math.hpp>

#include <algorithm>
#include <cmath>
#include <print>

int main()
{
    using nxx::math::root_epsilon;

    // clang-format off
// <noise-estimate-use>
// f(Re) from a solver converged to the tolerance 1e-10
auto const f = [](double re) {
    return book::friction(re, 1e-10);
};
double const Re = 1e5;
double const h = root_epsilon<3, double>()
                 * std::max(std::abs(Re), 1.0);
auto const e = book::estimate_noise(f, Re, h); // 8 values
// e.sigma / f(Re): the relative noise; e.k: the order
// </noise-estimate-use>
    // clang-format on
    if (e.k > 0) std::println("noise of f: {:.1e} relative (k = {})", e.sigma / std::abs(f(Re)), e.k);

    // The noisy exp of E17: a uniform relative noise of half-width eta has the root mean square eta / sqrt(3).
    double const eta = 1e-8;
    auto const g = [eta](double x) { return book::noisy(std::exp(x), x, eta); };
    double const x = 1.0;
    double const hx = root_epsilon<3, double>() * std::max(std::abs(x), 1.0);
    auto const n = book::estimate_noise(g, x, hx);
    if (e.k == 0 || n.k == 0) return 1;
    double const rel_f = e.sigma / std::abs(f(Re));
    double const rel_g = n.sigma / std::abs(g(x));
    std::println("noisy exp: {:.1e} relative (k = {}), against eta / sqrt(3) = {:.1e}", rel_g, n.k,
        eta / std::sqrt(3.0));

    // The friction factor's estimate within a factor 10 of the solver's tolerance (DS22 measures 4.3e-11 against a
    // root mean square deviation from the exact root of 4.9e-11); the noisy exp's within a factor 3 of
    // eta / sqrt(3); and a difference order of at least 2, so that the smooth part of f has dropped out.
    bool const ok = rel_f >= 1e-11 && rel_f <= 1e-9 && rel_g >= eta / std::sqrt(3.0) / 3
                    && rel_g <= 3 * eta / std::sqrt(3.0) && e.k >= 2 && n.k >= 2;
    return ok ? 0 : 1;
}

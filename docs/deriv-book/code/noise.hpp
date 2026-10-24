// The noisy function of Chapter 3 (E17, DS19), and the recipe that measures a noise level (E21, DS22). The generator
// and the examples include this one file, so that both perturb a value, and estimate a noise level, in exactly the
// same way.
//
// book::noisy: a function whose values carry a deterministic relative noise, what a solver output converged to a
// tolerance looks like to a finite difference.
//
// book::estimate_noise: the difference-table estimate of the noise level, a simplified form of the ECnoise procedure
// of J. J. More and S. M. Wild, "Estimating computational noise", SIAM J. Sci. Comput. 33(3), 1292-1314, 2011
// (the published procedure also tests the signs of the differences and adjusts the spacing). The k-th differences
// of the smooth part of f vanish like h^k while those of the noise do not, so sigma_k = sqrt(k!^2 / (2k)! times the
// mean square of the k-th differences) settles at the noise once k is large enough; k!^2 / (2k)! removes the growth
// of the differences of independent noise. The estimate is accepted at the first k whose sigma_k, sigma_(k+1) and
// sigma_(k+2) agree within a factor 4.

#pragma once

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

namespace book {

// clang-format off
// Snippet lines stay within 60 characters (code/README.md).
// <noise>
// f (1 + eta u) with u uniform in [-1, 1], drawn from a
// generator seeded with the bits of x: the same x always
// gives the same value, as a solver converged to a
// tolerance eta would.
inline double noisy(double f, double x, double eta)
{
    auto const bits = std::bit_cast<std::uint64_t>(x);
    auto const lo = static_cast<std::uint32_t>(bits);
    auto const hi = static_cast<std::uint32_t>(bits >> 32);
    std::mt19937 gen(lo ^ hi);
    double const u = double(gen()) / 4294967296.0;
    return f * (1.0 + eta * (2.0 * u - 1.0));
}
// </noise>

// <noise-estimate>
// The noise of fn near x from a difference table of its
// values at x, x + h, ..., x + m h (More and Wild).
struct NoiseEstimate
{
    double sigma; // the absolute noise; NaN if none settled
    int k;        // the difference order used; 0 if none
};
template<typename F>
NoiseEstimate estimate_noise(F fn, double x, double h,
                             std::size_t m = 7)
{
    std::vector<double> t(m + 1);
    for (std::size_t i = 0; i <= m; ++i)
        t[i] = fn(x + double(i) * h);
    std::vector<double> rms(m + 1, 0.0); // sigma_k
    double gamma = 1.0; // k!^2 / (2k)!
    for (std::size_t k = 1; k <= m; ++k) {
        for (std::size_t i = 0; i + k <= m; ++i)
            t[i] = t[i + 1] - t[i];
        gamma *= double(k) / double(2 * (2 * k - 1));
        double sum = 0.0;
        for (std::size_t i = 0; i + k <= m; ++i)
            sum += t[i] * t[i];
        rms[k] = std::sqrt(gamma * sum / double(m - k + 1));
    }
    for (std::size_t k = 1; k + 2 <= m; ++k) {
        auto const [lo, hi] = std::minmax(
            { rms[k], rms[k + 1], rms[k + 2] });
        if (hi <= 4 * lo) return { rms[k], int(k) };
    }
    return { std::numeric_limits<double>::quiet_NaN(), 0 };
}
// </noise-estimate>
// clang-format on

} // namespace book

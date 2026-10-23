// A function whose values carry a deterministic relative noise, the "noisy function" of Chapter 3 (E17, DS19): what
// a solver output converged to a tolerance, or a fitted table, looks like to a finite difference. The generator and
// the example include this one file, so that both perturb a value in exactly the same way.

#pragma once

#include <bit>
#include <cstdint>
#include <random>

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
// clang-format on

} // namespace book

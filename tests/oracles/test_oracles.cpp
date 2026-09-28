// Phase 0: the Boost.Math oracles are reachable. Later phases compare Numerixx's root finders, minimisers and
// quadrature rules against them (DESIGN §9.1). Errors go through a no-throw policy, as in the oracle tests to come.
#include <numerixx/core.hpp>

#include <boost/math/policies/policy.hpp>
#include <boost/math/tools/roots.hpp>
#include <doctest/doctest.h>

#include <cmath>
#include <cstdint>

TEST_SUITE("oracles")
{
    TEST_CASE("Boost.Math TOMS748 is available as an oracle")
    {
        namespace bmp  = boost::math::policies;
        using no_throw = bmp::policy<bmp::evaluation_error<bmp::errno_on_error>>;

        const auto     f          = [](double x) { return x * x - 2.0; };
        const auto     tolerance  = [](double a, double b) { return std::abs(b - a) < 1e-14; };
        std::uintmax_t iterations = 100;

        const auto [lo, hi] = boost::math::tools::toms748_solve(f, 1.0, 2.0, tolerance, iterations, no_throw {});
        CHECK(std::abs(0.5 * (lo + hi) - std::sqrt(2.0)) < 1e-13);
        CHECK(iterations < 100);
    }
}

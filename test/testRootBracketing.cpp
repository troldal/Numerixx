// ================================================================================================
// Catch2 test file for the root bracketing solvers.
// ================================================================================================

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <Roots.hpp>

#include <cmath>
#include <functional>
#include <vector>

TEST_CASE("nxx::roots - Bracketing Solver Tests", "[roots]")
{
    using namespace nxx::roots;

    // Test functions
    std::vector< std::function< double(double) > > functions {
        [](double x) { return std::sin(x) - x / 2.0; },
        [](double x) { return std::exp(x) - 3 * x; },
        [](double x) { return std::tan(x) - x; },
        [](double x) { return std::log(x) + x; },
        [](double x) { return std::cos(x) - std::pow(x, 3); },
        [](double x) { return std::sqrt(x) - std::cos(x); },
        [](double x) { return std::pow(x, 1.0 / 3.0) + std::pow(x, 1.0 / 5.0) - 1; },
        nxx::poly::Polynomial({-5.0, 0.0, 1.0})
    };

    // Test roots
    std::vector< double > roots {
        1.8954942670339812,
        0.6190612867359450,
        4.4934094579090642,
        0.5671432904097838,
        0.8654740331016144,
        0.6417143708728827,
        0.0700977093863724,
        2.2360679774997898
    };

    // Test brackets
    std::vector< std::pair< double, double > > brackets { { 1.0, 3.0 },
                                                          { 0.0, 1.0 },
                                                          { 4.0, 4.5 },
                                                          { 0.5, 1.0 },
                                                          { 0.5, 1.5 },
                                                          { 0.0, 1.0 },
                                                          { 0.0, 0.2 },
                                                          { 0.0, 2.5 } };
    // Test for the bisection solver
    SECTION("Bisection Solver")
    {
        size_t count = 0;
        for (const auto& func : functions) {
            INFO("Function " << count);
            auto solver = Bisection(func);
            auto root = *fsolve(solver, brackets[count], 1.0E-15);
            REQUIRE_THAT(std::abs(root - roots[count]), Catch::Matchers::WithinAbs(0.0, 0.000001));
            REQUIRE_THAT(solver.evaluate(root), Catch::Matchers::WithinAbs(0.0, 1.0E-8));
            ++count;
        }
    }

    // Test for the Ridder's solver
    SECTION("Ridders Solver Creation")
    {
        size_t count = 0;
        for (const auto& func : functions) {
            INFO("Function " << count);
            auto solver = Ridder(func);
            auto root = *fsolve(solver, brackets[count], 1.0E-15);
            REQUIRE_THAT(std::abs(root - roots[count]), Catch::Matchers::WithinAbs(0.0, 0.000001));
            REQUIRE_THAT(solver.evaluate(root), Catch::Matchers::WithinAbs(0.0, 1.0E-8));
            ++count;
        }
    }
}
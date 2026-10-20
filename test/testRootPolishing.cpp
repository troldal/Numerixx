// ================================================================================================
// Catch2 test file for the root polishing solvers.
// ================================================================================================

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <roots/Roots.hpp>

#include <cmath>
#include <functional>
#include <vector>

TEST_CASE("nxx::roots - Polishing Solver Creation", "[roots]")
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
        [](double x) { return std::pow(x, 1.0 / 3) + std::pow(x, 1.0 / 5) - 1; },
        nxx::poly::Polynomial({-5.0, 0.0, 1.0})
    };

    // Test derivatives
    std::vector< std::function< double(double) > > derivatives {
        [](double x) { return std::cos(x) - 0.5; },
        [](double x) { return std::exp(x) - 3; },
        [](double x) { return std::pow(1.0 / std::cos(x), 2) - 1; },
        [](double x) { return 1.0 / x + 1; },
        [](double x) { return -std::sin(x) - 3 * std::pow(x, 2); },
        [](double x) { return 1.0 / (2 * std::sqrt(x)) + std::sin(x); },
        [](double x) { return 1.0 / (3 * std::pow(x, 2.0 / 3)) + 1.0 / (5 * std::pow(x, 4.0 / 5)); },
        derivativeOf(nxx::poly::Polynomial({-5.0, 0.0, 1.0}))
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
                                                          { 4.4, 4.5 },
                                                          { 0.5, 1.0 },
                                                          { 0.5, 1.5 },
                                                          { 0.0, 1.0 },
                                                          { 0.0, 0.2 },
                                                          { 0.0, 2.5} };

    // Test the Newton solver
    SECTION("Newton Solver")
    {
        size_t count = 0;
        for (const auto& func : functions) {
            INFO("Function " << count);
            auto solver = Newton(func, derivatives[count]);
            auto root = *fdfsolve(solver, (brackets[count].first + brackets[count].second)/2.0, 1.0E-15);
            REQUIRE_THAT(std::abs(root - roots[count]), Catch::Matchers::WithinAbs(0.0, 0.000001));
            REQUIRE_THAT(solver.evaluate(root), Catch::Matchers::WithinAbs(0.0, 1.0E-8));
            ++count;
        }
    }

    // Test the Discrete Newton solver
    SECTION("Discrete Newton Creation")
    {
        size_t count = 0;
        for (const auto& func : functions) {
            INFO("Function " << count);
            auto solver = DNewton(func);
            auto root = *fdfsolve(solver, (brackets[count].first + brackets[count].second)/2.0, 1.0E-15);
            REQUIRE_THAT(std::abs(root - roots[count]), Catch::Matchers::WithinAbs(0.0, 0.000001));
            REQUIRE_THAT(solver.evaluate(root), Catch::Matchers::WithinAbs(0.0, 1.0E-8));
            ++count;
        }
    }
}
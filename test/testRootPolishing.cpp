//
// Created by Kenneth Balslev on 24/02/2023.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <roots/Roots.hpp>

#include <cmath>
#include <functional>
#include <vector>

TEST_CASE("Polishing Solver Creation", "[roots]")
{
    using namespace nxx::roots;

    std::vector< std::function< double(double) > > functions {
        [](double x) { return std::sin(x) - x / 2.0; },
        [](double x) { return std::exp(x) - 3 * x; },
        [](double x) { return std::tan(x) - x; },
        [](double x) { return std::log(x) + x; },
        [](double x) { return std::cos(x) - std::pow(x, 3); },
        [](double x) { return std::sqrt(x) - std::cos(x); },
        [](double x) { return std::pow(x, 1.0 / 3) + std::pow(x, 1.0 / 5) - 1; },
    };

    std::vector< std::function< double(double) > > derivatives {
        [](double x) { return std::cos(x) - 0.5; },
        [](double x) { return std::exp(x) - 3; },
        [](double x) { return std::pow(1.0 / std::cos(x), 2) - 1; },
        [](double x) { return 1.0 / x + 1; },
        [](double x) { return -std::sin(x) - 3 * std::pow(x, 2); },
        [](double x) { return 1.0 / (2 * std::sqrt(x)) + std::sin(x); },
        [](double x) { return 1.0 / (3 * std::pow(x, 2.0 / 3)) + 1.0 / (5 * std::pow(x, 4.0 / 5)); },
    };

    std::vector< double > roots {
        1.8954942670339812, 0.6190612867359450, 4.4934094579090642, 0.5671432904097838,
        0.8654740331016144, 0.6417143708728827, 0.0700977093863724,
    };

    std::vector< std::pair< double, double > > brackets { { 1.0, 3.0 }, { 0.0, 1.0 }, { 4.0, 4.5 }, { 0.5, 1.0 },
                                                          { 0.5, 1.5 }, { 0.0, 1.0 }, { 0.0, 0.2 } };

    size_t i = GENERATE(0, 1, 2, 3, 4, 5, 6);

    SECTION("Newton Solver Creation")
    {
        REQUIRE_THAT(std::abs(Newton(functions[i], derivatives[i]).evaluate(roots[i])), Catch::Matchers::WithinAbs(0.0, 0.000001));
    }

    SECTION("DNewton Solver Creation")
    {
        REQUIRE_THAT(std::abs(DNewton(functions[i]).evaluate(roots[i])), Catch::Matchers::WithinAbs(0.0, 0.000001));
    }
}
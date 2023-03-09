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

TEST_CASE("Bracketing Solver Tests", "[roots]")
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

    std::vector< double > roots {
        1.8954942670339812, 0.6190612867359450, 4.4934094579090642, 0.5671432904097838,
        0.8654740331016144, 0.6417143708728827, 0.0700977093863724,
    };

    std::vector< std::pair< double, double > > brackets { { 1.0, 3.0 }, { 0.0, 1.0 }, { 4.0, 4.5 }, { 0.5, 1.0 },
                                                          { 0.5, 1.5 }, { 0.0, 1.0 }, { 0.0, 0.2 } };

    size_t i = GENERATE(0, 1, 2, 3, 4, 5, 6);

    SECTION("Bisection Solver")
    {
        auto solver = Bisection(functions[i]);
        REQUIRE_THAT(std::abs(solver.evaluate(roots[i])), Catch::Matchers::WithinAbs(0.0, 0.000001));
        REQUIRE(solver.bounds() == std::make_pair(0.0, 0.0));

        solver.init({ -1.0, 1.0 });
        REQUIRE(solver.bounds() == std::make_pair(-1.0, 1.0));
    }

    SECTION("Ridders Solver Creation")
    {
        auto solver = Ridders(functions[i]);
        REQUIRE_THAT(std::abs(solver.evaluate(roots[i])), Catch::Matchers::WithinAbs(0.0, 0.000001));
        REQUIRE(solver.bounds() == std::make_pair(0.0, 0.0));

        solver.init({ -1.0, 1.0 });
        REQUIRE(solver.bounds() == std::make_pair(-1.0, 1.0));
    }
}
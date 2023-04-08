//
// Created by Kenneth Balslev on 25/03/2023.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <poly/Polynomial.hpp>
#include <poly/Polyroots.hpp>

#include <cmath>
#include <functional>
#include <numbers>
#include <vector>

TEST_CASE("nxx::poly - Polynomials with real coefficients and real roots", "[polynomials]")
{
    using namespace nxx::poly;
    using namespace std::complex_literals;

    auto poly = Polynomial({6.0, -5.0, 1.0});
    auto derivative = derivativeOf(poly);

    auto evaluations = std::vector<double>{6.0, 2.0, 0.0, 0.0, 2.0, 6.0, 12.0};
    auto derivatives = std::vector<double>{-5.0, -3.0, -1.0, 1.0, 3.0, 5.0, 7.0};

    SECTION("Polynomial evaluation")
    {
        size_t counter = 0;
        for (auto elem : evaluations)
            REQUIRE(poly(static_cast<double>(counter++)) == elem);
    }

    SECTION("Polynomial derivative")
    {
        size_t counter = 0;
        for (auto elem : derivatives)
            REQUIRE(derivative(static_cast<double>(counter++)) == elem);
    }

    SECTION("Polynomial roots")
    {
        auto roots = polysolve(poly);

        REQUIRE(roots.size() == 2);
        REQUIRE(roots[0] == 2.0);
        REQUIRE(roots[1] == 3.0);
        REQUIRE_THAT(poly(roots[0]), Catch::Matchers::WithinAbs(0.0, 1.0e-12));
        REQUIRE_THAT(poly(roots[1]), Catch::Matchers::WithinAbs(0.0, 1.0e-12));
    }


}
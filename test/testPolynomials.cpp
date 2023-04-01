//
// Created by Kenneth Balslev on 25/03/2023.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <poly/Polynomial.hpp>

#include <cmath>
#include <functional>
#include <numbers>
#include <vector>

TEST_CASE("Creation of polynomials with real coefficients", "[polynomials]")
{
    using namespace nxx::poly;
    using namespace std::complex_literals;

    std::vector< std::vector< double > > coefficients { { 6.0, -5.0, 1.0 },
                                                        { -6.0, 11.0, -6.0, 1.0 },
                                                        { 24.0, -50.0, 35.0, -10.0, 1.0 },
                                                        { -120.0, 274.0, -225.0, 85.0, -15.0, 1.0 },
                                                        { 720.0, -1764.0, 1624.0, -735.0, 175.0, -21.0, 1.0 },

                                                        { 2.0, 2.0, 1.0 },
                                                        {-1.0, 0.0, 0.0, 0.0, 0.0, 1.0} };

    SECTION("Check that the coefficients are stored correctly")
    {
        for (auto const& coeff : coefficients) {
            Polynomial< double > poly(coeff);

            REQUIRE(poly.coefficients() == coeff);
            REQUIRE(poly.order() == coeff.size() - 1);
        }
    }

    std::vector< std::vector< double > > rroots {
        {2.0,3.0},
        {1.0,2.0,3.0},
        {1.0,2.0,3.0,4.0},
        {1.0,2.0,3.0,4.0,5.0},
        {1.0,2.0,3.0,4.0,5.0,6.0},
    };

    std::vector< std::vector< std::complex<double> > > croots {
        {-1.0+1i,-1.0-1i},
        {-0.809017 +0.5878i, -0.809017 -0.5878i, +0.309017 +0.9510565i, +0.309017 -0.9510565i, 1.0}
    };

    SECTION("Check Polynomial evaluation with real roots") {
        for (int i = 0; i < 5; ++i) {
            Polynomial< double > poly(coefficients[i]);
                for (auto const& root : rroots[i]) {
                    REQUIRE_THAT(*poly(root), Catch::Matchers::WithinAbs(0.0, 1.0e-4));
                }
        }
    }

    SECTION("Check Polynomial evaluation with complex roots")
    {
        for (int i = 5; i < 7; ++i) {
                Polynomial< double > poly(coefficients[i]);
                for (auto const& root : croots[i - 5]) {
                    std::stringstream ss;
                    ss << "Root = " << root;
                    INFO(ss.str());
                    REQUIRE_THAT(abs(*poly(root)), Catch::Matchers::WithinAbs(0.0, 1.0e-4));
                }
        }
    }
}
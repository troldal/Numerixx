//
// Created by Kenneth Balslev on 25/03/2023.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <Poly.hpp>


#include <cmath>
#include <functional>
#include <numbers>
#include <vector>
#include <deque>

//TEST_CASE("nxx::poly - Polynomials with real coefficients and real roots", "[polynomials]")
//{
//    using namespace nxx::poly;
//    using namespace std::complex_literals;
//
//    auto poly = Polynomial({6.0, -5.0, 1.0});
//    auto derivative = derivativeOf(poly);
//
//    auto evaluations = std::vector<double>{6.0, 2.0, 0.0, 0.0, 2.0, 6.0, 12.0};
//    auto derivatives = std::vector<double>{-5.0, -3.0, -1.0, 1.0, 3.0, 5.0, 7.0};
//
//    SECTION("Polynomial evaluation")
//    {
//        size_t counter = 0;
//        for (auto elem : evaluations)
//            REQUIRE(poly(static_cast<double>(counter++)) == elem);
//    }
//
//    SECTION("Polynomial derivative")
//    {
//        size_t counter = 0;
//        for (auto elem : derivatives)
//            REQUIRE(derivative(static_cast<double>(counter++)) == elem);
//    }
//
//    SECTION("Polynomial roots")
//    {
//        auto roots = polysolve(poly);
//
//        REQUIRE(roots.size() == 2);
//        REQUIRE(roots[0] == 2.0);
//        REQUIRE(roots[1] == 3.0);
//        REQUIRE_THAT(poly(roots[0]), Catch::Matchers::WithinAbs(0.0, 1.0e-12));
//        REQUIRE_THAT(poly(roots[1]), Catch::Matchers::WithinAbs(0.0, 1.0e-12));
//    }
//
//
//}

TEST_CASE("Polynomial class tests", "[Polynomial]")
{
    using namespace nxx::poly;
    using namespace nxx::error;
    using namespace std::complex_literals;

    SECTION("Constructor tests")
    {
        // Test using std::complex<double> as coefficients
        Polynomial<std::complex<double>> p1({{1, 2}, {2, -1}, {-3, -4}});
        REQUIRE(p1.order() == 2);
        REQUIRE(p1.coefficients() == std::vector<std::complex<double>>{{1, 2}, {2, -1}, {-3, -4}});
        REQUIRE(p1.coefficients<std::deque<std::complex<double>>>() == std::deque<std::complex<double>>{{1, 2}, {2, -1}, {-3, -4}});
        auto c1 = std::vector<std::complex<double>>(p1.begin(), p1.end());
        REQUIRE(c1 == std::vector<std::complex<double>>{{1, 2}, {2, -1}, {-3, -4}});



        // Test using double as coefficients
        Polynomial p2({6.0, -5.0, 1.0});
        REQUIRE(p2.order() == 2);
        REQUIRE(p2.coefficients() == std::vector<double>{6.0, -5.0, 1.0});
        REQUIRE(p2.coefficients<std::deque<double>>() == std::deque<double>{6.0, -5.0, 1.0});
        auto c2 = std::vector<double>(p2.begin(), p2.end());
        REQUIRE(c2 == std::vector<double>{6.0, -5.0, 1.0});
    }

    SECTION("Evaluation tests")
    {
        Polynomial<std::complex<double>> p1({{1, 2}, {2, -1}, {-3, -4}});

        // Test using double as argument
        REQUIRE_THAT(p1(0).real(), Catch::Matchers::WithinAbs(1.0, 1.0e-12));
        REQUIRE_THAT(p1(0).imag(), Catch::Matchers::WithinAbs(2.0, 1.0e-12));
        REQUIRE_THAT(p1(1).real(), Catch::Matchers::WithinAbs(0.0, 1.0e-12));
        REQUIRE_THAT(p1(1).imag(), Catch::Matchers::WithinAbs(-3.0, 1.0e-12));
        REQUIRE_THAT(p1(-1).real(), Catch::Matchers::WithinAbs(-4.0, 1.0e-12));
        REQUIRE_THAT(p1(-1).imag(), Catch::Matchers::WithinAbs(-1.0, 1.0e-12));

        // Test using std::complex<double> as argument
        REQUIRE_THAT(p1(0.0 + 1.0i).real(), Catch::Matchers::WithinAbs(5.0, 1.0e-12));
        REQUIRE_THAT(p1(0.0 + 1.0i).imag(), Catch::Matchers::WithinAbs(8.0, 1.0e-12));
        REQUIRE_THAT(p1(1.0 + 1.0i).real(), Catch::Matchers::WithinAbs(12.0, 1.0e-12));
        REQUIRE_THAT(p1(1.0 + 1.0i).imag(), Catch::Matchers::WithinAbs(-3.0, 1.0e-12));
        REQUIRE_THAT(p1(-1.0 + 1.0i).real(), Catch::Matchers::WithinAbs(-8.0, 1.0e-12));
        REQUIRE_THAT(p1(-1.0 + 1.0i).imag(), Catch::Matchers::WithinAbs(11.0, 1.0e-12));

        Polynomial p2({6.0, -5.0, 1.0});

        // Test using double as argument
        REQUIRE_THAT(p2(0), Catch::Matchers::WithinAbs(6.0, 1.0e-12));
        REQUIRE_THAT(p2(1), Catch::Matchers::WithinAbs(2.0, 1.0e-12));
        REQUIRE_THAT(p2(-1), Catch::Matchers::WithinAbs(12.0, 1.0e-12));

        // Test using std::complex<double> as argument
        REQUIRE_THAT(p2(0.0 + 1.0i).real(), Catch::Matchers::WithinAbs(5.0, 1.0e-12));
        REQUIRE_THAT(p2(0.0 + 1.0i).imag(), Catch::Matchers::WithinAbs(-5.0, 1.0e-12));
        REQUIRE_THAT(p2(1.0 + 1.0i).real(), Catch::Matchers::WithinAbs(1.0, 1.0e-12));
        REQUIRE_THAT(p2(1.0 + 1.0i).imag(), Catch::Matchers::WithinAbs(-3.0, 1.0e-12));
        REQUIRE_THAT(p2(-1.0 + 1.0i).real(), Catch::Matchers::WithinAbs(11.0, 1.0e-12));
        REQUIRE_THAT(p2(-1.0 + 1.0i).imag(), Catch::Matchers::WithinAbs(-7.0, 1.0e-12));
    }

    SECTION("Arithmetic Operations tests")
    {
        Polynomial<double> p1({1, 2, 3});
        Polynomial<double> p2({4, 5, 6});
        Polynomial<double> p3({5, 6, 7, 8});

        auto p4 = p1 + p2;
        REQUIRE(p4.coefficients() == std::vector<double>{5, 7, 9});
        p4 = p2;
        p4 += p3;
        REQUIRE(p4.coefficients() == std::vector<double>{9, 11, 13, 8});

        auto p5 = p1 - p2;
        REQUIRE(p5.coefficients() == std::vector<double>{-3, -3, -3});
        p5 = p2;
        p5 -= p3;
        REQUIRE(p5.coefficients() == std::vector<double>{-1, -1, -1, 8});

        auto p6 = p1 * p2;
        REQUIRE(p6.coefficients() == std::vector<double>{4, 13, 28, 27, 18});
        p6 = p2;
        p6 *= p3;
        REQUIRE(p6.coefficients() == std::vector<double>{20, 49, 88, 103, 82, 48});

        auto p7 = p1 / p2;
        REQUIRE(p7.coefficients() == std::vector<double>{0.5});
        p7 = p1;
        p7 /= p2;
        REQUIRE(p7.coefficients() == std::vector<double>{0.5});

        auto p8 = p1 % p2;
        REQUIRE(p8.coefficients() == std::vector<double>{-1, -0.5});

        // todo: Test with other comlex types
        // todo: Test with cross-type operations
    }

    SECTION("Order and Coefficient Tests")
    {
        // Test using std::complex<double> as coefficients
        Polynomial<std::complex<double> > p1({1.0+1i, 2.0+1i, 3.0+1i, 0.0+0i, 0.0+0i});
        REQUIRE(p1.order() == 2);
        REQUIRE(p1.coefficients() == std::vector<std::complex<double>>{1.0+1i, 2.0+1i, 3.0+1i});

        // Test using double as coefficients
        Polynomial p2({1.0, 2.0, 3.0, 0.0, 0.0});
        REQUIRE(p2.order() == 2);
        REQUIRE(p2.coefficients() == std::vector<double>{1.0, 2.0, 3.0});

        // Test with zero coefficients (a constant term of 0.0 is added implicitly)
        Polynomial p3({});
        REQUIRE(p3.order() == 0);
        REQUIRE(p3.coefficients() == std::vector<double>{0.0});

        // Test equality and inequality
        REQUIRE(p1 == p1);
        REQUIRE(p2 == p2);
        REQUIRE(p3 == p3);
        REQUIRE(p2 != p3);
    }

    SECTION("Derivative Tests")
    {
        Polynomial<double> p1({1, 3, 3});
        auto p2 = derivativeOf(p1);
        REQUIRE(p2.coefficients() == std::vector<double>{3, 6});

        // todo: Test with other degrees
        // todo: Test with complex coefficients
    }

    SECTION("String Representation Tests")
    {
        Polynomial<double> p1({1, 2, 3});
        REQUIRE(p1.asString() == "1 + 2x + 3x^2");
    }

    SECTION("Boundary Tests")
    {
        Polynomial<double> p1({0.0, 0.0, 0.0});
        REQUIRE(p1.order() == 0);
        REQUIRE(p1.coefficients() == std::vector<double>{0.0});
    }

    SECTION("Error Handling Tests")
    {
        Polynomial<double> p1({1, 2, 3});
        Polynomial<double> p2({0.0, 0.0, 0.0});

        REQUIRE_THROWS_AS(p1 / p2, PolynomialError);
    }
}

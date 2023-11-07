//
// Created by Kenneth Balslev on 25/03/2023.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <Poly.hpp>

#include <cmath>
#include <deque>
#include <sstream>
#include <vector>

constexpr double EPS = 1E-5;

TEST_CASE("Polynomial class tests", "[Polynomial]")
{
    using namespace nxx::poly;
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

        // Test using double as coefficients
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

        // Test using std::complex<double> as coefficients
        Polynomial<std::complex<double>> c1({1.0+0i, 2.0+0i, 3.0+0i});
        Polynomial<std::complex<double>> c2({4.0+0i, 5.0+0i, 6.0+0i});
        Polynomial<std::complex<double>> c3({5.0+0i, 6.0+0i, 7.0+0i, 8.0+0i});

        auto c4 = c1 + c2;
        REQUIRE(c4.coefficients() == std::vector<std::complex<double>>{5.0+0i, 7.0+0i, 9.0+0i});
        c4 = c2;
        c4 += c3;
        REQUIRE(c4.coefficients() == std::vector<std::complex<double>>{9.0+0i, 11.0+0i, 13.0+0i, 8.0+0i});

        auto c5 = c1 - c2;
        REQUIRE(c5.coefficients() == std::vector<std::complex<double>>{-3.0+0i, -3.0+0i, -3.0+0i});
        c5 = c2;
        c5 -= c3;
        REQUIRE(c5.coefficients() == std::vector<std::complex<double>>{-1.0+0i, -1.0+0i, -1.0+0i, 8.0+0i});

        auto c6 = c1 * c2;
        REQUIRE(c6.coefficients() == std::vector<std::complex<double>>{4.0+0i, 13.0+0i, 28.0+0i, 27.0+0i, 18.0+0i});
        c6 = c2;
        c6 *= c3;
        REQUIRE(c6.coefficients() == std::vector<std::complex<double>>{20.0+0i, 49.0+0i, 88.0+0i, 103.0+0i, 82.0+0i, 48.0+0i});

        auto c7 = c1 / c2;
        REQUIRE(c7.coefficients() == std::vector<std::complex<double>>{0.5+0i});
        c7 = c1;
        c7 /= c2;
        REQUIRE(c7.coefficients() == std::vector<std::complex<double>>{0.5+0i});

        auto c8 = c1 % c2;
        REQUIRE(c8.coefficients() == std::vector<std::complex<double>>{-1.0+0i, -0.5+0i});

        // Test with cross-type operations
        auto t4 = p1 + c2;
        REQUIRE(t4.coefficients() == std::vector<std::complex<double>>{5.0+0i, 7.0+0i, 9.0+0i});

        auto t5 = p1 - c2;
        REQUIRE(t5.coefficients() == std::vector<std::complex<double>>{-3.0+0i, -3.0+0i, -3.0+0i});

        auto t6 = p1 * c2;
        REQUIRE(t6.coefficients() == std::vector<std::complex<double>>{4.0+0i, 13.0+0i, 28.0+0i, 27.0+0i, 18.0+0i});

        auto t7 = p1 / c2;
        REQUIRE(t7.coefficients() == std::vector<std::complex<double>>{0.5+0i});

        auto t8 = p1 % c2;
        REQUIRE(t8.coefficients() == std::vector<std::complex<double>>{-1.0+0i, -0.5+0i});

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
        auto d1 = derivativeOf(p1);
        REQUIRE(d1.coefficients() == std::vector<double>{3, 6});

        // test with 3rd degree polynomial
        Polynomial<double> p2({1, 2, 3, 4});
        auto d2 = derivativeOf(p2);
        REQUIRE(d2.coefficients() == std::vector<double>{2, 6, 12});

        // test with 2nd degree polynomial with complex coefficients
        Polynomial<std::complex<double>> p3({{1.0+0i, 3.0+0i, 3.0+0i}});
        auto d3 = derivativeOf(p3);
        REQUIRE(d3.coefficients() == std::vector<std::complex<double>>{{3.0+0i, 6.0+0i}});

        // test with 3rd degree polynomial with complex coefficients
        Polynomial<std::complex<double>> p4({{1.0+0i, 2.0+0i, 3.0+0i, 4.0+0i}});
        auto d4 = derivativeOf(p4);
        REQUIRE(d4.coefficients() == std::vector<std::complex<double>>{{2.0+0i, 6.0+0i, 12.0+0i}});
    }

    SECTION("String Representation Tests")
    {
        Polynomial<double> p1({1, 2, 3});
        REQUIRE(to_string(p1) == "1 + 2x + 3x^2");
    }

    SECTION("Boundary Tests")
    {
        Polynomial<double> p1({0.0, 0.0, 0.0});
        REQUIRE(p1.order() == 0);
        REQUIRE(p1.coefficients() == std::vector<double>{0.0});
    }

    //    SECTION("Error Handling Tests")
    //    {
    //        Polynomial<double> p1({1, 2, 3});
    //        Polynomial<double> p2({0.0, 0.0, 0.0});
    //
    //        REQUIRE_THROWS_AS(p1 / p2, PolynomialError);
    //    }

    SECTION("Other Tests")
    {
        Polynomial p1({1.0, 0.5, 0.3});
        REQUIRE_THAT(p1(0.5), Catch::Matchers::WithinAbs(1.0 + 0.5 * 0.5 + 0.3 * 0.5 * 0.5, 1.0E-12));

        Polynomial p2({1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1});
        REQUIRE_THAT(p2(1.0), Catch::Matchers::WithinAbs(1.0, 1.0E-12));

        Polynomial p3({0.3});
        REQUIRE_THAT(p3(0.75+1.2i).real(), Catch::Matchers::WithinAbs(0.3, 1.0E-12));
        REQUIRE_THAT(p3(0.75+1.2i).imag(), Catch::Matchers::WithinAbs(0.0, 1.0E-12));

        Polynomial p4({2.1, -1.34, 0.76, 0.45});
        REQUIRE_THAT(p4(0.49+0.95i).real(), Catch::Matchers::WithinAbs(0.3959143, 1.0E-5));
        REQUIRE_THAT(p4(0.49+0.95i).imag(), Catch::Matchers::WithinAbs(-0.643330, 1.0E-5));

        Polynomial p5({-2.31+0.44i, 4.21-3.19i, 0.93+1.04i, -0.42+0.68i});
        REQUIRE_THAT(p5(0.49+0.95i).real(), Catch::Matchers::WithinAbs(1.8246201, 1.0E-5));
        REQUIRE_THAT(p5(0.49+0.95i).imag(), Catch::Matchers::WithinAbs(2.30389412, 1.0E-5));
    }
}

TEST_CASE("Polynomial roots tests", "[Polynomial]")
{
    using namespace nxx::poly;
    using namespace std::complex_literals;

    SECTION("Quadratics")
    {
        Polynomial p1({26.0, -20.0, 4.0});
        auto rroots1 = polysolve(p1);
        REQUIRE(rroots1.value().empty());
        auto croots1 = polysolve<std::complex<double> >(p1);
        REQUIRE(croots1.value().size() == 2);
        REQUIRE_THAT(croots1.value()[0].real(), Catch::Matchers::WithinAbs( 2.5, EPS));
        REQUIRE_THAT(croots1.value()[0].imag(), Catch::Matchers::WithinAbs(-0.5, EPS));
        REQUIRE_THAT(croots1.value()[1].real(), Catch::Matchers::WithinAbs( 2.5, EPS));
        REQUIRE_THAT(croots1.value()[1].imag(), Catch::Matchers::WithinAbs( 0.5, EPS));

        Polynomial p2({25.0, -20.0, 4.0});
        auto rroots2 = polysolve(p2);
        REQUIRE(rroots2.value().size() == 2);
        REQUIRE_THAT(rroots2.value()[0], Catch::Matchers::WithinAbs(2.5, EPS));
        REQUIRE_THAT(rroots2.value()[1], Catch::Matchers::WithinAbs(2.5, EPS));

        auto croots2 = polysolve<std::complex<double> >(p2);
        REQUIRE(croots2.value().size() == 2);
        REQUIRE_THAT(croots2.value()[0].real(), Catch::Matchers::WithinAbs(2.5, EPS));
        REQUIRE_THAT(croots2.value()[0].imag(), Catch::Matchers::WithinAbs(0.0, EPS));
        REQUIRE_THAT(croots2.value()[1].real(), Catch::Matchers::WithinAbs(2.5, EPS));
        REQUIRE_THAT(croots2.value()[1].imag(), Catch::Matchers::WithinAbs(0.0, EPS));

        Polynomial p3({21.0, -20.0, 4.0});
        auto rroots3 = polysolve(p3);
        REQUIRE(rroots3.value().size() == 2);
        REQUIRE_THAT(rroots3.value()[0], Catch::Matchers::WithinAbs(1.5, EPS));
        REQUIRE_THAT(rroots3.value()[1], Catch::Matchers::WithinAbs(3.5, EPS));
        auto croots3 = polysolve<std::complex<double> >(p3);
        REQUIRE(croots3.value().size() == 2);
        REQUIRE_THAT(croots3.value()[0].real(), Catch::Matchers::WithinAbs(1.5, EPS));
        REQUIRE_THAT(croots3.value()[0].imag(), Catch::Matchers::WithinAbs(0.0, EPS));
        REQUIRE_THAT(croots3.value()[1].real(), Catch::Matchers::WithinAbs(3.5, EPS));
        REQUIRE_THAT(croots3.value()[1].imag(), Catch::Matchers::WithinAbs(0.0, EPS));

        Polynomial p4({0.0, 7.0, 4.0});
        auto rroots4 = polysolve(p4);
        REQUIRE(rroots4.value().size() == 2);
        REQUIRE_THAT(rroots4.value()[0], Catch::Matchers::WithinAbs(-1.75, EPS));
        REQUIRE_THAT(rroots4.value()[1], Catch::Matchers::WithinAbs( 0.00, EPS));
        auto croots4 = polysolve<std::complex<double> >(p4);
        REQUIRE(croots4.value().size() == 2);
        REQUIRE_THAT(croots4.value()[0].real(), Catch::Matchers::WithinAbs(-1.75, EPS));
        REQUIRE_THAT(croots4.value()[0].imag(), Catch::Matchers::WithinAbs( 0.00, EPS));
        REQUIRE_THAT(croots4.value()[1].real(), Catch::Matchers::WithinAbs( 0.00, EPS));
        REQUIRE_THAT(croots4.value()[1].imag(), Catch::Matchers::WithinAbs( 0.00, EPS));

        Polynomial p5({-20.0, 0.0, 5.0});
        auto rroots5 = polysolve(p5);
        REQUIRE(rroots5.value().size() == 2);
        REQUIRE_THAT(rroots5.value()[0], Catch::Matchers::WithinAbs(-2.0, EPS));
        REQUIRE_THAT(rroots5.value()[1], Catch::Matchers::WithinAbs( 2.0, EPS));
        auto croots5 = polysolve<std::complex<double> >(p5);
        REQUIRE(croots5.value().size() == 2);
        REQUIRE_THAT(croots5.value()[0].real(), Catch::Matchers::WithinAbs(-2.0, EPS));
        REQUIRE_THAT(croots5.value()[0].imag(), Catch::Matchers::WithinAbs( 0.0, EPS));
        REQUIRE_THAT(croots5.value()[1].real(), Catch::Matchers::WithinAbs( 2.0, EPS));
        REQUIRE_THAT(croots5.value()[1].imag(), Catch::Matchers::WithinAbs( 0.0, EPS));

        Polynomial p6({20.0, 0.0, 5.0});
        auto rroots6 = polysolve(p6);
        REQUIRE(rroots6.value().empty());
        auto croots6 = polysolve<std::complex<double> >(p6);
        REQUIRE(croots6.value().size() == 2);
        REQUIRE_THAT(croots6.value()[0].real(), Catch::Matchers::WithinAbs( 0.0, EPS));
        REQUIRE_THAT(croots6.value()[0].imag(), Catch::Matchers::WithinAbs(-2.0, EPS));
        REQUIRE_THAT(croots6.value()[1].real(), Catch::Matchers::WithinAbs( 0.0, EPS));
        REQUIRE_THAT(croots6.value()[1].imag(), Catch::Matchers::WithinAbs( 2.0, EPS));

        Polynomial p7({-21.0, 3.0, 0.0});
        auto rroots7 = polysolve(p7);
        REQUIRE(rroots7.value().size() == 1);
        REQUIRE_THAT(rroots7.value()[0], Catch::Matchers::WithinAbs(7.0, EPS));
        auto croots7 = polysolve<std::complex<double> >(p7);
        REQUIRE(croots7.value().size() == 1);
        REQUIRE_THAT(croots7.value()[0].real(), Catch::Matchers::WithinAbs(7.0, EPS));
        REQUIRE_THAT(croots7.value()[0].imag(), Catch::Matchers::WithinAbs(0.0, EPS));

    }

    SECTION("Cubics")
    {
        Polynomial p1({-27.0, 0.0, 0.0, 1.0});
        auto rroots1 = polysolve(p1);
        REQUIRE(rroots1.value().size() == 1);
        REQUIRE_THAT(rroots1.value()[0], Catch::Matchers::WithinAbs(3.0, EPS));
        auto croots1 = polysolve<std::complex<double> >(p1);
        REQUIRE(croots1.value().size() == 3);
        REQUIRE_THAT(croots1.value()[0].real(), Catch::Matchers::WithinAbs(-1.5, EPS));
        REQUIRE_THAT(croots1.value()[0].imag(), Catch::Matchers::WithinAbs(-1.5 * std::sqrt(3.0), EPS));
        REQUIRE_THAT(croots1.value()[1].real(), Catch::Matchers::WithinAbs(-1.5, EPS));
        REQUIRE_THAT(croots1.value()[1].imag(), Catch::Matchers::WithinAbs(1.5 * std::sqrt(3.0), EPS));
        REQUIRE_THAT(croots1.value()[2].real(), Catch::Matchers::WithinAbs(3.0, EPS));
        REQUIRE_THAT(croots1.value()[2].imag(), Catch::Matchers::WithinAbs(0.0, EPS));

        Polynomial p2({39.0, 1.0, -1.0, 1.0});
        auto rroots2 = polysolve(p2);
        REQUIRE(rroots2.value().size() == 1);
        REQUIRE_THAT(rroots2.value()[0], Catch::Matchers::WithinAbs(-3.0, EPS));
        auto croots2 = polysolve<std::complex<double> >(p2);
        REQUIRE(croots2.value().size() == 3);
        REQUIRE_THAT(croots2.value()[0].real(), Catch::Matchers::WithinAbs(-3.0, EPS));
        REQUIRE_THAT(croots2.value()[0].imag(), Catch::Matchers::WithinAbs(0.0, EPS));
        REQUIRE_THAT(croots2.value()[1].real(), Catch::Matchers::WithinAbs(2.0, EPS));
        REQUIRE_THAT(croots2.value()[1].imag(), Catch::Matchers::WithinAbs(-3.0, EPS));
        REQUIRE_THAT(croots2.value()[2].real(), Catch::Matchers::WithinAbs(2.0, EPS));
        REQUIRE_THAT(croots2.value()[2].imag(), Catch::Matchers::WithinAbs(3.0, EPS));

        Polynomial p3({-4913.0, 867.0, -51.0, 1.0});
        auto rroots3 = polysolve(p3);
        REQUIRE(rroots3.value().size() == 3);
        REQUIRE_THAT(rroots3.value()[0], Catch::Matchers::WithinAbs(17.0, EPS));
        REQUIRE_THAT(rroots3.value()[1], Catch::Matchers::WithinAbs(17.0, EPS));
        REQUIRE_THAT(rroots3.value()[2], Catch::Matchers::WithinAbs(17.0, EPS));
        auto croots3 = polysolve<std::complex<double> >(p3);
        REQUIRE(croots3.value().size() == 3);
        REQUIRE_THAT(croots3.value()[0].real(), Catch::Matchers::WithinAbs(17.0, EPS));
        REQUIRE_THAT(croots3.value()[0].imag(), Catch::Matchers::WithinAbs(0.0, EPS));
        REQUIRE_THAT(croots3.value()[1].real(), Catch::Matchers::WithinAbs(17.0, EPS));
        REQUIRE_THAT(croots3.value()[1].imag(), Catch::Matchers::WithinAbs(0.0, EPS));
        REQUIRE_THAT(croots3.value()[2].real(), Catch::Matchers::WithinAbs(17.0, EPS));
        REQUIRE_THAT(croots3.value()[2].imag(), Catch::Matchers::WithinAbs(0.0, EPS));

        Polynomial p4({-6647.0, 1071.0, -57.0, 1.0});
        auto rroots4 = polysolve(p4);
        REQUIRE(rroots4.value().size() == 3);
        REQUIRE_THAT(rroots4.value()[0], Catch::Matchers::WithinAbs(17.0, EPS));
        REQUIRE_THAT(rroots4.value()[1], Catch::Matchers::WithinAbs(17.0, EPS));
        REQUIRE_THAT(rroots4.value()[2], Catch::Matchers::WithinAbs(23.0, EPS));
        auto croots4 = polysolve<std::complex<double> >(p4);
        REQUIRE(croots4.value().size() == 3);
        REQUIRE_THAT(croots4.value()[0].real(), Catch::Matchers::WithinAbs(17.0, EPS));
        REQUIRE_THAT(croots4.value()[0].imag(), Catch::Matchers::WithinAbs(0.0, EPS));
        REQUIRE_THAT(croots4.value()[1].real(), Catch::Matchers::WithinAbs(17.0, EPS));
        REQUIRE_THAT(croots4.value()[1].imag(), Catch::Matchers::WithinAbs(0.0, EPS));
        REQUIRE_THAT(croots4.value()[2].real(), Catch::Matchers::WithinAbs(23.0, EPS));
        REQUIRE_THAT(croots4.value()[2].imag(), Catch::Matchers::WithinAbs(0.0, EPS));

        Polynomial p5({6647.0, -493.0, -11.0, 1.0});
        auto rroots5 = polysolve(p5);
        REQUIRE(rroots5.value().size() == 3);
        REQUIRE_THAT(rroots5.value()[0], Catch::Matchers::WithinAbs(-23.0, EPS));
        REQUIRE_THAT(rroots5.value()[1], Catch::Matchers::WithinAbs(17.0, EPS));
        REQUIRE_THAT(rroots5.value()[2], Catch::Matchers::WithinAbs(17.0, EPS));
        auto croots5 = polysolve<std::complex<double> >(p5);
        REQUIRE(croots5.value().size() == 3);
        REQUIRE_THAT(croots5.value()[0].real(), Catch::Matchers::WithinAbs(-23.0, EPS));
        REQUIRE_THAT(croots5.value()[0].imag(), Catch::Matchers::WithinAbs(0.0, EPS));
        REQUIRE_THAT(croots5.value()[1].real(), Catch::Matchers::WithinAbs(17.0, EPS));
        REQUIRE_THAT(croots5.value()[1].imag(), Catch::Matchers::WithinAbs(0.0, EPS));
        REQUIRE_THAT(croots5.value()[2].real(), Catch::Matchers::WithinAbs(17.0, EPS));
        REQUIRE_THAT(croots5.value()[2].imag(), Catch::Matchers::WithinAbs(0.0, EPS));

        Polynomial p6({-50065.0, 5087.0, -143.0, 1.0});
        auto rroots6 = polysolve(p6);
        REQUIRE(rroots6.value().size() == 3);
        REQUIRE_THAT(rroots6.value()[0], Catch::Matchers::WithinAbs(17.0, EPS));
        REQUIRE_THAT(rroots6.value()[1], Catch::Matchers::WithinAbs(31.0, EPS));
        REQUIRE_THAT(rroots6.value()[2], Catch::Matchers::WithinAbs(95.0, EPS));
        auto croots6 = polysolve<std::complex<double> >(p6);
        REQUIRE(croots6.value().size() == 3);
        REQUIRE_THAT(croots6.value()[0].real(), Catch::Matchers::WithinAbs(17.0, EPS));
        REQUIRE_THAT(croots6.value()[0].imag(), Catch::Matchers::WithinAbs(0.0, EPS));
        REQUIRE_THAT(croots6.value()[1].real(), Catch::Matchers::WithinAbs(31.0, EPS));
        REQUIRE_THAT(croots6.value()[1].imag(), Catch::Matchers::WithinAbs(0.0, EPS));
        REQUIRE_THAT(croots6.value()[2].real(), Catch::Matchers::WithinAbs(95.0, EPS));
        REQUIRE_THAT(croots6.value()[2].imag(), Catch::Matchers::WithinAbs(0.0, EPS));
    }

    SECTION("Higher-order")
    {
        Polynomial p1({-120, 274, -225, 85, -15, 1.0});
        auto rroots1 = polysolve(p1);
        REQUIRE(rroots1.value().size() == 5);
        REQUIRE_THAT(rroots1.value()[0], Catch::Matchers::WithinAbs(1.0, EPS));
        REQUIRE_THAT(rroots1.value()[1], Catch::Matchers::WithinAbs(2.0, EPS));
        REQUIRE_THAT(rroots1.value()[2], Catch::Matchers::WithinAbs(3.0, EPS));
        REQUIRE_THAT(rroots1.value()[3], Catch::Matchers::WithinAbs(4.0, EPS));
        REQUIRE_THAT(rroots1.value()[4], Catch::Matchers::WithinAbs(5.0, EPS));
        auto croots1 = polysolve<std::complex<double> >(p1);
        REQUIRE(croots1.value().size() == 5);
        REQUIRE_THAT(croots1.value()[0].real(), Catch::Matchers::WithinAbs(1.0, EPS));
        REQUIRE_THAT(croots1.value()[0].imag(), Catch::Matchers::WithinAbs(0.0, EPS));
        REQUIRE_THAT(croots1.value()[1].real(), Catch::Matchers::WithinAbs(2.0, EPS));
        REQUIRE_THAT(croots1.value()[1].imag(), Catch::Matchers::WithinAbs(0.0, EPS));
        REQUIRE_THAT(croots1.value()[2].real(), Catch::Matchers::WithinAbs(3.0, EPS));
        REQUIRE_THAT(croots1.value()[2].imag(), Catch::Matchers::WithinAbs(0.0, EPS));
        REQUIRE_THAT(croots1.value()[3].real(), Catch::Matchers::WithinAbs(4.0, EPS));
        REQUIRE_THAT(croots1.value()[3].imag(), Catch::Matchers::WithinAbs(0.0, EPS));
        REQUIRE_THAT(croots1.value()[4].real(), Catch::Matchers::WithinAbs(5.0, EPS));
        REQUIRE_THAT(croots1.value()[4].imag(), Catch::Matchers::WithinAbs(0.0, EPS));

        Polynomial p2({1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0});
        auto rroots2 = polysolve(p2);
        REQUIRE(rroots2.value().empty());
        auto croots2 = polysolve<std::complex<double> >(p2);
        REQUIRE(croots2.value().size() == 8);
        REQUIRE_THAT(croots2.value()[0].real(), Catch::Matchers::WithinAbs(-std::sqrt(3.0)/2.0, EPS));
        REQUIRE_THAT(croots2.value()[0].imag(), Catch::Matchers::WithinAbs(-0.5, EPS));
        REQUIRE_THAT(croots2.value()[1].real(), Catch::Matchers::WithinAbs(-std::sqrt(3.0)/2.0, EPS));
        REQUIRE_THAT(croots2.value()[1].imag(), Catch::Matchers::WithinAbs(0.5, EPS));
        REQUIRE_THAT(croots2.value()[2].real(), Catch::Matchers::WithinAbs(-0.5, EPS));
        REQUIRE_THAT(croots2.value()[2].imag(), Catch::Matchers::WithinAbs(-std::sqrt(3.0)/2.0, EPS));
        REQUIRE_THAT(croots2.value()[3].real(), Catch::Matchers::WithinAbs(-0.5, EPS));
        REQUIRE_THAT(croots2.value()[3].imag(), Catch::Matchers::WithinAbs(std::sqrt(3.0)/2.0, EPS));
        REQUIRE_THAT(croots2.value()[4].real(), Catch::Matchers::WithinAbs(0.5, EPS));
        REQUIRE_THAT(croots2.value()[4].imag(), Catch::Matchers::WithinAbs(-std::sqrt(3.0)/2.0, EPS));
        REQUIRE_THAT(croots2.value()[5].real(), Catch::Matchers::WithinAbs(0.5, EPS));
        REQUIRE_THAT(croots2.value()[5].imag(), Catch::Matchers::WithinAbs(std::sqrt(3.0)/2.0, EPS));
        REQUIRE_THAT(croots2.value()[6].real(), Catch::Matchers::WithinAbs(std::sqrt(3.0)/2.0, EPS));
        REQUIRE_THAT(croots2.value()[6].imag(), Catch::Matchers::WithinAbs(-0.5, EPS));
        REQUIRE_THAT(croots2.value()[7].real(), Catch::Matchers::WithinAbs(std::sqrt(3.0)/2.0, EPS));
        REQUIRE_THAT(croots2.value()[7].imag(), Catch::Matchers::WithinAbs(0.5, EPS));

        Polynomial p3({ 32.0, -48.0, -8.0, 28.0, -8.0, 16.0, -16.0, 12.0, -16.0, 6.0, 10.0, -17.0, 10.0, 2.0, -4.0, 1.0});
        auto rroots3 = polysolve(p3);
        REQUIRE(rroots3.value().size() == 7);
        REQUIRE_THAT(rroots3.value()[0], Catch::Matchers::WithinAbs(-1.6078107423472359, EPS));
        REQUIRE_THAT(rroots3.value()[1], Catch::Matchers::WithinAbs(-1.3066982484920768, EPS));
        REQUIRE_THAT(rroots3.value()[2], Catch::Matchers::WithinAbs(-1.0, EPS));
        REQUIRE_THAT(rroots3.value()[3], Catch::Matchers::WithinAbs(1.0, EPS));
        REQUIRE_THAT(rroots3.value()[4], Catch::Matchers::WithinAbs(1.0, EPS));
        REQUIRE_THAT(rroots3.value()[5], Catch::Matchers::WithinAbs(2.0, EPS));
        REQUIRE_THAT(rroots3.value()[6], Catch::Matchers::WithinAbs(2.0, EPS));
        auto croots3 = polysolve<std::complex<double> >(p3);
        REQUIRE(croots3.value().size() == 15);
        REQUIRE_THAT(croots3.value()[0].real(), Catch::Matchers::WithinAbs(-1.6078107423472359, EPS));
        REQUIRE_THAT(croots3.value()[0].imag(), Catch::Matchers::WithinAbs(0.0, EPS));
        REQUIRE_THAT(croots3.value()[1].real(), Catch::Matchers::WithinAbs(-1.3066982484920768, EPS));
        REQUIRE_THAT(croots3.value()[1].imag(), Catch::Matchers::WithinAbs(0.0, EPS));
        REQUIRE_THAT(croots3.value()[2].real(), Catch::Matchers::WithinAbs(-1.0, EPS));
        REQUIRE_THAT(croots3.value()[2].imag(), Catch::Matchers::WithinAbs(0.0, EPS));
        REQUIRE_THAT(croots3.value()[3].real(), Catch::Matchers::WithinAbs(-0.65893856175240950, EPS));
        REQUIRE_THAT(croots3.value()[3].imag(), Catch::Matchers::WithinAbs(-0.83459757287426684, EPS));
        REQUIRE_THAT(croots3.value()[4].real(), Catch::Matchers::WithinAbs(-0.65893856175240950, EPS));
        REQUIRE_THAT(croots3.value()[4].imag(), Catch::Matchers::WithinAbs(0.83459757287426684, EPS));
        REQUIRE_THAT(croots3.value()[5].real(), Catch::Matchers::WithinAbs(-0.070891117403341281, EPS));
        REQUIRE_THAT(croots3.value()[5].imag(), Catch::Matchers::WithinAbs(-1.1359249087587791, EPS));
        REQUIRE_THAT(croots3.value()[6].real(), Catch::Matchers::WithinAbs(-0.070891117403341281, EPS));
        REQUIRE_THAT(croots3.value()[6].imag(), Catch::Matchers::WithinAbs(1.1359249087587791, EPS));
        REQUIRE_THAT(croots3.value()[7].real(), Catch::Matchers::WithinAbs(0.57284747839410854, EPS));
        REQUIRE_THAT(croots3.value()[7].imag(), Catch::Matchers::WithinAbs(-1.1987808988289705, EPS));
        REQUIRE_THAT(croots3.value()[8].real(), Catch::Matchers::WithinAbs(0.57284747839410854, EPS));
        REQUIRE_THAT(croots3.value()[8].imag(), Catch::Matchers::WithinAbs(1.1987808988289705, EPS));
        REQUIRE_THAT(croots3.value()[9].real(), Catch::Matchers::WithinAbs(1.0, EPS));
        REQUIRE_THAT(croots3.value()[9].imag(), Catch::Matchers::WithinAbs(0.0, EPS));
        REQUIRE_THAT(croots3.value()[10].real(), Catch::Matchers::WithinAbs(1.0, EPS));
        REQUIRE_THAT(croots3.value()[10].imag(), Catch::Matchers::WithinAbs(0.0, EPS));
        REQUIRE_THAT(croots3.value()[11].real(), Catch::Matchers::WithinAbs(1.1142366961812986, EPS));
        REQUIRE_THAT(croots3.value()[11].imag(), Catch::Matchers::WithinAbs(-0.48083981203389980, EPS));
        REQUIRE_THAT(croots3.value()[12].real(), Catch::Matchers::WithinAbs(1.1142366961812986, EPS));
        REQUIRE_THAT(croots3.value()[12].imag(), Catch::Matchers::WithinAbs(0.48083981203389980, EPS));
        REQUIRE_THAT(croots3.value()[13].real(), Catch::Matchers::WithinAbs(2.0, EPS));
        REQUIRE_THAT(croots3.value()[13].imag(), Catch::Matchers::WithinAbs(0.0, EPS));
        REQUIRE_THAT(croots3.value()[14].real(), Catch::Matchers::WithinAbs(2.0, EPS));
        REQUIRE_THAT(croots3.value()[14].imag(), Catch::Matchers::WithinAbs(0.0, EPS));
    }
}

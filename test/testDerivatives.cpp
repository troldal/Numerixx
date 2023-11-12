//
// Created by Kenneth Balslev on 24/02/2023.
//

#include <Deriv.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <functional>
#include <numbers>
#include <vector>

TEST_CASE("nxx::deriv - Numerical Derivatives Test", "[derivatives]")
{
    using namespace nxx::deriv;

    std::vector< std::function< double(double) > > functions {
        [](double x) { return std::pow(x, 3) - 2 * x + 5; },
        [](double x) { return 2 * std::pow(x, 2) + 3 * x - 4; },
        [](double x) { return std::sin(x) + std::cos(x); },
        [](double x) { return std::log(x) + 2 * x; },
        [](double x) { return 4 * std::pow(x, 4) - 3 * std::pow(x, 3) + 2 * std::pow(x, 2) - x + 1; },
        [](double x) { return std::exp(x) + 3 * std::pow(x, 2); },
        [](double x) { return std::cos(x * x) - 2 * x; },
        [](double x) { return std::sqrt(x) + 2.0 / x; },
        [](double x) { return 3 * std::pow(x, 3) - 4 * std::pow(x, 2) + 5 * x - 6; },
        [](double x) { return 1.0 / (x + 1); },
        [](double x) { return std::exp(x); },
        [](double x) { return x * std::sqrt(x); },
        [](double x) { return std::sin(1.0 / x); },
        [](double x) { return std::exp(-x * x); },
        [](double x) { return x * x; },
        [](double x) { return 1.0 / x; }
    };

    std::vector< std::function< double(double) > > first_derivatives {
        [](double x) { return 3 * std::pow(x, 2) - 2; },
        [](double x) { return 4 * x + 3; },
        [](double x) { return std::cos(x) - std::sin(x); },
        [](double x) { return 1.0 / x + 2; },
        [](double x) { return 16 * std::pow(x, 3) - 9 * std::pow(x, 2) + 4 * x - 1; },
        [](double x) { return std::exp(x) + 6 * x; },
        [](double x) { return -2 * x * std::sin(x * x) - 2; },
        [](double x) { return 0.5 * std::pow(x, -0.5) - 2.0 / std::pow(x, 2); },
        [](double x) { return 9 * std::pow(x, 2) - 8 * x + 5; },
        [](double x) { return -1.0 / std::pow(x + 1, 2); },
        [](double x) { return std::exp(x); },
        [](double x) { return 1.5 * std::sqrt(x); },
        [](double x) { return -std::cos(1.0 / x) / (x * x); },
        [](double x) { return -2.0 * x * std::exp(-x * x); },
        [](double x) { return 2.0 * x; },
        [](double x) { return -1.0 / (x * x); }
    };

    std::vector< std::function< double(double) > > second_derivatives {
        [](double x) { return 6 * x; },
        [](double x) { return 4; },
        [](double x) { return -2 * std::sin(x); },
        [](double x) { return -1.0 / std::pow(x, 2); },
        [](double x) { return 48 * std::pow(x, 2) - 18 * x + 4; },
        [](double x) { return std::exp(x) + 6; },
        [](double x) { return -4 * x * x * std::cos(x * x) - 2 * std::sin(x * x); },
        [](double x) { return -0.25 * std::pow(x, -1.5) + 4.0 / std::pow(x, 3); },
        [](double x) { return 18 * x - 8; },
        [](double x) { return 2.0 / std::pow(x + 1, 3); },
        [](double x) { return std::exp(x); },
        [](double x) { return 0.75 / std::sqrt(x); },
        [](double x) { return 2.0 * std::cos(1 / x) / (x * x * x) - std::sin(1 / x) / (x * x * x * x); },
        [](double x) { return 4.0 * x * x * std::exp(-x * x) - 2.0 * std::exp(-x * x); },
        [](double x) { return 2.0; },
        [](double x) { return 2.0 / (x * x * x); }
    };

    std::vector< double > evals {
        2.0, 1.0, std::numbers::pi / 4, std::numbers::e, 0.0, 1.0, std::numbers::pi, 4.0, 2.0, 0.0, 1.0, 0.1, 0.45, 0.5, 0.0, 10.0
    };

    auto testDerivativeMethod = []< typename DerivativeMethod >(DerivativeMethod                                      method,
                                                                const std::vector< std::function< double(double) > >& functions,
                                                                const std::vector< std::function< double(double) > >& derivatives,
                                                                const std::vector< double >&                          evals,
                                                                double                                                tol) {
        for (size_t i = 0; i < functions.size(); ++i) {
            auto result = nxx::deriv::diff< DerivativeMethod >(functions[i], evals[i]);
            INFO("The function number is " << i);
            INFO("The result is " << *result);
            INFO("The expected result is " << derivatives[i](evals[i]));
            REQUIRE(result.has_value());
            REQUIRE_THAT(*result - derivatives[i](evals[i]), Catch::Matchers::WithinAbs(0.0, tol));
            REQUIRE_THROWS(nxx::deriv::diff< DerivativeMethod >(functions[i], evals[i], 0.0));
            REQUIRE(nxx::deriv::diff< DerivativeMethod >([](double x) { return std::sqrt(x); }, -1.0).has_value() == false);
            REQUIRE(nxx::deriv::diff< DerivativeMethod >([](double x) { return std::sqrt(x); }, 1.0).has_value() == true);
        }
    };

    // ============================================================================================
    // Computation of 1st order derivatives
    // ============================================================================================

    SECTION("nxx::deriv::central")
    { testDerivativeMethod(Order1CentralRichardson {}, functions, first_derivatives, evals, 1E-6); }

    SECTION("nxx::deriv::forward")
    { testDerivativeMethod(Order1ForwardRichardson {}, functions, first_derivatives, evals, 1E-6); }

    SECTION("nxx::deriv::backward")
    { testDerivativeMethod(Order1BackwardRichardson {}, functions, first_derivatives, evals, 1E-6); }

    SECTION("Order1CentralRichardson")
    { testDerivativeMethod(Order1CentralRichardson {}, functions, first_derivatives, evals, 1E-6); }

    SECTION("Order1Central3Point")
    { testDerivativeMethod(Order1Central3Point {}, functions, first_derivatives, evals, 1E-6); }

    SECTION("Order1Central5Point")
    { testDerivativeMethod(Order1Central5Point {}, functions, first_derivatives, evals, 1E-6); }

    SECTION("Order1ForwardRichardson")
    { testDerivativeMethod(Order1ForwardRichardson {}, functions, first_derivatives, evals, 1E-6); }

    SECTION("Order1Forward2Point")
    { testDerivativeMethod(Order1Forward2Point {}, functions, first_derivatives, evals, 1E-3); }

    SECTION("Order1Forward3Point")
    { testDerivativeMethod(Order1Forward3Point {}, functions, first_derivatives, evals, 1E-6); }

    SECTION("Order1BackwardRichardson")
    { testDerivativeMethod(Order1BackwardRichardson {}, functions, first_derivatives, evals, 1E-6); }

    SECTION("Order1Backward2Point")
    { testDerivativeMethod(Order1Backward2Point {}, functions, first_derivatives, evals, 1E-3); }

    SECTION("Order1Backward3Point")
    { testDerivativeMethod(Order1Backward3Point {}, functions, first_derivatives, evals, 1E-6); }

    // ============================================================================================
    // Computation of 2nd order derivatives
    // ============================================================================================

    SECTION("Order2Central3Point")
    { testDerivativeMethod(Order2Central3Point {}, functions, second_derivatives, evals, 1E-4); }

    SECTION("Order2Central5Point")
    { testDerivativeMethod(Order2Central5Point {}, functions, second_derivatives, evals, 1E-4); }

    SECTION("Order2Forward3Point")
    { testDerivativeMethod(Order2Forward3Point {}, functions, second_derivatives, evals, 1E-2); }

    SECTION("Order2Forward4Point")
    { testDerivativeMethod(Order2Forward4Point {}, functions, second_derivatives, evals, 1E-3); }

    SECTION("Order2Backward3Point")
    { testDerivativeMethod(Order2Backward3Point {}, functions, second_derivatives, evals, 1E-2); }

    SECTION("Order2Backward4Point")
    { testDerivativeMethod(Order2Backward4Point {}, functions, second_derivatives, evals, 1E-3); }
}
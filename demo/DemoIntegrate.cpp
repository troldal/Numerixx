////
//// Created by kenne on 12/11/2023.
////
//
//#include <Integration.hpp>
//#include <cmath>
//#include <fmt/format.h>
//#include <iomanip>
//#include <iostream>
//#include <numbers>
//
//int main()
//{
//    using namespace nxx::integrate;
//
//    double a = 0.0;
//    float  b = 10.0;
//
//    auto f      = [](double x) { return std::exp(x); };
////    auto result = integrate< Romberg >(f, a, b);
//
//    auto f_int = integralOf< Trapezoid >(f);
//
//     auto result = f_int(a, b);
//    auto check = std::exp(b) - std::exp(a);
//
//    std::cout << "Result: " << result << std::endl;
//    std::cout << "Check:  " << check << std::endl;
//
//    return 0;
//}

#include <Integ.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>

int main() {
    using namespace nxx::integrate;
    std::cout << std::fixed << std::setprecision(10);

    // Define some test functions
    auto f1 = [](std::floating_point auto x) { return x * x; }; // Integral from 0 to 1 is 1/3
    auto f2 = [](std::floating_point auto x) { return std::sin(x); }; // Integral from 0 to pi is 2
    auto f3 = [](std::floating_point auto x) { return std::exp(-x * x); }; // Integral from -inf to inf is sqrt(pi), but we'll approximate with a large interval

    // Compute the integrals using the different methods
    double integral_f1_romberg = integrate<Romberg>(f1, 0.0, 2.0);
    double integral_f1_simpson = integrate<Simpson>(f1, 0.0, 2.0);
    double integral_f1_trapezoid = integrate<Trapezoid>(f1, 0.0, 2.0);

    double integral_f2_romberg = integrate<Romberg>(f2, 0.0, std::numbers::pi);
    double integral_f2_simpson = integrate<Simpson>(f2, 0.0, std::numbers::pi);
    double integral_f2_trapezoid = integrate<Trapezoid>(f2, 0.0, std::numbers::pi);

    double integral_f3_romberg = integrate<Romberg>(f3, -10.0, 10.0);
    double integral_f3_simpson = integrate<Simpson>(f3, -10.0, 10.0);
    double integral_f3_trapezoid = integrate<Trapezoid>(f3, -10.0, 10.0);

    // Print the results
    std::cout << "Integral of x^2 from 0 to 1:\n";
    std::cout << "Romberg: " << integral_f1_romberg << "\n";
    std::cout << "Simpson: " << integral_f1_simpson << "\n";
    std::cout << "Trapezoid: " << integral_f1_trapezoid << "\n\n";

    std::cout << "Integral of sin(x) from 0 to pi:\n";
    std::cout << "Romberg: " << integral_f2_romberg << "\n";
    std::cout << "Simpson: " << integral_f2_simpson << "\n";
    std::cout << "Trapezoid: " << integral_f2_trapezoid << "\n\n";

    std::cout << "Integral of exp(-x^2) from -10 to 10 (approximation of sqrt(pi)):\n";
    std::cout << "Romberg: " << integral_f3_romberg << "\n";
    std::cout << "Simpson: " << integral_f3_simpson << "\n";
    std::cout << "Trapezoid: " << integral_f3_trapezoid << "\n\n";

    return 0;
}
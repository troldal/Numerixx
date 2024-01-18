////
//// Created by kenne on 12/11/2023.
////

#include <Integ.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <boost/multiprecision/cpp_bin_float.hpp>

using NXX_FLOAT = double; //boost::multiprecision::cpp_bin_float_50;

int main()
{
    using namespace nxx::integrate;
    std::cout << std::fixed << std::setprecision(50);

    // Define some test functions
    auto f1 = [](nxx::IsFloat auto x) { return x * x; };       // Integral from 0 to 1 is 1/3
    auto f2 = [](nxx::IsFloat auto x) { return sin(x); };      // Integral from 0 to pi is 2
    auto f3 = [](nxx::IsFloat auto x) { return exp(-x * x); }; // Integral from -inf to inf is sqrt(pi)

    struct bounds
    {
        NXX_FLOAT lower;
        NXX_FLOAT upper;
    };

    // Compute the integrals using the different methods
    auto f1_bounds = bounds{ 0.0, 2.0 };
    auto f2_bounds = bounds{ 0.0, std::numbers::pi };
    auto f3_bounds = bounds{ -10.0, 10.0 };

    std::cout << "Integral of x^2 from 0 to 2:\n";
    std::cout << "Romberg (bounds):        " << *integrate< Romberg >(f1, f1_bounds) << "\n";
    std::cout << "Romberg (init-list):     " << *integrate< Romberg >(f1, { 0.0, 2.0 }) << "\n";
    std::cout << "Simpson (bounds):        " << *integrate< Simpson >(f1, f1_bounds) << "\n";
    std::cout << "Simpson (init-list):     " << *integrate< Simpson >(f1, { 0.0, 2.0 }) << "\n";
    std::cout << "Trapezoid (bounds):      " << *integrate< Trapezoid >(f1, f1_bounds) << "\n";
    std::cout << "Trapezoid (init-list):   " << *integrate< Trapezoid >(f1, { 0.0, 2.0 }) << "\n\n";

    std::cout << "Integral of sin(x) from 0 to pi:\n";
    std::cout << "Romberg (bounds):        " << *integrate< Romberg >(f2, f2_bounds) << "\n";
    std::cout << "Romberg (init-list):     " << *integrate< Romberg >(f2, { 0.0, std::numbers::pi }) << "\n";
    std::cout << "Simpson (bounds):        " << *integrate< Simpson >(f2, f2_bounds) << "\n";
    std::cout << "Simpson (init-list):     " << *integrate< Simpson >(f2, { 0.0, std::numbers::pi }) << "\n";
    std::cout << "Trapezoid (bounds):      " << *integrate< Trapezoid >(f2, f2_bounds) << "\n";
    std::cout << "Trapezoid (init-list):   " << *integrate< Trapezoid >(f2, { 0.0, std::numbers::pi }) << "\n\n";

    std::cout << "Integral of exp(-x^2) from -10 to 10 (approximation of sqrt(pi)):\n";
    std::cout << "Romberg (bounds):        " << *integrate< Romberg >(f3, f3_bounds) << "\n";
    std::cout << "Romberg (init-list):     " << *integrate< Romberg >(f3, { -10.0, 10.0 }) << "\n";
    std::cout << "Simpson (bounds):        " << *integrate< Simpson >(f3, f3_bounds) << "\n";
    std::cout << "Simpson (init-list):     " << *integrate< Simpson >(f3, { -10.0, 10.0 }) << "\n";
    std::cout << "Trapezoid (bounds):      " << *integrate< Trapezoid >(f3, f3_bounds) << "\n";
    std::cout << "Trapezoid (init-list):   " << *integrate< Trapezoid >(f3, { -10.0, 10.0 }) << "\n\n";

    const auto integral_f1_romberg   = integralOf(f1);
    const auto integral_f1_simpson   = integralOf< Simpson >(f1);
    const auto integral_f1_trapezoid = integralOf< Trapezoid >(f1);

    std::cout << "Integration of x^2 from 0 to 2 using the integralOf function:\n";
    std::cout << "Romberg:     " << integral_f1_romberg({ 0.0, 2.0 }) << "\n";
    std::cout << "Simpson:     " << integral_f1_simpson({ 0.0, 2.0 }) << "\n";
    std::cout << "Trapezoid:   " << integral_f1_trapezoid({ 0.0, 2.0 }) << "\n\n";

    // TODO: This doesn't work on clang-cl
    // try {
    //     auto tmp = integralOf< Trapezoid >(f1);
    //     std::cout << tmp(f1_bounds, 1E-12, 5) << "\n";
    // }
    // catch (const nxx::NumerixxError& e) { std::cout << e.log() << "\n"; }

    auto manualIntegrate = [](auto solver, double tolerance = 1e-12, int maxIterations = 25) {
        auto result = solver.current();
        std::cout << "Manual:      " << result << "\n";
        for (auto i = 0; i < maxIterations; i++) {
            solver.iterate();
            if (abs(result - solver.current()) < tolerance)
                return solver.current();
            result = solver.current();
            std::cout << "Manual:      " << result << "\n";
        }

        return result;
    };

    auto rombergResult_1 = manualIntegrate(Trapezoid(f3, f3_bounds));
    auto rombergResult_2 = manualIntegrate(Trapezoid(f3, { -10.0, 10.0 }));


    return 0;
}

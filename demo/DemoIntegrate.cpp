//
// Created by kenne on 12/11/2023.
//

#include <Integration.hpp>
#include <cmath>
#include <fmt/format.h>
#include <iomanip>
#include <iostream>
#include <numbers>

int main()
{
    using namespace nxx::integrate;

    double a = 0.0;
    float  b = 2.0;

    auto f      = [](double x) { return std::exp(x); };
    auto result = integrate< Romberg >(f, a, b);

    auto f_int = integralOf< Trapezoid >(f);

    // auto result = f_int(a, b);
    auto check = std::exp(b) - std::exp(a);

    std::cout << "Result: " << result << std::endl;
    std::cout << "Check:  " << check << std::endl;

    return 0;
}
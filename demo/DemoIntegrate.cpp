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

    auto f      = [](double x) { return x * x; };
    auto result = integrate(RombergMethod {}, f, 0.0, 100.0);
    auto check  = (1.0 / 3.0) * std::pow(100.0, 3.0);

    std::cout << "Result: " << result << std::endl;
    std::cout << "Check:  " << check << std::endl;

    return 0;
}
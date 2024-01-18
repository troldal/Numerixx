#include <Optim.hpp>
#include <cmath>

#include <iomanip>
#include <iostream>
#include <numbers>

int main()
{
    using namespace nxx::optimize;

    auto myFunc = [](double x) { return x * x * x - 4 * x * x + x - 5; };    // Example function (quadratic)

    // Direct optimization
    GradientDescentOptimizer optimizer(0.01, GradientDescentOptimizer::Mode::Minimize);
    double                   result1 = optimize(optimizer, myFunc, 4);    // Starting guess is -5
    std::cout << "Optimized value: " << result1 << std::endl;

    // Using lambda function for optimization
    auto   optimizeFunc = optimizationOf(optimizer, myFunc);
    double result2      = optimizeFunc(4);    // Starting guess is 5
    std::cout << "Optimized value: " << result2 << std::endl;

    return 0;
}
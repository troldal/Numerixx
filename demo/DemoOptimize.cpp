#include <Deriv.hpp>
#include <Optim.hpp>
#include <cmath>

#include <iomanip>
#include <iostream>
#include <numbers>

int main()
{
    using namespace nxx::optim;
    using namespace nxx::deriv;

    auto myFunc = [](double x) { return x * x * x - 4 * x * x + x - 5; };    // Example function (quadratic)
        // auto myFunc = [](double x) { return -x*x; };    // Example function (quadratic)

    auto result  = optimize<GradientDescent, Minimize>(myFunc, derivativeOf(myFunc), 4.0, 1e-12, 10000);
    std::cout << "Optimized value: " << result << std::endl;

    return 0;
}
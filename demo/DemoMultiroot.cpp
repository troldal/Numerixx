//
// Created by Kenneth Balslev on 09/12/2022.
//

#include <deque>
#include <iostream>
#include <list>
#include <Multiroots.hpp>
#include <vector>
#include <cmath>
#include <numbers>

int main()
{
    using namespace nxx::deriv;
    using namespace nxx::multiroots;

    DynamicFunctionArray<double> functions = {
        [](const std::vector< double >& coeffs) { return 3 * coeffs[0] - std::cos(coeffs[1]*coeffs[2]) - 0.5; },
        [](const std::vector< double >& coeffs) { return coeffs[0] * coeffs[0] - 81*std::pow(coeffs[1] + 0.1,2)+std::sin(coeffs[2])+1.06; },
        [](const std::vector< double >& coeffs) { return std::exp(-coeffs[0]*coeffs[1]) + 20*coeffs[2] + (10*std::numbers::pi-3)/3; }
    };

    auto J2 = jacobian(functions, { 0.1, 0.1, -0.1 });
    std::cout << J2 << std::endl;


    std::cout << std::fixed;
    auto solver = DMultiNewton(functions, { 2.0, 2.0, 2.0 });
    auto result = multisolve(solver, { 2.0, 2.0, 2.0 }, 1.0e-10, 100);
    for (auto g : result) std::cout << g << std::endl;

    auto root = functions(result);
    for (auto r : root) std::cout << r << std::endl;

    return 0;
}

//
// Created by Kenneth Balslev on 09/12/2022.
//

#include <Multiroots.hpp>
#include <blaze/Blaze.h>
#include <cmath>
#include <deque>
#include <iomanip>
#include <iostream>
#include <list>
#include <numbers>
#include <span>
#include <vector>

int main()
{
    std::cout << std::fixed << std::setprecision(8);
    using namespace nxx::deriv;
    using namespace nxx::multiroots;

    auto f1 = [](std::span< double > coeffs) { return 3 * coeffs[0] - std::cos(coeffs[1] * coeffs[2]) - 0.5; };
    auto f2 = [](std::span< double > coeffs) {
        return coeffs[0] * coeffs[0] - 81 * std::pow(coeffs[1] + 0.1, 2) + std::sin(coeffs[2]) + 1.06;
    };
    auto f3 = [](std::span< double > coeffs) {
        return std::exp(-coeffs[0] * coeffs[1]) + 20 * coeffs[2] + (10 * std::numbers::pi - 3) / 3;
    };

    MultiFunctionArray functions { f1, f2, f3 };
    std::cout << functions.eval< blaze::DynamicVector >({ 0., 0., 0. }) << std::endl;


    auto J = jacobian(functions, blaze::DynamicVector{ 0., 0., 0. });
    std::cout << J << std::endl;

    // std::cout << std::fixed << std::setprecision(20);
    // auto solver = DMultiNewton(functions, { 0.0, 0.0, 0.0 });
    auto solver = SteepestDescent(functions, { 0.0, 0.0, 0.0 });
    auto result = multisolve(solver, { 2.0, 2.0, 2.0 }, 1.0e-2, 1000);
    std::cout << result << std::endl;

    std::cout << functions.eval< blaze::DynamicVector>(result) << std::endl;


    return 0;
}

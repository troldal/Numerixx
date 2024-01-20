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

    auto result1 = multisolve<SteepestDescent>(functions, { 2.0, 2.0, 2.0 });
    // auto result2 = multisolve<MultiNewton>(functions, *result1);

    std::cout << "Root:\n" << *result1 << std::endl;
    std::cout << "Result:\n" << functions(*result1) << std::endl;

    // auto               f1 = [](std::span< double > coeffs) { return 1 - coeffs[0]; };
    // auto               f2 = [](std::span< double > coeffs) { return 10 * (coeffs[1] - coeffs[0] * coeffs[0]); };
    // MultiFunctionArray functions { f1, f2 };
    //
    // auto result = multisolve< MultiNewton >(functions, { -10.0, -5.0 });
    // std::cout << "Root:\n" << *result << std::endl;
    // std::cout << "Result:\n" << functions(*result) << std::endl;

    // auto               f1 = [](std::span< double > coeffs) { return pow(coeffs[0], 2) + coeffs[0] * coeffs[1] - 10; };
    // auto               f2 = [](std::span< double > coeffs) { return coeffs[1] + 3 * coeffs[0] * pow(coeffs[1], 2) - 57; };
    // MultiFunctionArray functions { f1, f2 };
    //
    // auto result = multisolve< MultiNewton >(functions, { 10., 10. });
    // std::cout << "Root:\n" << *result << std::endl;
    // std::cout << "Result:\n" << functions(*result) << std::endl;

    return 0;
}

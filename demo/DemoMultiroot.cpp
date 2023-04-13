//
// Created by Kenneth Balslev on 09/12/2022.
//

#include <deque>
#include <iostream>
#include <list>
#include <numerixx.hpp>
#include <vector>
#include <cmath>

int main()
{
    using namespace nxx::deriv;
    using namespace nxx::multiroots;

//    std::vector< std::function< double(std::vector< double >) > > fns = {
//        [](const std::vector< double >& coeffs) { return coeffs[1] * coeffs[1] * (1.0 - coeffs[0]) - coeffs[0] * coeffs[0] * coeffs[0]; },
//        [](const std::vector< double >& coeffs) { return coeffs[0] * coeffs[0] + coeffs[1] * coeffs[1] - 1.0; }
//    };

//    std::vector< std::function< double(std::vector< double >) > > fns = {
//        [](const std::vector< double >& coeffs) { return 1 - coeffs[0]; },
//        [](const std::vector< double >& coeffs) { return 10 * (coeffs[1] - coeffs[0] * coeffs[0]);}
//    };

    std::vector< std::function< double(std::vector< double >) > > fns = {
        [](const std::vector< double >& coeffs) { return 3 * coeffs[0] - std::cos(coeffs[1]*coeffs[2]) - 0.5; },
        [](const std::vector< double >& coeffs) { return coeffs[0] * coeffs[0] - 81*std::pow(coeffs[1] + 0.1,2)+std::sin(coeffs[2])+1.06; },
        [](const std::vector< double >& coeffs) { return std::exp(-coeffs[0]*coeffs[1]) + 20*coeffs[2] + (10*M_PI-3)/3; }
    };

    std::cout << std::fixed;
    auto solver = nxx::multiroots::DMultiNewton(fns);
    auto result = multisolve(solver, { 2.0, 2.0, 2.0 }, 1.0e-6, 100);
    for (auto g : result) std::cout << g << std::endl;

    return 0;
}
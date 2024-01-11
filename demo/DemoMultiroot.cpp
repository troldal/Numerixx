//
// Created by Kenneth Balslev on 09/12/2022.
//

#include <deque>
#include <iostream>
#include <iomanip>
#include <list>
#include <Multiroots.hpp>
#include <vector>
#include <cmath>
#include <numbers>
#include <span>
#include <blaze/Blaze.h>

int main()
{
    std::cout << std::fixed << std::setprecision(8);
    using namespace nxx::deriv;
    using namespace nxx::multiroots;

    auto f1 = []( std::span< double > coeffs) { return 3 * coeffs[0] - std::cos(coeffs[1]*coeffs[2]) - 0.5; };
    auto f2 = []( std::span< double > coeffs) { return coeffs[0] * coeffs[0] - 81*std::pow(coeffs[1] + 0.1,2)+std::sin(coeffs[2])+1.06; };
    auto f3 = []( std::span< double > coeffs) { return std::exp(-coeffs[0]*coeffs[1]) + 20*coeffs[2] + (10*std::numbers::pi-3)/3; };

    auto fx = MultiFunction(f1);
    std::cout << fx({ 0.1, 0.1, -0.1 }) << std::endl;

    std::vector args{ 0.1, 0.1, -0.1 };
    std::cout << fx(args) << std::endl;

    blaze::DynamicVector<double> x{ 0.1, 0.1, -0.1 };
    std::cout << fx(x) << std::endl;

    auto dfx1 = partialdiff(fx, args);
    for (auto d : dfx1) std::cout << d << std::endl;

    auto dfx2 = partialdiff(fx, x);
    std::cout << dfx2 << std::endl;


    auto dfx3 = partialdiff<Order1Central3Point>(fx, { 0.1, 0.1, -0.1 });
    for (auto d : dfx3) std::cout << d << std::endl;

    auto dfx4 = partialdiff<std::deque, Order1Central3Point>(fx, { 0.1, 0.1, -0.1 });
    for (auto d : dfx4) std::cout << d << std::endl;
//    std::cout << dfx4 << std::endl;

    auto dfx5 = partialdiff<blaze::DynamicVector, Order1Central3Point>(fx, { 0.1, 0.1, -0.1 });
    for (auto d : dfx5) std::cout << d << std::endl;


    //    DynamicFunctionArray<double> functions { std::vector< std::function< double(const std::span< double >)> >{ f1, f2, f3 }};

//    DynamicFunctionArray<double> functions { std::vector< std::function< double(const std::vector< double >&)> >{
//        [](const std::vector< double >& coeffs) { return 3 * coeffs[0] - std::cos(coeffs[1]*coeffs[2]) - 0.5; },
//        [](const std::vector< double >& coeffs) { return coeffs[0] * coeffs[0] - 81*std::pow(coeffs[1] + 0.1,2)+std::sin(coeffs[2])+1.06; },
//        [](const std::vector< double >& coeffs) { return std::exp(-coeffs[0]*coeffs[1]) + 20*coeffs[2] + (10*std::numbers::pi-3)/3; }
//    }};

//    auto J2 = jacobian(functions, { 0.1, 0.1, -0.1 });
//    std::cout << J2 << std::endl;


//    std::cout << std::fixed << std::setprecision(20);
//    auto solver = DMultiNewton(functions, { 2.0, 2.0, 2.0 });
//    auto result = multisolve(solver, { 2.0, 2.0, 2.0 }, 1.0e-10, 100);
//    for (auto g : result) std::cout << g << std::endl;
//
//    auto root = functions(result);
//    for (auto r : root) std::cout << r << std::endl;

    return 0;
}

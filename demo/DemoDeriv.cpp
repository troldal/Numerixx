//
// Created by Troldal on 14/10/2022.
//

#include <calculus/derivatives.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <poly/Polynomial.hpp>
#include <polyroots/polyroots.hpp>

int main() {

    using namespace numerix;

    //numerix::poly::polynomial p({1.2, -0.25, -0.5, -0.15, -0.1});
    poly::Polynomial p({0.0, 0.0, 0.0, 0.0, 1.0});
    auto diff = deriv::Derivative(p);
    std::cout << std::fixed << std::setprecision(20);
    std::cout << diff(1.0) << std::endl;
    std::cout << p.derivative(1.0) << std::endl;

//    auto diff = calculus::Derivative([](long double x) { return std::sqrt(x); });
//    std::cout << std::fixed << std::setprecision(20);
//    std::cout << diff(0.0) << std::endl;

    return 0;
}
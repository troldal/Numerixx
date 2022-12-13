//
// Created by Troldal on 14/10/2022.
//

#include <deque>
#include <iomanip>
#include <iostream>
#include <poly/Polynomial.hpp>
#include <polyroots/Polyroots.hpp>

int main() {

    using namespace numerix::poly;
    using namespace numerix::polyroots;

    std::cout << "Solving f(x) = -5 + x^2\n";
    auto fun1 = Polynomial({-5.0, 0.0, 1.0});

    std::cout << "Coefficients: ";
    for (auto c : fun1.coefficients()) std::cout << c << " ";
    std::cout << std::endl;

    std::cout << "Roots: ";
    for (auto r : quadratic(fun1)) std::cout << r << " ";
    std::cout << std::endl;

    std::cout << "Derivatives at Roots: ";
    for (auto r : quadratic(fun1)) std::cout << fun1.derivative(r) << " ";
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "Solving f(x) = -2 - 1.5x + 0.75x^2 + 0.25x^3\n";
    auto fun2 = Polynomial({-2.0, -1.5, 0.75, 0.25});

    std::cout << "Coefficients: ";
    for (auto c : fun2.coefficients()) std::cout << c << " ";
    std::cout << std::endl;

    std::cout << "Roots: ";
    for (auto r : cubic(fun2)) std::cout << r << " ";
    std::cout << std::endl;

    std::cout << "Derivatives at Roots: ";
    for (auto r : cubic(fun2)) std::cout << fun2.derivative(r) << " ";
    std::cout << std::endl;
    std::cout << std::endl;

//    std::vector<double> c {-5, 0, 1};
//    auto fun = numerix::poly::polynomial(c);

//    std::cout << fun(2.23606797749978980505) << std::endl;
//    std::cout << fun(0) << std::endl;
//    std::cout << fun.derivative(-3.0) << std::endl;
//
//    //std::cout << std::fixed << std::setprecision(20);
//    std::cout << numerix::polyroots::quadratic(fun).front() << "\t" << numerix::polyroots::quadratic(fun).back() << std::endl;
//
//    auto fun2 = numerix::poly::Polynomial({-2.0, -1.5, 0.75, 0.25});
//    auto res = numerix::polyroots::cubic(fun2);
//
//    for (auto r : res)
//        std::cout << r << std::endl;
//
//    auto cof = fun.coefficients();
//    for (auto i : cof) std::cout << i << " ";

    return 0;
}
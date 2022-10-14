//
// Created by Troldal on 14/10/2022.
//

#include <iostream>
#include <iomanip>
#include <poly/polynomial.hpp>
#include <polyroots/polyroots.hpp>

int main() {

    auto fun = numerix::poly::polynomial({-5.0, 0.0, 1.0});

    std::cout << fun(2.23606797749978980505) << std::endl;
    std::cout << fun(0) << std::endl;
    std::cout << fun.derivative(-3.0) << std::endl;

    //std::cout << std::fixed << std::setprecision(20);
    std::cout << numerix::polyroots::quadratic(fun).front() << "\t" << numerix::polyroots::quadratic(fun).back() << std::endl;

    auto fun2 = numerix::poly::polynomial({-2.0, -1.5, 0.75, 0.25});
    auto res = numerix::polyroots::cubic(fun2);

    for (auto r : res)
        std::cout << r << std::endl;

    return 0;
}
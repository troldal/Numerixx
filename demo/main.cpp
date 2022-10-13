//
// Created by Troldal on 06/10/2022.
//


#include <iostream>
#include <iomanip>
#include <fmt/format.h>
#include <numerix.hpp>

auto fun(double x) {
    return x * x - 5;
}

int main() {

    using namespace numerix::roots;

    //auto solver = bisection([](double x) {return x * x - 5;});
    auto solver = dnewton(fun);
    auto guess = 10.0;
    solver.init(guess);

    std::cout << fmt::format("{:>10} | {:>25} | {:>25} ", "Iter", "Root", "Error") << std::endl;
    std::cout << "------------------------------------------------------------------\n";

    for (int i = 0; i <= 100; ++i) {
        std::cout << std::fixed << std::setprecision(10);

        std::cout << fmt::format("{:10} | {:25.20f} | {:25.20f} ", i, guess, abs(solver.evaluate(guess))) << std::endl;


        if (abs(solver.evaluate(guess)) < 1.0E-15) break;

        solver.iterate();
        guess = solver.result();
    }

    std::cout << "------------------------------------------------------------------\n";

    auto result1 = fsolve(ridders(fun), std::make_pair(0.0, 2.5), 1.0E-15);
    std::cout << std::setprecision(20) << result1 << std::endl;

    auto result2 = fsolve(bisection(fun), std::make_pair(0.0, 2.5), 1.0E-15);
    std::cout << std::setprecision(20) << result2 << std::endl;

    auto result3 = fdfsolve(dnewton(fun), 3.0, 1.0E-15);
    std::cout << std::setprecision(20) << result3 << std::endl;

    auto result4 = fdfsolve(newton(fun, [](double x){return 2 * x;}), 3.0, 1.0E-15);
    std::cout << std::setprecision(20) << result4 << std::endl;

    return 0;
}
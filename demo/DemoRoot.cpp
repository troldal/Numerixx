//
// Created by Troldal on 06/10/2022.
//


#include <iostream>
#include <iomanip>
#include <fmt/format.h>
#include <numerix.hpp>


auto fun = numerix::poly::Polynomial({-5.0, 0.0, 1.0});

template<typename S>
void print(S solver, std::pair<double, double> b) {

    auto bounds = b;
    solver.init(bounds);

    std::cout << "----------------------------------------------------------------------------------\n";
    std::cout << fmt::format("{:>10} | {:>15} | {:>15} | {:>15} | {:>15} ", "Iter", "Upper", "Lower", "Root", "Error") << std::endl;
    std::cout << "----------------------------------------------------------------------------------\n";

    for (int i = 0; i <= 100; ++i) {

        std::cout << fmt::format("{:10} | {:15.10f} | {:15.10f} | {:15.10f} | {:15.10f} "
                                 , i
                                 , bounds.first
                                 , bounds.second
                                 , (bounds.first + bounds.second) / 2.0
                                 , abs(solver.evaluate((bounds.first + bounds.second) / 2.0))) << std::endl;

        if (abs(solver.evaluate((bounds.first + bounds.second) / 2.0)) < 1.0E-15) break;

        solver.iterate();
        bounds = solver.bounds();

    }

    std::cout << std::fixed << std::setprecision(20 )<< "CONVERGED! Root found at: " << (bounds.first + bounds.second) / 2.0 << "\n";
    std::cout << "----------------------------------------------------------------------------------\n\n";

}

template<typename S>
void print(S solver, double g) {

    auto guess = g;
    solver.init(guess);

    std::cout << "------------------------------------------------------------------\n";
    std::cout << fmt::format("{:>10} | {:>25} | {:>25} ", "Iter", "Root", "Error") << std::endl;
    std::cout << "------------------------------------------------------------------\n";

    for (int i = 0; i <= 100; ++i) {
        std::cout << std::fixed << std::setprecision(10);

        std::cout << fmt::format("{:10} | {:25.20f} | {:25.20f} ", i, guess, abs(solver.evaluate(guess))) << std::endl;

        if (abs(solver.evaluate(guess)) < 1.0E-15) break;

        solver.iterate();
        guess = solver.result();
    }

    std::cout << std::fixed << std::setprecision(20 )<< "CONVERGED! Root found at: " << guess << "\n";
    std::cout << "------------------------------------------------------------------\n\n";
}

int main() {

    using numerix::roots::fsolve;
    using numerix::roots::fdfsolve;
    using numerix::roots::Ridders;
    using numerix::roots::Bisection;
    using numerix::roots::Newton;
    using numerix::roots::DNewton;

    std::cout << "RIDDERS:" << std::endl;
    print(Ridders(fun), {0.0, 2.5});

    std::cout << "BISECTION:" << std::endl;
    print(Bisection(fun), {0.0, 2.5});

    std::cout << "DISCRETE NEWTON:" << std::endl;
    print(DNewton(fun), 3.0);

    std::cout << "NEWTON:" << std::endl;
    print(Newton(fun, [&](double x){return fun.derivative(x);}), 1.25);

    std::cout << "Ridders:         " << fsolve(Ridders(fun), {0.0, 2.5}, 1.0E-15) << std::endl;
    std::cout << "Bisection:       " << fsolve(Bisection(fun), {0.0, 2.5}, 1.0E-15) << std::endl;
    std::cout << "Discrete Newton: " << fdfsolve(DNewton(fun), 1.25, 1.0E-15) << std::endl;
    std::cout << "Newton:          " << fdfsolve(Newton(fun, [&](double x){return fun.derivative(x);}), 1.25, 1.0E-15) << std::endl;

    std::cout << "Bisection2:      " << fsolve(Bisection([](double x){return x*x - 5;}), {0.0, 2.5}) << std::endl;
    std::cout << "Discrete Newton2:" << fdfsolve(DNewton([](double x){return x*x - 5;}), 1.25) << std::endl;
    std::cout << "Newton2:         " << fdfsolve(Newton([](double x){return x*x - 5;}, [](double x){return 2*x;}), 1.25) << std::endl;

    return 0;
}
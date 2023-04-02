//
// Created by Troldal on 06/10/2022.
//

#include <fmt/format.h>
#include <iomanip>
#include <iostream>
#include <numerixx.hpp>

auto fun = nxx::poly::Polynomial({-5.0, 0.0, 1.0});

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
    using nxx::roots::Ridders;
    using nxx::roots::Bisection;
    using nxx::roots::DNewton;
    using nxx::roots::Newton;
    using namespace nxx::deriv;
    using namespace nxx::poly;

    std::cout << "RIDDERS:" << std::endl;
    print(Ridders(fun), {0.0, 2.5});

    std::cout << "BISECTION:" << std::endl;
    print(Bisection(fun), {0.0, 2.5});

    std::cout << "DISCRETE NEWTON:" << std::endl;
    print(DNewton(fun), 3.0);

    std::cout << "NEWTON:" << std::endl;
    print(Newton(fun, derivativeOf(fun)), 1.25);

    std::cout << "Ridders:         " << *fsolve(Ridders(fun), {0.0, 2.5}, 1.0E-15) << std::endl;
    std::cout << "Bisection:       " << *fsolve(Bisection(fun), {0.0, 2.5}, 1.0E-15) << std::endl;
    std::cout << "Discrete Newton: " << *fdfsolve(DNewton(fun), 1.25, 1.0E-15) << std::endl;
    std::cout << "Newton:          " << *fdfsolve(Newton(fun, derivativeOf(fun)), 1.25, 1.0E-15) << std::endl;

    std::vector<std::function<double(double)>> functions {
        [](double x){return std::sin(x) - x/2.0;},
        [](double x){return std::exp(x) - 3*x;},
        [](double x){return std::tan(x) - x;},
        [](double x){return std::log(x) + x;},
        [](double x){return std::cos(x) - std::pow(x,3);},
        [](double x){return std::sqrt(x) - std::cos(x);},
        [](double x){return std::pow(x, 1.0/3) + std::pow(x, 1.0/5) - 1;},
    };

    std::vector<std::function<double(double)>> derivatives {
        [](double x){return std::cos(x) - 0.5;},
        [](double x){return std::exp(x) - 3;},
        [](double x){return std::pow(1.0/ std::cos(x), 2) - 1;},
        [](double x){return 1.0/x + 1;},
        [](double x){return -std::sin(x) - 3 * std::pow(x,2);},
        [](double x){return 1.0/(2*std::sqrt(x)) + std::sin(x);},
        [](double x){return 1.0/(3*std::pow(x, 2.0/3)) + 1.0/(5*std::pow(x, 4.0/5));},
    };

    std::vector<std::pair<double, double>> brackets {
        {1.0, 3.0},
        {0.0, 1.0},
        {4.0, 4.5},
        {0.5, 1.0},
        {0.5, 1.5},
        {0.0, 1.0},
        {0.0, 0.2}
    };

    for (size_t i = 0; i <= 6; ++i)
        std::cout << *fdfsolve(Newton(functions[i], derivatives[i]), (brackets[i].second + brackets[i].second)/2.0, 1.0E-15) << std::endl;

    return 0;
}
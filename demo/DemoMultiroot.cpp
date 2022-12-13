//
// Created by Kenneth Balslev on 09/12/2022.
//

#include <numerix.hpp>
#include <vector>
#include <deque>
#include <list>
#include <iostream>


int main() {

    auto func1 = [](std::vector< double> coeffs) {
        return coeffs[1] * coeffs[1] * (1.0 - coeffs[0]) - coeffs[0] * coeffs[0] * coeffs[0];
    };

    auto func2 = [](std::vector< double> coeffs) {
        return coeffs[0] * coeffs[0] + coeffs[1] * coeffs[1] - 1.0;
    };

    using function_type = std::function< double(std::vector< double>)>;
    std::list<function_type > fns;
    fns.emplace_back(func1);
    fns.emplace_back(func2);

    auto solver = numerix::multiroots::DMultiNewton(fns);
    solver.init({1.0, 1.0});

    for (int i = 0; i < 10; ++i)
        solver.iterate();

    auto res = solver.result();
    for(auto r : res) std::cout << r << std::endl;

}
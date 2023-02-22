//
// Created by Kenneth Balslev on 09/12/2022.
//

#include <deque>
#include <iostream>
#include <list>
#include <numerixx.hpp>
#include <vector>

int main() {

    using namespace nxx::multiroots;
    using nxx::linalg::Matrix;

    using function_type = std::function< double(Matrix< double>)>;

    function_type func1 = [](const Matrix< double>& coeffs) {
        return coeffs[1] * coeffs[1] * (1.0 - coeffs[0]) - coeffs[0] * coeffs[0] * coeffs[0];
    };

    function_type func2 = [](const Matrix< double>& coeffs) {
        return coeffs[0] * coeffs[0] + coeffs[1] * coeffs[1] - 1.0;
    };


//    auto solver = numerix::multiroots::DMultiNewton(fns);
//    solver.init({2.0, 2.0});
//
//    for (int i = 0; i < 10; ++i)
//        solver.iterate();
//
//    auto res = solver.result();
//    for(auto r : res) std::cout << r << std::endl;

    MultiFunction f({ func1, func2 });
    auto tmp = f({2.0, 2.0});
    std::cout << tmp << std::endl;

//    auto solver = numerix::multiroots::DMultiNewton(f);
//    solver.init({2.0, 2.0});
//
//    for (int i = 0; i < 10; ++i) solver.iterate();
//
//    auto res = solver.result();
//    for (auto r : res) std::cout << r << std::endl;
}
//
// Created by kenne on 30/08/2023.
//

#include <Poly.hpp>
#include <deque>
#include <iomanip>
#include <iostream>
#include <list>
#include <set>
#include <unordered_set>

int main() {
    using namespace nxx::poly;
    using namespace std::complex_literals;

    auto poly = Polynomial{ 32.0, -48.0, -8.0, 28.0, -8.0, 16.0, -16.0, 12.0, -16.0, 6.0, 10.0, -17.0, 10.0, 2.0, -4.0, 1.0};
    std::cout << poly.asString() << std::endl;

        auto rroots = polysolve(poly);

        std::cout << std::setprecision(8);
        for (const auto& root : rroots) {
            std::cout << root << std::endl;
        }

        auto croots = polysolve<std::complex<double>>(poly);
        for (const auto& root : croots) {
            std::cout << root << std::endl;
        }

    return 0;
}
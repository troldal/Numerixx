#include <Interp.hpp>

#include <cmath>

#include <iomanip>
#include <iostream>
#include <numbers>

int main()
{
    using namespace nxx::interp;

    auto curve = Steffen({ { 0, 0 }, { 1, 2 }, { 2, 3 }, { 3, 2 }, { 4, 0 } });
        for (double x = 0.0; x <= 4.0; x += 0.1) {
                std::cout << std::setw(10) << x << "\t" << std::setw(10) << curve(x) << std::endl;
        }

    return 0;
}
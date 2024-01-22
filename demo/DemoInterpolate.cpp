#include <Interp.hpp>

#include <cmath>

#include <iomanip>
#include <iostream>
#include <numbers>

int main()
{
    using namespace nxx::interp;

    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    //    auto curve = interpolationOf<Linear>({ { 0, 0 }, { 1, 2 }, { 2, 3 }, { 3, 2 }, { 4, 0 } });
    //    auto curve = interpolationOf<Lagrange>({ { 0, 0 }, { 1, 2 }, { 2, 3 }, { 3, 2 }, { 4, 0 } });
    //    auto curve = interpolationOf<Steffen>({ { 0, 0 }, { 1, 2 }, { 2, 3 }, { 3, 2 }, { 4, 0 } });
    //    auto curve = interpolationOf<Spline>({ { 0, 0 }, { 1, 2 }, { 2, 3 }, { 3, 2 }, { 4, 0 } });

    auto curve = interpolationOf< Steffen >({ { 0, 1 },
                                             { 1.420735492, 1.540302306 },
                                             { 2.454648713, 1.346356379 },
                                             { 3.070560004, 2.088869738 },
                                             { 3.621598752, 3.04234052 },
                                             { 4.520537863, 5.991202812 },
                                             { 5.860292251, 5.87203631 },
                                             { 7.328493299, 7.300592544 },
                                             { 8.494679123, 8.39185723 },
                                             { 9.206059243, 9.776685982 } });

    for (int start = 0; start <= 100; start++) std::cout << (start / 10.0) << "\t" << curve.extrapolate(start / 10.0) << std::endl;

    return 0;
}
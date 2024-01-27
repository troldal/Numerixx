#include "_external.hpp"

#include <Deriv.hpp>
#include <Optim.hpp>
#include <cmath>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <ranges>

int main()
{
    using namespace nxx::optim;
    using namespace nxx::deriv;
    using namespace sciplot;
    namespace rng = std::ranges;
    
    std::cout << std::setprecision(8) << std::fixed;

    auto myFunc = [](double x) { return -x * x * x + 4 * x * x - x - 5; };    // Example function (quadratic)
           //auto myFunc = [](double x) { return x*x; };    // Example function (quadratic)
     //auto myFunc = [](double x) { return pow(x,6) - 11*pow(x,3)+17*x*x - 7*x+1; };    // Example function (quadratic)

    auto result = optimize< GradientDescent, Minimize >(myFunc, derivativeOf(myFunc), 0.0, 1e-12, 10000);
    std::cout << "Optimized value: " << result << std::endl;
    std::cout << "Function value: " << myFunc(result) << std::endl;

    auto optimizer = Brent< decltype(myFunc), double, Minimize >(myFunc, { -100.0, 1.0 });

    for (int i = 0; i <= 100; ++i) {
        auto current = optimizer.current();
        std::cout << "Iteration " << i << ": " << current.first << "\t" << current.second << std::endl;
        // std::cout << "Iteration " << i << ": ";
        optimizer.iterate();
    }

    std::cout << myFunc(optimizer.current().first) << std::endl;
    std::cout << myFunc(optimizer.current().second) << std::endl;

    Vec x = linspace(-1.0, 1.0, 200);

    std::vector< double > y_lin;
    std::transform(begin(x), end(x), std::back_inserter(y_lin), [&](auto val) { return myFunc(val); });

    Plot2D plot1;
    plot1.palette("paired");
    plot1.drawCurve(x, y_lin).label("Linear").lineWidth(4);

    Figure figure = { { plot1 } };
    Canvas canvas = { { figure } };

    canvas.defaultPalette("set1");
    canvas.size(1600, 1200);
    canvas.show();

    return 0;
}

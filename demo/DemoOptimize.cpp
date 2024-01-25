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

    auto myFunc = [](double x) { return x * x * x - 4 * x * x + x - 5; };    // Example function (quadratic)
          // auto myFunc = [](double x) { return x*x; };    // Example function (quadratic)

    auto result  = optimize<GradientDescent, Minimize>(myFunc, derivativeOf(myFunc), 4.0, 1e-12, 10000);
    std::cout << "Optimized value: " << result << std::endl;

    auto optimizer = AutoSearch< decltype(myFunc), double, Minimize >(myFunc, { 4.0, 5.0 });

    for (int i = 0; i < 100; ++i) {
        auto current = optimizer.current();
        std::cout << "Iteration " << i << ": " << current.first << "\t" << current.second << std::endl;
        optimizer.iterate();
    }

    Vec x = linspace(-20.0, 5.0, 200);

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

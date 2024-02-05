#include "_external.hpp"

#include <Deriv.hpp>
#include <Optim.hpp>
#include <Roots.hpp>
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

    // std::cout << std::setprecision(8) << std::fixed;

    auto myFunc = [](double x) {
        return -x * x * x + 4 * x * x - x - 5;
    };    // Example function (quadratic)

    // auto myFunc = [](double x) { return x*x; };    // Example function (quadratic)
    // auto myFunc = [](double x) { return pow(x,6) - 11*pow(x,3)+17*x*x - 7*x+1; };

    std::cout << "Bracketing solver:\n";
    auto guess = fminimize< GoldenSearch >(myFunc, { -1.0, 1.0 }, [](const auto& data) {
        auto [iter, lower, guess, upper] = data;

        if (iter == 0) {
            std::cout << "----------------------------------------------------------------\n";
            std::cout << fmt::format("{:>10} | {:>15} | {:>15} | {:>15} ", "#", "Lower", "Guess", "Upper") << "\n";
            std::cout << "----------------------------------------------------------------\n";
        }

        std::cout << fmt::format("{:10} | {:15.10f} | {:15.10f} | {:15.10f} ", iter, lower, guess, upper) << "\n";

        BracketTerminator term;

        if (term(data)) {
            std::cout << "----------------------------------------------------------------\n";
            return true;
        }
        return false;
    });

    // auto guess = fminimize<Brent>(myFunc, { -1.0, 1.0 });

    std::cout << "Optimized value: " << guess << std::endl;
    std::cout << "Function value: " << myFunc(guess) << std::endl;

    auto result = fdfoptimize< Newton, Minimize >(myFunc, guess);
    std::cout << "Optimized value: " << result << std::endl;
    std::cout << "Function value: " << myFunc(result) << std::endl << std::endl;

    Vec x = linspace(-1.0, 1.0, 200);

    std::vector< double > y_lin;
    rng::transform(x, std::back_inserter(y_lin), [&](auto val) { return myFunc(val); });

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

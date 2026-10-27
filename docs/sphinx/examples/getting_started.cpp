// The program of docs/sphinx/getting-started.md: a root by Brent's method and a derivative by Ridders' method. The
// page shows the parts between the markers; main returns non-zero unless what the page says about the results holds,
// checked against closed forms by the bounds of the tests (tests/roots_check.hpp, tests/test_ridders.cpp).

// clang-format off
// <includes>
#include <numerixx/deriv.hpp>
#include <numerixx/roots.hpp>

#include <cmath>
#include <print>
// </includes>
// clang-format on

#include <limits>
#include <numbers>

int main()
{
    // clang-format off
    // <program>
    namespace rt = nxx::roots;
    auto const root = rt::solve<rt::Brent>([](double x) { return x * x - 2.0; }, rt::Bracket{1.0, 2.0});
    if (root) std::println("root  {:.17g} +- {:.1e}", root->value, root->abs_error);
    else std::println("no root: {}", root.error());

    namespace dv = nxx::deriv;
    auto const slope = dv::derivative<dv::Ridders>([](double x) { return std::exp(x); }, 1.0);
    if (slope) std::println("slope {:.17g} +- {:.1e}", slope->value, *slope->abs_error);
    else std::println("no slope: {}", slope.error());
    // </program>
    // clang-format on

    if (!root || !slope || !slope->abs_error) return 1;
    double const eps = std::numeric_limits<double>::epsilon();

    // Brent: within abs_error of sqrt(2), plus the shift of the sign change from fn's rounding: x * x rounds by at
    // most half an ulp of 2, the subtraction is exact, and |f'| = 2 sqrt(2), so 2.2e-16 / 2.83 < 8e-17. abs_error is
    // at most twice the default tolerance, 2 eps |x|.
    double const sqrt2 = std::numbers::sqrt2;
    if (std::abs(root->value - sqrt2) > root->abs_error + 8e-17) return 2;
    if (root->abs_error > 4 * eps * sqrt2) return 3;

    // Ridders: the rule of tests/test_ridders.cpp on the independent prototype's values for exp at 1
    // (docs/deriv/PHASE3_API.md): relative error 9.8e-16 and estimate 6.0e-14. The error is within 4 times the
    // prototype's and 6 times the estimate, and the estimate within twice the prototype's.
    double const e = std::numbers::e;
    double const error = std::abs(slope->value - e);
    if (error > 4 * 9.8e-16 * e || error > 6 * *slope->abs_error) return 4;
    if (*slope->abs_error > 2 * 6.0e-14 * e) return 5;
    return 0;
}

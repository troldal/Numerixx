// ================================================================================================
// This demo shows various ways to compute 1st and 2nd derivatives numerically using the methods
// in the calculus/Derivatives.hpp header.
// ================================================================================================

#include <Deriv.hpp>
#include <cmath>
#include <fmt/format.h>
#include <iomanip>
#include <iostream>
#include <numbers>

void print(const auto& testcases);

int main()
{
    using namespace nxx::deriv;
    std::cout << std::fixed << std::setprecision(8);
    auto func = [](std::floating_point auto x) { return std::log(x) + 2 * x; };

    // ============================================================================================
    // The most straightforward way of computing the numerical derivatives is to use the 'central'
    // function. This function will calculate the dericative of a function using a central difference
    // method. Alternatively, the 'forward' or 'backward' can be used to compute the derivative using
    // forward difference or backward difference methods, respectively. All three functions use
    // Richardson extrapolation for computing the derivatives.
    // The functions can take any callable function as the first argument, followed by the value at
    // which to compute the derivative. Optionally, the step size for the difference calculation
    // can be provided. If not provided, the cubic root of the machine epsilon will be used as default.
    //
    // Note that the return type is a tl::expected (developed by TartanLlama). When C++23 is generally
    // available, this will be replaced by std::expected. tl::expected contains either the return value
    // or an exception object. (For details, see https://en.cppreference.com/w/cpp/utility/expected)
    //
    // The following code shows how to compute the derivative for the lambda function stored in the
    // variable 'func':
    // ============================================================================================
    std::cout << "\nCompute the numerical derivative of the function ln(x) + 2x at the value 'e'\n";
    std::cout << "using central, forward, and backward difference calculation methods:\n";
    std::cout << "Central difference:  " << *central(func, std::numbers::e) << std::endl;
    std::cout << "Forward difference:  " << *forward(func, std::numbers::e) << std::endl;
    std::cout << "Backward difference: " << *backward(func, std::numbers::e) << std::endl;

    // ============================================================================================
    // The following code does the same as the previous, but the optional step size argument is also
    // provided:
    // ============================================================================================
    std::cout << "\nSame as above, but the optional 'stepsize' argument is provided:\n";
    std::cout << "Central difference:  " << *central(func, std::numbers::e, 1E-2) << std::endl;
    std::cout << "Forward difference:  " << *forward(func, std::numbers::e, 1E-2) << std::endl;
    std::cout << "Backward difference: " << *backward(func, std::numbers::e, 1E-2) << "\n\n";

    // ============================================================================================
    // If more fine-grained control is required, the 'derivative' template function can be used.
    // The 'derivative' function takes as a template argument, any callable type that can be used
    // as an algorithm for computing the derivative of an arbitrary function.
    // The function signature is otherwise identical to the 'central', 'forward', and 'backward'
    // functions; the only difference is that the algorithm for computing the derivative is
    // explicitly provided as a template parameter.
    //
    // A number of algorithms are provided:
    // ============================================================================================
    std::cout << "The following examples show how to use the nxx::deriv::derivative template function\n";
    std::cout << "to manually specify the algorithm used to compute the 1st derivative.\n";
    std::cout << "Order1CentralRichardson:   " << *diff<Order1CentralRichardson>(func, std::numbers::e) << std::endl;
    std::cout << "Order1Central3Point:       " << *diff<Order1Central3Point>(func, std::numbers::e) << std::endl;
    std::cout << "Order1Central5Point:       " << *diff<Order1Central5Point>(func, std::numbers::e) << "\n\n";

    std::cout << "Order1ForwardRichardson:   " << *diff< Order1ForwardRichardson >(func, std::numbers::e) << std::endl;
    std::cout << "Order1Forward2Point:       " << *diff<Order1Forward2Point>(func, std::numbers::e) << std::endl;
    std::cout << "Order1Forward3Point:       " << *diff<Order1Forward3Point>(func, std::numbers::e) << "\n\n";

    std::cout << "Order1BackwardRichardson:  " << *diff< Order1BackwardRichardson >(func, std::numbers::e) << std::endl;
    std::cout << "Order1Backward2Point:      " << *diff<Order1Backward2Point>(func, std::numbers::e) << std::endl;
    std::cout << "Order1Backward3Point:      " << *diff<Order1Backward3Point>(func, std::numbers::e) << "\n\n";

    // ============================================================================================
    // Similarly, algorithms for computing the 2nd derivatives are also provided:
    // ============================================================================================
    std::cout << "Similarly, the following examples show how to use the nxx::deriv::derivative template function\n";
    std::cout << "to manually specify the algorithm used to compute the 2nd derivative.\n";
    std::cout << "Order2Central3Point:      " << *diff<Order2Central3Point>(func, std::numbers::e) << std::endl;
    std::cout << "Order2Central5Point:      " << *diff<Order2Central5Point>(func, std::numbers::e) << "\n\n";

    std::cout << "Order2Forward3Point:      " << *diff<Order2Forward3Point>(func, std::numbers::e) << std::endl;
    std::cout << "Order2Forward4Point:      " << *diff<Order2Forward4Point>(func, std::numbers::e) << "\n\n";

    std::cout << "Order2Backward3Point:     " << *diff<Order2Backward3Point>(func, std::numbers::e) << std::endl;
    std::cout << "Order2Backward4Point:     " << *diff<Order2Backward4Point>(func, std::numbers::e) << "\n\n";

    // ============================================================================================
    // Finally, it is also possible to provide a custom algorithm, as long as it has the correct
    // signature. The following example shows how to use an algorithm implemented as a lambda:
    // ============================================================================================
    std::cout << "The following example show how to use the nxx::deriv::derivative template function\n";
    std::cout << "with a custom algorithm.\n";

    auto algo = [](auto function, double val, double stepsize = 1E-6) {
        return (4.0 * (function(val + stepsize) - function(val - stepsize)) -
                0.5 * (function(val + 2 * stepsize) - function(val - 2 * stepsize))) /
               (stepsize * 6);
    };

    std::cout << "Custom algorithm:         "<< *diff<decltype(algo)>(func, std::numbers::e) << "\n\n";

    // ============================================================================================
    // As a convenience function, the derivativeOf() template function can be used to create a
    // function object representing the derivative of the input function.
    // As a template argument, it is possible to pass any algorithm with the correct signature (see
    // above). If no template argument is given, the function will use the Order1CentralRichardson
    // algorithm as the default.
    // As above, it is possible to pass a custom algorithm and a custom stepsize. Also, the return
    // type of the function object created will be tl::expected (std::expected in the future).
    // The following code example shows how to create function objects for the 1st and 2nd derivatives.
    // ============================================================================================
    auto d1func = derivativeOf(func);
    auto d2func = derivativeOf<Order2Central5Point>(func);
    std::cout << "Derivative function objects using the .derivativeOf() function:\n";
    std::cout << "d1func:                   " << d1func(std::numbers::e) << "\n";
    std::cout << "d2func:                  " <<  d2func(std::numbers::e) << "\n\n";

    // ============================================================================================
    // The following code shows the results of computing the derivatives for 10 different functions
    // numerically, using the central, forward, and backward functions.
    // ============================================================================================
    auto f0 = [](double x) { return std::pow(x, 3) - 2 * x + 5; };
    auto f1 = [](double x) { return 2 * std::pow(x, 2) + 3 * x - 4; };
    auto f2 = [](double x) { return std::sin(x) + std::cos(x); };
    auto f3 = [](double x) { return std::log(x) + 2 * x; };
    auto f4 = [](double x) { return 4 * std::pow(x, 4) - 3 * std::pow(x, 3) + 2 * std::pow(x, 2) - x + 1; };
    auto f5 = [](double x) { return std::exp(x) + 3 * std::pow(x, 2); };
    auto f6 = [](double x) { return std::cos(x * x) - 2 * x; };
    auto f7 = [](double x) { return std::sqrt(x) + 2.0 / x; };
    auto f8 = [](double x) { return 3 * std::pow(x, 3) - 4 * std::pow(x, 2) + 5 * x - 6; };
    auto f9 = [](double x) { return 1.0 / (x + 1); };

    std::string s0 = "x^3 - 2*x + 5";
    std::string s1 = "2*x^2 + 3*x - 4";
    std::string s2 = "sin(x) + cos(x)";
    std::string s3 = "ln(x) + 2*x";
    std::string s4 = "4*x^4 - 3*x^3 + 2*x^2 - x + 1";
    std::string s5 = "exp(x) + 3x^2";
    std::string s6 = "cos(x^2) - 2*x";
    std::string s7 = "sqrt(x) + 2.0 / x";
    std::string s8 = "3*x^3 - 4*x^2 + 5*x - 6";
    std::string s9 = "1.0 / (x + 1)";

    double v0 = 2.0;
    double v1 = 1.0;
    double v2 = std::numbers::pi / 4;
    double v3 = std::numbers::e;
    double v4 = 0.0;
    double v5 = 1.0;
    double v6 = std::numbers::pi;
    double v7 = 4.0;
    double v8 = 2.0;
    double v9 = 0.0;

    double r0 = 10.0;
    double r1 = 7.0;
    double r2 = 0.0;
    double r3 = 2.367879441;
    double r4 = -1.0;
    double r5 = std::numbers::e + 6;
    double r6 = 0.703662284;
    double r7 = 0.125;
    double r8 = 25.0;
    double r9 = -1.0;

    std::cout << "CENTER DERIVATIVE USING RICHARDSON EXTRAPOLATION" << std::endl;
    std::vector< std::tuple< std::string, double, double, std::function<double(double)> > > test1;
    test1.emplace_back(s0, v0, r0, [&](double val){return *central(f0, val);});
    test1.emplace_back(s1, v1, r1, [&](double val){return *central(f1, val);});
    test1.emplace_back(s2, v2, r2, [&](double val){return *central(f2, val);});
    test1.emplace_back(s3, v3, r3, [&](double val){return *central(f3, val);});
    test1.emplace_back(s4, v4, r4, [&](double val){return *central(f4, val);});
    test1.emplace_back(s5, v5, r5, [&](double val){return *central(f5, val);});
    test1.emplace_back(s6, v6, r6, [&](double val){return *central(f6, val);});
    test1.emplace_back(s7, v7, r7, [&](double val){return *central(f7, val);});
    test1.emplace_back(s8, v8, r8, [&](double val){return *central(f8, val);});
    test1.emplace_back(s9, v9, r9, [&](double val){return *central(f9, val);});
    print(test1);

    std::cout << "FORWARD DERIVATIVE USING RICHARDSON EXTRAPOLATION" << std::endl;
    std::vector< std::tuple< std::string, double, double, std::function<double(double)> > > test2;
    test2.emplace_back(s0, v0, r0, [&](double val){return *forward(f0, val);});
    test2.emplace_back(s1, v1, r1, [&](double val){return *forward(f1, val);});
    test2.emplace_back(s2, v2, r2, [&](double val){return *forward(f2, val);});
    test2.emplace_back(s3, v3, r3, [&](double val){return *forward(f3, val);});
    test2.emplace_back(s4, v4, r4, [&](double val){return *forward(f4, val);});
    test2.emplace_back(s5, v5, r5, [&](double val){return *forward(f5, val);});
    test2.emplace_back(s6, v6, r6, [&](double val){return *forward(f6, val);});
    test2.emplace_back(s7, v7, r7, [&](double val){return *forward(f7, val);});
    test2.emplace_back(s8, v8, r8, [&](double val){return *forward(f8, val);});
    test2.emplace_back(s9, v9, r9, [&](double val){return *forward(f9, val);});
    print(test2);

    std::cout << "BACKWARD DERIVATIVE USING RICHARDSON EXTRAPOLATION" << std::endl;
    std::vector< std::tuple< std::string, double, double, std::function<double(double)> > > test3;
    test3.emplace_back(s0, v0, r0, [&](double val){return *backward(f0, val);});
    test3.emplace_back(s1, v1, r1, [&](double val){return *backward(f1, val);});
    test3.emplace_back(s2, v2, r2, [&](double val){return *backward(f2, val);});
    test3.emplace_back(s3, v3, r3, [&](double val){return *backward(f3, val);});
    test3.emplace_back(s4, v4, r4, [&](double val){return *backward(f4, val);});
    test3.emplace_back(s5, v5, r5, [&](double val){return *backward(f5, val);});
    test3.emplace_back(s6, v6, r6, [&](double val){return *backward(f6, val);});
    test3.emplace_back(s7, v7, r7, [&](double val){return *backward(f7, val);});
    test3.emplace_back(s8, v8, r8, [&](double val){return *backward(f8, val);});
    test3.emplace_back(s9, v9, r9, [&](double val){return *backward(f9, val);});
    print(test3);

    return 0;
}

void print(const auto& testcases)
{
    std::cout << "------------------------------------------------------------------------------------------------------\n";
    std::cout << fmt::format("{:>10} | {:>35} | {:>15} | {:>15} | {:>15} ", "#", "Function", "Evaluation Pt.", "True dF", "Calculated dF")
              << std::endl;
    std::cout << "------------------------------------------------------------------------------------------------------\n";

    size_t index = 1;
    for (auto item : testcases) {
        std::cout << fmt::format("{:10} | {:>35} | {:15.8f} | {:15.8f} | {:15.8f} ",
                                 index,
                                 std::get< 0 >(item),
                                 std::get< 1 >(item),
                                 std::get< 2 >(item),
                                 std::get< 3 >(item)(std::get< 1 >(item)))
                  << std::endl;
        ++index;
    }

    std::cout << "------------------------------------------------------------------------------------------------------\n\n";
}

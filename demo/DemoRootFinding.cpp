// ================================================================================================
// This demo shows how to use the nxx::roots namespace to find the roots of arbitrary functions.
// ================================================================================================

#include "_external.hpp"

#include <Deriv.hpp>
#include <Roots.hpp>
#include <algorithm>
#include <array>
#include <concepts>
#include <iomanip>
#include <iostream>

using NXX_FLOAT = double; // boost::multiprecision::cpp_bin_float_50;

// template<>
// struct fmt::formatter< NXX_FLOAT > : fmt::formatter< double >
// {
//     auto format(const NXX_FLOAT& d, fmt::format_context& ctx) const
//     {
//         return fmt::formatter< double >::format(static_cast< double >(d), ctx);
//     }
// };

int main()
{
    using namespace nxx::roots;
    using namespace boost::multiprecision;
    std::cout << std::fixed << std::setprecision(8);
    auto func = []<nxx::IsFloat VAL_T>(VAL_T x) { return x * x - decltype(x)(5.0); };
    constexpr std::pair<NXX_FLOAT, NXX_FLOAT> bounds = { 0.0, 2.5 };

    // ============================================================================================
    // The nxx::roots namespace contains a number of root-finding algorithms, for finding the roots
    // of arbitrary functions. The algorithms are implemented as objects that can be used to iterate
    // towards a solution. The algorithms are:
    //
    // 1. Ridders
    // 2. Bisection
    // 3. Discrete Newton
    // 4. Newton
    //
    // These algorithms can be used directly, by manually iterating towards a solution, or they can
    // be used indirectly, by using the fsolve() or fdfsolve() functions, which will automatically
    // iterate until a solution is found. The fsolve() function can be used with bracketing methods
    // (Ridders and Bisection), while the fdfsolve() function can be used with non-bracketing methods
    // (Discrete Newton and Newton) that require a derivative of the function.
    //
    // The easiest way to use the algorithms is to use the fsolve() and fdfsolve() functions, which
    // will automatically iterate until a solution is found. The following code shows how to use
    // the fsolve() and fdfsolve() functions to find the roots of the polynomial f(x) = x^2 - 5.
    // (Note that to find the roots of a polynomial, it is better to use the polysolve function
    // from the nxx::poly namespace, which is much faster and more accurate than the root-finding
    // algorithms in the nxx::roots namespace. However, as an example of how to use the root-finding
    // algorithms, it will work just fine.)
    //
    // It should be noted that both the fsolve() and fdfsolve() functions returns a tl::expected
    // object (std::expected when C++23 comes in widespread use), which can either contain a value
    // or an error. The error can be checked by calling the has_value() function, and the value can
    // be retrieved by calling the value() function. The value() function will throw an exception
    // if the expected object does not contain a value. The value can also be retrieved by using
    // the * operator, which will return a reference to the value. If the expected object does not
    // contain a value, the result of using the * operator is undefined.
    // ============================================================================================

    auto token = [](const PolishingIterData<size_t, double> &data) {
        auto [iter, guess, previous] = data;
        auto function = []<nxx::IsFloat VAL_T>(VAL_T x) { return x * x - decltype(x)(5.0); };
        auto eval = function(guess);

        if (iter == 0) {
            std::cout << "----------------------------------------------------------------------------------\n";
            std::cout << fmt::format("{:>10} | {:>15} | {:>15} ", "#", "Guess", "Eval") << "\n";
            std::cout << "----------------------------------------------------------------------------------\n";
        }

        std::cout << fmt::format("{:10} | {:15.10f} | {:15.10f} ", iter, guess, eval) << "\n";

        // PolishingStopToken term;

        if (iter >= 5) {
            std::cout << "----------------------------------------------------------------------------------\n";
            return true;
        }
        return false;
    };

    auto outputter = [](const auto &data) -> tl::expected<double, std::string> {
        auto [iter, guess, previous] = data;
        using expected = tl::expected<decltype(guess), std::string>;

        if (iter >= 5) return tl::make_unexpected(std::string("Too many iterations"));
        return expected(guess);
    };

    using Expected = decltype(outputter);
    using Token = decltype(token);

    std::cout << "\nCompute the root of the polynomial f(x) = x^2 - 5 using bracketing methods:\n";
    std::cout << "Bisection Method:         " << fsolve<Bisection>(func, { 0.0, 2.5 }).result() << std::endl;
    std::cout << "Ridder's Method:          " << fsolve<Ridder>(func, bounds).result() << std::endl;
    std::cout << "Regula Falsi Method:      " << fsolve<RegulaFalsi>(func, bounds, 1e-6, 100).result() << std::endl
              << std::endl;

    std::cout << "\nCompute the root of the polynomial f(x) = x^2 - 5 using polishing methods:\n";
    std::cout << "Newton's Method:          " << fdfsolve<Newton>(func, 1.25).result() << std::endl;
    std::cout << "Secant Method:            " << fdfsolve<Secant>(func, 1.25).result() << std::endl;
    std::cout << "Steffensen's Method:      " << fdfsolve<Steffensen>(func, 1.25).result() << std::endl << std::endl;

    std::cout << "Newton's Method:          \n"
              << fdfsolve<Newton>(func, 1.25, token).result<Expected>().error() << std::endl;
    std::cout << "Newton's Method:          \n"
              << fdfsolve<Newton, Token>(func, 1.25).result<Expected>().error() << std::endl;
    // std::cout << "Secant Method:            \n"
    //           << fdfsolve<Secant, Token>(func, 1.25).result<Expected>().error() << std::endl;
    // std::cout << "Steffensen's Method:      \n"
    //           << fdfsolve<Steffensen, Token>(func, 1.25).result<Expected>().error() << std::endl
    //           << std::endl;

    // ============================================================================================
    // If more fine-grained control is needed, the algorithms can be used directly. Both the bracketing
    // solvers and the polishing solvers behaves similarly. The evaluate() member function can be used
    // to evaluate the function at a given point (to check if convergence has been reached), and the
    // iterate() member function can be used to proceed with one iteration. The result() member function
    // can be used to retrieve the current bracketing interval or the current guess, depending on the
    // type of solver. The init() member function can be used to initialize the solver with a bracketing
    // interval or a guess, depending on the type of solver.
    // The following code shows how to use the bracketing solvers and the polishing solvers directly.
    // ============================================================================================



    return 0;
}

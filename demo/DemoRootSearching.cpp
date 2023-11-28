// ================================================================================================
// This demo shows how to use the nxx::roots namespace to search for a bracketing interval
// for a root of a function.
// ================================================================================================

#include <fmt/format.h>
#include <iostream>
#include <iomanip>
//#include <numerixx.hpp>
#include <Roots.hpp>

int main()
{
    using namespace nxx::roots;
    std::cout << std::fixed << std::setprecision(8);
    auto func = nxx::poly::Polynomial({ -5.0, 0.0, 1.0 });

    // ============================================================================================
    // The nxx::roots namespace contains a number of search algorithms, for finding a bracketing
    // interval for a root of a function. The algorithms are implemented as objects that can be
    // used to iterate towards a solution, until a bracketing interval is found. The algorithms
    // are:
    //
    // 1. BracketExpandUp
    // 2. BracketExpandDown
    // 3. BracketExpandOut
    // 4. BracketSearchUp
    // 5. BracketSearchDown
    // 6. BracketSubdivide
    //
    // The BracketExpandXX algorithms will expand the bracketing interval in the positive or
    // negative direction (or both), until a root is found. The BracketSearchXX algorithms will search
    // for a bracketing interval by moving the initial bracket in the positive or negative direction,
    // until a root is found. The BracketSubdivide algorithm will subdivide the bracketing interval
    // until a root is found. The BracketSubdivide is particularly useful when the bracketing
    // interval is known to contain a root, but the function is not known to be monotonic, and may
    // have multiple roots in the interval.
    //
    // The easiest way to use the algorithms is to use the search() and providing it with an
    // algorithm object, and the initial bracketing interval. The search() function will return
    // a tl::expected object, which will contain the bracketing interval if a root was found, or
    // an error code if no root was found. When C++23 is available, the tl::expected object will
    // be replaced with std::expected. The bracketing interval will be returned as a std::pair.
    // To get the bracketing interval from the tl::expected object, use the * operator, or the
    // value() method.
    //
    // The following code shows how to use the search algorithms to find a bracketing interval
    // for a root of a function.
    // ============================================================================================

    auto PrintBounds = [](auto b) { return std::string("(" + std::to_string(b.first) + ", " + std::to_string(b.second) + ")"); };

    std::cout << "\nIdentify the brackets around the root of the polynomial: " << to_string(func) << "\n";
    std::cout << "BracketExpandUp Method:         " << PrintBounds(*search< BracketExpandUp >(func, { 1.0, 1.1 })) << std::endl;
    std::cout << "BracketSearchUp Method:         " << PrintBounds(*search< BracketSearchUp >(func, { 1.0, 1.1 })) << std::endl;
    std::cout << "BracketExpandDown Method:       " << PrintBounds(*search< BracketExpandDown >(func, { 4.9, 5.0 })) << std::endl;
    std::cout << "BracketSearchDown Method:       " << PrintBounds(*search< BracketSearchDown >(func, { 4.9, 5.0 })) << std::endl;
    std::cout << "BracketExpandOut Method:        " << PrintBounds(*search< BracketExpandOut >(func, { 1.0, 1.1 })) << std::endl;
    std::cout << "BracketSubdivide Method:        " << PrintBounds(*search< BracketSubdivide >(func, { -5.0, 10.0 })) << std::endl;

    // The given polynomial has two roots, one at x = -2.23606798, and the other at x = 2.23606798. The examples above finds the
    // bracketing interval for the root at x = 2.23606798, except for the BracketSubdivide method, which finds the bracketing
    // interval for the root at x = -2.23606798. The reason is that the BracketSubdivide method always will find the lowest
    // root first.

    // The search() function also has two optional parameters, searchFactor, and maxiter. The searchFactor parameter is used
    // to control the rate at which the bracketing interval is expanded or subdivided. The default value is the golden ratio
    // (phi = ~1.6), but can take any value >= 1.0. The maxiter parameter is used to limit the number of iterations. The default
    // value is 100, but can take any value >= 1. The following code shows how to use the searchFactor and maxiter parameters.
    // The following code uses a searchFactor of 2.0, and a maxiter of 10.

    std::cout << "\nIdentify the brackets around the root of the polynomial, using a searchFactor of 2.0 and maxiter = 10: "
        << "\n";
    std::cout << "BracketExpandUp Method:         " << PrintBounds(*search< BracketExpandUp >(func, { 1.0, 1.1 }, 2.0, 10)) << std::endl;

    // ============================================================================================
    // As mentioned above, the search() function will return a tl::expected object, which will
    // contain the bracketing interval if a root was found, or an error code if no root was found.
    // This is useful if the user only wants to know if a root was found, or not.
    // ============================================================================================

    std::cout << "\nFind the brackets around the root of the function f(x) = log(x):\n\n";
    auto print_error = [](auto error) {
        std::cout << "Error Description: " << error.what() << std::endl;
        std::cout << "Error Type:        " << error.typeAsString() << std::endl;
        std::cout << "Last Value:        [" << error.value().first << ", " << error.value().second << "]" << std::endl;
        std::cout << "Iterations:        " << error.iterations() << std::endl << std::endl;
    };

    std::cout << "Initial Bracket:   [5.0, 10.0] (expanding down)\n";
    auto root = search< BracketExpandDown >([](double x) { return std::log(x); }, { 5.0, 10.0 });
    if (!root) print_error(root.error());

    std::cout << "Initial Bracket:   [5.0, 10.0] (expanding up)\n";
    root = search< BracketExpandUp >([](double x) { return std::log(x); }, { 5.0, 10.0 }, 1.0, 10);
    if (!root) print_error(root.error());

    // The error object is a subclass of the RootError class, which is a subclass of the std::runtime_error
    // class. This means that the error object can be used in a try-catch block, as shown below.
    try {
        throw root.error();
    }
    catch (const RootError& e) {
        std::cout << "Exception caught: " << e.what() << std::endl << std::endl;
    }

    // ============================================================================================
    // If more fine-grained control is needed, the algorithms can be used directly. The algorithms
    // are implemented as objects that can be used to iterate towards a solution, until a bracketing
    // interval is found.
    //
    // The following code shows how to use the algorithms directly.
    // ============================================================================================

    // Ad hoc lambda function to print the current iteration:
    auto find_bracket = [](auto solver) {
        // Print the header:
        std::cout << "----------------------------------------------------------------\n";
        std::cout << fmt::format("{:>10} | {:>15} | {:>15} | {:>15} ", "Iter", "Lower", "Upper", "Factor") << std::endl;
        std::cout << "----------------------------------------------------------------\n";

        // Iterate until convergence (or until 100 iterations have been performed):
        for (int i = 0; i <= 100; ++i) {
            // Print the current iteration:
            std::cout << fmt::format("{:10} | {:15.10f} | {:15.10f} | {:15.10f} ",
                                     i,                       // Iteration number
                                     solver.current().first,  // Lower bound
                                     solver.current().second, // Upper bound
                                     solver.ratio())          // Expansion factor
                << std::endl;

            // Check if a root is in the current interval:
            if (solver.evaluate(solver.current().first) * solver.evaluate(solver.current().second) < 0.0) break;

            // Perform one iteration:
            solver.iterate();
        }

        // Print the final result:
        std::cout << std::fixed << std::setprecision(10)
            << "SUCCESS! Bounds found at: (" << solver.current().first << ", "
            << solver.current().second << ")\n";
        std::cout << "----------------------------------------------------------------\n\n";
    };

    std::cout << "\nIdentify the brackets around the root of the polynomial, using the algorithms directly: "
              << "\n\n";

    std::cout << "BracketExpandUp Method:         " << std::endl;
    find_bracket(BracketExpandUp(func, { 1.0, 1.1 }));

    std::cout << "BracketSearchUp Method:         " << std::endl;
    find_bracket(BracketSearchUp(func, { 1.0, 1.1 }));

    std::cout << "BracketExpandDown Method:         " << std::endl;
    find_bracket(BracketExpandDown(func, { 4.9, 5.0 }));

    std::cout << "BracketSearchDown Method:         " << std::endl;
    find_bracket(BracketSearchDown(func, { 4.9, 5.0 }));

    std::cout << "BracketExpandOut Method:         " << std::endl;
    find_bracket(BracketExpandOut(func, { 1.0, 1.1 }));

    std::cout << "BracketSubdivide Method:         " << std::endl;
    find_bracket(BracketSubdivide(func, { -5.0, 10.0 }));

    return 0;
}
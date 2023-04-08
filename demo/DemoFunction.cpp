// ================================================================================================
// This demo shows how to use the Function class to evaluate a function object.
// ================================================================================================

#include <cmath>
#include <numerixx.hpp>

int main()
{
    using namespace nxx;
    using namespace std::complex_literals;
    auto func = [](auto x) { return std::log(x); };

    // ============================================================================================
    // The Function class is a wrapper class for functions that accepts and returns either floating point
    // values or a std::complex. The Function class can be used to evaluate a function object for a given
    // input value. This class has the same interface as any other function object, so one might wonder
    // why it is needed in the first place. The main reason is that the Function class also contains the
    // logic to handle errors that might occur during the evaluation of the function object. For example,
    // if the function object is not defined for a given input value, the Function class will return an
    // error instead of crashing the program. For this reason, the Function class is mostly used as an
    // internal class in the library, but it can also be used directly by the user.
    // The following example shows how to use the Function class to evaluate a function object. In this case,
    // the function object is a lambda function that computes the natural logarithm of its input, but
    // any callable object can be used.
    // ============================================================================================

    //  Create a function object
    auto f = func::Function(func);

    //  Evaluate the function object for a given input floating point value
    auto x       = 0.0;
    auto result1 = f.evaluate(x);

    //  Check if the evaluation was successful
    if (result1) {
        //  Print the result
        std::cout << "Result: " << result1.value() << std::endl;
    }
    else {
        //  Print the error message
        std::cout << "Error: " << result1.error().what() << std::endl;
    }

    //  Evaluate the function object for a given input complex value
    auto z       = 1.0 + 1.0i;
    auto result2 = f.evaluate(z);

    //  Check if the evaluation was successful
    if (result2) {
        //  Print the result
        std::cout << "Result: " << result2.value() << std::endl;
    }
    else {
        //  Print the error message
        std::cout << "Error: " << result2.error().what() << std::endl;
    }

    return 0;
}
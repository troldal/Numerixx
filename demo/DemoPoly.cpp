// ================================================================================================
// This demo shows how to use the nxx:poly::Polynomial class to create polynomial function objects.
// ================================================================================================

#include <Poly.hpp>
#include <deque>
#include <iomanip>
#include <iostream>
#include <list>
#include <set>
#include <unordered_set>
#include <fmt/format.h>

void printpoly(auto coefficients) {

    using namespace nxx::poly;
    auto func = Polynomial(coefficients);

    std::cout << func(-2) << "\t";
    std::cout << func(-1) << "\t";
    std::cout << func(0) << "\t";
    std::cout << func(1) << "\t";
    std::cout << func(2) << "\n\n";
}

int main() {

    using namespace nxx::poly;
    using namespace std::complex_literals;

    // ============================================================================================
    // Creating a polynomial function object is done by passing a list of the polynomial coefficients
    // to the Polynomial constructor. The list can be any container supporting bidirectional
    // iterators, as well as std::initializer_list. The value type must be a floating point type,
    // or a std::complex object with a floating point value type.
    //
    // The creted polynomial will be in the form c0 + c1*x + c2*x^2 + ... + cn*x^n, meaning
    // that the coefficients must be in the order of increasing power. Note that any trailing
    // zero-value coefficients will be ignored.
    //
    // The resulting object will be a function object, meaning that it has an overloaded function
    // call operator. As an alternative, the .evaluate() function can be used to evaluate the
    // polynomial at a given point. Note, however, that the .evaluate() function is slightly
    // different, in that it will return a tl::expected object (or std::expected when C++23 comes
    // in common use), which can be used to check for errors. The get the actual value, the
    // .value() function or the * operator must be called on the returned object.
    //
    // As a convenience, the Polynomial class also has a .asString() function, for creating
    // a std::string with the polynomial in textual form.
    //
    // In the following example, we create a polynomial function object using real coefficients,
    // but creating a polynomial with complex coefficients is just as easy.
    //
    // The return type will depend on the value type of the coefficients, as well as the type of
    // the argument passed to the function call operator. If the coefficients are of floating-point
    // type, and the argument is of floating-point type, the return type will also be of floating-point
    // type. If the coefficients are of floating-point type, and the argument is of std::complex type,
    // the return type will be of std::complex type. If the coefficients are of std::complex type,
    // the return type will always be of std::complex type.
    // ============================================================================================
    auto func1 = Polynomial({ 1.0, 2.0, 3.0, 4.0 });
    std::cout << "Created the polynomial f(x) = " << func1.asString() << std::endl << std::endl;

    std::cout << "Evaluation using function call operator:" << std::endl;
    std::cout << "Evaluation at -1.0: f(x) = " << func1(-1.0) << std::endl;
    std::cout << "Evaluation at 0.0:  f(x) = " << func1(0.0) << std::endl;
    std::cout << "Evaluation at 1.0:  f(x) = " << func1(1.0) << std::endl << std::endl;

    std::cout << "Evaluation using the .evaluate() function:" << std::endl;
    std::cout << "Evaluation at -1.0: f(x) = " << *func1.evaluate(-1.0) << std::endl;
    std::cout << "Evaluation at 0.0:  f(x) = " << *func1.evaluate(0.0) << std::endl;
    std::cout << "Evaluation at 1.0:  f(x) = " << *func1.evaluate(1.0) << std::endl << std::endl;

    // ============================================================================================
    // While the coefficients may be of floating-point type, the function object can be called
    // with an argument that is either a floating-point type or a std::complex object.
    // ============================================================================================
    std::cout << "Evaluation using function call operator, with std::complex argument:" << std::endl;
    std::cout << "Evaluation at (-0.072 - 0.638i): f(x) = " << func1(-0.0720852-0.638327i) << "\n\n";

    std::cout << "Evaluation using the .evaluate() function, with std::complex argument:" << std::endl;
    std::cout << "Evaluation at (-0.072 - 0.638i): f(x) = " << *func1.evaluate(-0.0720852-0.638327i) << "\n\n";

    // ============================================================================================
    // The coefficients of the polynomial can be retrieved using the .coefficients() function.
    // This function will return a const reference to a std::vector with the coefficients.
    // It is also possible to store the coefficients in a different (bidirectional) container
    // by passing the type as a template parameter to the function.
    // ============================================================================================
    std::cout << "Getting the polynomial coefficients as a std::vector:" << std::endl;
    auto coeffVector = func1.coefficients();
    for (auto c : coeffVector) std::cout << c << " ";
    std::cout << "\n\n";

    std::cout << "Getting the polynomial coefficients as a std::deque:" << std::endl;
    auto coeffDeque = func1.coefficients<std::deque<double>>();
    for (auto c : coeffDeque) std::cout << c << " ";
    std::cout << "\n\n";

    // ============================================================================================
    // As differentiation of polynomials is straightforward, this operation is also supported by the
    // Polynomial class. The derivativeOf function will return a new Polynomial function object
    // that represents the derivative function.
    // ============================================================================================
    auto d1func1 = derivativeOf(func1);
    std::cout << "Derivative of the function: f'(x) = " << d1func1.asString() << std::endl;
    std::cout << "\n";

    // ============================================================================================
    // The Polynomial class also supports addition, subtraction, multiplication, and division of
    // polynomials. The result is a new Polynomial object, which can be used in the same way as
    // the original polynomial.
    // While addition, subtraction, and multiplication operations are straightforward, division
    // is a bit more complex, as polynomial division may result in a non-zero remainder. The
    // overloaded division operator will return the quotient, whereas the modulo operator will
    // return the remainder. If both the quotient and remainder are needed, the divide() function
    // can be used, which will return a std::pair object with both the quotient and remainder.
    // ============================================================================================
    auto dividend = Polynomial({ 5.0, -3.0, 4.0, -1.0, 2.0 });
    auto divisor = Polynomial({ 1.0, 0.0, 1.0 });
    auto quotient = dividend / divisor;
    auto remainder = dividend % divisor;
    auto temp = quotient * divisor + remainder;

    std::cout << "Dividend: " << dividend.asString() << std::endl;
    std::cout << "Divisor: " << divisor.asString() << std::endl;
    std::cout << "Quotient: " << quotient.asString() << std::endl;
    std::cout << "Remainder: " << remainder.asString() << std::endl;
    std::cout << "Quotient * Divisor + Remainder: " << temp.asString() << std::endl << std::endl;

    // ============================================================================================
    // Finally, the nxx::poly namespace contains function for finding the roots of a polynomial.
    // The polyroots() function will return a tl::expected (std::expected) object containint a
    // std::vector with the roots of the polynomial.
    // This is the easiest way to find the roots of a polynomial, as it will find all the roots
    // of the polynomial, including complex roots, regardless of the order of the polynomial. In
    // addition, the quadratic() and cubic() functions can be used to find the roots of quadratic
    // and cubic polynomials, respectively. These functions will also return a tl::expected (std::expected)
    //
    // All three functions accepts an optional template argument that specifies the type of the
    // roots. The default type is the same type as the coefficients, but it is possible explicitly
    // specify the type, e.g. double or std::complex<double>. If the type is specified as a floating-
    // point type, the roots comlex roots will be ignored.
    //
    // Note that when getting the std::vector contained in the tl::expected object (using the .value()
    // function or the * operator), the vector will be returned by reference. This can cause problems
    // if the tl::expected object is a temporary object, as the vector will be destroyed when the
    // temporary object goes out of scope. This can happen when iterating over the roots using a
    // range-based for loop. To address this, the range-based for loop with initializer can be used,
    // as this will extend the lifetime of the temporary object until the end of the loop. The example
    // below shows how this can be done.
    // ============================================================================================
    auto poly1 = Polynomial({-1.0, 0.0, 0.0, 0.0, 0.0, 1.0}); // Polynomial with real coefficients and complex roots

    std::cout << "Real roots of the polynomial (automatically deduced): " << poly1.asString() << std::endl;
    for (auto res = polysolve(poly1); auto root : *res) std::cout << root << "\n"; // Will print the real roots (same type as coefficients)

    std::cout << "\nAll roots of the polynomial (specified by template parameter): " << poly1.asString() << std::endl;
    for (auto res = polysolve<std::complex<double> >(poly1); auto root : *res) std::cout << root << "\n"; // Will print all roots (complex and real)

    std::cout << "\nReal roots of the polynomial (specified by template parameter): " << poly1.asString() << std::endl;
    for (auto res = polysolve<double>(poly1); auto root : *res) std::cout << root << "\n"; // Will print the real roots (as specified by the template argument)

    auto poly2 = Polynomial({-1.0+0i, 0.0+0i, 0.0+0i, 0.0+0i, 0.0+0i, 1.0+0i}); // Polynomial with complex coefficients and complex roots

    std::cout << "\nAll roots of the complex polynomial (automatically deduced): " << poly2.asString() << std::endl;
    for (auto res = polysolve(poly2); auto root : *res) std::cout << root << "\n"; // Will print all the roots (same type as coefficients)

    std::cout << "\nAll roots of the complex polynomial (specified by template parameter): " << poly2.asString() << std::endl;
    for (auto res = polysolve<std::complex<double>>(poly2); auto root : *res) std::cout << root << "\n"; // Will print all roots (as specified by the template argument

    std::cout << "\nReal roots of the complex polynomial (specified by template parameter): " << poly2.asString() << std::endl;
    for (auto res = polysolve<double>(poly2); auto root : *res) std::cout << root << "\n"; // Will print the real roots (complex roots will be ignored)

    try {
        nlohmann::ordered_json data;
        data["Description"] = "Polynomial evaluation failed; non-finite result";
        data["Argument"] = std::nan("");
        data["Coefficients"] = poly2;

        throw nxx::NumerixxError("Integer error", nxx::NumerixxErrorType::General, data.dump());
    }
    catch (const nxx::NumerixxError& e) {
        std::cout << e.log() << std::endl;
    }

    return 0;
}

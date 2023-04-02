// ================================================================================================
// This demo shows how to use the nxx:poly::Polynomial class to create polynomial function objects.
// ================================================================================================

#include "poly/Polyroots.hpp"
#include <complex>
#include <deque>
#include <iomanip>
#include <iostream>
#include <list>
#include <poly/Polynomial.hpp>
#include <set>
#include <unordered_set>

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

    // ============================================================================================
    // Creating a polynomial function object is done by passing a list of the polynomial coefficients
    // to the Polynomial constructor. The list can be any container supporting bi-directional
    // iterators, as well as std::initializer_list. The value type must be a floating point type.
    //
    // The creted polynomial will be in the form c0 + c1*x + c2*x^2 + ... + cn*x^n, meaning
    // that the coefficients must be in the order of increasing power. Note that any trailing
    // zero-value coefficients will be ignored.
    //
    // The resulting object will be a function object, meaning that it has an overloaded function
    // call operator. Alternatively, the .evaluate() function can be used.
    //
    // As a convenience, the Polynomial class also has a .asString() function, for creating
    // a std::string with the polynomial in textual form.
    // ============================================================================================
    auto func1 = Polynomial({ 1.0, 2.0, 3.0, 4.0 });
    std::cout << "Created the polynomial f(x) = " << func1.asString() << std::endl << std::endl;

    std::cout << "Evaluation using function call operator:" << std::endl;
    std::cout << "Evaluation at -1.0: f(x) = " << func1(-1.0) << std::endl;
    std::cout << "Evaluation at 0.0:  f(x) = " << func1(0.0) << std::endl;
    std::cout << "Evaluation at 1.0:  f(x) = " << func1(1.0) << std::endl << std::endl;

    std::cout << "Evaluation using the .evaluate() function:" << std::endl;
    std::cout << "Evaluation at -1.0: f(x) = " << func1.evaluate(-1.0) << std::endl;
    std::cout << "Evaluation at 0.0:  f(x) = " << func1.evaluate(0.0) << std::endl;
    std::cout << "Evaluation at 1.0:  f(x) = " << func1.evaluate(1.0) << std::endl << std::endl;

    // ============================================================================================
    // While the coefficients must be of floating-point type, the function object can be called
    // with an argument that is either a floeting-point type or a std::complex object. Note that
    // the coefficients cannot be std::complex objects.
    // ============================================================================================
    std::cout << "Evaluation using function call operator, with std::complex argument:" << std::endl;
    std::cout << "Evaluation at (-0.072 - 0.638i): f(x) = " << func1(std::complex<double>(-0.0720852, -0.638327)) << "\n\n";

    std::cout << "Evaluation using the .evaluate() function, with std::complex argument:" << std::endl;
    std::cout << "Evaluation at (-0.072 - 0.638i): f(x) = " << func1.evaluate(std::complex<double>(-0.0720852, -0.638327)) << "\n\n";


    // ============================================================================================
    // The coefficients of the polynomial can be retrieved using the .coefficients() function.
    // This function will return a const reference to a std::vector with the coefficients.
    // It is also possible to store the coefficients in a different (bi-directional) container
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
    // Polynomial class. The .derivative() member function will return a new Polynomial function
    // object that represents the derivative function.
    // ============================================================================================
    auto d1func1 = derivativeOf(func1);
    std::cout << "Derivative of the function: f'(x) = " << d1func1.asString() << std::endl;
    std::cout << "\n\n";

    using namespace std::complex_literals;

    auto a_1 = Polynomial({9.0, -2.0, 0.0, 5.0});
    auto b_1 = Polynomial({1.0, 0.0, -1.0, 2.0});
    auto c_1 = a_1 * b_1;
    for (auto item : c_1.coefficients()) std::cout << item << " ";
    std::cout << "\n\n";

    auto rpoly = Polynomial({ 1.0, 2.0, 3.0, 4.0 });
    auto cpoly = Polynomial({ 1.0+0i, 2.0+0i, 3.0+0i, 4.0+0i });

    std::cout << "RPoly(1.0) = " << rpoly(1.0) << std::endl;
    std::cout << "RPoly(1.0 + 0.0i) = " << rpoly(std::complex{1.0, 0.0}) << std::endl;
    std::cout << "CPoly(1.0) = " << cpoly(1.0) << std::endl;
    std::cout << "CPoly(1.0 + 0.0i) = " << cpoly(std::complex{1.0, 0.0}) << std::endl;

    auto a_2 = Polynomial({-1.0, 2.0, -4.0, 3.0, 1.0});
    auto b_2 = Polynomial({1.0, -1.0, 1.0});
    auto c_2 = a_2 / b_2;
    for (auto item : c_2.coefficients()) std::cout << item << " ";
    std::cout << "\n";
    auto c_3 = a_2 % b_2;
    for (auto item : c_3.coefficients()) std::cout << item << " ";
    std::cout << "\n\n";


    using CoefficientList = std::vector<double>;

    auto coefficients = std::vector< CoefficientList > { { 5.7, 1.6, -3.2, 2.5 },
                                                         { -1.2, 0.8, -2.7, 4.0 },
                                                         { 2.4, -1.8, 3.2, -2.1, 5.3 },
                                                         { 3.9, -0.8, 2.6, -1.2, 3.7 },
                                                         { -0.9, 0.8, -2.4, 1.6, -3.1, 2.8 },
                                                         { 1.7, 0.6, -2.6, 3.2, -1.9, 4.5 },
                                                         { 0.7, -1.1, 2.4, -3.2, 1.8, -4.5, 6.2 },
                                                         { 0.6, -0.9, 1.8, -2.3, 1.2, -3.7, 2.1 },
                                                         { -0.6, 0.8, -1.2, 2.7, -3.5, 1.8, -4.2, 5.5 } };

    auto func0 = Polynomial({-1.0, 0.0, 0.0, 0.0, 0.0, 1.0});
    for (auto root : polysolve(func0)) std::cout << root << "\n";
    for (auto root : polysolve<std::complex<double> >(func0)) std::cout << root << "\n";
    for (auto root : polysolve<double>(func0)) std::cout << root << "\n";

    auto func2 = Polynomial({-1.0+0i, 0.0+0i, 0.0+0i, 0.0+0i, 0.0+0i, 1.0+0i});
    for (auto root : polysolve(func2)) std::cout << root << "\n";
    for (auto root : polysolve<std::complex<double> >(func2)) std::cout << root << "\n";
    for (auto root : polysolve<double>(func2)) std::cout << root << "\n";


//
//    std::cout << guess << std::endl;
//    std::cout << a << std::endl;
//    std::cout << *func(guess) << std::endl;
//
//    for (auto res : polysolve(Polynomial({5.7, 1.6, -3.2}))) std::cout << res << " ";
//    std::cout << std::endl;

//    std::cout << "Solving f(x) = -5 + x^2\n";
//    auto fun1 = Polynomial({-5.0, 0.0, 1.0});
//
//    std::cout << "Coefficients: ";
//    for (auto c : fun1.coefficients()) std::cout << c << " ";
//    std::cout << std::endl;
//
//    std::cout << "Roots: ";
//    for (auto r : quadratic(fun1)) std::cout << r << " ";
//    std::cout << std::endl;
//
//    std::cout << "Derivatives at Roots: ";
//    for (auto r : quadratic(fun1)) std::cout << fun1.derivative(r) << " ";
//    std::cout << std::endl;
//    std::cout << std::endl;
//
//    std::cout << "Solving f(x) = -2 - 1.5x + 0.75x^2 + 0.25x^3\n";
//    auto fun2 = Polynomial({-2.0, -1.5, 0.75, 0.25});
//
//    std::cout << "Coefficients: ";
//    for (auto c : fun2.coefficients()) std::cout << c << " ";
//    std::cout << std::endl;
//
//    std::cout << "Roots: ";
//    for (auto r : cubic(fun2)) std::cout << r << " ";
//    std::cout << std::endl;
//
//    std::cout << "Derivatives at Roots: ";
//    for (auto r : cubic(fun2)) std::cout << fun2.derivative(r) << " ";
//    std::cout << std::endl;
//    std::cout << std::endl;
//
//    auto cont = fun1.coefficients<std::deque<double>>();
//    for (auto c : cont) std::cout << c << std::endl;

//    std::vector<double> c {-5, 0, 1};
//    auto fun = numerix::poly::polynomial(c);

//    std::cout << fun(2.23606797749978980505) << std::endl;
//    std::cout << fun(0) << std::endl;
//    std::cout << fun.derivative(-3.0) << std::endl;
//
//    //std::cout << std::fixed << std::setprecision(20);
//    std::cout << numerix::polyroots::quadratic(fun).front() << "\t" << numerix::polyroots::quadratic(fun).back() << std::endl;
//
//    auto fun2 = numerix::poly::Polynomial({-2.0, -1.5, 0.75, 0.25});
//    auto res = numerix::polyroots::cubic(fun2);
//
//    for (auto r : res)
//        std::cout << r << std::endl;
//
//    auto cof = fun.coefficients();
//    for (auto i : cof) std::cout << i << " ";

    return 0;
}
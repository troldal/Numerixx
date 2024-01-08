//
// Created by Kenneth Balslev on 09/12/2022.
//

#include <deque>
#include <iostream>
#include <list>
#include <Multiroots.hpp>
#include <vector>
#include <cmath>
#include <numbers>
#include <span>
#include <blaze/Blaze.h>

/**
 * @brief Primary template for deducing function traits of callables.
 *
 * This structure template is used to deduce the return type and argument type
 * of a callable object, particularly focusing on single-parameter functions.
 * It is specialized for function pointers and member function pointers.
 *
 * @tparam Func The callable object type.
 */
template <typename Func>
struct function_traits : function_traits<decltype(&Func::operator())> {};

/**
 * @brief Specialization of function_traits for function pointers with a single parameter.
 *
 * This specialization deduces the return type and argument type of a single-parameter
 * function pointer.
 *
 * @tparam R Return type of the function.
 * @tparam Arg Argument type of the function.
 */
template <typename R, typename Arg>
struct function_traits<R(*)(Arg)> {
    using return_type = R;
    using argument_type = Arg;
};

/**
 * @brief Specialization of function_traits for member function pointers (including lambdas) with a single parameter.
 *
 * This specialization deduces the return type and argument type of a member function
 * pointer or a lambda function with a single parameter.
 *
 * @tparam R Return type of the function.
 * @tparam C Class type of the member function pointer.
 * @tparam Arg Argument type of the function.
 */
template <typename R, typename C, typename Arg>
struct function_traits<R(C::*)(Arg) const> {
    using return_type = R;
    using argument_type = Arg;
};

/**
 * @brief Creates a std::function from a callable with specific return and argument types.
 *
 * This function template creates a std::function object from a given callable.
 * It statically asserts that the return type must be a floating-point type (float, double, or long double)
 * and the argument must be a std::span of a floating-point type.
 *
 * @tparam Func The callable object type.
 * @param f The callable object.
 * @throws std::static_assert if the return type or argument type does not meet the requirements.
 * @return std::function with deduced return and argument types.
 */
template <typename Func>
auto make_function(Func f) {
    using traits = function_traits<Func>;
    using return_type = typename traits::return_type;
    using argument_type = typename traits::argument_type;

    static_assert(std::is_same_v<return_type, float> || std::is_same_v<return_type, double> || std::is_same_v<return_type, long double>,
                  "Return type must be float, double, or long double");

    static_assert(std::is_same_v<argument_type, std::span<float>> || std::is_same_v<argument_type, std::span<double>> || std::is_same_v<argument_type, std::span<long double>>,
                  "Argument must be a std::span of float, double, or long double");

    return std::function<return_type(argument_type)>{f};
}

int main()
{
    using namespace nxx::deriv;
    using namespace nxx::multiroots;

    auto f1 = []( std::span< double > coeffs) { return 3 * coeffs[0] - std::cos(coeffs[1]*coeffs[2]) - 0.5; };
    auto f2 = []( std::span< double > coeffs) { return coeffs[0] * coeffs[0] - 81*std::pow(coeffs[1] + 0.1,2)+std::sin(coeffs[2])+1.06; };
    auto f3 = []( std::span< double > coeffs) { return std::exp(-coeffs[0]*coeffs[1]) + 20*coeffs[2] + (10*std::numbers::pi-3)/3; };

    //std::function< double(const std::span< double >)> fx = f1;
//    auto fx = make_function(f1);


//    if constexpr (std::is_same_v<typename decltype(fx)::result_type, double> )
//        std::cout << "double" << std::endl;
//    else
//        std::cout << "not double" << std::endl;
//
//    std::vector args = { 0.1f, 0.1f, -0.1f };
//    std::cout << fx(args) << std::endl;

    auto fx = MultiFunction(f1);
    std::cout << fx({ 0.1, 0.1, -0.1 }) << std::endl;

    std::vector args{ 0.1f, 0.1f, -0.1f };
    std::cout << fx(args) << std::endl;

    blaze::DynamicVector<float> x{ 0.1f, 0.1f, -0.1f };
        std::cout << fx(x) << std::endl;



//    DynamicFunctionArray<double> functions { std::vector< std::function< double(const std::span< double >)> >{ f1, f2, f3 }};

//    DynamicFunctionArray<double> functions { std::vector< std::function< double(const std::vector< double >&)> >{
//        [](const std::vector< double >& coeffs) { return 3 * coeffs[0] - std::cos(coeffs[1]*coeffs[2]) - 0.5; },
//        [](const std::vector< double >& coeffs) { return coeffs[0] * coeffs[0] - 81*std::pow(coeffs[1] + 0.1,2)+std::sin(coeffs[2])+1.06; },
//        [](const std::vector< double >& coeffs) { return std::exp(-coeffs[0]*coeffs[1]) + 20*coeffs[2] + (10*std::numbers::pi-3)/3; }
//    }};

//    auto J2 = jacobian(functions, { 0.1, 0.1, -0.1 });
//    std::cout << J2 << std::endl;


//    std::cout << std::fixed << std::setprecision(20);
//    auto solver = DMultiNewton(functions, { 2.0, 2.0, 2.0 });
//    auto result = multisolve(solver, { 2.0, 2.0, 2.0 }, 1.0e-10, 100);
//    for (auto g : result) std::cout << g << std::endl;
//
//    auto root = functions(result);
//    for (auto r : root) std::cout << r << std::endl;

    return 0;
}

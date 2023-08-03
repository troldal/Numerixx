/*
    o.     O O       o Oo      oO o.OOoOoo `OooOOo.  ooOoOOo o      O o      O
    Oo     o o       O O O    o o  O        o     `o    O     O    o   O    o
    O O    O O       o o  o  O  O  o        O      O    o      o  O     o  O
    O  o   o o       o O   Oo   O  ooOO     o     .O    O       oO       oO
    O   o  O o       O O        o  O        OOooOO'     o       Oo       Oo
    o    O O O       O o        O  o        o    o      O      o  o     o  o
    o     Oo `o     Oo o        O  O        O     O     O     O    O   O    O
    O     `o  `OoooO'O O        o ooOooOoO  O      o ooOOoOo O      o O      o

    Copyright © 2023 Kenneth Troldal Balslev

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the “Software”), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is furnished
    to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef NUMERIXX_MULTIDERIVATIVES_HPP
#define NUMERIXX_MULTIDERIVATIVES_HPP

#include "Derivatives.hpp"

#include <blaze/Blaze.h>

namespace nxx::deriv
{

    namespace impl
    {
        /**
         * @brief This structure serves as a traits class to provide additional information
         * about vectors.
         *
         * @tparam typename... variadic template parameters
         */
        template< typename... >
        struct VectorTraits;

        /**
         * @brief Specialization of the VectorTraits class for DynamicVector<T>.
         *
         * @tparam T Type of the elements in the vector.
         */
        template< typename T >
        struct VectorTraits< blaze::DynamicVector< T > >
        {
            using value_type = T;
        };

        /**
         * @brief This structure serves as a traits class to provide additional information
         * about matrices.
         *
         * @tparam typename... variadic template parameters
         */
        template< typename... >
        struct MatrixTraits;

        /**
         * @brief Specialization of the MatrixTraits class for DynamicMatrix<T>.
         *
         * @tparam T Type of the elements in the matrix.
         */
        template< typename T >
        struct MatrixTraits< blaze::DynamicMatrix< T > >
        {
            using value_type = T;
        };

    }    // namespace impl

    /**
     * @brief This concept checks whether a given function object (FN) is a multivariate function.
     *
     * @tparam FN Type of the function object.
     */
    template< typename FN >
    concept IsMultiFunction = requires(FN fn) {
                                  // clang-format off
                                  { std::floating_point< std::invoke_result_t< FN, std::vector< double > > >};
                                  // clang-format on
                              };

    /**
     * @brief This concept checks whether a given array (ARR) contains multi-argument functions.
     *
     * @tparam ARR Type of the array.
     */
    template< typename ARR >
    concept IsMultiFunctionArray = requires(ARR arr) {
                                       // clang-format off
                                       { arr.begin() } -> std::input_iterator;
                                       { arr.end() } -> std::input_iterator;
                                       { IsMultiFunction< typename ARR::value_type> };
                                       // clang-format on
                                   };

    /**
     * @brief Helper template alias that deduces the return type of the multi-argument function.
     *
     * @tparam T Type of the function.
     */
    template< typename T >
    using MultiReturnType = std::invoke_result_t< T, std::vector< double > >;

    /**
     * @brief Computes the partial derivative of a multivariate function at a specified index.
     *
     * This function computes the partial derivative of a function with multiple arguments,
     * i.e., a function of the form f(x1, x2, ..., xn). The function `partialdiff` uses
     * the provided algorithm template parameter to compute the derivative.
     *
     * @tparam ALGO Algorithm to be used for differentiation. This template parameter can be
     * any type that satisfies the requirements for a differentiation algorithm.
     *
     * @param function Multivariate function to differentiate. This is a callable object that
     * takes a vector of doubles as argument and returns a floating-point number.
     *
     * @param args Arguments to the function. This parameter can be of any type that can be
     * converted into a vector of floating-point numbers.
     *
     * @param index Index at which the partial derivative should be computed. This parameter
     * specifies the variable with respect to which the derivative is taken.
     *
     * @param stepsize (Optional) The step size for the derivative. This parameter determines
     * the accuracy of the derivative computation. By default, it uses a function `StepSize`
     * to compute an appropriate step size.
     *
     * @return The partial derivative of the function with respect to the variable at the specified index.
     *
     * @note The returned derivative is evaluated at the provided argument values.
     */
    template< typename ALGO = Order1CentralRichardson >
    inline auto partialdiff(IsMultiFunction auto                  function,
                            auto                                  args,
                            size_t                                index,
                            MultiReturnType< decltype(function) > stepsize = StepSize< MultiReturnType< decltype(function) > >)
       // requires std::convertible_to< typename impl::VectorTraits< decltype(args) >::value_type, MultiReturnType< decltype(function) > > ||
       //          std::convertible_to< typename decltype(args)::value_type, MultiReturnType< decltype(function) > >
    {
        std::vector< MultiReturnType< decltype(function) > > argvector(args.begin(), args.end());

        // Create a lambda function for a multi variable function,
        // that takes a single argument and returns the function value
        auto f = [&](double value) {
            argvector[index] = value;
            return function(argvector);
        };

        // Return the partial derivative
        return diff< ALGO >(f, argvector[index], stepsize);
    }

    /**
     * @brief Computes the partial derivative of a multivariate function at a specified index.
     *
     * This function computes the partial derivative of a function with multiple arguments,
     * i.e., a function of the form f(x1, x2, ..., xn). The function `partialdiff` uses
     * the provided algorithm template parameter to compute the derivative.
     *
     * @tparam ALGOLOWER Algorithm to be used for differentiation. This template parameter can be
     * any type that satisfies the requirements for a differentiation algorithm. This is used
     * to compute the partial derivative when the point of evaluation is near the lower boundary.
     *
     * @tparam ALGOCENTER Algorithm to be used for differentiation. This template parameter can be
     * any type that satisfies the requirements for a differentiation algorithm.
     *
     * @tparam ALGOUPPER Algorithm to be used for differentiation. This template parameter can be
     * any type that satisfies the requirements for a differentiation algorithm. This is used
     * to compute the partial derivative when the point of evaluation is near the upper boundary.
     *
     * @param function Multivariate function to differentiate. This is a callable object that
     * takes a vector of doubles as argument and returns a floating-point number.
     *
     * @param args Arguments to the function. This parameter can be of any type that can be
     * converted into a vector of floating-point numbers.
     *
     * @param index Index at which the partial derivative should be computed. This parameter
     * specifies the variable with respect to which the derivative is taken.
     *
     * @param limits Pair of lower and upper limits for the variable at the specified index.
     *
     * @param stepsize (Optional) The step size for the derivative. This parameter determines
     * the accuracy of the derivative computation. By default, it uses a function `StepSize`
     * to compute an appropriate step size.
     *
     * @return The partial derivative of the function with respect to the variable at the specified index.
     *
     * @note The returned derivative is evaluated at the provided argument values.
     */
    template< typename ALGOLOWER   = Order1ForwardRichardson,
              typename ALGOCENTRAL = Order1CentralRichardson,
              typename ALGOUPPER   = Order1BackwardRichardson >
    inline auto partialdiff(IsMultiFunction auto                                                                      function,
                            auto                                                                                      args,
                            size_t                                                                                    index,
                            std::pair< MultiReturnType< decltype(function) >, MultiReturnType< decltype(function) > > limits,
                            MultiReturnType< decltype(function) > stepsize = StepSize< MultiReturnType< decltype(function) > >)
       // requires std::convertible_to< typename impl::VectorTraits< decltype(args) >::value_type, MultiReturnType< decltype(function) > > ||
       //          std::convertible_to< typename decltype(args)::value_type, MultiReturnType< decltype(function) > >
    {
        std::vector< MultiReturnType< decltype(function) > > argvector(args.begin(), args.end());

        // Create a lambda function for a multi variable function,
        // that takes a single argument and returns the function value
        auto f = [&](double value) {
            argvector[index] = value;
            return function(argvector);
        };

        if (argvector[index] < limits.first) throw DerivativeError { "Partial derivative undefined. Value is below lower limit." };
        if (argvector[index] > limits.second) throw DerivativeError { "Partial derivative undefined. Value is above upper limit." };

        // Return the partial derivative
        if ((argvector[index] - stepsize) < limits.first) return diff< ALGOLOWER >(f, argvector[index], stepsize);
        if ((argvector[index] + stepsize) > limits.second) return diff< ALGOUPPER >(f, argvector[index], stepsize);

        return diff< ALGOCENTRAL >(f, argvector[index], stepsize);
    }

    /**
     * @brief Computes the derivative of a multivariate function.
     *
     * This function computes the derivative of a function with multiple arguments, i.e.,
     * a function of the form f(x1, x2, ..., xn). The function `multidiff` uses the provided
     * algorithm template parameter to compute the derivative.
     *
     * @tparam ALGO Algorithm to be used for differentiation. This template parameter can be
     * any type that satisfies the requirements for a differentiation algorithm.
     *
     * @param function Multivariate function to differentiate. This is a callable object that
     * takes a vector of doubles as argument and returns a floating-point number.
     *
     * @param args Arguments to the function. This parameter can be of any type that can be
     * converted into a vector of floating-point numbers.
     *
     * @param stepsize (Optional) The step size for the derivative. This parameter determines
     * the accuracy of the derivative computation. By default, it uses a function `StepSize`
     * to compute an appropriate step size.
     *
     * @return A vector containing the derivatives of the function with respect to each argument.
     *
     * @note The returned derivatives are evaluated at the provided argument values.
     */
    template< typename ALGO = Order1CentralRichardson >
    inline auto multidiff(IsMultiFunction auto                  function,
                          auto                                  args,
                          MultiReturnType< decltype(function) > stepsize = StepSize< MultiReturnType< decltype(function) > >)
        //requires std::convertible_to< typename impl::VectorTraits< decltype(args) >::value_type, MultiReturnType< decltype(function) > > ||
        //         std::convertible_to< typename decltype(args)::value_type, MultiReturnType< decltype(function) > >
    {
        std::vector< MultiReturnType< decltype(function) > > argvector(args.begin(), args.end());

        size_t index = 0;
        for (auto arg : argvector) {
            args[index++] = *partialdiff< ALGO >(function, argvector, index, stepsize);
        }

        return args;
    }

    /**
     * @brief Computes the derivative of a multivariate function.
     *
     * This function computes the derivative of a function with multiple arguments, i.e.,
     * a function of the form f(x1, x2, ..., xn). The function `multidiff` uses the provided
     * algorithm template parameter to compute the derivative.
     *
     * @tparam ALGOLOWER Algorithm to be used for differentiation. This template parameter can be
     * any type that satisfies the requirements for a differentiation algorithm. This algorithm
     * is used to compute the partial derivative when the point of evaluation is near the lower
     * limit.
     *
     * @tparam ALGOCENTRAL Algorithm to be used for differentiation. This template parameter can be
     * any type that satisfies the requirements for a differentiation algorithm.
     *
     * @tparam ALGOUPPER Algorithm to be used for differentiation. This template parameter can be
     * any type that satisfies the requirements for a differentiation algorithm. This algorithm
     * is used to compute the partial derivative when the point of evaluation is near the upper
     * limit.
     *
     * @param function Multivariate function to differentiate. This is a callable object that
     * takes a vector of doubles as argument and returns a floating-point number.
     *
     * @param args Arguments to the function. This parameter can be of any type that can be
     * converted into a vector of floating-point numbers.
     *
     * @param limits An array of pairs of values that specify the lower and upper limits for each argument.
     *
     * @param stepsize (Optional) The step size for the derivative. This parameter determines
     * the accuracy of the derivative computation. By default, it uses a function `StepSize`
     * to compute an appropriate step size.
     *
     * @return A vector containing the derivatives of the function with respect to each argument.
     *
     * @note The returned derivatives are evaluated at the provided argument values.
     */
    template< typename ALGOLOWER   = Order1ForwardRichardson,
              typename ALGOCENTRAL = Order1CentralRichardson,
              typename ALGOUPPER   = Order1BackwardRichardson >
    inline auto multidiff(IsMultiFunction auto                                                                                     function,
                          auto                                                                                                     args,
                          std::vector< std::pair< MultiReturnType< decltype(function) >, MultiReturnType< decltype(function) > > > limits,
                          MultiReturnType< decltype(function) > stepsize = StepSize< MultiReturnType< decltype(function) > >)
        requires std::convertible_to< typename impl::VectorTraits< decltype(args) >::value_type, MultiReturnType< decltype(function) > > ||
                 std::convertible_to< typename decltype(args)::value_type, MultiReturnType< decltype(function) > >
    {
        if (args.size() != limits.size()) throw DerivativeError { "Number of limits does not match number of arguments." };

        std::vector< MultiReturnType< decltype(function) > > argvector(args.begin(), args.end());

        size_t index = 0;
        for (auto arg : argvector) {
            args[index++] = *partialdiff< ALGOLOWER, ALGOCENTRAL, ALGOUPPER >(function, argvector, index, limits[index], stepsize);
        }

        return args;
    }

    /**
     * @brief Computes the Jacobian matrix for an array of multivariate functions.
     *
     * The Jacobian matrix is a matrix that represents the first-order partial derivatives
     * of a vector-valued function. This function `jacobian` computes the Jacobian matrix
     * for an array of functions, each of which takes a vector of doubles as argument and
     * returns a floating-point number.
     *
     * @tparam ALGO Algorithm to be used for differentiation. This template parameter can be
     * any type that satisfies the requirements for a differentiation algorithm.
     *
     * @param functions Array of multivariate functions to compute the Jacobian for.
     * Each function in the array is a callable object that takes a vector of doubles
     * as argument and returns a floating-point number.
     *
     * @param args Arguments to the functions. This parameter can be of any type that can be
     * converted into a vector of floating-point numbers.
     *
     * @param stepsize (Optional) The step size for the derivative. This parameter determines
     * the accuracy of the Jacobian computation. By default, it uses a function `StepSize`
     * to compute an appropriate step size.
     *
     * @return The Jacobian matrix.
     *
     * @note The returned Jacobian matrix is evaluated at the provided argument values.
     */
    template< typename ALGO = Order1CentralRichardson >
    inline auto jacobian(IsMultiFunctionArray auto                                   functions,
                         auto                                                        args,
                         MultiReturnType< typename decltype(functions)::value_type > stepsize =
                             StepSize< MultiReturnType< typename decltype(functions)::value_type > >)
        requires std::convertible_to< typename impl::VectorTraits< decltype(args) >::value_type,
                                      MultiReturnType< typename decltype(functions)::value_type > > ||
                 std::convertible_to< typename decltype(args)::value_type, MultiReturnType< typename decltype(functions)::value_type > >
    {
        using namespace nxx::deriv;
        using namespace blaze;

        using return_type = MultiReturnType< typename decltype(functions)::value_type >;

        blaze::DynamicMatrix< return_type > J(functions.size(), functions.size());

        size_t index = 0;
        for (auto& fn : functions) {
            auto pdiffs = multidiff< ALGO >(fn, args, stepsize);
            auto row    = blaze::row(J, index);
            auto dst    = row.begin();
            for (auto src = pdiffs.begin(); src != pdiffs.end(); ++src, ++dst) *dst = *src;
            ++index;
        }

        return J;
    }


    /**
     * @brief Computes the Jacobian matrix for an array of multivariate functions.
     *
     * The Jacobian matrix is a matrix that represents the first-order partial derivatives
     * of a vector-valued function. This function `jacobian` computes the Jacobian matrix
     * for an array of functions, each of which takes a vector of doubles as argument and
     * returns a floating-point number.
     *
     * @tparam ALGOLOWER Algorithm to be used for differentiation. This template parameter can be
     * any type that satisfies the requirements for a differentiation algorithm. This algorithm
     * is used to compute the partial derivative when the point of evaluation is near the lower
     * limit.
     *
     * @tparam ALGOCENTRAL Algorithm to be used for differentiation. This template parameter can be
     * any type that satisfies the requirements for a differentiation algorithm.
     *
     * @tparam ALGOUPPER Algorithm to be used for differentiation. This template parameter can be
     * any type that satisfies the requirements for a differentiation algorithm. This algorithm
     * is used to compute the partial derivative when the point of evaluation is near the upper
     * limit.
     *
     * @param functions Array of multivariate functions to compute the Jacobian for.
     * Each function in the array is a callable object that takes a vector of doubles
     * as argument and returns a floating-point number.
     *
     * @param args Arguments to the functions. This parameter can be of any type that can be
     * converted into a vector of floating-point numbers.
     *
     * @param limits Array of pairs of lower and upper limits for each argument.
     *
     * @param stepsize (Optional) The step size for the derivative. This parameter determines
     * the accuracy of the Jacobian computation. By default, it uses a function `StepSize`
     * to compute an appropriate step size.
     *
     * @return The Jacobian matrix.
     *
     * @note The returned Jacobian matrix is evaluated at the provided argument values.
     */
    template< typename ALGOLOWER   = Order1ForwardRichardson,
              typename ALGOCENTRAL = Order1CentralRichardson,
              typename ALGOUPPER   = Order1BackwardRichardson >
    inline auto jacobian(IsMultiFunctionArray auto                                                               functions,
                         auto                                                                                    args,
                         std::vector< std::pair< MultiReturnType< typename decltype(functions)::value_type >,
                                                 MultiReturnType< typename decltype(functions)::value_type > > > limits,

                         MultiReturnType< typename decltype(functions)::value_type > stepsize =
                             StepSize< MultiReturnType< typename decltype(functions)::value_type > >)
        requires std::convertible_to< typename impl::VectorTraits< decltype(args) >::value_type,
                                      MultiReturnType< typename decltype(functions)::value_type > > ||
                 std::convertible_to< typename decltype(args)::value_type, MultiReturnType< typename decltype(functions)::value_type > >
    {
        using namespace nxx::deriv;
        using namespace blaze;

        if (limits.size() != functions.size()) throw DerivativeError { "The number of limits must match the number of functions." };

        using return_type = MultiReturnType< typename decltype(functions)::value_type >;

        blaze::DynamicMatrix< return_type > J(functions.size(), functions.size());

        size_t index = 0;
        for (auto& fn : functions) {
            auto pdiffs = multidiff< ALGOLOWER, ALGOCENTRAL, ALGOUPPER >(fn, args, limits, stepsize);
            auto row    = blaze::row(J, index);
            auto dst    = row.begin();
            for (auto src = pdiffs.begin(); src != pdiffs.end(); ++src, ++dst) *dst = *src;
            ++index;
        }

        return J;
    }

}    // namespace nxx::deriv

#endif    // NUMERIXX_MULTIDERIVATIVES_HPP

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

#include "calculus/Derivatives.hpp"

#include "blaze/Blaze.h"

namespace nxx::deriv
{

    namespace impl
    {
        /*
         * Forward declaration of the VectorTraits class.
         */
        template< typename... >
        struct VectorTraits;

        /*
         * Specialization of the VectorTraits class for DynamicVector<T>
         */
        template< typename T >
        struct VectorTraits< blaze::DynamicVector< T > >
        {
            using value_type = T;
        };

        /*
         * Specialization of the VectorTraits class for STL containers
         */
        template< typename T >
        struct VectorTraits< T >
        {
            using value_type = typename T::value_type;
        };

        /*
         * Forward declaration of the MatrixTraits class.
         */
        template< typename... >
        struct MatrixTraits;

        /*
         * Specialization of the MatrixTraits class for DynamicMatrix<T>
         */
        template< typename T >
        struct MatrixTraits< blaze::DynamicMatrix< T > >
        {
            using value_type = T;
        };

    }    // namespace impl

    /**
     * @brief Represents the concept of a function that can be applied to a vector of doubles and returns a floating-point value.
     *
     * The IsMultiFunction concept is satisfied by any function that meets the following criteria:
     * - Accepts a single argument of type `FN`, which represents the function to be checked.
     * - The function `FN` can be invoked with a vector of doubles as its argument.
     * - The return type of the function `FN` is a floating-point value.
     *
     * @param fn The function to be checked.
     * @return `true` if the given function satisfies the IsMultiFunction concept, otherwise `false`.
     *
     * @note This concept uses the `std::invoke_result_t` template to check the return type of the function.
     * @note The `std::floating_point` template is used to check that the return type is a floating-point value.
     * @note The `std::vector` template is used to represent the argument type of the function.
     */
    template< typename FN >
    concept IsMultiFunction = requires(FN fn) {
                                  // clang-format off
                                  { std::floating_point< std::invoke_result_t< FN, std::vector< double > > >};
                                  // clang-format on
                              };

    /**
     * @brief Concept specifying requirements for a multi-function array.
     *
     * This concept defines the requirements for a type `ARR` to be considered a
     * multi-function array. A multi-function array is an array-like container
     * where each element can be evaluated with a function and meets the
     * requirements of the `IsMultiFunction` concept.
     *
     * @tparam ARR The type of the array to be checked.
     */
    template< typename ARR >
    concept IsMultiFunctionArray = requires(ARR arr) {
                                       // clang-format off
                                       { arr.begin() } -> std::input_iterator;
                                       { arr.end() } -> std::input_iterator;
                                       { IsMultiFunction< typename ARR::value_type> };
                                       // clang-format on
                                   };

    template< typename T >
    using MultiReturnType = std::invoke_result_t< T, std::vector< double > >;

    /**
     * @brief Calculate the partial derivative of a multi-function.
     *
     * This function calculates the partial derivative of a multi-function with respect to the
     * argument at the specified index. The user can optionally provide a step size for the
     * numerical approximation of the derivative.
     *
     * @tparam IsMultiFunction   Boolean flag indicating whether the provided function is a multi-function.
     * @param function          The multi-function for which the partial derivative needs to be calculated.
     * @param args              The arguments at which the function is evaluated.
     * @param index             The index of the argument with respect to which the derivative is computed.
     * @param stepsize          Optional step size for numerical approximation. Default is the
     *                          step size obtained from StepSize<MultiReturnType<decltype(function)>>.
     *
     * @return The partial derivative of the multi-function.
     *
     * @note This function requires the argument type to be convertible to the return type of the function.
     * @see StepSize
     */
    template< typename ALGO >
    inline auto partialdiff(IsMultiFunction auto                  function,
                            auto                                  args,
                            size_t                                index,
                            MultiReturnType< decltype(function) > stepsize = StepSize< MultiReturnType< decltype(function) > >)
        requires std::convertible_to< typename impl::VectorTraits< decltype(args) >::value_type, MultiReturnType< decltype(function) > >
    {
        std::vector< MultiReturnType< decltype(function) > > argvector(args.begin(), args.end());

        // Create a lambda function for a multi variable function,
        // that takes a single argument and returns the function value
        auto f = [&](double value) {
            argvector[index] = value;
            return function(argvector);
        };

        // Return the partial derivative
        return diff< ALGO >(f, argvector[index], std::max(stepsize, stepsize * argvector[index]));
    }

    /**
     * @brief Calculates the multidimensional difference for the given function and arguments.
     *
     * The multidiff function calculates the multidimensional difference for the given function and arguments.
     * It takes in three parameters:
     *     - function: The function for which the difference is to be calculated.
     *     - args: The arguments for the function.
     *     - stepsize: The step size to be used for the difference calculation (default: StepSize<MultiReturnType<decltype(function)>>).
     *
     * The function returns the multidimensional difference.
     *
     * @param function The function for which the difference is to be calculated.
     * @param args The arguments for the function.
     * @param stepsize The step size to be used for the difference calculation (default: StepSize<MultiReturnType<decltype(function)>>).
     * @return The multidimensional difference.
     *
     * @tparam function The type of function for which the difference is to be calculated.
     * @tparam args The type of arguments for the function.
     * @tparam stepsize The type of step size to be used for the difference calculation.
     *                 It must be convertible to MultiReturnType<decltype(function)>.
     *
     * @requires std::convertible_to<typename impl::VectorTraits<decltype(args)>::value_type, MultiReturnType<decltype(function)>>.
     */
    template< typename ALGO >
    inline auto multidiff(IsMultiFunction auto                  function,
                          auto                                  args,
                          MultiReturnType< decltype(function) > stepsize = StepSize< MultiReturnType< decltype(function) > >)
        requires std::convertible_to< typename impl::VectorTraits< decltype(args) >::value_type, MultiReturnType< decltype(function) > >
    {
        std::vector< MultiReturnType< decltype(function) > > argvector(args.begin(), args.end());

        size_t index = 0;
        for (auto arg : argvector) {
            args[index++] = *partialdiff< ALGO >(function, argvector, index, std::max(stepsize, stepsize * arg));
        }

        return args;
    }

    /**
     * \brief Calculate the Jacobian matrix for a given set of functions and arguments.
     *
     * This function calculates the Jacobian matrix for a given set of functions and arguments.
     * The Jacobian matrix represents the partial derivatives of each function with respect to each argument.
     *
     * \tparam IsMultiFunctionArray A type that represents an array or range of functions.
     * \tparam args A type that represents the arguments for the functions.
     * \param functions An array or range of functions.
     * \param args The arguments for the functions.
     * \param stepsize The step size used for numerical differentiation (optional).
     *                If not provided, the default step size based on the return type of the functions will be used.
     *
     * \return The Jacobian matrix as a matrix-like data structure.
     *
     * \requirements
     * - The `functions` must be an array or range of functions.
     * - The `args` must be convertible to the return type of the functions.
     *
     * \note
     * - The `functions` and `args` must have the same number of elements, otherwise undefined behavior may occur.
     * - The `impl::VectorTraits` type provides information about the vector-like properties of the `args` type.
     *
     * \see StepSize
     * \see MultiReturnType
     * \see impl::VectorTraits
     * \see JacobianMatrix
     */
    template< typename ALGO >
    inline auto jacobian(IsMultiFunctionArray auto                                   functions,
                         auto                                                        args,
                         MultiReturnType< typename decltype(functions)::value_type > stepsize = StepSize< MultiReturnType< typename decltype(functions)::value_type > >)
        requires std::convertible_to< typename impl::VectorTraits< decltype(args) >::value_type, MultiReturnType< typename decltype(functions)::value_type > >
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

}    // namespace nxx::deriv

#endif    // NUMERIXX_MULTIDERIVATIVES_HPP

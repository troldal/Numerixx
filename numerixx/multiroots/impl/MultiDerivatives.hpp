/*
    888b      88  88        88  88b           d88  88888888888  88888888ba   88  8b        d8  8b        d8
    8888b     88  88        88  888b         d888  88           88      "8b  88   Y8,    ,8P    Y8,    ,8P
    88 `8b    88  88        88  88`8b       d8'88  88           88      ,8P  88    `8b  d8'      `8b  d8'
    88  `8b   88  88        88  88 `8b     d8' 88  88aaaaa      88aaaaaa8P'  88      Y88P          Y88P
    88   `8b  88  88        88  88  `8b   d8'  88  88"""""      88""""88'    88      d88b          d88b
    88    `8b 88  88        88  88   `8b d8'   88  88           88    `8b    88    ,8P  Y8,      ,8P  Y8,
    88     `8888  Y8a.    .a8P  88    `888'    88  88           88     `8b   88   d8'    `8b    d8'    `8b
    88      `888   `"Y8888Y"'   88     `8'     88  88888888888  88      `8b  88  8P        Y8  8P        Y8

    Copyright © 2022 Kenneth Troldal Balslev

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

// ===== Numerixx Includes
#include <Deriv.hpp>
#include "MultiFunctionArray.hpp"

// ===== External Includes
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


    template<typename ALGO = Order1CentralRichardson, typename T, typename Func>
    std::vector<T> partialdiff(const Func& func, const std::vector<T>& point) {
        std::vector<T> derivatives(point.size(), T(0));

        for (size_t i = 0; i < point.size(); ++i) {
            // Create a single-variable function by fixing all variables except the i-th variable
            auto singleVarFunc = [&, i](T x) {
                std::vector<T> tempPoint = point;
                tempPoint[i] = x;
                return func(tempPoint);
            };

            // Compute the derivative using the diff function
            derivatives[i] = *diff<ALGO>(singleVarFunc, point[i]);
        }

        return derivatives;
    }

    template<typename T>
    blaze::DynamicMatrix<T> jacobian(const multiroots::DynamicFunctionArray<T>& functions, const std::vector<T>& point) {
        // Determine the size of the Jacobian matrix
        size_t numRows = std::distance(functions.begin(), functions.end());
        size_t numCols = point.size();

        // Create a matrix to hold the Jacobian
        blaze::DynamicMatrix<T> J(numRows, numCols);

        // Compute the partial derivatives for each function
        size_t row = 0;
        for (const auto& func : functions) {
            std::vector<T> partials = partialdiff(func, point);
            for (size_t col = 0; col < numCols; ++col) {
                J(row, col) = partials[col];
            }
            ++row;
        }

        return J;
    }

    template< typename T >
    blaze::DynamicMatrix< T > jacobian(const multiroots::DynamicFunctionArray< T >& functions,
                                              const blaze::DynamicVector< T >&             point)
    {
        return jacobian(functions, std::vector< T >(point.begin(), point.end()));
    }

    template< typename T >
    blaze::DynamicMatrix< T > jacobian(const multiroots::DynamicFunctionArray< T >& functions,
                                              const std::initializer_list< T >&            point)
    {
        return jacobian(functions, std::vector< T >(point));
    }

}    // namespace nxx::deriv

#endif    // NUMERIXX_MULTIDERIVATIVES_HPP

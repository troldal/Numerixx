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
#include "ContainerTraits.hpp"
#include "MultiFunctionArray.hpp"
#include <Deriv.hpp>

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

    /**
     * @brief Partial derivative computation implementation.
     *
     * @tparam ALGO The algorithm used for computing derivatives.
     * @tparam RET_T The return type of the function whose derivative is being computed.
     * @param func A function object representing the function to differentiate.
     * @param point A container representing the point at which the derivative is computed.
     * @param derivatives A container to store the computed derivatives.
     *
     * @details
     * Computes the partial derivatives of a multi-variable function at a given point.
     * The function `func` is differentiated with respect to each variable in `point` separately.
     * The derivative computation is performed using the specified algorithm `ALGO`.
     *
     * @note Allocating `singleVarFunc` inside the loop is not optimal and can be improved.
     *
     * @throws std::runtime_error If the algorithm fails to compute the derivative.
     */
    template< typename ALGO, typename RET_T >
    void partialdiff_impl(const auto& func, const auto& point, auto& derivatives)
    {
        using CONT_T = std::remove_cvref_t< decltype(point) >;
        using ARG_T  = traits::ContainerValueType_t< CONT_T >;
        static_assert(sizeof(ARG_T) <= sizeof(RET_T), "The precision of the argument types exceeds that of the return type.");

        for (size_t i = 0; i < point.size(); ++i) {
            // TODO: allocating inside the loop is not optimal
            auto singleVarFunc = [&, i](ARG_T x) {
                auto tempPoint = point;
                tempPoint[i]   = x;
                return func(std::span< ARG_T >(tempPoint.data(), tempPoint.size()));
            };

            // Compute the derivative using the diff function
            derivatives[i] = *diff< ALGO >(singleVarFunc, point[i]);
        }
    }

    /**
     * @brief Computes partial derivatives of a multi-variable function.
     *
     * @tparam ALGO The algorithm used for computing derivatives. Defaults to Order1CentralRichardson.
     * @tparam RET_T The return type of the function whose derivative is being computed.
     * @tparam PARAM_T The parameter type of the function.
     * @tparam CONT_T The container type for the function arguments and derivatives.
     * @param func A `multiroots::MultiFunction` object representing the function to differentiate.
     * @param point A container representing the point at which the derivative is computed.
     * @return A container of the same type as `point` containing the computed derivatives.
     *
     * @details
     * This function template provides a convenient interface for computing partial derivatives of a
     * multi-variable function represented by a `multiroots::MultiFunction` object. The derivatives
     * are computed at the specified `point`.
     *
     * @note Requires `ALGO` to be a valid differential solver and the elements of `CONT_T` to be floating-point types.
     * @note The CONT_T container must be a contiguous container, and support the interface of the standard containers.
     */
    template< typename ALGO = Order1CentralRichardson, typename RET_T, typename PARAM_T, typename CONT_T >
    requires ALGO::IsDiffSolver && nxx::IsFloat< traits::ContainerValueType_t< CONT_T > >
    auto partialdiff(const multiroots::MultiFunction< RET_T, PARAM_T >& func, const CONT_T& point)
    {
        using ARG_T = traits::ContainerValueType_t< CONT_T >;
        CONT_T derivatives(point.size(), ARG_T {});

        partialdiff_impl< ALGO, RET_T >(func, point, derivatives);

        return derivatives;
    }

    /**
     * @brief Computes partial derivatives of a multi-variable function, outputting results in a generalized container type.
     *
     * @tparam OUT_T The template class of the output container for storing derivatives, with a flexible template parameter list.
     * @tparam ALGO The algorithm used for computing derivatives. Defaults to Order1CentralRichardson.
     * @tparam RET_T The return type of the function whose derivative is being computed.
     * @tparam PARAM_T The parameter type of the function.
     * @tparam CONT_T The container type for the function arguments.
     * @param func A `multiroots::MultiFunction` object representing the function to differentiate.
     * @param point A container representing the point at which the derivative is computed.
     * @return An `OUT_T` container containing the computed derivatives.
     *
     * @details
     * This template function computes the partial derivatives of a multi-variable function at a specified point,
     * storing the results in a container of type `OUT_T`. The `OUT_T` template is designed to accommodate a variety
     * of container types with different template parameter lists, offering greater flexibility in the choice of output
     * container for the derivatives.
     *
     * The function `func` is differentiated at each variable in `point` separately, using the specified algorithm `ALGO`.
     * The derivatives are stored in an `OUT_T<ARG_T>` container, where `ARG_T` is the value type of the input container `CONT_T`.
     *
     * @note Requires `ALGO` to be a valid differential solver and the elements of `CONT_T` to be floating-point types.
     * @note The CONT_T container must be a contiguous container, and support the interface of the standard containers.
     * @note The OUT_T container must be a random access container, and support the interface of the standard containers.
     *
     * @throws std::runtime_error If the derivative computation algorithm fails.
     */
    template< template< typename... > class OUT_T,
              typename ALGO = Order1CentralRichardson,
              typename RET_T,
              typename PARAM_T,
              typename CONT_T >
    requires ALGO::IsDiffSolver && nxx::IsFloat< traits::ContainerValueType_t< CONT_T > >
    auto partialdiff(const multiroots::MultiFunction< RET_T, PARAM_T >& func, const CONT_T& point)
    {
        using ARG_T = traits::ContainerValueType_t< CONT_T >;
        OUT_T< ARG_T > derivatives(point.size(), ARG_T {});

        partialdiff_impl< ALGO, RET_T >(func, point, derivatives);

        return derivatives;
    }

    /**
     * @brief Computes partial derivatives of a multi-variable function, outputting results in a specified container type.
     *
     * @tparam OUT_T The template class of the output container for storing derivatives.
     * @tparam ALGO The algorithm used for computing derivatives. Defaults to Order1CentralRichardson.
     * @tparam RET_T The return type of the function whose derivative is being computed.
     * @tparam PARAM_T The parameter type of the function.
     * @tparam CONT_T The container type for the function arguments.
     * @param func A `multiroots::MultiFunction` object representing the function to differentiate.
     * @param point A container representing the point at which the derivative is computed.
     * @return An `OUT_T` container containing the computed derivatives.
     *
     * @details
     * This template function computes the partial derivatives of a given multi-variable function
     * at a specified point. It stores the computed derivatives in an `OUT_T` type container.
     * This allows for flexibility in the choice of container type for the derivatives.
     *
     * The function `func` is represented by a `multiroots::MultiFunction` object and is differentiated
     * at each variable in `point` separately. The derivative computation is carried out using the algorithm
     * specified by `ALGO`.
     *
     * @note Requires `ALGO` to be a valid differential solver and the elements of `CONT_T` to be floating-point types.
     * @note The CONT_T container must be a contiguous container, and support the interface of the standard containers.
     * @note This particular overload has the purpose of supporting blaze::DynamicVector as the output container.
     *
     * @throws std::runtime_error If the algorithm fails to compute the derivative.
     */
    template< template< typename, bool, typename... > class OUT_T,
              typename ALGO = Order1CentralRichardson,
              typename RET_T,
              typename PARAM_T,
              typename CONT_T >
    requires ALGO::IsDiffSolver && nxx::IsFloat< traits::ContainerValueType_t< CONT_T > >
    auto partialdiff(const multiroots::MultiFunction< RET_T, PARAM_T >& func, const CONT_T& point)
    {
        using ARG_T = traits::ContainerValueType_t< CONT_T >;
        OUT_T< ARG_T, false > derivatives(point.size(), ARG_T {});

        partialdiff_impl< ALGO, RET_T >(func, point, derivatives);

        return derivatives;
    }

    /**
     * @brief Computes partial derivatives of a multi-variable function using an initializer list for the point of differentiation.
     *
     * @tparam ALGO The algorithm used for computing derivatives. Defaults to Order1CentralRichardson.
     * @tparam RET_T The return type of the function whose derivative is being computed.
     * @tparam PARAM_T The parameter type of the function.
     * @tparam ARG_T The type of the elements in the initializer list.
     * @param func A `multiroots::MultiFunction` object representing the function to differentiate.
     * @param point An `std::initializer_list` representing the point at which the derivative is computed.
     * @return A container (determined by the default behavior of `partialdiff`) containing the computed derivatives.
     *
     * @details
     * This template function provides a convenient way to compute the partial derivatives of a multi-variable function
     * at a specified point, using an `std::initializer_list` to define the point. This allows for a more concise and
     * readable syntax when the dimensions of the point are fixed and known at compile time.
     *
     * The function internally converts the initializer list into a `std::vector` and then delegates the computation
     * to the `partialdiff` function template. The type of the output container for the derivatives is determined by
     * the default behavior of the called `partialdiff` overload.
     *
     * @note Requires `ALGO` to be a valid differential solver.
     *
     * @throws std::runtime_error If the derivative computation algorithm fails.
     */
    template< typename ALGO = Order1CentralRichardson, typename RET_T, typename PARAM_T, typename ARG_T >
    requires ALGO::IsDiffSolver
    auto partialdiff(const multiroots::MultiFunction< RET_T, PARAM_T >& func, const std::initializer_list< ARG_T >& point)
    {
        return partialdiff< ALGO >(func, std::vector< ARG_T >(point));
    }

    /**
     * @brief Computes partial derivatives of a multi-variable function using an initializer list, with a specified output container type.
     *
     * @tparam OUT_T The template class of the output container for storing derivatives.
     * @tparam ALGO The algorithm used for computing derivatives. Defaults to Order1CentralRichardson.
     * @tparam RET_T The return type of the function whose derivative is being computed.
     * @tparam PARAM_T The parameter type of the function.
     * @tparam ARG_T The type of the elements in the initializer list.
     * @param func A `multiroots::MultiFunction` object representing the function to differentiate.
     * @param point An `std::initializer_list` representing the point at which the derivative is computed.
     * @return An `OUT_T` container containing the computed derivatives.
     *
     * @details
     * This template function extends the functionality of `partialdiff` to allow the use of an initializer list for specifying
     * the point of differentiation. This offers a concise way to define points, especially useful in cases with known
     * dimensions at compile time. The computed derivatives are stored in a container specified by the `OUT_T` template parameter.
     *
     * Internally, the function converts the initializer list into a `std::vector`, then delegates the computation to another
     * `partialdiff` overload that handles vector inputs. This overload allows the user to specify the type of the output container
     * for the derivatives, adding flexibility in how the results are stored and used.
     *
     * @note Requires `ALGO` to be a valid differential solver.
     * @note The OUT_T container must be a random access container, and support the interface of the standard containers.
     *
     * @throws std::runtime_error If the derivative computation algorithm fails.
     */
    template< template< typename... > class OUT_T,
              typename ALGO = Order1CentralRichardson,
              typename RET_T,
              typename PARAM_T,
              typename ARG_T >
    requires ALGO::IsDiffSolver
    auto partialdiff(const multiroots::MultiFunction< RET_T, PARAM_T >& func, const std::initializer_list< ARG_T >& point)
    {
        return partialdiff< OUT_T, ALGO >(func, std::vector< ARG_T >(point));
    }

    /**
     * @brief Computes partial derivatives of a multi-variable function using an initializer list, with a specified output container type.
     *
     * @tparam OUT_T The template class of the output container for storing derivatives, with a specific template parameter format.
     * @tparam ALGO The algorithm used for computing derivatives. Defaults to Order1CentralRichardson.
     * @tparam RET_T The return type of the function whose derivative is being computed.
     * @tparam PARAM_T The parameter type of the function.
     * @tparam ARG_T The type of the elements in the initializer list.
     * @param func A `multiroots::MultiFunction` object representing the function to differentiate.
     * @param point An `std::initializer_list` representing the point at which the derivative is computed.
     * @return An `OUT_T` container containing the computed derivatives.
     *
     * @details
     * This template function allows for the computation of partial derivatives of a multi-variable function at a point
     * specified by an initializer list. It provides the flexibility to define the output container type through the `OUT_T`
     * template parameter. This variant of `OUT_T` expects a specific format for its template parameters, enabling more control
     * over the type of the output container.
     *
     * The function internally converts the initializer list to a `std::vector` of `ARG_T` type and then utilizes the
     * `partialdiff` overload that takes a vector as input. The use of an initializer list simplifies the syntax for defining
     * the point of differentiation, particularly when its dimensions are known at compile time.
     *
     * @note Requires `ALGO` to be a valid differential solver.
     * @note This particular overload has the purpose of supporting blaze::DynamicVector as the output container.
     *
     * @throws std::runtime_error If the derivative computation algorithm fails.
     */
    template< template< typename, bool, typename... > class OUT_T,
              typename ALGO = Order1CentralRichardson,
              typename RET_T,
              typename PARAM_T,
              typename ARG_T >
    requires ALGO::IsDiffSolver
    auto partialdiff(const multiroots::MultiFunction< RET_T, PARAM_T >& func, const std::initializer_list< ARG_T >& point)
    {
        return partialdiff< OUT_T, ALGO >(func, std::vector< ARG_T >(point));
    }

    template< typename T >
    blaze::DynamicMatrix< T > jacobian(const multiroots::DynamicFunctionArray< T >& functions, const std::vector< T >& point)
    {
        // Determine the size of the Jacobian matrix
        size_t numRows = std::distance(functions.begin(), functions.end());
        size_t numCols = point.size();

        // Create a matrix to hold the Jacobian
        blaze::DynamicMatrix< T > J(numRows, numCols);

        // Compute the partial derivatives for each function
        size_t row = 0;
        for (const auto& func : functions) {
            std::vector< T > partials = partialdiff(func, point);
            for (size_t col = 0; col < numCols; ++col) {
                J(row, col) = partials[col];
            }
            ++row;
        }

        return J;
    }

    template< typename T >
    blaze::DynamicMatrix< T > jacobian(const multiroots::DynamicFunctionArray< T >& functions, const blaze::DynamicVector< T >& point)
    {
        return jacobian(functions, std::vector< T >(point.begin(), point.end()));
    }

    template< typename T >
    blaze::DynamicMatrix< T > jacobian(const multiroots::DynamicFunctionArray< T >& functions, const std::initializer_list< T >& point)
    {
        return jacobian(functions, std::vector< T >(point));
    }

}    // namespace nxx::deriv

#endif    // NUMERIXX_MULTIDERIVATIVES_HPP

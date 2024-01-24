/*
    888b      88  88        88  88b           d88  88888888888  88888888ba   88  8b        d8  8b        d8
    8888b     88  88        88  888b         d888  88           88      "8b  88   Y8,    ,8P    Y8,    ,8P
    88 `8b    88  88        88  88`8b       d8'88  88           88      ,8P  88    `8b  d8'      `8b  d8'
    88  `8b   88  88        88  88 `8b     d8' 88  88aaaaa      88aaaaaa8P'  88      Y88P          Y88P
    88   `8b  88  88        88  88  `8b   d8'  88  88"""""      88""""88'    88      d88b          d88b
    88    `8b 88  88        88  88   `8b d8'   88  88           88    `8b    88    ,8P  Y8,      ,8P  Y8,
    88     `8888  Y8a.    .a8P  88    `888'    88  88           88     `8b   88   d8'    `8b    d8'    `8b
    88      `888   `"Y8888Y"'   88     `8'     88  88888888888  88      `8b  88  8P        Y8  8P        Y8

    Copyright © 2024 Kenneth Troldal Balslev

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

/**
 * @file partialdiff.hpp
 * @brief Provides functionality for computing partial derivatives of multi-variable functions,
 * as well as matrices of derivatives such as the Jacobian and Hessian matrices.
 *
 * This file contains templates for calculating the partial derivatives of functions that
 * take multiple arguments. It utilizes the algorithm specified by the `ALGO` template parameter
 * to compute these derivatives.
 */

#pragma once

// ===== Numerixx Includes
#include "ContainerTraits.hpp"
#include "MultiFunctionArray.hpp"
#include <Deriv.hpp>

// ===== External Includes
#include "_external.hpp"

namespace nxx::deriv
{
    namespace detail
    {
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
    }    // namespace detail

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

        detail::partialdiff_impl< ALGO, RET_T >(func, point, derivatives);

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

        detail::partialdiff_impl< ALGO, RET_T >(func, point, derivatives);

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

        detail::partialdiff_impl< ALGO, RET_T >(func, point, derivatives);

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

    /**
     * @brief Computes a matrix of derivatives (such as the Jacobian) for a set of multi-variable functions.
     *
     * @tparam ALGO The algorithm used for computing derivatives.
     * @tparam RES_T The result type of the functions.
     * @tparam PARAM_T The parameter type of the functions.
     * @tparam CONTAINER_T The container type for the function arguments.
     * @param functions A `multiroots::MultiFunctionArray` representing the set of functions.
     * @param point A container representing the point at which the derivatives are computed.
     * @return A `blaze::DynamicMatrix` containing the derivatives, structured as a matrix.
     *
     * @details
     * This template function computes a matrix of derivatives, such as the Jacobian matrix, for a set of multi-variable functions.
     * Each function is represented within a `multiroots::MultiFunctionArray`. The function calculates the partial derivatives of each
     * function at a specified point, arranging these derivatives into a matrix.
     *
     * The number of rows in the resulting matrix corresponds to the number of functions in the array (`functions.size()`),
     * and the number of columns corresponds to the number of variables in each function (`point.size()`). Each row in the matrix
     * represents the gradient of a single function at the specified point.
     *
     * The derivative computation for each function is performed using the specified algorithm `ALGO`. This flexibility allows for
     * the use of different numerical methods to suit the specific requirements of the derivative calculation.
     *
     * @throws std::runtime_error If the derivative computation algorithm fails.
     */
    template< typename ALGO, typename RES_T, typename PARAM_T, typename CONTAINER_T >
    blaze::DynamicMatrix< RES_T > multidiff(const multiroots::MultiFunctionArray< RES_T, PARAM_T >& functions, const CONTAINER_T& point)
    {
        // Determine the size of the Jacobian matrix
        size_t numRows = functions.size();
        size_t numCols = point.size();

        // Create a matrix to hold the Jacobian
        blaze::DynamicMatrix< RES_T > J(numRows, numCols);

        // Compute the partial derivatives for each function
        size_t row = 0;
        for (const auto& func : functions) {
            auto partials = partialdiff< ALGO >(func, point);
            for (size_t col = 0; col < numCols; ++col) {
                J(row, col) = partials[col];
            }
            ++row;
        }

        return J;
    }

    /**
     * @brief Computes a matrix of derivatives for a set of multi-variable functions using an initializer list to specify the point.
     *
     * @tparam ALGO The algorithm used for computing derivatives.
     * @tparam RES_T The result type of the functions.
     * @tparam PARAM_T The parameter type of the functions.
     * @param functions A `multiroots::MultiFunctionArray` representing the set of functions.
     * @param point An `std::initializer_list` representing the point at which the derivatives are computed.
     * @return A `blaze::DynamicMatrix` containing the derivatives, structured as a matrix.
     *
     * @details
     * This template function extends `multidiff` to allow the use of an initializer list for specifying the point
     * of differentiation. It is particularly useful when the dimensions of the point are fixed and known at compile time,
     * allowing for a more concise and readable syntax.
     *
     * The function internally converts the initializer list into a `std::vector` and then delegates the computation
     * to the `multidiff` function template that takes a container as input. This allows for the same flexibility in
     * choosing the derivative computation algorithm (`ALGO`) and supports the same output format as the primary `multidiff` function.
     *
     * @throws std::runtime_error If the derivative computation algorithm fails.
     */
    template< typename ALGO, typename RES_T, typename PARAM_T >
    blaze::DynamicMatrix< RES_T > multidiff(const multiroots::MultiFunctionArray< RES_T, PARAM_T >& functions,
                                            const std::initializer_list< RES_T >&                   point)
    {
        return multidiff< ALGO >(functions, std::vector< RES_T >(point));
    }

    /**
     * @brief Computes the Jacobian matrix for a set of multi-variable functions.
     *
     * @tparam RES_T The result type of the functions.
     * @tparam PARAM_T The parameter type of the functions.
     * @tparam CONTAINER_T The container type for the function arguments.
     * @param functions A `multiroots::MultiFunctionArray` representing the set of functions.
     * @param point A container representing the point at which the Jacobian is computed.
     * @return A `blaze::DynamicMatrix` containing the Jacobian matrix.
     *
     * @details
     * This function calculates the Jacobian matrix of a set of multi-variable functions represented by a
     * `multiroots::MultiFunctionArray`. The Jacobian matrix is a matrix of all first-order partial derivatives
     * of a vector-valued function. Each row of the Jacobian matrix corresponds to the gradient of one function
     * in the array with respect to the variables specified in `point`.
     *
     * The derivative computation utilizes the `Order1CentralRichardson` algorithm, which is suitable for first-order
     * derivative calculations. This algorithm choice provides a balance between computational efficiency and
     * accuracy for most use cases.
     *
     * @throws std::runtime_error If the derivative computation algorithm fails.
     */
    template< typename RES_T, typename PARAM_T, typename CONTAINER_T >
    blaze::DynamicMatrix< RES_T > jacobian(const multiroots::MultiFunctionArray< RES_T, PARAM_T >& functions, const CONTAINER_T& point)
    {
        return multidiff< Order1CentralRichardson >(functions, point);
    }

    /**
     * @brief Computes the Jacobian matrix for a set of multi-variable functions using an initializer list to specify the point.
     *
     * @tparam RES_T The result type of the functions.
     * @tparam PARAM_T The parameter type of the functions.
     * @param functions A `multiroots::MultiFunctionArray` representing the set of functions.
     * @param point An `std::initializer_list` representing the point at which the Jacobian is computed.
     * @return A `blaze::DynamicMatrix` containing the Jacobian matrix.
     *
     * @details
     * This function extends the functionality of `jacobian` to allow the specification of the differentiation point using an
     * initializer list. This is particularly useful for cases where the point's dimensions are fixed and known at compile time,
     * offering a more concise and readable way to define the point.
     *
     * The function internally converts the initializer list into a `std::vector` and then delegates the Jacobian matrix computation
     * to the `multidiff` function template, using the `Order1CentralRichardson` algorithm. This approach maintains the consistency
     * of using a first-order central difference method (Richardson extrapolation) for calculating the Jacobian, balancing accuracy
     * and computational efficiency.
     *
     * @throws std::runtime_error If the derivative computation algorithm fails.
     */
    template< typename RES_T, typename PARAM_T >
    blaze::DynamicMatrix< RES_T > jacobian(const multiroots::MultiFunctionArray< RES_T, PARAM_T >& functions,
                                           const std::initializer_list< RES_T >&                   point)
    {
        return multidiff< Order1CentralRichardson >(functions, std::vector< RES_T >(point));
    }

    /**
     * @brief Computes the Hessian matrix for a set of multi-variable functions.
     *
     * @tparam RES_T The result type of the functions.
     * @tparam PARAM_T The parameter type of the functions.
     * @tparam CONTAINER_T The container type for the function arguments.
     * @param functions A `multiroots::MultiFunctionArray` representing the set of functions.
     * @param point A container representing the point at which the Hessian is computed.
     * @return A `blaze::DynamicMatrix` containing the Hessian matrix.
     *
     * @details
     * This function calculates the Hessian matrix of a set of multi-variable functions represented by a
     * `multiroots::MultiFunctionArray`. The Hessian matrix is a square matrix of second-order partial derivatives
     * of a vector-valued function. It provides critical insights into the curvature of each function in the array
     * at the specified `point`.
     *
     * The derivative computation utilizes the `Order2Central5Point` algorithm, specifically designed for second-order
     * derivative calculations. This algorithm offers enhanced accuracy for computing second-order derivatives, which is
     * crucial for the correct estimation of the Hessian matrix.
     *
     * @throws std::runtime_error If the derivative computation algorithm fails.
     */
    template< typename RES_T, typename PARAM_T, typename CONTAINER_T >
    blaze::DynamicMatrix< RES_T > hessian(const multiroots::MultiFunctionArray< RES_T, PARAM_T >& functions, const CONTAINER_T& point)
    {
        return multidiff< Order2Central5Point >(functions, point);
    }

    /**
     * @brief Computes the Hessian matrix for a set of multi-variable functions using an initializer list to specify the point.
     *
     * @tparam RES_T The result type of the functions.
     * @tparam PARAM_T The parameter type of the functions.
     * @param functions A `multiroots::MultiFunctionArray` representing the set of functions.
     * @param point An `std::initializer_list` representing the point at which the Hessian is computed.
     * @return A `blaze::DynamicMatrix` containing the Hessian matrix.
     *
     * @details
     * This function extends the functionality of `hessian` to allow for the specification of the differentiation point
     * using an initializer list. This approach is particularly useful when the dimensions of the point are fixed and
     * known at compile time, providing a concise and readable way to define the point.
     *
     * Internally, the function converts the initializer list into a `std::vector` and then delegates the Hessian matrix
     * computation to the `multidiff` function template, utilizing the `Order2Central5Point` algorithm. This algorithm
     * is specifically tailored for accurate second-order derivative calculations, which are essential for the correct
     * estimation of the Hessian matrix.
     *
     * @throws std::runtime_error If the derivative computation algorithm fails.
     */
    template< typename RES_T, typename PARAM_T >
    blaze::DynamicMatrix< RES_T > hessian(const multiroots::MultiFunctionArray< RES_T, PARAM_T >& functions,
                                          const std::initializer_list< RES_T >&                   point)
    {
        return multidiff< Order2Central5Point >(functions, std::vector< RES_T >(point));
    }

}    // namespace nxx::deriv


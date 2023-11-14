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

#ifndef NUMERIXX_DIFFERENTIATION_HPP
#define NUMERIXX_DIFFERENTIATION_HPP

// ===== Numerixx Includes
#include <Concepts.hpp>
#include <Error.hpp>

// ===== External Includes
#include <gcem.hpp>
#include <tl/expected.hpp>

// ===== Standard Library Includes
#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>

namespace nxx::deriv
{

    namespace detail
    {
        template< typename T >
        struct DerivErrorData
        {
            T x;
            T h;
            T f;
            T df;

            friend std::ostream& operator<<(std::ostream& os, const DerivErrorData& data)
            {
                os << "x = " << data.x << ", h = " << data.h << ", f = " << data.f << ", df = " << data.df;
                return os;
            }
        };
    }    // namespace detail

    /**
     * @brief Concept checking whether a type is a callable function object that returns a floating point type.
     * @tparam FN The type to check.
     */
    template< typename FN >
    concept IsFunction = requires(FN fn) {
                             {
                                 fn(0.0)
                             } -> std::floating_point;
                         };

    /**
     * @brief Alias template for determining the step size for a given type when computing numerical derivatives.
     * @tparam T The type for which to compute the step size.
     */
    template< typename T >
    requires std::floating_point< T >
    struct StepSizeHelper
    {
        static constexpr T value = gcem::pow(std::numeric_limits< T >::epsilon(), 1.0 / 3.0);
    };

    /**
     * @brief Constant expression for the step size used when computing numerical derivatives.
     * @tparam T The type for which to compute the step size.
     */
    template< typename T >
    requires std::floating_point< T >
    inline constexpr T StepSize = StepSizeHelper< T >::value;

    // ====================================================================
    // Central finite difference formulas
    // ====================================================================

    template< std::floating_point T >
    void validateStepSize(T stepsize, T minStepSize)
    {
        if (stepsize < minStepSize) throw NumerixxError("Step size is too low.");
    }

    /**
     * @brief A class defining a function object for computing the derivative of an arbitrary function.
     *
     * @tparam FN The type of function object to use for computing the derivative.
     */
    template< typename FN >
    class SolverTemplate
    {
    public:
        /**
         * @brief Function call operator.
         *
         * @param function The function for which to compute the derivative. The function can be any callable type taking a
         * floating point type as an argument, and returns a value of the same type.
         * @param val The value at which to compute the derivative.
         * @param stepsize The finite difference used for computing the derivative.
         *
         * @return The derivative. The return type is the same as the return type of the provided function.
         *
         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
         * or one of the convenience functions.
         *
         * @throws NumerixxError if stepsize is invalid.
         */
        inline auto operator()(IsFunction auto                  function, std::floating_point auto val, std::floating_point auto stepsize) const
        {
            using RETURN_T = std::invoke_result_t< decltype(function), decltype(val) >;
            static_assert(std::floating_point< RETURN_T >, "The return type of the provided function must be a floating point type.");
            validateStepSize(stepsize, std::sqrt(std::numeric_limits< RETURN_T >::epsilon()));
            return FN {}(function, val, stepsize);
        }
    };

    /**
     * @brief A class defining a function object for computing the 1st order derivative of an arbitrary function,
     * using a centered Richardson extrapolation method.
     *
     * @param function The function for which to compute the derivative. The function can be any callable type taking a
     * floating point type as an argument, and returns a value of the same type.
     * @param val The value at which to compute the derivative.
     * @param stepsize The finite difference used for computing the derivative.
     *
     * @return The derivative. The return type is the same as the return type of the provided function.
     *
     * @note This function is not intended to be used directly. Instead, use the \c diff template function,
     * or one of the convenience functions.
     *
     * @details This function uses the Richardson extrapolation method for computing the derivative.
     * It approximates the derivative by taking the difference of the function evaluated at `val + stepsize` and `val - stepsize`,
     * and scaling by the reciprocal of `stepsize * 6`. This method provides a more accurate approximation of the derivative
     * than simple finite difference methods. See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra,
     * for details.
     *
     * @throws NumerixxError if stepsize is invalid.
     */
    constexpr auto Order1CentralRichardsonLambda =
        [](IsFunction auto function, std::floating_point auto val, std::floating_point auto stepsize) {
            return (4.0 * (function(val + stepsize) - function(val - stepsize)) -
                    0.5 * (function(val + 2 * stepsize) - function(val - 2 * stepsize))) /
                   (stepsize * 6);
        };

    using Order1CentralRichardson = SolverTemplate< decltype(Order1CentralRichardsonLambda) >;

    /**
     * @brief A class defining a function object for computing the 1st order derivative of an arbitrary function,
     * using a 3-point centered finite divided-difference method.
     */
    //    class Order1Central3Point
    //    {
    //    public:
    //        /**
    //         * @brief Function call operator.
    //         *
    //         * @param function The function for which to compute the derivative. The function can be any callable type taking a
    //         * floating point type as an argument, and returns a value of the same type.
    //         * @param val The value at which to compute the derivative.
    //         * @param stepsize The finite difference used for computing the derivative.
    //         *
    //         * @return The derivative. The return type is the same as the return type of the provided function.
    //         *
    //         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
    //         * or one of the convenience functions.
    //         *
    //         * @details This function uses the 3-point centered finite divided-difference method for computing the derivative.
    //         * It approximates the derivative by taking the difference of the function evaluated at `val + stepsize` and `val - stepsize`,
    //         * and scaling by the reciprocal of `2 * stepsize`. See chapter 23 in "Numerical Methods for Engineers", 8th Edition by
    //         * Steven C. Chapra, for details.
    //         *
    //         * @throws NumerixxError if stepsize is invalid.
    //         */
    //        inline auto operator()(IsFunction auto                  function,
    //                               std::floating_point auto val,
    //                               std::floating_point auto stepsize) const
    //        {
    //            using RETURN_T = std::invoke_result_t< decltype(function), decltype(val) >;
    //            static_assert(std::floating_point< RETURN_T >, "The return type of the provided function must be a floating point type.");
    //            validateStepSize(stepsize, std::sqrt(std::numeric_limits< RETURN_T >::epsilon()));
    //            return (function(val + stepsize) - function(val - stepsize)) / (2 * stepsize);
    //        }
    //    };

    constexpr auto Order1Central3PointLambda =
        [](IsFunction auto function, std::floating_point auto val, std::floating_point auto stepsize) {
            return (function(val + stepsize) - function(val - stepsize)) / (2 * stepsize);
        };

    using Order1Central3Point = SolverTemplate< decltype(Order1Central3PointLambda) >;

    /**
     * @brief A class defining a function object for computing the 1st order derivative of an arbitrary function,
     * using a 5-point centered finite divided-difference method.
     */
    //    class Order1Central5Point
    //    {
    //    public:
    //        /**
    //         * @brief Function call operator.
    //         *
    //         * @param function The function for which to compute the derivative. The function can be any callable type taking a
    //         * floating point type as an argument, and returns a value of the same type.
    //         * @param val The value at which to compute the derivative.
    //         * @param stepsize The finite difference used for computing the derivative.
    //         *
    //         * @return The derivative. The return type is the same as the return type of the provided function.
    //         *
    //         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
    //         * or one of the convenience functions.
    //         *
    //         * @details This function uses the 5-point centered finite divided-difference method for computing the derivative.
    //         * It approximates the derivative by taking the difference of the function evaluated at `val + 2 * stepsize` and `val - 2 *
    //         * stepsize`, and scaling by the reciprocal of `12 * stepsize`. See chapter 23 in "Numerical Methods for Engineers", 8th
    //         Edition
    //         * by Steven C. Chapra, for details.
    //         *
    //         * @throws NumerixxError if stepsize is invalid.
    //         */
    //        inline auto operator()(IsFunction auto function, std::floating_point auto val, std::floating_point auto stepsize) const
    //        {
    //            using RETURN_T = std::invoke_result_t< decltype(function), decltype(val) >;
    //            static_assert(std::floating_point< RETURN_T >, "The return type of the provided function must be a floating point type.");
    //            validateStepSize(stepsize, std::sqrt(std::numeric_limits< RETURN_T >::epsilon()));
    //            return (-function(val + 2 * stepsize) + 8 * function(val + stepsize) - 8 * function(val - stepsize) +
    //                    function(val - 2 * stepsize)) /
    //                   (12 * stepsize);
    //        }
    //    };

    constexpr auto Order1Central5PointLambda =
        [](IsFunction auto function, std::floating_point auto val, std::floating_point auto stepsize) {
            return (-function(val + 2 * stepsize) + 8 * function(val + stepsize) - 8 * function(val - stepsize) +
                    function(val - 2 * stepsize)) /
                   (12 * stepsize);
        };

    using Order1Central5Point = SolverTemplate< decltype(Order1Central5PointLambda) >;

    /**
     * @brief A class defining a function object for computing the 2nd order derivative of an arbitrary function,
     * using a 3-point centered finite divided-difference method.
     */
    //    class Order2Central3Point
    //    {
    //    public:
    //        /**
    //         * @brief Function call operator.
    //         *
    //         * @param function The function for which to compute the derivative. The function can be any callable type taking a
    //         * floating point type as an argument, and returns a value of the same type.
    //         * @param val The value at which to compute the derivative.
    //         * @param stepsize The finite difference used for computing the derivative.
    //         *
    //         * @return The derivative. The return type is the same as the return type of the provided function.
    //         *
    //         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
    //         * or one of the convenience functions.
    //         *
    //         * @details This function uses the 3-point centered finite divided-difference method for computing the derivative.
    //         * It approximates the derivative by taking the difference of the function evaluated at `val + stepsize` and `val - stepsize`,
    //         * and scaling by the reciprocal of `2 * stepsize`. See chapter 23 in "Numerical Methods for Engineers", 8th Edition by
    //         * Steven C. Chapra, for details.
    //         *
    //         * @throws NumerixxError if stepsize is invalid.
    //         */
    //        inline auto operator()(IsFunction auto                  function,
    //                               std::floating_point auto val,
    //                               std::floating_point auto stepsize) const
    //        {
    //            using RETURN_T = std::invoke_result_t< decltype(function), decltype(val) >;
    //            static_assert(std::floating_point< RETURN_T >, "The return type of the provided function must be a floating point type.");
    //            validateStepSize(stepsize, std::cbrt(std::numeric_limits< RETURN_T >::epsilon()));
    //            return (function(val + stepsize) - 2 * function(val) + function(val - stepsize)) / std::pow(stepsize, 2);
    //        }
    //    };

    constexpr auto Order2Central3PointLambda =
        [](IsFunction auto function, std::floating_point auto val, std::floating_point auto stepsize) {
            return (function(val + stepsize) - 2 * function(val) + function(val - stepsize)) / std::pow(stepsize, 2);
        };

    using Order2Central3Point = SolverTemplate< decltype(Order2Central3PointLambda) >;

    /**
     * @brief A class defining a function object for computing the 2nd order derivative of an arbitrary function,
     * using a 5-point centered finite divided-difference method.
     */
    //    class Order2Central5Point
    //    {
    //    public:
    //        /**
    //         * @brief Function call operator.
    //         *
    //         * @param function The function for which to compute the derivative. The function can be any callable type taking a
    //         * floating point type as an argument, and returns a value of the same type.
    //         * @param val The value at which to compute the derivative.
    //         * @param stepsize The finite difference used for computing the derivative.
    //         *
    //         * @return The derivative. The return type is the same as the return type of the provided function.
    //         *
    //         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
    //         * or one of the convenience functions.
    //         *
    //         * @details This function uses the 5-point centered finite divided-difference method for computing the derivative.
    //         * It approximates the derivative by taking the difference of the function evaluated at `val + 2 * stepsize` and `val - 2 *
    //         * stepsize`, and scaling by the reciprocal of `12 * stepsize`. See chapter 23 in "Numerical Methods for Engineers", 8th
    //         Edition
    //         * by Steven C. Chapra, for details.
    //         *
    //         * @throws NumerixxError if stepsize is invalid.
    //         */
    //        inline auto operator()(IsFunction auto                  function,
    //                               std::floating_point auto val,
    //                               std::floating_point auto stepsize) const
    //        {
    //            using RETURN_T = std::invoke_result_t< decltype(function), decltype(val) >;
    //            static_assert(std::floating_point< RETURN_T >, "The return type of the provided function must be a floating point type.");
    //            validateStepSize(stepsize, std::cbrt(std::numeric_limits< RETURN_T >::epsilon()));
    //            return (-function(val + 2 * stepsize) + 16 * function(val + stepsize) - 30 * function(val) + 16 * function(val - stepsize)
    //            -
    //                    function(val - 2 * stepsize)) /
    //                   (12 * std::pow(stepsize, 2));
    //        }
    //    };

    constexpr auto Order2Central5PointLambda =
        [](IsFunction auto function, std::floating_point auto val, std::floating_point auto stepsize) {
            return (-function(val + 2 * stepsize) + 16 * function(val + stepsize) - 30 * function(val) + 16 * function(val - stepsize) -
                    function(val - 2 * stepsize)) /
                   (12 * std::pow(stepsize, 2));
        };

    using Order2Central5Point = SolverTemplate< decltype(Order2Central5PointLambda) >;

    // ====================================================================
    // Forward finite difference formulas
    // ====================================================================

    /**
     * @brief A class defining a function object for computing the 1st order derivative of an arbitrary function,
     * using a forward Richardson extrapolation method.
     */
    //    class Order1ForwardHighOrder
    //    {
    //    public:
    //        /**
    //         * @brief Function call operator.
    //         *
    //         * @param function The function for which to compute the derivative. The function can be any callable type taking a
    //         * floating point type as an argument, and returns a value of the same type.
    //         * @param val The value at which to compute the derivative.
    //         * @param stepsize The finite difference used for computing the derivative.
    //         *
    //         * @return The derivative. The return type is the same as the return type of the provided function.
    //         *
    //         * * @note This function is not intended to be used directly. Instead, use the \c diff template function,
    //         * or one of the convenience functions.
    //         *
    //         * @details This function uses the Richardson extrapolation method for computing the derivative.
    //         * It approximates the derivative by taking the difference of the function evaluated at `val + stepsize` and `val - stepsize`,
    //         * and scaling by the reciprocal of `stepsize * 6`. This method provides a more accurate approximation of the derivative
    //         * than simple finite difference methods. See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C.
    //         Chapra,
    //         * for details.
    //         *
    //         * @throws NumerixxError if stepsize is invalid.
    //         */
    //        inline auto operator()(IsFunction auto                  function,
    //                               std::floating_point auto val,
    //                               std::floating_point auto stepsize) const
    //        {
    //            using RETURN_T = std::invoke_result_t< decltype(function), decltype(val) >;
    //            static_assert(std::floating_point< RETURN_T >, "The return type of the provided function must be a floating point type.");
    //
    //            validateStepSize(stepsize, std::sqrt(std::numeric_limits< RETURN_T >::epsilon()));
    //
    //            const auto diff1 = function(val + stepsize);
    //            const auto diff2 = function(val + stepsize * 2);
    //            const auto diff3 = function(val + stepsize * 3);
    //            const auto diff4 = function(val + stepsize * 4);
    //
    //            return (22.0 * (diff4 - diff3) - 62.0 * (diff3 - diff2) + 52.0 * (diff2 - diff1)) / (stepsize * 12);
    //        }
    //    };

    constexpr auto Order1ForwardRichardsonLambda =
        [](IsFunction auto function, std::floating_point auto val, std::floating_point auto stepsize) {
            const auto diff1 = function(val + stepsize);
            const auto diff2 = function(val + stepsize * 2);
            const auto diff3 = function(val + stepsize * 3);
            const auto diff4 = function(val + stepsize * 4);

            return (22.0 * (diff4 - diff3) - 62.0 * (diff3 - diff2) + 52.0 * (diff2 - diff1)) / (stepsize * 12);
        };

    using Order1ForwardRichardson = SolverTemplate< decltype(Order1ForwardRichardsonLambda) >;

    /**
     * @brief A class defining a function object for computing the 1st order derivative of an arbitrary function,
     * using a 2-point forward finite divided-difference method.
     */
    //    class Order1Forward2Point
    //    {
    //    public:
    //        /**
    //         * @brief Function call operator.
    //         *
    //         * @param function The function for which to compute the derivative. The function can be any callable type taking a
    //         * floating point type as an argument, and returns a value of the same type.
    //         * @param val The value at which to compute the derivative.
    //         * @param stepsize The finite difference used for computing the derivative.
    //         *
    //         * @return The derivative. The return type is the same as the return type of the provided function.
    //         *
    //         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
    //         * or one of the convenience functions.
    //         *
    //         * @details This function uses the 2-point forward finite divided-difference method for computing the derivative.
    //         * It approximates the derivative by taking the difference of the function evaluated at `val + stepsize` and `val`,
    //         * and scaling by the reciprocal of `stepsize`. See chapter 23 in "Numerical Methods for Engineers", 8th Edition by
    //         * Steven C. Chapra, for details.
    //         *
    //         * @throws NumerixxError if stepsize is invalid.
    //         */
    //        inline auto operator()(IsFunction auto                  function,
    //                               std::floating_point auto val,
    //                               std::floating_point auto stepsize) const
    //        {
    //            using RETURN_T = std::invoke_result_t< decltype(function), decltype(val) >;
    //            static_assert(std::floating_point< RETURN_T >, "The return type of the provided function must be a floating point type.");
    //            validateStepSize(stepsize, std::sqrt(std::numeric_limits< RETURN_T >::epsilon()));
    //            return (function(val + stepsize) - function(val)) / stepsize;
    //        }
    //    };

    constexpr auto Order1Forward2PointLambda =
        [](IsFunction auto function, std::floating_point auto val, std::floating_point auto stepsize) {
            return (function(val + stepsize) - function(val)) / stepsize;
        };

    using Order1Forward2Point = SolverTemplate< decltype(Order1Forward2PointLambda) >;

    /**
     * @brief A class defining a function object for computing the 1st order derivative of an arbitrary function,
     * using a 3-point forward finite divided-difference method.
     */
    //    class Order1Forward3Point
    //    {
    //    public:
    //        /**
    //         * @brief Function call operator.
    //         *
    //         * @param function The function for which to compute the derivative. The function can be any callable type taking a
    //         * floating point type as an argument, and returns a value of the same type.
    //         * @param val The value at which to compute the derivative.
    //         * @param stepsize The finite difference used for computing the derivative.
    //         *
    //         * @return The derivative. The return type is the same as the return type of the provided function.
    //         *
    //         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
    //         * or one of the convenience functions.
    //         *
    //         * @details This function uses the 3-point forward finite divided-difference method for computing the derivative.
    //         * It approximates the derivative by taking the difference of the function evaluated at `val + 2 * stepsize` and `val`,
    //         * and scaling by the reciprocal of `2 * stepsize`. See chapter 23 in "Numerical Methods for Engineers", 8th Edition by
    //         * Steven C. Chapra, for details.
    //         *
    //         * @throws NumerixxError if stepsize is invalid.
    //         */
    //        inline auto operator()(IsFunction auto                  function,
    //                               std::floating_point auto val,
    //                               std::floating_point auto stepsize) const
    //        {
    //            using RETURN_T = std::invoke_result_t< decltype(function), decltype(val) >;
    //            static_assert(std::floating_point< RETURN_T >, "The return type of the provided function must be a floating point type.");
    //            validateStepSize(stepsize, std::sqrt(std::numeric_limits< RETURN_T >::epsilon()));
    //            return (-function(val + 2 * stepsize) + 4 * function(val + stepsize) - 3 * function(val)) / (2 * stepsize);
    //        }
    //    };

    constexpr auto Order1Forward3PointLambda =
        [](IsFunction auto function, std::floating_point auto val, std::floating_point auto stepsize) {
            return (-function(val + 2 * stepsize) + 4 * function(val + stepsize) - 3 * function(val)) / (2 * stepsize);
        };

    using Order1Forward3Point = SolverTemplate< decltype(Order1Forward3PointLambda) >;

    /**
     * @brief A class defining a function object for computing the 2nd order derivative of an arbitrary function,
     * using a 3-point forward finite divided-difference method.
     */
    //    class Order2Forward3Point
    //    {
    //    public:
    //        /**
    //         * @brief Function call operator.
    //         *
    //         * @param function The function for which to compute the derivative. The function can be any callable type taking a
    //         * floating point type as an argument, and returns a value of the same type.
    //         * @param val The value at which to compute the derivative.
    //         * @param stepsize The finite difference used for computing the derivative.
    //         *
    //         * @return The derivative. The return type is the same as the return type of the provided function.
    //         *
    //         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
    //         * or one of the convenience functions.
    //         *
    //         * @details This function uses the 3-point forward finite divided-difference method for computing the derivative.
    //         * It approximates the derivative by taking the difference of the function evaluated at `val + 2 * stepsize` and `val`,
    //         * and scaling by the reciprocal of `stepsize^2`. See chapter 23 in "Numerical Methods for Engineers", 8th Edition by
    //         * Steven C. Chapra, for details.
    //         *
    //         * @throws NumerixxError if stepsize is invalid.
    //         */
    //        inline auto operator()(IsFunction auto                  function,
    //                               std::floating_point auto val,
    //                               std::floating_point auto stepsize) const
    //        {
    //            using RETURN_T = std::invoke_result_t< decltype(function), decltype(val) >;
    //            static_assert(std::floating_point< RETURN_T >, "The return type of the provided function must be a floating point type.");
    //            validateStepSize(stepsize, std::cbrt(std::numeric_limits< RETURN_T >::epsilon()));
    //            return (function(val + 2 * stepsize) - 2 * function(val + stepsize) + function(val)) / std::pow(stepsize, 2);
    //        }
    //    };

    constexpr auto Order2Forward3PointLambda =
        [](IsFunction auto function, std::floating_point auto val, std::floating_point auto stepsize) {
            return (function(val + 2 * stepsize) - 2 * function(val + stepsize) + function(val)) / std::pow(stepsize, 2);
        };

    using Order2Forward3Point = SolverTemplate< decltype(Order2Forward3PointLambda) >;

    /**
     * @brief A class defining a function object for computing the 2nd order derivative of an arbitrary function,
     * using a 4-point forward finite divided-difference method.
     */
    //    class Order2Forward4Point
    //    {
    //    public:
    //        /**
    //         * @brief Function call operator.
    //         *
    //         * @param function The function for which to compute the derivative. The function can be any callable type taking a
    //         * floating point type as an argument, and returns a value of the same type.
    //         * @param val The value at which to compute the derivative.
    //         * @param stepsize The finite difference used for computing the derivative.
    //         *
    //         * @return The derivative. The return type is the same as the return type of the provided function.
    //         *
    //         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
    //         * or one of the convenience functions.
    //         *
    //         * @details This function uses the 4-point forward finite divided-difference method for computing the derivative.
    //         * It approximates the derivative by taking the difference of the function evaluated at `val + 3 * stepsize` and `val`,
    //         * and scaling by the reciprocal of `stepsize^2`. See chapter 23 in "Numerical Methods for Engineers", 8th Edition by
    //         * Steven C. Chapra, for details.
    //         *
    //         * @throws NumerixxError if stepsize is invalid.
    //         */
    //        inline auto operator()(IsFunction auto                  function,
    //                               std::floating_point auto val,
    //                               std::floating_point auto stepsize) const
    //        {
    //            using RETURN_T = std::invoke_result_t< decltype(function), decltype(val) >;
    //            static_assert(std::floating_point< RETURN_T >, "The return type of the provided function must be a floating point type.");
    //            validateStepSize(stepsize, std::cbrt(std::numeric_limits< RETURN_T >::epsilon()));
    //            return (-function(val + 3 * stepsize) + 4 * function(val + 2 * stepsize) - 5 * function(val + stepsize) + 2 *
    //            function(val)) /
    //                   std::pow(stepsize, 2);
    //        }
    //    };

    constexpr auto Order2Forward4PointLambda =
        [](IsFunction auto function, std::floating_point auto val, std::floating_point auto stepsize) {
            return (-function(val + 3 * stepsize) + 4 * function(val + 2 * stepsize) - 5 * function(val + stepsize) + 2 * function(val)) /
                   std::pow(stepsize, 2);
        };

    using Order2Forward4Point = SolverTemplate< decltype(Order2Forward4PointLambda) >;

    // ====================================================================
    // Backward finite difference formulas
    // ====================================================================

    /**
     * @brief A class defining a function object for computing the 1st order derivative of an arbitrary function,
     * using a backward Richardson extrapolation method.
     */
    //    class Order1BackwardHighOrder
    //    {
    //    public:
    //        /**
    //         * @brief Function call operator.
    //         *
    //         * @param function The function for which to compute the derivative. The function can be any callable type taking a
    //         * floating point type as an argument, and returns a value of the same type.
    //         * @param val The value at which to compute the derivative.
    //         * @param stepsize The finite difference used for computing the derivative.
    //         *
    //         * @return The derivative. The return type is the same as the return type of the provided function.
    //         *
    //         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
    //         * or one of the convenience functions.
    //         *
    //         * @details This function uses the Richardson extrapolation method for computing the derivative.
    //         * It approximates the derivative by taking the difference of the function evaluated at `val - stepsize` and `val - 4 *
    //         stepsize`,
    //         * and scaling by the reciprocal of `-stepsize * 12`. This method provides a more accurate approximation of the derivative
    //         * than simple finite difference methods. See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C.
    //         Chapra,
    //         * for details.
    //         *
    //         * @throws NumerixxError if stepsize is invalid.
    //         */
    //        inline auto operator()(IsFunction auto                  function,
    //                               std::floating_point auto val,
    //                               std::floating_point auto stepsize) const
    //        {
    //            using RETURN_T = std::invoke_result_t< decltype(function), decltype(val) >;
    //            static_assert(std::floating_point< RETURN_T >, "The return type of the provided function must be a floating point type.");
    //
    //            validateStepSize(stepsize, std::sqrt(std::numeric_limits< RETURN_T >::epsilon()));
    //
    //            const auto diff1 = function(val - stepsize);
    //            const auto diff2 = function(val - stepsize * 2);
    //            const auto diff3 = function(val - stepsize * 3);
    //            const auto diff4 = function(val - stepsize * 4);
    //
    //            return (22.0 * (diff4 - diff3) - 62.0 * (diff3 - diff2) + 52.0 * (diff2 - diff1)) / -(stepsize * 12);
    //        }
    //    };

    constexpr auto Order1BackwardRichardsonLambda =
        [](IsFunction auto function, std::floating_point auto val, std::floating_point auto stepsize) {
            const auto diff1 = function(val - stepsize);
            const auto diff2 = function(val - stepsize * 2);
            const auto diff3 = function(val - stepsize * 3);
            const auto diff4 = function(val - stepsize * 4);

            return (22.0 * (diff4 - diff3) - 62.0 * (diff3 - diff2) + 52.0 * (diff2 - diff1)) / -(stepsize * 12);
        };

    using Order1BackwardRichardson = SolverTemplate< decltype(Order1BackwardRichardsonLambda) >;

    /**
     * @brief A class defining a function object for computing the 1st order derivative of an arbitrary function,
     * using a 2-point backward finite divided-difference method.
     */
    //    class Order1Backward2Point
    //    {
    //    public:
    //        /**
    //         * @brief Function call operator.
    //         *
    //         * @param function The function for which to compute the derivative. The function can be any callable type taking a
    //         * floating point type as an argument, and returns a value of the same type.
    //         * @param val The value at which to compute the derivative.
    //         * @param stepsize The finite difference used for computing the derivative.
    //         *
    //         * @return The derivative. The return type is the same as the return type of the provided function.
    //         *
    //         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
    //         * or one of the convenience functions.
    //         *
    //         * @details This function uses the 2-point backward finite divided-difference method for computing the derivative.
    //         * It approximates the derivative by taking the difference of the function evaluated at `val` and `val - stepsize`,
    //         * and scaling by the reciprocal of `stepsize`. See chapter 23 in "Numerical Methods for Engineers", 8th Edition by
    //         * Steven C. Chapra, for details.
    //         *
    //         * @throws NumerixxError if stepsize is invalid.
    //         */
    //        inline auto operator()(IsFunction auto                  function,
    //                               std::floating_point auto val,
    //                               std::floating_point auto stepsize) const
    //        {
    //            using RETURN_T = std::invoke_result_t< decltype(function), decltype(val) >;
    //            static_assert(std::floating_point< RETURN_T >, "The return type of the provided function must be a floating point type.");
    //
    //            validateStepSize(stepsize, std::sqrt(std::numeric_limits< RETURN_T >::epsilon()));
    //            return (function(val) - function(val - stepsize)) / stepsize;
    //        }
    //    };

    constexpr auto Order1Backward2PointLambda =
        [](IsFunction auto function, std::floating_point auto val, std::floating_point auto stepsize) {
            return (function(val) - function(val - stepsize)) / stepsize;
        };

    using Order1Backward2Point = SolverTemplate< decltype(Order1Backward2PointLambda) >;

    /**
     * @brief A class defining a function object for computing the 1st order derivative of an arbitrary function,
     * using a 3-point backward finite divided-difference method.
     */
    //    class Order1Backward3Point
    //    {
    //    public:
    //        /**
    //         * @brief Function call operator.
    //         *
    //         * @param function The function for which to compute the derivative. The function can be any callable type taking a
    //         * floating point type as an argument, and returns a value of the same type.
    //         * @param val The value at which to compute the derivative.
    //         * @param stepsize The finite difference used for computing the derivative.
    //         *
    //         * @return The derivative. The return type is the same as the return type of the provided function.
    //         *
    //         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
    //         * or one of the convenience functions.
    //         *
    //         * @details This function uses the 3-point backward finite divided-difference method for computing the derivative.
    //         * It approximates the derivative by taking the difference of the function evaluated at `val` and `val - 2 * stepsize`,
    //         * and scaling by the reciprocal of `2 * stepsize`. See chapter 23 in "Numerical Methods for Engineers", 8th Edition by
    //         * Steven C. Chapra, for details.
    //         *
    //         * @throws NumerixxError if stepsize is invalid.
    //         */
    //        inline auto operator()(IsFunction auto                  function,
    //                               std::floating_point auto val,
    //                               std::floating_point auto stepsize) const
    //        {
    //            using RETURN_T = std::invoke_result_t< decltype(function), decltype(val) >;
    //            static_assert(std::floating_point< RETURN_T >, "The return type of the provided function must be a floating point type.");
    //            validateStepSize(stepsize, std::sqrt(std::numeric_limits< RETURN_T >::epsilon()));
    //            return (3 * function(val) - 4 * function(val - stepsize) + function(val - 2 * stepsize)) / (2 * stepsize);
    //        }
    //    };

    constexpr auto Order1Backward3PointLambda =
        [](IsFunction auto function, std::floating_point auto val, std::floating_point auto stepsize) {
            return (3 * function(val) - 4 * function(val - stepsize) + function(val - 2 * stepsize)) / (2 * stepsize);
        };

    using Order1Backward3Point = SolverTemplate< decltype(Order1Backward3PointLambda) >;

    /**
     * @brief A class defining a function object for computing the 2nd order derivative of an arbitrary function,
     * using a 3-point backward finite divided-difference method.
     */
    //    class Order2Backward3Point
    //    {
    //    public:
    //        /**
    //         * @brief Function call operator.
    //         *
    //         * @param function The function for which to compute the derivative. The function can be any callable type taking a
    //         * floating point type as an argument, and returns a value of the same type.
    //         * @param val The value at which to compute the derivative.
    //         * @param stepsize The finite difference used for computing the derivative.
    //         *
    //         * @return The derivative. The return type is the same as the return type of the provided function.
    //         *
    //         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
    //         * or one of the convenience functions.
    //         *
    //         * @details This function uses the 3-point backward finite divided-difference method for computing the derivative.
    //         * It approximates the derivative by taking the difference of the function evaluated at `val` and `val - 2 * stepsize`,
    //         * and scaling by the reciprocal of `stepsize^2`. See chapter 23 in "Numerical Methods for Engineers", 8th Edition by
    //         * Steven C. Chapra, for details.
    //         *
    //         * @throws NumerixxError if stepsize is invalid.
    //         */
    //        inline auto operator()(IsFunction auto                  function,
    //                               std::floating_point auto val,
    //                               std::floating_point auto stepsize) const
    //        {
    //            using RETURN_T = std::invoke_result_t< decltype(function), decltype(val) >;
    //            static_assert(std::floating_point< RETURN_T >, "The return type of the provided function must be a floating point type.");
    //            validateStepSize(stepsize, std::cbrt(std::numeric_limits< RETURN_T >::epsilon()));
    //            return (function(val) - 2 * function(val - stepsize) + function(val - 2 * stepsize)) / std::pow(stepsize, 2);
    //        }
    //    };

    constexpr auto Order2Backward3PointLambda =
        [](IsFunction auto function, std::floating_point auto val, std::floating_point auto stepsize) {
            return (function(val) - 2 * function(val - stepsize) + function(val - 2 * stepsize)) / std::pow(stepsize, 2);
        };

    using Order2Backward3Point = SolverTemplate< decltype(Order2Backward3PointLambda) >;

    /**
     * @brief A class defining a function object for computing the 2nd order derivative of an arbitrary function,
     * using a 4-point backward finite divided-difference method.
     */
    //    class Order2Backward4Point
    //    {
    //    public:
    //        /**
    //         * @brief Function call operator.
    //         *
    //         * @param function The function for which to compute the derivative. The function can be any callable type taking a
    //         * floating point type as an argument, and returns a value of the same type.
    //         * @param val The value at which to compute the derivative.
    //         * @param stepsize The finite difference used for computing the derivative.
    //         *
    //         * @return The derivative. The return type is the same as the return type of the provided function.
    //         *
    //         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
    //         * or one of the convenience functions.
    //         *
    //         * @details This function uses the 4-point backward finite divided-difference method for computing the derivative.
    //         * It approximates the derivative by taking the difference of the function evaluated at `val` and `val - 3 * stepsize`,
    //         * and scaling by the reciprocal of `stepsize^2`. See chapter 23 in "Numerical Methods for Engineers", 8th Edition by
    //         * Steven C. Chapra, for details.
    //         *
    //         * @throws NumerixxError if stepsize is invalid.
    //         */
    //        inline auto operator()(IsFunction auto                  function,
    //                               std::floating_point auto val,
    //                               std::floating_point auto stepsize) const
    //        {
    //            using RETURN_T = std::invoke_result_t< decltype(function), decltype(val) >;
    //            static_assert(std::floating_point< RETURN_T >, "The return type of the provided function must be a floating point type.");
    //            validateStepSize(stepsize, std::cbrt(std::numeric_limits< RETURN_T >::epsilon()));
    //            return (2 * function(val) - 5 * function(val - stepsize) + 4 * function(val - 2 * stepsize) - function(val - 3 *
    //            stepsize)) /
    //                   std::pow(stepsize, 2);
    //        }
    //    };

    constexpr auto Order2Backward4PointLambda =
        [](IsFunction auto function, std::floating_point auto val, std::floating_point auto stepsize) {
            return (2 * function(val) - 5 * function(val - stepsize) + 4 * function(val - 2 * stepsize) - function(val - 3 * stepsize)) /
                   std::pow(stepsize, 2);
        };

    using Order2Backward4Point = SolverTemplate< decltype(Order2Backward4PointLambda) >;

    /**
     * @brief Compute the derivative of a function, using the specified algorithm. This function checks the result
     * for errors and will return a tl::expected object that will contain the result or an error object.
     *
     * @tparam ALGO The algorithm for computing the derivative. This can be one of the provided function objects, or
     * a custom function object type conforming to the required interface.
     * @param function The function for which to compute the derivative. The function can be any callable type taking a
     * floating point type as an argument, and returns a value of the same type.
     * @param val The value at which to compute the derivative.
     * @param stepsize (Optional) The finite difference used to compute the derivative. If an argument for this parameter is not provided,
     * a default value vill be used. The default value is the cubic root of the machine epsilon for the function return type.
     *
     * @return A \c tl::expected (\c std::expected) containing the (approximated) derivative of the function, or (in case of an error)
     * a \c DerivativeError exception object describing the error.
     *
     * @throws NumerixxError if function is not callable.
     */
    template< typename ALGO >
    inline auto diff(IsFunction auto                  function,
                     std::floating_point auto                                  val,
                     std::invoke_result_t< decltype(function), decltype(val) > stepsize =
                         StepSize< std::invoke_result_t< decltype(function), decltype(val) > >)
    {
        if (!function) throw NumerixxError("Function object is invalid.");

        using RETURN_T = std::invoke_result_t< decltype(function), decltype(val) >;
        static_assert(std::floating_point< RETURN_T >, "The return type of the provided function must be a floating point type.");
        using DerivError = Error< detail::DerivErrorData< decltype(val) > >;
        using EXPECTED_T = tl::expected< RETURN_T, DerivError >;

        auto deriv = ALGO {}(function, val, std::max(stepsize, stepsize * val));

        if (std::isfinite(deriv))
            return EXPECTED_T(deriv);
        else
            return EXPECTED_T(tl::make_unexpected(DerivError("Computation of derivative gave non-finite result.",
                                                             NumerixxErrorType::Deriv,
                                                             { .x = val, .h = stepsize, .f = function(val), .df = deriv })));
    }

    /**
     * @brief A convenience function for computing the derivative using a centered divided difference method.
     *
     * @param function The function for which to compute the derivative.
     * @param val The value at which to compute the derivative.
     * @param stepsize (Optional) The finite difference used to compute the derivative. If an argument for this parameter is not provided,
     * a default value vill be used. The default value is the cubic root of the machine epsilon for the function return type.
     *
     * @return A \c tl::expected (\c std::expected) containing the (approximated) derivative of the function, or (in case of an error)
     * a \c DerivativeError exception object describing the error.
     */
    template< typename FN >
    requires IsFunction< FN >
    inline auto central(FN                               function,
                        std::floating_point auto                                  val,
                        std::invoke_result_t< decltype(function), decltype(val) > stepsize =
                            StepSize< std::invoke_result_t< decltype(function), decltype(val) > >)
    {
        return diff< Order1CentralRichardson >(function, val, stepsize);
    }

    /**
     * @brief A convenience function for computing the derivative using a forward divided difference method.
     *
     * @param function The function for which to compute the derivative.
     * @param val The value at which to compute the derivative.
     * @param stepsize (Optional) The finite difference used to compute the derivative. If an argument for this parameter is not provided,
     * a default value vill be used. The default value is the cubic root of the machine epsilon for the function return type.
     *
     * @return A \c tl::expected (\c std::expected) containing the (approximated) derivative of the function, or (in case of an error)
     * a \c DerivativeError exception object describing the error.
     */
    template< typename FN >
    requires IsFunction< FN >
    inline auto forward(FN                               function,
                        std::floating_point auto                                  val,
                        std::invoke_result_t< decltype(function), decltype(val) > stepsize =
                            StepSize< std::invoke_result_t< decltype(function), decltype(val) > >)
    {
        return diff< Order1ForwardRichardson >(function, val, stepsize);
    }

    /**
     * @brief A convenience function for computing the derivative using a backward divided difference method.
     *
     * @param function The function for which to compute the derivative.
     * @param val The value at which to compute the derivative.
     * @param stepsize (Optional) The finite difference used to compute the derivative. If an argument for this parameter is not provided,
     * a default value vill be used. The default value is the cubic root of the machine epsilon for the function return type.
     *
     * @return A \c tl::expected (\c std::expected) containing the (approximated) derivative of the function, or (in case of an error)
     * a \c DerivativeError exception object describing the error.
     */
    template< typename FN >
    requires IsFunction< FN >
    inline auto backward(FN                               function,
                         std::floating_point auto                                  val,
                         std::invoke_result_t< decltype(function), decltype(val) > stepsize =
                             StepSize< std::invoke_result_t< decltype(function), decltype(val) > >)
    {
        return diff< Order1BackwardRichardson >(function, val, stepsize);
    }

    /**
     * @brief Create a function object representing the derivative of the input function, using numerical differentiation.
     *
     * @tparam ALGO The algorithm type for computing the derivative. This can be any algorithm with the right interface,
     * and can be both 1st order and higher order derivatives. The default is \c Order1CentralRichardson.
     * @param function The function to compute the derivative of.
     * @param stepsize (Optional) The step size to use in the computation. The default is the cubic root of the machine epsilon.
     *
     * @return A lambda function representing the derivative of the input function. The lambda will take one floating point
     * argument, and return the (approximated) derivative of the function.
     *
     * @note For objects of the \c Polynomial class, an overload of the \c derivativeOf function is provided for computing the
     * derivative function analytically.
     *
     * @note The returned function object will not check the result for errors. If you want to check the result for errors,
     * use the \c diff function instead, or create a custom function object (e.g., a lambda) that checks the result for errors.
     *
     * @throws NumerixxError if function is not callable.
     */
    template< typename ALGO = Order1CentralRichardson >
    inline auto derivativeOf(IsFunction auto                  function,
        std::invoke_result_t< decltype(function), double > stepsize = StepSize< std::invoke_result_t< decltype(function), double > >)
        requires(!poly::IsPolynomial< decltype(function) >)
    {
        if (!function) throw NumerixxError("Function object is invalid.");

        return [=](std::floating_point auto val) -> std::invoke_result_t< decltype(function), decltype(val) > {
            return ALGO {}(function, val, stepsize);
        };
    }

}    // namespace nxx::deriv

#endif    // NUMERIXX_DIFFERENTIATION_HPP

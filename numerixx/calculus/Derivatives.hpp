/*
    ##    ## ##     ## ##     ## ######## ########  #### ##     ## ##     ##
    ###   ## ##     ## ###   ### ##       ##     ##  ##   ##   ##   ##   ##
    ####  ## ##     ## #### #### ##       ##     ##  ##    ## ##     ## ##
    ## ## ## ##     ## ## ### ## ######   ########   ##     ###       ###
    ##  #### ##     ## ##     ## ##       ##   ##    ##    ## ##     ## ##
    ##   ### ##     ## ##     ## ##       ##    ##   ##   ##   ##   ##   ##
    ##    ##  #######  ##     ## ######## ##     ## #### ##     ## ##     ##

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

#include "DerivativeError.hpp"
#include "../.dependencies/expected/expected.hpp"

#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <tuple>
#include <utility>
#include <type_traits>

namespace nxx::deriv
{
    template< typename FN >
    concept IsFunction = requires(FN fn) {
                             {
                                 fn(0.0)
                             } -> std::floating_point;
                         };

    template< typename SOLVER >
    concept IsSolver = std::invocable< SOLVER, std::function< double(double) >, double, double >;

    template< typename T >
    using ReturnType = std::invoke_result_t< T, double >;


    // ====================================================================
    // Central finite difference formulas
    // ====================================================================

    class Order1CentralRichardson
    {
    public:
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            return (4.0 * (function(val + stepsize) - function(val - stepsize)) -
                    0.5 * (function(val + 2 * stepsize) - function(val - 2 * stepsize))) /
                   (stepsize * 6);
        }
    };

    /**
     * @brief Compute the 1st order derivative of the provided function, using a 3-point centered finite
     * divided-difference formula.
     * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
     * @tparam Fn The type of the function. Can be any callable object, e.g. a lambda.
     * @param function The function for which to compute the derivative.
     * @param val The value at which to compute the derivative.
     * @param stepsize The step size to use when using the divided-difference formula.
     * @return The approximated derivative. The type is the same as the return type of the provided function.
     * @note This function is not intended to be used directly. Instead, use the Derivative class, or one of the
     * convenience functions.
     */
    class Order1Central3Point
    {
    public:
        /**
         * @brief
         * @tparam FN
         * @param function
         * @param val
         * @param stepsize
         * @return
         */
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            return (function(val + stepsize) - function(val - stepsize)) / (2 * stepsize);
        }
    };

    /**
     * @brief Compute the 1st order derivative of the provided function, using a 5-point centered finite
     * divided-difference formula.
     * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
     * @tparam Fn The type of the function. Can be any callable object, e.g. a lambda.
     * @param function The function for which to compute the derivative.
     * @param val The value at which to compute the derivative.
     * @param stepsize The step size to use when using the divided-difference formula.
     * @return The approximated derivative. The type is the same as the return type of the provided function.
     * @note This function is not intended to be used directly. Instead, use the Derivative class, or one of the
     * convenience functions.
     */
    class Order1Central5Point
    {
    public:
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            return (-function(val + 2 * stepsize) + 8 * function(val + stepsize) - 8 * function(val - stepsize) +
                    function(val - 2 * stepsize)) /
                   (12 * stepsize);
        }
    };

    /**
     * @brief Compute the 2nd order derivative of the provided function, using a 3-point centered finite
     * divided-difference formula.
     * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
     * @tparam Fn The type of the function. Can be any callable object, e.g. a lambda.
     * @param function The function for which to compute the derivative.
     * @param val The value at which to compute the derivative.
     * @param stepsize The step size to use when using the divided-difference formula.
     * @return The approximated derivative. The type is the same as the return type of the provided function.
     * @note This function is not intended to be used directly. Instead, use the Derivative class, or one of the
     * convenience functions.
     */
    class Order2Central3Point
    {
    public:
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            return (function(val + stepsize) - 2 * function(val) + function(val - stepsize)) / std::pow(stepsize, 2);
        }
    };

    /**
     * @brief Compute the 2nd order derivative of the provided function, using a 5-point centered finite
     * divided-difference formula.
     * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
     * @tparam Fn The type of the function. Can be any callable object, e.g. a lambda.
     * @param function The function for which to compute the derivative.
     * @param val The value at which to compute the derivative.
     * @param stepsize The step size to use when using the divided-difference formula.
     * @return The approximated derivative. The type is the same as the return type of the provided function.
     * @note This function is not intended to be used directly. Instead, use the Derivative class, or one of the
     * convenience functions.
     */
    class Order2Central5Point
    {
    public:
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            return (-function(val + 2 * stepsize) + 16 * function(val + stepsize) - 30 * function(val) + 16 * function(val - stepsize) -
                    function(val - 2 * stepsize)) /
                   (12 * std::pow(stepsize, 2));
        }
    };

    // ====================================================================
    // Forward finite difference formulas
    // ====================================================================

    class Order1ForwardRichardson
    {
    public:
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            const auto diff1 = function(val + stepsize);
            const auto diff2 = function(val + stepsize * 2);
            const auto diff3 = function(val + stepsize * 3);
            const auto diff4 = function(val + stepsize * 4);

            return (22.0 * (diff4 - diff3) - 62.0 * (diff3 - diff2) + 52.0 * (diff2 - diff1)) / (stepsize * 12);
        }
    };

    /**
     * @brief Compute the 1st order derivative of the provided function, using a 2-point forward finite
     * divided-difference formula.
     * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
     * @tparam Fn The type of the function. Can be any callable object, e.g. a lambda.
     * @param function The function for which to compute the derivative.
     * @param val The value at which to compute the derivative.
     * @param stepsize The step size to use when using the divided-difference formula.
     * @return The approximated derivative. The type is the same as the return type of the provided function.
     * @note This function is not intended to be used directly. Instead, use the Derivative class, or one of the
     * convenience functions.
     */
    class Order1Forward2Point
    {
    public:
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            return (function(val + stepsize) - function(val)) / stepsize;
        }
    };

    /**
     * @brief Compute the 1st order derivative of the provided function, using a 3-point forward finite
     * divided-difference formula.
     * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
     * @tparam Fn The type of the function. Can be any callable object, e.g. a lambda.
     * @param function The function for which to compute the derivative.
     * @param val The value at which to compute the derivative.
     * @param stepsize The step size to use when using the divided-difference formula.
     * @return The approximated derivative. The type is the same as the return type of the provided function.
     * @note This function is not intended to be used directly. Instead, use the Derivative class, or one of the
     * convenience functions.
     */
    class Order1Forward3Point
    {
    public:
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            return (-function(val + 2 * stepsize) + 4 * function(val + stepsize) - 3 * function(val)) / (2 * stepsize);
        }
    };

    /**
     * @brief Compute the 2nd order derivative of the provided function, using a 3-point forward finite
     * divided-difference formula.
     * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
     * @tparam Fn The type of the function. Can be any callable object, e.g. a lambda.
     * @param function The function for which to compute the derivative.
     * @param val The value at which to compute the derivative.
     * @param stepsize The step size to use when using the divided-difference formula.
     * @return The approximated derivative. The type is the same as the return type of the provided function.
     * @note This function is not intended to be used directly. Instead, use the Derivative class, or one of the
     * convenience functions.
     */
    class Order2Forward3Point
    {
    public:
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            return (function(val + 2 * stepsize) - 2 * function(val + stepsize) + function(val)) / std::pow(stepsize, 2);
        }
    };

    /**
     * @brief Compute the 2nd order derivative of the provided function, using a 4-point forward finite
     * divided-difference formula.
     * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
     * @tparam Fn The type of the function. Can be any callable object, e.g. a lambda.
     * @param function The function for which to compute the derivative.
     * @param val The value at which to compute the derivative.
     * @param stepsize The step size to use when using the divided-difference formula.
     * @return The approximated derivative. The type is the same as the return type of the provided function.
     * @note This function is not intended to be used directly. Instead, use the Derivative class, or one of the
     * convenience functions.
     */
    class Order2Forward4Point
    {
    public:
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            return (-function(val + 3 * stepsize) + 4 * function(val + 2 * stepsize) - 5 * function(val + stepsize) + 2 * function(val)) /
                   std::pow(stepsize, 2);
        }
    };

    // ====================================================================
    // Backward finite difference formulas
    // ====================================================================

    class Order1BackwardRichardson
    {
    public:
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            const auto diff1 = function(val - stepsize);
            const auto diff2 = function(val - stepsize * 2);
            const auto diff3 = function(val - stepsize * 3);
            const auto diff4 = function(val - stepsize * 4);

            return (22.0 * (diff4 - diff3) - 62.0 * (diff3 - diff2) + 52.0 * (diff2 - diff1)) / -(stepsize * 12);
        }
    };

    /**
     * @brief Compute the 1st order derivative of the provided function, using a 2-point backward finite
     * divided-difference formula.
     * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
     * @tparam Fn The type of the function. Can be any callable object, e.g. a lambda.
     * @param function The function for which to compute the derivative.
     * @param val The value at which to compute the derivative.
     * @param stepsize The step size to use when using the divided-difference formula.
     * @return The approximated derivative. The type is the same as the return type of the provided function.
     * @note This function is not intended to be used directly. Instead, use the Derivative class, or one of the
     * convenience functions.
     */
    class Order1Backward2Point
    {
    public:
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            return (function(val) - function(val - stepsize)) / stepsize;
        }
    };

    /**
     * @brief Compute the 1st order derivative of the provided function, using a 3-point backward finite
     * divided-difference formula.
     * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
     * @tparam Fn The type of the function. Can be any callable object, e.g. a lambda.
     * @param function The function for which to compute the derivative.
     * @param val The value at which to compute the derivative.
     * @param stepsize The step size to use when using the divided-difference formula.
     * @return The approximated derivative. The type is the same as the return type of the provided function.
     * @note This function is not intended to be used directly. Instead, use the Derivative class, or one of the
     * convenience functions.
     */
    class Order1Backward3Point
    {
    public:
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            return (3 * function(val) - 4 * function(val - stepsize) + function(val - 2 * stepsize)) / (2 * stepsize);
        }
    };

    /**
     * @brief Compute the 2nd order derivative of the provided function, using a 3-point backward finite
     * divided-difference formula.
     * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
     * @tparam Fn The type of the function. Can be any callable object, e.g. a lambda.
     * @param function The function for which to compute the derivative.
     * @param val The value at which to compute the derivative.
     * @param stepsize The step size to use when using the divided-difference formula.
     * @return The approximated derivative. The type is the same as the return type of the provided function.
     * @note This function is not intended to be used directly. Instead, use the Derivative class, or one of the
     * convenience functions.
     */
    class Order2Backward3Point
    {
    public:
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            return (function(val) - 2 * function(val - stepsize) + function(val - 2 * stepsize)) / std::pow(stepsize, 2);
        }
    };

    /**
     * @brief Compute the 2nd order derivative of the provided function, using a 4-point backward finite
     * divided-difference formula.
     * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
     * @tparam Fn The type of the function. Can be any callable object, e.g. a lambda.
     * @param function The function for which to compute the derivative.
     * @param val The value at which to compute the derivative.
     * @param stepsize The step size to use when using the divided-difference formula.
     * @return The approximated derivative. The type is the same as the return type of the provided function.
     * @note This function is not intended to be used directly. Instead, use the Derivative class, or one of the
     * convenience functions.
     */
    class Order2Backward4Point
    {
    public:
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            return (2 * function(val) - 5 * function(val - stepsize) + 4 * function(val - 2 * stepsize) - function(val - 3 * stepsize)) /
                   std::pow(stepsize, 2);
        }
    };

    /**
     * @brief
     * @tparam ALGO
     * @param function
     * @param val
     * @param stepsize
     * @return
     */
    template< typename ALGO >
    inline auto derivative(
        IsFunction auto                  function,
        ReturnType< decltype(function) > val,
        ReturnType< decltype(function) > stepsize = std::cbrt(std::numeric_limits< ReturnType< decltype(function) > >::epsilon()))
    {
        tl::expected< ReturnType< decltype(function) >, error::DerivativeError > result;

        auto solver = ALGO {};
        auto deriv  = solver(function, val, stepsize);

        if (std::isfinite(deriv))
            result = deriv;
        else
            result = tl::make_unexpected(error::DerivativeError { "Computation of derivative gave non-finite result." });

        return result;
    }

    /**
     * @brief
     * @tparam Fn
     * @param function
     * @param val
     * @param order
     * @return
     */
    inline auto
        central(IsFunction auto                  function,
                ReturnType< decltype(function) > val,
                ReturnType< decltype(function) > stepsize = std::cbrt(std::numeric_limits< ReturnType< decltype(function) > >::epsilon()))
    {
        return derivative< Order1CentralRichardson >(function, val, stepsize);
    }

    /**
     * @brief
     * @tparam Fn
     * @param function
     * @param val
     * @param order
     * @return
     */
    inline auto
        forward(IsFunction auto                  function,
                ReturnType< decltype(function) > val,
                ReturnType< decltype(function) > stepsize = std::cbrt(std::numeric_limits< ReturnType< decltype(function) > >::epsilon()))
    {
        return derivative< Order1ForwardRichardson >(function, val, stepsize);
    }

    /**
     * @brief
     * @tparam Fn
     * @param function
     * @param val
     * @param order
     * @return
     */
    inline auto
        backward(IsFunction auto                  function,
                 ReturnType< decltype(function) > val,
                 ReturnType< decltype(function) > stepsize = std::cbrt(std::numeric_limits< ReturnType< decltype(function) > >::epsilon()))
    {
        return derivative< Order1BackwardRichardson >(function, val, stepsize);
    }

}    // namespace nxx::deriv

#endif    // NUMERIXX_DIFFERENTIATION_HPP

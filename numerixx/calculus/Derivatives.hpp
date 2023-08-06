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

#ifndef NUMERIXX_DIFFERENTIATION_HPP
#define NUMERIXX_DIFFERENTIATION_HPP

// ===== Numerixx Includes
#include "../poly/Polynomial.hpp"

// ===== External Includes
#include <tl/expected.hpp>
#include <gcem.hpp>

// ===== Standard Library Includes
#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <tuple>
#include <type_traits>
#include <utility>

namespace nxx::deriv
{

    class DerivativeError : public std::runtime_error
    {
    public:
        explicit DerivativeError(const char* msg)
            : std::runtime_error(msg) {}
    };

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

    template< typename T >
        requires std::floating_point< T >
    struct StepSizeHelper
    {
        static constexpr T value = gcem::pow(std::numeric_limits< T >::epsilon(), 1.0 / 3.0);
    };

    template< typename T >
        requires std::floating_point< T >
    inline constexpr T StepSize = StepSizeHelper< T >::value;

    // ====================================================================
    // Central finite difference formulas
    // ====================================================================

    /**
     * @brief A class defining a function object for computing the 1st order derivative of an arbitrary function,
     * using a centered Richardson extrapolation method.
     */
    class Order1CentralRichardson
    {
    public:
        /**
         * @brief Function call operator.
         * @param function The function for which to compute the derivative. The function can be any callable type taking a
         * floating point type as an argument, and returns a value of the same type.
         * @param val The value at which to compute the derivative.
         * @param stepsize The finite difference used for computing the derivative.
         * @return The derivative. The return type is the same as the retrun type of the provided function.
         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
         * or one of the convenience functions.
         * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
         */
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            return (4.0 * (function(val + stepsize) - function(val - stepsize)) -
                    0.5 * (function(val + 2 * stepsize) - function(val - 2 * stepsize))) /
                   (stepsize * 6);
        }
    };

    /**
     * @brief A class defining a function object for computing the 1st order derivative of an arbitrary function,
     * using a 3-point centered finite divided-difference method.
     */
    class Order1Central3Point
    {
    public:
        /**
         * @brief Function call operator.
         * @param function The function for which to compute the derivative. The function can be any callable type taking a
         * floating point type as an argument, and returns a value of the same type.
         * @param val The value at which to compute the derivative.
         * @param stepsize The finite difference used for computing the derivative.
         * @return The derivative. The return type is the same as the retrun type of the provided function.
         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
         * or one of the convenience functions.
         * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
         */
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            return (function(val + stepsize) - function(val - stepsize)) / (2 * stepsize);
        }
    };

    /**
     * @brief A class defining a function object for computing the 1st order derivative of an arbitrary function,
     * using a 5-point centered finite divided-difference method.
     */
    class Order1Central5Point
    {
    public:
        /**
         * @brief Function call operator.
         * @param function The function for which to compute the derivative. The function can be any callable type taking a
         * floating point type as an argument, and returns a value of the same type.
         * @param val The value at which to compute the derivative.
         * @param stepsize The finite difference used for computing the derivative.
         * @return The derivative. The return type is the same as the retrun type of the provided function.
         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
         * or one of the convenience functions.
         * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
         */
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            return (-function(val + 2 * stepsize) + 8 * function(val + stepsize) - 8 * function(val - stepsize) +
                    function(val - 2 * stepsize)) /
                   (12 * stepsize);
        }
    };

    /**
     * @brief A class defining a function object for computing the 2nd order derivative of an arbitrary function,
     * using a 3-point centered finite divided-difference method.
     */
    class Order2Central3Point
    {
    public:
        /**
         * @brief Function call operator.
         * @param function The function for which to compute the derivative. The function can be any callable type taking a
         * floating point type as an argument, and returns a value of the same type.
         * @param val The value at which to compute the derivative.
         * @param stepsize The finite difference used for computing the derivative.
         * @return The derivative. The return type is the same as the retrun type of the provided function.
         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
         * or one of the convenience functions.
         * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
         */
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            return (function(val + stepsize) - 2 * function(val) + function(val - stepsize)) / std::pow(stepsize, 2);
        }
    };

    /**
     * @brief A class defining a function object for computing the 2nd order derivative of an arbitrary function,
     * using a 5-point centered finite divided-difference method.
     */
    class Order2Central5Point
    {
    public:
        /**
         * @brief Function call operator.
         * @param function The function for which to compute the derivative. The function can be any callable type taking a
         * floating point type as an argument, and returns a value of the same type.
         * @param val The value at which to compute the derivative.
         * @param stepsize The finite difference used for computing the derivative.
         * @return The derivative. The return type is the same as the retrun type of the provided function.
         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
         * or one of the convenience functions.
         * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
         */
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

    /**
     * @brief A class defining a function object for computing the 1st order derivative of an arbitrary function,
     * using a forward Richardson extrapolation method.
     */
    class Order1ForwardRichardson
    {
    public:
        /**
         * @brief Function call operator.
         * @param function The function for which to compute the derivative. The function can be any callable type taking a
         * floating point type as an argument, and returns a value of the same type.
         * @param val The value at which to compute the derivative.
         * @param stepsize The finite difference used for computing the derivative.
         * @return The derivative. The return type is the same as the retrun type of the provided function.
         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
         * or one of the convenience functions.
         * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
         */
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
     * @brief A class defining a function object for computing the 1st order derivative of an arbitrary function,
     * using a 2-point forward finite divided-difference method.
     */
    class Order1Forward2Point
    {
    public:
        /**
         * @brief Function call operator.
         * @param function The function for which to compute the derivative. The function can be any callable type taking a
         * floating point type as an argument, and returns a value of the same type.
         * @param val The value at which to compute the derivative.
         * @param stepsize The finite difference used for computing the derivative.
         * @return The derivative. The return type is the same as the retrun type of the provided function.
         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
         * or one of the convenience functions.
         * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
         */
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            return (function(val + stepsize) - function(val)) / stepsize;
        }
    };

    /**
     * @brief A class defining a function object for computing the 1st order derivative of an arbitrary function,
     * using a 3-point forward finite divided-difference method.
     */
    class Order1Forward3Point
    {
    public:
        /**
         * @brief Function call operator.
         * @param function The function for which to compute the derivative. The function can be any callable type taking a
         * floating point type as an argument, and returns a value of the same type.
         * @param val The value at which to compute the derivative.
         * @param stepsize The finite difference used for computing the derivative.
         * @return The derivative. The return type is the same as the retrun type of the provided function.
         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
         * or one of the convenience functions.
         * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
         */
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            return (-function(val + 2 * stepsize) + 4 * function(val + stepsize) - 3 * function(val)) / (2 * stepsize);
        }
    };

    /**
     * @brief A class defining a function object for computing the 2nd order derivative of an arbitrary function,
     * using a 3-point forward finite divided-difference method.
     */
    class Order2Forward3Point
    {
    public:
        /**
         * @brief Function call operator.
         * @param function The function for which to compute the derivative. The function can be any callable type taking a
         * floating point type as an argument, and returns a value of the same type.
         * @param val The value at which to compute the derivative.
         * @param stepsize The finite difference used for computing the derivative.
         * @return The derivative. The return type is the same as the retrun type of the provided function.
         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
         * or one of the convenience functions.
         * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
         */
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            return (function(val + 2 * stepsize) - 2 * function(val + stepsize) + function(val)) / std::pow(stepsize, 2);
        }
    };

    /**
     * @brief A class defining a function object for computing the 2nd order derivative of an arbitrary function,
     * using a 4-point forward finite divided-difference method.
     */
    class Order2Forward4Point
    {
    public:
        /**
         * @brief Function call operator.
         * @param function The function for which to compute the derivative. The function can be any callable type taking a
         * floating point type as an argument, and returns a value of the same type.
         * @param val The value at which to compute the derivative.
         * @param stepsize The finite difference used for computing the derivative.
         * @return The derivative. The return type is the same as the retrun type of the provided function.
         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
         * or one of the convenience functions.
         * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
         */
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

    /**
     * @brief A class defining a function object for computing the 1st order derivative of an arbitrary function,
     * using a backward Richardson extrapolation method.
     */
    class Order1BackwardRichardson
    {
    public:
        /**
         * @brief Function call operator.
         * @param function The function for which to compute the derivative. The function can be any callable type taking a
         * floating point type as an argument, and returns a value of the same type.
         * @param val The value at which to compute the derivative.
         * @param stepsize The finite difference used for computing the derivative.
         * @return The derivative. The return type is the same as the retrun type of the provided function.
         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
         * or one of the convenience functions.
         * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
         */
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
     * @brief A class defining a function object for computing the 1st order derivative of an arbitrary function,
     * using a 2-point backward finite divided-difference method.
     */
    class Order1Backward2Point
    {
    public:
        /**
         * @brief Function call operator.
         * @param function The function for which to compute the derivative. The function can be any callable type taking a
         * floating point type as an argument, and returns a value of the same type.
         * @param val The value at which to compute the derivative.
         * @param stepsize The finite difference used for computing the derivative.
         * @return The derivative. The return type is the same as the retrun type of the provided function.
         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
         * or one of the convenience functions.
         * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
         */
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            return (function(val) - function(val - stepsize)) / stepsize;
        }
    };

    /**
     * @brief A class defining a function object for computing the 1st order derivative of an arbitrary function,
     * using a 3-point backward finite divided-difference method.
     */
    class Order1Backward3Point
    {
    public:
        /**
         * @brief Function call operator.
         * @param function The function for which to compute the derivative. The function can be any callable type taking a
         * floating point type as an argument, and returns a value of the same type.
         * @param val The value at which to compute the derivative.
         * @param stepsize The finite difference used for computing the derivative.
         * @return The derivative. The return type is the same as the retrun type of the provided function.
         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
         * or one of the convenience functions.
         * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
         */
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            return (3 * function(val) - 4 * function(val - stepsize) + function(val - 2 * stepsize)) / (2 * stepsize);
        }
    };

    /**
     * @brief A class defining a function object for computing the 2nd order derivative of an arbitrary function,
     * using a 3-point backward finite divided-difference method.
     */
    class Order2Backward3Point
    {
    public:
        /**
         * @brief Function call operator.
         * @param function The function for which to compute the derivative. The function can be any callable type taking a
         * floating point type as an argument, and returns a value of the same type.
         * @param val The value at which to compute the derivative.
         * @param stepsize The finite difference used for computing the derivative.
         * @return The derivative. The return type is the same as the retrun type of the provided function.
         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
         * or one of the convenience functions.
         * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
         */
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            return (function(val) - 2 * function(val - stepsize) + function(val - 2 * stepsize)) / std::pow(stepsize, 2);
        }
    };

    /**
     * @brief A class defining a function object for computing the 2nd order derivative of an arbitrary function,
     * using a 4-point backward finite divided-difference method.
     */
    class Order2Backward4Point
    {
    public:
        /**
         * @brief Function call operator.
         * @param function The function for which to compute the derivative. The function can be any callable type taking a
         * floating point type as an argument, and returns a value of the same type.
         * @param val The value at which to compute the derivative.
         * @param stepsize The finite difference used for computing the derivative.
         * @return The derivative. The return type is the same as the retrun type of the provided function.
         * @note This function is not intended to be used directly. Instead, use the \c diff template function,
         * or one of the convenience functions.
         * @details See chapter 23 in "Numerical Methods for Engineers", 8th Edition by Steven C. Chapra, for details.
         */
        inline auto
            operator()(IsFunction auto function, ReturnType< decltype(function) > val, ReturnType< decltype(function) > stepsize) const
        {
            return (2 * function(val) - 5 * function(val - stepsize) + 4 * function(val - 2 * stepsize) - function(val - 3 * stepsize)) /
                   std::pow(stepsize, 2);
        }
    };

    /**
     * @brief Compute the derivative of a function, using the specified algorithm. This function checks the result
     * for errors and will return a tl::expected object that will contain the result or an error object.
     * @tparam ALGO The algorithm for computing the derivative. This can be one of the provided function objects, or
     * a custom function object type conforming to the required interface.
     * @param function The function for which to compute the derivative.
     * @param val The value at which to compute the derivative.
     * @param stepsize (Optional) The finite difference used to compute the derivative. If an argument for this parameter is not provided,
     * a default value vill be used. The default value is the cubic root of the machine epsilon for the function return type.
     * @return A \c tl::expected (\c std::expected) containing the (approximated) derivative of the function, or (in case of an error)
     * a \c DerivativeError exception object describing the error.
     */
    template< typename ALGO >
    inline auto diff(IsFunction auto                  function,
                     ReturnType< decltype(function) > val,
                     ReturnType< decltype(function) > stepsize = StepSize< ReturnType< decltype(function) > >)
    {
        tl::expected< ReturnType< decltype(function) >, DerivativeError > result;
        auto                                                              deriv = ALGO {}(function, val, std::max(stepsize, stepsize * val));

        if (std::isfinite(deriv))
            result = deriv;
        else
            result = tl::make_unexpected(DerivativeError { "Computation of derivative gave non-finite result." });
        return result;
    }

    /**
     * @brief A convenience function for computing the derivative using a centered divided difference method.
     * @param function The function for which to compute the derivative.
     * @param val The value at which to compute the derivative.
     * @param stepsize (Optional) The finite difference used to compute the derivative. If an argument for this parameter is not provided,
     * a default value vill be used. The default value is the cubic root of the machine epsilon for the function return type.
     * @return A \c tl::expected (\c std::expected) containing the (approximated) derivative of the function, or (in case of an error)
     * a \c DerivativeError exception object describing the error.
     */
    template< typename FN >
        requires IsFunction< FN >
    inline auto central(FN                               function,
                        ReturnType< decltype(function) > val,
                        ReturnType< decltype(function) > stepsize = StepSize< ReturnType< decltype(function) > >)
    {
        return diff< Order1CentralRichardson >(function, val, stepsize);
    }

    /**
     * @brief A convenience function for computing the derivative using a forward divided difference method.
     * @param function The function for which to compute the derivative.
     * @param val The value at which to compute the derivative.
     * @param stepsize (Optional) The finite difference used to compute the derivative. If an argument for this parameter is not provided,
     * a default value vill be used. The default value is the cubic root of the machine epsilon for the function return type.
     * @return A \c tl::expected (\c std::expected) containing the (approximated) derivative of the function, or (in case of an error)
     * a \c DerivativeError exception object describing the error.
     */
    template< typename FN >
    requires IsFunction< FN >
    inline auto forward(FN                               function,
                        ReturnType< decltype(function) > val,
                        ReturnType< decltype(function) > stepsize = StepSize< ReturnType< decltype(function) > >)
    {
        return diff< Order1ForwardRichardson >(function, val, stepsize);
    }

    /**
     * @brief A convenience function for computing the derivative using a backward divided difference method.
     * @param function The function for which to compute the derivative.
     * @param val The value at which to compute the derivative.
     * @param stepsize (Optional) The finite difference used to compute the derivative. If an argument for this parameter is not provided,
     * a default value vill be used. The default value is the cubic root of the machine epsilon for the function return type.
     * @return A \c tl::expected (\c std::expected) containing the (approximated) derivative of the function, or (in case of an error)
     * a \c DerivativeError exception object describing the error.
     */
    template< typename FN >
    requires IsFunction< FN >
    inline auto backward(FN                               function,
                         ReturnType< decltype(function) > val,
                         ReturnType< decltype(function) > stepsize = StepSize< ReturnType< decltype(function) > >)
    {
        return diff< Order1BackwardRichardson >(function, val, stepsize);
    }

    /**
     * @brief Create a function object representing the derivative of the input function, using numerical differentiation.
     * @tparam ALGO The algorithm type for computing the derivative. This can be any algorithm with the right interface,
     * and can be both 1st order and higher order derivatives. The default is \c Order1CentralRichardson.
     * @param function The function to compute the derivative of.
     * @param stepsize (Optional) The step size to use in the computation. The default is the cubic root of the machine epsilon.
     * @return A lambda function representing the derivative of the input function. The lambda will take one floating point
     * argument, and return the (approximated) derivative of the function.
     * @note For objects of the \c Polynomial class, an overload of the \c derivativeOf function is provided for computing the
     * derivative function analytically.
     * @note The returned function object will not check the result for errors. If you want to check the result for errors,
     * use the \c diff function instead, or create a custom function object (e.g., a lambda) that checks the result for errors.
     */
    template< typename ALGO = Order1CentralRichardson >
    inline auto derivativeOf(IsFunction auto                  function,
                             ReturnType< decltype(function) > stepsize = StepSize< ReturnType< decltype(function) > >)
        requires(!poly::IsPolynomial< decltype(function) >)
    {
        using RT = ReturnType< decltype(function) >;
        return [=](RT val) { return ALGO {}(function, val, stepsize); };
    }

}    // namespace nxx::deriv

#endif    // NUMERIXX_DIFFERENTIATION_HPP

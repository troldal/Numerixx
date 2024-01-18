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

#ifndef NUMERIXX_FUNCTION_HPP
#define NUMERIXX_FUNCTION_HPP

// ===== Numerixx Includes
#include <Concepts.hpp>

// ===== External Includes
#include <tl/expected.hpp>

namespace nxx::func
{

    /**
     * @brief A custom error class for the Function class.
     */
    class FunctionError : public std::runtime_error
    {
    public:
        explicit FunctionError(const char* msg)
            : std::runtime_error(msg)
        {}
    };

    /**
     * @brief A wrapper class for functions that accepts and returns either double or std::complex<double>.
     *
     * @tparam FN The type of the function that is wrapped by the Function class.
     * @requires FN should be invocable with double or std::complex<double> as its argument type.
     */
    template< typename FN >
    requires std::invocable< FN, double > || std::invocable< FN, std::complex< double > >
    class Function
    {
    public:
        /**
         * @brief Constructs a Function object.
         *
         * @param fn The function to be wrapped by the Function object.
         */
        explicit Function(FN fn)
            : m_fn(fn)
        {}

        /**
         * @brief Constructs a Function object from an rvalue reference.
         *
         * @param fn The function to be wrapped by the Function object.
         */
        explicit Function(FN&& fn)
            : m_fn(std::forward< FN >(fn))
        {}

        /**
         * @brief Function call operator.
         *
         * @tparam T The type of the input argument.
         * @param x The input value for which the function will be evaluated.
         * @return The result of the function evaluation.
         */
        template< typename T >
        requires(std::floating_point< T > && std::floating_point< std::invoke_result_t< FN, T > >) ||
                (IsComplex< T > && IsComplex< std::invoke_result_t< FN, T > >)
        auto operator()(T x) const
        {
            return m_fn(std::forward< T >(x));
        }

        /**
         * @brief Evaluates the wrapped function for the given input and returns the result or an error.
         *
         * @tparam T The type of the input argument.
         * @param x The input value for which the function will be evaluated.
         * @return A tl::expected object containing either the result of the function evaluation or a FunctionError.
         */
        template< typename T >
        requires(std::floating_point< T > && std::floating_point< std::invoke_result_t< FN, T > >) ||
                (IsComplex< T > && IsComplex< std::invoke_result_t< FN, T > >)
        auto evaluate(T x) const
        {
            using RT = std::invoke_result_t< FN, T >;
            tl::expected< RT, FunctionError > result;
            auto                              eval = m_fn(std::forward< T >(x));

            if (std::isfinite(abs(eval)))
                result = eval;
            else
                result = tl::make_unexpected(FunctionError { "Function evaluation gave non-finite result." });
            return result;
        }

    private:
        FN m_fn; /**< The wrapped function object. */
    };

}    // namespace nxx::func

#endif    // NUMERIXX_FUNCTION_HPP

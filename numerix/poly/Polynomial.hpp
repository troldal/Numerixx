/*
    888b      88  88        88  88b           d88  88888888888  88888888ba   88  8b        d8
    8888b     88  88        88  888b         d888  88           88      "8b  88   Y8,    ,8P
    88 `8b    88  88        88  88`8b       d8'88  88           88      ,8P  88    `8b  d8'
    88  `8b   88  88        88  88 `8b     d8' 88  88aaaaa      88aaaaaa8P'  88      Y88P
    88   `8b  88  88        88  88  `8b   d8'  88  88"""""      88""""88'    88      d88b
    88    `8b 88  88        88  88   `8b d8'   88  88           88    `8b    88    ,8P  Y8,
    88     `8888  Y8a.    .a8P  88    `888'    88  88           88     `8b   88   d8'    `8b
    88      `888   `"Y8888Y"'   88     `8'     88  88888888888  88      `8b  88  8P        Y8

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

#ifndef NUMERIX_POLYNOMIAL_HPP
#define NUMERIX_POLYNOMIAL_HPP

#include <algorithm>
#include <numeric>
#include <type_traits>
#include <vector>

namespace numerix::poly
{

    /**
     * @brief The Polynomial class implements the functionality of a polynomial with arbitrary number of coefficients.
     * @tparam T The type of the polynomial coefficients (must be of floating-point type)
     */
    template<typename T>
        requires std::is_floating_point_v<T>
    class Polynomial final
    {
        std::vector<T> m_coefficients; /**< The internal store of polynomial coefficients. */

    public:
        /*
         * Alias declarations.
         */
        using value_type = T;

        /**
         * @brief Constructor taking a container of coefficients as an argument.
         * @param coefficients The container of coefficients. The value type must be of floating-point type, and must
         * support iterators.
         */
        template<typename CONTAINER>
        explicit Polynomial(CONTAINER coefficients) : m_coefficients { coefficients.begin(), coefficients.end() }
        {}

        /**
         * @brief Constructor taking an std::initializer_list with the polynomial coefficients as an argument.
         * @param coefficients The list of coefficients of the polynomial.
         */
        Polynomial(std::initializer_list<T> coefficients) : m_coefficients { coefficients } {}

        /**
         * @brief Function call operator, for evaluating the polynomial at a given value.
         * @param value The value at which to evaluate the polynomial.
         * @return The result of the evaluation.
         */
        inline value_type operator()(value_type value) const { return evaluate(value); }

        /**
         * @brief Function for evaluating the polynomial at a given value.
         * @param value The value at which to evaluate the polynomial.
         * @return The result of the evaluation.
         */
        inline value_type evaluate(value_type value) const
        {
            // ===== Horner's method implemented in terms of std::accummulate.
            return std::accumulate(m_coefficients.rbegin() + 1,
                                   m_coefficients.rend(),
                                   *m_coefficients.rbegin(),
                                   [value](value_type curr, value_type coeff) { return curr * value + coeff; });
        }

        /**
         * @brief Evaluate the first derivative at a given value.
         * @param value The value at which to evaluate the derivative.
         * @return The result of the evaluation.
         */
        inline value_type derivative(value_type value) const
        {
            std::vector<value_type> newCoeff = m_coefficients;
            int                     n        = newCoeff.size() - 1;

            // ===== Convert the coefficients, by multiplying by the exponent
            std::transform(newCoeff.crbegin(), newCoeff.crend(), newCoeff.rbegin(), [&n](value_type elem) { return elem * n--; });

            // ===== The order of the resulting polynomial is one less than the original.
            // ===== This handled by rotating the elements around the 2nd element and setting the last to zero.
            std::rotate(newCoeff.begin(), newCoeff.begin() + 1, newCoeff.end());
            newCoeff.back() = 0.0;

            return Polynomial(newCoeff).evaluate(value);
        }

        /**
         * @brief Get the vector of polynomial coefficients.
         * @return A const reference to the vector of coefficients.
         */
        const std::vector<value_type>& coefficients() const { return m_coefficients; }

        /**
         * @brief Get the polynomial coefficients in an arbitrary container.
         * @tparam CONTAINER The type of the container. The value type must be floating type, and it must support iterators.
         * @return A CONTAINER object with the polynomial coefficients.
         */
        template<typename CONTAINER>
            requires std::is_floating_point_v<typename CONTAINER::value_type>
        CONTAINER coefficients() const
        {
            CONTAINER result { m_coefficients.begin(), m_coefficients.end() };
            return result;
        }
    };

    /*
     * Deduction guide for using arbitrary containers for coefficient input.
     */
    template<typename CONTAINER>
    Polynomial(CONTAINER coefficients) -> Polynomial<typename CONTAINER::value_type>;
}    // namespace numerix::poly

#endif    // NUMERIX_POLYNOMIAL_HPP

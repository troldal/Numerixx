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
     * @brief
     * @tparam Cont
     */
    template<typename Cont = std::vector<double>>
    requires std::is_floating_point_v<typename Cont::value_type>
    class polynomial final {

        using RT = typename Cont::value_type;

        Cont m_coefficients;

    public:

        using value_type = typename Cont::value_type;
        using container_type = Cont;

        /**
         * @brief
         * @param coefficients
         */
        explicit polynomial(Cont coefficients) : m_coefficients{coefficients} {}

        /**
         * @brief
         * @param value
         * @return
         */
        inline RT operator()(RT value) const {
            return evaluate(value);
        }

        /**
         * @brief
         * @param value
         * @return
         */
        inline RT evaluate(RT value) const {

            // ===== Horner's method implemented in terms of std::accummulate.
            return std::accumulate(m_coefficients.rbegin() + 1, m_coefficients.rend(), *m_coefficients.rbegin(), [value](RT curr, RT coeff) {
                return curr * value + coeff;
            });
        }

        /**
         * @brief
         * @param value
         * @return
         */
        inline RT derivative(RT value) const {
            Cont newCoeff = m_coefficients;
            int n = newCoeff.size() - 1;

            // ===== Convert the coefficients, by multiplying by the exponent
            std::transform(newCoeff.crbegin(), newCoeff.crend(), newCoeff.rbegin(), [&n](RT elem){
                return elem * n--; });

            // ===== The order of the resulting polynomial is one less than the original.
            // ===== This handled by rotating the elements around the 2nd element and setting the last to zero.
            std::rotate(newCoeff.begin(), newCoeff.begin() + 1, newCoeff.end());
            newCoeff.back() = 0.0;

            return polynomial(newCoeff).evaluate(value);
        }

        /**
         * @brief
         * @return
         */
        const Cont& coefficients() const {
            return m_coefficients;
        }

    };
}

#endif    // NUMERIX_POLYNOMIAL_HPP

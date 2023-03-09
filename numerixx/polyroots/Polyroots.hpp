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

#ifndef NUMERIXX_POLYROOTS_HPP
#define NUMERIXX_POLYROOTS_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numbers>
#include <vector>

namespace nxx::polyroots
{

    /**
     * @brief Compute the real roots of a quadratic polynomial.
     * @tparam POLYTYPE The type of the polynomial (must be poly::Polynomial, with value_type = floating point).
     * @param poly The Polynomial object.
     * @return The roots of the polynomial.
     */
    template<typename POLYTYPE>
        requires std::is_floating_point_v<typename POLYTYPE::value_type> &&
                 std::same_as<POLYTYPE, poly::Polynomial<typename POLYTYPE::value_type>>
    inline auto quadratic(POLYTYPE poly)
    {
        using result_vector = std::vector<typename POLYTYPE::value_type>;
        auto& coeff         = poly.coefficients();

        // ===== Check that the polynomial is quadratic.
        if (coeff.size() != 3) throw std::invalid_argument("Polynomial Error: Polynomial is not quadratic.");

        // ===== Compute the discriminant. If it is negative, there are no real roots.
        auto sq_num = coeff[1] * coeff[1] - 4 * coeff[2] * coeff[0];
        if (sq_num < 0.0) return result_vector {};

        // ===== Compute the roots, sort them lowest to highest, and return them.
        auto roots = result_vector { (-coeff[1] + std::sqrt(sq_num)) / (2 * coeff[2]), (-coeff[1] - std::sqrt(sq_num)) / (2 * coeff[2]) };
        std::sort(roots.begin(), roots.end());
        return roots;
    }

    /**
     * @brief Compute the real roots of a cubic polynomial.
     * @tparam POLYTYPE The type of the polynomial (must be poly::Polynomial, with value_type = floating point).
     * @param poly The Polynomial object.
     * @return The roots of the polynomial.
     */
    template<typename POLYTYPE>
        requires std::is_floating_point_v<typename POLYTYPE::value_type> &&
                 std::same_as<POLYTYPE, poly::Polynomial<typename POLYTYPE::value_type>>
    inline auto cubic(POLYTYPE poly)
    {
        using std::acos;
        using std::cbrt;
        using std::cos;
        using std::sqrt;
        using namespace std::numbers;

        using result_vector = std::vector<typename POLYTYPE::value_type>;
        using value_type    = typename POLYTYPE::value_type;

        auto coeff = poly.coefficients();
        if (coeff.size() != 4) throw std::invalid_argument("Polynomial Error: Polynomial is not cubic.");

        std::transform(coeff.cbegin(), coeff.cend(), coeff.begin(), [&coeff](value_type elem) { return elem / coeff.back(); });

        auto& a_0 = coeff[0];
        auto& a_1 = coeff[1];
        auto& a_2 = coeff[2];

        // ===== Compute the constants required for an analytic solution.
        auto p = (1.0 / 3.0) * (3 * a_1 - pow(a_2, 2));
        auto q = (1.0 / 27.0) * (2 * pow(a_2, 3) - 9 * a_2 * a_1 + 27 * a_0);
        auto R = (pow(q, 2) / 4.0) + (pow(p, 3) / 27.0);

        // ===== If R <= 0, there are three real roots
        if (R <= 0.0) {
            auto m     = 2 * sqrt(-p / 3);
            auto theta = acos(3 * q / (p * m)) / 3.0;

            result_vector roots { m * cos(theta) - a_2 / 3,
                                  m * cos(theta + 2 * pi_v<value_type> / 3) - a_2 / 3,
                                  m * cos(theta + 4 * pi_v<value_type> / 3) - a_2 / 3 };
            std::sort(roots.begin(), roots.end());

            return roots;
        }

        // ===== If R > 0, there is one real root
        auto P = cbrt(-q / 2.0 + sqrt(R));
        auto Q = cbrt(-q / 2.0 - sqrt(R));

        return result_vector { P + Q - a_2 / 3.0 };
    }

}    // namespace numerix::polyroots

#endif    // NUMERIXX_POLYROOTS_HPP

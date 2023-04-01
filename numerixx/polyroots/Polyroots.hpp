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

#ifndef NUMERIXX_POLYROOTS_HPP
#define NUMERIXX_POLYROOTS_HPP

#include "../poly/Polynomial.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numbers>
#include <vector>

namespace nxx::polyroots
{

    /**
     * @brief Finds the root of a linear polynomial (monomial) using basic algebraic manipulation.
     *
     * This function accepts a linear polynomial (monomial) as input and calculates its root using
     * basic algebraic manipulation. It returns the root as a single value of the same type as the
     * polynomial coefficients.
     *
     * @param poly A linear polynomial, which should satisfy the poly::IsPolynomial concept.
     *
     * @return The root of the linear polynomial.
     *
     * @throws std::invalid_argument If the input polynomial is not a monomial.
     *
     * @note The function checks if the input polynomial is linear by verifying its order.
     * If the order is not equal to 1, an std::invalid_argument exception is thrown.
     */
    inline auto linear(poly::IsPolynomial auto poly) {
        if (poly.order() != 1) throw std::invalid_argument("Polynomial Error: Input is not a monomial.");
        return -poly.coefficients().front() / poly.coefficients().back() ;
    }

    /**
     * @brief Finds the roots of a quadratic polynomial using the quadratic formula, returning either complex
     * or real roots depending on the RT template parameter.
     *
     * This function accepts a quadratic polynomial as input and calculates its roots using the quadratic formula.
     * It returns a vector containing the roots, either real or complex, depending on the discriminant of the
     * polynomial and the RT template parameter.
     *
     * @tparam RT The desired return type for the roots. Defaults to void, which will return the same type as the
     * polynomial coefficients. If specified, the roots will be of type RT.
     * @param poly A quadratic polynomial, which should satisfy the poly::IsPolynomial concept.
     *
     * @return A vector containing the roots of the quadratic polynomial. If RT is void, the return type will be
     * a vector of complex numbers if the input polynomial has complex coefficients, or a vector of real numbers
     * if the input polynomial has real coefficients. If RT is specified, the return type will be a vector of RT.
     *
     * @throws std::invalid_argument If the input polynomial is not quadratic.
     *
     * @note The function checks if the input polynomial is quadratic by verifying the size of its coefficients
     * vector. If the size is not equal to 3, an std::invalid_argument exception is thrown.
     * @note The discriminant is used to determine the nature of the roots. If the discriminant is negative, the
     * roots will be complex; otherwise, they will be real.
     * @note The returned roots are either complex or real, depending on the provided RT template parameter. If
     * the return type is complex, all roots will be returned. If the return type is real, only roots with
     * imaginary parts smaller than a specified tolerance (1e-6) will be returned.
     */
    template<typename RT = void>
    inline auto quadratic(poly::IsPolynomial auto poly)
    {
        using VALUE_TYPE = typename poly::PolynomialTraits< decltype(poly) >::value_type;
        using FLOAT_TYPE = typename poly::PolynomialTraits< decltype(poly) >::fundamental_type;
        using COMPLEX_TYPE = std::complex< FLOAT_TYPE >;

        auto& coeff         = poly.coefficients();

        // ===== Check that the polynomial is quadratic.
        if (coeff.size() != 3) throw std::invalid_argument("Polynomial Error: Polynomial is not quadratic.");

        // ===== Compute the discriminant. If it is negative, there are no real roots.
        auto sq_num = coeff[1] * coeff[1] - 4.0 * coeff[2] * coeff[0];

        // ===== Compute the roots, sort them lowest to highest, and return them.
        auto roots = std::vector<COMPLEX_TYPE> { (-coeff[1] + std::sqrt(sq_num)) / (2.0 * coeff[2]), (-coeff[1] - std::sqrt(sq_num)) / (2.0 * coeff[2]) };
        std::sort(roots.begin(), roots.end(), [](COMPLEX_TYPE lhs, COMPLEX_TYPE rhs) { return lhs.real() < rhs.real(); });

        // If the return type is complex, return the roots as is.
        if constexpr (utils::IsComplex< RT > || (std::same_as< RT, void > && utils::IsComplex< VALUE_TYPE >)) {
            return roots;

            // Otherwise, return only the real roots.
        } else {
            using RET = std::conditional_t< std::same_as< RT, void >, VALUE_TYPE, RT >;
            auto realroots = std::vector< RET > {};
            std::for_each(roots.begin(), roots.end(), [&realroots](COMPLEX_TYPE elem) { if (abs(elem.imag()) < 1e-6) realroots.push_back(elem.real()); });
            return realroots;
        }
    }

    /**
     * @brief Finds the roots of a cubic polynomial using an analytic solution, returning either
     * complex or real roots depending on the RT template parameter.
     *
     * This function accepts a cubic polynomial as input and calculates its roots using an
     * analytic solution. It returns a vector containing the roots, either real or complex,
     * depending on the discriminant of the polynomial and the RT template parameter.
     *
     * @tparam RT The desired return type for the roots. Defaults to void, which will return the
     * same type as the polynomial coefficients. If specified, the roots will be of type RT.
     * @param poly A cubic polynomial, which should satisfy the poly::IsPolynomial concept.
     *
     * @return A vector containing the roots of the cubic polynomial. If RT is void, the return
     * type will be a vector of complex numbers if the input polynomial has complex coefficients,
     * or a vector of real numbers if the input polynomial has real coefficients. If RT is specified,
     * the return type will be a vector of RT.
     *
     * @throws std::invalid_argument If the input polynomial is not cubic.
     * @throws static_assert If the input polynomial has a value type that is not floating point.
     *
     * @note The function checks if the input polynomial is cubic by verifying the size of its
     * coefficients vector. If the size is not equal to 4, an std::invalid_argument exception is thrown.
     * @note The roots are calculated using an analytic solution that differentiates between the cases
     * of having three real roots and one real root with two complex conjugate roots.
     * @note The returned roots are either complex or real, depending on the provided RT template
     * parameter. If the return type is complex, all roots will be returned. If the return type is real,
     * only roots with imaginary parts smaller than a specified tolerance (1e-6) will be returned.
     * @note Polynomials with complex coefficients are not supported. If the input polynomial has complex
     * coefficients, a static_assert will be triggered.
     */
    template<typename RT = void>
    inline auto cubic(poly::IsPolynomial auto poly)
    {
        static_assert(std::is_floating_point_v< typename poly::PolynomialTraits< decltype(poly) >::value_type >,
                      "Polynomial Error: Polynomial value type must be floating point.");

        using VALUE_TYPE = typename poly::PolynomialTraits< decltype(poly) >::value_type;
        using FLOAT_TYPE = typename poly::PolynomialTraits< decltype(poly) >::fundamental_type;
        using COMPLEX_TYPE = std::complex< FLOAT_TYPE >;

        using std::acos;
        using std::cbrt;
        using std::cos;
        using std::sqrt;
        using namespace std::numbers;

        auto coeff = poly.coefficients();
        if (coeff.size() != 4) throw std::invalid_argument("Polynomial Error: Polynomial is not cubic.");

        std::transform(coeff.cbegin(), coeff.cend(), coeff.begin(), [&coeff](FLOAT_TYPE elem) { return elem / coeff.back(); });

        auto& a_0 = coeff[0];
        auto& a_1 = coeff[1];
        auto& a_2 = coeff[2];

        // ===== Compute the constants required for an analytic solution.
        auto p = (1.0 / 3.0) * (3 * a_1 - pow(a_2, 2));
        auto q = (1.0 / 27.0) * (2 * pow(a_2, 3) - 9 * a_2 * a_1 + 27 * a_0);
        auto R = (pow(q, 2) / 4.0) + (pow(p, 3) / 27.0);

        std::vector< COMPLEX_TYPE > roots {};

        // ===== If R <= 0, there are three real roots
        if (R <= 0.0) {
            auto m     = 2 * sqrt(-p / 3);
            auto theta = acos(3 * q / (p * m)) / 3.0;

            roots = { m * cos(theta) - a_2 / 3,
                      m * cos(theta + 2 * pi_v< FLOAT_TYPE > / 3) - a_2 / 3,
                      m * cos(theta + 4 * pi_v< FLOAT_TYPE > / 3) - a_2 / 3 };
        } else {
            // ===== If R > 0, there is one real root and two complex conjugate roots.
            auto P = cbrt(-q / 2.0 + sqrt(R));
            auto Q = cbrt(-q / 2.0 - sqrt(R));

            roots = { P + Q - a_2 / 3.0,
                      { -0.5 * (P + Q) - a_2 / 3.0, 0.5 * sqrt(3.0) * (P - Q) },
                      { -0.5 * (P + Q) - a_2 / 3.0, -0.5 * sqrt(3.0) * (P - Q) } };
        }

        std::sort(roots.begin(), roots.end(), [](COMPLEX_TYPE a, COMPLEX_TYPE b) { return a.real() < b.real(); });
        // If the return type is complex, return the roots as is.
        if constexpr (utils::IsComplex< RT > || (std::same_as< RT, void > && utils::IsComplex< VALUE_TYPE >)) {
            return roots;

            // Otherwise, return only the real roots.
        } else {
            using RET = std::conditional_t< std::same_as< RT, void >, VALUE_TYPE, RT >;
            auto realroots = std::vector< RET > {};
            std::for_each(roots.begin(), roots.end(), [&realroots](COMPLEX_TYPE elem) { if (abs(elem.imag()) < 1e-6) realroots.push_back(elem.real()); });
            return realroots;
        }
    }

    /**
     * @brief Finds an approximate root of a polynomial using Laguerre's method, given an initial guess.
     *
     * This function utilizes Laguerre's method to find an approximate root of a polynomial equation.
     * It takes a polynomial and an optional initial guess for a root as input, and returns the
     * approximate root after a fixed number of iterations (default: 100) or when the difference between
     * iterations is below a specified tolerance (1E-12).
     *
     * @param poly A polynomial, which should satisfy the poly::IsPolynomial concept.
     * @param guess An optional initial guess for a root of the polynomial. Defaults to 2.0.
     *
     * @return An approximate root of the polynomial as a std::complex. Even if the polynomial is real,
     * the root may be complex due to the nature of the Laguerre's method.
     *
     * @note The function calculates the first and second derivatives of the input polynomial to perform the
     * Laguerre's method iterations.
     * @note The Laguerre's method is an iterative root-finding technique that converges rapidly for most
     * polynomials. However, it may fail to converge for certain ill-conditioned polynomials. In such cases,
     * using a different root-finding method may be necessary.
     */
    inline auto laguerre(poly::IsPolynomial auto poly, std::complex<typename poly::PolynomialTraits< decltype(poly) >::fundamental_type> guess = 2.0)
    {
        using FLOAT_TYPE = poly::PolynomialTraits< decltype(poly) >::fundamental_type;
        using COMPLEX_TYPE = std::complex< FLOAT_TYPE >;

        // ===== Lambda function for computing the Laguerre step.
        auto laguerrestep = [] (COMPLEX_TYPE g_param, COMPLEX_TYPE h_param){
            auto temp = std::sqrt(2.0 * (3.0 * h_param - g_param*g_param));
            return 3.0 / (abs(g_param + temp) > abs(g_param - temp) ? (g_param + temp) : (g_param - temp));
        };

        COMPLEX_TYPE root = guess;
        COMPLEX_TYPE G;
        COMPLEX_TYPE H;
        COMPLEX_TYPE step;

        auto d1poly = derivativeOf(poly);
        auto d2poly = derivativeOf(d1poly);

        for (int i = 0; i < 100; ++i) {
            G = *d1poly(root) / *poly(root);
            H = G * G - *d2poly(root) / *poly(root);
            step = laguerrestep(G, H);
            if (abs(step) < 1E-12) break;
            root = root - step;
        }

        return root;
    }

    /**
     * @brief Solves a polynomial equation using Laguerre's method and the quadratic formula, returning
     * either complex or real roots depending on the RT template parameter.
     *
     * This function accepts a polynomial as input and solves it using a combination of Laguerre's method
     * and the quadratic formula. If the polynomial is of degree higher than 2, the Laguerre method is used
     * to find the roots. For quadratics, the quadratic formula is used. The roots can be returned as complex
     * or real numbers depending on the RT template parameter.
     *
     * @tparam RT The desired return type for the roots. Defaults to void, which will return the same type as
     * the polynomial coefficients. If specified, the roots will be of type RT.
     * @param poly A polynomial, which should satisfy the poly::IsPolynomial concept. The input polynomial
     * can have real or complex coefficients.
     *
     * @return A vector containing the roots of the polynomial. If RT is void, the return type will be a
     * vector of complex numbers if the input polynomial has complex coefficients, or a vector of real numbers
     * if the input polynomial has real coefficients. If RT is specified, the return type will be a vector of RT.
     *
     * @note The implementation uses Laguerre's method for polynomials of degree higher than 2 and the quadratic
     * formula for quadratics. The Laguerre method is an iterative root-finding technique that converges rapidly
     * for most polynomials. A polishing step is performed after finding each root using the Laguerre method to
     * improve the accuracy of the root.
     * @note The returned roots are either complex or real, depending on the provided RT template parameter. If
     * the return type is complex, all roots will be returned. If the return type is real, only roots with
     * imaginary parts smaller than a specified tolerance (1e-6) will be returned.
     */
    template<typename RT = void>
    inline auto polysolve(poly::IsPolynomial auto poly)
    {
        using VALUE_TYPE = typename poly::PolynomialTraits< decltype(poly) >::value_type;
        using FLOAT_TYPE = typename poly::PolynomialTraits< decltype(poly) >::fundamental_type;
        using COMPLEX_TYPE = std::complex< FLOAT_TYPE >;

        auto polynomial = poly::Polynomial< COMPLEX_TYPE >(std::vector<COMPLEX_TYPE>{ poly.begin(), poly.end() });
        auto roots      = std::vector< COMPLEX_TYPE > {};

        // If the polynomial is higher than quadratic, use the Laguerre method.
        while (polynomial.order() > 2) {
            roots.emplace_back(laguerre(polynomial));
            roots.back() = laguerre(poly, roots.back()); // Polishing step
            polynomial = polynomial / poly::Polynomial< COMPLEX_TYPE > { -roots.back(), 1.0 };
        }

        // If the polynomial is quadratic, use the quadratic formula.
        auto quadroots = quadratic(polynomial);
        roots.insert(roots.end(), quadroots.begin(), quadroots.end());
        std::sort(roots.begin(), roots.end(), [](COMPLEX_TYPE a, COMPLEX_TYPE b) { return a.real() < b.real(); });

        // If the return type is complex, return the roots as is.
        if constexpr (utils::IsComplex< RT > || (std::same_as< RT, void > && utils::IsComplex< VALUE_TYPE >)) {
            return roots;

        // Otherwise, return only the real roots.
        } else {
            using RET = std::conditional_t< std::same_as< RT, void >, VALUE_TYPE, RT >;
            auto realroots = std::vector< RET > {};
            std::for_each(roots.begin(), roots.end(), [&realroots](COMPLEX_TYPE elem) { if (abs(elem.imag()) < 1e-6) realroots.push_back(elem.real()); });
            return realroots;
        }
    }
}    // namespace numerix::polyroots

#endif    // NUMERIXX_POLYROOTS_HPP

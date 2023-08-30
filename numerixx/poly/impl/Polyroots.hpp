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

// ===== Numerixx Includes
#include "Polynomial.hpp"
#include "PolynomialError.hpp"
#include <Constants.hpp>
#include <Roots.hpp>

// ===== Standard Library Includes
#include <algorithm>
#include <cassert>
#include <cmath>
#include <numbers>
#include <random>
#include <vector>

namespace nxx::poly
{
    namespace impl
    {

        /**
         * @brief Sorts a vector of complex polynomial roots and returns them as the specified return type.
         *
         * This function takes a vector of complex polynomial roots, sorts them based on their real parts,
         * and returns them as the specified return type. If the return type is a floating point, only the
         * real roots will be returned.
         *
         * @tparam RT The desired return type. Must be either a floating point type or a complex type.
         * @param roots A vector of complex polynomial roots.
         * @return A vector of roots in the specified return type. If RT is a floating point type, only
         * real roots will be returned.
         *
         * @note The input roots must be of a complex type.
         * @note The input roots vector must be of type std::vector.
         */
        template< typename RT >
        requires(std::floating_point< RT > || IsComplex< RT >)
        inline auto sortAndReturn(auto roots, auto tolerance)
            requires IsComplex< typename decltype(roots)::value_type > &&
                     std::same_as< decltype(roots), std::vector< typename decltype(roots)::value_type > >
        {
            // Sorting function
            auto sortingFunc = [tolerance](auto const& a, auto const& b) {
                // NOTE: Without sqrt of the tolerance, the sorting may fail.
                return std::abs(b.real() - a.real()) < std::sqrt(tolerance) ? a.imag() < b.imag() : a.real() < b.real();
            };

            // Sort the roots
            std::sort(roots.begin(), roots.end(), sortingFunc);

            // If the return type is complex, return the roots as is.
            if constexpr (IsComplex< RT >) {
                return roots;
            }
            // Otherwise, return only the real roots.
            else {
                auto realroots = std::vector< RT > {};
                std::for_each(roots.begin(), roots.end(), [&realroots, tolerance](auto elem) {
                    // NOTE: Without sqrt of the tolerance, the test fails.
                    if (std::abs(elem.imag()) < std::sqrt(tolerance)) realroots.push_back(elem.real());
                });
                return realroots;
            }
        }
    }    // namespace impl

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
     * @throws error::PolynomialError If the input polynomial is not a monomial.
     *
     * @note The function checks if the input polynomial is linear by verifying its order.
     * If the order is not equal to 1, an std::invalid_argument exception is thrown.
     */
    inline auto linear(IsPolynomial auto poly)
    {
        if (poly.order() != 1) throw error::PolynomialError("Polynomial Error: Input is not a monomial.");
        return -poly.coefficients().front() / poly.coefficients().back();
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
     * @throws error::PolynomialError If the input polynomial is not quadratic.
     *
     * @note The function checks if the input polynomial is quadratic by verifying the size of its coefficients
     * vector. If the size is not equal to 3, an std::invalid_argument exception is thrown.
     * @note The discriminant is used to determine the nature of the roots. If the discriminant is negative, the
     * roots will be complex; otherwise, they will be real.
     * @note The returned roots are either complex or real, depending on the provided RT template parameter. If
     * the return type is complex, all roots will be returned. If the return type is real, only roots with
     * imaginary parts smaller than a specified tolerance (1e-6) will be returned.
     */
    template< typename RT = void >
    inline auto quadratic(IsPolynomial auto poly, typename PolynomialTraits< decltype(poly) >::fundamental_type tolerance = nxx::EPS)
    {
        // ===== Check that the polynomial is quadratic.
        if (poly.order() != 2) throw error::PolynomialError("Quadratic Error: Polynomial is not quadratic.");

        using POLY_TYPE    = PolynomialTraits< decltype(poly) >;
        using VALUE_TYPE   = typename POLY_TYPE::value_type;
        using FLOAT_TYPE   = typename POLY_TYPE::fundamental_type;
        using COMPLEX_TYPE = std::complex< FLOAT_TYPE >;

        const auto& coeffs = poly.coefficients();
        const auto& a = coeffs[2];
        const auto& b = coeffs[1];
        const auto& c = coeffs[0];

        COMPLEX_TYPE discriminant = sqrt(b * b - 4.0 * a * c);
        COMPLEX_TYPE sqrt_component = std::conj(b) * discriminant;

        COMPLEX_TYPE q = -0.5 * (b + (sqrt_component.real() >= 0.0 ? discriminant : -discriminant));

        if (std::abs(q) < tolerance || std::abs(a) < tolerance) throw error::PolynomialError("Quadratic polynomial is ill formed.");

        std::vector< COMPLEX_TYPE > roots = { q / a, c / q };

        using RET = std::conditional_t< std::same_as< RT, void >, VALUE_TYPE, RT >;
        return impl::sortAndReturn< RET >(roots, tolerance);
    }

    /**
     * @brief Finds the roots of a cubic polynomial using an analytic solution, returning either
     * complex or real roots depending on the RT template parameter.
     *
     * This function accepts a cubic polynomial as input and calculates its roots using an
     * analytic solution. It returns a vector containing the roots, either real or complex,
     * depending on the discriminant of the polynomial and the RT template parameter.
     *
     * @tparam RT The desired return type. If not specified, the function will return roots with
     * the same type as the polynomial value type.
     * @param poly The input polynomial. Must be a cubic polynomial.
     * @return A vector of roots in the specified return type. If RT is a floating point type,
     * only real roots will be returned.
     *
     * @throws error::PolynomialError if the input polynomial is not cubic.
     *
     * @note The function checks if the input polynomial is cubic by verifying the size of its coefficients
     * vector. If the size is not equal to 4, an error::PolynomialError exception is thrown.
     * @note The returned roots are either complex or real, depending on the provided RT template parameter. If
     * the return type is complex, all roots will be returned. If the return type is real, only roots with
     * imaginary parts smaller than a specified tolerance (1e-6) will be returned.
     */
    template< typename RT = void >
    inline auto cubic(IsPolynomial auto poly, typename PolynomialTraits< decltype(poly) >::fundamental_type tolerance = nxx::EPS)
    {
        if (poly.order() != 3) throw error::PolynomialError("Cubic Error: Polynomial is not cubic.");

        using POLY_TYPE    = PolynomialTraits< decltype(poly) >;
        using VALUE_TYPE   = typename POLY_TYPE::value_type;
        using FLOAT_TYPE   = typename POLY_TYPE::fundamental_type;
        using COMPLEX_TYPE = std::complex< FLOAT_TYPE >;

        using std::sqrt;
        using namespace std::numbers;
        using namespace std::complex_literals;

        // ===== Ad hoc lambda function to calculate the cube root of a complex number.
        auto cbrt = [](COMPLEX_TYPE x) { return std::pow(x, 1.0 / 3.0); };

        auto coeff = poly.coefficients();
        std::transform(coeff.cbegin(), coeff.cend(), coeff.begin(), [&coeff](auto elem) { return elem / coeff.back(); });

        const auto& a = coeff[2];
        const auto& b = coeff[1];
        const auto& c = coeff[0];

        COMPLEX_TYPE Q = (a * a - 3.0 * b) / 9.0;
        COMPLEX_TYPE R = (2.0 * a * a * a - 9.0 * a * b + 27.0 * c) / 54.0;
        COMPLEX_TYPE A =
            -cbrt(R + ((std::conj(R) * sqrt(R * R - Q * Q * Q)).real() >= 0.0 ? sqrt(R * R - Q * Q * Q) : -sqrt(R * R - Q * Q * Q)));
        COMPLEX_TYPE B = (abs(A) == 0.0 ? 0.0 : Q / A);

        std::vector< COMPLEX_TYPE > roots = { A + B - a / 3.0,
                                              -0.5 * (A + B) - a / 3.0 + 0.5 * sqrt(3.0) * (A - B) * 1.0i,
                                              -0.5 * (A + B) - a / 3.0 - 0.5 * sqrt(3.0) * (A - B) * 1.0i };

        using RET = std::conditional_t< std::same_as< RT, void >, VALUE_TYPE, RT >;
        return impl::sortAndReturn< RET >(roots, tolerance);
    }

    /**
     * @brief Finds an approximate root of a polynomial using Laguerre's method, given an initial guess.
     *
     * This function utilizes Laguerre's method to find an approximate root of a polynomial equation.
     * It takes a polynomial and an optional initial guess for a root as input, and returns the
     * approximate root after a fixed number of iterations or when the difference between
     * iterations is below a specified tolerance.
     *
     * @param poly A polynomial, which should satisfy the poly::IsPolynomial concept.
     * @param guess An optional initial guess for a root of the polynomial. Defaults to 1.0.
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
    template< typename POLY >
    requires IsPolynomial< POLY >
    inline auto laguerre(POLY                                                                          poly,
                         std::complex< typename PolynomialTraits< decltype(poly) >::fundamental_type > guess          = 1.0,
                         typename PolynomialTraits< decltype(poly) >::fundamental_type                 tolerance      = nxx::EPS,
                         int                                                                           max_iterations = nxx::MAXITER)
    {
        if (poly.order() <= 3) throw error::PolynomialError("Laguerre Error: Polynomial is cubic or lower.");

        using POLY_TYPE    = PolynomialTraits< decltype(poly) >;
        using FLOAT_TYPE   = typename POLY_TYPE::fundamental_type;
        using COMPLEX_TYPE = std::complex< FLOAT_TYPE >;

        // ===== Lambda function for computing the Laguerre step.
        auto laguerrestep = [](COMPLEX_TYPE g_param, COMPLEX_TYPE h_param) {
            COMPLEX_TYPE temp = std::sqrt(2.0 * (3.0 * h_param - g_param * g_param));
            return 3.0 / (abs(g_param + temp) > abs(g_param - temp) ? (g_param + temp) : (g_param - temp));
        };

        if (abs(poly(guess)) < tolerance) return guess;

        COMPLEX_TYPE root = guess;
        COMPLEX_TYPE G;
        COMPLEX_TYPE H;
        COMPLEX_TYPE step;

        auto d1poly = derivativeOf(poly);
        auto d2poly = derivativeOf(d1poly);

        // ===== Use a random number generator to perturb the step size every 10 iterations.
        std::random_device                       rd;
        std::mt19937                             mt(rd());
        std::uniform_real_distribution< double > dist(0.9, 1.1);

        // TODO: Return a std::expected if the Laguerre method fails to converge.
        // ===== Perform the Laguerre iterations.
        for (int i = 0; i < max_iterations; ++i) {
            G    = d1poly(root) / poly(root);
            H    = G * G - d2poly(root) / poly(root);
            step = laguerrestep(G, H);
            if (!std::isfinite(abs(step))) step = 0.1;    // ===== If the step is not finite, use a small value.
            if (abs(step) < tolerance) break;             // ===== If the step is below the tolerance, stop.
            if (i % 10 == 0) step *= dist(mt);            // ===== Perturb the step size every 10 iterations.
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
     * @param poly A polynomial, which should satisfy the IsPolynomial concept. The input polynomial can have
     * real or complex coefficients.
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
     * imaginary parts smaller than the specified tolerance will be returned.
     */
    template< typename RT = void >
    inline auto polysolve(IsPolynomial auto                                             poly,
                          typename PolynomialTraits< decltype(poly) >::fundamental_type tolerance      = nxx::EPS,
                          int                                                           max_iterations = nxx::MAXITER)
    {
        if (max_iterations < 1) throw error::PolynomialError("Maximum number of iterations must be greater than zero.");
        if (poly.order() < 1) throw error::PolynomialError("Polynomial must have at least two coefficients (a monomial).");

        using POLY_TYPE    = PolynomialTraits< decltype(poly) >;
        using VALUE_TYPE   = typename POLY_TYPE::value_type;
        using FLOAT_TYPE   = typename POLY_TYPE::fundamental_type;
        using COMPLEX_TYPE = std::complex< FLOAT_TYPE >;

        // Create a polynomial object and a vector for the roots. The value type of the polynomial is converted to
        // COMPLEX_TYPE if it is not already complex.
        auto polynomial = Polynomial< COMPLEX_TYPE >(std::vector< COMPLEX_TYPE > { poly.begin(), poly.end() });
        auto original   = Polynomial< COMPLEX_TYPE >(std::vector< COMPLEX_TYPE > { poly.begin(), poly.end() });
        auto roots      = std::vector< COMPLEX_TYPE > {};

        // A lambda for computing the root of the polynomial using Newton's method:
        auto newt = [=](auto f, auto x) {
            auto df = derivativeOf(f);
            for (int i = 0; i < max_iterations; ++i) {
                auto dx = f(x) / df(x);
                x -= dx;
                auto fval = f(x);
                if (std::abs(fval.real()) < tolerance &&
                    std::abs(fval.imag()) < tolerance &&
                    std::abs(dx.real()) < tolerance &&
                    std::abs(dx.imag()) < tolerance)
                    break;
            }
            return x;
        };

        switch (poly.order()) {
            case 1: {    // ===== Linear equation
                roots.emplace_back(linear(polynomial));
                break;
            }

            case 2: {    // ===== Quadratic equation
                auto quadroots = quadratic<COMPLEX_TYPE>(polynomial);
                roots.insert(roots.end(), quadroots.begin(), quadroots.end());
                break;
            }

            case 3: {    // ===== Cubic equation
                auto cubicroots = cubic<COMPLEX_TYPE>(polynomial);
                roots.insert(roots.end(), cubicroots.begin(), cubicroots.end());
                break;
            }

            default: {    // ===== Higher order equation; use Laguerre's method
                while (polynomial.order() > 3) {
                    using namespace nxx::roots;

                    // ===== Find a root using Laguerre's method
                    roots.emplace_back(laguerre(polynomial, 0.0, tolerance, max_iterations));

                    // ===== Polish the root on the original polynomial using Newton's method
                    //roots.back() = *fdfsolve(Newton(polynomial, derivativeOf(polynomial)), roots.back(), tolerance, max_iterations);
                    roots.back() = newt(original, roots.back()); //TODO: Use the roots::fdfsolve function instead

                    // ===== Deflate the polynomial by the root to reduce the order
                    polynomial /= Polynomial< COMPLEX_TYPE > { -roots.back(), 1.0 };
                }

                // ===== Solve the remaining cubic equation
                auto cuberoots = cubic<COMPLEX_TYPE>(polynomial);
                roots.insert(roots.end(), cuberoots.begin(), cuberoots.end());
                break;
            }
        }

        using RET = std::conditional_t< std::same_as< RT, void >, VALUE_TYPE, RT >;
        return impl::sortAndReturn< RET >(roots, tolerance);
    }
}    // namespace nxx::poly

#endif    // NUMERIXX_POLYROOTS_HPP

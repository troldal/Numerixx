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
#include <Constants.hpp>
#include <Roots.hpp>

// ===== Standard Library Includes
#include <algorithm>
#include <cassert>
#include <cmath>
#include <numbers>
#include <optional>
#include <random>
#include <vector>

namespace nxx::poly
{
    namespace impl
    {

        void validateTolerance(std::floating_point auto tolerance)
        {
            if (tolerance <= 0)
                throw NumerixxError("Invalid tolerance value: " + std::to_string(tolerance) + ". Tolerance must be a positive number.");
        }

        void validateMaxIterations(std::integral auto max_iterations)
        {
            if (max_iterations < 1)
                throw NumerixxError("Invalid maximum number of iterations: " + std::to_string(max_iterations) +
                                    ". Maximum number of iterations must be greater than zero.");
        }

        void validatePolynomialOrder(std::integral auto order, std::integral auto min_order)
        {
            if (order < min_order) throw NumerixxError("Polynomial must have order of at least " + std::to_string(min_order) + ".");
        }

        /**
         * @brief Sorts a vector of roots either real or complex based on their values.
         *
         * This function sorts a given vector of roots which can be either real numbers or complex numbers.
         * If the type RT is not complex, it filters out roots with an imaginary part greater or equal to
         * the square root of the given tolerance and then returns a vector of the real parts.
         * If the type RT is complex, it sorts the complex roots based on their norm and returns them.
         *
         * @tparam RT The data type of the roots to be sorted. Should be either floating-point or complex.
         * @param roots The vector of roots to sort.
         * @param tolerance The tolerance used to determine the threshold for filtering out non-real roots.
         *        It must be a positive value. This is not used for sorting complex types.
         * @return A vector of sorted roots. The type of the vector is the same as the input type RT.
         *
         * @pre The roots vector should not be empty.
         * @pre The tolerance should be a positive number.
         *
         * @throw NumerixxError if the tolerance is less than or equal to zero.
         *
         * @note The sorting is stable for complex roots and is based on the norm of the complex numbers.
         *       For real roots, the sort is based on the natural ordering of the numbers.
         */
        template< typename RT >
        requires(std::floating_point< RT > || IsComplex< RT >)
        inline auto sortRoots(auto roots, auto tolerance)
            requires IsComplex< typename decltype(roots)::value_type > &&
                     std::same_as< decltype(roots), std::vector< typename decltype(roots)::value_type > >
        {
            validateTolerance(tolerance);

            // Calculate the square root of tolerance once
            auto toleranceSqrt = std::sqrt(tolerance);

            // If the type RT is not complex, filter out roots with an imaginary part
            // greater or equal to the square root of the given tolerance
            if constexpr (!IsComplex< RT >) {
                std::erase_if(roots, [toleranceSqrt](const auto& elem) { return std::abs(elem.imag()) >= toleranceSqrt; });
            }

            // Sorting function
            auto sortingFunc = [toleranceSqrt](const auto& root1, const auto& root2) {
                return std::abs(root2.real() - root1.real()) < toleranceSqrt ? root1.imag() < root2.imag() : root1.real() < root2.real();
            };

            // Sort the roots
            std::sort(roots.begin(), roots.end(), sortingFunc);

            // If the type RT is complex, return the roots as they are. If not, return only the real parts of the roots.
            if constexpr (IsComplex< RT >) {
                return roots;
            }
            else {
                std::vector< RT > realroots;
                std::transform(roots.begin(), roots.end(), std::back_inserter(realroots), [](const auto& elem) { return elem.real(); });
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
     * @param tolerance The tolerance used to determine if the imaginary part of a complex root is
     * sufficiently small to be considered real. Only relevant if the input polynomial has complex
     * coefficients. Defaults to nxx::EPS.
     *
     * @return The root of the linear polynomial.
     *
     * @throws error::PolynomialError If the input polynomial is not a monomial.
     *
     * @note The function checks if the input polynomial is linear by verifying its order.
     * If the order is not equal to 1, an std::invalid_argument exception is thrown.
     */
    template< typename RT = void >
    inline auto linear(IsPolynomial auto poly, typename PolynomialTraits< decltype(poly) >::fundamental_type tolerance = nxx::EPS)
    {
        impl::validateTolerance(tolerance);
        impl::validatePolynomialOrder(poly.order(), 1);

        // Define type aliases for readability
        using POLY_T    = PolynomialTraits< decltype(poly) >;
        using VALUE_T   = typename POLY_T::value_type;
        using FLOAT_T   = typename POLY_T::fundamental_type;
        using COMPLEX_T = std::complex< FLOAT_T >;

        // Calculate the root of the linear polynomial
        std::vector< COMPLEX_T > root;
        root.emplace_back(-poly.coefficients().front() / poly.coefficients().back());

        // Define the return type based on the template parameter RT
        using RETURN_T = std::conditional_t< std::same_as< RT, void >, VALUE_T, RT >;

        // Sort the root and return it
        tl::expected< std::vector< RETURN_T >, NumerixxError > result = impl::sortRoots< RETURN_T >(root, tolerance);
        return result;
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
        impl::validateTolerance(tolerance);
        impl::validatePolynomialOrder(poly.order(), 2);

        // Define type aliases for readability
        using POLY_T     = PolynomialTraits< decltype(poly) >;
        using VALUE_T    = typename POLY_T::value_type;
        using FLOAT_T    = typename POLY_T::fundamental_type;
        using COMPLEX_T  = std::complex< FLOAT_T >;
        using RETURN_T   = std::conditional_t< std::same_as< RT, void >, VALUE_T, RT >;
        using EXPECTED_T = tl::expected< std::vector< RETURN_T >, NumerixxError >;

        // Extract the coefficients of the polynomial
        const auto& coeffs = poly.coefficients();
        const auto& a      = coeffs[2];
        const auto& b      = coeffs[1];
        const auto& c      = coeffs[0];

        // Calculate the discriminant
        const COMPLEX_T discriminant   = sqrt(b * b - 4.0 * a * c);
        const COMPLEX_T sqrt_component = std::conj(b) * discriminant;

        // Calculate the roots of the quadratic polynomial
        const COMPLEX_T q = -0.5 * (b + (sqrt_component.real() >= 0.0 ? discriminant : -discriminant));

        // Check if the discriminant or the coefficient 'a' is less than the tolerance
        if (std::abs(q) < tolerance || std::abs(a) < tolerance)
            return EXPECTED_T(tl::unexpected(NumerixxError("Quadratic polynomial is ill formed.")));

        // Calculate the roots
        std::vector< COMPLEX_T > roots = { q / a, c / q };

        // Sort the roots and return them
        EXPECTED_T result = impl::sortRoots< RETURN_T >(roots, tolerance);
        return result;
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
        impl::validateTolerance(tolerance);
        impl::validatePolynomialOrder(poly.order(), 3);

        using POLY_T     = PolynomialTraits< decltype(poly) >;
        using VALUE_T    = typename POLY_T::value_type;
        using FLOAT_T    = typename POLY_T::fundamental_type;
        using COMPLEX_T  = std::complex< FLOAT_T >;
        using RETURN_T   = std::conditional_t< std::same_as< RT, void >, VALUE_T, RT >;
        using EXPECTED_T = tl::expected< std::vector< RETURN_T >, NumerixxError >;

        using std::sqrt;
        using namespace std::numbers;
        using namespace std::complex_literals;

        // ===== Ad hoc lambda function to calculate the cube root of a complex number.
        auto cbrt = [](COMPLEX_T x) { return std::pow(x, 1.0 / 3.0); };

        auto coeff = poly.coefficients();
        std::transform(coeff.cbegin(), coeff.cend(), coeff.begin(), [&coeff](auto elem) { return elem / coeff.back(); });

        const auto& a = coeff[2];
        const auto& b = coeff[1];
        const auto& c = coeff[0];

        const COMPLEX_T Q = (a * a - 3.0 * b) / 9.0;
        const COMPLEX_T R = (2.0 * a * a * a - 9.0 * a * b + 27.0 * c) / 54.0;
        const COMPLEX_T A =
            -cbrt(R + ((std::conj(R) * sqrt(R * R - Q * Q * Q)).real() >= 0.0 ? sqrt(R * R - Q * Q * Q) : -sqrt(R * R - Q * Q * Q)));
        const COMPLEX_T B = (abs(A) == 0.0 ? 0.0 : Q / A);

        std::vector< COMPLEX_T > roots = { A + B - a / 3.0,
                                           -0.5 * (A + B) - a / 3.0 + 0.5 * sqrt(3.0) * (A - B) * 1.0i,
                                           -0.5 * (A + B) - a / 3.0 - 0.5 * sqrt(3.0) * (A - B) * 1.0i };

        return EXPECTED_T(impl::sortRoots< RETURN_T >(roots, tolerance));
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
                         std::complex< typename PolynomialTraits< decltype(poly) >::fundamental_type > guess          = 0.0,
                         typename PolynomialTraits< decltype(poly) >::fundamental_type                 tolerance      = nxx::EPS,
                         int                                                                           max_iterations = nxx::MAXITER)
    {
        impl::validateTolerance(tolerance);
        impl::validateMaxIterations(max_iterations);
        impl::validatePolynomialOrder(poly.order(), 4);

        // Define type aliases for readability
        using POLY_T     = PolynomialTraits< decltype(poly) >;
        using FLOAT_T    = typename POLY_T::fundamental_type;
        using COMPLEX_T  = std::complex< FLOAT_T >;
        using EXPECTED_T = tl::expected< std::vector< COMPLEX_T >, NumerixxError >;
        using OPTIONAL_T = std::optional< COMPLEX_T >;

        const COMPLEX_T& order = poly.order();

        // Define a lambda function for computing the Laguerre step.
        auto laguerrestep = [&](COMPLEX_T g_param, COMPLEX_T h_param) -> OPTIONAL_T {
            const COMPLEX_T arg = std::sqrt((order - COMPLEX_T(1)) * (order * h_param - g_param * g_param));
            const COMPLEX_T den = (abs(g_param + arg) > abs(g_param - arg) ? (g_param + arg) : (g_param - arg));
            return (abs(den) < nxx::EPS ? OPTIONAL_T(std::nullopt) : OPTIONAL_T(order / den));
        };

        // A lambda for computing the root of the polynomial using Newton's method:
        auto newt = [=](auto f, COMPLEX_T x) -> OPTIONAL_T {
            auto df = derivativeOf(f);    // Compute the derivative of the polynomial.

            int i = 0;
            while (true) {
                if (i >= max_iterations) return std::nullopt;
                if (abs(df(x)) < nxx::EPS) return std::nullopt;

                auto dx = f(x) / df(x);    // Calculate the step.
                x -= dx;                   // Update the root estimate.

                // Check for convergence.
                auto fval = f(x);
                if (std::abs(fval.real()) < tolerance && std::abs(fval.imag()) < tolerance && std::abs(dx.real()) < tolerance &&
                    std::abs(dx.imag()) < tolerance)
                    break;

                ++i;
            }
            return OPTIONAL_T(x);
        };

        // Initialize the root with the guess
        COMPLEX_T  root = guess;
        COMPLEX_T  G;
        COMPLEX_T  H;
        OPTIONAL_T step;

        // Get function objects for the first and second derivatives of the polynomial
        auto d1poly = derivativeOf(poly);
        auto d2poly = derivativeOf(d1poly);

        // Initialize a random number generator to perturb the step size every 10 iterations.
        std::random_device                        rd;
        std::mt19937                              mt(rd());
        std::uniform_real_distribution< FLOAT_T > dist(0.0, 1.0);

        // Perform the Laguerre iterations.
        int i = 0;
        while (true) {
            // If the absolute value of the polynomial evaluated at the root is less than the tolerance, return the root.
            if (abs(poly(root)) < tolerance) break;

            // Return an error if the maximum number of iterations is reached.
            if (i >= max_iterations) return EXPECTED_T(tl::unexpected(NumerixxError("Maximum number of iterations reached.")));

            // Calculate G and H for the Laguerre step
            G    = d1poly(root) / poly(root);
            H    = G * G - d2poly(root) / poly(root);
            step = laguerrestep(G, H);

            // If the step is invalid, use a small value.
            if (!step) step = OPTIONAL_T(root * 0.1);

            // If the step is below the tolerance, stop.
            if (abs(*step) < tolerance) break;

            // Perturb the step size every 10 iterations.
            if (i % 10 == 0) *step = dist(mt);

            // Update the root
            root = root - *step;

            // Increment the iteration counter
            ++i;
        }

        // ===== Polish the root on the original polynomial using Newton's method
        // roots.back() = *fdfsolve(Newton(polynomial, derivativeOf(polynomial)), roots.back(), tolerance, max_iterations);
        const auto polished_root = newt(poly, root);
        if (polished_root) root = *polished_root;    // TODO: Use the roots::fdfsolve function instead

        // Return the root
        return EXPECTED_T(std::vector { root });
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
    // Template function to solve polynomial equations of various orders.
    template< typename RT = void >
    inline auto polysolve(IsPolynomial auto                                             poly,
                          typename PolynomialTraits< decltype(poly) >::fundamental_type tolerance      = nxx::EPS,
                          int                                                           max_iterations = nxx::MAXITER)
    {
        // Validate input parameters: tolerance, max_iterations, and polynomial order.
        impl::validateTolerance(tolerance);
        impl::validateMaxIterations(max_iterations);
        impl::validatePolynomialOrder(poly.order(), 1);

        // Define types for readability and flexibility.
        using POLY_T     = PolynomialTraits< decltype(poly) >;                             // Traits of the input polynomial.
        using VALUE_T    = typename POLY_T::value_type;                                    // Value type of the polynomial.
        using FLOAT_T    = typename POLY_T::fundamental_type;                              // Fundamental type, used for complex numbers.
        using COMPLEX_T  = std::complex< FLOAT_T >;                                        // Complex type based on the fundamental type.
        using RETURN_T   = std::conditional_t< std::same_as< RT, void >, VALUE_T, RT >;    // Return type.
        using EXPECTED_T = tl::expected< std::vector< RETURN_T >, NumerixxError >;         // Expected return type.

        // Convert input polynomial to complex type and initialize a vector for roots.
        auto polynomial = Polynomial< COMPLEX_T >(std::vector< COMPLEX_T > { poly.begin(), poly.end() });
        auto roots      = std::vector< COMPLEX_T > {};

        // Lambda function to find roots based on the order of the polynomial.
        auto findRoots = [&](Polynomial< COMPLEX_T > p) {
            switch (p.order()) {
                case 1:
                    return linear< COMPLEX_T >(p);    // Linear polynomial.
                case 2:
                    return quadratic< COMPLEX_T >(p);    // Quadratic polynomial.
                case 3:
                    return cubic< COMPLEX_T >(p);    // Cubic polynomial.
                default:
                    return laguerre(p, 1.0, tolerance, max_iterations);    // Higher-order polynomial.
            }
        };

        // Loop to solve and deflate the polynomial until its order is reduced to 3 or less.
        int order;
        do {
            order                  = polynomial.order();       // Current order of the polynomial.
            const auto roots_found = findRoots(polynomial);    // Find roots.

            // Check if roots are found, return an error if not.
            if (!roots_found) [[unlikely]]
                return EXPECTED_T(tl::unexpected(NumerixxError("Error: Root-finding failed.")));

            // Insert found roots into the roots vector.
            roots.insert(roots.end(), (*roots_found).begin(), (*roots_found).end());

            // Deflate the polynomial if its order is greater than 3.
            if (polynomial.order() > 3) polynomial /= Polynomial< COMPLEX_T > { -roots.back(), 1.0 };
        }
        while (order > 3);    // Continue until the polynomial is of order 3 or less.

        // Sort the roots and return them as the expected return type.
        return EXPECTED_T(impl::sortRoots< RETURN_T >(roots, tolerance));
    }

}    // namespace nxx::poly

#endif    // NUMERIXX_POLYROOTS_HPP

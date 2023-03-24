/*
    o.     O O       o Oo      oO o.OOoOoo `OooOOo.  ooOoOOo o      O o      O
    Oo     o o       O O O    o o  O        o     `o    O     O    o   O    o
    O O    O O       o o  o  O  O  o        O      O    o      o  O     o  O
    O  o   o o       o O   Oo   O  ooOO     o     .O    O       oO       oO
    O   o  O o       O O        o  O        OOooOO'     o       Oo       Oo
    o    O O O       O o        O  o        o    o      O      o  o     o  o
    o     Oo `o     Oo o        O  O        O     O     O     O    O   O    O
    O     `o  `OoooO'O O        o ooOooOoO  O      o ooOOoOo O      o O      o

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

#ifndef NUMERIXX_POLYNOMIAL_HPP
#define NUMERIXX_POLYNOMIAL_HPP

#include "../utils/Concepts.hpp"

#include <algorithm>
#include <complex>
#include <iterator>
#include <numeric>
#include <sstream>
#include <type_traits>
#include <vector>

namespace nxx::poly
{

    template< typename CONTAINER >
    concept IsCoefficientContainer = (std::floating_point< typename CONTAINER::value_type > ||
                                      utils::IsComplex< typename CONTAINER::value_type >) &&
                                     (std::bidirectional_iterator< typename CONTAINER::iterator >);

    /**
     * @brief The Polynomial class implements the functionality of a polynomial with arbitrary number of coefficients.
     * @tparam T The type of the polynomial coefficients (must be of floating-point type)
     */
    template< typename T >
        requires std::floating_point< T > || utils::IsComplex< T >
    class Polynomial final
    {
        std::vector< T > m_coefficients; /**< The internal store of polynomial coefficients. */

    public:
        /*
         * Alias declarations.
         */
        using value_type = T;

        /**
         * @brief Constructor taking a container of coefficients as an argument.
         * @param coefficients The container of coefficients. The value type must be of floating-point type, and must
         * support iterators. The coefficients must be in order of increasing power, starting with the constant (power of zero).
         * @note Trailing zeros will be trimmed away.
         */
        explicit Polynomial(const IsCoefficientContainer auto& coefficients)
            : m_coefficients { coefficients.cbegin(),
                               std::find_if(coefficients.crbegin(), coefficients.crend(), [](T val) { return val != 0.0; }).base() }
        {}

        /**
         * @brief Constructor taking an std::initializer_list with the polynomial coefficients as an argument.
         * @param coefficients The list of coefficients of the polynomial. The coefficients must be in order of increasing power,
         * starting with the constant (power of zero).
         * @note This constructor calls the previous constructor, by converting the std::initializer_list to a std::vector.
         */
        Polynomial(std::initializer_list< T > coefficients) : Polynomial(std::vector< T >(coefficients)) {}

        /**
         * @brief Function call operator, for evaluating the polynomial at a given value.
         * @param value The floating point value at which to evaluate the polynomial.
         * @return A floating point value with the result of the evaluation.
         */
        inline auto operator()(value_type value) const { return evaluate(value); }

        /**
         * @brief Function call operator, for evaluating the polynomial at a given value.
         * @param value The std::complex value at which to evaluate the polynomial.
         * @return A std::complex value with the result of the evaluation.
         */
        inline auto operator()(std::complex< value_type > value) const { return evaluate(value); }

        /**
         * @brief Function for evaluating the polynomial at a given value.
         * @param value The floating point value at which to evaluate the polynomial.
         * @return A floating point value with the result of the evaluation.
         */
        [[nodiscard]] inline auto evaluate(value_type value) const
        {
            // ===== Horner's method implemented in terms of std::accummulate.
            return std::accumulate(m_coefficients.crbegin() + 1,
                                   m_coefficients.crend(),
                                   m_coefficients.back(),
                                   [value](value_type curr, value_type coeff) { return curr * value + coeff; });
        }

        /**
         * @brief Function for evaluating the polynomial at a given value.
         * @param value The std::complex value at which to evaluate the polynomial.
         * @return A std::complex value with the result of the evaluation.
         */
        [[nodiscard]] inline auto evaluate(std::complex< value_type > value) const
        {
            // ===== Horner's method implemented in terms of std::accummulate.
            return std::accumulate(m_coefficients.crbegin() + 1,
                                   m_coefficients.crend(),
                                   std::complex< value_type >(m_coefficients.back()),
                                   [value](std::complex< value_type > curr, value_type coeff) { return curr * value + coeff; });
        }

        /**
         * @brief Get the vector of polynomial coefficients.
         * @return A const reference to the vector of coefficients.
         */
        [[nodiscard]] const auto& coefficients() const { return m_coefficients; }

        /**
         * @brief Get the polynomial coefficients in an arbitrary container.
         * @tparam CONTAINER The type of the container. The value type must be floating type, and it must support iterators.
         * @return A CONTAINER object with the polynomial coefficients.
         */
        template< typename CONTAINER >
            requires IsCoefficientContainer< CONTAINER >
        [[nodiscard]] CONTAINER coefficients() const
        {
            return CONTAINER { m_coefficients.cbegin(), m_coefficients.cend() };
        }

        /**
         * @brief Create a string representation of the polynomial.
         * @return A std::string with the polynomial in text format.
         */
        [[nodiscard]] std::string asString() const
        {
            std::stringstream ss;

            ss << m_coefficients.front();

            int power = 1;
            std::for_each(m_coefficients.begin() + 1, m_coefficients.end(), [&](T value) {
                if (value == 0.0) {
                    ++power;
                    return;
                }

                if (value < 0.0)
                    ss << " - ";
                else
                    ss << " + ";

                ss << std::abs(value);
                ss << "x";

                if (power >= 2) ss << "^" << power;

                ++power;
            });

            return ss.str();
        }

        /**
         * @brief Get the order of the polynomial, i.e. the highest highest power of the polynomial.
         * @return The power of the polynomial as an int.
         */
        [[nodiscard]] auto order() const { return m_coefficients.size() - 1; }

        /**
         * @brief
         * @param rhs
         * @return
         */
        auto& operator+=(Polynomial< T > const& rhs)
        {
            auto temp = *this + rhs;
            std::swap(*this, temp);
        }

        /**
         * @brief
         * @param rhs
         * @return
         */
        auto& operator-=(Polynomial< T > const& rhs)
        {
            auto temp = *this - rhs;
            std::swap(*this, temp);
        }
    };

    /*
     * Deduction guide for using arbitrary containers for coefficient input.
     */
    Polynomial(IsCoefficientContainer auto coefficients) -> Polynomial< typename decltype(coefficients)::value_type >;

    /**
     * @brief A concept that checks if the given type is a specialization of the Polynomial class template.
     * This concept verifies if the given type is a specialization of the Polynomial class template with the same
     * value type as the type itself.
     * @tparam POLY The type to be checked if it is a specialization of the Polynomial class template.
     * @return std::true_type if the type is a specialization of the Polynomial class template; otherwise std::false_type.
     */
    template<typename POLY>
    concept IsPolynomial = std::same_as< POLY, Polynomial< typename POLY::value_type> >;

    /**
     * @brief Create a function object representing the derivative of the input function.
     * @tparam T The exact type of the polynomial coefficients (will be auto-deduced from the argument)
     * @param func The function object to get the derivative of.
     * @return A Polynomial object representing the derivative of the input function.
     * @note This function is similar to the .derivativeOf() function in the nxx::deriv namespace. This
     * function works only on polynomials, where the derivative will always be another polynomial.
     */
    template< typename T >
    inline auto derivativeOf(const Polynomial< T >& func)
    {
        using value_type = T;
        std::vector< value_type > coefficients { func.coefficients().cbegin() + 1, func.coefficients().cend() };

        int n = 1;
        std::transform(coefficients.cbegin(), coefficients.cend(), coefficients.begin(), [&n](value_type elem) { return elem * n++; });
        return Polynomial(coefficients);
    }

    /**
     * @brief Adds two polynomials.
     * @tparam T Type of coefficients of the first polynomial.
     * @tparam U Type of coefficients of the second polynomial.
     * @param lhs The first polynomial to add.
     * @param rhs The second polynomial to add.
     * @return The polynomial resulting from adding rhs to lhs.
     */
    template< typename T, typename U >
    auto operator+(Polynomial< T > const& lhs, Polynomial< U > const& rhs)
    {
        using TYPE = std::common_type_t< T, U >;

        if (lhs.order() < rhs.order()) {
            auto coeffs = rhs.coefficients();
            std::transform(coeffs.cbegin(), coeffs.cend(), lhs.coefficients().begin(), coeffs.begin(), std::plus< TYPE >());
            return Polynomial(coeffs);
        }
        else {
            auto coeffs = lhs.coefficients();
            std::transform(coeffs.cbegin(), coeffs.cend(), rhs.coefficients().begin(), coeffs.begin(), std::plus< TYPE >());
            return Polynomial(coeffs);
        }
    }

    /**
     * @brief Subtracts two polynomials.
     * @tparam T Type of coefficients of the first polynomial.
     * @tparam U Type of coefficients of the second polynomial.
     * @param lhs The first polynomial to subtract.
     * @param rhs The second polynomial to subtract.
     * @return The polynomial resulting from subtracting rhs from lhs.
     */
    template< typename T, typename U >
    auto operator-(Polynomial< T > const& lhs, Polynomial< U > const& rhs)
    {
        using TYPE = std::common_type_t< T, U >;

        if (lhs.order() < rhs.order()) {
            auto coeffs = rhs.coefficients();
            std::transform(coeffs.cbegin(), coeffs.cend(), lhs.coefficients().begin(), coeffs.begin(), [](TYPE a, TYPE b) { return b - a; });
            return Polynomial(coeffs);
        }
        else {
            auto coeffs = lhs.coefficients();
            std::transform(coeffs.cbegin(), coeffs.cend(), rhs.coefficients().begin(), coeffs.begin(), [](TYPE a, TYPE b) { return a - b; });
            return Polynomial(coeffs);
        }
    }

    /**
     * @brief Multiply two polynomials using the `*` operator.
     * This operator multiplies the polynomial `lhs` by the polynomial `rhs` and returns the result as a new polynomial.
     * @tparam T The type of coefficients in the polynomial `lhs`.
     * @tparam U The type of coefficients in the polynomial `rhs`.
     * @param lhs The first polynomial to multiply.
     * @param rhs The second polynomial to multiply.
     * @return The product of `lhs` and `rhs`.
     * @note The function assumes that the coefficient type `T` is convertible to the coefficient type `U`, or vice versa.
     * @see Polynomial
     */
    template< typename T, typename U >
    auto operator*(Polynomial< T > const& lhs, Polynomial< U > const& rhs)
    {
        using TYPE = std::common_type_t< T, U >;

        std::vector< TYPE > result(lhs.order() + rhs.order() + 1, 0.0);

        int counter = 0;
        for (auto elem : rhs.coefficients()) {
            std::vector< TYPE > subtractor(lhs.coefficients().cbegin(), lhs.coefficients().cend());
            std::transform(subtractor.cbegin(), subtractor.cend(), subtractor.begin(), [=](TYPE val) { return val * elem; });
            std::transform(subtractor.cbegin(), subtractor.cend(), result.begin() + counter, result.begin() + counter, std::plus< TYPE >());
            ++counter;
        }

        return Polynomial(result);
    }

    /**
     *   @brief Divide one polynomial by another.
     *   This function divides the polynomial lhs by the polynomial rhs and returns
     *   the quotient and remainder as polynomials.
     *   @tparam T The type of coefficients in the polynomial lhs.
     *   @tparam U The type of coefficients in the polynomial rhs.
     *   @param lhs The polynomial to be divided.
     *   @param rhs The polynomial to divide by.
     *   @return A pair of polynomials (quotient and remainder).
     *   @note The function assumes that the coefficient type T is convertible to the coefficient type U, or vice versa.
     *   @note The degree of the remainder polynomial is less than the degree of the divisor polynomial.
     *   @see Polynomial
     */
    template< typename T, typename U >
    auto divide(Polynomial< T > const& lhs, Polynomial< U > const& rhs)
    {
        using TYPE = std::common_type_t< T, U >;

        const auto& lhs_coeffs = lhs.coefficients();
        const auto& rhs_coeffs = rhs.coefficients();

        std::vector< TYPE > quotient(lhs.order() - rhs.order() + 1, 0.0);
        std::vector< TYPE > remainder(lhs_coeffs.begin(), lhs_coeffs.end());

        // ===== Consider simplifying this loop =====
        auto odiff  = lhs.order() - rhs.order();
        for (auto q_iter = quotient.rbegin(), r_iter = remainder.rbegin(); q_iter != quotient.rend(); ++q_iter, ++r_iter) {
            *q_iter = *r_iter / rhs_coeffs.back();

            auto counter = lhs.order() - 1;
            for (auto iter = r_iter + 1; iter != r_iter + odiff + 1; ++iter) {
                *iter -= *q_iter * rhs_coeffs[(counter--) - odiff];
            }
        }

        remainder.erase(remainder.begin() + rhs.order(), remainder.end());
        return std::make_pair(Polynomial(quotient), Polynomial(remainder));
    }

    /**
     * @brief Divide one polynomial by another using the `/` operator.
     * This operator divides the polynomial `lhs` by the polynomial `rhs` using the `divide` function
     * and returns the quotient as a polynomial.
     * @param lhs The polynomial to be divided.
     * @param rhs The polynomial to divide by.
     * @return The quotient of dividing `lhs` by `rhs`.
     * @see Polynomial
     * @see divide
     */
    auto operator/(const IsPolynomial auto& lhs, const IsPolynomial auto& rhs) { return divide(lhs, rhs).first; }

    /**
     * @brief Compute the remainder of dividing one polynomial by another using the % operator.
     * This operator divides the polynomial lhs by the polynomial rhs using the divide function
     * and returns the remainder as a polynomial.
     * @param lhs The polynomial to be divided.
     * @param rhs The polynomial to divide by.
     * @return The remainder of dividing `lhs` by `rhs`.
     * @see Polynomial
     * @see divide
     */
    auto operator%(const IsPolynomial auto& lhs, const IsPolynomial auto& rhs) { return divide(lhs, rhs).second; }

}    // namespace nxx::poly

#endif    // NUMERIXX_POLYNOMIAL_HPP

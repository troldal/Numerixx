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

#ifndef NUMERIXX_POLYNOMIAL_HPP
#define NUMERIXX_POLYNOMIAL_HPP

// ===== Numerixx Includes
#include "PolynomialError.hpp"

// ===== External Includes
#include <tl/expected.hpp>
#include "../.utils/Concepts.hpp"

// ===== Standard Library Includes
#include <algorithm>
#include <complex>
#include <iterator>
#include <numeric>
#include <sstream>
#include <type_traits>
#include <vector>

namespace nxx::poly
{

    /**
     * @brief A concept that checks whether a container is suitable for storing polynomial coefficients.
     *
     * This concept checks whether a container has a value type that is a floating-point type or a
     * complex type, and whether its iterator type is bidirectional.
     *
     * @tparam CONTAINER The container to check.
     */
    template< typename CONTAINER >
    concept IsCoefficientContainer = (std::floating_point< typename CONTAINER::value_type > ||
                                      IsComplex< typename CONTAINER::value_type >) &&
                                     (std::bidirectional_iterator< typename CONTAINER::iterator >);

    /*
     * Forward declaration of the Polynomial class.
     */
    template< typename T >
        requires std::floating_point< T > || IsComplex< T >
    class Polynomial;

    /*
     * Forward declaration of the PolynomialTraits class.
     */
    template<typename...>
    struct PolynomialTraits;

    /*
     * Specialization of the PolynomialTraits class for Polynomial objects with floating point coefficients.
     */
    template<typename T>
//        requires std::floating_point< T >
    struct PolynomialTraits<Polynomial<T>>
    {
        using value_type = T;
        using fundamental_type = T;
    };

    /*
     * Specialization of the PolynomialTraits class for Polynomial objects with complex coefficients.
     */
    template<typename T>
    struct PolynomialTraits<Polynomial<std::complex<T>>>
    {
        using value_type = std::complex<T>;
        using fundamental_type = T;
    };

    /**
     * @brief A class representing a polynomial with coefficients of type T.
     *
     * This class provides functionality to evaluate a polynomial at a point,
     * retrieve the polynomial's coefficients, get the polynomial's order, and
     * output the polynomial as a string. Additionally, it provides operators
     * for addition and subtraction of polynomials.
     *
     * @tparam T The type of the polynomial coefficients. This must be a floating
     * point type or a type that satisfies the `utils::IsComplex` concept.
     */
    template< typename T >
        requires std::floating_point< T > || IsComplex< T >
    class Polynomial final
    {
        std::vector< T > m_coefficients; /**< The internal store of polynomial coefficients. */

    public:
        /**
         * @brief The type of the polynomial coefficients.
         */
        using value_type = T;

        /**
         * @brief Constructs a polynomial from a range of coefficients.
         *
         * The range must be a container whose value type is convertible to T. The
         * constructed polynomial will have a degree equal to the number of
         * non-zero coefficients minus one.
         *
         * @tparam IsCoefficientContainer A type trait that verifies that the
         * template argument is a container of coefficients.
         * @param coefficients The range of coefficients.
         */
        explicit Polynomial(const IsCoefficientContainer auto& coefficients)
            : m_coefficients { coefficients.cbegin(),
                               std::find_if(coefficients.crbegin(), coefficients.crend(), [](T val) { return val != 0.0; }).base() } {
                if (m_coefficients.empty()) {
                    m_coefficients.push_back(0.0);
                }
        }

        /**
         * @brief Constructs a polynomial from a list of coefficients.
         *
         * The constructed polynomial will have a degree equal to the number of
         * non-zero coefficients minus one.
         *
         * @param coefficients The list of coefficients.
         */
        Polynomial(std::initializer_list< T > coefficients) : Polynomial(std::vector< T >(coefficients)) {}

        /**
         * @brief Evaluates the polynomial at a given value.
         *
         * This function evaluates the polynomial at a given value using Horner's method and returns the result.
         * The function template parameter `U` specifies the type of the input value, which can be either a real or
         * complex number.
         *
         * @param value The value to evaluate the polynomial at.
         * @return The value of the polynomial at the specified input value.
         */
        inline auto operator()(auto value) const { return *evaluate(value); }

        /**
         * @brief Evaluates the polynomial at a given value.
         *
         * This function evaluates the polynomial at a given value using Horner's method and returns the result.
         * The function template parameter `U` specifies the type of the input value, which can be either a real or
         * complex number.
         *
         * @tparam U The type of the input value.
         * @param value The value to evaluate the polynomial at.
         * @return The value of the polynomial at the specified input value.
         */
        template< typename U >
            requires std::floating_point< U > || IsComplex< U >
        [[nodiscard]] inline auto evaluate(U value) const
        {
                using TYPE = std::common_type_t< T, U >;
                tl::expected< TYPE, error::PolynomialError > result;

                // ===== Horner's method implemented in terms of std::accummulate.
                TYPE eval = std::accumulate(m_coefficients.crbegin() + 1,
                                            m_coefficients.crend(),
                                            static_cast< TYPE >(m_coefficients.back()),
                                            [value](TYPE curr, TYPE coeff) { return curr * value + coeff; });

                if (std::isfinite(abs(eval)))
                    result = eval;
                else
                    result = tl::make_unexpected(error::PolynomialError { "Computation of polynomial gave non-finite result." });
                return result;
        }

        /**
         * @brief Gets the coefficients of the polynomial.
         *
         * @return A constant reference to the vector of coefficients.
         */
        [[nodiscard]] const auto& coefficients() const { return m_coefficients; }

        /**
         * @brief Gets the coefficients of the polynomial as a container of a different type.
         *
         * The container must be constructible from a range of values convertible to T.
         *
         * @tparam CONTAINER The type of the container to use. Must be a container of coefficients.
         * @param coefficients The container in which to store the coefficients.
         * @return The container of coefficients.
         */
        template< typename CONTAINER >
            requires IsCoefficientContainer< CONTAINER >
        [[nodiscard]] CONTAINER coefficients() const
        {
            return CONTAINER { m_coefficients.cbegin(), m_coefficients.cend() };
        }

        /**
         * @brief Returns a string representation of the polynomial.
         *
         * The polynomial is output in the form `a_0 + a_1 x + a_2 x^2 + ... + a_{n-1} x^{n-1} + a_n x^n`.
         * If a coefficient is zero, it is omitted. If a coefficient is negative, it is
         * shown with a negative sign.
         *
         * @return A string representation of the polynomial.
         */
        [[nodiscard]] std::string asString() const
        {
            std::ostringstream oss;

            // Print highest degree term with sign, if non-zero
            if (!m_coefficients.empty()) {
                oss << m_coefficients.front();
                for (size_t i = 1; i < m_coefficients.size(); ++i) {
                    auto coeff = m_coefficients[i];
                    if constexpr (!IsComplex< T > ) {
                        if (coeff == 0.0) { continue; }
                        else if (coeff > 0.0) { oss << " + "; }
                        else { oss << " - "; }
                        oss << abs(coeff) << "x";
                    }
                    else {
                        if (abs(coeff) == 0.0) { continue; }
                        else { oss << " + "; }
                        oss << coeff << "x";
                    }
                    if (i >= 2) { oss << "^" << i; }
                }
            }
            return oss.str();
        }

        /**
         * @brief Returns the order of the polynomial.
         *
         * The order of the polynomial is one less than the number of coefficients.
         *
         * @return The order of the polynomial.
         */
        [[nodiscard]] auto order() const { return m_coefficients.size() - 1; }

        /**
         * @brief Adds another polynomial to this polynomial.
         *
         * This operator adds the given polynomial to this polynomial and returns
         * the result as a new polynomial. The degree of the result will be equal
         * to the maximum degree of the two polynomials.
         *
         * @param rhs The polynomial to add to this polynomial.
         * @return The sum of the two polynomials.
         */
        auto& operator+=(Polynomial< T > const& rhs)
        {
            auto temp = *this + rhs;
            std::swap(*this, temp);
            return *this;
        }

        /**
         * @brief Subtracts another polynomial from this polynomial.
         *
         * This operator subtracts the given polynomial from this polynomial and returns
         * the result as a new polynomial. The degree of the result will be equal
         * to the maximum degree of the two polynomials.
         *
         * @param rhs The polynomial to subtract from this polynomial.
         * @return The difference between the two polynomials.
         */
        auto& operator-=(Polynomial< T > const& rhs)
        {
            auto temp = *this - rhs;
            std::swap(*this, temp);
            return *this;
        }

        /**
         * @brief Multiplies the polynomial by another polynomial and returns the result.
         *
         * This operator multiplies the polynomial by another polynomial and assigns the resulting polynomial
         * to the current polynomial object. The function uses the free operator*() function to perform
         * the multiplication.
         *
         * @param rhs The other polynomial to multiply by.
         * @return A reference to the modified polynomial object.
         */
        auto& operator*=(Polynomial< T > const& rhs)
        {
            auto temp = *this * rhs;
            std::swap(*this, temp);
            return *this;
        }

        /**
         * @brief Divides the polynomial by another polynomial and returns the result.
         *
         * This operator divides the polynomial by another polynomial and assigns the resulting polynomial
         * to the current polynomial object. The function uses the operator/() function to perform the division.
         *
         * @param rhs The other polynomial to divide by.
         * @return A reference to the modified polynomial object.
         */
        auto& operator/=(Polynomial< T > const& rhs)
        {
            auto temp = *this / rhs;
            std::swap(*this, temp);
            return *this;
        }

        /**
         * @brief Returns a const iterator to the beginning of the coefficients container of the Polynomial.
         *
         * This member function provides a const iterator pointing to the first coefficient of the Polynomial.
         * It can be used for traversing the coefficients in a read-only manner.
         *
         * @return A const iterator pointing to the first coefficient of the Polynomial.
         *
         * @note The returned iterator is compatible with standard C++ library algorithms and functions that
         * expect const_iterators as input.
         */
        auto begin() const { return m_coefficients.cbegin(); }
        auto cbegin() const { return m_coefficients.cbegin(); }


        /**
         * @brief Returns a const iterator to the end of the coefficients container of the Polynomial.
         *
         * This member function provides a const iterator pointing to one past the last coefficient of the
         * Polynomial. It can be used for traversing the coefficients in a read-only manner.
         *
         * @return A const iterator pointing to one past the last coefficient of the Polynomial.
         *
         * @note The returned iterator is compatible with standard C++ library algorithms and functions that
         * expect const_iterators as input.
         */
        auto end() const { return m_coefficients.cend(); }
        auto cend() const { return m_coefficients.cend(); }
    };

    /*
     * Deduction guide for using arbitrary containers for coefficient input.
     */
    template < typename CONTAINER >
        requires IsCoefficientContainer< CONTAINER >
    Polynomial(CONTAINER coefficients) -> Polynomial< typename decltype(coefficients)::value_type >;

    /**
     * @brief A concept that checks if the given type is a Polynomial of some type T.
     *
     * This concept checks if the given type is a Polynomial of some type T. It requires
     * that the value_type of the Polynomial be the same as the given type.
     *
     * @tparam POLY The type to check.
     */
    template< typename POLY >
    concept IsPolynomial = std::same_as< POLY, Polynomial< typename POLY::value_type > >;

    /**
     * @brief Computes the derivative of a polynomial.
     *
     * This function computes the derivative of the given polynomial and returns the result
     * as a new polynomial. The derivative of a polynomial f(x) of degree n is given by
     * f'(x) = a_1 + 2 a_2 x + ... + n a_n x^{n-1}.
     *
     * @param func The polynomial to differentiate.
     * @return The derivative of the polynomial.
     */
    inline auto derivativeOf(IsPolynomial auto func)
    {
        if (func.order() == 0) throw error::PolynomialError("Cannot differentiate a constant polynomial.");

        using VALUE_TYPE   = typename PolynomialTraits< decltype(func) >::value_type;
        using FLOAT_TYPE   = typename PolynomialTraits< decltype(func) >::fundamental_type;

        std::vector< VALUE_TYPE > coefficients { func.coefficients().cbegin() + 1, func.coefficients().cend() };

        int n = 1;
        std::transform(coefficients.cbegin(), coefficients.cend(), coefficients.begin(), [&n](VALUE_TYPE elem) { return elem * static_cast<FLOAT_TYPE>(n++); });
        return Polynomial(coefficients);
    }

    /**
     * @brief Adds two polynomials.
     *
     * This operator adds the given two polynomials and returns the result as a new polynomial.
     * The degree of the result will be equal to the maximum degree of the two polynomials.
     *
     * @tparam T The type of the coefficients of the first polynomial.
     * @tparam U The type of the coefficients of the second polynomial.
     * @param lhs The first polynomial to add.
     * @param rhs The second polynomial to add.
     * @return The sum of the two polynomials.
     */
    template< typename T, typename U >
    auto operator+(Polynomial< T > const& lhs, Polynomial< U > const& rhs)
    {
        using TYPE = std::common_type_t< T, U >;

        if (lhs.order() < rhs.order()) {
            auto coeffs = rhs.coefficients();
            std::transform(lhs.coefficients().cbegin(), lhs.coefficients().cend(), coeffs.cbegin(), coeffs.begin(), std::plus< TYPE >());
            return Polynomial(coeffs);
        }
        else {
            auto coeffs = lhs.coefficients();
            std::transform(rhs.coefficients().cbegin(), rhs.coefficients().cend(), coeffs.cbegin(), coeffs.begin(), std::plus< TYPE >());
            return Polynomial(coeffs);
        }
    }

    /**
     * @brief Subtracts two polynomials.
     *
     * This operator subtracts the second polynomial from the first polynomial and returns
     * the result as a new polynomial. The degree of the result will be equal to the maximum
     * degree of the two polynomials.
     *
     * @tparam T The type of the coefficients of the first polynomial.
     * @tparam U The type of the coefficients of the second polynomial.
     * @param lhs The polynomial to subtract from.
     * @param rhs The polynomial to subtract.
     * @return The difference between the two polynomials.
     */
    template< typename T, typename U >
    auto operator-(Polynomial< T > const& lhs, Polynomial< U > const& rhs)
    {
        using TYPE = std::common_type_t< T, U >;

        if (lhs.order() < rhs.order()) {
            auto coeffs = rhs.coefficients();
            std::transform(lhs.coefficients().cbegin(), lhs.coefficients().cend(), coeffs.cbegin(), coeffs.begin(), [](TYPE a, TYPE b) {
                return b - a;
            });
            return Polynomial(coeffs);
        }
        else {
            auto coeffs = lhs.coefficients();
            std::transform(rhs.coefficients().cbegin(), rhs.coefficients().cend(), coeffs.cbegin(), coeffs.begin(), [](TYPE a, TYPE b) {
                return a - b;
            });
            return Polynomial(coeffs);
        }
    }

    /**
     * @brief Multiplies two polynomials.
     *
     * This operator multiplies the given two polynomials and returns the result as a new
     * polynomial. The degree of the result will be equal to the sum of the degrees of the
     * two polynomials.
     *
     * @tparam T The type of the coefficients of the first polynomial.
     * @tparam U The type of the coefficients of the second polynomial.
     * @param lhs The first polynomial to multiply.
     * @param rhs The second polynomial to multiply.
     * @return The product of the two polynomials.
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
     * @brief Divides two polynomials and returns the quotient and remainder.
     *
     * This function divides the first polynomial by the second polynomial and returns the
     * quotient and remainder as a pair of polynomials. The degree of the remainder will be
     * less than the degree of the second polynomial.
     *
     * @tparam T The type of the coefficients of the first polynomial.
     * @tparam U The type of the coefficients of the second polynomial.
     * @param lhs The polynomial to divide.
     * @param rhs The polynomial to divide by.
     * @return The quotient and remainder of the two polynomials.
     */
    template< typename T, typename U >
    auto divide(Polynomial< T > const& lhs, Polynomial< U > const& rhs)
    {
        using TYPE = std::common_type_t< T, U >;

        std::vector<TYPE> dividend = lhs.coefficients();
        std::vector<TYPE> divisor = rhs.coefficients();
        std::vector<TYPE> remainder {};

        if (divisor.empty() || divisor.back() == 0.0 || rhs.order() > lhs.order())
            throw error::PolynomialError("Invalid divisor polynomial");

        std::vector<TYPE> quotient(lhs.order() - rhs.order() + 1, 0);
        remainder = dividend;

        for (std::size_t i = lhs.order(); i >= rhs.order() && i < dividend.size(); --i) {
            TYPE coef = remainder[i] / divisor.back();
            quotient[i - rhs.order()] = coef;

            for (std::size_t j = 0; j <= rhs.order(); ++j) {
                remainder[i - j] -= coef * divisor[rhs.order() - j];
            }
        }

        remainder.erase(std::find_if(remainder.crbegin(), remainder.crend(), [](TYPE val) { return val != 0.0; }).base() + 1, remainder.end());

        return std::make_pair(Polynomial(quotient), Polynomial(remainder));
    }

    /**
     * @brief Divides two polynomials and returns the quotient.
     *
     * This operator divides the first polynomial by the second polynomial and returns the
     * quotient as a new polynomial. The degree of the quotient will be less than or equal to
     * the degree of the first polynomial minus the degree of the second polynomial.
     *
     * @param lhs The polynomial to divide.
     * @param rhs The polynomial to divide by.
     * @return The quotient of the two polynomials.
     */
    auto operator/(const IsPolynomial auto& lhs, const IsPolynomial auto& rhs) { return divide(lhs, rhs).first; }

    /**
     * @brief Divides two polynomials and returns the remainder.
     *
     * This operator divides the first polynomial by the second polynomial and returns the
     * remainder as a new polynomial. The degree of the remainder will be less than the
     * degree of the second polynomial.
     *
     * @param lhs The polynomial to divide.
     * @param rhs The polynomial to divide by.
     * @return The remainder of the two polynomials.
     */
    auto operator%(const IsPolynomial auto& lhs, const IsPolynomial auto& rhs) { return divide(lhs, rhs).second; }

}    // namespace nxx::poly

#endif    // NUMERIXX_POLYNOMIAL_HPP

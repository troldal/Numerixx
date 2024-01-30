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

#pragma once

// ===== Numerixx Includes
#include <Concepts.hpp>
#include <Error.hpp>

// ===== External Includes
#include "_external.hpp"

// ===== Standard Library Includes
#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <iterator>
#include <numeric>
#include <optional>
#include <sstream>
#include <type_traits>
#include <vector>

namespace nxx::poly
{

    namespace detail
    {

        template< typename T >
        struct PolyErrorData
        {
            std::string        details;
            std::vector< T >   coefficients;
            std::optional< T > arg    = std::nullopt;
            std::optional< T > result = std::nullopt;

            friend std::ostream& operator<<(std::ostream& os, const PolyErrorData& data)
            {
                os << data.details << "\n";
                os << "Coefficients: ";
                for (size_t i = 0; i < data.coefficients.size(); ++i) {
                    os << data.coefficients[i];
                    if (i < data.coefficients.size() - 1) {
                        os << ", ";
                    }
                }
                if (data.arg) {
                    os << "\nArgument: " << *data.arg;
                }
                if (data.result) {
                    os << "\nResult: " << *data.result;
                }
                return os;
            }
        };

    }    // namespace detail

    /**
     * @concept IsCoefficientContainer
     * @brief Checks if a container type is suitable for storing polynomial coefficients.
     *
     * This concept checks if a container type can be used to store the coefficients of a polynomial.
     * It enforces that the container must have bidirectional iterators and its value type must be
     * either a floating point or a complex number. Additionally, it must provide begin and end
     * member functions that return forward iterators, ensuring the container can be iterated over
     * in a range-based for loop.
     *
     * @tparam CONTAINER The container type to be checked against the concept.
     *
     * @requirements
     * - CONTAINER must provide begin and end functions that return forward iterators.
     * - CONTAINER::value_type must be a floating point or a complex number type.
     * - CONTAINER::iterator must be a bidirectional iterator.
     */
    template< typename CONTAINER >
    concept IsCoefficientContainer =
        requires(CONTAINER a) {
            {
                begin(a)
            } -> std::forward_iterator;    // Check if 'begin' gives a forward iterator
            {
                end(a)
            } -> std::forward_iterator;    // Check if 'end' gives a forward iterator
        } && (nxx::IsFloat< typename CONTAINER::value_type > ||
              IsComplex< typename CONTAINER::value_type >)&&(std::bidirectional_iterator< typename CONTAINER::iterator >);

    /*
     * Forward declaration of the PolynomialTraits class.
     */
    template< typename... >
    struct PolynomialTraits;

    /*
     * Specialization of the PolynomialTraits class for Polynomial objects with floating point coefficients.
     */
    template< typename T >
        requires nxx::IsFloat< T >
    struct PolynomialTraits< Polynomial< T > >
    {
        using value_type       = T;
        using fundamental_type = T;
    };

    /*
     * Specialization of the PolynomialTraits class for Polynomial objects with complex coefficients.
     */
    template<typename T>
        requires nxx::IsFloat< T >
    struct PolynomialTraits< Polynomial< std::complex< T > > >
    {
        using value_type       = std::complex< T >;
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
    template<typename T = double>
        requires nxx::IsFloat< T > || IsComplex< T >
    class Polynomial final
    {
        std::vector< T > m_coefficients; /**< The internal store of polynomial coefficients. */

    public:
        /**
         * @brief The type of the polynomial coefficients.
         */
        using value_type = T;

        /**
         * @brief Constructs a zero polynomial.
         *
         * This constructor initializes the polynomial with a single zero coefficient.
         */
        Polynomial()
            : Polynomial({ 0 })
        {}

        /**
         * @brief Constructs a Polynomial with the given coefficients and an optional custom serializer.
         *
         * This constructor initializes a Polynomial with the given coefficients, trimming any trailing zeros.
         * It also sets a custom serializer function for representing the Polynomial as a string, if provided.
         * The coefficients are processed in reverse to remove any trailing zeros, ensuring that the
         * polynomial is stored in its minimal form. If no non-zero coefficients are provided, a zero polynomial
         * is created with a single zero coefficient.
         *
         * @tparam IsCoefficientContainer Checks if the input is a container with a valid
         *                                value type for the coefficients.
         * @param coefficients An container of coefficients representing the polynomial.
         *                     The coefficients are expected to be provided in increasing order
         *                     of degree, i.e., {a0, a1, a2, ...}, where a0 is the constant term.
         *
         * @throws NumerixxError if the serializer function is not provided.
         */
        explicit Polynomial(const IsCoefficientContainer auto& coefficients)
        {
            // Lambda function to determine if a coefficient is near zero, considering both floating-point and complex numbers.
            auto is_near_zero = [&](typename std::remove_cvref_t< decltype(coefficients) >::value_type val) -> bool {
                // Use a different epsilon value based on whether the type is complex or not.
                if constexpr (IsComplex< decltype(val) >) {
                    // Calculate epsilon for complex numbers.
                    auto epsilon = 1e-8;//std::numeric_limits< typename decltype(val)::value_type >::epsilon();
                    // Check if the norm of the complex number is within the tolerance defined by epsilon.
                    return std::norm(val) <= sqrt(epsilon);
                }
                else {
                    // Calculate epsilon for floating-point numbers.
                    auto epsilon = std::numeric_limits< decltype(val) >::epsilon();
                    // Check if the value is within the tolerance defined by epsilon.
                    return std::norm(val) <= sqrt(epsilon);
                }
            };

            // Find the iterator to the first non-zero coefficient when traversing the container in reverse.
            auto rev_it = std::find_if_not(coefficients.crbegin(), coefficients.crend(), is_near_zero);

            // If all coefficients are near zero, initialize the polynomial as a zero polynomial.
            if (rev_it == coefficients.crend()) {
                m_coefficients.push_back(T {});    // Ensure at least one coefficient for a zero polynomial
            }
            else {
                // Assign coefficients from the start to the first non-zero coefficient found, excluding trailing zeros.
                m_coefficients.assign(coefficients.cbegin(), rev_it.base());
            }
        }

        /**
         * @brief Constructs a polynomial from an initializer list of coefficients.
         *
         * This constructor allows for creating a Polynomial object using brace-enclosed
         * initializer list for coefficients. It will automatically remove any trailing zeros
         * from the polynomial representation. A serializer function can be provided to
         * customize the serialization of the polynomial. If no serializer is provided, a
         * default one is used which serializes the coefficients in a human-readable string format.
         *
         * @param coefficients An initializer list of coefficients representing the polynomial.
         *                     The coefficients are expected to be provided in increasing order
         *                     of degree, i.e., {a0, a1, a2, ...}, where a0 is the constant term.
         *
         * @note If the initializer list is empty, the constructed polynomial will be the zero
         *       polynomial. Also, if the serializer function provided is empty, an exception will
         *       be thrown.
         *
         * @throw NumerixxError if the serializer function is empty.
         */
        Polynomial(std::initializer_list< T > coefficients)
            : Polynomial(std::vector< T >(coefficients))
        {}

        /**
         * @brief Constructs a polynomial with a single coefficient.
         *
         * This constructor is explicit to prevent unintentional type conversions.
         *
         * @param coefficient The single coefficient of the polynomial.
         */
        explicit Polynomial(T coefficient)
            : Polynomial({ coefficient })
        {}

        /**
         * @brief Constructs a polynomial from a range of coefficients.
         *
         * This constructor template allows the creation of a polynomial from iterators
         * pointing to the beginning and end of a range containing coefficients.
         *
         * @tparam InputIterator The type of the input iterators. Must satisfy the requirements
         * of InputIterator.
         * @param first The beginning iterator of the coefficient range.
         * @param last The ending iterator of the coefficient range (one-past-the-end).
         */
        template< typename InputIterator >
        Polynomial(InputIterator first, InputIterator last)
            : Polynomial(std::vector< typename std::iterator_traits< InputIterator >::value_type > { first, last })
        {}

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
         * @brief Evaluates the polynomial at a given point using Horner's method.
         *
         * This function evaluates the polynomial for a given value of `U`. It checks for
         * the polynomial's order and handles edge cases, such as when the polynomial
         * order is zero or there are no coefficients. It also includes error handling for
         * non-finite arguments or results. Horner's method is used for the evaluation
         * which is efficient and numerically stable.
         *
         * @tparam U The type of the value at which the polynomial is evaluated. It must
         *           be convertible to `T` or be a floating-point or complex type.
         * @param value The point at which to evaluate the polynomial.
         *
         * @return tl::expected<std::common_type_t<T, U>, nxx::NumerixxError>
         *         The result of the polynomial evaluation. It returns a `tl::expected` object
         *         that contains the evaluation result or an error. The evaluation result is of
         *         type `std::common_type_t<T, U>`, which is the common type between `T` and `U`.
         *         If an error occurs during the evaluation, `tl::unexpected` is returned with
         *         an instance of `nxx::NumerixxError` containing the error details.
         *
         * @note This function uses `std::accumulate` for its implementation. The
         *       accumulation starts from the second to last coefficient to the first,
         *       ensuring that the iterator does not surpass the bounds.
         *
         * @exception nxx::NumerixxError Thrown if the polynomial has no coefficients, or if
         *            the argument or result of the evaluation is non-finite.
         */
        template<typename U>
            requires std::convertible_to< U, T > || nxx::IsFloat< U > || IsComplex< U >
        [[nodiscard]]
        inline auto evaluate(U value) const
            -> tl::expected< std::common_type_t< T, U >, Error< detail::PolyErrorData< std::common_type_t< T, U > > > >
        {
            using TYPE      = std::common_type_t< T, U >;
            using PolyError = Error< detail::PolyErrorData< TYPE > >;

            if (order() == 0) [[unlikely]]
                return m_coefficients.front();

            auto coeff_begin = m_coefficients.crbegin();
            auto coeff_end   = m_coefficients.crend();

            if (m_coefficients.size() <= 1 || coeff_begin == coeff_end) [[unlikely]]

                return tl::unexpected(PolyError("Polynomial error",
                                                nxx::NumerixxErrorType::Poly,
                                                { .details      = "Polynomial evaluation failed; no coefficients.",
                                                  .coefficients = { m_coefficients.begin(), m_coefficients.end() } }));

            // ===== Horner's method implemented in terms of std::accummulate.
            TYPE result =
                std::accumulate(coeff_begin + 1, coeff_end, static_cast< TYPE >(m_coefficients.back()), [value](TYPE curr, TYPE coeff) {
                    return curr * static_cast< TYPE >(value) + coeff;
            });

            if (!std::isfinite(std::abs(result))) [[unlikely]]
                return tl::unexpected(PolyError("Polynomial error",
                                                nxx::NumerixxErrorType::Poly,
                                                { .details      = "Polynomial evaluation failed; non-finite result.",
                                                  .coefficients = { m_coefficients.begin(), m_coefficients.end() },
                                                  .arg          = value,
                                                  .result       = result }));

            return result;
        }

        /**
         * @brief Gets the coefficients of the polynomial.
         *
         * @return A constant reference to the vector of coefficients.
         */
        [[nodiscard]]
        const std::vector< T >& coefficients() const
        {
            return m_coefficients;
        }

        /**
         * @brief Gets the coefficients of the polynomial as a container of a different type.
         *
         * The container must be constructible from a range of values convertible to T.
         *
         * @tparam CONTAINER The type of the container to use. Must be a container of coefficients.
         * @return The container of coefficients.
         */
        template< typename CONTAINER >
        requires IsCoefficientContainer< CONTAINER >
        [[nodiscard]]
        CONTAINER coefficients() const
        {
            return CONTAINER { m_coefficients.cbegin(), m_coefficients.cend() };
        }

        /**
         * @brief Outputs the polynomial to an output stream.
         *
         * The polynomial is output in the form `a_0 + a_1 x + a_2 x^2 + ... + a_{n-1} x^{n-1} + a_n x^n`.
         * If a coefficient is zero, it is omitted. If a coefficient is negative, it is
         * shown with a negative sign.
         *
         * @return A string representation of the polynomial.
         */
        friend std::ostream& operator<<(std::ostream& os, Polynomial const& p)
        {
            if (p.m_coefficients.empty()) {
                return os << 0;    // output zero if the polynomial is empty
            }

            // Iterator for the coefficients
            auto coeff_it = p.m_coefficients.cbegin();
            // Handle the constant term (0th degree) separately
            os << *coeff_it++;

            // Now handle the rest of the coefficients
            int degree = 1;
            for (; coeff_it != p.m_coefficients.cend(); ++coeff_it, ++degree) {
                // We use 'auto' to work with different types (e.g., float, double, complex)
                const auto& coeff = *coeff_it;

                // Skip zero coefficients
                if (std::abs(coeff) < std::abs(std::sqrt(std::numeric_limits< T >::epsilon()))) continue;

                // Use different logic for non-complex and complex types
                if constexpr (!IsComplex< T >) {
                    os << (coeff > T {} ? " + " : " - ") << std::abs(coeff);
                }
                else {
                    // Complex numbers might need a different approach for sign and abs
                    os << " + " << coeff;
                }

                // For non-zero coefficients, print 'x' and maybe the exponent
                os << "x";
                if (degree > 1) {
                    os << "^" << degree;
                }
            }

            return os;

            // Old code:
            //            std::ostringstream oss;
            //
            //            // Print highest degree term with sign, if non-zero
            //            if (!p.m_coefficients.empty()) {
            //                oss << p.m_coefficients.front();
            //                for (size_t i = 1; i < p.m_coefficients.size(); ++i) {
            //                    auto coeff = p.m_coefficients[i];
            //                    if constexpr (!IsComplex< T >) {
            //                        if (coeff == 0.0) {
            //                            continue;
            //                        }
            //                        else if (coeff > 0.0) {
            //                            oss << " + ";
            //                        }
            //                        else {
            //                            oss << " - ";
            //                        }
            //                        oss << abs(coeff) << "x";
            //                    }
            //                    else {
            //                        if (abs(coeff) == 0.0) {
            //                            continue;
            //                        }
            //                        else {
            //                            oss << " + ";
            //                        }
            //                        oss << coeff << "x";
            //                    }
            //                    if (i >= 2) {
            //                        oss << "^" << i;
            //                    }
            //                }
            //            }
            //
            //            os << oss.str();
            //            return os;
        }

        /**
         * @brief Returns the order of the polynomial.
         *
         * The order of the polynomial is one less than the number of coefficients.
         *
         * @return The order of the polynomial.
         */
        [[nodiscard]]
        auto order() const
        {
            return m_coefficients.size() - 1;
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
        //        [[nodiscard]]
        //        std::string serialize() const
        //        {
        //            return m_serializer(m_coefficients);
        //        }

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
        template<typename U>
            requires nxx::IsFloat< U > || (IsComplex< T > && IsComplex< U >)
        Polynomial< T >& operator+=(Polynomial< U > const& rhs)
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
        template<typename U>
            requires nxx::IsFloat< U > || (IsComplex< T > && IsComplex< U >)
        Polynomial< T >& operator-=(Polynomial< U > const& rhs)
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
        template<typename U>
            requires nxx::IsFloat< U > || (IsComplex< T > && IsComplex< U >)
        Polynomial< T >& operator*=(Polynomial< U > const& rhs)
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
        template<typename U>
            requires nxx::IsFloat< U > || (IsComplex< T > && IsComplex< U >)
        Polynomial< T >& operator/=(Polynomial< U > const& rhs)
        {
            auto temp = *this / rhs;
            std::swap(*this, temp);
            return *this;
        }

        /**
         * @brief Equality operator for polynomials.
         *
         * This operator compares two polynomials for equality. Two polynomials are equal if all their
         * coefficients are equal.
         *
         * @param rhs The second polynomial to compare.
         * @return True if the two polynomials are equal, false otherwise.
         */
        template<typename U>
            requires(nxx::IsFloat< T > && nxx::IsFloat< U >) || (IsComplex< T > && IsComplex< U >)
        bool operator==(Polynomial< U > const& rhs) const
        {
            return m_coefficients == rhs.m_coefficients;
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
    template< typename CONTAINER >
    requires IsCoefficientContainer< CONTAINER >
    Polynomial(CONTAINER coefficients) -> Polynomial< typename CONTAINER::value_type >;

    template< typename CONTAINER, typename FUNC >
    requires IsCoefficientContainer< CONTAINER >
    Polynomial(CONTAINER coefficients, FUNC f) -> Polynomial< typename CONTAINER::value_type >;

    template<typename T>
        requires nxx::IsFloat< T > || IsComplex< T >
    Polynomial(std::initializer_list< T > coefficients) -> Polynomial< T >;

    template<typename T, typename FUNC>
        requires nxx::IsFloat< T > || IsComplex< T >
    Polynomial(std::initializer_list< T > coefficients, FUNC f) -> Polynomial< T >;

    /**
     * @brief Creates a polynomial from a given set of roots.
     *
     * Constructs a polynomial whose roots are specified in the input container. This function
     * template is constrained to accept only containers that fulfill the IsCoefficientContainer
     * concept, ensuring the container supports operations required to iterate over its elements.
     *
     * @tparam CONTAINER A container type that satisfies the IsCoefficientContainer concept.
     *                   The container should hold elements of a numeric type.
     * @param roots A container holding the roots of the polynomial.
     * @return A Polynomial object whose coefficients are calculated to have the specified roots.
     *
     * @note The function assumes the underlying numeric type of the container can form a polynomial,
     *       which includes arithmetic operations and value initialization.
     */
    template< typename CONTAINER >
    requires IsCoefficientContainer< CONTAINER >
    Polynomial< typename CONTAINER::value_type > createPolynomialFromRoots(const CONTAINER& roots)
    {
        using ValueType             = typename CONTAINER::value_type;
        using CoefficientsContainer = std::vector< ValueType >;

        // Start with a constant polynomial p(x) = 1
        CoefficientsContainer coefficients = { ValueType { 1 } };

        // Iterate over each root to construct the polynomial
        for (const auto& root : roots) {
            // Initialize a new vector for the new coefficients with one extra degree
            CoefficientsContainer newCoefficients(coefficients.size() + 1, ValueType { 0 });

            // Compute the new coefficients
            for (size_t i = 0; i < coefficients.size(); ++i) {
                // No out-of-bounds access should occur here, but we'll add an assert as a sanity check
                assert(i < newCoefficients.size() - 1);    // Ensure we do not go out of bounds

                newCoefficients[i + 1] += coefficients[i];       // x * coefficient
                newCoefficients[i] -= root * coefficients[i];    // -root * coefficient
            }

            // Update the coefficients for the next iteration
            coefficients.swap(newCoefficients);
        }

        // Construct and return the Polynomial with the final coefficients
        return Polynomial< ValueType >(coefficients);
    }

    /**
     * @brief Creates a polynomial with the given roots specified in an initializer list.
     *
     * This overload allows for convenient creation of a polynomial from an initializer list.
     * It internally converts the initializer list to a vector and calls the vector-based
     * createPolynomialFromRoots function template.
     *
     * @tparam T The numeric type of the polynomial coefficients (e.g., double, std::complex<double>).
     * @param roots An initializer list containing the roots of the polynomial.
     * @return A Polynomial object with coefficients calculated based on the provided roots.
     */
    template<typename T>
        requires nxx::IsFloat< T > || IsComplex< T >
    Polynomial< T > createPolynomialFromRoots(const std::initializer_list< T >& roots)
    {
        // Convert initializer_list to a vector
        std::vector< T > rootsVector(roots.begin(), roots.end());

        // Use the previously defined function template
        return createPolynomialFromRoots(rootsVector);
    }

    /**
     * @brief A concept that checks if the given type is a Polynomial of some type T.
     *
     * This concept checks if the given type is a Polynomial of some type T. It requires
     * that the value_type of the Polynomial be the same as the given type.
     *
     * @tparam POLY The type to check.
     */
    //    template< typename POLY >
    //    concept IsPolynomial = std::same_as< POLY, Polynomial< typename POLY::value_type > >;

    /**
     * @brief Computes the derivative of a given polynomial function.
     *
     * This function takes as an input a polynomial (assumed by the "IsPolynomial auto" parameter type),
     * and computes its derivative. This is done using the power rule of differentiation,
     * which states that d/dx[aX^n] = a*n*X^(n-1), where a is a coefficient and n is the power of each term in the polynomial.
     * If the input function is a constant (i.e., it has an order of 0), it throws an error as constants do not have meaningful derivatives.
     *
     * @param func The polynomial function to compute the derivative of.
     * @return A polynomial function, which is the derivative of the input function.
     * @throws error::PolynomialError if the input polynomial is constant.
     */
    inline auto derivativeOf(IsPolynomial auto func)
    {
        // Throw an error if the order of the polynomial function is 0 (i.e., the function is constant).
        if (func.order() == 0) throw std::runtime_error("Cannot differentiate a constant polynomial.");

        using VALUE_TYPE =
            typename PolynomialTraits< decltype(func) >::value_type;    // This is the type of the coefficients of the polynomial.
        using FLOAT_TYPE = typename PolynomialTraits< decltype(func) >::fundamental_type;    // This is the type of the exponents and other
                                                                                             // arithmetic operations in the polynomial.

        // Get a list of the coefficients of the input polynomial, skipping the first one (which corresponds to the constant term).
        std::vector< VALUE_TYPE > coefficients { func.coefficients().cbegin() + 1, func.coefficients().cend() };

        int n = 1;
        // Compute the derivative of each term in the polynomial using the power rule, and store them back in the 'coefficients' vector.
        std::transform(coefficients.cbegin(), coefficients.cend(), coefficients.begin(), [&n](VALUE_TYPE elem) {
            return elem * static_cast< FLOAT_TYPE >(n++);    // d/dx[aX^n] = a*n*X^(n-1)
        });

        // Return a new polynomial that represents the derivative of the input function.
        return Polynomial(coefficients);
    }

    /**
     * @brief Converts a Polynomial object to its string representation.
     *
     * This function takes a Polynomial of a generic type and converts it into a
     * string by streaming it into a std::stringstream. It assumes that the
     * Polynomial class has an overload of the stream insertion operator (`operator<<`)
     * defined for std::ostream that dictates how the Polynomial is converted to a string.
     *
     * @tparam T The type of the coefficients in the Polynomial.
     * @param p The Polynomial object to convert to a string.
     * @return A std::string representing the Polynomial.
     */
    template< typename T >
    std::string to_string(Polynomial< T > const& p)
    {
        std::stringstream ss;
        ss << p;
        return ss.str();
    }

    /**
     * @brief Adds two polynomials.
     *
     * This operator adds the given two polynomials and returns the result as a new polynomial.
     * The degree of the result will be equal to the maximum degree of the two polynomials.
     *
     * @tparam T The type of the first operand. Can be any type.
     * @tparam U The type of the second operand. Can be any type.
     * @param lhs The first operand, a Polynomial.
     * @param rhs The second operand, a Polynomial.
     *
     * @returns An object of type Polynomial that represents the sum of lhs and rhs.
     */
    template<typename T, typename U>
        requires(nxx::IsFloat< T > || IsComplex< T >) && (nxx::IsFloat< U > || IsComplex< U >)
    auto operator+(Polynomial< T > const& lhs, Polynomial< U > const& rhs)
    {
        // Determine the common type between T and U
        using TYPE = std::common_type_t< T, U >;

        // Check if lhs's degree is less than rhs's degree
        if (lhs.order() < rhs.order()) {
            // Copy rhs's coefficients to add to
            std::vector< TYPE > coeffs(rhs.coefficients().cbegin(), rhs.coefficients().cend());

            // Add lhs's coefficients to copied rhs's coefficients one by one
            std::transform(lhs.coefficients().cbegin(), lhs.coefficients().cend(), coeffs.cbegin(), coeffs.begin(), std::plus< TYPE >());

            // Return a Polynomial constructed from the resulting coefficients
            return Polynomial< TYPE >(coeffs);
        }
        else {
            // Copy lhs's coefficients to add to
            std::vector< TYPE > coeffs(lhs.coefficients().cbegin(), lhs.coefficients().cend());

            // Add rhs's coefficients to copied lhs's coefficients one by one
            std::transform(rhs.coefficients().cbegin(), rhs.coefficients().cend(), coeffs.cbegin(), coeffs.begin(), std::plus< TYPE >());

            // Return a Polynomial constructed from the resulting coefficients
            return Polynomial< TYPE >(coeffs);
        }
    }

    /**
     * @brief Subtracts two polynomials.
     *
     * This operator subtracts the second polynomial from the first polynomial and returns
     * the result as a new polynomial. The degree of the result will be equal to the maximum
     * degree of the two polynomials.
     *
     * @tparam T The type of the first operand. Can be any type.
     * @tparam U The type of the second operand. Can be any type.
     * @param lhs The first operand, a Polynomial.
     * @param rhs The second operand, a Polynomial.
     *
     * @returns An object of type Polynomial that represents the difference of lhs and rhs.
     */
    template<typename T, typename U>
        requires(nxx::IsFloat< T > || IsComplex< T >) && (nxx::IsFloat< U > || IsComplex< U >)
    auto operator-(Polynomial< T > const& lhs, Polynomial< U > const& rhs)
    {
        // Determine the common type between T and U
        using TYPE = std::common_type_t< T, U >;

        // Check if lhs's degree is less than rhs's degree
        if (lhs.order() < rhs.order()) {
            // Copy rhs's coefficients to subtract from
            std::vector< TYPE > coeffs(rhs.coefficients().cbegin(), rhs.coefficients().cend());

            // Subtract lhs's coefficients from copied rhs's coefficients one by one
            std::transform(lhs.coefficients().cbegin(), lhs.coefficients().cend(), coeffs.cbegin(), coeffs.begin(), [](TYPE a, TYPE b) {
                return a - b;
            });

            // Return a Polynomial constructed from the resulting coefficients
            return Polynomial(coeffs);
        }
        else {
            // Copy lhs's coefficients to subtract from
            std::vector< TYPE > coeffs(lhs.coefficients().cbegin(), lhs.coefficients().cend());

            // Subtract rhs's coefficients from copied lhs's coefficients one by one
            std::transform(rhs.coefficients().cbegin(), rhs.coefficients().cend(), coeffs.cbegin(), coeffs.begin(), [](TYPE a, TYPE b) {
                return b - a;
            });

            // Return a Polynomial constructed from the resulting coefficients
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
     * @tparam T The type of the first operand. Can be any type.
     * @tparam U The type of the second operand. Can be any type.
     * @param lhs The first operand, a Polynomial.
     * @param rhs The second operand, a Polynomial.
     *
     * @returns An object of type Polynomial that represents the product of lhs and rhs.
     */
    template<typename T, typename U>
        requires(nxx::IsFloat< T > || IsComplex< T >) && (nxx::IsFloat< U > || IsComplex< U >)
    auto operator*(Polynomial< T > const& lhs, Polynomial< U > const& rhs)
    {
        // Determine the common type between T and U
        using TYPE = std::common_type_t< T, U >;

        // Initialize a structure to contain the result
        // The degree of the product is the sum of degrees of multipliers (plus one for the constant term)
        std::vector< TYPE > result(lhs.order() + rhs.order() + 1, 0.0);

        // Counter for the number of terms already processed
        int counter = 0;

        // Iterate over every term in the RHS polynomial
        for (auto elem : rhs.coefficients()) {
            // Get a copy of the LHS coefficients
            std::vector< TYPE > subtractor(lhs.coefficients().cbegin(), lhs.coefficients().cend());

            // Multiply each term in LHS by the current term in RHS
            std::transform(subtractor.cbegin(), subtractor.cend(), subtractor.begin(), [=](TYPE val) { return val * elem; });

            // Add the resulting terms to the result polynomial, taking into consideration their respective powers
            // (represented by the 'counter' variable)
            std::transform(subtractor.cbegin(), subtractor.cend(), result.begin() + counter, result.begin() + counter, std::plus< TYPE >());

            // Advance the counter
            ++counter;
        }

        // Return a Polynomial constructed from the resulting coefficients
        return Polynomial(result);
    }

    /**
     * @brief Divides two polynomials and returns the quotient and remainder.
     *
     * This function divides the first polynomial by the second polynomial and returns the
     * quotient and remainder as a pair of polynomials. The degree of the remainder will be
     * less than the degree of the second polynomial.
     *
     * @tparam T The type of the first operand. Can be any type.
     * @tparam U The type of the second operand. Can be any type.
     * @param lhs The first operand, the dividend (polynomial to be divided).
     * @param rhs The second operand, the divisor (polynomial to divide by).
     *
     * @return a pair of polynomials for quotient and remainder respectively.
     *
     * @throws error::PolynomialError when invalid divisor polynomial is passed i.e it either handles singleton case when
     * the divisor doesn't have coefficients or the polynomial order of divisor is larger than the dividend.
     */
    template<typename T, typename U>
        requires(nxx::IsFloat< T > || IsComplex< T >) && (nxx::IsFloat< U > || IsComplex< U >)
    auto divide(Polynomial< T > const& lhs, Polynomial< U > const& rhs)
    {
        // Determine the type of polynomial which is common between T and U
        using TYPE = std::common_type_t< T, U >;

        // The coefficients of the polynomials
        std::vector< TYPE > dividend(lhs.cbegin(), lhs.cend());
        std::vector< TYPE > divisor(rhs.cbegin(), rhs.cend());

        // Handle case when divisor doesn't have coefficients or
        // when polynomial order of divisor is larger than the dividend
        if (divisor.empty() || divisor.back() == 0.0 || rhs.order() > lhs.order())
            throw NumerixxError("Divisor polynomial cannot be empty or have a higher degree than the dividend.");

        // Initialize the quotient polynomial coefficients with 0
        assert(lhs.order() - rhs.order() + 1 > 0);    // Ensure we do not go out of bounds
        std::vector< TYPE > quotient(lhs.order() - rhs.order() + 1, TYPE { 0 });
        std::vector< TYPE > remainder = dividend;

        // Iteratively compute the coefficients of the quotient polynomial
        for (std::size_t i = lhs.order(); i >= rhs.order() && i < dividend.size(); --i) {
            assert(i < remainder.size());                 // Ensure we do not go out of bounds
            assert(i - rhs.order() < quotient.size());    // Ensure we do not go out of bounds
            TYPE coef                 = remainder[i] / divisor.back();
            quotient[i - rhs.order()] = coef;

            // For each coefficient in the divisor subtract coef*divisor from the dividend (remainder)
            for (std::size_t j = 0; j <= rhs.order(); ++j) {
                assert(i - j < remainder.size());            // Ensure we do not go out of bounds
                assert(rhs.order() - j < divisor.size());    // Ensure we do not go out of bounds
                remainder[i - j] -= coef * divisor[rhs.order() - j];
            }
        }

        // Trim the leading zeros in the remainder
        remainder.erase(std::find_if(remainder.crbegin(), remainder.crend(), [](TYPE val) { return val != 0.0; }).base() + 1,
                        remainder.end());

        // Return the quotient and the remainder as a pair
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


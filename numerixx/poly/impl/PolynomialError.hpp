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

#ifndef NUMERIXX_POLYNOMIALERROR_HPP
#define NUMERIXX_POLYNOMIALERROR_HPP

//===== Standard Library Includes
#include <stdexcept>
#include <vector>

namespace nxx::error
{
    /**
     * @brief The RootErrorType enum is an enum class for the different types of root-finding errors.
     */
    enum class PolyErrorType { RootNotFound, MaxIterationsExceeded, NumericalError };

    /**
     * @brief The PolynomialError class represents an exception that is thrown when there is an error related to polynomials.
     *
     * This class is derived from the std::runtime_error class.
     */
    class PolynomialError : public std::runtime_error
    {
    public:
        // Takes a lvalue reference (std::string) and just copy it to base class
        explicit PolynomialError(const std::string& msg)
            : std::runtime_error(msg)
        {}

        // Takes a rvalue reference (std::string) and move it to base class
        explicit PolynomialError(std::string&& msg)
            : std::runtime_error(msg)
        {}
    };

    namespace impl
    { /**
       * @brief The RootErrorImpl class is a template class for root-finding errors.
       * @tparam T The type of the root value.
       */
        template< typename T >
        class PolyErrorImpl : public PolynomialError
        {
            PolyErrorType   m_type;           /*< The type of the error. */
            std::vector<T>  m_roots;          /*< The value of the root before the error was thrown. */

        public:
            /**
             * @brief Constructor.
             * @param msg The error message.
             * @param type The type of the error.
             * @param value The value of the root before the error was thrown.
             * @param iter The number of iterations performed before the error was thrown.
             */
            explicit PolyErrorImpl(const char* msg, PolyErrorType type, std::vector<T> roots)
                : PolynomialError(msg),
                  m_type(type),
                  m_roots(value)
            {}

            /**
             * @brief Returns the type of the error.
             * @return The type of the error.
             */
            [[nodiscard]]
            PolyErrorType type() const
            {
                return m_type;
            }

            /**
             * @brief Returns a string representation of the error type.
             * @return A string representation of the error type.
             */
            [[nodiscard]]
            auto typeAsString() const
            {
                switch (m_type) {
                    case PolyErrorType::RootNotFound:
                        return "Root missing; number of roots does not match degree of polynomial";
                    case PolyErrorType::MaxIterationsExceeded:
                        return "Max iterations exceeded";
                    case PolyErrorType::NumericalError:
                        return "Numerical error";
                }

                return "Unknown error";
            }

            /**
             * @brief Returns the value of the root before the error was thrown.
             * @return The value of the root.
             */
            [[nodiscard]]
            auto value() const
            {
                return m_roots;
            }
        };
    }    // namespace impl

}    // namespace nxx::error

#endif    // NUMERIXX_POLYNOMIALERROR_HPP

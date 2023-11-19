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

#ifndef NUMERIXX_ROOTCOMMON_HPP
#define NUMERIXX_ROOTCOMMON_HPP

#include <Concepts.hpp>

// ===== Standard Library Includes
#include <stdexcept>

namespace nxx::roots
{

    // ========================================================================
    // ROOT ERROR/EXCEPTION CLASSES
    // ========================================================================

    /**
     * @brief The RootErrorType enum is an enum class for the different types of root-finding errors.
     */
    enum class RootErrorType { NoRootInBracket, MaxIterationsExceeded, NumericalError };

    /**
     * @brief The RootError class is a base class for root-finding errors.
     */
    class RootError : public std::runtime_error
    {
    public:
        /**
         * @brief Constructor.
         * @param msg The error message.
         */
        explicit RootError(const char* msg)
            : std::runtime_error(msg)
        {}
    };

    namespace impl
    { /**
       * @brief The RootErrorImpl class is a template class for root-finding errors.
       * @tparam T The type of the root value.
       */
        template< typename T >
        class RootErrorImpl : public RootError
        {
            RootErrorType m_type;           /*< The type of the error. */
            T             m_value;          /*< The value of the root before the error was thrown. */
            int           m_iterations = 0; /*< The number of iterations performed before the error was thrown. */

        public:
            /**
             * @brief Constructor.
             * @param msg The error message.
             * @param type The type of the error.
             * @param value The value of the root before the error was thrown.
             * @param iter The number of iterations performed before the error was thrown.
             */
            explicit RootErrorImpl(const char* msg, RootErrorType type, T value, int iter = 0)
                : RootError(msg),
                  m_type(type),
                  m_value(value),
                  m_iterations(iter)
            {}

            /**
             * @brief Returns the type of the error.
             * @return The type of the error.
             */
            [[nodiscard]]
            RootErrorType type() const
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
                    case RootErrorType::NoRootInBracket:
                        return "No root in bracket";
                    case RootErrorType::MaxIterationsExceeded:
                        return "Max iterations exceeded";
                    case RootErrorType::NumericalError:
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
                return m_value;
            }

            /**
             * @brief Returns the number of iterations performed before the error was thrown.
             * @return The number of iterations performed before the error was thrown.
             */
            [[nodiscard]]
            auto iterations() const
            {
                return m_iterations;
            }
        };
    }    // namespace impl

    // ========================================================================
    // TRAITS CLASSES FOR ROOT-FINDING WITHOUT DERIVATIVES
    // ========================================================================

    /*
     * Forward declaration of the Ridders class.
     */
    template<typename FN, typename ARG_T>
        requires IsFloatInvocable< FN >
    class Ridder;

    /*
     * Forward declaration of the Bisection class.
     */
    template<typename FN, typename ARG_T>
        requires IsFloatInvocable< FN >
    class Bisection;

    /*
     * Forward declaration of the RegulaFalsi class.
     */
    template<typename FN, typename ARG_T>
        requires IsFloatInvocable< FN >
    class RegulaFalsi;

    /*
     * Private implementation details.
     */
    namespace impl
    {
        /*
         * Forward declaration of the BracketingTraits class.
         */
        template< typename FN >
        struct BracketingTraits;

        /*
         * Specialization of the BracketingTraits class for Ridders<FN>
         */
        template<typename FN, typename T>
        struct BracketingTraits< Ridder< FN, T > >
        {
            using FUNCTION_T = FN;
            using ARG_T = T;
            using RETURN_T = std::invoke_result_t< FN, ARG_T >;
        };

        /*
         * Specialization of the BracketingTraits class for Bisection<FN>
         */
        template<typename FN, typename T>
        struct BracketingTraits< Bisection< FN, T > >
        {
            using FUNCTION_T = FN;
            using ARG_T = T;
            using RETURN_T = std::invoke_result_t< FN, ARG_T >;
        };

        /*
         * Specialization of the BracketingTraits class for RegulaFalsi<FN>
         */
        template<typename FN, typename T>
        struct BracketingTraits< RegulaFalsi< FN, T > >
        {
            using FUNCTION_T = FN;
            using ARG_T = T;
            using RETURN_T = std::invoke_result_t< FN, ARG_T >;
        };
    } // namespace impl

    // ========================================================================
    // TRAITS CLASSES FOR ROOT-FINDING WITH DERIVATIVES
    // ========================================================================

    /*
     * Forward declaration of the DNewton class.
     */
    template< typename FN, typename DFN >
    requires IsFloatInvocable< FN > && IsFloatInvocable< DFN >
    class DNewton;

    /*
     * Forward declaration of the Newton class.
     */
    template< typename FN, typename DFN >
    requires(IsFloatInvocable< FN > && IsFloatInvocable< DFN >) ||
            (IsComplex< std::invoke_result_t< FN, double > > && IsComplex< std::invoke_result_t< DFN, double > >)
    class Newton;

    namespace impl
    {
        /*
         * Forward declaration of the PolishingTraits class.
         */
        template< typename... ARGS >
        struct PolishingTraits;

        /*
         * Specialization of the PolishingTraits class for Newton<FN, DFN>
         */
        template< typename FN, typename DFN >
        struct PolishingTraits< Newton< FN, DFN > >
        {
            using FUNCTION_T        = FN;
            using DERIV_T           = DFN;
            using FUNCTION_RETURN_T = std::invoke_result_t< FN, double >;
            using DERIV_RETURN_T    = std::invoke_result_t< DFN, double >;
        };

        /*
         * Specialization of the PolishingTraits class for DNewton<FN, DFN>
         */
        template< typename FN, typename DFN >
        struct PolishingTraits< DNewton< FN, DFN > >
        {
            using FUNCTION_T        = FN;
            using DERIV_T           = DFN;
            using FUNCTION_RETURN_T = std::invoke_result_t< FN, double >;
            using DERIV_RETURN_T    = std::invoke_result_t< DFN, double >;
        };

    }    // namespace impl

    template< typename FN >
    requires IsFloatInvocable< FN >
    class BracketSearchUp;

    template< typename FN >
    requires IsFloatInvocable< FN >
    class BracketSearchDown;

    template< typename FN >
    requires IsFloatInvocable< FN >
    class BracketExpandUp;

    template< typename FN >
    requires IsFloatInvocable< FN >
    class BracketExpandDown;

    template< typename FN >
    requires IsFloatInvocable< FN >
    class BracketExpandOut;

    template< typename FN >
    requires IsFloatInvocable< FN >
    class BracketSubdivide;

    /*
     * Private implementation details.
     */
    namespace impl
    {
        /*
         * Forward declaration of the BracketingTraits class.
         */
        template< typename FN >
        struct SearchingTraits;

        /*
         * Specialization of the BracketingTraits class for Ridders<FN>
         */
        template< typename FN >
        struct SearchingTraits< BracketSearchUp< FN > >
        {
            using FUNCTION_T = FN;
            using RETURN_T   = std::invoke_result_t< FN, double >;
        };

        template< typename FN >
        struct SearchingTraits< BracketSearchDown< FN > >
        {
            using FUNCTION_T = FN;
            using RETURN_T   = std::invoke_result_t< FN, double >;
        };

        template< typename FN >
        struct SearchingTraits< BracketExpandUp< FN > >
        {
            using FUNCTION_T = FN;
            using RETURN_T   = std::invoke_result_t< FN, double >;
        };

        template< typename FN >
        struct SearchingTraits< BracketExpandDown< FN > >
        {
            using FUNCTION_T = FN;
            using RETURN_T   = std::invoke_result_t< FN, double >;
        };

        template< typename FN >
        struct SearchingTraits< BracketExpandOut< FN > >
        {
            using FUNCTION_T = FN;
            using RETURN_T   = std::invoke_result_t< FN, double >;
        };

        template< typename FN >
        struct SearchingTraits< BracketSubdivide< FN > >
        {
            using FUNCTION_T = FN;
            using RETURN_T   = std::invoke_result_t< FN, double >;
        };

    }    // namespace impl

}    // namespace nxx::roots

#endif    // NUMERIXX_ROOTCOMMON_HPP

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

#ifndef NUMERIXX_ROOTPOLISHING_HPP
#define NUMERIXX_ROOTPOLISHING_HPP

// ===== Numerixx Includes
#include "RootCommon.hpp"
#include <Constants.hpp>
#include <Deriv.hpp>
#include <Poly.hpp>

// ===== Standard Library Includes
#include <functional>
#include <stdexcept>
#include <utility>

namespace nxx::roots
{
    using nxx::deriv::derivativeOf;
    using nxx::poly::derivativeOf;

    // ========================================================================
    // ROOT-FINDING WITH DERIVATIVES
    // ========================================================================

    /*
     * Private implementation details.
     */
    namespace impl
    {

        /**
         * @brief A base class for root polishing methods.
         *
         * This class serves as a base for various root polishing methods, providing common functionality
         * and interface. The POLICY template parameter should provide specific algorithm implementation.
         *
         * @tparam POLICY The policy class defining the specific root polishing method.
         * @requires The POLICY class should have a specialization of PolishingTraits that defines:
         *           function_type, deriv_type, function_return_type, and deriv_return_type.
         */
        template< typename POLICY >
        requires(std::floating_point< typename PolishingTraits< POLICY >::FUNCTION_RETURN_T > &&
                 std::floating_point< typename PolishingTraits< POLICY >::DERIV_RETURN_T >) ||
                (IsComplex< typename PolishingTraits< POLICY >::FUNCTION_RETURN_T > &&
                 IsComplex< typename PolishingTraits< POLICY >::DERIV_RETURN_T >)
        class PolishingBase
        {
            /*
             * Friend declarations.
             */
            friend POLICY;

        public:
            using FUNCTION_T = /**< The type of the function object. */
                typename impl::PolishingTraits< POLICY >::FUNCTION_T;
            using DERIV_T = /**< The type of the derivative function object. */
                typename impl::PolishingTraits< POLICY >::DERIV_T;
            using FUNCTION_RETURN_T = /**< The return type of the function object. */
                typename impl::PolishingTraits< POLICY >::FUNCTION_RETURN_T;
            using DERIV_RETURN_T = /**< The return type of the derivative function object. */
                typename impl::PolishingTraits< POLICY >::DERIV_RETURN_T;

        protected:
            /**
             * @brief Default constructor.
             */
            ~PolishingBase() = default;

        private:
            FUNCTION_T        m_func {};               /**< The function object to find the root for. */
            DERIV_T           m_deriv {};              /**< The function object for the derivative. */
            FUNCTION_RETURN_T m_guess;                 /**< The current root estimate. */
            bool              m_initialized { false }; /**< Flag indicating whether the object has been initialized. */

            /**
             * @brief Constructor, taking function objects for the function and its derivative as arguments.
             * @param objective The function object to find the root for.
             * @param derivative The function object for the derivative
             * @note Constructor is private to avoid direct usage by clients.
             */
            explicit PolishingBase(FUNCTION_T objective, DERIV_T derivative)
                : m_func { objective },
                  m_deriv { derivative }
            {}

        public:
            /**
             * @brief Copy constructor.
             *
             * @param other Another PolishingBase object to be copied.
             */
            PolishingBase(const PolishingBase& other) = default;

            /**
             * @brief Move constructor.
             *
             * @param other Another PolishingBase object to be moved.
             */
            PolishingBase(PolishingBase&& other) noexcept = default;

            /**
             * @brief Copy assignment operator.
             *
             * @param other Another PolishingBase object to be copied.
             * @return A reference to the assigned object.
             */
            PolishingBase& operator=(const PolishingBase& other) = default;

            /**
             * @brief Move assignment operator.
             *
             * @param other Another PolishingBase object to be moved.
             * @return A reference to the assigned object.
             */
            PolishingBase& operator=(PolishingBase&& other) noexcept = default;

            /**
             * @brief Initialize the solver by setting the initial bounds around the root.
             * @tparam T The type of the root estimate.
             * @param guess The root estimate.
             */
            template< typename T >
            requires std::convertible_to< T, FUNCTION_RETURN_T >
            void init(T guess)
            {
                m_initialized = true;
                m_guess       = guess;
            }

            /**
             * @brief Reset the solver to its initial state.
             */
            void reset() { m_initialized = false; }

            /**
             * @brief Evaluate the function object at a given point.
             * @tparam T The type of the argument (must be floating point type)
             * @param value The value at which to evaluate the function.
             * @return The result of the evaluation.
             */
            template< typename T >
            requires std::convertible_to< T, FUNCTION_RETURN_T >
            auto evaluate(T value)
            {
                return m_func(value);
            }

            /**
             * @brief Evaluate the derivative at a given point.
             * @tparam T The type of the argument (must be floating point type)
             * @param value The value at which to evaluate the derivative.
             * @return The result of the evaluation.
             */
            template< typename T >
            requires std::convertible_to< T, FUNCTION_RETURN_T >
            auto derivative(T value)
            {
                return m_deriv(value);
            }

            /**
             * @brief Get the current estimate of the root.
             * @return A const reference to the root.
             */
            auto result() const
            {
                if (!m_initialized) throw std::runtime_error("Solver has not been initialized.");
                return m_guess;
            }
        };
    }    // namespace impl

    /**
     * @brief A class implementing the Discrete Newton method for root finding.
     *
     * This class is derived from the PolishingBase class and implements the Discrete Newton method
     * for finding the roots of a given function. The function and its derivative are provided
     * as template arguments.
     *
     * @tparam FN The type of the function object for which to find the root.
     * @tparam DFN The type of the derivative function object.
     * @requires FN and DFN should be callable with double and return a floating point type.
     */
    template< typename FN, typename DFN >
        requires IsFloatInvocable< FN > && IsFloatInvocable< DFN >
    class DNewton final : public impl::PolishingBase< DNewton< FN, DFN > >
    {
        /*
         * Private alias declarations.
         */
        using BASE = impl::PolishingBase< DNewton< FN, DFN > >;

    public:
        /*
         * Public alias declarations.
         */
        using function_type = FN;
        using deriv_type    = DFN;

        /**
         * @brief Constructor, taking the function object as an argument.
         * @param objective The function object for which to find the root.
         */
        explicit DNewton(FN objective)
            : BASE(objective, derivativeOf(objective))
        {}

        /**
         * @brief Perform one iteration. This is the main algorithm of the DNewton method.
         */
        void iterate() { BASE::m_guess = BASE::m_guess - BASE::evaluate(BASE::m_guess) / BASE::derivative(BASE::m_guess); }
    };

    /**
     * @brief Deduction guide for the DNewton algorithm
     */
    template< typename FN >
    DNewton(FN objective) -> DNewton< FN, std::function< decltype(objective(0.0))(decltype(objective(0.0))) > >;

    /**
     * @brief A class implementing the Newton-Raphson method for root finding.
     *
     * This class is derived from the PolishingBase class and implements the Newton-Raphson method
     * for finding the roots of a given function. The function and its derivative are provided
     * as template arguments.
     *
     * @tparam FN The type of the function object for which to find the root.
     * @tparam DFN The type of the derivative function object.
     * @requires FN and DFN should be callable with double and return a floating point type.
     */
    template<typename FN, typename DFN>
        requires(IsFloatInvocable< FN > && IsFloatInvocable< DFN >) ||
                (IsComplex< std::invoke_result_t< FN, double > > && IsComplex< std::invoke_result_t< DFN, double > >)
    class Newton final : public impl::PolishingBase< Newton< FN, DFN > >
    {
        /*
         * Private alias declarations.
         */
        using BASE = impl::PolishingBase< Newton< FN, DFN > >;

    public:
        /*
         * Public alias declarations
         */
        using function_type = FN;
        using deriv_type    = DFN;

        /**
         * @brief Constructor, taking function objects for the function and its derivative.
         * @param objective The function object for the function.
         * @param deriv The function object for the derivative.
         */
        explicit Newton(FN objective, DFN deriv)
            : BASE(objective, deriv)
        {}

        /**
         * @brief Perform one iteration. This is the main algorithm of the Newton method.
         */
        void iterate() { BASE::m_guess = BASE::m_guess - BASE::evaluate(BASE::m_guess) / BASE::derivative(BASE::m_guess); }
    };

    template< typename FN, typename DERIV >
    requires std::invocable< FN, double > && std::invocable< DERIV, double >
    Newton(FN func, DERIV deriv) -> Newton< decltype(func), decltype(deriv) >;

    /**
     * @brief The fdfsolve function is a convenience function for running a polishing solver (i.e. with derivative),
     * without dealing with low level details. If fine grained control is needed, such as advanced search stopping
     * criteria or running each iteration manually, please see the documentation for the solver classes.
     * @tparam SOLVER The type of the solver. This could be the Newton or DNewton solvers, but any solver with the
     * correct interface can be used.
     * @param solver The actual solver object.
     * @param guess The initial guess of the root. The guess must be reasonably close to the actual root.
     * @param eps The max. allowed error.
     * @param maxiter The max. number of allowed iterations.
     * @return The root estimate.
     */
    template< typename SOLVER >
    requires requires(SOLVER solver, typename SOLVER::FUNCTION_RETURN_T guess) {
                 // clang-format off
                 { solver.evaluate(0.0) } -> std::same_as< typename SOLVER::FUNCTION_RETURN_T >;
                 { solver.init(guess) };
                 { solver.iterate() };
                 // clang-format on
    }
    auto fdfsolve(SOLVER                             solver,
                  typename SOLVER::FUNCTION_RETURN_T guess,
                  std::floating_point auto           eps,
                  //     = nxx::EPS,
                  int maxiter) // = nxx::MAXITER)
    {
        using EXPECTED_T = impl::RootErrorImpl< typename SOLVER::FUNCTION_RETURN_T >;
        using RETURN_T   = tl::expected< typename SOLVER::FUNCTION_RETURN_T, EXPECTED_T >;

        solver.init(guess);
        RETURN_T result = solver.result();

        // Check for NaN or Inf.
        if (!std::isfinite(abs(solver.evaluate(result.value())))) {
            result = tl::make_unexpected(EXPECTED_T("Invalid initial guess!", RootErrorType::NumericalError, result.value()));
            return result;
        }

        // Begin iteration loop.
        int iter = 1;
        while (true) {
            result = solver.result();

            // Check for NaN or Inf
            if (!std::isfinite(abs(result.value()))) {
                result = tl::make_unexpected(EXPECTED_T("Non-finite result!", RootErrorType::NumericalError, result.value(), iter));
                break;
            }

            // Check for convergence
            if (abs(solver.evaluate(solver.result())) < eps) break;

            // Check for max. iterations
            if (iter >= maxiter) {
                result = tl::make_unexpected(
                    EXPECTED_T("Maximum number of iterations exceeded!", RootErrorType::MaxIterationsExceeded, result.value(), iter));
                break;
            }

            // Perform one iteration
            ++iter;
            solver.iterate();
        }

        return result;
    }
}    // namespace nxx::roots

#endif    // NUMERIXX_ROOTPOLISHING_HPP

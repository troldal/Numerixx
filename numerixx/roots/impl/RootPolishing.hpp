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
#include <Poly.hpp>

// ===== Standard Library Includes

namespace nxx::roots
{
    // =================================================================================================================
    //
    //  88888888ba               88  88             88           88
    //  88      "8b              88  ""             88           ""
    //  88      ,8P              88                 88
    //  88aaaaaa8P'  ,adPPYba,   88  88  ,adPPYba,  88,dPPYba,   88  8b,dPPYba,    ,adPPYb,d8
    //  88""""""'   a8"     "8a  88  88  I8[    ""  88P'    "8a  88  88P'   `"8a  a8"    `Y88
    //  88          8b       d8  88  88   `"Y8ba,   88       88  88  88       88  8b       88
    //  88          "8a,   ,a8"  88  88  aa    ]8I  88       88  88  88       88  "8a,   ,d88
    //  88           `"YbbdP"'   88  88  `"YbbdP"'  88       88  88  88       88   `"YbbdP"Y8
    //                                                                             aa,    ,88
    //                                                                              "Y8bbdP"
    // =================================================================================================================

    /*
     * Private implementation details.
     */
    namespace detail
    {
        /**
         * @brief Provides a base class template for root polishing algorithms.
         *
         * The PolishingBase class template serves as a foundational component for
         * algorithms that refine or 'polish' roots of a given function. It encapsulates
         * common functionalities such as storing the objective function, its derivative,
         * and the current guess of the root. This class enforces certain type constraints
         * on the template parameters to ensure compatibility with root polishing algorithms.
         *
         * @tparam SUBCLASS The subclass inheriting from PolishingBase.
         * @tparam FUNCTION_T The type of the function for which the root is being polished.
         * @tparam DERIV_T The type of the derivative function of FUNCTION_T.
         * @tparam ARG_T The type of the argument to the function and its derivative.
         */
        template<typename SUBCLASS, typename FUNCTION_T, typename DERIV_T, typename ARG_T>
            requires std::same_as< typename PolishingTraits< SUBCLASS >::FUNCTION_T, FUNCTION_T > &&
                     std::same_as< typename PolishingTraits< SUBCLASS >::DERIV_T, DERIV_T > &&
                     IsFloatOrComplexInvocable< FUNCTION_T > &&
                     IsFloatOrComplexInvocable< DERIV_T > &&
                     IsFloatOrComplex< ARG_T >
        class PolishingBase
        {
            friend SUBCLASS;

        public:
            static constexpr bool IsPolishingSolver = true; /**< Flag indicating the class is a polishing solver. */

            using FUNCT_RES_T = std::invoke_result_t< FUNCTION_T, ARG_T >;   /**< Result type of the function. */
            using DERIV_RES_T = std::invoke_result_t< DERIV_T, ARG_T >;      /**< Result type of the derivative function. */
            using RESULT_T = std::common_type_t< FUNCT_RES_T, DERIV_RES_T >; /**< Common type for results of function and derivative. */

        protected:
            ~PolishingBase() = default; /**< Protected destructor to prevent direct instantiation. */

        private:
            FUNCTION_T m_func{};               /**< The function object to find the root for. */
            DERIV_T    m_deriv{};              /**< The function object for the derivative. */
            RESULT_T   m_guess;                /**< The current root estimate. */

        public:

            /**
             * @brief Constructs the PolishingBase with a function, its derivative, and an initial guess.
             * @param objective The function for which the root is being refined.
             * @param derivative The derivative of the objective function.
             * @param guess Initial guess for the root.
             */
            PolishingBase(FUNCTION_T objective, DERIV_T derivative, ARG_T guess)
                : m_func{ objective },
                  m_deriv{ derivative },
                  m_guess { guess }
            {}

            PolishingBase(const PolishingBase& other)                = default; /**< Default copy constructor. */
            PolishingBase(PolishingBase&& other) noexcept            = default; /**< Default move constructor. */
            PolishingBase& operator=(const PolishingBase& other)     = default; /**< Default copy assignment operator. */
            PolishingBase& operator=(PolishingBase&& other) noexcept = default; /**< Default move assignment operator. */

            /**
             * @brief Evaluates the function with the given value.
             * @tparam T Type of the value, must be a float or complex type.
             * @param value The value to evaluate the function at.
             * @return The result of the function evaluation.
             */
            template<typename T>
                requires nxx::IsFloatOrComplex< T >
            auto evaluate(T value) { return m_func(value); }

            /**
             * @brief Evaluates the derivative function with the given value.
             * @tparam T Type of the value, must be a float or complex type.
             * @param value The value to evaluate the derivative at.
             * @return The result of the derivative function evaluation.
             */
            template<typename T>
                requires nxx::IsFloatOrComplex< T >
            auto derivative(T value) { return std::invoke(m_deriv, value); }

            /**
             * @brief Returns the current result of the solver.
             * @details This function returns the current root estimate.
             *          It throws an exception if the solver has not been initialized.
             * @throws NumerixxError If the solver has not been initialized.
             * @return The current root estimate.
             */
            RESULT_T current() const
            {
                return m_guess;
            }

            void iterate() { std::invoke(static_cast< SUBCLASS& >(*this)); }
        };
    } // namespace impl

    // =================================================================================================================
    //
    //  888b      88
    //  8888b     88                                    ,d
    //  88 `8b    88                                    88
    //  88  `8b   88   ,adPPYba,  8b      db      d8  MM88MMM  ,adPPYba,   8b,dPPYba,
    //  88   `8b  88  a8P_____88  `8b    d88b    d8'    88    a8"     "8a  88P'   `"8a
    //  88    `8b 88  8PP"""""""   `8b  d8'`8b  d8'     88    8b       d8  88       88
    //  88     `8888  "8b,   ,aa    `8bd8'  `8bd8'      88,   "8a,   ,a8"  88       88
    //  88      `888   `"Ybbd8"'      YP      YP        "Y888  `"YbbdP"'   88       88
    //
    // =================================================================================================================


    /**
     * @brief Defines the Newton class for performing Newton's method root polishing.
     *
     * The Newton class template is a specialized implementation of the root polishing
     * algorithm known as Newton's method (or the Newton-Raphson method). It inherits
     * from a base class that provides common functionalities for root polishing algorithms,
     * and adds the specific iteration logic for Newton's method. This class is templated
     * to accept a function, its derivative, and an optional argument type.
     *
     * @tparam FN Type of the function for which the root is being refined.
     * @tparam DFN Type of the derivative function of FN.
     * @tparam ARG_T Type of the argument to the function and its derivative, defaults to double.
     */
    template<IsFloatOrComplexInvocable FN, IsFloatOrComplexInvocable DFN, IsFloatOrComplex ARG_T = double>
    class Newton final : public detail::PolishingBase< Newton< FN, DFN, ARG_T >, FN, DFN, ARG_T >
    {
        using BASE = detail::PolishingBase< Newton< FN, DFN, ARG_T >, FN, DFN, ARG_T >; /**< Base class alias for readability. */

    public:
        using BASE::BASE; /**< Inherits constructors from PolishingBase. */

        /**
         * @brief Performs a single iteration of Newton's method.
         * @details This method updates the root estimate using the Newton-Raphson formula.
         *          It assumes the class has been properly initialized.
         */
        void operator()() { BASE::m_guess = BASE::m_guess - BASE::evaluate(BASE::m_guess) / BASE::derivative(BASE::m_guess); }
    };

    /**
     * @brief Deduction guides for Newton class.
     * Allows the type of Newton class to be deduced from the constructor parameters.
     */
    template<typename FN, typename DERIV, typename ARG_T>
        requires IsFloatOrComplexInvocable< FN > && IsFloatOrComplexInvocable< DERIV > && IsFloatOrComplex< ARG_T >
    Newton(FN, DERIV, ARG_T) -> Newton< FN, DERIV, ARG_T >;


    // =================================================================================================================
    //
    //  ad88888ba
    // d8"     "8b                                                    ,d
    // Y8,                                                            88
    // `Y8aaaaa,     ,adPPYba,   ,adPPYba,  ,adPPYYba,  8b,dPPYba,  MM88MMM
    //   `"""""8b,  a8P_____88  a8"     ""  ""     `Y8  88P'   `"8a   88
    //         `8b  8PP"""""""  8b          ,adPPPPP88  88       88   88
    // Y8a     a8P  "8b,   ,aa  "8a,   ,aa  88,    ,88  88       88   88,
    //  "Y88888P"    `"Ybbd8"'   `"Ybbd8"'  `"8bbdP"Y8  88       88   "Y888
    //
    // =================================================================================================================

    /**
     * @brief Defines the Secant class for performing the secant method of root polishing.
     *
     * The Secant class template implements the secant method, an iterative root finding algorithm
     * that uses a succession of roots of secant lines to approximate a root of a function. This class
     * template inherits from a base class that provides common functionalities for root polishing algorithms,
     * and adds specific logic for the secant method iterations. It is designed to handle an initial guess
     * and a subsequent series of approximations to converge on the root. The secant method is a derivative-free
     * alternative to Newton's method and is particularly useful when the derivative of the function is not
     * readily available or is expensive to compute.
     *
     * @tparam FN The type of the function for which the root is being polished.
     * @tparam DFN The type of the derivative function of FN, used in the initial step.
     * @tparam ARG_T The type of the argument to the function, defaults to double.
     */
    template<IsFloatOrComplexInvocable FN, IsFloatOrComplexInvocable DFN, IsFloatOrComplex ARG_T = double>
    class Secant final : public detail::PolishingBase< Secant< FN, DFN, ARG_T >, FN, DFN, ARG_T >
    {
        using BASE = detail::PolishingBase< Secant< FN, DFN, ARG_T >, FN, DFN, ARG_T >; /**< Base class alias for readability. */

        ARG_T m_prevGuess;             /**< Stores the previous guess for the root. */
        bool  m_hasPrevGuess{ false }; /**< Flag to indicate whether a previous guess is available. */
        bool  m_firstStep{ true };     /**< Flag to indicate whether the first step is to be taken. */

    public:
        using BASE::BASE; /**< Inherits constructors from PolishingBase. */

        /**
         * @brief Performs a single iteration of the secant method.
         * @details This method switches between a Newton-Raphson step for the first iteration
         *          and a Secant step for subsequent iterations.
         * @throws std::runtime_error If a division by near-zero occurs.
         */
        void operator()()
        {
            if (m_firstStep) {
                ARG_T f_x       = BASE::evaluate(BASE::m_guess);
                ARG_T f_prime_x = BASE::derivative(BASE::m_guess);

                if (abs(f_prime_x) < std::numeric_limits< ARG_T >::epsilon()) {
                    throw std::runtime_error("Division by near-zero in Newton-Raphson step.");
                    // TODO: Return a tl::expected instead of throwing an exception.
                }

                m_prevGuess = BASE::m_guess;
                BASE::m_guess -= f_x / f_prime_x;
                m_firstStep    = false;
                m_hasPrevGuess = true;
            }
            else {
                ARG_T f_x      = BASE::evaluate(BASE::m_guess);
                ARG_T f_x_prev = BASE::evaluate(m_prevGuess);

                if (abs(f_x - f_x_prev) < std::numeric_limits< ARG_T >::epsilon()) {
                    throw std::runtime_error("Division by near-zero in Secant method.");
                    // TODO: Return a tl::expected instead of throwing an exception.
                }

                ARG_T newGuess = BASE::m_guess - f_x * (BASE::m_guess - m_prevGuess) / (f_x - f_x_prev);
                m_prevGuess    = BASE::m_guess;
                BASE::m_guess  = newGuess;
            }
        }
    };

    /**
     * @brief Deduction guides for Secant class.
     * Allows the type of Secant class to be deduced from the constructor parameters.
     */
    template<typename FN, typename DFN, typename ARG_T>
        requires IsFloatOrComplexInvocable< FN > && IsFloatOrComplexInvocable< DFN > && IsFloatOrComplex< ARG_T >
    Secant(FN, DFN, ARG_T) -> Secant< FN, DFN, ARG_T >;


    // =================================================================================================================
    //
    //  ad88888ba                        ad88     ad88
    // d8"     "8b  ,d                  d8"      d8"
    // Y8,          88                  88       88
    // `Y8aaaaa,  MM88MMM  ,adPPYba,  MM88MMM  MM88MMM  ,adPPYba,  8b,dPPYba,   ,adPPYba,   ,adPPYba,  8b,dPPYba,
    //   `"""""8b,  88    a8P_____88    88       88    a8P_____88  88P'   `"8a  I8[    ""  a8P_____88  88P'   `"8a
    //         `8b  88    8PP"""""""    88       88    8PP"""""""  88       88   `"Y8ba,   8PP"""""""  88       88
    // Y8a     a8P  88,   "8b,   ,aa    88       88    "8b,   ,aa  88       88  aa    ]8I  "8b,   ,aa  88       88
    //  "Y88888P"   "Y888  `"Ybbd8"'    88       88     `"Ybbd8"'  88       88  `"YbbdP"'   `"Ybbd8"'  88       88
    //
    // =================================================================================================================

    /**
     * @brief Defines the Steffensen class for performing Steffensen's method of root polishing.
     *
     * Steffensen's method is an iterative root finding algorithm that improves upon the simple fixed-point
     * iteration by incorporating a form of Aitken's Δ² process. This class template, `Steffensen`, inherits
     * from a base class that provides common functionalities for root polishing algorithms, and adds the
     * specific logic for Steffensen's method iterations. It starts with a Newton-Raphson step and then
     * switches to Steffensen's method for subsequent iterations. This method is particularly effective
     * for functions where the derivative is difficult to compute or is not readily available.
     *
     * @tparam FN The type of the function for which the root is being polished.
     * @tparam DFN The type of the derivative function of FN, used in the initial step.
     * @tparam ARG_T The type of the argument to the function, defaults to double.
     */
    template<IsFloatOrComplexInvocable FN, IsFloatOrComplexInvocable DFN, IsFloatOrComplex ARG_T = double>
    class Steffensen final : public detail::PolishingBase< Steffensen< FN, DFN, ARG_T >, FN, DFN, ARG_T >
    {
        using BASE = detail::PolishingBase< Steffensen< FN, DFN, ARG_T >, FN, DFN, ARG_T >; /**< Base class alias for readability. */

        bool m_firstStep{ true }; /**< Flag to indicate whether the first step is to be taken. */

    public:
        using BASE::BASE; /**< Inherits constructors from PolishingBase. */

        /**
         * @brief Performs a single iteration of the hybrid Steffensen method.
         * @details Uses Newton-Raphson for the first iteration and Steffensen's method subsequently.
         * @throws std::runtime_error If a division by near-zero occurs.
         */
        void operator()()
        {
            if (m_firstStep) {
                // Perform a Newton-Raphson step for the first iteration.
                ARG_T f_x       = BASE::evaluate(BASE::m_guess);
                ARG_T f_prime_x = BASE::derivative(BASE::m_guess);

                if (abs(f_prime_x) < std::numeric_limits< ARG_T >::epsilon()) {
                    throw std::runtime_error("Division by near-zero in Newton-Raphson step.");
                    // TODO: Return a tl::expected instead of throwing an exception.
                }

                BASE::m_guess -= f_x / f_prime_x;
                m_firstStep = false;
            }
            else {
                // Perform a Steffensen's method step for subsequent iterations.
                ARG_T x  = BASE::m_guess;
                ARG_T fx = BASE::evaluate(x);

                ARG_T x1  = x + fx;
                ARG_T fx1 = BASE::evaluate(x1);
                ARG_T x2  = x1 + fx1;

                ARG_T denominator = fx1 - fx;
                if (abs(denominator) < std::numeric_limits< ARG_T >::epsilon()) {
                    throw std::runtime_error("Division by near-zero in Steffensen's method.");
                    // TODO: Return a tl::expected instead of throwing an exception.
                }

                BASE::m_guess = x - (fx * fx) / denominator;
            }
        }
    };

    /**
     * @brief Deduction guides for Steffensen class.
     * Allows the type of Steffensen class to be deduced from the constructor parameters.
     */
    template<typename FN, typename DFN, typename ARG_T>
        requires IsFloatOrComplexInvocable< FN > && IsFloatOrComplexInvocable< DFN > && IsFloatOrComplex< ARG_T >
    Steffensen(FN, DFN, ARG_T) -> Steffensen< FN, DFN, ARG_T >;


    // =================================================================================================================
    //
    //    ad88          88     ad88                          88
    //   d8"            88    d8"                            88
    //   88             88    88                             88
    // MM88MMM  ,adPPYb,88  MM88MMM  ,adPPYba,   ,adPPYba,   88  8b       d8   ,adPPYba,
    //   88    a8"    `Y88    88     I8[    ""  a8"     "8a  88  `8b     d8'  a8P_____88
    //   88    8b       88    88      `"Y8ba,   8b       d8  88   `8b   d8'   8PP"""""""
    //   88    "8a,   ,d88    88     aa    ]8I  "8a,   ,a8"  88    `8b,d8'    "8b,   ,aa
    //   88     `"8bbdP"Y8    88     `"YbbdP"'   `"YbbdP"'   88      "8"       `"Ybbd8"'
    //
    // =================================================================================================================


    namespace detail
    {
        /**
         * @brief Implements a generic root finding solver function template for polishing solvers.
         *
         * This function template, `fdfsolve_impl`, provides a generic implementation for root finding
         * algorithms that utilize polishing solvers. It is designed to work with solvers that conform
         * to the requirements of polishing solvers, such as having a defined `IsPolishingSolver` static member,
         * initialization, and iteration methods. The function handles initialization, iteration, and
         * convergence checking, returning the result along with any potential errors encountered during
         * the solving process.
         *
         * @tparam SOLVER The type of the solver to be used in root finding. Must conform to the polishing solver concept.
         */
        template<typename SOLVER>
            requires SOLVER::IsPolishingSolver
        auto fdfsolve_impl(SOLVER                solver,
                           IsFloat auto          eps,
                           std::integral auto    maxiter)
        {
            using ERROR_T = detail::RootErrorImpl< typename SOLVER::RESULT_T >;  /**< Type for error handling. */
            using RETURN_T = tl::expected< typename SOLVER::RESULT_T, ERROR_T >; /**< Type for the function return value. */

            RETURN_T result = solver.current();

            // Check for NaN or Inf in the initial guess evaluation.
            using std::isfinite;
            if (!isfinite(abs(solver.evaluate(result.value())))) {
                result = tl::make_unexpected(ERROR_T("Invalid initial guess!", RootErrorType::NumericalError, result.value()));
                return result;
            }

            // Begin the iteration loop.
            int iter = 1;
            while (true) {
                result = solver.current();

                // Check for non-finite results.
                if (!isfinite(abs(result.value()))) {
                    result = tl::make_unexpected(ERROR_T("Non-finite result!", RootErrorType::NumericalError, result.value(), iter));
                    break;
                }

                // Check for convergence.
                if (abs(solver.evaluate(solver.current())) < eps) break;

                // Check for exceeding the maximum number of iterations.
                if (iter >= maxiter) {
                    result = tl::make_unexpected(
                        ERROR_T("Maximum number of iterations exceeded!", RootErrorType::MaxIterationsExceeded, result.value(), iter));
                    break;
                }

                // Perform one iteration.
                ++iter;
                solver.iterate();
            }

            return result;
        }
    }

    /**
     * @brief Defines a high-level root finding function template `fdfsolve` with default parameters.
     *
     * The `fdfsolve` function template is an extension of the root finding interface that provides
     * default values for the tolerance (epsilon) and the maximum number of iterations. This overloaded
     * version allows users to specify only the function, its derivative, and an initial guess, while
     * using default values for epsilon and maximum iterations based on the type of the guess. It supports
     * different polishing solvers, making it versatile for various root finding needs.
     *
     * @tparam SOLVER_T The template class of the solver to be used. Must be a valid polishing solver type.
     * @tparam FN_T The type of the function for which the root is being refined.
     * @tparam DERIV_T The type of the derivative function of FN_T.
     * @tparam GUESS_T The type of the initial guess for the root.
     * @tparam EPS_T The type of the epsilon value for convergence check, defaulted based on GUESS_T.
     * @tparam ITER_T The type of the maximum iterations count, defaulted to int.
     */
    template<template< typename, typename, typename > class SOLVER_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloatOrComplexInvocable DERIV_T,
        IsFloatOrComplex GUESS_T,
        IsFloat EPS_T,
        std::integral ITER_T>
    auto fdfsolve(FN_T    function,
                  DERIV_T derivative,
                  GUESS_T guess,
                  EPS_T   eps, ITER_T maxiter)
    {
        // Instantiates the solver with the given function, its derivative, and types.
        auto solver = SOLVER_T< FN_T, DERIV_T, GUESS_T >(function, derivative, guess);

        // Delegates the solving process to fdfsolve_impl, passing in the solver and other parameters.
        return detail::fdfsolve_impl(solver, eps, maxiter);
    }
} // namespace nxx::roots

#endif    // NUMERIXX_ROOTPOLISHING_HPP

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

#include <algorithm>
#include <array>
#include <cassert>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <tuple>
#include <utility>

#include "calculus/Derivatives.hpp"
#include "poly/Polynomial.hpp"

namespace nxx::roots
{
    using nxx::poly::derivativeOf;
    using nxx::deriv::derivativeOf;


    // ========================================================================
    // ROOT-FINDING WITH DERIVATIVES
    // ========================================================================

    /*
     * Forward declaration of the DNewton class.
     */
    template< typename FN, typename DFN >
        requires std::invocable< FN, double > && std::invocable< DFN, double >
    class DNewton;

    /*
     * Forward declaration of the Newton class.
     */
    template< typename FN, typename DFN >
        requires std::invocable< FN, double > && std::invocable< DFN, double >
    class Newton;

    /*
     * Private implementation details.
     */
    namespace impl
    {
        /*
         * Forward declaration of the PolishingTraits class.
         */
        template< typename... SOLVER >
        struct PolishingTraits;

        /*
         * Specialization of the PolishingTraits class for Newton<FN, DFN>
         */
        template< typename FN, typename DFN >
        struct PolishingTraits< Newton< FN, DFN > >
        {
            using function_type = FN;
            using deriv_type    = DFN;
        };

        /*
         * Specialization of the PolishingTraits class for DNewton<FN, DFN>
         */
        template< typename FN, typename DFN >
        struct PolishingTraits< DNewton< FN, DFN > >
        {
            using function_type = FN;
            using deriv_type    = DFN;
        };

        /**
         * @brief A CRTP base class for root-finding algorithms using derivatives.
         *
         * This class serves as a CRTP base class for root-finding algorithms that require the use of derivatives.
         * It requires that the derived class provides the necessary function and derivative callable objects,
         * which should be invocable with a double argument.
         *
         * @tparam POLICY The derived class, which should inherit from PolishingBase.
         *
         * @note This class should not be used directly by clients; it is intended for use as a base class.
         * @note The derived class should satisfy the requirements for the PolishingTraits type trait.
         */
        template< typename POLICY >
            requires std::invocable< typename impl::PolishingTraits< POLICY >::function_type, double > &&
                     std::invocable< typename impl::PolishingTraits< POLICY >::deriv_type, double >
        class PolishingBase
        {
            /*
             * Friend declarations.
             */
            friend POLICY;

        private:
            using function_type = typename impl::PolishingTraits< POLICY >::function_type; /**< The type of the function object. */
            function_type m_func {};                                                       /**< The function object to find the root for. */

            using deriv_type = typename impl::PolishingTraits< POLICY >::deriv_type; /**< The type of the derivative function object. */
            deriv_type m_deriv {};                                                   /**< The function object for the derivative. */

            using RT = std::invoke_result_t< function_type, double >;                /**< The return type of the function object. */
            RT m_guess;                                                              /**< The current root estimate. */

        public:
            /**
             * @brief Constructor, taking function objects for the function and its derivative as arguments.
             * @param objective The function object to find the root for.
             * @param derivative The function object for the derivative
             * @note Constructor is private to avoid direct usage by clients.
             */
            explicit PolishingBase(function_type objective, deriv_type derivative)
                : m_func { objective },
                  m_deriv { derivative }
            {}

        public:
            /**
             * @brief Initialize the solver by setting the initial bounds around the root.
             * @tparam T The type of the root estimate.
             * @param guess The root estimate.
             */
            template< typename T >
                requires std::is_floating_point_v< T >
            void init(T guess)
            {
                m_guess = guess;
            }

            /**
             * @brief Evaluate the function object at a given point.
             * @tparam T The type of the argument (must be floating point type)
             * @param value The value at which to evaluate the function.
             * @return The result of the evaluation.
             */
            template< typename T >
                requires std::is_floating_point_v< T >
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
                requires std::is_floating_point_v< T >
            auto derivative(T value)
            {
                return m_deriv(value);
            }

            /**
             * @brief Get the current estimate of the root.
             * @return A const reference to the root.
             */
            const auto& result() const { return m_guess; }
        };
    }    // namespace impl

    /**
     * @brief The DNewton class is a derived class of the PolishingBase CRTP base class.
     * It implements Newton's method, but uses numerical differentation to evaluate the derivatives,
     * i.e. discrete Newton's method.
     * @tparam FN The type of the function object for the function for which to find the root.
     * @tparam DFN The type of the function object for the derivative (deducet automatically).
     */
    template< typename FN, typename DFN >
        requires std::invocable< FN, double > && std::invocable< DFN, double >
    class DNewton final : public impl::PolishingBase< DNewton< FN, DFN > >
    {
        /*
         * Private alias declarations.
         */
        using Base = impl::PolishingBase< DNewton< FN, DFN > >;

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
        explicit DNewton(FN objective) : Base(objective, derivativeOf(objective)) {}


        /**
         * @brief Perform one iteration. This is the main algorithm of the DNewton method.
         */
        void iterate()
        {
            Base::m_guess = Base::m_guess - Base::evaluate(Base::m_guess) / Base::derivative(Base::m_guess);
        }
    };

    /**
     * @brief Deduction guide for the DNewton algorithm
     */
    template< typename FN >
    DNewton(FN objective) -> DNewton< FN, std::function< decltype(objective(0.0))(decltype(objective(0.0))) > >;

    /**
     * @brief The Newton class is a derived class of the PolishingBase CRTP base class.
     * It implements Newton's method for root finding with derivatives.
     * @tparam FN The type of the function object for the function for which to find the root.
     * @tparam DFN The type of the function object for the derivative.
     */
    template< typename FN, typename DFN >
        requires std::invocable< FN, double > && std::invocable< DFN, double >
    class Newton final : public impl::PolishingBase< Newton< FN, DFN > >
    {
        /*
         * Private alias declarations.
         */
        using Base = impl::PolishingBase< Newton< FN, DFN > >;

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
        explicit Newton(FN objective, DFN deriv) : Base(objective, deriv) {}

        /**
         * @brief Perform one iteration. This is the main algorithm of the Newton method.
         */
        void iterate()
        {
            Base::m_guess = Base::m_guess - Base::evaluate(Base::m_guess) / Base::derivative(Base::m_guess);
        }
    };

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
    inline auto fdfsolve(SOLVER solver, double guess, double eps = 1.0E-6, int maxiter = 100)
        -> tl::expected<double, error::RootError>
    {
        // using RT = decltype(solver.evaluate(0.0));

        solver.init(guess);

        int iter = 1;
        while (true) {
            if (abs(solver.evaluate(solver.result())) < eps) break;
            solver.iterate();
            if (!std::isfinite(solver.result())) return tl::make_unexpected(error::RootError("Root Error!"));;

            ++iter;
            if (iter > maxiter) break;
        }

        return solver.result();
    }

}    // namespace nxx::roots

#endif    // NUMERIXX_ROOTPOLISHING_HPP

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

#ifndef NUMERIXX_MULTIROOTS_HPP
#define NUMERIXX_MULTIROOTS_HPP

#include "../calculus/Jacobian.hpp"
#include "../linalg/FactorizeGJ.hpp"
#include "../linalg/Matrix.hpp"
#include "MultiFunction.hpp"
#include <functional>

namespace nxx::multiroots
{
    template< typename TMultiFunc >
    // requires std::same_as<typename FUNCARR::value_type, std::function<typename FUNCARR::value_type::result_type(std::vector<typename
    // FUNCARR::value_type::result_type>)> >
    class DMultiNewton;

    /*
     * Forward declaration of the PolishingTraits class.
     */
    template< typename SOLVER >
    struct MultirootsSolverTraits;

    /*
     * Specialization of the PolishingTraits class for Newton<FN, DFN>
     */
    template< typename TMultiFunc >
    struct MultirootsSolverTraits< DMultiNewton< TMultiFunc > >
    {
        using function_array = typename TMultiFunc::function_array;
        using function_type = typename TMultiFunc::function_type;
        using return_type   = typename TMultiFunc::return_type;
        using argument_type = typename TMultiFunc::argument_type;
    };

    // ========================================================================
    // MULTIROOT-FINDING WITH DERIVATIVES
    // ========================================================================

    template< typename SOLVER >
    // requires std::invocable<typename impl::PolishingTraits<SOLVER>::function_type, double> &&
    //          std::invocable<typename impl::PolishingTraits<SOLVER>::deriv_type, double>
    class MultirootBase
    {
        /*
         * Friend declarations.
         */
        friend SOLVER;

    public:
        using function_array = typename MultirootsSolverTraits< SOLVER >::function_array;
        using function_type = typename MultirootsSolverTraits< SOLVER >::function_type;
        using return_type   = typename MultirootsSolverTraits< SOLVER >::return_type;
        using argument_type = typename MultirootsSolverTraits< SOLVER >::argument_type;

        MultiFunction< function_array > m_functions {}; /**< The function object to find the root for. */

        using RT = nxx::linalg::Matrix< return_type >;
        RT m_guess; /**< The current root estimate. */

    public:
        /**
         * @brief Constructor, taking function objects for the function and its derivative as arguments.
         * @param objective The function object to find the root for.
         * @param derivative The function object for the derivative
         * @note Constructor is private to avoid direct usage by clients.
         */
        explicit MultirootBase(const MultiFunction< function_array >& functions) : m_functions { functions }, m_guess(functions.size(), 1) {}

    public:
        /**
         * @brief
         * @tparam ARR
         * @param guess
         */
        template< typename ARR >
            requires std::convertible_to< typename ARR::value_type, return_type >
        void init(const ARR& guess)
        {
            m_guess = guess;
        }

        /**
         * @brief
         * @tparam T
         * @param guess
         */
        template< typename T >
            requires std::convertible_to< T, return_type >
        void init(std::initializer_list< T > guess)
        {
            m_guess = guess;
        }

        /**
         * @brief
         * @tparam ARR
         * @param values
         * @return
         */
        template< typename ARR >
            requires std::convertible_to< typename ARR::value_type, return_type >
        auto evaluate(const ARR& values)
        {

            return m_functions(values);
//            using numerix::linalg::Matrix;
//
//            std::vector< return_type > vals = values;
//
//            std::vector< return_type > result;
//            result.reserve(m_functions.size());
//            for (auto& f : m_functions) result.emplace_back(f(vals));
//
//            Matrix< return_type > res(result.size(), 1);
//            res = result;
//
//            return res;
        }

        /**
         * @brief Evaluate the derivative at a given point.
         * @tparam T The type of the argument (must be floating point type)
         * @param value The value at which to evaluate the derivative.
         * @return The result of the evaluation.
         */
        // template<typename T>
        //     requires std::is_floating_point_v<T>
        // auto derivative(T values)
        //{
        //     return m_deriv(value);
        // }

        /**
         * @brief Get the current estimate of the root.
         * @return A const reference to the root.
         */
        const auto& result() const { return m_guess; }
    };

    /**
     * @brief
     * @tparam FUNCARR
     */
    template< typename TMultiFunc >
    // requires std::same_as<typename FUNCARR::value_type, std::function<typename FUNCARR::value_type::result_type(std::vector<typename
    // FUNCARR::value_type::result_type>)> >
    class DMultiNewton final : public MultirootBase< DMultiNewton< TMultiFunc > >
    {
        /*
         * Private alias declarations.
         */
        using Base = MultirootBase< DMultiNewton< TMultiFunc > >;

    public:
        /*
         * Public alias declarations.
         */
        using function_type = typename TMultiFunc::function_type;
        using return_type   = typename TMultiFunc::return_type;
        using argument_type = typename TMultiFunc::argument_type;

        /**
         * @brief Constructor, taking the function object as an argument.
         * @param objective The function object for which to find the root.
         */
        explicit DMultiNewton(const TMultiFunc& funcArray) : Base(funcArray) {}

        /**
         * @brief Perform one iteration. This is the main algorithm of the DMultiNewton method.
         */
        void iterate()
        {
            using nxx::deriv::computeJacobian;
            using nxx::linalg::FactorizeGJ;
            using nxx::linalg::Matrix;

            // ===== Create vector of function evaluations, using vector of guesses
            Matrix< return_type > functionEvals(Base::m_functions.size(), 1);
            functionEvals = Base::evaluate(Base::m_guess) * (-1);

            // ===== Solve the linear system using the jacobian matrix, and estimate new guesses
            auto jacobian = computeJacobian(Base::m_functions.functionArray(), Base::m_guess);
            auto [inv, x] = FactorizeGJ(jacobian, functionEvals);
            Base::m_guess += x;
        }
    };
}    // namespace numerix::multiroots

#endif    // NUMERIXX_MULTIROOTS_HPP

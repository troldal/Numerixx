/*
    888b      88  88        88  88b           d88  88888888888  88888888ba   88  8b        d8
    8888b     88  88        88  888b         d888  88           88      "8b  88   Y8,    ,8P
    88 `8b    88  88        88  88`8b       d8'88  88           88      ,8P  88    `8b  d8'
    88  `8b   88  88        88  88 `8b     d8' 88  88aaaaa      88aaaaaa8P'  88      Y88P
    88   `8b  88  88        88  88  `8b   d8'  88  88"""""      88""""88'    88      d88b
    88    `8b 88  88        88  88   `8b d8'   88  88           88    `8b    88    ,8P  Y8,
    88     `8888  Y8a.    .a8P  88    `888'    88  88           88     `8b   88   d8'    `8b
    88      `888   `"Y8888Y"'   88     `8'     88  88888888888  88      `8b  88  8P        Y8

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

#ifndef NUMERIX_MULTIROOTS_HPP
#define NUMERIX_MULTIROOTS_HPP

#include "../linalg/Matrix.hpp"
#include "../linalg/FactorizeGJ.hpp"
#include "../calculus/Jacobian.hpp"
#include <functional>

namespace numerix::multiroots
{
    template<typename FUNCARR>
        requires std::same_as<typename FUNCARR::value_type, std::function<typename FUNCARR::value_type::result_type(std::vector<typename FUNCARR::value_type::result_type>)> >
    class DMultiNewton;

    /*
         * Forward declaration of the PolishingTraits class.
     */
    template<typename FUNCARR>
    struct MultirootsTraits;

    /*
         * Specialization of the PolishingTraits class for Newton<FN, DFN>
     */
    template<typename FUNCARR>
    struct MultirootsTraits<DMultiNewton<FUNCARR>>
    {
        using function_type = typename FUNCARR::value_type;
        using result_type = typename function_type::result_type;
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

    private:
        using function_type = typename MultirootsTraits<SOLVER>::function_type;//std::function< double(std::vector< double >) >;
        using result_type = typename MultirootsTraits<SOLVER>::result_type;

        std::vector< function_type > m_functions {}; /**< The function object to find the root for. */

        // using deriv_type = typename impl::PolishingTraits<SOLVER>::deriv_type;
        // deriv_type m_deriv {}; /**< The function object for the derivative. */

        using RT = numerix::linalg::Matrix< result_type >;
        RT m_guess; /**< The current root estimate. */

    public:
        /**
         * @brief Constructor, taking function objects for the function and its derivative as arguments.
         * @param objective The function object to find the root for.
         * @param derivative The function object for the derivative
         * @note Constructor is private to avoid direct usage by clients.
         */
        explicit MultirootBase(const std::vector< function_type >& functions) : m_functions { functions },
                                                                                m_guess(functions.size(), 1) {}

    public:

        /**
         * @brief
         * @tparam ARR
         * @param guess
         */
        template<typename ARR>
            requires std::convertible_to<typename ARR::value_type, result_type>
        void init(const ARR& guess) { m_guess = guess; }

        /**
         * @brief
         * @tparam T
         * @param guess
         */
        template<typename T>
            requires std::convertible_to<T, result_type>
        void init(std::initializer_list<T> guess) { m_guess = guess; }

        /**
         * @brief
         * @tparam ARR
         * @param values
         * @return
         */
        template<typename ARR>
            requires std::convertible_to<typename ARR::value_type, result_type>
        auto evaluate(const ARR& values) {

            using numerix::linalg::Matrix;

            std::vector< result_type > vals = values;

            std::vector< result_type > result;
            result.reserve(m_functions.size());
            for (auto& f : m_functions) result.emplace_back(f(vals));

            Matrix< result_type > res(result.size(), 1);
            res = result;

            return res;

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
    template<typename FUNCARR>
        requires std::same_as<typename FUNCARR::value_type, std::function<typename FUNCARR::value_type::result_type(std::vector<typename FUNCARR::value_type::result_type>)> >
    class DMultiNewton final : public MultirootBase< DMultiNewton<FUNCARR> >
    {
        /*
         * Private alias declarations.
         */
        using Base = MultirootBase< DMultiNewton<FUNCARR> >;

    public:
        /*
         * Public alias declarations.
         */
        using function_type = typename FUNCARR::value_type;
        using result_type = typename function_type::result_type;

        /**
         * @brief Constructor, taking the function object as an argument.
         * @param objective The function object for which to find the root.
         */
        explicit DMultiNewton(const FUNCARR& functions) : Base(std::vector(functions.begin(), functions.end())) {}

        /**
         * @brief Perform one iteration. This is the main algorithm of the DMultiNewton method.
         */
        void iterate()
        {
            using numerix::linalg::Matrix;
            using numerix::deriv::computeJacobian;
            using numerix::linalg::FactorizeGJ;

            // ===== Create vector of function evaluations, using vector of guesses
            Matrix< result_type >functionEvals(Base::m_functions.size(), 1);
            functionEvals = Base::evaluate(Base::m_guess) * (-1);

            // ===== Solve the linear system using the jacobian matrix, and estimate new guesses
            auto jacobian = computeJacobian(Base::m_functions, Base::m_guess);
            auto [inv, x] = FactorizeGJ(jacobian, functionEvals);
            Base::m_guess += x;
        }
    };
}    // namespace numerix::multiroots

#endif    // NUMERIX_MULTIROOTS_HPP

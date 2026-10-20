//
// Created by kenne on 13/11/2023.
//

#pragma once

#include "OptimCommon.hpp"
#include "_external.hpp"

#include <cmath>
#include <functional>

#include <Concepts.hpp>
#include <Deriv.hpp>
#include <Interp.hpp>
#include <Poly.hpp>

namespace nxx::optim
{
    namespace detail
    {

        // =================================================================================================================
        //
        //   ,ad8888ba,                         88888888ba               88  88             88
        //  d8"'    `"8b                 ,d     88      "8b              88  ""             88
        // d8'        `8b                88     88      ,8P              88                 88
        // 88          88  8b,dPPYba,  MM88MMM  88aaaaaa8P'  ,adPPYba,   88  88  ,adPPYba,  88,dPPYba,
        // 88          88  88P'    "8a   88     88""""""'   a8"     "8a  88  88  I8[    ""  88P'    "8a
        // Y8,        ,8P  88       d8   88     88          8b       d8  88  88   `"Y8ba,   88       88
        //  Y8a.    .a8P   88b,   ,a8"   88,    88          "8a,   ,a8"  88  88  aa    ]8I  88       88
        //   `"Y8888Y"'    88`YbbdP"'    "Y888  88           `"YbbdP"'   88  88  `"YbbdP"'  88       88
        //                 88
        //                 88
        // =================================================================================================================

        template< typename DERIVED, typename FUNCTION_T, typename DERIV_T, typename ARG_T, typename MODE_T >
        class OptimPolishBase
        {
            friend DERIVED;

        public:
            static constexpr bool IsDerivativeOptimizer = true;

            using FUNCT_RES_T = std::invoke_result_t< FUNCTION_T, ARG_T >;
            using DERIV_RES_T = std::invoke_result_t< DERIV_T, ARG_T >;
            using RESULT_T    = std::common_type_t< FUNCT_RES_T, DERIV_RES_T >;

        protected:
            ~OptimPolishBase() = default;

        private:
            FUNCTION_T m_func {};
            DERIV_T    m_deriv {};
            RESULT_T   m_guess;
            MODE_T     m_mode;

        public:
            OptimPolishBase(FUNCTION_T objective, DERIV_T derivative, ARG_T guess, MODE_T mode = {})
                : m_func { objective },
                  m_deriv { derivative },
                  m_guess { guess },
                  m_mode(mode)
            {}

            OptimPolishBase(const OptimPolishBase& other)     = default; /**< Default copy constructor. */
            OptimPolishBase(OptimPolishBase&& other) noexcept = default; /**< Default move constructor. */

            OptimPolishBase& operator=(const OptimPolishBase& other)     = default; /**< Default copy assignment operator. */
            OptimPolishBase& operator=(OptimPolishBase&& other) noexcept = default; /**< Default move assignment operator. */

            template< typename T >
            requires nxx::IsFloatOrComplex< T >
            auto evaluate(T value)
            {
                return m_func(value);
            }

            template< typename T >
            requires nxx::IsFloatOrComplex< T >
            auto derivative(T value)
            {
                return std::invoke(m_deriv, value);
            }

            RESULT_T current() const { return m_guess; }

            void iterate() { std::invoke(static_cast< DERIVED& >(*this)); }
        };
    }    // namespace detail

    // =================================================================================================================
    //
    //   ,ad8888ba,                                    88  88
    //  d8"'    `"8b                                   88  ""                            ,d
    // d8'                                             88                                88
    // 88             8b,dPPYba,  ,adPPYYba,   ,adPPYb,88  88   ,adPPYba,  8b,dPPYba,  MM88MMM
    // 88      88888  88P'   "Y8  ""     `Y8  a8"    `Y88  88  a8P_____88  88P'   `"8a   88
    // Y8,        88  88          ,adPPPPP88  8b       88  88  8PP"""""""  88       88   88
    //  Y8a.    .a88  88          88,    ,88  "8a,   ,d88  88  "8b,   ,aa  88       88   88,
    //   `"Y88888P"   88          `"8bbdP"Y8   `"8bbdP"Y8  88   `"Ybbd8"'  88       88   "Y888
    //
    // =================================================================================================================

    template< IsFloatOrComplexInvocable FN, IsFloatOrComplexInvocable DFN, IsFloatOrComplex ARG_T, typename MODE_T = Minimize >
    class GradientDescent final : public detail::OptimPolishBase< GradientDescent< FN, DFN, ARG_T, MODE_T >, FN, DFN, ARG_T, MODE_T >
    {
        using BASE = detail::OptimPolishBase< GradientDescent< FN, DFN, ARG_T, MODE_T >, FN, DFN, ARG_T, MODE_T >; /**< Base class alias for
                                                                                                                      readability. */

    public:
        using BASE::BASE;

        void operator()()
        {
            ARG_T gradient;    // Should be RESULT_T. Requires traits classes to be defined.
            if constexpr (std::is_same_v< MODE_T, Minimize >) {
                gradient = BASE::derivative(BASE::m_guess);
            }
            else {
                gradient = -BASE::derivative(BASE::m_guess);
            }

            std::vector< std::pair< double, double > > points;
            points.emplace_back(0.0, BASE::evaluate(BASE::m_guess));
            points.emplace_back(0.5, BASE::evaluate(BASE::m_guess - gradient * 0.5));
            points.emplace_back(1.0, BASE::evaluate(BASE::m_guess - gradient));
            const auto interp = nxx::interp::makepoly(points);

            if (interp.order() < 2) return;

            auto stepsize = nxx::poly::polysolve(derivativeOf(interp))->front();

            if (stepsize > 1.0) {
                BASE::m_guess -= gradient;
                return;
            }

            if (stepsize < 0.0) {
                BASE::m_guess -= gradient * 0.01;
                return;
            }

            BASE::m_guess -= gradient * stepsize;
        }
    };

    /**
     * @brief Deduction guides for Newton class.
     * Allows the type of Newton class to be deduced from the constructor parameters.
     */
    template< typename FN, typename DERIV, typename ARG_T, typename MODE_T >
    requires IsFloatOrComplexInvocable< FN > && IsFloatOrComplexInvocable< DERIV > && IsFloatOrComplex< ARG_T >
    GradientDescent(FN, DERIV, ARG_T, MODE_T) -> GradientDescent< FN, DERIV, ARG_T, MODE_T >;

    template< typename FN, typename DERIV, typename ARG_T, typename MODE_T = Minimize >
    requires IsFloatOrComplexInvocable< FN > && IsFloatOrComplexInvocable< DERIV > && IsFloatOrComplex< ARG_T >
    GradientDescent(FN, DERIV, ARG_T) -> GradientDescent< FN, DERIV, ARG_T, MODE_T >;

    // =================================================================================================================
    //
    // 888b      88
    // 8888b     88                                    ,d
    // 88 `8b    88                                    88
    // 88  `8b   88   ,adPPYba,  8b      db      d8  MM88MMM  ,adPPYba,   8b,dPPYba,
    // 88   `8b  88  a8P_____88  `8b    d88b    d8'    88    a8"     "8a  88P'   `"8a
    // 88    `8b 88  8PP"""""""   `8b  d8'`8b  d8'     88    8b       d8  88       88
    // 88     `8888  "8b,   ,aa    `8bd8'  `8bd8'      88,   "8a,   ,a8"  88       88
    // 88      `888   `"Ybbd8"'      YP      YP        "Y888  `"YbbdP"'   88       88
    //
    // =================================================================================================================

    template< IsFloatOrComplexInvocable FN, IsFloatOrComplexInvocable DFN, IsFloatOrComplex ARG_T, typename MODE_T = Minimize >
    class Newton final : public detail::OptimPolishBase< Newton< FN, DFN, ARG_T, MODE_T >, FN, DFN, ARG_T, MODE_T >
    {
        using BASE =
            detail::OptimPolishBase< Newton< FN, DFN, ARG_T, MODE_T >, FN, DFN, ARG_T, MODE_T >; /**< Base class alias for readability. */

    public:
        using BASE::BASE;

        void operator()()
        {
            using nxx::deriv::derivativeOf;
            using nxx::poly::derivativeOf;

            auto firstDerivative  = [&](ARG_T x) { return BASE::derivative(x); };
            auto secondDerivative = [&](ARG_T x) { return derivativeOf(firstDerivative)(x); };

            BASE::m_guess = BASE::m_guess - firstDerivative(BASE::current()) / secondDerivative(BASE::current());
        }
    };

    /**
     * @brief Deduction guides for Newton class.
     * Allows the type of Newton class to be deduced from the constructor parameters.
     */
    template< typename FN, typename DERIV, typename ARG_T, typename MODE_T >
    requires IsFloatOrComplexInvocable< FN > && IsFloatOrComplexInvocable< DERIV > && IsFloatOrComplex< ARG_T >
    Newton(FN, DERIV, ARG_T, MODE_T) -> Newton< FN, DERIV, ARG_T, MODE_T >;

    template< typename FN, typename DERIV, typename ARG_T, typename MODE_T = Minimize >
    requires IsFloatOrComplexInvocable< FN > && IsFloatOrComplexInvocable< DERIV > && IsFloatOrComplex< ARG_T >
    Newton(FN, DERIV, ARG_T) -> Newton< FN, DERIV, ARG_T, MODE_T >;

    // =================================================================================================================
    //
    //    ad88          88     ad88                                   88                      88
    //   d8"            88    d8"                              ,d     ""                      ""
    //   88             88    88                               88
    // MM88MMM  ,adPPYb,88  MM88MMM  ,adPPYba,   8b,dPPYba,  MM88MMM  88  88,dPYba,,adPYba,   88  888888888   ,adPPYba,
    //   88    a8"    `Y88    88    a8"     "8a  88P'    "8a   88     88  88P'   "88"    "8a  88       a8P"  a8P_____88
    //   88    8b       88    88    8b       d8  88       d8   88     88  88      88      88  88    ,d8P'    8PP"""""""
    //   88    "8a,   ,d88    88    "8a,   ,a8"  88b,   ,a8"   88,    88  88      88      88  88  ,d8"       "8b,   ,aa
    //   88     `"8bbdP"Y8    88     `"YbbdP"'   88`YbbdP"'    "Y888  88  88      88      88  88  888888888   `"Ybbd8"'
    //                                           88
    //                                           88
    //
    // =================================================================================================================

    template< typename SOLVER >
    auto fdfoptimize_impl(SOLVER solver, double eps, int maxiter)
    {
        auto result = solver.current();
        int  iter   = 1;
        while (true) {
            result = solver.current();
            // std::cout << "Iteration " << iter << ": " << result << std::endl;
            if (abs(solver.derivative(result)) < eps * result + eps / 2) break;
            if (iter >= maxiter) break;
            ++iter;
            solver.iterate();
            auto new_guess = solver.current();
            if (abs(new_guess - result) < eps) break;
        }
        return result;
    }

    template< template< typename, typename, typename, typename > class SOLVER_T,
              typename MODE_T = Minimize,
              IsFloatOrComplexInvocable FN_T,
              IsFloatOrComplexInvocable DERIV_T,
              IsFloatOrComplex          GUESS_T,
              IsFloat                   EPS_T  = GUESS_T,
              std::integral             ITER_T = int >
    auto fdfoptimize(FN_T func, DERIV_T deriv, GUESS_T guess, EPS_T eps = epsilon< GUESS_T >(), ITER_T maxiter = iterations< GUESS_T >())
    {
        return fdfoptimize_impl(SOLVER_T(func, deriv, guess, MODE_T {}), eps, maxiter);
    }

    template< template< typename, typename, typename, typename > class SOLVER_T,
              typename MODE_T = Minimize,
              IsFloatOrComplexInvocable FN_T,
              IsFloatOrComplex          GUESS_T,
              IsFloat                   EPS_T  = GUESS_T,
              std::integral             ITER_T = int >
    auto fdfoptimize(FN_T func, GUESS_T guess, EPS_T eps = epsilon< GUESS_T >(), ITER_T maxiter = iterations< GUESS_T >())
    {
        using nxx::deriv::derivativeOf;
        using nxx::poly::derivativeOf;
        return fdfoptimize_impl(SOLVER_T(func, derivativeOf(func), guess, MODE_T {}), eps, maxiter);
    }

}    // namespace nxx::optim

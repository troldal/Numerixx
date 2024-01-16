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

#ifndef NUMERIXX_MULTIROOTS_IMPL_HPP
#define NUMERIXX_MULTIROOTS_IMPL_HPP

// ===== Numerixx Includes
#include "MultiDerivatives.hpp"
#include <Constants.hpp>
#include <Deriv.hpp>
#include <Poly.hpp>
#include <Roots.hpp>

// ===== External Includes
#include <blaze/Blaze.h>

// ===== Standard Library Includes
#include <iostream>
#include <stdexcept>

namespace nxx::multiroots
{
    template< typename RES_T, typename PARAM_T >
    class DMultiNewton;

    /*
     * Forward declaration of the PolishingTraits class.
     */
    template< typename SOLVER >
    struct MultirootsSolverTraits;

    /*
     * Specialization of the PolishingTraits class for Newton<FN, DFN>
     */
    template< typename RES_T, typename PARAM_T >
    struct MultirootsSolverTraits< DMultiNewton< RES_T, PARAM_T > >
    {
        using return_type = RES_T;
        using param_type  = PARAM_T;
    };

    // ========================================================================
    // MULTIROOT-FINDING WITH DERIVATIVES
    // ========================================================================

    template< typename DERIVED, typename RES_T, typename PARAM_T >
    // requires IsMultirootFunction< typename MultirootsSolverTraits< SOLVER >::function_type >
    class MultirootBase
    {
        friend DERIVED;

    public:
        static constexpr bool IsMultirootSolver = true;

        using return_type = RES_T;
        using param_type  = PARAM_T;

    protected:
        ~MultirootBase() = default;

    private:
        MultiFunctionArray< RES_T, PARAM_T > m_functions {};

        using RT = blaze::DynamicVector< RES_T >;
        RT m_guess; /**< The current root estimate. */

        auto size() const { return m_functions.size(); }

    public:
        template< typename ARR >
        explicit MultirootBase(const MultiFunctionArray< RES_T, PARAM_T >& functions, const ARR& guess)
            : m_functions { functions },
              m_guess { RT(guess.size()) }
        {
            std::copy(guess.begin(), guess.end(), m_guess.begin());
        }

        template< typename ARG_T >
        explicit MultirootBase(const MultiFunctionArray< RES_T, PARAM_T >& functions, std::initializer_list< ARG_T > guess)
            : m_functions { functions },
              m_guess { RT(guess.size()) }
        {
            std::copy(guess.begin(), guess.end(), m_guess.begin());
        }

                       MultirootBase(const MultirootBase&)     = default; /**< Default copy constructor. */
                       MultirootBase(MultirootBase&&) noexcept = default; /**< Default move constructor. */
        MultirootBase& operator=(const MultirootBase&)         = default; /**< Default copy assignment operator. */
        MultirootBase& operator=(MultirootBase&&) noexcept     = default; /**< Default move assignment operator. */

        template< typename ARR >
        //        requires std::floating_point< typename nxx::deriv::impl::VectorTraits< ARR >::value_type > ||
        //                 std::floating_point< typename ARR::value_type >
        auto evaluate(ARR values)
        {
            return m_functions.template eval< blaze::DynamicVector >(values);
        }

        template< template< typename... > typename OUT_T >
        [[nodiscard]]
        OUT_T< RES_T > evaluate() const
        {
            return m_functions.template eval< OUT_T >(m_guess);
        }

        [[nodiscard]]
        auto current() const
        {
            return m_guess;
        }

        template< template< typename... > typename OUT_T >
        [[nodiscard]]
        OUT_T< RES_T > current() const
        {
            return OUT_T< RES_T >(m_guess.begin(), m_guess.end());
        }

        void iterate() { std::invoke(static_cast< DERIVED& >(*this)); }
    };

    template< typename RES_T, typename PARAM_T >
    // requires IsMultirootFunction< FunctionType >
    class DMultiNewton final : public MultirootBase< DMultiNewton< RES_T, PARAM_T >, RES_T, PARAM_T >
    {
        using BASE = MultirootBase< DMultiNewton< RES_T, PARAM_T >, RES_T, PARAM_T >;

    public:
        using BASE::BASE;

        void operator()()
        {
            using namespace blaze;
            using namespace nxx::deriv;

            // Solve the linear system J * dx = -f(x) (solving for dx) and update the root estimate (x_new = x_old + dx).
            BASE::m_guess += solve(jacobian(BASE::m_functions, BASE::m_guess), -BASE::evaluate(BASE::m_guess));
        }
    };

    // deduction guide:
    template< typename RES_T, typename PARAM_T, typename ARR_T >
    DMultiNewton(MultiFunctionArray< RES_T, PARAM_T >, ARR_T) -> DMultiNewton< RES_T, PARAM_T >;

    template< typename RES_T, typename PARAM_T, typename ARG_T >
    DMultiNewton(MultiFunctionArray< RES_T, PARAM_T >, std::initializer_list< ARG_T >) -> DMultiNewton< RES_T, PARAM_T >;


    template< typename RES_T, typename PARAM_T >
    class SteepestDescent final : public MultirootBase< SteepestDescent< RES_T, PARAM_T >, RES_T, PARAM_T >
    {
        using BASE = MultirootBase< SteepestDescent< RES_T, PARAM_T >, RES_T, PARAM_T >;

    public:
        using BASE::BASE;

        void operator()()
        {
            auto gradient  = computeGradient(BASE::m_guess);
            auto direction = 1.0 / norm(gradient) * gradient;
            auto stepsize  = computeStepSize(direction);
            BASE::m_guess  = computeGuess(BASE::m_guess, direction, stepsize);
        }

    private:
        auto computeGFunction(const auto& point)
        {
            auto  f_eval = BASE::evaluate(point);
            RES_T result = std::accumulate(f_eval.begin(), f_eval.end(), RES_T {}, [](RES_T sum, const auto& f) { return sum + f * f; });
            return result;
        }

        blaze::DynamicVector< RES_T > computeGradient(const blaze::DynamicVector< RES_T >& point)
        {
            using namespace blaze;
            using namespace nxx::deriv;

            auto J = jacobian(BASE::m_functions, point);
            auto F = BASE::evaluate(point);
            return 2 * trans(J) * F;
        }

        RES_T computeStepSize(const blaze::DynamicVector< RES_T >& direction)
        {
            using namespace blaze;
            using namespace nxx::poly;

            const RES_T a1 = 0;
            const RES_T a2 = 0.5;
            const RES_T a3 = 1.0;

            const DynamicVector< RES_T > arg1 = BASE::m_guess; // - a1 * direction;
            const DynamicVector< RES_T > arg2 = BASE::m_guess - a2 * direction;
            const DynamicVector< RES_T > arg3 = BASE::m_guess - direction; // - a3 * direction;

            const RES_T g1 = computeGFunction(arg1);
            const RES_T g2 = computeGFunction(arg2);
            const RES_T g3 = computeGFunction(arg3);

            RES_T stepsize = 0.001;    // A default small step size

            if (g3 < g1) {

                RES_T               h1 = (g2 - g1) / (a2 - a1);
                RES_T               h2 = (g3 - g2) / (a3 - a2);
                RES_T               h3 = (h2 - h1) / (a3 - a1);
                Polynomial< RES_T > P  = { g1, (h1 - 0.5 * h3), h3 };

                auto P_prime = derivativeOf(P);
                auto roots   = poly::polysolve(P_prime);

                if (!roots->empty()) {
                    stepsize = roots->front();    // Take the first root as the optimal step size
                    if (stepsize < 0 || stepsize > 1) {
                        stepsize = 0.001;    // If the root is not in the expected range, revert to default
                    }
                }
            }

            return stepsize;
        }

        blaze::DynamicVector< RES_T > computeGuess(const blaze::DynamicVector< RES_T >& point, const blaze::DynamicVector< RES_T >& direction, RES_T stepsize)
        {
            using namespace blaze;
            DynamicVector< RES_T > new_guess = point - stepsize * direction;
            auto                   g_new     = computeGFunction(new_guess);
            auto                   g_current = computeGFunction(point);

            if (g_new >= g_current) {
                return point;
                // throw std::runtime_error("Steepest descent failed to improve the solution.");
            }

            return new_guess;
        }
    };

    // Deduction guides (similar to those for DMultiNewton)
    template< typename RES_T, typename PARAM_T, typename ARR_T >
    SteepestDescent(MultiFunctionArray< RES_T, PARAM_T >, ARR_T) -> SteepestDescent< RES_T, PARAM_T >;

    template< typename RES_T, typename PARAM_T, typename ARG_T >
    SteepestDescent(MultiFunctionArray< RES_T, PARAM_T >, std::initializer_list< ARG_T >) -> SteepestDescent< RES_T, PARAM_T >;

    template< typename SOLVER >
    //    requires requires(SOLVER solver, typename SOLVER::function_return_type guess) {
    //                 // clang-format off
    //                 { solver.evaluate(0.0) } -> std::floating_point;
    //                 { solver.init(guess) };
    //                 { solver.iterate() };
    //                 // clang-format on
    //             }
    inline auto multisolve(SOLVER                                      solver,
                           std::vector< typename SOLVER::return_type > guess,
                           typename SOLVER::return_type                eps     = nxx::EPS,
                           int                                         maxiter = nxx::MAXITER)
    {
        // using ET = std::runtime_error;
        // using RT = tl::expected< typename SOLVER::return_type, ET >;
        using RT = std::vector< typename SOLVER::return_type >;

        auto result = solver.current();

        // Check for NaN or Inf.
        if (!std::isfinite(norm(solver.evaluate(result)))) {
            throw std::runtime_error("Non-finite result!");
            // result = tl::make_unexpected(ET("Invalid initial guess!", RootErrorType::NumericalError, result.value()));
            // return result;
        }

        // Begin iteration loop.
        int iter = 1;
        while (true) {
            result = solver.current();

            // Check for NaN or Inf
            if (!std::isfinite(norm(solver.evaluate(result)))) {
                throw std::runtime_error("Non-finite result!");
                // result = tl::make_unexpected(ET("Non-finite result!", RootErrorType::NumericalError, result.value(), iter));
                // break;
            }

            // Check for convergence
            if (norm(solver.evaluate(result)) < eps) {
                std::cout << "Converged after " << iter << " iterations." << std::endl;
                break;
            }

            if (iter >= maxiter) {
                std::cout << "Maximum number of iterations exceeded!" << std::endl;
                break;
            }

            // Perform one iteration
            ++iter;
            solver.iterate();
        }

        return result;
    }

}    // namespace nxx::multiroots

#endif    // NUMERIXX_MULTIROOTS_IMPL_HPP

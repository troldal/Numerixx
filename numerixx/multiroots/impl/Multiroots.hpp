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
#include <Interp.hpp>
#include <Optim.hpp>

// ===== External Includes
#include <blaze/Blaze.h>

// ===== Standard Library Includes
#include <iostream>
#include <stdexcept>

namespace nxx::multiroots
{
    template< IsFloat RES_T, IsFloat PARAM_T, IsFloat ARG_T >
    class MultiNewton;

    template< IsFloat RES_T, IsFloat PARAM_T, IsFloat ARG_T >
    class SteepestDescent;

    /*
     * Forward declaration of the PolishingTraits class.
     */
    template< typename SOLVER >
    struct MultirootsSolverTraits;

    /*
     * Specialization of the PolishingTraits class for Newton<FN, DFN>
     */
    template< IsFloat RES_T, IsFloat PARAM_T, IsFloat ARG_T >
    struct MultirootsSolverTraits< MultiNewton< RES_T, PARAM_T, ARG_T > >
    {
        using return_type = RES_T;
        using param_type  = PARAM_T;
        using arg_type    = ARG_T;
    };

    template< IsFloat RES_T, IsFloat PARAM_T, IsFloat ARG_T >
    struct MultirootsSolverTraits< SteepestDescent< RES_T, PARAM_T, ARG_T > >
    {
        using return_type = RES_T;
        using param_type  = PARAM_T;
        using arg_type    = ARG_T;
    };

    // =================================================================================================================
    //
    // 88b           d88               88           88  88888888ba
    // 888b         d888               88    ,d     ""  88      "8b                            ,d
    // 88`8b       d8'88               88    88         88      ,8P                            88
    // 88 `8b     d8' 88  88       88  88  MM88MMM  88  88aaaaaa8P'  ,adPPYba,    ,adPPYba,  MM88MMM  ,adPPYba,
    // 88  `8b   d8'  88  88       88  88    88     88  88""""88'   a8"     "8a  a8"     "8a   88     I8[    ""
    // 88   `8b d8'   88  88       88  88    88     88  88    `8b   8b       d8  8b       d8   88      `"Y8ba,
    // 88    `888'    88  "8a,   ,a88  88    88,    88  88     `8b  "8a,   ,a8"  "8a,   ,a8"   88,    aa    ]8I
    // 88     `8'     88   `"YbbdP'Y8  88    "Y888  88  88      `8b  `"YbbdP"'    `"YbbdP"'    "Y888  `"YbbdP"'
    //
    // =================================================================================================================

    namespace detail
    {
        /**
         * @class MultirootBase
         * @brief Template class for base multi-root solver.
         *
         * @tparam DERIVED The derived class that inherits from MultirootBase.
         * @tparam RES_T The return type of the functions used in the solver.
         * @tparam PARAM_T The parameter type of the functions used in the solver.
         * @tparam ARG_T The type of the arguments used for initial guess.
         */
        template< typename DERIVED, IsFloat RES_T, IsFloat PARAM_T, IsFloat ARG_T >
        class MultirootBase
        {
            friend DERIVED;

        public:
            static constexpr bool IsMultirootSolver = true;    ///< Flag indicating this is a multi-root solver.

            using return_type = RES_T;                                   ///< Alias for return type.
            using param_type  = PARAM_T;                                 ///< Alias for parameter type.
            using arg_type    = ARG_T;                                   ///< Alias for argument type.
            using RETURN_T    = blaze::DynamicVector< RES_T >;           ///< Alias for Blaze dynamic vector.
            using FUNCTION_T  = MultiFunctionArray< RES_T, PARAM_T >;    ///< Alias for function array type.

        protected:
            /**
             * @brief Protected destructor to prevent direct instantiation of base class.
             */
            ~MultirootBase() = default;

        private:
            FUNCTION_T m_functions {};    ///< Storage for function objects.
            RETURN_T   m_guess {};        ///< Initial guess for the roots.

        public:
            /**
             * @brief Constructor initializing the multi-root solver with functions and an initial guess.
             * @tparam ARR Container type for the initial guess.
             * @param functions A MultiFunctionArray object containing the functions for root solving.
             * @param guess Initial guess for the roots.
             */
            template< typename ARR >    // Use IsContainer instead??
            explicit MultirootBase(const MultiFunctionArray< RES_T, PARAM_T >& functions, const ARR& guess)
                : m_functions { functions },
                  m_guess { RETURN_T(guess.size()) }
            {
                std::copy(guess.begin(), guess.end(), m_guess.begin());
            }

            /**
             * @brief Constructor initializing the multi-root solver with functions and an initializer list as the guess.
             * @param functions A MultiFunctionArray object containing the functions for root solving.
             * @param guess Initial guess for the roots provided as an initializer list.
             */
            explicit MultirootBase(const MultiFunctionArray< RES_T, PARAM_T >& functions, std::initializer_list< ARG_T > guess)
                : m_functions { functions },
                  m_guess { RETURN_T(guess.size()) }
            {
                std::copy(guess.begin(), guess.end(), m_guess.begin());
            }

            // Rule of five: Default copy/move constructors and assignment operators
            MultirootBase(const MultirootBase&)     = default; /**< Default copy constructor. */
            MultirootBase(MultirootBase&&) noexcept = default; /**< Default move constructor. */

            MultirootBase& operator=(const MultirootBase&)     = default; /**< Default copy assignment operator. */
            MultirootBase& operator=(MultirootBase&&) noexcept = default; /**< Default move assignment operator. */

            /**
             * @brief Evaluates the functions at given values.
             * @tparam ARR Container type for the input values.
             * @param values Container of values for function evaluation.
             * @return A Blaze dynamic vector containing the results of function evaluations.
             */
            template< typename ARR >    // Use IsContainer instead??
            auto evaluate(ARR values)
            {
                return m_functions.template eval< blaze::DynamicVector >(values);
            }

            /**
             * @brief Evaluates the functions using the current guess.
             * @tparam OUT_T Template template parameter for the output container type.
             * @return OUT_T<RES_T> Container of the specified type containing the results of function evaluations.
             */
            template< template< typename... > typename OUT_T >
            OUT_T< RES_T > evaluate() const
            {
                return m_functions.template eval< OUT_T >(m_guess);
            }

            /**
             * @brief Retrieves the current guess.
             * @return A Blaze dynamic vector representing the current guess.
             */
            auto current() const { return m_guess; }

            /**
             * @brief Retrieves the current guess in a specified container type.
             * @tparam OUT_T Template template parameter for the output container type.
             * @return OUT_T<RES_T> Container of the specified type containing the current guess.
             */
            template< template< typename... > typename OUT_T >
            OUT_T< RES_T > current() const
            {
                return OUT_T< RES_T >(m_guess.begin(), m_guess.end());
            }

            /**
             * @brief Iterates the solver once using the derived class's implementation.
             */
            void iterate() { std::invoke(static_cast< DERIVED& >(*this)); }
        };
    }    // namespace detail

    // =================================================================================================================
    //
    // 88b           d88  888b      88
    // 888b         d888  8888b     88                                    ,d
    // 88`8b       d8'88  88 `8b    88                                    88
    // 88 `8b     d8' 88  88  `8b   88   ,adPPYba,  8b      db      d8  MM88MMM  ,adPPYba,   8b,dPPYba,
    // 88  `8b   d8'  88  88   `8b  88  a8P_____88  `8b    d88b    d8'    88    a8"     "8a  88P'   `"8a
    // 88   `8b d8'   88  88    `8b 88  8PP"""""""   `8b  d8'`8b  d8'     88    8b       d8  88       88
    // 88    `888'    88  88     `8888  "8b,   ,aa    `8bd8'  `8bd8'      88,   "8a,   ,a8"  88       88
    // 88     `8'     88  88      `888   `"Ybbd8"'      YP      YP        "Y888  `"YbbdP"'   88       88
    //
    // =================================================================================================================

    /**
     * @class MultiNewton
     * @brief Template class implementing the Newton-Raphson method for multi-root solving.
     *
     * @tparam RES_T The return type of the functions used in the solver.
     * @tparam PARAM_T The parameter type of the functions used in the solver.
     * @tparam ARG_T The type of the arguments used for initial guess.
     */
    template< IsFloat RES_T, IsFloat PARAM_T, IsFloat ARG_T >
    class MultiNewton final : public detail::MultirootBase< MultiNewton< RES_T, PARAM_T, ARG_T >, RES_T, PARAM_T, ARG_T >
    {
        using BASE = detail::MultirootBase< MultiNewton< RES_T, PARAM_T, ARG_T >, RES_T, PARAM_T, ARG_T >;

    public:
        /**
         * @brief Inherits constructors from the base class MultirootBase.
         */
        using BASE::BASE;

        /**
         * @brief Overloaded function call operator implementing the Newton-Raphson iteration step.
         * @details This method updates the current guess by solving the linear system J * dx = -f(x) and updating the root estimate.
         */
        void operator()()
        {
            using namespace blaze;
            using namespace nxx::deriv;

            // Solve the linear system J * dx = -f(x) (solving for dx) and update the root estimate (x_new = x_old + dx).
            BASE::m_guess += solve(jacobian(BASE::m_functions, BASE::m_guess), -BASE::evaluate(BASE::m_guess));
        }
    };

    /**
     * @brief Deduction guide for MultiNewton with a MultiFunctionArray and an arbitrary container type.
     * @tparam RES_T The return type of the functions.
     * @tparam PARAM_T The parameter type of the functions.
     * @tparam ARR_T The container type for the initial guess.
     */
    template< IsFloat RES_T, IsFloat PARAM_T, typename ARR_T >
    MultiNewton(MultiFunctionArray< RES_T, PARAM_T >, ARR_T) -> MultiNewton< RES_T, PARAM_T, traits::ContainerValueType_t< ARR_T > >;

    /**
     * @brief Deduction guide for MultiNewton with a MultiFunctionArray and an initializer list.
     * @tparam RES_T The return type of the functions.
     * @tparam PARAM_T The parameter type of the functions.
     * @tparam ARG_T The type of the arguments in the initializer list.
     */
    template< IsFloat RES_T, IsFloat PARAM_T, IsFloat ARG_T >
    MultiNewton(MultiFunctionArray< RES_T, PARAM_T >, std::initializer_list< ARG_T >) -> MultiNewton< RES_T, PARAM_T, ARG_T >;

    // =================================================================================================================
    //
    //
    //  ad88888ba   88888888ba,
    // d8"     "8b  88      `"8b                                                                ,d
    // Y8,          88        `8b                                                               88
    // `Y8aaaaa,    88         88   ,adPPYba,  ,adPPYba,   ,adPPYba,   ,adPPYba,  8b,dPPYba,  MM88MMM
    //   `"""""8b,  88         88  a8P_____88  I8[    ""  a8"     ""  a8P_____88  88P'   `"8a   88
    //         `8b  88         8P  8PP"""""""   `"Y8ba,   8b          8PP"""""""  88       88   88
    // Y8a     a8P  88      .a8P   "8b,   ,aa  aa    ]8I  "8a,   ,aa  "8b,   ,aa  88       88   88,
    //  "Y88888P"   88888888Y"'     `"Ybbd8"'  `"YbbdP"'   `"Ybbd8"'   `"Ybbd8"'  88       88   "Y888
    //
    // =================================================================================================================

    /**
     * @class SteepestDescent
     * @brief Template class implementing the Steepest Descent method for multi-root solving.
     *
     * @tparam RES_T The return type of the functions used in the solver.
     * @tparam PARAM_T The parameter type of the functions used in the solver.
     * @tparam ARG_T The type of the arguments used for initial guess.
     */
    template< IsFloat RES_T, IsFloat PARAM_T, IsFloat ARG_T >
    class SteepestDescent final : public detail::MultirootBase< SteepestDescent< RES_T, PARAM_T, ARG_T >, RES_T, PARAM_T, ARG_T >
    {
        using BASE = detail::MultirootBase< SteepestDescent< RES_T, PARAM_T, ARG_T >, RES_T, PARAM_T, ARG_T >;

    public:
        /**
         * @brief Inherits constructors from the base class MultirootBase.
         */
        using BASE::BASE;

        /**
         * @brief Overloaded function call operator implementing the Steepest Descent iteration step.
         * @details This method updates the current guess by computing the gradient, direction, and step size,
         *          and then updating the root estimate accordingly.
         */
        void operator()()
        {
            auto gradient  = computeGradient(BASE::m_guess);
            auto direction = 1.0 / norm(gradient) * gradient;
            auto stepsize  = computeStepSize(direction);
            BASE::m_guess  = computeGuess(BASE::m_guess, direction, stepsize);
        }

    private:
        /**
         * @brief Computes the G function for a given point.
         * @param point The point at which to compute the G function.
         * @return The computed G function value.
         */
        auto computeGFunction(const auto& point)
        {
            auto  f_eval = BASE::evaluate(point);
            RES_T result = std::accumulate(f_eval.begin(), f_eval.end(), RES_T {}, [](RES_T sum, const auto& f) { return sum + f * f; });
            return result;
        }

        /**
         * @brief Computes the gradient vector at a given point.
         * @param point The point at which to compute the gradient.
         * @return The gradient vector at the given point.
         */
        blaze::DynamicVector< RES_T > computeGradient(const blaze::DynamicVector< RES_T >& point)
        {
            using namespace blaze;
            using namespace nxx::deriv;

            auto J = jacobian(BASE::m_functions, point);
            auto F = BASE::evaluate(point);
            return 2 * trans(J) * F;
        }

        /**
         * @brief Computes the step size for the steepest descent.
         * @param direction The descent direction.
         * @return The computed step size.
         */
        RES_T computeStepSize(const blaze::DynamicVector< RES_T >& direction)
        {
            using namespace blaze;
            using namespace nxx::poly;

            const RES_T a1 = 0;
            const RES_T a2 = 0.5;
            const RES_T a3 = 1.0;

            const DynamicVector< RES_T > arg1 = BASE::m_guess;    // - a1 * direction;
            const DynamicVector< RES_T > arg2 = BASE::m_guess - a2 * direction;
            const DynamicVector< RES_T > arg3 = BASE::m_guess - direction;    // - a3 * direction;

            const RES_T g1 = computeGFunction(arg1);
            const RES_T g2 = computeGFunction(arg2);
            const RES_T g3 = computeGFunction(arg3);

            // auto func = [&](RES_T a) {return computeGFunction(DynamicVector< RES_T >(BASE::m_guess - a * direction)); };

            // RES_T stepsize = 0.001;    // A default small step size
            //
            // if (g3 < g1) {
            //     RES_T               h1 = (g2 - g1) / (a2 - a1);
            //     RES_T               h2 = (g3 - g2) / (a3 - a2);
            //     RES_T               h3 = (h2 - h1) / (a3 - a1);
            //     Polynomial< RES_T > P  = { g1, (h1 - 0.5 * h3), h3 };
            //
            //     auto P_prime = derivativeOf(P);
            //     if (auto roots = poly::polysolve(P_prime); !roots->empty()) {
            //         stepsize = roots->front();    // Take the first root as the optimal step size
            //         if (stepsize < 0 || stepsize > 1) {
            //             stepsize = 0.001;    // If the root is not in the expected range, revert to default
            //         }
            //     }
            // }


            std::vector<std::pair<double, double>> points;
            points.emplace_back(0.0, g1);
            points.emplace_back(0.5, g2);
            points.emplace_back(1.0, g3);
            auto interp = nxx::interp::makepoly(points);

            if (interp.order() < 2) return 0.0;

            auto stepsize = nxx::poly::polysolve(derivativeOf(interp))->front();

            if (stepsize > 1.0) return 1.0;

            if (stepsize < 0.0) return 0.001;

            // using namespace nxx::optim;
            // RES_T stepsize = optimize<GradientDescent, Minimize>(func, nxx::deriv::derivativeOf(func), 0.5);

            return stepsize;
        }

        /**
         * @brief Computes the next guess in the steepest descent method.
         * @param point The current point.
         * @param direction The descent direction.
         * @param stepsize The step size.
         * @return The next guess point.
         */
        blaze::DynamicVector< RES_T >
            computeGuess(const blaze::DynamicVector< RES_T >& point, const blaze::DynamicVector< RES_T >& direction, RES_T stepsize)
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

    /**
     * @brief Deduction guide for SteepestDescent with a MultiFunctionArray and an arbitrary container type.
     * @tparam RES_T The return type of the functions.
     * @tparam PARAM_T The parameter type of the functions.
     * @tparam ARR_T The container type for the initial guess.
     */
    template< IsFloat RES_T, IsFloat PARAM_T, typename ARR_T >
    SteepestDescent(MultiFunctionArray< RES_T, PARAM_T >, ARR_T)
        -> SteepestDescent< RES_T, PARAM_T, traits::ContainerValueType_t< ARR_T > >;

    /**
     * @brief Deduction guide for SteepestDescent with a MultiFunctionArray and an initializer list.
     * @tparam RES_T The return type of the functions.
     * @tparam PARAM_T The parameter type of the functions.
     * @tparam ARG_T The type of the arguments in the initializer list.
     */
    template< IsFloat RES_T, IsFloat PARAM_T, IsFloat ARG_T >
    SteepestDescent(MultiFunctionArray< RES_T, PARAM_T >, std::initializer_list< ARG_T >) -> SteepestDescent< RES_T, PARAM_T, ARG_T >;

    // =================================================================================================================
    //
    //                                  88           88                          88
    //                                  88    ,d     ""                          88
    //                                  88    88                                 88
    // 88,dPYba,,adPYba,   88       88  88  MM88MMM  88  ,adPPYba,   ,adPPYba,   88  8b       d8   ,adPPYba,
    // 88P'   "88"    "8a  88       88  88    88     88  I8[    ""  a8"     "8a  88  `8b     d8'  a8P_____88
    // 88      88      88  88       88  88    88     88   `"Y8ba,   8b       d8  88   `8b   d8'   8PP"""""""
    // 88      88      88  "8a,   ,a88  88    88,    88  aa    ]8I  "8a,   ,a8"  88    `8b,d8'    "8b,   ,aa
    // 88      88      88   `"YbbdP'Y8  88    "Y888  88  `"YbbdP"'   `"YbbdP"'   88      "8"       `"Ybbd8"'
    //
    //
    // =================================================================================================================

    namespace detail
    {

        /**
         * @brief Solves multi-root problems using a given solver.
         * @details This function iteratively applies the solver to find the roots of a system of equations,
         *          and checks for convergence or non-finite results.
         * @tparam SOLVER The type of the solver, which must conform to the MultirootSolver interface.
         * @param solver An instance of the solver.
         * @param eps The convergence tolerance.
         * @param maxiter The maximum number of iterations allowed.
         * @return tl::expected<RESULT_T, ERROR_T> A tl::expected object containing the result vector
         *         or an error if the solution did not converge or if a non-finite result was encountered.
         */
        template< typename SOLVER >
        requires SOLVER::IsMultirootSolver
        auto multisolve_impl(SOLVER solver, IsFloat auto eps, std::integral auto maxiter)
        {
            using ERROR_T  = std::runtime_error;
            using RESULT_T = blaze::DynamicVector< typename SOLVER::return_type >;
            using RETURN_T = tl::expected< RESULT_T, ERROR_T >;

            RETURN_T result = solver.current();

            // Check for NaN or Inf.
            if (!std::isfinite(norm(solver.evaluate(*result)))) {
                result = tl::make_unexpected(ERROR_T("Non-finite result!"));
                return result;
            }

            // Begin iteration loop.
            int iter = 1;
            while (true) {
                result = solver.current();

                // Check for NaN or Inf
                if (!std::isfinite(norm(solver.evaluate(*result)))) {
                    result = tl::make_unexpected(ERROR_T("Non-finite result!"));
                    return result;
                }

                // Check for convergence
                if (norm(solver.evaluate(*result)) < eps) {
                    std::cout << "Converged after " << iter << " iterations." << std::endl;
                    break;
                }

                if (iter >= maxiter) {
                    std::cout << "Maximum number of iterations exceeded: " << maxiter << std::endl;
                    break;
                }

                // Perform one iteration
                ++iter;
                solver.iterate();
            }

            return result;
        }
    }    // namespace detail

    /**
     * @brief Solves multi-root problems using a specified solver template, functions, and an initial guess.
     * @details This function template creates a solver instance and calls the multisolve_impl function to find roots.
     * @tparam SOLVER_T Template class of the solver.
     * @tparam RES_T The return type of the functions used in the solver.
     * @tparam PARAM_T The parameter type of the functions used in the solver.
     * @tparam ARR_T Container type for the initial guess.
     * @tparam EPS_T Floating point type for the convergence tolerance.
     * @tparam ITER_T Integral type for the maximum number of iterations.
     * @param functions A MultiFunctionArray object containing the functions for root solving.
     * @param guess Initial guess for the roots.
     * @param eps Convergence tolerance.
     * @param maxiter Maximum number of iterations.
     * @return An instance of tl::expected containing the result or an error.
     */
    template< template< typename, typename, typename > class SOLVER_T,
              IsFloat RES_T,
              IsFloat PARAM_T,
              typename ARR_T,    // change to IsContainer??
              IsFloat       EPS_T  = traits::ContainerValueType_t< ARR_T >,
              std::integral ITER_T = int >
    auto multisolve(MultiFunctionArray< RES_T, PARAM_T > functions,
                    ARR_T                                guess,
                    EPS_T                                eps     = epsilon< traits::ContainerValueType_t< ARR_T > >(),
                    ITER_T                               maxiter = iterations< traits::ContainerValueType_t< ARR_T > >())
    {
        using SOLVER = SOLVER_T< RES_T, PARAM_T, traits::ContainerValueType_t< ARR_T > >;
        return detail::multisolve_impl(SOLVER(functions, guess), eps, maxiter);
    }

    /**
     * @brief Solves multi-root problems using a specified solver template, functions, and an initializer list as the guess.
     * @details This function template creates a solver instance and calls the multisolve_impl function to find roots.
     * @tparam SOLVER_T Template class of the solver.
     * @tparam RES_T The return type of the functions used in the solver.
     * @tparam PARAM_T The parameter type of the functions used in the solver.
     * @tparam ARG_T The type of the arguments in the initializer list.
     * @tparam EPS_T Floating point type for the convergence tolerance.
     * @tparam ITER_T Integral type for the maximum number of iterations.
     * @param functions A MultiFunctionArray object containing the functions for root solving.
     * @param guess Initial guess for the roots provided as an initializer list.
     * @param eps Convergence tolerance.
     * @param maxiter Maximum number of iterations.
     * @return An instance of tl::expected containing the result or an error.
     */
    template< template< typename, typename, typename > class SOLVER_T,
              IsFloat       RES_T,
              IsFloat       PARAM_T,
              IsFloat       ARG_T,
              IsFloat       EPS_T  = ARG_T,
              std::integral ITER_T = int >
    auto multisolve(MultiFunctionArray< RES_T, PARAM_T > functions,
                    std::initializer_list< ARG_T >       guess,
                    EPS_T                                eps     = epsilon< ARG_T >(),
                    ITER_T                               maxiter = iterations< ARG_T >())
    {
        using SOLVER = SOLVER_T< RES_T, PARAM_T, ARG_T >;
        return detail::multisolve_impl(SOLVER(functions, guess), eps, maxiter);
    }

}    // namespace nxx::multiroots

#endif    // NUMERIXX_MULTIROOTS_IMPL_HPP

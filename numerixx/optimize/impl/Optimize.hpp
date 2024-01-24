//
// Created by kenne on 13/11/2023.
//

#pragma once

#include <cmath>
#include <functional>

#include <Concepts.hpp>
#include <Deriv.hpp>
#include <Poly.hpp>
#include <Interp.hpp>
#include <Roots.hpp>

namespace nxx::optim
{

    struct Minimize
    {
    };
    struct Maximize
    {
    };

    template< typename DERIVED, typename FUNCTION_T, typename DERIV_T, typename ARG_T, typename MODE_T >
    class OptimDerivBase
    {
        friend DERIVED;

    public:
        static constexpr bool IsDerivativeOptimizer = true;

        using FUNCT_RES_T = std::invoke_result_t< FUNCTION_T, ARG_T >;
        using DERIV_RES_T = std::invoke_result_t< DERIV_T, ARG_T >;
        using RESULT_T    = std::common_type_t< FUNCT_RES_T, DERIV_RES_T >;

    protected:
        ~OptimDerivBase() = default;

    private:
        FUNCTION_T m_func {};
        DERIV_T    m_deriv {};
        RESULT_T   m_guess;
        MODE_T     m_mode;

    public:
        OptimDerivBase(FUNCTION_T objective, DERIV_T derivative, ARG_T guess, MODE_T mode = {})
            : m_func { objective },
              m_deriv { derivative },
              m_guess { guess },
              m_mode(mode)
        {}



        OptimDerivBase(const OptimDerivBase& other)     = default; /**< Default copy constructor. */
        OptimDerivBase(OptimDerivBase&& other) noexcept = default; /**< Default move constructor. */

        OptimDerivBase& operator=(const OptimDerivBase& other)     = default; /**< Default copy assignment operator. */
        OptimDerivBase& operator=(OptimDerivBase&& other) noexcept = default; /**< Default move assignment operator. */

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

    template< IsFloatOrComplexInvocable FN, IsFloatOrComplexInvocable DFN, IsFloatOrComplex ARG_T, typename MODE_T = Minimize >
    class GradientDescent final : public OptimDerivBase< GradientDescent< FN, DFN, ARG_T, MODE_T >, FN, DFN, ARG_T, MODE_T >
    {
        using BASE =
            OptimDerivBase< GradientDescent< FN, DFN, ARG_T, MODE_T >, FN, DFN, ARG_T, MODE_T >; /**< Base class alias for readability. */

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

            std::vector<std::pair<double, double>> points;
            points.emplace_back(0.0, BASE::evaluate(BASE::m_guess));
            points.emplace_back(0.5, BASE::evaluate(BASE::m_guess - gradient * 0.5));
            points.emplace_back(1.0, BASE::evaluate(BASE::m_guess - gradient));
            auto interp = nxx::interp::makepoly(points);

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

    // template< typename ALGO, typename Function >
    // double optimize(ALGO algorithm, Function func, double initialGuess)
    // {
    //     return algorithm(func, initialGuess);
    // }
    //
    // template< typename ALGO, typename Function >
    // auto optimizationOf(ALGO algorithm, Function func)
    // {
    //     return [=](double initialGuess) { return algorithm(func, initialGuess); };
    // }

    template< typename SOLVER >
    auto optimize_impl(SOLVER solver, double eps, int maxiter)
    {
        auto result = solver.current();
        int  iter   = 1;
        while (true) {

            result = solver.current();
            std::cout << "Iteration " << iter << ": " << result << std::endl;
            if (abs(solver.derivative(result)) < eps) break;
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
              IsFloat                   EPS_T = GUESS_T,
              std::integral             ITER_T = int >
    auto optimize(FN_T func, DERIV_T deriv, GUESS_T guess, EPS_T eps = epsilon< GUESS_T >(), ITER_T maxiter = iterations< GUESS_T >())
    {
        return optimize_impl(SOLVER_T(func, deriv, guess, MODE_T{}), eps, maxiter);
    }

}    // namespace nxx::optim


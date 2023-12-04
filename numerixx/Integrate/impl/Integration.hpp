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

#ifndef NUMERIXX_INTEGRATION_HPP
#define NUMERIXX_INTEGRATION_HPP

// ===== Local Includes
#include "IntegrationTraits.hpp"
#include "IntegrationError.hpp"
#include "IntegrationValidation.hpp"

// ===== Numerixx Includes
#include <Concepts.hpp>
#include <Constants.hpp>
#include <Error.hpp>

// ===== External Includes
#include <tl/expected.hpp>
#include <boost/multi_array.hpp>

// ===== Standard Library Includes
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <random>
#include <type_traits>
#include <span>

namespace nxx::integrate
{
    namespace detail
    {
        // =============================================================================================================
        //
        // 88888888ba
        // 88      "8b
        // 88      ,8P
        // 88aaaaaa8P'  ,adPPYYba,  ,adPPYba,   ,adPPYba,
        // 88""""""8b,  ""     `Y8  I8[    ""  a8P_____88
        // 88      `8b  ,adPPPPP88   `"Y8ba,   8PP"""""""
        // 88      a8P  88,    ,88  aa    ]8I  "8b,   ,aa
        // 88888888P"   `"8bbdP"Y8  `"YbbdP"'   `"Ybbd8"'
        //
        // =============================================================================================================

        template<typename DERIVED, typename FUNCTION_T, typename ARG_T>
            requires std::same_as< typename IntegrationTraits< DERIVED >::FUNCTION_T, FUNCTION_T > &&
                     nxx::IsFloatInvocable< FUNCTION_T > &&
                     nxx::IsFloat< ARG_T > &&
                     nxx::IsFloat< typename IntegrationTraits< DERIVED >::RETURN_T >
        class IntegrationBase
        {
            friend DERIVED;

        public:
            static constexpr bool IsIntegrationSolver = true; /**< Flag indicating the class is a bracketing solver. */

            using RESULT_T = std::invoke_result_t< FUNCTION_T, ARG_T >; /**< Result type of the function. */
            using BOUNDS_T = std::pair< ARG_T, ARG_T >;                 /**< Type for representing the bounds around the root. */

        protected:
            ~IntegrationBase() = default; /**< Protected destructor to prevent direct instantiation. */

        private:
            FUNCTION_T m_func{};     /**< The function object to find the root for. */
            BOUNDS_T   m_bounds{};   /**< Holds the current bounds around the root. */
            RESULT_T   m_estimate{}; /**< Holds the current estimate of the root. */
            ARG_T      m_interval{}; /**< Holds the current interval size. */

        public:
            IntegrationBase(FUNCTION_T objective, const IsFloatStruct auto& bounds)
                : m_func{ std::move(objective) } { init(bounds); }

            template<size_t N> requires (N == 2)
            IntegrationBase(FUNCTION_T objective, const ARG_T (&bounds)[N])
                : m_func{ std::move(objective) }
            {
                auto bnds = std::span(bounds, N);
                init(std::pair{ bnds.front(), bnds.back() });
            }

            IntegrationBase(const IntegrationBase& other)                = default; /**< Default copy constructor. */
            IntegrationBase(IntegrationBase&& other) noexcept            = default; /**< Default move constructor. */
            IntegrationBase& operator=(const IntegrationBase& other)     = default; /**< Default copy assignment operator. */
            IntegrationBase& operator=(IntegrationBase&& other) noexcept = default; /**< Default move assignment operator. */

            template<IsFloatStruct STRUCT_T>
            void init(STRUCT_T bounds)
            {
                auto [lower, upper] = bounds;
                m_bounds            = BOUNDS_T{ lower, upper };
                m_interval          = upper - lower;
                m_estimate          = m_interval * (evaluate(lower) + evaluate(upper)) / 2;
            }

            RESULT_T evaluate(ARG_T value) { return m_func(value); }

            RESULT_T current() const { return m_estimate; }

            void iterate() { std::invoke(static_cast< DERIVED& >(*this)); }
        };
    } // namespace detail

    // =================================================================================================================
    //
    // 888888888888                                                                      88           88
    //      88                                                                           ""           88
    //      88                                                                                        88
    //      88  8b,dPPYba,  ,adPPYYba,  8b,dPPYba,    ,adPPYba,  888888888   ,adPPYba,   88   ,adPPYb,88
    //      88  88P'   "Y8  ""     `Y8  88P'    "8a  a8P_____88       a8P"  a8"     "8a  88  a8"    `Y88
    //      88  88          ,adPPPPP88  88       d8  8PP"""""""    ,d8P'    8b       d8  88  8b       88
    //      88  88          88,    ,88  88b,   ,a8"  "8b,   ,aa  ,d8"       "8a,   ,a8"  88  "8a,   ,d88
    //      88  88          `"8bbdP"Y8  88`YbbdP"'    `"Ybbd8"'  888888888   `"YbbdP"'   88   `"8bbdP"Y8
    //                                  88
    //                                  88
    // =================================================================================================================


    template<IsFloatInvocable FN, IsFloat ARG_T = double>
    class Trapezoid final : public detail::IntegrationBase< Trapezoid< FN, ARG_T >, FN, ARG_T >
    {
        using BASE = detail::IntegrationBase< Trapezoid< FN, ARG_T >, FN, ARG_T >; /**< Base class alias for readability. */

    public:
        using BASE::BASE; /**< Inherits constructors from BracketingBase. */
        using RESULT_T = std::invoke_result_t< FN, ARG_T >;

        static const inline std::string SolverName = "Trapezoid";

        int m_iter{ 1 };

        void operator()()
        {
            const auto& [lower, upper] = BASE::m_bounds;

            BASE::m_interval /= 2;                              // Halve the step size in each iteration
            const int numOfMidpoints = std::pow(2, m_iter - 1); // Calculate the number of midpoints for this iteration

            // Create a vector of midpoints
            std::vector< RESULT_T > midPoints(numOfMidpoints);
            std::generate(midPoints.begin(),
                          midPoints.end(),
                          [n = 1, lower, h = BASE::m_interval]() mutable { return lower + (2 * n++ - 1) * h; });

            // Calculate the sum of function values at the midpoints
            const RESULT_T sum =
                std::transform_reduce(midPoints.begin(),
                                      midPoints.end(),
                                      ARG_T(0.0),
                                      std::plus(),
                                      [&](RESULT_T x) { return BASE::evaluate(x); });

            // Update the integral estimate
            BASE::m_estimate = BASE::m_estimate / 2 + BASE::m_interval * sum;
            m_iter++;
        }
    };

    template<typename FN, typename BOUNDS_T>
        requires IsFloatInvocable< FN > && IsFloatStruct< BOUNDS_T >
    Trapezoid(FN, BOUNDS_T) -> Trapezoid< FN, StructCommonType_t< BOUNDS_T > >;

    template<typename FN, typename ARG_T>
        requires IsFloatInvocable< FN > && IsFloat< ARG_T >
    Trapezoid(FN, std::initializer_list< ARG_T >) -> Trapezoid< FN, ARG_T >;


    // =================================================================================================================
    // 88888888ba                                   88
    // 88      "8b                                  88
    // 88      ,8P                                  88
    // 88aaaaaa8P'  ,adPPYba,   88,dPYba,,adPYba,   88,dPPYba,    ,adPPYba,  8b,dPPYba,   ,adPPYb,d8
    // 88""""88'   a8"     "8a  88P'   "88"    "8a  88P'    "8a  a8P_____88  88P'   "Y8  a8"    `Y88
    // 88    `8b   8b       d8  88      88      88  88       d8  8PP"""""""  88          8b       88
    // 88     `8b  "8a,   ,a8"  88      88      88  88b,   ,a8"  "8b,   ,aa  88          "8a,   ,d88
    // 88      `8b  `"YbbdP"'   88      88      88  8Y"Ybbd8"'    `"Ybbd8"'  88           `"YbbdP"Y8
    //                                                                                    aa,    ,88
    //                                                                                     "Y8bbdP"
    // =================================================================================================================

    template<IsFloatInvocable FN, IsFloat ARG_T = double>
    class Romberg final : public detail::IntegrationBase< Romberg< FN, ARG_T >, FN, ARG_T >
    {
        using BASE = detail::IntegrationBase< Romberg< FN, ARG_T >, FN, ARG_T >; /**< Base class alias for readability. */

    public:
        using BASE::BASE; /**< Inherits constructors from BracketingBase. */
        using RESULT_T = std::invoke_result_t< FN, ARG_T >;
        using ARRAY_T = boost::multi_array< RESULT_T, 2 >;

        static const inline std::string SolverName = "Romberg";

        int     m_iter{ 1 };
        ARRAY_T R;

        void operator()()
        {
            using boost::extents;
            using boost::multi_array;

            R.resize(extents[m_iter + 1][m_iter + 1]);

            const auto& [lower, upper] = BASE::m_bounds;

            // Initialize the first element with the basic trapezoidal rule
            if (m_iter == 1) { R[0][0] = (upper - lower) * (BASE::evaluate(lower) + BASE::evaluate(upper)) / 2; }

            // Calculate the step size for this iteration
            const RESULT_T h = (upper - lower) / std::pow(2, m_iter);

            // Trapezoidal rule: Calculate the sum of the function values at the midpoints
            RESULT_T sum = 0.0;
            for (int k = 1; k <= std::pow(2, m_iter - 1); ++k)
                sum += BASE::evaluate(lower + (2 * k - 1) * h);

            // Update the Romberg table for the first column (trapezoidal rule)
            R[m_iter][0] = R[m_iter - 1][0] / 2 + h * sum;

            // Apply the Romberg integration formula
            for (int j = 1; j <= m_iter; ++j)
                // Use the previously computed values to extrapolate to higher orders
                R[m_iter][j] = R[m_iter][j - 1] + (R[m_iter][j - 1] - R[m_iter - 1][j - 1]) / (std::pow(4, j) - 1);

            BASE::m_estimate = R[m_iter][m_iter];
            m_iter++;
        }
    };

    template<typename FN, typename BOUNDS_T>
        requires IsFloatInvocable< FN > && IsFloatStruct< BOUNDS_T >
    Romberg(FN, BOUNDS_T) -> Romberg< FN, StructCommonType_t< BOUNDS_T > >;


    template<typename FN, typename ARG_T>
        requires IsFloatInvocable< FN > && IsFloat< ARG_T >
    Romberg(FN, std::initializer_list< ARG_T >) -> Romberg< FN, ARG_T >;


    // =================================================================================================================
    //  ad88888ba   88
    // d8"     "8b  ""
    // Y8,
    // `Y8aaaaa,    88  88,dPYba,,adPYba,   8b,dPPYba,   ,adPPYba,   ,adPPYba,   8b,dPPYba,
    //   `"""""8b,  88  88P'   "88"    "8a  88P'    "8a  I8[    ""  a8"     "8a  88P'   `"8a
    //         `8b  88  88      88      88  88       d8   `"Y8ba,   8b       d8  88       88
    // Y8a     a8P  88  88      88      88  88b,   ,a8"  aa    ]8I  "8a,   ,a8"  88       88
    //  "Y88888P"   88  88      88      88  88`YbbdP"'   `"YbbdP"'   `"YbbdP"'   88       88
    //                                      88
    //                                      88
    // =================================================================================================================

    template<IsFloatInvocable FN, IsFloat ARG_T = double>
    class Simpson final : public detail::IntegrationBase< Simpson< FN, ARG_T >, FN, ARG_T >
    {
        using BASE = detail::IntegrationBase< Simpson< FN, ARG_T >, FN, ARG_T >; // Base class alias for readability.

    public:
        using BASE::BASE; // Inherits constructors from IntegrationBase.
        using RESULT_T = std::invoke_result_t< FN, ARG_T >;

        static const inline std::string SolverName = "Simpson";

        int m_iter{ 1 };

        void operator()()
        {
            const auto& [lower, upper] = BASE::m_bounds;

            BASE::m_interval /= 2;                     // Halve the step size in each iteration
            int numOfPoints = 1 + std::pow(2, m_iter); // Calculate the number of points for this iteration

            // Create a vector of points
            std::vector< RESULT_T > points(numOfPoints);
            std::generate(points.begin(),
                          points.end(),
                          [n = 0, lower, h = BASE::m_interval]() mutable { return lower + n++ * h; });

            // Calculate the sum of function values at the points, with alternating coefficients 4 and 2
            RESULT_T sum = 0.0;
            for (size_t i = 1; i < numOfPoints - 1; ++i) { sum += BASE::evaluate(points[i]) * (i % 2 == 0 ? 2 : 4); }

            // Update the integral estimate
            BASE::m_estimate = BASE::m_interval / 3 * (BASE::evaluate(lower) + BASE::evaluate(upper) + sum);
            m_iter++;
        }
    };

    // Deduction guide similar to the one for Trapezoid
    template<typename FN, typename BOUNDS_T>
        requires IsFloatInvocable< FN > && IsFloatStruct< BOUNDS_T >
    Simpson(FN, BOUNDS_T) -> Simpson< FN, StructCommonType_t< BOUNDS_T > >;

    template<typename FN, typename ARG_T>
        requires IsFloatInvocable< FN > && IsFloat< ARG_T >
    Simpson(FN, std::initializer_list< ARG_T >) -> Simpson< FN, ARG_T >;


    // =================================================================================================================
    // 88
    // ""                ,d                                                       ,d
    //                   88                                                       88
    // 88  8b,dPPYba,  MM88MMM  ,adPPYba,   ,adPPYb,d8  8b,dPPYba,  ,adPPYYba,  MM88MMM  ,adPPYba,
    // 88  88P'   `"8a   88    a8P_____88  a8"    `Y88  88P'   "Y8  ""     `Y8    88    a8P_____88
    // 88  88       88   88    8PP"""""""  8b       88  88          ,adPPPPP88    88    8PP"""""""
    // 88  88       88   88,   "8b,   ,aa  "8a,   ,d88  88          88,    ,88    88,   "8b,   ,aa
    // 88  88       88   "Y888  `"Ybbd8"'   `"YbbdP"Y8  88          `"8bbdP"Y8    "Y888  `"Ybbd8"'
    //                                      aa,    ,88
    //                                       "Y8bbdP"
    // =================================================================================================================

    template<template< typename, typename > class SOLVER_T,
        IsFloatInvocable FN, IsFloatStruct STRUCT_T, IsFloat TOL_T = StructCommonType_t< STRUCT_T >, std::integral ITER_T = int>
    auto integrate(FN       function,
                   STRUCT_T bounds,
                   TOL_T    tolerance     = epsilon< StructCommonType_t< STRUCT_T > >(),
                   ITER_T   maxIterations = 25)
    {
        using RESULT_T = StructCommonType_t< STRUCT_T >;
        using ERROR_T = Error< detail::IntegrationErrorData< RESULT_T, ITER_T > >; /**< Type for error handling. */
        using RETURN_T = tl::expected< RESULT_T, ERROR_T >;                        /**< Type for the function return value. */
        using std::isfinite;

        const auto& [lower, upper] = bounds;
        auto        solver         = SOLVER_T< FN, RESULT_T >(function, bounds);

        auto result = solver.current();
        if (!isfinite(result)) {
            return RETURN_T(tl::make_unexpected(ERROR_T(
                decltype(solver)::SolverName + " integration failed: Initial estimate is not finite.",
                NumerixxErrorType::Integral,
                { .value = result, .eabs = 0.0, .erel = 0.0, .iterations = 0 })));
        }

        RESULT_T eabs;
        RESULT_T erel;

        for (auto i = 0; i < maxIterations; i++) {
            solver.iterate();

            if (!isfinite(solver.current())) {
                return RETURN_T(tl::make_unexpected(ERROR_T(decltype(solver)::SolverName + " integration failed: Result is not finite.",
                                                            NumerixxErrorType::Integral,
                                                            { .value = result, .eabs = 0.0, .erel = 0.0, .iterations = 0 })));
            }

            eabs = abs(solver.current() - result);
            erel = 1.0 - abs(solver.current() / result);

            if (eabs < tolerance)
                return RETURN_T(solver.current());
            result = solver.current();
        }

        return RETURN_T(tl::make_unexpected(ERROR_T(
            decltype(solver)::SolverName + " integration failed: Maximum number of iterations reached.",
            NumerixxErrorType::Integral,
            { .value = result, .eabs = eabs, .erel = erel, .iterations = maxIterations })));
    }

    template<template< typename, typename > class SOLVER_T,
             IsFloatInvocable FN, IsFloat ARG_T, IsFloat TOL_T = ARG_T, std::integral ITER_T = int,
             size_t N>
        requires (N == 2)
    auto integrate(FN            function,
                   const ARG_T (&bounds)[N],
                   TOL_T         tolerance     = epsilon< ARG_T >(),
                   ITER_T        maxIterations = 25)
    {
        auto bnds = std::span(bounds, N);
        return integrate< SOLVER_T >(function, std::pair(bnds.front(), bnds.back()), tolerance, maxIterations);
    }


    // =================================================================================================================
    //
    // 88                                                                       88    ,ad8888ba,       ad88
    // ""                ,d                                                     88   d8"'    `"8b     d8"
    //                   88                                                     88  d8'        `8b    88
    // 88  8b,dPPYba,  MM88MMM  ,adPPYba,   ,adPPYb,d8  8b,dPPYba,  ,adPPYYba,  88  88          88  MM88MMM
    // 88  88P'   `"8a   88    a8P_____88  a8"    `Y88  88P'   "Y8  ""     `Y8  88  88          88    88
    // 88  88       88   88    8PP"""""""  8b       88  88          ,adPPPPP88  88  Y8,        ,8P    88
    // 88  88       88   88,   "8b,   ,aa  "8a,   ,d88  88          88,    ,88  88   Y8a.    .a8P     88
    // 88  88       88   "Y888  `"Ybbd8"'   `"YbbdP"Y8  88          `"8bbdP"Y8  88    `"Y8888Y"'      88
    //                                      aa,    ,88
    //                                       "Y8bbdP"
    // =================================================================================================================

    namespace detail
    {
        template<template< typename, typename > class ALGO, IsFloatInvocable FN>
        class IntegrationFunctor
        {
            FN m_function{};

        public:
            template<IsFloatStruct STRUCT_T, IsFloat TOL_T = StructCommonType_t< STRUCT_T >>
            auto operator()(STRUCT_T bounds,
                            TOL_T    tol  = epsilon< StructCommonType_t< STRUCT_T > >(),
                            int      iter = 25) const
            {
                auto result = integrate< ALGO >(m_function, bounds, tol, iter);
                if (result) return result.value();
                throw result.error();
            }

            template<IsFloat ARG_T, IsFloat TOL_T = ARG_T, size_t N>
            auto operator()(const ARG_T (&bounds)[N],
                            TOL_T         tol  = epsilon< ARG_T >(),
                            int           iter = 25) const
                requires (N == 2)
            {
                auto bnds   = std::span(bounds, N);
                auto result = integrate< ALGO >(m_function, std::pair(bnds.front(), bnds.back()), tol, iter);
                if (result) return result.value();
                throw result.error();
            }
        };
    } // namespace detail

    template<template< typename, typename > class ALGO_T = Romberg, IsFloatInvocable FN>
    auto integralOf(FN) { return detail::IntegrationFunctor< ALGO_T, FN >(); }
} // namespace nxx::integrate

#endif    // NUMERIXX_INTEGRATION_HPP

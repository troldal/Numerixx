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

/**
 * @file Integration.hpp
 * @brief Comprehensive header file for numerical integration utilities.
 *
 * This header file, Integration.hpp, provides a full suite of tools and templates for performing
 * numerical integration in C++. It includes a range of classes and functions that facilitate
 * numerical integration using various algorithms like Trapezoid, Simpson, Romberg, and others.
 * The file is structured to provide a cohesive and extensible framework for integration tasks,
 * suitable for both educational and professional applications in scientific computing,
 * engineering, and related fields.
 *
 * The file contains:
 * - Base class templates for integration solvers (IntegrationBase).
 * - Specific integration algorithm implementations (Trapezoid, Simpson, Romberg).
 * - Generic integration function templates (integrate) that work with various solver classes.
 * - Overloads of the integrate function to support different types of bounds (e.g., arrays, structures).
 * - A functor class template (IntegrationFunctor) for encapsulating integration algorithms as callable objects.
 * - A factory function template (integralOf) for easy creation of integration functors.
 *
 * @note
 * This file is designed to be included as a single header, providing all necessary tools for
 * numerical integration without the need for multiple includes. It relies on several external
 * dependencies, including the standard library, Boost (for multi_array), and tl::expected for
 * error handling.
 *
 * @warning
 * The user must ensure that the function to be integrated and the integration bounds are
 * compatible with the algorithms used. Improper use or configuration may result in inaccurate
 * results or runtime errors.
 *
 * @todo Consider approaches to improve memory usage, especially for the Trapezoid and Simpson classes.
 * @todo Consider implementing a parallel version of the integration algorithms.
 * @todo Implement performance tests and benchmarks.
 */

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

        /**
         * @class IntegrationBase
         * @brief Template base class for numerical integration solvers.
         *
         * @details This class template serves as a base for creating various numerical integration solvers.
         *          It enforces constraints on the types of function and arguments used in the integration process.
         *          It provides the fundamental interface and members that all derived integration solver classes should adhere to.
         *
         * @tparam DERIVED The derived class following the Curiously Recurring Template Pattern (CRTP).
         * @tparam FUNCTION_T The type of function to be integrated.
         * @tparam ARG_T The type of the argument for the function.
         *
         * @note This class uses the CRTP design pattern to enable static polymorphism.
         */
        template<typename DERIVED, typename FUNCTION_T, typename ARG_T>
            requires std::same_as< typename IntegrationTraits< DERIVED >::FUNCTION_T, FUNCTION_T > &&
                     nxx::IsFloatInvocable< FUNCTION_T > &&
                     nxx::IsFloat< ARG_T > &&
                     nxx::IsFloat< typename IntegrationTraits< DERIVED >::RETURN_T >
        class IntegrationBase
        {
            friend DERIVED; // Granting access to derived classes.

        public:
            static constexpr bool IsIntegrationSolver = true; /**< Flag indicating the class is a bracketing solver. */

            using RESULT_T = std::invoke_result_t< FUNCTION_T, ARG_T >; /**< Result type of the function. */
            using BOUNDS_T = std::pair< ARG_T, ARG_T >;                 /**< Type for representing the bounds around the root. */

        protected:
            ~IntegrationBase() = default; /**< Protected destructor to prevent direct instantiation. */

        private:
            FUNCTION_T m_func{};     /**< The function object to integrate. */
            BOUNDS_T   m_bounds{};   /**< Holds the current bounds of the integration. */
            RESULT_T   m_estimate{}; /**< Holds the current estimate of the integral. */
            ARG_T      m_interval{}; /**< Holds the current interval size. */

        public:
            /**
             * @brief Constructor initializing the integration solver with a function and bounds.
             * @param objective The function to integrate.
             * @param bounds The bounds of integration.
             */
            IntegrationBase(FUNCTION_T objective, const IsFloatStruct auto& bounds)
                : m_func{ std::move(objective) } { init(bounds); }

            /**
             * @brief Constructor initializing the integration solver with a function and bounds array.
             * @param objective The function to integrate.
             * @param bounds Array representing the bounds of integration.
             * @tparam N Size of the bounds array, must be 2.
             */
            template<size_t N> requires (N == 2)
            IntegrationBase(FUNCTION_T objective, const ARG_T (&bounds)[N])
                : m_func{ std::move(objective) } { init(std::pair{ bounds[0], bounds[1] }); }

            // Rule of Five for proper management of resources.
            IntegrationBase(const IntegrationBase& other)                = default; /**< Default copy constructor. */
            IntegrationBase(IntegrationBase&& other) noexcept            = default; /**< Default move constructor. */
            IntegrationBase& operator=(const IntegrationBase& other)     = default; /**< Default copy assignment operator. */
            IntegrationBase& operator=(IntegrationBase&& other) noexcept = default; /**< Default move assignment operator. */

            /**
             * @brief Initializes the integration solver with bounds.
             * @param bounds The bounds of integration.
             * @tparam STRUCT_T The type of the bounds, must satisfy the IsFloatStruct concept.
             */
            template<IsFloatStruct STRUCT_T>
            void init(STRUCT_T bounds)
            {
                auto [lower, upper] = bounds;
                validateRange(lower, upper);
                m_bounds   = BOUNDS_T{ lower, upper };
                m_interval = upper - lower;
                m_estimate = m_interval * (evaluate(lower) + evaluate(upper)) / 2;
            }

            /**
             * @brief Evaluates the function at a given value.
             * @param value The argument value to evaluate the function at.
             * @return The result of the function evaluation.
             */
            RESULT_T evaluate(ARG_T value) const { return m_func(value); }

            /**
             * @brief Retrieves the current estimate of the integral.
             * @return The current integral estimate.
             */
            RESULT_T current() const { return m_estimate; }

            /**
             * @brief Performs an iteration of the integration algorithm.
             *
             * @details This function triggers the derived class's specific integration algorithm.
             *          This enforces that the derived classes has an overloaded function call operator.
             *          The derived class uses the current state of this base class (like bounds, function, etc.)
             *          to perform an iteration of the integration process.
             */
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

    /**
     * @class Trapezoid
     * @brief Implementation of numerical integration using the trapezoidal rule.
     *
     * @details This class implements the numerical integration using the trapezoidal rule.
     *          It is derived from IntegrationBase and overrides the operator() to implement
     *          the specific logic for the trapezoidal rule.
     *
     *          The Trapezoid class is designed for use in single-threaded environments.
     *          It is not thread-safe and should not be accessed concurrently from multiple
     *          threads. If multithreaded usage is required, it is the responsibility of the
     *          user to ensure proper synchronization mechanisms are in place to prevent
     *          concurrent access and modification.
     *
     *          The class focuses on performance and efficiency in single-threaded scenarios.
     *          Using this class in a multithreaded context without adequate synchronization
     *          can lead to undefined behavior and incorrect computation results.
     *
     * @tparam FN The type of function to be integrated.
     * @tparam ARG_T The type of the argument for the function, defaulted to double.
     */
    template<IsFloatInvocable FN, IsFloat ARG_T = double>
    class Trapezoid final : public detail::IntegrationBase< Trapezoid< FN, ARG_T >, FN, ARG_T >
    {
        using BASE = detail::IntegrationBase< Trapezoid< FN, ARG_T >, FN, ARG_T >; /**< Base class alias for readability. */

    public:
        using BASE::BASE; /**< Inherits constructors from IntegrationBase. */
        using RESULT_T = std::invoke_result_t< FN, ARG_T >;

        static const inline std::string SolverName = "Trapezoid"; /**< Name of the solver. */

        int                     m_iter{ 1 };          /**< Iteration counter. */
        uint64_t                m_midpointCount{ 1 }; // Initialize numOfMidpoints
        std::vector< RESULT_T > m_midpoints;          // Member variable for midpoints

        /**
         * @brief Overloaded function call operator that performs a single iteration of the trapezoidal rule.
         *
         * @details This method calculates the next approximation of the integral using the trapezoidal rule.
         *          It divides the integration interval into smaller sub-intervals and calculates the area
         *          of trapezoids formed under the curve of the function.
         */
        void operator()()
        {
            const auto& [lower, upper] = BASE::m_bounds;

            // Recalculate the interval based on the iteration number
            uint64_t divisor = 1 << m_iter;
            BASE::m_interval = (upper - lower) / divisor;

            // Double the number of midpoints for this iteration
            if (m_iter > 1) m_midpointCount *= 2;

            // Resize midPoints vector if necessary
            m_midpoints.resize(m_midpointCount);
            std::generate(m_midpoints.begin(),
                          m_midpoints.end(),
                          [n = 1, lower, h = BASE::m_interval]() mutable { return lower + (2 * n++ - 1) * h; });

            // Calculate the sum of function values at the midpoints
            const RESULT_T sum =
                std::transform_reduce(m_midpoints.begin(),
                                      m_midpoints.end(),
                                      ARG_T(0.0),
                                      std::plus(),
                                      [&](RESULT_T x) { return BASE::evaluate(x); });

            // Update the integral estimate
            BASE::m_estimate = BASE::m_estimate / 2 + BASE::m_interval * sum;
            m_iter++;
        }
    };

    /**
     * @brief Deduction guides for Trapezoid class.
     * Allows the type of Trapezoid class to be deduced from the constructor parameters.
     */
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

    /**
     * @class Romberg
     * @brief Implementation of numerical integration using the Romberg integration method.
     *
     * @details This class implements the numerical integration using the Romberg integration method.
     *          It is derived from IntegrationBase and overrides the operator() to implement
     *          the specific logic for Romberg integration.
     *
     *          The Romberg class is designed for use in single-threaded environments.
     *          It is not thread-safe and should not be accessed concurrently from multiple
     *          threads. If multithreaded usage is required, it is the responsibility of the
     *          user to ensure proper synchronization mechanisms are in place to prevent
     *          concurrent access and modification.
     *
     *          The class focuses on performance and efficiency in single-threaded scenarios.
     *          Using this class in a multithreaded context without adequate synchronization
     *          can lead to undefined behavior and incorrect computation results.
     *
     * @tparam FN The type of function to be integrated.
     * @tparam ARG_T The type of the argument for the function, defaulted to double.
     */
    template<IsFloatInvocable FN, IsFloat ARG_T = double>
    class Romberg final : public detail::IntegrationBase< Romberg< FN, ARG_T >, FN, ARG_T >
    {
        using BASE = detail::IntegrationBase< Romberg< FN, ARG_T >, FN, ARG_T >; /**< Base class alias for readability. */

    public:
        using BASE::BASE; /**< Inherits constructors from BracketingBase. */
        using RESULT_T = std::invoke_result_t< FN, ARG_T >;
        using ARRAY_T = boost::multi_array< RESULT_T, 2 >; /**< Type for the Romberg integration table. */

        static const inline std::string SolverName = "Romberg"; /**< Name of the solver. */

        int     m_iter{ 1 }; /**< Iteration counter. */
        ARRAY_T R;           /**< Romberg integration table. */

        /**
         * @brief Overloaded function call operator that performs a single iteration of Romberg integration.
         *
         * @details This method calculates the next approximation of the integral using the Romberg integration method.
         *          It updates the Romberg table with new values computed in each iteration and uses these values to
         *          extrapolate to higher order estimates of the integral.
         */
        void operator()()
        {
            using boost::extents;
            R.resize(extents[m_iter + 1][m_iter + 1]);

            const auto& [lower, upper] = BASE::m_bounds;

            // Initialize the first element with the basic trapezoidal rule
            if (m_iter == 1) { R[0][0] = BASE::m_interval * (BASE::evaluate(lower) + BASE::evaluate(upper)) / 2; }

            // Halve the step size for this iteration
            uint64_t divisor = 1 << m_iter;
            BASE::m_interval = (upper - lower) / divisor;

            // Trapezoidal rule: Calculate the sum of the function values at the midpoints
            RESULT_T       sum          = 0.0;
            const uint64_t numMidpoints = 1 << (m_iter - 1); // 2^(m_iter - 1) using bitwise shift
            for (uint64_t k = 1; k <= numMidpoints; ++k)
                sum += BASE::evaluate(lower + (2 * k - 1) * BASE::m_interval);

            // Update the Romberg table for the first column (trapezoidal rule)
            R[m_iter][0] = R[m_iter - 1][0] / 2 + BASE::m_interval * sum;

            // Apply the Romberg integration formula
            RESULT_T four_pow_j = 1; // Initial value for 4^j
            for (int j = 1; j <= m_iter; ++j) {
                four_pow_j *= 4; // Multiply by 4 at each step
                R[m_iter][j] = R[m_iter][j - 1] + (R[m_iter][j - 1] - R[m_iter - 1][j - 1]) / (four_pow_j - 1);
            }

            BASE::m_estimate = R[m_iter][m_iter];
            m_iter++;
        }
    };

    /**
     * @brief Deduction guides for Romberg class.
     * Allows the type of Romberg class to be deduced from the constructor parameters.
     */
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

    /**
     * @class Simpson
     * @brief Implementation of numerical integration using Simpson's rule.
     *
     * @details This class implements the numerical integration using Simpson's rule.
     *          It is derived from IntegrationBase and overrides the operator() to implement
     *          the specific logic for Simpson's rule integration.
     *
     *          The Simpson class is designed for use in single-threaded environments.
     *          It is not thread-safe and should not be accessed concurrently from multiple
     *          threads. If multithreaded usage is required, it is the responsibility of the
     *          user to ensure proper synchronization mechanisms are in place to prevent
     *          concurrent access and modification.
     *
     *          The class focuses on performance and efficiency in single-threaded scenarios.
     *          Using this class in a multithreaded context without adequate synchronization
     *          can lead to undefined behavior and incorrect computation results.
     *
     * @tparam FN The type of function to be integrated.
     * @tparam ARG_T The type of the argument for the function, defaulted to double.
     */
    template<IsFloatInvocable FN, IsFloat ARG_T = double>
    class Simpson final : public detail::IntegrationBase< Simpson< FN, ARG_T >, FN, ARG_T >
    {
        using BASE = detail::IntegrationBase< Simpson< FN, ARG_T >, FN, ARG_T >; // Base class alias for readability.

    public:
        using BASE::BASE; // Inherits constructors from IntegrationBase.
        using RESULT_T = std::invoke_result_t< FN, ARG_T >;

        static const inline std::string SolverName = "Simpson"; /**< Name of the solver. */

        int                     m_iter{ 1 };
        int                     numOfIntervals{ 1 }; // Initialize numOfIntervals
        std::vector< RESULT_T > points;              // Member variable for points

        /**
         * @brief Overloaded function call operator that performs a single iteration of Simpson's rule.
         *
         * @details This method calculates the next approximation of the integral using Simpson's rule.
         *          It divides the integration interval into smaller sub-intervals and calculates the area
         *          under the curve of the function using Simpson's rule.
         */
        void operator()()
        {
            const auto& [lower, upper] = BASE::m_bounds;

            // Halve the step size for this iteration
            uint64_t divisor = 1 << m_iter;
            BASE::m_interval = (upper - lower) / divisor;

            // Double the number of intervals using arithmetic operation
            numOfIntervals *= 2;

            // Calculate the total number of points
            int numOfPoints = numOfIntervals + 1;

            // Resize points vector if necessary
            points.resize(numOfPoints);
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

    /**
     * @brief Deduction guides for Romberg class.
     * Allows the type of Romberg class to be deduced from the constructor parameters.
     */
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

    /**
     * @brief Performs numerical integration using a specified solver.
     *
     * @details This function template provides a generic interface for numerical integration. It constructs
     *          an instance of the specified solver class and iteratively refines the estimate of the integral
     *          within given bounds and tolerance limits.
     *
     * @tparam SOLVER_T Template template parameter specifying the solver class.
     * @tparam FN The type of function to be integrated.
     * @tparam STRUCT_T The type of the structure representing bounds.
     * @tparam TOL_T The type of the tolerance value, defaulted to the common type of STRUCT_T.
     * @tparam ITER_T The type of the iteration count, defaulted to int.
     * @param function The function to be integrated.
     * @param bounds The bounds of integration.
     * @param tolerance The tolerance for the result accuracy, defaults to machine epsilon for the result type.
     * @param maxIterations The maximum number of iterations allowed, defaults to 25.
     * @return A tl::expected object containing either the result of the integration or an error.
     */
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

    /**
     * @brief Overload of the integrate function accepting bounds as an array.
     *
     * @details This function template is an overload of the primary `integrate` function. It accepts
     *          the bounds of integration as an array and delegates to the primary `integrate` function
     *          by converting these bounds into a pair.
     *
     * @tparam SOLVER_T Template template parameter specifying the solver class.
     * @tparam FN The type of function to be integrated.
     * @tparam ARG_T The type of the argument for the bounds.
     * @tparam TOL_T The type of the tolerance value, defaulted to ARG_T.
     * @tparam ITER_T The type of the iteration count, defaulted to int.
     * @tparam N The size of the bounds array, must be 2.
     * @param function The function to be integrated.
     * @param bounds Array representing the bounds of integration.
     * @param tolerance The tolerance for the result accuracy, defaults to machine epsilon for the argument type.
     * @param maxIterations The maximum number of iterations allowed, defaults to 25.
     * @return A tl::expected object containing either the result of the integration or an error.
     */
    template<template< typename, typename > class SOLVER_T,
             IsFloatInvocable FN, IsFloat ARG_T, IsFloat TOL_T = ARG_T, std::integral ITER_T = int,
             size_t N>
        requires (N == 2)
    auto integrate(FN            function,
                   const ARG_T (&bounds)[N],
                   TOL_T         tolerance     = epsilon< ARG_T >(),
                   ITER_T        maxIterations = 25)
    {
        return integrate< SOLVER_T >(function, std::pair(bounds[0], bounds[1]), tolerance, maxIterations);
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
        /**
         * @class IntegrationFunctor
         * @brief Functor class for numerical integration algorithms.
         *
         * @details This class encapsulates a function to be integrated and provides a convenient functor
         *          interface for performing numerical integration. It utilizes the integrate function and
         *          handles the result and error scenarios.
         *
         * @tparam ALGO Template template parameter specifying the integration algorithm class.
         * @tparam FN The type of function to be integrated.
         */
        template<template< typename, typename > class ALGO, IsFloatInvocable FN>
        class IntegrationFunctor
        {
            FN m_function{}; /**< The function to be integrated. */

        public:
            /**
             * @brief Functor operator for integration with structure-based bounds.
             *
             * @tparam STRUCT_T The type of the structure representing bounds.
             * @tparam TOL_T The type of the tolerance value, defaulted to the common type of STRUCT_T.
             * @param bounds The bounds of integration.
             * @param tol The tolerance for the result accuracy, defaults to machine epsilon for the result type.
             * @param iter The maximum number of iterations allowed, defaults to 25.
             * @return The result of the integration.
             * @throws Throws an error if the integration result is an unexpected value.
             */
            template<IsFloatStruct STRUCT_T, IsFloat TOL_T = StructCommonType_t< STRUCT_T >>
            auto operator()(STRUCT_T bounds,
                            TOL_T    tol  = epsilon< StructCommonType_t< STRUCT_T > >(),
                            int      iter = 25) const
            {
                auto result = integrate< ALGO >(m_function, bounds, tol, iter);
                if (result) return result.value();
                throw result.error();
            }

            /**
             * @brief Functor operator for integration with array-based bounds.
             *
             * @tparam ARG_T The type of the argument for the bounds.
             * @tparam TOL_T The type of the tolerance value, defaulted to ARG_T.
             * @tparam N The size of the bounds array, must be 2.
             * @param bounds Array representing the bounds of integration.
             * @param tol The tolerance for the result accuracy, defaults to machine epsilon for the argument type.
             * @param iter The maximum number of iterations allowed, defaults to 25.
             * @return The result of the integration.
             * @throws Throws an error if the integration result is an unexpected value.
             */
            template<IsFloat ARG_T, IsFloat TOL_T = ARG_T, size_t N>
            auto operator()(const ARG_T (&bounds)[N],
                            TOL_T         tol  = epsilon< ARG_T >(),
                            int           iter = 25) const
                requires (N == 2)
            {
                auto result = integrate< ALGO >(m_function, std::pair(bounds[0], bounds[1]), tol, iter);
                if (result) return result.value();
                throw result.error();
            }
        };
    } // namespace detail

    /**
     * @brief Factory function to create an IntegrationFunctor for a given function.
     *
     * @details This function template provides a convenient way to create an IntegrationFunctor for a
     *          specific function and a numerical integration algorithm. It abstracts the creation process
     *          and allows users to focus on specifying the function and optionally the integration algorithm.
     *
     * @tparam ALGO_T Template template parameter specifying the integration algorithm class, defaulted to Romberg.
     * @tparam FN The type of function to be integrated.
     * @param function The function for which the integral functor is to be created.
     * @return An instance of IntegrationFunctor configured with the given function and algorithm.
     */
    template<template< typename, typename > class ALGO_T = Romberg, IsFloatInvocable FN>
    auto integralOf(FN function) { return detail::IntegrationFunctor< ALGO_T, FN >(); }
} // namespace nxx::integrate

#endif    // NUMERIXX_INTEGRATION_HPP

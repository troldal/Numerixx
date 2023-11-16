//
// Created by kenne on 12/11/2023.
//

#ifndef NUMERIXX_INTEGRATION_HPP
#define NUMERIXX_INTEGRATION_HPP

// ===== Numerixx Includes
#include <Concepts.hpp>
#include <Error.hpp>

// ===== External Includes
#include <gcem.hpp>
#include <tl/expected.hpp>
#include <boost/multi_array.hpp>

// ===== Standard Library Includes
#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>

namespace nxx::integrate
{
    namespace detail
    {
        /**
         * @brief Alias template for determining the step size for a given type when computing numerical derivatives.
         * @tparam T The type for which to compute the step size.
         */
        template<typename T>
            requires std::floating_point< T >
        struct StepSizeHelper
        {
        private:
            static constexpr T exponent = T(1.0 / 3.0);

        public:
            static constexpr T value = gcem::pow(std::numeric_limits< T >::epsilon(), exponent);
        };

        /**
         * @brief Constant expression for the step size used when computing numerical derivatives.
         * @tparam T The type for which to compute the step size.
         */
        template<typename T>
            requires std::floating_point< T >
        inline constexpr T StepSize = StepSizeHelper< T >::value;

        template<typename ALGO>
        class IntSolverTemplate
        {
        public:
            static constexpr auto IsIntSolver = true;

            /**
             * @brief Function call operator.
             *
             * @param function The function for which to compute the derivative. The function can be any callable type taking a
             * floating point type as an argument, and returns a value of the same type.
             * @param val The value at which to compute the derivative.
             * @param stepsize The finite difference used for computing the derivative.
             *
             * @return The derivative. The return type is the same as the return type of the provided function.
             *
             * @note This function is not intended to be used directly. Instead, use the \c diff template function,
             * or one of the convenience functions.
             *
             * @throws NumerixxError if stepsize is invalid.
             */
            template<std::floating_point TOL_T = double, std::integral ITER_T = int>
            inline auto operator()(IsFloatInvocable auto    function,
                                   std::floating_point auto lower,
                                   std::floating_point auto upper,
                                   TOL_T                    tolerance,
                                   ITER_T                   maxIterations) const
            {
                using ARG_T = std::common_type_t< decltype(lower), decltype(upper) >;
                using RETURN_T = std::invoke_result_t< decltype(function), ARG_T >;
                static_assert(std::floating_point< RETURN_T >, "The return type of the provided function must be a floating point type.");
                return ALGO{}(function, lower, upper, tolerance, maxIterations);
            }
        };

        /**
         * @brief Class implementing the Trapezoid Integration Method.
         *
         * The Trapezoid Method is a numerical technique for approximating the definite integral of a function.
         */
        class TrapezoidFunctor
        {
        public:
            /**
             * @brief Function call operator that performs integration using the Trapezoid method.
             *
             * This function calculates the integral of a provided function over a specified interval [a, b],
             * using the Trapezoid integration method with adaptive step sizing.
             *
             * @tparam Function Type of the function to integrate.
             * @param function The function to integrate.
             * @param a The lower limit of integration.
             * @param b The upper limit of integration.
             * @param tolerance The tolerance for convergence of the method.
             * @param maxIterations The maximum number of iterations allowed for convergence.
             * @return The computed integral of the function over [a, b].
             * @throws std::runtime_error If the method does not converge within the specified number of iterations.
             */
            template<std::floating_point TOL_T = double, std::integral ITER_T = int>
            auto operator()(IsFloatInvocable auto    function,
                            std::floating_point auto lower,
                            std::floating_point auto upper,
                            TOL_T                    tolerance     = 1e-6,
                            ITER_T                   maxIterations = 100) const
            {
                // Initial step size and initial trapezoidal estimate

                using ARG_T = std::common_type_t< decltype(lower), decltype(upper) >;
                using RETURN_T = std::invoke_result_t< decltype(function), ARG_T >;

                ARG_T    h    = (upper - lower);
                RETURN_T Iold = h * (function(lower) + function(upper)) / 2;

                // Main loop for adaptive trapezoidal rule
                for (ITER_T i = 1; i <= maxIterations; ++i) {
                    h /= 2;                                     // Halve the step size in each iteration
                    ITER_T numOfMidpoints = std::pow(2, i - 1); // Calculate the number of midpoints for this iteration

                    // Create a vector of midpoints
                    std::vector< RETURN_T > midPoints(numOfMidpoints);
                    std::generate(midPoints.begin(),
                                  midPoints.end(),
                                  [n = ITER_T{ 1 }, lower, h]() mutable { return lower + (2 * n++ - 1) * h; });

                    // Calculate the sum of function values at the midpoints
                    RETURN_T sum =
                        std::transform_reduce(midPoints.begin(),
                                              midPoints.end(),
                                              0.0,
                                              std::plus< >(),
                                              [&](ARG_T x) { return function(x); });

                    // Update the integral estimate
                    RETURN_T Inew = Iold / 2 + h * sum;

                    // Check for convergence
                    if (i > 1 && std::abs(Inew - Iold) < tolerance) return Inew; // Return the new estimate if converged

                    Iold = Inew; // Update the old estimate for the next iteration
                }

                // Throw an exception if the method does not converge within the maximum number of iterations
                throw std::runtime_error("AdaptiveTrapezoidMethod did not converge within the specified number of iterations.");
            }
        };

        /**
         * @brief Class implementing the Adaptive Simpson's Integration Method.
         *
         * Adaptive Simpson's Method is a numerical technique to approximate the definite integral of a function,
         * which adapts the size of the intervals based on the function's behavior.
         */
        class SimpsonFunctor
        {
        public:
            /**
             * @brief Function call operator that performs integration using the Adaptive Simpson's method.
             *
             * This function calculates the integral of a provided function over a specified interval [a, b],
             * using an adaptive approach to improve the accuracy of Simpson's method.
             *
             * @tparam Function Type of the function to integrate.
             * @param function The function to integrate.
             * @param a The lower limit of integration.
             * @param b The upper limit of integration.
             * @param tolerance The tolerance for convergence of the method.
             * @param maxDepth The maximum recursion depth to prevent infinite recursion.
             * @return The computed integral of the function over [a, b].
             */
            template<IsFloatInvocable FN, std::floating_point TOL_T = double, std::integral ITER_T = int>
            auto operator()(FN                       function,
                            std::floating_point auto lower,
                            std::floating_point auto upper,
                            TOL_T                    tolerance = 1e-6,
                            ITER_T                   maxDepth  = 100) const
            {
                // Begin the recursive adaptive Simpson's method
                return adaptiveSimpson(function, lower, upper, tolerance, maxDepth, 0);
            }

        private:
            /**
             * @brief Private helper function for recursive adaptive Simpson's integration.
             *
             * This function recursively applies Simpson's method to smaller and smaller intervals until the desired tolerance is achieved.
             *
             * @tparam Function Type of the function to integrate.
             * @param function The function to integrate.
             * @param a The lower limit of the current interval.
             * @param b The upper limit of the current interval.
             * @param tolerance The tolerance for convergence in the current recursive call.
             * @param maxDepth The maximum allowed recursion depth.
             * @param depth The current depth of recursion.
             * @return The approximate integral over the interval [a, b].
             */
            template<IsFloatInvocable FN, std::floating_point TOL_T = double, std::integral ITER_T = int>
            auto adaptiveSimpson(FN                       function,
                                 std::floating_point auto lower,
                                 std::floating_point auto upper,
                                 TOL_T                    tolerance,
                                 ITER_T                   maxDepth,
                                 ITER_T                   depth) const
            {
                double m  = (lower + upper) / 2; // Midpoint of the interval
                double h  = (upper - lower) / 2; // Half the interval length
                double fa = function(lower), fb = function(upper), fm = function(m);

                // Standard Simpson's rule
                double S = (h / 3) * (fa + 4 * fm + fb);

                // Simpson's rule on the left and right subintervals
                double Sleft  = (h / 6) * (fa + 4 * function(lower + h / 2) + fm);
                double Sright = (h / 6) * (fm + 4 * function(m + h / 2) + fb);
                double S2     = Sleft + Sright;

                // Check if the current approximation is within the tolerance
                if (depth >= maxDepth || std::abs(S2 - S) < 15 * tolerance) return S2;

                // Recursively apply the method to each half
                return adaptiveSimpson(function, lower, m, tolerance / 2, maxDepth, depth + 1) +
                       adaptiveSimpson(function, m, upper, tolerance / 2, maxDepth, depth + 1);
            }
        };

        /**
         * @brief Class implementing the Romberg Integration Method.
         *
         * The Romberg Method is used for numerical integration, which is an application of Richardson extrapolation
         * to evaluate the definite integral of a function.
         */
        class RombergFunctor
        {
        public:
            /**
             * @brief Function call operator that performs Romberg integration.
             *
             * This template function calculates the integral of a provided function over a specified interval [a, b],
             * using the Romberg integration method.
             *
             * @param function The function to integrate.
             * @param lower The lower limit of integration.
             * @param upper The upper limit of integration.
             * @param tolerance The tolerance for convergence of the method.
             * @param maxIterations The maximum number of iterations allowed for convergence.
             * @return The computed integral of the function over [a, b].
             * @throws std::runtime_error If the method does not converge within the specified number of iterations.
             */
            template<std::floating_point TOL_T = double, std::integral ITER_T = int>
            auto operator()(IsFloatInvocable auto    function,
                            std::floating_point auto lower,
                            std::floating_point auto upper,
                            TOL_T                    tolerance     = 1e-6,
                            ITER_T                   maxIterations = 1000) const
            {
                using ARG_T = std::common_type_t< decltype(lower), decltype(upper) >;
                using RETURN_T = std::invoke_result_t< decltype(function), ARG_T >;
                using TYPE_T = std::common_type_t< RETURN_T, ARG_T >;

                using boost::extents;
                using boost::multi_array;

                multi_array< TYPE_T, 2 > R(extents[maxIterations][maxIterations]);

                // Initialize the first element with the basic trapezoidal rule
                R[0][0] = (upper - lower) * (function(lower) + function(upper)) / 2;

                for (ITER_T i = 1; i < maxIterations; ++i) {
                    // Calculate the step size for this iteration
                    TYPE_T h = (upper - lower) / std::pow(2, i);

                    // Trapezoidal rule: Calculate the sum of the function values at the midpoints
                    TYPE_T sum = 0.0;
                    for (ITER_T k = 1; k <= std::pow(2, i - 1); ++k)
                        sum += function(lower + (2 * k - 1) * h);

                    // Update the Romberg table for the first column (trapezoidal rule)
                    R[i][0] = R[i - 1][0] / 2 + h * sum;

                    // Apply the Romberg integration formula
                    for (ITER_T j = 1; j <= i; ++j)
                        // Use the previously computed values to extrapolate to higher orders
                        R[i][j] = R[i][j - 1] + (R[i][j - 1] - R[i - 1][j - 1]) / (std::pow(4, j) - 1);

                    // Check for convergence: if the difference between the last two diagonal elements is within tolerance
                    if (i > 1 && std::abs(R[i][i] - R[i - 1][i - 1]) < tolerance)
                        return R[i][i]; // Return the converged value
                }

                // If max iterations are reached without converging, throw an exception
                throw std::runtime_error("RombergMethod did not converge within the specified number of iterations.");
            }
        };
    } // namespace detail

    using Romberg = detail::IntSolverTemplate< detail::RombergFunctor >;
    using Simpson = detail::IntSolverTemplate< detail::SimpsonFunctor >;
    using Trapezoid = detail::IntSolverTemplate< detail::TrapezoidFunctor >;

    template< typename ALGO = Romberg, IsFloatInvocable FN, std::floating_point TOL_T = double, std::integral ITER_T = int >
    auto integrate(FN                       function,
                   std::floating_point auto lower,
                   std::floating_point auto upper,
                   TOL_T                    tolerance     = 1e-6,
                   ITER_T                   maxIterations = 100)
    {
        // if (!function) throw std::runtime_error("Function object is invalid.");
        return ALGO{}(function, lower, upper, tolerance, maxIterations);
    }

    template<typename ALGO = Romberg >
    inline auto integralOf(IsFloatInvocable auto function)
    {
        // if (!function) throw std::runtime_error("Function object is invalid.");
        return [=](std::floating_point auto a, std::floating_point auto b) { return ALGO {}(function, a, b); };
    }
}    // namespace nxx::integrate

#endif    // NUMERIXX_INTEGRATION_HPP
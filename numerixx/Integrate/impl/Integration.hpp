//
// Created by kenne on 12/11/2023.
//

#ifndef NUMERIXX_INTEGRATION_IMP_HPP
#define NUMERIXX_INTEGRATION_IMP_HPP

// ===== Numerixx Includes
#include <Concepts.hpp>
#include <Error.hpp>

// ===== External Includes
#include <gcem.hpp>
#include <tl/expected.hpp>

// ===== Standard Library Includes
#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>

namespace nxx::integrate
{

    /**
     * @brief Alias template for the return type of a function.
     */
    template< typename T >
    using ReturnType = std::invoke_result_t< T, double >;

    /**
     * @brief Alias template for determining the step size for a given type when computing numerical derivatives.
     * @tparam T The type for which to compute the step size.
     */
    template< typename T >
    requires std::floating_point< T >
    struct StepSizeHelper
    {
        static constexpr T value = gcem::pow(std::numeric_limits< T >::epsilon(), 1.0 / 3.0);
    };

    /**
     * @brief Constant expression for the step size used when computing numerical derivatives.
     * @tparam T The type for which to compute the step size.
     */
    template< typename T >
    requires std::floating_point< T >
    inline constexpr T StepSize = StepSizeHelper< T >::value;

    template< typename FN >
    concept IsFunction = requires(FN fn) {
                             {
                                 fn(0.0)
                             } -> std::floating_point;
                         };

    /**
     * @brief Class implementing the Romberg Integration Method.
     *
     * The Romberg Method is used for numerical integration, which is an application of Richardson extrapolation
     * to evaluate the definite integral of a function.
     */
    class Romberg
    {
    public:
        /**
         * @brief Function call operator that performs Romberg integration.
         *
         * This template function calculates the integral of a provided function over a specified interval [a, b],
         * using the Romberg integration method.
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
        template< IsFunction FN, std::floating_point TOL_T = double, std::integral ITER_T = int >
        auto operator()(FN                       function,
                        std::floating_point auto lower,
                        std::floating_point auto upper,
                        TOL_T                    tolerance     = 1e-6,
                        ITER_T                   maxIterations = 100) const
        {
            using ReturnType = std::invoke_result_t< FN, std::common_type_t< decltype(lower), decltype(upper) > >;
            std::vector< std::vector< ReturnType > > R(
                1,
                std::vector< ReturnType >(1, 0.5 * (upper - lower) * (function(lower) + function(upper))));

            for (int i = 1; i <= maxIterations; ++i) {
                // Enlarge the R vector for the next iteration.
                R.emplace_back(i + 1, 0);
                double h = (upper - lower) / std::pow(2.0, i);

                // Generate midpoints for the trapezoidal rule.
                std::vector< double > midPoints;
                midPoints.reserve(std::pow(2, i - 1));
                std::generate_n(std::back_inserter(midPoints), std::pow(2, i - 1), [n = 1, lower, h]() mutable {
                    return lower + (2 * n++ - 1) * h;
                });

                // Calculate the trapezoidal rule for the current iteration.
                R.at(i).at(0) =
                    0.5 * R.at(i - 1).at(0) +
                    h * std::accumulate(midPoints.begin(), midPoints.end(), 0.0, [&](double acc, double x) { return acc + function(x); });

                // Apply the Romberg extrapolation for the current iteration.
                for (int j = 1; j <= i; ++j) {
                    R.at(i).at(j) = R.at(i).at(j - 1) + (R.at(i).at(j - 1) - R.at(i - 1).at(j - 1)) / (std::pow(4.0, j) - 1);
                }

                // Check for convergence.
                if (i > 1 && std::abs(R.at(i).at(i) - R.at(i - 1).at(i - 1)) < tolerance) {
                    std::cout << "Iterations: " << i << std::endl;
                    return R.at(i).at(i);
                }
            }
            // Throw an exception if the method does not converge within the maximum number of iterations.
            throw std::runtime_error("RombergMethod did not converge within the specified number of iterations.");
        }
    };

    /**
     * @brief Class implementing the Adaptive Simpson's Integration Method.
     *
     * Adaptive Simpson's Method is a numerical technique to approximate the definite integral of a function,
     * which adapts the size of the intervals based on the function's behavior.
     */
    class Simpson
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
        template< IsFunction FN, std::floating_point TOL_T = double, std::integral ITER_T = int >
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
        template< IsFunction FN, std::floating_point TOL_T = double, std::integral ITER_T = int >
        auto adaptiveSimpson(FN                       function,
                             std::floating_point auto lower,
                             std::floating_point auto upper,
                             TOL_T                    tolerance,
                             ITER_T                   maxDepth,
                             ITER_T                   depth) const
        {
            double m  = (lower + upper) / 2;    // Midpoint of the interval
            double h  = (upper - lower) / 2;    // Half the interval length
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
     * @brief Class implementing the Trapezoid Integration Method.
     *
     * The Trapezoid Method is a numerical technique for approximating the definite integral of a function.
     */
    class Trapezoid
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
        template< IsFunction FN, std::floating_point TOL_T = double, std::integral ITER_T = int >
        auto operator()(FN                       function,
                        std::floating_point auto lower,
                        std::floating_point auto upper,
                        TOL_T                    tolerance     = 1e-6,
                        ITER_T                   maxIterations = 100) const
        {
            // Initial step size and initial trapezoidal estimate
            ReturnType< decltype(function) > h    = (upper - lower);
            ReturnType< decltype(function) > Iold = 0.5 * h * (function(lower) + function(upper));

            // Main loop for adaptive trapezoidal rule
            for (int i = 1; i <= maxIterations; ++i) {
                h *= 0.5;                                   // Halve the step size in each iteration
                int numOfMidpoints = std::pow(2, i - 1);    // Calculate the number of midpoints for this iteration

                // Create a vector of midpoints
                std::vector< double > midPoints(numOfMidpoints);
                std::generate(midPoints.begin(), midPoints.end(), [n = 1, lower, h]() mutable { return lower + (2 * n++ - 1) * h; });

                // Calculate the sum of function values at the midpoints
                ReturnType< decltype(function) > sum =
                    std::transform_reduce(midPoints.begin(), midPoints.end(), 0.0, std::plus<>(), [&](double x) { return function(x); });

                // Update the integral estimate
                ReturnType< decltype(function) > Inew = 0.5 * Iold + h * sum;

                // Check for convergence
                if (i > 1 && std::abs(Inew - Iold) < tolerance) {
                    std::cout << "Iterations: " << i << std::endl;    // Output the number of iterations for debugging
                    return Inew;                                      // Return the new estimate if converged
                }
                Iold = Inew;    // Update the old estimate for the next iteration
            }

            // Throw an exception if the method does not converge within the maximum number of iterations
            throw std::runtime_error("AdaptiveTrapezoidMethod did not converge within the specified number of iterations.");
        }
    };

    template< typename ALGO = Romberg, IsFunction FN, std::floating_point TOL_T = double, std::integral ITER_T = int >
    auto integrate(FN                       function,
                   std::floating_point auto lower,
                   std::floating_point auto upper,
                   TOL_T                    tolerance     = 1e-6,
                   ITER_T                   maxIterations = 100)
    {
        if (!function) throw std::runtime_error("Function object is invalid.");
        return ALGO {}(function, lower, upper, tolerance, maxIterations);
    }

    template< typename ALGO = Romberg >
    inline auto integralOf(IsFunction auto function)
    {
        if (!function) throw std::runtime_error("Function object is invalid.");

        using RT = ReturnType< decltype(function) >;
        return [=](std::floating_point auto a, std::floating_point auto b) { return ALGO {}(function, a, b); };
    }
}    // namespace nxx::integrate

#endif    // NUMERIXX_INTEGRATION_IMP_HPP
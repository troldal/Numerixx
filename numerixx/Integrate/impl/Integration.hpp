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
#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
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

    class RombergMethod
    {
    public:
        inline auto operator()(IsFunction auto                  function,
                               ReturnType< decltype(function) > a,
                               ReturnType< decltype(function) > b,
                               double                           tolerance     = 1e-6,
                               int                              maxIterations = 1000) const -> ReturnType< decltype(function) >
        {
            std::vector< std::vector< ReturnType< decltype(function) > > > R(
                1,
                std::vector< ReturnType< decltype(function) > >(1, 0.5 * (b - a) * (function(a) + function(b))));

            for (int i = 1; i <= maxIterations; ++i) {
                R.push_back(std::vector< ReturnType< decltype(function) > >(i + 1, 0));
                double h = (b - a) / std::pow(2.0, i);

                // Corrected trapezoidal rule calculation for subsequent iterations
                double sum = 0.0;
                for (int k = 1; k <= std::pow(2, i - 1); ++k) {
                    sum += function(a + (2 * k - 1) * h);
                }
                R[i][0] = 0.5 * R[i - 1][0] + h * sum;

                for (int j = 1; j <= i; ++j) {
                    R[i][j] = R[i][j - 1] + (R[i][j - 1] - R[i - 1][j - 1]) / (std::pow(4.0, j) - 1);
                }

                if (i > 1 && std::abs(R[i][i] - R[i - 1][i - 1]) < tolerance) {
                    return R[i][i];
                }
            }
            throw NumerixxError("RombergMethod did not converge within the specified number of iterations.");
        }
    };

    class SimpsonsMethod
    {
    public:
        inline auto operator()(IsFunction auto function, ReturnType< decltype(function) > a, ReturnType< decltype(function) > b) const
            -> ReturnType< decltype(function) >
        {
            ReturnType< decltype(function) > h = (b - a) / 2;
            return (h / 3) * (function(a) + 4 * function(a + h) + function(b));
        }
    };

    class TrapezoidMethod
    {
    public:
        inline auto operator()(IsFunction auto                  function,
                               ReturnType< decltype(function) > a,
                               ReturnType< decltype(function) > b,
                               double                           tolerance     = 1e-6,
                               int                              maxIterations = 1000) const -> ReturnType< decltype(function) >
        {
            ReturnType< decltype(function) > h    = (b - a);
            ReturnType< decltype(function) > Iold = 0.5 * h * (function(a) + function(b));
            for (int i = 1; i <= maxIterations; ++i) {
                h *= 0.5;
                ReturnType< decltype(function) > sum = 0;
                for (int j = 1; j <= std::pow(2, i - 1); ++j) {
                    sum += function(a + (2 * j - 1) * h);
                }
                ReturnType< decltype(function) > Inew = 0.5 * Iold + h * sum;
                if (i > 1 && std::abs(Inew - Iold) < tolerance) {
                    return Inew;
                }
                Iold = Inew;
            }
            throw NumerixxError("AdaptiveTrapezoidMethod did not converge within the specified number of iterations.");
        }
    };

    template< typename ALGO >
    inline auto integrate(ALGO algorithm, IsFunction auto function, ReturnType< decltype(function) > a, ReturnType< decltype(function) > b)
    {
        if (!function) throw NumerixxError("Function object is invalid.");
        return algorithm(function, a, b);
    }

    template< typename ALGO = SimpsonsMethod >
    inline auto integralOf(IsFunction auto function, ReturnType< decltype(function) > a)
    {
        if (!function) throw NumerixxError("Function object is invalid.");

        using RT = ReturnType< decltype(function) >;
        return [=](RT b) { return ALGO {}(function, a, b); };
    }
}    // namespace nxx::integrate

#endif    // NUMERIXX_INTEGRATION_IMP_HPP

/*
    888b      88  88        88  88b           d88  88888888888  88888888ba   88  8b        d8  8b        d8
    8888b     88  88        88  888b         d888  88           88      "8b  88   Y8,    ,8P    Y8,    ,8P
    88 `8b    88  88        88  88`8b       d8'88  88           88      ,8P  88    `8b  d8'      `8b  d8'
    88  `8b   88  88        88  88 `8b     d8' 88  88aaaaa      88aaaaaa8P'  88      Y88P          Y88P
    88   `8b  88  88        88  88  `8b   d8'  88  88"""""      88""""88'    88      d88b          d88b
    88    `8b 88  88        88  88   `8b d8'   88  88           88    `8b    88    ,8P  Y8,      ,8P  Y8,
    88     `8888  Y8a.    .a8P  88    `888'    88  88           88     `8b   88   d8'    `8b    d8'    `8b
    88      `888   `"Y8888Y"'   88     `8'     88  88888888888  88      `8b  88  8P        Y8  8P        Y8

    Copyright © 2024 Kenneth Troldal Balslev

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

#pragma once

// ===== Numerixx Includes
#include <Concepts.hpp>
#include <Deriv.hpp>
#include <Poly.hpp>

// ===== External Includes
#include <blaze/Blaze.h>

// ===== Standard Library Includes
#include <algorithm>
#include <cassert>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <vector>

/**
 * @file Interpolate.hpp
 * @brief Comprehensive interpolation library for various algorithms and utilities.
 *
 * This header file includes a suite of classes and functions designed for interpolation tasks.
 * It encompasses various interpolation methods, utility functions for creating interpolation
 * objects, and a function for generating a polynomial that fits a set of points.
 *
 * The following components are included:
 *
 * 1. Base class for interpolation algorithms (InterpBase) using the CRTP pattern.
 * 2. Derived classes for specific interpolation methods:
 *    - Linear interpolation (Linear)
 *    - Lagrange polynomial interpolation (Lagrange)
 *    - Steffen interpolation (Steffen)
 *    - Cubic spline interpolation (Spline)
 * 3. Utility functions for interpolation:
 *    - interpolate: Functions to interpolate at a single point using a specified algorithm.
 *    - interpolationOf: Functions to return interpolation function objects based on given points and an algorithm.
 *    - makepoly: Function to create a polynomial passing through a set of points.
 *
 * Each class and function is documented with details on usage, parameters, and return types.
 * This library requires the Blaze library for some of the matrix and vector operations.
 *
 * @note This library is designed to be flexible and extensible, allowing for easy addition of new interpolation methods.
 */

namespace nxx::interp
{
    namespace detail
    {

        // =================================================================================================================
        //
        //    88                                                           88888888ba
        //    88                ,d                                         88      "8b
        //    88                88                                         88      ,8P
        //    88  8b,dPPYba,  MM88MMM  ,adPPYba,  8b,dPPYba,  8b,dPPYba,   88aaaaaa8P'  ,adPPYYba,  ,adPPYba,   ,adPPYba,
        //    88  88P'   `"8a   88    a8P_____88  88P'   "Y8  88P'    "8a  88""""""8b,  ""     `Y8  I8[    ""  a8P_____88
        //    88  88       88   88    8PP"""""""  88          88       d8  88      `8b  ,adPPPPP88   `"Y8ba,   8PP"""""""
        //    88  88       88   88,   "8b,   ,aa  88          88b,   ,a8"  88      a8P  88,    ,88  aa    ]8I  "8b,   ,aa
        //    88  88       88   "Y888  `"Ybbd8"'  88          88`YbbdP"'   88888888P"   `"8bbdP"Y8  `"YbbdP"'   `"Ybbd8"'
        //                                                    88
        //                                                    88
        //
        // =================================================================================================================

        /**
         * @brief Provides a CRTP base class for interpolation of a series of points.
         *
         * The InterpBase class template is designed to be used as a base class for various interpolation
         * techniques. It uses the Curiously Recurring Template Pattern (CRTP) to allow derived classes
         * to specify their interpolation method.
         *
         * @tparam DERIVED The derived class implementing the interpolation logic.
         * @tparam CONTAINER_T The type of container used to store the points.
         * @tparam POINT_T The type of points stored in the container. Must be a floating point structure.
         *
         * @note Requires that POINT_T is a floating point structure.
         */
        template< typename DERIVED, IsFloatStruct POINT_T >
        class InterpBase
        {
            friend DERIVED;

            std::vector< POINT_T > m_points;    ///< Container holding the points for interpolation.

        public:
            static constexpr bool IsInterpolator = true;    ///< Flag to identify interpolator classes.

        protected:
            ~InterpBase() = default;                          ///< Default destructor.
            using VALUE_T = StructCommonType_t< POINT_T >;    ///< Common value type derived from POINT_T.

        public:
            InterpBase(const InterpBase& other)     = default;    ///< Default copy constructor.
            InterpBase(InterpBase&& other) noexcept = default;    ///< Default move constructor.

            InterpBase& operator=(const InterpBase& other)     = default;    ///< Default copy assignment operator.
            InterpBase& operator=(InterpBase&& other) noexcept = default;    ///< Default move assignment operator.

            /**
             * @brief Constructs an interpolation base object with the given points.
             *
             * @param points The container of points to be used for interpolation.
             * @throws std::runtime_error If the number of points is less than 2.
             */
            template <template< typename... > class CONTAINER_T >
            explicit InterpBase(const CONTAINER_T< POINT_T >& points)
                : m_points(points.begin(), points.end())
            {
                if (m_points.size() < 2) throw std::runtime_error("Interpolation requires at least two points.");
                std::sort(m_points.begin(), m_points.end(), [](const auto& p1, const auto& p2) {
                    auto [x1, y1] = p1;
                    auto [x2, y2] = p2;
                    return x1 < x2;
                });
            }

            /**
             * @brief Constructs an interpolation base object with the given points.
             *
             * @tparam N The number of points in the initializer list.
             * @param points The initializer list of points to be used for interpolation.
             * @throws std::runtime_error If the number of points is less than 2.
             */
            template< size_t N >
            requires(N >= 2)
            explicit InterpBase(const POINT_T (&points)[N])
                : m_points(points, points + N)
            {
                if (m_points.size() < 2) throw std::runtime_error("Interpolation requires at least two points.");
                std::sort(m_points.begin(), m_points.end(), [](const auto& p1, const auto& p2) {
                    auto [x1, y1] = p1;
                    auto [x2, y2] = p2;
                    return x1 < x2;
                });
            }

            /**
             * @brief Constructs an interpolation base object with the given x- and y-coordinates.
             *
             * @tparam XCONT_T The type of container used to store the x-coordinates.
             * @tparam YCONT_T The type of container used to store the y-coordinates.
             * @tparam XVAL_T The type of values stored in the x-container.
             * @tparam YVAL_T The type of values stored in the y-container.
             * @param x The container of x-coordinates to be used for interpolation.
             * @param y The container of y-coordinates to be used for interpolation.
             * @throws std::runtime_error If the number of points is less than 2.
             * @throws std::runtime_error If the number of x- and y-coordinates are not equal.
             */
            template< template< typename... > class XCONT_T, template< typename... > class YCONT_T, typename XVAL_T, typename YVAL_T >
            InterpBase(const XCONT_T< XVAL_T >& x, const YCONT_T< YVAL_T >& y)
            {
                if (x.size() != y.size()) throw std::runtime_error("Interpolation requires equal number of x and y points.");
                if (x.size() < 2) throw std::runtime_error("Interpolation requires at least two points.");
                for (size_t i = 0; i < x.size(); ++i) {
                    m_points.push_back({ x[i], y[i] });
                }
                std::sort(m_points.begin(), m_points.end(), [](const auto& p1, const auto& p2) {
                    auto [x1, y1] = p1;
                    auto [x2, y2] = p2;
                    return x1 < x2;
                });
            }

            /**
             * @brief Constructs an interpolation base object with the given x- and y-coordinates.
             *
             * @tparam XVAL_T The type of values stored in the x-container.
             * @tparam YVAL_T The type of values stored in the y-container.
             * @tparam N The number of points in the initializer list.
             * @param x The initializer list of x-coordinates to be used for interpolation.
             * @param y The initializer list of y-coordinates to be used for interpolation.
             * @throws std::runtime_error If the number of points is less than 2.
             * @throws std::runtime_error If the number of x- and y-coordinates are not equal.
             */
            template< typename XVAL_T, size_t N, typename YVAL_T, size_t M >
            requires(N == M && N >= 2)
            InterpBase(const XVAL_T (&x)[N], const YVAL_T (&y)[M])
            {
                for (size_t i = 0; i < N; ++i) {
                    m_points.push_back({ x[i], y[i] });
                }
                std::sort(m_points.begin(), m_points.end(), [](const auto& p1, const auto& p2) {
                    auto [x1, y1] = p1;
                    auto [x2, y2] = p2;
                    return x1 < x2;
                });
            }

            /**
             * @brief Function call operator to perform interpolation.
             *
             * This method delegates the actual interpolation task to the derived class by calling
             * its interpolate method.
             *
             * @param x The x-coordinate value for which interpolation is to be performed.
             * @return VALUE_T The interpolated value at the specified x-coordinate.
             * @throws std::runtime_error If the interpolation point is out of bounds.
             */
            VALUE_T operator()(VALUE_T x) const
            {
                if (auto [xa, _] = m_points.front(); x < xa) throw std::runtime_error("Interpolation point is out of bounds.");
                if (auto [xb, _] = m_points.back(); x > xb) throw std::runtime_error("Interpolation point is out of bounds.");
                return static_cast< const DERIVED& >(*this).interpolate(x);
            }
        };
    }    // namespace detail

    // =================================================================================================================
    //
    //    88           88
    //    88           ""
    //    88
    //    88           88  8b,dPPYba,    ,adPPYba,  ,adPPYYba,  8b,dPPYba,
    //    88           88  88P'   `"8a  a8P_____88  ""     `Y8  88P'   "Y8
    //    88           88  88       88  8PP"""""""  ,adPPPPP88  88
    //    88           88  88       88  "8b,   ,aa  88,    ,88  88
    //    88888888888  88  88       88   `"Ybbd8"'  `"8bbdP"Y8  88
    //
    // =================================================================================================================

    /**
     * @brief Provides a derived class for linear interpolation using the CRTP pattern.
     *
     * The Linear class template extends the InterpBase class for performing linear interpolation
     * and extrapolation of a series of points.
     *
     * @tparam CONTAINER_T The type of container used to store the points.
     * @tparam POINT_T The type of points stored in the container.
     */
    template< typename POINT_T >
    class Linear : public detail::InterpBase< Linear< POINT_T >, POINT_T >
    {
        using BASE = detail::InterpBase< Linear< POINT_T >, POINT_T >;

    public:
        using BASE::BASE;                          ///< Inherits constructors from the base class.
        using VALUE_T = typename BASE::VALUE_T;    ///< Type for value representation, inherited from base class.

        /**
         * @brief Performs linear interpolation for a given x-value.
         *
         * This method finds the interval that the x-value falls into and performs linear interpolation
         * using the points defining that interval.
         *
         * @param x The x-coordinate value for which interpolation is to be performed.
         * @return VALUE_T The interpolated value at the specified x-coordinate.
         * @throws std::out_of_range If x is out of the range of the points.
         */
        VALUE_T interpolate(VALUE_T x) const
        {
            // Find the interval that x falls into
            auto it = std::lower_bound(BASE::m_points.begin(), BASE::m_points.end(), x, [](const auto& p, double x) {
                auto [x1, _] = p;
                return x1 <= x;
            });

            // Check if x is equal to the x-value of the last point
            if (it == BASE::m_points.end()) {
                auto [_, y] = BASE::m_points.back();
                return y;    // Return the y-value of the last point
            }

            // Assert that the interval is valid
            assert(std::distance(BASE::m_points.begin(), it) >= 1);
            assert(std::distance(it, BASE::m_points.end()) >= 1);

            // Perform linear interpolation
            auto [x1, y1] = *(it - 1);
            auto [x2, y2] = *it;
            return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
        }

        /**
         * @brief Performs linear extrapolation for a given x-value.
         *
         * This method handles cases where the x-value is outside the range of the points,
         * using linear extrapolation based on the slope of the closest interval.
         *
         * @param x The x-coordinate value for which extrapolation is to be performed.
         * @return VALUE_T The extrapolated value at the specified x-coordinate.
         * @throws std::runtime_error If the number of points is less than 2.
         * @throws std::out_of_range If x is out of the range of the points and extrapolation is not possible.
         */
        VALUE_T extrapolate(VALUE_T x) const
        {
            // Handle the case where x is before the first point
            if (x < BASE::m_points.front().first) {
                auto [x1, y1] = BASE::m_points[0];
                auto [x2, y2] = BASE::m_points[1];
                double slope  = (y2 - y1) / (x2 - x1);
                return y1 + slope * (x - x1);
            }

            // Handle the case where x is after the last point
            if (x > BASE::m_points.back().first) {
                auto [x1, y1] = BASE::m_points[BASE::m_points.size() - 2];
                auto [x2, y2] = BASE::m_points.back();
                double slope  = (y2 - y1) / (x2 - x1);
                return y1 + slope * (x - x1);
            }

            // If x is within the range, use linear interpolation
            return interpolate(x);
        }
    };

    /*
     * Deduction guides for Linear class with a container of POINT_T.
     */
    template< template< typename... > class CONTAINER_T, typename POINT_T >
    Linear(CONTAINER_T< POINT_T >) -> Linear< POINT_T >;

    template< size_t N >
    requires(N >= 2)
    Linear(const std::pair< double, double > (&)[N]) -> Linear< std::pair< double, double > >;

    template< template< typename... > class XCONT_T, template< typename... > class YCONT_T, typename XVAL_T, typename YVAL_T >
    Linear(const XCONT_T< XVAL_T >& x, const YCONT_T< YVAL_T >& y) -> Linear< std::pair< XVAL_T, YVAL_T > >;

    template< typename XVAL_T, size_t N, typename YVAL_T, size_t M >
    requires(N == M && N >= 2)
    Linear(const XVAL_T (&)[N], const YVAL_T (&)[M]) -> Linear< std::pair< XVAL_T, YVAL_T > >;

    // =================================================================================================================
    //
    //    88
    //    88
    //    88
    //    88           ,adPPYYba,   ,adPPYb,d8  8b,dPPYba,  ,adPPYYba,  8b,dPPYba,    ,adPPYb,d8   ,adPPYba,
    //    88           ""     `Y8  a8"    `Y88  88P'   "Y8  ""     `Y8  88P'   `"8a  a8"    `Y88  a8P_____88
    //    88           ,adPPPPP88  8b       88  88          ,adPPPPP88  88       88  8b       88  8PP"""""""
    //    88           88,    ,88  "8a,   ,d88  88          88,    ,88  88       88  "8a,   ,d88  "8b,   ,aa
    //    88888888888  `"8bbdP"Y8   `"YbbdP"Y8  88          `"8bbdP"Y8  88       88   `"YbbdP"Y8   `"Ybbd8"'
    //                              aa,    ,88                                        aa,    ,88
    //                               "Y8bbdP"                                          "Y8bbdP"
    //
    // =================================================================================================================

    /**
     * @brief Provides a derived class for Lagrange polynomial interpolation using the CRTP pattern.
     *
     * The Lagrange class template extends the InterpBase class for performing Lagrange polynomial
     * interpolation and extrapolation of a series of points.
     *
     * @tparam CONTAINER_T The type of container used to store the points.
     * @tparam POINT_T The type of points stored in the container.
     */
    template< typename POINT_T >
    class Lagrange : public detail::InterpBase< Lagrange< POINT_T >, POINT_T >
    {
        using BASE    = detail::InterpBase< Lagrange< POINT_T >, POINT_T >;
        using VALUE_T = typename BASE::VALUE_T;

        /**
         * @brief Implements the Lagrange interpolation algorithm.
         *
         * This method calculates the Lagrange polynomial based on the provided points
         * and evaluates it at the given x-coordinate.
         *
         * @param x The x-coordinate value for which interpolation is to be performed.
         * @return VALUE_T The interpolated value at the specified x-coordinate.
         * @throws std::runtime_error If no points are provided for interpolation.
         */
        VALUE_T implementation(VALUE_T x) const
        {
            // Extrapolate using the slopes between the outer points:
            //            // Check if x is outside the bounds of the data points
            //            if (x < BASE::m_points.front().first) {
            //                // Linear extrapolation using the slope of the first interval
            //                auto [x0, y0] = BASE::m_points[0];
            //                auto [x1, y1] = BASE::m_points[1];
            //                VALUE_T slope = (y1 - y0) / (x1 - x0);
            //                return y0 + slope * (x - x0);
            //            } else if (x > BASE::m_points.back().first) {
            //                // Linear extrapolation using the slope of the last interval
            //                size_t n = BASE::m_points.size();
            //                auto [xn_1, yn_1] = BASE::m_points[n - 2];
            //                auto [xn, yn] = BASE::m_points[n - 1];
            //                VALUE_T slope = (yn - yn_1) / (xn - xn_1);
            //                return yn + slope * (x - xn);
            //            } else {
            //                // Lagrange interpolation
            //                VALUE_T result = 0;
            //                for (size_t j = 0; j < BASE::m_points.size(); ++j) {
            //                    VALUE_T term = BASE::m_points[j].second;
            //                    for (size_t m = 0; m < BASE::m_points.size(); ++m) {
            //                        if (m != j) {
            //                            term *= (x - BASE::m_points[m].first) / (BASE::m_points[j].first - BASE::m_points[m].first);
            //                        }
            //                    }
            //                    result += term;
            //                }
            //                return result;
            //            }

            // Extrapolate using the slopes of the Lagrange polynomial at the outer points:
            // Check if x is outside the bounds of the data points
            if (x < BASE::m_points.front().first) {
                // Linear extrapolation using the slope of the Lagrange polynomial at the first point
                VALUE_T x0    = BASE::m_points.front().first;
                VALUE_T y0    = BASE::m_points.front().second;
                VALUE_T slope = *nxx::deriv::forward(*this, x0);
                return y0 + slope * (x - x0);
            }
            else if (x > BASE::m_points.back().first) {
                // Linear extrapolation using the slope of the Lagrange polynomial at the last point
                VALUE_T xn    = BASE::m_points.back().first;
                VALUE_T yn    = BASE::m_points.back().second;
                VALUE_T slope = *nxx::deriv::backward(*this, xn);
                return yn + slope * (x - xn);
            }
            else {
                // Regular Lagrange interpolation
                VALUE_T result = 0;
                for (size_t j = 0; j < BASE::m_points.size(); ++j) {
                    VALUE_T term = BASE::m_points[j].second;
                    for (size_t m = 0; m < BASE::m_points.size(); ++m) {
                        if (m != j) {
                            term *= (x - BASE::m_points[m].first) / (BASE::m_points[j].first - BASE::m_points[m].first);
                        }
                    }
                    result += term;
                }
                return result;
            }
        }

    public:
        using BASE::BASE;    ///< Inherits constructors from the base class.

        /**
         * @brief Performs Lagrange interpolation for a given x-value.
         *
         * This method checks if the x-value is within the bounds of the points and then
         * uses the implementation method to perform the interpolation.
         *
         * @param x The x-coordinate value for which interpolation is to be performed.
         * @return VALUE_T The interpolated value at the specified x-coordinate.
         * @throws std::runtime_error If the interpolation point is out of bounds.
         */
        VALUE_T interpolate(VALUE_T x) const { return implementation(x); }

        /**
         * @brief Performs Lagrange extrapolation for a given x-value.
         *
         * This method uses the implementation method to perform extrapolation irrespective
         * of whether the x-value is within the bounds of the points.
         *
         * @param x The x-coordinate value for which extrapolation is to be performed.
         * @return VALUE_T The extrapolated value at the specified x-coordinate.
         */
        VALUE_T extrapolate(VALUE_T x) const { return implementation(x); }
    };

    /*
     * Deduction guides for Lagrange class with a container of POINT_T.
     */
    template< template< typename... > class CONTAINER_T, typename POINT_T >
    Lagrange(CONTAINER_T< POINT_T >) -> Lagrange< POINT_T >;

    template< size_t N >
    requires(N >= 2)
    Lagrange(const std::pair< double, double > (&)[N]) -> Lagrange< std::pair< double, double > >;

    template< template< typename... > class XCONT_T, template< typename... > class YCONT_T, typename XVAL_T, typename YVAL_T >
    Lagrange(const XCONT_T< XVAL_T >& x, const YCONT_T< YVAL_T >& y) -> Lagrange< std::pair< XVAL_T, YVAL_T > >;

    template< typename XVAL_T, size_t N, typename YVAL_T, size_t M >
    requires(N == M && N >= 2)
    Lagrange(const XVAL_T (&)[N], const YVAL_T (&)[M]) -> Lagrange< std::pair< XVAL_T, YVAL_T > >;

    // =================================================================================================================
    //
    //     ad88888ba                        ad88     ad88
    //    d8"     "8b  ,d                  d8"      d8"
    //    Y8,          88                  88       88
    //    `Y8aaaaa,  MM88MMM  ,adPPYba,  MM88MMM  MM88MMM  ,adPPYba,  8b,dPPYba,
    //      `"""""8b,  88    a8P_____88    88       88    a8P_____88  88P'   `"8a
    //            `8b  88    8PP"""""""    88       88    8PP"""""""  88       88
    //    Y8a     a8P  88,   "8b,   ,aa    88       88    "8b,   ,aa  88       88
    //     "Y88888P"   "Y888  `"Ybbd8"'    88       88     `"Ybbd8"'  88       88
    //
    // =================================================================================================================

    /**
     * @brief Provides a derived class for interpolation using the Steffen method.
     *
     * The Steffen class template extends the InterpBase class for performing Steffen interpolation
     * and extrapolation of a series of points. The Steffen method uses a Hermite interpolation
     * approach to ensure smoothness at the data points.
     *
     * @tparam CONTAINER_T The type of container used to store the points.
     * @tparam POINT_T The type of points stored in the container.
     */
    template< typename POINT_T >
    class Steffen : public detail::InterpBase< Steffen< POINT_T >, POINT_T >
    {
        using BASE    = detail::InterpBase< Steffen< POINT_T >, POINT_T >;
        using VALUE_T = typename BASE::VALUE_T;

        mutable std::optional< std::vector< VALUE_T > > m_slopes;    ///< Slopes at each point for Hermite interpolation.

    public:
        using BASE::BASE;    ///< Inherits constructors from the base class.

        /**
         * @brief Calculates slopes for Hermite interpolation at each point.
         *
         * This method computes the slopes used in the Hermite interpolation formula,
         * ensuring a smooth transition at each point.
         */
        void calculateSlopes() const
        {
            m_slopes     = std::vector< VALUE_T > {};
            auto& slopes = m_slopes.value();

            size_t n = BASE::m_points.size();
            slopes.resize(n);

            // Calculate slopes
            for (size_t i = 0; i < n; ++i) {
                if (i == 0 || i == n - 1) {    // Endpoint slopes
                    slopes[i] = (BASE::m_points[i == 0 ? 1 : n - 1].second - BASE::m_points[i == 0 ? 0 : n - 2].second) /
                                (BASE::m_points[i == 0 ? 1 : n - 1].first - BASE::m_points[i == 0 ? 0 : n - 2].first);
                }
                else {    // Intermediate slopes
                    VALUE_T slope1 =
                        (BASE::m_points[i].second - BASE::m_points[i - 1].second) / (BASE::m_points[i].first - BASE::m_points[i - 1].first);
                    VALUE_T slope2 =
                        (BASE::m_points[i + 1].second - BASE::m_points[i].second) / (BASE::m_points[i + 1].first - BASE::m_points[i].first);
                    slopes[i] = (slope1 * slope2 <= 0) ? 0 : 2 / (1 / slope1 + 1 / slope2);
                }
            }
        }

        /**
         * @brief Performs Steffen interpolation for a given x-value.
         *
         * This method uses Hermite interpolation based on the slopes calculated at each point
         * to interpolate a value at the given x-coordinate.
         *
         * @param x The x-coordinate value for which interpolation is to be performed.
         * @return VALUE_T The interpolated value at the specified x-coordinate.
         * @throws std::runtime_error If the interpolation point is out of bounds.
         */
        VALUE_T interpolate(VALUE_T x) const
        {
            if (!m_slopes) calculateSlopes();
            auto& slopes = m_slopes.value();

            // Find the interval that x falls into
            auto it =
                std::upper_bound(BASE::m_points.begin(), BASE::m_points.end(), x, [](double x, const auto& p) { return x < p.first; });

            auto [x1, y1]  = *(it - 1);
            auto [x2, y2]  = *it;
            VALUE_T slope1 = slopes[std::distance(BASE::m_points.begin(), it) - 1];
            VALUE_T slope2 = slopes[std::distance(BASE::m_points.begin(), it)];

            // Hermite interpolation
            VALUE_T t   = (x - x1) / (x2 - x1);
            VALUE_T h00 = (1 + 2 * t) * (1 - t) * (1 - t);
            VALUE_T h10 = t * (1 - t) * (1 - t);
            VALUE_T h01 = t * t * (3 - 2 * t);
            VALUE_T h11 = t * t * (t - 1);

            return h00 * y1 + h10 * slope1 * (x2 - x1) + h01 * y2 + h11 * slope2 * (x2 - x1);
        }

        /**
         * @brief Performs Steffen extrapolation for a given x-value.
         *
         * This method handles cases where the x-value is outside the range of the points,
         * using linear extrapolation based on the slope at the closest endpoint.
         *
         * @param x The x-coordinate value for which extrapolation is to be performed.
         * @return VALUE_T The extrapolated value at the specified x-coordinate.
         */
        VALUE_T extrapolate(VALUE_T x) const
        {
            size_t n = BASE::m_points.size();

            // Extrapolate at the beginning
            if (x <= BASE::m_points.front().first) {
                auto [x0, y0] = BASE::m_points[0];
                VALUE_T slope = *nxx::deriv::forward(*this, x0);
                return y0 + slope * (x - x0);
            }

            // Extrapolate at the end
            if (x >= BASE::m_points.back().first) {
                auto [xn, yn] = BASE::m_points[n - 1];
                VALUE_T slope = *nxx::deriv::backward(*this, xn);
                return yn + slope * (x - xn);
            }

            // If x is within the range, use interpolation
            return interpolate(x);
        }
    };

    /*
     * Deduction guides for Steffen class with a container of POINT_T.
     */
    template< template< typename... > class CONTAINER_T, typename POINT_T >
    Steffen(CONTAINER_T< POINT_T >) -> Steffen< POINT_T >;

    template< size_t N >
    requires(N >= 2)
    Steffen(const std::pair< double, double > (&)[N]) -> Steffen< std::pair< double, double > >;

    template< template< typename... > class XCONT_T, template< typename... > class YCONT_T, typename XVAL_T, typename YVAL_T >
    Steffen(const XCONT_T< XVAL_T >& x, const YCONT_T< YVAL_T >& y) -> Steffen< std::pair< XVAL_T, YVAL_T > >;

    template< typename XVAL_T, size_t N, typename YVAL_T, size_t M >
    requires(N == M && N >= 2)
    Steffen(const XVAL_T (&)[N], const YVAL_T (&)[M]) -> Steffen< std::pair< XVAL_T, YVAL_T > >;

    // =================================================================================================================
    //
    //     ad88888ba                88  88
    //    d8"     "8b               88  ""
    //    Y8,                       88
    //    `Y8aaaaa,    8b,dPPYba,   88  88  8b,dPPYba,    ,adPPYba,
    //      `"""""8b,  88P'    "8a  88  88  88P'   `"8a  a8P_____88
    //            `8b  88       d8  88  88  88       88  8PP"""""""
    //    Y8a     a8P  88b,   ,a8"  88  88  88       88  "8b,   ,aa
    //     "Y88888P"   88`YbbdP"'   88  88  88       88   `"Ybbd8"'
    //                 88
    //                 88
    //
    // =================================================================================================================

    /**
     * @brief Provides a derived class for cubic spline interpolation.
     *
     * The Spline class template extends the InterpBase class for performing cubic spline
     * interpolation of a series of points. It computes and stores coefficients for each cubic
     * spline segment and uses them to interpolate values at arbitrary points.
     *
     * @tparam CONTAINER_T The type of container used to store the points.
     * @tparam POINT_T The type of points stored in the container.
     */
    template< typename POINT_T >
    class Spline : public detail::InterpBase< Spline< POINT_T >, POINT_T >
    {
        using BASE    = detail::InterpBase< Spline< POINT_T >, POINT_T >;
        using VALUE_T = typename BASE::VALUE_T;

        /**
         * @struct SplineCoefficients
         * @brief Structure to store the coefficients of the cubic spline segments.
         */
        struct SplineCoefficients
        {
            std::vector< VALUE_T > a;    ///< Coefficients of the zeroth degree.
            std::vector< VALUE_T > b;    ///< Coefficients of the first degree.
            std::vector< VALUE_T > c;    ///< Coefficients of the second degree.
            std::vector< VALUE_T > d;    ///< Coefficients of the third degree.
        };

        mutable std::optional< SplineCoefficients > m_coefficients;    ///< Instance variable to hold the coefficients for the cubic spline.

    public:
        using BASE::BASE;    ///< Inherits constructors from the base class.

        /**
         * @brief Calculates coefficients for the cubic spline.
         *
         * This method computes the coefficients for each cubic spline segment based on the provided points.
         *
         * @throws std::runtime_error If there are less than three points, which is insufficient for cubic spline interpolation.
         */
        void calculateSplineCoefficients() const
        {
            m_coefficients     = SplineCoefficients {};
            auto& coefficients = m_coefficients.value();

            size_t n = BASE::m_points.size() - 1;    // Number of intervals

            // Vectors for spline calculation
            blaze::DynamicVector< VALUE_T > a(n + 1), b(n), d(n + 1), h(n);

            // Initializing 'a' with y-values and calculating 'h' (intervals)
            for (size_t i = 0; i < n; ++i) {
                a[i] = BASE::m_points[i].second;
                h[i] = BASE::m_points[i + 1].first - BASE::m_points[i].first;
            }
            a[n] = BASE::m_points[n].second;

            // Calculating the alpha vector used in the spline calculation
            blaze::DynamicVector< VALUE_T > alpha(n);
            for (size_t i = 1; i < n; ++i) {
                alpha[i] = 3.0 / h[i] * (a[i + 1] - a[i]) - 3.0 / h[i - 1] * (a[i] - a[i - 1]);
            }

            // Vectors used in solving the tridiagonal system for 'c'
            blaze::DynamicVector< VALUE_T > c(n + 1), l(n + 1), mu(n + 1), z(n + 1);
            l[0]  = 1.0;
            mu[0] = z[0] = 0.0;

            // Forward sweep to solve for 'c'
            for (size_t i = 1; i < n; ++i) {
                l[i]  = 2.0 * (BASE::m_points[i + 1].first - BASE::m_points[i - 1].first) - h[i - 1] * mu[i - 1];
                mu[i] = h[i] / l[i];
                z[i]  = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
            }

            l[n] = 1.0;
            z[n] = c[n] = 0.0;

            // Backward sweep for 'c', and solving for 'b' and 'd'
            for (int j = n - 1; j >= 0; --j) { // NOLINT
                c[j] = z[j] - mu[j] * c[j + 1];
                b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
                d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
            }

            // Storing the calculated coefficients
            coefficients.a = std::vector< VALUE_T >(a.begin(), a.end());
            coefficients.b = std::vector< VALUE_T >(b.begin(), b.end());
            coefficients.c = std::vector< VALUE_T >(c.begin(), c.end());
            coefficients.d = std::vector< VALUE_T >(d.begin(), d.end());
        }

        /**
         * @brief Interpolates a value at a given x-coordinate using the spline.
         *
         * This method finds the appropriate spline segment for the given x-coordinate and evaluates
         * the cubic polynomial to find the interpolated value.
         *
         * @param x The x-coordinate value for which interpolation is to be performed.
         * @return VALUE_T The interpolated value at the specified x-coordinate.
         */
        VALUE_T interpolate(VALUE_T x) const
        {
            if (!m_coefficients) calculateSplineCoefficients();

            auto& a = m_coefficients->a;
            auto& b = m_coefficients->b;
            auto& c = m_coefficients->c;
            auto& d = m_coefficients->d;

            // Find the right interval for x
            size_t interval = BASE::m_points.size() - 2;    // Default to the last interval
            for (size_t i = 0; i < BASE::m_points.size() - 1; ++i) {
                if (x <= BASE::m_points[i + 1].first) {
                    interval = i;
                    break;
                }
            }

            // Calculate the difference between x and the starting x of the interval
            VALUE_T dx = x - BASE::m_points[interval].first;

            VALUE_T result = a[interval] + b[interval] * dx + c[interval] * dx * dx + d[interval] * dx * dx * dx;
            return result;
        }

        /**
         * @brief Extrapolates a value at a given x-coordinate using the spline.
         *
         * Currently, this method behaves the same as interpolate, but it could be extended
         * for specific extrapolation behavior in the future.
         *
         * @param x The x-coordinate value for which extrapolation is to be performed.
         * @return VALUE_T The extrapolated value at the specified x-coordinate.
         */
        VALUE_T extrapolate(VALUE_T x) const
        {
            size_t n = BASE::m_points.size();

            // Extrapolate at the beginning
            if (x <= BASE::m_points.front().first) {
                auto [x0, y0] = BASE::m_points[0];
                VALUE_T slope = *nxx::deriv::forward(*this, x0);
                return y0 + slope * (x - x0);
            }

            // Extrapolate at the end
            if (x >= BASE::m_points.back().first) {
                auto [xn, yn] = BASE::m_points[n - 1];
                VALUE_T slope = *nxx::deriv::backward(*this, xn);
                return yn + slope * (x - xn);
            }

            return interpolate(x);
        }
    };

    /*
     * Deduction guides for Spline class with a container of POINT_T.
     */
    template< template< typename... > class CONTAINER_T, typename POINT_T >
    Spline(CONTAINER_T< POINT_T >) -> Spline< POINT_T >;

    template< size_t N >
    requires(N >= 2)
    Spline(const std::pair< double, double > (&)[N]) -> Spline< std::pair< double, double > >;

    template< template< typename... > class XCONT_T, template< typename... > class YCONT_T, typename XVAL_T, typename YVAL_T >
    Spline(const XCONT_T< XVAL_T >& x, const YCONT_T< YVAL_T >& y) -> Spline< std::pair< XVAL_T, YVAL_T > >;

    template< typename XVAL_T, size_t N, typename YVAL_T, size_t M >
    requires(N == M && N >= 2)
    Spline(const XVAL_T (&)[N], const YVAL_T (&)[M]) -> Spline< std::pair< XVAL_T, YVAL_T > >;

    // =================================================================================================================
    //
    //    88                                                                        88
    //    ""                ,d                                                      88                ,d
    //                      88                                                      88                88
    //    88  8b,dPPYba,  MM88MMM  ,adPPYba,  8b,dPPYba,  8b,dPPYba,    ,adPPYba,   88  ,adPPYYba,  MM88MMM  ,adPPYba,
    //    88  88P'   `"8a   88    a8P_____88  88P'   "Y8  88P'    "8a  a8"     "8a  88  ""     `Y8    88    a8P_____88
    //    88  88       88   88    8PP"""""""  88          88       d8  8b       d8  88  ,adPPPPP88    88    8PP"""""""
    //    88  88       88   88,   "8b,   ,aa  88          88b,   ,a8"  "8a,   ,a8"  88  88,    ,88    88,   "8b,   ,aa
    //    88  88       88   "Y888  `"Ybbd8"'  88          88`YbbdP"'    `"YbbdP"'   88  `"8bbdP"Y8    "Y888  `"Ybbd8"'
    //                                                    88
    //                                                    88
    //
    // =================================================================================================================

    /**
     * @brief Interpolates at a single point using a specified interpolation algorithm.
     *
     * This function template creates an instance of the specified interpolation algorithm and
     * uses it to interpolate at the given x-coordinate.
     *
     * @tparam ALGO The interpolation algorithm class template.
     * @tparam CONTAINER_T The type of container used to store the points.
     * @tparam POINT_T The type of points stored in the container.
     * @param points The container of points to be used for interpolation.
     * @param x The x-coordinate value for which interpolation is to be performed.
     * @return The interpolated value at the specified x-coordinate.
     *
     * @note Requires that ALGO< CONTAINER_T, POINT_T > has a static member IsInterpolator set to true.
     */
    template< template< typename > class ALGO, template< typename... > class CONTAINER_T, typename POINT_T >
    requires ALGO< POINT_T >::IsInterpolator
    auto interpolate(const CONTAINER_T< POINT_T >& points, double x)
    {
        auto algo = ALGO(points);
        return algo(x);
    }

    /**
     * @brief Interpolates at a single point using a specified interpolation algorithm, with points provided as an initializer list.
     *
     * This function template is a specialization for the case where points are provided as an initializer list.
     * It creates an instance of the specified interpolation algorithm and uses it to interpolate at the given x-coordinate.
     *
     * @tparam ALGO The interpolation algorithm class template.
     * @tparam N The number of points provided in the initializer list.
     * @param points An initializer list of points to be used for interpolation.
     * @param x The x-coordinate value for which interpolation is to be performed.
     * @return The interpolated value at the specified x-coordinate.
     *
     * @note Requires that ALGO< std::vector, std::pair< double, double > > has a static member IsInterpolator set to true.
     */
    template< template< typename > class ALGO, size_t N >
    requires ALGO< std::pair< double, double > >::IsInterpolator
    auto interpolate(const std::pair< double, double > (&points)[N], double x)
    {
        auto algo = ALGO(points);
        return algo(x);
    }

    /**
     * @brief Provides a function template for interpolation with separate x and y containers.
     *
     * This function template is designed for cases where x and y values are stored in separate containers.
     * It enables interpolation at a specified x-coordinate using a chosen interpolation algorithm. It combines
     * the x and y values into pairs, creates an instance of the specified interpolation algorithm with these
     * pairs, and then performs interpolation at the specified x-coordinate.
     *
     * @tparam ALGO The interpolation algorithm class template.
     * @tparam XCONT_T The type of container used to store the x-coordinate values.
     * @tparam YCONT_T The type of container used to store the y-coordinate values.
     * @tparam XVAL_T The type of values stored in the x-coordinates container.
     * @tparam YVAL_T The type of values stored in the y-coordinates container.
     * @param x A container of x-coordinate values.
     * @param y A container of y-coordinate values corresponding to the x-coordinates.
     * @param xval The x-coordinate value for which interpolation is to be performed.
     * @return The interpolated value at the specified x-coordinate.
     *
     * @note Requires that ALGO< std::vector, std::pair< XVAL_T, YVAL_T > > has a static member IsInterpolator set to true.
     */
    template< template< typename > class ALGO,
              template< typename... >
              class XCONT_T,
              template< typename... >
              class YCONT_T,
              typename XVAL_T,
              typename YVAL_T >
    requires ALGO< std::pair< XVAL_T, YVAL_T > >::IsInterpolator
    auto interpolate(const XCONT_T< XVAL_T >& x, const YCONT_T< YVAL_T >& y, double xval)
    {
        auto algo = ALGO(x, y);
        return algo(xval);
    }

    /**
     * @brief Provides a function template for interpolation with separate x and y arrays.
     *
     * This function template is a specialization for cases where x and y values are stored in separate arrays.
     * It enables interpolation at a specified x-coordinate using a chosen interpolation algorithm. It combines
     * the x and y values into pairs, creates an instance of the specified interpolation algorithm with these
     * pairs, and then performs interpolation at the specified x-coordinate.
     *
     * @tparam ALGO The interpolation algorithm class template.
     * @tparam XVAL_T The type of values stored in the x-coordinates array.
     * @tparam N The number of elements in the x-coordinates array.
     * @tparam YVAL_T The type of values stored in the y-coordinates array.
     * @tparam M The number of elements in the y-coordinates array.
     * @param x An array of x-coordinate values.
     * @param y An array of y-coordinate values corresponding to the x-coordinates.
     * @param xval The x-coordinate value for which interpolation is to be performed.
     * @return The interpolated value at the specified x-coordinate.
     *
     * @note Requires that ALGO< std::vector, std::pair< XVAL_T, YVAL_T > > has a static member IsInterpolator set to true.
     */
    template< template< typename > class ALGO, typename XVAL_T, size_t N, typename YVAL_T, size_t M >
    requires ALGO< std::pair< XVAL_T, YVAL_T > >::IsInterpolator
    auto interpolate(const XVAL_T (&x)[N], const YVAL_T (&y)[M], double xval)
    {
        auto algo = ALGO(x, y);
        return algo(xval);
    }

    // =================================================================================================================
    //
    //    88                                                             ,ad8888ba,       ad88
    //    ""                ,d                                          d8"'    `"8b     d8"
    //                      88                                         d8'        `8b    88
    //    88  8b,dPPYba,  MM88MMM  ,adPPYba,  8b,dPPYba,  8b,dPPYba,   88          88  MM88MMM
    //    88  88P'   `"8a   88    a8P_____88  88P'   "Y8  88P'    "8a  88          88    88
    //    88  88       88   88    8PP"""""""  88          88       d8  Y8,        ,8P    88
    //    88  88       88   88,   "8b,   ,aa  88          88b,   ,a8"   Y8a.    .a8P     88
    //    88  88       88   "Y888  `"Ybbd8"'  88          88`YbbdP"'     `"Y8888Y"'      88
    //                                                    88
    //                                                    88
    //
    // =================================================================================================================

    /**
     * @brief Returns an interpolation function object based on the given points and algorithm.
     *
     * This function template creates an instance of the specified interpolation algorithm with the provided
     * points. The returned object can be used to perform interpolation at different x-coordinate values.
     *
     * @tparam ALGO The interpolation algorithm class template.
     * @tparam CONTAINER_T The type of container used to store the points.
     * @tparam POINT_T The type of points stored in the container.
     * @param points The container of points to be used for creating the interpolation function object.
     * @return An instance of the specified interpolation algorithm initialized with the given points.
     *
     * @note Requires that ALGO< CONTAINER_T, POINT_T > has a static member IsInterpolator set to true.
     */
    template< template< typename > class ALGO, template< typename... > class CONTAINER_T, typename POINT_T >
    requires ALGO< POINT_T >::IsInterpolator
    auto interpolationOf(const CONTAINER_T< POINT_T >& points)
    {
        return ALGO(points);
    }

    /**
     * @brief Returns an interpolation function object based on the given points and algorithm, with points provided as an initializer list.
     *
     * This function template is a specialization for cases where points are provided as an initializer list.
     * It creates an instance of the specified interpolation algorithm with the provided points.
     * The returned object can be used to perform interpolation at different x-coordinate values.
     *
     * @tparam ALGO The interpolation algorithm class template.
     * @tparam N The number of points provided in the initializer list.
     * @param points An initializer list of points to be used for creating the interpolation function object.
     * @return An instance of the specified interpolation algorithm initialized with the given points.
     *
     * @note Requires that ALGO< std::vector, std::pair< double, double > > has a static member IsInterpolator set to true.
     */
    template< template< typename > class ALGO, size_t N >
    requires ALGO< std::pair< double, double > >::IsInterpolator
    auto interpolationOf(const std::pair< double, double > (&points)[N])
    {
        return ALGO(points);
    }

    /**
     * @brief Returns an interpolation function object based on the given points and algorithm, with points provided as separate x and y
     * containers.
     *
     * This function template is designed for cases where x and y values are stored in separate containers.
     * It creates an instance of the specified interpolation algorithm with the provided points.
     * The returned object can be used to perform interpolation at different x-coordinate values.
     *
     * @tparam ALGO The interpolation algorithm class template.
     * @tparam XCONT_T The type of container used to store the x-coordinate values.
     * @tparam YCONT_T The type of container used to store the y-coordinate values.
     * @tparam XVAL_T The type of values stored in the x-coordinates container.
     * @tparam YVAL_T The type of values stored in the y-coordinates container.
     * @param x A container of x-coordinate values.
     * @param y A container of y-coordinate values corresponding to the x-coordinates.
     * @return An instance of the specified interpolation algorithm initialized with the given points.
     *
     * @note Requires that ALGO< std::vector, std::pair< XVAL_T, YVAL_T > > has a static member IsInterpolator set to true.
     */
    template< template< typename > class ALGO,
              template< typename... >
              class XCONT_T,
              template< typename... >
              class YCONT_T,
              typename XVAL_T,
              typename YVAL_T >
    requires ALGO< std::pair< XVAL_T, YVAL_T > >::IsInterpolator
    auto interpolationOf(const XCONT_T< XVAL_T >& x, const YCONT_T< YVAL_T >& y)
    {
        return ALGO(x, y);
    }

    /**
     * @brief Returns an interpolation function object based on the given points and algorithm, with points provided as separate x and y
     * arrays.
     *
     * This function template is a specialization for cases where x and y values are stored in separate arrays.
     * It creates an instance of the specified interpolation algorithm with the provided points.
     * The returned object can be used to perform interpolation at different x-coordinate values.
     *
     * @tparam ALGO The interpolation algorithm class template.
     * @tparam XVAL_T The type of values stored in the x-coordinates array.
     * @tparam N The number of elements in the x-coordinates array.
     * @tparam YVAL_T The type of values stored in the y-coordinates array.
     * @tparam M The number of elements in the y-coordinates array.
     * @param x An array of x-coordinate values.
     * @param y An array of y-coordinate values corresponding to the x-coordinates.
     * @return An instance of the specified interpolation algorithm initialized with the given points.
     *
     * @note Requires that ALGO< std::vector, std::pair< XVAL_T, YVAL_T > > has a static member IsInterpolator set to true.
     */
    template< template< typename > class ALGO, typename XVAL_T, size_t N, typename YVAL_T, size_t M >
    requires ALGO< std::pair< XVAL_T, YVAL_T > >::IsInterpolator
    auto interpolationOf(const XVAL_T (&x)[N], const YVAL_T (&y)[M])
    {
        return ALGO(x, y);
    }

    // =================================================================================================================
    //
    //                                  88                                              88
    //                                  88                                              88
    //                                  88                                              88
    //  88,dPYba,,adPYba,   ,adPPYYba,  88   ,d8   ,adPPYba,  8b,dPPYba,    ,adPPYba,   88  8b       d8
    //  88P'   "88"    "8a  ""     `Y8  88 ,a8"   a8P_____88  88P'    "8a  a8"     "8a  88  `8b     d8'
    //  88      88      88  ,adPPPPP88  8888[     8PP"""""""  88       d8  8b       d8  88   `8b   d8'
    //  88      88      88  88,    ,88  88`"Yba,  "8b,   ,aa  88b,   ,a8"  "8a,   ,a8"  88    `8b,d8'
    //  88      88      88  `"8bbdP"Y8  88   `Y8a  `"Ybbd8"'  88`YbbdP"'    `"YbbdP"'   88      Y88'
    //                                                        88                                d8'
    //                                                        88                               d8'
    //
    // =================================================================================================================

    /**
     * @brief Creates a polynomial that passes through all the given points.
     *
     * This function sets up and solves a linear system to find the coefficients of a polynomial
     * that passes through the given set of points. It constructs a polynomial using these coefficients.
     *
     * @param points A vector of pairs, where each pair represents a point (x, y).
     * @return An instance of the Polynomial class representing the polynomial that fits all the given points.
     *
     * @note This function requires the blaze library for matrix and vector operations.
     */
    inline auto makepoly(const std::vector< std::pair< double, double > >& points)
    {
        size_t                         n = points.size();
        blaze::DynamicMatrix< double > A(n, n);
        blaze::DynamicVector< double > b(n);
        blaze::DynamicVector< double > x(n);

        // Setting up the system Ax = b
        for (size_t i = 0; i < n; ++i) {
            double xi = points[i].first;
            b[i]      = points[i].second;
            for (size_t j = 0; j < n; ++j) {
                A(i, j) = std::pow(xi, j);
            }
        }

        // Solving the system
        x = blaze::solve(A, b);

        // Extracting coefficients
        std::vector< double > coefficients(n);
        for (size_t i = 0; i < n; ++i) {
            coefficients[i] = x[i];
        }

        return nxx::poly::Polynomial(coefficients);
    }

}    // namespace nxx::interp

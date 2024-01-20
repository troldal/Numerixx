//
// Created by kenne on 13/11/2023.
//

#ifndef NUMERIXX_INTERPOLATE_HPP
#define NUMERIXX_INTERPOLATE_HPP

//
// Created by I22696 on 13-11-2023.
//

#include <Concepts.hpp>
#include <Poly.hpp>
#include <algorithm>
#include <blaze/Blaze.h>
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace nxx::interp
{

    template< typename DERIVED, template< typename... > class CONTAINER_T, typename POINT_T >
    requires IsFloatStruct< POINT_T >
    class InterpolationBase
    {
        friend DERIVED;

        CONTAINER_T< POINT_T > m_points;

    public:
        static constexpr bool IsInterpolator = true;

    protected:
        ~InterpolationBase() = default;
        using VALUE_T        = StructCommonType_t< POINT_T >;

    public:
        explicit InterpolationBase(const CONTAINER_T< POINT_T >& points)
            : m_points(points)
        {
            if (m_points.size() < 2) throw std::runtime_error("Linear interpolation requires at least two points.");
            std::sort(m_points.begin(), m_points.end(), [](const auto& p1, const auto& p2) {
                auto [x1, y1] = p1;
                auto [x2, y2] = p2;
                return x1 < x2;
            });
        }

        VALUE_T operator()(VALUE_T x) const
        {
            if (auto [xa, _] = m_points.front(); x < xa) throw std::runtime_error("Interpolation point is out of bounds.");
            if (auto [xb, _] = m_points.back(); x > xb) throw std::runtime_error("Interpolation point is out of bounds.");
            return static_cast< const DERIVED& >(*this).interpolate(x);
        }
    };

    template< template< typename... > class CONTAINER_T, typename POINT_T >
    class Linear : public InterpolationBase< Linear< CONTAINER_T, POINT_T >, CONTAINER_T, POINT_T >
    {
        using BASE = InterpolationBase< Linear< CONTAINER_T, POINT_T >, CONTAINER_T, POINT_T >;

    public:
        using BASE::BASE;
        using VALUE_T = typename BASE::VALUE_T;

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

            assert(std::distance(BASE::m_points.begin(), it) >= 1);
            assert(std::distance(it, BASE::m_points.end()) >= 1);

            // Perform linear interpolation
            auto [x1, y1] = *(it - 1);
            auto [x2, y2] = *it;
            return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
        }

        VALUE_T extrapolate(VALUE_T x) const
        {
            // Check if the points array is too small
            if (BASE::m_points.size() < 2) {
                throw std::runtime_error("Insufficient points for extrapolation");
            }

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

    template< template< typename... > class CONTAINER_T, typename POINT_T >
    Linear(CONTAINER_T< POINT_T >) -> Linear< CONTAINER_T, POINT_T >;

    Linear(std::initializer_list< std::pair< double, double > >) -> Linear< std::vector, std::pair< double, double > >;

    template< template< typename... > class CONTAINER_T, typename POINT_T >
    class Lagrange : public InterpolationBase< Lagrange< CONTAINER_T, POINT_T >, CONTAINER_T, POINT_T >
    {
        using BASE    = InterpolationBase< Lagrange< CONTAINER_T, POINT_T >, CONTAINER_T, POINT_T >;
        using VALUE_T = typename BASE::VALUE_T;

        VALUE_T implementation(VALUE_T x) const
        {
            if (BASE::m_points.empty()) {
                throw std::runtime_error("No points provided for interpolation.");
            }

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

    public:
        using BASE::BASE;

        VALUE_T interpolate(VALUE_T x) const
        {
            auto [x1, y1] = BASE::m_points.front();
            auto [x2, y2] = BASE::m_points.back();
            if (x < x1 || x > x2) throw std::runtime_error("Interpolation point is out of bounds.");
            return implementation(x);
        }

        VALUE_T extrapolate(VALUE_T x) const { return implementation(x); }
    };

    template< template< typename... > class CONTAINER_T, typename POINT_T >
    Lagrange(CONTAINER_T< POINT_T >) -> Lagrange< CONTAINER_T, POINT_T >;

    Lagrange(std::initializer_list< std::pair< double, double > >) -> Lagrange< std::vector, std::pair< double, double > >;

    template< template< typename... > class CONTAINER_T, typename POINT_T >
    class Steffen : public InterpolationBase< Steffen< CONTAINER_T, POINT_T >, CONTAINER_T, POINT_T >
    {
        using BASE    = InterpolationBase< Steffen< CONTAINER_T, POINT_T >, CONTAINER_T, POINT_T >;
        using VALUE_T = typename BASE::VALUE_T;

        mutable std::vector< VALUE_T > slopes;

    public:
        using BASE::BASE;

        // Constructor to initialize the slopes
        explicit Steffen(const CONTAINER_T< POINT_T >& points)
            : BASE(points)
        {
            calculateSlopes();
        }

        void calculateSlopes() const
        {
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

        VALUE_T interpolate(VALUE_T x) const
        {
            // Find the interval that x falls into
            auto it =
                std::upper_bound(BASE::m_points.begin(), BASE::m_points.end(), x, [](double x, const auto& p) { return x < p.first; });

            if (it == BASE::m_points.begin() || it == BASE::m_points.end()) {
                throw std::runtime_error("Interpolation point is out of bounds.");
            }

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

        VALUE_T extrapolate(VALUE_T x) const
        {
            size_t n = BASE::m_points.size();

            // Extrapolate at the beginning
            if (x < BASE::m_points.front().first) {
                auto [x0, y0] = BASE::m_points[0];
                auto [x1, y1] = BASE::m_points[1];
                VALUE_T slope = (y1 - y0) / (x1 - x0);
                return y0 + slope * (x - x0);
            }

            // Extrapolate at the end
            if (x > BASE::m_points.back().first) {
                auto [xn_1, yn_1] = BASE::m_points[n - 2];
                auto [xn, yn]     = BASE::m_points[n - 1];
                VALUE_T slope     = (yn - yn_1) / (xn - xn_1);
                return yn + slope * (x - xn);
            }

            // If x is within the range, use interpolation
            return interpolate(x);
        }
    };

    template< template< typename... > class CONTAINER_T, typename POINT_T >
    Steffen(CONTAINER_T< POINT_T >) -> Steffen< CONTAINER_T, POINT_T >;

    Steffen(std::initializer_list< std::pair< double, double > >) -> Steffen< std::vector, std::pair< double, double > >;

template<template<typename...> class CONTAINER_T, typename POINT_T>
class Spline : public InterpolationBase<Spline<CONTAINER_T, POINT_T>, CONTAINER_T, POINT_T>
{
    using BASE = InterpolationBase<Spline<CONTAINER_T, POINT_T>, CONTAINER_T, POINT_T>;
    using VALUE_T = typename BASE::VALUE_T;

    struct SplineSegment {
        VALUE_T a, b, c, d, x;
    };

    mutable std::vector<SplineSegment> segments;

public:
    using BASE::BASE;

    Spline(const CONTAINER_T<POINT_T>& points) : BASE(points) {
        calculateSplineCoefficients();
    }

    void calculateSplineCoefficients() const {
        size_t n = BASE::m_points.size() - 1;
        if (n < 1) {
            throw std::runtime_error("Not enough points for spline interpolation");
        }

        blaze::DynamicVector<VALUE_T> h(n), alpha(n), l(n+1), mu(n), z(n+1), c(n+1), b(n), d(n);

        // Compute h and alpha
        for (size_t i = 0; i < n; ++i) {
            h[i] = BASE::m_points[i+1].first - BASE::m_points[i].first;
            alpha[i] = (3 / h[i]) * (BASE::m_points[i+1].second - BASE::m_points[i].second)
                     - (3 / h[i-1]) * (BASE::m_points[i].second - BASE::m_points[i-1].second);
        }

        // Compute l, mu, z
        l[0] = 1;
        mu[0] = z[0] = 0;
        for (size_t i = 1; i < n; ++i) {
            l[i] = 2 * (BASE::m_points[i+1].first - BASE::m_points[i-1].first) - h[i-1] * mu[i-1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i-1] * z[i-1]) / l[i];
        }
        l[n] = 1;
        z[n] = 0;

        // Back substitution loop for the cubic coefficients
        c[n] = 0; // Natural spline: second derivative of n is zero
        for (size_t j = n-1; j-- > 0; ) {
            c[j] = z[j] - mu[j] * c[j+1];
            b[j] = (BASE::m_points[j+1].second - BASE::m_points[j].second) / h[j]
                   - h[j] * (c[j+1] + 2 * c[j]) / 3;
            d[j] = (c[j+1] - c[j]) / (3 * h[j]);
        }

        // Store coefficients
        for (size_t i = 0; i < n; ++i) {
            segments[i] = {BASE::m_points[i].second, b[i], c[i], d[i], BASE::m_points[i].first};
        }
    }

    VALUE_T interpolate(VALUE_T x) const {
        // Find the right segment
        auto it = std::lower_bound(segments.begin(), segments.end(), x,
                                   [](const SplineSegment& seg, VALUE_T x) { return seg.x < x; });
        if (it != segments.begin()) {
            --it;
        }

        // Compute the interpolated value
        VALUE_T dx = x - it->x;
        return it->a + it->b * dx + it->c * dx * dx + it->d * dx * dx * dx;
    }
};

    template< template< typename... > class CONTAINER_T, typename POINT_T >
    Spline(CONTAINER_T< POINT_T >) -> Spline< CONTAINER_T, POINT_T >;

    Spline(std::initializer_list< std::pair< double, double > >) -> Spline< std::vector, std::pair< double, double > >;

    struct SplineCoefficients
    {
        std::vector< double > a, b, c, d;
    };

    class CubicSplineInterp
    {
    public:
        template< typename Point >
        auto operator()(const std::vector< Point >& points) const
        {
            if (points.size() < 3) throw std::runtime_error("Cubic spline interpolation requires at least three points.");

            size_t                         n = points.size() - 1;
            blaze::DynamicVector< double > a(n + 1), b(n), d(n + 1), h(n);

            for (size_t i = 0; i < n; ++i) {
                a[i] = points[i].second;
                h[i] = points[i + 1].first - points[i].first;
            }
            a[n] = points[n].second;

            blaze::DynamicVector< double > alpha(n);
            for (size_t i = 1; i < n; ++i) {
                alpha[i] = 3.0 / h[i] * (a[i + 1] - a[i]) - 3.0 / h[i - 1] * (a[i] - a[i - 1]);
            }

            blaze::DynamicVector< double > c(n + 1), l(n + 1), mu(n + 1), z(n + 1);
            l[0]  = 1.0;
            mu[0] = z[0] = 0.0;

            for (size_t i = 1; i < n; ++i) {
                l[i]  = 2.0 * (points[i + 1].first - points[i - 1].first) - h[i - 1] * mu[i - 1];
                mu[i] = h[i] / l[i];
                z[i]  = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
            }

            l[n] = 1.0;
            z[n] = c[n] = 0.0;

            for (int j = n - 1; j >= 0; --j) {
                c[j] = z[j] - mu[j] * c[j + 1];
                b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
                d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
            }

            // The coefficients of the cubic splines are in a, b, c, d
            // For simplicity, only returning 'a', modify as needed
            // std::vector<double> coefficients(a.begin(), a.end());

            SplineCoefficients coefficients;

            coefficients.a = std::vector< double >(a.begin(), a.end());
            coefficients.b = std::vector< double >(b.begin(), b.end());
            coefficients.c = std::vector< double >(c.begin(), c.end());
            coefficients.d = std::vector< double >(d.begin(), d.end());

            return coefficients;
        }
    };

    double evaluateSpline(const std::vector< std::pair< double, double > >& points,
                          const std::vector< double >&                      a,
                          const std::vector< double >&                      b,
                          const std::vector< double >&                      c,
                          const std::vector< double >&                      d,
                          double                                            x)
    {
        if (points.size() < 2) throw std::runtime_error("Insufficient points for spline evaluation.");

        // Find the right interval for x
        size_t interval = points.size() - 2;    // Default to the last interval
        for (size_t i = 0; i < points.size() - 1; ++i) {
            if (x <= points[i + 1].first) {
                interval = i;
                break;
            }
        }

        // Calculate the difference between x and the starting x of the interval
        double dx = x - points[interval].first;

        // Evaluate the cubic polynomial
        double result = a[interval] + b[interval] * dx + c[interval] * dx * dx + d[interval] * dx * dx * dx;
        return result;
    }

    template< typename ALGO, typename Point >
    double interpolate(ALGO algorithm, const std::vector< Point >& points, double x)
    {
        return algorithm(points, x);
    }

    template< typename ALGO, typename Point >
auto interpolationOf(ALGO algorithm, const std::vector< Point >& points)
    {
        return [=](double x) { return algorithm(points, x); };
    }

    // =================================================================================================================
    //
    //                                 88                                              88
    //                                 88                                              88
    //                                 88                                              88
    // 88,dPYba,,adPYba,   ,adPPYYba,  88   ,d8   ,adPPYba,  8b,dPPYba,    ,adPPYba,   88  8b       d8
    // 88P'   "88"    "8a  ""     `Y8  88 ,a8"   a8P_____88  88P'    "8a  a8"     "8a  88  `8b     d8'
    // 88      88      88  ,adPPPPP88  8888[     8PP"""""""  88       d8  8b       d8  88   `8b   d8'
    // 88      88      88  88,    ,88  88`"Yba,  "8b,   ,aa  88b,   ,a8"  "8a,   ,a8"  88    `8b,d8'
    // 88      88      88  `"8bbdP"Y8  88   `Y8a  `"Ybbd8"'  88`YbbdP"'    `"YbbdP"'   88      Y88'
    //                                                       88                                d8'
    //                                                       88                               d8'
    //
    // =================================================================================================================

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

#endif    // NUMERIXX_INTERPOLATE_HPP

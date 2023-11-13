//
// Created by kenne on 13/11/2023.
//

#ifndef NUMERIXX_INTERPOLATE_HPP
#define NUMERIXX_INTERPOLATE_HPP

//
// Created by I22696 on 13-11-2023.
//

#include <algorithm>
#include <blaze/Blaze.h>
#include <stdexcept>
#include <vector>

namespace nxx::interp
{

    class LinearInterp
    {
    public:
        // Assumes points is a vector of pairs, where each pair is (x, y)
        template< typename Point >
        double operator()(const std::vector< Point >& points, double x) const
        {
            if (points.size() < 2) throw std::runtime_error("Linear interpolation requires at least two points.");

            // Find the interval that x falls into
            auto it = std::lower_bound(points.begin(), points.end(), x, [](const Point& p, double x) { return p.first < x; });

            // Handle edge cases: x is out of bounds
            if (it == points.begin() || it == points.end()) throw std::runtime_error("Interpolation point is out of bounds.");

            // Perform linear interpolation
            auto [x1, y1] = *(it - 1);
            auto [x2, y2] = *it;
            return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
        }
    };

    class PolynomialInterp
    {
    public:
        template< typename Point >
        double operator()(const std::vector< Point >& points, double x) const
        {
            if (points.empty()) throw std::runtime_error("No points provided for interpolation.");

            double result = 0;
            for (size_t i = 0; i < points.size(); ++i) {
                double term = points[i].second;
                for (size_t j = 0; j < points.size(); ++j) {
                    if (j != i) {
                        term *= (x - points[j].first) / (points[i].first - points[j].first);
                    }
                }
                result += term;
            }
            return result;
        }
    };

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

}    // namespace nxx::interp

#endif    // NUMERIXX_INTERPOLATE_HPP

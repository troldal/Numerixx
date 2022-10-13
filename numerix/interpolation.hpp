//
// Created by Kenneth Balslev on 22/01/2021.
//

#ifndef PCPROPS_INTERPOLATION_HPP
#define PCPROPS_INTERPOLATION_HPP

#include <algorithm>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace numeric
{
    class Interpolator
    {
        std::vector<std::pair<double, double> > m_points {};

    public:

        Interpolator(std::vector<std::pair<double, double> > points) : m_points(points) {
            std::sort(m_points.begin(), m_points.end(), [](const auto& a, const auto& b) {return a.first < b.first;} );
        }

        Interpolator(std::vector<double> x, std::vector<double> y) {
            if (x.size() != y.size()) throw std::invalid_argument("Vectors of x's and y's must be of same length.");
            for (auto [xi, yi] = std::tuple {x.begin(), y.begin()}; xi != x.end(); ++xi, ++yi)
                m_points.emplace_back(std::make_pair(*xi, *yi));
            std::sort(m_points.begin(), m_points.end(), [](const auto& a, const auto& b) {return a.first < b.first;} );
        }

        Interpolator(const Interpolator& other) = default;

        Interpolator(Interpolator&& other) noexcept = default;

        ~Interpolator() = default;

        Interpolator& operator=(const Interpolator& other) = default;

        Interpolator& operator=(Interpolator&& other) noexcept = default;

        double operator()(double arg) {
            if (arg < m_points[0].first || arg > m_points.back().first) throw std::invalid_argument("Argument outside valid range");
            auto upper = std::upper_bound(m_points.begin(), m_points.end(), arg, [](double a, const std::pair<double, double>& b) {return a < b.first;} );
            auto lower = upper - 1;

            auto slope = (upper->second - lower->second) / (upper->first - lower->first);

            return lower->second + slope * (arg - lower->first);
        }

    };
}    // namespace numeric

#endif    // PCPROPS_INTERPOLATION_HPP

//
// Created by Kenneth Balslev on 16/01/2021.
//

#ifndef PCPROPS_INTEGRATION_HPP
#define PCPROPS_INTEGRATION_HPP

#include <cmath>
#include <functional>

namespace numeric {

    /**
     * @brief
     * @param func
     * @param x1
     * @param x2
     * @param precision
     * @return
     */
    inline double integrate(const std::function<double(double)>& func, double x1, double x2, double precision = 1E-6)
    {
        using std::abs;
        using std::min;
        using std::pow;
        using Point = std::pair<double, double>;

        // ===== Internal recursive lambda.
        // ===== If the precision is reached, the segment integral (area) is returned. Otherwise, the function is called on two sub-segments.
        std::function<double(const Point&, const Point&, const Point&)> calcSegmentIntegral = [&] (const Point& lower, const Point& upper, const Point& mid)
        {
               if (abs(func(mid.first) - mid.second) / func(mid.first) > (precision == 0.0 ? 1E-6 : precision)) {
                   auto real_mid = std::make_pair(mid.first, func(mid.first));
                   return calcSegmentIntegral(lower, real_mid, std::make_pair((lower.first + real_mid.first) / 2, (lower.second + real_mid.second) / 2)) +
                          calcSegmentIntegral(real_mid, upper, std::make_pair((real_mid.first + upper.first) / 2, (real_mid.second + upper.second) / 2));
               }

               return 0.5 * (upper.first - lower.first) * (lower.second + upper.second);
        };

        // ===== Ensure that the function is subdivided to at least 10 segments.
        auto step = (x2 - x1) / 10;
        double result = 0.0;
        for (int i = 0; i < 10; ++i) {
            auto lower = std::make_pair(x1 + i * step, func(x1 + i * step));
            auto upper = std::make_pair(x1 + (i + 1) * step, func(x1 + (i + 1) * step));
            auto mid   = std::make_pair((lower.first + upper.first) / 2, (lower.second + upper.second) / 2);
            result += calcSegmentIntegral(lower, upper, mid);
        }

        return result;
    }

} // namespace numeric


#endif    // PCPROPS_INTEGRATION_HPP

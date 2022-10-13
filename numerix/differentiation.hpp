//
// Created by Kenneth Balslev on 16/01/2021.
//

#ifndef PCPROPS_DIFFERENTIATION_HPP
#define PCPROPS_DIFFERENTIATION_HPP

#include <functional>
#include <cmath>

namespace numeric {

    /**
     * @brief
     * @param func
     * @param x
     * @return
     */
    inline double diff_central(const std::function<double(double)>& func, double x)
    {
        using std::sqrt;

        auto h = sqrt(std::numeric_limits<double>::epsilon());

        return (func(x + h) - func(x - h)) / (2 * h);
    }

    /**
     * @brief
     * @param func
     * @param x
     * @return
     */
    inline double diff_backward(const std::function<double(double)>& func, double x)
    {
        using std::sqrt;

        auto h = sqrt(std::numeric_limits<double>::epsilon());

        return (func(x) - func(x - h)) / h;
    }

} // namespace numeric

#endif    // PCPROPS_DIFFERENTIATION_HPP

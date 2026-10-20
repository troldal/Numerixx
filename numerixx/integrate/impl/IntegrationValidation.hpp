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

#pragma once

// ===== Numerixx Includes
#include <Error.hpp>

/**
 * @file IntegrationValidation.hpp
 * @brief Header file defining the validateRange function for checking integration bounds.
 *
 * This file contains the definition of the validateRange function, which is used to validate
 * the range of integration bounds. It ensures that the lower bound is less than the upper bound,
 * throwing an exception if this condition is not met. This function is crucial for ensuring
 * the correctness and stability of numerical integration processes.
 */
namespace nxx::integrate::detail
{
    /**
     * @brief Validates the range of integration bounds.
     *
     * @details This function checks if the lower bound is less than the upper bound for an integration range.
     *          If the lower bound is not less than the upper bound, it throws a NumerixxError.
     *
     * @tparam T The type of the bounds, constrained to floating-point types.
     * @param lower The lower bound of the range.
     * @param upper The upper bound of the range.
     * @throw NumerixxError if the lower bound is not less than the upper bound.
     */
    void validateRange(IsFloat auto lower, IsFloat auto upper)
    {
        if (lower >= upper) throw NumerixxError("The lower bound must be less than the upper bound.");
    }
} // namespace nxx::integrate::detail



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

#ifndef INTEGRATIONERROR_HPP
#define INTEGRATIONERROR_HPP

/**
 * @file IntegrationError.hpp
 * @brief Header file defining the IntegrationErrorData structure.
 *
 * This file contains the definition of the IntegrationErrorData structure, which is used to
 * encapsulate detailed error information related to numerical integration processes.
 * It provides a comprehensive view of the error state, including the last computed value,
 * absolute and relative errors, and the number of iterations performed.
 *
 * This structure is particularly useful for diagnosing issues in numerical integration routines,
 * allowing users to assess the quality of the result and the convergence behavior of the algorithm.
 */
namespace nxx::integrate::detail
{
    /**
     * @struct IntegrationErrorData
     * @brief Struct to hold detailed error information for integration processes.
     *
     * @details This structure encapsulates error-related information for numerical integration routines.
     *          It stores the last computed value, absolute and relative errors, and the number of iterations.
     *
     * @note This structure is used internally as a template parameter to the nxx::Error class.
     *
     * @tparam T The data type of the value and error measurements (typically a floating-point type).
     * @tparam ITER_T The data type for counting iterations (typically an integral type).
     */
    template<typename T, typename ITER_T>
    struct IntegrationErrorData
    {
        T      value;      /**< The last computed value of the integration process. */
        T      eabs;       /**< The absolute error of the computed value. */
        T      erel;       /**< The relative error of the computed value. */
        ITER_T iterations; /**< The total number of iterations performed. */

        /**
         * @brief Overloads the output stream operator for IntegrationErrorData.
         *
         * @details This function allows IntegrationErrorData instances to be easily printed to an output stream.
         *          It formats the data with high precision and labels for clarity.
         *
         * @param os The output stream to write to.
         * @param data The IntegrationErrorData instance to output.
         * @return Returns the modified output stream.
         */
        friend std::ostream& operator<<(std::ostream& os, const IntegrationErrorData& data)
        {
            os << std::setprecision(16)
                << "Value: " << data.value << "\n"
                << "Abs. Error: " << data.eabs << "\n"
                << "Rel. Error: " << data.erel << "\n"
                << "Iterations: " << data.iterations << "\n";
            return os;
        }
    };
} // namespace nxx::integrate::detail


#endif //INTEGRATIONERROR_HPP

/*
    888b      88  88        88  88b           d88  88888888888  88888888ba   88  8b        d8
    8888b     88  88        88  888b         d888  88           88      "8b  88   Y8,    ,8P
    88 `8b    88  88        88  88`8b       d8'88  88           88      ,8P  88    `8b  d8'
    88  `8b   88  88        88  88 `8b     d8' 88  88aaaaa      88aaaaaa8P'  88      Y88P
    88   `8b  88  88        88  88  `8b   d8'  88  88"""""      88""""88'    88      d88b
    88    `8b 88  88        88  88   `8b d8'   88  88           88    `8b    88    ,8P  Y8,
    88     `8888  Y8a.    .a8P  88    `888'    88  88           88     `8b   88   d8'    `8b
    88      `888   `"Y8888Y"'   88     `8'     88  88888888888  88      `8b  88  8P        Y8

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

#ifndef NUMERIXX_FACTORIZEGJ_HPP
#define NUMERIXX_FACTORIZEGJ_HPP

#include "Matrix.hpp"

namespace nxx::linalg
{

    /**
     * @brief Gauss-Jordan elimination with back-substitution. This function determines the inverse of the
     * coefficient matrix, while also determining the solution vector.
     * @param coefficients The Matrix of coefficients.
     * @param results The vector of results (the b in Ax=b). The contents of the results vector is replaced by the
     * solution vector.
     * @return A std::pair holding the inverse to the coefficient matrix, and the solution vector.
     */
    auto FactorizeGJ(is_matrix auto coefficients, is_matrix auto results)
    {
        // ===== Check that the dimensions of the input matrix and vector match.
        if (coefficients.colCount() != results.rowCount() || results.colCount() != 1)
            throw std::logic_error("Gauss Jordan error: Dimensions of input matrix/vector does not match.");

        // ===== Check that the coefficient matrix has the same extents in each direction.
        if (!coefficients.isSquare()) throw std::logic_error("Gauss Jordan error: The coefficient matrix must be square.");

        // ===== Create the identity matrix, which will be modified and become the inverse of the coefficient matrix.
        auto inverse = std::remove_cvref_t<decltype(coefficients)>::CreateIdentityMatrix(coefficients.colCount());

        // ===== Elimination
        // TODO: Implement pivoting.
        // TODO: Implement scaling.
        for (int i = 0; i < coefficients.rowCount(); ++i) {
            // ===== Select the i'th row of the coefficient and identity matrices.
            // ===== The pivot element is the i'th element of the i'th row of the coefficient matrix
            auto row   = coefficients.row(i);
            auto inv   = inverse.row(i);
            auto pivot = row(0, i);

            // ===== Divide all elements of the i'th row of the coefficient matrix and identity matrices.
            // ===== Do the same to the i'th element of the results vector.
            // ===== This will set the pivot element equal to 1.
            row /= pivot;
            inv /= pivot;
            results(i, 0) /= pivot;

            // ===== For each of the rows below the i'th row, subtract the right amount of the i'th row
            // ===== to set set all elements below the pivot element to zero.
            for (int j = i + 1; j < coefficients.rowCount(); ++j) {
                // ===== Select the j'th row (below the i'th row) of the coefficient and identity matrices.
                // ===== Select the element below the pivot element in the coefficient matrix.
                auto row2 = coefficients.row(j);
                auto inv2 = inverse.row(j);
                auto elem = row2(0, i);

                // ===== Multiply the i'th row by the element below the pivot, and subtract from the j'th row,
                // ===== for both the coefficient and identity matrices.
                // ===== Do the same to the j'th element of the results vector.
                row2          = row2 - row * elem;
                inv2          = inv2 - inv * elem;
                results(j, 0) = results(j, 0) - results(i, 0) * elem;
            }
        }

        // =====Back-substitution
        for (int i = coefficients.rowCount() - 1; i >= 0; --i) {
            // ===== Select the i'th row (starting from the bottom) of the coefficient and identity matrices.
            auto row = coefficients.row(i);
            auto inv = inverse.row(i);

            // ===== For each of the rows above the i'th row, subtract the right amount of the i'th row
            for (int j = i - 1; j >= 0; --j) {
                // ===== Select the j'th row (above the i'th row) of the coefficient and identity matrices.
                // ===== Select the element above the pivot element in the coefficient matrix.
                auto row2 = coefficients.row(j);
                auto inv2 = inverse.row(j);
                auto elem = row2(0, i);

                // ===== Multiply the i'th row by the element above the pivot, and subtract from the j'th row,
                // ===== for both the coefficient and identity matrices.
                // ===== Do the same to the j'th element of the results vector.
                row2          = row2 - row * elem;
                inv2          = inv2 - inv * elem;
                results(j, 0) = results(j, 0) - results(i, 0) * elem;
            }
        }

        return std::make_pair( inverse, results );
    }
}    // namespace numerix::linalg

#endif    // NUMERIXX_FACTORIZEGJ_HPP

//
// Created by Kenneth Balslev on 19/10/2022.
//

#ifndef NUMERIX_FACTORIZEGJ_H
#define NUMERIX_FACTORIZEGJ_H

#include "Matrix.hpp"

namespace numerix::linalg {

    /**
     * @brief
     * @tparam T1
     * @tparam T2
     * @param coefficients
     * @param results
     * @return
     */
    template<typename T1, typename T2>
    requires is_matrix<T1> && is_matrix<T2>
    std::pair<T1, T2> FactorizeGJ(T1 coefficients, T2 results) {

        // ===== Check that the dimensions of the input matrix and vector match.
        if (coefficients.colCount() != results.rowCount() || results.colCount() != 1)
            throw std::logic_error("Gauss Jordan error: Dimensions of input matrix/vector does not match.");

        // ===== Check that the coefficient matrix has the same extents in each direction.
        if (!coefficients.isSquare())
            throw std::logic_error("Gauss Jordan error: The coefficient matrix must be square.");

        // ===== Create the identity matrix, which will be modified and become the inverse of the coefficient matrix.
        T1 inverse = T1::CreateIdentityMatrix(coefficients.colCount());

        // ===== Elimination
        // TODO: Implement pivoting.
        // TODO: Implement scaling.
        for (size_t i = 0; i < coefficients.rowCount(); ++i){

            // ===== Select the i'th row of the coefficient and identity matrices.
            // ===== The pivot element is the i'th element of the i'th row of the coefficient matrix
            auto row = coefficients.row(i);
            auto inv = inverse.row(i);
            auto pivot = row(0,i);

            // ===== Divide all elements of the i'th row of the coefficient matrix and identity matrices.
            // ===== Do the same to the i'th element of the results vector.
            // ===== This will set the pivot element equal to 1.
            row /= pivot;
            inv /= pivot;
            results(i, 0) /= pivot;

            // ===== For each of the rows below the i'th row, subtract the right amount of the i'th row
            // ===== to set set all elements below the pivot element to zero.
            for (size_t j = i + 1; j < coefficients.rowCount(); ++j) {

                // ===== Select the j'th row (below the i'th row) of the coefficient and identity matrices.
                // ===== Select the element below the pivot element in the coefficient matrix.
                auto row2 = coefficients.row(j);
                auto inv2 = inverse.row(j);
                auto elem = row2(0,i);

                // ===== Multiply the i'th row by the element below the pivot, and subtract from the j'th row,
                // ===== for both the coefficient and identity matrices.
                // ===== Do the same to the j'th element of the results vector.
                row2 = row2 - row * elem;
                inv2 = inv2 - inv * elem;
                results(j,0) = results(j,0) - results(i, 0) * elem;
            }
        }

        // =====Back-substitution
        for (int64_t i = coefficients.rowCount() - 1; i >= 0; --i){

            // ===== Select the i'th row (starting from the bottom) of the coefficient and identity matrices.
            auto row = coefficients.row(i);
            auto inv = inverse.row(i);

            // ===== For each of the rows above the i'th row, subtract the right amount of the i'th row
            for (int64_t j = i - 1; j >= 0; --j) {

                // ===== Select the j'th row (above the i'th row) of the coefficient and identity matrices.
                // ===== Select the element above the pivot element in the coefficient matrix.
                auto row2 = coefficients.row(j);
                auto inv2 = inverse.row(j);
                auto elem = row2(0,i);

                // ===== Multiply the i'th row by the element above the pivot, and subtract from the j'th row,
                // ===== for both the coefficient and identity matrices.
                // ===== Do the same to the j'th element of the results vector.
                row2 = row2 - row * elem;
                inv2 = inv2 - inv * elem;
                results(j,0) = results(j,0) - results(i, 0) * elem;
            }
        }

        return {inverse, results};
    }
}


#endif    // NUMERIX_FACTORIZEGJ_H

//
// Created by Kenneth Balslev on 19/10/2022.
//

#ifndef NUMERIX_GAUSSJORDAN_H
#define NUMERIX_GAUSSJORDAN_H

#include "Matrix.hpp"

namespace numerix::linalg {

    template<typename T1, typename T2>
    requires is_matrix<T1> && is_matrix<T2>
    Matrix<std::common_type<typename T1::value_type, typename T2::value_type>> GaussJordan(T1 mat, T2 vec) {

        using result_t = std::common_type<typename T1::value_type, typename T2::value_type>;




        mat.augment(vec);

        for (int i = 0, j = 0; i < mat.rowCount(); ++i, ++j) {

            for (int k = mat.colCount() - 1; k >= j; --k)
                mat[i][k] /= mat[i][j];

            for (int k = i + 1; k < mat.rowCount(); ++k) {
                auto pivot = mat[k][j];
                for (int l = 0; l < mat.colCount(); ++l){
                    mat[k][l] -= mat[i][l] * pivot;
                }
            }
        }

        int k = 2;
        for (int i = mat.rowCount() - 1; i >= 0; --i) {
            auto sol = mat[i][mat.colCount() - 1];
            for (int j = i - 1; j >= 0; --j) {
                mat[j][mat.colCount() - 1] -= mat[j][mat.colCount() - k] * sol;
                mat[j][mat.colCount() - k] = 0.0;
            }
            ++k;
        }

        result_t result(mat.rowCount(), 1);

        for (int i = 0; i < mat.rowCount(); ++i)
            result[i][0] = mat[i][mat.colCount() - 1];

        return result;

    }

}


#endif    // NUMERIX_GAUSSJORDAN_H

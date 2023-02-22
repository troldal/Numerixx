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

#ifndef NUMERIXX_MATRIXFUNCTIONS_HPP
#define NUMERIXX_MATRIXFUNCTIONS_HPP

#include "MatrixCommon.hpp"

namespace nxx::linalg
{

    /**
     * @brief Stream output operator.
     * @param out The output stream.
     * @param mat The matrix to output to stream.
     * @return A reference to the output stream.
     */
    std::ostream& operator<<(std::ostream& out, const is_matrix auto& mat)
    {
        auto printRow = [&](int row) {
            std::cout << "{ ";
            std::cout << std::fixed << std::setw(2);
            for (auto col = 0; col < mat.colCount(); ++col) std::cout << mat(row, col) << " ";
            std::cout << "}";
        };

        for (int i = 0; i < mat.rowCount(); ++i) {
            printRow(i);
            if (i < mat.rowCount() - 1) std::cout << "\n";
        }

        return out;
    }

    /**
     * @brief Matrix addition operator. Adds the elements of two (identically sized) matrices.
     * @param mat1 Matrix 1
     * @param mat2 Matrix 2
     * @return A new matrix with the added numbers.
     */
    inline auto operator+(const is_matrix auto& mat1, const is_matrix auto& mat2)
    {
        using value_t  = std::common_type_t<typename std::remove_cvref_t<decltype(mat1)>::value_type, typename std::remove_cvref_t<decltype(mat2)>::value_type>;
        using result_t = Matrix<value_t>;

        if (mat1.rowCount() != mat2.rowCount()) throw std::invalid_argument("Matrix Addition Error: Row count must be identical.");
        if (mat1.colCount() != mat2.colCount()) throw std::invalid_argument("Matrix Addition Error: Column count must be identical.");

        result_t result { mat1.rowCount(), mat1.colCount() };
        std::transform(mat1.begin(), mat1.end(), mat2.begin(), result.begin(), std::plus());
        return result;
    }

    /**
     * @brief Add a scalar to each element of a matrix.
     * @param mat The Matrix object.
     * @param scalar The scalar to add to each element.
     * @return A new Matrix object with the new values.
     */
    inline auto operator+(const is_matrix auto& mat, is_number auto scalar)
    {
        using value_t  = std::common_type_t<typename std::remove_cvref_t<decltype(mat)>::value_type, std::remove_cvref_t<decltype(scalar)>>;
        using result_t = Matrix<value_t>;

        result_t result { mat.rowCount(), mat.colCount() };
        std::transform(mat.begin(), mat.end(), result.begin(), [&](const auto& v) { return static_cast<result_t>(v) + scalar; });
        return result;
    }

    /**
     * @brief Matrix subtraction operator. Subtracts the elements of two (identically sized) matrices.
     * @param mat1 Matrix 1
     * @param mat2 Matrix 2
     * @return A new matrix with the subtracted numbers.
     */
    inline auto operator-(const is_matrix auto& mat1, const is_matrix auto& mat2)
    {
        using value_t  = std::common_type_t<typename std::remove_cvref_t<decltype(mat1)>::value_type, typename std::remove_cvref_t<decltype(mat2)>::value_type>;
        using result_t = Matrix<value_t>;

        if (mat1.rowCount() != mat2.rowCount()) throw std::invalid_argument("Matrix Subtraction Error: Row count must be identical.");
        if (mat1.colCount() != mat2.colCount()) throw std::invalid_argument("Matrix Subtraction Error: Column count must be identical.");

        result_t result { mat1.rowCount(), mat1.colCount() };
        std::transform(mat1.begin(), mat1.end(), mat2.begin(), result.begin(), [&](const auto& a, const auto& b) {
            return static_cast<value_t>(a) - static_cast<value_t>(b);
        });

        return result;
    }

    /**
     * @brief Subtract a scalar from each element of a matrix.
     * @param mat The Matrix object.
     * @param scalar The scalar to subtract from each element.
     * @return A new Matrix object with the new values.
     */
    inline auto operator-(const is_matrix auto& mat, is_number auto scalar)
    {
        using value_t  = std::common_type_t<typename std::remove_cvref_t<decltype(mat)>::value_type, std::remove_cvref_t<decltype(scalar)>>;
        using result_t = Matrix<value_t>;

        result_t result { mat.rowCount(), mat.colCount() };
        std::transform(mat.begin(), mat.end(), result.begin(), [&](const auto& v) { return static_cast<result_t>(v) - scalar; });

        return result;
    }

    /**
     * @brief Matrix multiplication operator. Multiplies two matrices. The number of columns of the first Matrix must be equal to the
     * number of rows of the second Matrix.
     * @param mat1 Matrix 1
     * @param mat2 Matrix 2
     * @return A new matrix with the result of multiplication.
     */
    inline auto operator*(const is_matrix auto& mat1, const is_matrix auto& mat2)
    {
        using value_t  = std::common_type_t<typename std::remove_cvref_t<decltype(mat1)>::value_type, typename std::remove_cvref_t<decltype(mat2)>::value_type>;
        using result_t = Matrix<value_t>;

        if (mat1.colCount() != mat2.rowCount()) throw std::invalid_argument("Matrix Multiplication Error: Matrix 1 row count must be equal to Matrix 2 column count.");

        result_t result { mat1.rowCount(), mat2.colCount() };
        for (int i = 0; i < result.rowCount(); ++i)
            for (int j = 0; j < result.colCount(); ++j)
                result(i, j) = std::inner_product(mat1.row(i).begin(), mat1.row(i).end(), mat2.col(j).begin(), static_cast<value_t>(0.0));

        return result;
    }

    /**
     * @brief Multiply each element of a Matrix by a scalar.
     * @param mat The Matrix object.
     * @param scalar The scalar to multiply by.
     * @return A new Matrix object with the new values.
     */
    inline auto operator*(const is_matrix auto& mat, is_number auto scalar)
    {
        using value_t  = std::common_type_t<typename std::remove_cvref_t<decltype(mat)>::value_type, std::remove_cvref_t<decltype(scalar)>>;
        using result_t = Matrix<value_t>;

        result_t result { mat.rowCount(), mat.colCount() };
        std::transform(mat.begin(), mat.end(), result.begin(), [&](const auto& v) {
            return static_cast<typename result_t::value_type>(v) * scalar;
        });

        return result;
    }

    /**
     * @brief Multiply each element of a Matrix by a scalar.
     * @param mat The Matrix object.
     * @param scalar The scalar to multiply by.
     * @return A new Matrix object with the new values.
     */
    inline auto operator*(const is_number auto scalar, const is_matrix auto& mat)
    {
        return (mat * scalar);
    }

    /**
     * @brief Divide each element of a Matrix by a scalar.
     * @param mat The Matrix object.
     * @param scalar The scalar to divide by.
     * @return A new Matrix object with the new values.
     */
    inline auto operator/(const is_matrix auto& mat, is_number auto scalar)
    {
        using value_t  = std::common_type_t<typename std::remove_cvref_t<decltype(mat)>::value_type, std::remove_cvref_t<decltype(scalar)>>;
        using result_t = Matrix<value_t>;

        result_t result { mat.rowCount(), mat.colCount() };
        std::transform(mat.begin(), mat.end(), result.begin(), [&](const auto& v) { return static_cast<value_t>(v) / scalar; });

        return result;
    }

    /**
     * @brief Compute the transpose of the input Matrix (or MatrixProxy.
     * @param mat The input Matrix
     * @return A new Matrix object, with the transpose of the input.
     * @note This function does not modify the input. It merely produces a new Matrix object with the transpose.
     */
    inline auto transpose(const is_matrix auto& mat)
    {
        using value_t  = typename std::remove_cvref_t<decltype(mat)>::value_type;
        using result_t = Matrix<value_t>;

        result_t result(mat.colCount(), mat.rowCount());
        for (int i = 0; i < result.rowCount(); ++i) {
            if (result.colCount() == 1)    // TODO: Not sure why this requires special treatment. See if there's a better way.
                result.row(i)(0, 0) = mat.col(i)(0, 0);
            std::copy(mat.col(i).begin(), mat.col(i).end(), result.row(i).begin());
        }

        return result;
    }

}    // namespace numerix::linalg

#endif    // NUMERIXX_MATRIXFUNCTIONS_HPP

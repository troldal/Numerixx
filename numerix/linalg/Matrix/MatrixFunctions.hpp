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

#ifndef NUMERIX_MATRIXFUNCTIONS_HPP
#define NUMERIX_MATRIXFUNCTIONS_HPP

#include "MatrixCommon.hpp"

namespace numerix::linalg
{

    /**
     * @brief
     * @tparam T
     * @param out
     * @param mat
     * @return
     */
    template<typename T>
        requires std::is_same_v<T, Matrix<typename T::value_type>> || std::is_same_v<T, MatrixView<typename T::value_type>> ||
                 std::is_same_v<T, MatrixViewConst<typename T::value_type>>
    std::ostream& operator<<(std::ostream& out, const T& mat)
    {
        auto printRow = [&](int row) {
            std::cout << "{ ";
            std::cout << std::fixed << std::setw(2);
            for (auto col = 0; col < mat.colCount(); ++col) std::cout << mat(row, col) << " ";
            std::cout << "}\n";
        };

        for (int i = 0; i < mat.rowCount(); ++i) printRow(i);

        return out;
    }

    /**
     * @brief
     * @tparam T
     * @param a
     * @param b
     * @return
     */
    template<typename T, typename U>
        requires is_matrix<T> && is_matrix<U>
    inline auto operator+(const T& mat1, const U& mat2)
    {
        using value_t  = std::common_type_t<typename T::value_type, typename U::value_type>;
        using result_t = Matrix<value_t>;

        result_t result { mat1.rowCount(), mat1.colCount() };
        std::transform(mat1.begin(), mat1.end(), mat2.begin(), result.begin(), std::plus());
        return result;
    }

    /**
     * @brief
     * @tparam T
     * @tparam U
     * @param a
     * @param b
     * @return
     */
    template<typename T, typename U>
        requires is_matrix<T> && is_number<U>
    inline auto operator+(const T& mat, U scalar)
    {
        using value_t  = std::common_type_t<typename T::value_type, U>;
        using result_t = Matrix<value_t>;

        result_t result { mat.rowCount(), mat.colCount() };
        std::transform(mat.begin(), mat.end(), result.begin(), [&](const auto& v) { return static_cast<result_t>(v) + scalar; });
        return result;
    }

    /**
     * @brief
     * @tparam T
     * @param a
     * @param b
     * @return
     */
    template<typename T, typename U>
        requires is_matrix<T> && is_matrix<U>
    inline auto operator-(const T& mat1, const U& mat2)
    {
        using value_t  = std::common_type_t<typename T::value_type, typename U::value_type>;
        using result_t = Matrix<value_t>;

        result_t result { mat1.rowCount(), mat1.colCount() };
        std::transform(mat1.begin(), mat1.end(), mat2.begin(), result.begin(), [&](const auto& a, const auto& b) {
            return static_cast<value_t>(a) - static_cast<value_t>(b);
        });

        return result;
    }

    /**
     * @brief
     * @tparam T
     * @tparam U
     * @param a
     * @param b
     * @return
     */
    template<typename T, typename U>
        requires is_matrix<T> && is_number<U>
    inline auto operator-(const T& mat, U scalar)
    {
        using value_t  = std::common_type_t<typename T::value_type, U>;
        using result_t = Matrix<value_t>;

        result_t result { mat.rowCount(), mat.colCount() };
        std::transform(mat.begin(), mat.end(), result.begin(), [&](const auto& v) { return static_cast<result_t>(v) - scalar; });

        return result;
    }

    /**
     * @brief
     * @tparam T
     * @tparam U
     * @param a
     * @param b
     * @return
     */
    template<typename T, typename U>
        requires is_matrix<T> && is_number<U>
    inline auto operator*(const T& mat, U scalar)
    {
        using value_t  = std::common_type_t<typename T::value_type, U>;
        using result_t = Matrix<value_t>;

        result_t result { mat.rowCount(), mat.colCount() };
        std::transform(mat.begin(), mat.end(), result.begin(), [&](const auto& v) {
            return static_cast<typename result_t::value_type>(v) * scalar;
        });

        return result;
    }

    /**
     * @brief
     * @tparam T
     * @tparam U
     * @param a
     * @param b
     * @return
     */
    template<typename T, typename U>
        requires is_matrix<U> && is_number<T>
    inline auto operator*(const T scalar, U& mat)
    {
        return (mat * scalar);
    }

    /**
     * @brief
     * @tparam T
     * @tparam U
     * @param a
     * @param b
     * @return
     */
    template<typename T, typename U>
        requires is_matrix<U> && is_matrix<T>
    inline auto operator*(const T& mat1, const U& mat2)
    {
        using value_t  = std::common_type_t<typename T::value_type, typename U::value_type>;
        using result_t = Matrix<value_t>;

        result_t result { mat1.rowCount(), mat2.colCount() };
        for (int i = 0; i < result.rowCount(); ++i)
            for (int j = 0; j < result.colCount(); ++j)
                result(i, j) = std::inner_product(mat1.row(i).begin(), mat1.row(i).end(), mat2.col(j).begin(), static_cast<value_t>(0.0));

        return result;
    }

    /**
     * @brief
     * @tparam T
     * @tparam U
     * @param a
     * @param b
     * @return
     */
    template<typename T, typename U>
        requires is_matrix<T> && is_number<U>
    inline auto operator/(const T& mat, U scalar)
    {
        using value_t  = std::common_type_t<typename T::value_type, U>;
        using result_t = Matrix<value_t>;

        result_t result { mat.rowCount(), mat.colCount() };
        std::transform(mat.begin(), mat.end(), result.begin(), [&](const auto& v) { return static_cast<value_t>(v) / scalar; });

        return result;
    }

    /**
     * @brief Compute the transpose of the input Matrix (or MatrixProxy.
     * @tparam T The type of the input Matrix. Must be of Matrix or MatrixProxy types.
     * @param mat The input Matrix
     * @return A new Matrix object, with the transpose of the input.
     * @note This function does not modify the input. It merely produces a new Matrix object with the transpose.
     */
    template<typename T>
        requires is_matrix<T>
    inline Matrix<typename T::value_type> transpose(const T& mat)
    {
        using value_t  = typename T::value_type;
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

// template<typename T>
//     requires std::same_as<T, numerix::linalg::MatrixElements<typename numerix::linalg::impl::MatrixTraits<T>::value_type>> ||
//              std::same_as<T, numerix::linalg::MatrixElementsConst<typename numerix::linalg::impl::MatrixTraits<T>::value_type>>
// inline void swap(T lhs, T rhs)
//{
//     lhs.swap(rhs);
// }

#endif    // NUMERIX_MATRIXFUNCTIONS_HPP

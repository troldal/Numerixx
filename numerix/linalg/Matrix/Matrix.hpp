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

#ifndef NUMERIX_MATRIX_IMPL_HPP
#define NUMERIX_MATRIX_IMPL_HPP

#include "MatrixBase.hpp"
#include "MatrixCommon.hpp"

namespace numerix::linalg
{

    /**
     * @brief The Matrix class is the main abstraction for matrices. It derives from the MatrixBase class using CRTP.
     * The Matrix class owns the underlying data, unlike the MatrixView class which simply is a view into a Matrix object.
     * @tparam T The type of the underlying matrix elements. Default is double, but can be any kind of floating point or integer number.
     */
    template<typename T = double>
        requires is_number<T>
    class Matrix : public impl::MatrixBase<Matrix<T>>
    {
        /**
         * Friend declarations. Necessary for the base class to access elements of Matrix using CRTP.
         */
        friend impl::MatrixBase<Matrix<T>>;
        friend MatrixView<T>;
        friend MatrixViewConst<T>;

        /**
         * Private alias declatations.
         */
        using parent = impl::MatrixBase<Matrix<T>>;

    private:
        std::vector<T> m_data;     /**< The underlying array of matrix elements. */
        Slice          m_rowSlice; /**< The Slice describing the rows. Required to provide a common interface with MatrixView. */
        Slice          m_colSlice; /**< The Slice describing the columns. Required to provide a common interface with MatrixView. */

        /**
         * @brief
         * @return
         */
        auto gslice() const
        {
            auto start = m_rowSlice.start() * extents().first + m_colSlice.start();
            return impl::GSlice(start, { m_rowSlice.length(), m_colSlice.length() }, { m_rowSlice.stride(), m_colSlice.stride() });
        }

        /**
         * @brief Get the extents in each of the dimensions (rows and columns)
         * @return A std::pair with the row and column extents.
         */
        auto extents() const { return std::make_pair(parent::rowCount(), parent::colCount()); }

    public:
        /**
         * Public alias declatations. To be consistant with standard library containers, and to provide access to
         * non-standard assignment operator.
         */
        using value_type = T;
        using parent::operator=;

        /**
         * @brief Constructor taking the number of rows and columns.
         * @param rows Number of rows.
         * @param cols Number of cols.
         */
        Matrix(int rows, int cols) : m_data(rows * cols), m_rowSlice(0, rows, cols), m_colSlice(0, cols, 1)
        {
            if (rows <= 0) throw std::invalid_argument("Invalid Matrix Extents: A Matrix object must have at least one row.");
            if (cols <= 0) throw std::invalid_argument("Invalid Matrix Extents: A Matrix object must have at least one column.");
        }

        /**
         * @brief Copy constructor.
         * @param other Object to be copied.
         */
        Matrix(const Matrix& other) = default;

        /**
         * @brief Move constructor.
         * @param other Object to be moved.
         */
        Matrix(Matrix&& other) noexcept = default;

        /**
         * @brief Destructor.
         */
        ~Matrix() = default;

        /**
         * @brief Copy assignment operator.
         * @param other Object to be copied.
         * @return The copied-to object.
         */
        Matrix& operator=(const Matrix& other) = default;

        /**
         * @brief Move assignment operator.
         * @param other Object to be moved.
         * @return The moved-to object.
         */
        Matrix& operator=(Matrix&& other) noexcept = default;

        /**
         * @brief Access the raw array of Matrix elements.
         * @return A pointer to the first element.
         */
        T* data() { return m_data.data(); }

        /**
         * @brief Access the raw array of Matrix elements.
         * @return A const pointer to the first element.
         */
        const T* data() const { return m_data.data(); }
    };

}    // namespace numerix::linalg

#endif    // NUMERIX_MATRIX_HPP

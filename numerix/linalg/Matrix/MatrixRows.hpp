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

#ifndef NUMERIX_MATRIXROWS_HPP
#define NUMERIX_MATRIXROWS_HPP

#include "MatrixCommon.h"

namespace numerix::linalg
{

    /**
     * @brief Collection of the rows of a Matrix or MatrixView object. This can be used to iterate through all rows of a Matrix.
     * @tparam MATRIX The type of the Matrix object
     * @tparam IsConst Flag indicating if the object should be treated as const (non-mutable)
     */
    template<typename MATRIX, bool IsConst>
        requires std::same_as<MATRIX, MatrixView<typename MATRIX::value_type>> ||
                 std::same_as<MATRIX, MatrixViewConst<typename MATRIX::value_type>>
    class MatrixRowsConcept
    {
        /*
         * Friend declarations. Required to allow use of private constructor.
         */
        friend impl::MatrixBase<Matrix<typename MATRIX::value_type>>;
        friend impl::MatrixBase<MatrixView<typename MATRIX::value_type>>;
        friend impl::MatrixBase<MatrixViewConst<typename MATRIX::value_type>>;

        /**
         * Private alias declatations.
         */
        using matrix_t = std::conditional_t<IsConst,
                                            MatrixViewConst<typename std::remove_cvref_t<MATRIX>::value_type>,
                                            MatrixView<typename std::remove_cvref_t<MATRIX>::value_type>>;

    private:
        matrix_t m_matrix; /**< A MatrixView object containing the columns. */

        /**
         * @brief Constructor taking a MatrixView as an argument.
         * @param data The MatrixView representing the collection of rows.
         * @note Constructor is private to avoid manual creation by client.
         */
        explicit MatrixRowsConcept(matrix_t data) : m_matrix(data) {}

    public:
        /**
         * Public alias declatations. To be consistant with standard library containers.
         */
        using matrix_type = matrix_t;
        using value_type  = typename matrix_t::value_type;

        /**
         * @brief Copy constructor.
         * @param other Object to be copied.
         */
        MatrixRowsConcept(const MatrixRowsConcept& other) = default;

        /**
         * @brief Move constructor.
         * @param other Object to be moved.
         */
        MatrixRowsConcept(MatrixRowsConcept&& other) noexcept = default;

        /**
         * @brief Destructor.
         */
        ~MatrixRowsConcept() = default;

        /**
         * @brief Copy assignment operator.
         * @param other Object to be copied.
         * @return The copied-to object.
         */
        MatrixRowsConcept& operator=(const MatrixRowsConcept& other) = default;

        /**
         * @brief Move assignment operator.
         * @param other Object to be moved.
         * @return The moved-to object.
         */
        MatrixRowsConcept& operator=(MatrixRowsConcept&& other) noexcept = default;

        /**
         * @brief Function call operator, taking the index of the row to return, as an argument.
         * @param index The index of the row to return.
         * @return The row at the given index.
         * @note This version is disabled for non-mutable objects.
         */
        auto operator()(size_t index)
            requires(!std::same_as<MATRIX, MatrixViewConst<typename MATRIX::value_type>>)
        {
            return m_matrix.row(index);
        }

        /**
         * @brief Function call operator, taking the index of the row to return, as an argument.
         * @param index The index of the row to return.
         * @return The row at the given index.
         */
        auto operator()(size_t index) const { return m_matrix.row(index); }

        /**
         * @brief Subscript operator, taking the index of the row to return, as an argument.
         * @param index The index of the row to return.
         * @return The row at the given index.
         * @note This version is disabled for non-mutable objects.
         */
        auto operator[](size_t index)
            requires(!std::same_as<MATRIX, MatrixViewConst<typename MATRIX::value_type>>)
        {
            return m_matrix.row(index);
        }

        /**
         * @brief Subscript operator, taking the index of the row to return, as an argument.
         * @param index The index of the row to return.
         * @return The row at the given index.
         */
        auto operator[](size_t index) const { return m_matrix.row(index); }

        /**
         * @brief Get the size, i.e. the count of rows.
         * @return The row count.
         */
        size_t size() const { return m_matrix.rowCount(); }

        /**
         * @brief Get an iterator to the first row in the collection.
         * @return An iterator to the first row.
         * @note This version is disabled for non-mutable objects.
         */
        auto begin()
            requires(!std::same_as<MATRIX, MatrixViewConst<typename MATRIX::value_type>>)
        {
            return MatrixRowIter<typename std::remove_reference_t<decltype(*this)>>(*this, 0);
        }

        /**
         * @brief Get a const iterator to the first row in the collection.
         * @return A const iterator to the first row.
         */
        auto begin() const { return MatrixRowIterConst<typename std::remove_cvref_t<decltype(*this)>>(*this, 0); }

        /**
         * @brief Get a const iterator to the first row in the collection.
         * @return A const iterator to the first row.
         */
        auto cbegin() const { return MatrixRowIterConst<typename std::remove_cvref_t<decltype(*this)>>(*this, 0); }

        /**
         * @brief Get an iterator to one past the last row in the collection.
         * @return An iterator to one past the last row.
         * @note This version is disabled for non-mutable objects.
         */
        auto end()
            requires(!std::same_as<MATRIX, MatrixViewConst<typename MATRIX::value_type>>)
        {
            return MatrixRowIter<typename std::remove_reference_t<decltype(*this)>>(*this, size());
        }

        /**
         * @brief Get a const iterator to one past the last row in the collection.
         * @return A const iterator to one past the last row.
         */
        auto end() const { return MatrixRowIterConst<typename std::remove_cvref_t<decltype(*this)>>(*this, size()); }

        /**
         * @brief Get a const iterator to one past the last row in the collection.
         * @return A const iterator to one past the last row.
         */
        auto cend() const { return MatrixRowIterConst<typename std::remove_cvref_t<decltype(*this)>>(*this, size()); }

        /**
         * @brief Get the first row in the collection.
         * @return The first row.
         * @note This version is disabled for non-mutable objects.
         */
        auto front()
            requires(!std::same_as<MATRIX, MatrixViewConst<typename MATRIX::value_type>>)
        {
            return (*this)[0];
        }

        /**
         * @brief Get the first row in the collection.
         * @return The first row.
         */
        auto front() const { return (*this)[0]; }

        /**
         * @brief Get the last row in the collection.
         * @return The last row.
         * @note This version is disabled for non-mutable objects.
         */
        auto back()
            requires(!std::same_as<MATRIX, MatrixViewConst<typename MATRIX::value_type>>)
        {
            return (*this)[m_matrix.rowCount() - 1];
        }

        /**
         * @brief Get the last row in the collection.
         * @return The last row.
         */
        auto back() const { return (*this)[m_matrix.rowCount() - 1]; }
    };

}    // namespace numerix::linalg

#endif    // NUMERIX_MATRIXROWS_HPP

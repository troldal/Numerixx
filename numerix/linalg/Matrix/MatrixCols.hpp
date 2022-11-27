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

#ifndef NUMERIX_MATRIXCOLS_HPP
#define NUMERIX_MATRIXCOLS_HPP

#include "MatrixCommon.h"

namespace numerix::linalg
{

    /**
     * @brief
     * @tparam T
     * @tparam IsConst
     */
    template<typename T, bool IsConst>
        requires std::same_as<T, MatrixView<typename T::value_type>> || std::same_as<T, MatrixViewConst<typename T::value_type>>
    class MatrixColsConcept
    {
        /*
         *
         */
        using matrix_t = std::conditional_t<IsConst,
                                            MatrixViewConst<typename std::remove_reference_t<T>::value_type>,
                                            MatrixView<typename std::remove_reference_t<T>::value_type>>;

        matrix_t m_matrix; /**< A MatrixView object containing the columns. */

    public:
        /*
         *
         */
        using value_type = matrix_t;

        /**
         * @brief
         * @param data
         */
        explicit MatrixColsConcept(matrix_t data) : m_matrix(data) {}

        /**
         * @brief
         * @param other
         */
        MatrixColsConcept(const MatrixColsConcept& other) = default;

        /**
         * @brief
         * @param other
         */
        MatrixColsConcept(MatrixColsConcept&& other) noexcept = default;

        /**
         * @brief
         * @param other
         * @return
         */
        MatrixColsConcept& operator=(const MatrixColsConcept& other) = default;

        /**
         * @brief
         * @param other
         * @return
         */
        MatrixColsConcept& operator=(MatrixColsConcept&& other) noexcept = default;

        /**
         * @brief
         * @param index
         * @return
         */
        auto operator()(size_t index) { return m_matrix.col(index); }

        /**
         * @brief
         * @param index
         * @return
         */
        auto operator()(size_t index) const { return m_matrix.col(index); }

        /**
         * @brief
         * @param index
         * @return
         */
        auto operator[](size_t index) { return m_matrix.col(index); }

        /**
         * @brief
         * @param index
         * @return
         */
        auto operator[](size_t index) const { return m_matrix.col(index); }

        /**
         * @brief
         * @return
         */
        size_t size() const { return m_matrix.colCount(); }

        /**
         * @brief
         * @return
         */
        auto begin()
            requires(!std::same_as<T, MatrixViewConst<typename T::value_type>>)
        {
            return MatrixColIter<typename std::remove_reference_t<decltype(*this)>>(*this, 0);
        }

        /**
         * @brief
         * @return
         */
        auto begin() const { return MatrixColIterConst<typename std::remove_cvref_t<decltype(*this)>>(*this, 0); }

        /**
         * @brief
         * @return
         */
        auto cbegin() const { return MatrixColIterConst<typename std::remove_cvref_t<decltype(*this)>>(*this, 0); }

        /**
         * @brief
         * @return
         */
        auto end()
            requires(!std::same_as<T, MatrixViewConst<typename T::value_type>>)
        {
            return MatrixColIter<typename std::remove_reference_t<decltype(*this)>>(*this, size());
        }

        /**
         * @brief
         * @return
         */
        auto end() const { return MatrixColIterConst<typename std::remove_cvref_t<decltype(*this)>>(*this, size()); }

        /**
         * @brief
         * @return
         */
        auto cend() const { return MatrixColIterConst<typename std::remove_cvref_t<decltype(*this)>>(*this, size()); }

        auto front() { return (*this)[0]; }

        auto front() const { return (*this)[0]; }

        auto back() { return (*this)[m_matrix.colCount() - 1]; }

        auto back() const { return (*this)[m_matrix.colCount() - 1]; }
    };

}

#endif    // NUMERIX_MATRIXCOLS_HPP

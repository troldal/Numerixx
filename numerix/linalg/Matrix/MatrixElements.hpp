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

#ifndef NUMERIX_MATRIXELEMENTS_HPP
#define NUMERIX_MATRIXELEMENTS_HPP

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
    class MatrixElementsConcept
    {
    private:
        /*
         *
         */
        using matrix_t = std::conditional_t<IsConst,
                                            MatrixViewConst<typename std::remove_reference_t<T>::value_type>,
                                            MatrixView<typename std::remove_reference_t<T>::value_type>>;

        matrix_t m_matrix; /**< A MatrixView object containing the columns. */

        friend class impl::MatrixBase<Matrix<typename matrix_t::value_type>>;
        friend class impl::MatrixBase<MatrixViewConcept<typename matrix_t::value_type, false>>;
        friend class impl::MatrixBase<MatrixViewConcept<typename matrix_t::value_type, true>>;
        friend void ::swap(MatrixElementsConcept lhs, MatrixElementsConcept rhs);

        /**
         * @brief
         * @param data
         */
        explicit MatrixElementsConcept(matrix_t data) : m_matrix(data) {}

        /**
         * @brief
         * @param other
         */
        MatrixElementsConcept(const MatrixElementsConcept& other) = default;

        /**
         * @brief
         * @param other
         */
        MatrixElementsConcept(MatrixElementsConcept&& other) noexcept = default;

    public:
        /*
         *
         */
        using value_type = matrix_t;

        /**
         * @brief
         * @param other
         * @return
         */
        MatrixElementsConcept& operator=(const MatrixElementsConcept& other)
        {
            if (other.m_matrix.rowCount() != m_matrix.rowCount() || other.m_matrix.colCount() != m_matrix.colCount())
                throw std::out_of_range("Matrices have different sizes.");

            std::copy(other.begin(), other.end(), begin());

            return *this;
        }

        /**
         * @brief
         * @param other
         * @return
         */
        MatrixElementsConcept& operator=(MatrixElementsConcept&& other) noexcept
        {
            if (other.m_matrix.rowCount() != m_matrix.rowCount() || other.m_matrix.colCount() != m_matrix.colCount())
                throw std::out_of_range("Matrices have different sizes.");

            std::copy(other.begin(), other.end(), begin());

            return *this;
        }

        /**
         * @brief
         * @tparam C
         * @param other
         * @return
         */
        template<typename C>
            requires std::same_as<typename matrix_t::value_type, typename C::value_type>
        MatrixElementsConcept& operator=(const C& other)
        {
            if (other.size() != m_matrix.size()) throw std::out_of_range("Matrices have different sizes.");

            std::copy(other.begin(), other.end(), begin());

            return *this;
        }

        /**
         * @brief
         * @tparam C
         * @return
         */
        template<typename C>
            requires (!is_matrix<C>) && std::same_as<typename matrix_t::value_type, typename C::value_type>
        operator C() const
        {
            C result(m_matrix.size());
            std::copy(begin(), end(), result.begin());
            return result;
        }

        /**
         * @brief
         * @return
         */
        auto begin()
            requires(!std::same_as<T, MatrixViewConst<typename T::value_type>>) && (!std::is_const_v<T>)
        {
            return m_matrix.begin();
        }

        /**
         * @brief
         * @return
         */
        auto begin() const { return m_matrix.begin(); }

        /**
         * @brief
         * @return
         */
        auto cbegin() const { return m_matrix.begin(); }

        /**
         * @brief
         * @return
         */
        auto end()
            requires(!std::same_as<T, MatrixViewConst<typename T::value_type>>) && (!std::is_const_v<T>)
        {
            return m_matrix.end();
        }

        /**
         * @brief
         * @return
         */
        auto end() const { return m_matrix.end(); }

        /**
         * @brief
         * @return
         */
        auto cend() const { return m_matrix.end(); }

        /**
         * @brief
         * @param other
         */
        void swap(MatrixElementsConcept other)
        {
            std::vector<typename matrix_t::value_type> temp = *this;
            *this                                           = other;
            other                                           = temp;
        }
    };

}    // namespace numerix::linalg


#endif    // NUMERIX_MATRIXELEMENTS_HPP

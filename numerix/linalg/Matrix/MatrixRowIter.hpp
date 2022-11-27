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

#ifndef NUMERIX_MATRIXROWITER_HPP
#define NUMERIX_MATRIXROWITER_HPP

#include "MatrixCommon.h"

namespace numerix::linalg
{

    /**
     * @brief
     * @tparam T
     * @tparam IsConst
     */
    template<typename T, bool IsConst>
        requires std::same_as<T, MatrixRows<typename T::value_type>> || std::same_as<T, MatrixRowsConst<typename T::value_type>>
    class MatrixRowIterConcept
    {
        /*
         *
         */
        using rows_t = std::conditional_t<IsConst, MatrixRowsConst<typename T::value_type>, MatrixRows<typename T::value_type>>;
        using row_t  = std::
            conditional_t<IsConst, MatrixViewConst<typename T::value_type::value_type>, MatrixView<typename T::value_type::value_type>>;

        rows_t                 m_rows;             /**< A pointer to the matrix element array. */
        size_t                 m_current;          /**< The current index. */
        std::unique_ptr<row_t> m_currow = nullptr; /**< */

    public:
        /*
         * Alias declarations.
         */
        using iterator_category = std::forward_iterator_tag;
        using value_type        = row_t;
        using difference_type   = size_t;
        using pointer           = row_t*;
        using reference         = row_t&;

        /**
         * @brief
         * @param data
         * @param slice
         * @param pos
         */
        MatrixRowIterConcept(rows_t data, size_t pos = 0) : m_rows(data), m_current(pos) {}

        /**
         * @brief
         * @return
         */
        MatrixRowIterConcept end() const { return { m_rows, m_rows.size() }; }

        /**
         * @brief
         * @return
         */
        MatrixRowIterConcept& operator++()
        {
            ++m_current;
            return *this;
        }

        /**
         * @brief
         * @return
         */
        MatrixRowIterConcept operator++(int)
        {
            MatrixColIterConcept slice = *this;
            ++m_current;
            return slice;
        }

        /**
         * @brief
         * @return
         */
        row_t& operator*()
        {
            m_currow = std::make_unique<row_t>(m_rows(m_current));
            return *m_currow;
        }

        /**
         * @brief
         * @return
         */
        row_t* operator->()
        {
            m_currow = std::make_unique<row_t>(m_rows(m_current));
            return m_currow.get();
        }

        /**
         * @brief
         * @param other
         * @return
         */
        bool operator==(const MatrixRowIterConcept& other) const { return m_current == other.m_current; }

        /**
         * @brief
         * @param other
         * @return
         */
        bool operator!=(const MatrixRowIterConcept& other) const { return !(*this == other); }

        /**
         * @brief
         * @param other
         * @return
         */
        bool operator<(const MatrixRowIterConcept& other) const { return m_current < other.m_current; }
    };
}


#endif    // NUMERIX_MATRIXROWITER_HPP

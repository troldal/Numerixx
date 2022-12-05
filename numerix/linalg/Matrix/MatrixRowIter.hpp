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
     * @brief Row iterator. Used to iterate over the rows of a Matrix or MatrixView object.
     * @tparam ROWCOLL The type of collection to iterate over. This can be a MatrixRows or MatrixRowsConst object.
     * @tparam IsConst Flag indicating if objects of the class should be non-mutable.
     */
    template<typename ROWCOLL, bool IsConst>
        requires std::same_as<ROWCOLL, MatrixRows<typename ROWCOLL::matrix_type>> ||
                 std::same_as<ROWCOLL, MatrixRowsConst<typename ROWCOLL::matrix_type>>
    class MatrixRowIterConcept
    {
        /*
         * Friend declarations. Required to allow use of private constructor.
         */
        friend MatrixRows<typename ROWCOLL::matrix_type>;
        friend MatrixRowsConst<typename ROWCOLL::matrix_type>;

        /**
         * Alias for the row container (MatrixRows or MatrixRowsConst).
         */
        using rows_t =
            std::conditional_t<IsConst, MatrixRowsConst<typename ROWCOLL::matrix_type>, MatrixRows<typename ROWCOLL::matrix_type>>;

        /**
         * Alias for the row type (MatrixView or MatrixViewConst).
         */
        using row_t = std::conditional_t<IsConst,
                                         MatrixViewConst<typename ROWCOLL::matrix_type::value_type>,
                                         MatrixView<typename ROWCOLL::matrix_type::value_type>>;

    private:
        rows_t                 m_rows;             /**< A pointer to the row container. */
        size_t                 m_current;          /**< The current index. */
        std::unique_ptr<row_t> m_currow = nullptr; /**< A unique_ptr to the current row (can be nullptr, e.g. for the end()-iterator. */

        /**
         * @brief Constructor taking the row container and starting position (default 0) as arguments.
         * @param data The row container to iterate over (can be MatrixRows or MatrixRowsConst).
         * @param pos The starting position (default 0).
         * @note The constructor is private to avoid direct usage by clients.
         */
        MatrixRowIterConcept(rows_t data, size_t pos = 0) : m_rows(data), m_current(pos) {}

    public:
        /*
         * Alias declarations. Required to conform to the interface of standard iterators.
         */
        using iterator_category = std::forward_iterator_tag;
        using value_type        = row_t;
        using difference_type   = size_t;
        using pointer           = row_t*;
        using reference         = row_t&;
        using matrix_type       = row_t;

        /**
         * @brief Increment operator. Moves the iterator one step forward.
         * @return A reference to the incremented iterator.
         */
        MatrixRowIterConcept& operator++()
        {
            ++m_current;
            return *this;
        }

        /**
         * @brief Post-increment operator. Moves the iterator one step forward.
         * @return The iterator prior to incrementing it.
         */
        MatrixRowIterConcept operator++(int)
        {
            MatrixColIterConcept slice = *this;
            ++m_current;
            return slice;
        }

        /**
         * @brief Indirection operator.
         * @return A reference to the pointed-to object.
         */
        row_t& operator*()
        {
            m_currow = std::make_unique<row_t>(m_rows(m_current));
            return *m_currow;
        }

        /**
         * @brief Arrow operator.
         * @return A pointer to the pointed-to object.
         */
        row_t* operator->()
        {
            m_currow = std::make_unique<row_t>(m_rows(m_current));
            return m_currow.get();
        }

        /**
         * @brief Equality operator. Check if the argument is equal to the iterator (i.e. the position is the same).
         * @param other The iterator to compare to.
         * @return If they are equal, true; otherwise false.
         */
        bool operator==(const MatrixRowIterConcept& other) const { return m_current == other.m_current; }

        /**
         * @brief Inequality operator. Check if the argument is not equal to the iterator (i.e. the position is not the same).
         * @param other The iterator to compare to.
         * @return If they are not equal, true; otherwise false.
         */
        bool operator!=(const MatrixRowIterConcept& other) const { return !(*this == other); }

        /**
         * @brief Less-than operator. Check if the argument is less than the iterator (i.e. the position is lower than for the iterator).
         * @param other The iterator to compare to.
         * @return If the argument is less than the current, true; otherwise false.
         */
        bool operator<(const MatrixRowIterConcept& other) const { return m_current < other.m_current; }
    };
}    // namespace numerix::linalg

#endif    // NUMERIX_MATRIXROWITER_HPP

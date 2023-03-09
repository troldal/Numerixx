/*
    888b      88  88        88  88b           d88  88888888888  88888888ba   88  8b        d8  8b        d8
    8888b     88  88        88  888b         d888  88           88      "8b  88   Y8,    ,8P    Y8,    ,8P
    88 `8b    88  88        88  88`8b       d8'88  88           88      ,8P  88    `8b  d8'      `8b  d8'
    88  `8b   88  88        88  88 `8b     d8' 88  88aaaaa      88aaaaaa8P'  88      Y88P          Y88P
    88   `8b  88  88        88  88  `8b   d8'  88  88"""""      88""""88'    88      d88b          d88b
    88    `8b 88  88        88  88   `8b d8'   88  88           88    `8b    88    ,8P  Y8,      ,8P  Y8,
    88     `8888  Y8a.    .a8P  88    `888'    88  88           88     `8b   88   d8'    `8b    d8'    `8b
    88      `888   `"Y8888Y"'   88     `8'     88  88888888888  88      `8b  88  8P        Y8  8P        Y8

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

#ifndef NUMERIXX_MATRIXELEMITER_HPP
#define NUMERIXX_MATRIXELEMITER_HPP

#include <utility>

#include "MatrixCommon.hpp"
#include "MatrixSlice.hpp"

namespace nxx::linalg
{
    /**
     * @brief
     * @tparam T
     * @tparam IsConst
     */
    template<typename T, bool IsConst>
    class MatrixElemIterConcept
    {
        /*
         * Friend declarations. Required to allow use of private constructor.
         */
        friend impl::MatrixBase<Matrix<T>>;
        friend impl::MatrixBase<MatrixView<T>>;
        friend impl::MatrixBase<MatrixViewConst<T>>;

        /**
         * Private alias declatations.
         */
        using iter_t = std::conditional_t<IsConst, const T, T>;

    private:
        iter_t* m_data;    /**< A pointer to the matrix element array. */
        impl::GSlice  m_slice;   /**< The generalized slice for the data array. */
        int  m_current; /**< The current index. */

        /**
         * @brief Constructor taking a raw array of elements, a gslice, and the starting position (default 0) as arguments.
         * @param data The raw array of elements (pointer to the first element).
         * @param slice The generalized slice, describing the subset of elements to iterate over.
         * @param pos The starting position (default 0)
         * @note The constructor is private to prevent direct usage by clients.
         */
        MatrixElemIterConcept(iter_t* data, impl::GSlice slice, int pos = 0) : m_data(data), m_slice(std::move(slice)), m_current(pos) {}

        /**
         * @brief Get a copy of the end iterator
         * @return The end iterator (one past the last element).
         */
        MatrixElemIterConcept end() const { return { m_data, m_slice, m_slice.size() }; }

    public:
        /*
         * Alias declarations. Required to conform to the interface of standard iterators.
         */
        using iterator_category = std::forward_iterator_tag;
        using value_type        = iter_t;
        using difference_type   = size_t;
        using pointer           = iter_t*;
        using reference         = iter_t&;
        
        /**
         * @brief Increment operator. Moves the iterator one step forward.
         * @return A reference to the incremented iterator.
         */
        MatrixElemIterConcept& operator++()
        {
            ++m_current;
            return *this;
        }

        /**
         * @brief Post-increment operator. Moves the iterator one step forward.
         * @return The iterator prior to incrementing it.
         */
        MatrixElemIterConcept operator++(int) // NOLINT
        {
            MatrixElemIterConcept slice = *this;
            ++m_current;
            return slice;
        }

        /**
         * @brief Indirection operator.
         * @return A reference to the pointed-to object.
         */
        iter_t& operator*() { return m_data[m_slice(m_current)]; }

        /**
         * @brief Arrow operator.
         * @return A pointer to the pointed-to object.
         */
        iter_t* operator->() { return &m_data[m_slice(m_current)]; }

        /**
         * @brief Equality operator. Check if the argument is equal to the iterator (i.e. the position is the same).
         * @param other The iterator to compare to.
         * @return If they are equal, true; otherwise false.
         */
        bool operator==(const MatrixElemIterConcept& other) const { return m_current == other.m_current; }

        /**
         * @brief Inequality operator. Check if the argument is not equal to the iterator (i.e. the position is the not same).
         * @param other The iterator to compare to.
         * @return If they are not equal, true; otherwise false.
         */
        bool operator!=(const MatrixElemIterConcept& other) const { return !(*this == other); } // NOLINT

        /**
         * @brief Less-than operator. Check if the argument is less than the iterator (i.e. the position is lower than for the iterator).
         * @param other The iterator to compare to.
         * @return If the argument is less than the current, true; otherwise false.
         */
        bool operator<(const MatrixElemIterConcept& other) const { return m_current < other.m_current; }
    };

}    // namespace numerix::linalg

#endif    // NUMERIXX_MATRIXELEMITER_HPP

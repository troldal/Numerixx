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

#ifndef NUMERIX_MATRIXELEMENTITER_HPP
#define NUMERIX_MATRIXELEMENTITER_HPP

#include "MatrixCommon.h"
#include "MatrixSlice.hpp"


namespace numerix::linalg
{
    /**
     * @brief
     * @tparam T
     * @tparam IsConst
     */
    template<typename T, bool IsConst>
    class MatrixElementIterConcept
    {
        using iter_t = std::conditional_t<IsConst, const T, T>;

        iter_t* m_data;    /**< A pointer to the matrix element array. */
        GSlice  m_slice;   /**< The generalized slice for the data array. */
        size_t  m_current; /**< The current index. */

    public:
        /*
         * Alias declarations.
         */
        using iterator_category = std::forward_iterator_tag;
        using value_type        = iter_t;
        using difference_type   = size_t;
        using pointer           = iter_t*;
        using reference         = iter_t&;

        /**
         * @brief
         * @param data
         * @param slice
         * @param pos
         */
        MatrixElementIterConcept(iter_t* data, GSlice slice, size_t pos = 0) : m_data(data), m_slice(slice), m_current(pos) {}

        /**
         * @brief
         * @return
         */
        MatrixElementIterConcept end() const { return { m_data, m_slice, m_slice.size() }; }

        /**
         * @brief
         * @return
         */
        MatrixElementIterConcept& operator++()
        {
            ++m_current;
            return *this;
        }

        /**
         * @brief
         * @return
         */
        MatrixElementIterConcept operator++(int)
        {
            MatrixElementIterConcept slice = *this;
            ++m_current;
            return slice;
        }

        /**
         * @brief
         * @return
         */
        iter_t& operator*() { return m_data[m_slice(m_current)]; }

        /**
         * @brief
         * @return
         */
        iter_t* operator->() { return &m_data[m_slice(m_current)]; }

        /**
         * @brief
         * @param other
         * @return
         */
        bool operator==(const MatrixElementIterConcept& other) const { return m_current == other.m_current; }

        /**
         * @brief
         * @param other
         * @return
         */
        bool operator!=(const MatrixElementIterConcept& other) const { return !(*this == other); }

        /**
         * @brief
         * @param other
         * @return
         */
        bool operator<(const MatrixElementIterConcept& other) const { return m_current < other.m_current; }
    };

}

#endif    // NUMERIX_MATRIXELEMENTITER_HPP

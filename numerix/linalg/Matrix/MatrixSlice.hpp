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

#ifndef NUMERIX_MATRIXSLICE_HPP
#define NUMERIX_MATRIXSLICE_HPP

#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <vector>

namespace numerix::linalg
{

    /**
     * @brief
     */
    class Slice
    {
    public:
        /**
         * @brief
         */
        Slice() : m_start(0), m_length(0), m_stride(1) {}

        /**
         * @brief
         * @param start
         */
        explicit Slice(size_t start) : m_start(start), m_length(0), m_stride(1) {}

        /**
         * @brief
         * @param start
         * @param length
         * @param stride
         */
        Slice(size_t start, size_t length, size_t stride = 1) : m_start(start), m_length(length), m_stride(stride) {}

        /**
         * @brief
         * @param index
         * @return
         */
        size_t operator()(size_t index) const
        {
            if (index > m_length) throw std::out_of_range("Index out of bounds.");
            return m_start + index * m_stride;
        }

        /**
         * @brief
         * @return
         */
        size_t start() const { return m_start; }

        /**
         * @brief
         * @return
         */
        size_t length() const { return m_length; }

        /**
         * @brief
         * @return
         */
        size_t stride() const { return m_stride; }

    private:
        size_t m_start;  /**< */
        size_t m_length; /**< */
        size_t m_stride; /**< */
    };

    /**
     * @brief
     */
    class GSlice
    {
    public:
        /**
         * @brief
         */
        GSlice() = default;

        /**
         * @brief
         * @param start
         * @param extents
         */
        GSlice(size_t start, std::initializer_list<size_t> extents)
            : m_size(extents.size() > 1 ? (*extents.begin()) * (*extents.end() - 1) : (*extents.begin())),
              m_start(start),
              m_extents(extents)
        {
            m_size = m_extents.size() > 1 ? m_extents.front() * m_extents.back() : m_extents.front();
            if (extents.size() > 2) throw std::invalid_argument("Only 2-dimensional matrices are supported.");
        }

        /**
         * @brief
         * @param start
         * @param extents
         * @param strides
         */
        GSlice(size_t start, std::initializer_list<size_t> extents, std::initializer_list<size_t> strides)
            : m_size(extents.size() > 1 ? (*extents.begin()) * (*extents.end() - 1) : (*extents.begin())),
              m_start(start),
              m_extents(extents),
              m_strides(strides)
        {
            m_size = m_extents.size() > 1 ? m_extents.front() * m_extents.back() : m_extents.front();
            if (extents.size() > 2 || strides.size() > 2) throw std::invalid_argument("Only 2-dimensional matrices are supported.");
        }

        /**
         * @brief
         * @param row
         * @param col
         * @return
         */
        size_t operator()(size_t row, size_t col) const
        {
            if (row > m_extents.front() - 1) throw std::invalid_argument("Invalid row number.");
            if (col > m_extents.back() - 1) throw std::invalid_argument("Invalid column number.");
            return m_start + row * m_strides.front() + col * m_strides.back();
        }

        /**
         * @brief
         * @param index
         * @return
         */
        size_t operator()(size_t index) const
        {
            size_t row = index / colCount();
            size_t col = index % colCount();
            return (*this)(row, col);
        }

        /**
         * @brief
         * @return
         */
        size_t rowCount() const { return m_extents.front(); }

        /**
         * @brief
         * @return
         */
        size_t colCount() const { return m_extents.back(); }

        /**
         * @brief
         * @return
         */
        size_t size() const { return m_size; }

        /**
         * @brief
         * @return
         */
        size_t start() const { return m_start; }

    private:
        size_t              m_size;    /**< */
        size_t              m_start;   /**< */
        std::vector<size_t> m_extents; /**< */
        std::vector<size_t> m_strides; /**< */
    };

}

#endif    // NUMERIX_MATRIXSLICE_HPP

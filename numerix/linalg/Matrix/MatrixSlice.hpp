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

#include "MatrixCommon.hpp"

namespace numerix::linalg
{

    /**
     * @brief The Slice class is used to slice up a Matrix or MatrixView in each direction. A Slice object contains a start
     * position, and (optionally) a length and a stride (i.e. step or jump size). For example, Slice{0,4,2} will being at the
     * first element (element 0), have a length of four elements, and select every 2nd element in the parent Matrix or
     * MatrixView.
     */
    class Slice
    {
    public:
        /**
         * @brief Default constructor with no arguments (members initialized with default values).
         */
        Slice() : m_start(0), m_length(0), m_stride(1) {}

        /**
         * @brief Constructor, taking the starting element as argument. The others have default values.
         * @param start The index of the first element.
         */
        explicit Slice(int start) : m_start(start), m_length(0), m_stride(1) {}

        /**
         * @brief Constructor, taking start, length, and stride as arguments.
         * @param start The index of the first element.
         * @param length The length (number of elements in one direction)
         * @param stride The step size in the parent Matrix or MatrixView.
         */
        Slice(int start, int length, int stride = 1) : m_start(start), m_length(length), m_stride(stride) {}

        /**
         * @brief Function call operator for converting local index, to index of the parent array.
         * @param index Index relative to the slice.
         * @return An index relative to the parent array.
         */
        auto operator()(int index) const
        {
            if (index > m_length) throw std::out_of_range("Slice Bounds Error: Index out of bounds.");
            if (index < 0) throw std::out_of_range("Slice Bounds Error: Index is negative.");
            return m_start + index * m_stride;
        }

        /**
         * @brief Get the start index of the Slice object.
         * @return The start index.
         */
        auto start() const { return m_start; }

        /**
         * @brief Get the length of the Slice object.
         * @return The length.
         */
        auto length() const { return m_length; }

        /**
         * @brief Get the stride of the Slice object.
         * @return The stride.
         */
        auto stride() const { return m_stride; }

    private:
        int m_start;  /**< Index of the start element. */
        int m_length; /**< The length of the Slice object. */
        int m_stride; /**< The stride of the Slice element. */
    };

    namespace impl
    {
        /**
         * @brief A GSlice, or Generalized Slice, is a multi-dimensional slice (although, only two dimensions are supported
         * in the Matrix and MatrixView classes). This class is only used internally; it is not supposed to be used
         * directly by clients.
         */
        class GSlice
        {
        public:
            /**
             * @brief Default constructor, with no arguments.
             */
            GSlice() = default;

            /**
             * @brief Constructor, taking the starting element, and a list of extents as arguments.
             * @param start The index of the starting element.
             * @param extents List of extents in each direction.
             */
            GSlice(int start, std::initializer_list<int> extents)
                : m_size(extents.size() > 1 ? (*extents.begin()) * (*extents.end() - 1) : (*extents.begin())),
                  m_start(start),
                  m_extents(extents)
            {
                m_size = m_extents.size() > 1 ? m_extents.front() * m_extents.back() : m_extents.front();
                if (extents.size() > 2) throw std::invalid_argument("Only 2-dimensional matrices are supported.");
            }

            /**
             * @brief Constructor, taking the starting element, a list of extents, and a list of strides as arguments.
             * @param start The index of the starting element.
             * @param extents List of extents in each direction.
             * @param strides List of strides for each direction.
             */
            GSlice(int start, std::initializer_list<int> extents, std::initializer_list<int> strides)
                : m_size(extents.size() > 1 ? (*extents.begin()) * (*extents.end() - 1) : (*extents.begin())),
                  m_start(start),
                  m_extents(extents),
                  m_strides(strides)
            {
                m_size = m_extents.size() > 1 ? m_extents.front() * m_extents.back() : m_extents.front();
                if (extents.size() > 2 || strides.size() > 2) throw std::invalid_argument("Only 2-dimensional matrices are supported.");
            }

            /**
             * @brief Function call operator for converting coordinates to index of the parent array.
             * @param row Row index, relative to the slice.
             * @param col Column index, relative to the slice.
             * @return Index relative to the parent array.
             */
            int operator()(int row, int col) const
            {
                if (row > m_extents.front() - 1) throw std::invalid_argument("GSlice Bounds Error: Invalid row number.");
                if (col > m_extents.back() - 1) throw std::invalid_argument("GSlice Bounds Error: Invalid column number.");
                return m_start + row * m_strides.front() + col * m_strides.back();
            }

            /**
             * @brief Function call operator for converting local index, to index of the parent array.
             * @param index Index relative to the slice.
             * @return An index relative to the parent array.
             */
            auto operator()(int index) const
            {
                if (index > size() - 1) throw std::out_of_range("GSlice Bounds Error: Index out of bounds.");
                if (index < 0) throw std::out_of_range("GSlice Bounds Error: Index is negative.");
                auto row = index / colCount();
                auto col = index % colCount();
                return (*this)(row, col);
            }

            /**
             * @brief Get the number of rows in the slice.
             * @return The number of rows.
             */
            int rowCount() const { return m_extents.front(); }

            /**
             * @brief Get the number of columns in the slice.
             * @return The number of columns.
             */
            int colCount() const { return m_extents.back(); }

            /**
             * @brief Get the size (number of elements) of the slice.
             * @return The number of elements.
             */
            int size() const { return m_size; }

            /**
             * @brief Get the index of the starting element.
             * @return Index of the start element.
             */
            auto start() const { return m_start; }

        private:
            int              m_size;    /**< The size of the slice. */
            int              m_start;   /**< Index of the starting element. */
            std::vector<int> m_extents; /**< Vector of extents in each direction. */
            std::vector<int> m_strides; /**< Vector of strides for each direction. */
        };
    }    // namespace impl
}    // namespace numerix::linalg

#endif    // NUMERIX_MATRIXSLICE_HPP

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

#ifndef NUMERIX_MATRIXBASE_HPP
#define NUMERIX_MATRIXBASE_HPP

#include <cstddef>

#include "MatrixCommon.hpp"
#include "MatrixSlice.hpp"

namespace numerix::linalg::impl
{

    /**
     * @brief The CRTP base class for the Matrix and MatrixProxy classes.
     * This class provides basic functionality common to both the Matrix and MatrixProxy classes.
     * @tparam DERIVED The derived class used for CRTP
     */
    template<typename DERIVED>
    class MatrixBase
    {
        friend Matrix<typename MatrixTraits<DERIVED>::value_type>;
        friend MatrixView<typename MatrixTraits<DERIVED>::value_type>;
        friend MatrixViewConst<typename MatrixTraits<DERIVED>::value_type>;

    private:
        // ========================================================================================
        // Private member functions
        // ========================================================================================

        /**
         * @brief The index() function converts the row/column indeces to an index of the item in the raw array.
         * @param row The row index of the element.
         * @param col The column index of the element.
         * @return The index of the element in the raw array.
         */
        size_t index(size_t row, size_t col) const
        {
            const DERIVED& derived = static_cast<const DERIVED&>(*this);
            if (row >= derived.m_rowSlice.length()) throw std::out_of_range("Bounds Error: Row index out of bounds.");
            if (col >= derived.m_colSlice.length()) throw std::out_of_range("Bounds Error: Column index out of bounds.");

            size_t start = derived.m_rowSlice.start() * derived.extents().second + derived.m_colSlice.start();
            return start + row * derived.m_rowSlice.stride() + col * derived.m_colSlice.stride();
        }

        /**
         * @brief Check the slice bounds, and modify them to be relative to the parent Matrix, if necessary.
         * @param rowSlice The Slice object representing the rows.
         * @param colSlice The Slice object representing the columns.
         * @return New, bounds checked slices for rows and columns.
         */
        std::pair<Slice, Slice> checkSliceBounds(Slice rowSlice, Slice colSlice) const
        {
            const DERIVED& derived = static_cast<const DERIVED&>(*this);

            // ===== If the length() member is zero, replace it with the distance to the last element.
            auto rSlice =
                Slice(rowSlice.start(), (rowSlice.length() == 0 ? rowCount() - rowSlice.start() : rowSlice.length()), rowSlice.stride());
            auto cSlice =
                Slice(colSlice.start(), (colSlice.length() == 0 ? colCount() - colSlice.start() : colSlice.length()), colSlice.stride());

            // ===== Check that the slices are inside the bounds of the matrix object.
            if (rSlice.length() * rSlice.stride() + rSlice.start() - rSlice.stride() > rowCount() - 1)
                throw std::out_of_range("Slice Error: Row slice out of bounds.");
            if (cSlice.length() * cSlice.stride() + cSlice.start() - cSlice.stride() > colCount() - 1)
                throw std::out_of_range("SliceError: Column slice out of bounds.");

            // ===== If the DERIVED type is a Matrix, modify the row slice accordingly.
            if constexpr (std::same_as<DERIVED, Matrix<typename MatrixTraits<DERIVED>::value_type>>) {
                rSlice = Slice(rSlice.start(), rSlice.length(), rSlice.stride() * derived.m_rowSlice.stride());
            }

            // ===== If the DERIVED type is a MatrixView, convert the slices to represent the parent Matrix object.
            if constexpr (std::same_as<DERIVED, MatrixView<typename MatrixTraits<DERIVED>::value_type>> ||
                          std::same_as<DERIVED, MatrixViewConst<typename MatrixTraits<DERIVED>::value_type>>)
            {
                rSlice = Slice(rSlice.start() * (derived.m_rowSlice.stride() / derived.extents().second) + derived.m_rowSlice.start(),
                               rSlice.length(),
                               rSlice.stride() * derived.m_rowSlice.stride());
                cSlice = Slice(cSlice.start() * derived.m_colSlice.stride() + derived.m_colSlice.start(),
                               cSlice.length(),
                               cSlice.stride() * derived.m_colSlice.stride());
            }

            return std::pair<Slice, Slice>(rSlice, cSlice);
        }

        /**
         * @brief Default constructor
         * @note Constructor made private to avoid non-matrix classes to use it as a base.
         */
        MatrixBase() = default;

    public:
        /**
         * Public alias declatations. To be consistant with standard library containers.
         */
        using value_type = typename MatrixTraits<DERIVED>::value_type;

        // ========================================================================================
        // Special Member Functions (i.e. copy constructors and assignment operators)
        // ========================================================================================

        /**
         * @brief Copy constructor.
         * @param other Object to be copied.
         */
        MatrixBase(const MatrixBase& other) = default;

        /**
         * @brief Move constructor.
         * @param other Object to be moved.
         */
        MatrixBase(MatrixBase&& other) noexcept = default;

        /**
         * @brief Destructor.
         */
        ~MatrixBase() = default;

        /**
         * @brief Copy assignment operator.
         * @param other Object to be copied.
         * @return The copied-to object.
         */
        MatrixBase& operator=(const MatrixBase& other) = default;

        /**
         * @brief Move assignment operator.
         * @param other Object to be moved.
         * @return The moved-to object.
         */
        MatrixBase& operator=(MatrixBase&& other) noexcept = default;

        // ========================================================================================
        // Iterator access (i.e. begin/end iterators) and row/column ranges
        // ========================================================================================

        /**
         * @brief Get a begin-iterator, i.e. an iterator pointing to the first (upper left) element.
         * @return An iterator to the first element
         * @note This function will be disabled for const objects.
         */
        MatrixElemIter<typename MatrixTraits<DERIVED>::value_type> begin()
            requires(!std::same_as<DERIVED, MatrixViewConst<typename MatrixTraits<DERIVED>::value_type>>) && (!std::is_const_v<DERIVED>)
        {
            DERIVED& derived = static_cast<DERIVED&>(*this);
            return MatrixElemIter<value_type>(derived.data(), derived.gslice());
        }

        /**
         * @brief Get a const begin-iterator, i.e. an iterator pointing to the first (upper left) element.
         * @return A const iterator to the first element
         */
        MatrixElemIterConst<typename MatrixTraits<DERIVED>::value_type> begin() const
        {
            const DERIVED& derived = static_cast<const DERIVED&>(*this);
            return MatrixElemIterConst<value_type>(derived.data(), derived.gslice());
        }

        /**
         * @brief Get a const begin-iterator, i.e. an iterator pointing to the first (upper left) element.
         * @return A const iterator to the first element
         */
        MatrixElemIterConst<typename MatrixTraits<DERIVED>::value_type> cbegin() const
        {
            const DERIVED& derived = static_cast<const DERIVED&>(*this);
            return MatrixElemIterConst<value_type>(derived.data(), derived.gslice());
        }

        /**
         * @brief Get an end-iterator, i.e. an iterator pointing to one past the last (lower right) element.
         * @return An iterator to one past the last element.
         * @note This function will be disabled for const objects.
         */
        MatrixElemIter<typename MatrixTraits<DERIVED>::value_type> end()
            requires(!std::same_as<DERIVED, MatrixViewConst<typename MatrixTraits<DERIVED>::value_type>>) && (!std::is_const_v<DERIVED>)
        {
            DERIVED& derived = static_cast<DERIVED&>(*this);
            return MatrixElemIter<value_type>(derived.data(), derived.gslice()).end();
        }

        /**
         * @brief Get a const end-iterator, i.e. an iterator pointing to one past the last (lower right) element.
         * @return A const iterator to one past the last element.
         */
        MatrixElemIterConst<typename MatrixTraits<DERIVED>::value_type> end() const
        {
            const DERIVED& derived = static_cast<const DERIVED&>(*this);
            return MatrixElemIterConst<value_type>(derived.data(), derived.gslice()).end();
        }

        /**
         * @brief Get a const end-iterator, i.e. an iterator pointing to one past the last (lower right) element.
         * @return A const iterator to one past the last element.
         */
        MatrixElemIterConst<typename MatrixTraits<DERIVED>::value_type> cend() const
        {
            const DERIVED& derived = static_cast<const DERIVED&>(*this);
            return MatrixElemIterConst<value_type>(derived.data(), derived.gslice()).end();
        }

        /**
         * @brief Get a collection of all columns in a Matrix of MatrixView object.
         * @return A MatrixCols object with all the columns.
         * @note This function will be disabled for const objects.
         */
        auto cols()
            requires(!std::same_as<DERIVED, MatrixViewConst<value_type>>) && (!std::is_const_v<DERIVED>)
        {
            auto view = (*this)({ 0, rowCount(), 1 }, { 0, colCount(), 1 });
            return MatrixCols<decltype(view)>(view);
        }

        /**
         * @brief Get a collection of all columns in a Matrix of MatrixView object.
         * @return A MatrixColsConst object with all the columns.
         */
        auto cols() const
        {
            auto view = (*this)({ 0, rowCount(), 1 }, { 0, colCount(), 1 });
            return MatrixColsConst<decltype(view)>(view);
        }

        /**
         * @brief Get a collection of all rows in a Matrix of MatrixView object.
         * @return A MatrixCols object with all the rows.
         * @note This function will be disabled for const objects.
         */
        auto rows()
            requires(!std::same_as<DERIVED, MatrixViewConst<value_type>>) && (!std::is_const_v<DERIVED>)
        {
            auto view = (*this)({ 0, rowCount(), 1 }, { 0, colCount(), 1 });
            return MatrixRows<decltype(view)>(view);
        }

        /**
         * @brief Get a collection of all rows in a Matrix of MatrixView object.
         * @return A MatrixCols object with all the rows.
         */
        auto rows() const
        {
            auto view = (*this)({ 0, rowCount(), 1 }, { 0, colCount(), 1 });
            return MatrixRowsConst<decltype(view)>(view);
        }

        // ========================================================================================
        // View and direct element access
        // ========================================================================================

        /**
         * @brief Function call operator overload, for providing Fortran-like element access.
         * @param row The row index.
         * @param col The column index.
         * @return A reference to the matrix element.
         * @note Bounds-checking is done in the index function.
         */
        value_type& operator()(size_t row, size_t col)
            requires(!std::same_as<DERIVED, MatrixViewConst<value_type>>) && (!std::is_const_v<DERIVED>)
        {
            DERIVED& derived = static_cast<DERIVED&>(*this);
            return derived.data()[index(row, col)];
        }

        /**
         * @brief Const function call operator overload, for providing Fortran-like element access.
         * @param row The row index.
         * @param col The column index.
         * @return A const reference to the matrix element.
         * @note Bounds-checking is done in the index function.
         */
        const value_type& operator()(size_t row, size_t col) const
        {
            const DERIVED& derived = static_cast<const DERIVED&>(*this);
            return derived.data()[index(row, col)];
        }

        /**
         * @brief Function call operator overload, for getting a MatrixView object for a subset of the elements.
         * @param rowSlice The Slice object defining the rows.
         * @param colSlice The Slice object defining the columns.
         * @return A MatrixView object with the subset of matrix elements.
         * @note Bounds-checking is done in the checkSliceBounds function.
         */
        auto operator()(const Slice& rowSlice, const Slice& colSlice)
            requires(!std::same_as<DERIVED, MatrixViewConst<value_type>>) && (!std::is_const_v<DERIVED>)
        {
            DERIVED& derived = static_cast<DERIVED&>(*this);

            // ===== Check and modify the slices as required.
            auto [rSlice, cSlice] = checkSliceBounds(rowSlice, colSlice);

            // ===== For Matrix objects, create a MatrixView object using *this as an argument.
            if constexpr (std::same_as<DERIVED, Matrix<typename MatrixTraits<DERIVED>::value_type>> && (!std::is_const_v<DERIVED>))
                return MatrixView<typename MatrixTraits<DERIVED>::value_type>(rSlice, cSlice, static_cast<DERIVED*>(this));

            // ===== For MatrixView objects, create a MatrixView object using the parent Matrix as an argument.
            if constexpr (std::same_as<DERIVED, MatrixView<typename MatrixTraits<DERIVED>::value_type>>)
                return MatrixView<typename MatrixTraits<DERIVED>::value_type>(rSlice, cSlice, derived.m_matrix);
        }

        /**
         * @brief Function call operator overload (const), for getting a MatrixView object for a subset of the elements.
         * @param rowSlice The Slice object defining the rows.
         * @param colSlice The Slice object defining the columns.
         * @return A MatrixView object with the subset of matrix elements.
         * Bounds-checking is done in the checkSliceBounds function.
         */
        auto operator()(const Slice& rowSlice, const Slice& colSlice) const
        {
            const DERIVED& derived = static_cast<const DERIVED&>(*this);

            // ===== Check and modify the slices as required.
            auto [rSlice, cSlice] = checkSliceBounds(rowSlice, colSlice);

            // ===== For Matrix objects, create a MatrixViewConst object using *this as an argument.
            if constexpr (std::same_as<DERIVED, Matrix<typename MatrixTraits<DERIVED>::value_type>>)
                return MatrixViewConst<typename MatrixTraits<DERIVED>::value_type>(rSlice, cSlice, static_cast<const DERIVED*>(this));

            // ===== For MatrixView and MatrixViewConst objects, create a MatrixViewConst object using the parent Matrix as an argument.
            if constexpr (std::same_as<DERIVED, MatrixViewConst<typename MatrixTraits<DERIVED>::value_type>> ||
                          std::same_as<DERIVED, MatrixView<typename MatrixTraits<DERIVED>::value_type>>)
                return MatrixViewConst<typename MatrixTraits<DERIVED>::value_type>(rSlice, cSlice, derived.m_matrix);
        }

        /**
         * @brief Get the row at the given index.
         * @param index The index of the row to get.
         * @return A MatrixView object representing the row.
         * @note This function will be disabled for const objects.
         */
        auto row(size_t index)
            requires(!std::same_as<DERIVED, MatrixViewConst<value_type>>) && (!std::is_const_v<DERIVED>)
        {
            auto a = rowCount();
            if (index >= rowCount()) throw std::out_of_range("Bounds Error: Row index out of bounds.");
            return (*this)({ index, 1, 1 }, { 0, colCount(), 1 });
        }

        /**
         * @brief Get the row at the given index.
         * @param index The index of the row to get.
         * @return A MatrixViewConst object representing the row.
         */
        auto row(size_t index) const
        {
            if (index >= rowCount()) throw std::out_of_range("Bounds Error: Row index out of bounds.");
            return (*this)({ index, 1, 1 }, { 0, colCount(), 1 });
        }

        /**
         * @brief Get the column at the given index.
         * @param index The index of the column to get.
         * @return A MatrixView object representing the column.
         * @note This function will be disabled for const objects.
         */
        auto col(size_t index)
            requires(!std::same_as<DERIVED, MatrixViewConst<value_type>>) && (!std::is_const_v<DERIVED>)
        {
            if (index >= colCount()) throw std::out_of_range("Bounds Error: Column index out of bounds.");
            return (*this)({ 0, rowCount(), 1 }, { index, 1, 1 });
        }

        /**
         * @brief Get the column at the given index.
         * @param index The index of the column to get.
         * @return A MatrixViewConst object representing the column.
         */
        auto col(size_t index) const
        {
            if (index >= colCount()) throw std::out_of_range("Bounds Error: Column index out of bounds.");
            return (*this)({ 0, rowCount(), 1 }, { index, 1, 1 });
        }

        // ========================================================================================
        // Meta data
        // ========================================================================================

        /**
         * @brief Get the number of rows in the Matrix (or MatrixProxy) object.
         * @return The number of rows.
         */
        size_t rowCount() const
        {
            const DERIVED& derived = static_cast<const DERIVED&>(*this);
            return derived.m_rowSlice.length();
        }

        /**
         * @brief Get the number of columns in the Matrix (or MatrixProxy) object.
         * @return The number of columns.
         */
        size_t colCount() const
        {
            const DERIVED& derived = static_cast<const DERIVED&>(*this);
            return derived.m_colSlice.length();
        }

        /**
         * @brief Get the total number of elements of the Matrix or MatrixView object.
         * @return The number of elements.
         */
        size_t size() const { return rowCount() * colCount(); }

        /**
         * @brief Determine if the Matrix (or MatrixView) object is square (i.e. rowCount == colCount).
         * @return If yes, true. Otherwise false.
         */
        bool isSquare() const { return rowCount() == colCount(); }

        // ========================================================================================
        // Arithmetic and non-standard assignment operators.
        // ========================================================================================

        /**
         * @brief Assignment operator
         * @param other The Matrix or MatrixView object from which the elements should be added to the current.
         * @return A reference to the current object.
         * @pre The two matrix objects must have same extents in each dimension.
         * @todo What should happen when the extents do not match?
         */
        template<typename U>
            requires is_matrix<U> && (!std::is_same_v<U, DERIVED>)
        DERIVED& operator=(const U& other)
        {
            // assert(other.rowCount() == rowCount() && other.colCount() == colCount());
            DERIVED& derived = static_cast<DERIVED&>(*this);
            std::copy(other.begin(), other.end(), derived.begin());
            return derived;
        }

        /**
         * @brief Add-assign a scalar value to each matrix element.
         * @param value The scalar value to add.
         * @return A reference to the object.
         */
        DERIVED& operator+=(is_number auto value)
        {
            DERIVED& derived = static_cast<DERIVED&>(*this);
            derived          = (derived + static_cast<value_type>(value));
            return derived;
        }

        /**
         * @brief Add-assign the elements of a matrix to the elements of the current object.
         * @param mat The matrix object with the elements to add.
         * @return A reference to the object.
         */
        DERIVED& operator+=(is_matrix auto mat)
        {
            DERIVED& derived = static_cast<DERIVED&>(*this);
            derived          = (derived + mat);
            return derived;
        }

        /**
         * @brief Subtract-assign a scalar value from each matrix element.
         * @param value The scalar value to subtract.
         * @return A reference to the object.
         */
        DERIVED& operator-=(is_number auto value)
        {
            DERIVED& derived = static_cast<DERIVED&>(*this);
            derived          = (derived - static_cast<value_type>(value));
            return derived;
        }

        /**
         * @brief Subtract-assign the elements of a matrix from the elements of the current object.
         * @param mat The matrix object with the elements to subtract.
         * @return A reference to the object.
         */
        DERIVED& operator-=(is_matrix auto mat)
        {
            DERIVED& derived = static_cast<DERIVED&>(*this);
            derived          = (derived - mat);
            return derived;
        }

        /**
         * @brief Divide-assign each matrix element with a scalar value.
         * @param value The scalar value to divide by.
         * @return A reference to the object.
         */
        DERIVED& operator/=(is_number auto value)
        {
            DERIVED& derived = static_cast<DERIVED&>(*this);
            derived          = (derived / static_cast<value_type>(value));
            return derived;
        }

        /**
         * @brief Multiply-assign each matrix element by a scalar value.
         * @param value The scalar value to multiply by.
         * @return A reference to the object.
         */
        DERIVED& operator*=(is_number auto value)
        {
            DERIVED& derived = static_cast<DERIVED&>(*this);
            derived          = derived * static_cast<value_type>(value);
            return derived;
        }

        // ========================================================================================
        // Static Member Functions
        // ========================================================================================

        /**
         * @brief Create an identity Matrix object (with ones (1) down the diagonal).
         * @param extents The extents of the identity matrix (number of rows and columns must be identical).
         * @return A Matrix object representing the identity matrix.
         */
        static Matrix<typename MatrixTraits<DERIVED>::value_type> CreateIdentityMatrix(size_t extents)
        {
            using value_t = typename MatrixTraits<DERIVED>::value_type;

            Matrix<typename MatrixTraits<DERIVED>::value_type> result(extents, extents);
            for (auto& elem : result) elem = static_cast<value_t>(0.0);
            for (int i = 0; i < extents; ++i) result(i, i) = static_cast<value_t>(1.0);
            return result;
        }
    };

}    // namespace numerix::linalg::impl

#endif    // NUMERIX_MATRIXBASE_HPP

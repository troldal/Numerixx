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

#ifndef NUMERIX_MATRIX_HPP
#define NUMERIX_MATRIX_HPP

#include <cassert>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <span>
#include <exception>

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
        Slice() : m_start(0),
                  m_length(-1),
                  m_stride(1) {}

        /**
         * @brief
         * @param start
         */
        explicit Slice(size_t start) : m_start(start),
                                       m_length(-1),
                                       m_stride(1) {}

        /**
         * @brief
         * @param start
         * @param length
         * @param stride
         */
        Slice(size_t start, size_t length, size_t stride = 1) : m_start(start),
                                                                m_length(length),
                                                                m_stride(stride) {}

        /**
         * @brief
         * @param index
         * @return
         */
        size_t operator()(size_t index) const {
            if (index > m_length) throw std::out_of_range("Index out of bounds.");
            return m_start + index * m_stride;
        }

        /**
         * @brief
         * @return
         */
        size_t start() const {return m_start;}

        /**
         * @brief
         * @return
         */
        size_t length() const {return m_length;}

        /**
         * @brief
         * @return
         */
        size_t stride() const {return m_stride;}

    private:
        size_t m_start; /**< */
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
        GSlice(size_t start, std::initializer_list<size_t> extents) : m_size(extents.size() > 1 ? (*extents.begin()) * (*extents.end() - 1) : (*extents.begin())),
                                                                           m_start(start),
                                                                           m_extents(extents) {
            m_size = m_extents.size() > 1 ? m_extents.front() * m_extents.back() : m_extents.front();
            if (extents.size() > 2) throw std::invalid_argument("Only 2-dimensional matrices are supported.");
        }

        /**
         * @brief
         * @param start
         * @param extents
         * @param strides
         */
        GSlice(size_t start, std::initializer_list<size_t> extents, std::initializer_list<size_t> strides) : m_size(extents.size() > 1 ? (*extents.begin()) * (*extents.end() - 1) : (*extents.begin())),
                                                                                                                  m_start(start),
                                                                                                                  m_extents(extents),
                                                                                                                  m_strides(strides) {
            m_size = m_extents.size() > 1 ? m_extents.front() * m_extents.back() : m_extents.front();
            if (extents.size() > 2 || strides.size() > 2)  throw std::invalid_argument("Only 2-dimensional matrices are supported.");
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
        size_t operator()(size_t index) const {
            size_t row = index / colCount();
            size_t col = index % colCount();
            return (*this)(row, col);
        }

        /**
         * @brief
         * @return
         */
        size_t rowCount() const {
            return m_extents.front();
        }

        /**
         * @brief
         * @return
         */
        size_t colCount() const {
            return m_extents.back();
        }

        /**
         * @brief
         * @return
         */
        size_t size() const {
            return m_size;
        }

        /**
         * @brief
         * @return
         */
        size_t start() const {
            return m_start;
        }

    private:
        size_t                m_size; /**< */
        size_t                m_start; /**< */
        std::vector<size_t> m_extents; /**< */
        std::vector<size_t> m_strides; /**< */
    };

    /**
     * @brief
     * @tparam T
     * @tparam IsConst
     */
    template<typename T, bool IsConst>
    class MatrixIterConcept
    {
        using iter_t = std::conditional_t<IsConst, const T, T>;

        iter_t* m_data; /**< A pointer to the matrix element array. */
        GSlice m_slice; /**< The generalized slice for the data array. */
        size_t m_current; /**< The current index. */

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
        MatrixIterConcept(iter_t* data, GSlice slice, size_t pos = 0) : m_data(data),
                                                           m_slice(slice),
                                                           m_current(pos){}

        /**
         * @brief
         * @return
         */
        MatrixIterConcept end() const {
            return {m_data, m_slice, m_slice.size()};
        }

        /**
         * @brief
         * @return
         */
        MatrixIterConcept& operator++() {
            ++m_current;
            return *this;
        }

        /**
         * @brief
         * @return
         */
        MatrixIterConcept operator++(int) {
            MatrixIterConcept slice = *this;
            ++m_current;
            return slice;
        }

        /**
         * @brief
         * @return
         */
        iter_t& operator*() {
            return m_data[m_slice(m_current)];
        }

        /**
         * @brief
         * @return
         */
        iter_t* operator->() {
            return &m_data[m_slice(m_current)];
        }

        /**
         * @brief
         * @param other
         * @return
         */
        bool operator==(const MatrixIterConcept& other) const {
            return m_current == other.m_current;
        }

        /**
         * @brief
         * @param other
         * @return
         */
        bool operator!=(const MatrixIterConcept& other) const {
            return !(*this == other);
        }

        /**
         * @brief
         * @param other
         * @return
         */
        bool operator<(const MatrixIterConcept& other) const {
            return m_current < other.m_current;
        }
    };

    /**
     * @brief
     */
    template<typename T>
    using MatrixIter = MatrixIterConcept<T, false>;

    /**
     * @brief
     */
    template<typename T>
    using MatrixIterConst = MatrixIterConcept<T, true>;

    template<typename T, bool IsConst>
        requires std::is_arithmetic_v<T> && (!std::is_same_v<T, bool>) && (!std::is_same_v<T, char>)
    class MatrixViewConcept;

    template<typename T>
    using MatrixView = MatrixViewConcept<T, false>;

    template<typename T>
    using MatrixViewConst = MatrixViewConcept<T, true>;

    /*
     *
     */
//    template<typename T>
//        requires std::is_arithmetic_v<T> && (!std::is_same_v<T, bool>) && (!std::is_same_v<T, char>)
//    class MatrixView;

    /*
     *
     */
//    template<typename T>
//        requires std::is_arithmetic_v<T> && (!std::is_same_v<T, bool>) && (!std::is_same_v<T, char>)
//    class MatrixViewConst;

    /*
     *
     */
    template<typename T>
        requires std::is_arithmetic_v<T> && (!std::is_same_v<T, bool>) && (!std::is_same_v<T, char>)
    class Matrix;

    /*
     *
     */
    template<typename T>
    concept is_matrix_type = std::is_same_v<T, Matrix<typename T::value_type>> || std::is_same_v<T, MatrixView<typename T::value_type>> || std::is_same_v<T, MatrixViewConst<typename T::value_type>>;

    /*
     *
     */
    template<typename T>
    concept is_number_type = std::is_arithmetic_v<T> && (!std::is_same_v<T, bool>) && (!std::is_same_v<T, char>);

    namespace impl {

        /*
         *
         */
        template<typename T>
        struct MatrixTraits;

        /**
         * @brief The CRTP base class for the Matrix and MatrixProxy classes.
         * This class provides basic functionality common to both the Matrix and MatrixProxy classes.
         * @tparam DERIVED The derived class used for CRTP
         */
        template<typename DERIVED>
        class MatrixBase {

            /*
             * Alias declarations
             */
            using value_type = typename MatrixTraits<DERIVED>::value_type;

            /**
             * @brief The index() function converts the row/column indeces to an index of the item in the raw array.
             * @param row The row index of the element.
             * @param col The column index of the element.
             * @return The index of the element in the raw array.
             */
            size_t index(size_t row, size_t col) const {
                const DERIVED& derived = static_cast<const DERIVED&>(*this);
                if (row >= derived.m_rowSlice.length() || col >= derived.m_colSlice.length())
                    throw std::out_of_range("Index out of bounds.");

                size_t start = derived.m_rowSlice.start() * derived.dimensions().second + derived.m_colSlice.start();
                return start + row * derived.m_rowSlice.stride() + col * derived.m_colSlice.stride();
            }

        public:

            /**
             * @brief Function call operator overload, for providing Fortran-like element access.
             * @param row The row index.
             * @param col The column index.
             * @return A reference to the matrix element.
             */
            value_type& operator()(size_t row, size_t col) requires (!std::is_same_v<DERIVED, MatrixViewConst<value_type>>) {
                DERIVED& derived = static_cast<DERIVED&>(*this);
                return derived.data()[index(row, col)];
            }

            /**
             * @brief Const function call operator overload, for providing Fortran-like element access.
             * @param row The row index.
             * @param col The column index.
             * @return A const reference to the matrix element.
             */
            const value_type& operator()(size_t row, size_t col) const {
                const DERIVED& derived = static_cast<const DERIVED&>(*this);
                return derived.data()[index(row, col)];
            }

            /**
             * @brief Function call operator overload, for getting a MatrixProxy object for a subset of the elements.
             * @param rowSlice The Slice object defining the rows.
             * @param colSlice The Slice object defining the columns.
             * @return A MatrixProxy object with the subset of matrix elements.
             */
            auto operator()(const Slice& rowSlice, const Slice& colSlice) requires (!std::is_same_v<DERIVED, MatrixViewConst<value_type>>) {
                DERIVED& derived = static_cast<DERIVED&>(*this);
                return derived.slice(rowSlice, colSlice);
            }

            /**
             * @brief
             * @param rowSlice
             * @param colSlice
             * @return
             */
            auto operator()(const Slice& rowSlice, const Slice& colSlice) const {
                const DERIVED& derived = static_cast<const DERIVED&>(*this);
                return derived.slice(rowSlice, colSlice);
            }

            /**
             * @brief Get the number of rows in the Matrix (or MatrixProxy) object.
             * @return The number of rows.
             */
            size_t rowCount() const {
                const DERIVED& derived = static_cast<const DERIVED&>(*this);
                return derived.m_rowSlice.length();
            }

            /**
             * @brief Get the number of columns in the Matrix (or MatrixProxy) object.
             * @return The number of columns.
             */
            size_t colCount() const {
                const DERIVED& derived = static_cast<const DERIVED&>(*this);
                return derived.m_colSlice.length();
            }

            /**
             * @brief Determine if the Matrix (or MatrixProxy) object is square (i.e. rowCount == colCount).
             * @return If yes, true. Otherwise false.
             */
            bool isSquare() const {
                return rowCount() == colCount();
            }

            /**
             * @brief Get a begin-iterator, i.e. an iterator pointing to the first (upper left) element.
             * @return An iterator to the first element
             */
            MatrixIter<value_type> begin() requires (!std::is_same_v<DERIVED, MatrixViewConst<value_type>>) {
                DERIVED& derived = static_cast<DERIVED&>(*this);
                return MatrixIter<value_type>(derived.data(), derived.gslice());
            }

            /**
             * @brief Get a const begin-iterator, i.e. an iterator pointing to the first (upper left) element.
             * @return A const iterator to the first element
             */
            MatrixIterConst<value_type> begin() const {
                const DERIVED& derived = static_cast<const DERIVED&>(*this);
                return MatrixIterConst<value_type>(derived.data(), derived.gslice());
            }

            /**
             * @brief
             * @return
             */
            MatrixIterConst<value_type> cbegin() const {
                const DERIVED& derived = static_cast<const DERIVED&>(*this);
                return MatrixIterConst<value_type>(derived.data(), derived.gslice());
            }

            /**
             * @brief Get an end-iterator, i.e. an iterator pointing to one past the last (lower right) element.
             * @return An iterator to one past the last element.
             */
            MatrixIter<value_type> end() requires (!std::is_same_v<DERIVED, MatrixViewConst<value_type>>) {
                DERIVED& derived = static_cast<DERIVED&>(*this);
                return MatrixIter<value_type>(derived.data(), derived.gslice()).end();
            }

            /**
             * @brief Get a const end-iterator, i.e. an iterator pointing to one past the last (lower right) element.
             * @return A const iterator to one past the last element.
             */
            MatrixIterConst<value_type> end() const {
                const DERIVED& derived = static_cast<const DERIVED&>(*this);
                return MatrixIterConst<value_type>(derived.data(), derived.gslice()).end();
            }

            /**
             * @brief Get a const end-iterator, i.e. an iterator pointing to one past the last (lower right) element.
             * @return A const iterator to one past the last element.
             */
            MatrixIterConst<value_type> cend() const {
                const DERIVED& derived = static_cast<const DERIVED&>(*this);
                return MatrixIterConst<value_type>(derived.data(), derived.gslice()).end();
            }

            /**
             * @brief
             * @param index
             * @return
             */
            auto row(size_t index) requires (!std::is_same_v<DERIVED, MatrixViewConst<value_type>>) {
                return (*this)({index, 1, 1}, {0, colCount(), 1});
            }

            /**
             * @brief
             * @param index
             * @return
             */
            auto row(size_t index) const {
                return (*this)({index, 1, 1}, {0, colCount(), 1});
            }

            /**
             * @brief
             * @param index
             * @return
             */
            auto col(size_t index) requires (!std::is_same_v<DERIVED, MatrixViewConst<value_type>>) {
                return (*this)({0, rowCount(), 1}, {index, 1, 1});
            }

            /**
             * @brief
             * @param index
             * @return
             */
            auto col(size_t index) const {
                return (*this)({0, rowCount(), 1}, {index, 1, 1});
            }

            /**
             * @brief
             * @tparam T
             * @param other
             * @return
             */
            template<typename T>
            requires is_matrix_type<T>
            DERIVED& operator=(const T& other) {
                //assert(other.rowCount() == rowCount() && other.colCount() == colCount());
                DERIVED& derived = static_cast<DERIVED&>(*this);
                std::copy(other.begin(), other.end(), derived.begin());
                return derived;
            }

            /**
             * @brief Matrix/Matrix add assignment operator. I.e. adds two Matrix of MatrixProxy objects, with the same size.
             * @tparam T The type of matrix to be added, i.e. Matrix or MatrixProxy.
             * @param other The Matrix of MatrixProxy to be added.
             * @return A reference to the current object, after addition.
             */
            template<typename T>
            requires is_matrix_type<T>
            DERIVED& operator+=(const T& other) {
                assert(other.rowCount() == rowCount() && other.colCount() == colCount());
                DERIVED& derived = static_cast<DERIVED&>(*this);
                std::transform(other.begin(), other.end(), derived.begin(), derived.begin(), std::plus<value_type>());
                return derived;
            }

            /**
             * @brief Scalar add assignment operator. I.e. adds a scalar value to all elements of the Matrix.
             * @tparam T The type of scalar to add. Can be any number type, except char and bool.
             * @param value The scalar value.
             * @return A reference to the current object, after addition.
             */
            template<typename T>
            requires std::is_arithmetic_v<T> && (!std::is_same_v<T, bool>) && (!std::is_same_v<T, char>)
            DERIVED& operator+=(T value) {
                DERIVED& derived = static_cast<DERIVED&>(*this);
                derived = derived + static_cast<value_type>(value);
                return derived;
            }

            /**
             * @brief
             * @tparam T
             * @param other
             * @return
             */
            template<typename T>
                requires is_matrix_type<T>
            DERIVED& operator-=(const T& other) {
                assert(other.rowCount() == rowCount() && other.colCount() == colCount());
                DERIVED& derived = static_cast<DERIVED&>(*this);
                std::transform(other.begin(), other.end(), derived.begin(), derived.begin(), std::minus<value_type>());
                return derived;
            }

            /**
             * @brief
             * @tparam T
             * @param value
             * @return
             */
            template<typename T>
                requires std::is_arithmetic_v<T> && (!std::is_same_v<T, bool>) && (!std::is_same_v<T, char>)
            DERIVED& operator-=(T value) {
                DERIVED& derived = static_cast<DERIVED&>(*this);
                derived = derived - static_cast<value_type>(value);
                return derived;
            }

            /**
             * @brief
             * @tparam T
             * @param value
             * @return
             */
            template<typename T>
                requires std::is_arithmetic_v<T> && (!std::is_same_v<T, bool>) && (!std::is_same_v<T, char>)
            DERIVED& operator/=(T value) {
                DERIVED& derived = static_cast<DERIVED&>(*this);
                derived = derived / static_cast<value_type>(value);
                return derived;
            }

            /**
             * @brief
             * @tparam T
             * @param value
             * @return
             */
            template<typename T>
                requires std::is_arithmetic_v<T> && (!std::is_same_v<T, bool>) && (!std::is_same_v<T, char>)
            DERIVED& operator*=(T value) {
                DERIVED& derived = static_cast<DERIVED&>(*this);
                derived = derived * static_cast<value_type>(value);
                return derived;
            }
        };
    }

    /**
     * @brief The Matrix class is the main abstraction for matrices. It derives from the MatrixBase class using CRTP.
     * The Matrix class owns the underlying data, unlike the MatrixProxy class which simply is a view into a Matrix object.
     * @tparam T The type of the underlying matrix elements. Default is double, but can be any kind of floating point or integer number.
     */
    template<typename T = double>
        requires std::is_arithmetic_v<T> && (!std::is_same_v<T, bool>) && (!std::is_same_v<T, char>)
    class Matrix : public impl::MatrixBase<Matrix<T>>
    {
    private:

        /**
         * Friend declarations. Necessary for the base class to access elements of Matrix using CRTP.
         */
        friend impl::MatrixBase<Matrix<T>>;
        friend MatrixView<T>;
        friend MatrixViewConst<T>;

        /**
         * Private alias declatations.
         */
        using parent = impl::MatrixBase<Matrix<T>>;

        // TODO: Consider using an Extents class, rather than row and col count.
        std::vector<T> m_data; /**< The underlying array of matrix elements. */
        // TODO: Delete m_rowCount and m_colCount from the internal datastructure.
        int m_rowCount; /**< The number of rows in the matrix. */
        int m_colCount; /**< The number of columns in the matrix. */
        Slice m_rowSlice; /**< The Slice describing the rows. Required to provide a common interface with MatrixProxy. */
        Slice m_colSlice; /**< The Slice describing the columns. Required to provide a common interface with MatrixProxy. */

        /**
         * @brief Create a MatrixProxy object, based on the input row and column slices.
         * @param rowSlice The Slice object describing the row elements.
         * @param colSlice The Slice object describing the column elements.
         * @return A MatrixProxy object for the subset of Matrix elements included in the slices.
         */
        auto slice(const Slice& rowSlice, const Slice& colSlice); // NOTE: Implementation can be found after the MatrixProxy class.

        /**
         * @brief
         * @param rowSlice
         * @param colSlice
         * @return
         */
        auto slice(const Slice& rowSlice, const Slice& colSlice) const;

        /**
         * @brief
         * @return
         */
        auto gslice() const {
            auto start = m_rowSlice.start() * dimensions().first + m_colSlice.start();
            return GSlice(start, {m_rowSlice.length(), m_colSlice.length()}, {m_rowSlice.stride(), m_colSlice.stride()});
        }

    public:

        /**
         * Public alias declatations. To be consistant with standard library containers.
         */
        using value_type = T;

        /**
         * @brief Constructor taking the number of rows and columns.
         * @param rows Number of rows.
         * @param cols Number of cols.
         */
        Matrix(int rows, int cols) : m_rowCount{rows},
                                     m_colCount{cols},
                                     m_data(rows * cols),
                                     m_rowSlice(0, rows, cols),
                                     m_colSlice(0, cols, 1)
        {}

        /**
         * @brief Copy constructor.
         * @param other Object to be copied.
         */
        Matrix(const Matrix& other) = default;

        /**
         * @brief Move constructor.
         * @param other Object to be moved.
         */
        Matrix(Matrix&& other) noexcept = default;

        /**
         * @brief Destructor.
         */
        ~Matrix() = default;

        /**
         * @brief Copy assignment operator.
         * @param other Object to be copied.
         * @return The copied-to object.
         */
        //Matrix& operator=(const Matrix& other) { return parent::operator=(other); }
        //Matrix& operator=(const MatrixProxy<T>& other) { return parent::operator=(other); }


        /**
         * @brief Move assignment operator.
         * @param other Object to be moved.
         * @return The moved-to object.
         */
        Matrix& operator=(Matrix&& other) noexcept = default;

        /**
         * @brief
         * @return
         */
        auto dimensions() const {
            return std::make_pair(parent::rowCount(), parent::colCount());
        }

        /**
         * @brief
         * @param i
         * @return
         * @todo Needs new implementation.
         */
        auto operator[](const int i) {
            return std::span (&m_data.data()[i * m_colSlice.length()], m_colSlice.length());
        }

        /**
         * @brief
         * @param i
         * @return
         * @todo Needs new implementation.
         */
        const auto operator[](int i) const {
            return std::span (&m_data.data()[i * m_colSlice.length()], m_colSlice.length());
        }

        /**
         * @brief Access the raw array of Matrix elements.
         * @return A pointer to the first element.
         */
        T* data() {
            return m_data.data();
        }

        /**
         * @brief Access the raw array of Matrix elements.
         * @return A const pointer to the first element.
         */
        const T* data() const {
            return m_data.data();
        }

        /**
         * @brief
         * @param other
         * @return
         */
//        Matrix& operator*=(const Matrix& other) {
//            assert(m_colSlice.length() == other.rowCount());
//
//            std::vector<T> result;
//            for (int r = 0; r < m_rowSlice.length(); ++r)
//                for (int c = 0; c < other.colCount(); ++c)
//                    result.push_back(std::inner_product(parent::row(r).begin(), parent::row(r).end(), other.col(c).begin(), 0.0));
//
//            m_data = result;
//            m_colCount = other.colCount();
//            return *this;
//        }
    };

    template<typename T, bool IsConst>
        requires std::is_arithmetic_v<T> && (!std::is_same_v<T, bool>) && (!std::is_same_v<T, char>)
    class MatrixViewConcept : public impl::MatrixBase<MatrixViewConcept<T, IsConst>>
    {
    private:

        /*
         * Friend declarations. Necessary for the base class to access elements of MatrixProxy using CRTP.
         */
        friend impl::MatrixBase<MatrixViewConcept<T, IsConst>>;

        /**
         * Private alias declatations.
         */
        using parent = impl::MatrixBase<MatrixViewConcept<T, IsConst>>;
        using data_t = std::conditional_t<IsConst, const T, T>;
        using matrix_t = std::conditional_t<IsConst, const Matrix<T>, Matrix<T>>;


        /**
         * @brief
         * @return
         */
        auto gslice() const {
            auto start = m_rowSlice.start() * m_matrix->dimensions().first + m_colSlice.start();
            return GSlice(start, {m_rowSlice.length(), m_colSlice.length()}, {m_rowSlice.stride(), m_colSlice.stride()});
        }

    public:

        /**
         * Alias declatations. To be consistant with standard library containers.
         */
        using value_type = T;

        /**
         * @brief
         * @param rowSlice
         * @param colSlice
         * @param data
         */
        MatrixViewConcept(const Slice& rowSlice, const Slice& colSlice, matrix_t* matrix) : m_rowSlice(rowSlice),
                                                                                      m_colSlice(colSlice),
                                                                                      m_matrix(matrix) {}

        /**
         * @brief
         * @param other
         * @return
         */
        //MatrixProxy& operator=(const MatrixProxy& other) { return parent::operator=(other); }
        //MatrixProxy& operator=(const Matrix<T>& other) { return parent::operator=(other); }

        /**
         * @brief
         * @param rowSlice
         * @param colSlice
         * @return
         */
        auto slice(const Slice& rowSlice, const Slice& colSlice) {

            //if (rowSlice.start() + rowSlice.stride() * (rowSlice.length() - 1) >= m_slice.rowCount())
            //    throw std::out_of_range("Row index out of bounds.");

            //if (colSlice.start() + colSlice.stride() * (colSlice.length() - 1) >= m_slice.colCount())
            //    throw std::out_of_range("Column index out of bounds.");
            auto rSlice = Slice(rowSlice.start() * (m_rowSlice.stride()/dimensions().second) + m_rowSlice.start(), rowSlice.length(), rowSlice.stride() * m_rowSlice.stride());
            auto cSlice = Slice(colSlice.start() * m_colSlice.stride() + m_colSlice.start(), colSlice.length(), colSlice.stride() * m_colSlice.stride());

            return MatrixViewConcept(rSlice, cSlice, m_matrix);
        }

        auto slice(const Slice& rowSlice, const Slice& colSlice) const {

            //if (rowSlice.start() + rowSlice.stride() * (rowSlice.length() - 1) >= m_slice.rowCount())
            //    throw std::out_of_range("Row index out of bounds.");

            //if (colSlice.start() + colSlice.stride() * (colSlice.length() - 1) >= m_slice.colCount())
            //    throw std::out_of_range("Column index out of bounds.");
            auto rSlice = Slice(rowSlice.start() * (m_rowSlice.stride()/dimensions().second) + m_rowSlice.start(), rowSlice.length(), rowSlice.stride() * m_rowSlice.stride());
            auto cSlice = Slice(colSlice.start() * m_colSlice.stride() + m_colSlice.start(), colSlice.length(), colSlice.stride() * m_colSlice.stride());

            return MatrixViewConcept(rSlice, cSlice, m_matrix);
        }

        /**
         * @brief
         * @return
         */
        auto dimensions() const {
            return m_matrix->dimensions();
        }

    private:

        /**
         * @brief
         * @return
         */
        data_t* data() {
            return m_matrix->data();
        }

        /**
         * @brief
         * @return
         */
        const T* data() const {
            return m_matrix->data();
        }

        Slice m_rowSlice; /**< */
        Slice m_colSlice; /**< */
        matrix_t* m_matrix; /**< */
    };


    /**
     * @brief The MatrixProxy class is a view into a subset of the elements of a Matrix. It derives from the MatrixBase class using CRTP.
     * An object of the MatrixProxy class does not own the data; the data is owned by the corresponding Matrix object.
     * @tparam T The type of the underlying matrix elements. Default is double, but can be any kind of floating point or integer number.
     * @todo Consider putting in the impl namespace, as MatrixProxy objects should only be created from a Matrix object.
     */
//    template<typename T>
//    requires std::is_arithmetic_v<T> && (!std::is_same_v<T, bool>) && (!std::is_same_v<T, char>)
//    class MatrixView : public impl::MatrixBase<MatrixView<T>>
//    {
//    private:
//
//        /*
//         * Friend declarations. Necessary for the base class to access elements of MatrixProxy using CRTP.
//         */
//        friend impl::MatrixBase<MatrixView<T>>;
//
//        /**
//         * Private alias declatations.
//         */
//        using parent = impl::MatrixBase<MatrixView<T>>;
//
//        /**
//         * @brief
//         * @return
//         */
//        auto gslice() const {
//            auto start = m_rowSlice.start() * m_matrix->dimensions().first + m_colSlice.start();
//            return GSlice(start, {m_rowSlice.length(), m_colSlice.length()}, {m_rowSlice.stride(), m_colSlice.stride()});
//        }
//
//    public:
//
//        /**
//         * Alias declatations. To be consistant with standard library containers.
//         */
//        using value_type = T;
//
//        /**
//         * @brief
//         * @param rowSlice
//         * @param colSlice
//         * @param data
//         */
//        MatrixView(const Slice& rowSlice, const Slice& colSlice, Matrix<T>* matrix) : m_rowSlice(rowSlice),
//                                                                                       m_colSlice(colSlice),
//                                                                                       m_matrix(matrix) {}
//
//        /**
//         * @brief
//         * @param other
//         * @return
//         */
//        //MatrixProxy& operator=(const MatrixProxy& other) { return parent::operator=(other); }
//        //MatrixProxy& operator=(const Matrix<T>& other) { return parent::operator=(other); }
//
//        /**
//         * @brief
//         * @param rowSlice
//         * @param colSlice
//         * @return
//         */
//        auto slice(const Slice& rowSlice, const Slice& colSlice) {
//
//            //if (rowSlice.start() + rowSlice.stride() * (rowSlice.length() - 1) >= m_slice.rowCount())
//            //    throw std::out_of_range("Row index out of bounds.");
//
//            //if (colSlice.start() + colSlice.stride() * (colSlice.length() - 1) >= m_slice.colCount())
//            //    throw std::out_of_range("Column index out of bounds.");
//            auto rSlice = Slice(rowSlice.start() * (m_rowSlice.stride()/dimensions().second) + m_rowSlice.start(), rowSlice.length(), rowSlice.stride() * m_rowSlice.stride());
//            auto cSlice = Slice(colSlice.start() * m_colSlice.stride() + m_colSlice.start(), colSlice.length(), colSlice.stride() * m_colSlice.stride());
//
//            return MatrixView(rSlice, cSlice, m_matrix);
//        }
//
//        auto slice(const Slice& rowSlice, const Slice& colSlice) const {
//
//            //if (rowSlice.start() + rowSlice.stride() * (rowSlice.length() - 1) >= m_slice.rowCount())
//            //    throw std::out_of_range("Row index out of bounds.");
//
//            //if (colSlice.start() + colSlice.stride() * (colSlice.length() - 1) >= m_slice.colCount())
//            //    throw std::out_of_range("Column index out of bounds.");
//            auto rSlice = Slice(rowSlice.start() * (m_rowSlice.stride()/dimensions().second) + m_rowSlice.start(), rowSlice.length(), rowSlice.stride() * m_rowSlice.stride());
//            auto cSlice = Slice(colSlice.start() * m_colSlice.stride() + m_colSlice.start(), colSlice.length(), colSlice.stride() * m_colSlice.stride());
//
//            return MatrixViewConst(rSlice, cSlice, m_matrix);
//        }
//
//        /**
//         * @brief
//         * @return
//         */
//        auto dimensions() const {
//            return m_matrix->dimensions();
//        }
//
//    private:
//
//        /**
//         * @brief
//         * @return
//         */
//        T* data() {
//            return m_matrix->data();
//        }
//
//        /**
//         * @brief
//         * @return
//         */
//        const T* data() const {
//            return m_matrix->data();
//        }
//
//        Slice m_rowSlice; /**< */
//        Slice m_colSlice; /**< */
//        Matrix<T>* m_matrix; /**< */
//    };

    /**
     * @brief The MatrixProxyConst class is a (const) view into a subset of the elements of a Matrix. It derives from the MatrixBase class using CRTP.
     * An object of the MatrixProxyConst class does not own the data; the data is owned by the corresponding Matrix object. The underlying data cannot be manipulated
     * through the MatrixProxyConst interface.
     * @tparam T The type of the underlying matrix elements. Default is double, but can be any kind of floating point or integer number.
     * @todo Consider putting in the impl namespace, as MatrixProxy objects should only be created from a Matrix object.
     */
//    template<typename T>
//        requires std::is_arithmetic_v<T> && (!std::is_same_v<T, bool>) && (!std::is_same_v<T, char>)
//    class MatrixViewConst : public impl::MatrixBase<MatrixViewConst<T>>
//    {
//    private:
//
//        /*
//         * Friend declarations. Necessary for the base class to access elements of MatrixProxy using CRTP.
//         */
//        friend impl::MatrixBase<MatrixViewConst<T>>;
//
//        /**
//         * Private alias declatations.
//         */
//        using parent = impl::MatrixBase<MatrixViewConst<T>>;
//
//        /**
//         * @brief
//         * @return
//         */
//        auto gslice() const {
//            auto start = m_rowSlice.start() * m_matrix->dimensions().first + m_colSlice.start();
//            return GSlice(start, {m_rowSlice.length(), m_colSlice.length()}, {m_rowSlice.stride(), m_colSlice.stride()});
//        }
//
//    public:
//
//        /**
//         * Alias declatations. To be consistant with standard library containers.
//         */
//        using value_type = T;
//
//        /**
//         * @brief
//         * @param rowSlice
//         * @param colSlice
//         * @param data
//         */
//        MatrixViewConst(const Slice& rowSlice, const Slice& colSlice, const Matrix<T>* matrix) : m_rowSlice(rowSlice),
//                                                                                       m_colSlice(colSlice),
//                                                                                       m_matrix(matrix) {}
//
//        /**
//         * @brief
//         * @param other
//         * @return
//         */
//        //MatrixProxy& operator=(const MatrixProxy& other) { return parent::operator=(other); }
//        //MatrixProxy& operator=(const Matrix<T>& other) { return parent::operator=(other); }
//
//        /**
//         * @brief
//         * @param rowSlice
//         * @param colSlice
//         * @return
//         */
//        auto slice(const Slice& rowSlice, const Slice& colSlice) const {
//
//            //if (rowSlice.start() + rowSlice.stride() * (rowSlice.length() - 1) >= m_slice.rowCount())
//            //    throw std::out_of_range("Row index out of bounds.");
//
//            //if (colSlice.start() + colSlice.stride() * (colSlice.length() - 1) >= m_slice.colCount())
//            //    throw std::out_of_range("Column index out of bounds.");
//            auto rSlice = Slice(rowSlice.start() * (m_rowSlice.stride()/dimensions().second) + m_rowSlice.start(), rowSlice.length(), rowSlice.stride() * m_rowSlice.stride());
//            auto cSlice = Slice(colSlice.start() * m_colSlice.stride() + m_colSlice.start(), colSlice.length(), colSlice.stride() * m_colSlice.stride());
//
//            return MatrixViewConst(rSlice, cSlice, m_matrix);
//        }
//
//        /**
//         * @brief
//         * @return
//         */
//        auto dimensions() const {
//            return m_matrix->dimensions();
//        }
//
//    private:
//
//        /**
//         * @brief
//         * @return
//         */
//        const T* data() const {
//            return m_matrix->data();
//        }
//
//        Slice m_rowSlice; /**< */
//        Slice m_colSlice; /**< */
//        const Matrix<T>* m_matrix; /**< */
//    };


    /*
     * Implementation if the (non-const) slice function for the Matrix class.
     */
    template<typename T>
    requires std::is_arithmetic_v<T> && (!std::is_same_v<T, bool>) && (!std::is_same_v<T, char>)
    auto Matrix<T>::slice(const Slice& rowSlice, const Slice& colSlice)
    {

        //            if (rowSlice.start() + rowSlice.stride() * (rowSlice.length() - 1) >= m_slice.rowCount())
        if (rowSlice.length() > m_rowSlice.length())
            throw std::out_of_range("Row index out of bounds.");

        //            if (colSlice.start() + colSlice.stride() * (colSlice.length() - 1) >= m_slice.colCount())
        if (colSlice.length() > m_colSlice.length())
            throw std::out_of_range("Column index out of bounds.");

        auto rSlice = Slice(rowSlice.start(), rowSlice.length(), rowSlice.stride() * m_rowSlice.stride());

        return MatrixView<T>(rSlice, colSlice, this);
    }

    /*
     * Implementation if the (const) slice function for the Matrix class.
     */
    template<typename T>
    requires std::is_arithmetic_v<T> && (!std::is_same_v<T, bool>) && (!std::is_same_v<T, char>)
    auto Matrix<T>::slice(const Slice& rowSlice, const Slice& colSlice) const
    {

        //            if (rowSlice.start() + rowSlice.stride() * (rowSlice.length() - 1) >= m_slice.rowCount())
        if (rowSlice.length() > m_rowSlice.length())
            throw std::out_of_range("Row index out of bounds.");

        //            if (colSlice.start() + colSlice.stride() * (colSlice.length() - 1) >= m_slice.colCount())
        if (colSlice.length() > m_colSlice.length())
            throw std::out_of_range("Column index out of bounds.");

        auto rSlice = Slice(rowSlice.start(), rowSlice.length(), rowSlice.stride() * m_rowSlice.stride());

        return MatrixViewConst<T>(rSlice, colSlice, this);
    }

    namespace impl {

        /**
         * @brief
         * @tparam T
         */
        template<typename T>
        struct MatrixTraits<Matrix<T>> {
            using value_type = T;
        };

        /**
         * @brief
         * @tparam T
         */
        template<typename T>
        struct MatrixTraits<MatrixView<T>> {
            using value_type = T;
        };

        /**
         * @brief
         * @tparam T
         */
        template<typename T>
        struct MatrixTraits<MatrixViewConst<T>> {
            using value_type = T;
        };
    }

    /**
     * @brief
     * @tparam T
     * @param out
     * @param mat
     * @return
     */
    template <typename T>
    requires std::is_same_v<T, Matrix<typename T::value_type>> || std::is_same_v<T, MatrixView<typename T::value_type>> || std::is_same_v<T, MatrixViewConst<typename T::value_type>>
    std::ostream& operator<<(std::ostream& out, const T& mat) {
        auto printRow = [&](int row) {
            std::cout << "{ ";
            std::cout << std::fixed << std::setw(2);
            for (auto col = 0; col < mat.colCount(); ++col) std::cout << mat(row, col) << " ";
            std::cout << "}\n";
        };

        for (int i = 0; i < mat.rowCount(); ++i)
            printRow(i);

        return out;
    }

    /**
     * @brief
     * @tparam T
     * @param a
     * @param b
     * @return
     */
    template<typename T, typename U>
    requires is_matrix_type<T> && is_matrix_type<U>
    inline Matrix<typename std::common_type<typename T::value_type, typename U::value_type>::type> operator+(const T& a, const U& b) {

        auto result = a;
        result += b;
        return result;
    }

    /**
     * @brief
     * @tparam T
     * @tparam U
     * @param a
     * @param b
     * @return
     */
    template<typename T, typename U>
    requires is_matrix_type<T> && is_number_type<U>
    inline Matrix<typename std::common_type<typename T::value_type, U>::type> operator+(const T& a, U b) {

        using result_t = Matrix<typename std::common_type<typename T::value_type, U>::type>;

        result_t result {a.rowCount(), a.colCount()};
        std::transform(a.begin(), a.end(), result.begin(), [&](const auto& v) { return static_cast<typename result_t::value_type>(v) + b; });
        return result;
    }

    /**
     * @brief
     * @tparam T
     * @param a
     * @param b
     * @return
     */
    template<typename T, typename U>
    requires is_matrix_type<T> && is_matrix_type<U>
    inline Matrix<typename std::common_type<typename T::value_type, typename U::value_type>::type> operator-(const T& a, const U& b) {

        auto result = a;
        result -= b;
        return result;
    }

    /**
     * @brief
     * @tparam T
     * @tparam U
     * @param a
     * @param b
     * @return
     */
    template<typename T, typename U>
    requires is_matrix_type<T> && is_number_type<U>
    inline Matrix<typename std::common_type<typename T::value_type, U>::type> operator-(const T& a, U b)  {

        using result_t = Matrix<typename std::common_type<typename T::value_type, U>::type>;

        result_t result {a.rowCount(), a.colCount()};
        std::transform(a.begin(), a.end(), result.begin(), [&](const auto& v) { return static_cast<typename result_t::value_type>(v) - b; });
        return result;
    }

    /**
     * @brief Compute the result of Matrix-Matrix multiplication.
     * @tparam T The type of the first Matrix (Matrix or MatrixProxy)
     * @tparam U The type of the second Matrix (Matrix or MatrixProxy)
     * @param a The first Matrix
     * @param b The second Matrix
     * @return A Metrix object with the result of the Matrix multiplication
     */
    template<typename T, typename U>
        requires is_matrix_type<T> && is_matrix_type<U>
    inline Matrix<typename std::common_type_t<typename T::value_type, typename U::value_type>> operator*(const T& a, const U& b) {

        using result_t = Matrix<typename std::common_type_t<typename T::value_type, typename U::value_type>>;

        result_t result {a.rowCount(), b.colCount()};
        for (int i = 0; i < result.rowCount(); ++i)
            for (int j = 0; j < result.colCount(); ++j)
                result(i, j) = std::inner_product(a.row(i).begin(), a.row(i).end(), b.col(j).begin(), static_cast<typename result_t::value_type>(0.0));

        return result;
    }

    /**
     * @brief
     * @tparam T
     * @tparam U
     * @param a
     * @param b
     * @return
     */
    template<typename T, typename U>
        requires is_matrix_type<T> && is_number_type<U>
    inline Matrix<typename std::common_type<typename T::value_type, U>::type> operator*(const T& a, U b)  {

        using result_t = Matrix<typename std::common_type<typename T::value_type, U>::type>;

        result_t result {a.rowCount(), a.colCount()};
        std::transform(a.begin(), a.end(), result.begin(), [&](const auto& v) { return static_cast<typename result_t::value_type>(v) * b; });
        return result;
    }

    /**
     * @brief
     * @tparam T
     * @tparam U
     * @param a
     * @param b
     * @return
     */
    template<typename T, typename U>
        requires is_matrix_type<U> && is_number_type<T>
    inline Matrix<typename std::common_type<typename U::value_type, T>::type> operator*(const T a, U& b)  {

        using result_t = Matrix<typename std::common_type<typename U::value_type, T>::type>;

        result_t result {b.rowCount(), b.colCount()};
        std::transform(b.begin(), b.end(), result.begin(), [&](const auto& v) { return static_cast<typename result_t::value_type>(v) * a; });
        return result;
    }

    /**
     * @brief
     * @tparam T
     * @tparam U
     * @param a
     * @param b
     * @return
     */
    template<typename T, typename U>
    requires is_matrix_type<T> && is_number_type<U>
    inline Matrix<typename std::common_type<typename T::value_type, U>::type> operator/(const T& a, U b)  {

        using result_t = Matrix<typename std::common_type<typename T::value_type, U>::type>;

        result_t result {a.rowCount(), a.colCount()};
        std::transform(a.begin(), a.end(), result.begin(), [&](const auto& v) { return static_cast<typename result_t::value_type>(v) / b; });
        return result;
    }

    /**
     * @brief Compute the transpose of the input Matrix (or MatrixProxy.
     * @tparam T The type of the input Matrix. Must be of Matrix or MatrixProxy types.
     * @param mat The input Matrix
     * @return A new Matrix object, with the transpose of the input.
     * @note This function does not modify the input. It merely produces a new Matrix object with the transpose.
     */
    template<typename T>
    requires is_matrix_type<T>
    inline Matrix<typename T::value_type> transpose(const T& mat) {
        Matrix<typename T::value_type> result(mat.colCount(), mat.rowCount());
        for (int i = 0; i < mat.rowCount(); ++i)
            std::copy(mat.row(i).begin(), mat.row(i).end(), result.col(i).begin());
        return result;
    }
}

#endif    // NUMERIX_MATRIX_HPP

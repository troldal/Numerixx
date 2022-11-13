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
        Slice() : m_start(0),
                  m_length(-1),
                  m_stride(1) {}

        explicit Slice(size_t start) : m_start(start),
                                       m_length(-1),
                                       m_stride(1) {}

        Slice(size_t start, size_t length, size_t stride = 1) : m_start(start),
                                                                m_length(length),
                                                                m_stride(stride) {}

        size_t operator()(size_t index) const {
            if (index > m_length) throw std::out_of_range("Index out of bounds.");
            return m_start + index * m_stride;
        }

        size_t start() const {return m_start;}

        size_t length() const {return m_length;}

        size_t stride() const {return m_stride;}


    private:
        size_t m_start;
        size_t m_length;
        size_t m_stride;
    };

    class MatrixExtents {
    public:
        MatrixExtents(size_t rows, size_t cols) : m_extents(std::make_pair(rows, cols)) {}

        size_t rowCount() const { return m_extents.first; }

        size_t colCount() const { return m_extents.second; }

    private:
        std::pair<size_t, size_t> m_extents;
    };

    /**
     * @brief
     */
    class MatrixSlice
    {
    public:

        MatrixSlice() = default;

        MatrixSlice(size_t start, std::initializer_list<size_t> extents) : m_size(extents.size() > 1 ? (*extents.begin()) * (*extents.end() - 1) : (*extents.begin())),
                                                                           m_start(start),
                                                                           m_extents(extents) {
            m_size = m_extents.size() > 1 ? m_extents.front() * m_extents.back() : m_extents.front();
            if (extents.size() > 2) throw std::invalid_argument("Only 2-dimensional matrices are supported.");
        }

        MatrixSlice(size_t start, std::initializer_list<size_t> extents, std::initializer_list<size_t> strides) : m_size(extents.size() > 1 ? (*extents.begin()) * (*extents.end() - 1) : (*extents.begin())),
                                                                                                                  m_start(start),
                                                                                                                  m_extents(extents),
                                                                                                                  m_strides(strides) {
            m_size = m_extents.size() > 1 ? m_extents.front() * m_extents.back() : m_extents.front();
            if (extents.size() > 2 || strides.size() > 2)  throw std::invalid_argument("Only 2-dimensional matrices are supported.");
        }


        size_t operator()(size_t row, size_t col) const
        {
            if (row > m_extents.front() - 1) throw std::invalid_argument("Invalid row number.");
            if (col > m_extents.back() - 1) throw std::invalid_argument("Invalid column number.");
            return m_start + row * m_strides.front() + col * m_strides.back();
        }

        size_t operator()(size_t index) const {
            size_t row = index / colCount();
            size_t col = index % colCount();
            return (*this)(row, col);
        }

        size_t rowCount() const {
            return m_extents.front();
        }

        size_t colCount() const {
            return m_extents.back();
        }

        size_t size() const {
            return m_size;
        }

        size_t start() const {
            return m_start;
        }

    private:
        size_t                m_size;
        size_t                m_start;
        std::vector<size_t> m_extents;
        std::vector<size_t> m_strides;
    };

    template<typename T>
    class SliceIter
    {
        T* m_data;
        MatrixSlice m_slice;
        size_t m_current;

    public:
        SliceIter(T* data, MatrixSlice slice, size_t pos = 0) : m_data(data),
                                                                m_slice(slice),
                                                                m_current(pos){}

        SliceIter end() const {
            return {m_data, m_slice, m_slice.size()};
        }

        SliceIter& operator++() {
            ++m_current;
            return *this;
        }

        SliceIter operator++(int) {
            SliceIter slice = *this;
            ++m_current;
            return slice;
        }

        T& operator*() {
            return m_data[m_slice(m_current)];
        }

        bool operator==(const SliceIter& other) const {
            return m_current == other.m_current;
        }

        bool operator!=(const SliceIter& other) const {
            return !(*this == other);
        }

        bool operator<(const SliceIter& other) const {
            return m_current < other.m_current;
        }

    };


    template<typename T>
    class SliceIterConst
    {
        const T* m_data;
        MatrixSlice m_slice;
        size_t m_current;

    public:
        SliceIterConst(const T* data, MatrixSlice slice, size_t pos = 0) : m_data(data),
                                                                m_slice(slice),
                                                                m_current(pos){}

        SliceIterConst end() const {
            return {m_data, m_slice, m_slice.size()};
        }

        SliceIterConst& operator++() {
            ++m_current;
            return *this;
        }

        SliceIterConst operator++(int) {
            SliceIter slice = *this;
            ++m_current;
            return slice;
        }

        const T& operator*() const {
            return m_data[m_slice(m_current)];
        }

        bool operator==(const SliceIterConst& other) const {
            return m_current == other.m_current;
        }

        bool operator!=(const SliceIterConst& other) const {
            return !(*this == other);
        }

        bool operator<(const SliceIterConst& other) const {
            return m_current < other.m_current;
        }

    };

    template<typename T>
        requires std::is_arithmetic_v<T> && (!std::is_same_v<T, bool>) && (!std::is_same_v<T, char>)
    class MatrixProxy;

    template<typename T>
        requires std::is_arithmetic_v<T> && (!std::is_same_v<T, bool>) && (!std::is_same_v<T, char>)
    class Matrix;

    namespace impl {

        template<typename T>
        struct MatrixTraits;



        /**
         * @brief The CRTP base class for the Matrix and MatrixProxy classes.
         * @tparam DERIVED The derived class used for CRTP
         */
        template<typename DERIVED>
        class MatrixBase {

            /*
             * Alias declarations
             */
            using value_type = typename MatrixTraits<DERIVED>::value_type;

            /**
             * @brief
             * @param row
             * @param col
             * @return
             */
            size_t index(size_t row, size_t col) const {
                const DERIVED& derived = static_cast<const DERIVED&>(*this);
                if (row >= derived.m_slice.rowCount() || col >= derived.m_slice.colCount())
                    throw std::out_of_range("Index out of bounds.");

                return derived.m_slice(row, col);
            }

            MatrixExtents extents() const {
                const DERIVED& derived = static_cast<const DERIVED&>(*this);
                return derived.extents();
            }

        public:

            /**
             * @brief
             * @param row
             * @param col
             * @return
             */
            value_type& operator()(size_t row, size_t col) {
                DERIVED& derived = static_cast<DERIVED&>(*this);
                return derived.m_data[index(row, col)];
            }

            /**
             * @brief
             * @param row
             * @param col
             * @return
             */
            value_type operator()(size_t row, size_t col) const {
                const DERIVED& derived = static_cast<const DERIVED&>(*this);

                return derived.m_data[index(row, col)];
            }

            /**
             * @brief
             * @param rowSlice
             * @param colSlice
             * @return
             */
            auto operator()(const Slice& rowSlice, const Slice& colSlice) {
                DERIVED& derived = static_cast<DERIVED&>(*this);
                return derived.slice(rowSlice, colSlice);
            }

            size_t rowCount() const {
                const DERIVED& derived = static_cast<const DERIVED&>(*this);
                return derived.m_slice.rowCount();
            }

            size_t colCount() const {
                const DERIVED& derived = static_cast<const DERIVED&>(*this);
                return derived.m_slice.colCount();
            }

            bool isSquare() const {
                return rowCount() == colCount();
            }

            SliceIter<value_type> begin() {
                DERIVED& derived = static_cast<DERIVED&>(*this);
                return SliceIter<value_type>(derived.data(), derived.m_slice);
            }

            SliceIterConst<value_type> begin() const {
                const DERIVED& derived = static_cast<const DERIVED&>(*this);
                return SliceIterConst<value_type>(derived.data(), derived.m_slice);
            }

            SliceIter<value_type> end() {
                DERIVED& derived = static_cast<DERIVED&>(*this);
                return SliceIter<value_type>(derived.data(), derived.m_slice).end();
            }

            SliceIterConst<value_type> end() const {
                const DERIVED& derived = static_cast<const DERIVED&>(*this);
                return SliceIterConst<value_type>(derived.data(), derived.m_slice).end();
            }

            template<typename T>
            requires (std::is_same_v<T, Matrix<value_type>> || std::is_same_v<T, MatrixProxy<value_type>>)
            DERIVED& operator+=(const T& other) {
                assert(other.rowCount() == rowCount() && other.colCount() == colCount());
                DERIVED& derived = static_cast<DERIVED&>(*this);
                std::transform(other.begin(), other.end(), derived.begin(), derived.begin(), std::plus<value_type>());
                return derived;
            }

            template<typename T>
                requires std::is_arithmetic_v<T> && (!std::is_same_v<T, bool>) && (!std::is_same_v<T, char>)
            DERIVED& operator+=(T value) {
                DERIVED& derived = static_cast<DERIVED&>(*this);
                for (auto& item : derived) item += value;
                return derived;
            }

        };
    }

    /**
     * @brief
     * @tparam T
     */
    template<typename T>
        requires std::is_arithmetic_v<T> && (!std::is_same_v<T, bool>) && (!std::is_same_v<T, char>)
    class MatrixProxy : public impl::MatrixBase<MatrixProxy<T>>
    {
        /*
         * Friend declarations
         */
        friend impl::MatrixBase<MatrixProxy<T>>;

        MatrixExtents extents() const {
            return m_extents;
        }

    public:

        using value_type = T;

        MatrixProxy(const MatrixExtents& extents, const MatrixSlice& slice, T* data) : m_extents(extents),
                                                                                       m_slice(slice),
                                                                                       m_data(data) {}

        auto slice(const Slice& rowSlice, const Slice& colSlice) {

            //if (rowSlice.start() + rowSlice.stride() * (rowSlice.length() - 1) >= m_slice.rowCount())
            //    throw std::out_of_range("Row index out of bounds.");

            //if (colSlice.start() + colSlice.stride() * (colSlice.length() - 1) >= m_slice.colCount())
            //    throw std::out_of_range("Column index out of bounds.");

            auto delta = MatrixExtents(extents().rowCount() - m_slice.rowCount(), extents().colCount() - m_slice.colCount());
            auto rSlice = Slice(rowSlice.start() + delta.rowCount(), rowSlice.length(), rowSlice.stride() + delta.rowCount());
            auto cSlice = Slice(colSlice.start() + delta.colCount(), colSlice.length(), colSlice.stride());


            auto start = rSlice.start() * m_slice.colCount() + cSlice.start();
            auto slice = MatrixSlice(start, {rSlice.length(), cSlice.length()}, {rSlice.stride() * m_slice.colCount(), cSlice.stride()});

            return MatrixProxy(m_extents, slice, data());
        }

    private:

        T* data() {
            return m_data;
        }

        const T* data() const {
            return m_data;
        }

        Slice m_rowSlice;
        Slice m_colSlice;
        //MatrixSlice m_slice;
        T* m_data;
        //MatrixExtents m_extents;
    };




    /**
     * @brief
     * @tparam T
     */
    template<typename T>
        requires std::is_arithmetic_v<T> && (!std::is_same_v<T, bool>) && (!std::is_same_v<T, char>)
    class Matrix : public impl::MatrixBase<Matrix<T>>
    {
        /*
         * Friend declarations
         */
        friend impl::MatrixBase<Matrix<T>>;

    private:
        std::vector<T> m_data; /**< */
        int m_rowCount; /**< */
        int m_colCount; /**< */
        //MatrixSlice m_slice; /**< */
        Slice m_rowSlice;
        Slice m_colSlice;

        /**
         * @brief
         * @param rowSlice
         * @param colSlice
         * @return
         */
        auto slice(const Slice& rowSlice, const Slice& colSlice) {

//            if (rowSlice.start() + rowSlice.stride() * (rowSlice.length() - 1) >= m_slice.rowCount())
            if (rowSlice.length() > m_slice.rowCount())
                throw std::out_of_range("Row index out of bounds.");

//            if (colSlice.start() + colSlice.stride() * (colSlice.length() - 1) >= m_slice.colCount())
            if (colSlice.length() > m_slice.colCount())
                throw std::out_of_range("Column index out of bounds.");

            auto start = rowSlice.start() * m_slice.colCount() + colSlice.start();
            auto slice = MatrixSlice(start, {rowSlice.length(), colSlice.length()}, {rowSlice.stride() * m_slice.colCount(), colSlice.stride()});

            return MatrixProxy(extents(), slice, data());
        }

        MatrixExtents extents() const {
            return MatrixExtents(m_slice.rowCount(), m_slice.colCount());
        }

    public:

        /*
         * Alias declatations
         */
        using value_type = T;

        /**
         * @brief
         * @param rows
         * @param columns
         */
        Matrix(int rows, int columns) : m_rowCount{rows},
                                        m_colCount{columns},
                                        m_data(rows * columns),
                                        //m_slice(0, {rows, columns}, {columns, 1})
                                        m_rowSlice(0, rows, columns),
                                        m_colSlice(0, columns, 1)
        {}

        /**
         * @brief
         * @param other
         */
        Matrix(const Matrix& other) = default;

        /**
         * @brief
         * @param other
         */
        Matrix(Matrix&& other) noexcept = default;

        /**
         * @brief
         */
        ~Matrix() = default;

        /**
         * @brief
         * @return
         */
        Matrix& operator=(const Matrix&) = default;

        /**
         * @brief
         * @param other
         * @return
         */
        Matrix& operator=(Matrix&& other) noexcept = default;

        /**
         * @brief
         * @param i
         * @return
         */
        auto operator[](const int i) {
            return std::span (&m_data.data()[i * m_colSlice.length()], m_colSlice.length());
        }

        /**
         * @brief
         * @param i
         * @return
         */
        const auto operator[](int i) const {
            return std::span (&m_data.data()[i * m_colSlice.length()], m_colSlice.length());
        }

        /**
         * @brief
         * @param i
         * @return
         */
        const auto row(int i) const {
            return (*this)[i];
        }

        /**
         * @brief
         * @param i
         * @return
         */
        auto row(int i) {
            return (*this)[i];
        }

        /**
         * @brief
         * @param i
         * @return
         */
        auto col(int i) const {

            std::vector<T> result;
            for (int row = 0; row < m_rowSlice.length(); ++row)
                result.push_back(m_data[row * m_colSlice.length() + i]);

            return result;
        }

        /**
         * @brief
         * @return
         */
        T* data() {
            return m_data.data();
        }

        /**
         * @brief
         * @return
         */
        const T* data() const {
            return m_data.data();
        }

        /**
         * @brief
         * @param other
         * @return
         */
        Matrix& operator-=(const Matrix& other) {
            assert(other.rowCount() == m_rowSlice.length() && other.colCount() == m_colSlice.length());
            std::transform(other.m_data.begin(), other.m_data.end(), m_data.begin(), m_data.begin(), std::minus<T>());
            return *this;
        }

        /**
         * @brief
         * @tparam U
         * @param value
         * @return
         */
        template<typename U>
        requires std::is_arithmetic_v<T> && (!std::is_same_v<T, bool>) && (!std::is_same_v<T, char>)
        Matrix& operator-=(U value) {
            std::transform(m_data.begin(), m_data.end(), m_data.begin(), [&](T elem) {return elem - value;} );
            return *this;
        }

        /**
         * @brief
         * @tparam U
         * @param value
         * @return
         */
        template<typename U>
        requires std::is_arithmetic_v<T> && (!std::is_same_v<T, bool>) && (!std::is_same_v<T, char>)
        Matrix& operator*=(U value) {
            std::transform(m_data.begin(), m_data.end(), m_data.begin(), [&](T elem) {return elem * value;} );
            return *this;
        }

        /**
         * @brief
         * @tparam U
         * @param value
         * @return
         */
        template<typename U>
        requires std::is_arithmetic_v<T> && (!std::is_same_v<T, bool>) && (!std::is_same_v<T, char>)
        Matrix& operator/=(U value) {
            std::transform(m_data.begin(), m_data.end(), m_data.begin(), [&](T elem) {return elem / value;} );
            return *this;
        }

        Matrix& operator*=(const Matrix& other) {
            assert(m_colSlice.length() == other.rowCount());
            
            std::vector<T> result;
            for (int r = 0; r < m_rowSlice.length(); ++r)
                for (int c = 0; c < other.colCount(); ++c)
                    result.push_back(std::inner_product(row(r).begin(), row(r).end(), other.col(c).begin(), 0.0));

            m_data = result;
            m_colCount = other.colCount();
            return *this;
        }

        void augment(Matrix<T>& mat) {
            assert(mat.colCount() == 1 && mat.rowCount() == m_rowSlice.length());

            //Matrix<T> result (rowCount(), colCount() + 1);


            int i = 0;
            for (auto& r : m_data) {
                r.push_back(mat[i][0]);
                ++i;
            }
            m_colCount++;
        }
    };

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
        struct MatrixTraits<MatrixProxy<T>> {
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
    requires std::is_same_v<T, Matrix<typename T::value_type>> || std::is_same_v<T, MatrixProxy<typename T::value_type>>
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
    template<typename T>
    inline Matrix<T> operator+(const Matrix<T>& a, const Matrix<T>& b) {

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
    inline Matrix<T> operator+(const Matrix<T>& a, U b) {

        auto result = a;
        result += b;
        return result;
    }

    /**
     * @brief
     * @tparam T
     * @param a
     * @param b
     * @return
     */
    template<typename T>
    inline Matrix<T> operator-(const Matrix<T>& a, const Matrix<T>& b) {

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
    inline Matrix<T> operator-(const Matrix<T>& a, U b) {

        auto result = a;
        result -= b;
        return result;
    }

    /**
     * @brief
     * @tparam T
     * @param a
     * @param b
     * @return
     */
    template<typename T>
    inline Matrix<T> operator*(const Matrix<T>& a, const Matrix<T>& b) {

        auto result = a;
        result *= b;
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
    inline Matrix<T> operator*(const Matrix<T>& a, U b) {

        auto result = a;
        result *= b;
        return result;
    }

    /**
     * @brief
     * @tparam T
     * @param a
     * @param b
     * @return
     */
    template<typename T>
    inline Matrix<T> operator/(const Matrix<T>& a, const T& b) {

        auto result = a;
        result /= b;
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
    inline Matrix<T> operator/(const Matrix<T>& a, U b) {

        auto result = a;
        result /= b;
        return result;
    }

    template<typename T>
    inline Matrix<T> transpose(const Matrix<T>& mat) {
        Matrix<T> result(mat.colCount(), mat.rowCount());

        for (int r = 0; r < result.rowCount(); ++r)
            for (int c = 0; c < result.colCount(); ++c)
                result[r][c] = mat[c][r];

        return result;
    }

}





// MATRIX CLASS FROM STROUSTRUP TCPL 4th EDITION

//#include <algorithm>
//#include <array>
//#include <cassert>
//#include <numeric>
//#include <type_traits>
//#include <vector>
//
//namespace matrix_impl {
//
//    template<typename T, size_t N>
//    struct matrix_init {
//        using type = std::initializer_list<typename matrix_init<T,N-1>::type>;
//    };
//
//    template<typename T>
//    struct matrix_init<T,1> {
//        using type = std::initializer_list<T>;
//    };
//
//    template<typename T>
//    struct matrix_init<T,0>;
//
//    template<size_t N, typename LIST>
//    std::array<size_t, N> derive_extents(const LIST& list);
//
//    template<size_t N, typename I, typename LIST>
//    std::enable_if<(N>1),void> add_extents(I& first, const LIST& list);
//
//    template<size_t N, typename I, typename LIST>
//    std::enable_if<(N==1),void> add_extents(I& first, const LIST& list);
//
//    template<typename LIST>
//    bool check_non_jagged(const LIST& list);
//
//    template<typename T, typename VEC>
//    void insert_flat(std::initializer_list<T> list, VEC& vec);
//
//    template<typename T, typename VEC>
//    void add_list(const std::initializer_list<T>* first, const std::initializer_list<T>* last, VEC& vec);
//
//    template<typename T, typename VEC>
//    void add_list(const T* first, const T* last, VEC& vec);
//
//}
//
//template<size_t N>
//struct matrix_slice {
//
//    matrix_slice() = default;
//
//    matrix_slice(size_t s, std::initializer_list<size_t> exts);
//    matrix_slice(size_t s, std::initializer_list<size_t> exts, std::initializer_list<size_t> strs);
//
//    template<typename... DIMS>
//    matrix_slice(DIMS... dims) {
//        static_assert(sizeof...(DIMS)==N,"");
//
//        size_t args[N] {size_t(dims)...};
//        return std::inner_product(args, args+N, strides.begin(), size_t(0));
//    }
//
//    template<typename... DIMS, typename = std::enable_if<All(Convertible<DIMS, size_t>()...)> >
//    size_t operator()(DIMS... dims) const;
//
//    size_t size;
//    size_t start;
//    std::array<size_t,N> extents;
//    std::array<size_t,N> strides;
//};
//
//template<typename T, size_t N>
//class matrix_ref {
//
//public:
//    matrix_ref(const matrix_slice<N>& s, T* p) : desc{s}, ptr{p} {}
//
//
//private:
//    matrix_slice<N> desc;
//    T* ptr;
//
//};
//
//template<typename T, size_t N>
//using matrix_initializer = typename matrix_impl::matrix_init<T,N>::type;
//
//
//namespace matrix_impl {
//
//    template<size_t N, typename... DIMS>
//    bool check_bounds(const matrix_slice<N>& slice, DIMS... dims);
//
//    constexpr bool All() {return true;}
//
//    template<typename... Args>
//    constexpr bool All(bool b, Args... args) {
//        return b && All(args...);
//    }
//
//    template<typename... Args>
//    constexpr bool requesting_element() {
//        return All(std::is_convertible_v<Args, size_t>()...);
//    }
//
//    template<typename... Args>
//    constexpr bool requesting_slice() {
//        return All((std::is_convertible_v<Args, size_t>() || std::is_same_v<Args, slice>())...)
//        // WHAT!!!???
//    }
//
//}
//
//template<typename T, size_t N>
//class matrix {
//
//public:
//    //static constexpr size_t order = N;
//    using value_type = T;
//    using iterator = typename std::vector<T>::iterator;
//    using const_iterator = typename std::vector<T>::const_iterator;
//
//    matrix() = default;
//    matrix(matrix&& other) noexcept = default;
//    matrix& operator=(matrix&& other) noexcept = default;
//    matrix(const matrix& other) = default;
//    matrix& operator=(const matrix& other) = default;
//    ~matrix() = default;
//
//    template<typename U>
//    matrix(const matrix_ref<U,N>& other) : desc{other.desc},
//                                           elems{other.begin(), other.end()}
//    {
//        static_assert(std::is_convertible_v<U,T>(), "matrix constructor: incompatibel element types");
//    }
//
//    template<typename U>
//    matrix& operator=(const matrix_ref<U, N>& other) {
//        static_assert(std::is_convertible_v<U,T>(), "matrix =: incompatibel element types");
//
//        desc = other.desc;
//        elems.assign(other.begin(), other.end());
//        return *this;
//    }
//
//    template<typename...EXTS>
//    explicit matrix(EXTS... exts) : desc{exts...},
//                                    elems{desc.size}
//    {}
//
//    matrix(matrix_initializer<T,N> init) {
//        matrix_impl::derive_extents(init, desc.extents);
//        elems.reserve(desc.size);
//        matrix_impl::insert_flat(init, elems);
//        assert(elems.size() == desc.size);
//    }
//
//    matrix& operator=(matrix_initializer<T,N> list);
//
//    template<typename U>
//    matrix(std::initializer_list<U>) = delete;
//
//    template<typename U>
//    matrix& operator=(std::initializer_list<U>) = delete;
//
//    static constexpr size_t order() {return N;}
//    size_t extent(size_t n) const {return desc.extents[n];}
//    size_t size() const {return elems.size();}
//    const matrix_slice<N>& descriptor() const {return desc;}
//
//    T* data() {return elems.data();}
//    const T* data() const {return elems.data();}
//
//    template<typename... Args>
//    std::enable_if<matrix_impl::requesting_element<Args...>(), T&> operator()(Args... args) {
//        assert(matrix_impl::check_bounds(desc, args...));
//        return *(data() + desc(args...));
//    }
//
//    template<typename... Args>
//    std::enable_if<matrix_impl::requesting_element<Args...>(), const T&> operator()(Args... args) const {
//        assert(matrix_impl::check_bounds(desc, args...));
//        return *(data() + desc(args...));
//    }
//
//    template<typename... Args>
//    std::enable_if<matrix_impl::requesting_element<Args...>(), matrix_ref<T,N>>
//        operator()(const Args&... args);
//
//    template<typename... Args>
//    std::enable_if<matrix_impl::requesting_element<Args...>(), matrix_ref<const T,N>>
//        operator()(const Args&... args) const;
//
//    matrix_ref<T,N-1> operator[](size_t i) {return row(i); }
//    matrix_ref<const T, N-1> operator[](size_t i) const {return row(i); }
//
//    matrix_ref<T,N-1> row(size_t n) {
//        assert(n<rows());
//        matrix_slice<N-1> row;
//        matrix_impl::slice_dim<0>(n, desc, row);
//        return {row, data()};
//    }
//
//    matrix_ref<const T,N-1> row(size_t n) const {
//        assert(n<rows());
//        matrix_slice<N-1> row;
//        matrix_impl::slice_dim<0>(n, desc, row);
//        return {row, data()};
//    }
//
//    // Specializations for  N == 1 and N == 0 missing
//
//    matrix_ref<T,N-1> col(size_t n) {
//        assert(n<cols());
//        matrix_slice<N-1> col;
//        matrix_impl::slice_dim<1>(n, desc, col);
//        return {col, data()};
//    }
//
//    matrix_ref<const T,N-1> col(size_t n) const {
//        assert(n<cols());
//        matrix_slice<N-1> col;
//        matrix_impl::slice_dim<1>(n, desc, col);
//        return {col, data()};
//    }
//
//    // Specializations for  N == 1 and N == 0 missing
//
//
//    template<typename F>
//    matrix& apply(F f) {
//        for (auto& x : elems) f(x);
//        return *this;
//    }
//
//    template<typename M, typename F>
//    matrix& apply(const M& m, F f) {
//        assert(same_extents(desc, m.descriptor()));
//        for(auto i = begin(), j = m.begin(); i!=end(); ++i, ++j)
//            f(*i, *j);
//
//        return *this;
//    }
//
//    matrix& operator=(const T& value);
//
//    matrix& operator+=(const T& value) {
//        return apply([&](T& a){a+=value;});
//    }
//
//    matrix& operator-=(const T& value) {
//        return apply([&](T& a){a-=value;});
//    }
//
//    matrix& operator*=(const T& value) {
//        return apply([&](T& a){a*=value;});
//    }
//
//    matrix& operator/=(const T& value) {
//        return apply([&](T& a){a/=value;});
//    }
//
//    matrix& operator%=(const T& value) {
//        return apply([&](T& a){a%=value;});
//    }
//
//    template<typename M>
//    std::enable_if<matrix_type<M>, matrix<T,N>&> operator+=(const M& x) {
//        static_assert(m.order == N, "+=: mismatched matrix dimensions");
//        assert(same_extents(desc, m.descriptor()));
//
//        return apply(m, [](T& a, value_type<M>&b){a+=b;});
//    }
//
//    template<typename M>
//    std::enable_if<matrix_type<M>, matrix<T,N>&> operator-=(const M& x) {
//        static_assert(m.order == N, "-=: mismatched matrix dimensions");
//        assert(same extents(desc, m.descriptor()));
//
//        return apply(m, [](T& a, value_type<M>&b){a-=b;});
//    }
//
//
//private:
//    matrix_slice<N> desc;
//    std::vector<T> elems;
//
//
//};
//
//namespace matrix_impl
//{
//
//    template<typename T, size_t N>
//    inline matrix<T, N> operator+(const matrix<T, N>& m, const T& val)
//    {
//        matrix<T, N> res = m;
//        res += val;
//        return res;
//    }
//
//    template<typename T,
//             typename T2,
//             size_t                                                                                        N,
//             typename RT = matrix<std::common_type<value_type<T>, value_type<T2>>, N> inline matrix<RT, N> operator+(const matrix<T, N>&  a,
//                                                                                                                     const matrix<T2, N>& b)
//    {
//        matrix<RT, N> res = a;
//        res += b;
//        return res;
//    }
//
//    template<typename T, size_t N>
//    inline matrix<T, N> operator+(const matrix_ref<T, N>& x, const T& n)
//    {
//        matrix<T, N> res = x;
//        res += n;
//        return res;
//    }
//
//    template<typename T>
//    inline matrix<T, 2> operator*(const matrix<T, 1>& u, const matrix<T, 1>& v)
//    {
//        const size_t n = u.extent(0);
//        const size_t m = v.extent(0);
//        matrix<T, 2> res(n, m);
//
//        for (size_t i = 0; i != n; ++i)
//            for (size_t j = 0; j != m; ++j) res(i, j) = u[i] * v[j];
//
//        return res;
//    }
//
//    template<typename T>
//    inline matrix<T, 1> operator*(const matrix<T, 2>& m, const matrix<T, 1>& v)
//    {
//        assert(m.extent(1) == v.extent(0));
//
//        const size_t n = m.extent(0);
//        matrix<T, 1> res(n);
//        for (size_t i = 0; i != n; ++i)
//            for (size_t j = 0; j != n; ++j) res(i) += m(i, j) * v(j);
//    }
//
//    template<typename T>
//    inline matrix<T, 2> operator*(const matrix<T, 2>& m1, const matrix<T, 2>& m2)
//    {
//        const size_t n = m1.extent(0);
//        const size_t m = m1.extent(1);
//        assert(m == m2.extent(0));
//
//        const size_t p = m2.extent(1);
//        matrix<T, 2> res(n, p);
//        for (size_t i = 0; i != n; ++i)
//            for (size_t j = 0; j != m; ++j)
//                for (size_t k = 0; k != p; ++k) res(i, j) = m1(i, k) * m2(k, j);
//
//        return res;
//    }
//
//    template<size_t N, typename LIST>
//    std::array<size_t, N> derive_extents(const LIST& list)
//    {
//        std::array<size_t, N> a;
//        auto f = a.begin();
//        add_extents<N>(f, list);
//        return a;
//    }
//
//    template<size_t N, typename I, typename LIST>
//    std::enable_if<(N > 1), void> add_extents(I& first, const LIST& list)
//    {
//        assert(check_non_jagged(list));
//        *first = list.size();
//        add_extents<N-1>(++first, *list.begin());
//    }
//
//    template<size_t N, typename I, typename LIST>
//    std::enable_if<(N == 1), void> add_extents(I& first, const LIST& list)
//    {
//        *first++ = list.size();
//    }
//
//    template<typename LIST>
//    bool check_non_jagged(const LIST& list)
//    {
//        auto i = list.begin();
//        for (auto j = i+1; j!=list.end(); ++j)
//            if (i->size() != j->size())
//                return false;
//
//        return true;
//    }
//
//    template<typename T, typename VEC>
//    void insert_flat(std::initializer_list<T> list, VEC& vec)
//    {
//        add_list(list.begin(), list.end(), vec);
//    }
//
//    template<typename T, typename VEC>
//    void add_list(const std::initializer_list<T>* first, const std::initializer_list<T>* last, VEC& vec)
//    {
//        for (;first!=last; ++first)
//            add_list(first->begin(), first->end(), vec);
//    }
//
//    template<typename T, typename VEC>
//    void add_list(const T* first, const T* last, VEC& vec) {
//        vec.insert(vec.end(), first, last);
//    }
//
//    template<size_t N, typename... DIMS>
//    bool check_bounds(const matrix_slice<N>& slice, DIMS... dims)
//    {
//        size_t indexes[N]{size_t(dims...)};
//        return std::equal(indexes, indexes+N, slice.extents, std::less<size_t>{});
//    }
//}
//
#endif    // NUMERIX_MATRIX_HPP

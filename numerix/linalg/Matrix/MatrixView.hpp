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

#ifndef NUMERIX_MATRIXVIEW_HPP
#define NUMERIX_MATRIXVIEW_HPP

#include "MatrixCommon.h"

namespace numerix::linalg
{

    /**
     * @brief The MatrixViewConcept class is a genereic implementation of a view into a subset of the elements of a Matrix.
     * Similar to the Matrix class, it derives from the MatrixBase class using CRTP. Specializations of this class are the
     * MatrixView and MatrixViewConst classes, for mutable and const views, respectively.
     * An object of the MatrixProxy class does not own the data; the data is owned by the corresponding Matrix object.
     * @tparam T The type of the underlying matrix elements. Default is double, but can be any kind of floating point or integer number.
     * @tparam IsConst
     * @todo Consider putting in the impl namespace, as MatrixProxy objects should only be created from a Matrix object.
     */
    template<typename T, bool IsConst>
        requires is_number<T>
    class MatrixViewConcept : public impl::MatrixBase<MatrixViewConcept<T, IsConst>>
    {
        /*
         * Friend declarations. Necessary for the base class to access elements of MatrixProxy using CRTP.
         */
        friend impl::MatrixBase<MatrixViewConcept<T, IsConst>>;

        /**
         * Private alias declatations.
         */
        using parent   = impl::MatrixBase<MatrixViewConcept<T, IsConst>>;
        using matrix_t = std::conditional_t<IsConst, const Matrix<T>, Matrix<T>>;

    private:
        Slice     m_rowSlice; /**< */
        Slice     m_colSlice; /**< */
        matrix_t* m_matrix;   /**< */

        /**
         * @brief
         * @return
         */
        T* data()
            requires(!std::is_const_v<T>);

        /**
         * @brief
         * @return
         */
        const T* data() const;

        /**
         * @brief
         * @return
         */
        auto gslice() const;

    public:
        using value_type = T;

        /**
         * @brief
         * @param rowSlice
         * @param colSlice
         * @param data
         */
        MatrixViewConcept(const Slice& rowSlice, const Slice& colSlice, matrix_t* matrix)
            : m_rowSlice(rowSlice),
              m_colSlice(colSlice),
              m_matrix(matrix)
        {}

        /**
         * @brief
         * @param other
         */
        MatrixViewConcept(const MatrixViewConcept& other) = default;

        /**
         * @brief
         * @param other
         */
        MatrixViewConcept(MatrixViewConcept&& other) noexcept = default;

        /**
         * @brief
         */
        ~MatrixViewConcept() = default;

        /**
         * @brief
         * @param other
         * @return
         */
        MatrixViewConcept& operator=(const MatrixViewConcept& other) = default;

        /**
         * @brief
         * @param other
         * @return
         */
        MatrixViewConcept& operator=(MatrixViewConcept&& other) noexcept = default;

        /**
         * @brief
         * @return
         */
        auto extents() const;
    };

    /**
     * @details
     */
    template<typename T, bool IsConst>
        requires is_number<T>
    T* MatrixViewConcept<T, IsConst>::data()
        requires(!std::is_const_v<T>)
    {
        return m_matrix->data();
    }

    /**
     * @details
     */
    template<typename T, bool IsConst>
        requires is_number<T>
    const T* MatrixViewConcept<T, IsConst>::data() const
    {
        return m_matrix->data();
    }

    /**
     * @details
     */
    template<typename T, bool IsConst>
        requires is_number<T>
    auto MatrixViewConcept<T, IsConst>::gslice() const
    {
        auto start = m_rowSlice.start() * m_matrix->extents().first + m_colSlice.start();
        return GSlice(start, { m_rowSlice.length(), m_colSlice.length() }, { m_rowSlice.stride(), m_colSlice.stride() });
    }

    /**
     * @details
     */
    template<typename T, bool IsConst>
        requires is_number<T>
    auto MatrixViewConcept<T, IsConst>::extents() const
    {
        return m_matrix->extents();
    }

}    // namespace numerix::linalg


#endif    // NUMERIX_MATRIXVIEW_HPP

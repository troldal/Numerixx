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

#include "MatrixCommon.hpp"

namespace numerix::linalg
{

    /**
     * @brief The MatrixViewConcept class is a generic implementation of a view into a subset of the elements of a Matrix.
     * Similar to the Matrix class, it derives from the MatrixBase class using CRTP. Specializations of this class are the
     * MatrixView and MatrixViewConst classes, for mutable and const views, respectively.
     * An object of the MatrixProxy class does not own the data; the data is owned by the corresponding Matrix object.
     * @tparam T The type of the underlying matrix elements. Default is double, but can be any kind of floating point or integer number.
     * @tparam IsConst A flag indicating if the object should be treated as const (i.e. non-mutable, similar to a const_iterator).
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
        friend impl::MatrixBase<Matrix<T>>;

        /**
         * Private alias declatations.
         */
        using parent   = impl::MatrixBase<MatrixViewConcept<T, IsConst>>;
        using matrix_t = std::conditional_t<IsConst, const Matrix<T>, Matrix<T>>;

    private:
        matrix_t* m_matrix;   /**< A raw pointer (i.e. non-owning) to the underlying Matrix object. */
        Slice     m_rowSlice; /**< The Slice describing the rows. */
        Slice     m_colSlice; /**< The Slice describing the columns. */

        /**
         * @brief Access the raw array of Matrix elements.
         * @return A pointer to the first element.
         * @note Unlike in the Matrix class, this is a private member. The reason is that it provides access to
         * the entire array of matrix elements, not only the ones visible in the view. Direct access by clients
         * would result in unexpected behaviour.
         * @note The non-const version is only available for mutable types.
         */
        T* data()
            requires(!IsConst) && (!std::is_const_v<T>)
        {
            return m_matrix->data();
        }

        /**
         * @brief Access the raw array of Matrix elements.
         * @return A pointer to the first element.
         * @note Unlike in the Matrix class, this is a private member. The reason is that it provides access to
         * the entire array of matrix elements, not only the ones visible in the view. Direct access by clients
         * would result in unexpected behaviour.
         */
        const T* data() const { return m_matrix->data(); }

        /**
         * @brief
         * @return
         */
        auto gslice() const
        {
            auto start = m_rowSlice.start() * m_matrix->extents().first + m_colSlice.start();
            return impl::GSlice(start, { m_rowSlice.length(), m_colSlice.length() }, { m_rowSlice.stride(), m_colSlice.stride() });
        }

        /**
         * @brief Get the extents in each of the dimensions (rows and columns), of the parent Matrix.
         * @return A std::pair with the row and column extents.
         */
        auto extents() const { return m_matrix->extents(); }

        /**
         * @brief Constructor taking row/column slices and the parent Matrix object as arguments.
         * @param rowSlice Slice object representing the rows.
         * @param colSlice Slice object representing the columns.
         * @param mat The parent Matrix object
         * @note This constructor is private, as clients should not be allowed to create MatrixViews directly.
         */
        MatrixViewConcept(const Slice& rowSlice, const Slice& colSlice, matrix_t* mat)
            : m_matrix(mat),
              m_rowSlice(rowSlice),
              m_colSlice(colSlice)

        {}

    public:
        /**
         * Public alias declatations. To be consistant with standard library containers, and to provide access to
         * non-standard assignment operator.
         */
        using value_type = T;
        using parent::operator=;

        /**
         * @brief Copy constructor.
         * @param other Object to be copied.
         */
        MatrixViewConcept(const MatrixViewConcept& other) = default;

        /**
         * @brief Move constructor.
         * @param other Object to be moved.
         */
        MatrixViewConcept(MatrixViewConcept&& other) noexcept = default;

        /**
         * @brief Destructor.
         */
        ~MatrixViewConcept() = default;

        /**
         * @brief Copy assignment operator.
         * @param other Object to be copied.
         * @return The copied-to object.
         */
        MatrixViewConcept& operator=(const MatrixViewConcept& other) = default;

        /**
         * @brief Move assignment operator.
         * @param other Object to be moved.
         * @return The moved-to object.
         */
        MatrixViewConcept& operator=(MatrixViewConcept&& other) noexcept = default;
    };
}    // namespace numerix::linalg

#endif    // NUMERIX_MATRIXVIEW_HPP

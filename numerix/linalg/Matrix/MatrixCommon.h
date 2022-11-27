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

#ifndef NUMERIX_MATRIXCOMMON_H
#define NUMERIX_MATRIXCOMMON_H

#include <cassert>
#include <cstddef>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace numerix::linalg
{

    /*
     *
     */
    template<typename T>
    concept is_number = (std::integral<T> || std::floating_point<T>) && (!std::same_as<T, bool>) && (!std::same_as<T, char>);


    /*
     * Forward declaration of the Matrix class.
     */
    template<typename T>
        requires is_number<T>
    class Matrix;

    /*
     * Forward declaration of the MatrixViewConcept class.
     */
    template<typename T, bool IsConst>
        requires is_number<T>
    class MatrixViewConcept;

    /**
     * @brief Declaration of the spcialization of the (mutable) MatrixView class, based on the MatrixViewConcept class.
     */
    template<typename T>
    using MatrixView = MatrixViewConcept<T, false>;

    /**
     * @brief Declaration of the spcialization of the (const) MatrixViewConst class, based on the MatrixViewConcept class.
     */
    template<typename T>
    using MatrixViewConst = MatrixViewConcept<T, true>;


    /**
     * @brief
     * @tparam T
     * @tparam IsConst
     */
    template<typename T, bool IsConst>
    class MatrixElementIterConcept;

    /**
     * @brief
     */
    template<typename T>
    using MatrixElementIter = MatrixElementIterConcept<T, false>;

    /**
     * @brief
     */
    template<typename T>
    using MatrixElementIterConst = MatrixElementIterConcept<T, true>;

    /**
     * @brief The concept of a matrix. Can be a Matrix, MatrixView or MatrixViewConst object.
     * @tparam T The type of the matrix object.
     */
    template<typename T>
    concept is_matrix = std::same_as<T, Matrix<typename T::value_type>> || std::same_as<T, MatrixView<typename T::value_type>> ||
                        std::same_as<T, MatrixViewConst<typename T::value_type>>;

    namespace impl
    {
        /*
         * Forward declaration of the MatrixTraits class.
         */
        template<typename T>
        struct MatrixTraits;

        template<typename T>
        struct MatrixTraits<Matrix<T>>
        {
            using value_type = T;
        };

        template<typename T>
        struct MatrixTraits<MatrixView<T>>
        {
            using value_type = T;
        };

        template<typename T>
        struct MatrixTraits<MatrixViewConst<T>>
        {
            using value_type = T;
        };

    }    // namespace impl

    /*
     *
     */
    template<typename T, bool IsConst>
        requires std::same_as<T, MatrixView<typename T::value_type>> || std::same_as<T, MatrixViewConst<typename T::value_type>>
    class MatrixColsConcept;

    /**
     * @brief
     */
    template<typename T>
    using MatrixCols = MatrixColsConcept<T, false>;

    /**
     * @brief
     */
    template<typename T>
    using MatrixColsConst = MatrixColsConcept<T, true>;

    /*
     *
     */
    template<typename T, bool IsConst>
        requires std::same_as<T, MatrixView<typename T::value_type>> || std::same_as<T, MatrixViewConst<typename T::value_type>>
    class MatrixRowsConcept;

    /**
     * @brief
     */
    template<typename T>
    using MatrixRows = MatrixRowsConcept<T, false>;

    /**
     * @brief
     */
    template<typename T>
    using MatrixRowsConst = MatrixRowsConcept<T, true>;

    /*
     *
     */
    template<typename T, bool IsConst>
        requires std::same_as<T, MatrixView<typename impl::MatrixTraits<T>::value_type>> || std::same_as<T, MatrixViewConst<typename impl::MatrixTraits<T>::value_type>>
    class MatrixElementsConcept;

    /**
     * @brief
     */
    template<typename T>
    using MatrixElements = MatrixElementsConcept<T, false>;

    /**
     * @brief
     */
    template<typename T>
    using MatrixElementsConst = MatrixElementsConcept<T, true>;

    namespace impl
    {

        template<typename T>
        struct MatrixTraits<MatrixElements<T>>
        {
            using value_type = T;
        };

        template<typename T>
        struct MatrixTraits<MatrixElementsConst<T>>
        {
            using value_type = T;
        };


    }    // namespace impl

    /*
     *
     */
    template<typename T, bool IsConst>
        requires std::same_as<T, MatrixCols<typename T::value_type>> || std::same_as<T, MatrixColsConst<typename T::value_type>>
    class MatrixColIterConcept;

    /**
     * @brief
     */
    template<typename T>
    using MatrixColIter = MatrixColIterConcept<T, false>;

    /**
     * @brief
     */
    template<typename T>
    using MatrixColIterConst = MatrixColIterConcept<T, true>;

    /*
     *
     */
    template<typename T, bool IsConst>
        requires std::same_as<T, MatrixRows<typename T::value_type>> || std::same_as<T, MatrixRowsConst<typename T::value_type>>
    class MatrixRowIterConcept;

    /**
     * @brief
     */
    template<typename T>
    using MatrixRowIter = MatrixRowIterConcept<T, false>;

    /**
     * @brief
     */
    template<typename T>
    using MatrixRowIterConst = MatrixRowIterConcept<T, true>;



}

template<typename T>
    requires std::same_as<T, numerix::linalg::MatrixElements<typename numerix::linalg::impl::MatrixTraits<T>::value_type>> ||
             std::same_as<T, numerix::linalg::MatrixElementsConst<typename numerix::linalg::impl::MatrixTraits<T>::value_type>>
inline void swap(T lhs, T rhs);

#endif    // NUMERIX_MATRIXCOMMON_H

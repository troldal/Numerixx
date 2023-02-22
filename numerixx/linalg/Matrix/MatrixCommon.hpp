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

#ifndef NUMERIXX_MATRIXCOMMON_HPP
#define NUMERIXX_MATRIXCOMMON_HPP

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <vector>
#include <iterator>

namespace nxx::linalg
{

    /**
     * @brief Concept defining an integral or floating-point number to be used in the Matrix classes
     * @tparam T
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
     * @brief Declaration of the specialization of the (mutable) MatrixView class, based on the MatrixViewConcept class.
     */
    template<typename T>
    using MatrixView = MatrixViewConcept<T, false>;

    /**
     * @brief Declaration of the specialization of the (const) MatrixViewConst class, based on the MatrixViewConcept class.
     */
    template<typename T>
    using MatrixViewConst = MatrixViewConcept<T, true>;

    /*
     * Forward declaration of the MatrixElementIterConcept class.
     */
    template<typename T, bool IsConst>
    class MatrixElemIterConcept;

    /**
     * @brief Declaration of the specialization of the (mutable) MatrixElementIter class, based on the MatrixElementIterConcept class.
     */
    template<typename T>
    using MatrixElemIter = MatrixElemIterConcept<T, false>;

    /**
     * @brief Declaration of the specialization of the (const) MatrixElementIterConst class, based on the MatrixElementIterConcept class.
     */
    template<typename T>
    using MatrixElemIterConst = MatrixElemIterConcept<T, true>;

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

        /*
         * Specialization of the MatrixTraits class for Matrix<T>
         */
        template<typename T>
        struct MatrixTraits<Matrix<T>>
        {
            using value_type = T;
        };

        /*
         * Specialization of the MatrixTraits class for MatrixView<T>
         */
        template<typename T>
        struct MatrixTraits<MatrixView<T>>
        {
            using value_type = T;
        };

        /*
         * Specialization of the MatrixTraits class for MatrixViewConst<T>
         */
        template<typename T>
        struct MatrixTraits<MatrixViewConst<T>>
        {
            using value_type = T;
        };

    }    // namespace impl

    /*
     * Forward declaration of the MatrixColsConcept class.
     */
    template<typename T, bool IsConst>
        requires std::same_as<T, MatrixView<typename T::value_type>> || std::same_as<T, MatrixViewConst<typename T::value_type>>
    class MatrixColsConcept;

    /**
     * @brief Declaration of the specialization of the (mutable) MatrixCols class, based on the MatrixColsConcept class.
     */
    template<typename T>
    using MatrixCols = MatrixColsConcept<T, false>;

    /**
     * @brief Declaration of the specialization of the (const) MatrixColsConst class, based on the MatrixColsConcept class.
     */
    template<typename T>
    using MatrixColsConst = MatrixColsConcept<T, true>;

    /*
     * Forward declaration of the MatrixColsConcept class.
     */
    template<typename T, bool IsConst>
        requires std::same_as<T, MatrixView<typename T::value_type>> || std::same_as<T, MatrixViewConst<typename T::value_type>>
    class MatrixRowsConcept;

    /**
     * @brief Declaration of the specialization of the (mutable) MatrixRows class, based on the MatrixRowsConcept class.
     */
    template<typename T>
    using MatrixRows = MatrixRowsConcept<T, false>;

    /**
     * @brief Declaration of the specialization of the (const) MatrixRowsConst class, based on the MatrixRowsConcept class.
     */
    template<typename T>
    using MatrixRowsConst = MatrixRowsConcept<T, true>;

    /*
     * Forward declaration of the MatrixColIterConcept class.
     */
    template<typename T, bool IsConst>
        requires std::same_as<T, MatrixCols<typename T::matrix_type>> || std::same_as<T, MatrixColsConst<typename T::matrix_type>>
    class MatrixColIterConcept;

    /**
     * @brief Declaration of the specialization of the (mutable) MatrixColIter class, based on the MatrixColIterConcept class.
     */
    template<typename T>
    using MatrixColIter = MatrixColIterConcept<T, false>;

    /**
     * @brief Declaration of the specialization of the (const) MatrixColIterConst class, based on the MatrixColIterConcept class.
     */
    template<typename T>
    using MatrixColIterConst = MatrixColIterConcept<T, true>;

    /*
     * Forward declaration of the MatrixRowIterConcept class.
     */
    template<typename T, bool IsConst>
        requires std::same_as<T, MatrixRows<typename T::matrix_type>> || std::same_as<T, MatrixRowsConst<typename T::matrix_type>>
    class MatrixRowIterConcept;

    /**
     * @brief Declaration of the specialization of the (mutable) MatrixRowIter class, based on the MatrixRowIterConcept class.
     */
    template<typename T>
    using MatrixRowIter = MatrixRowIterConcept<T, false>;

    /**
     * @brief Declaration of the specialization of the (const) MatrixRowIterConst class, based on the MatrixRowIterConcept class.
     */
    template<typename T>
    using MatrixRowIterConst = MatrixRowIterConcept<T, true>;
}    // namespace numerix::linalg

#endif    // NUMERIXX_MATRIXCOMMON_HPP

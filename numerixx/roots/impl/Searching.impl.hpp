/*
    888b      88  88        88  88b           d88  88888888888  88888888ba   88  8b        d8  8b        d8
    8888b     88  88        88  888b         d888  88           88      "8b  88   Y8,    ,8P    Y8,    ,8P
    88 `8b    88  88        88  88`8b       d8'88  88           88      ,8P  88    `8b  d8'      `8b  d8'
    88  `8b   88  88        88  88 `8b     d8' 88  88aaaaa      88aaaaaa8P'  88      Y88P          Y88P
    88   `8b  88  88        88  88  `8b   d8'  88  88"""""      88""""88'    88      d88b          d88b
    88    `8b 88  88        88  88   `8b d8'   88  88           88    `8b    88    ,8P  Y8,      ,8P  Y8,
    88     `8888  Y8a.    .a8P  88    `888'    88  88           88     `8b   88   d8'    `8b    d8'    `8b
    88      `888   `"Y8888Y"'   88     `8'     88  88888888888  88      `8b  88  8P        Y8  8P        Y8

    Copyright © 2023 Kenneth Troldal Balslev

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

#pragma once

// ===== Numerixx Includes
#include "Common.impl.hpp"
#include <Constants.hpp>
#include <Functions.hpp>

// ===== Standard Library Includes
#include <numbers>

/**
 * @file RootSearching.hpp
 * @brief This file contains various search algorithms for finding roots of a function.
 *
 * The search algorithms provided in this file include bracketing searchers and bracket expansion searchers.
 * Bracketing searchers incrementally expand the search bounds either upwards or downwards to find a bracket
 * where a root exists. Bracket expansion searchers expand both the lower and upper bounds symmetrically outwards
 * or expand one bound while keeping the other fixed. These algorithms are useful for finding a bracket where a root
 * exists when the initial guess is not near the actual root.
 *
 * @author Kenneth Troldal Balslev
 * @copyright Copyright (c) 2023 Kenneth Troldal Balslev
 * @license MIT License
 */
namespace nxx::roots {
    // =================================================================================================================
    //
    //  ad88888ba                                                   88           88
    // d8"     "8b                                                  88           ""
    // Y8,                                                          88
    // `Y8aaaaa,     ,adPPYba,  ,adPPYYba,  8b,dPPYba,   ,adPPYba,  88,dPPYba,   88  8b,dPPYba,    ,adPPYb,d8
    //   `"""""8b,  a8P_____88  ""     `Y8  88P'   "Y8  a8"     ""  88P'    "8a  88  88P'   `"8a  a8"    `Y88
    //         `8b  8PP"""""""  ,adPPPPP88  88          8b          88       88  88  88       88  8b       88
    // Y8a     a8P  "8b,   ,aa  88,    ,88  88          "8a,   ,aa  88       88  88  88       88  "8a,   ,d88
    //  "Y88888P"    `"Ybbd8"'  `"8bbdP"Y8  88           `"Ybbd8"'  88       88  88  88       88   `"YbbdP"Y8
    //                                                                                             aa,    ,88
    //                                                                                              "Y8bbdP"
    // =================================================================================================================

    namespace detail {

        template<typename>
        struct SearchingTraits; // Forward declaration with variadic template parameters

        // Generic specialization of SearchingTraits
        template<template<typename, typename> class SearchMethod, // Template template parameter for the search method
            typename FN,
            typename ARG_T>
        struct SearchingTraits<SearchMethod<FN, ARG_T>>
        {
            using FUNCTION_T = FN;
            using RETURN_T = std::invoke_result_t<FN, double>;
        };

        /**
         * @brief Provides a base class template for search-based algorithms.
         *
         * The SearchBase class template serves as a foundational component for
         * algorithms that conduct searches within a specified bound, utilizing a factor for
         * controlling the search behavior. It encapsulates common functionalities such as storing the
         * objective function, maintaining the current search bounds, and a factor influencing the search process.
         * This class enforces certain type constraints on the template parameters to ensure compatibility with
         * search-based algorithms.
         *
         * @tparam SUBCLASS The subclass inheriting from SearchBase.
         * @tparam FUNCTION_T The type of the function involved in the search process.
         * @tparam ARG_T The type of the argument to the function.
         *
         * @note This class template uses SFINAE to enforce type constraints on its template parameters.
         */
        template<typename SUBCLASS, IsFloatInvocable FUNCTION_T, IsFloat ARG_T>
        requires std::same_as<typename SearchingTraits<SUBCLASS>::FUNCTION_T, FUNCTION_T>
                 && nxx::IsFloatInvocable<FUNCTION_T> && nxx::IsFloat<ARG_T>
                 && nxx::IsFloat<typename SearchingTraits<SUBCLASS>::RETURN_T>
        class SearchBase
        {
            friend SUBCLASS;

          public:
            static constexpr bool IsBracketingSearcher =
                true; /**< Flag indicating the class is a bracketing searcher. */

            using RESULT_T = std::invoke_result_t<FUNCTION_T, ARG_T>; /**< Result type of the function. */
            using BOUNDS_T = std::pair<ARG_T, ARG_T>; /**< Type for representing the search bounds. */
            using RATIO_T = ARG_T; /**< Type for representing the search adjustment ratio. */

          protected:
            ~SearchBase() = default; /**< Protected destructor to prevent direct instantiation. */

          private:
            FUNCTION_T m_objective{}; /**< The objective function for the search. */
            BOUNDS_T m_bounds{}; /**< Holds the current search bounds. */
            RATIO_T m_ratio{}; /**< The factor influencing the search process. */

          public:
            /**
             * @brief Constructs the SearchBase with an objective function, bounds from a float struct, and an optional
             * factor.
             * @param objective The function involved in the search process.
             * @param bounds Struct with exactly two members representing the search bounds.
             * @param factor The factor influencing the search process, defaults to the golden ratio.
             * @details This constructor initializes the search algorithm with a specific objective
             *          function, search bounds defined by a struct, and an optional factor.
             */
            SearchBase(FUNCTION_T objective, IsFloatStruct auto bounds, RATIO_T factor = std::numbers::phi)
              : m_objective{ std::move(objective) }, m_bounds(toPair(bounds)), m_ratio{ factor }
            {
                validateBounds(m_bounds);
            }

            /**
             * @brief Constructor for the SearchBase class.
             * @tparam N The size of the bounds array. Must be 2.
             * @param objective The objective function for the search process.
             * @param bounds An array of two elements representing the search bounds.
             * @param factor The factor influencing the search process. Defaults to the golden ratio.
             */
            template<size_t N>
            requires(N == 2)
            SearchBase(FUNCTION_T objective, const ARG_T (&bounds)[N], RATIO_T factor = std::numbers::phi)
              : m_objective{ std::move(objective) }, m_bounds(std::pair{ bounds[0], bounds[1] }), m_ratio{ factor }
            {
                validateBounds(m_bounds);
            }

            /**
             * @brief Sets the search bounds.
             * @param bounds The new bounds to be set for the search, represented as a pair of values.
             * @throws NumerixxError If the search algorithm has not been initialized or if bounds are invalid.
             */
            void setBounds(const BOUNDS_T &bounds)
            {
                m_bounds = toPair(bounds);
                validateBounds(m_bounds);
            }

            /**
             * @brief Sets the search adjustment ratio.
             * @param factor The new ratio to be set for the search process.
             * @throws NumerixxError If the search algorithm has not been initialized or if ratio is invalid.
             */
            void setRatio(RATIO_T factor)
            {
                if (factor < 1.0) throw NumerixxError("Invalid factor.");
                m_ratio = factor;
            }

            SearchBase(const SearchBase &other) = default; /**< Default copy constructor. */
            SearchBase(SearchBase &&other) noexcept = default; /**< Default move constructor. */
            SearchBase &operator=(const SearchBase &other) = default; /**< Default copy assignment operator. */
            SearchBase &operator=(SearchBase &&other) noexcept = default; /**< Default move assignment operator. */

            /**
             * @brief Evaluates the objective function at a given value.
             * @param value The value at which the function is to be evaluated.
             * @return The result of evaluating the function at the specified value.
             */
            [[nodiscard]] RESULT_T evaluate(ARG_T value) const { return m_objective(value); }

            /**
             * @brief Returns the current search bounds.
             * @details This method returns the current bounds being used by the search algorithm.
             *          It throws an exception if the search has not been initialized.
             * @throws NumerixxError If the search algorithm has not been initialized.
             * @return The current search bounds.
             */
            [[nodiscard]] const BOUNDS_T &current() const { return m_bounds; }

            /**
             * @brief Returns the current search adjustment ratio.
             * @details This method returns the ratio influencing the search process.
             *          It throws an exception if the search has not been initialized.
             * @throws NumerixxError If the search algorithm has not been initialized.
             * @return The current search adjustment ratio.
             */
            [[nodiscard]] RATIO_T ratio() const { return m_ratio; }

            void iterate() { std::invoke(static_cast<SUBCLASS &>(*this)); }
        };
    } // namespace detail

    // =================================================================================================================
    //
    //  ad88888ba                                                   88           88        88
    // d8"     "8b                                                  88           88        88
    // Y8,                                                          88           88        88
    // `Y8aaaaa,     ,adPPYba,  ,adPPYYba,  8b,dPPYba,   ,adPPYba,  88,dPPYba,   88        88  8b,dPPYba,
    //   `"""""8b,  a8P_____88  ""     `Y8  88P'   "Y8  a8"     ""  88P'    "8a  88        88  88P'    "8a
    //         `8b  8PP"""""""  ,adPPPPP88  88          8b          88       88  88        88  88       d8
    // Y8a     a8P  "8b,   ,aa  88,    ,88  88          "8a,   ,aa  88       88  Y8a.    .a8P  88b,   ,a8"
    //  "Y88888P"    `"Ybbd8"'  `"8bbdP"Y8  88           `"Ybbd8"'  88       88   `"Y8888Y"'   88`YbbdP"'
    //                                                                                         88
    //                                                                                         88
    // =================================================================================================================

    /**
     * @brief Defines the BracketSearchUp class for performing upward bracketing search.
     *
     * The BracketSearchUp class template is a specialized search algorithm designed to incrementally
     * expand the search bounds upwards (increasing values) to find a bracket where a root exists.
     * It inherits from a base class that provides common functionalities for search-based algorithms,
     * and adds the specific logic for upward bracketing. This class is templated to accept a function
     * and an optional argument type, along with an optional factor that controls the expansion of the bounds.
     *
     * @tparam FN The type of the function for which the bracket is being searched.
     * @tparam ARG_T The type of the argument to the function, defaults to double.
     */

    template<IsFloatInvocable FN, IsFloat ARG_T = double>
    class BracketSearchUp final : public detail::SearchBase<BracketSearchUp<FN, ARG_T>, FN, ARG_T>
    {
        using BASE =
            detail::SearchBase<BracketSearchUp<FN, ARG_T>, FN, ARG_T>; /**< Base class alias for readability. */

      public:
        using BASE::BASE; /**< Inherits constructors from SearchBase. */

        /**
         * @brief Performs a single iteration of the upward bracketing search.
         * @details This method expands the search bounds upwards if the current bounds do not bracket a root.
         *          The expansion factor controls the rate at which the bounds are expanded.
         */
        void operator()()
        {
            const auto &bounds = BASE::current();

            // Check if current bounds already bracket a root
            if (BASE::evaluate(bounds.first) * BASE::evaluate(bounds.second) < 0.0) return;

            // Expand the bounds upwards
            auto newBounds = bounds;
            newBounds.first = bounds.second;
            newBounds.second = bounds.second + (bounds.second - bounds.first) * BASE::ratio();
            BASE::setBounds(newBounds);
        }
    };

    /**
     * @brief Deduction guides for BracketSearchUp class.
     * Allows the type of BracketSearchUp class to be deduced from the constructor parameters.
     */
    template<typename FN, typename ARG_T, typename FACTOR_T = ARG_T>
    requires IsFloatInvocable<FN> && IsFloat<ARG_T> && IsFloat<FACTOR_T>
    BracketSearchUp(FN, std::initializer_list<ARG_T>, FACTOR_T factor = std::numbers::phi)
        -> BracketSearchUp<FN, ARG_T>;

    template<typename FN, typename BOUNDS_T, typename FACTOR_T = StructCommonType_t<BOUNDS_T>>
    requires IsFloatInvocable<FN> && IsFloatStruct<BOUNDS_T>
    BracketSearchUp(FN, BOUNDS_T, FACTOR_T factor = std::numbers::phi)
        -> BracketSearchUp<FN, StructCommonType_t<BOUNDS_T>>;

    // =================================================================================================================
    //
    //  ad88888ba                                                   88           88888888ba,
    // d8"     "8b                                                  88           88      `"8b
    // Y8,                                                          88           88        `8b
    // `Y8aaaaa,     ,adPPYba,  ,adPPYYba,  8b,dPPYba,   ,adPPYba,  88,dPPYba,   88         88   ,adPPYba,   8b      db
    // d8  8b,dPPYba,
    //   `"""""8b,  a8P_____88  ""     `Y8  88P'   "Y8  a8"     ""  88P'    "8a  88         88  a8"     "8a  `8b    d88b
    //   d8'  88P'   `"8a
    //         `8b  8PP"""""""  ,adPPPPP88  88          8b          88       88  88         8P  8b       d8   `8b d8'`8b
    //         d8'   88       88
    // Y8a     a8P  "8b,   ,aa  88,    ,88  88          "8a,   ,aa  88       88  88      .a8P   "8a,   ,a8"    `8bd8'
    // `8bd8'    88       88
    //  "Y88888P"    `"Ybbd8"'  `"8bbdP"Y8  88           `"Ybbd8"'  88       88  88888888Y"'     `"YbbdP"'       YP YP
    //  88       88
    //
    //
    // =================================================================================================================

    /**
     * @brief Defines the BracketSearchDown class for performing downward bracketing search.
     *
     * The BracketSearchDown class template is a specialized search algorithm designed to incrementally
     * expand the search bounds downwards (decreasing values) to find a bracket where a root exists.
     * It inherits from a base class that provides common functionalities for search-based algorithms,
     * and adds the specific logic for downward bracketing. This class is templated to accept a function
     * and an optional argument type, along with an optional factor that controls the contraction of the bounds.
     *
     * @tparam FN The type of the function for which the bracket is being searched.
     * @tparam ARG_T The type of the argument to the function, defaults to double.
     */
    template<IsFloatInvocable FN, IsFloat ARG_T = double>
    class BracketSearchDown final : public detail::SearchBase<BracketSearchDown<FN, ARG_T>, FN, ARG_T>
    {
        using BASE =
            detail::SearchBase<BracketSearchDown<FN, ARG_T>, FN, ARG_T>; /**< Base class alias for readability. */

      public:
        using BASE::BASE; /**< Inherits constructors from SearchBase. */

        /**
         * @brief Performs a single iteration of the downward bracketing search.
         * @details This method expands the search bounds downwards if the current bounds do not bracket a root.
         *          The expansion factor controls the rate at which the bounds are contracted.
         */
        void operator()()
        {
            const auto &bounds = BASE::current();

            // Check if current bounds already bracket a root
            if (BASE::evaluate(bounds.first) * BASE::evaluate(bounds.second) < 0.0) return;

            // Expand the bounds downwards
            auto newBounds = bounds;
            newBounds.second = bounds.first;
            newBounds.first = bounds.first - (bounds.second - bounds.first) * BASE::ratio();
            BASE::setBounds(newBounds);
        }
    };

    /**
     * @brief Deduction guides for BracketSearchDown class.
     * Allows the type of BracketSearchDown class to be deduced from the constructor parameters.
     */
    template<typename FN, typename ARG_T, typename FACTOR_T = ARG_T>
    requires IsFloatInvocable<FN> && IsFloat<ARG_T> && IsFloat<FACTOR_T>
    BracketSearchDown(FN, std::initializer_list<ARG_T>, FACTOR_T factor = std::numbers::phi)
        -> BracketSearchDown<FN, ARG_T>;

    template<typename FN, typename BOUNDS_T, typename FACTOR_T = StructCommonType_t<BOUNDS_T>>
    requires IsFloatInvocable<FN> && IsFloatStruct<BOUNDS_T>
    BracketSearchDown(FN, BOUNDS_T, FACTOR_T factor = std::numbers::phi)
        -> BracketSearchDown<FN, StructCommonType_t<BOUNDS_T>>;

    // =================================================================================================================
    //
    // 88888888888                                                              88  88        88
    // 88                                                                       88  88        88
    // 88                                                                       88  88        88
    // 88aaaaa      8b,     ,d8  8b,dPPYba,   ,adPPYYba,  8b,dPPYba,    ,adPPYb,88  88        88  8b,dPPYba,
    // 88"""""       `Y8, ,8P'   88P'    "8a  ""     `Y8  88P'   `"8a  a8"    `Y88  88        88  88P'    "8a
    // 88              )888(     88       d8  ,adPPPPP88  88       88  8b       88  88        88  88       d8
    // 88            ,d8" "8b,   88b,   ,a8"  88,    ,88  88       88  "8a,   ,d88  Y8a.    .a8P  88b,   ,a8"
    // 88888888888  8P'     `Y8  88`YbbdP"'   `"8bbdP"Y8  88       88   `"8bbdP"Y8   `"Y8888Y"'   88`YbbdP"'
    //                           88                                                               88
    //                           88                                                               88
    // =================================================================================================================

    /**
     * @brief Defines the BracketExpandUp class for performing upward bracket expansion.
     *
     * The BracketExpandUp class template is a specialized search algorithm designed to incrementally
     * expand the upper bound upwards (increasing values) while keeping the lower bound fixed. This is useful
     * for finding a bracket where a root exists when the initial guess is lower than the actual root. It inherits
     * from a base class that provides common functionalities for search-based algorithms and adds the specific
     * logic for upward bracket expansion. This class is templated to accept a function and an optional argument type,
     * along with an optional factor that controls the expansion of the upper bound.
     *
     * @tparam FN The type of the function for which the bracket is being expanded.
     * @tparam ARG_T The type of the argument to the function, defaults to double.
     */
    template<IsFloatInvocable FN, IsFloat ARG_T = double>
    class BracketExpandUp final : public detail::SearchBase<BracketExpandUp<FN, ARG_T>, FN, ARG_T>
    {
        using BASE =
            detail::SearchBase<BracketExpandUp<FN, ARG_T>, FN, ARG_T>; /**< Base class alias for readability. */

      public:
        using BASE::BASE; /**< Inherits constructors from SearchBase. */

        /**
         * @brief Performs a single iteration of the upward bracket expansion.
         * @details This method expands the upper bound upwards if the current bounds do not bracket a root.
         *          The expansion factor controls the rate at which the upper bound is expanded.
         */
        void operator()()
        {
            const auto &bounds = BASE::current();

            // Check if current bounds already bracket a root
            if (BASE::evaluate(bounds.first) * BASE::evaluate(bounds.second) < 0.0) return;

            // Expand the upper bound upwards
            auto newBounds = bounds;
            newBounds.second = bounds.second + (bounds.second - bounds.first) * BASE::ratio();
            BASE::setBounds(newBounds);
        }
    };

    /**
     * @brief Deduction guides for BracketExpandUp class.
     * Allows the type of BracketExpandUp class to be deduced from the constructor parameters.
     */
    template<typename FN, typename ARG_T, typename FACTOR_T = ARG_T>
    requires IsFloatInvocable<FN> && IsFloat<ARG_T> && IsFloat<FACTOR_T>
    BracketExpandUp(FN, std::initializer_list<ARG_T>, FACTOR_T factor = std::numbers::phi)
        -> BracketExpandUp<FN, ARG_T>;

    template<typename FN, typename BOUNDS_T, typename FACTOR_T = StructCommonType_t<BOUNDS_T>>
    requires IsFloatInvocable<FN> && IsFloatStruct<BOUNDS_T>
    BracketExpandUp(FN, BOUNDS_T, FACTOR_T factor = std::numbers::phi)
        -> BracketExpandUp<FN, StructCommonType_t<BOUNDS_T>>;

    // =================================================================================================================
    //
    // 88888888888                                                              88  88888888ba,
    // 88                                                                       88  88      `"8b
    // 88                                                                       88  88        `8b
    // 88aaaaa      8b,     ,d8  8b,dPPYba,   ,adPPYYba,  8b,dPPYba,    ,adPPYb,88  88         88   ,adPPYba,   8b db d8
    // 8b,dPPYba, 88"""""       `Y8, ,8P'   88P'    "8a  ""     `Y8  88P'   `"8a  a8"    `Y88  88         88  a8" "8a
    // `8b    d88b d8' 88P'   `"8a 88              )888(     88       d8  ,adPPPPP88  88       88  8b       88  88 8P 8b
    // d8   `8b  d8'`8b d8' 88       88 88            ,d8" "8b,   88b,   ,a8"  88,    ,88  88       88  "8a,   ,d88  88
    // .a8P   "8a,   ,a8"    `8bd8'  `8bd8' 88       88 88888888888  8P'     `Y8  88`YbbdP"'   `"8bbdP"Y8  88       88
    // `"8bbdP"Y8  88888888Y"'     `"YbbdP"'       YP      YP 88       88
    //                           88
    //                           88
    // =================================================================================================================

    /**
     * @brief Defines the BracketExpandDown class for performing downward bracket expansion.
     *
     * The BracketExpandDown class template is a specialized search algorithm designed to incrementally
     * expand the lower bound downwards (decreasing values) while keeping the upper bound fixed. This is useful
     * for finding a bracket where a root exists when the initial guess is higher than the actual root. It inherits
     * from a base class that provides common functionalities for search-based algorithms and adds the specific
     * logic for downward bracket expansion. This class is templated to accept a function and an optional argument type,
     * along with an optional factor that controls the expansion of the lower bound.
     *
     * @tparam FN The type of the function for which the bracket is being expanded.
     * @tparam ARG_T The type of the argument to the function, defaults to double.
     */
    template<IsFloatInvocable FN, IsFloat ARG_T = double>
    class BracketExpandDown final : public detail::SearchBase<BracketExpandDown<FN, ARG_T>, FN, ARG_T>
    {
        using BASE =
            detail::SearchBase<BracketExpandDown<FN, ARG_T>, FN, ARG_T>; /**< Base class alias for readability. */

      public:
        using BASE::BASE; /**< Inherits constructors from SearchBase. */

        /**
         * @brief Performs a single iteration of the downward bracket expansion.
         * @details This method expands the lower bound downwards if the current bounds do not bracket a root.
         *          The expansion factor controls the rate at which the lower bound is expanded.
         */
        void operator()()
        {
            const auto &bounds = BASE::current();

            // Check if current bounds already bracket a root
            if (BASE::evaluate(bounds.first) * BASE::evaluate(bounds.second) < 0.0) return;

            // Expand the lower bound downwards
            auto newBounds = bounds;
            newBounds.first = bounds.first - (bounds.second - bounds.first) * BASE::ratio();
            BASE::setBounds(newBounds);
        }
    };

    /**
     * @brief Deduction guides for BracketExpandDown class.
     * Allows the type of BracketExpandDown class to be deduced from the constructor parameters.
     */
    template<typename FN, typename ARG_T, typename FACTOR_T = ARG_T>
    requires IsFloatInvocable<FN> && IsFloat<ARG_T> && IsFloat<FACTOR_T>
    BracketExpandDown(FN, std::initializer_list<ARG_T>, FACTOR_T factor = std::numbers::phi)
        -> BracketExpandDown<FN, ARG_T>;

    template<typename FN, typename BOUNDS_T, typename FACTOR_T = StructCommonType_t<BOUNDS_T>>
    requires IsFloatInvocable<FN> && IsFloatStruct<BOUNDS_T>
    BracketExpandDown(FN, BOUNDS_T, FACTOR_T factor = std::numbers::phi)
        -> BracketExpandDown<FN, StructCommonType_t<BOUNDS_T>>;

    // =================================================================================================================
    //
    // 88888888888                                                              88    ,ad8888ba,
    // 88                                                                       88   d8"'    `"8b                  ,d
    // 88                                                                       88  d8'        `8b                 88
    // 88aaaaa      8b,     ,d8  8b,dPPYba,   ,adPPYYba,  8b,dPPYba,    ,adPPYb,88  88          88  88       88  MM88MMM
    // 88"""""       `Y8, ,8P'   88P'    "8a  ""     `Y8  88P'   `"8a  a8"    `Y88  88          88  88       88    88
    // 88              )888(     88       d8  ,adPPPPP88  88       88  8b       88  Y8,        ,8P  88       88    88
    // 88            ,d8" "8b,   88b,   ,a8"  88,    ,88  88       88  "8a,   ,d88   Y8a.    .a8P   "8a,   ,a88    88,
    // 88888888888  8P'     `Y8  88`YbbdP"'   `"8bbdP"Y8  88       88   `"8bbdP"Y8    `"Y8888Y"'     `"YbbdP'Y8    "Y888
    //                           88
    //                           88
    // =================================================================================================================

    /**
     * @brief Defines the BracketExpandOut class for performing outward bracket expansion.
     *
     * The BracketExpandOut class template is a specialized search algorithm designed to incrementally
     * expand both the lower and upper bounds symmetrically outwards. This is useful for finding a bracket
     * where a root exists when the initial guess is not near the actual root. It inherits from a base class
     * that provides common functionalities for search-based algorithms and adds the specific logic for
     * outward bracket expansion. This class is templated to accept a function and an optional argument type,
     * along with an optional factor that controls the symmetric expansion of the bounds.
     *
     * @tparam FN The type of the function for which the bracket is being expanded.
     * @tparam ARG_T The type of the argument to the function, defaults to double.
     */
    template<IsFloatInvocable FN, IsFloat ARG_T = double>
    class BracketExpandOut final : public detail::SearchBase<BracketExpandOut<FN, ARG_T>, FN, ARG_T>
    {
        using BASE =
            detail::SearchBase<BracketExpandOut<FN, ARG_T>, FN, ARG_T>; /**< Base class alias for readability. */

      public:
        using BASE::BASE; /**< Inherits constructors from SearchBase. */

        /**
         * @brief Performs a single iteration of the outward bracket expansion.
         * @details This method expands both the lower and upper bounds outwards symmetrically if the current
         *          bounds do not bracket a root. The expansion factor controls the rate at which the bounds are
         * expanded.
         */
        void operator()()
        {
            const auto &bounds = BASE::current();

            // Check if current bounds already bracket a root
            if (BASE::evaluate(bounds.first) * BASE::evaluate(bounds.second) < 0.0) return;

            // Expand the bounds symmetrically outwards
            auto newBounds = bounds;
            newBounds.first = bounds.first - (bounds.second - bounds.first) * BASE::ratio() / 2.0;
            newBounds.second = bounds.second + (bounds.second - bounds.first) * BASE::ratio() / 2.0;
            BASE::setBounds(newBounds);
        }
    };

    /**
     * @brief Deduction guides for BracketExpandOut class.
     * Allows the type of BracketExpandOut class to be deduced from the constructor parameters.
     */
    template<typename FN, typename ARG_T, typename FACTOR_T = ARG_T>
    requires IsFloatInvocable<FN> && IsFloat<ARG_T> && IsFloat<FACTOR_T>
    BracketExpandOut(FN, std::initializer_list<ARG_T>, FACTOR_T factor = std::numbers::phi)
        -> BracketExpandOut<FN, ARG_T>;

    template<typename FN, typename BOUNDS_T, typename FACTOR_T = StructCommonType_t<BOUNDS_T>>
    requires IsFloatInvocable<FN> && IsFloatStruct<BOUNDS_T>
    BracketExpandOut(FN, BOUNDS_T, FACTOR_T factor = std::numbers::phi)
        -> BracketExpandOut<FN, StructCommonType_t<BOUNDS_T>>;

    // =================================================================================================================
    //
    //  ad88888ba                88                    88  88               88           88
    // d8"     "8b               88                    88  ""               ""           88
    // Y8,                       88                    88                                88
    // `Y8aaaaa,    88       88  88,dPPYba,    ,adPPYb,88  88  8b       d8  88   ,adPPYb,88   ,adPPYba,
    //   `"""""8b,  88       88  88P'    "8a  a8"    `Y88  88  `8b     d8'  88  a8"    `Y88  a8P_____88
    //         `8b  88       88  88       d8  8b       88  88   `8b   d8'   88  8b       88  8PP"""""""
    // Y8a     a8P  "8a,   ,a88  88b,   ,a8"  "8a,   ,d88  88    `8b,d8'    88  "8a,   ,d88  "8b,   ,aa
    //  "Y88888P"    `"YbbdP'Y8  8Y"Ybbd8"'    `"8bbdP"Y8  88      "8"      88   `"8bbdP"Y8   `"Ybbd8"'
    //
    //
    // =================================================================================================================

    /**
     * @brief Defines the BracketSubdivide class for performing bracket subdivision search.
     *
     * The BracketSubdivide class template is a specialized search algorithm designed to subdivide
     * the current search bounds into smaller segments in an attempt to find a bracket where a root exists.
     * It inherits from a base class that provides common functionalities for search-based algorithms,
     * and adds the specific logic for subdividing the search bounds. This class is templated to accept
     * a function and an optional argument type, along with an optional factor that controls the subdivision process.
     *
     * @tparam FN The type of the function for which the bracket is being searched.
     * @tparam ARG_T The type of the argument to the function, defaults to double.
     */
    template<IsFloatInvocable FN, IsFloat ARG_T = double>
    class BracketSubdivide final : public detail::SearchBase<BracketSubdivide<FN, ARG_T>, FN, ARG_T>
    {
        using BASE =
            detail::SearchBase<BracketSubdivide<FN, ARG_T>, FN, ARG_T>; /**< Base class alias for readability. */

      public:
        using BASE::BASE; /**< Inherits constructors from SearchBase. */

        /**
         * @brief Performs a single iteration of the bracket subdivision search.
         * @details This method subdivides the current bounds into smaller segments based on the factor,
         *          attempting to find a segment where the function changes sign, indicating a bracket.
         *          If no such segment is found, the factor is doubled to increase the search range.
         */
        void operator()()
        {
            const auto &bounds = BASE::current();
            if (BASE::evaluate(bounds.first) * BASE::evaluate(bounds.second) < 0.0) return;

            size_t factor = std::ceil(BASE::ratio());
            auto diff = (bounds.second - bounds.first) / factor;
            auto lower = bounds.first;
            // auto   upper  = std::min(bounds.first + diff, bounds.second);
            auto upper = std::min(bounds.first + diff, bounds.second);
            for (size_t i = 0; i < factor; ++i) {
                if (BASE::evaluate(lower) * BASE::evaluate(upper) < 0.0) {
                    BASE::setBounds({ lower, upper });
                    return;
                }
                lower = upper;
                upper = std::min(upper + diff, bounds.second);
            }

            BASE::setRatio(BASE::ratio() * 2.0); // Increase the factor to expand the search range
        }
    };

    /**
     * @brief Deduction guides for BracketSubdivide class.
     * Allows the type of BracketSubdivide class to be deduced from the constructor parameters.
     */
    template<typename FN, typename ARG_T, typename FACTOR_T = ARG_T>
    requires IsFloatInvocable<FN> && IsFloat<ARG_T> && IsFloat<FACTOR_T>
    BracketSubdivide(FN, std::initializer_list<ARG_T>, FACTOR_T factor = std::numbers::phi)
        -> BracketSubdivide<FN, ARG_T>;

    template<typename FN, typename BOUNDS_T, typename FACTOR_T = StructCommonType_t<BOUNDS_T>>
    requires IsFloatInvocable<FN> && IsFloatStruct<BOUNDS_T>
    BracketSubdivide(FN, BOUNDS_T, FACTOR_T factor = std::numbers::phi)
        -> BracketSubdivide<FN, StructCommonType_t<BOUNDS_T>>;

    // =================================================================================================================
    //                                                            88
    //                                                            88
    //                                                            88
    // ,adPPYba,   ,adPPYba,  ,adPPYYba,  8b,dPPYba,   ,adPPYba,  88,dPPYba,
    // I8[    ""  a8P_____88  ""     `Y8  88P'   "Y8  a8"     ""  88P'    "8a
    //  `"Y8ba,   8PP"""""""  ,adPPPPP88  88          8b          88       88
    // aa    ]8I  "8b,   ,aa  88,    ,88  88          "8a,   ,aa  88       88
    // `"YbbdP"'   `"Ybbd8"'  `"8bbdP"Y8  88           `"Ybbd8"'  88       88
    //
    // =================================================================================================================

    template<std::integral ITER_T, IsFloat RESULT_T>
    using SearchIterData = std::tuple<ITER_T, RESULT_T, RESULT_T>;

    template<IsFloatInvocable FN_T, IsFloat RATIO_T, std::integral ITER_T>
    class SearchStopToken
    {
        FN_T m_fn; /**< The function for which the root is being searched. */
        RATIO_T m_ratio; /**< The epsilon value for the termination condition. */
        ITER_T m_maxiter; /**< The maximum iteration count for the termination condition. */

      public:
        /**
         * @brief Default constructor. Initializes the epsilon value and maximum iteration count to default values.
         */
        explicit SearchStopToken() : m_ratio(std::numbers::phi), m_maxiter(iterations<double>()) {}

        /**
         * @brief Constructor. Initializes the epsilon value and maximum iteration count to the specified values.
         * @param eps The epsilon value for the termination condition.
         * @param maxiter The maximum iteration count for the termination condition.
         */
        explicit SearchStopToken(RATIO_T ratio, ITER_T maxiter) : m_ratio(ratio), m_maxiter(maxiter) {}

        /**
         * @brief Constructor. Initializes the epsilon value and maximum iteration count to the specified values.
         * @param maxiter The maximum iteration count for the termination condition.
         * @param eps The epsilon value for the termination condition.
         */
        explicit SearchStopToken(ITER_T maxiter, RATIO_T ratio) : m_ratio(ratio), m_maxiter(maxiter) {}

        /**
         * @brief Constructor. Initializes the epsilon value to the specified value and the maximum iteration count to a
         * default value.
         * @param eps The epsilon value for the termination condition.
         */
        explicit SearchStopToken(RATIO_T ratio) : m_ratio(ratio), m_maxiter(iterations<double>()) {}

        /**
         * @brief Constructor. Initializes the maximum iteration count to the specified value and the epsilon value to a
         * default value.
         * @param maxiter The maximum iteration count for the termination condition.
         */
        explicit SearchStopToken(ITER_T maxiter) : m_ratio(std::numbers::phi), m_maxiter(maxiter) {}

        /**
         * @brief Checks if the termination condition is met.
         * @param data The current iteration data.
         * @return true if the termination condition is met, false otherwise.
         */
        bool operator()(const auto &data) const
        {
            const auto &[iter, lower, upper] = data;

            if (m_fn(lower) * m_fn(upper) <= 0.0) return true;
            if (iter >= m_maxiter) return true;

            return false;
        }
    };

    /**
     * @brief Deduction guides for the BracketStopToken class template.
     * Allows the type of BracketStopToken to be deduced from the constructor parameters.
     */
    //    SearchStopToken() -> SearchStopToken<double, size_t>;
    //    SearchStopToken(IsFloat auto eps) -> SearchStopToken<decltype(eps), size_t>;
    //    SearchStopToken(std::integral auto maxiter) -> SearchStopToken<double, decltype(maxiter)>;
    //    SearchStopToken(IsFloat auto eps, std::integral auto maxiter)
    //        -> SearchStopToken<decltype(eps), decltype(maxiter)>;
    //    SearchStopToken(std::integral auto maxiter, IsFloat auto eps)
    //        -> SearchStopToken<decltype(eps), decltype(maxiter)>;


    namespace detail {

        template<std::integral ITER_T, IsFloat RESULT_T>
        class SearchResult
        {
            SearchIterData<ITER_T, RESULT_T>
                m_iterData; /**< The IterData object holding the result of the root-finding problem. */

          public:
            /**
             * @brief Constructs the BracketSolverResult with an IterData object.
             * @param iterData The IterData object holding the result of the root-finding problem.
             */
            explicit SearchResult(SearchIterData<ITER_T, RESULT_T> iterData) : m_iterData(iterData) {}

            SearchResult(const SearchResult &) = delete; // No copy constructor
            SearchResult(SearchResult &&) = delete; // No move constructor

            SearchResult &operator=(const SearchResult &) = delete; // No copy assignment
            SearchResult &operator=(SearchResult &&) = delete; // No move assignment

            /**
             * @brief Returns the result of the root-finding problem.
             * @tparam OUTPUT_T The type of the output. Defaults to the type of the result of the root-finding problem.
             * @return The result of the root-finding problem. If OUTPUT_T is a class, it is constructed with the
             * IterData object. Otherwise, the guess from the IterData object is returned.
             * @note This method is only available for rvalue references.
             */
            template<typename OUTPUT_T = RESULT_T>
            auto result() &&
            {
                if constexpr (std::is_class_v<OUTPUT_T>)
                    return OUTPUT_T{}(m_iterData);
                else
                    //                    return m_iterData.guess;
                    return std::make_pair(std::get<1>(m_iterData), std::get<2>(m_iterData));
            }

            /**
             * @brief Returns the result of the root-finding problem using a specified outputter.
             * @tparam OUTPUTTER_T The type of the outputter. Must be a callable object that accepts an IterData object.
             * @param outputter The outputter to use to format the result.
             * @return The result of the root-finding problem, formatted by the outputter.
             * @note This method is only available for rvalue references.
             */
            template<typename OUTPUTTER_T>
            auto result(OUTPUTTER_T outputter) &&
            {
                return outputter(m_iterData);
            }
        };

        template<typename ITER_T, typename RESULT_T>
        SearchResult(BracketIterData<ITER_T, RESULT_T>) -> SearchResult<ITER_T, RESULT_T>;

        template<typename SOLVER, typename TOKEN_T>
        requires SOLVER::IsBracketingSearcher
        auto search_impl(const SOLVER &solver, const TOKEN_T &terminator)
        {
            SOLVER _solver = solver;
            TOKEN_T _terminator = terminator;
            using ARG_T = typename SOLVER::RESULT_T;

            const auto &[lower, upper] = _solver.current();
            size_t iter = 0;

            SearchIterData<size_t, ARG_T> iterData{ iter, lower, upper };
            auto &[_iter, _lower, _upper] = iterData;

            while (true) {
                _iter = iter;
                _lower = lower;
                _upper = upper;

                if (_terminator(iterData)) break;
                _solver.iterate();
                ++iter;
            }

            return SearchResult(iterData);
        }


        template<template<typename, typename> class SOLVER_T, typename FN_T, typename STRUCT_T, typename... Args>
        requires(sizeof...(Args) <= 2)
        auto search_common(FN_T func, const STRUCT_T &bounds, const Args &...args)
        {
            using TUPLE_T = std::tuple<Args...>;
            using SOLVER = SOLVER_T<FN_T, StructCommonType_t<STRUCT_T>>;
            using ITERDATA_T = BracketIterData<size_t, StructCommonType_t<STRUCT_T>>;

            TUPLE_T args_tuple = std::make_tuple(args...);

            // Zero arguments are passed...
            if constexpr (std::tuple_size_v<TUPLE_T> == 0)
                return detail::search_impl(SOLVER(func, bounds), SearchStopToken<FN_T, double, size_t>{});

            // One argument is passed...
            else if constexpr (std::tuple_size_v<TUPLE_T> == 1) {
                using ArgType = std::tuple_element_t<0, TUPLE_T>;

                // If the argument is a floating point or integral type, use it as maxiter/eps
                if constexpr (std::is_floating_point_v<ArgType> || std::is_integral_v<ArgType>)
                    return detail::search_impl(SOLVER(func, bounds), SearchStopToken(std::get<0>(args_tuple)));

                // If the argument is a callable, use as a stop token
                // TODO: Should be able to accept functors that take any king of argument, not just references.
                else if constexpr (std::same_as<std::invoke_result_t<ArgType, ITERDATA_T &>, bool>)
                    return detail::search_impl(SOLVER(func, bounds), std::get<0>(args_tuple));

                else
                    std::invoke(
                        []<bool flag = false>() { static_assert(flag, "Invalid argument passed to search_common"); });
            }

            // Two arguments are passed...
            else if constexpr (std::tuple_size_v<TUPLE_T> == 2) {
                // Unpack and use the two arguments
                using ArgType1 = std::tuple_element_t<0, decltype(args_tuple)>;
                using ArgType2 = std::tuple_element_t<1, decltype(args_tuple)>;
                static_assert(std::is_floating_point_v<ArgType1> != std::is_floating_point_v<ArgType2>,
                    "Two arguments must be one floating point and "
                    "one integral type");
                return detail::search_impl(
                    SOLVER(func, bounds), SearchStopToken(std::get<0>(args_tuple), std::get<1>(args_tuple)));
            } else
                std::invoke(
                    []<bool flag = false>() { static_assert(flag, "Invalid argument passed to fsolve_common"); });
        }
    } // namespace detail


    template<template<typename, typename> class SOLVER_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloatStruct STRUCT_T,
        typename... ARGS>
    auto search(FN_T func, STRUCT_T bounds, ARGS... args)
    {
        return detail::search_common<SOLVER_T>(func, bounds, args...);
    }

    template<template<typename, typename> class SOLVER_T,
        typename TOKEN_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloatStruct STRUCT_T>
    // TODO: Should be able to accept functors that take any kind of argument, not just references.
    requires std::invocable<TOKEN_T, SearchIterData<size_t, StructCommonType_t<STRUCT_T>> &>
    auto search(FN_T func, STRUCT_T bounds)
    {
        return detail::search_common<SOLVER_T>(func, bounds, TOKEN_T{});
    }

    template<template<typename, typename> class SOLVER_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloat ARG_T,
        size_t N,
        typename... ARGS>
    requires(N == 2)
    auto search(FN_T func, const ARG_T (&bounds)[N], ARGS... args)
    {
        return detail::search_common<SOLVER_T>(func, std::pair{ bounds[0], bounds[1] }, args...);
    }

    template<template<typename, typename> class SOLVER_T,
        typename TOKEN_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloat ARG_T,
        size_t N>
    // TODO: Should be able to accept functors that take any kingd of argument, not just references.
    requires(N == 2) && std::invocable<TOKEN_T, SearchIterData<size_t, ARG_T> &>
    auto search(FN_T func, const ARG_T (&bounds)[N])
    {
        auto bounds_ = std::span(bounds); // Mostly needed to suppress compiler warning.
        return detail::search_common<SOLVER_T>(func, std::pair{ bounds_.front(), bounds_.back() }, TOKEN_T{});
    }


} // namespace nxx::roots

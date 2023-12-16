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

#ifndef NUMERIXX_ROOTSEARCHING_HPP
#define NUMERIXX_ROOTSEARCHING_HPP

// ===== Numerixx Includes
#include "RootCommon.hpp"
#include <Constants.hpp>

// ===== External Includes
#include <tl/expected.hpp>

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
namespace nxx::roots
{
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

    namespace detail
    {
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
        template< typename SUBCLASS, IsFloatInvocable FUNCTION_T, IsFloat ARG_T >
        requires std::same_as< typename SearchingTraits< SUBCLASS >::FUNCTION_T, FUNCTION_T > && nxx::IsFloatInvocable< FUNCTION_T > &&
                 nxx::IsFloat< ARG_T > && nxx::IsFloat< typename SearchingTraits< SUBCLASS >::RETURN_T >
        class SearchBase
        {
            friend SUBCLASS;

        public:
            static constexpr bool IsBracketingSearcher = true; /**< Flag indicating the class is a bracketing searcher. */

            using RESULT_T = std::invoke_result_t< FUNCTION_T, ARG_T >; /**< Result type of the function. */
            using BOUNDS_T = std::pair< ARG_T, ARG_T >;                 /**< Type for representing the search bounds. */
            using RATIO_T  = ARG_T;                                     /**< Type for representing the search adjustment ratio. */

        protected:
            ~SearchBase() = default; /**< Protected destructor to prevent direct instantiation. */

        private:
            FUNCTION_T m_objective {}; /**< The objective function for the search. */
            BOUNDS_T   m_bounds {};    /**< Holds the current search bounds. */
            RATIO_T    m_ratio {};     /**< The factor influencing the search process. */

        public:
            /**
             * @brief Constructs the SearchBase with an objective function, bounds from a float struct, and an optional factor.
             * @param objective The function involved in the search process.
             * @param bounds Struct with exactly two members representing the search bounds.
             * @param factor The factor influencing the search process, defaults to the golden ratio.
             * @details This constructor initializes the search algorithm with a specific objective
             *          function, search bounds defined by a struct, and an optional factor.
             */
            SearchBase(FUNCTION_T objective, IsFloatStruct auto bounds, RATIO_T factor = std::numbers::phi)
                : m_objective { std::move(objective) },
                  m_ratio { factor }
            {
                init(bounds, factor);
            }

            template< size_t N >
            requires(N == 2)
            SearchBase(FUNCTION_T objective, const ARG_T (&bounds)[N], RATIO_T factor = std::numbers::phi)
                : m_objective { std::move(objective) },
                  m_ratio { factor }
            {
                init(std::pair { bounds[0], bounds[1] }, factor);
            }

            /**
             * @brief Sets the search bounds.
             * @param bounds The new bounds to be set for the search, represented as a pair of values.
             * @throws NumerixxError If the search algorithm has not been initialized or if bounds are invalid.
             */
            void setBounds(const BOUNDS_T& bounds)
            {
                auto [lower, upper] = bounds;
                if (lower == upper) throw NumerixxError("Invalid bounds.");

                if (lower > upper)
                    m_bounds = BOUNDS_T { upper, lower };
                else
                    m_bounds = BOUNDS_T { lower, upper };
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

            SearchBase(const SearchBase& other)                = default; /**< Default copy constructor. */
            SearchBase(SearchBase&& other) noexcept            = default; /**< Default move constructor. */
            SearchBase& operator=(const SearchBase& other)     = default; /**< Default copy assignment operator. */
            SearchBase& operator=(SearchBase&& other) noexcept = default; /**< Default move assignment operator. */

            /**
             * @brief Initializes the search with bounds from a float struct and an optional factor.
             * @param bounds Struct with exactly two members representing the search bounds.
             * @param factor The factor influencing the search process, defaults to the golden ratio.
             * @throws NumerixxError If the factor is invalid.
             */
            void init(IsFloatStruct auto bounds, RATIO_T factor = std::numbers::phi)
            {
                auto [lower, upper] = bounds;
                setBounds(BOUNDS_T { lower, upper });
                setRatio(factor);
            }

            /**
             * @brief Evaluates the objective function at a given value.
             * @param value The value at which the function is to be evaluated.
             * @return The result of evaluating the function at the specified value.
             */
            [[nodiscard]]
            RESULT_T evaluate(ARG_T value) const
            {
                return m_objective(value);
            }

            /**
             * @brief Returns the current search bounds.
             * @details This method returns the current bounds being used by the search algorithm.
             *          It throws an exception if the search has not been initialized.
             * @throws NumerixxError If the search algorithm has not been initialized.
             * @return The current search bounds.
             */
            [[nodiscard]]
            const BOUNDS_T& current() const
            {
                return m_bounds;
            }

            /**
             * @brief Returns the current search adjustment ratio.
             * @details This method returns the ratio influencing the search process.
             *          It throws an exception if the search has not been initialized.
             * @throws NumerixxError If the search algorithm has not been initialized.
             * @return The current search adjustment ratio.
             */
            [[nodiscard]]
            RATIO_T ratio() const
            {
                return m_ratio;
            }
        };
    }    // namespace detail

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

    template< IsFloatInvocable FN, IsFloat ARG_T = double >
    class BracketSearchUp final : public detail::SearchBase< BracketSearchUp< FN, ARG_T >, FN, ARG_T >
    {
        using BASE = detail::SearchBase< BracketSearchUp< FN, ARG_T >, FN, ARG_T >; /**< Base class alias for readability. */

    public:
        using BASE::BASE; /**< Inherits constructors from SearchBase. */

        /**
         * @brief Performs a single iteration of the upward bracketing search.
         * @details This method expands the search bounds upwards if the current bounds do not bracket a root.
         *          The expansion factor controls the rate at which the bounds are expanded.
         */
        void iterate()
        {
            const auto& bounds = BASE::current();

            // Check if current bounds already bracket a root
            if (BASE::evaluate(bounds.first) * BASE::evaluate(bounds.second) < 0.0) return;

            // Expand the bounds upwards
            auto newBounds   = bounds;
            newBounds.first  = bounds.second;
            newBounds.second = bounds.second + (bounds.second - bounds.first) * BASE::ratio();
            BASE::setBounds(newBounds);
        }
    };

    /**
     * @brief Deduction guides for BracketSearchUp class.
     * Allows the type of BracketSearchUp class to be deduced from the constructor parameters.
     */
    template< typename FN, typename ARG_T, typename FACTOR_T = ARG_T >
    requires IsFloatInvocable< FN > && IsFloat< ARG_T > && IsFloat< FACTOR_T >
    BracketSearchUp(FN, std::initializer_list< ARG_T >, FACTOR_T factor = std::numbers::phi) -> BracketSearchUp< FN, ARG_T >;

    template< typename FN, typename BOUNDS_T, typename FACTOR_T = StructCommonType_t< BOUNDS_T > >
    requires IsFloatInvocable< FN > && IsFloatStruct< BOUNDS_T >
    BracketSearchUp(FN, BOUNDS_T, FACTOR_T factor = std::numbers::phi) -> BracketSearchUp< FN, StructCommonType_t< BOUNDS_T > >;

    // =================================================================================================================
    //
    //  ad88888ba                                                   88           88888888ba,
    // d8"     "8b                                                  88           88      `"8b
    // Y8,                                                          88           88        `8b
    // `Y8aaaaa,     ,adPPYba,  ,adPPYYba,  8b,dPPYba,   ,adPPYba,  88,dPPYba,   88         88   ,adPPYba,   8b      db      d8  8b,dPPYba,
    //   `"""""8b,  a8P_____88  ""     `Y8  88P'   "Y8  a8"     ""  88P'    "8a  88         88  a8"     "8a  `8b    d88b    d8'  88P'   `"8a
    //         `8b  8PP"""""""  ,adPPPPP88  88          8b          88       88  88         8P  8b       d8   `8b  d8'`8b  d8'   88       88
    // Y8a     a8P  "8b,   ,aa  88,    ,88  88          "8a,   ,aa  88       88  88      .a8P   "8a,   ,a8"    `8bd8'  `8bd8'    88       88
    //  "Y88888P"    `"Ybbd8"'  `"8bbdP"Y8  88           `"Ybbd8"'  88       88  88888888Y"'     `"YbbdP"'       YP      YP      88       88
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
    template< IsFloatInvocable FN, IsFloat ARG_T = double >
    class BracketSearchDown final : public detail::SearchBase< BracketSearchDown< FN, ARG_T >, FN, ARG_T >
    {
        using BASE = detail::SearchBase< BracketSearchDown< FN, ARG_T >, FN, ARG_T >; /**< Base class alias for readability. */

    public:
        using BASE::BASE; /**< Inherits constructors from SearchBase. */

        /**
         * @brief Performs a single iteration of the downward bracketing search.
         * @details This method expands the search bounds downwards if the current bounds do not bracket a root.
         *          The expansion factor controls the rate at which the bounds are contracted.
         */
        void iterate()
        {
            const auto& bounds = BASE::current();

            // Check if current bounds already bracket a root
            if (BASE::evaluate(bounds.first) * BASE::evaluate(bounds.second) < 0.0) return;

            // Expand the bounds downwards
            auto newBounds   = bounds;
            newBounds.second = bounds.first;
            newBounds.first  = bounds.first - (bounds.second - bounds.first) * BASE::ratio();
            BASE::setBounds(newBounds);
        }
    };

    /**
     * @brief Deduction guides for BracketSearchDown class.
     * Allows the type of BracketSearchDown class to be deduced from the constructor parameters.
     */
    template< typename FN, typename ARG_T, typename FACTOR_T = ARG_T >
    requires IsFloatInvocable< FN > && IsFloat< ARG_T > && IsFloat< FACTOR_T >
    BracketSearchDown(FN, std::initializer_list< ARG_T >, FACTOR_T factor = std::numbers::phi) -> BracketSearchDown< FN, ARG_T >;

    template< typename FN, typename BOUNDS_T, typename FACTOR_T = StructCommonType_t< BOUNDS_T > >
    requires IsFloatInvocable< FN > && IsFloatStruct< BOUNDS_T >
    BracketSearchDown(FN, BOUNDS_T, FACTOR_T factor = std::numbers::phi) -> BracketSearchDown< FN, StructCommonType_t< BOUNDS_T > >;

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
    template< IsFloatInvocable FN, IsFloat ARG_T = double >
    class BracketExpandUp final : public detail::SearchBase< BracketExpandUp< FN, ARG_T >, FN, ARG_T >
    {
        using BASE = detail::SearchBase< BracketExpandUp< FN, ARG_T >, FN, ARG_T >; /**< Base class alias for readability. */

    public:
        using BASE::BASE; /**< Inherits constructors from SearchBase. */

        /**
         * @brief Performs a single iteration of the upward bracket expansion.
         * @details This method expands the upper bound upwards if the current bounds do not bracket a root.
         *          The expansion factor controls the rate at which the upper bound is expanded.
         */
        void iterate()
        {
            const auto& bounds = BASE::current();

            // Check if current bounds already bracket a root
            if (BASE::evaluate(bounds.first) * BASE::evaluate(bounds.second) < 0.0) return;

            // Expand the upper bound upwards
            auto newBounds   = bounds;
            newBounds.second = bounds.second + (bounds.second - bounds.first) * BASE::ratio();
            BASE::setBounds(newBounds);
        }
    };

    /**
     * @brief Deduction guides for BracketExpandUp class.
     * Allows the type of BracketExpandUp class to be deduced from the constructor parameters.
     */
    template< typename FN, typename ARG_T, typename FACTOR_T = ARG_T >
    requires IsFloatInvocable< FN > && IsFloat< ARG_T > && IsFloat< FACTOR_T >
    BracketExpandUp(FN, std::initializer_list< ARG_T >, FACTOR_T factor = std::numbers::phi) -> BracketExpandUp< FN, ARG_T >;

    template< typename FN, typename BOUNDS_T, typename FACTOR_T = StructCommonType_t< BOUNDS_T > >
    requires IsFloatInvocable< FN > && IsFloatStruct< BOUNDS_T >
    BracketExpandUp(FN, BOUNDS_T, FACTOR_T factor = std::numbers::phi) -> BracketExpandUp< FN, StructCommonType_t< BOUNDS_T > >;

    // =================================================================================================================
    //
    // 88888888888                                                              88  88888888ba,
    // 88                                                                       88  88      `"8b
    // 88                                                                       88  88        `8b
    // 88aaaaa      8b,     ,d8  8b,dPPYba,   ,adPPYYba,  8b,dPPYba,    ,adPPYb,88  88         88   ,adPPYba,   8b      db      d8
    // 8b,dPPYba, 88"""""       `Y8, ,8P'   88P'    "8a  ""     `Y8  88P'   `"8a  a8"    `Y88  88         88  a8"     "8a  `8b    d88b d8'
    // 88P'   `"8a 88              )888(     88       d8  ,adPPPPP88  88       88  8b       88  88         8P  8b       d8   `8b  d8'`8b d8'
    // 88       88 88            ,d8" "8b,   88b,   ,a8"  88,    ,88  88       88  "8a,   ,d88  88      .a8P   "8a,   ,a8"    `8bd8'  `8bd8'
    // 88       88 88888888888  8P'     `Y8  88`YbbdP"'   `"8bbdP"Y8  88       88   `"8bbdP"Y8  88888888Y"'     `"YbbdP"'       YP      YP
    // 88       88
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
    template< IsFloatInvocable FN, IsFloat ARG_T = double >
    class BracketExpandDown final : public detail::SearchBase< BracketExpandDown< FN, ARG_T >, FN, ARG_T >
    {
        using BASE = detail::SearchBase< BracketExpandDown< FN, ARG_T >, FN, ARG_T >; /**< Base class alias for readability. */

    public:
        using BASE::BASE; /**< Inherits constructors from SearchBase. */

        /**
         * @brief Performs a single iteration of the downward bracket expansion.
         * @details This method expands the lower bound downwards if the current bounds do not bracket a root.
         *          The expansion factor controls the rate at which the lower bound is expanded.
         */
        void iterate()
        {
            const auto& bounds = BASE::current();

            // Check if current bounds already bracket a root
            if (BASE::evaluate(bounds.first) * BASE::evaluate(bounds.second) < 0.0) return;

            // Expand the lower bound downwards
            auto newBounds  = bounds;
            newBounds.first = bounds.first - (bounds.second - bounds.first) * BASE::ratio();
            BASE::setBounds(newBounds);
        }
    };

    /**
     * @brief Deduction guides for BracketExpandDown class.
     * Allows the type of BracketExpandDown class to be deduced from the constructor parameters.
     */
    template< typename FN, typename ARG_T, typename FACTOR_T = ARG_T >
    requires IsFloatInvocable< FN > && IsFloat< ARG_T > && IsFloat< FACTOR_T >
    BracketExpandDown(FN, std::initializer_list< ARG_T >, FACTOR_T factor = std::numbers::phi) -> BracketExpandDown< FN, ARG_T >;

    template< typename FN, typename BOUNDS_T, typename FACTOR_T = StructCommonType_t< BOUNDS_T > >
    requires IsFloatInvocable< FN > && IsFloatStruct< BOUNDS_T >
    BracketExpandDown(FN, BOUNDS_T, FACTOR_T factor = std::numbers::phi) -> BracketExpandDown< FN, StructCommonType_t< BOUNDS_T > >;

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
    template< IsFloatInvocable FN, IsFloat ARG_T = double >
    class BracketExpandOut final : public detail::SearchBase< BracketExpandOut< FN, ARG_T >, FN, ARG_T >
    {
        using BASE = detail::SearchBase< BracketExpandOut< FN, ARG_T >, FN, ARG_T >; /**< Base class alias for readability. */

    public:
        using BASE::BASE; /**< Inherits constructors from SearchBase. */

        /**
         * @brief Performs a single iteration of the outward bracket expansion.
         * @details This method expands both the lower and upper bounds outwards symmetrically if the current
         *          bounds do not bracket a root. The expansion factor controls the rate at which the bounds are expanded.
         */
        void iterate()
        {
            const auto& bounds = BASE::current();

            // Check if current bounds already bracket a root
            if (BASE::evaluate(bounds.first) * BASE::evaluate(bounds.second) < 0.0) return;

            // Expand the bounds symmetrically outwards
            auto newBounds   = bounds;
            newBounds.first  = bounds.first - (bounds.second - bounds.first) * BASE::ratio() / 2.0;
            newBounds.second = bounds.second + (bounds.second - bounds.first) * BASE::ratio() / 2.0;
            BASE::setBounds(newBounds);
        }
    };

    /**
     * @brief Deduction guides for BracketExpandOut class.
     * Allows the type of BracketExpandOut class to be deduced from the constructor parameters.
     */
    template< typename FN, typename ARG_T, typename FACTOR_T = ARG_T >
    requires IsFloatInvocable< FN > && IsFloat< ARG_T > && IsFloat< FACTOR_T >
    BracketExpandOut(FN, std::initializer_list< ARG_T >, FACTOR_T factor = std::numbers::phi) -> BracketExpandOut< FN, ARG_T >;

    template< typename FN, typename BOUNDS_T, typename FACTOR_T = StructCommonType_t< BOUNDS_T > >
    requires IsFloatInvocable< FN > && IsFloatStruct< BOUNDS_T >
    BracketExpandOut(FN, BOUNDS_T, FACTOR_T factor = std::numbers::phi) -> BracketExpandOut< FN, StructCommonType_t< BOUNDS_T > >;

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
    template< IsFloatInvocable FN, IsFloat ARG_T = double >
    class BracketSubdivide final : public detail::SearchBase< BracketSubdivide< FN, ARG_T >, FN, ARG_T >
    {
        using BASE = detail::SearchBase< BracketSubdivide< FN, ARG_T >, FN, ARG_T >; /**< Base class alias for readability. */

    public:
        using BASE::BASE; /**< Inherits constructors from SearchBase. */

        /**
         * @brief Performs a single iteration of the bracket subdivision search.
         * @details This method subdivides the current bounds into smaller segments based on the factor,
         *          attempting to find a segment where the function changes sign, indicating a bracket.
         *          If no such segment is found, the factor is doubled to increase the search range.
         */
        void iterate()
        {
            const auto& bounds = BASE::current();
            if (BASE::evaluate(bounds.first) * BASE::evaluate(bounds.second) < 0.0) return;

            size_t factor = std::ceil(BASE::ratio());
            auto   diff   = (bounds.second - bounds.first) / factor;
            auto   lower  = bounds.first;
            auto   upper  = bounds.first + diff;
            for (size_t i = 0; i < factor; ++i) {
                if (BASE::evaluate(lower) * BASE::evaluate(upper) < 0.0) {
                    BASE::setBounds({ lower, upper });
                    return;
                }
                lower = upper;
                upper += diff;
            }

            BASE::setRatio(BASE::ratio() * 2.0);    // Increase the factor to expand the search range
        }
    };

    /**
     * @brief Deduction guides for BracketSubdivide class.
     * Allows the type of BracketSubdivide class to be deduced from the constructor parameters.
     */
    template< typename FN, typename ARG_T, typename FACTOR_T = ARG_T >
    requires IsFloatInvocable< FN > && IsFloat< ARG_T > && IsFloat< FACTOR_T >
    BracketSubdivide(FN, std::initializer_list< ARG_T >, FACTOR_T factor = std::numbers::phi) -> BracketSubdivide< FN, ARG_T >;

    template< typename FN, typename BOUNDS_T, typename FACTOR_T = StructCommonType_t< BOUNDS_T > >
    requires IsFloatInvocable< FN > && IsFloatStruct< BOUNDS_T >
    BracketSubdivide(FN, BOUNDS_T, FACTOR_T factor = std::numbers::phi) -> BracketSubdivide< FN, StructCommonType_t< BOUNDS_T > >;

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

    namespace detail
    {
        /**
         * @brief Implements a generic search solver function template for bracketing searchers.
         *
         * This function template, `search_impl`, provides a generic implementation for search operations
         * using various bracketing searcher algorithms. It is designed to work with solvers that conform
         * to the requirements of bracketing searchers, such as having a defined `IsBracketingSearcher` static member,
         * initialization, and iteration methods. The function handles initialization, iteration, and
         * checking for a successful bracketing of the root, returning the result along with any potential errors encountered
         * during the searching process.
         *
         * @tparam SOLVER The type of the solver to be used in search operations. Must conform to the bracketing searcher concept.
         */
        template< typename SOLVER >
        requires SOLVER::IsBracketingSearcher
        auto search_impl(SOLVER solver, IsFloatStruct auto bounds, IsFloat auto ratio, std::integral auto maxiter)
        {
            using ET = RootErrorImpl< decltype(bounds) >;    /**< Type for error handling. */
            using RT = tl::expected< decltype(bounds), ET >; /**< Type for the function return value. */

            solver.init(bounds, ratio);
            auto                      curBounds = solver.current();
            RT                        result    = curBounds;
            typename SOLVER::RESULT_T eval_lower;
            typename SOLVER::RESULT_T eval_upper;

            // Check for NaN or Inf in the initial bounds.
            if (!std::isfinite(solver.evaluate(curBounds.first)) || !std::isfinite(solver.evaluate(curBounds.second))) {
                result = tl::make_unexpected(ET("Invalid initial brackets!", RootErrorType::NumericalError, result.value()));
                return result;
            }

            int iter = 1;
            while (true) {
                curBounds  = solver.current();
                eval_lower = solver.evaluate(curBounds.first);
                eval_upper = solver.evaluate(curBounds.second);

                // Check for non-finite results.
                if (!std::isfinite(curBounds.first) || !std::isfinite(curBounds.second) || !std::isfinite(eval_lower) ||
                    !std::isfinite(eval_upper))
                {
                    result = tl::make_unexpected(ET("Non-finite result!", RootErrorType::NumericalError, curBounds, iter));
                    break;
                }

                // Check if the root is bracketed by the bounds. If yes, return the bounds.
                if (eval_lower * eval_upper < 0.0) {
                    result = curBounds;
                    break;
                }

                // Check for maximum number of iterations.
                if (iter >= maxiter) {
                    result = tl::make_unexpected(
                        ET("Maximum number of iterations exceeded!", RootErrorType::MaxIterationsExceeded, curBounds, iter));
                    break;
                }

                // Otherwise, iterate the solver
                solver.iterate();
                ++iter;
            }

            return result;
        }
    }    // namespace detail

    /**
     * @brief Defines a high-level search function template `search` using bracketing searchers.
     *
     * The `search` function template provides a convenient interface for performing search operations
     * using various bracketing searcher algorithms. It abstracts the creation and configuration of the
     * solver instance and then delegates the actual search process to `search_impl`. This function is templated
     * to accept a solver type, the function, bounds for the search, a search factor for controlling the search
     * process, and a maximum number of iterations. It supports different types of bracketing searchers, making it
     * versatile for various search needs.
     *
     * @tparam SOLVER_T The template class of the solver to be used. Must be a valid bracketing searcher type.
     * @tparam FN_T The type of the function for which the search is being conducted.
     * @tparam STRUCT_T The struct type holding the bounds for the search.
     * @tparam FACTOR_T The type of the search factor, defaulted based on STRUCT_T.
     * @tparam ITER_T The type of the maximum iterations count, defaulted to int.
     */
    template< template< typename, typename > class SOLVER_T,
              IsFloatInvocable FN_T,
              IsFloatStruct    STRUCT_T,
              IsFloat          FACTOR_T = StructCommonType_t< STRUCT_T >,
              std::integral    ITER_T   = int >
    auto search(FN_T     function,
                STRUCT_T bounds,
                FACTOR_T ratio   = std::numbers::phi,
                ITER_T   maxiter = iterations< StructCommonType_t< STRUCT_T > >())
    {
        auto [lo, hi] = bounds; /**< Extract lower and upper bounds from the struct. */

        using ARG_T = StructCommonType_t< STRUCT_T >;            /**< Common type for the bounds. */
        auto solver = SOLVER_T< FN_T, ARG_T >(function, bounds); /**< Instantiates the solver with the given function. */

        // Delegates the search process to search_impl, passing in the solver and other parameters.
        return detail::search_impl(solver, std::pair< ARG_T, ARG_T > { lo, hi }, ratio, maxiter);
    }

    /**
     * @brief Extends the high-level search function template `search` for initializer list bounds.
     *
     * This version of the `search` function template allows for specifying the bounds using an initializer list.
     * It is particularly useful when the bounds are known at compile time or for concise inline specifications.
     * The function checks the size of the initializer list to ensure exactly two elements are provided for the
     * bounds. It then creates a solver instance and delegates the search process to `search_impl`.
     *
     * @tparam SOLVER_T The template class of the solver to be used. Must be a valid bracketing searcher type.
     * @tparam FN_T The type of the function for which the search is being conducted.
     * @tparam ARG_T The type of the bounds and the argument to the function.
     * @tparam FACTOR_T The type of the search factor, defaulted based on ARG_T.
     * @tparam ITER_T The type of the maximum iterations count, defaulted to int.
     */
    template< template< typename, typename > class SOLVER_T,
              IsFloatInvocable FN_T,
              IsFloat          ARG_T,
              IsFloat          FACTOR_T = ARG_T,
              std::integral    ITER_T   = int,
              size_t           N >
    requires(N == 2)
    auto search(FN_T function, const ARG_T (&bounds)[N], FACTOR_T ratio = std::numbers::phi, ITER_T maxiter = iterations< ARG_T >())
    {
        auto solver =
            SOLVER_T< FN_T, ARG_T >(function, std::pair(bounds[0], bounds[1])); /**< Instantiates the solver with the given function. */

        // Delegates the search process to search_impl, passing in the solver and other parameters.
        return detail::search_impl(solver, std::pair(bounds[0], bounds[1]), ratio, maxiter);
    }
}    // namespace nxx::roots

#endif    // NUMERIXX_ROOTSEARCHING_HPP

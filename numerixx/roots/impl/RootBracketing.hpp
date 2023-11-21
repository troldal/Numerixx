/*
    888b      88  88        88  88b           d88  88888888888  88888888ba   88  8b        d8  8b        d8
    8888b     88  88        88  888b         d888  88           88      "8b  88   Y8,    ,8P    Y8,    ,8P
    88 `8b    88  88        88  88`8b       d8'88  88           88      ,8P  88    `8b  d8'      `8b  d8'
    88  `8b   88  88        88  88 `8b     d8' 88  88aaaaa      88aaaaaa8P'  88      Y88P          Y88P
    88   `8b  88  88        88  88  `8b   d8'  88  88"""""      88""""88'    88      d88b          d88b
    88    `8b 88  88        88  88   `8b d8'   88  88           88    `8b    88    ,8P  Y8,      ,8P  Y8,
    88     `8888  Y8a.    .a8P  88    `888'    88  88           88     `8b   88   d8'    `8b    d8'    `8b
    88      `888   `"Y8888Y"'   88     `8'     88  88888888888  88      `8b  88  8P        Y8  8P        Y8

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

#ifndef NUMERIXX_ROOTBRACKETING_HPP
#define NUMERIXX_ROOTBRACKETING_HPP

// ===== Numerixx Includes
#include "RootCommon.hpp"
#include <Constants.hpp>
#include <Deriv.hpp>

// ===== Standard Library Includes
#include <algorithm>
#include <array>
#include <span>

namespace nxx::roots
{
    // =================================================================================================================
    //
    //  88888888ba                                       88                             88
    //  88      "8b                                      88                      ,d     ""
    //  88      ,8P                                      88                      88
    //  88aaaaaa8P'  8b,dPPYba,  ,adPPYYba,   ,adPPYba,  88   ,d8   ,adPPYba,  MM88MMM  88  8b,dPPYba,    ,adPPYb,d8
    //  88""""""8b,  88P'   "Y8  ""     `Y8  a8"     ""  88 ,a8"   a8P_____88    88     88  88P'   `"8a  a8"    `Y88
    //  88      `8b  88          ,adPPPPP88  8b          8888[     8PP"""""""    88     88  88       88  8b       88
    //  88      a8P  88          88,    ,88  "8a,   ,aa  88`"Yba,  "8b,   ,aa    88,    88  88       88  "8a,   ,d88
    //  88888888P"   88          `"8bbdP"Y8   `"Ybbd8"'  88   `Y8a  `"Ybbd8"'    "Y888  88  88       88   `"YbbdP"Y8
    //                                                                                                    aa,    ,88
    //                                                                                                     "Y8bbdP"
    //
    // =================================================================================================================

    namespace impl
    {
        /**
         * @brief A base class for bracketing root-finding algorithms.
         *
         * This class provides a generic interface for bracketing root-finding algorithms. The actual algorithm
         * is implemented in a derived class, which is passed as a template argument (POLICY).
         *
         * @tparam DERIVED The derived class implementing the specific root-finding algorithm.
         * @requires POLICY must be invocable with a floating point type as its argument type.
         */
        template<typename DERIVED, typename FUNCTION_T, typename ARG_T>
            requires std::same_as< typename BracketingTraits< DERIVED >::FUNCTION_T, FUNCTION_T > &&
                     nxx::IsFloatInvocable< FUNCTION_T > &&
                     nxx::FloatingPoint< ARG_T > &&
                     nxx::FloatingPoint< typename BracketingTraits< DERIVED >::RETURN_T >
        class BracketingBase
        {
            /*
             * Friend declarations.
             */
            friend DERIVED;

        public:
            using RESULT_T = std::invoke_result_t< FUNCTION_T, ARG_T >;
            using BOUNDS_T = std::pair< ARG_T, ARG_T >;

        protected:
            /**
             * @brief Default destructor.
             */
            ~BracketingBase() = default;

        public:
            FUNCTION_T m_func{};                 /**< The function object to find the root for. */
            BOUNDS_T   m_bounds{};               /**< Holds the current bounds around the root. */
            bool       m_isInitialized{ false }; /**< Indicates whether the solver has been initialized. */

            /**
             * @brief Constructor, taking a function object as an argument.
             *
             * @param objective The function object to find the root for.
             * @note Constructor is private to avoid direct usage by clients.
             */
            explicit BracketingBase(FUNCTION_T objective)
                : m_func{ std::move(objective) } {}

            /**
             * @brief Constructor, taking a function object and std::initialilzer_list with the bounds as arguments.
             *
             * @param objective The function object to find the root for.
             * @param bounds An std::initializer_list object holding the initial bounds around the root. The root must be contained
             * inside these bounds. The list must contain exactly two elements, which will be interpreted as the lower and
             * upper bounds, respectively.
             * @note Constructor is private to avoid direct usage by clients.
             */
            template<typename T>
                requires nxx::FloatingPoint< T >
            BracketingBase(FUNCTION_T objective, std::initializer_list< T > bounds)
                : m_func{ std::move(objective) } { init(bounds); }

            /**
             * @brief Constructor, taking a function object and a container with the bounds as arguments.
             *
             * @param objective The function object to find the root for.
             * @param bounds A container holding the initial bounds around the root. The root must be contained
             * inside these bounds. The container must contain exactly two elements, which will be interpreted as the lower and
             * upper bounds, respectively.
             * @note Constructor is private to avoid direct usage by clients.
             */
            template<IsContainer CONT_T>
                requires nxx::FloatingPoint< typename CONT_T::value_type >
            BracketingBase(FUNCTION_T objective, CONT_T bounds)
                : m_func{ std::move(objective) } { init(bounds); }

            /**
             * @brief Constructor, taking a function object and a struct with the bounds as arguments.
             *
             * @param objective The function object to find the root for.
             * @param bounds A struct holding the initial bounds around the root. The root must be contained
             * inside these bounds. The struct must support structured bindings to provide two values: lower and upper bounds.
             * Examples of supported types include pairs, tuples, or custom structs with structured bindings support.
             * @note Constructor is private to avoid direct usage by clients.
             */
            template<IsFloatStruct STRUCT_T>
            BracketingBase(FUNCTION_T objective, STRUCT_T bounds)
                : m_func{ std::move(objective) } { init(bounds); }

            /**
             * @brief Sets the current bounds around the root.
             *
             * @param bounds A std::pair object holding the lower and upper bounds, respectively.
             */
            void setBounds(const BOUNDS_T& bounds)
            {
                if (!m_isInitialized) throw NumerixxError("Solver has not been initialized!");
                auto [lower, upper] = bounds;
                static_assert(nxx::FloatingPoint< decltype(lower) >);
                m_bounds = BOUNDS_T{ lower, upper };
            }

            /**
             * @brief Copy constructor.
             *
             * @param other Another BracketingBase object to be copied.
             */
            BracketingBase(const BracketingBase& other) = default;

            /**
             * @brief Move constructor.
             *
             * @param other Another BracketingBase object to be moved.
             */
            BracketingBase(BracketingBase&& other) noexcept = default;

            /**
             * @brief Copy assignment operator.
             *
             * @param other Another BracketingBase object to be copied.
             * @return A reference to the assigned object.
             */
            BracketingBase& operator=(const BracketingBase& other) = default;

            /**
             * @brief Move assignment operator.
             *
             * @param other Another BracketingBase object to be moved.
             * @return A reference to the assigned object.
             */
            BracketingBase& operator=(BracketingBase&& other) noexcept = default;

            /**
             * @brief Initializes the solver with the initial bounds.
             *
             * @param bounds An std::initializer_list object holding the initial bounds around the root. The root must be contained
             * inside these bounds. The list must contain exactly two elements, which will be interpreted as the lower and
             * upper bounds, respectively.
             */
            template<typename T>
                requires nxx::FloatingPoint< T >
            void init(std::initializer_list< T > bounds)
            {
                m_isInitialized = true;
                if (bounds.size() != 2) throw NumerixxError("Container must contain exactly two elements!");
                auto bnds = std::span(bounds.begin(), bounds.end());
                setBounds(BOUNDS_T{ bnds.front(), bnds.back() });
            }

            /**
             * @brief Initializes the solver with the initial bounds.
             *
             * @param bounds A container holding the initial bounds around the root. The root must be contained
             * inside these bounds. The container must contain exactly two elements, which will be interpreted as the lower and
             * upper bounds, respectively.
             */
            template<IsContainer CONT_T>
                requires nxx::FloatingPoint< typename CONT_T::value_type >
            void init(CONT_T bounds)
            {
                m_isInitialized = true;
                if (bounds.size() != 2) throw NumerixxError("Container must contain exactly two elements!");
                setBounds(BOUNDS_T{ bounds.front(), bounds.back() });
            }

            /**
             * @brief Initializes the solver with the initial bounds.
             *
             * @param bounds A struct holding the initial bounds around the root. The root must be contained
             * inside these bounds. The struct must support structured bindings to provide two values: lower and upper bounds.
             * Examples of supported types include pairs, tuples, or custom structs with structured bindings support.
             */
            template<IsFloatStruct STRUCT_T>
            void init(STRUCT_T bounds)
            {
                m_isInitialized     = true;
                auto [lower, upper] = bounds;
                setBounds(BOUNDS_T{ lower, upper });
            }

            /**
             * @brief Resets the solver to its initial state. To be called before reusing the solver.
             *
             * @note After calling this function, the solver must be initialized again before it can be used.
             * @warning If the solver is used without initializing, the behavior is undefined.
             */
            void reset() { m_isInitialized = false; }

            /**
             * @brief Evaluates the function to solve at a given point.
             *
             * Passes the given argument to the function object to solve and returns the result of the evaluation.
             * The return type will be the same as the return type of the given function object.
             *
             * @param value The value at which to evaluate the function.
             * @return The result of the evaluation.
             */
            RESULT_T evaluate(ARG_T value) { return m_func(value); }

            /**
             * @brief Returns the current bounds around the root.
             *
             * Every time an iteration is executed, the bounds will narrow. This function returns a const reference
             * to the current bounds as a std::pair. The value type of the bounds is the same as the return type
             * of the function object.
             *
             * @return A const reference to the current bounds.
             */
            const BOUNDS_T& bounds() const
            {
                if (!m_isInitialized) throw NumerixxError("Solver has not been initialized!");
                return m_bounds;
            }
        };
    } // namespace impl


    // =================================================================================================================
    //
    //  88888888ba   88           88           88
    //  88      "8b  ""           88           88
    //  88      ,8P               88           88
    //  88aaaaaa8P'  88   ,adPPYb,88   ,adPPYb,88   ,adPPYba,  8b,dPPYba,  ,adPPYba,
    //  88""""88'    88  a8"    `Y88  a8"    `Y88  a8P_____88  88P'   "Y8  I8[    ""
    //  88    `8b    88  8b       88  8b       88  8PP"""""""  88           `"Y8ba,
    //  88     `8b   88  "8a,   ,d88  "8a,   ,d88  "8b,   ,aa  88          aa    ]8I
    //  88      `8b  88   `"8bbdP"Y8   `"8bbdP"Y8   `"Ybbd8"'  88          `"YbbdP"'
    //
    // =================================================================================================================

    /**
     * @brief Implements Ridder's method for root-finding.
     *
     * This class implements Ridder's method, a bracketing root-finding algorithm, as a derived class of
     * impl::BracketingBase. It inherits the base functionality from impl::BracketingBase and adds the
     * specific algorithm implementation for Ridder's method.
     *
     * @tparam FN The function object type for which to find the root.
     * @requires FN must be invocable with a double as its argument type.
     */
    template<IsFloatInvocable FN, FloatingPoint ARG_T = double>
    class Ridder final : public impl::BracketingBase< Ridder< FN, ARG_T >, FN, ARG_T >
    {
        /*
         * Private alias declarations.
         */
        using BASE = impl::BracketingBase< Ridder< FN, ARG_T >, FN, ARG_T >;

    public:
        using BASE::BASE;

        /**
         * @brief Perform one iteration of Ridder's method.
         *
         * This function implements the main algorithm of Ridder's method for root-finding. It updates the
         * bounds around the root during each iteration, gradually narrowing the search interval.
         */
        void iterate()
        {
            const auto& bounds = BASE::bounds();
            using RT = std::invoke_result_t< FN, decltype(bounds.first) >;

            using std::abs;
            using std::pow;
            using std::sqrt;

            const RT& x_lo = bounds.first;
            const RT& x_hi = bounds.second;
            RT        f_lo = BASE::evaluate(x_lo);
            RT        f_hi = BASE::evaluate(x_hi);

            RT x_mid;
            RT f_mid;

            RT x_new;
            RT f_new;

            // ===== Calculate new bounds
            x_mid    = (x_lo + x_hi) / 2.0;
            f_mid    = BASE::evaluate(x_mid);
            int sign = ((f_lo - f_hi) < 0.0 ? -1 : 1);
            x_new    = x_mid + (x_mid - x_lo) * ((sign * f_mid) / sqrt(f_mid * f_mid - f_lo * f_hi));
            f_new    = BASE::evaluate(x_new);

            // ===== General case: The root is between x_mid and x_new
            if (f_mid * f_new < 0.0) {
                if (x_mid < x_new)
                    BASE::setBounds({ x_mid, x_new });
                else
                    BASE::setBounds({ x_new, x_mid });
            }

            // ===== Degenerate cases: The root is between x_new and either x_lo or x_hi
            if (f_hi * f_new < 0.0) {
                if (x_hi < x_new)
                    BASE::setBounds({ x_hi, x_new });
                else
                    BASE::setBounds({ x_new, x_hi });
            }

            else {
                if (x_lo < x_new)
                    BASE::setBounds({ x_lo, x_new });
                else
                    BASE::setBounds({ x_new, x_lo });
            }
        }
    };

    /**
     * @brief Deduction guide for the Ridder class.
     * @tparam FN The type of the function object for which to find the root. The function must be invocable
     * with a double argument.
     */
    template<IsFloatInvocable FN>
    Ridder(FN func) -> Ridder< decltype(func) >;

    /**
     * @brief Deduction guide for the Ridder class.
     * @tparam FN The type of the function object for which to find the root. The function must be invocable
     * with a double argument.
     * @tparam ARG_T The type of the bounds. Must be a floating point type.
     */
    template<IsFloatInvocable FN, FloatingPoint ARG_T>
    Ridder(FN func, std::initializer_list< ARG_T > bounds) -> Ridder< decltype(func), ARG_T >;

    /**
     * @brief Deduction guide for the Ridder class.
     * @tparam FN The type of the function object for which to find the root. The function must be invocable
     * with a double argument.
     * @tparam CONT_T The type of the container holding the bounds. Must be a container of floating point types.
     */
    template<IsFloatInvocable FN, IsContainer CONT_T>
    Ridder(FN func, CONT_T bounds) -> Ridder< decltype(func), typename CONT_T::value_type >;

    /**
     * @brief Deduction guide for the Ridder class.
     * @tparam FN The type of the function object for which to find the root. The function must be invocable
     * with a double argument.
     * @tparam BOUNDS_T The type of the struct holding the bounds. Must be a struct with structured bindings support
     * and with floating point types as its members.
     */
    template<IsFloatInvocable FN, IsFloatStruct BOUNDS_T>
    Ridder(FN func, BOUNDS_T bounds) -> Ridder< decltype(func), StructCommonType_t< BOUNDS_T > >;


    // =================================================================================================================
    //
    //  88888888ba   88                                              88
    //  88      "8b  ""                                       ,d     ""
    //  88      ,8P                                           88
    //  88aaaaaa8P'  88  ,adPPYba,   ,adPPYba,   ,adPPYba,  MM88MMM  88   ,adPPYba,   8b,dPPYba,
    //  88""""""8b,  88  I8[    ""  a8P_____88  a8"     ""    88     88  a8"     "8a  88P'   `"8a
    //  88      `8b  88   `"Y8ba,   8PP"""""""  8b            88     88  8b       d8  88       88
    //  88      a8P  88  aa    ]8I  "8b,   ,aa  "8a,   ,aa    88,    88  "8a,   ,a8"  88       88
    //  88888888P"   88  `"YbbdP"'   `"Ybbd8"'   `"Ybbd8"'    "Y888  88   `"YbbdP"'   88       88
    //
    // =================================================================================================================

    /**
     * @brief Implements the bisection method for root-finding.
     *
     * This class implements the bisection method, a bracketing root-finding algorithm, as a derived
     * class of impl::BracketingBase. It inherits the base functionality from impl::BracketingBase and
     * adds the specific algorithm implementation for the bisection method.
     *
     * @tparam FN The function object type for which to find the root.
     * @requires FN must be invocable with a double as its argument type.
     */
    template<IsFloatInvocable FN, FloatingPoint ARG_T = double>
    class Bisection final : public impl::BracketingBase< Bisection< FN, ARG_T >, FN, ARG_T >
    {
        /*
         * Private alias declarations.
         */
        using BASE = impl::BracketingBase< Bisection< FN, ARG_T >, FN, ARG_T >;

    public:
        using BASE::BASE;

        /**
         * @brief Perform one iteration of the bisection method.
         *
         * This function implements the main algorithm of the bisection method for root-finding. It
         * updates the bounds around the root during each iteration, gradually narrowing the search
         * interval.
         */
        void iterate()
        {
            const auto& bounds = BASE::bounds();
            using RT = std::invoke_result_t< FN, decltype(bounds.first) >;

            if (RT root = (bounds.first + bounds.second) / 2.0; BASE::evaluate(bounds.first) * BASE::evaluate(root) < 0.0)
                BASE::setBounds({ bounds.first, root });
            else
                BASE::setBounds({ root, bounds.second });
        }
    };

    /**
     * @brief Deduction guide for the Bisection class.
     * @tparam FN The type of the function object for which to find the root. The function must be invocable
     * with a double argument.
     */
    template<IsFloatInvocable FN>
    Bisection(FN func) -> Bisection< decltype(func) >;

    /**
     * @brief Deduction guide for the Bisection class.
     * @tparam FN The type of the function object for which to find the root. The function must be invocable
     * with a double argument.
     * @tparam ARG_T The type of the bounds. Must be a floating point type.
     */
    template<IsFloatInvocable FN, FloatingPoint ARG_T>
    Bisection(FN func, std::initializer_list< ARG_T > bounds) -> Bisection< decltype(func), ARG_T >;

    /**
     * @brief Deduction guide for the Bisection class.
     * @tparam FN The type of the function object for which to find the root. The function must be invocable
     * with a double argument.
     * @tparam CONT_T The type of the container holding the bounds. Must be a container of floating point types.
     */
    template<IsFloatInvocable FN, IsContainer CONT_T>
    Bisection(FN func, CONT_T bounds) -> Bisection< decltype(func), typename CONT_T::value_type >;

    /**
     * @brief Deduction guide for the Bisection class.
     * @tparam FN The type of the function object for which to find the root. The function must be invocable
     * with a double argument.
     * @tparam BOUNDS_T The type of the struct holding the bounds. Must be a struct with structured bindings support
     * and with floating point types as its members.
     */
    template<IsFloatInvocable FN, IsFloatStruct BOUNDS_T>
    Bisection(FN func, BOUNDS_T bounds) -> Bisection< decltype(func), StructCommonType_t< BOUNDS_T > >;


    // =================================================================================================================
    //
    //  88888888ba                                        88              88888888888          88             88
    //  88      "8b                                       88              88                   88             ""
    //  88      ,8P                                       88              88                   88
    //  88aaaaaa8P'  ,adPPYba,   ,adPPYb,d8  88       88  88  ,adPPYYba,  88aaaaa  ,adPPYYba,  88  ,adPPYba,  88
    //  88""""88'   a8P_____88  a8"    `Y88  88       88  88  ""     `Y8  88"""""  ""     `Y8  88  I8[    ""  88
    //  88    `8b   8PP"""""""  8b       88  88       88  88  ,adPPPPP88  88       ,adPPPPP88  88   `"Y8ba,   88
    //  88     `8b  "8b,   ,aa  "8a,   ,d88  "8a,   ,a88  88  88,    ,88  88       88,    ,88  88  aa    ]8I  88
    //  88      `8b  `"Ybbd8"'   `"YbbdP"Y8   `"YbbdP'Y8  88  `"8bbdP"Y8  88       `"8bbdP"Y8  88  `"YbbdP"'  88
    //                           aa,    ,88
    //                            "Y8bbdP"
    //
    // =================================================================================================================

    /**
     * @brief Regula Falsi (False Position) method for root-finding.
     *
     * This class implements the Regula Falsi algorithm, also known as the False Position method, for
     * finding the root of a given function. It inherits from the BracketingBase class and provides
     * the specific implementation for the Regula Falsi method.
     *
     * @tparam FN The type of the function object for which to find the root. The function must be invocable
     * with a double argument.
     */
    template<IsFloatInvocable FN, FloatingPoint ARG_T = double>
    class RegulaFalsi final : public impl::BracketingBase< RegulaFalsi< FN, ARG_T >, FN, ARG_T >
    {
        /*
         * Private alias declarations.
         */
        using BASE = impl::BracketingBase< RegulaFalsi< FN, ARG_T >, FN, ARG_T >;

    public:
        using BASE::BASE;

        /**
         * @brief Perform one iteration of the Regula Falsi algorithm.
         *
         * This function implements the main algorithm of the Regula Falsi method for root-finding.
         * It updates the bounds around the root during each iteration, refining the search interval.
         */
        void iterate()
        {
            const auto& bounds = BASE::bounds();
            using RT = std::invoke_result_t< FN, decltype(bounds.first) >;

            RT f_lo = BASE::evaluate(bounds.first);
            RT f_hi = BASE::evaluate(bounds.second);

            RT root   = (bounds.first * f_hi - bounds.second * f_lo) / (f_hi - f_lo);
            RT f_root = BASE::evaluate(root);

            if (f_lo * f_root < 0.0) { BASE::setBounds({ bounds.first, root }); }
            else { BASE::setBounds({ root, bounds.second }); }
        }
    };

    /**
     * @brief Deduction guide for the RegulaFalsi class.
     * @tparam FN The type of the function object for which to find the root. The function must be invocable
     * with a double argument.
     */
    template<IsFloatInvocable FN>
    RegulaFalsi(FN func) -> RegulaFalsi< decltype(func) >;

    /**
     * @brief Deduction guide for the RegulaFalsi class.
     * @tparam FN The type of the function object for which to find the root. The function must be invocable
     * with a double argument.
     * @tparam ARG_T The type of the bounds. Must be a floating point type.
     */
    template<IsFloatInvocable FN, FloatingPoint ARG_T>
    RegulaFalsi(FN func, std::initializer_list< ARG_T > bounds) -> RegulaFalsi< decltype(func), ARG_T >;

    /**
     * @brief Deduction guide for the RegulaFalsi class.
     * @tparam FN The type of the function object for which to find the root. The function must be invocable
     * with a double argument.
     * @tparam CONT_T The type of the container holding the bounds. Must be a container of floating point types.
     */
    template<IsFloatInvocable FN, IsContainer CONT_T>
    RegulaFalsi(FN func, CONT_T bounds) -> RegulaFalsi< decltype(func), typename CONT_T::value_type >;

    /**
     * @brief Deduction guide for the RegulaFalsi class.
     * @tparam FN The type of the function object for which to find the root. The function must be invocable
     * with a double argument.
     * @tparam BOUNDS_T The type of the struct holding the bounds. Must be a struct with structured bindings support
     * and with floating point types as its members.
     */
    template<IsFloatInvocable FN, IsFloatStruct BOUNDS_T>
    RegulaFalsi(FN func, BOUNDS_T bounds) -> RegulaFalsi< decltype(func), StructCommonType_t< BOUNDS_T > >;


    // =================================================================================================================
    //
    //     ad88                          88
    //    d8"                            88
    //    88                             88
    //  MM88MMM  ,adPPYba,   ,adPPYba,   88  8b       d8   ,adPPYba,
    //    88     I8[    ""  a8"     "8a  88  `8b     d8'  a8P_____88
    //    88      `"Y8ba,   8b       d8  88   `8b   d8'   8PP"""""""
    //    88     aa    ]8I  "8a,   ,a8"  88    `8b,d8'    "8b,   ,aa
    //    88     `"YbbdP"'   `"YbbdP"'   88      "8"       `"Ybbd8"'
    //
    // =================================================================================================================

    namespace impl
    {
        /**
         * @brief Implementation function for the fsolve functions.
         *
         * This function template takes a solver object, initial bounds, an optional convergence
         * tolerance (epsilon), and an optional maximum number of iterations. It attempts to find the
         * root of the function within the given bounds using the solver's algorithm.
         *
         * @tparam SOLVER The solver type, which must implement the required interface (e.g., evaluate(), init(), iterate()).
         * @param solver The solver object configured with the function for which to find the root.
         * @param bounds A std::pair containing the initial lower and upper bounds for the search interval.
         * @param eps The convergence tolerance (optional, default is 1.0E-6).
         * @param maxiter The maximum number of iterations allowed (optional, default is 100).
         * @return A tl::expected<double, RootError> object, which contains the root on success, or a RootError on failure.
         * @note The solver must implement a compatible interface with the required member functions,
         *       such as evaluate(), init(), and iterate().
         */
        template<typename SOLVER, typename EPS_T, typename ITER_T>
            requires nxx::FloatingPoint< typename SOLVER::RESULT_T > &&
                     std::convertible_to< EPS_T, typename SOLVER::RESULT_T > &&
                     requires(SOLVER solver, std::pair< typename SOLVER::RESULT_T, typename SOLVER::RESULT_T > bounds)
                     {
                         { solver.evaluate(std::declval< double >()) } -> nxx::FloatingPoint;
                         { solver.init(bounds) };
                         { solver.iterate() };
                     }
        auto fsolve_impl(SOLVER                                                            solver,
                         std::pair< typename SOLVER::RESULT_T, typename SOLVER::RESULT_T > bounds,
                         EPS_T                                                             eps,
                         ITER_T                                                            maxiter)
        {
            using ET = RootErrorImpl< typename SOLVER::RESULT_T >;
            using RT = tl::expected< typename SOLVER::RESULT_T, ET >;
            using std::isfinite;

            solver.init(bounds);

            // Declare variables for use in the iteration loop.
            auto curBounds = solver.bounds();
            RT result = (curBounds.first + curBounds.second) / 2.0;
            std::array< std::pair< typename SOLVER::RESULT_T, typename SOLVER::RESULT_T >, 2 > roots{};
            decltype(roots.begin()) min;

            // Check for NaN or Inf.
            if (!isfinite(solver.evaluate(curBounds.first)) || !isfinite(solver.evaluate(curBounds.second))) {
                result = tl::make_unexpected(ET("Invalid initial brackets!", RootErrorType::NumericalError, result.value()));
                return result;
            }

            // Check for a root in the initial bracket.
            if (solver.evaluate(curBounds.first) * solver.evaluate(curBounds.second) > 0.0) {
                result = tl::make_unexpected(ET("Root not bracketed!", RootErrorType::NoRootInBracket, result.value()));
                return result;
            }

            // Begin the iteration loop.
            int iter = 1;
            while (true) {
                curBounds = solver.bounds();
                roots     = { std::make_pair(curBounds.first, abs(solver.evaluate(curBounds.first))),
                              std::make_pair(curBounds.second, abs(solver.evaluate(curBounds.second))) };

                // Check for NaN or Inf.
                if (std::any_of(roots.begin(), roots.end(), [](const auto& r) { return !isfinite(r.second); })) {
                    result = tl::make_unexpected(ET("Non-finite result!", RootErrorType::NumericalError, result.value(), iter));
                    break;
                }

                // Check for convergence.
                min = std::min_element(roots.begin(), roots.end(), [](const auto& a, const auto& b) { return a.second < b.second; });
                if (min->second < eps) {
                    result = min->first;
                    break;
                }

                // Check for maximum number of iterations.
                if (iter >= maxiter) {
                    result = tl::make_unexpected(ET("Max. iterations exceeded!", RootErrorType::MaxIterationsExceeded, min->first, iter));
                    break;
                }

                // Perform one iteration.
                solver.iterate();
                ++iter;
            }

            return result;
        }
    } // namespace impl

    /**
     * @brief Overload of fsolve function that accepts a struct-like object for bounds, e.g., a pair or tuple.
     *
     * @tparam SOLVER_T The solver type used to find the root of the function.
     * @tparam FN_T The type of the function object for which to find the root.
     * @tparam STRUCT_T The type of the bounds object.
     * @tparam EPS_T The type of the tolerance for stopping the algorithm.
     * @tparam ITER_T The type of the maximum number of iterations allowed.
     * @param function The function object for which to find the root.
     * @param bounds Any object that supports structured bindings and provides two values for the lower and upper bounds.
     * @param eps The tolerance for stopping the algorithm.
     * @param maxiter The maximum number of iterations allowed.
     * @return tl::expected object containing either the root of the function or an error.
     */
    template<template< typename, typename > class SOLVER_T,
        IsFloatInvocable FN_T,
        IsFloatStruct STRUCT_T,
        FloatingPoint EPS_T = StructCommonType_t< STRUCT_T >,
        std::integral ITER_T = int>
    auto fsolve(FN_T     function,
                STRUCT_T bounds,
                EPS_T    eps     = epsilon< StructCommonType_t< STRUCT_T > >(),
                ITER_T   maxiter = iterations< StructCommonType_t< STRUCT_T > >())
    {
        auto [lo, hi] = bounds;

        using ARG_T = std::common_type_t< decltype(lo), decltype(hi) >;
        auto solver = SOLVER_T< FN_T, ARG_T >(function);

        return impl::fsolve_impl(solver, std::pair< ARG_T, ARG_T >{ lo, hi }, eps, maxiter);
    }

    /**
     * @brief Overload of fsolve function that accepts an initializer list for bounds.
     *
     * @tparam SOLVER_T The solver type used to find the root of the function.
     * @tparam FN_T The type of the function object for which to find the root.
     * @tparam ARG_T The type of the argument to the function object.
     * @tparam EPS_T The type of the tolerance for stopping the algorithm.
     * @tparam ITER_T The type of the maximum number of iterations allowed.
     * @param function The function object for which to find the root.
     * @param bounds An initializer list containing exactly two elements representing the lower and upper bounds.
     * @param eps The tolerance for stopping the algorithm.
     * @param maxiter The maximum number of iterations allowed.
     * @return tl::expected object containing either the root of the function or an error.
     * @throws NumerixxError if the initializer list does not contain exactly two elements.
     */
    template<template< typename, typename > class SOLVER_T,
        IsFloatInvocable FN_T,
        FloatingPoint ARG_T,
        FloatingPoint EPS_T = ARG_T,
        std::integral ITER_T = int>
    auto fsolve(FN_T                           function,
                std::initializer_list< ARG_T > bounds,
                EPS_T                          eps     = epsilon< ARG_T >(),
                ITER_T                         maxiter = iterations< ARG_T >())
    {
        if (bounds.size() != 2) throw NumerixxError("Initializer list must contain exactly two elements!");
        auto solver = SOLVER_T< FN_T, ARG_T >(function);
        auto bnds   = std::span(bounds.begin(), bounds.end());
        return impl::fsolve_impl(solver, { bnds.front(), bnds.back() }, eps, maxiter);
    }

    /**
     * @brief Overload of fsolve function that accepts a container for bounds.
     *
     * @tparam SOLVER_T The solver type used to find the root of the function.
     * @tparam FN_T The type of the function object for which to find the root.
     * @tparam CONT_T The type of the container for the bounds.
     * @tparam EPS_T The type of the tolerance for stopping the algorithm.
     * @tparam ITER_T The type of the maximum number of iterations allowed.
     * @param function The function object for which to find the root.
     * @param bounds A container containing exactly two elements representing the lower and upper bounds.
     * @param eps The tolerance for stopping the algorithm.
     * @param maxiter The maximum number of iterations allowed.
     * @return tl::expected object containing either the root of the function or an error.
     * @throws NumerixxError if the container does not contain exactly two elements.
     */
    template<template< typename, typename > class SOLVER_T,
        IsFloatInvocable FN_T,
        IsContainer CONT_T,
        FloatingPoint EPS_T = typename CONT_T::value_type,
        std::integral ITER_T = int>
        requires nxx::FloatingPoint< typename CONT_T::value_type >
    auto fsolve(FN_T          function,
                const CONT_T& bounds,
                EPS_T         eps     = epsilon< typename CONT_T::value_type >(),
                ITER_T        maxiter = iterations< typename CONT_T::value_type >())
    {
        if (bounds.size() != 2) throw NumerixxError("Container must contain exactly two elements!");

        using ARG_T = typename CONT_T::value_type;

        auto solver = SOLVER_T< FN_T, ARG_T >(function);
        return impl::fsolve_impl(solver, { bounds.front(), bounds.back() }, eps, maxiter);
    }
} // namespace nxx::roots

#endif    // NUMERIXX_ROOTBRACKETING_HPP

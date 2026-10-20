/*
    888b      88  88        88  88b           d88  88888888888  88888888ba   88
   8b        d8  8b        d8 8888b     88  88        88  888b         d888  88
   88      "8b  88   Y8,    ,8P    Y8,    ,8P 88 `8b    88  88        88  88`8b
   d8'88  88           88      ,8P  88    `8b  d8'      `8b  d8' 88  `8b   88 88
   88  88 `8b     d8' 88  88aaaaa      88aaaaaa8P'  88      Y88P          Y88P
    88   `8b  88  88        88  88  `8b   d8'  88  88"""""      88""""88'    88
   d88b          d88b 88    `8b 88  88        88  88   `8b d8'   88  88 88 `8b
   88    ,8P  Y8,      ,8P  Y8, 88     `8888  Y8a.    .a8P  88    `888'    88 88
   88     `8b   88   d8'    `8b    d8'    `8b 88      `888   `"Y8888Y"'   88 `8'
   88  88888888888  88      `8b  88  8P        Y8  8P        Y8

    Copyright © 2022 Kenneth Troldal Balslev

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the “Software”), to
   deal in the Software without restriction, including without limitation the
   rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
   sell copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#pragma once

// ===== Numerixx Includes
#include "Common.impl.hpp"
#include <Constants.hpp>
#include <Deriv.hpp>

// ===== Standard Library Includes
#include <Functions.hpp>
#include <algorithm>
#include <array>

namespace nxx::roots {
    // =================================================================================================================
    //
    //  88888888ba                                       88 88 88      "8b 88 ,d
    //  "" 88      ,8P                                      88 88 88aaaaaa8P'
    //  8b,dPPYba,  ,adPPYYba,   ,adPPYba,  88   ,d8   ,adPPYba,  MM88MMM  88
    //  8b,dPPYba,    ,adPPYb,d8 88""""""8b,  88P'   "Y8  ""     `Y8  a8"     ""
    //  88 ,a8"   a8P_____88    88     88  88P'   `"8a  a8"    `Y88 88      `8b
    //  88          ,adPPPPP88  8b          8888[     8PP"""""""    88     88 88
    //  88  8b       88 88      a8P  88          88,    ,88  "8a,   ,aa 88`"Yba,
    //  "8b,   ,aa    88,    88  88       88  "8a,   ,d88 88888888P"   88
    //  `"8bbdP"Y8   `"Ybbd8"'  88   `Y8a  `"Ybbd8"'    "Y888  88  88       88
    //  `"YbbdP"Y8
    //                                                                                                    aa,
    //                                                                                                    ,88
    //                                                                                                     "Y8bbdP"
    //
    // =================================================================================================================

    namespace detail {

        template<typename>
        struct BracketingTraits; // Forward declaration remains the same

        // Generic specialization of BracketingTraits
        template<template<typename, typename> class Solver, // Template template parameter
            typename FN,
            typename T>
        struct BracketingTraits<Solver<FN, T>> // Partial specialization
        {
            using FUNCTION_T = FN;
            using ARG_T = T;
            using RETURN_T = std::invoke_result_t<FN, ARG_T>;
        };

        /**
         * @brief Provides a base class template for root bracketing algorithms.
         *
         * The BracketingBase class template serves as a foundational component
         * for algorithms that bracket roots of a given function. It
         * encapsulates common functionalities such as storing the objective
         * function and maintaining the current bounds around the root. This
         * class enforces certain type constraints on the template parameters to
         * ensure compatibility with root bracketing algorithms.
         *
         * @tparam DERIVED The subclass inheriting from BracketingBase.
         * @tparam FUNCTION_T The type of the function for which the root is
         * being bracketed.
         * @tparam ARG_T The type of the argument to the function.
         */
        template<typename DERIVED, typename FUNCTION_T, typename ARG_T>
        requires std::same_as<typename BracketingTraits<DERIVED>::FUNCTION_T, FUNCTION_T>
                 && nxx::IsFloatInvocable<FUNCTION_T> && nxx::IsFloat<ARG_T>
                 && nxx::IsFloat<typename BracketingTraits<DERIVED>::RETURN_T>
        class BracketingBase
        {
            friend DERIVED;

          public:
            static constexpr bool IsBracketingSolver = true; /**< Flag indicating the class is a bracketing solver. */

            using RESULT_T = std::invoke_result_t<FUNCTION_T, ARG_T>; /**< Result type of the function. */
            using BOUNDS_T = std::pair<ARG_T, ARG_T>; /**< Type for representing the bounds around the root. */
            using RETURN_T = std::tuple<ARG_T, ARG_T, ARG_T>;

          protected:
            ~BracketingBase() = default; /**< Protected destructor to prevent direct instantiation. */

          private:
            FUNCTION_T m_func{}; /**< The function object to find the root for. */
            BOUNDS_T m_bounds{}; /**< Holds the current bounds around the root. */
            RETURN_T m_result{};

          public:
            /**
             * @brief Constructs the BracketingBase with a function and bounds
             * from a float struct.
             * @tparam STRUCT_T The struct type holding the bounds.
             * @param objective The function for which the root is being
             * bracketed.
             * @param bounds Struct with the initial bounds.
             */
            BracketingBase(FUNCTION_T objective, IsFloatStruct auto bounds)
              : m_func{ std::move(objective) }, m_bounds{ toPair(bounds) },
                m_result{ m_bounds.first, 0.0, m_bounds.second }
            {
                validateBounds(m_bounds);
            }

            /**
             * @brief Constructs the BracketingBase with a function and bounds
             * from an array.
             * @tparam N The size of the array. Must be 2.
             * @param objective The function for which the root is being
             * bracketed.
             * @param bounds Array with the initial bounds.
             */
            template<size_t N>
            requires(N == 2)
            BracketingBase(FUNCTION_T objective, const ARG_T (&bounds)[N])
              : m_func{ std::move(objective) }, m_bounds{ std::pair{ bounds[0], bounds[1] } },
                m_result{ m_bounds.first, 0.0, m_bounds.second }
            {
                validateBounds(m_bounds);
            }

            /**
             * @brief Sets the bounds for the bracketing solver.
             * @param bounds The new bounds to be set, represented as a pair of
             * values.
             * @throws NumerixxError If the solver has not been initialized.
             * @note This method assumes that the bounds are provided in the
             * correct order (lower, upper).
             */
            void setBounds(const BOUNDS_T &bounds)
            {
                // if (!m_isInitialized) throw NumerixxError("Solver has not
                // been initialized!"); auto [lower, upper] = bounds;
                // static_assert(nxx::IsFloat< decltype(lower) >);
                // m_bounds = BOUNDS_T { lower, upper };
                m_bounds = toPair(bounds);
                std::get<0>(m_result) = m_bounds.first;
                std::get<2>(m_result) = m_bounds.second;
                validateBounds(m_bounds);
            }

            void setBounds(const RETURN_T &range)
            {
                // if (!m_isInitialized) throw NumerixxError("Solver has not
                // been initialized!"); auto [lower, upper] = bounds;
                // static_assert(nxx::IsFloat< decltype(lower) >);
                // m_bounds = BOUNDS_T { lower, upper };
                m_result = range;
                m_bounds = { std::get<0>(m_result), std::get<2>(m_result) };
                validateBounds(m_bounds);
            }

            BracketingBase(const BracketingBase &other) = default; /**< Default copy constructor. */
            BracketingBase(BracketingBase &&other) noexcept = default; /**< Default move constructor. */

            BracketingBase &operator=(const BracketingBase &other) = default; /**< Default copy assignment operator. */
            BracketingBase &operator=(
                BracketingBase &&other) noexcept = default; /**< Default move assignment operator. */

            /**
             * @brief Evaluates the function at a given value.
             * @param value The value at which the function is to be evaluated.
             * @return The result of evaluating the function at the specified
             * value.
             */
            RESULT_T evaluate(IsFloat auto value) { return static_cast<RESULT_T>(m_func(value)); }
            /**
             * @brief Returns the current bounds of the solver.
             * @details This method returns the current bounds being used by the
             * solver. It throws an exception if the solver has not been
             * initialized.
             * @throws NumerixxError If the solver has not been initialized.
             * @return The current bounds of the solver.
             */
            const RETURN_T &current() const { return m_result; }

            /**
             * @brief Performs a single iteration of the algorithm.
             * @details This method must be implemented in the derived class,
             * and performs a single iteration of the root bracketing algorithm.
             */
            void iterate() { std::invoke(static_cast<DERIVED &>(*this)); }
        };
    } // namespace detail

    // =================================================================================================================
    //
    //  88888888ba   88           88           88
    //  88      "8b  ""           88           88
    //  88      ,8P               88           88
    //  88aaaaaa8P'  88   ,adPPYb,88   ,adPPYb,88   ,adPPYba,  8b,dPPYba,
    //  ,adPPYba, 88""""88'    88  a8"    `Y88  a8"    `Y88  a8P_____88  88P'
    //  "Y8  I8[    "" 88    `8b    88  8b       88  8b       88  8PP"""""""  88
    //  `"Y8ba, 88     `8b   88  "8a,   ,d88  "8a,   ,d88  "8b,   ,aa  88 aa ]8I
    //  88      `8b  88   `"8bbdP"Y8   `"8bbdP"Y8   `"Ybbd8"'  88 `"YbbdP"'
    //
    // =================================================================================================================

    /**
     * @brief Defines the Ridder class for performing Ridder's method of root
     * bracketing.
     *
     * Ridder's method is a root-finding algorithm that provides a more robust
     * and often faster convergence than simple bisection. This class template,
     * `Ridder`, inherits from a base class that provides common functionalities
     * for root bracketing algorithms, and adds the specific iteration logic for
     * Ridder's method. It is templated to accept a function and an optional
     * argument type.
     *
     * @tparam FN The type of the function for which the root is being
     * bracketed.
     * @tparam ARG_T The type of the argument to the function, defaults to
     * double.
     */
    template<IsFloatInvocable FN, IsFloat ARG_T = double>
    class Ridder final : public detail::BracketingBase<Ridder<FN, ARG_T>, FN, ARG_T>
    {
        using BASE = detail::BracketingBase<Ridder<FN, ARG_T>, FN, ARG_T>; /**< Base class alias for
                                                                              readability. */

      public:
        using BASE::BASE; /**< Inherits constructors from BracketingBase. */

        /**
         * @brief Performs a single iteration of Ridder's method.
         * @details This method updates the bounds using Ridder's algorithm. It
         * calculates a new estimate for the root and adjusts the bounds
         * accordingly.
         */
        void operator()()
        {
            const auto &bounds = BASE::current();
            using RT = std::invoke_result_t<FN, decltype(std::get<0>(bounds))>;

            const RT &x_lo = std::get<0>(bounds);
            const RT &x_hi = std::get<2>(bounds);
            RT f_lo = BASE::evaluate(x_lo);
            RT f_hi = BASE::evaluate(x_hi);

            RT x_mid = (x_lo + x_hi) / 2.0;
            RT f_mid = BASE::evaluate(x_mid);

            int sign = ((f_lo - f_hi) < 0.0 ? -1 : 1);
            RT x_new = x_mid + (x_mid - x_lo) * ((sign * f_mid) / sqrt(f_mid * f_mid - f_lo * f_hi));
            RT f_new = BASE::evaluate(x_new);

            if (f_mid * f_new < 0.0)
                BASE::setBounds(
                    x_mid < x_new ? std::make_tuple(x_mid, x_new, x_new) : std::make_tuple(x_new, x_new, x_mid));
            else if (f_hi * f_new < 0.0)
                BASE::setBounds(
                    x_hi < x_new ? std::make_tuple(x_hi, x_new, x_new) : std::make_tuple(x_new, x_new, x_hi));
            else
                BASE::setBounds(
                    x_lo < x_new ? std::make_tuple(x_lo, x_new, x_new) : std::make_tuple(x_new, x_new, x_lo));
        }
    };

    /**
     * @brief Deduction guides for Ridder class.
     * Allows the type of Ridder class to be deduced from the constructor
     * parameters.
     */
    template<typename FN, typename ARG_T>
    requires IsFloatInvocable<FN> && IsFloat<ARG_T>
    Ridder(FN, std::initializer_list<ARG_T>) -> Ridder<FN, ARG_T>;

    template<typename FN, typename BOUNDS_T>
    requires IsFloatInvocable<FN> && IsFloatStruct<BOUNDS_T>
    Ridder(FN, BOUNDS_T) -> Ridder<FN, StructCommonType_t<BOUNDS_T>>;

    // =================================================================================================================
    //
    //  88888888ba   88                                              88
    //  88      "8b  ""                                       ,d     ""
    //  88      ,8P                                           88
    //  88aaaaaa8P'  88  ,adPPYba,   ,adPPYba,   ,adPPYba,  MM88MMM  88
    //  ,adPPYba,   8b,dPPYba, 88""""""8b,  88  I8[    ""  a8P_____88  a8" "" 88
    //  88  a8"     "8a  88P'   `"8a 88      `8b  88   `"Y8ba,   8PP"""""""  8b
    //  88     88  8b       d8  88       88 88      a8P  88  aa    ]8I  "8b, ,aa
    //  "8a,   ,aa    88,    88  "8a,   ,a8"  88       88 88888888P"   88
    //  `"YbbdP"'   `"Ybbd8"'   `"Ybbd8"'    "Y888  88   `"YbbdP"'   88       88
    //
    // =================================================================================================================

    /**
     * @brief Defines the Bisection class for performing the bisection method of
     * root bracketing.
     *
     * The Bisection class template is an implementation of the classic
     * bisection method for root finding. This method is a bracketing algorithm
     * that repeatedly bisects an interval and then selects a subinterval in
     * which a root must lie for further processing. It inherits from a base
     * class that provides common functionalities for root bracketing
     * algorithms, and adds the specific iteration logic for the bisection
     * method. The class is templated to accept a function and an optional
     * argument type.
     *
     * @tparam FN The type of the function for which the root is being
     * bracketed.
     * @tparam ARG_T The type of the argument to the function, defaults to
     * double.
     */
    template<IsFloatInvocable FN, IsFloat ARG_T = double>
    class Bisection final : public detail::BracketingBase<Bisection<FN, ARG_T>, FN, ARG_T>
    {
        using BASE = detail::BracketingBase<Bisection<FN, ARG_T>, FN, ARG_T>; /**< Base class alias for
                                                                                 readability. */

      public:
        using BASE::BASE; /**< Inherits constructors from BracketingBase. */

        /**
         * @brief Performs a single iteration of the bisection method.
         * @details This method updates the bounds by bisecting the current
         * interval and choosing the subinterval where the sign of the function
         * changes.
         */
        void operator()()
        {
            const auto &bounds = BASE::current();
            using RT = std::invoke_result_t<FN, decltype(std::get<0>(bounds))>;

            if (RT root = (std::get<0>(bounds) + std::get<2>(bounds)) / 2.0;
                BASE::evaluate(std::get<0>(bounds)) * BASE::evaluate(root) < 0.0)
                BASE::setBounds({ std::make_tuple(std::get<0>(bounds), (std::get<0>(bounds) + root) / 2, root) });
            else
                BASE::setBounds({ std::make_tuple(root, (root + std::get<2>(bounds)) / 2, std::get<2>(bounds)) });
        }
    };

    /**
     * @brief Deduction guides for Bisection class.
     * Allows the type of Bisection class to be deduced from the constructor
     * parameters.
     */
    template<typename FN, typename ARG_T>
    requires IsFloatInvocable<FN> && IsFloat<ARG_T>
    Bisection(FN, std::initializer_list<ARG_T>) -> Bisection<FN, ARG_T>;

    template<typename FN, typename BOUNDS_T>
    requires IsFloatInvocable<FN> && IsFloatStruct<BOUNDS_T>
    Bisection(FN, BOUNDS_T) -> Bisection<FN, StructCommonType_t<BOUNDS_T>>;

    // =================================================================================================================
    //
    //  88888888ba                                        88 88888888888 88 88
    //  88      "8b                                       88              88 88
    //  "" 88      ,8P                                       88              88
    //  88 88aaaaaa8P'  ,adPPYba,   ,adPPYb,d8  88       88  88  ,adPPYYba,
    //  88aaaaa  ,adPPYYba,  88  ,adPPYba,  88 88""""88'   a8P_____88  a8" `Y88
    //  88       88  88  ""     `Y8  88"""""  ""     `Y8  88  I8[    ""  88 88
    //  `8b   8PP"""""""  8b       88  88       88  88  ,adPPPPP88  88
    //  ,adPPPPP88  88   `"Y8ba,   88 88     `8b  "8b,   ,aa  "8a,   ,d88  "8a,
    //  ,a88  88  88,    ,88  88       88,    ,88  88  aa    ]8I  88 88      `8b
    //  `"Ybbd8"'   `"YbbdP"Y8   `"YbbdP'Y8  88  `"8bbdP"Y8  88       `"8bbdP"Y8
    //  88  `"YbbdP"'  88
    //                           aa,    ,88
    //                            "Y8bbdP"
    //
    // =================================================================================================================

    /**
     * @brief Defines the RegulaFalsi class for performing the regula falsi
     * (false position) method of root bracketing.
     *
     * The RegulaFalsi class template implements the regula falsi method, also
     * known as the false position method, for root finding. This method is a
     * bracketing algorithm similar to the bisection method but, instead of
     * bisecting the interval, it uses a linear approximation to guess the root.
     * This class inherits from a base class that provides common
     * functionalities for root bracketing algorithms and adds the specific
     * iteration logic for the regula falsi method. It is templated to accept a
     * function and an optional argument type.
     *
     * @tparam FN The type of the function for which the root is being
     * bracketed.
     * @tparam ARG_T The type of the argument to the function, defaults to
     * double.
     */
    template<IsFloatInvocable FN, IsFloat ARG_T = double>
    class RegulaFalsi final : public detail::BracketingBase<RegulaFalsi<FN, ARG_T>, FN, ARG_T>
    {
        using BASE = detail::BracketingBase<RegulaFalsi<FN, ARG_T>, FN, ARG_T>; /**< Base class alias for
                                                                                   readability. */

      public:
        using BASE::BASE; /**< Inherits constructors from BracketingBase. */

        /**
         * @brief Performs a single iteration of the regula falsi method.
         * @details This method updates the bounds by applying the regula falsi
         * formula to find a new estimate for the root, then adjusts the bounds
         * accordingly.
         */
        void operator()()
        {
            const auto &bounds = BASE::current();
            using RT = std::invoke_result_t<FN, decltype(std::get<0>(bounds))>;

            RT f_lo = BASE::evaluate(std::get<0>(bounds));
            RT f_hi = BASE::evaluate(std::get<2>(bounds));

            RT root = (std::get<0>(bounds) * f_hi - std::get<2>(bounds) * f_lo) / (f_hi - f_lo);
            RT f_root = BASE::evaluate(root);

            if (f_lo * f_root < 0.0)
                BASE::setBounds({ std::make_tuple(std::get<0>(bounds), root, root) });
            else
                BASE::setBounds({ std::make_tuple(root, root, std::get<2>(bounds)) });
        }
    };

    /**
     * @brief Deduction guides for RegulaFalsi class.
     * Allows the type of RegulaFalsi class to be deduced from the constructor
     * parameters.
     */
    template<typename FN, typename ARG_T>
    requires IsFloatInvocable<FN> && IsFloat<ARG_T>
    RegulaFalsi(FN, std::initializer_list<ARG_T>) -> RegulaFalsi<FN, ARG_T>;

    template<typename FN, typename BOUNDS_T>
    requires IsFloatInvocable<FN> && IsFloatStruct<BOUNDS_T>
    RegulaFalsi(FN, BOUNDS_T) -> RegulaFalsi<FN, StructCommonType_t<BOUNDS_T>>;

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

    template<std::integral ITER_T, IsFloat RESULT_T>
    using BracketIterData = IterData<ITER_T, RESULT_T, RESULT_T, RESULT_T>;

    struct BracketBehavior
    {
        template<typename ITER_T, typename EPS_T>
        bool operator()(BracketIterData<ITER_T, EPS_T> iterData, std::integral auto maxiter, IsFloat auto eps) const
        {
            const auto &[iter, lower, x, upper] = iterData;

            if ((upper - lower) <= eps * x + eps / 2) return true;
            if (iter >= maxiter) return true;

            return false;
        }
    };

    template<typename... Args>
    using BracketStopToken = StopToken<BracketBehavior, Args...>;

    namespace detail {

        /**
         * @brief A function template that implements the solver loop for root-finding problems.
         *
         * This function template is used by the fsolve_common function to perform the iterative process of a
         * root-finding problem. It accepts a solver and a termination condition, and performs iterations until the
         * termination condition is met.
         *
         * @tparam SOLVER The type of the solver to use for the root-finding problem. Must be a class that has
         * IsBracketingSolver as true.
         * @tparam TOKEN_T The type of the callable object that specifies the termination condition for the solver.
         * @param solver The solver to use for the root-finding problem.
         * @param terminator The termination condition for the solver.
         * @return The result of the root-finding problem, encapsulated in a BracketSolverResult object.
         * @note This function requires that the termination condition callable is invocable with an IterData object.
         */
        template<typename SOLVER, typename TOKEN_T>
        requires SOLVER::IsBracketingSolver
        auto fsolve_impl(const SOLVER &solver, const TOKEN_T &terminator)
        {
            SOLVER _solver = solver;
            TOKEN_T _terminator = terminator;
            using ARG_T = typename SOLVER::RESULT_T;

            const auto &[lower, x, upper] = _solver.current();
            size_t iter = 0;

            BracketIterData<size_t, ARG_T> iterData{ iter, lower, x, upper };

            auto &[_iter, _lower, _guess, _upper] = iterData;

            while (true) {
                _iter = iter;
                _lower = lower;
                _guess = x;
                _upper = upper;

                if (_terminator(iterData)) break;
                _solver.iterate();
                ++iter;
            }

            return ResultProxy<BracketIterData<size_t, ARG_T>, 0, 2>(iterData);
        }

    } // namespace detail

    /**
     * @brief A function template that solves a root-finding problem using a specified solver and optionally the
     *required epsilon value and/or maximum iteration count.
     *
     * This function template is an overload of the fsolve function that accepts a structure for the bounds, such as a
     * std::pair, a std::tuple, or a custom structure with two floating point elements. It uses a specified solver to
     * find the root of a given function within the provided bounds. The required epsilon value and/or maximum iteration
     * count can be optionally passed as additional arguments. If no additional arguments are passed, the default values
     * for the epsilon value and maximum iteration count are used.
     *
     * @tparam SOLVER_T The type of the solver to use for the root-finding problem. Must be one of the pre-defined
     * solver classes, such as Bisection, Ridder, or RegulaFalsi, or a custom solver class that meets the requirements.
     * @tparam FN_T The type of the function for which the root is being found. Must be invocable with a floating point
     * or complex number. This can be a function pointer, a function object, or a lambda function.
     * @tparam STRUCT_T The type of the structure for the bounds. This can be a std::pair, std::tuple, or a custom
     * structure. The structure must have exactly two floating point elements.
     * @tparam ARGS The type of the additional arguments. These arguments are passed to the fsolve_common function.
     *
     * @param func The function for which the root is being found.
     * @param bounds The bounds within which the root is being found. Must be a structure, such as a std::pair or a
     * std::tuple, that contains exactly two floating point elements.
     * @param args The additional arguments passed to the fsolve_common function.
     * @return The result of the root-finding problem.
     *
     * @pre The function must be continuous and have opposite signs at the bounds to ensure that a root exists within
     * the bounds.
     * @post The result of the root-finding problem is returned as a BracketSolverResult object, which can be used to
     * retrieve the result in a specified format.
     */
    template<template<typename, typename> class SOLVER_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloatStruct STRUCT_T,
        typename... ARGS>
    auto fsolve(FN_T func, STRUCT_T bounds, ARGS... args)
    {
        using SOLVER = SOLVER_T<FN_T, StructCommonType_t<STRUCT_T>>;
        return detail::fsolve_impl(SOLVER(func, bounds), makeToken<BracketStopToken>(args...));
    }

    /**
     * @brief A function template that solves a root-finding problem using a specified solver and termination condition.
     *
     * This function template is an overload of the fsolve function that accepts a structure for the bounds, such as a
     * std::pair, a std::tuple, or a custom structure with two floating point elements. It uses a specified solver to
     * find the root of a given function within the provided bounds. The termination condition for the solver is
     * specified by a callable object, that is called for each iteration of the solver. The termination condition
     * callable must accept an IterData object, which contains the current iteration count and the current bounds around
     * the root.
     *
     * @tparam SOLVER_T The type of the solver to use for the root-finding problem. Must be a template class that
     * accepts two type parameters.
     * @tparam TOKEN_T The type of the callable object that specifies the termination condition for the solver.
     * @tparam FN_T The type of the function for which the root is being found. Must be invocable with a floating point
     * or complex number.
     * @tparam STRUCT_T The type of the structure for the bounds.
     * @param func The function for which the root is being found.
     * @param bounds The bounds within which the root is being found.
     * @return A proxy object containing the result of the root-finding problem.
     * @note This function requires that the termination condition callable is invocable with an IterData object.
     */
    template<template<typename, typename> class SOLVER_T,
        typename TOKEN_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloatStruct STRUCT_T>
    // TODO: Should be able to accept functors that take any king of argument, not just references.
    requires std::invocable<TOKEN_T, BracketIterData<size_t, StructCommonType_t<STRUCT_T>> &>
    auto fsolve(FN_T func, STRUCT_T bounds)
    {
        using SOLVER = SOLVER_T<FN_T, StructCommonType_t<STRUCT_T>>;
        return detail::fsolve_impl(SOLVER(func, bounds), TOKEN_T{});
    }

    /**
     * @brief A function template that solves a root-finding problem using a specified solver and termination condition.
     *
     * This function template is an overload of the fsolve function that accepts a fixed-size array for the bounds.
     * It uses a specified solver to find the root of a given function within the provided bounds.
     * The termination condition for the solver is specified by a callable object.
     *
     * @tparam SOLVER_T The type of the solver to use for the root-finding problem. Must be a template class that
     * accepts two type parameters.
     * @tparam FN_T The type of the function for which the root is being found. Must be invocable with a floating point
     * or complex number.
     * @tparam ARG_T The type of the argument to the function. Must be a floating point number.
     * @tparam N The size of the array for the bounds. Must be 2.
     * @tparam ARGS The type of the additional arguments. These arguments are passed to the fsolve_common function.
     * @param func The function for which the root is being found.
     * @param bounds The bounds within which the root is being found. Must be a fixed-size array of 2 elements.
     * @param args The additional arguments passed to the fsolve_common function.
     * @return The result of the root-finding problem.
     * @note This function requires that the termination condition callable is invocable with an IterData object.
     */
    template<template<typename, typename> class SOLVER_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloat ARG_T,
        size_t N,
        typename... ARGS>
    requires(N == 2)
    auto fsolve(FN_T func, const ARG_T (&bounds)[N], ARGS &&...args)
    {
        using SOLVER = SOLVER_T<FN_T, ARG_T>;
        return detail::fsolve_impl(
            SOLVER(func, std::pair{ bounds[0], bounds[1] }), makeToken<BracketStopToken>(std::forward<ARGS>(args)...));
    }

    /**
     * @brief A function template that solves a root-finding problem using a specified solver and termination condition.
     *
     * This function template is an overload of the fsolve function that accepts a fixed-size array for the bounds.
     * It uses a specified solver to find the root of a given function within the provided bounds.
     * The termination condition for the solver is specified by a callable object.
     *
     * @tparam SOLVER_T The type of the solver to use for the root-finding problem. Must be a template class that
     * accepts two type parameters.
     * @tparam TOKEN_T The type of the callable object that specifies the termination condition for the solver.
     * @tparam FN_T The type of the function for which the root is being found. Must be invocable with a floating point
     * or complex number.
     * @tparam ARG_T The type of the argument to the function. Must be a floating point number.
     * @tparam N The size of the array for the bounds. Must be 2.
     * @param func The function for which the root is being found.
     * @param bounds The bounds within which the root is being found. Must be a fixed-size array of 2 elements.
     * @return The result of the root-finding problem.
     * @note This function requires that the termination condition callable is invocable with an IterData object.
     */
    template<template<typename, typename> class SOLVER_T,
        typename TOKEN_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloat ARG_T,
        size_t N>
    // TODO: Should be able to accept functors that take any king of argument, not just references.
    requires(N == 2) && std::invocable<TOKEN_T, BracketIterData<size_t, ARG_T> &>
    auto fsolve(FN_T func, const ARG_T (&bounds)[N])
    {
        using SOLVER = SOLVER_T<FN_T, ARG_T>;
        return detail::fsolve_impl(SOLVER(func, std::pair{ bounds[0], bounds[1] }), TOKEN_T{});
    }

} // namespace nxx::roots

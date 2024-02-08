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

#pragma once

// ===== Numerixx Includes
#include "RootCommon.hpp"
#include <Constants.hpp>
#include <Deriv.hpp>

// ===== Standard Library Includes
#include <Functions.hpp>
#include <algorithm>
#include <array>

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

    namespace detail
    {
        /**
         * @brief Provides a base class template for root bracketing algorithms.
         *
         * The BracketingBase class template serves as a foundational component for
         * algorithms that bracket roots of a given function. It encapsulates common
         * functionalities such as storing the objective function and maintaining the
         * current bounds around the root. This class enforces certain type constraints
         * on the template parameters to ensure compatibility with root bracketing algorithms.
         *
         * @tparam DERIVED The subclass inheriting from BracketingBase.
         * @tparam FUNCTION_T The type of the function for which the root is being bracketed.
         * @tparam ARG_T The type of the argument to the function.
         */
        template< typename DERIVED, typename FUNCTION_T, typename ARG_T >
        requires std::same_as< typename BracketingTraits< DERIVED >::FUNCTION_T, FUNCTION_T > && nxx::IsFloatInvocable< FUNCTION_T > &&
                 nxx::IsFloat< ARG_T > && nxx::IsFloat< typename BracketingTraits< DERIVED >::RETURN_T >
        class BracketingBase
        {
            friend DERIVED;

        public:
            static constexpr bool IsBracketingSolver = true; /**< Flag indicating the class is a bracketing solver. */

            using RESULT_T = std::invoke_result_t< FUNCTION_T, ARG_T >; /**< Result type of the function. */
            using BOUNDS_T = std::pair< ARG_T, ARG_T >;                 /**< Type for representing the bounds around the root. */
            using RETURN_T = std::tuple< ARG_T, ARG_T, ARG_T >;

        protected:
            ~BracketingBase() = default; /**< Protected destructor to prevent direct instantiation. */

        private:
            FUNCTION_T m_func {};   /**< The function object to find the root for. */
            BOUNDS_T   m_bounds {}; /**< Holds the current bounds around the root. */
            RETURN_T   m_result {};

        public:
            /**
             * @brief Constructs the BracketingBase with a function and bounds from a float struct.
             * @tparam STRUCT_T The struct type holding the bounds.
             * @param objective The function for which the root is being bracketed.
             * @param bounds Struct with the initial bounds.
             */
            BracketingBase(FUNCTION_T objective, IsFloatStruct auto bounds)
                : m_func { std::move(objective) },
                  m_bounds { toPair(bounds) },
                  m_result { m_bounds.first, 0.0, m_bounds.second }
            {
                validateBounds(m_bounds);
            }

            /**
             * @brief Constructs the BracketingBase with a function and bounds from an array.
             * @tparam N The size of the array. Must be 2.
             * @param objective The function for which the root is being bracketed.
             * @param bounds Array with the initial bounds.
             */
            template< size_t N >
            requires(N == 2)
            BracketingBase(FUNCTION_T objective, const ARG_T (&bounds)[N])
                : m_func { std::move(objective) },
                  m_bounds { std::pair { bounds[0], bounds[1] } },
                  m_result { m_bounds.first, 0.0, m_bounds.second }
            {
                validateBounds(m_bounds);
            }

            /**
             * @brief Sets the bounds for the bracketing solver.
             * @param bounds The new bounds to be set, represented as a pair of values.
             * @throws NumerixxError If the solver has not been initialized.
             * @note This method assumes that the bounds are provided in the correct order (lower, upper).
             */
            void setBounds(const BOUNDS_T& bounds)
            {
                // if (!m_isInitialized) throw NumerixxError("Solver has not been initialized!");
                // auto [lower, upper] = bounds;
                // static_assert(nxx::IsFloat< decltype(lower) >);
                // m_bounds = BOUNDS_T { lower, upper };
                m_bounds                = toPair(bounds);
                std::get< 0 >(m_result) = m_bounds.first;
                std::get< 2 >(m_result) = m_bounds.second;
                validateBounds(m_bounds);
            }

            void setBounds(const RETURN_T& range)
            {
                // if (!m_isInitialized) throw NumerixxError("Solver has not been initialized!");
                // auto [lower, upper] = bounds;
                // static_assert(nxx::IsFloat< decltype(lower) >);
                // m_bounds = BOUNDS_T { lower, upper };
                m_result = range;
                m_bounds = { std::get< 0 >(m_result), std::get< 2 >(m_result) };
                validateBounds(m_bounds);
            }

            BracketingBase(const BracketingBase& other)     = default; /**< Default copy constructor. */
            BracketingBase(BracketingBase&& other) noexcept = default; /**< Default move constructor. */

            BracketingBase& operator=(const BracketingBase& other)     = default; /**< Default copy assignment operator. */
            BracketingBase& operator=(BracketingBase&& other) noexcept = default; /**< Default move assignment operator. */

            /**
             * @brief Evaluates the function at a given value.
             * @param value The value at which the function is to be evaluated.
             * @return The result of evaluating the function at the specified value.
             */
            RESULT_T evaluate(IsFloat auto value) { return static_cast< RESULT_T >(m_func(value)); }
            /**
             * @brief Returns the current bounds of the solver.
             * @details This method returns the current bounds being used by the solver.
             *          It throws an exception if the solver has not been initialized.
             * @throws NumerixxError If the solver has not been initialized.
             * @return The current bounds of the solver.
             */
            const RETURN_T& current() const { return m_result; }

            /**
             * @brief Performs a single iteration of the algorithm.
             * @details This method must be implemented in the derived class, and
             *          performs a single iteration of the root bracketing algorithm.
             */
            void iterate() { std::invoke(static_cast< DERIVED& >(*this)); }
        };
    }    // namespace detail

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
     * @brief Defines the Ridder class for performing Ridder's method of root bracketing.
     *
     * Ridder's method is a root-finding algorithm that provides a more robust and often faster
     * convergence than simple bisection. This class template, `Ridder`, inherits from a base class
     * that provides common functionalities for root bracketing algorithms, and adds the specific
     * iteration logic for Ridder's method. It is templated to accept a function and an optional
     * argument type.
     *
     * @tparam FN The type of the function for which the root is being bracketed.
     * @tparam ARG_T The type of the argument to the function, defaults to double.
     */
    template< IsFloatInvocable FN, IsFloat ARG_T = double >
    class Ridder final : public detail::BracketingBase< Ridder< FN, ARG_T >, FN, ARG_T >
    {
        using BASE = detail::BracketingBase< Ridder< FN, ARG_T >, FN, ARG_T >; /**< Base class alias for readability. */

    public:
        using BASE::BASE; /**< Inherits constructors from BracketingBase. */

        /**
         * @brief Performs a single iteration of Ridder's method.
         * @details This method updates the bounds using Ridder's algorithm. It calculates
         *          a new estimate for the root and adjusts the bounds accordingly.
         */
        void operator()()
        {
            const auto& bounds = BASE::current();
            using RT           = std::invoke_result_t< FN, decltype(std::get< 0 >(bounds)) >;

            const RT& x_lo = std::get< 0 >(bounds);
            const RT& x_hi = std::get< 2 >(bounds);
            RT        f_lo = BASE::evaluate(x_lo);
            RT        f_hi = BASE::evaluate(x_hi);

            RT x_mid = (x_lo + x_hi) / 2.0;
            RT f_mid = BASE::evaluate(x_mid);

            int sign  = ((f_lo - f_hi) < 0.0 ? -1 : 1);
            RT  x_new = x_mid + (x_mid - x_lo) * ((sign * f_mid) / sqrt(f_mid * f_mid - f_lo * f_hi));
            RT  f_new = BASE::evaluate(x_new);

            // Update bounds based on the results of Ridder's method.
            // if (f_mid * f_new < 0.0)
            //     BASE::setBounds(x_mid < x_new ? std::make_pair(x_mid, x_new) : std::make_pair(x_new, x_mid));
            // else if (f_hi * f_new < 0.0)
            //     BASE::setBounds(x_hi < x_new ? std::make_pair(x_hi, x_new) : std::make_pair(x_new, x_hi));
            // else
            //     BASE::setBounds(x_lo < x_new ? std::make_pair(x_lo, x_new) : std::make_pair(x_new, x_lo));
            if (f_mid * f_new < 0.0)
                BASE::setBounds(x_mid < x_new ? std::make_tuple(x_mid, x_new, x_new) : std::make_tuple(x_new, x_new, x_mid));
            else if (f_hi * f_new < 0.0)
                BASE::setBounds(x_hi < x_new ? std::make_tuple(x_hi, x_new, x_new) : std::make_tuple(x_new, x_new, x_hi));
            else
                BASE::setBounds(x_lo < x_new ? std::make_tuple(x_lo, x_new, x_new) : std::make_tuple(x_new, x_new, x_lo));
        }
    };

    /**
     * @brief Deduction guides for Ridder class.
     * Allows the type of Ridder class to be deduced from the constructor parameters.
     */
    template< typename FN, typename ARG_T >
    requires IsFloatInvocable< FN > && IsFloat< ARG_T >
    Ridder(FN, std::initializer_list< ARG_T >) -> Ridder< FN, ARG_T >;

    template< typename FN, typename BOUNDS_T >
    requires IsFloatInvocable< FN > && IsFloatStruct< BOUNDS_T >
    Ridder(FN, BOUNDS_T) -> Ridder< FN, StructCommonType_t< BOUNDS_T > >;

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
     * @brief Defines the Bisection class for performing the bisection method of root bracketing.
     *
     * The Bisection class template is an implementation of the classic bisection method for root finding.
     * This method is a bracketing algorithm that repeatedly bisects an interval and then selects a subinterval
     * in which a root must lie for further processing. It inherits from a base class that provides common
     * functionalities for root bracketing algorithms, and adds the specific iteration logic for the bisection method.
     * The class is templated to accept a function and an optional argument type.
     *
     * @tparam FN The type of the function for which the root is being bracketed.
     * @tparam ARG_T The type of the argument to the function, defaults to double.
     */
    template< IsFloatInvocable FN, IsFloat ARG_T = double >
    class Bisection final : public detail::BracketingBase< Bisection< FN, ARG_T >, FN, ARG_T >
    {
        using BASE = detail::BracketingBase< Bisection< FN, ARG_T >, FN, ARG_T >; /**< Base class alias for readability. */

    public:
        using BASE::BASE; /**< Inherits constructors from BracketingBase. */

        /**
         * @brief Performs a single iteration of the bisection method.
         * @details This method updates the bounds by bisecting the current interval and
         *          choosing the subinterval where the sign of the function changes.
         */
        void operator()()
        {
            const auto& bounds = BASE::current();
            using RT           = std::invoke_result_t< FN, decltype(std::get< 0 >(bounds)) >;

            if (RT root = (std::get< 0 >(bounds) + std::get< 2 >(bounds)) / 2.0;
                BASE::evaluate(std::get< 0 >(bounds)) * BASE::evaluate(root) < 0.0)
                BASE::setBounds({ std::make_tuple(std::get< 0 >(bounds), (std::get< 0 >(bounds) + root) / 2, root) });
            else
                BASE::setBounds({ std::make_tuple(root, (root + std::get< 2 >(bounds)) / 2, std::get< 2 >(bounds)) });
        }
    };

    /**
     * @brief Deduction guides for Bisection class.
     * Allows the type of Bisection class to be deduced from the constructor parameters.
     */
    template< typename FN, typename ARG_T >
    requires IsFloatInvocable< FN > && IsFloat< ARG_T >
    Bisection(FN, std::initializer_list< ARG_T >) -> Bisection< FN, ARG_T >;

    template< typename FN, typename BOUNDS_T >
    requires IsFloatInvocable< FN > && IsFloatStruct< BOUNDS_T >
    Bisection(FN, BOUNDS_T) -> Bisection< FN, StructCommonType_t< BOUNDS_T > >;

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
     * @brief Defines the RegulaFalsi class for performing the regula falsi (false position) method of root bracketing.
     *
     * The RegulaFalsi class template implements the regula falsi method, also known as the false position method,
     * for root finding. This method is a bracketing algorithm similar to the bisection method but, instead of bisecting
     * the interval, it uses a linear approximation to guess the root. This class inherits from a base class that provides
     * common functionalities for root bracketing algorithms and adds the specific iteration logic for the regula falsi method.
     * It is templated to accept a function and an optional argument type.
     *
     * @tparam FN The type of the function for which the root is being bracketed.
     * @tparam ARG_T The type of the argument to the function, defaults to double.
     */
    template< IsFloatInvocable FN, IsFloat ARG_T = double >
    class RegulaFalsi final : public detail::BracketingBase< RegulaFalsi< FN, ARG_T >, FN, ARG_T >
    {
        using BASE = detail::BracketingBase< RegulaFalsi< FN, ARG_T >, FN, ARG_T >; /**< Base class alias for readability. */

    public:
        using BASE::BASE; /**< Inherits constructors from BracketingBase. */

        /**
         * @brief Performs a single iteration of the regula falsi method.
         * @details This method updates the bounds by applying the regula falsi formula
         *          to find a new estimate for the root, then adjusts the bounds accordingly.
         */
        void operator()()
        {
            const auto& bounds = BASE::current();
            using RT           = std::invoke_result_t< FN, decltype(std::get< 0 >(bounds)) >;

            RT f_lo = BASE::evaluate(std::get< 0 >(bounds));
            RT f_hi = BASE::evaluate(std::get< 2 >(bounds));

            RT root   = (std::get< 0 >(bounds) * f_hi - std::get< 2 >(bounds) * f_lo) / (f_hi - f_lo);
            RT f_root = BASE::evaluate(root);

            if (f_lo * f_root < 0.0)
                BASE::setBounds({ std::make_tuple(std::get< 0 >(bounds), root, root) });
            else
                BASE::setBounds({ std::make_tuple(root, root, std::get< 2 >(bounds)) });
        }
    };

    /**
     * @brief Deduction guides for RegulaFalsi class.
     * Allows the type of RegulaFalsi class to be deduced from the constructor parameters.
     */
    template< typename FN, typename ARG_T >
    requires IsFloatInvocable< FN > && IsFloat< ARG_T >
    RegulaFalsi(FN, std::initializer_list< ARG_T >) -> RegulaFalsi< FN, ARG_T >;

    template< typename FN, typename BOUNDS_T >
    requires IsFloatInvocable< FN > && IsFloatStruct< BOUNDS_T >
    RegulaFalsi(FN, BOUNDS_T) -> RegulaFalsi< FN, StructCommonType_t< BOUNDS_T > >;

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

    template< std::integral ITER_T, IsFloat RESULT_T >
    struct IterData
    {
        ITER_T   iter;
        RESULT_T lower;
        RESULT_T guess;
        RESULT_T upper;
    };

    template< IsFloat EPS_T, std::integral ITER_T >
    class BracketTerminator
    {
        EPS_T  m_eps;
        ITER_T m_maxiter;

    public:
        explicit BracketTerminator()
            : m_eps(epsilon< double >()),
              m_maxiter(iterations< double >())
        {}

        explicit BracketTerminator(EPS_T eps, ITER_T maxiter)
            : m_eps(eps),
              m_maxiter(maxiter)
        {}

        explicit BracketTerminator(ITER_T maxiter, EPS_T eps)
            : m_eps(eps),
              m_maxiter(maxiter)
        {}

        explicit BracketTerminator(EPS_T eps)
            : m_eps(eps),
              m_maxiter(iterations< double >())
        {}

        explicit BracketTerminator(ITER_T maxiter)
            : m_eps(epsilon< double >()),
              m_maxiter(maxiter)
        {}

        bool operator()(const auto& data) const
        {
            const auto& [iter, lower, x, upper] = data;

            if ((upper - lower) <= m_eps * x + m_eps / 2) return true;
            if (iter >= m_maxiter) return true;

            return false;
        }
    };

    BracketTerminator() -> BracketTerminator< double, size_t >;
    BracketTerminator(IsFloat auto eps) -> BracketTerminator< decltype(eps), size_t >;
    BracketTerminator(std::integral auto maxiter) -> BracketTerminator< double, decltype(maxiter) >;
    BracketTerminator(IsFloat auto eps, std::integral auto maxiter) -> BracketTerminator< decltype(eps), decltype(maxiter) >;
    BracketTerminator(std::integral auto maxiter, IsFloat auto eps) -> BracketTerminator< decltype(eps), decltype(maxiter) >;

    namespace detail
    {
        // Concept for checking if a type is float or integral
        template< typename T >
        concept FloatOrIntegral = std::is_floating_point_v< T > || std::is_integral_v< T >;

        // Concept to check the validity of arguments
        template< typename... Args >
        concept ValidArgs =
            requires {
                requires(sizeof...(Args) <= 2);
                requires(sizeof...(Args) == 0) || (sizeof...(Args) == 1 && (FloatOrIntegral< Args > && ...)) ||
                            (sizeof...(Args) == 2 && ((FloatOrIntegral< Args > && ...) && (std::is_floating_point_v< Args > != ...)));
            };

        /**
         * @brief Implements a generic root finding solver function template for bracketing solvers.
         *
         * This function template, `fsolve_impl`, provides a generic implementation for root finding
         * using various bracketing solver algorithms. It is designed to work with solvers that conform
         * to the requirements of bracketing solvers, such as having a defined `IsBracketingSolver` static member,
         * initialization, and iteration methods. The function handles initialization, iteration, and
         * convergence checking, returning the result along with any potential errors encountered during
         * the solving process.
         *
         * @tparam SOLVER The type of the solver to be used in root finding. Must conform to the bracketing solver concept.
         */

        // template< typename SOLVER >
        // requires SOLVER::IsBracketingSolver
        // auto fsolve_impl(SOLVER solver, IsFloat auto eps, std::integral auto maxiter)
        // {
        //     using ET = RootErrorImpl< typename SOLVER::RESULT_T >;    /**< Type for error handling. */
        //     using RT = tl::expected< typename SOLVER::RESULT_T, ET >; /**< Type for the function return value. */
        //     using std::isfinite;
        //
        //     // Declare variables for use in the iteration loop.
        //     auto curBounds = solver.current();
        //     RT   result    = (std::get< 0 >(curBounds) + std::get< 2 >(curBounds)) / 2.0;
        //     std::array< std::pair< typename SOLVER::RESULT_T, typename SOLVER::RESULT_T >, 2 > roots {};
        //     decltype(roots.begin())                                                            min;
        //
        //     // Check for NaN or Inf in the initial bounds.
        //     if (!isfinite(solver.evaluate(std::get< 0 >(curBounds))) || !isfinite(solver.evaluate(std::get< 2 >(curBounds)))) {
        //         result = tl::make_unexpected(ET("Invalid initial brackets!", RootErrorType::NumericalError, result.value()));
        //         return result;
        //     }
        //
        //     // Check for a root in the initial bracket.
        //     if (solver.evaluate(std::get< 0 >(curBounds)) * solver.evaluate(std::get< 2 >(curBounds)) > 0.0) {
        //         result = tl::make_unexpected(ET("Root not bracketed!", RootErrorType::NoRootInBracket, result.value()));
        //         return result;
        //     }
        //
        //     // Begin the iteration loop.
        //     int iter = 1;
        //     while (true) {
        //         curBounds = solver.current();
        //         roots     = { std::make_pair(std::get< 0 >(curBounds), abs(solver.evaluate(std::get< 0 >(curBounds)))),
        //                       std::make_pair(std::get< 2 >(curBounds), abs(solver.evaluate(std::get< 2 >(curBounds)))) };
        //
        //         // Check for NaN or Inf.
        //         if (std::any_of(roots.begin(), roots.end(), [](const auto& r) { return !isfinite(r.second); })) {
        //             result = tl::make_unexpected(ET("Non-finite result!", RootErrorType::NumericalError, result.value(), iter));
        //             break;
        //         }
        //
        //         // Check for convergence.
        //         min = std::min_element(roots.begin(), roots.end(), [](const auto& a, const auto& b) { return a.second < b.second; });
        //         if (min->second < eps) {
        //             result = min->first;
        //             break;
        //         }
        //
        //         // Check for maximum number of iterations.
        //         if (iter >= maxiter) {
        //             result = tl::make_unexpected(ET("Max. iterations exceeded!", RootErrorType::MaxIterationsExceeded, min->first,
        //             iter)); break;
        //         }
        //
        //         // Perform one iteration.
        //         solver.iterate();
        //         ++iter;
        //     }
        //
        //     return result;
        // }

        template< std::integral ITER_T, IsFloat RESULT_T >
        class BracketSolverResult
        {
            IterData< ITER_T, RESULT_T > m_iterData;

        public:
            explicit BracketSolverResult(IterData< ITER_T, RESULT_T > iterData)
                : m_iterData(iterData)
            {}

            template< typename OUTPUT_T = RESULT_T >
            OUTPUT_T result()
            {
                return m_iterData.guess;
            }
        };

        template< typename SOLVER, typename TERMINATOR >
        // requires SOLVER::IsBracketOptimizer
        auto fsolve_impl(const SOLVER& solver, const TERMINATOR& terminator)
        {
            SOLVER     _solver     = solver;
            TERMINATOR _terminator = terminator;

            const auto& [lower, x, upper] = _solver.current();
            size_t iter                   = 0;

            IterData< size_t, double > iterData { iter, lower, x, upper };

            while (true) {
                iterData.iter  = iter;
                iterData.lower = lower;
                iterData.guess = x;
                iterData.upper = upper;

                if (_terminator(iterData)) break;
                _solver.iterate();
                ++iter;
            }

            // return std::get< 1 >(_solver.current());
            return BracketSolverResult(iterData);
        }

        template< template< typename, typename > class SOLVER_T, typename FN_T, typename STRUCT_T, typename... Args >
        auto fsolve_common(FN_T func, const STRUCT_T& bounds, const Args&... args)
        {
            using TUPLE_T = std::tuple< Args... >;
            using SOLVER  = SOLVER_T< FN_T, StructCommonType_t< STRUCT_T > >;

            TUPLE_T args_tuple = std::make_tuple(args...);
            static_assert(std::tuple_size_v< TUPLE_T > <= 2, "Too many arguments passed to fsolve_common");

            // Zero arguments are passed...
            if constexpr (std::tuple_size_v< TUPLE_T > == 0) return detail::fsolve_impl(SOLVER(func, bounds), BracketTerminator {});

            // One argument is passed...
            else if constexpr (std::tuple_size_v< TUPLE_T > == 1) {
                using ArgType = std::tuple_element_t< 0, TUPLE_T >;

                if constexpr (std::is_floating_point_v< ArgType > || std::is_integral_v< ArgType >) {
                    return detail::fsolve_impl(SOLVER(func, bounds), BracketTerminator(std::get< 0 >(args_tuple)));
                }

                else if constexpr (std::is_same_v< std::invoke_result_t< ArgType, IterData< size_t, StructCommonType_t< STRUCT_T > > >,
                                                   bool >)
                {
                    return detail::fsolve_impl(SOLVER(func, bounds), std::get< 0 >(args_tuple));
                }
                else
                    []< bool flag = false >() { static_assert(flag, "Invalid argument passed to fsolve_common"); }
                ();
            }

            // Two arguments are passed...
            else if constexpr (std::tuple_size_v< TUPLE_T > == 2) {
                // Unpack and use the two arguments
                using ArgType1 = std::tuple_element_t< 0, decltype(args_tuple) >;
                using ArgType2 = std::tuple_element_t< 1, decltype(args_tuple) >;
                static_assert(std::is_floating_point_v< ArgType1 > != std::is_floating_point_v< ArgType2 >,
                              "Two arguments must be one floating point and one integral type");
                return detail::fsolve_impl(SOLVER(func, bounds), BracketTerminator(std::get< 0 >(args_tuple), std::get< 1 >(args_tuple)));
            }
            else
                []< bool flag = false >() { static_assert(flag, "Invalid argument passed to fsolve_common"); }
            ();
        }

    }    // namespace detail

    /**
     * @brief Defines a high-level root finding function template `fsolve` using bracketing solvers.
     *
     * The `fsolve` function template provides a convenient interface for performing root finding
     * using various bracketing solver algorithms. It abstracts the creation and configuration of the
     * solver instance and then delegates the actual root finding process to `fsolve_impl`. This
     * function is templated to accept a solver type, the function, bounds for the root, a tolerance
     * for convergence, and a maximum number of iterations. It supports different types of bracketing
     * solvers, making it versatile for various root finding needs.
     *
     * @tparam SOLVER_T The template class of the solver to be used. Must be a valid bracketing solver type.
     * @tparam FN_T The type of the function for which the root is being bracketed.
     * @tparam STRUCT_T The struct type holding the bounds for the root.
     * @tparam EPS_T The type of the epsilon value for convergence check, defaulted based on STRUCT_T.
     * @tparam ITER_T The type of the maximum iterations count, defaulted to int.
     *
     * @param function The function object for which to find the root.
     * @param bounds A struct with two members representing the lower and upper bounds.
     * @param eps The tolerance for stopping the algorithm.
     * @param maxiter The maximum number of iterations allowed.
     */
    // template< template< typename, typename > class SOLVER_T,
    //           IsFloatInvocable FN_T,
    //           IsFloatStruct    STRUCT_T,
    //           IsFloat          EPS_T  = StructCommonType_t< STRUCT_T >,
    //           std::integral    ITER_T = int >
    // auto fsolve(FN_T     function,
    //             STRUCT_T bounds,
    //             EPS_T    eps     = epsilon< StructCommonType_t< STRUCT_T > >(),
    //             ITER_T   maxiter = iterations< StructCommonType_t< STRUCT_T > >())
    // {
    //     using ARG_T = StructCommonType_t< STRUCT_T >;            /**< Common type for the bounds. */
    //     auto solver = SOLVER_T< FN_T, ARG_T >(function, bounds); /**< Instantiates the solver with the given function. */
    //
    //     // Delegates the solving process to fsolve_impl, passing in the solver and other parameters.
    //     return detail::fsolve_impl(solver, eps, maxiter);
    // }

    /**
     * @brief Extends the high-level root finding function template `fsolve` for initializer list bounds.
     *
     * This version of `fsolve` function template allows for specifying the bounds using an initializer list.
     * It is particularly useful when the bounds are known at compile time or for concise inline specifications.
     * The function checks the size of the initializer list to ensure exactly two elements are provided for the
     * bounds. It then creates a solver instance and delegates the root finding process to `fsolve_impl`.
     *
     * @tparam SOLVER_T The template class of the solver to be used. Must be a valid bracketing solver type.
     * @tparam FN_T The type of the function for which the root is being bracketed.
     * @tparam ARG_T The type of the bounds and the argument to the function.
     * @tparam EPS_T The type of the epsilon value for convergence check, defaulted based on ARG_T.
     * @tparam ITER_T The type of the maximum iterations count, defaulted to int.
     *
     * @param function The function object for which to find the root.
     * @param bounds Array representing the initial bounds around the root.
     * @param eps The tolerance for stopping the algorithm.
     * @param maxiter The maximum number of iterations allowed.
     */
    // template< template< typename, typename > class SOLVER_T,
    //           IsFloatInvocable FN_T,
    //           IsFloat          ARG_T,
    //           IsFloat          EPS_T  = ARG_T,
    //           std::integral    ITER_T = int,
    //           size_t           N >
    // requires(N == 2)
    // auto fsolve(FN_T function, const ARG_T (&bounds)[N], EPS_T eps = epsilon< ARG_T >(), ITER_T maxiter = iterations< ARG_T >())
    // {
    //     auto solver =
    //         SOLVER_T< FN_T, ARG_T >(function, std::pair(bounds[0], bounds[1])); /**< Instantiates the solver with the given function. */
    //
    //     // Delegates the solving process to fsolve_impl, passing in the solver and other parameters.
    //     return detail::fsolve_impl(solver, eps, maxiter);
    // }

    template< template< typename, typename > class SOLVER_T, IsFloatOrComplexInvocable FN_T, IsFloatStruct STRUCT_T, typename... ARGS >
    auto fsolve(FN_T func, STRUCT_T bounds, ARGS... args)
    {
        return detail::fsolve_common< SOLVER_T >(func, bounds, args...);
    }

    template< template< typename, typename > class SOLVER_T, IsFloatOrComplexInvocable FN_T, IsFloat ARG_T, size_t N, typename... ARGS >
    requires(N == 2)
    auto fsolve(FN_T func, const ARG_T (&bounds)[N], ARGS... args)
    {
        return detail::fsolve_common< SOLVER_T >(func, std::pair { bounds[0], bounds[1] }, args...);
    }

}    // namespace nxx::roots

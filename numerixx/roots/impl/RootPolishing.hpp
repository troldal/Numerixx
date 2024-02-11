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

#include <impl/Polynomial.hpp>

// ===== Standard Library Includes

namespace nxx::roots {
    // =================================================================================================================
    //
    //  88888888ba               88  88             88           88
    //  88      "8b              88  ""             88           ""
    //  88      ,8P              88                 88
    //  88aaaaaa8P'  ,adPPYba,   88  88  ,adPPYba,  88,dPPYba,   88  8b,dPPYba,    ,adPPYb,d8
    //  88""""""'   a8"     "8a  88  88  I8[    ""  88P'    "8a  88  88P'   `"8a  a8"    `Y88
    //  88          8b       d8  88  88   `"Y8ba,   88       88  88  88       88  8b       88
    //  88          "8a,   ,a8"  88  88  aa    ]8I  88       88  88  88       88  "8a,   ,d88
    //  88           `"YbbdP"'   88  88  `"YbbdP"'  88       88  88  88       88   `"YbbdP"Y8
    //                                                                             aa,    ,88
    //                                                                              "Y8bbdP"
    // =================================================================================================================

    /*
     * Private implementation details.
     */
    namespace detail {
        /**
         * @brief Provides a base class template for root polishing algorithms.
         *
         * The PolishingBase class template serves as a foundational component for
         * algorithms that refine or 'polish' roots of a given function. It encapsulates
         * common functionalities such as storing the objective function, its derivative,
         * and the current guess of the root. This class enforces certain type constraints
         * on the template parameters to ensure compatibility with root polishing algorithms.
         *
         * @tparam SUBCLASS The subclass inheriting from PolishingBase.
         * @tparam FUNCTION_T The type of the function for which the root is being polished.
         * @tparam DERIV_T The type of the derivative function of FUNCTION_T.
         * @tparam ARG_T The type of the argument to the function and its derivative.
         */
        template<typename SUBCLASS, typename FUNCTION_T, typename DERIV_T, typename ARG_T>
        requires std::same_as<typename PolishingTraits<SUBCLASS>::FUNCTION_T, FUNCTION_T>
                 && std::same_as<typename PolishingTraits<SUBCLASS>::DERIV_T, DERIV_T>
                 && IsFloatOrComplexInvocable<FUNCTION_T> && IsFloatOrComplexInvocable<DERIV_T>
                 && IsFloatOrComplex<ARG_T>
        class PolishingBase
        {
            friend SUBCLASS;

          public:
            static constexpr bool IsPolishingSolver = true; /**< Flag indicating the class is a polishing solver. */

            using FUNCT_RES_T = std::invoke_result_t<FUNCTION_T, ARG_T>; /**< Result type of the function. */
            using DERIV_RES_T = std::invoke_result_t<DERIV_T, ARG_T>; /**< Result type of the derivative function. */
            using RESULT_T = std::common_type_t<FUNCT_RES_T,
                DERIV_RES_T>; /**< Common type for results of function and derivative. */

          protected:
            ~PolishingBase() = default; /**< Protected destructor to prevent direct instantiation. */

          private:
            FUNCTION_T m_func{}; /**< The function object to find the root for. */
            DERIV_T m_deriv{}; /**< The function object for the derivative. */
            RESULT_T m_guess; /**< The current root estimate. */

          public:
            /**
             * @brief Constructs the PolishingBase with a function, its derivative, and an initial guess.
             * @param objective The function for which the root is being refined.
             * @param derivative The derivative of the objective function.
             * @param guess Initial guess for the root.
             */
            PolishingBase(FUNCTION_T objective, DERIV_T derivative, ARG_T guess)
              : m_func{ objective }, m_deriv{ derivative }, m_guess{ guess }
            {}

            PolishingBase(const PolishingBase &other) = default; /**< Default copy constructor. */
            PolishingBase(PolishingBase &&other) noexcept = default; /**< Default move constructor. */
            PolishingBase &operator=(const PolishingBase &other) = default; /**< Default copy assignment operator. */
            PolishingBase &operator=(
                PolishingBase &&other) noexcept = default; /**< Default move assignment operator. */

            /**
             * @brief Evaluates the function with the given value.
             * @tparam T Type of the value, must be a float or complex type.
             * @param value The value to evaluate the function at.
             * @return The result of the function evaluation.
             */
            template<typename T>
            requires nxx::IsFloatOrComplex<T>
            auto evaluate(T value)
            {
                return m_func(value);
            }

            /**
             * @brief Evaluates the derivative function with the given value.
             * @tparam T Type of the value, must be a float or complex type.
             * @param value The value to evaluate the derivative at.
             * @return The result of the derivative function evaluation.
             */
            template<typename T>
            requires nxx::IsFloatOrComplex<T>
            auto derivative(T value)
            {
                return std::invoke(m_deriv, value);
            }

            /**
             * @brief Returns the current result of the solver.
             * @details This function returns the current root estimate.
             *          It throws an exception if the solver has not been initialized.
             * @throws NumerixxError If the solver has not been initialized.
             * @return The current root estimate.
             */
            const RESULT_T &current() const { return m_guess; }

            void iterate() { std::invoke(static_cast<SUBCLASS &>(*this)); }
        };
    } // namespace detail

    // =================================================================================================================
    //
    //  888b      88
    //  8888b     88                                    ,d
    //  88 `8b    88                                    88
    //  88  `8b   88   ,adPPYba,  8b      db      d8  MM88MMM  ,adPPYba,   8b,dPPYba,
    //  88   `8b  88  a8P_____88  `8b    d88b    d8'    88    a8"     "8a  88P'   `"8a
    //  88    `8b 88  8PP"""""""   `8b  d8'`8b  d8'     88    8b       d8  88       88
    //  88     `8888  "8b,   ,aa    `8bd8'  `8bd8'      88,   "8a,   ,a8"  88       88
    //  88      `888   `"Ybbd8"'      YP      YP        "Y888  `"YbbdP"'   88       88
    //
    // =================================================================================================================


    /**
     * @brief Defines the Newton class for performing Newton's method root polishing.
     *
     * The Newton class template is a specialized implementation of the root polishing
     * algorithm known as Newton's method (or the Newton-Raphson method). It inherits
     * from a base class that provides common functionalities for root polishing algorithms,
     * and adds the specific iteration logic for Newton's method. This class is templated
     * to accept a function, its derivative, and an optional argument type.
     *
     * @tparam FN Type of the function for which the root is being refined.
     * @tparam DFN Type of the derivative function of FN.
     * @tparam ARG_T Type of the argument to the function and its derivative, defaults to double.
     */
    template<IsFloatOrComplexInvocable FN, IsFloatOrComplexInvocable DFN, IsFloatOrComplex ARG_T = double>
    class Newton final : public detail::PolishingBase<Newton<FN, DFN, ARG_T>, FN, DFN, ARG_T>
    {
        using BASE =
            detail::PolishingBase<Newton<FN, DFN, ARG_T>, FN, DFN, ARG_T>; /**< Base class alias for readability. */

      public:
        using BASE::BASE; /**< Inherits constructors from PolishingBase. */

        /**
         * @brief Performs a single iteration of Newton's method.
         * @details This method updates the root estimate using the Newton-Raphson formula.
         *          It assumes the class has been properly initialized.
         */
        void operator()()
        {
            BASE::m_guess = BASE::m_guess - BASE::evaluate(BASE::m_guess) / BASE::derivative(BASE::m_guess);
        }
    };

    /**
     * @brief Deduction guides for Newton class.
     * Allows the type of Newton class to be deduced from the constructor parameters.
     */
    template<typename FN, typename DERIV, typename ARG_T>
    requires IsFloatOrComplexInvocable<FN> && IsFloatOrComplexInvocable<DERIV> && IsFloatOrComplex<ARG_T>
    Newton(FN, DERIV, ARG_T) -> Newton<FN, DERIV, ARG_T>;


    // =================================================================================================================
    //
    //  ad88888ba
    // d8"     "8b                                                    ,d
    // Y8,                                                            88
    // `Y8aaaaa,     ,adPPYba,   ,adPPYba,  ,adPPYYba,  8b,dPPYba,  MM88MMM
    //   `"""""8b,  a8P_____88  a8"     ""  ""     `Y8  88P'   `"8a   88
    //         `8b  8PP"""""""  8b          ,adPPPPP88  88       88   88
    // Y8a     a8P  "8b,   ,aa  "8a,   ,aa  88,    ,88  88       88   88,
    //  "Y88888P"    `"Ybbd8"'   `"Ybbd8"'  `"8bbdP"Y8  88       88   "Y888
    //
    // =================================================================================================================

    /**
     * @brief Defines the Secant class for performing the secant method of root polishing.
     *
     * The Secant class template implements the secant method, an iterative root finding algorithm
     * that uses a succession of roots of secant lines to approximate a root of a function. This class
     * template inherits from a base class that provides common functionalities for root polishing algorithms,
     * and adds specific logic for the secant method iterations. It is designed to handle an initial guess
     * and a subsequent series of approximations to converge on the root. The secant method is a derivative-free
     * alternative to Newton's method and is particularly useful when the derivative of the function is not
     * readily available or is expensive to compute.
     *
     * @tparam FN The type of the function for which the root is being polished.
     * @tparam DFN The type of the derivative function of FN, used in the initial step.
     * @tparam ARG_T The type of the argument to the function, defaults to double.
     */
    template<IsFloatOrComplexInvocable FN, IsFloatOrComplexInvocable DFN, IsFloatOrComplex ARG_T = double>
    class Secant final : public detail::PolishingBase<Secant<FN, DFN, ARG_T>, FN, DFN, ARG_T>
    {
        using BASE =
            detail::PolishingBase<Secant<FN, DFN, ARG_T>, FN, DFN, ARG_T>; /**< Base class alias for readability. */

        ARG_T m_prevGuess; /**< Stores the previous guess for the root. */
        bool m_hasPrevGuess{ false }; /**< Flag to indicate whether a previous guess is available. */
        bool m_firstStep{ true }; /**< Flag to indicate whether the first step is to be taken. */

      public:
        using BASE::BASE; /**< Inherits constructors from PolishingBase. */

        /**
         * @brief Performs a single iteration of the secant method.
         * @details This method switches between a Newton-Raphson step for the first iteration
         *          and a Secant step for subsequent iterations.
         */
        void operator()()
        {
            if (m_firstStep) {
                ARG_T f_x = BASE::evaluate(BASE::m_guess);
                ARG_T f_prime_x = BASE::derivative(BASE::m_guess);

                if (abs(f_prime_x) < std::numeric_limits<ARG_T>::epsilon()) return;

                m_prevGuess = BASE::m_guess;
                BASE::m_guess -= f_x / f_prime_x;
                m_firstStep = false;
                m_hasPrevGuess = true;
            } else {
                ARG_T f_x = BASE::evaluate(BASE::m_guess);
                ARG_T f_x_prev = BASE::evaluate(m_prevGuess);

                if (abs(f_x - f_x_prev) < std::numeric_limits<ARG_T>::epsilon()) return;

                ARG_T newGuess = BASE::m_guess - f_x * (BASE::m_guess - m_prevGuess) / (f_x - f_x_prev);
                m_prevGuess = BASE::m_guess;
                BASE::m_guess = newGuess;
            }
        }
    };

    /**
     * @brief Deduction guides for Secant class.
     * Allows the type of Secant class to be deduced from the constructor parameters.
     */
    template<typename FN, typename DFN, typename ARG_T>
    requires IsFloatOrComplexInvocable<FN> && IsFloatOrComplexInvocable<DFN> && IsFloatOrComplex<ARG_T>
    Secant(FN, DFN, ARG_T) -> Secant<FN, DFN, ARG_T>;


    // =================================================================================================================
    //
    //  ad88888ba                        ad88     ad88
    // d8"     "8b  ,d                  d8"      d8"
    // Y8,          88                  88       88
    // `Y8aaaaa,  MM88MMM  ,adPPYba,  MM88MMM  MM88MMM  ,adPPYba,  8b,dPPYba,   ,adPPYba,   ,adPPYba,  8b,dPPYba,
    //   `"""""8b,  88    a8P_____88    88       88    a8P_____88  88P'   `"8a  I8[    ""  a8P_____88  88P'   `"8a
    //         `8b  88    8PP"""""""    88       88    8PP"""""""  88       88   `"Y8ba,   8PP"""""""  88       88
    // Y8a     a8P  88,   "8b,   ,aa    88       88    "8b,   ,aa  88       88  aa    ]8I  "8b,   ,aa  88       88
    //  "Y88888P"   "Y888  `"Ybbd8"'    88       88     `"Ybbd8"'  88       88  `"YbbdP"'   `"Ybbd8"'  88       88
    //
    // =================================================================================================================

    /**
     * @brief Defines the Steffensen class for performing Steffensen's method of root polishing.
     *
     * Steffensen's method is an iterative root finding algorithm that improves upon the simple fixed-point
     * iteration by incorporating a form of Aitken's Δ² process. This class template, `Steffensen`, inherits
     * from a base class that provides common functionalities for root polishing algorithms, and adds the
     * specific logic for Steffensen's method iterations. It starts with a Newton-Raphson step and then
     * switches to Steffensen's method for subsequent iterations. This method is particularly effective
     * for functions where the derivative is difficult to compute or is not readily available.
     *
     * @tparam FN The type of the function for which the root is being polished.
     * @tparam DFN The type of the derivative function of FN, used in the initial step.
     * @tparam ARG_T The type of the argument to the function, defaults to double.
     */
    template<IsFloatOrComplexInvocable FN, IsFloatOrComplexInvocable DFN, IsFloatOrComplex ARG_T = double>
    class Steffensen final : public detail::PolishingBase<Steffensen<FN, DFN, ARG_T>, FN, DFN, ARG_T>
    {
        using BASE =
            detail::PolishingBase<Steffensen<FN, DFN, ARG_T>, FN, DFN, ARG_T>; /**< Base class alias for readability. */

        bool m_firstStep{ true }; /**< Flag to indicate whether the first step is to be taken. */

      public:
        using BASE::BASE; /**< Inherits constructors from PolishingBase. */

        /**
         * @brief Performs a single iteration of the hybrid Steffensen method.
         * @details Uses Newton-Raphson for the first iteration and Steffensen's method subsequently.
         * @throws std::runtime_error If a division by near-zero occurs.
         */
        void operator()()
        {
            if (m_firstStep) {
                // Perform a Newton-Raphson step for the first iteration.
                ARG_T f_x = BASE::evaluate(BASE::m_guess);
                ARG_T f_prime_x = BASE::derivative(BASE::m_guess);

                if (abs(f_prime_x) < std::numeric_limits<ARG_T>::epsilon()) {
                    throw std::runtime_error("Division by near-zero in Newton-Raphson step.");
                    // TODO: Return a tl::expected instead of throwing an exception.
                }

                BASE::m_guess -= f_x / f_prime_x;
                m_firstStep = false;
            } else {
                // Perform a Steffensen's method step for subsequent iterations.
                ARG_T x = BASE::m_guess;
                ARG_T fx = BASE::evaluate(x);

                ARG_T x1 = x + fx;
                ARG_T fx1 = BASE::evaluate(x1);
                ARG_T x2 = x1 + fx1;

                ARG_T denominator = fx1 - fx;
                if (abs(denominator) < std::numeric_limits<ARG_T>::epsilon()) {
                    throw std::runtime_error("Division by near-zero in Steffensen's method.");
                    // TODO: Return a tl::expected instead of throwing an exception.
                }

                BASE::m_guess = x - (fx * fx) / denominator;
            }
        }
    };

    /**
     * @brief Deduction guides for Steffensen class.
     * Allows the type of Steffensen class to be deduced from the constructor parameters.
     */
    template<typename FN, typename DFN, typename ARG_T>
    requires IsFloatOrComplexInvocable<FN> && IsFloatOrComplexInvocable<DFN> && IsFloatOrComplex<ARG_T>
    Steffensen(FN, DFN, ARG_T) -> Steffensen<FN, DFN, ARG_T>;


    // =================================================================================================================
    //
    //    ad88          88     ad88                          88
    //   d8"            88    d8"                            88
    //   88             88    88                             88
    // MM88MMM  ,adPPYb,88  MM88MMM  ,adPPYba,   ,adPPYba,   88  8b       d8   ,adPPYba,
    //   88    a8"    `Y88    88     I8[    ""  a8"     "8a  88  `8b     d8'  a8P_____88
    //   88    8b       88    88      `"Y8ba,   8b       d8  88   `8b   d8'   8PP"""""""
    //   88    "8a,   ,d88    88     aa    ]8I  "8a,   ,a8"  88    `8b,d8'    "8b,   ,aa
    //   88     `"8bbdP"Y8    88     `"YbbdP"'   `"YbbdP"'   88      "8"       `"Ybbd8"'
    //
    // =================================================================================================================

    /**
     * @brief A struct template that holds the iteration data for a root finding problem.
     *
     * The PolishingIterData struct template is designed to store the current iteration count
     * and the current guess for the root in a root finding problem. It is used in conjunction
     * with root polishing algorithms to keep track of the current state of the algorithm.
     *
     * @tparam ITER_T The type of the iteration count. Must be an integral type.
     * @tparam RESULT_T The type of the root estimate. Must be a floating point type.
     */
    template<std::integral ITER_T, IsFloatOrComplex RESULT_T>
    struct PolishingIterData
    {
        ITER_T iter; /**< The current iteration count. */
        RESULT_T guess; /**< The current guess for the root. */
        std::vector<RESULT_T> previous{};
    };

    /**
     * @brief A class template that defines a termination condition for root polishing algorithms.
     *
     * The PolishingStopToken class template provides a mechanism to define when a root polishing
     * algorithm should stop iterating. It encapsulates an epsilon value and a maximum iteration count,
     * and provides a callable interface to check if the termination condition is met. The termination
     * condition is met when the current iteration count exceeds the maximum iteration count.
     *
     * @tparam EPS_T The type of the epsilon value. Must be a floating point type.
     * @tparam ITER_T The type of the iteration count. Must be an integral type.
     */
    template<IsFloat EPS_T, std::integral ITER_T>
    class PolishingStopToken
    {
        EPS_T m_eps; /**< The epsilon value for the termination condition. */
        ITER_T m_maxiter; /**< The maximum iteration count for the termination condition. */

      public:
        /**
         * @brief Default constructor. Initializes the epsilon value and maximum iteration count to default values.
         */
        explicit PolishingStopToken() : m_eps(epsilon<double>()), m_maxiter(iterations<double>()) {}

        /**
         * @brief Constructor. Initializes the epsilon value and maximum iteration count to the specified values.
         * @param eps The epsilon value for the termination condition.
         * @param maxiter The maximum iteration count for the termination condition.
         */
        explicit PolishingStopToken(EPS_T eps, ITER_T maxiter) : m_eps(eps), m_maxiter(maxiter) {}

        /**
         * @brief Constructor. Initializes the epsilon value and maximum iteration count to the specified values.
         * @param maxiter The maximum iteration count for the termination condition.
         * @param eps The epsilon value for the termination condition.
         */
        explicit PolishingStopToken(ITER_T maxiter, EPS_T eps) : m_eps(eps), m_maxiter(maxiter) {}

        /**
         * @brief Constructor. Initializes the epsilon value to the specified value and the maximum iteration count to a
         * default value.
         * @param eps The epsilon value for the termination condition.
         */
        explicit PolishingStopToken(EPS_T eps) : m_eps(eps), m_maxiter(iterations<double>()) {}

        /**
         * @brief Constructor. Initializes the maximum iteration count to the specified value and the epsilon value to a
         * default value.
         * @param maxiter The maximum iteration count for the termination condition.
         */
        explicit PolishingStopToken(ITER_T maxiter) : m_eps(epsilon<double>()), m_maxiter(maxiter) {}

        /**
         * @brief Checks if the termination condition is met.
         * @param data The current iteration data.
         * @return true if the termination condition is met, false otherwise.
         */
        bool operator()(const auto &data) const
        {
            const auto &[iter, guess, previous] = data;

            if (!previous.empty() && abs(guess - previous.back()) <= m_eps * abs(guess) + m_eps / 2) return true;
            if (iter >= m_maxiter) return true;

            return false;
        }
    };

    // Deduction guides for PolishingStopToken class.
    PolishingStopToken() -> PolishingStopToken<double, size_t>;
    PolishingStopToken(IsFloat auto eps) -> PolishingStopToken<decltype(eps), size_t>;
    PolishingStopToken(std::integral auto maxiter) -> PolishingStopToken<double, decltype(maxiter)>;
    PolishingStopToken(IsFloat auto eps, std::integral auto maxiter)
        -> PolishingStopToken<decltype(eps), decltype(maxiter)>;
    PolishingStopToken(std::integral auto maxiter, IsFloat auto eps)
        -> PolishingStopToken<decltype(eps), decltype(maxiter)>;

    namespace detail {

        /**
         * @brief A class template that encapsulates the result of a root finding problem.
         *
         * The PolishingSolverResult class template is designed to hold the result of a root finding problem
         * solved using a polishing solver. It stores an IterData object that contains the final iteration count
         * and the final root estimate. The class provides methods to retrieve the result in different formats,
         * depending on the needs of the user. The class is non-copyable and non-movable to ensure the integrity
         * of the result data.
         *
         * @tparam ITER_T The type of the iteration count. Must be an integral type.
         * @tparam RESULT_T The type of the root estimate. Must be a floating point type.
         */
        template<std::integral ITER_T, IsFloatOrComplex RESULT_T>
        class PolishingSolverResult
        {
            PolishingIterData<ITER_T, RESULT_T>
                m_iterData; /**< The IterData object holding the result of the root-finding problem. */

          public:
            /**
             * @brief Constructor that initializes the IterData object with the provided data.
             * @param iterData The IterData object that contains the final iteration count and root estimate.
             */
            explicit PolishingSolverResult(PolishingIterData<ITER_T, RESULT_T> iterData) : m_iterData(iterData) {}

            PolishingSolverResult(const PolishingSolverResult &) = delete; // No copy constructor
            PolishingSolverResult(PolishingSolverResult &&) = delete; // No move constructor

            PolishingSolverResult &operator=(const PolishingSolverResult &) = delete; // No copy assignment
            PolishingSolverResult &operator=(PolishingSolverResult &&) = delete; // No move assignment

            /**
             * @brief Retrieves the result of the root finding problem.
             * @tparam OUTPUT_T The type in which to return the result. Defaults to the type of the root estimate.
             * @return The root estimate. If OUTPUT_T is a class type, an object of that type is constructed with the
             * IterData.
             */
            template<typename OUTPUT_T = RESULT_T>
            auto result() &&
            {
                if constexpr (!std::constructible_from<OUTPUT_T, RESULT_T>)
                    return OUTPUT_T{}(m_iterData);
                else
                    return m_iterData.guess;
            }

            /**
             * @brief Retrieves the result of the root finding problem using a provided outputter function.
             * @tparam OUTPUTTER_T The type of the outputter function.
             * @param outputter The function to use to output the result.
             * @return The result of the root finding problem, outputted using the provided function.
             */
            template<typename OUTPUTTER_T>
            auto result(OUTPUTTER_T outputter) &&
            {
                return outputter(m_iterData);
            }
        };

        /**
         * @brief Deduction guide for PolishingSolverResult class.
         * Allows the type of PolishingSolverResult class to be deduced from the constructor parameters.
         */
        template<typename ITER_T, typename RESULT_T>
        PolishingSolverResult(PolishingSolverResult<ITER_T, RESULT_T>) -> PolishingSolverResult<ITER_T, RESULT_T>;

        /**
         * @brief Implements the root finding process for a given solver and termination token.
         *
         * This function template, `fdfsolve_impl`, provides a generic implementation for root finding
         * algorithms that utilize polishing solvers. It is designed to work with solvers that conform
         * to the requirements of polishing solvers. The function handles initialization, iteration, and
         * convergence checking, returning the result along with any potential errors encountered during
         * the solving process.
         *
         * @tparam SOLVER The type of the solver to be used in root finding. Must conform to the polishing solver
         *                concept.
         * @tparam TOKEN_T The type of the termination token. Must be a callable object that accepts an IterData object.
         *
         * @param solver The solver to be used in root finding.
         * @param terminator The termination token that determines when to stop the iteration process.
         *
         * @return The result of the root finding process, encapsulated in a PolishingSolverResult object.
         *
         * @note This function requires the termination token to be a callable object that accepts an IterData object.
         */
        template<typename SOLVER, typename TOKEN_T>
        requires SOLVER::IsPolishingSolver
        auto fdfsolve_impl(const SOLVER &solver, const TOKEN_T &terminator)
        {
            SOLVER _solver = solver;
            TOKEN_T _terminator = terminator;
            using ARG_T = typename SOLVER::RESULT_T;

            const auto &x = _solver.current();
            size_t iter = 0;

            PolishingIterData<size_t, ARG_T> iterData{ iter, x };

            while (true) {
                iterData.iter = iter;
                iterData.guess = x;


                if (_terminator(iterData)) break;
                // if (iter >= 100) break;
                iterData.previous.push_back(x);
                _solver.iterate();
                ++iter;
            }

            return PolishingSolverResult(iterData);
        }

        /**
         * @brief A function template that provides a common implementation for root finding algorithms.
         *
         * This function template, `fdfsolve_common`, provides a generic implementation for root finding
         * algorithms that utilize polishing solvers. It is designed to work with solvers that conform
         * to the requirements of polishing solvers. The function handles initialization, iteration, and
         * convergence checking, returning the result along with any potential errors encountered during
         * the solving process.
         *
         * @tparam SOLVER_T The template class of the solver to be used. Must be a valid polishing solver type.
         * @tparam FN_T The type of the function for which the root is being refined.
         * @tparam DERIV_T The type of the derivative function of FN_T.
         * @tparam ARG_T The type of the initial guess for the root.
         * @tparam Args The type of additional arguments passed to the function.
         *
         * @param func The function for which the root is being refined.
         * @param derivative The derivative of the function.
         * @param point The initial guess for the root.
         * @param args Additional arguments passed to the function.
         *
         * @return The result of the root finding process.
         *
         * @note This function requires the termination token to be a callable object that accepts an IterData object.
         */
        template<template<typename, typename, typename> class SOLVER_T,
            typename FN_T,
            typename DERIV_T,
            typename ARG_T,
            typename... Args>
        requires(sizeof...(Args) <= 2)
        auto fdfsolve_common(FN_T func, DERIV_T derivative, const ARG_T &point, const Args &...args)
        {
            using TUPLE_T = std::tuple<Args...>;
            using SOLVER = SOLVER_T<FN_T, DERIV_T, ARG_T>;
            using ITERDATA_T = PolishingIterData<size_t, ARG_T>;

            TUPLE_T args_tuple = std::make_tuple(args...);

            // Zero arguments are passed...
            if constexpr (std::tuple_size_v<TUPLE_T> == 0)
                return detail::fdfsolve_impl(SOLVER(func, derivative, point), PolishingStopToken{});

            // One argument is passed...
            else if constexpr (std::tuple_size_v<TUPLE_T> == 1) {
                using ArgType = std::tuple_element_t<0, TUPLE_T>;

                // If the argument is a floating point or integral type, use it as maxiter/eps
                if constexpr (std::is_floating_point_v<ArgType> || std::is_integral_v<ArgType>)
                    return detail::fdfsolve_impl(
                        SOLVER(func, derivative, point), PolishingStopToken(std::get<0>(args_tuple)));

                // If the argument is a callable, use as a stop token
                // TODO: Should be able to accept functors that take any king of argument, not just references.
                else if constexpr (std::same_as<std::invoke_result_t<ArgType, ITERDATA_T &>, bool>)
                    return detail::fdfsolve_impl(SOLVER(func, derivative, point), std::get<0>(args_tuple));

                else
                    std::invoke(
                        []<bool flag = false>() { static_assert(flag, "Invalid argument passed to fsolve_common"); });
            }

            // Two arguments are passed...
            else if constexpr (std::tuple_size_v<TUPLE_T> == 2) {
                // Unpack and use the two arguments
                using ArgType1 = std::tuple_element_t<0, decltype(args_tuple)>;
                using ArgType2 = std::tuple_element_t<1, decltype(args_tuple)>;
                static_assert(std::is_floating_point_v<ArgType1> != std::is_floating_point_v<ArgType2>,
                    "Two arguments must be one floating point and "
                    "one integral type");
                return detail::fdfsolve_impl(SOLVER(func, derivative, point),
                    PolishingStopToken(std::get<0>(args_tuple), std::get<1>(args_tuple)));
            } else
                std::invoke(
                    []<bool flag = false>() { static_assert(flag, "Invalid argument passed to fsolve_common"); });
        }
    } // namespace detail

    /**
     * @brief A function template that solves a root finding problem using a specified solver.
     *
     * This function template, `fdfsolve`, provides a generic implementation for root finding
     * algorithms that utilize polishing solvers. It is designed to work with solvers that conform
     * to the requirements of polishing solvers, such as having a defined `IsPolishingSolver` static member,
     * initialization, and iteration methods. The function handles initialization, iteration, and
     * convergence checking, returning the result along with any potential errors encountered during
     * the solving process.
     *
     * @tparam SOLVER_T The template class of the solver to be used. Must be a valid polishing solver type.
     * @tparam FN_T The type of the function for which the root is being refined.
     * @tparam DERIV_T The type of the derivative function of FN_T.
     * @tparam GUESS_T The type of the initial guess for the root.
     * @tparam ARGS The type of additional arguments passed to the function.
     *
     * @param func The function for which the root is being refined.
     * @param derivative The derivative of the function.
     * @param guess The initial guess for the root.
     * @param args Additional arguments passed to the function.
     *
     * @return The result of the root finding process.
     */
    template<template<typename, typename, typename> class SOLVER_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloatOrComplexInvocable DERIV_T,
        IsFloatOrComplex GUESS_T,
        typename... ARGS>
    auto fdfsolve(FN_T func, DERIV_T derivative, GUESS_T guess, ARGS... args)
    {
        return detail::fdfsolve_common<SOLVER_T>(func, derivative, guess, args...);
    }

    template<template<typename, typename, typename> class SOLVER_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloatOrComplex GUESS_T,
        typename... ARGS>
    auto fdfsolve(FN_T func, GUESS_T guess, ARGS... args)
    {
        using nxx::deriv::derivativeOf;
        using nxx::poly::derivativeOf;

        return detail::fdfsolve_common<SOLVER_T>(func, derivativeOf(func), guess, args...);
    }


    /**
     * @brief A function template that solves a root finding problem using a specified solver and termination token.
     *
     * This function template, `fdfsolve`, provides a generic implementation for root finding
     * algorithms that utilize polishing solvers. It is designed to work with solvers that conform
     * to the requirements of polishing solvers, such as having a defined `IsPolishingSolver` static member,
     * initialization, and iteration methods. The function handles initialization, iteration, and
     * convergence checking, returning the result along with any potential errors encountered during
     * the solving process.
     *
     * @tparam SOLVER_T The template class of the solver to be used. Must be a valid polishing solver type.
     * @tparam TOKEN_T The type of the termination token. Must be a callable object that accepts an IterData object.
     * @tparam FN_T The type of the function for which the root is being refined.
     * @tparam DERIV_T The type of the derivative function of FN_T.
     * @tparam GUESS_T The type of the initial guess for the root.
     *
     * @param func The function for which the root is being refined.
     * @param derivative The derivative of the function.
     * @param guess The initial guess for the root.
     *
     * @return The result of the root finding process.
     *
     * @note This function requires the termination token to be a callable object that accepts an IterData object.
     */
    template<template<typename, typename, typename> class SOLVER_T,
        typename TOKEN_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloatOrComplexInvocable DERIV_T,
        IsFloatOrComplex GUESS_T>
    // TODO: Should be able to accept functors that take any king of argument, not just references.
    requires std::invocable<TOKEN_T, PolishingIterData<size_t, GUESS_T> &>
    auto fdfsolve(FN_T func, DERIV_T derivative, GUESS_T guess)
    {
        return detail::fdfsolve_common<SOLVER_T>(func, derivative, guess, TOKEN_T{});
    }

    template<template<typename, typename, typename> class SOLVER_T,
        typename TOKEN_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloatOrComplex GUESS_T>
    // TODO: Should be able to accept functors that take any king of argument, not just references.
    requires std::invocable<TOKEN_T, PolishingIterData<size_t, GUESS_T> &>
    auto fdfsolve(FN_T func, GUESS_T guess)
    {
        using nxx::deriv::derivativeOf;
        using nxx::poly::derivativeOf;

        return detail::fdfsolve_common<SOLVER_T>(func, derivativeOf(func), guess, TOKEN_T{});
    }


} // namespace nxx::roots

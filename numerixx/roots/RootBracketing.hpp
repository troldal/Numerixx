/*
    o.     O O       o Oo      oO o.OOoOoo `OooOOo.  ooOoOOo o      O o      O
    Oo     o o       O O O    o o  O        o     `o    O     O    o   O    o
    O O    O O       o o  o  O  O  o        O      O    o      o  O     o  O
    O  o   o o       o O   Oo   O  ooOO     o     .O    O       oO       oO
    O   o  O o       O O        o  O        OOooOO'     o       Oo       Oo
    o    O O O       O o        O  o        o    o      O      o  o     o  o
    o     Oo `o     Oo o        O  O        O     O     O     O    O   O    O
    O     `o  `OoooO'O O        o ooOooOoO  O      o ooOOoOo O      o O      o

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

#ifndef NUMERIXX_ROOTBRACKETING_HPP
#define NUMERIXX_ROOTBRACKETING_HPP

// ===== Numerixx Includes
#include "RootCommon.hpp"
#include "calculus/Derivatives.hpp"

// ===== External Includes
#include "../.utils/Constants.hpp"

// ===== Standard Library Includes
#include <array>

namespace nxx::roots
{
    // ========================================================================
    // ROOT-FINDING WITHOUT DERIVATIVES (BRACKETING)
    // ========================================================================

    namespace impl
    {
        /**
         * @brief A base class for bracketing root-finding algorithms.
         *
         * This class provides a generic interface for bracketing root-finding algorithms. The actual algorithm
         * is implemented in a derived class, which is passed as a template argument (POLICY).
         *
         * @tparam POLICY The derived class implementing the specific root-finding algorithm.
         * @requires POLICY must be invocable with a floating point type as its argument type.
         */
        template< typename POLICY >
        requires std::floating_point< typename BracketingTraits< POLICY >::return_type >
        class BracketingBase
        {
            /*
             * Friend declarations.
             */
            friend POLICY;

        public:
            using function_type = typename BracketingTraits< POLICY >::function_type;
            using return_type   = typename BracketingTraits< POLICY >::return_type;

        protected:
            /**
             * @brief Default destructor.
             */
            ~BracketingBase() = default;

        private:
            function_type                         m_func {};                 /**< The function object to find the root for. */
            std::pair< return_type, return_type > m_bounds { 0.0, 0.0 };     /**< Holds the current bounds around the root. */
            bool                                  m_isInitialized { false }; /**< Indicates whether the solver has been initialized. */

            /**
             * @brief Constructor, taking a function object as an argument.
             *
             * @param objective The function object to find the root for.
             * @note Constructor is private to avoid direct usage by clients.
             */
            explicit BracketingBase(function_type objective)
                : m_func { std::move(objective) }
            {}

            /**
             * @brief Sets the current bounds around the root.
             *
             * @param bounds An object holding the new bounds around the root.
             */
            void setBounds(auto bounds)
            {
                if (!m_isInitialized) throw std::logic_error("Solver has not been initialized!");
                auto [lower, upper] = bounds;
                static_assert(std::floating_point< decltype(lower) >);
                m_bounds = std::pair< return_type, return_type > { lower, upper };
            }

            /**
             * @brief Sets the current bounds around the root.
             *
             * @param bounds std::initializer_list holding the new bounds around the root.
             */
            template< typename T >
            requires std::floating_point< T >
            constexpr void setBounds(std::initializer_list< T > bounds)
            {
                if (bounds.size() != 2) throw std::logic_error("Initializer list must contain exactly two elements!");
                setBounds(std::pair< return_type, return_type > { *bounds.begin(), *(bounds.begin() + 1) });
            }

        public:
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
             * @brief Initializes the solver by setting the initial bounds around the root.
             *
             * @param bounds An object holding the initial bounds around the root. The root must be contained
             * inside these bounds. The object must support structured bindings to provide two values: lower and upper bounds.
             * Examples of supported types include pairs, tuples, or custom structs with structured bindings support.
             * @warning If the solver is used without initializing, the behavior is undefined.
             */
            void init(auto bounds)
            {
                m_isInitialized     = true;
                auto [lower, upper] = bounds;
                setBounds(std::pair< return_type, return_type > { lower, upper });
            }

            /**
             * @brief Initializes the solver by setting the initial bounds around the root.
             *
             * @param bounds std::initializer_list holding the initial bounds around the root. The root must be contained
             * inside these bounds. The list must contain exactly two elements, which will be interpreted as the lower and
             * upper bounds, respectively.
             * @warning If the solver is used without initializing, the behavior is undefined.
             */
            template< typename T >
            requires std::floating_point< T >
            void init(std::initializer_list< T > bounds)
            {
                m_isInitialized = true;
                setBounds(bounds);
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
            auto evaluate(return_type value) { return m_func(value); }

            /**
             * @brief Returns the current bounds around the root.
             *
             * Every time an iteration is executed, the bounds will narrow. This function returns a const reference
             * to the current bounds as a std::pair. The value type of the bounds is the same as the return type
             * of the function object.
             *
             * @return A const reference to the current bounds.
             */
            const auto& bounds() const
            {
                if (!m_isInitialized) throw std::logic_error("Solver has not been initialized!");
                return m_bounds;
            }
        };
    }    // namespace impl

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
    template< typename FN >
    requires std::floating_point< std::invoke_result_t< FN, double > >
    class Ridder final : public impl::BracketingBase< Ridder< FN > >
    {
        /*
         * Private alias declarations.
         */
        using BASE = impl::BracketingBase< Ridder< FN > >;

    public:
        /*
         * Public alias declarations.
         */
        using function_type = FN;

        /**
         * @brief Constructor, taking the function object as an argument.
         *
         * @param objective The function object for which to find the root.
         * @note This constructor must call the BracketingBase constructor.
         */
        explicit Ridder(FN objective)
            : BASE { objective }
        {}

        /**
         * @brief Perform one iteration of Ridder's method.
         *
         * This function implements the main algorithm of Ridder's method for root-finding. It updates the
         * bounds around the root during each iteration, gradually narrowing the search interval.
         */
        void iterate()
        {
            const auto& bounds = BASE::bounds();
            using RT           = std::invoke_result_t< FN, decltype(bounds.first) >;

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

    Ridder(std::invocable auto func) -> Ridder< decltype(func) >;

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
    template< typename FN >
    requires std::floating_point< std::invoke_result_t< FN, double > >
    class Bisection final : public impl::BracketingBase< Bisection< FN > >
    {
        /*
         * Private alias declarations.
         */
        using BASE = impl::BracketingBase< Bisection< FN > >;

    public:
        /*
         * Public alias declarations.
         */
        using function_type = FN;

        /**
         * @brief Constructor, taking the function object as an argument.
         * @param objective The function object for which to find the root.
         * @note This constructor must call the BracketingBase constructor.
         */
        explicit Bisection(FN objective)
            : BASE { objective }
        {}

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
            using RT           = std::invoke_result_t< FN, decltype(bounds.first) >;

            RT root = (bounds.first + bounds.second) / 2.0;

            if (BASE::evaluate(bounds.first) * BASE::evaluate(root) < 0.0)
                BASE::setBounds({ bounds.first, root });
            else
                BASE::setBounds({ root, bounds.second });
        }
    };

    Bisection(std::invocable auto func) -> Bisection< decltype(func) >;

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
    template< typename FN >
    requires std::floating_point< std::invoke_result_t< FN, double > >
    class RegulaFalsi final : public impl::BracketingBase< RegulaFalsi< FN > >
    {
        /*
         * Private alias declarations.
         */
        using BASE = impl::BracketingBase< RegulaFalsi< FN > >;

    public:
        /*
         * Public alias declarations.
         */
        using function_type = FN;

        /**
         * @brief Constructor, taking the function object as an argument.
         *
         * @param objective The function object for which to find the root.
         * @note This constructor must call the BracketingBase constructor.
         */
        explicit RegulaFalsi(FN objective)
            : BASE { objective }
        {}

        /**
         * @brief Perform one iteration of the Regula Falsi algorithm.
         *
         * This function implements the main algorithm of the Regula Falsi method for root-finding.
         * It updates the bounds around the root during each iteration, refining the search interval.
         */
        void iterate()
        {
            const auto& bounds = BASE::bounds();
            using RT           = std::invoke_result_t< FN, decltype(bounds.first) >;

            RT f_lo = BASE::evaluate(bounds.first);
            RT f_hi = BASE::evaluate(bounds.second);

            RT root   = (bounds.first * f_hi - bounds.second * f_lo) / (f_hi - f_lo);
            RT f_root = BASE::evaluate(root);

            if (f_lo * f_root < 0.0) {
                BASE::setBounds({ bounds.first, root });
            }
            else {
                BASE::setBounds({ root, bounds.second });
            }
        }
    };

    RegulaFalsi(std::invocable auto func) -> RegulaFalsi< decltype(func) >;

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
        template< typename SOLVER >
        requires requires(SOLVER solver, std::pair< typename SOLVER::return_type, typename SOLVER::return_type > bounds) {
                     // clang-format off
                     { solver.evaluate(0.0) } -> std::floating_point;
                     { solver.init(bounds) };
                     { solver.iterate() };
                     // clang-format on
                 }
        inline auto fsolve_impl(SOLVER                                                                  solver,
                                std::pair< typename SOLVER::return_type, typename SOLVER::return_type > bounds,
                                typename SOLVER::return_type                                            eps     = nxx::EPS,
                                int                                                                     maxiter = nxx::MAXITER)
        {
            using ET = RootErrorImpl< typename SOLVER::return_type >;
            using RT = tl::expected< typename SOLVER::return_type, ET >;

            solver.init(bounds);

            // Declare variables for use in the iteration loop.
            auto curBounds = solver.bounds();
            RT   result    = (curBounds.first + curBounds.second) / 2.0;
            std::array< std::pair< typename SOLVER::return_type, typename SOLVER::return_type >, 2 > roots {};
            decltype(roots.begin())                                                                  min;

            // Check for NaN or Inf.
            if (!std::isfinite(solver.evaluate(curBounds.first)) || !std::isfinite(solver.evaluate(curBounds.second))) {
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
                if (std::any_of(roots.begin(), roots.end(), [](const auto& r) { return !std::isfinite(r.second); })) {
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
    }    // namespace impl

    /**
     * @brief Main fsolve function to find the root of the function with specified bounds.
     *
     * @tparam SOLVER The solver type used to find the root of the function.
     * @param solver The solver instance used to find the root of the function.
     * @param bounds Any object that supports structured bindings and provides two values for the lower and upper bounds.
     * @param eps The tolerance for stopping the algorithm.
     * @param maxiter The maximum number of iterations allowed.
     * @return tl::expected object containing either the root of the function or an error.
     */
    template< typename SOLVER >
    inline auto fsolve(SOLVER solver, auto bounds, typename SOLVER::return_type eps = nxx::EPS, int maxiter = nxx::MAXITER)
    {
        using RT      = typename SOLVER::return_type;
        auto [lo, hi] = bounds;
        return impl::fsolve_impl(solver, std::pair< RT, RT > { lo, hi }, eps, maxiter);
    }

    /**
     * @brief Overload of fsolve function that accepts an initializer list for bounds.
     *
     * @tparam SOLVER The solver type used to find the root of the function.
     * @tparam T The type of elements in the initializer list.
     * @param solver The solver instance used to find the root of the function.
     * @param bounds An initializer list containing exactly two elements representing the lower and upper bounds.
     * @param eps The tolerance for stopping the algorithm.
     * @param maxiter The maximum number of iterations allowed.
     * @return tl::expected object containing either the root of the function or an error.
     * @throws std::logic_error if the initializer list does not contain exactly two elements.
     */
    template< typename SOLVER, typename T >
    inline auto
        fsolve(SOLVER solver, std::initializer_list< T > bounds, typename SOLVER::return_type eps = nxx::EPS, int maxiter = nxx::MAXITER)
    {
        using RT = typename SOLVER::return_type;
        if (bounds.size() != 2) throw std::logic_error("Initializer list must contain exactly two elements!");
        return impl::fsolve_impl(solver, std::pair< RT, RT > { *bounds.begin(), *(bounds.begin() + 1) }, eps, maxiter);
    }

}    // namespace nxx::roots

#endif    // NUMERIXX_ROOTBRACKETING_HPP

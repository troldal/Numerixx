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

#ifndef NUMERIXX_ROOTSEARCHING_HPP
#define NUMERIXX_ROOTSEARCHING_HPP

// ===== Numerixx Includes
#include "RootCommon.hpp"

// ===== External Includes
#include "../.dependencies/expected/expected.hpp"
#include "../.utils/Constants.hpp"

// ===== Standard Library Includes
#include <numbers>
#include <optional>

namespace nxx::roots
{

    // ========================================================================
    // BRACKET SEARCHING
    // ========================================================================
    namespace impl
    {
        /**
         * @brief A template base class for search algorithms.
         *
         * This class serves as the base for various search algorithms, providing common
         * functionality and data members. It is designed to work with a specific search
         * policy provided as a template parameter. The policy should define the traits
         * used for proper instantiation of this class.
         *
         * @tparam POLICY The search policy implementation.
         */
        template< typename POLICY >
        requires std::floating_point< typename SearchingTraits< POLICY >::return_type >
        class SearchBase
        {
            /*
             * Friend declarations.
             */
            friend POLICY;

        public:
            using function_type = typename SearchingTraits< POLICY >::function_type;
            using return_type   = typename SearchingTraits< POLICY >::return_type;

        protected:
            /**
             * @brief Default constructor.
             */
            ~SearchBase() = default;

        private:
            std::pair< return_type, return_type > m_bounds { 0.0, 1.0 };
            std::optional< decltype(m_bounds) >   m_limits = std::nullopt;
            function_type                         m_objective {};
            return_type                           m_factor {};
            bool                                  m_isInitialized { false };

            /**
             * @brief Constructor with objective function.
             *
             * Initializes the search algorithm with the provided objective function.
             * The function is stored internally and used for evaluation during the search process.
             *
             * @param objective The objective function to optimize.
             */
            explicit SearchBase(function_type objective, std::optional< decltype(m_bounds) > limits, return_type factor)
                : m_objective(std::move(objective)),
                  m_limits(std::move(limits)),
                  m_factor(factor)
            {
                if (m_factor < 1.0) throw std::invalid_argument("Invalid factor.");
            }

            /**
             * @brief Set the search bounds.
             *
             * This method sets the search bounds, represented by a pair of lower and upper
             * bounds. If the lower bound is greater than the upper bound, the values will be
             * swapped internally to ensure a valid search space.
             *
             * @param bounds A pair of lower and upper bounds.
             * @throws std::invalid_argument If the lower and upper bounds are equal.
             */
            void setBounds(auto bounds)
            {
                if (!m_isInitialized) throw std::logic_error("Search algorithm not initialized!");
                auto [lower, upper] = bounds;
                static_assert(std::floating_point< decltype(lower) >);
                if (lower == upper) throw std::invalid_argument("Invalid bounds.");
                if (lower > upper) std::swap(lower, upper);

                if (m_limits.has_value()) {
                    auto [lowerLimit, upperLimit] = m_limits.value();
                    if (lower < lowerLimit) lower = lowerLimit;
                    if (upper > upperLimit) upper = upperLimit;
                }
                m_bounds = std::pair< return_type, return_type > { lower, upper };
            }

            /**
             * @brief Set the search bounds using an initializer list.
             *
             * This method sets the search bounds using an initializer list containing two elements.
             * The elements are treated as the lower and upper bounds for the search space.
             * If the first element is greater than the second element, the values will be
             * swapped internally to ensure a valid search space.
             *
             * @tparam T The floating-point type.
             * @param bounds An initializer list with exactly two elements.
             * @throws std::logic_error If the initializer list does not contain exactly two elements.
             */
            template< typename T >
            requires std::floating_point< T >
            void setBounds(std::initializer_list< T > bounds)
            {
                if (bounds.size() != 2) throw std::logic_error("Initializer list must contain exactly two elements!");
                setBounds(std::pair< return_type, return_type > { *bounds.begin(), *(bounds.begin() + 1) });
            }

            /**
             * @brief Set the search factor.
             *
             * This method sets the search factor, which is used to determine the next search
             * point. The search factor is multiplied with the current search point to obtain
             * the next search point. The search factor must be greater than 1.0.
             *
             * @param factor The search factor.
             * @throws std::invalid_argument If the search factor is less than or equal to 1.0.
             */
            void setFactor(auto factor) requires std::floating_point< decltype(factor) >
            {
                if (!m_isInitialized) throw std::logic_error("Search algorithm not initialized!");
                if (factor < 1.0) throw std::invalid_argument("Invalid factor.");
                m_factor = factor;
            }

        public:
            /**
             * @brief Copy constructor.
             *
             * @param other Another SearchBase object to be copied.
             */
            SearchBase(const SearchBase& other) = default;

            /**
             * @brief Move constructor.
             *
             * @param other Another SearchBase object to be moved.
             */
            SearchBase(SearchBase&& other) noexcept = default;

            /**
             * @brief Copy assignment operator.
             *
             * @param other Another SearchBase object to be copied.
             * @return A reference to this object.
             */
            SearchBase& operator=(const SearchBase& other) = default;

            /**
             * @brief Move assignment operator.
             *
             * @param other Another SearchBase object to be moved.
             * @return A reference to this object.
             */
            SearchBase& operator=(SearchBase&& other) noexcept = default;

            /**
             * @brief Returns the current search bounds.
             *
             * @return A const reference to a pair of lower and upper bounds.
             */
            [[nodiscard]]
            std::pair< double, double > bounds() const
            {
                if (!m_isInitialized) throw std::logic_error("Search algorithm not initialized!");
                return m_bounds;
            }

            /**
             * @brief Get the factor used in the search algorithm.
             *
             * This method returns the factor used by the search algorithm for certain
             * operations, such as adjusting the search bounds. The factor may be
             * different for different search policies and can affect the behavior of
             * the algorithm.
             *
             * @return The factor used in the search algorithm.
             */
            return_type factor() const
            {
                if (!m_isInitialized) throw std::logic_error("Search algorithm not initialized!");
                return m_factor;
            }

            /**
             * @brief Initialize the search bounds.
             *
             * This method is a shorthand for calling setBounds with a pair of lower and upper bounds.
             *
             * @param bounds A pair of lower and upper bounds.
             */
            void init(auto bounds, return_type factor = 0.0)
            {
                m_isInitialized = true;
                setBounds(bounds);
                if (factor >= 1.0) setFactor(factor);
            }

            /**
             * @brief Initialize the search bounds using an initializer list.
             *
             * This method is a shorthand for calling setBounds with an initializer list containing two elements.
             *
             * @tparam T The floating-point type.
             * @param bounds An initializer list with exactly two elements.
             */
            template< typename T >
            requires std::floating_point< T >
            void init(std::initializer_list< T > bounds, return_type factor = 0.0)
            {
                m_isInitialized = true;
                setBounds(bounds);
                if (factor >= 1.0) setFactor(factor);
            }

            /**
             * @brief Reset the search algorithm.
             *
             * This method resets the search algorithm to its initial state.
             */
            void reset() { m_isInitialized = false; }

            /**
             * @brief Evaluate the objective function at a given point.
             *
             * This method evaluates the stored objective function at the specified point `x` and
             * returns the result.
             *
             * @param value The point at which to evaluate the objective function.
             * @return The value of the objective function at the specified point.
             */
            [[nodiscard]]
            auto evaluate(double value) const
            {
                return m_objective(value);
            }
        };
    }    // namespace impl

    /**
     * @brief A class for bracket search up algorithm.
     *
     * The BracketSearchUp class implements a bracket search algorithm that moves the search
     * bounds upwards. The algorithm starts with an initial bracket and iteratively expands the
     * bracket until a root is found. The objective function is expected to be a continuous
     * function for which the algorithm tries to find a root.
     *
     * @tparam FN The objective function type.
     */
    template< typename FN >
    requires std::floating_point< std::invoke_result_t< FN, double > >
    class BracketSearchUp final : public impl::SearchBase< BracketSearchUp< FN > >
    {
        using BASE = impl::SearchBase< BracketSearchUp< FN > >;

    public:
        using function_type = FN;
        using return_type   = std::invoke_result_t< FN, double >;
        using limit_type = std::pair< return_type, return_type >;

        /**
         * @brief Constructor with objective function.
         *
         * Initializes the BracketSearchUp algorithm with the provided objective function.
         * The function is stored internally and used for evaluation during the search process.
         *
         * @param objective The objective function to optimize.
         */
        explicit BracketSearchUp(function_type objective, std::optional<limit_type> limits = std::nullopt, return_type factor = std::numbers::phi_v< return_type >)
            : BASE(std::move(objective), limits, factor)
        {}

        /**
         * @brief Perform a single iteration of the search algorithm.
         *
         * During each iteration, the search bounds are expanded upwards using the golden ratio.
         * The iteration stops when a root is found, i.e., when the objective function values
         * at the bounds have opposite signs.
         */
        void iterate()
        {
            const auto& bounds = BASE::bounds();

            if (BASE::evaluate(bounds.first) * BASE::evaluate(bounds.second) < 0.0) return;

            auto newBounds   = bounds;
            newBounds.first  = bounds.second;
            newBounds.second = bounds.second + (bounds.second - bounds.first) * BASE::factor();
            BASE::setBounds(newBounds);
        }
    };

    /**
     * @brief A class for bracket search down algorithm.
     *
     * The BracketSearchDown class implements a bracket search algorithm that moves the search
     * bounds downwards. The algorithm starts with an initial bracket and iteratively expands the
     * bracket until a root is found. The objective function is expected to be a continuous
     * function for which the algorithm tries to find a root.
     *
     * @tparam FN The objective function type.
     */
    template< typename FN >
    requires std::floating_point< std::invoke_result_t< FN, double > >
    class BracketSearchDown final : public impl::SearchBase< BracketSearchDown< FN > >
    {
        using BASE = impl::SearchBase< BracketSearchDown< FN > >;

    public:
        using function_type = FN;
        using return_type   = std::invoke_result_t< FN, double >;
        using limit_type = std::pair< return_type, return_type >;

        /**
         * @brief Constructor with objective function.
         *
         * Initializes the BracketSearchDown algorithm with the provided objective function.
         * The function is stored internally and used for evaluation during the search process.
         *
         * @param objective The objective function to optimize.
         */
        explicit BracketSearchDown(function_type objective, std::optional<limit_type> limits = std::nullopt, return_type factor = std::numbers::phi_v< return_type >)
            : BASE(std::move(objective), limits, factor)
        {}

        /**
         * @brief Perform a single iteration of the search algorithm.
         *
         * During each iteration, the search bounds are expanded downwards using the golden ratio.
         * The iteration stops when a root is found, i.e., when the objective function values
         * at the bounds have opposite signs.
         */
        void iterate()
        {
            const auto& bounds = BASE::bounds();

            if (BASE::evaluate(bounds.first) * BASE::evaluate(bounds.second) < 0.0) return;

            auto newBounds   = bounds;
            newBounds.second = bounds.first;
            newBounds.first  = bounds.first - (bounds.second - bounds.first) * BASE::factor();
            BASE::setBounds(newBounds);
        }
    };

    /**
     * @brief A class for bracket expand up algorithm.
     *
     * The BracketExpandUp class implements a bracket expansion algorithm that expands the search
     * bounds upwards. The algorithm starts with an initial bracket and iteratively expands the
     * bracket until a root is found. The objective function is expected to be a continuous
     * function for which the algorithm tries to find a root.
     *
     * @tparam FN The objective function type.
     */
    template< typename FN >
    requires std::floating_point< std::invoke_result_t< FN, double > >
    class BracketExpandUp final : public impl::SearchBase< BracketExpandUp< FN > >
    {
        using BASE = impl::SearchBase< BracketExpandUp< FN > >;

    public:
        using function_type = FN;
        using return_type   = std::invoke_result_t< FN, double >;
        using limit_type = std::pair< return_type, return_type >;

        /**
         * @brief Constructor with objective function.
         *
         * Initializes the BracketExpandUp algorithm with the provided objective function.
         * The function is stored internally and used for evaluation during the search process.
         *
         * @param objective The objective function to optimize.
         */
        explicit BracketExpandUp(function_type objective, std::optional<limit_type> limits = std::nullopt, return_type factor = std::numbers::phi_v< return_type >)
            : BASE(std::move(objective), limits, factor)
        {}

        /**
         * @brief Perform a single iteration of the search algorithm.
         *
         * During each iteration, the search bounds are expanded upwards using the golden ratio.
         * The iteration stops when a root is found, i.e., when the objective function values
         * at the bounds have opposite signs.
         */
        void iterate()
        {
            const auto& bounds = BASE::bounds();

            if (BASE::evaluate(bounds.first) * BASE::evaluate(bounds.second) < 0.0) return;

            auto newBounds   = bounds;
            newBounds.second = bounds.second + (bounds.second - bounds.first) * BASE::factor();
            BASE::setBounds(newBounds);
        }
    };

    /**
     * @brief A class for bracket expand down algorithm.
     *
     * The BracketExpandDown class implements a bracket expansion algorithm that expands the search
     * bounds downwards. The algorithm starts with an initial bracket and iteratively expands the
     * bracket until a root is found. The objective function is expected to be a continuous
     * function for which the algorithm tries to find a root.
     *
     * @tparam FN The objective function type.
     */
    template< typename FN >
    requires std::floating_point< std::invoke_result_t< FN, double > >
    class BracketExpandDown final : public impl::SearchBase< BracketExpandDown< FN > >
    {
        using BASE = impl::SearchBase< BracketExpandDown< FN > >;

    public:
        using function_type = FN;
        using return_type   = std::invoke_result_t< FN, double >;
        using limit_type = std::pair< return_type, return_type >;

        /**
         * @brief Constructor with objective function.
         *
         * Initializes the BracketExpandDown algorithm with the provided objective function.
         * The function is stored internally and used for evaluation during the search process.
         *
         * @param objective The objective function to optimize.
         */
        explicit BracketExpandDown(function_type objective, std::optional<limit_type> limits = std::nullopt, return_type factor = std::numbers::phi_v< return_type >)
            : BASE(std::move(objective), limits, factor)
        {}

        /**
         * @brief Perform a single iteration of the search algorithm.
         *
         * During each iteration, the search bounds are expanded downwards using the golden ratio.
         * The iteration stops when a root is found, i.e., when the objective function values
         * at the bounds have opposite signs.
         */
        void iterate()
        {
            const auto& bounds = BASE::bounds();

            if (BASE::evaluate(bounds.first) * BASE::evaluate(bounds.second) < 0.0) return;

            auto newBounds  = bounds;
            newBounds.first = bounds.first - (bounds.second - bounds.first) * BASE::factor();
            BASE::setBounds(newBounds);
        }
    };

    template< typename FN >
    requires std::floating_point< std::invoke_result_t< FN, double > >
    class BracketExpandOut final : public impl::SearchBase< BracketExpandOut< FN > >
    {
        using BASE = impl::SearchBase< BracketExpandOut< FN > >;

    public:
        using function_type = FN;
        using return_type   = std::invoke_result_t< FN, double >;
        using limit_type = std::pair< return_type, return_type >;

        /**
         * @brief Constructor with objective function.
         *
         * Initializes the BracketExpandDown algorithm with the provided objective function.
         * The function is stored internally and used for evaluation during the search process.
         *
         * @param objective The objective function to optimize.
         */
        explicit BracketExpandOut(function_type objective, std::optional<limit_type> limits = std::nullopt, return_type factor = std::numbers::phi_v< return_type >)
            : BASE(std::move(objective), limits, factor)
        {}

        /**
         * @brief Perform a single iteration of the search algorithm.
         *
         * During each iteration, the search bounds are expanded downwards using the golden ratio.
         * The iteration stops when a root is found, i.e., when the objective function values
         * at the bounds have opposite signs.
         */
        void iterate()
        {
            const auto& bounds = BASE::bounds();

            if (BASE::evaluate(bounds.first) * BASE::evaluate(bounds.second) < 0.0) return;

            auto newBounds   = bounds;
            newBounds.first  = bounds.first - (bounds.second - bounds.first) * BASE::factor() / 2.0;
            newBounds.second = bounds.second + (bounds.second - bounds.first) * BASE::factor() / 2.0;
            BASE::setBounds(newBounds);
        }
    };

    /**
     * @brief A class for root finding of a function within an interval using a subdivide and search approach.
     *
     * This class is designed for finding the root of a function within an interval by iteratively
     * subdividing the interval and searching for the subinterval where the root lies. The process
     * continues until the root is found or the number of subdivisions exceeds a specified threshold.
     * The class inherits from the SearchBase class.
     *
     * @tparam FN A function object type representing the objective function to be optimized. The
     *            function must have a signature compatible with `double function(double)`.
     */
    template< typename FN >
    requires std::floating_point< std::invoke_result_t< FN, double > >
    class BracketSubdivide final : public impl::SearchBase< BracketSubdivide< FN > >
    {
        using BASE = impl::SearchBase< BracketSubdivide< FN > >;

    public:
        using function_type = FN;
        using return_type   = std::invoke_result_t< FN, double >;

        /**
         * @brief Constructor with objective function and optional factor.
         *
         * Initializes the search algorithm with the provided objective function and an optional
         * factor determining the number of subdivisions of the search interval.
         *
         * @param objective The objective function to optimize.
         * @param factor The initial number of subdivisions of the search interval (default is 1.0).
         */
        explicit BracketSubdivide(function_type objective, return_type factor = std::numbers::phi_v< return_type >)
            : BASE(std::move(objective), std::nullopt, factor)
        {}

        /**
         * @brief Iterates the search algorithm.
         *
         * This method subdivides the current search interval into smaller subintervals based on the
         * current factor value. It then searches for the subinterval where the root lies, updates the
         * bounds accordingly, and returns. If the root is not found in any of the subintervals, the
         * factor is doubled, and the search continues in the next iteration.
         */
        void iterate()
        {
            const auto& bounds = BASE::bounds();
            if (BASE::evaluate(bounds.first) * BASE::evaluate(bounds.second) < 0.0) return;

            size_t factor = std::ceil(BASE::factor());
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

            BASE::setFactor(BASE::factor() * 2.0);
        }
    };

    namespace impl
    {
        /**
         * @brief The main search implementation function.
         *
         * This function template implements the search algorithm for a root finding solver.
         * The SOLVER type must meet specific requirements defined by the requires expression.
         *
         * @tparam SOLVER The solver type to use for root finding.
         * @param solver The solver instance used to perform root finding.
         * @param bounds The initial bounds of the search interval.
         * @param searchFactor An optional parameter specifying the search factor for the solver (default is the golden ratio).
         * @param maxiter The maximum number of iterations allowed before the search terminates (default is nxx::MAXITER).
         * @return An expected object containing either the bounds of the found root or an error.
         */
        template< typename SOLVER >
        requires requires(SOLVER solver, std::pair< typename SOLVER::return_type, typename SOLVER::return_type > bounds) {
                     // clang-format off
                     { solver.evaluate(0.0) } -> std::floating_point;
                     { solver.init(bounds) };
                     { solver.iterate() };
                     { solver.bounds() };
                     { solver.factor() };
                     // clang-format on
                 }
        inline auto search_impl(SOLVER                                                                  solver,
                                std::pair< typename SOLVER::return_type, typename SOLVER::return_type > bounds,
                                typename SOLVER::return_type searchFactor = std::numbers::phi_v< typename SOLVER::return_type >,
                                int                          maxiter      = nxx::MAXITER)
        {
            using ET = RootErrorImpl< decltype(bounds) >;
            using RT = tl::expected< decltype(bounds), ET >;

            solver.init(bounds, searchFactor);
            auto                         curBounds = solver.bounds();
            RT                           result    = curBounds;
            typename SOLVER::return_type eval_lower;
            typename SOLVER::return_type eval_upper;

            // Check for NaN or Inf.
            if (!std::isfinite(solver.evaluate(curBounds.first)) || !std::isfinite(solver.evaluate(curBounds.second))) {
                result = tl::make_unexpected(ET("Invalid initial brackets!", RootErrorType::NumericalError, result.value()));
                return result;
            }

            int iter = 1;
            while (true) {
                curBounds  = solver.bounds();
                eval_lower = solver.evaluate(curBounds.first);
                eval_upper = solver.evaluate(curBounds.second);

                // Check for NaN or Inf.
                if (!std::isfinite(curBounds.first) || !std::isfinite(curBounds.second) || !std::isfinite(eval_lower) ||
                    !std::isfinite(eval_upper))
                {
                    result = tl::make_unexpected(ET("Non-finite result!", RootErrorType::NumericalError, curBounds, iter));
                    break;
                }

                // Check if the root is bracketed by the bounds. If yes, return the bounds.
                if (eval_lower * eval_upper <= 0.0) {
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
    }    // namespace impl

    /**
     * @brief A search function template for root finding using a solver.
     *
     * This function template wraps the search_impl function with a more user-friendly interface.
     *
     * @tparam SOLVER The solver type to use for root finding.
     * @param solver The solver instance used to perform root finding.
     * @param bounds The initial bounds of the search interval as a pair.
     * @param searchFactor An optional parameter specifying the search factor for the solver (default is the golden ratio).
     * @param maxiter The maximum number of iterations allowed before the search terminates (default is nxx::MAXITER).
     * @return An expected object containing either the bounds of the found root or an error.
     */
    template< typename SOLVER >
    inline auto search(SOLVER                       solver,
                       auto                         bounds,
                       typename SOLVER::return_type searchFactor = std::numbers::phi_v< typename SOLVER::return_type >,
                       int                          maxiter      = nxx::MAXITER)
    {
        using RT      = typename SOLVER::return_type;
        auto [lo, hi] = bounds;
        return search_impl(solver, std::pair< RT, RT > { lo, hi }, searchFactor, maxiter);
    }

    /**
     * @brief A search function template for root finding using a solver with initializer_list bounds.
     *
     * This function template wraps the search_impl function with a more user-friendly interface that
     * accepts an initializer list for the bounds.
     *
     * @tparam SOLVER The solver type to use for root finding.
     * @tparam T The type of the bounds in the initializer list.
     * @param solver The solver instance used to perform root finding.
     * @param bounds The initial bounds of the search interval as an initializer_list.
     * @param searchFactor An optional parameter specifying the search factor for the solver (default is the golden ratio).
     * @param maxiter The maximum number of iterations allowed before the search terminates (default is nxx::MAXITER).
     * @return An expected object containing either the bounds of the found root or an error.
     * @throws std::logic_error if the initializer list contains a number of elements different from 2.
     */
    template< typename SOLVER, typename T >
    inline auto search(SOLVER                       solver,
                       std::initializer_list< T >   bounds,
                       typename SOLVER::return_type searchFactor = std::numbers::phi_v< typename SOLVER::return_type >,
                       int                          maxiter      = nxx::MAXITER)
    {
        using RT = typename SOLVER::return_type;
        if (bounds.size() != 2) throw std::logic_error("Initializer list must contain exactly two elements!");
        return search_impl(solver, std::pair< RT, RT > { *bounds.begin(), *(bounds.begin() + 1) }, searchFactor, maxiter);
    }

}    // namespace nxx::roots

#endif    // NUMERIXX_ROOTSEARCHING_HPP

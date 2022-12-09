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

#ifndef NUMERIX_ROOTBRACKETING_HPP
#define NUMERIX_ROOTBRACKETING_HPP

#include "calculus/Derivatives.hpp"

namespace numerix::roots
{

    // ========================================================================
    // BRACKET SEARCHING
    // ========================================================================

    /**
     * @brief
     * @param objective
     * @param lower
     * @param upper
     * @return
     */
    inline std::pair<double, double>
        bracket_search_up(const std::function<double(double)>& objective, double lower, double upper, double max_iter = 100)
    {
        if (upper <= lower) throw std::logic_error("Upper value must be higher than the lower value!");

        auto diff = upper - lower;
        for (int i = 0; i < max_iter; ++i) {
            if (objective(lower) * objective(upper) < 0.0) return std::make_pair(lower, upper);
            lower = upper;
            upper += diff;
        }

        throw std::logic_error("Bracket not found!");
    }

    // ========================================================================
    // ROOT-FINDING WITHOUT DERIVATIVES
    // ========================================================================

    /*
     * Forward declaration of the Ridders class.
     */
    template<typename FN>
        requires std::invocable<FN, double>
    class Ridders;

    /*
     * Forward declaration of the Bisection class.
     */
    template<typename FN>
        requires std::invocable<FN, double>
    class Bisection;

    /*
     * Private implementation details.
     */
    namespace impl
    {
        /*
         * Forward declaration of the BracketingTraits class.
         */
        template<typename FN>
        struct BracketingTraits;

        /*
         * Specialization of the BracketingTraits class for Ridders<FN>
         */
        template<typename FN>
        struct BracketingTraits<Ridders<FN>>
        {
            using function_type = FN;
        };

        /*
         * Specialization of the BracketingTraits class for Bisection<FN>
         */
        template<typename FN>
        struct BracketingTraits<Bisection<FN>>
        {
            using function_type = FN;
        };

        /**
         * @brief The BracketingBase class is a CRTP base class for bracketing solvers.
         * @tparam SOLVER The type of the (derived) solver class.
         */
        template<typename SOLVER>
            requires std::invocable<typename BracketingTraits<SOLVER>::function_type, double>
        class BracketingBase
        {
            /*
             * Friend declarations.
             */
            friend SOLVER;

        private:
            using function_type = typename BracketingTraits<SOLVER>::function_type;
            function_type m_func {}; /**< The function object to find the root for. */

            using RT = decltype(m_func(0.0));
            std::pair<RT, RT> m_bounds; /**< The bounds around the root. */

            /**
             * @brief Constructor, taking a function object as an argument.
             * @param objective The function object to find the root for.
             * @note Constructor is private to avoid direct usage by clients.
             */
            explicit BracketingBase(function_type objective) : m_func { objective } {}

        public:
            /**
             * @brief Initialize the solver by setting the initial bounds around the root.
             * @tparam T1 The type of the first bracket.
             * @tparam T2 The type of the second bracket.
             * @param bounds A std::pair with the bounds.
             */
            template<typename T1, typename T2>
                requires std::is_floating_point_v<T1> && std::is_floating_point_v<T2>
            void init(std::pair<T1, T2> bounds)
            {
                m_bounds = bounds;
            }

            /**
             * @brief Evaluate the function object at a given point.
             * @tparam T The type of the argument (must be floating point type)
             * @param value The value at which to evaluate the function.
             * @return The result of the evaluation.
             */
            template<typename T>
                requires std::is_floating_point_v<T>
            auto evaluate(T value)
            {
                return m_func(value);
            }

            /**
             * @brief Get the current bounds around the root.
             * @return A const reference to the bounds pair.
             */
            const auto& result() const { return m_bounds; }
        };
    }    // namespace impl

    /**
     * @brief The Ridders class is a derived class of the BracketingBase CRTP base class.
     * It implements Ridder's method for root finding without derivatives.
     * @tparam FN The type of the function object for which to find the root.
     */
    template<typename FN>
        requires std::invocable<FN, double>
    class Ridders final : public impl::BracketingBase<Ridders<FN>>
    {
        /*
         * Private alias declarations.
         */
        using Base = impl::BracketingBase<Ridders<FN>>;

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
        explicit Ridders(FN objective) : Base { objective } {}

        /**
         * @brief Perform one iteration. This is the main algorithm of Ridder's method.
         */
        void iterate()
        {
            using RT = decltype(Base::evaluate(Base::m_bounds.first));

            using std::abs;
            using std::pow;
            using std::sqrt;

            RT& x_lo = Base::m_bounds.first;
            RT& x_hi = Base::m_bounds.second;
            RT  f_lo = Base::evaluate(x_lo);
            RT  f_hi = Base::evaluate(x_hi);

            RT x_mid;
            RT f_mid;

            RT x_new;
            RT f_new;

            // ===== Calculate new bounds
            x_mid    = (x_lo + x_hi) / 2.0;
            f_mid    = Base::evaluate(x_mid);
            int sign = ((f_lo - f_hi) < 0.0 ? -1 : 1);
            x_new    = x_mid + (x_mid - x_lo) * ((sign * f_mid) / sqrt(f_mid * f_mid - f_lo * f_hi));
            f_new    = Base::evaluate(x_new);

            // ===== If x_new is NaN (i.e. the expression in the sqrt is negative), then return the input bounds.
            if (std::isnan(x_new)) Base::m_bounds = std::make_pair(x_lo, x_hi);

            // ===== General case: The root is between x_mid and x_new
            if (f_mid * f_new < 0.0) {
                if (x_mid < x_new)
                    Base::m_bounds = std::make_pair(x_mid, x_new);
                else
                    Base::m_bounds = std::make_pair(x_new, x_mid);
            }

            // ===== Degenerate cases: The root is between x_new and either x_lo or x_hi
            if (f_hi * f_new < 0.0) {
                if (x_hi < x_new)
                    Base::m_bounds = std::make_pair(x_hi, x_new);
                else
                    Base::m_bounds = std::make_pair(x_new, x_hi);
            }

            else {
                if (x_lo < x_new)
                    Base::m_bounds = std::make_pair(x_lo, x_new);
                else
                    Base::m_bounds = std::make_pair(x_new, x_lo);
            }
        }
    };

    /**
     * @brief The Bisection class is a derived class of the BracketingBase CRTP base class.
     * It implements the bisection method for root finding without derivatives.
     * @tparam FN The type of the function object for which to find the root.
     */
    template<typename FN>
        requires std::invocable<FN, double>
    class Bisection final : public impl::BracketingBase<Bisection<FN>>
    {
        /*
         * Private alias declarations.
         */
        using Base = impl::BracketingBase<Bisection<FN>>;

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
        explicit Bisection(FN objective) : Base { objective } {}

        /**
         * @brief Perform one iteration. This is the main algorithm of the bisection method.
         */
        void iterate()
        {
            using RT = decltype(Base::evaluate(Base::m_bounds.first));

            RT root = (Base::m_bounds.first + Base::m_bounds.second) / 2.0;

            if (Base::evaluate(Base::m_bounds.first) * Base::evaluate(root) < 0.0)
                Base::m_bounds = std::pair<RT, RT>(Base::m_bounds.first, root);
            else
                Base::m_bounds = std::pair<RT, RT>(root, Base::m_bounds.second);
        }
    };

    /**
     * @brief Convenience function for running a bracketing solver (i.e. without derivative), in one go.
     * @tparam S The type of the solver.
     * @param solver The solver object.
     * @param bounds The initial bounds around the root. A root must exist between the brackets.
     * @param eps The max. allowed error.
     * @param maxiter The max. number of allowed iterations.
     * @return The root estimate (mid-point between brackets).
     */
    template<typename S>
    inline auto fsolve(S solver, std::pair<double, double> bounds, double eps = 1.0E-6, int maxiter = 100)
    {
        using RT = decltype(solver.evaluate(0.0));

        solver.init(bounds);
        RT result;

        int iter = 1;
        while (true) {
            result = (solver.result().first + solver.result().second) / 2.0;
            if (abs(solver.result().first - solver.result().second) < eps || abs(solver.evaluate(result)) < eps) break;
            solver.iterate();

            ++iter;
            if (iter > maxiter) break;
        }

        return result;
    }

}    // namespace numerix::roots

#endif    // NUMERIX_ROOTBRACKETING_HPP

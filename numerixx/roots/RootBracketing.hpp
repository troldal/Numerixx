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


#include "RootError.hpp"
#include "calculus/Derivatives.hpp"

namespace nxx::roots
{

    // ========================================================================
    // BRACKET SEARCHING
    // ========================================================================

    /**
     * @brief HOLD
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
         * @brief The BracketingBase class serves as a base class for all bracketing solvers.
         * It functions as the CRTP base class for concrete bracketing solvers, and it therefore
         * does not contain any virtual functions.
         * @tparam POLICY The type of the (derived) solver class, e.g. Bisection or Ridders.
         */
        template<typename POLICY>
            requires std::invocable<typename BracketingTraits<POLICY>::function_type, double>
        class BracketingBase
        {
            /*
             * Friend declarations.
             */
            friend POLICY;

        private:
            using function_type = typename BracketingTraits<POLICY>::function_type;
            function_type m_func {}; /**< The function object to find the root for. */

            using return_type = decltype(m_func(0.0));
            std::pair<return_type, return_type> m_bounds {0.0, 0.0}; /**< Holds the current bounds around the root. */

            /**
             * @brief Constructor, taking a function object as an argument.
             * @param objective The function object to find the root for.
             * @note Constructor is private to avoid direct usage by clients.
             */
            explicit BracketingBase(function_type objective) : m_func { objective } {}

            /**
             * @brief
             * @param bounds
             */
            void setBounds(std::pair<return_type, return_type> bounds) { m_bounds = bounds; }

        public:
            /**
             * @brief The \ref init function initialises the solver by setting the initial bounds around the root.
             * At the point of initialisation, the solver has already been constructed and provided with the
             * function to solve. The purpose of the initialisation is only to set the initial bounds.
             * @param bounds A std::pair holding the initial bounds around the root. The root must be contained
             * inside these bounds.
             * @warning If the solver is used without initialising, the behaviour is undefined.
             */
            void init(std::pair< return_type, return_type > bounds) { setBounds(bounds); }

            /**
             * @brief The \ref evaluate function evaluates the function to solve, at a given point. This is
             * done simply by passing the given argument to the function object to solve.
             * @param value The value at which to evaluate the function.
             * @return The result of the evaluation. The return type will be the same as the return type of the
             * given function object.
             */
            auto evaluate(return_type value) { return m_func(value); }

            /**
             * @brief The result function returns the current bounds around the root. Every time an iteration
             * is executed, the bounds will narrow
             * @return A const reference to the current bounds. The bounds is returned as std::pair.
             * The value type of the bounds are the same as the return type of the function object.
             */
            const auto& bounds() const { return m_bounds; }

            /**
             * @brief
             */
            void iterate() { static_cast< POLICY& >(*this).iterate();}
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
            const auto& bounds = Base::bounds();
            using RT = decltype(Base::evaluate(bounds.first));

            using std::abs;
            using std::pow;
            using std::sqrt;

            const RT& x_lo = bounds.first;
            const RT& x_hi = bounds.second;
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
//            if (std::isnan(x_new)) Base::setBounds({x_lo, x_hi});

            // ===== General case: The root is between x_mid and x_new
            if (f_mid * f_new < 0.0) {
                if (x_mid < x_new)
                    Base::setBounds({x_mid, x_new});
                else
                    Base::setBounds({x_new, x_mid});
            }

            // ===== Degenerate cases: The root is between x_new and either x_lo or x_hi
            if (f_hi * f_new < 0.0) {
                if (x_hi < x_new)
                    Base::setBounds({x_hi, x_new});
                else
                    Base::setBounds({x_new, x_hi});
            }

            else {
                if (x_lo < x_new)
                    Base::setBounds({x_lo, x_new});
                else
                    Base::setBounds({x_new, x_lo});
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
            const auto& bounds = Base::bounds();
            using RT = decltype(Base::evaluate(bounds.first));

            RT root = (bounds.first + bounds.second) / 2.0;

            if (Base::evaluate(bounds.first) * Base::evaluate(root) < 0.0)
                Base::setBounds({bounds.first, root});
            else
                Base::setBounds({root, bounds.second});
        }
    };

    /**
     * @brief The fsolve function is a convenience function for running a bracketing solver (i.e. without derivative), without
     * dealing with low level details. If fine grained control is needed, such as advanced search stopping criteria or running each
     * iteration manually, please see the documentation for the solver classes.
     * @tparam SOLVER The type of the solver. This could be the Bisection or the Ridders solvers, but any solver with the correct interface can be used.
     * @param solver The actual solver object.
     * @param bounds The initial bounds around the root. A root must exist between the brackets.
     * @param eps The max. allowed error.
     * @param maxiter The max. number of allowed iterations.
     * @return The root estimate (mid-point between brackets).
     * @note The function returns only a single estimate of the root, i.e. the mid-point between the brackets after the final iteration.
     * If the actual brackets are needed, please use the solver object directly.
     */
    template<typename SOLVER >
    inline auto fsolve(SOLVER solver, std::pair<double, double> bounds, double eps = 1.0E-6, int maxiter = 100)
        -> tl::expected<double, error::RootError>
    {
        using RT = decltype(solver.evaluate(0.0));
        const auto& curBounds = solver.bounds();

        solver.init(bounds);
        RT result;

        int iter = 1;
        while (true) {
            if (!std::isfinite(curBounds.first) || !std::isfinite(curBounds.second)) return tl::make_unexpected(error::RootError("Root Error!"));
            result = (curBounds.first + curBounds.second) / 2.0;
            if (abs(curBounds.first - curBounds.second) < eps || abs(solver.evaluate(result)) < eps) break;
            solver.iterate();

            ++iter;
            if (iter > maxiter) break;
        }

        return result;
    }

}    // namespace numerix::roots

#endif    // NUMERIXX_ROOTBRACKETING_HPP

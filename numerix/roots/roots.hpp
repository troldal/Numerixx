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

#pragma clang diagnostic push
#pragma ide diagnostic ignored "HidingNonVirtualFunction"
#ifndef NUMERIX_ROOTS_HPP
#define NUMERIX_ROOTS_HPP

#include <array>
#include <cassert>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <algorithm>
#include <tuple>

#include "calculus/derivatives.hpp"

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

    /**
     * @brief
     * @tparam S
     * @tparam T1
     * @tparam T2
     * @tparam T3
     * @param solver
     * @param bounds
     * @param eps
     * @param maxiter
     * @return
     */
    template<typename S, typename T1, typename T2, typename T3 = double>
    requires std::is_floating_point_v<T1> && std::is_floating_point_v<T2> && std::is_floating_point_v<T3>
    inline auto fsolve(S solver, std::pair<T1, T2> bounds, T3 eps = 1.0E-6, int maxiter = 100) {

        using RT = decltype(solver.result().first);

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

    /**
     * @brief
     * @tparam Algo
     * @tparam Fn
     */
    template<typename Algo, typename Fn>
    requires std::is_invocable_v<Fn, double>
    class BracketingBase {

    private:

        Fn m_func {};

    protected:

        using RT = decltype(m_func(0.0));
        std::pair<RT, RT> m_bounds;

        /**
         * @brief
         * @param objective
         */
        explicit BracketingBase(Fn objective) : m_func{ objective }{}

    public:

        /**
         * @brief
         * @tparam T1
         * @tparam T2
         * @param bounds
         */
        template<typename T1, typename T2>
        requires std::is_floating_point_v<T1> && std::is_floating_point_v<T2>
        void init(std::pair<T1, T2> bounds) {
            m_bounds = bounds;
        }

        /**
         * @brief
         * @tparam T
         * @param value
         * @return
         */
        template<typename T>
        requires std::is_floating_point_v<T>
        auto evaluate(T value) {
            return m_func(value);
        }

        /**
         * @brief
         * @return
         */
        void iterate() {
            return static_cast<Algo*>(this)->iterate();
        }

        const auto& result() const {
            return m_bounds;
        }
    };

    /**
     * @brief
     */
    template<typename Fn>
    requires std::is_invocable_v<Fn, double>
    class ridders final : public BracketingBase<ridders<Fn>, Fn>
    {
        using Base = BracketingBase<ridders<Fn>, Fn>;

    public:

        /**
         * @brief
         * @param objective
         */
        explicit ridders(Fn objective) : Base(objective){}

        /**
         * @brief
         * @tparam T1
         * @tparam T2
         * @param bounds
         * @return
         */
        void iterate() {

            using RT = decltype(Base::evaluate(Base::m_bounds.first));

            using std::abs;
            using std::pow;
            using std::sqrt;

            RT& x_lo = Base::m_bounds.first;
            RT& x_hi = Base::m_bounds.second;
            RT f_lo = Base::evaluate(x_lo);
            RT f_hi = Base::evaluate(x_hi);

            RT x_mid;
            RT f_mid;

            RT x_new;
            RT f_new;

            // ===== Calculate new bounds
            x_mid = (x_lo + x_hi) / 2.0;
            f_mid = Base::evaluate(x_mid);
            int sign  = ((f_lo - f_hi) < 0.0 ? -1 : 1);
            x_new = x_mid + (x_mid - x_lo) * ((sign * f_mid) / sqrt(f_mid * f_mid - f_lo * f_hi));
            f_new = Base::evaluate(x_new);

            // ===== If x_new is NaN (i.e. the expression in the sqrt is negative), then return the input bounds.
            if (std::isnan(x_new))
                Base::m_bounds = std::make_pair(x_lo, x_hi);

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
     * @brief
     * @tparam Fn
     */
    template<typename Fn>
    requires std::is_invocable_v<Fn, double>
    class bisection final : public BracketingBase<bisection<Fn>, Fn>
    {
        using Base = BracketingBase<bisection<Fn>, Fn>;

    public:

        /**
         * @brief
         * @param objective
         */
        explicit bisection(Fn objective) : Base(objective){}

        /**
         * @brief
         * @tparam T1
         * @tparam T2
         * @param bounds
         * @return
         */
        //template<typename T1, typename T2>
        //requires std::is_floating_point_v<T1> && std::is_floating_point_v<T2>
        void iterate() {

            using RT = decltype(Base::evaluate(Base::m_bounds.first));

            RT root = (Base::m_bounds.first + Base::m_bounds.second) / 2.0;

            if (Base::evaluate(Base::m_bounds.first) * Base::evaluate(root) < 0.0)
                Base::m_bounds = std::pair<RT, RT>(Base::m_bounds.first, root);
            else
                Base::m_bounds = std::pair<RT, RT>(root, Base::m_bounds.second);
        }
    };



    // ========================================================================
    // ROOT-FINDING WITH DERIVATIVES
    // ========================================================================

    /**
     * @brief
     * @tparam S
     * @tparam T1
     * @tparam T2
     * @param solver
     * @param guess
     * @param eps
     * @param maxiter
     * @return
     */
    template<typename S, typename T1, typename T2>
    requires std::is_floating_point_v<T1> && std::is_floating_point_v<T2>
    inline auto fdfsolve(S solver, T1 guess, T2 eps = 1.0E-6, int maxiter = 100) {

        using RT = decltype(solver.result());

        solver.init(guess);

        int iter = 1;
        while (true) {

            if (abs(solver.evaluate(solver.result())) < eps) break;
            solver.iterate();

            ++iter;
            if (iter > maxiter) break;
        }

        return solver.result();
    }

    /**
     * @brief
     * @tparam Algo
     * @tparam Fn
     * @tparam DFn
     */
    template<typename Algo, typename Fn, typename DFn>
    requires std::is_invocable_v<Fn, float> && std::is_invocable_v<DFn, float>
    class PolishingBase {

    private:

        Fn m_func {};
        DFn m_deriv {};

    protected:

        using RT = decltype(m_func(0.0));
        RT m_guess;

        /**
         * @brief
         * @param objective
         * @param derivative
         */
        explicit PolishingBase(Fn objective, DFn derivative) : m_func{ objective },
                                                               m_deriv { derivative }
        {}

    public:

        /**
         * @brief
         * @tparam T
         * @param guess
         */
        template<typename T>
        requires std::is_floating_point_v<T>
        void init(T guess) {
            m_guess = guess;
        }

        /**
         * @brief
         * @tparam T
         * @param value
         * @return
         */
        template<typename T>
        requires std::is_floating_point_v<T>
        auto evaluate(T value) {
            return m_func(value);
        }

        /**
         * @brief
         * @tparam T
         * @param value
         * @return
         */
        template<typename T>
        requires std::is_floating_point_v<T>
        auto derivative(T value) {
            return m_deriv(value);
        }

        /**
         * @brief
         * @return
         */
        void iterate() {
            return static_cast<Algo*>(this)->iterate();
        }

        /**
         * @brief
         * @return
         */
        const auto& result() const {
            return m_guess;
        }
    };

    /**
     * @brief
     * @tparam Fn
     * @tparam DFn
     */
    template<typename Fn, typename DFn>
    requires std::is_invocable_v<Fn, float>
    class dnewton final : public PolishingBase<dnewton<Fn, DFn>, Fn, DFn>
    {
        using Base = PolishingBase<dnewton<Fn, DFn>, Fn, DFn>;

    public:

        /**
         * @brief
         * @param objective
         */
        explicit dnewton(Fn objective) : Base(objective, [=](double x){return numerix::deriv::central(objective, x);}){}

        /**
         * @brief
         */
        void iterate() {

            Base::m_guess = Base::m_guess - Base::evaluate(Base::m_guess) / Base::derivative(Base::m_guess);
        }
    };

    /**
     * @brief Deduction guide for the dnewton algorithm
     */
    template<typename Fn> dnewton(Fn objective) -> dnewton<Fn, std::function<decltype(objective(0.0))(decltype(objective(0.0)))>>;

    /**
     * @brief
     * @tparam Fn
     * @tparam DFn
     */
    template<typename Fn, typename DFn>
    requires std::is_invocable_v<Fn, float>
    class newton final : public PolishingBase<dnewton<Fn, DFn>, Fn, DFn>
    {
        using Base = PolishingBase<dnewton<Fn, DFn>, Fn, DFn>;

    public:

        /**
         * @brief
         * @param objective
         * @param deriv
         */
        explicit newton(Fn objective, DFn deriv) : Base(objective, deriv){}

        /**
         * @brief
         */
        void iterate() {
            Base::m_guess = Base::m_guess - Base::evaluate(Base::m_guess) / Base::derivative(Base::m_guess);
        }
    };


}    // namespace numerix

#endif    // NUMERIX_ROOTS_HPP

#pragma clang diagnostic pop
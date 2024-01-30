//
// Created by kenne on 27/01/2024.
//

#pragma once

#include "OptimCommon.hpp"

#include <Functions.hpp>

#include <algorithm>
#include <array>
#include <numbers>
#include <optional>
#include <tuple>
#include <utility>

namespace nxx::optim
{
    namespace detail
    {

        // =============================================================================================================
        //
        //   ,ad8888ba,                         88888888ba                                       88
        //  d8"'    `"8b                 ,d     88      "8b                                      88                      ,d
        // d8'        `8b                88     88      ,8P                                      88                      88
        // 88          88  8b,dPPYba,  MM88MMM  88aaaaaa8P'  8b,dPPYba,  ,adPPYYba,   ,adPPYba,  88   ,d8   ,adPPYba,  MM88MMM
        // 88          88  88P'    "8a   88     88""""""8b,  88P'   "Y8  ""     `Y8  a8"     ""  88 ,a8"   a8P_____88    88
        // Y8,        ,8P  88       d8   88     88      `8b  88          ,adPPPPP88  8b          8888[     8PP"""""""    88
        //  Y8a.    .a8P   88b,   ,a8"   88,    88      a8P  88          88,    ,88  "8a,   ,aa  88`"Yba,  "8b,   ,aa    88,
        //   `"Y8888Y"'    88`YbbdP"'    "Y888  88888888P"   88          `"8bbdP"Y8   `"Ybbd8"'  88   `Y8a  `"Ybbd8"'    "Y888
        //                 88
        //                 88
        //
        // =============================================================================================================

        template< typename DERIVED, IsFloatInvocable FUNCTION_T, IsFloat ARG_T, typename MODE_T >
        // requires...
        class OptimBracketBase
        {
            friend DERIVED;

        public:
            static constexpr bool IsBracketOptimizer = true;

            using RESULT_T = std::invoke_result_t< FUNCTION_T, ARG_T >;
            using BOUNDS_T = std::pair< ARG_T, ARG_T >;
            using RETURN_T = std::tuple< ARG_T, ARG_T, ARG_T >;

        protected:
            ~OptimBracketBase() = default;

        private:
            FUNCTION_T m_func {};
            BOUNDS_T   m_bounds {};
            MODE_T     m_mode;
            RETURN_T   m_result {};

        public:
            OptimBracketBase(FUNCTION_T objective, const IsFloatStruct auto& bounds)    // std::numbers::phi)
                : m_func { objective },
                  m_bounds { toPair(bounds) },
                  m_mode({}),
                  m_result { m_bounds.first, 0.0, m_bounds.second }
            {
                validateBounds(m_bounds);
            }

            template< size_t N >
            requires(N == 2)
            OptimBracketBase(FUNCTION_T objective, const ARG_T (&bounds)[N])    // std::numbers::phi)
                : m_func { objective },
                  m_bounds { std::pair { bounds[0], bounds[1] } },
                  m_mode({}),
                  m_result { m_bounds.first, 0.0, m_bounds.second }
            {
                validateBounds(m_bounds);
            }

            void setBounds(const BOUNDS_T& bounds)
            {
                m_bounds                = toPair(bounds);
                std::get< 0 >(m_result) = m_bounds.first;
                std::get< 2 >(m_result) = m_bounds.second;
                validateBounds(m_bounds);
            }

            void setBounds(const RETURN_T& range)
            {
                m_result = range;
                m_bounds = { std::get< 0 >(m_result), std::get< 2 >(m_result) };
                validateBounds(m_bounds);
            }

            // Default move/copy constructors/assignment operators:
            OptimBracketBase(const OptimBracketBase& other)     = default; /**< Default copy constructor. */
            OptimBracketBase(OptimBracketBase&& other) noexcept = default; /**< Default move constructor. */

            OptimBracketBase& operator=(const OptimBracketBase& other)     = default; /**< Default copy assignment operator. */
            OptimBracketBase& operator=(OptimBracketBase&& other) noexcept = default; /**< Default move assignment operator. */

            [[nodiscard]]
            RESULT_T evaluate(ARG_T value) const
            {
                if constexpr (std::same_as< MODE_T, Minimize >)
                    return m_func(value);
                else
                    return -m_func(value);
            }

            [[nodiscard]]
            const RETURN_T& current() const
            {
                return m_result;
            }

            void iterate() { std::invoke(static_cast< DERIVED& >(*this)); }
        };
    }    // namespace detail

    // =================================================================================================================
    //
    //   ,ad8888ba,                88           88                            ad88888ba
    //  d8"'    `"8b               88           88                           d8"     "8b
    // d8'                         88           88                           Y8,
    // 88              ,adPPYba,   88   ,adPPYb,88   ,adPPYba,  8b,dPPYba,   `Y8aaaaa,    8b,dPPYba,   ,adPPYba,
    // 88      88888  a8"     "8a  88  a8"    `Y88  a8P_____88  88P'   `"8a    `"""""8b,  88P'   "Y8  a8"     ""
    // Y8,        88  8b       d8  88  8b       88  8PP"""""""  88       88          `8b  88          8b
    //  Y8a.    .a88  "8a,   ,a8"  88  "8a,   ,d88  "8b,   ,aa  88       88  Y8a     a8P  88          "8a,   ,aa
    //   `"Y88888P"    `"YbbdP"'   88   `"8bbdP"Y8   `"Ybbd8"'  88       88   "Y88888P"   88           `"Ybbd8"'
    //
    // =================================================================================================================

    template< IsFloatInvocable FN, IsFloat ARG_T = double, typename MODE_T = Minimize >
    class GoldenSearch final : public detail::OptimBracketBase< GoldenSearch< FN, ARG_T, MODE_T >, FN, ARG_T, MODE_T >
    {
        using BASE =
            detail::OptimBracketBase< GoldenSearch< FN, ARG_T, MODE_T >, FN, ARG_T, MODE_T >; /**< Base class alias for readability. */
        using POINT_T = std::pair< ARG_T, ARG_T >;
        using RANGE_T = std::array< POINT_T, 4 >;

        std::optional< RANGE_T > m_range {};

    public:
        using BASE::BASE;

        void operator()()
        {
            if (!m_range) {
                initializeRange();
            }

            iterate_();
        }

    private:
        static constexpr size_t A  = 0;
        static constexpr size_t X1 = 1;
        static constexpr size_t X2 = 2;
        static constexpr size_t B  = 3;

        POINT_T calcPoint(ARG_T x) const { return { x, BASE::evaluate(x) }; }

        void initializeRange()
        {
            using std::numbers::phi;

            auto  bounds = BASE::current();
            ARG_T a      = std::get< 0 >(bounds);
            ARG_T d      = std::get< 2 >(bounds);
            ARG_T b      = d - (d - a) / phi;
            ARG_T c      = a + (d - a) / phi;

            m_range = RANGE_T { { calcPoint(a), calcPoint(b), calcPoint(c), calcPoint(d) } };
        }

        void iterate_()
        {
            using std::numbers::phi;

            auto& range = *m_range;

            // Comparing the function values at the intermediate points
            if (range[X1].second <= range[X2].second) {
                std::rotate(range.begin() + 1, range.begin() + 3, range.end());
                range[X1] = calcPoint(range[B].first - (range[B].first - range[A].first) / phi);
            }
            else {
                std::rotate(range.rbegin() + 1, range.rbegin() + 3, range.rend());
                range[X2] = calcPoint(range[A].first + (range[B].first - range[A].first) / phi);
            }

            BASE::setBounds({ std::make_tuple(range[A].first, range[X1].first, range[B].first) });
        }
    };

    // =================================================================================================================
    //
    //    88888888ba
    //    88      "8b                                        ,d
    //    88      ,8P                                        88
    //    88aaaaaa8P'  8b,dPPYba,   ,adPPYba,  8b,dPPYba,  MM88MMM
    //    88""""""8b,  88P'   "Y8  a8P_____88  88P'   `"8a   88
    //    88      `8b  88          8PP"""""""  88       88   88
    //    88      a8P  88          "8b,   ,aa  88       88   88,
    //    88888888P"   88           `"Ybbd8"'  88       88   "Y888
    //
    // =================================================================================================================

    template< IsFloatInvocable FN, IsFloat ARG_T = double, typename MODE_T = Minimize >
    class Brent final : public detail::OptimBracketBase< Brent< FN, ARG_T, MODE_T >, FN, ARG_T, MODE_T >
    {
        using BASE    = detail::OptimBracketBase< Brent< FN, ARG_T, MODE_T >, FN, ARG_T, MODE_T >;
        using POINT_T = std::pair< ARG_T, ARG_T >;
        using RANGE_T = std::array< POINT_T, 3 >;    // Three points for parabolic interpolation

        ARG_T tolerance = 1e-8;    // = static_cast<T>(ldexp(1.0, 1-bits));
        ARG_T x;                   // minima so far
        ARG_T w;                   // second best point
        ARG_T v;                   // previous value of w
        ARG_T u;                   // most recent evaluation point
        ARG_T delta;               // The distance moved in the last step
        ARG_T delta2;              // The distance moved in the step before last
        ARG_T fu, fv, fw, fx;      // function evaluations at u, v, w, x
        ARG_T mid;                 // midpoint of min and max
        ARG_T fract1, fract2;      // minimal relative movement in x

        static constexpr ARG_T golden = 0.3819660f;    // golden ratio, don't need too much precision here!

        bool m_isInitialised = false;

    public:
        using BASE::BASE;

        void operator()()
        {
            if (!m_isInitialised) {
                // x = w = v = BASE::current().second;
                x = w = v = std::get< 2 >(BASE::current());
                fw = fv = fx = BASE::evaluate(x);
                delta2 = delta  = 0;
                m_isInitialised = true;
            }

            iterate_();
        }

    private:
        POINT_T calcPoint(ARG_T x) const { return { x, BASE::evaluate(x) }; }

        void iterate_()
        {
            using std::abs;

            auto lower = std::get< 0 >(BASE::current());
            auto upper = std::get< 2 >(BASE::current());

            // The following implementation is taken from Boost.Math, *almost* verbatim.

            // get midpoint
            mid = (lower + upper) / 2;
            // // work out if we're done already:
            // fract1 = tolerance * tolerance * abs(x) + tolerance * tolerance / 4;
            // fract2 = 2 * fract1;
            fract1 = std::numeric_limits< ARG_T >::epsilon() * 2;
            fract2 = 2 * fract1;
            // if (abs(x - mid) <= (fract2 - (upper - lower) / 2)) return;

            if (abs(delta2) > fract1) {
                // try and construct a parabolic fit:
                ARG_T r = (x - w) * (fx - fv);
                ARG_T q = (x - v) * (fx - fw);
                ARG_T p = (x - v) * q - (x - w) * r;
                q       = 2 * (q - r);
                if (q > 0) p = -p;
                q        = abs(q);
                ARG_T td = delta2;
                delta2   = delta;
                // determine whether a parabolic step is acceptable or not:
                if ((abs(p) >= abs(q * td / 2)) || (p <= q * (lower - x)) || (p >= q * (upper - x))) {
                    // nope, try golden section instead
                    delta2 = (x >= mid) ? lower - x : upper - x;
                    delta  = golden * delta2;
                }
                else {
                    // whew, parabolic fit:
                    delta = p / q;
                    u     = x + delta;
                    if (((u - lower) < fract2) || ((upper - u) < fract2))
                        delta = (mid - x) < 0 ? static_cast< ARG_T >(-abs(fract1)) : static_cast< ARG_T >(abs(fract1));
                }
            }
            else {
                // golden section:
                delta2 = (x >= mid) ? lower - x : upper - x;
                delta  = golden * delta2;
            }
            // update current position:
            u  = (abs(delta) >= fract1) ? ARG_T(x + delta) : (delta > 0 ? ARG_T(x + abs(fract1)) : ARG_T(x - abs(fract1)));
            fu = BASE::evaluate(u);
            if (fu <= fx) {
                // good new point is an improvement!
                // update brackets:
                u >= x ? lower = x : upper = x;
                // update control points:
                v  = w;
                w  = x;
                x  = u;
                fv = fw;
                fw = fx;
                fx = fu;
            }
            else {
                // Oh dear, point u is worse than what we have already,
                // even so it *must* be better than one of our endpoints:
                u < x ? lower = u : upper = u;
                if ((fu <= fw) || (w == x)) {
                    // however it is at least second best:
                    v  = w;
                    w  = u;
                    fv = fw;
                    fw = fu;
                }
                else if ((fu <= fv) || (v == x) || (v == w)) {
                    // third best:
                    v  = u;
                    fv = fu;
                }
            }

            BASE::setBounds({ std::make_tuple(lower, x, upper) });
        }
    };

    template< typename FN, typename ARG_T, typename MODE_T = Minimize >
    requires IsFloatInvocable< FN > && IsFloat< ARG_T >
    Brent(FN, std::initializer_list< ARG_T >) -> Brent< FN, ARG_T, MODE_T >;

    template< typename FN, typename BOUNDS_T, typename MODE_T = Minimize >
    requires IsFloatInvocable< FN > && IsFloatStruct< BOUNDS_T >
    Brent(FN, const BOUNDS_T&) -> Brent< FN, StructCommonType_t< BOUNDS_T >, MODE_T >;

    // =================================================================================================================
    //
    //      ad88                                   88                      88
    //     d8"                              ,d     ""                      ""
    //     88                               88
    //   MM88MMM  ,adPPYba,   8b,dPPYba,  MM88MMM  88  88,dPYba,,adPYba,   88  888888888   ,adPPYba,
    //     88    a8"     "8a  88P'    "8a   88     88  88P'   "88"    "8a  88       a8P"  a8P_____88
    //     88    8b       d8  88       d8   88     88  88      88      88  88    ,d8P'    8PP"""""""
    //     88    "8a,   ,a8"  88b,   ,a8"   88,    88  88      88      88  88  ,d8"       "8b,   ,aa
    //     88     `"YbbdP"'   88`YbbdP"'    "Y888  88  88      88      88  88  888888888   `"Ybbd8"'
    //                        88
    //                        88
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

        template< typename SOLVER, typename TERMINATOR >
        requires SOLVER::IsBracketOptimizer
        auto foptimize_impl(const SOLVER& solver, const TERMINATOR& terminator)
        {
            SOLVER     _solver     = solver;
            TERMINATOR _terminator = terminator;

            const auto& [lower, x, upper] = _solver.current();
            size_t iter                   = 0;

            IterData< size_t, double > iterData { iter, lower, x, upper };

            while (true) {
                // std::cout << std::setprecision(10) << std::fixed;
                // std::cout << "Iteration " << iter << ": " << lower << " " << x << " " << upper << std::endl;
                // if (iter >= maxiter) break;

                // Termination logic from Boost.Math
                // auto mid = (lower + upper) / 2;
                // auto fract1 = eps * abs(x) + eps / 4;
                // auto fract2 = 2 * fract1;
                // if (abs(x - mid) <= (fract2 - (upper - lower) / 2)) break;

                // Termination logic from Numerical Recipes
                // if ((upper - lower) <= eps * x + eps / 2) break;

                iterData.iter  = iter;
                iterData.lower = lower;
                iterData.guess = x;
                iterData.upper = upper;

                if (_terminator(iterData)) break;
                _solver.iterate();
                ++iter;
            }

            return std::get< 1 >(_solver.current());
        }

        template< template< typename, typename, typename > class SOLVER_T,
                  typename MODE_T,
                  typename FN_T,
                  typename STRUCT_T,
                  typename... Args >
        auto foptimize_common(FN_T func, const STRUCT_T& bounds, const Args&... args)
        {
            using TUPLE_T = std::tuple< Args... >;
            using SOLVER = SOLVER_T< FN_T, StructCommonType_t< STRUCT_T >, MODE_T >;

            TUPLE_T args_tuple = std::make_tuple(args...);
            static_assert(std::tuple_size_v< TUPLE_T > <= 2, "Too many arguments passed to foptimize_common");

            // Zero arguments are passed...
            if constexpr (std::tuple_size_v< TUPLE_T > == 0)
                return detail::foptimize_impl(SOLVER(func, bounds), BracketTerminator {});


            // One argument is passed...
            else if constexpr (std::tuple_size_v< TUPLE_T > == 1) {
                using ArgType = std::tuple_element_t< 0, TUPLE_T >;

                if constexpr (std::is_floating_point_v< ArgType > || std::is_integral_v< ArgType >) {
                    return detail::foptimize_impl(SOLVER(func, bounds), BracketTerminator(std::get< 0 >(args_tuple)));
                }

                else if constexpr(std::is_same_v< std::invoke_result_t< ArgType, IterData< size_t, StructCommonType_t< STRUCT_T > > >, bool >) {
                    return detail::foptimize_impl(SOLVER(func, bounds), std::get< 0 >(args_tuple));
                }
                else
                    []<bool flag = false>(){static_assert(flag, "Invalid argument passed to foptimize_common");}();
            }

            // Two arguments are passed...
            else if constexpr (std::tuple_size_v< TUPLE_T > == 2) {
                // Unpack and use the two arguments
                using ArgType1 = std::tuple_element_t< 0, decltype(args_tuple) >;
                using ArgType2 = std::tuple_element_t< 1, decltype(args_tuple) >;
                static_assert(std::is_floating_point_v< ArgType1 > != std::is_floating_point_v< ArgType2 >,
                              "Two arguments must be one floating point and one integral type");
                return detail::foptimize_impl(SOLVER(func, bounds), BracketTerminator(std::get< 0 >(args_tuple), std::get< 1 >(args_tuple)));
            }
            else
                []<bool flag = false>(){static_assert(flag, "Invalid argument passed to foptimize_common");}();
        }

    }    // namespace detail

    // =============================================================//
    // ===== Specialization functions for general optimization =====//

    template< template< typename, typename, typename > class SOLVER_T,
              typename MODE_T,
              IsFloatOrComplexInvocable FN_T,
              IsFloatStruct             STRUCT_T,
              typename... ARGS >
    auto foptimize(FN_T func, STRUCT_T bounds, ARGS... args)
    {
        return detail::foptimize_common< SOLVER_T, MODE_T >(func, bounds, args...);
    }

    template< template< typename, typename, typename > class SOLVER_T,
              typename MODE_T,
              IsFloatOrComplexInvocable FN_T,
              IsFloat                   ARG_T,
              size_t                    N,
              typename... ARGS >
    requires(N == 2)
    auto foptimize(FN_T func, const ARG_T (&bounds)[N], ARGS... args)
    {
        return detail::foptimize_common< SOLVER_T, MODE_T >(func, std::pair { bounds[0], bounds[1] }, args...);
    }

    // =====================================================//
    // ===== Specialization functions for minimization =====//

    template< template< typename, typename, typename > class SOLVER_T,
              IsFloatOrComplexInvocable FN_T,
              IsFloatStruct             STRUCT_T,
              typename... ARGS >
    auto fminimize(FN_T func, const STRUCT_T& bounds, ARGS... args)
    {
        return foptimize< SOLVER_T, Minimize >(func, bounds, args...);
    }

    template< template< typename, typename, typename > class SOLVER_T,
              IsFloatOrComplexInvocable FN_T,
              IsFloat                   ARG_T,
              size_t                    N,
              typename... ARGS >
    requires(N == 2)
    auto fminimize(FN_T func, const ARG_T (&bounds)[N], ARGS... args)
    {
        return foptimize< SOLVER_T, Minimize >(func, std::pair { bounds[0], bounds[1] }, args...);
    }

    // =====================================================//
    // ===== Specialization functions for maximization =====//

    template< template< typename, typename, typename > class SOLVER_T,
              IsFloatOrComplexInvocable FN_T,
              IsFloatStruct             STRUCT_T,
              typename... ARGS >
    auto fmaximize(FN_T func, const STRUCT_T& bounds, ARGS... args)
    {
        return foptimize< SOLVER_T, Maximize >(func, bounds, args...);
    }

    template< template< typename, typename, typename > class SOLVER_T,
              IsFloatOrComplexInvocable FN_T,
              IsFloat                   ARG_T,
              size_t                    N,
              typename... ARGS >
    auto fmaximize(FN_T func, const ARG_T (&bounds)[N], ARGS... args)
    {
        return foptimize< SOLVER_T, Maximize >(func, std::pair { bounds[0], bounds[1] }, args...);
    }

}    // namespace nxx::optim
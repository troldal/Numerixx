//
// Created by kenne on 27/01/2024.
//

#pragma once

#include "OptimCommon.hpp"

#include <algorithm>
#include <array>
#include <numbers>
#include <optional>
#include <utility>

namespace nxx::optim
{

    template< typename DERIVED, IsFloatInvocable FUNCTION_T, IsFloat ARG_T, typename MODE_T >
    // requires...
    class OptimBracketBase
    {
        friend DERIVED;

    public:
        static constexpr bool IsSearchOptimizer = true;

        using RESULT_T = std::invoke_result_t< FUNCTION_T, ARG_T >;
        using BOUNDS_T = std::pair< ARG_T, ARG_T >;

    protected:
        ~OptimBracketBase() = default;

    private:
        FUNCTION_T m_func {};
        BOUNDS_T   m_bounds {};
        MODE_T     m_mode;

    public:
        OptimBracketBase(FUNCTION_T objective, IsFloatStruct auto bounds)    // std::numbers::phi)
            : m_func { objective },
              m_bounds { toPair(bounds) },
              m_mode({})
        {
            validateBounds(m_bounds);
        }

        template< size_t N >
        requires(N == 2)
        OptimBracketBase(FUNCTION_T objective, const ARG_T (&bounds)[N])    // std::numbers::phi)
            : m_func { objective },
              m_bounds { std::pair { bounds[0], bounds[1] } },
              m_mode({})
        {
            validateBounds(m_bounds);
        }

        void setBounds(const BOUNDS_T& bounds)
        {
            m_bounds = toPair(bounds);
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
            return m_func(value);
        }

        [[nodiscard]]
        const BOUNDS_T& current() const
        {
            return m_bounds;
        }

        void iterate() { std::invoke(static_cast< DERIVED& >(*this)); }
    };

    template< IsFloatInvocable FN, IsFloat ARG_T = double, typename MODE_T = Minimize >
    class GoldenSearch final : public OptimSearchBase< GoldenSearch< FN, ARG_T, MODE_T >, FN, ARG_T, MODE_T >
    {
        using BASE    = OptimSearchBase< GoldenSearch< FN, ARG_T, MODE_T >, FN, ARG_T, MODE_T >; /**< Base class alias for readability. */
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
            ARG_T a      = bounds.first;
            ARG_T d      = bounds.second;
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

            BASE::setBounds({ range[A].first, range[B].first });
        }
    };

    template< IsFloatInvocable FN, IsFloat ARG_T = double, typename MODE_T = Minimize >
    class Parabolic final : public OptimSearchBase< Parabolic< FN, ARG_T, MODE_T >, FN, ARG_T, MODE_T >
    {
        using BASE    = OptimSearchBase< Parabolic< FN, ARG_T, MODE_T >, FN, ARG_T, MODE_T >;
        using POINT_T = std::pair< ARG_T, ARG_T >;
        using RANGE_T = std::array< POINT_T, 3 >;    // Three points for parabolic interpolation

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
        POINT_T calcPoint(ARG_T x) const { return { x, BASE::evaluate(x) }; }

        void initializeRange()
        {
            using std::numbers::phi;
            auto  bounds = BASE::current();
            ARG_T a      = bounds.first;
            ARG_T b      = bounds.first + (bounds.second - bounds.first) / phi;    // Golden ratio
            ARG_T c      = bounds.second;

            m_range = RANGE_T { { calcPoint(a), calcPoint(b), calcPoint(c) } };
        }

        static constexpr size_t R = 0;
        static constexpr size_t S = 1;
        static constexpr size_t T = 2;

        void iterate_()
        {
            auto& range = *m_range;
            using std::numbers::phi;
            namespace rng = std::ranges;

            // Fit a parabola to the points and find the vertex
            ARG_T x_vertex = parabolicVertex(range[R], range[S], range[T]);

            x_vertex <= range[S].first ? range[T] = range[S] : range[R] = range[S];
            range[S] = calcPoint(x_vertex);
            rng::sort(range, [](const POINT_T& a, const POINT_T& b) { return a.first < b.first; });

            // Update bounds in the base class
            BASE::setBounds({ range[R].first, range[T].first });
        }

        ARG_T parabolicVertex(const POINT_T& p0, const POINT_T& p1, const POINT_T& p2) const
        {
            auto& [x0, f0] = p0;
            auto& [x1, f1] = p1;
            auto& [x2, f2] = p2;

            const ARG_T eps = std::sqrt(std::numeric_limits< ARG_T >::epsilon());

            const ARG_T quotient  = f0 * (x1 * x1 - x2 * x2) + f1 * (x2 * x2 - x0 * x0) + f2 * (x0 * x0 - x1 * x1);
            const ARG_T remainder = 2.0 * (f0 * (x1 - x2) + f1 * (x2 - x0) + f2 * (x0 - x1));
            // const ARG_T guessPoly = quotient / std::copysign(std::max(std::abs(remainder), eps), remainder);
            const ARG_T guessPoly = quotient / remainder;

            return guessPoly;
        }
    };

    template< IsFloatInvocable FN, IsFloat ARG_T = double, typename MODE_T = Minimize >
    class Brent final : public OptimSearchBase< Brent< FN, ARG_T, MODE_T >, FN, ARG_T, MODE_T >
    {
        using BASE    = OptimSearchBase< Brent< FN, ARG_T, MODE_T >, FN, ARG_T, MODE_T >;
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
                x = w = v = BASE::current().second;
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

            auto lower = BASE::current().first;
            auto upper = BASE::current().second;

            // The following implementation is taken from Boost.Math, *almost* verbatim.

            // get midpoint
            mid = (lower + upper) / 2;
            // work out if we're done already:
            fract1 = tolerance * abs(x) + tolerance / 4;
            fract2 = 2 * fract1;
            if (abs(x - mid) <= (fract2 - (upper - lower) / 2)) return;

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

            BASE::setBounds({ lower, upper });
        }
    };
}    // namespace nxx::optim
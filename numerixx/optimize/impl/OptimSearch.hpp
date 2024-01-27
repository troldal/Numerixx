//
// Created by kenne on 25/01/2024.
//

#pragma once

#include "OptimCommon.hpp"
#include "_external.hpp"

#include <Concepts.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <optional>

namespace nxx::optim
{



    template< typename DERIVED, IsFloatInvocable FUNCTION_T, IsFloat ARG_T, typename MODE_T >
    // requires...
    class OptimSearchBase
    {
        friend DERIVED;

    public:
        static constexpr bool IsSearchOptimizer = true;

        using RESULT_T = std::invoke_result_t< FUNCTION_T, ARG_T >;
        using BOUNDS_T = std::pair< ARG_T, ARG_T >;
        using RATIO_T  = ARG_T;

    protected:
        ~OptimSearchBase() = default;

    private:
        FUNCTION_T m_func {};
        BOUNDS_T   m_bounds {};
        RATIO_T    m_ratio {};
        MODE_T     m_mode;

    public:
        OptimSearchBase(FUNCTION_T objective, IsFloatStruct auto bounds, RATIO_T ratio = 1.0)    // std::numbers::phi)
            : m_func { objective },
              m_bounds { toPair(bounds) },
              m_ratio(ratio),
              m_mode({})
        {
            validateBounds(m_bounds);
        }

        template< size_t N >
        requires(N == 2)
        OptimSearchBase(FUNCTION_T objective, const ARG_T (&bounds)[N], RATIO_T ratio = 1.0)    // std::numbers::phi)
            : m_func { objective },
              m_bounds { std::pair { bounds[0], bounds[1] } },
              m_ratio(ratio),
              m_mode({})
        {
            validateBounds(m_bounds);
        }

        void setBounds(const BOUNDS_T& bounds)
        {
            m_bounds = toPair(bounds);
            validateBounds(m_bounds);
        }

        void setRatio(RATIO_T factor)
        {
            if (factor < 1.0) throw NumerixxError("Invalid factor.");
            m_ratio = factor;
        }

        // Default move/copy constructors/assignment operators:
        OptimSearchBase(const OptimSearchBase& other)     = default; /**< Default copy constructor. */
        OptimSearchBase(OptimSearchBase&& other) noexcept = default; /**< Default move constructor. */

        OptimSearchBase& operator=(const OptimSearchBase& other)     = default; /**< Default copy assignment operator. */
        OptimSearchBase& operator=(OptimSearchBase&& other) noexcept = default; /**< Default move assignment operator. */

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

        [[nodiscard]]
        RATIO_T ratio() const
        {
            return m_ratio;
        }

        void iterate() { std::invoke(static_cast< DERIVED& >(*this)); }
    };

    template< IsFloatInvocable FN, IsFloat ARG_T = double, typename MODE_T = Minimize >
    class AutoSearch final : public OptimSearchBase< AutoSearch< FN, ARG_T, MODE_T >, FN, ARG_T, MODE_T >
    {
        using BASE    = OptimSearchBase< AutoSearch< FN, ARG_T, MODE_T >, FN, ARG_T, MODE_T >; /**< Base class alias for readability. */
        using POINT_T = std::pair< ARG_T, ARG_T >;
        using RANGE_T = std::array< POINT_T, 3 >;

        std::optional< RANGE_T > m_range {};

        static constexpr size_t LOW  = 0;
        static constexpr size_t MID  = 1;
        static constexpr size_t HIGH = 2;

        POINT_T calcPoint(ARG_T x) const { return { x, BASE::evaluate(x) }; }

    public:
        using BASE::BASE;

        void operator()()
        {
            namespace rng = std::ranges;

            const ARG_T eps      = std::sqrt(std::numeric_limits< ARG_T >::epsilon());
            auto bounds = BASE::current();

            if (!m_range) {

                auto lowEval = BASE::evaluate(bounds.first);
                auto highEval = BASE::evaluate(bounds.second);
                if (lowEval < highEval) std::swap(bounds.first, bounds.second);

                m_range          = RANGE_T {};
                (*m_range)[LOW]  = calcPoint(bounds.first);
                (*m_range)[MID] = calcPoint(bounds.second);
                (*m_range)[HIGH]  = calcPoint(2 * bounds.second - bounds.first);
            }

            auto& range = *m_range;

            auto& x0 = range[LOW].first;
            auto& x1 = range[MID].first;
            auto& x2 = range[HIGH].first;
            auto& f0 = range[LOW].second;
            auto& f1 = range[MID].second;
            auto& f2 = range[HIGH].second;

            // If the mid point is the lowest, the minimum is already bracketed.
            const auto minimum = rng::min_element(range, [](const POINT_T& a, const POINT_T& b) { return a.second < b.second; });
            if (minimum == range.begin() + 1) return;

            const ARG_T maxStepSize = (range[MID].first - range[LOW].first) * BASE::ratio();
            const ARG_T guessStep = x2 + maxStepSize * 2;

            const ARG_T quotient = f0*(x1*x1 - x2*x2) + f1*(x2*x2 - x0*x0) + f2*(x0*x0 - x1*x1);
            const ARG_T remainder = 2.0 * (f0*(x1 - x2) + f1*(x2 - x0) + f2*(x0 - x1));
            const ARG_T guessPoly = quotient / std::copysign(std::max(std::abs(remainder), eps), remainder);

            const ARG_T guess = std::min(guessStep, guessPoly + eps);

            std::rotate(range.begin(), range.begin() + 1, range.end());
            range.back() = calcPoint(guess);
            BASE::setBounds({ range[LOW].first, range[MID].first });
        }

        // void operator()()
        // {
        //     // Constants used in the algorithm
        //
        //     // if (!m_stepSize) m_stepSize = BASE::m_bounds.second - BASE::m_bounds.first;
        //
        //     const ARG_T maxRatio = 1.5;
        //     const ARG_T eps      = std::sqrt(std::numeric_limits< ARG_T >::epsilon());    // Small number to avoid division by zero
        //
        //     ARG_T leftBound  = BASE::m_bounds.first;
        //     ARG_T rightBound = BASE::m_bounds.second;
        //     ARG_T newPoint   = 0.0;
        //     ARG_T trialPoint = 0.0;
        //
        //     ARG_T leftEval  = BASE::evaluate(leftBound);
        //     ARG_T rightEval = BASE::evaluate(rightBound);
        //     ARG_T newEval   = 0.0;
        //     ARG_T trialEval = 0.0;
        //
        //     // Ensure we are moving towards decreasing function values
        //     if (rightEval > leftEval) {
        //         std::swap(leftBound, rightBound);
        //         std::swap(leftEval, rightEval);
        //     }
        //
        //     // First expansion using the golden ratio
        //     newPoint = rightBound + BASE::m_ratio * (rightBound - leftBound);
        //     newEval  = BASE::evaluate(newPoint);
        //
        //     // Main bracketing logic
        //     if (rightEval > newEval || (rightEval < newEval && rightEval < leftEval)) {
        //
        //         ARG_T trialLimit = rightBound + maxRatio * (newPoint - rightBound);
        //
        //         trialPoint =
        //             (leftEval * (rightBound * rightBound - newPoint * newPoint) +
        //              rightEval * (newPoint * newPoint - leftBound * leftBound) +
        //              newEval * (leftBound * leftBound - rightBound * rightBound)) /
        //             (2.0 * (leftEval * (rightBound - newPoint) + rightEval * (newPoint - leftBound) + newEval * (leftBound -
        //             rightBound)));
        //
        //         // If the trial point is between the current right bound and newPoint...
        //         if ((rightBound - trialPoint) * (trialPoint - newPoint) > 0.0) {
        //             trialEval = BASE::evaluate(trialPoint);
        //
        //             // The condition trialEval < newEval checks if the function value at trialPoint
        //             // is less than the function value at newPoint. If true, this implies that
        //             // trialPoint is a better candidate for the minimum than newPoint. Consequently,
        //             // the algorithm updates the bounds to [rightBound, newPoint], narrowing down
        //             // the search interval around the region where the lower function value (trialEval)
        //             // was found. This step is essential in refining the search for the minimum,
        //             // ensuring that subsequent iterations focus on the most promising part of the
        //             // search space.
        //             if (trialEval < newEval) {
        //                 BASE::setBounds({ rightBound, newPoint });
        //                 return;
        //             }
        //
        //             // If trialEval is greater than rightEval, it indicates that the function's value
        //             // is increasing as we move from rightBound to trialPoint. Under the assumption
        //             // that the function is unimodal within the current interval, this increase
        //             // suggests that we have passed the minimum. Therefore, the minimum is likely
        //             // located in the interval [leftBound, trialPoint]. The algorithm thus updates
        //             // the bounds to focus the search within this narrowed range, where the minimum
        //             // is more probable.
        //             if (trialEval > rightEval) {
        //                 BASE::setBounds({ leftBound, trialPoint });
        //                 return;
        //             }
        //
        //             // Default Case: Reached when trialPoint is between rightBound and newPoint,
        //             // but does not indicate a better minimum candidate than newPoint nor suggests
        //             // the minimum lies between leftBound and trialPoint. In this scenario, the
        //             // algorithm opts to continue searching by further expanding the interval.
        //             // This is done by calculating a new trialPoint beyond newPoint, using
        //             // BASE::m_ratio (typically the golden ratio or a similar factor) to guide the
        //             // expansion. This step is crucial for exploring more of the search space
        //             // in pursuit of the minimum, especially when the current interval has not
        //             // yielded a conclusive direction for where the minimum might be.
        //             BASE::setBounds({ rightBound, newPoint });
        //             return;
        //         }
        //
        //         // This condition checks if the trialPoint lies between newPoint and trialLimit.
        //         // If true, it implies that the trialPoint is within a reasonable extension beyond
        //         // the current bounds, avoiding excessive expansion. The algorithm then evaluates
        //         // the function value at trialPoint (trialEval). If trialEval is found to be less
        //         // than the function value at newPoint (newEval), it indicates that a potentially
        //         // better minimum exists between newPoint and trialPoint. Consequently, the
        //         // algorithm updates the bounds to this narrower interval [newPoint, trialPoint],
        //         // focusing future search efforts in an area where the minimum is more likely to be found.
        //         if ((newPoint - trialPoint) * (trialPoint - trialLimit) > 0.0) {
        //             trialEval = BASE::evaluate(trialPoint);
        //             if (trialEval < newEval) {
        //                 BASE::setBounds({ newPoint, trialPoint });
        //                 return;
        //             }
        //
        //             // If trialEval is greater than newEval, it indicates that the function's value
        //             // is increasing as we move from newPoint to trialPoint. Under the assumption
        //             // that the function is unimodal within the current interval, this increase
        //             // suggests that we have passed the minimum. Therefore, the minimum is likely
        //             // located in the interval [rightBound, trialPoint]. The algorithm thus updates
        //             // the bounds to focus the search within this narrowed range, where the minimum
        //             // is more probable.
        //             BASE::setBounds({ rightBound, trialPoint });
        //             return;
        //         }
        //
        //         // This condition evaluates whether the trialPoint has reached or exceeded the
        //         // trialLimit relative to newPoint. Specifically, it checks if trialPoint is not
        //         // between newPoint and trialLimit. If true, this implies that the trialPoint is
        //         // either exactly at the trialLimit or has gone beyond it. In such cases, the
        //         // algorithm recognizes that further exploration in this direction might not be
        //         // productive or could lead to overshooting the minimum. Therefore, it sets the
        //         // new bounds to [newPoint, trialLimit], limiting the search to this range and
        //         // preventing further expansion beyond the trialLimit. This approach ensures that
        //         // the search remains focused and within a reasonable range, avoiding unnecessary
        //         // exploration of areas unlikely to contain the minimum.
        //         if ((trialPoint - trialLimit) * (trialLimit - newPoint) >= 0.0) {
        //             BASE::setBounds({ rightBound, newPoint });
        //             return;
        //         }
        //
        //         // In this default case, the algorithm has determined that the minimum has not
        //         // been effectively bracketed within the current bounds. This situation arises
        //         // when trialPoint does not yield a better minimum than newPoint and is not
        //         // indicative of the minimum being between leftBound and trialPoint.
        //         // Consequently, the algorithm decides to explore further by extending the
        //         // search interval beyond the current newPoint. The extension is calculated
        //         // using BASE::m_ratio (commonly the golden ratio or another factor) applied
        //         // to the distance between newPoint and rightBound. This approach allows the
        //         // algorithm to systematically and proportionally expand the search area,
        //         // moving trialPoint further in the direction where the minimum might be
        //         // located based on the current understanding of the function's behavior.
        //         BASE::setBounds({ rightBound, newPoint });
        //         return;
        //     }
        //
        //     return;
        // }
    };

}    // namespace nxx::optim
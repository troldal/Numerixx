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

namespace nxx::optim
{

    template< IsFloat ARG1, IsFloat ARG2 >
    void validateBounds(std::pair< ARG1, ARG2 >& bounds)
    {
        auto& [lower, upper] = bounds;
        if (lower == upper) throw NumerixxError("Invalid bounds.");
        if (lower > upper) std::swap(lower, upper);
    }

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
        OptimSearchBase(FUNCTION_T objective, IsFloatStruct auto bounds, RATIO_T ratio = std::numbers::phi)
            : m_func { objective },
              m_bounds { toPair(bounds) },
              m_ratio(ratio),
              m_mode({})
        {
            validateBounds(m_bounds);
        }

        template< size_t N >
        requires(N == 2)
        OptimSearchBase(FUNCTION_T objective, const ARG_T (&bounds)[N], RATIO_T ratio = std::numbers::phi)
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
        using BASE = OptimSearchBase< AutoSearch< FN, ARG_T, MODE_T >, FN, ARG_T, MODE_T >; /**< Base class alias for readability. */

    public:
        using BASE::BASE;

        void operator()()
        {
            // Constants used in the algorithm
            const ARG_T maxRatio = 100.0;
            const ARG_T tiny     = std::sqrt(std::numeric_limits< ARG_T >::epsilon());    // Small number to avoid division by zero

            std::array< ARG_T, 4 > points {};
            std::array< ARG_T, 4 > evals {};
            std::array< ARG_T, 4 > vars {};

            constexpr int LEFT  = 0;
            constexpr int RIGHT = 1;
            constexpr int NEW   = 2;
            constexpr int TRIAL = 3;

            points[LEFT]  = BASE::m_bounds.first;
            points[RIGHT] = BASE::m_bounds.second;
            points[NEW]   = 0.0;
            points[TRIAL] = 0.0;

            auto& leftBound  = points[LEFT];
            auto& rightBound = points[RIGHT];
            auto& newPoint   = points[NEW];
            auto& trialPoint = points[TRIAL];

            evals[LEFT]  = BASE::evaluate(leftBound);
            evals[RIGHT] = BASE::evaluate(rightBound);
            evals[NEW]   = 0.0;
            evals[TRIAL] = 0.0;

            auto& leftEval  = evals[LEFT];
            auto& rightEval = evals[RIGHT];
            auto& newEval   = evals[NEW];
            auto& trialEval = evals[TRIAL];

            // Ensure we are moving towards decreasing function values
            if (rightEval > leftEval) {
                std::swap(leftBound, rightBound);
                std::swap(leftEval, rightEval);
            }

            // First expansion using the golden ratio
            newPoint = rightBound + BASE::m_ratio * (rightBound - leftBound);
            newEval  = BASE::evaluate(newPoint);

            // Main bracketing logic
            if (rightEval > newEval) {
                // Parabolic fit calculations
                ARG_T r          = (rightBound - leftBound) * (rightEval - newEval);
                ARG_T q          = (rightBound - newPoint) * (rightEval - leftEval);
                ARG_T trialLimit = rightBound + maxRatio * (newPoint - rightBound);    // Limit for u based on GLIMIT
                trialPoint       = rightBound - ((rightBound - newPoint) * q - (rightBound - leftBound) * r) /
                                              (2.0 * std::max(std::abs(q - r), tiny) * ((q - r) > 0 ? 1 : -1));

                // If the trial point is between the current right bound and newPoint...
                if ((rightBound - trialPoint) * (trialPoint - newPoint) > 0.0) {
                    trialEval = BASE::evaluate(trialPoint);

                    // The condition trialEval < newEval checks if the function value at trialPoint
                    // is less than the function value at newPoint. If true, this implies that
                    // trialPoint is a better candidate for the minimum than newPoint. Consequently,
                    // the algorithm updates the bounds to [rightBound, newPoint], narrowing down
                    // the search interval around the region where the lower function value (trialEval)
                    // was found. This step is essential in refining the search for the minimum,
                    // ensuring that subsequent iterations focus on the most promising part of the
                    // search space.
                    if (trialEval < newEval) {
                        BASE::setBounds({ rightBound, newPoint });
                        return;
                    }

                    // If trialEval is greater than rightEval, it indicates that the function's value
                    // is increasing as we move from rightBound to trialPoint. Under the assumption
                    // that the function is unimodal within the current interval, this increase
                    // suggests that we have passed the minimum. Therefore, the minimum is likely
                    // located in the interval [leftBound, trialPoint]. The algorithm thus updates
                    // the bounds to focus the search within this narrowed range, where the minimum
                    // is more probable.
                    if (trialEval > rightEval) {
                        BASE::setBounds({ leftBound, trialPoint });
                        return;
                    }

                    // Default Case: Reached when trialPoint is between rightBound and newPoint,
                    // but does not indicate a better minimum candidate than newPoint nor suggests
                    // the minimum lies between leftBound and trialPoint. In this scenario, the
                    // algorithm opts to continue searching by further expanding the interval.
                    // This is done by calculating a new trialPoint beyond newPoint, using
                    // BASE::m_ratio (typically the golden ratio or a similar factor) to guide the
                    // expansion. This step is crucial for exploring more of the search space
                    // in pursuit of the minimum, especially when the current interval has not
                    // yielded a conclusive direction for where the minimum might be.
                    BASE::setBounds({ rightBound, newPoint });
                    return;
                }

                // This condition checks if the trialPoint lies between newPoint and trialLimit.
                // If true, it implies that the trialPoint is within a reasonable extension beyond
                // the current bounds, avoiding excessive expansion. The algorithm then evaluates
                // the function value at trialPoint (trialEval). If trialEval is found to be less
                // than the function value at newPoint (newEval), it indicates that a potentially
                // better minimum exists between newPoint and trialPoint. Consequently, the
                // algorithm updates the bounds to this narrower interval [newPoint, trialPoint],
                // focusing future search efforts in an area where the minimum is more likely to be found.
                else if ((newPoint - trialPoint) * (trialPoint - trialLimit) > 0.0) {
                    trialEval = BASE::evaluate(trialPoint);
                    if (trialEval < newEval) {
                        BASE::setBounds({ newPoint, trialPoint });
                        return;
                    }
                }

                // This condition evaluates whether the trialPoint has reached or exceeded the
                // trialLimit relative to newPoint. Specifically, it checks if trialPoint is not
                // between newPoint and trialLimit. If true, this implies that the trialPoint is
                // either exactly at the trialLimit or has gone beyond it. In such cases, the
                // algorithm recognizes that further exploration in this direction might not be
                // productive or could lead to overshooting the minimum. Therefore, it sets the
                // new bounds to [newPoint, trialLimit], limiting the search to this range and
                // preventing further expansion beyond the trialLimit. This approach ensures that
                // the search remains focused and within a reasonable range, avoiding unnecessary
                // exploration of areas unlikely to contain the minimum.
                else if ((trialPoint - trialLimit) * (trialLimit - newPoint) >= 0.0) {
                    BASE::setBounds({ newPoint, trialLimit });
                    return;
                }

                // In this default case, the algorithm has determined that the minimum has not
                // been effectively bracketed within the current bounds. This situation arises
                // when trialPoint does not yield a better minimum than newPoint and is not
                // indicative of the minimum being between leftBound and trialPoint.
                // Consequently, the algorithm decides to explore further by extending the
                // search interval beyond the current newPoint. The extension is calculated
                // using BASE::m_ratio (commonly the golden ratio or another factor) applied
                // to the distance between newPoint and rightBound. This approach allows the
                // algorithm to systematically and proportionally expand the search area,
                // moving trialPoint further in the direction where the minimum might be
                // located based on the current understanding of the function's behavior.
                else {
                    BASE::setBounds({ rightBound, newPoint });
                    return;
                }
            }

            return;
        }
    };

}    // namespace nxx::optim
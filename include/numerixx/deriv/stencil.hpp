// Finite-difference stencils as integer data (DESIGN §7.1): exact for every scalar type. A stencil<Order, Accuracy,
// Points> computes the Order-th derivative with error O(h^Accuracy) from sum(weight[i] * f(x + offset[i] h)) /
// (denominator h^Order).
//
// The spike provides the central stencils and the first-order one-sided pair; phase 2 adds forward/backward 1_2, 1_3,
// 2_1 and 2_2.
#pragma once

#include <numerixx/config.hpp>

#include <array>
#include <cstddef>
#include <cstdint>

NXX_BEGIN_HEADER

namespace nxx::deriv
{
    template<int Order, int Accuracy, std::size_t Points>
    struct stencil
    {
        static constexpr int order    = Order;
        static constexpr int accuracy = Accuracy;

        std::array<int, Points> offset;
        std::array<int, Points> weight;
        int                     denominator;

        // Points with a non-zero weight: the evaluations one derivative costs.
        constexpr std::uint32_t nonzero_points() const noexcept
        {
            std::uint32_t n = 0;
            for (const int w : weight)
                if (w != 0) ++n;
            return n;
        }
    };

    inline constexpr stencil<1, 2, 2> central_1_2 { { -1, 1 }, { -1, 1 }, 2 };
    inline constexpr stencil<1, 4, 4> central_1_4 { { -2, -1, 1, 2 }, { 1, -8, 8, -1 }, 12 };
    inline constexpr stencil<2, 2, 3> central_2_2 { { -1, 0, 1 }, { 1, -2, 1 }, 1 };
    inline constexpr stencil<2, 4, 5> central_2_4 { { -2, -1, 0, 1, 2 }, { -1, 16, -30, 16, -1 }, 12 };
    inline constexpr stencil<1, 1, 2> forward_1_1 { { 0, 1 }, { -1, 1 }, 1 };
    inline constexpr stencil<1, 1, 2> backward_1_1 { { -1, 0 }, { -1, 1 }, 1 };
}    // namespace nxx::deriv

NXX_END_HEADER

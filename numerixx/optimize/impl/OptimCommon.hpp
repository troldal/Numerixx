//
// Created by kenne on 25/01/2024.
//

#pragma once

namespace nxx::optim
{

    struct Minimize
    {
    };
    struct Maximize
    {
    };

    template< IsFloat ARG1, IsFloat ARG2 >
    void validateBounds(std::pair< ARG1, ARG2 >& bounds)
    {
        auto& [lower, upper] = bounds;
        if (lower == upper) throw NumerixxError("Invalid bounds.");
        if (lower > upper) std::swap(lower, upper);
    }

}    // namespace nxx::optim
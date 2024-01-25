//
// Created by kenne on 25/01/2024.
//

#pragma once

#include "Concepts.hpp"

namespace nxx
{

    auto toPair(const nxx::IsFloatStruct auto& data)
    {
        auto [val1, val2] = data;
        return std::pair { val1, val2 };
    }

    template< nxx::IsFloatStruct STRUCT_T, nxx::IsFloat VAL_1, nxx::IsFloat VAL_2 >
    STRUCT_T toStruct(const std::pair< VAL_1, VAL_2 >& data)
    {
        STRUCT_T result {};
        result.val1 = data.first;
        result.val2 = data.second;
        return result;
    }

}    // namespace nxx

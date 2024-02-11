//
// Created by kenne on 11/02/2024.
//

#pragma once

#include <concepts>

namespace nxx::roots {

    template<typename DERIVED, typename FUNCTION_T, typename ARG_T>
    requires std::same_as<typename BracketingTraits<DERIVED>::FUNCTION_T, FUNCTION_T>
             && nxx::IsFloatInvocable<FUNCTION_T> && nxx::IsFloat<ARG_T>
             && nxx::IsFloat<typename BracketingTraits<DERIVED>::RETURN_T>
    class BracketingBase;

}

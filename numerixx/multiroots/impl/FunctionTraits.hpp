//
// Created by kenne on 08/01/2024.
//

#pragma once

namespace nxx::traits
{
    /**
     * @brief Trait class to deduce the return type and argument type of a single-parameter callable.
     *
     * @details This primary template deduces the function traits of a callable object by
     *          decomposing the type of its call operator. It inherits from a specialization
     *          of FunctionTraits for the type of the callable's operator().
     *
     * @tparam FUNC_T The type of the callable object.
     */
    template< typename FUNC_T >
    struct FunctionTraits : public FunctionTraits< decltype(&FUNC_T::operator()) >
    {
    };

    /**
     * @brief Specialization of FunctionTraits for function pointers with a single parameter.
     *
     * @details This specialization deduces the return type and argument type of a single-parameter
     *          function pointer.
     *
     * @tparam RES_T The return type of the function.
     * @tparam PARAM_T The type of the single argument of the function.
     */
    template< typename RES_T, typename PARAM_T >
    struct FunctionTraits< RES_T (*)(PARAM_T) >
    {
        using return_type   = RES_T;
        using argument_type = PARAM_T;
    };

    /**
     * @brief Specialization of FunctionTraits for member function pointers (including lambdas)
     *        with a single parameter.
     *
     * @details This specialization handles member function pointers and captures the return type
     *          and argument type for member functions (or lambdas) that take a single parameter.
     *
     * @tparam RES_T The return type of the function.
     * @tparam CLASS_T The class type of the member function.
     * @tparam PARAM_T The type of the single argument of the function.
     */
    template< typename RES_T, typename CLASS_T, typename PARAM_T >
    struct FunctionTraits< RES_T (CLASS_T::*)(PARAM_T) const >
    {
        using return_type   = RES_T;
        using argument_type = PARAM_T;
    };
}


//
// Created by kenne on 08/01/2024.
//

#ifndef NUMERIXX_FUNCTIONTRAITS_HPP
#define NUMERIXX_FUNCTIONTRAITS_HPP

namespace nxx::traits
{
    /**
     * @brief Trait class to deduce the return type and argument type of a single-parameter callable.
     *
     * @details This primary template deduces the function traits of a callable object by
     *          decomposing the type of its call operator. It inherits from a specialization
     *          of FunctionTraits for the type of the callable's operator().
     *
     * @tparam Func The type of the callable object.
     */
    template< typename Func >
    struct FunctionTraits : public FunctionTraits< decltype(&Func::operator()) >
    {
    };

    /**
     * @brief Specialization of FunctionTraits for function pointers with a single parameter.
     *
     * @details This specialization deduces the return type and argument type of a single-parameter
     *          function pointer.
     *
     * @tparam R The return type of the function.
     * @tparam Arg The type of the single argument of the function.
     */
    template< typename R, typename Arg >
    struct FunctionTraits< R (*)(Arg) >
    {
        using return_type   = R;
        using argument_type = Arg;
    };

    /**
     * @brief Specialization of FunctionTraits for member function pointers (including lambdas)
     *        with a single parameter.
     *
     * @details This specialization handles member function pointers and captures the return type
     *          and argument type for member functions (or lambdas) that take a single parameter.
     *
     * @tparam R The return type of the function.
     * @tparam C The class type of the member function.
     * @tparam Arg The type of the single argument of the function.
     */
    template< typename R, typename C, typename Arg >
    struct FunctionTraits< R (C::*)(Arg) const >
    {
        using return_type   = R;
        using argument_type = Arg;
    };
}

#endif    // NUMERIXX_FUNCTIONTRAITS_HPP

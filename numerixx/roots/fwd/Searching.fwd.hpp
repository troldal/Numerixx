//
// Created by kenne on 11/02/2024.
//

#pragma once

namespace nxx::roots {

    /*
     * Forward declaration of the BracketSearchUp class.
     */
    template<IsFloatInvocable FN, IsFloat ARG_T>
    class BracketSearchUp;

    /*
     * Forward declaration of the BracketSearchDown class.
     */
    template<IsFloatInvocable FN, IsFloat ARG_T>
    class BracketSearchDown;

    /*
     * Forward declaration of the BracketExpandUp class.
     */
    template<IsFloatInvocable FN, IsFloat ARG_T>
    class BracketExpandUp;

    /*
     * Forward declaration of the BracketExpandDown class.
     */
    template<IsFloatInvocable FN, IsFloat ARG_T>
    class BracketExpandDown;

    /*
     * Forward declaration of the BracketExpandOut class.
     */
    template<IsFloatInvocable FN, IsFloat ARG_T>
    class BracketExpandOut;

    /*
     * Forward declaration of the BracketSubdivide class.
     */
    template<IsFloatInvocable FN, IsFloat ARG_T>
    class BracketSubdivide;

    /*
     * Forward declaration of the BracketIterData structure
     */
    template<std::integral ITER_T, IsFloat RESULT_T>
    using SearchIterData = std::tuple<ITER_T, RESULT_T, RESULT_T>;

    /*
     * Forward declaration of the search functions.
     */
    template<template<typename, typename> class SOLVER_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloatStruct STRUCT_T,
        typename... ARGS>
    auto search(FN_T func, STRUCT_T bounds, ARGS... args);

    template<template<typename, typename> class SOLVER_T,
        typename TOKEN_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloatStruct STRUCT_T>
    requires std::invocable<TOKEN_T, SearchIterData<size_t, StructCommonType_t<STRUCT_T>> &>
    auto search(FN_T func, STRUCT_T bounds);

    template<template<typename, typename> class SOLVER_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloat ARG_T,
        size_t N,
        typename... ARGS>
    requires(N == 2)
    auto search(FN_T func, const ARG_T (&bounds)[N], ARGS... args);

    template<template<typename, typename> class SOLVER_T,
        typename TOKEN_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloat ARG_T,
        size_t N>
    requires(N == 2) && std::invocable<TOKEN_T, SearchIterData<size_t, ARG_T> &>
    auto search(FN_T func, const ARG_T (&bounds)[N]);

} // namespace nxx::roots

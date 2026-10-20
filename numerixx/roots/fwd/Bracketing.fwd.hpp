//
// Created by kenne on 11/02/2024.
//

#pragma once

#include "Common.fwd.hpp"

#include <Concepts.hpp>

namespace nxx::roots {

    /*
     * Forward declaration of the Ridders class.
     */
    template<IsFloatInvocable FN, nxx::IsFloat ARG_T>
    class Ridder;

    /*
     * Forward declaration of the Bisection class.
     */
    template<IsFloatInvocable FN, nxx::IsFloat ARG_T>
    class Bisection;

    /*
     * Forward declaration of the RegulaFalsi class.
     */
    template<IsFloatInvocable FN, nxx::IsFloat ARG_T>
    class RegulaFalsi;

    /*
     * Forward declaration of the BracketIterData structure
     */
    template<std::integral ITER_T, IsFloat RESULT_T>
    using BracketIterData = IterData<ITER_T, RESULT_T, RESULT_T, RESULT_T>;

    /*
     * Forward declaration of the fsolve functions.
     */
    template<template<typename, typename> class SOLVER_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloatStruct STRUCT_T,
        typename... ARGS>
    auto fsolve(FN_T func, STRUCT_T bounds, ARGS... args);

    template<template<typename, typename> class SOLVER_T,
        typename TOKEN_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloatStruct STRUCT_T>
    requires std::invocable<TOKEN_T, BracketIterData<size_t, StructCommonType_t<STRUCT_T>> &>
    auto fsolve(FN_T func, STRUCT_T bounds);

    template<template<typename, typename> class SOLVER_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloat ARG_T,
        size_t N,
        typename... ARGS>
    requires(N == 2)
    auto fsolve(FN_T func, const ARG_T (&bounds)[N], ARGS &&...args);

    template<template<typename, typename> class SOLVER_T,
        typename TOKEN_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloat ARG_T,
        size_t N>
    requires(N == 2) && std::invocable<TOKEN_T, BracketIterData<size_t, ARG_T> &>
    auto fsolve(FN_T func, const ARG_T (&bounds)[N]);


} // namespace nxx::roots

//
// Created by kenne on 11/02/2024.
//

#pragma once

namespace nxx::roots {

    /*
     * Forward declaration of the Newton class.
     */
    template<IsFloatOrComplexInvocable FN, IsFloatOrComplexInvocable DFN, IsFloatOrComplex ARG_T>
    class Newton;

    /*
     * Forward declaration of the Secant class.
     */
    template<IsFloatOrComplexInvocable FN, IsFloatOrComplexInvocable DFN, IsFloatOrComplex ARG_T>
    class Secant;

    /*
     * Forward declaration of the Steffensen class.
     */
    template<IsFloatOrComplexInvocable FN, IsFloatOrComplexInvocable DFN, IsFloatOrComplex ARG_T>
    class Steffensen;

    /*
     * Forward declaration of the PolishingIterData structure
     */
    template<std::integral ITER_T, IsFloatOrComplex RESULT_T>
    using PolishingIterData = std::tuple<ITER_T, RESULT_T, std::vector<RESULT_T>>;

    /*
     * Forward declaration of the fdfsolve functions.
     */
    template<template<typename, typename, typename> class SOLVER_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloatOrComplexInvocable DERIV_T,
        IsFloatOrComplex GUESS_T,
        typename... ARGS>
    auto fdfsolve(FN_T func, DERIV_T derivative, GUESS_T guess, ARGS... args);

    template<template<typename, typename, typename> class SOLVER_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloatOrComplex GUESS_T,
        typename... ARGS>
    auto fdfsolve(FN_T func, GUESS_T guess, ARGS... args);

    template<template<typename, typename, typename> class SOLVER_T,
        typename TOKEN_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloatOrComplexInvocable DERIV_T,
        IsFloatOrComplex GUESS_T>
    requires std::invocable<TOKEN_T, PolishingIterData<size_t, GUESS_T> &>
    auto fdfsolve(FN_T func, DERIV_T derivative, GUESS_T guess);

    template<template<typename, typename, typename> class SOLVER_T,
        typename TOKEN_T,
        IsFloatOrComplexInvocable FN_T,
        IsFloatOrComplex GUESS_T>
    requires std::invocable<TOKEN_T, PolishingIterData<size_t, GUESS_T> &>
    auto fdfsolve(FN_T func, GUESS_T guess);

} // namespace nxx::roots

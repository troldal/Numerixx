/*
    888b      88  88        88  88b           d88  88888888888  88888888ba   88  8b        d8  8b        d8
    8888b     88  88        88  888b         d888  88           88      "8b  88   Y8,    ,8P    Y8,    ,8P
    88 `8b    88  88        88  88`8b       d8'88  88           88      ,8P  88    `8b  d8'      `8b  d8'
    88  `8b   88  88        88  88 `8b     d8' 88  88aaaaa      88aaaaaa8P'  88      Y88P          Y88P
    88   `8b  88  88        88  88  `8b   d8'  88  88"""""      88""""88'    88      d88b          d88b
    88    `8b 88  88        88  88   `8b d8'   88  88           88    `8b    88    ,8P  Y8,      ,8P  Y8,
    88     `8888  Y8a.    .a8P  88    `888'    88  88           88     `8b   88   d8'    `8b    d8'    `8b
    88      `888   `"Y8888Y"'   88     `8'     88  88888888888  88      `8b  88  8P        Y8  8P        Y8

    Copyright © 2022 Kenneth Troldal Balslev

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the “Software”), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is furnished
    to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#pragma once

#include <Concepts.hpp>
#include <Error.hpp>

// ===== Standard Library Includes
#include <stdexcept>

namespace nxx::roots {

    template<IsFloat ARG1, IsFloat ARG2>
    void validateBounds(std::pair<ARG1, ARG2> &bounds)
    {
        auto &[lower, upper] = bounds;
        if (lower == upper) throw NumerixxError("Invalid bounds.");
        if (lower > upper) std::swap(lower, upper);
    }

    // Helper struct to assist with deducing types from arguments
    template<typename... Args>
    struct ArgTypes;

    // Specialize for zero arguments
    template<>
    struct ArgTypes<>
    {
        using EPS_T = double;
        using ITER_T = size_t;
    };

    // Specialize for one argument (float or integral)
    template<typename Arg>
    struct ArgTypes<Arg>
    {
        using EPS_T = std::conditional_t<IsFloat<Arg>, Arg, double>;
        using ITER_T = std::conditional_t<std::integral<Arg>, Arg, size_t>;
    };

    // Specialize for two arguments
    template<typename Arg1, typename Arg2>
    struct ArgTypes<Arg1, Arg2>
    {
        using EPS_T = std::conditional_t<IsFloat<Arg1>, Arg1, Arg2>; // Assumes at least one is float
        using ITER_T = std::conditional_t<std::integral<Arg1>, Arg1, Arg2>; // Assumes at least one is integral
    };

    template<typename IMPL_T, typename... Args>
    class StopToken
    {
        using EPS_T = typename ArgTypes<Args...>::EPS_T;
        using ITER_T = typename ArgTypes<Args...>::ITER_T;

        IMPL_T m_impl{};
        EPS_T m_eps; /**< The epsilon value for the termination condition. */
        ITER_T m_maxiter; /**< The maximum iteration count for the termination condition. */

      public:
        explicit StopToken() : m_eps(epsilon<double>()), m_maxiter(iterations<double>()) {}

        constexpr explicit StopToken(Args... args)
        requires(sizeof...(Args) <= 2)
          : m_eps(std::tuple_size_v<decltype(std::make_tuple(args...))> == 2
                      ? std::get<EPS_T>(std::make_tuple(args...))
                      : (std::same_as<decltype(std::get<0>(std::make_tuple(args...))), EPS_T>
                                ? std::get<0>(std::make_tuple(args...))
                                : epsilon<double>())),
            m_maxiter(std::tuple_size_v<decltype(std::make_tuple(args...))> == 2
                          ? std::get<ITER_T>(std::make_tuple(args...))
                          : (std::same_as<decltype(std::get<0>(std::make_tuple(args...))), ITER_T>
                                    ? std::get<0>(std::make_tuple(args...))
                                    : iterations<double>()))
        {}

        /**
         * @brief Checks if the termination condition is met.
         * @param data The current iteration data.
         * @return true if the termination condition is met, false otherwise.
         */
        bool operator()(const auto &data) const { return m_impl(data, m_maxiter, m_eps); }

        EPS_T eps() const { return m_eps; }
        ITER_T maxiter() const { return m_maxiter; }
    };

    template<template<typename...> class TOKEN_T, typename... Args>
    auto makeToken(Args... args)
    {
        using TUPLE_T = std::tuple<Args...>;
        TUPLE_T args_tuple = std::make_tuple(args...);

        if constexpr (sizeof...(Args) == 1 && !IsFloat<decltype(std::get<0>(args_tuple))>
                      && !std::integral<decltype(std::get<0>(args_tuple))>)
            return std::get<0>(args_tuple);
        else
            return TOKEN_T<Args...>(args...);
    }

    template<typename ITERDATA_T, size_t IterIndex, size_t ResultIndex>
    class ResultProxy
    {
        ITERDATA_T
        m_iterData; /**< The IterData object holding the result of the root-finding problem. */

        using ITER_T = decltype(std::get<IterIndex>(m_iterData));
        using RESULT_T = decltype(std::get<ResultIndex>(m_iterData));

      public:
        /**
         * @brief Constructs the BracketSolverResult with an IterData object.
         * @param iterData The IterData object holding the result of the root-finding problem.
         */
        explicit ResultProxy(ITERDATA_T iterData) : m_iterData(iterData) {}

        ResultProxy(const ResultProxy &) = delete; // No copy constructor
        ResultProxy(ResultProxy &&) = delete; // No move constructor

        ResultProxy &operator=(const ResultProxy &) = delete; // No copy assignment
        ResultProxy &operator=(ResultProxy &&) = delete; // No move assignment

        /**
         * @brief Returns the result of the root-finding problem.
         * @tparam OUTPUT_T The type of the output. Defaults to the type of the result of the root-finding problem.
         * @return The result of the root-finding problem. If OUTPUT_T is a class, it is constructed with the
         * IterData object. Otherwise, the guess from the IterData object is returned.
         * @note This method is only available for rvalue references.
         */
        template<typename OUTPUT_T = RESULT_T>
        auto result() &&
        {
            if constexpr (std::is_class_v<OUTPUT_T>)
                return OUTPUT_T{}(m_iterData);
            else
                return std::get<ResultIndex>(m_iterData);
        }

        /**
         * @brief Returns the result of the root-finding problem using a specified outputter.
         * @tparam OUTPUTTER_T The type of the outputter. Must be a callable object that accepts an IterData object.
         * @param outputter The outputter to use to format the result.
         * @return The result of the root-finding problem, formatted by the outputter.
         * @note This method is only available for rvalue references.
         */
        template<typename OUTPUTTER_T>
        auto result(OUTPUTTER_T outputter) &&
        {
            return outputter(m_iterData);
        }
    };

} // namespace nxx::roots

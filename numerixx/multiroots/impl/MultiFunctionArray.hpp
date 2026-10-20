/*
    888b      88  88        88  88b           d88  88888888888  88888888ba   88  8b        d8  8b        d8
    8888b     88  88        88  888b         d888  88           88      "8b  88   Y8,    ,8P    Y8,    ,8P
    88 `8b    88  88        88  88`8b       d8'88  88           88      ,8P  88    `8b  d8'      `8b  d8'
    88  `8b   88  88        88  88 `8b     d8' 88  88aaaaa      88aaaaaa8P'  88      Y88P          Y88P
    88   `8b  88  88        88  88  `8b   d8'  88  88"""""      88""""88'    88      d88b          d88b
    88    `8b 88  88        88  88   `8b d8'   88  88           88    `8b    88    ,8P  Y8,      ,8P  Y8,
    88     `8888  Y8a.    .a8P  88    `888'    88  88           88     `8b   88   d8'    `8b    d8'    `8b
    88      `888   `"Y8888Y"'   88     `8'     88  88888888888  88      `8b  88  8P        Y8  8P        Y8

    Copyright © 2024 Kenneth Troldal Balslev

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

/**
 * @file MultiFunctionArray.hpp
 * @brief This file contains the declaration of the MultiFunctionArray class.
 *
 * The MultiFunctionArray class is a template class designed to store and manage a collection
 * of function objects, where each function follows a specified signature. It provides
 * functionalities to add functions to the array, invoke all functions with a set of parameters,
 * and access individual functions.
 */

#pragma once

#include <concepts>
#include <initializer_list>
#include <stdexcept>
#include <vector>

#include <Concepts.hpp>

namespace nxx::multiroots
{
    /**
     * @class MultiFunctionArray
     * @brief Template class for managing a collection of function objects.
     *
     * @tparam RET_T The return type of the functions in the array.
     * @tparam PARAM_T The parameter type of the functions in the array.
     */
    template< IsFloat RET_T, IsFloat PARAM_T >
    class MultiFunctionArray
    {
    public:
        using FUNC_T         = MultiFunction< RET_T, PARAM_T >;                   ///< Alias for function type.
        using const_iterator = typename std::vector< FUNC_T >::const_iterator;    ///< Iterator type.

        /**
         * @brief Default constructor for MultiFunctionArray.
         */
        MultiFunctionArray() = default;

        /**
         * @brief Constructor for initializing with a list of functions.
         * @details This constructor takes a variadic list of function objects and adds them to the array.
         * @tparam FIRSTFUNC_T Type of the first function in the list.
         * @tparam RESTFUNCS_T Types of the remaining functions in the list.
         * @param firstFunc The first function to be added to the array.
         * @param restFuncs The remaining functions to be added to the array.
         */
        template< typename FIRSTFUNC_T, typename... RESTFUNCS_T >
        requires IsInvocable< FIRSTFUNC_T >
        explicit MultiFunctionArray(FIRSTFUNC_T firstFunc, RESTFUNCS_T... restFuncs)
        {
            // Use addFunction for the first function
            addFunction(firstFunc);

            // Use addFunction for the rest of the functions
            (addFunction(restFuncs), ...);
        }

        /**
         * @brief Constructor for initializing with a container of functions.
         * @details This constructor takes a container with function objects and adds them to the array.
         * @tparam CONTAINER_T Type of the container.
         * @param container The container holding the functions to be added.
         */
        template< typename CONTAINER_T >
        requires IsInvocable< typename CONTAINER_T::value_type >
        explicit MultiFunctionArray(const CONTAINER_T& container)
        {
            // Iterate over the container and add each function to the array
            for (const auto& func : container) {
                if constexpr (std::is_same_v< MultiFunction< RET_T, PARAM_T >, typename CONTAINER_T::value_type >)
                    functions.push_back(func);
                else
                    addFunction(func);
            }
        }

        /**
         * @brief Adds a function to the array.
         * @details The function to be added must match the signature defined by RES_T and PARAM_T.
         * @tparam CALLABLE_T Type of the function to be added.
         * @param func The function to be added.
         */
        template< typename CALLABLE_T >
        requires std::same_as< RET_T, typename traits::FunctionTraits< std::decay_t< CALLABLE_T > >::return_type > &&
                 std::same_as< PARAM_T, typename traits::FunctionTraits< std::decay_t< CALLABLE_T > >::argument_type::value_type >
        void addFunction(const CALLABLE_T& func)
        {
            functions.push_back(FUNC_T(func));
        }

        /**
         * @brief Applies all functions to the input container and returns a container of the same type.
         * @tparam CONTAINER_T Type of the input container.
         * @param input The input container.
         * @return CONTAINER_T A container of the same type as input, containing the results of function applications.
         */
        template< typename CONTAINER_T >
        CONTAINER_T operator()(const CONTAINER_T& input) const
        {
            return evaluate< CONTAINER_T >(input);
        }

        /**
         * @brief Applies all functions to the initializer list and returns a vector of results.
         * @param input An initializer list of input values.
         * @return std::vector<RET_T> A vector containing the results of each function invocation.
         */
        std::vector< RET_T > operator()(std::initializer_list< RET_T > input) const { return evaluate< std::vector< RET_T > >(input); }

        /**
         * @brief Evaluates the input container using the stored functions, returning a container of the same type.
         * @tparam CONTAINER_T Type of the input container.
         * @param input The input container.
         * @return CONTAINER_T A container of the same type as input, containing the results of function applications.
         */
        template< typename CONTAINER_T >
        CONTAINER_T eval(const CONTAINER_T& input) const
        {
            return operator()< CONTAINER_T >(input);
        }

        /**
         * @brief Evaluates the initializer list using the stored functions, returning a vector of results.
         * @param input An initializer list of input values.
         * @return std::vector<RET_T> A vector containing the results of each function invocation.
         */
        std::vector< RET_T > eval(std::initializer_list< RET_T > input) const { return operator()< std::vector< RET_T > >(input); }

        /**
         * @brief Applies all functions to the input container and returns a container of the specified output type.
         * @tparam OUT_T Template template parameter specifying the type of the output container.
         * @tparam CONTAINER_T Type of the input container.
         * @param input The input container.
         * @return OUT_T<RET_T> An output container of the specified type containing the results.
         */
        template< template< typename... > class OUT_T, typename CONTAINER_T >
        OUT_T< RET_T > operator()(const CONTAINER_T& input) const
        {
            return evaluate< OUT_T< RET_T > >(input);
        }

        /**
         * @brief Applies all functions to the initializer list and returns a container of the specified output type.
         * @tparam OUT_T Template template parameter specifying the type of the output container.
         * @param input An initializer list of input values.
         * @return OUT_T<RET_T> An output container of the specified type containing the results.
         */
        template< template< typename... > class OUT_T >
        OUT_T< RET_T > operator()(std::initializer_list< RET_T > input) const
        {
            return evaluate< OUT_T< RET_T > >(input);
        }

        /**
         * @brief Evaluates the input container using the stored functions, returning a container of the specified output type.
         * @tparam OUT_T Template template parameter specifying the type of the output container.
         * @tparam CONTAINER_T Type of the input container.
         * @param input The input container.
         * @return OUT_T<RET_T> An output container of the specified type containing the results.
         */
        template< template< typename... > class OUT_T, typename CONTAINER_T >
        OUT_T< RET_T > eval(const CONTAINER_T& input) const
        {
            return operator()< OUT_T >(input);
        }

        /**
         * @brief Evaluates the initializer list using the stored functions, returning a container of the specified output type.
         * @tparam OUT_T Template template parameter specifying the type of the output container.
         * @param input An initializer list of input values.
         * @return OUT_T<RET_T> An output container of the specified type containing the results.
         */
        template< template< typename... > class OUT_T >
        OUT_T< RET_T > eval(std::initializer_list< RET_T > input) const
        {
            return operator()< OUT_T >(input);
        }

        /**
         * @brief Applies all functions to the input container and returns an output container with additional template parameters.
         * @tparam OUT_T Template template parameter specifying the type of the output container, with additional template parameters.
         * @tparam CONTAINER_T Type of the input container.
         * @param input The input container.
         * @return OUT_T<RET_T, false> An output container of the specified type containing the results.
         */
        template< template< typename, bool, typename... > class OUT_T, typename CONTAINER_T >
        OUT_T< RET_T, false > operator()(const CONTAINER_T& input) const
        {
            return evaluate< OUT_T< RET_T, false > >(input);
        }

        /**
         * @brief Applies all functions to the initializer list and returns an output container with additional template parameters.
         * @tparam OUT_T Template template parameter specifying the type of the output container, with additional template parameters.
         * @param input An initializer list of input values.
         * @return OUT_T<RET_T, false> An output container of the specified type containing the results.
         */
        template< template< typename, bool, typename... > class OUT_T >
        OUT_T< RET_T, false > operator()(std::initializer_list< RET_T > input) const
        {
            return evaluate< OUT_T< RET_T, false > >(input);
        }

        /**
         * @brief Evaluates the input container using the stored functions, returning an output container with additional template
         * parameters.
         * @tparam OUT_T Template template parameter specifying the type of the output container, with additional template parameters.
         * @tparam CONTAINER_T Type of the input container.
         * @param input The input container.
         * @return OUT_T<RET_T, false> An output container of the specified type containing the results.
         */
        template< template< typename, bool, typename... > class OUT_T, typename CONTAINER_T >
        OUT_T< RET_T, false > eval(const CONTAINER_T& input) const
        {
            return operator()< OUT_T >(input);
        }

        /**
         * @brief Evaluates the initializer list using the stored functions, returning an output container with additional template
         * parameters.
         * @tparam OUT_T Template template parameter specifying the type of the output container, with additional template parameters.
         * @param input An initializer list of input values.
         * @return OUT_T<RET_T, false> An output container of the specified type containing the results.
         */
        template< template< typename, bool, typename... > class OUT_T >
        OUT_T< RET_T, false > eval(std::initializer_list< RET_T > input) const
        {
            return operator()< OUT_T >(input);
        }

        /**
         * @brief Accesses a function at a specified index.
         * @param index The index of the function to access.
         * @return const FUNC_T& A reference to the function at the specified index.
         */
        const FUNC_T& operator[](size_t index) const
        {
            if (index >= functions.size()) {
                throw std::out_of_range("Index out of range");
            }
            return functions[index];
        }

        /**
         * @brief Returns an iterator to the beginning of the function array.
         * @return const_iterator An iterator to the first function in the array.
         */
        const_iterator begin() const { return functions.cbegin(); }
        const_iterator cbegin() const { return functions.cbegin(); }

        const_iterator end() const { return functions.cend(); }
        const_iterator cend() const { return functions.cend(); }

        /**
         * @brief Returns the number of functions in the array.
         * @return auto The number of functions stored in the array.
         */
        auto size() const { return functions.size(); }

    private:
        std::vector< FUNC_T > functions;    ///< Internal storage for function objects.

        /**
         * @brief Evaluates the input container using the stored functions and returns an output container of the specified type.
         * @tparam OUT_T Type of     the output container.
         * @tparam CONT_T Type of the input container.
         * @param input The input container.
         * @return OUT_T An output container of the specified type containing the results.
         * @details This method transforms each function's result into an output container of type OUT_T.
         *          It requires that OUT_T is constructible with the size of the input and provides an iterator interface.
         */
        template< typename OUT_T, typename CONT_T >
        OUT_T evaluate(const CONT_T& input) const
        {
            OUT_T                result(functions.size());
            std::vector< RET_T > args(input.begin(), input.end());

            std::transform(functions.begin(), functions.end(), result.begin(), [&args](const auto& func) {
                return std::invoke(func, args);
            });

            return result;
        }
    };

    /**
     *@brief Deduction guide for MultiFunctionArray to deduce template parameters from function types.
     */
    template< typename FirstFunc, typename... RestFuncs >
    MultiFunctionArray(FirstFunc&& firstFunc, RestFuncs&&... restFuncs)
        -> MultiFunctionArray< typename traits::FunctionTraits< std::decay_t< FirstFunc > >::return_type,
                               typename traits::FunctionTraits< std::decay_t< FirstFunc > >::argument_type::value_type >;

}    // namespace nxx::multiroots


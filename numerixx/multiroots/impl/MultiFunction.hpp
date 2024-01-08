//
// Created by kenne on 08/01/2024.
//

#ifndef NUMERIXX_MULTIFUNCTION_HPP
#define NUMERIXX_MULTIFUNCTION_HPP

#include "ContainerTraits.hpp"
#include "FunctionTraits.hpp"
#include <Concepts.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <span>
#include <type_traits>

namespace nxx::multiroots
{

    /**
     * @brief A class template to encapsulate a multi-dimensional function for root finding.
     *
     * @details This class wraps a callable object, allowing it to be invoked with different types
     *          of inputs, such as initializer lists and containers. The class ensures that the
     *          types of the callable's return value and arguments meet certain requirements.
     *
     * @tparam T The return type of the function. Must be a floating-point type as per nxx::IsFloat.
     * @tparam U The type of the elements in the input container. Must be a floating-point type after
     *           removing const/volatile qualifiers and references.
     */
    template< typename T, typename U >
    requires nxx::IsFloat< T > && nxx::IsFloat< std::remove_cvref_t< U > >
    class MultiFunction
    {
    public:
        using FUNC_T = std::function< T(std::span< U >) >;    ///< Type of the encapsulated function.

        /**
         * @brief Constructs a MultiFunction object with a given callable.
         *
         * @details This constructor wraps a given callable, ensuring it has a compatible signature.
         *
         * @tparam CALLABLE_T The type of the callable to be wrapped.
         * @param  f The callable object to be wrapped.
         */
        template< typename CALLABLE_T >
        requires(
            std::is_same_v< typename traits::FunctionTraits< CALLABLE_T >::argument_type, std::span< std::remove_cvref_t< U > > > ||
            std::is_same_v< typename traits::FunctionTraits< CALLABLE_T >::argument_type, std::span< const std::remove_cvref_t< U > > >)
        explicit MultiFunction(CALLABLE_T f)
            : function(f)
        {}

        /**
         * @brief Invokes the function using a std::initializer_list, specifically for const U types.
         *
         * @details This operator allows the function to be called with an initializer list.
         *          It is enabled only when U is a const type. The initializer list is converted
         *          to a std::span before invoking the encapsulated function.
         *
         * @param  list An initializer list of elements of type U.
         * @return The return value of the encapsulated function.
         */
        T operator()(std::initializer_list< U > list) const requires std::is_const_v< U >
        {
            std::span< U > span(data(list), list.size());
            return function(span);
        }

        /**
         * @brief Invokes the function using a std::initializer_list, specifically for non-const U types.
         *
         * @details This operator allows the function to be called with an initializer list.
         *          It is enabled only when U is a non-const type. The initializer list is first
         *          converted to a std::vector, then to a std::span, before invoking the encapsulated function.
         *
         * @param  list An initializer list of elements of type U.
         * @return The return value of the encapsulated function.
         */
        T operator()(std::initializer_list< U > list) const requires(!std::is_const_v< U >)
        {
            std::vector< U > tempVec(list);    // Create a temporary vector from the initializer list
            std::span< U >   span(tempVec.data(), tempVec.size());
            return function(span);
        }

        /**
         * @brief Invokes the function using a container, where the container's value type matches U.
         *
         * @details This operator allows the function to be called with a container,
         *          provided the container's value type exactly matches U (after removing cv-ref qualifiers).
         *          It converts the container to a std::span before invoking the function.
         *          If U is a const type, the container is used directly. Otherwise, a copy is made.
         *
         * @tparam CONTAINER_T The type of the input container.
         * @param  container The container holding elements of type U.
         * @return The return value of the encapsulated function.
         */
        template< typename CONTAINER_T >
        requires std::is_same_v< std::remove_cvref_t< U >, traits::ContainerValueType_t< CONTAINER_T > >
        T operator()(const CONTAINER_T& container) const
        {
            if constexpr (std::is_const_v< U >) {
                std::span< U > span(container.data(), container.size());
                return function(span);
            }
            else {
                CONTAINER_T    tempVec(container);    // Create a temporary vector from the container
                std::span< U > span(tempVec.data(), tempVec.size());
                return function(span);
            }
        }

        /**
         * @brief Invokes the function using a container, where the container's value type is convertible to U.
         *
         * @details This operator allows the function to be called with a container,
         *          provided the container's value type is convertible to U, and not exactly U.
         *          It creates a temporary std::vector of U, copies and converts the elements
         *          of the container into it, and then uses this vector to invoke the function.
         *
         * @tparam CONTAINER_T The type of the input container.
         * @param  container The container holding elements convertible to U.
         * @return The return value of the encapsulated function.
         */
        template< typename CONTAINER_T >
        requires std::convertible_to< traits::ContainerValueType_t< CONTAINER_T >, std::remove_cvref_t< U > > &&
                 (!std::is_same_v< std::remove_cvref_t< U >, traits::ContainerValueType_t< CONTAINER_T > >)
        T operator()(const CONTAINER_T& container) const
        {
            std::vector< U > tempVec(container.size());    // Create a temporary vector from the container
            std::transform(container.cbegin(), container.cend(), tempVec.begin(), [](float value) { return static_cast< double >(value); });

            std::span< U > span(tempVec.data(), tempVec.size());
            return function(span);
        }

    private:
        FUNC_T function;    ///< The internal function object.
    };

    /**
     * @brief Deduction guide for MultiFunction.
     *
     * @details Provides a deduction guide for the MultiFunction template, allowing the compiler
     *          to deduce the template arguments from the type of the function object passed to
     *          the constructor.
     *
     * @tparam Func The type of the function object.
     */
    template< typename FUNC_T >
    MultiFunction(FUNC_T) -> MultiFunction< typename traits::FunctionTraits< FUNC_T >::return_type,
                                            typename traits::FunctionTraits< FUNC_T >::argument_type::value_type >;

}    // namespace nxx::multiroots

#endif    // NUMERIXX_MULTIFUNCTION_HPP

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

// ===== External Includes
#include "_external.hpp"

// ===== Standard Library Includes
#include <complex>
#include <span>

namespace nxx
{
    template< typename T >
    concept IsFloat = (!std::integral<T>)&&(
        std::floating_point<T> || boost::multiprecision::is_number<T>::value || std::convertible_to<T, double>);

    //    template<typename T, typename F>
    //    concept IsFloatFunction = requires(F f, const std::vector<T>& v) {
    //                                        { f(v) } -> std::convertible_to<T>;
    //                                    };

    /**
     * @brief Concept to check if a type is a complex number.
     * @tparam T The type to check.
     */
    template< typename T >
    concept IsComplex = requires(T x) {
                            typename T::value_type;
                            {
                                std::real(x)
                            } -> nxx::IsFloat;
                            {
                                std::imag(x)
                            } -> nxx::IsFloat;
                        };

    template< typename T >
    concept IsFloatOrComplex = IsComplex< T > || nxx::IsFloat< T >;

    namespace detail
    {
        template< typename T, typename = void >
        struct is_callable : std::false_type
        {
        };

        template< typename T >
        struct is_callable< T, std::void_t< decltype(&T::operator()) > > : std::true_type
        {
        };

        // Helper variable template
        template< typename T >
        inline constexpr bool is_callable_v = is_callable< T >::value;
    }    // namespace detail

    template< typename Func >
    concept IsInvocable = detail::is_callable_v< Func >;

    /**
     * @brief Concept checking whether a type is a callable function object that returns a floating point type.
     * @tparam Func The type to check.
     */
    template< typename Func >
    concept IsFloatInvocable = requires(Func f) {
                                   requires nxx::IsFloatOrComplex< std::invoke_result_t< Func, float > >;
                                   requires nxx::IsFloatOrComplex< std::invoke_result_t< Func, double > >;
                                   requires nxx::IsFloatOrComplex< std::invoke_result_t< Func, long double > >;
                               };

    template< typename Func >
    concept IsComplexInvocable = requires(Func f) {
                                     requires nxx::IsComplex< std::invoke_result_t< Func, std::complex< float > > >;
                                     requires nxx::IsComplex< std::invoke_result_t< Func, std::complex< double > > >;
                                     requires nxx::IsComplex< std::invoke_result_t< Func, std::complex< long double > > >;
                                 };

    template< typename Func >
    concept IsFloatOrComplexInvocable = IsFloatInvocable< Func > || IsComplexInvocable< Func >;

    template<typename Func, typename Arg>
    struct is_invocable_with_span {
        static constexpr bool value = requires(Func f, std::span<Arg> arg) {
                                          { f(arg) } -> nxx::IsFloatOrComplex;
                                      };
    };

    template<typename Func>
    concept IsSpanInvocable =
        is_invocable_with_span<Func, float>::value ||
        is_invocable_with_span<Func, double>::value ||
        is_invocable_with_span<Func, long double>::value;


//    template< typename Func >
//    concept IsSpanInvocable = requires(Func f) {
//                                  requires nxx::IsFloatOrComplex< std::invoke_result_t< Func, std::span< float > > >;
//                                  requires nxx::IsFloatOrComplex< std::invoke_result_t< Func, std::span< double > > >;
//                                  requires nxx::IsFloatOrComplex< std::invoke_result_t< Func, std::span< long double > > >;
//                              };

    template< typename T >
    concept IsContainer = requires(T a) {
                              typename T::value_type;
                              typename T::iterator;
                              typename T::const_iterator;
                              {
                                  a.begin()
                              } -> std::same_as< typename T::iterator >;
                              {
                                  a.end()
                              } -> std::same_as< typename T::iterator >;
                              {
                                  a.size()
                              } -> std::convertible_to< std::size_t >;
                          };

    template< typename T >
    concept IsFloatContainer = IsContainer< T > && IsFloat< typename T::value_type >;

    template< typename S >
    requires(!IsContainer< S > && !std::is_pointer_v< S > && !std::is_array_v< S >)
    struct StructTraits
    {
    private:
        static constexpr auto commonType = [](S s) {
            auto [first, second] = s;
            using RT             = std::common_type_t< decltype(first), decltype(second) >;
            return RT { (first + second) };
        };

        static constexpr auto firstType = [](S s) {
            auto [first, second] = s;
            return first;
        };

        static constexpr auto secondType = [](S s) {
            auto [first, second] = s;
            return second;
        };

    public:
        using common_type = std::invoke_result_t< decltype(commonType), S >;
        using first_type  = std::invoke_result_t< decltype(firstType), S >;
        using second_type = std::invoke_result_t< decltype(secondType), S >;
    };

    template< typename S >
    using StructCommonType_t = typename StructTraits< S >::common_type;

    template< typename S >
    using StructFirstType_t = typename StructTraits< S >::first_type;

    template< typename S >
    using StructSecondType_t = typename StructTraits< S >::second_type;

    template< typename S >
    concept IsFloatStruct = !std::is_array_v< std::remove_cvref_t< S > > && nxx::IsFloat< typename StructTraits< S >::first_type > &&
                            nxx::IsFloat< typename StructTraits< S >::second_type >;

}    // namespace nxx

namespace nxx::poly
{
    template< typename T >
    requires nxx::IsFloat< T > || IsComplex< T >
    class Polynomial;

    template< typename POLY >
    concept IsPolynomial = std::same_as< POLY, Polynomial< typename POLY::value_type > >;
}    // namespace nxx::poly

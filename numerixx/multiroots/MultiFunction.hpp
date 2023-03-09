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

#ifndef NUMERIXX_MULTIFUNCTION_HPP
#define NUMERIXX_MULTIFUNCTION_HPP

#include <concepts>
#include <functional>

#include "../linalg/Matrix.hpp"

namespace nxx::multiroots
{
    namespace impl
    {
        template< typename... >
        struct FunctionTraits;

        template< typename RET, typename ARG >
        struct FunctionTraits< std::function< RET(ARG) > >
        {
            using ret_type = RET;
            using arg_type = ARG;
        };

    }    // namespace impl

    /**
     * @brief
     * @tparam FUNCARR
     */
    template< typename FUNCARR >
        requires std::floating_point< typename impl::FunctionTraits< typename FUNCARR::value_type >::ret_type >
    class MultiFunction
    {
    public:
        /*
         * Public alias declarations
         */
        using function_array = FUNCARR;
        using function_type  = typename FUNCARR::value_type;
        using return_type    = typename impl::FunctionTraits< function_type >::ret_type;
        using argument_type  = typename impl::FunctionTraits< function_type >::arg_type;

        /**
         * @brief
         * @param funcArray
         */
        explicit MultiFunction(const FUNCARR& funcArray) : m_functionArray(funcArray) {}

        /**
         * @brief
         * @tparam TFunc
         * @param funcArray
         */
        template< typename TFunc >
        MultiFunction(std::initializer_list< TFunc > funcArray) : m_functionArray(funcArray.begin(), funcArray.end())
        {}

        /**
         * @brief
         * @param arg
         * @return
         */
        argument_type operator()(const argument_type& arg) { return evaluate(arg); }

        /**
         * @brief
         * @param arg
         * @return
         */
        argument_type operator()(std::initializer_list< typename argument_type::value_type > arg) { return evaluate(arg); }

        /**
         * @brief
         * @return
         */
        auto size() const
        {
            return m_functionArray.size();
        }

        /**
         * @brief
         * @return
         */
        const auto& functionArray() const
        {
            return m_functionArray;
        }

    private:
        /**
         * @brief
         * @tparam TArg
         * @param arg
         * @return
         */
        template< typename TArg >
        argument_type evaluate(const TArg& arg)
        {
            using namespace nxx::linalg;

            if (arg.size() != m_functionArray.size())
                throw std::invalid_argument("MultiFunction Evaluation Error: number of arguments does not match number of equations.");

            argument_type args(arg.begin(), arg.end());

            std::vector< return_type > temp_result;
            temp_result.reserve(m_functionArray.size());
            for (auto& f : m_functionArray) temp_result.emplace_back(f(args));
            argument_type result(temp_result.begin(), temp_result.end());

            return result;
        }

        FUNCARR m_functionArray;
    };

    /*
     * MultiFunction deduction guide
     */
    template< typename TFunc >
    MultiFunction(std::initializer_list< TFunc > funcArray) -> MultiFunction< std::vector< TFunc > >;

    /**
     * @brief
     */
    class MultiDerivative
    {
    private:
    public:
    };

}    // namespace numerix::multiroots

#endif    // NUMERIXX_MULTIFUNCTION_HPP

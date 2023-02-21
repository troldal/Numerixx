//
// Created by Kenneth Balslev on 13/12/2022.
//

#ifndef NUMERIX_MULTIFUNCTION_HPP
#define NUMERIX_MULTIFUNCTION_HPP

#include <concepts>
#include <functional>

#include "../linalg/Matrix.hpp"

namespace numerix::multiroots
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
            using namespace numerix::linalg;

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

#endif    // NUMERIX_MULTIFUNCTION_HPP

//
// Created by Kenneth Balslev on 09/12/2022.
//

#ifndef NUMERIXX_JACOBIAN_HPP
#define NUMERIXX_JACOBIAN_HPP

#include <vector>
#include "Derivatives.hpp"

namespace nxx::deriv
{

    namespace impl
    {
        template< typename FN >
            requires std::invocable< FN, std::vector< std::invoke_result_t< FN, std::vector< typename FN::result_type > > > >
        auto getSingleVariableFunction(FN                                                                      multiVarFunc,
                                       const std::vector< std::invoke_result_t< FN, std::vector< typename FN::result_type > > >& coeffs,
                                       size_t                                                                  item)
        {
            using RT = std::invoke_result_t< FN, std::vector< typename FN::result_type > >;

            auto singleVarFunction = [=](RT value) {
                auto c  = coeffs;
                c[item] = value;
                return multiVarFunc(c);
            };

            return singleVarFunction;
        }

        template< typename FN >
            requires std::invocable< FN, std::vector< std::invoke_result_t< FN, std::vector< typename FN::result_type > > > >
        auto computePartialDerivs(FN multiVarFunc,
                                   const std::vector< std::invoke_result_t< FN, std::vector< typename FN::result_type > > >& coeffs)
        {

            using namespace nxx::deriv;
            std::vector<double> results;
            results.reserve(coeffs.size());

            size_t index = 0;
            for (auto item : coeffs) {
                auto fun = getSingleVariableFunction(multiVarFunc, coeffs, index);
                results.emplace_back(central(fun, coeffs[index]));
                ++index;
            }
            return results;
        }

    }    // namespace impl

    template< typename FNs, typename COEFF >
    auto computeJacobian(FNs functions, const COEFF& coeffs)
    {

        std::vector<typename COEFF::value_type> c(coeffs.begin(), coeffs.end());
        linalg::Matrix<typename COEFF::value_type> results(functions.size(), coeffs.size());
        size_t rowindex = 0;
        for(auto& fn : functions) {
            auto rowdata = impl::computePartialDerivs(fn, c);
            auto row = results.row(rowindex);

            size_t colindex = 0;
            for (auto elem : rowdata) {
                row(0, colindex) = elem;
                ++colindex;
            }
            ++rowindex;
        }

        return results;
    }



}    // namespace numerix::deriv

#endif    // NUMERIXX_JACOBIAN_HPP

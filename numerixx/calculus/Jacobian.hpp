/*
    o.     O O       o Oo      oO o.OOoOoo `OooOOo.  ooOoOOo o      O o      O
    Oo     o o       O O O    o o  O        o     `o    O     O    o   O    o
    O O    O O       o o  o  O  O  o        O      O    o      o  O     o  O
    O  o   o o       o O   Oo   O  ooOO     o     .O    O       oO       oO
    O   o  O o       O O        o  O        OOooOO'     o       Oo       Oo
    o    O O O       O o        O  o        o    o      O      o  o     o  o
    o     Oo `o     Oo o        O  O        O     O     O     O    O   O    O
    O     `o  `OoooO'O O        o ooOooOoO  O      o ooOOoOo O      o O      o

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

/*
    o.     O O       o Oo      oO o.OOoOoo `OooOOo.  ooOoOOo o      O o      O
    Oo     o o       O O O    o o  O        o     `o    O     O    o   O    o
    O O    O O       o o  o  O  O  o        O      O    o      o  O     o  O
    O  o   o o       o O   Oo   O  ooOO     o     .O    O       oO       oO
    O   o  O o       O O        o  O        OOooOO'     o       Oo       Oo
    o    O O O       O o        O  o        o    o      O      o  o     o  o
    o     Oo `o     Oo o        O  O        O     O     O     O    O   O    O
    O     `o  `OoooO'O O        o ooOooOoO  O      o ooOOoOo O      o O      o

    Copyright © 2023 Kenneth Troldal Balslev

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

#ifndef NUMERIXX_MULTIDERIVATIVES_HPP
#define NUMERIXX_MULTIDERIVATIVES_HPP

#include "Derivatives.hpp"

namespace nxx::deriv
{

    namespace impl
    {
        /*
         * Forward declaration of the VectorTraits class.
         */
        template< typename... >
        struct VectorTraits;

        /*
         * Specialization of the VectorTraits class for DynamicVector<T>
         */
        template< typename T >
        struct VectorTraits< blaze::DynamicVector< T > >
        {
            using value_type = T;
        };

        /*
         * Forward declaration of the MatrixTraits class.
         */
        template< typename... >
        struct MatrixTraits;

        /*
         * Specialization of the MatrixTraits class for DynamicMatrix<T>
         */
        template< typename T >
        struct MatrixTraits< blaze::DynamicMatrix< T > >
        {
            using value_type = T;
        };

    }    // namespace impl

    template< typename FN >
    concept IsMultiFunction = requires(FN fn) {
                                  // clang-format off
                                  { std::floating_point< std::invoke_result_t< FN, std::vector< double > > >};
                                  // clang-format on
                              };

    template< typename ARR >
    concept IsMultiFunctionArray = requires(ARR arr) {
                                       // clang-format off
                                       { arr.begin() } -> std::input_iterator;
                                       { arr.end() } -> std::input_iterator;
                                       { IsMultiFunction< typename ARR::value_type> };
                                       // clang-format on
                                   };

    template< typename T >
    using MultiReturnType = std::invoke_result_t< T, std::vector< double > >;

    template< typename ALGO >
    inline auto partialdiff(IsMultiFunction auto                  function,
                            auto                                  args,
                            size_t                                index,
                            MultiReturnType< decltype(function) > stepsize = StepSize< MultiReturnType< decltype(function) > >)
        requires std::convertible_to< typename impl::VectorTraits< decltype(args) >::value_type, MultiReturnType< decltype(function) > > ||
                 std::convertible_to< typename decltype(args)::value_type, MultiReturnType< decltype(function) > >
    {
        std::vector< MultiReturnType< decltype(function) > > argvector(args.begin(), args.end());

        // Create a lambda function for a multi variable function,
        // that takes a single argument and returns the function value
        auto f = [&](double value) {
            argvector[index] = value;
            return function(argvector);
        };

        // Return the partial derivative
        return diff< ALGO >(f, argvector[index], std::max(stepsize, stepsize * argvector[index]));
    }

    template< typename ALGO >
    inline auto multidiff(IsMultiFunction auto                  function,
                          auto                                  args,
                          MultiReturnType< decltype(function) > stepsize = StepSize< MultiReturnType< decltype(function) > >)
        requires std::convertible_to< typename impl::VectorTraits< decltype(args) >::value_type, MultiReturnType< decltype(function) > > ||
                 std::convertible_to< typename decltype(args)::value_type, MultiReturnType< decltype(function) > >
    {
        std::vector< MultiReturnType< decltype(function) > > argvector(args.begin(), args.end());

        size_t index = 0;
        for (auto arg : argvector) {
            args[index++] = *partialdiff< ALGO >(function, argvector, index, std::max(stepsize, stepsize * arg));
        }

        return args;
    }

    template< typename ALGO >
    inline auto jacobian(IsMultiFunctionArray auto                                   functions,
                         auto                                                        args,
                         MultiReturnType< typename decltype(functions)::value_type > stepsize = StepSize< decltype(stepsize) >)
        requires std::convertible_to< typename impl::VectorTraits< decltype(args) >::value_type, MultiReturnType< typename decltype(functions)::value_type > > ||
                 std::convertible_to< typename decltype(args)::value_type, MultiReturnType< typename decltype(functions)::value_type > >
    {
        using namespace nxx::deriv;
        using namespace blaze;

        using return_type = MultiReturnType< typename decltype(functions)::value_type >;

        blaze::DynamicMatrix< return_type > J(functions.size(), functions.size());

        size_t index = 0;
        for (auto& fn : functions) {
            auto pdiffs = multidiff< ALGO >(fn, args, stepsize);
            auto row    = blaze::row(J, index);
            auto dst    = row.begin();
            for (auto src = pdiffs.begin(); src != pdiffs.end(); ++src, ++dst) *dst = *src;
            ++index;
        }

        return J;
    }

}    // namespace nxx::deriv

#endif    // NUMERIXX_MULTIDERIVATIVES_HPP

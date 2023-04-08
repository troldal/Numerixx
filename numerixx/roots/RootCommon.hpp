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

#ifndef NUMERIXX_ROOTCOMMON_HPP
#define NUMERIXX_ROOTCOMMON_HPP

#include <stdexcept>

namespace nxx::roots
{

    // ========================================================================
    // ROOT ERROR/EXCEPTION CLASSES
    // ========================================================================

    class RootError : public std::runtime_error
    {
    public:
        explicit RootError(const char* msg) : std::runtime_error(msg) {};
    };

    // ========================================================================
    // TRAITS CLASSES FOR ROOT-FINDING WITHOUT DERIVATIVES
    // ========================================================================

    /*
     * Forward declaration of the Ridders class.
     */
    template< typename FN >
        requires std::floating_point< std::invoke_result_t< FN, double > >
    class Ridder;

    /*
     * Forward declaration of the Bisection class.
     */
    template< typename FN >
        requires std::floating_point< std::invoke_result_t< FN, double > >
    class Bisection;

    /*
     * Forward declaration of the RegulaFalsi class.
     */
    template< typename FN >
        requires std::floating_point< std::invoke_result_t< FN, double > >
    class RegulaFalsi;

    /*
     * Private implementation details.
     */
    namespace impl
    {
        /*
         * Forward declaration of the BracketingTraits class.
         */
        template< typename FN >
        struct BracketingTraits;

        /*
         * Specialization of the BracketingTraits class for Ridders<FN>
         */
        template< typename FN >
        struct BracketingTraits< Ridder< FN > >
        {
            using function_type = FN;
            using return_type   = std::invoke_result_t< FN, double >;
        };

        /*
         * Specialization of the BracketingTraits class for Bisection<FN>
         */
        template< typename FN >
        struct BracketingTraits< Bisection< FN > >
        {
            using function_type = FN;
            using return_type   = std::invoke_result_t< FN, double >;
        };

        /*
         * Specialization of the BracketingTraits class for RegulaFalsi<FN>
         */
        template< typename FN >
        struct BracketingTraits< RegulaFalsi< FN > >
        {
            using function_type = FN;
            using return_type   = std::invoke_result_t< FN, double >;
        };

    }    // namespace impl

    // ========================================================================
    // TRAITS CLASSES FOR ROOT-FINDING WITH DERIVATIVES
    // ========================================================================

    /*
     * Forward declaration of the DNewton class.
     */
    template< typename FN, typename DFN >
        requires std::floating_point< std::invoke_result_t< FN, double > > && std::floating_point< std::invoke_result_t< DFN, double > >
    class DNewton;

    /*
     * Forward declaration of the Newton class.
     */
    template< typename FN, typename DFN >
        requires std::floating_point< std::invoke_result_t< FN, double > > && std::floating_point< std::invoke_result_t< DFN, double > >
    class Newton;

    namespace impl
    {
        /*
         * Forward declaration of the PolishingTraits class.
         */
        template< typename... ARGS >
        struct PolishingTraits;

        /*
         * Specialization of the PolishingTraits class for Newton<FN, DFN>
         */
        template< typename FN, typename DFN >
        struct PolishingTraits< Newton< FN, DFN > >
        {
            using function_type        = FN;
            using deriv_type           = DFN;
            using function_return_type = std::invoke_result_t< FN, double >;
            using deriv_return_type    = std::invoke_result_t< DFN, double >;
        };

        /*
         * Specialization of the PolishingTraits class for DNewton<FN, DFN>
         */
        template< typename FN, typename DFN >
        struct PolishingTraits< DNewton< FN, DFN > >
        {
            using function_type        = FN;
            using deriv_type           = DFN;
            using function_return_type = std::invoke_result_t< FN, double >;
            using deriv_return_type    = std::invoke_result_t< DFN, double >;
        };

    }    // namespace impl

    template< typename FN >
        requires std::floating_point< std::invoke_result_t< FN, double > >
    class BracketSearchUp;

    template< typename FN >
        requires std::floating_point< std::invoke_result_t< FN, double > >
    class BracketSearchDown;

    template< typename FN >
    requires std::floating_point< std::invoke_result_t< FN, double > >
    class BracketExpandUp;

    template< typename FN >
    requires std::floating_point< std::invoke_result_t< FN, double > >
    class BracketExpandDown;

    template< typename FN >
    requires std::floating_point< std::invoke_result_t< FN, double > >
    class BracketExpandOut;

    template< typename FN >
    requires std::floating_point< std::invoke_result_t< FN, double > >
    class BracketSubdivide;

    /*
     * Private implementation details.
     */
    namespace impl
    {
        /*
         * Forward declaration of the BracketingTraits class.
         */
        template< typename FN >
        struct SearchingTraits;

        /*
         * Specialization of the BracketingTraits class for Ridders<FN>
         */
        template< typename FN >
        struct SearchingTraits< BracketSearchUp< FN > >
        {
            using function_type = FN;
            using return_type   = std::invoke_result_t< FN, double >;
        };

        template< typename FN >
        struct SearchingTraits< BracketSearchDown< FN > >
        {
            using function_type = FN;
            using return_type   = std::invoke_result_t< FN, double >;
        };

        template< typename FN >
        struct SearchingTraits< BracketExpandUp< FN > >
        {
            using function_type = FN;
            using return_type   = std::invoke_result_t< FN, double >;
        };

        template< typename FN >
        struct SearchingTraits< BracketExpandDown< FN > >
        {
            using function_type = FN;
            using return_type   = std::invoke_result_t< FN, double >;
        };

        template< typename FN >
        struct SearchingTraits< BracketExpandOut< FN > >
        {
            using function_type = FN;
            using return_type   = std::invoke_result_t< FN, double >;
        };

        template< typename FN >
        struct SearchingTraits< BracketSubdivide< FN > >
        {
            using function_type = FN;
            using return_type   = std::invoke_result_t< FN, double >;
        };

    }

}    // namespace nxx::roots

#endif    // NUMERIXX_ROOTCOMMON_HPP

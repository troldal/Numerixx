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


#ifndef NUMERIXX_POLYNOMIALERROR_HPP
#define NUMERIXX_POLYNOMIALERROR_HPP

#include <stdexcept>

namespace nxx::error
{
    class PolynomialError : public std::runtime_error
    {
    public:
        explicit PolynomialError(const char* msg) : std::runtime_error(msg) {};
    };

}    // namespace numerix::error

#endif    // NUMERIXX_POLYNOMIALERROR_HPP

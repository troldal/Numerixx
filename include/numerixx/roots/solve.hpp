// The one-call facade (DESIGN §6.13): a composition of the core, not a second implementation.
//
//     nxx::roots::solve(f, {lo, hi});   // the bracketing default (Brent, provisionally; phase 3 chooses by corpus counts)
//
// Deferred to phase 3: solve(f, x0) = then(expand from a window around x0, brent) and solve(f, df, x0) with rtsafe.
#pragma once

#include <numerixx/roots/bracket.hpp>
#include <numerixx/roots/brent.hpp>

NXX_BEGIN_HEADER

namespace nxx::roots
{
    template<class F, class In>
        requires detail::bracket_input_v<In>
    constexpr auto solve(const F& fn, const In& in)
    { return brent {}(fn, in); }

    template<class F, real T>
    constexpr auto solve(const F& fn, const T (&lo_hi)[2])
    {
        return brent {}(fn, std::pair<T, T> { lo_hi[0], lo_hi[1] });    // a pair, not the array: cl cannot order the overloads
    }
}    // namespace nxx::roots

NXX_END_HEADER

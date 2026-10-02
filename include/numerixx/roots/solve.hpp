// The one-call facade (DESIGN §6.13): a composition of the core, not a second implementation.
//
//     nxx::roots::solve(f, {lo, hi});   // the bracketing default (Brent, provisionally; phase 3 chooses by corpus counts)
//
// Deferred to phase 3: solve(f, x0) = then(expand from a window around x0, brent) and solve(f, df, x0) with rtsafe.
#pragma once

#include <numerixx/roots/bracket.hpp>
#include <numerixx/roots/brent.hpp>

#include <cstddef>
#include <type_traits>
#include <utility>

NXX_BEGIN_HEADER

namespace nxx::roots
{
    // Each overload is constrained on brent accepting the call, so that std::is_invocable_v is false for a function that
    // cannot take the bracket's scalar type; a deleted sibling states the reason, as in bracketing_facade (DESIGN §6.6).
    template<class F, class In>
        requires(detail::bracket_input_v<In> && std::is_invocable_v<const brent<>&, const F&, const In&>)
    constexpr auto solve(const F& fn, const In& in)
    { return brent {}(fn, in); }

    template<class F, class In>
        requires(detail::bracket_input_v<In> && !std::is_invocable_v<const brent<>&, const F&, const In&>)
    void solve(const F&, const In&) NXX_DELETE("the function cannot be called with the scalar type of the bracket");

    // A braced list or a C array: N is deduced, as in bracketing_facade, so {x} is not taken as {x, 0}.
    template<class F, real T, std::size_t N>
        requires(N == 2 && std::is_invocable_v<const brent<>&, const F&, std::pair<T, T>>)
    constexpr auto solve(const F& fn, const T (&lo_hi)[N])
    {
        return brent {}(fn, std::pair<T, T> { lo_hi[0], lo_hi[1] });    // a pair, not the array: cl cannot order the overloads
    }

    template<class F, real T, std::size_t N>
        requires(N == 2 && !std::is_invocable_v<const brent<>&, const F&, std::pair<T, T>>)
    void solve(const F&, const T (&)[N]) NXX_DELETE("the function cannot be called with the scalar type of the bracket");

    template<class F, class T, std::size_t N>
        requires(N != 2 || !real<T>)
    void solve(const F&, const T (&)[N]) NXX_DELETE("a bracket has two ends of a real type: write {lo, hi}, for example "
                                                    "{1.0, 2.0}");

    // Anything else (a scalar guess, a pair of integers, a pointer) is not a bracket; arrays are left to the overloads
    // above, as rejected_v does in bracketing_facade.
    template<class F, class In>
        requires(!detail::bracket_input_v<In> && !std::is_array_v<std::remove_cvref_t<In>>)
    void solve(const F&, const In&) NXX_DELETE("bracketing solvers need a bracket: pass {lo, hi}, nxx::bracket<T>::make(a, b), "
                                               "or a search result");
}    // namespace nxx::roots

NXX_END_HEADER

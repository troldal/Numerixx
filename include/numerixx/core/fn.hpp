// Function adaptors (DESIGN §6.12). The spike provides fn::counted, used to check that every solver's evaluation count
// equals the number of calls of the user's function; the other adaptors (negate, shift, extend_linearly, catching,
// value_of) arrive with the modules that need them.
#pragma once

#include <numerixx/core/callable.hpp>

#include <cstdint>
#include <functional>
#include <utility>

NXX_BEGIN_HEADER

namespace nxx::fn
{
    // Wraps fn and counts its calls in an external counter (stateful objects go in by reference, DESIGN §3.2).
    template<class F>
    class counted_fn
    {
        nxx::detail::copyable_box<F>          fn_;
        std::reference_wrapper<std::uint32_t> calls_;

    public:
        constexpr counted_fn(F fn, std::uint32_t& calls) : fn_(std::move(fn)), calls_(calls) {}

        template<class X>
        constexpr decltype(auto) operator()(const X& x) const
        {
            ++calls_.get();
            return std::invoke(*fn_, x);
        }

        constexpr std::uint32_t evaluation_cost() const noexcept { return nxx::cost_of(*fn_); }
    };

    template<class F>
    constexpr auto counted(F fn, std::uint32_t& calls)
    { return counted_fn<F> { std::move(fn), calls }; }
}    // namespace nxx::fn

NXX_END_HEADER
